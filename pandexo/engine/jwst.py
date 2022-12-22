import os
import sys
import json
import numpy as np
import pandas as pd
from copy import deepcopy 
from astropy.io import fits
from pandeia.engine.instrument_factory import InstrumentFactory
from pandeia.engine.perform_calculation import perform_calculation
from . import create_input as create
from .compute_noise import ExtractSpec
import astropy.units as u
import pickle
from pandeia.engine.calc_utils import build_default_calc, build_default_source

#constant parameters.. consider putting these into json file 
#max groups in integration
max_ngroup = 65536.0 
#minimum number of integrations
min_nint_trans = 1

#refdata directory
default_refdata_directory = os.environ.get("pandeia_refdata")

def compute_full_sim(dictinput,verbose=False): 
    """Top level function to set up exoplanet obs. for JW
    
    Function to set up explanet observations for JWST only and 
    compute simulated spectrum. It uses STScI's Pandeia to compute 
    instrument throughputs and WebbPSF to compute PSFs. 
    
    Parameters
    ----------
    dictinput : dict
        dictionary containing instrument parameters and exoplanet specific 
        parameters. {"pandeia_input":dict1, "pandexo_input":dict1}
    verbose : bool 
        (Optional) prints out check points throughout code 
    
    Returns
    -------
    dict
        large dictionary with 1d, 2d simualtions, timing info, instrument info, warnings
    
    Examples
    --------  
      
    >>> from .pandexo.engine.jwst import compute_full_sim 
    >>> from .pandexo.engine.justplotit import jwst_1d_spec
    >>> a = compute_full_sim({"pandeia_input": pandeiadict, "pandexo_input":exodict})
    >>> jwst_1d_spec(a)
    .. image:: 1d_spec.png
    
    Notes
    -----
    It is much easier to run simulations through either **run_online** or **justdoit**. **justdoit** contains functions to create input dictionaries and **run_online** contains web forms to create input dictionaries.
    
    See Also
    -------- 
        pandexo.engine.justdoit.run_pandexo : Best function for running pandexo runs
        pandexo.engine.run_online : Allows running functions through online interface
    """
    pandeia_input = dictinput['pandeia_input']
    pandexo_input = dictinput['pandexo_input']    	
	
    #define the calculation we'll be doing 
    if pandexo_input['planet']['w_unit'] == 'sec':
        calculation = 'phase_spec'
    else: 
        calculation = pandexo_input['calculation'].lower()

    #which instrument 
    instrument = pandeia_input['configuration']['instrument']['instrument']
    conf = pandeia_input['configuration']
    
    #if optimize is in the ngroups section, this will throw an error 
    #so create temp conf with 2 groups 
    if 'optimize' in str(conf['detector']['ngroup']): 
        conf_temp = deepcopy(conf) 
        conf_temp['detector']['ngroup'] = 2
    else: 
        conf_temp = conf


    i = InstrumentFactory(config=conf_temp)
    
    #detector parameters
    det_pars = i.read_detector_pars()
    fullwell = det_pars['fullwell']
    rn = det_pars['rn']
    mingroups = det_pars['mingroups']
        
    #exposure parameters 
    exp_pars = i.the_detector.exposure_spec
    tframe =exp_pars.tframe
    nframe = exp_pars.nframe
    nskip = exp_pars.nsample_skip

    sat_unit = pandexo_input['observation']['sat_unit']

    if sat_unit =='%':
        sat_level = pandexo_input['observation']['sat_level']/100.0*fullwell
    elif sat_unit =='e':
        sat_level = pandexo_input['observation']['sat_level']
    else: 
        raise Exception("Saturation Level Needs Units: % fullwell or Electrons ")
    

    
    #parameteres needed from exo_input
    mag = pandexo_input['star']['mag']
    
    
    noccultations = pandexo_input['observation']['noccultations']
    R = pandexo_input['observation']['R']

    noise_floor = pandexo_input['observation']['noise_floor']

    
    #get stellar spectrum and in transit spec
    star_spec = create.outTrans(pandexo_input['star'])
    #get rstar if user calling from grid 
    both_spec = create.bothTrans(star_spec, pandexo_input['planet'], star=pandexo_input['star'])
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])
    
    #get transit duration from phase curve or from input 
    if calculation == 'phase_spec': 
        transit_duration = max(both_spec['time']) - min(both_spec['time'])
    else: 
        #convert to seconds, then remove quantity and convert back to float 
        transit_duration = float((pandexo_input['planet']['transit_duration']*u.Unit(pandexo_input['planet']['td_unit'])).to(u.second)/u.second)

    #amount of exposure time out-of-occultation, as a fraction of in-occ time 
    try:
        expfact_out = pandexo_input['observation']['fraction'] 
        print("WARNING: key input fraction has been replaced with new 'baseline option'. See notebook example")
        pandexo_input['observation']['baseline'] = pandexo_input['observation']['fraction'] 
        pandexo_input['observation']['baseline_unit'] ='frac'
    except:
        if pandexo_input['observation']['baseline_unit'] =='frac':
            expfact_out = pandexo_input['observation']['baseline'] 
        elif pandexo_input['observation']['baseline_unit'] =='total':
            expfact_out = transit_duration/( pandexo_input['observation']['baseline'] - transit_duration)
        elif pandexo_input['observation']['baseline_unit'] =='total_hrs':
            expfact_out = transit_duration/( pandexo_input['observation']['baseline']*3600.0 - transit_duration)
        else: 
            raise Exception("Wrong units for baseine: either 'frac' or 'total' or 'total_hrs' accepted")

    #add to pandeia input 
    pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = out_spectrum
    
    if isinstance(pandeia_input["configuration"]["detector"]["ngroup"], (float,int)):
        m = {"ngroup":int(pandeia_input["configuration"]["detector"]["ngroup"]), "tframe":tframe,
            "nframe":nframe,"mingroups":mingroups,"nskip":nskip}
    else:
        #run pandeia once to determine max exposure time per int and get exposure params
        if verbose: print("Optimization Reqested: Computing Duty Cycle")
        m = {"maxexptime_per_int":compute_maxexptime_per_int(pandeia_input, sat_level) , 
            "tframe":tframe,"nframe":nframe,"mingroups":mingroups,"nskip":nskip}
        if verbose: print("Finished Duty Cycle Calc")

    #calculate all timing info
    timing, flags = compute_timing(m,transit_duration,expfact_out,noccultations)
    
    #Simulate out trans and in transit
    if verbose: print("Starting Out of Transit Simulation")
    out = perform_out(pandeia_input, pandexo_input,timing, both_spec)
    
    #extract extraction area before dict conversion
    extraction_area = out.extraction_area
    out = out.as_dict()
    out.pop('3d')
    if verbose: print("End out of Transit")

    #Remove effects of Quantum Yield from shot noise 
    out = remove_QY(out, instrument)

    #this kind of redundant going to compute inn from out instead 
    #keep perform_in but change inputs to (out, timing, both_spec)
    if verbose: print("Starting In Transit Simulation")
    inn = perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation)
    if verbose: print("End In Transit")
    


    #compute warning flags for timing info 
    warnings = add_warnings(out, timing, sat_level/fullwell, flags, instrument) 

    compNoise = ExtractSpec(inn, out, rn, extraction_area, timing)
    
    #slope method is pandeia's pure noise calculation (taken from SNR)
    #contains correlated noise, RN, dark current, sky, 
    #uses MULTIACCUM formula so we deviated from this. 
    #could eventually come back to this if Pandeia adopts First-Last formula
    if calculation == 'slope method': 
        #Extract relevant info from pandeia output (1d curves and wavelength) 
        #extracted flux in units of electron/s
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_slope_method()

    #derives noise from 2d postage stamps. Doing this results in a higher 
    #1d flux rate than the Pandeia gets from extracting its own. 
    #this should be used to benchmark Pandeia's 1d extraction  
    elif calculation == '2d extract':
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_2d_extract()
    
    #this is the noise calculation that PandExo uses online. It derives 
    #its own calculation of readnoise and does not use MULTIACUMM 
    #noise formula  
    elif calculation == 'fml':
        w = out['1d']['extracted_flux'][0]
        result = compNoise.run_f_minus_l()
    
    elif calculation == 'phase_spec':
        result = compNoise.run_phase_spec()
        w = result['time']
    else:
        result = None
        raise Exception('WARNING: Calculation method not found.')
        
    varin = result['var_in_1d']
    varout = result['var_out_1d']
    extracted_flux_out = result['photon_out_1d']
    extracted_flux_inn = result['photon_in_1d']

        
    #bin the data according to user input 
    if R != None: 
        wbin = bin_wave_to_R(w, R)

        photon_out_bin = uniform_tophat_sum(wbin, w,extracted_flux_out)
        photon_in_bin = uniform_tophat_sum(wbin,w, extracted_flux_inn)
        var_in_bin = uniform_tophat_sum(wbin, w,varin)
        var_out_bin = uniform_tophat_sum(wbin,w, varout)

        wbin = wbin[photon_out_bin > 0 ]
        photon_in_bin = photon_in_bin[photon_out_bin > 0 ]
        var_in_bin = var_in_bin[photon_out_bin > 0 ]
        var_out_bin = var_out_bin[photon_out_bin > 0 ]
        photon_out_bin = photon_out_bin[photon_out_bin>0]
    else: 
        wbin = w
        photon_out_bin = extracted_flux_out
        wbin = wbin[photon_out_bin>0]
        photon_in_bin = extracted_flux_inn
        photon_in_bin = photon_in_bin[photon_out_bin>0]
        var_in_bin = varin
        var_in_bin = var_in_bin[photon_out_bin>0]
        var_out_bin = varout
        var_out_bin = var_out_bin[photon_out_bin>0]
        photon_out_bin = photon_out_bin[photon_out_bin>0]
        
    
    if calculation == 'phase_spec':
        to = (timing["APT: Num Groups per Integration"]-1)*tframe
        ti = (timing["APT: Num Groups per Integration"]-1)*tframe
    else: 
        #otherwise error propagation and account for different 
        #times in transit and out 
        to = result['on_source_out']
        ti = result['on_source_in']
        
    var_tot = (to/ti/photon_out_bin)**2.0 * var_in_bin + (photon_in_bin*to/ti/photon_out_bin**2.0)**2.0 * var_out_bin
    error_spec = np.sqrt(var_tot)
        
    #factor in occultations to noise 
    nocc = timing['Number of Transits']
    error_spec = error_spec / np.sqrt(nocc)
        
    #Add in user specified noise floor 
    error_spec_nfloor = add_noise_floor(noise_floor, wbin, error_spec) 

    #add in random noise for the simulated spectrum 
    np.random.seed()
    rand_noise= error_spec_nfloor*(np.random.randn(len(wbin)))
    raw_spec = (photon_out_bin/to-photon_in_bin/ti)/(photon_out_bin/to)
    sim_spec = raw_spec + rand_noise 
    
    #if secondary tranist, multiply spectra by -1 
    if pandexo_input['planet']['f_unit'] == 'fp/f*':
        sim_spec = -1.0*sim_spec
        raw_spec = -1.0*raw_spec    
   
    #package processed data
    finalspec = {'wave':wbin,
              'spectrum': raw_spec,
              'spectrum_w_rand':sim_spec,
              'error_w_floor':error_spec_nfloor}
    
    rawstuff = {
                'electrons_out':photon_out_bin*nocc, 
                'electrons_in':photon_in_bin*nocc,
                'var_in':var_in_bin*nocc, 
                'var_out':var_out_bin*nocc,
                'e_rate_out':photon_out_bin/to,
                'e_rate_in':photon_in_bin/ti,
                'wave':wbin,
                'error_no_floor':error_spec, 
                'rn[out,in]':result['rn[out,in]'],
                'bkg[out,in]':result['bkg[out,in]']
                }
 
    result_dict = as_dict(out,both_spec ,finalspec, 
                timing, mag, sat_level, warnings,
                pandexo_input['planet']['f_unit'], rawstuff,calculation)

    return result_dict 
    
def compute_maxexptime_per_int(pandeia_input, sat_level):
    """Computes optimal maximum exposure time per integration
    
    Function to simulate 2d jwst image with 2 groups, 1 integration, 1 exposure 
    and return the maximum time 
    for one integration before saturation occurs. If saturation has 
    already occured, returns maxexptime_per_int as np.nan. This then 
    tells Pandexo to set min number of groups (ngroups =2). This avoids 
    error if saturation occurs. This routine assumes that min ngroups is 2. 
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia dictionary input 
    sat_level : int or float
        user defined saturation level in units of electrons
    
    Returns
    ------- 
    float 
        Maximum exposure time per integration before specified saturation level
    
    Examples
    --------
    
    >>> max_time = compute_maxexptime_per_int(pandeia_input, 50000.0)
    >>> print(max_time)
    12.0
    """
    
    #run once to get 2d rate image 
    pandeia_input['configuration']['detector']['ngroup'] = int(2 )
    pandeia_input['configuration']['detector']['nint'] = 1 
    pandeia_input['configuration']['detector']['nexp'] = 1
    
    report = perform_calculation(pandeia_input, dict_report=False)
    report_dict = report.as_dict() 
    
    # count rate on the detector in e-/second/pixel
    #det = report_dict['2d']['detector']
    det = report.signal.rate_plus_bg_list[0]['fp_pix']

    timeinfo = report_dict['information']['exposure_specification']
    #totaltime = timeinfo['tgroup']*timeinfo['ngroup']*timeinfo['nint']
    
    maxdetvalue = np.max(det)


    #maximum time before saturation per integration 
    #based on user specified saturation level

    try:
        maxexptime_per_int = sat_level/maxdetvalue
    except: 
        maxexptime_per_int = np.nan
    
    return maxexptime_per_int
        
def compute_timing(m,transit_duration,expfact_out,noccultations): 
    """Computes all timing info for observation
    
    Computes all JWST specific timing info for observation including. Some pertinent 
    JWST terminology:

        - frame: The result of sequentially clocking and digitizing all pixels in a rectangular area of an SCA. **Full-fame readout** means to digitize all pixels in an SCA, including reference pixels. **Frame** also applies to the result of clocking and digitizing a subarray on an SCA.
        - group: One or more consecutively read frames. There are no intervening resets. Frames may be averaged to form a group but for exoplanets the read out scheme is always 1 frame = 1 group
        - integration: The end result of resetting the detector and then non-destructively sampling it one or more times over a finite period of time before resetting the detector again. This is a unit of data for which signal is proportional to intensity, and it consists of one or more GROUPS.
        - exposure: The end result of one or more INTEGRATIONS over a finite period of time.  EXPOSURE defines the contents of a single FITS file.
    
    Parameters
    ---------
    m : dict 
        Dictionary output from **compute_maxexptime_per_int**
    transit_duration : float or int 
        transit duration in seconds 
    expfact_out : float or int 
        fraction of time spent in transit versus out of transit 
    noccultations : int 
        number of transits 
    
    Returns
    ------- 
    timing : dict
        All timing info
    warningflag : dict 
        Warning flags
    
    Examples
    --------
    >>> timing, flags = compute_timing(m, 2*60.0*60.0, 1.0, 1.0)
    >>> print((list(timing.keys())))
    ['Number of Transits', 'Num Integrations Out of Transit', 'Num Integrations In Transit', 
    'APT: Num Groups per Integration', 'Seconds per Frame', 'Observing Efficiency (%)', 'On Source Time(sec)', 
    'Exposure Time Per Integration (secs)', 'Reset time Plus 30 min TA time (hrs)',
    'APT: Num Integrations per Occultation', 'Transit Duration']
    """
    tframe = m['tframe']
    nframe = m['nframe']
    nskip = m['nskip']
    mingroups = m['mingroups']
    overhead_per_int = tframe #overhead time added per integration 
    try: 
        #are we starting with a exposure time ?
        maxexptime_per_int = m['maxexptime_per_int']
    except:
        #or a pre defined number of groups specified by user
        ngroups_per_int = m['ngroup']
        
    flag_default = "All good"
    flag_high = "All good"
    if 'maxexptime_per_int' in locals():
        #Frist, if maxexptime_per_int has been defined (from above), compute ngroups_per_int
        
        #number of frames in one integration is the maximum time beofre exposure 
        #divided by the time it takes for one frame. Note this does not include 
        #reset frames 

        nframes_per_int = np.floor(maxexptime_per_int/tframe)
    
        #for exoplanets nframe =1 an nskip always = 0 so ngroups_per_int 
        #and nframes_per_int area always the same 
        ngroups_per_int = np.floor(nframes_per_int/(nframe + nskip)) 
    
        #put restriction on number of groups 
        #there is a hard limit to the maximum number groups. 
        #if you exceed that limit, set it to the maximum value instead.
        #also set another check for saturation

        if ngroups_per_int > max_ngroup:
            ngroups_per_int = max_ngroup
            flag_high = "Groups/int > max num of allowed groups"
 
        if (ngroups_per_int < mingroups) | np.isnan(ngroups_per_int):
            ngroups_per_int = mingroups  
            nframes_per_int = mingroups
            flag_default = "NGROUPS<"+str(mingroups)+"SET TO NGROUPS="+str(mingroups)

    elif 'ngroups_per_int' in locals(): 
        #if it maxexptime_per_int been defined then set nframes per int 
        nframes_per_int = ngroups_per_int*(nframe+nskip)
        
        #if that didn't work its because maxexptime_per_int is nan .. run calc with mingroups
    else:
        #if maxexptime_per_int is nan then just ngroups and nframe to 2 
        #for the sake of not returning error
        ngroups_per_int = mingroups
        nframes_per_int = mingroups
        flag_default = "Something went wrong. SET TO NGROUPS="+str(mingroups)

          
    #the integration time is related to the number of groups and the time of each 
    #group 
    exptime_per_int = ngroups_per_int*tframe
    
    #clock time includes the reset frame 
    clocktime_per_int = (ngroups_per_int+1.0)*tframe
    
    #observing efficiency (i.e. what percentage of total time is spent on soure)
    eff = (ngroups_per_int - 1.0)/(ngroups_per_int + 1.0)
    
    #this says "per occultation" but this is just the in transit frames.. See below
    # transit duration / ((ngroups + reset)*frame time)
    nint_per_occultation =  transit_duration/((ngroups_per_int+1.0)*tframe)
    
    #figure out how many integrations are in transit and how many are out of transit 
    nint_in = np.ceil(nint_per_occultation)
    nint_out = np.ceil(nint_in/expfact_out)
    
    #you would never want a single integration in transit. 
    #here we assume that for very dim things, you would want at least 
    #3 integrations in transit 
    if nint_in < min_nint_trans:
        ngroups_per_int = np.floor(ngroups_per_int/min_nint_trans)
        exptime_per_int = (ngroups_per_int)*tframe
        clocktime_per_int = ngroups_per_int*tframe
        eff = (ngroups_per_int - 1.0)/(ngroups_per_int + 1.0)
        nint_per_occultation =  transit_duration/((ngroups_per_int+1.0)*tframe)
        nint_in = np.ceil(nint_per_occultation)
        nint_out = np.ceil(nint_in/expfact_out)
        
    if nint_out < min_nint_trans:
        nint_out = min_nint_trans
   
    timing = {
        "Transit Duration" : (transit_duration)/60.0/60.0,
        "Seconds per Frame" : tframe,
        "Time/Integration incl reset (sec)":clocktime_per_int,
        "APT: Num Groups per Integration" :int(ngroups_per_int), 
        "Num Integrations Out of Transit":int(nint_out),
        "Num Integrations In Transit":int(nint_in),
        "APT: Num Integrations per Occultation":int(nint_out+nint_in),
        "Observing Efficiency (%)": eff*100.0,
        "Transit+Baseline, no overhead (hrs)": (nint_out+nint_in)*clocktime_per_int/60.0/60.0, 
        "Number of Transits": noccultations
        }      
        
    return timing, {'flag_default':flag_default,'flag_high':flag_high}

def remove_QY(pandeia_dict, instrument):
    """Removes Quantum Yield from Pandeia Fluxes. Place Holder. 
    
    Parameters
    ----------
    pandeia_dict : dict 
        pandeia output dictionary
    instrument : str
        instrument running
    
    Returns
    -------
    dict 
        same exact dictionary with extracted_flux = extracted_flux/QY
    """
    if instrument == 'niriss':
        try:
            qy = fits.open(os.path.join(default_refdata_directory,'jwst', instrument,'qe' ,'jwst_niriss_h2rg_qe_20221003172003.fits'))
        except: 
            raise Exception('PANDEIA REFERENCE DATA NEEDS TO BE UPDATED')

        x_grid = pandeia_dict['1d']['extracted_flux'][0]
        qy_on_grid = np.interp(x_grid, qy[1].data['WAVELENGTH'], qy[1].data['CONVERSION'])
    elif instrument == 'nirspec':
        try:
            qy = fits.open(os.path.join(default_refdata_directory,'jwst', instrument,'qe' ,'jwst_nirspec_qe_20160902193401.fits'))
        except: 
            raise Exception('PANDEIA REFERENCE DATA NEEDS TO BE UPDATED')        
        x_grid = pandeia_dict['1d']['extracted_flux'][0]
        qy_on_grid = np.interp(x_grid, qy[1].data['WAVELENGTH'], qy[1].data['CONVERSION'])
    else:
        #nircam and miri currently have no qy effects
        qy_on_grid = 1.0

    pandeia_dict['1d']['extracted_flux'][1] = pandeia_dict['1d']['extracted_flux'][1]/qy_on_grid
    return pandeia_dict

def perform_out(pandeia_input, pandexo_input,timing, both_spec):
    """Runs pandeia for the out of transit data
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia specific input info 
    pandexo_input : dict 
        exoplanet specific observation info 
    timing : dict 
        timing dictionary from **compute_timing** 
    both_spec : dict 
        dictionary transit spectra computed from **createInput.bothTrans** 

    Returns
    -------
    dict 
        pandeia output dictionary for out of transit data 
    """
    #pandeia inputs, simulate one integration at a time 
    pandeia_input['configuration']['detector']['ngroup'] = int(timing['APT: Num Groups per Integration'])
    pandeia_input['configuration']['detector']['nint'] = int(timing['Num Integrations Out of Transit'])
    pandeia_input['configuration']['detector']['nexp'] = 1 

    report_out = perform_calculation(pandeia_input, dict_report=False)
    
    return report_out

    
def perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation): 
    """Computes in transit data 
    
    Runs Pandeia for the in transit data or computes the in transit simulation 
    from the out of transit pandeia run 
    
    Parameters
    ----------
    pandeia_input : dict 
        pandeia specific input info 
    pandexo_input : dict 
        exoplanet specific observation info 
    timing : dict 
        timing dictionary from **compute_timing** 
    both_spec : dict 
        dictionary transit spectra computed from **createInput.bothTrans** 
    out : dict 
        out of transit dictionary from **perform_in**
    calculation : str
        key which speficies the kind of noise calcualtion 
        (2d extract, slope method, fml, phase_spec). 
        Recommended for transit transmisstion spectra = fml

    Returns
    -------
    dict
        pandeia output dictionary 
    """
    
    #function to run pandeia for in transit
    if calculation == 'phase_spec':
        #return the phase curve since it's all we need 
        report_in = {'time': both_spec['time'],'planet_phase': both_spec['planet_phase']}
    elif calculation == 'fml':
        #for FML method, we only use the flux rate calculated in pandeia so 
        #can compute in transit flux rate without running pandeia a third time     
        report_in = deepcopy(out)
        
        transit_depth = np.interp(report_in['1d']['extracted_flux'][0],
                                    both_spec['wave'], both_spec['frac'])
        report_in['1d']['extracted_flux'][1] = report_in['1d']['extracted_flux'][1]*transit_depth
    else: 
        #only run pandeia a third time if doing slope method and need accurate run for the 
        #nint and timing
        pandeia_input['configuration']['detector']['ngroup'] = int(timing['APT: Num Groups per Integration'])
        pandeia_input['configuration']['detector']['nint'] = int(timing['Num Integrations In Transit'])
        pandeia_input['configuration']['detector']['nexp'] = 1
  
        in_transit_spec = np.array([both_spec['wave'], both_spec['flux_in_trans']])
    
        pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = in_transit_spec

        report_in = perform_calculation(pandeia_input, dict_report=True)
        instrument = pandeia_input['configuration']['instrument']['instrument']
        #remove QY effects 
        report_in = remove_QY(report_in, instrument)
        report_in.pop('3d')
    
    return report_in
          
def add_warnings(pand_dict, timing, sat_level, flags,instrument): 
    """Add warnings for front end 
    
    Adds in necessary warning flags for a JWST observation usually associated with 
    too few or too many groups or saturation. Alerts user if saturation level is higher 
    than 80 percent and if the number of groups is less than 5. Or, if the full well is 
    greater than 80. These warnings are currently very arbitrary. Will be updated as 
    better JWST recommendations are made. 
    
    Parameters
    ----------
    pand_dict : 
        output from pandeia run 
    timing : dict 
        output from **compute_timing** 
    sat_level : int or float 
        user specified saturation level in fractional (00/100)
    flags : dict 
        warning flags taken from output of **compute_timing**
    instrument : str 
        Only allowable strings are: "nirspec", "niriss", "nircam", "miri"
    
    Returns
    -------
    dict 
        all warnings 

    Notes
    -----    
    These are warnings are just suggestions and are not yet required.   
    """

    ngroups_per_int = timing['APT: Num Groups per Integration']
  
    #check for saturation 
    try:  
        flag_nonl = pand_dict['warnings']['partial_saturated']
    except: 
        flag_nonl = "All good"    
    try: 
        flag_sat = pand_dict['warnings']['full_saturated']
    except: 
        flag_sat = "All good"
        
    #check for too small number of groups
    flag_low = "All good"
    flag_perc = "All good"

    if (sat_level > .80) & (ngroups_per_int <3):
        flag_low = "% full well>80% & only " + str(ngroups_per_int) + " groups"
    if (sat_level > .80): 
        flag_perc = "% full well>80%"

     
    warnings = {
            "Group Number Too Low?" : flag_low,
            "Group Number Too High?": flags["flag_high"],
            "Non linear?" : flag_nonl,
            "Saturated?" : flag_sat,
            "% full well high?": flag_perc, 
            "Num Groups Reset?": flags["flag_default"]
    }

    return warnings     
    
def add_noise_floor(noise_floor, wave_bin, error_spec):
    """Add in noise floor 
    
    This adds in a user speficied noise floor. Does not add the noise floor in quadrature 
    isntead it sets error[error<noise_floor] = noise_floor. If a wavelength dependent 
    noise floor is given and the wavelength ranges are off, it interpolates the out of 
    range noise floor. 
    
    Parameters
    ----------
    noise_floor : str or int 
        file with two column [wavelength, noise(ppm)] or single number with constant noise floor in ppm 
    wave_bin : array of float
        final binned wavelength grid from simulation 
    error_spec : array of float 
        final computed error on the planet spectrum in units of rp^2/r*^2 or fp/f*
    
    Returns
    -------
    array of float 
        error_spec-- new error
    
    Examples
    --------
    
    >>> import numpy as np
    >>> wave = np.linspace(1,2.7,10)
    >>> error = np.zeros(10)+1e-6
    >>> newerror = add_noise_floor(20, wave, error)
    >>> print(newerror)
    [  2.00000000e-05   2.00000000e-05   2.00000000e-05   2.00000000e-05
       2.00000000e-05   2.00000000e-05   2.00000000e-05   2.00000000e-05
       2.00000000e-05   2.00000000e-05]
    """
    #add user specified noise floor 
    if (type(noise_floor)==float) | (type(noise_floor) == int):
        error_spec[error_spec<noise_floor*1e-6] = noise_floor*1e-6
    elif (type(noise_floor)==str):
        read_noise = np.genfromtxt(noise_floor, dtype=(float, float), names='w, n')
        w_overlap = (wave_bin>=min(read_noise['w'])) & (wave_bin<=max(read_noise['w'])) 
        wnoise = wave_bin[w_overlap]
        noise = np.zeros(len(wave_bin))
        noise[w_overlap] = np.interp(wnoise , read_noise['w'], read_noise['n'])
        noise[(wave_bin>max(read_noise['w']))] = read_noise['n'][read_noise['w'] == max(read_noise['w'])]
        noise[(wave_bin<min(read_noise['w']))] = read_noise['n'][read_noise['w'] == min(read_noise['w'])]
        error_spec[error_spec<noise*1e-6] = noise[error_spec<noise*1e-6]*1e-6
    else: 
        raise ValueError('Noise Floor added was not integer or file')
    return error_spec

def bin_wave_to_R(w, R):
    """Creates new wavelength axis at specified resolution
    
    Parameters
    ----------
    w : list of float or numpy array of float
        Wavelength axis to be rebinned 
    R : float or int 
        Resolution to bin axis to 
    
    Returns
    -------
    list of float
        New wavelength axis at specified resolution
    
    Examples
    --------
    
    >>> newwave = bin_wave_to_R(np.linspace(1,2,1000), 10)
    >>> print((len(newwave)))
    11
    """
    wave = []
    tracker = min(w)
    i = 1 
    ind= 0
    firsttime = True
    while(tracker<max(w)):
        if i <len(w)-1:
            dlambda = w[i]-w[ind]
            newR = w[i]/dlambda
            if (newR < R) & (firsttime):
                tracker = w[ind]
                wave += [tracker]
                ind += 1
                i += 1 
                firsttime = True
            elif newR < R:
                tracker = w[ind]+dlambda/2.0
                wave +=[tracker]
                ind = (np.abs(w-tracker)).argmin()
                i = ind+1
                firsttime = True
            else:
                firsttime = False            
                i+=1    
        else:
            tracker = max(w)
            wave += [tracker]
    return np.array(wave)
    
def uniform_tophat_sum(newx,x, y):
    """Adapted from Mike R. Line to rebin spectra
    
    Sums groups of points in certain wave bin 
    
    Parameters
    ----------
    newx : list of float or numpy array of float
        New wavelength grid to rebin to 
    x : list of float or numpy array of float 
        Old wavelength grid to get rid of 
    y : list of float or numpy array of float 
        New rebinned y axis 
    
    Returns
    -------
    array of floats 
        new wavelength grid 
        
    Examples 
    --------
        
    >>> from .pandexo.engine.jwst import uniform_tophat_sum
    >>> oldgrid = np.linspace(1,3,100)
    >>> y = np.zeros(100)+10.0
    >>> newy = uniform_tophat_sum(np.linspace(2,3,3), oldgrid, y)
    >>> newy
    array([ 240.,  250.,  130.])
    """
    newx = np.array(newx)
    szmod=newx.shape[0]
    delta=np.zeros(szmod)
    ynew=np.zeros(szmod)
    delta[0:-1]=newx[1:]-newx[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((x >= newx[i]-0.5*delta[i-1]) & (x < newx[i]+0.5*delta[i]))
        ynew[i]=np.sum(y[loc])
    loc=np.where((x > newx[0]-0.5*delta[0]) & (x < newx[0]+0.5*delta[0]))
    ynew[0]=np.sum(y[loc])
    return ynew

def target_acq(instrument, both_spec, warning): 
    """Contains functionality to compute optimal TA strategy 

    Takes pandexo normalized flux from create_input and checks for saturation, or 
    if SNR is below the minimum requirement for each. Then adds warnings and 2d displays 
    and target acq info to final output dict 

    Parameters
    ----------
    instrument : str 
        possible options are niriss, nirspec, miri and nircam 
    both_spec : dict
        output dictionary from **create_input** 
    warning : dict 
        output dictionary from **add_warnings** 

    Retruns
    -------

    """
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])

    #this automatically builds a default calculation 
    #I got reasonable answers for everything so all you should need to do here is swap out (instrument = 'niriss', 'nirspec','miri' or 'nircam')
    c = build_default_calc(telescope='jwst', instrument=instrument, mode='target_acq', method='taphot')
    c['scene'][0]['spectrum']['sed'] = {'sed_type':'input','spectrum':out_spectrum}
    c['scene'][0]['spectrum']['normalization']['type'] = 'none'
    rphot = perform_calculation(c, dict_report=True)

    #check warnings (pandeia doesn't return values for these warnings, so try will fail if all good)
    try: 
        warnings['TA Satruated?'] = rphot['warnings']['saturated']
    except:
        warnings['TA Satruated?'] = 'All good'

    try: 
        warnings['TA SNR Threshold'] = rphot['warnings']['ta_snr_threshold']
    except:
        warnings['TA SNR Threshold'] = 'All good'

    #build TA dict 
    ta = {'sn':rphot['scalar']['sn'],
            'ngroup': rphot['input']['configuration']['detector']['ngroup'], 
            'saturation':rphot['2d']['saturation']}

def as_dict(out, both_spec ,binned, timing, mag, sat_level, warnings, punit, unbinned,calculation): 
    """Format dictionary for output data 
    
    Takes all output from jwst run and converts it to simple dictionary 
    
    Parameters
    ----------
    out : dict 
        output dictionary from **compute_out**
    both_spec : dict 
        output dictionary from **createInput.bothTrans**
    binned : dict 
        dictionary from **wrapper** 
    timing : dict 
        dictionary from **compute_timing**
    mag : dict 
        magnitude of system 
    sat_level : float or int
        saturation level in electrons 
    warnings : dict 
        warning dictionary from **add_warnings**
    punit : "fp/f*" or "rp^2/r*^2"
        unit of supplied spectra options are: only options are fp/f* or rp^2/r*^2
    unbinned : dict 
        unbinned raw data from **wrapper**
    calculation : str 
        noise calculation type
    
    Returns
    -------
    dict 
        compressed dictionary 
    """
    #for emission spectrum

    p=1.0
    if punit == 'fp/f*': p = -1.0

    timing_div = pd.DataFrame.from_dict(timing, orient='index')
    timing_div.columns = ['Value']
    timing_div = timing_div.to_html()
    timing_div = '<table class="table table-striped"> \n' + timing_div[36:len(timing_div)] 
    timing_div = timing_div.encode()

    warnings_div = pd.DataFrame.from_dict(warnings, orient='index')
    warnings_div.columns = ['Value']
    warnings_div = warnings_div.to_html()
    warnings_div = '<table class="table table-striped"> \n' + warnings_div[36:len(warnings_div)]
    warnings_div = warnings_div.encode()
       
    input_dict = {
   	 "Target Mag": mag , 
   	 "Saturation Level (electons)": sat_level, 
   	 "Instrument": out['input']['configuration']['instrument']['instrument'], 
   	 "Mode": out['input']['configuration']['instrument']['mode'], 
   	 "Aperture": out['input']['configuration']['instrument']['aperture'], 
   	 "Disperser": out['input']['configuration']['instrument']['disperser'], 
   	 "Subarray": out['input']['configuration']['detector']['subarray'], 
   	 "Readmode": out['input']['configuration']['detector']['readmode'], 
 	 "Filter": out['input']['configuration']['instrument']['filter'],
 	 "Primary/Secondary": punit
    }
    
    input_div = pd.DataFrame.from_dict(input_dict, orient='index')
    input_div.columns = ['Value']
    input_div = input_div.to_html()
    input_div = '<table class="table table-striped"> \n' + input_div[36:len(input_div)]
    input_div = input_div.encode()
    
    #add calc type to input dict (doing it here so it doesn't output on webpage
    input_dict["Calculation Type"]= calculation
    
    final_dict = {
    'OriginalInput': {'model_spec':both_spec['model_spec'],
                     'model_wave' : both_spec['model_wave'],
                     'star_spec': both_spec['flux_out_trans']},
    'RawData': unbinned,
    'FinalSpectrum': binned,
        
    #pic output 
    'PandeiaOutTrans': out, 

    #all timing info 
    'timing': timing,
    'warning':warnings,
    'input':input_dict,
    
    #divs for html rendering    
    'timing_div':timing_div, 
    'input_div':input_div,
    'warnings_div':warnings_div,
    }
    return final_dict

    

