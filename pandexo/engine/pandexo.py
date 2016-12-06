import sys
import json
import numpy as np
import pandas as pd
from copy import deepcopy 

from pandeia.engine.instrument_factory import InstrumentFactory
from pandeia.engine.perform_calculation import perform_calculation
import create_input as create
from compute_noise import ExtractSpec

#max groups in integration
max_ngroup = 65536.0 
#minimum number of integrations
min_nint_trans = 1
#electron capacity full well 


def wrapper(dictinput):
    """
	Function to set up exoplanet observation and compute simulated spectrum. 

	Parameters
	----------
	input: 
	    Dictionary containing pandeia output 
        Dictionary containing extra exoplanet spcific parameters
    
	output: 
	    Dictionary of Pandeia and Pandexo products: 
	        out: out of transit pandeia dict 
	        inn: in transit pandeia dict 
	        timing: timing info for exoplanet calc 
	        binning: output from binned up spectrum
	
	Attributes
	----------
	    compute_maxexpttime_per_int
	    compute_timing 
	    perform_in
	    perform_out
	    add_warnings 
	    bin_data
	    add_noise_floor
	    as_dict
	"""
    #constant parameters.. consider putting these into json file 
    pandeia_input = dictinput['pandeia_input']
    pandexo_input = dictinput['pandexo_input']    	
	
    #define the calculation we'll be doing 
    if pandexo_input['planet']['w_unit'] == 'sec':
        calculation = 'phase_spec'
    else: 
        calculation = pandexo_input['calculation'].lower()
    
    if calculation == 'scale':
        import HST_TExoNS as hst
        reload(hst)
        hmag            = pandexo_input['star']['hmag']
        trdur           = pandexo_input['observation']['transit_duration']
        numTr           = pandexo_input['observation']['noccultations']
        nchan           = pandexo_input['observation']['nchan']
        scanDirection   = pandexo_input['observation']['scanDirection']
        norbits         = pandexo_input['observation']['norbits']
        schedulability  = pandexo_input['observation']['schedulability']
        disperser       = pandeia_input['configuration']['instrument']['disperser'].lower()
        subarray        = pandeia_input['configuration']['detector']['subarray'].lower()
        nsamp           = pandeia_input['configuration']['detector']['nsamp']
        samp_seq        = pandeia_input['configuration']['detector']['samp_seq']
        deptherr        = hst.wfc3_TExoNS(hmag, trdur, numTr, nchan, disperser, scanDirection, subarray, nsamp, samp_seq, norbits)
        
        return deptherr
    
    #which instrument 
    instrument = pandeia_input['configuration']['instrument']['instrument']
    conf = {'instrument': pandeia_input['configuration']['instrument']}
    i = InstrumentFactory(config=conf)
    det_pars = i.get_detector_pars()
    fullwell = det_pars['fullwell']
    rn = det_pars['rn']
    pix_size = det_pars['pix_size']*1e-3 #convert from miliarcsec to arcsec
    sat_level = pandexo_input['observation']['sat_level']/100.0*fullwell

    #parameteres needed from exo_input
    mag = pandexo_input['star']['mag']
    
    
    noccultations = pandexo_input['observation']['noccultations']
    wave_bin = pandexo_input['observation']['wave_bin']
    #amount of exposure time out-of-occultation, as a fraction of in-occ time 
    expfact_out = pandexo_input['observation']['fraction'] 
    noise_floor = pandexo_input['observation']['noise_floor']


    #get stellar spectrum and in transit spec
    star_spec = create.outTrans(pandexo_input['star'])
    both_spec = create.bothTrans(star_spec, pandexo_input['planet'])
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])
        
    #get transit duration from phase curve or from input 
    if calculation == 'phase_spec': 
        transit_duration = max(both_spec['time']) - min(both_spec['time'])
    else: 
        transit_duration = pandexo_input['observation']['transit_duration']

    #add to pandeia input 
    pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = out_spectrum
    
    #run pandeia once to determine max exposure time per int and get exposure params
    print "Computing Duty Cycle"
    m = compute_maxexptime_per_int(pandeia_input, sat_level) 
    print "Finished Duty Cucle Calc"

    #calculate all timing info
    timing, flags = compute_timing(m,transit_duration,expfact_out,noccultations)
    
    #Simulate out trans and in transit
    print "Starting Out of Transit Simulation"
    out = perform_out(pandeia_input, pandexo_input,timing, both_spec)
    print "End out of Transit"

    #this kind of redundant going to compute inn from out instead 
    #keep perform_in but change inputs to (out, timing, both_spec)
    print "Starting In Transit Simulation"
    inn = perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation)
    print "End In Transit" 

    #compute warning flags for timing info 
    warnings = add_warnings(out, timing, sat_level, flags, instrument) 

    compNoise = ExtractSpec(inn, out, rn, pix_size, timing)
    
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
    wbin, photon_out_bin = bin_data(w, extracted_flux_out, wave_bin)
    wbin, photon_in_bin = bin_data(w, extracted_flux_inn, wave_bin)
    wbin, var_in_bin = bin_data(w, varin, wave_bin)
    wbin, var_out_bin = bin_data(w, varout, wave_bin)
    
    #calculate total variance
    var_tot = var_in_bin + var_out_bin
    error = np.sqrt(var_tot)
    
    #calculate error on spectrum
    error_spec = error/photon_out_bin
   
    #Add in user specified noise floor 
    error_spec_nfloor = add_noise_floor(noise_floor, wbin, error_spec) 

    
    #add in random noise for the simulated spectrum 
    rand_noise= np.sqrt((var_in_bin+var_out_bin))*(np.random.randn(len(wbin)))
    raw_spec = (photon_out_bin-photon_in_bin)/photon_out_bin
    sim_spec = (photon_out_bin-photon_in_bin + rand_noise)/photon_out_bin 
    
    #if secondary tranist, multiply spectra by -1 
    if pandexo_input['planet']['f_unit'] == 'fp/f*':
        sim_spec = -1.0*sim_spec
        raw_spec = -1.0*raw_spec
    
   
    #package binned data
    binned = {'wave':wbin,
              'spectrum': raw_spec,
              'spectrum_w_rand':sim_spec,
              'error_w_floor':error_spec_nfloor}
    
    unbinned = {
                'flux_out':extracted_flux_out, 
                'flux_in':extracted_flux_inn,
                'var_in':varin, 
                'var_out':varout, 
                'wave':w,
                'error_no_floor':np.sqrt(varin+varout)/extracted_flux_out
                }
 
    result_dict = as_dict(out,both_spec ,binned, 
                timing, mag, sat_level, warnings,
                pandexo_input['planet']['f_unit'], unbinned,calculation)

    return result_dict 



def compute_maxexptime_per_int(pandeia_input, sat_level):
    """ 
    Function to take in first 2d image and return the maximum time 
    for one integration before saturation occurs. If saturation has 
    already occured, returns maxexptime_per_int as np.nan. This then 
    tells Pandexo to set min number of groups (ngroups =2). This avoids 
    error if saturation occurs. 
    
    Parameters
	----------
    inputs:     
        -pandeia_input: Pandeia dictionary input
        -sat_level: user defined saturaton level
    
    output: 
        Dictionary output with 
            -maxexptime_per_int: time needed to reach the saturation level 
            -nframe: number of frames per group 
            -nskip: number of skipped frames 
            -exptime_per_frame: time per frame 
    """
    
    #run once to get 2d rate image 
    pandeia_input['configuration']['detector']['ngroup'] = 2 
    pandeia_input['configuration']['detector']['nint'] = 1 
    pandeia_input['configuration']['detector']['nexp'] = 1
    
    report = perform_calculation(pandeia_input, dict_report=False)
    report_dict = report.as_dict() 

            
    #check for hard saturation 
    if 'saturated' in report_dict['warnings'].keys(): 
        if report_dict['warnings']['saturated'][0:4] == 'Hard':
            print('Hard saturation with minimum number of groups')
    
    # count rate on the detector in e-/second 
    det = report_dict['2d']['detector']
    
    timeinfo = report_dict['information']['exposure_specification']
    #totaltime = timeinfo['tgroup']*timeinfo['ngroup']*timeinfo['nint']
    
    maxdetvalue = np.max(det)
    #maximum time before saturation per integration 
    #based on user specified saturation level
    try:
        maxexptime_per_int = sat_level/maxdetvalue
    except: 
        maxexptime_per_int = np.nan
        
    exptime_per_frame = report_dict['information']['exposure_specification']['tframe']
    nframe = report_dict['information']['exposure_specification']['nframe']
    nskip = report_dict['information']['exposure_specification']['nskip']
    return {'maxexptime_per_int':maxexptime_per_int, 'nframe':nframe, 'nskip':nskip, 'exptime_per_frame': exptime_per_frame}
    
def compute_timing(m,transit_duration,expfact_out,noccultations): 
    """
    Computes timing info for observation:  
	
	Parameters
	----------    
    Input: 
        - output dictionary from compute_maxexptime_per_int()
        - transit_duration 
        - expfact_out: fraction of time spent in versus out of transit
        - noccultations: number of transits
        
    Output: 
        -Dictionary of all timing info 
            frame time 
            exposure per integration 
            groups per integration
            integrations out of transit 
            integrations in transit 
            integrations per occultation 
            total observing hours 
            reset time plust TA time (set at 30 minutes)
            observing efficiency 
            number of transits  
            """
    exptime_per_frame = m['exptime_per_frame']
    nframe = m['nframe']
    nskip = m['nskip']
    overhead_per_int = exptime_per_frame #overhead time added per integration 
    maxexptime_per_int = m['maxexptime_per_int']

    flag_default = "All good"
    flag_high = "All good"
    try:
        #number of frames in one integration is the maximum time beofre exposure 
        #divided by the time it takes for one frame. Note this does not include 
        #reset frames 

        nframes_per_int = long(maxexptime_per_int/exptime_per_frame)
    
        #for exoplanets nframe =1 an nskip always = 1 so ngroups_per_int 
        #and nframes_per_int area always the same 
        ngroups_per_int = long(nframes_per_int/(nframe + nskip)) 
    
        #put restriction on number of groups 
        #there is a hard limit to the maximum number groups. 
        #if you exceed that limit, set it to the maximum value instead.
        #also set another check for saturation
    
        if ngroups_per_int > max_ngroup:
            ngroups_per_int = max_ngroup
            print("Num of groups per int exceeded max num of allowed groups"+str(ngroups_per_int))
            print("Setting number of groups to max value = 65536.0")
            flag_high = "Groups/int > max num of allowed groups"
 
        if ngroups_per_int < 2:
            ngroups_per_int = 2.0  
            nframes_per_int = 2
            print("Hard saturation during first group. Check Pandeia Warnings.")
            flag_default = "NGROUPS<2, SET TO NGROUPS=2 BY DEFAULT"
    except: 
        #if maxexptime_per_int is nan then just ngroups and nframe to 2 
        #for the sake of not returning error
        nframes_per_int = 2
        ngroups_per_int = 2
        flag_default = "NGROUPS<2, SET TO NGROUPS=2 BY DEFAULT"
                
    #the integration time is related to the number of groups and the time of each 
    #group 
    exptime_per_int = ngroups_per_int*exptime_per_frame
    
    #clock time includes the reset frame 
    clocktime_per_int = ngroups_per_int*exptime_per_frame
    
    #observing efficiency (i.e. what percentage of total time is spent on soure)
    eff = (ngroups_per_int - 1.0)/(ngroups_per_int + 1.0)
    
    #this says "per occultation" but this is just the in transit frames.. See below
    #nframes_per_occultation = long(transit_duration/exptime_per_frame)
    #ngroups_per_occultation = long(nframes_per_occultation/(nframe + nskip))
    nint_per_occultation =  transit_duration*eff/exptime_per_int
    
    #figure out how many integrations are in transit and how many are out of transit 
    nint_in = np.ceil(nint_per_occultation)
    nint_out = np.ceil(nint_in/expfact_out)
    
    #you would never want a single integration in transit. 
    #here we assume that for very dim things, you would want at least 
    #3 integrations in transit 
    if nint_in < min_nint_trans:
        ngroups_per_int = np.floor(ngroups_per_int/3.0)
        exptime_per_int = (ngroups_per_int-1.)*exptime_per_frame
        clocktime_per_int = ngroups_per_int*exptime_per_frame
        eff = (ngroups_per_int - 1.0)/(ngroups_per_int + 1.0)
        nint_per_occultation =  transit_duration*eff/exptime_per_int
        nint_in = np.ceil(nint_per_occultation)
        nint_out = np.ceil(nint_in/expfact_out)
        
    if nint_out < min_nint_trans:
        nint_out = min_nint_trans
   
    timing = {
        "Transit Duration" : transit_duration/60.0/60.0,
        "Seconds per Frame" : exptime_per_frame,
        "Exposure Time Per Integration (secs)":exptime_per_int,
        "Num Groups per Integration" :ngroups_per_int, 
        "Num Integrations Out of Transit":nint_out,
        "Num Integrations In Transit":nint_in,
        "Num Integrations per Occultation":nint_out+nint_in,
        "On Source Time(sec)": noccultations*clocktime_per_int*(nint_out+nint_in),
        "Reset time Plus 30 min TA time (hrs)": overhead_per_int*(nint_in + nint_out)/60.0/60.0 + 0.5,
        "Observing Efficiency (%)": eff*100.0,
        "Number of Transits": noccultations
        }      
        
    return timing, {'flag_default':flag_default,'flag_high':flag_high}

def perform_in(pandeia_input, pandexo_input,timing, both_spec, out, calculation): 
    """
    Runs Pandeia for the in transit data or computs in transit from out of transit 
    Pandeia run

	Parameters
	----------
    Input: 
        -pandeia dictionary input
        -pandexo dictionary input
        -timing info calculated from compute_timing
        -spectra calculated from createInput()
        
    Ouput: 
        -pandeia dictionary output 
    
    Attributes
    ---------- 
        - perform_calculation
        - bin_data_wave 
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
        pandeia_input['configuration']['detector']['ngroup'] = timing['Num Groups per Integration']
        pandeia_input['configuration']['detector']['nint'] = timing['Num Integrations In Transit']
        pandeia_input['configuration']['detector']['nexp'] = 1
  
        in_transit_spec = np.array([both_spec['wave'], both_spec['flux_in_trans']])
    
        pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = in_transit_spec

        report_in = perform_calculation(pandeia_input, dict_report=True)
        report_in.pop('3d')
    
    return report_in
          
	
def perform_out(pandeia_input, pandexo_input,timing, both_spec):
    """
    Runs Pandeia for the out of transit data 

	Parameters
	----------
    Input: 
        -pandeia dictionary input
        -pandexo dictionary input
        -timing info calculated from compute_timing
        -spectra calculated from createInput()
        
    Ouput: 
        -pandeia dictionary output 
    
    Attributes
    ---------- 
        - perform_calculation    
    """
    #pandeia inputs, simulate one integration at a time 
    pandeia_input['configuration']['detector']['ngroup'] = timing['Num Groups per Integration']
    pandeia_input['configuration']['detector']['nint'] = timing['Num Integrations Out of Transit']
    pandeia_input['configuration']['detector']['nexp'] = 1 

    report_out = perform_calculation(pandeia_input, dict_report=True)
    report_out.pop('3d')

    return report_out
    
def add_warnings(pand_dict, timing, sat_level, flags,instrument): 
    """
    Adds in necessary warning flags. 

	Parameters
	----------               
    Input:
        -pand_dict: dictionary output from pandeia run 
        -timing: timing calculated from "compute_timing()"
        -sat_level: user specified saturaiton level 
        -flags: two flags which are calculated in "compute_timing()"
        -instrumnet: what instrument defines # electrons permitted in well
    Output: 
        -Dictionary of warning values:
            Group Number Too Low? : groups calculated are less than 5
            Group Number Too High?: Group calculated are more than 65535 
            Non Linear?: Taken from Pandeia output for "soft saturation" 
            Saturated?: Taken from Pandeia output for "hard saturation"
            % full well high? : user inputed value above 60 for percent full well 
            Num groups reset?: User put in parameters which resulted in ngroups <2 so 
                                pandexo set ngroups=2
        
    """
    ngroups_per_int = timing['Num Groups per Integration']
  
    #check for saturation 
    try:  
        flag_nonl = pand_dict['warnings']['nonlinear']
    except: 
        flag_nonl = "All good"    
    try: 
        flag_sat = pand_dict['warnings']['saturated']
    except: 
        flag_sat = "All good"
        
    #check for too small number of groups
    flag_low = "All good"
    flag_perc = "All good"

    if (sat_level > 80) & (ngroups_per_int <5):
        flag_low = "% full well>60% & only " + str(ngroups_per_int) + " groups"
    if (sat_level > 80): 
        flag_perc = "% full well>60%"

     
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
    """
    Adds in User specified noise floor
    
    Does not add noise floor in quadrature. Instead: 
    error_spec(error_spec<noise_floor) = noise_floor
    
    If wavelength dependent noisefloor given, interpolation is done 

	Parameters
	---------- 
    Inputs: 
        -noise_floor: file or single number 
        -wave_bin: wavelength axis 
        -error_spec: noise curve calculated from wrapper
    Outputs: 
        -error_spec: new noise curve with noise floor added
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
        
def bin_data_wave(wlgrid,old_wave, spec):
    """
    Instead of binning to a specific wave bin, this bins to a specific 
    wavelength. 
    Used for computing Pandeia In transit dict wihtout running pandeia a third time 
    
    Parameters:
    -----------
        Inputs:
            wlgrid: wave grid you want to bin to
            old_wave: old wavelength you want to change 
            spec: spectra that needs to be binned to wlgrid
        Outputs: 
            Fint: new spectra on new wlgrid
    """
    wno = 1e4/old_wave
    Fp = spec[::-1]
    szmod=wlgrid.shape[0]

    delta=np.zeros(szmod)
    Fint=np.zeros(szmod)
    delta[0:-1]=wlgrid[1:]-wlgrid[:-1]  
    delta[szmod-1]=delta[szmod-2] 
    #pdb.set_trace()
    for i in range(szmod-1):
        i=i+1
        loc=np.where((1E4/wno >= wlgrid[i]-0.5*delta[i-1]) & (1E4/wno < wlgrid[i]+0.5*delta[i]))
        Fint[i]=np.mean(Fp[loc])
        print wlgrid[i]-0.5*delta[i-1], wlgrid[i]+0.5*delta[i-1], np.mean(Fp[loc])
    loc=np.where((1E4/wno > wlgrid[0]-0.5*delta[0]) & (1E4/wno < wlgrid[0]+0.5*delta[0]))
    Fint[0]=np.mean(Fp[loc])
    
    return Fint

def bin_data(x,y,wlength):
    """ 
    Takes 2 arrays x and y and bins them into groups of blength.
    
    Parameters
	----------
        Inputs:     
            -x, y:                   1D lists or numpy arrays
        Outputs:    
            - xout, yout, yerrout,noise:    1D numpy arrays
    """

    ii = 0
    start = 0
    ind = []
    for i in range(0, len(x)-1):
        if x[i+1] - x[start] >= wlength:
            ind.append(i+1)
            start = i 
    
    if ind[len(ind)-1] != (len(x)):
        ind.append(len(x))
    
    # convert to arrays if necessary
    x = np.array(x)
    y = np.array(y)


    xout,yout= [],[]
    first = 0
    for i in ind:
        xout.append(sum(x[first:i])/len(x[first:i]))
        yout.append(sum(y[first:i]))
        first = i 
    return np.array(xout),np.array(yout)

def as_dict(out, both_spec ,binned, timing, mag, sat_level, warnings, punit, unbinned,calculation): 
    """
    computes full dictionary output for either output to user 
    or condenses the output for online pandexo interface
    """
    #for emission spectrum

    p=1.0
    if punit == 'fp/f*': p = -1.0
    
    timing_div = pd.DataFrame(timing.items(), columns=['Timing Info', 'Values']).to_html().encode()
    timing_div = '<table class="table table-striped"> \n' + timing_div[36:len(timing_div)]
    
    warnings_div = pd.DataFrame(warnings.items(), columns=['Check', 'Status']).to_html().encode()
    warnings_div = '<table class="table table-striped"> \n' + warnings_div[36:len(warnings_div)]
       
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
    
    input_div = pd.DataFrame(input_dict.items(), columns=['Component', 'Values']).to_html().encode()
    input_div = '<table class="table table-striped"> \n' + input_div[36:len(input_div)]
    
    #add calc type to input dict (doing it here so it doesn't output on webpage
    input_dict["Calculation Type"]= calculation
    
    final_dict = {
    'OriginalInput': {'og_spec':both_spec['og_spec'],
                     'og_wave' : both_spec['og_wave']},
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
