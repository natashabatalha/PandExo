import sys
import json
import copy
import numpy as np
import matplotlib.pyplot as plt
from pandeia.engine.perform_calculation import perform_calculation
import createInput as create
import matplotlib.pyplot as plt
import pandas as pd
import warnings 

#max groups in integration
max_ngroup = 65536.0 
#electron capacity full well 
fullwell = {
'niriss' : 85500.0, 
'nircam' : 90000.0,
'miri' : 250000.0, 
'nirspec' : 55000.0}

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

	#which instrument 
    instrument = pandeia_input['configuration']['instrument']['instrument']
    
    sat_level = pandexo_input['observation']['sat_level']/100.0*fullwell[instrument]

    #parameteres needed from exo_input
    mag = pandexo_input['star']['mag']
    
    transit_duration = pandexo_input['observation']['transit_duration']
    noccultations = pandexo_input['observation']['noccultations']
    wave_bin = pandexo_input['observation']['wave_bin']
    #amount of exposure time out-of-occultation, as a fraction of in-occ time 
    expfact_out = pandexo_input['observation']['fraction'] 
    noise_floor = pandexo_input['observation']['noise_floor']

    #do you want condensed input? or all the dictionary output parameters 
    online_true = pandexo_input['online']

    #get stellar spectrum and in transit spec
    star_spec = create.outTrans(pandexo_input['star'])
    both_spec = create.bothTrans(star_spec, pandexo_input['planet'])
    out_spectrum = np.array([both_spec['wave'], both_spec['flux_out_trans']])

    #add to pandeia input 
    pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = out_spectrum
    
    print "Starting compute maxexptime"    
    #run pandeia once to determine max exposure time per int and get exposure params
    m = compute_maxexptime_per_int(pandeia_input, sat_level) 
    print 'end compute maxexptime'
    #calculate all timing info
    timing, flags = compute_timing(m,transit_duration,expfact_out,noccultations)

    #Simulate out trans and in transit
    print "start out"
    out = perform_out(pandeia_input, pandexo_input,timing, both_spec)
    print "finish out"
    inn = perform_in(pandeia_input, pandexo_input,timing, both_spec)
    print "finish in"
    #compute warning flags for timing info 
    warnings = add_warnings(inn, timing, sat_level, flags, instrument) 
        
    #Extract relevant info from pandeia output (1d curves and wavelength) 
    w = out.curves['extracted_flux'][0]
    curves_out = out.curves
    curves_inn = inn.curves

    #In the following the SN is changed to incorporate number of occultations 
    #i.e. multiply by sqrt(n) 
    sn_in = curves_inn['sn'][1]*np.sqrt(noccultations)
    sn_out = curves_out['sn'][1]*np.sqrt(noccultations)
    
    extracted_flux_inn = curves_inn['extracted_flux'][1]
    extracted_noise_inn = curves_inn['extracted_flux'][1]/(sn_in)

    extracted_flux_out = curves_out['extracted_flux'][1]
    extracted_noise_out = curves_out['extracted_flux'][1]/(sn_out)
    
    varin = (extracted_noise_inn)**2.0
    varout = (extracted_noise_out)**2.0
    
    #bin the data according to user input 
    wbin, photon_out_bin = bin_data(w, extracted_flux_out, wave_bin)
    wbin, photon_in_bin = bin_data(w, extracted_flux_inn, wave_bin)
    wbin, var_in_bin = bin_data(w, varin, wave_bin)
    wbin, var_out_bin = bin_data(w, varout, wave_bin)
    
    #calculate total variance
    var_tot = var_in_bin + var_out_bin
    error = np.sqrt(var_tot)
    
    #add in random noise for the simulated spectrum 
    rand_noise= np.sqrt((var_in_bin+var_out_bin))*(np.random.randn(len(wbin)))
    raw_spec = (photon_out_bin-photon_in_bin)/photon_out_bin
    sim_spec = (photon_out_bin-photon_in_bin + rand_noise)/photon_out_bin 
    
    #if secondary tranist, multiply spectra by -1 
    if pandexo_input['planet']['f_unit'] == 'fp/f*':
        sim_spec = -1.0*sim_spec
        raw_spec = -1.0*raw_spec
    
    error_spec = error/photon_out_bin
   
    #Add in user specified noise floor 
    error_spec_nfloor = add_noise_floor(noise_floor, wbin, error_spec) 
   
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
    #condense output if online interface is running
    if online_true: 
        result_dict = as_dict(out,inn, both_spec ,binned, 
                timing, mag, sat_level, warnings,pandexo_input['planet']['f_unit'], unbinned)
    else: 
        result_dict = {'out': out, 'in' : inn, 'both_spec':both_spec, 
                        'binned':binned,'unbinned':unbinned ,'timing': timing}

    
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
    
    #maximum time before saturation per integration 
    #based on user specified saturation level
    try:
        maxexptime_per_int = sat_level/np.max(det)
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
    exptime_per_int = (ngroups_per_int-1.)*exptime_per_frame
    
    #clock time includes the reset frame 
    clocktime_per_int = ngroups_per_int*exptime_per_frame
    
    #observing efficiency (i.e. what percentage of total time is spent on soure)
    eff = exptime_per_int / (clocktime_per_int+overhead_per_int)
    
    #this says "per occultation" but this is just the in transit frames.. See below
    nframes_per_occultation = long(transit_duration/exptime_per_frame)
    ngroups_per_occultation = long(nframes_per_occultation/(nframe + nskip))
    nint_per_occultation = ngroups_per_occultation/(ngroups_per_int + 2)
    
    #figure out how many integrations are in transit and how many are out of transit 
    nint_in = long(nint_per_occultation)
    nint_out = long(nint_in/expfact_out)

    if nint_in == 0:
        nint_in = 1.0
        
    if nint_out == 0.0:
        nint_out = 1.0
   
    timing = {
        "Seconds per Frame" : exptime_per_frame,
        "Exposure Time Per Integration (secs)":exptime_per_int,
        "Num Groups per Integration" :ngroups_per_int, 
        "Num Integrations Out of Transit":nint_out,
        "Num Integrations In Transit":nint_in,
        "Num Integrations per Occultation":nint_out+nint_in,
        "On Source Time": noccultations*clocktime_per_int*(nint_out+nint_in)/60.0/60.0,
        "Reset time Plus TA time (hrs)": overhead_per_int*(nint_in + nint_out)/60.0/60.0 + 0.5,
        "Observing Efficiency (%)": eff*100.0,
        "Number of Transits": noccultations
        }      
        
    return timing, {'flag_default':flag_default,'flag_high':flag_high}

def perform_in(pandeia_input, pandexo_input,timing, both_spec): 
    """
    Runs Pandeia for the in transit data 

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
    #function to run pandeia for in transit
    pandeia_input['configuration']['detector']['ngroup'] = timing['Num Groups per Integration']
    pandeia_input['configuration']['detector']['nint'] = timing['Num Integrations In Transit']
    pandeia_input['configuration']['detector']['nexp'] = 1
  
    in_transit_spec = np.array([both_spec['wave'], both_spec['flux_in_trans']])
    
    pandeia_input['scene'][0]['spectrum']['sed']['spectrum'] = in_transit_spec

    report_in = perform_calculation(pandeia_input, dict_report=False)
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

    report_out = perform_calculation(pandeia_input, dict_report=False)

    return report_out
    
def add_warnings(inn, timing, sat_level, flags,instrument): 
    """
    Adds in necessary warning flags. 

	Parameters
	----------               
    Input:
        -inn: dictionary output from in transit data 
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
        flag_nonl = inn.as_dict()['warnings']['nonlinear']
    except: 
        flag_nonl = "All good"    
    try: 
        flag_sat = inn.as_dict()['warnings']['saturated']
    except: 
        flag_sat = "All good"
        
    #check for too small number of groups
    flag_low = "All good"
    flag_perc = "All good"

    if (sat_level > 60) & (ngroups_per_int <5):
        flag_low = "% full well>60% & only " + str(ngroups_per_int) + " groups"
    if (sat_level > 60): 
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

def as_dict(out, inn, both_spec ,binned, timing, mag, sat_level, warnings, punit, unbinned): 
    """
    computes full dictionary output for either output to user 
    or condenses the output for online pandexo interface
    """
    #for emission spectrum

    p=1.0
    if punit == 'fp/f*': p = -1.0
    
    out = out.as_dict()
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
    
    final_dict = {
    'OriginalInput': {'og_spec':p*(both_spec['flux_out_trans']-both_spec['flux_in_trans'])/both_spec['flux_out_trans'],
                     'og_wave' : both_spec['wave']},
    'RawData': unbinned,
    'FinalSpectrum': binned,
        
    #pic output 
    'PandeiaOutTrans': out, 
    #'inn':inn.as_dict(),

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
