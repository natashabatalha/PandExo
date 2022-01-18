from .justdoit import run_pandexo,load_exo_dict,load_mode_dict
import numpy as np 
def test_nircam(): 
    exo_dict = load_exo_dict()
    exo_dict['observation']['sat_level'] = 100    #saturation level in percent of full well 
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits 
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                                 #without having to redo the calcualtion
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out 
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 

    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs
    exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit 

    inst = 'NIRCam F322W2'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 34
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 164
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'
    
    print(f'{inst} passed all tests')

    inst = 'NIRCam F444W'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 69
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 203
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'

    print(f'{inst} passed all tests')


def test_niriss(): 
    exo_dict = load_exo_dict()
    exo_dict['observation']['sat_level'] = 100    #saturation level in percent of full well 
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits 
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                                 #without having to redo the calcualtion
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out 
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 

    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs
    exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit 

    # order 1
    inst = 'NIRISS SOSS'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 3
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 73
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'
    
    print(f'{inst} Order 1 passed all tests')

    inst = 'NIRISS SOSS'
    inst_dict = load_mode_dict(inst)
    inst_dict['strategy']['order'] = 2
    inst_dict['configuration']['detector']['subarray'] = 'substrip256'
    result = run_pandexo(exo_dict,inst_dict, verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 2
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 158
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'

    print(f'{inst} Order 2 passed all tests')



def test_nirspec(): 
    exo_dict = load_exo_dict()
    exo_dict['observation']['sat_level'] = 100    #saturation level in percent of full well 
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits 
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                                 #without having to redo the calcualtion
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out 
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 

    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs
    exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit 

    inst = 'NIRSpec G395H'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 10
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 168
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'
    print(f'{inst} passed all tests')

    exo_dict['star']['mag'] = 12.0
    inst = 'NIRSpec Prism'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 5
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 76
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'

    print(f'{inst} passed all tests')  



def test_miri():
    exo_dict = load_exo_dict()
    exo_dict['observation']['sat_level'] = 100    #saturation level in percent of full well 
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits 
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                                 #without having to redo the calcualtion
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out 
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 

    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs
    exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit 
    
    inst = 'MIRI LRS'
    result = run_pandexo(exo_dict,[inst], verbose=False)
    ngroup = result['timing']['APT: Num Groups per Integration']
    last = 23
    assert np.isclose(ngroup,last,atol=1), f'{inst} failed ngroups test with {ngroup} vs. {last}'
    min_err = min(result['FinalSpectrum']['error_w_floor'])*1e6
    last = 53
    assert np.isclose(min_err,last,atol=5), f'{inst} surpassed 5 ppm difference with {min_err} vs. {last}'

    print(f'{inst} passed all tests')  

def run_all():
    try:
        test_miri()
    except: 
        print('MIRI Failed Tests')
    try:
        test_nircam()
    except: 
        print('NIRCam Failed Tests')
    try:
        test_niriss()
    except:
        print('NIRISS Failed Tests')
    try: 
        test_nirspec()
    except:
        print('NIRSpec Failed Tests')
