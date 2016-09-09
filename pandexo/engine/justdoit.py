import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory
from pandexo import wrapper
import os 
import pickle as pkl
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

ALL = {"MIRI LRS":False,
       "NIRISS SOSS":False,
       "NIRSpec G140M":False,
       "NIRSpec G140H":False,
       "NIRSpec G235M":False,
       "NIRSpec G235H":False,
       "NIRSpec G395M":False,
       "NIRSpec G395H":False,
       "NIRSpec Prism":False,
       "NIRCam F322W2":False,
       "NIRCam F444W":False}


def print_instruments():
    print "Choose from the following:"
    print ALL.keys()
    return

def load_exo_dict():
    pandexo_input = {
        "star":{
                "type" : "user or phoenix", 
                "starpath" : "file path",
                "w_unit": "Angs,cm,um,cm or Hz",
                "f_unit": "W/m2/um, FLAM, Jy, or erg/s/cm2/Hz",
                "mag": "magnitude",
                "ref_wave": "corresponding ref wave",
                "temp": "Only if phoenix, (in K)",
                "metal": "in log Fe/H",
                "logg": "cgs"
            },

    "planet" :{
	        "type": "user",
            "exopath" : "file path",
            "w_unit": "Angs,cm,um,cm or Hz",
            "f_unit": "rp/r* or fp/f*"
            },

    "observation": {
        "sat_level": "in % sat level",
        "transit_duration": "in seconds",
        "noccultations": "num transits",
        "wave_bin": "in micron",
        "fraction": "time in/out", 
        "noise_floor":"constant number or file name"
        }
    }
    print "Replace all inputs before feeding to run_modes:"
    return pandexo_input   
    
def load_mode_dict(inst):
    
    return

def get_thruput(instrument):
    obsmode = {
               'instrument': instrument,
               'mode': mode,
               'filter': config['filter'],
               'aperture': config['aperture'],
               'disperser': config['disperser']
               }
                             
    conf = {'instrument': obsmode}

    i = InstrumentFactory(config=conf)
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
    pce = i.get_total_eff(wave)

    return wave,pce


def run_pandexo(exo=None, inst=None, param_space = None, param_range = None,
                                            output_path=os.getcwd()):
        
    #single instrument mode with dictionary input OR single planet 
    if type(inst) == dict: 
        print "Running Single Case w/ User Instrument Dict"
        result =wrapper({"pandeia_input": inst , "pandexo_input":exo})
        return result 
 
    #single instrument mode and single planet OR several planets   
    if (type(inst) == str) or ((type(inst)==list or type(inst)==np.ndarray) and len(inst)==1): 
        if param_space==None or param_range==None:
            print "Running Single Case for: " + inst
            inst_dict = load_mode_dict(inst)
            result =wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})
            return result
         
        print "Running through parameters in parallel: " + param_space 
        def run_param_space(i): 
            #break up parameter space to two separate dictionary keys
            key1 = param_space[0:param_space.find('+')]
            key2 = param_space[param_space.find('+')+1:len(param_space)]
            exo[key1][key2] = i 
            #load in correct dict format
            inst_dict = load_mode_dict(inst)
            name = os.path.split(i)[1]
            return {name: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})}
        
        #run the above function in parallel 
        results = Parallel(n_jobs=num_cores)(delayed(run_param_space)(i) for i in param_range)
        
        #dump all results [an array of dictionaries] into single file
        filename = param_space + '.p'
        pkl.dump(results, open(os.path.join(output_path,filename),'w'))
        return results
        
    #run several different instrument modes and single planet
    if type(inst) == list or type(inst) == np.ndarray:
        return  
        
    