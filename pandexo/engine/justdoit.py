import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory
from pandexo import wrapper
import load_modes as lm
import os 
import pickle as pkl
from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

ALL = {"WFC3 G141":False,
       "MIRI LRS":False,
       "NIRISS SOSS_Or1":False,
       "NIRISS SOSS_Or2":False,
       "NIRSpec G140M":False,
       "NIRSpec G140H":False,
       "NIRSpec G235M":False,
       "NIRSpec G235H":False,
       "NIRSpec G395M":False,
       "NIRSpec G395H":False,
       "NIRSpec Prism":False,
       "NIRCam F322W":False,
       "NIRCam F444W":False}


def print_instruments():
    print "Choose from the following:"
    print ALL.keys()
    return

def load_exo_dict():
    """
    Function loads in empty exoplanet dictionary for pandexo input
    """
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
    '''
    Small function to pull in correct instrument dictionary 
    See load_modes.py 
    '''
    return lm.SetDefaultModes(inst).pick()

def get_thruput(inst):
    """
    Pulls complete instrument wave length solution and photon to 
    electron conversion efficiency (PCE) based on instrument key input 
    
    input: 
        one of the instrument keys above (ALL) 
    attributes: 
        calls load_mode_dict
    returns dictionary with wave and pce
    """
    
    #pull correct dictionary
    input_dict = lm.SetDefaultModes(inst).pick()
                             
    conf = {'instrument': input_dict['configuration']['instrument']}
    i = InstrumentFactory(config=conf)
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
    pce = i.get_total_eff(wave)
    return {'wave':wave,'pce':pce}

def run_param_space(i,exo,inst,param_space): 
    """
    Function to run through Parallel for running multiple 
    planet types with a single isntrument 
    """
    #break up parameter space to two separate dictionary keys
    key1 = param_space[0:param_space.find('+')]
    key2 = param_space[param_space.find('+')+1:len(param_space)]
    exo[key1][key2] = i 
    #load in correct dict format
    inst_dict = load_mode_dict(inst)
    name = os.path.split(str(i))[1]
    return {name: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})}

def run_inst_space(inst,exo): 
    """
    Function to run through Parallel for running multiple 
    instrument modes with a single planet 
    """
    #load in correct dict format
    inst_dict = load_mode_dict(inst)
    return {inst: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})}


def run_pandexo(exo, inst, param_space = 0, param_range = 0,save_file = True,
                            output_path=os.getcwd(), output_file = ''):  
    """
    Functionality: program contains functionality for running single or 
    multiple runs of PandExo 
    
    Inputs: 
    exo: exoplanet input dictionary 
    inst: instrument input dictionary OR LIST of keys (for allowable keys 
            see print_instruments()
    param_space: Default is 0, no exoplanet parameter space, to run through a parameter 
                need to specify which one need to specify two keys
                from exo dict with + in between. 
                i.e. observation+fraction
                    star+temp
                    planet+exopath
    param_range: supply array or list over which to run the parameters space.. 
                    i.e. array of temperatures if running through stellar temp or 
                    array of files if running through planet models
    save_file: default saves file output but you can spefiy no file output
    output_path: path for output files. default will print out files in current working directory 
    output_file: file name, default = single.p for single run 
                            param_space.p for parameter run 
                            instrument_run.p for instrument parameter space 
    
    Outputs: 
    result: For single run output will just be a single PandExo output dictionary 
            https://github.com/natashabatalha/PandExo/wiki/PandExo-Output
            For multiple runs the output will be organized into a list with each 
            a dictionary named by whatever you are looping through 
            i.e. [{'First temp': PandExoDict}, {'Second temp': PandExoDict}, etc..]
            
    Attributes: 
    load_mode_dict: to call instrument dictionaries based off keys 
    run_inst_space: to run instrument parameters space in parallel 
    run_param_space: to run exoplanet parameter space in parallel
    """

    #single instrument mode with dictionary input OR single planet 
    if type(inst) == dict: 
        print "Running Single Case w/ User Instrument Dict"
        results =wrapper({"pandeia_input": inst , "pandexo_input":exo})
        if output_file == '':
            output_file = 'singlerun.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'w'))
        return results

    #make sure inst is in list format.. makes my life so much easier
    try:
        if type(inst) != list: 
            raise ValueError
    except ValueError:
        print 'Instrument input is not dict so must be list'
        print 'Enter in format [bla] or [bla,bla]' 
        return    
         
    #single instrument mode and single planet OR several planets  
     
    if len(inst)==1 and inst[0] != 'RUN ALL': 
        
        #start case of no parameter space run 
        if param_space==0 or param_range==0:
            print "Running Single Case for: " + inst[0]
            inst_dict = load_mode_dict(inst[0])
            results =wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})
            if output_file == '':
                output_file = 'singlerun.p'
            if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'w'))
            return results
         
        #if there are parameters to cycle through this will run
        print "Running through exo parameters in parallel: " + param_space 
        #run the above function in parallel 
        results = Parallel(n_jobs=num_cores)(delayed(run_param_space)(i,exo,inst[0],param_space) for i in param_range)
        
        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file = param_space + '.p'

        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'w'))
        return results
        
    #run several different instrument modes and single planet
    print "Running select instruments" 
    if len(inst)>1:
        
        results = Parallel(n_jobs=num_cores)(delayed(run_inst_space)(i, exo) for i in inst)

        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file =  'instrument_run.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'w'))
        return results
            
    #cycle through all options  
    elif inst[0].lower() == 'run all':
        print "Running through all instruments"  
        results = Parallel(n_jobs=num_cores)(delayed(run_inst_space)(i, exo) for i in ALL.keys())
    
        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file =  'instrument_run.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'w'))
        return results
