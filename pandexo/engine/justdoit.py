import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory
from pandexo import wrapper
import load_modes as lm
import os 
import pickle as pkl
from joblib import Parallel, delayed
import multiprocessing
import json 

num_cores = multiprocessing.cpu_count()

ALL = {"WFC3 G141":False,
       "MIRI LRS":False,
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
    """Prints a list of the possible instrument templates to load
    """
    print "Choose from the following:"
    print ALL.keys()
    return

def load_exo_dict():
    """Loads in empty exoplanet dictionary for pandexo input
    
    Loads in empty exoplanet dictionary so that the user can manually edit different planet 
    parameters. Must be done for every bash run. User must edit each keys within the dictionary. 
    
    Returns
    -------
    dict 
        Empty dictionary to be filled out by the user before running PandExo
    
    Example
    -------
    >>> exo_dict = load_exo_dict()
    >>> exo_dict['planet']['transit_duration'] = 2*60*60 #2 hours
    """
    with open(os.path.join(os.path.dirname(__file__), "reference",
                               "exo_input.json")) as data_file:
        pandexo_input = json.load(data_file)
    return pandexo_input
    
def load_mode_dict(inst):
    """Function to pull in correct instrument dictionary 
    
    This is the instrument counterpart to load_exo_dict. It loads in a template instrument 
    dictionary for a specific instrument mode. 
    
    Parameters
    ----------
    inst : str
        One of the allowable instrument keys. To see which instrument modes are available 
        use `print_instruments()`
    
    Returns
    -------
    dict 
        Filled out template of instrument dictionary, which can be editted before 
        running to PandExo (not required). 
    
    Example
    -------
    >>> inst_dict = load_mode_dict('MIRI LRS')
    >>> inst_dict['configuration']['instrument']['aperture'] = 'lrsslit'
    """
    return lm.SetDefaultModes(inst).pick()

def get_thruput(inst):
    """Returns complete instrument photon to electron conversion efficiency
    Pulls complete instrument photon to electron conversion efficiency 
    (PCE) based on instrument key input 
    
    Parameters
    ----------
    inst : str
        One of the instrument keys in `print_instruments`
    
    Returns
    ------- 
    dict 
        Dictionary with wave solution and PCE
    
    Example
    -------
    >>> thru_dict = get_thruput('NIRISS SOSS_Or1')
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
    """Changes exo dictionary and submits run  
    
    This function is used to reset the exo dictionary based on what parameter
    is being looped over and submits run to `wrapper` so that all the jobs 
    are run in parallel 
    
    Parameters
    ----------
    i : str or float
        Can be either a str or float based on what you are looping through (str for 
        filenames, float for stellar temps, float for magnitudes, etc)
    exo : dict
        Exoplanet dictionary which can be loaded in and editted through `load_exo_dict`
    inst : str 
        Key which indicates with instrument 
    param_space : str 
        Set of keys within exo_dict to indicate which parameter to loop through. 
        Should be in the format of "first level of dict"+"second level of dict". 
        For example, for stellar temp `param_space` would be "star+temp"
    
    Returns
    -------
    dict
        Dictionary with output of pandexo. Key is the value of the parameter that was 
        looped through.  
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
    """Changes inst dictionary and submits run  
    
    This function is used to reset the instrument dictionary. 
    
    Parameters
    ----------
    exo : dict
        Exoplanet dictionary which can be loaded in and editted through `load_exo_dict`
    inst : str 
        Key which indicates with instrument 

    Returns
    -------
    dict
        Dictionary with output of pandexo. Key is the value of the parameter that was 
        looped through.  
    """
    #load in correct dict format
    inst_dict = load_mode_dict(inst)
    return {inst: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo})}


def run_pandexo(exo, inst, param_space = 0, param_range = 0,save_file = True,
                            output_path=os.getcwd(), output_file = ''):  
    """Submits multiple runs of pandexo in parallel. 
    
    Functionality: program contains functionality for running single or 
    multiple runs of PandExo 
    
    Parameters
    ---------- 
    exo : dict 
        exoplanet input dictionary 
    inst : dict or str or list of str
        instrument input dictionary OR LIST of keys (for allowable keys see `print_instruments()`
    param_space : str or 0 
        (Optional) Default is 0 = no exoplanet parameter space. To run through a parameter 
        specify which one need to specify two keys from exo dict with + in between. 
        i.e. observation+fraction
        star+temp
        planet+exopath
    param_range : list of str or list of float
        (Optional) Default = 0 An array or list over which to run the parameters space.
        i.e. array of temperatures if running through stellar temp or 
        array of files if running through planet models. Must specify param_space 
        if using this. 
    save_file : bool
        (Optional) Default = True saves file, False does not 
    output_path : str
        (Optional) Defaults to current working directory
    output_file : str 
        (Optional) Default is "singlerun.p" for single runs, "param_space.p" for exo parameter runs 
        or "instrument_run.p" for instrument parameter space runs. 
    
    Returns
    ------- 
    dict 
        For single run output will just be a single PandExo output dictionary 
        https://github.com/natashabatalha/PandExo/wiki/PandExo-Output
        For multiple runs the output will be organized into a list with each 
        a dictionary named by whatever you are looping through 
        i.e. [{'First temp': PandExoDict}, {'Second temp': PandExoDict}, etc..]
    
    Example
    -------
    For single run: 
    
    >>> a = run_pandexo(exo_dict, ['MIRI LRS'])
    
    For multiple instruments:
    
    >>> a = run_pandexo(exo_dict, ['MIRI LRS','NIRSpec G395H']
    
    Loop through a exoplanet parameter (stellar magnitude): 
    
    >>> a = run_pandexo(exo_dict, ['NIRSpec G395M'], 
            param_space ='star+mag',param_range = np.linspace(6,10,5))
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
        print 'Enter in format ["NIRSpec G140M"] or ["NIRISS SOSS","MIRI LRS"]' 
        return    
         
    #single instrument mode and single planet OR several planets  
     
    if len(inst)==1 and inst[0] != 'RUN ALL': 
        
        #start case of no parameter space run 
        if isinstance(param_space, (float, int)) or isinstance(param_range, (float, int)):
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
