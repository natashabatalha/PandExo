import numpy as np
from pandeia.engine.instrument_factory import InstrumentFactory
from .pandexo import wrapper
from .load_modes import SetDefaultModes
import os
import pickle as pkl
from joblib import Parallel, delayed
import multiprocessing
import json
from .exomast import get_target_data
from astroquery.simbad import Simbad
import astropy.units as u
import copy

user_cores = multiprocessing.cpu_count()

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


def print_instruments(verbose=True):
    """Prints a list of the possible instrument templates to load
    """
    if verbose: print("Choose from the following:")
    return ALL.keys()

def getStarName(planet_name):
    """
    Given a string with a (supposed) planet name, this function returns the star name. For example:

    - If `planet_name` is 'HATS-5b' this returns 'HATS-5'.
    - If `planet_name` is 'Kepler-12Ab' this returns 'Kepler-12A'.
    
    It also handles the corner case in which `planet_name` is *not* a planet name, but a star name itself, e.g.:

    - If `planet_name` is 'HAT-P-1' it returns 'HAT-P-1'.
    - If `planet_name` is 'HAT-P-1  ' it returns 'HAT-P-1'.
    """

    star_name = planet_name.strip()

    # Check if last character is a letter:
    if str.isalpha(star_name[-1]):
        if star_name[-1] == star_name[-1].lower():
            star_name = star_name[:-1]
            
    # Return trimmed string:
    return star_name.strip()

def load_exo_dict(planet_name=None,pl_kwargs={}):
    """Loads in empty exoplanet dictionary for pandexo input

    Loads in empty exoplanet dictionary so that the user can manually edit different planet
    parameters. Must be done for every bash run. User must edit each keys within the dictionary.
    
    Parameters
    ----------
    planet_name : str
        (Optional) Planet name e.g. 'HD 189733 b' or 'HD189733b'
    pl_kwargs : dict 
        (Optional) : dict
        if you get an error that there is a missing field you can enter it in dictionary form using this
        e.g. pl_kwargs = {"Jmag":7}

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

    if not isinstance(planet_name, type(None)):
        planet_data = get_target_data(planet_name)[0]

        pandexo_input['star']['type'] = 'phoenix'
        pandexo_input['star']['temp'] = planet_data['Teff']
        pandexo_input['star']['metal'] = planet_data['Fe/H']
        pandexo_input['star']['logg'] = planet_data['stellar_gravity']
        Simbad.add_votable_fields('flux(H)')
        Simbad.add_votable_fields('flux(J)')

        star_name = getStarName(planet_name)


        if 'Jmag' in planet_data.keys():
            jmag = planet_data['Jmag']
        else: 
            try:
                jmag = Simbad.query_object(star_name)['FLUX_J'][0]

                #removing for how as blind str cutoffs are bug prone
                #if np.ma.is_masked(jmag):
                #    # Remove 'A' from star_name for systems with binary stars (e.g., WASP-77A)
                #    star_name = star_name[:-1]
                #    jmag = Simbad.query_object(star_name)['FLUX_J'][0]
            except: 
                jmag = pl_kwargs.get('Jmag',0)
                if jmag==0: 
                    raise Exception("Uh oh. Exo.MAST/simbad do not have a Jmag. Please enter it with pl_kwargs. E.g. pl_wargs={'Jmag':8} ")

        if 'Hmag' in planet_data.keys():
            hmag = planet_data['Hmag']
        else: 
            try:
                hmag = Simbad.query_object(star_name)['FLUX_H'][0]
                
                #removing for how as blind str cutoffs are bug prone
                #if np.ma.is_masked(hmag):
                #    # Remove 'A' from star_name for systems with binary stars (e.g., WASP-77A)
                #    star_name = star_name[:-1]
                #    hmag = Simbad.query_object(star_name)['FLUX_H'][0]
            except: 
                hmag = pl_kwargs.get('Hmag',0)
                if hmag==0: 
                    raise Exception("Uh oh. Exo.MAST/simbad do not have a Hmag. Please enter it with pl_kwargs. E.g. pl_wargs={'Hmag':8} ")


        pandexo_input["star"]["mag"] = jmag
        pandexo_input["star"]["ref_wave"] = 1.25
        pandexo_input["star"]["jmag"] = jmag
        pandexo_input["star"]["hmag"] = hmag
        #optinoal star radius
        pandexo_input["star"]["radius"] = planet_data['Rs']
        pandexo_input["star"]["r_unit"] = planet_data['Rs_unit'][0]+ planet_data['Rs_unit'][1:].lower()

       #optional planet radius/mass
        pandexo_input["planet"]["radius"] = planet_data['Rp']
        pandexo_input["planet"]["r_unit"] = planet_data['Rp_unit'][0]+ planet_data['Rp_unit'][1:].lower()
        try: 
            pandexo_input["planet"]["mass"] = planet_data['Mp']
            pandexo_input["planet"]["m_unit"] = planet_data['Mp_unit'][0]+ planet_data['Mp_unit'][1:].lower()

        except: 
            print("No mass found. Setting mass to np.nan. Mass is only required for model grids. If you want to enter one please enter it with pl_kwargs. E.g. pl_wargs={'Mp':1,'Mp_unit':'M_jupiter'}")
            pandexo_input["planet"]["mass"] = pl_kwargs.get('mass',np.nan)
            pandexo_input["planet"]["m_unit"] = pl_kwargs.get('m_unit',np.nan)

        pandexo_input["planet"]["transit_duration"] = planet_data['transit_duration']
        pandexo_input["planet"]["td_unit"] = planet_data['transit_duration_unit']
        depth = pandexo_input["planet"]["radius"]**2 / ((pandexo_input["star"]["radius"]
                                                    *u.Unit(pandexo_input["star"]["r_unit"]) )
                                                        .to(u.Unit(pandexo_input["planet"]["r_unit"]))).value**2
        pandexo_input["planet"]["depth"]      = depth
        if planet_data['inclination'] == None:
            inc = 90
        else:
            inc = planet_data['inclination']

        pandexo_input["planet"]["i"]          = inc
        pandexo_input["planet"]["ars"]        = planet_data['a/Rs']
        period = planet_data['orbital_period']
        period_unit = planet_data['orbital_period_unit']
        pandexo_input["planet"]["period"]     = (period*u.Unit(period_unit)).to(u.Unit('day')).value
        pandexo_input["planet"]["ecc"]        = planet_data['eccentricity']
        pandexo_input["planet"]["ecc"]        = planet_data['eccentricity']
        try:
            pandexo_input["planet"]["w"]      = float(planet_data['omega'] )
        except:
            pandexo_input["planet"]["w"]      = 90.
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
    return SetDefaultModes(inst).pick()

def get_thruput(inst, niriss=1, nirspec='f100lp'):
    """Returns complete instrument photon to electron conversion efficiency
    Pulls complete instrument photon to electron conversion efficiency
    (PCE) based on instrument key input

    Parameters
    ----------
    inst : str
        One of the instrument keys in `print_instruments`
    niriss : int
        (Optional) defines which niriss order you want (1 or 2)
    nirspec : str
        (Optional) for NIRISS G140M/H there are two available filters (f100lp and f070lp)
        if you are selecting G140M or G140H, this allows you to pick which one
    Returns
    -------
    dict
        Dictionary with wave solution and PCE

    Example
    -------
    >>> thru_dict = get_thruput('NIRISS SOSS_Or1')
    """

    #pull correct dictionary
    input_dict =  SetDefaultModes(inst).pick()
    conf = input_dict['configuration']
    conf['detector']['ngroup'] = 2

    if conf['instrument']['instrument'].lower() =='niriss':
        #pandeia handles slit losses inside the 2d engine. So, you need to account for the
        #extra .663 here
        conf["instrument"]["disperser"] = conf["instrument"]["disperser"] +'_'+str(niriss)
        i = InstrumentFactory(config=conf)
        wr = i.get_wave_range()
        wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
        pce = i.get_total_eff(wave)
        return {'wave':wave,'pce':0.663*pce}
    elif (conf['instrument']['instrument'].lower() =='nirspec') and ('g140' in conf["instrument"]["disperser"]):
        conf["instrument"]["filter"] = nirspec

    i = InstrumentFactory(config=conf)
    wr = i.get_wave_range()
    wave = np.linspace(wr['wmin'], wr['wmax'], num=500)
    pce = i.get_total_eff(wave)
    return {'wave':wave,'pce':pce}

def run_param_space(i,exo,inst,param_space, verbose=False):
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
    verbose : bool 
        (Optional) prints out checkpoints. Assuming the user does not want a load of 
        print statements for parallel runs

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
    return {name: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo}, verbose=verbose)}

def run_inst_space(inst,exo, verbose=False):
    """Changes inst dictionary and submits run

    This function is used to reset the instrument dictionary.

    Parameters
    ----------
    exo : dict
        Exoplanet dictionary which can be loaded in and editted through `load_exo_dict`
    inst : str
        Key which indicates with instrument
    verbose : bool 
        (Optional) prints out checkpoints. Assuming the user does not want a load of 
        print statements for parallel runs

    Returns
    -------
    dict
        Dictionary with output of pandexo. Key is the value of the parameter that was
        looped through.
    """
    #load in correct dict format
    inst_dict = load_mode_dict(inst)
    return {inst: wrapper({"pandeia_input": inst_dict , "pandexo_input":exo}, verbose=verbose)}


def run_pandexo(exo, inst, param_space = 0, param_range = 0,save_file = True,
                            output_path=os.getcwd(), output_file = '',num_cores=user_cores, verbose=True):
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
    verbose : bool 
        (Optional) For single runs, if false, it turns off all print statements. For parameter space 
        runs it is defaulted to never print statements out.

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

    >>> a = run_pandexo(exo_dict, ['MIRI LRS','NIRSpec G395H'])

    Loop through a exoplanet parameter (stellar magnitude):

    >>> a = run_pandexo(exo_dict, ['NIRSpec G395M'],
            param_space ='star+mag',param_range = np.linspace(6,10,5))
    """

    #single instrument mode with dictionary input OR single planet
    if type(inst) == dict:
        if verbose: print("Running Single Case w/ User Instrument Dict")
        results =wrapper({"pandeia_input": inst , "pandexo_input":exo}, verbose=verbose)
        if output_file == '':
            output_file = 'singlerun.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'wb'))
        return results

    #make sure inst is in list format.. makes my life so much easier
    try:
        if type(inst) != list:
            raise ValueError
    except ValueError:
        print('Instrument input is not dict so must be list')
        print('Enter in format ["NIRSpec G140M"] or ["NIRISS SOSS","MIRI LRS"]')
        return

    #single instrument mode and single planet OR several planets

    if len(inst)==1 and inst[0] != 'RUN ALL':

        #start case of no parameter space run
        if isinstance(param_space, (float, int)) or isinstance(param_range, (float, int)):
            if verbose: print("Running Single Case for: " + inst[0])
            inst_dict = load_mode_dict(inst[0])
            results =wrapper({"pandeia_input": inst_dict , "pandexo_input":exo}, verbose=verbose)
            if output_file == '':
                output_file = 'singlerun.p'
            if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'wb'))
            return results

        #if there are parameters to cycle through this will run
        if verbose: print("Running through exo parameters in parallel: " + param_space)
        #run the above function in parallel
        results = Parallel(n_jobs=num_cores)(delayed(run_param_space)(i,exo,inst[0],param_space) for i in param_range)

        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file = param_space + '.p'

        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'wb'))
        return results

    #run several different instrument modes and single planet
    if verbose: print("Running select instruments")
    if len(inst)>1:

        results = Parallel(n_jobs=num_cores)(delayed(run_inst_space)(i, exo) for i in inst)

        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file =  'instrument_run.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'wb'))
        return results

    #cycle through all options
    elif inst[0].lower() == 'run all':
        if verbose: print("Running through all instruments")
        results = Parallel(n_jobs=num_cores)(delayed(run_inst_space)(i, exo) for i in ALL.keys())

        #Default dump all results [an array of dictionaries] into single file
        #and return results immediately to user
        if output_file == '':
            output_file =  'instrument_run.p'
        if save_file: pkl.dump(results, open(os.path.join(output_path,output_file),'wb'))
        return results

def subarrays(inst):
  """function to show availalble subarrays and their times (in secons)

  Parameters
  ----------
  inst : str
    string of either niriss, nirspec, miri or nircam

  Returns
  -------
  dict
    dictionary with name of subarray as keys and time in seconds as entry
  """
  print("Subarray field stored in inst_dict['configuration']['detector']['subarray']")

  if inst.lower() == 'niriss':
    return {'substrip96':2.2129,'substrip256':5.4913}
  elif inst.lower() == 'nirspec':
    return {'sub1024a':0.451,'sub1024b':0.451,'sub2048':0.90156,'sub512':0.22572,'sub512s':0.14392}
  elif inst.lower() == 'miri':
    return {'slitlessprism':0.159}
  elif inst.lower()  == 'nircam':
    return {"subgrism64":0.34, "subgrism128":0.67, "subgrism256":1.34,
      "subgrism64 (noutputs=1)":1.3, "subgrism128 (noutputs=1)":2.6, "subgrism256 (noutputs=1)":5.2}
  else:
    raise Exception("Only instruments are niriss, nirspec, miri, nircam. Pick one.")

def dispersers(inst):
  """function to show available dispersers

  Parameters
  ----------
  inst : str
    string of either niriss, nirspec, miri or nircam

  Returns
  -------
  list
    lsit with available dispersers
  """
  print("Dispersers field stored in inst_dict['configuration']['instrument']['disperser']")
  if inst.lower() == 'niriss':
    return ['gr700xd']
  elif inst.lower() == 'nirspec':
    return ['g140m','g140h','g235m','g235h','g395m','g395h','prism']
  elif inst.lower() == 'miri':
    return ['p750l']
  elif inst.lower()  == 'nircam':
    return ['grismr']
  else:
    raise Exception("Only instruments are niriss, nirspec, miri, nircam. Pick one.")

def filters(inst):
    """Function to show availalbe filters

  Parameters
  ----------
  inst : str
    string of either niriss, nirspec, miri or nircam

  Returns
  -------
  list
    list with availalbe filters
  """
    print("Filters field stored in inst_dict['configuration']['instrument']['filter']")

    if inst.lower() == 'niriss':
        return ["clear","f277w"]
    elif inst.lower() == 'nirspec':
        return ['f070lp','f100lp','f170lp','f290lp','clear']
    elif inst.lower() == 'miri':
        print("No filters for miri lrs, type None, or null in filter field")
        return [None]
    elif inst.lower()  == 'nircam':
        return ['f322w2','f444w']
    else:
        raise Exception("Only instruments are niriss, nirspec, miri, nircam. Pick one.")

def grid_options(grid = 'fortney'):
    """Function to show available grid options

    PandExo now supports various grid options. Currently, the only one that is availalbe
    is the Fortney grid for giant planets. We will be implementing others, as they
    become available. It will become increasingly difficult for users to see what
    options are availalbe to them. This function should guide users to choose
    the grid that best fits their needs.

    Parameters
    ----------
    grid : str
        (Optional) string which 'fortney'

    """
    return
