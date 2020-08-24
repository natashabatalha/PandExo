def wrapper(dictinput, verbose=False):
    """Tpo level function to call either jwst, hst, wfirst
    
    Top level function which calls either jwst, hst or wfirst noise simulation. 
    
    Parameters
    ----------
    dictinput : 
        dictionary containing instrument parameters and exoplanet 
        specific parameters. {"pandeia_input":dict1, "pandexo_input":dict1}
    verbose : bool 
        (Optional) prints out checkpoints throughout code

    Returns
    -------
    dict 
        output specific to observatory requested

    
    Notes
    -----
    You should not run simulations through this. It is much easier to run simulations through 
    either **run_online** or **justdoit**. **justdoit** contains functions 
    to create input dictionaries and **run_online** contains web forms to create input dictionaries.
    
    See Also
    --------
    pandexo.engine.justdoit : gives ability to simply submit runs 
    pandexo.engine.run_online : submit runs through user interface 
    """
    pandeia_input = dictinput['pandeia_input']    	
    telescope = pandeia_input['telescope']

    if telescope=='jwst':
        from .jwst import compute_full_sim
        return compute_full_sim(dictinput, verbose=verbose)
    elif telescope=='hst':
        from .hst import compute_sim_hst
        return compute_sim_hst(dictinput, verbose=verbose)
    elif telescope=='wfirst':
        return
    else:
        print("INVALID TELESCOPE. PandExo only accepts: jwst, hst, wfirst")
