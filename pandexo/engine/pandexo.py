def wrapper(dictinput):
    """
    Top level function which calls either jwst, hst or wfirst noise simulation. 
    
    :param dictinput: dictionary containing instrument parameters and exoplanet specific parameters. {"pandeia_input":dict1, "pandexo_input":dict1}
    :type inputdict: dict
    :returns: dict specific to observatory requested
    rtype: dict
    
    .. note:: You should not run simulations through this. It is much easier to run simulations through 
        either **run_online** or **justdoit**. **justdoit** contains functions 
        to create input dictionaries and **run_online** contains web forms to create input dictionaries.
    .. seealso:: 
        Module :justdoit:`justdoit`
            Documentation of the :mod:`justdoit` 
        Module :run_online:`run_online`
            Documentation of the :mod:`run_online`
    """
    pandexo_input = dictinput['pandexo_input']    	
    telescope = pandexo_input['telescope']
    
    if telescope=='jwst':
        from jwst import compute_full_sim
        return compute_full_sim(dictinput)
    elif telescope=='hst':
        from hst import wfc3_TExoNS
        return wfc3_TExoNS(dictinput)
    elif telescope=='wfirst':
        return     
    else:
        print "INVALID TELESCOPE. PandExo only accepts: jwst, hst, wfirst"