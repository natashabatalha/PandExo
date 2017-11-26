JWST Tutorial 
=============

This tutorial can be downloaded as a iPython notebook on the PandExo Github. You also always have the option of running PandExo UI by running in `bash`

.. code:: bash

    start_pandexo 

And then going to https://locahost:1111

In python though, start by loading in some universal packages, as well as PandExo.

.. code:: python

    import warnings
    warnings.filterwarnings('ignore')
    import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
    import numpy as np
    import os

Editting Input Dictionaries
---------------------------

Step 1) Load in a blank exoplanet dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To start, load in a blank exoplanet dictionary with empty keys. You will
fill these out for yourself in the next step.

.. code:: python

    exo_dict = jdi.load_exo_dict()

Edit exoplanet observation inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Editting each keys are annoying. But, do this carefully or it could
result in nonsense runs

.. code:: python

    exo_dict['observation']['sat_level'] = 80                  #saturation level in percent of full well 
    exo_dict['observation']['sat_unit'] = '%'                  # other option = 'e' for electrons
    exo_dict['observation']['noccultations'] = 2               #number of transits 
    exo_dict['observation']['R'] = None                        #fixed binning. I usually suggest ZERO binning.. you can always bin later 
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0        #time spent observing out of transit, make sure to speciy units
    exo_dict['observation']['baseline_unit'] = 'total'         #total obersving time, other option 'frac' = in/out 
    exo_dict['observation']['noise_floor'] = 0                 #this can be a fixed level or it can be a filepath 
                                                               #to a wavelength dependent noise floor solution (units are ppm)

Edit exoplanet host star inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note… If you select ‘phoenix’ you **do not** have to provide a starpath,
w\_unit or f\_unit, but you **do** have to provide a temp, metal and
logg. If you select ‘user’ you **do not** need to provide a temp, metal
and logg, but you **do** need to provide units and starpath.

.. code:: python

    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs

Edit exoplanet inputs using one of  three options: 1) user specified, 2) constant value, 3) select from grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1) Edit exoplanet planet inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    exo_dict['planet']['type'] ='user'
    exo_dict['planet']['exopath'] = 'wasp12b.txt'
    exo_dict['planet']['w_unit'] = 'cm'                      #other options include "um","nm" ,"Angs", "secs" (for phase curves)
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'               #other options are 'fp/f*' 
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's'                      #Any unit of time in accordance with astropy.units can be added

2) Users can also add in a constant temperature or a constant transit depth 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    exo_dict['planet']['type'] = 'constant' 
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit 
    #ORRRRR....
    #if you wanted to instead to secondary transit at constant temperature 
    exo_dict['planet']['f_unit'] = 'fp/f*' 
    exo_dict['planet']['temp'] = 1000

3) Select from grid
^^^^^^^^^^^^^^^^^^^
NOTE: Currently only the fortney grid for hot Jupiters from Fortney+2010 is supported. Holler though, if you want another grid supported   

.. code:: python 

    exo_dict['planet']['type'] = 'grid'
    exo_dict['planet']['temp'] = 1000                 #grid: 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500
    exo_dict['planet']['chem'] = 'noTiO'              #options: 'noTiO' and 'eqchem', noTiO is chemical eq. without TiO
    exo_dict['planet']['cloud'] = 'ray10'               #options: nothing: '0', 
    #                                                   Weak, medium, strong scattering: ray10,ray100, ray1000
    #                                                   Weak, medium, strong cloud: flat1,flat10, flat100
    exo_dict['planet']['mass'] = 1
    exo_dict['planet']['m_unit'] = 'M_jup'            #Any unit of mass in accordance with astropy.units can be added here
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration 
    exo_dict['planet']['td_unit'] = 's' 

Step 2) Load in instrument dictionary (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 2 is optional because PandExo has templates for the most 
common modes. Within those templates the subarrays selected are the ones 
with the lowest frame times, and the read modes selected are the ones with one frame
per 1 group (standard). If this suits you see Option 1 below. There is additional 
fine tuning that can be done within each of these observing modes.  
As a first pass, I'd suggest skipping this for now. Then come back and fine tune after your first run. 

    - NIRCam F444W 
    - NIRSpec Prism 
    - NIRSpec G395M 
    - NIRSpec G395H 
    - NIRSpec G235H 
    - NIRSpec G235M 
    - NIRCam F322W2 
    - NIRSpec G140M 
    - NIRSpec G140H 
    - MIRI LRS 
    - NIRISS SOSS

.. code:: python

    jdi.print_instruments()

::

    Choose from the following:
    ['NIRCam F444W', 'NIRSpec Prism', 'NIRSpec G395M', 'NIRCam F322W2', 'NIRSpec G395H', 'NIRSpec G235H', 'NIRSpec G235M', 'NIRSpec G140M', 'NIRSpec G140H', 'MIRI LRS', 'NIRISS SOSS', 'WFC3 G141']

.. code:: python

    inst_dict = jdi.load_mode_dict('NIRSpec Prism')

Change subarray: 

.. code:: python

    inst_dict["configuration"]["detector"]["subarray"] = 'sub512'

Don't know which options are available?? 

..code:: python 

    #subarrays
    print(jdi.subarrays('nirspec'))
    #filters
    print(jdi.filters('nircam'))
    #dispersers
    print(jdi.dispersers('nirspec'))

The last thing to note is that PandExo (by Default) optimizes the amount of groups 
that can fit into an integration. If you want to set your own number of groups you can 
change that too. 

..code:: python 

    inst_dict["configuration"]["detector"]["ngroup"] = 5

Running PandExo Command Line
----------------------------

You have **four options** for running PandExo. All of them are accessed
through attribute **jdi.run\_pandexo**. See examples below.

``jdi.run_pandexo(exo, inst, param_space = 0, param_range = 0,save_file = True, output_path=os.getcwd(), output_file = '')``

Option 1- Run single instrument mode, single planet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you forget which instruments are available run
**jdi.print\_instruments()** and pick one

.. code:: python

    jdi.print_instruments()

::

    Choose from the following:
    ['NIRCam F444W', 'NIRSpec Prism', 'NIRSpec G395M', 'NIRCam F322W2', 'NIRSpec G395H', 'NIRSpec G235H', 'NIRSpec G235M', 'NIRSpec G140M', 'NIRSpec G140H', 'MIRI LRS', 'NIRISS SOSS_Or1', 'NIRISS SOSS_Or2', 'WFC3 G141']

.. code:: python

    result = jdi.run_pandexo(exo_dict,['NIRCam F322W2'])

::

    Running Single Case for: NIRCam F322W2
    Computing Duty Cycle
    Finished Duty Cucle Calc
    Starting Out of Transit Simulation
    End out of Transit
    Starting In Transit Simulation
    End In Transit

Option 2- Run single instrument mode (with user dict), single planet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the same thing as option 1 but instead of feeding it a list of
keys, you can feed it a instrument dictionary (this is for users who
wanted to simulate something NOT pre defined within pandexo)

.. code:: python

    inst_dict = jdi.load_mode_dict('NIRSpec G140H')
    #personalize subarray
    inst_dict["configuration"]["detector"]["subarray"] = 'sub2048'
    result = jdi.run_pandexo(exo_dict, inst_dict)

::

    Running Single Case w/ User Instrument Dict
    Computing Duty Cycle
    Finished Duty Cucle Calc
    Starting Out of Transit Simulation
    End out of Transit
    Starting In Transit Simulation
    End In Transit

Option 3- Run several modes, single planet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use several modes from **print\_isntruments()** options.

.. code:: python

    #choose select 
    result = jdi.run_pandexo(exo_dict,['NIRCam F444W','NIRCam F322W2','MIRI LRS'],
                   output_path = '/Users/nbatalh1/Desktop/JWSTFUN')
    #run all 
    result = jdi.run_pandexo(exo_dict, ['RUN ALL'], save_file = False)

Option 4- Run single mode, several planet cases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a single modes from **print\_instruments()** options. But explore
parameter space with respect to **any** parameter in the exo dict. The
example below shows how to loop over several planet models

You can loop through anything in the exoplanet dictionary. It will be
planet, star or observation followed by whatever you want to loop
through in that set.

i.e. planet+exopath, star+temp, star+metal, star+logg,
observation+sat\_level.. etc

.. code:: python

    #looping over different exoplanet models 
    jdi.run_pandexo(exo_dict, ['NIRCam F444W'], param_space = 'planet+exopath',
                    param_range = os.listdir('/Users/nbatalha1/all_my_models_here'),
                   output_path = '/Users/nbatalh1/Desktop/JWSTFUN')

    #looping over different stellar temperatures 
    jdi.run_pandexo(exo_dict, ['NIRCam F444W'], param_space = 'star+temp',
                    param_range = np.linspace(5000,8000,2),
                   output_path = '/Users/nbatalh1/Desktop/JWSTFUN')

    #looping over different saturation levels
    jdi.run_pandexo(exo_dict, ['NIRCam F444W'], param_space = 'observation+sat_level',
                    param_range = np.linspace(.5,1,5),
                   output_path = '/Users/nbatalh1/Desktop/JWSTFUN')

Running PandExo GUI
-------------------
The same interface that is available online is also available for use on your machine. 
Using the GUI is very simple and good alternative if editing the input dictionaries is 
confusing. 

.. code:: python 

    import pandexo.engine.run_online as ro
    ro.main()

Then open up your favorite internet browser and go to: http://localhost:1111

.. note:: Some WebApp functions may not be available. For example, the precomputed noise simulations and transmission and emission modeling are only available online. 
                   
Analyzing Output
----------------

There are pre computed functions for analyzing most common outputs. You can also explore 
the dictionary structure yourself. 

.. code:: python

    import pandexo.engine.justplotit as jpi
    import pickle as pk

Plot 1D Data with Errorbars
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple plotting options exist within `jwst\_1d\_spec`

1. Plot a single run

.. code:: python

    #load in output from run
    out = pk.load(open('singlerun.p','r'))
    #for a single run 
    x,y, e = jpi.jwst_1d_spec(out, R=100, num_tran=10, model=False, x_range=[.8,1.28])

.. image:: jwst_1d_spec.png

2. Plot several runs from parameters space run 

.. code:: python

    #load in output from multiple runs
    multi = pk.load(open('three_nirspec_modes.p','r'))

    #get into list format 
    list_multi = [multi[0]['NIRSpec G140M'], multi[1]['NIRSpec G235M'], multi[2]['NIRSpec G395M']]

    x,y,e = jpi.jwst_1d_spec(list_multi, R=100, model=False, x_range=[1,5])

.. image:: jwst_1d_spec_multi.png

Plot Noise & More
~~~~~~~~~~~~~~~~~

Several functions exist to plot various outputs.

See also **jwst\_1d\_bkg**, **jwst\_1d\_snr**, **jwst\_1d\_flux**,

.. code:: python

    x,y = jpi.jwst_noise(out)

.. image:: jwst_noise.png

Plot 2D Detector Profile
~~~~~~~~~~~~~~~~~~~~~~~~

See also **jwst\_2d\_sat** to plot saturation profile

.. code:: python

    data = jpi.jwst_2d_det(out)

.. image:: jwst_1d_det.png