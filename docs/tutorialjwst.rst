.. AUTO-GENERATED FILE. DO NOT EDIT.
.. Generated from notebooks/JWST_Running_Pandexo.ipynb by docs/generate_tutorial_rst.py.

JWST Tutorial
=============

Getting Started
---------------

Before starting here, all the instructions on the `installation
page <https://natashabatalha.github.io/PandExo/installation.html>`__
should be completed!

Here you will learn how to:

- set planet properties
- set stellar properties
- run default instrument modes
- adjust instrument modes
- run pandexo

.. code:: ipython3

    import warnings
    warnings.filterwarnings('ignore')
    import pandexo.engine.justdoit as jdi
    import pandexo.engine.justplotit as jpi
    import numpy as np
    import os
    #pip install pandexo.engine --upgrade

Confirm that the Pandeia engine, reference data, and PSF library all
match the version required by PandExo.

.. code:: ipython3

    print(os.environ['pandeia_refdata'] )
    import pandeia.engine
    pandeia.engine.pandeia_version()

Load blank exo dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~

To start, load in a blank exoplanet dictionary with empty keys. You will
fill these out for yourself in the next step.

.. code:: ipython3

    exo_dict = jdi.load_exo_dict()
    print(exo_dict.keys())
    #print(exo_dict['star']['w_unit'])

Editing Input Dictionaries
~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit exoplanet observation inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Review each input value and its units carefully so the calculation is
reproducible.

.. code:: ipython3

    exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
    exo_dict['observation']['sat_unit'] = '%'
    exo_dict['observation']['noccultations'] = 1 #number of transits
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                                 #without having to redo the calculation
    exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                        #'frac' : fraction of time in transit versus out = in/out
                                                        #'total' : total observing time (seconds)
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)

    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                                 #to a wavelength dependent noise floor solution (units are ppm)

Edit exoplanet host star inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note... If you select ``phoenix`` you **do not** have to provide a
``starpath``, ``w_unit`` or ``f_unit``, but you **do** have to provide a
``temp``, ``metal`` and ``logg``. If you select ``user`` you **do not**
need to provide a ``temp``, ``metal`` and ``logg``, but you **do** need
to provide units and starpath.

Option 1) Grab stellar model from database
""""""""""""""""""""""""""""""""""""""""""

.. code:: ipython3

    #OPTION 1 get start from database
    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 10.0              #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs

Option 2) Input as dictionary or filename
"""""""""""""""""""""""""""""""""""""""""

.. code:: ipython3

    # Create an arbitrary user-supplied stellar spectral shape. PandExo
    # normalizes it to the magnitude and reference wavelength below.
    wl = np.linspace(0.5, 15, 3000)

    # A user spectrum can be supplied as a dictionary or a two-column file.
    starflux = {'f': (wl / 1.25)**-2, 'w': wl}

    exo_dict['star']['type'] = 'user'
    exo_dict['star']['mag'] = 10.0              #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25
    exo_dict['star']['starpath'] = starflux
    exo_dict['star']['w_unit'] = 'um'
    exo_dict['star']['f_unit'] = 'jy'

Edit exoplanet inputs using one of three options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1) user specified
2) constant value
3) select from grid

1) Edit exoplanet planet inputs if using your own model
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code:: ipython3

    exo_dict['planet']['type'] = 'user'                      # Upload a user spectrum.
    wavelength = np.linspace(0.6, 5.3, 500)
    spectrum = np.full_like(wavelength, 0.015)
    exo_dict['planet']['exopath'] = {'f': spectrum, 'w': wavelength}

    exo_dict['planet']['w_unit'] = 'um'                      # Other options include nm, Angs, and sec (phase curves).
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'               # The other option is fp/f*.
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration
    exo_dict['planet']['td_unit'] = 's'                      #Any unit of time in accordance with astropy.units can be added

2) Users can also add in a constant temperature or a constant transit depth
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code:: ipython3

    exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration
    exo_dict['planet']['td_unit'] = 's'
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'            #Any unit of distance in accordance with astropy.units can be added here
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'              #Same deal with astropy.units here
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit

    #ORRRRR....
    #if you wanted to instead to secondary transit at constant temperature
    #exo_dict['planet']['f_unit'] = 'fp/f*'
    #exo_dict['planet']['temp'] = 1000

3) Select from a grid
"""""""""""""""""""""

Currently, PandExo supports the Fortney et al. (2010) hot-Jupiter grid.

.. code:: ipython3

    exo_dict['planet']['type'] = 'grid'                #tells pandexo you want to pull from the grid
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0
    exo_dict['planet']['td_unit'] = 's'
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


Load in instrument dictionary (OPTIONAL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step 2 is optional because PandExo has the functionality to
automatically load in instrument dictionaries. Skip this if you plan on
observing with one of the following and want to use the subarray with
the smallest frame time and the readout mode with 1 frame/1 group
(standard): - NIRCam F444W - NIRSpec Prism - NIRSpec G395M - NIRSpec
G395H - NIRSpec G235H - NIRSpec G235M - NIRCam F322W2 - NIRCam DHS -
NIRSpec G140M - NIRSpec G140H - MIRI LRS - NIRISS SOSS

.. code:: ipython3

    jdi.print_instruments()

.. code:: ipython3

    inst_dict = jdi.load_mode_dict('NIRSpec G140H')

    # Loading an instrument dictionary lets you personalize template fields.
    # Templates use subarrays with the lowest frame times and readout modes
    # with one frame per group; change these fields when needed.

    #Try printing this out to get a feel for how it is structured:

    print(inst_dict['configuration'])

.. code:: ipython3

    #Another way to display this is to print out the keys
    inst_dict.keys()

Don't know what instrument options there are?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: ipython3

    print("SUBARRAYS")
    print(jdi.subarrays('nirspec'))

    print("FILTERS")
    print(jdi.filters('nircam'))

    print("DISPERSERS")
    print(jdi.dispersers('nirspec'))

.. code:: ipython3

    #you can try personalizing some of these fields

    inst_dict["configuration"]["detector"]["ngroup"] = 'optimize'   #running "optimize" will select the maximum
                                                                    #possible groups before saturation.
                                                                    #You can also write in any integer between 2-65536

    inst_dict["configuration"]["detector"]["subarray"] = 'sub2048'   #change the NIRSpec subarray



Adjusting the Background Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may want to think about adjusting the background level of your
observation, based on the position of your target. PandExo two options
and three levels for the position:

- ``ecliptic`` or ``minzodi``
- ``low``, ``medium``, ``high``

.. code:: ipython3

    inst_dict['background'] = 'ecliptic'
    inst_dict['background_level'] = 'high'

Running NIRISS SOSS Order 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^

PandExo only will extract a single order at a time. By default, it is
set to extract Order 1. Below you can see how to extract the second
order.

**NOTE!** Users should be careful with this calculation. Saturation will
be limited by the **first** order. Therefore, first run one calculation
with ``ngroup='optimize'`` for Order 1, then use that group count for
the Order 2 calculation.

.. code:: ipython3

    inst_dict = jdi.load_mode_dict('NIRISS SOSS')
    inst_dict['strategy']['order'] = 2
    inst_dict['configuration']['detector']['subarray'] = 'substrip256'
    ngroup_from_order1_run = 2
    inst_dict["configuration"]["detector"]["ngroup"] = ngroup_from_order1_run

Running PandExo
~~~~~~~~~~~~~~~

You have **four options** for running PandExo. All use
**jdi.run_pandexo**.

``jdi.run_pandexo(exo, inst, param_space=0, param_range=0, save_file=True, output_path=None, output_file='', num_cores=None, verbose=True)``

Option 1- Run single instrument mode, single planet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you forget which instruments are available run
**jdi.print_instruments()** and pick one

.. code:: ipython3

    jdi.print_instruments()

.. code:: ipython3

    # MIRI/LRS needs source and planet spectra that extend into the mid-IR.
    # The Fortney grid example above only extends to about 5 microns, so use
    # the constant-depth setup here.
    exo_dict['planet']['type'] = 'constant'
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0
    exo_dict['planet']['td_unit'] = 's'
    exo_dict['planet']['radius'] = 1
    exo_dict['planet']['r_unit'] = 'R_jup'
    exo_dict['star']['radius'] = 1
    exo_dict['star']['r_unit'] = 'R_sun'
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'

    result = jdi.run_pandexo(exo_dict, ['MIRI LRS'], save_file=False)

Note, you can turn off print statements with ``verbose=False``

Analyzing Output
^^^^^^^^^^^^^^^^

The output dictionary can be passed directly to ``justplotit`` to
inspect the simulated transit or eclipse spectrum and the expected
precision. Here the spectrum is binned to ``R=100`` for display.

.. code:: ipython3

    wave, spectrum, error = jpi.jwst_1d_spec(result, R=100, num_tran=1, model=False, x_range=[5, 14])

Option 2- Run single instrument mode (with user dict), single planet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the same as option 1, but it accepts an instrument dictionary so
you can customize a supported configuration.

.. code:: ipython3

    inst_dict = jdi.load_mode_dict('NIRSpec G140H')
    # Personalize the subarray.
    inst_dict["configuration"]["detector"]["subarray"] = 'sub2048'
    result = jdi.run_pandexo(exo_dict, inst_dict, save_file=False)

Option 3- Run several modes, single planet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use several modes from **print_instruments()** options.

.. code:: ipython3

    # Choose selected modes.
    result = jdi.run_pandexo(exo_dict,['NIRSpec G140M','NIRSpec G235M','NIRSpec G395M'],
                   save_file=False, verbose=True)
    #run all
    #result = jdi.run_pandexo(exo_dict, ['RUN ALL'], save_file = False)

Option 4- Run single mode, several planet cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use a single mode from **print_instruments()** and explore parameter
space for **any** value in the exoplanet dictionary. The example below
uses a few planet radii so the notebook remains runnable without
external model files.

Use ``planet``, ``star``, or ``observation`` followed by the key to
vary.

For example: ``planet+exopath``, ``star+temp``, ``star+metal``,
``star+logg``, or ``observation+sat_level``.

.. code:: ipython3

    # Loop over a small, self-contained planet-radius grid.
    # This keeps the example runnable without requiring external model files.
    exo_dict['planet']['type'] = 'constant'
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'

    radius_grid = np.linspace(0.8, 1.2, 3)
    radius_results = jdi.run_pandexo(exo_dict, ['NIRCam F444W'],
                                    param_space='planet+radius',
                                    param_range=radius_grid,
                                    save_file=False,
                                    verbose=True)
