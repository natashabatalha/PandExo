Py Dict Structure of JWST Output
================================

The PandExo output is organized into `Python Dictionaries`_. See
`example notebooks`_ for an explanation of how to manipulate these and
plot the most common outputs. Below is a breakdown of everything
contained in the PandExo output dictionary. 

FinalSpectrum (contains 4 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **spectrum\_w\_rand**: Planet spectrum with random noise added in
   units of *(Rp/Rs)2* or *(Fp/Fs)*
-  **spectrum**: Planet spectrum with no random noise
-  **error\_w\_floor**: Error with user defined noise floor. If no floor
   was specified, there is no floor.
-  **wave**: wavelength (microns)

.. code:: python 

    print dict['FinalSpectrum']['spectrum_w_rand'] 

OriginalInput (contains 2 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **model\_wave**: original wavelength input by user
-  **model\_spec**: original spectrum input by user

.. code:: python 

    print dict['OriginalInput']['model_wave']

Warning (contains 6 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Num Groups Reset?**: Before PandExo simulates in and out of transit
   observations, it computes a single integration with 2 groups in order
   to figure out how many additional groups it can add before hitting
   saturation. If it computes a number less than 2, it resets the number
   of groups to 2.
-  **Group Number Too Low?**: Prints out warning if number of groups is
   less than 5 and the saturation level less than 60%.
-  **Group Number Too High?**: Prints out warning if number of groups
   per integration exceeds 65536
-  **Saturated?**: This is an output directly taken from Pandeia’s “hard
   saturation” flag. If there are any saturated pixels, it will alert
   you here. You can also see the saturation profile.
-  **Non linear?**: This is an output directly taken from Pandeia’s
   “soft saturation” flag.
-  **% full well high?**: If you’ve set the saturation level over 80%,
   it will warn you.

.. code:: python 

    print dict['Warning']['Num Groups Reset?']

PandeiaOutTrans (contains 9 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the raw output of Pandeia’s simulation of the out of transit
observation. For a complete breakdown of these outputs go to `STScI’s
Pandeia Documentation`_. 

- **sub\_reports** 
- **information** 
- **warnings** 
- **transform** 
- **2d** 
- **scalar** 
- **1d** 
- **input**
- **3d**

.. code:: python 

    print dict['PandeiaOutTrans']['information']

RawData (contains 6 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~

-  **var\_in**: The variance of only the in transit data
-  **wave**: wavelength vector
-  **flux\_in**: Flux of the in transit data in units of e-/s
-  **flux\_out**: Flux of the out of transit data in units of e-/s
-  **error\_no\_floor**: Error without any noise floor
-  **var\_out**: The variance of only the out of transit data

.. code:: python
    
    print dict['RawData']['var_in']    

Timing (contains 10 keys)
~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Seconds per Frame**
-  **Number of Transits**
-  **Observing Efficiency (%)** = (num groups - 1)/(num groups + 1)
-  **Num Integrations Out of Transit**
-  **On Source Time**
-  **Exposure Time Per Integration (secs)**
-  **Reset time Plus TA time (hrs)**: Target acquisition time is assumed
   to be 30 minutes
-  **Num Integrations In Transit**
-  **Num Groups per Integration**
-  **Num Integrations per Occultation**

.. code:: python 
  
    print dict['Timing']['Seconds per Frame']

Input (contains 10 keys)
~~~~~~~~~~~~~~~~~~~~~~~~

-  **Target Mag**
-  **Readmode**
-  **Disperser**
-  **Filter**
-  **Instrument**
-  **Mode**
-  **Saturation Level (electons)**
-  **Aperture**
-  **Subarray**
-  **Primary/Secondary**

.. code:: python
 
    print dict['Input']['Target Mag']

Timing, Warning, and Input Divs (all contain html script)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here are the html scripts for the three tables that are rendered on the
website output page.

.. _Python Dictionaries: http://www.tutorialspoint.com/python/python_dictionary.htm
.. _example notebooks: https://github.com/natashabatalha/PandExo/wiki/Notebooks
.. _STScI’s Pandeia Documentation: https://jwst.etc.stsci.edu
