Py Dict Structure of JWST Output
================================

The PandExo output is organized into `Python Dictionaries`_. See
`example notebooks`_ for an explanation of how to manipulate these and
plot the most common outputs. Below is a breakdown of everything
contained in the PandExo output dictionary. 

FinalSpectrum
~~~~~~~~~~~~~

-  **spectrum\_w\_rand**: Planet spectrum with random noise added in
   units of *(Rp/Rs)2* or *(Fp/Fs)*
-  **spectrum**: Planet spectrum with no random noise
-  **error\_w\_floor**: Error with user defined noise floor. If no floor
   was specified, there is no floor.
-  **wave**: wavelength (microns)

.. code:: python 

    print(result['FinalSpectrum']['spectrum_w_rand'])

OriginalInput
~~~~~~~~~~~~~

-  **model\_wave**: original wavelength input by user
-  **model\_spec**: original spectrum input by user
-  **star\_spec**: out-of-transit stellar spectrum used by PandExo

.. code:: python 

    print(result['OriginalInput']['model_wave'])

warning
~~~~~~~

-  **Num Groups Reset?**: Reports whether optimization had to reset the group
   count to the detector minimum.
-  **Group Number Too Low?**: Reports a high requested full-well level with
   fewer than three groups, and cautions when a one-group ramp is used.
-  **Group Number Too High?**: Prints out warning if number of groups
   per integration exceeds 65536
-  **Saturated?**: This is an output directly taken from Pandeia's "hard
   saturation" flag. If there are any saturated pixels, it will alert
   you here. You can also see the saturation profile.
-  **Non linear?**: This is an output directly taken from Pandeia's
   "soft saturation" flag.
-  **% full well high?**: If you've set the saturation level over 80%,
   it will warn you.
-  **Minimum Integrations?**: Reports when the calculation has fewer than
   three in-transit integrations or when optimization reduces the group count
   to retain at least three.

Mode-specific calculations can add warnings for NIRCam readout optimization
and data excess, NIRSpec long exposures, MIRI slit-mode TSO use, or target
acquisition.

.. code:: python 

    print(result['warning']['Num Groups Reset?'])

PandeiaOutTrans
~~~~~~~~~~~~~~~

This is the raw output of Pandeia's simulation of the out of transit
observation. For a complete breakdown of these outputs go to `STScI's
Pandeia Documentation`_. 

- **sub\_reports** 
- **information** 
- **warnings** 
- **transform** 
- **2d** 
- **scalar** 
- **1d** 
- **input**

.. code:: python 

    print(result['PandeiaOutTrans']['information'])

RawData
~~~~~~~

-  **var\_in**: The variance of only the in transit data
-  **wave**: wavelength vector
-  **electrons\_in**, **electrons\_out**: Total in- and out-of-transit electrons
-  **e\_rate\_in**, **e\_rate\_out**: In- and out-of-transit electron rates
-  **electron\_per\_int** and **snr\_int**: Per-integration diagnostics
-  **error\_no\_floor**: Error without any noise floor
-  **var\_out**: The variance of only the out of transit data
-  **rn[out,in]** and **bkg[out,in]**: Read-noise and background diagnostics

.. code:: python
    
    print(result['RawData']['var_in'])

timing
~~~~~~

-  **Transit Duration** and **Number of Transits**
-  **Seconds per Frame** and **Time/Integration incl reset (sec)**
-  **Measurement Time per Integration (sec)**
-  **APT: Num Groups per Integration**
-  **Num Integrations In Transit** and **Num Integrations Out of Transit**
-  **APT: Exposures/Dith**, **APT: Num Integrations per Exposure**, and
   **APT: Num Integrations per Occultation**
-  **Observing Efficiency (%)**
-  **Transit+Baseline, no overhead (hrs)**

Multistripe and NIRCam modes add mode-specific stripe, on-source-time, and
data-excess fields.

.. code:: python 
  
    print(result['timing']['Seconds per Frame'])

input
~~~~~

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
 
    print(result['input']['Target Mag'])

HTML display fields
~~~~~~~~~~~~~~~~~~~

``apt_div``, ``calculation_div``, ``timing_div``, ``input_div``, and
``warnings_div`` contain HTML tables rendered on the website. Use the
corresponding dictionaries for programmatic analysis.

.. _Python Dictionaries: https://docs.python.org/3/tutorial/datastructures.html#dictionaries
.. _example notebooks: https://github.com/natashabatalha/PandExo/tree/master/notebooks
.. _STScI's Pandeia Documentation: https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial
