JWST Output Dictionary
======================

PandExo returns a nested `Python dictionary`_. See the `example notebooks`_
for common analysis and plotting patterns. The exact fields vary by instrument
mode and selected options, so inspect the result you receive rather than
assuming every calculation has every key:

.. code:: python

    print(result.keys())
    print(result['FinalSpectrum'].keys())

The arrays in ``FinalSpectrum`` and ``RawData`` are NumPy arrays. Wavelengths
are in microns unless a field name or the associated metadata states otherwise.

FinalSpectrum
~~~~~~~~~~~~~

-  **spectrum\_w\_rand**: Planet spectrum with random noise added, in
   ``(Rp/Rs)^2`` or ``Fp/Fs``.
-  **spectrum**: Planet spectrum without random noise.
-  **error\_w\_floor**: Error with user defined noise floor. If no floor
   was specified, there is no floor.
-  **wave**: Wavelength in microns.

.. code:: python 

    print(result['FinalSpectrum']['spectrum_w_rand'])

OriginalInput
~~~~~~~~~~~~~

-  **model\_wave**: Original wavelength array supplied by the user.
-  **model\_spec**: Original planet spectrum supplied by the user.
-  **star\_spec**: Out-of-transit stellar spectrum used by PandExo.

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
acquisition. Warning values and their availability vary by mode.

.. code:: python 

    print(result['warning']['Num Groups Reset?'])

PandeiaOutTrans
~~~~~~~~~~~~~~~

This is the raw output of Pandeia's out-of-transit simulation. For a complete
breakdown, see `STScI's Pandeia documentation`_. Its nested fields are
Pandeia-version and mode dependent.

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

-  **var\_in**: Variance of the in-transit data.
-  **wave**: Wavelength vector in microns.
-  **electrons\_in**, **electrons\_out**: Total in- and out-of-transit electrons.
-  **e\_rate\_in**, **e\_rate\_out**: In- and out-of-transit electron rates.
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
-  **Saturation Level (electrons)**
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

Saved ``.p`` results are Python pickle files. Only load files from trusted
sources, because unpickling an untrusted file can execute arbitrary code.

.. code:: python

    import pickle

    with open("singlerun.p", "rb") as handle:
        result = pickle.load(handle)  # Load only files from trusted sources.

.. _Python dictionary: https://docs.python.org/3/tutorial/datastructures.html#dictionaries
.. _example notebooks: https://github.com/natashabatalha/PandExo/tree/master/notebooks
.. _STScI's Pandeia Documentation: https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial
