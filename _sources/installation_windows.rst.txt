Windows Installing Guide for PandExo
====================================

Created by J. N. van Haastere (Physics and Astronomy bachelor's student),
University of Amsterdam and VU University Amsterdam.

Step 1 
``````

Create a dedicated conda environment; uninstalling other Python versions is
not necessary.

Step 2
``````

`Download most recent Anaconda version <https://www.anaconda.com/download/>`_.
PandExo supports Python 3.10 and newer.
    
Step 3
``````

Install Anaconda or Miniconda, then create and activate an environment:

.. code-block:: bat

    conda create -n pandexo python=3.11
    conda activate pandexo

.. note::
    If successfully installed, ``python --version`` reports Python 3.10 or newer.

Step 4
``````
Visual Studio build tools are not normally required. If pip must compile a
dependency, install the current `Microsoft C++ Build Tools
<https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_.
    
Step 5
``````

`Follow PandExo Installation guide in getting files <https://natashabatalha.github.io/PandExo/installation.html>`_

- get Pandeia data v2026p7
- get Pandeia PSFs v2026p7
- get stellar reference data for synphot/stsynphot (includes phoenix models and
  throughputs)


Step 6: Setting Reference Data
```````````````````````````````
Using an archiver (for example, 7-Zip), extract Pandeia, the Pandeia PSFs, and
the stellar reference data into separate folders. The 2026.7 JWST files are
available here:

- `Pandeia data v2026p7 JWST <https://stsci.app.box.com/v/pandeia-data-v2026p7-jwst>`_
- `Pandeia PSFs v2026p7 JWST <https://stsci.app.box.com/v/pandeia-psfs-v2026p7-jwst>`_

.. code-block:: bash

    C:/Users/USERNAME/pandeia_data-2026.7-jwst
    C:/Users/USERNAME/pandeia_psfs-2026.7-jwst
    C:/Users/USERNAME/stellar_reference_data/grp/redcat/trds

The variable names are ``pandeia_refdata``, ``PSF_DIR``, and ``PYSYN_CDBS``.

Set these in Windows environment variables / registry using SETX. Open CMD or Anaconda Prompt and run:   

.. code-block:: bash

    SETX pandeia_refdata "C:/Users/USERNAME/pandeia_data-2026.7-jwst"
    SETX PSF_DIR "C:/Users/USERNAME/pandeia_psfs-2026.7-jwst"
    SETX PYSYN_CDBS "C:/Users/USERNAME/stellar_reference_data/grp/redcat/trds"

The ``SETX`` command persists these variables for future terminals, but it
usually does not update the current open terminal. Close and reopen CMD,
PowerShell, or Anaconda Prompt before checking the values.

You can check this step by opening a new CMD and run > SET

Step 7: Installing needed Packages
``````````````````````````````````

Install PandExo Engine

.. code-block:: bash

    python -m pip install pandexo.engine

See the main :doc:`installation` guide if you run into problems.

Step 8: Run Test
````````````````
Navigate to the PandExo source checkout and run the smoke test.

.. code-block:: bash 

    cd '...'
    python -m pytest tests/test_run.py -q

Congrats on being persevering/stubborn enough to get PandExo to work on Windows!
There are some example Jupyter notebooks to try in ``PandExo/notebooks``.
