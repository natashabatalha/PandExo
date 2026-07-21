Windows Installation Guide for PandExo
======================================

Created by J. N. van Haastere (Physics and Astronomy bachelor's student),
University of Amsterdam and VU University Amsterdam.

Step 1: Install Conda
`````````````````````

Create a dedicated conda environment; uninstalling other Python versions is
not necessary.

`Download most recent Anaconda version <https://www.anaconda.com/download/>`_.
Miniconda is sufficient. PandExo requires Python 3.10 or later; continuous
integration currently tests Python 3.10 through 3.12, while Python 3.13 is
also exercised as a non-blocking compatibility check.

Step 2: Create an Environment
`````````````````````````````

Install Anaconda or Miniconda, then create and activate an environment:

.. code-block:: bat

    conda create -n pandexo python=3.11
    conda activate pandexo

.. note::
    Verify that ``python --version`` reports Python 3.10 or later.

Step 3: Optional Build Tools
````````````````````````````
Visual Studio build tools are not normally required. If pip must compile a
dependency, install the current `Microsoft C++ Build Tools
<https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_.

Step 4: Download Reference Data
```````````````````````````````
Using an archiver such as 7-Zip, extract the following archives. Keep Pandeia
reference data, PSFs, and synphot data in separate locations. The 2026.7 JWST
files are available here:

- `Pandeia data v2026p7 JWST <https://stsci.app.box.com/v/pandeia-data-v2026p7-jwst>`_
- `Pandeia PSFs v2026p7 JWST <https://stsci.app.box.com/v/pandeia-psfs-v2026p7-jwst>`_
- `PHOENIX stellar reference atlas <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar>`_
- `Normalization bandpasses <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar>`_

The stellar and normalization archives both extract a ``grp/redcat/trds``
tree. ``PYSYN_CDBS`` must name the PHOENIX archive's ``trds`` directory, not
its parent ``grp`` directory. Copy the normalization archive's ``comp`` and
``mtab`` *contents* into the matching directories under that same
``PYSYN_CDBS`` directory. Do not nest another ``grp/redcat/trds`` tree there.

Step 5: Set Reference-Data Variables
`````````````````````````````````````
The following PowerShell commands set each variable for the current terminal
and persist it for future terminals. Replace the example paths before running
them. Start a new PowerShell or Anaconda Prompt after using this step so the
persistent values are loaded automatically.

.. code-block:: powershell

    $env:pandeia_refdata = "$HOME\pandexo-data\pandeia_data-2026.7-jwst"
    $env:PSF_DIR = "$HOME\pandexo-data\pandeia_psfs-2026.7-jwst"
    $env:PYSYN_CDBS = "$HOME\pandexo-data\phoenix\grp\redcat\trds"

    [Environment]::SetEnvironmentVariable('pandeia_refdata', $env:pandeia_refdata, 'User')
    [Environment]::SetEnvironmentVariable('PSF_DIR', $env:PSF_DIR, 'User')
    [Environment]::SetEnvironmentVariable('PYSYN_CDBS', $env:PYSYN_CDBS, 'User')

Copy the normalization ``comp`` and ``mtab`` contents into ``PYSYN_CDBS``. In
this example, ``$normalizationTrds`` is the ``trds`` directory extracted from
the normalization archive:

.. code-block:: powershell

    $normalizationTrds = "$HOME\pandexo-data\normalization\grp\redcat\trds"
    Copy-Item -Recurse -Force "$normalizationTrds\comp\*" "$env:PYSYN_CDBS\comp\"
    Copy-Item -Recurse -Force "$normalizationTrds\mtab\*" "$env:PYSYN_CDBS\mtab\"

Verify the current terminal's variables and required files:

.. code-block:: powershell

    Get-ChildItem $env:PYSYN_CDBS
    Test-Path "$env:PYSYN_CDBS\comp\nonhst\bessell_j_003_syn.fits"
    Test-Path "$env:PYSYN_CDBS\grid\phoenix\catalog.fits"

For Command Prompt instead of PowerShell, set the current terminal and its
persistent replacement separately:

.. code-block:: bat

    set "pandeia_refdata=C:\Users\USERNAME\pandexo-data\pandeia_data-2026.7-jwst"
    setx pandeia_refdata "%pandeia_refdata%"

Repeat the pair of commands for ``PSF_DIR`` and ``PYSYN_CDBS``. ``setx`` only
affects future terminals, so retain ``set`` for the current Command Prompt.

Step 6: Install PandExo
```````````````````````

.. code-block:: powershell

    python -m pip install pandexo.engine

Confirm that PandExo can see the reference data:

.. code-block:: powershell

    python -c "import pandeia.engine; pandeia.engine.pandeia_version()"

Step 7: Validate the Installation
``````````````````````````````````

For a normal pip installation, run this lightweight installed-package smoke
test. It does not start the web application or run a simulation:

.. code-block:: powershell

    python -c "import pandexo.engine.justdoit as jdi; mode = jdi.load_mode_dict('NIRSpec Prism'); assert mode['configuration']['instrument']['disperser'] == 'prism'; print('PandExo configuration smoke test passed.')"

To run repository tests, first clone the source checkout and install its test
dependencies. This command applies only inside that checkout:

.. code-block:: powershell

    git clone https://github.com/natashabatalha/PandExo.git
    cd PandExo
    python -m pip install -e ".[test]"
    python -m pytest tests/test_import.py -q

Congrats on being persevering/stubborn enough to get PandExo to work on Windows!
There are some example Jupyter notebooks to try in ``PandExo/notebooks``.
