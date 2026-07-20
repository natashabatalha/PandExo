.. note::
    Before reading further, if you do not wish to install PandExo, there is an
    online version of the code here at \
    `STScI's ExoCTK <https://exoctk.stsci.edu/pandexo/>`_. 


Pre-installation Data Download
==============================

PandExo requires downloading reference data for 1) JWST instrument info,
2) JWST PSFs, 3) stellar SEDs, and 4) normalization bandpasses. It also
requires setting up the ``pandeia_refdata``, ``PSF_DIR``, and ``PYSYN_CDBS``
environment variables.

JWST Reference Data
````````````````````
JWST reference data has been updated to 2026.7.

.. warning::
    Reference data must match the software version. For example, Pandeia
    reference data 2026.7 should match the Pandeia software release 2026.7,
    and that should match PandExo 2026.7.

The 2026.7 JWST reference data is available here:

- `Pandeia data v2026p7 JWST <https://stsci.app.box.com/v/pandeia-data-v2026p7-jwst>`_
- `Pandeia PSFs v2026p7 JWST <https://stsci.app.box.com/v/pandeia-psfs-v2026p7-jwst>`_

More information is available in the `Pandeia Engine installation guide <https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia+Engine+Installation>`_.


After downloading the reference data, add the environment variables to your
shell startup file so they persist across terminal sessions. For Bash, use
``~/.bashrc`` or ``~/.bash_profile``; for Zsh, use ``~/.zshrc``. Add these lines
to the appropriate file, replacing the example paths:

.. code-block:: bash

    export pandeia_refdata=/path/to/pandeia_data-2026.7-jwst
    export PSF_DIR=/path/to/pandeia_psfs-2026.7-jwst

Source the startup file to make the changes available in your current terminal,
for example:

.. code-block:: bash

    source ~/.bashrc

Use the path to the file you edited, such as ``source ~/.bash_profile`` or
``source ~/.zshrc``. Opening a new terminal will also load the saved variables.

Stellar SEDs  
````````````
PandExo uses synphot/stsynphot for stellar spectra and PHOENIX model
interpolation. The PHOENIX reference atlas can be `downloaded here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar>`_.

Once untarred, the files will produce a directory tree of `grp/redcat/trds`. The pandeia.engine uses the contents of the `trds` directory.

**Environment variable: $PYSYN_CDBS must point to the trds directory (NOT grp)**

Add this environment variable to the same shell startup file and source that
file again:

.. code-block:: bash 

    export PYSYN_CDBS=/path/to/grp/redcat/trds

Normalization Files  
```````````````````
PandExo also needs the STScI/CDBS-style throughput files used for J/H/K
normalization bandpasses.

`Download the file here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar>`_

Once untarred, this archive also produces a ``grp/redcat/trds`` directory with
``comp`` and ``mtab`` subdirectories. Merge their contents into the existing
``comp`` and ``mtab`` directories under ``$PYSYN_CDBS``. Do not create a second
``grp/redcat/trds`` tree inside ``$PYSYN_CDBS``.

For example, if ``NORMALIZATION_TRDS`` is the ``trds`` directory extracted
from the normalization archive:

.. code-block:: bash

    NORMALIZATION_TRDS=/path/to/extracted/grp/redcat/trds
    rsync -a "$NORMALIZATION_TRDS/comp/" "$PYSYN_CDBS/comp/"
    rsync -a "$NORMALIZATION_TRDS/mtab/" "$PYSYN_CDBS/mtab/"

The ``$PYSYN_CDBS`` directory should now contain ``comp``, ``grid``, and
``mtab`` (among other synphot directories). Verify that layout and the files
PandExo uses for J/H/K normalization and PHOENIX spectra:

.. code-block:: bash

    ls "$PYSYN_CDBS"
    ls "$PYSYN_CDBS/comp/nonhst/bessell_j_003_syn.fits"
    ls "$PYSYN_CDBS/comp/nonhst/bessell_h_004_syn.fits"
    ls "$PYSYN_CDBS/comp/nonhst/bessell_k_003_syn.fits"
    ls "$PYSYN_CDBS/grid/phoenix/catalog.fits"

Verify the complete reference-data setup with:

.. code-block:: bash

    python -c "import pandeia.engine; pandeia.engine.pandeia_version()"

If properly installed and configured, it should show the matching versions and
stellar reference-data directory, like this:

.. code-block:: text

    Pandeia Engine version:  2026.7
    Pandeia RefData version: 2026.7
    Pandeia PSFs version:    2026.7
    Synphot Data:            /path/to/grp/redcat/trds

Fortney+ 2010 Planet Grid (Optional)
````````````````````````````````````
A user may want to install a grid of atmospheric models to simulate planet atmospheres. Some of the example notebooks use an atmospheric model grid. This grid can be obtained from:

The required ``fortney_models.db`` database is included in the `ExoCTK Fortney
data archive <https://data.science.stsci.edu/redirect/JWST/ExoCTK/compressed/fortney.tar.gz>`_.
Download and extract that archive, or install ExoCTK and run
``exoctk.utils.download_exoctk_data('fortney')``.

After downloading the Fortney files, add this optional environment variable to
the same shell startup file and source that file again:

.. code-block:: bash 

    export FORTGRID_DIR=/path/to/exoctk_data/fortney/fortney_models.db


Installation with Pip or Git
============================

Install with pip: 

.. code-block:: bash

    python -m pip install pandexo.engine


OR Download PandExo's repository via Github. The Github also has helpful notebooks for getting started!

.. code-block:: bash

    git clone https://github.com/natashabatalha/PandExo.git
    cd PandExo
    python -m pip install .



Final Test for Success
======================

For a normal pip installation, run this installed-package smoke test. It loads
PandExo's bundled NIRSpec configuration and does not download data, start the
web application, or run a simulation:

.. code-block:: bash

    python -c "import pandexo.engine.justdoit as jdi; mode = jdi.load_mode_dict('NIRSpec Prism'); assert mode['configuration']['instrument']['disperser'] == 'prism'; print('PandExo configuration smoke test passed.')"

If you cloned the PandExo source repository and want to run its regression
test, install the test dependency and run the test from the repository root:

.. code-block:: bash

    python -m pip install -e ".[test]"
    python -m pytest tests/test_run.py -q


The Importance of Upgrading PandExo
===================================

It is crucial that your version of PandExo is up to date. There were many critical updates in the reference files after launch and as a result of commissioning work. Updating PandExo requires three crucial steps:


1) Verify pandexo.engine is Current 
```````````````````````````````````

.. code-block:: bash 

    python -m pip install --upgrade pandexo.engine


2) Verify pandeia.engine version compatible
```````````````````````````````````````````

Currently PandExo requires pandeia.engine==2026.7.

.. code-block:: bash 

    python -m pip install pandeia.engine==2026.7

3) Grab pandeia.engine data 2026.7
``````````````````````````````````

The 2026.7 JWST reference data is available here:

- `Pandeia data v2026p7 JWST <https://stsci.app.box.com/v/pandeia-data-v2026p7-jwst>`_
- `Pandeia PSFs v2026p7 JWST <https://stsci.app.box.com/v/pandeia-psfs-v2026p7-jwst>`_
