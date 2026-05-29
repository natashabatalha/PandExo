.. note::
    Before reading further, if you do not wish to install PandExo,\
    there is an online version of the code here at \
    `STScI's ExoCTK <https://exoctk.stsci.edu/pandexo/>`_. 


Pre-installation Data Download
==============================

PandExo requires downloading reference data for 1) JWST instrument info,
2) JWST PSFs, 3) stellar SEDs, and 4) normalization bandpasses. It also
requires setting up the ``pandeia_refdata``, ``PSF_DIR``, and ``PYSYN_CDBS``
environment variables.

JWST Reference Data
````````````````````
JWST Reference data has been updated to 2026.2.

.. warning::
    Reference data must match the software version. For example, Pandeia
    reference data 2026.2 should match the Pandeia software release 2026.2,
    and that should match PandExo 2026.2.

The 2026.2 JWST reference data is available here:

- `Pandeia data v2026p2 JWST <https://stsci.box.com/v/pandeia-data-v2026p2-jwst>`_
- `Pandeia PSFs v2026p2 JWST <https://stsci.box.com/v/pandeia-psfs-v2026p2-jwst>`_

More information on `pandeia installation can be found here <https://outerspace.stsci.edu/display/PEN/Pandeia+Engine+Installation>`_


After you have downloaded the reference data, create environment variable (`more resources on how to create environment variables are located here <https://natashabatalha.github.io/picaso/installation.html#create-environment-variable>`_). 

You can verify your installation by opening up a terminal with access to the conda installation, and type

.. code-block:: bash 

    python -c "import pandeia.engine; pandeia.engine.pandeia_version()"

If properly installed and configured, it should show the refdata version and
stellar reference-data directory, like this:

.. code-block:: bash 

    Pandeia Engine version:  2026.2
    Pandeia RefData version:  2026.2
    Pandeia PSFs version:    2026.2
    Synphot Data:  /your/data/directory/synphot


.. code-block:: bash 

    export pandeia_refdata=/path/to/pandeia-data-v2026p2-jwst
    export PSF_DIR=/path/to/pandeia-psfs-v2026p2-jwst

These commands set the variables for the current shell session. To make them
persist, add them to your shell startup file, such as ``~/.bashrc``,
``~/.bash_profile``, ``~/.zshrc``, or the equivalent file for your shell.

Stellar SEDs  
````````````
PandExo uses synphot/stsynphot for stellar spectra and PHOENIX model
interpolation. The PHOENIX reference atlas can be `downloaded here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar>`_.

Once untarred, the files will produce a directory tree of `grp/redcat/trds`. The pandeia.engine uses the contents of the `trds` directory.

**Environment variable: $PYSYN_CDBS must point to the trds directory (NOT grp)**

Create your environment variable:

.. code-block:: bash 

    echo 'export PYSYN_CDBS="$USRDIR/grp/redcat/trds"' >>~/.bash_profile

Normalization Files  
````````````````````
PandExo also needs the STScI/CDBS-style throughput files used for J/H/K
normalization bandpasses.

`Download the file here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar>`_

Once untarred this will also produce a directory tree of `grp/redcat/trds` with two folders `comp` and `mtab`. Place these folders into the folder you created above `$USRDIR/grp/redcat/trds` 

.. code-block:: bash 

    >> ls $USRDIR/grp/redcat/trds
    comp grid mtab

Now you should have three folders in your `trds` folder. 

Fortney+ 20210  Planet Grid (Optional)
````````````````````````````````````````
A user may want to install a grid of atmospheric models to simulate planet atmospheres. Some of the example notebooks use an atmospheric model grid. This grid can be obtained from:

- `The ExoCTK Website <https://github.com/ExoCTK/exoctk#obtain-the-exoctk-data>`_.

After downloading the Fortney files, create an environmental variable to point to them.

.. code-block:: bash 

    echo 'export FORTGRID_DIR="$USRDIR/fortney_models.db"' >>~/.bash_profile


Installation with Pip or Git
============================

Install with pip: 

.. code-block:: bash

    pip install pandexo.engine


OR Download PandExo's repository via Github. The Github also has helpful notebooks for getting started!

.. code-block:: bash

    git clone --recursive https://github.com/natashabatalha/pandexo
    cd pandexo
    pip install .



Final Test for Success
======================
 
Run the smoke test to confirm that your code is working:

.. code-block:: bash 

    python -m pytest tests/test_run.py -q


Troubleshooting-Common Errors
=============================



The Importance of Upgrading PandExo
===================================

It is crucial that your verison of PandExo is up to date. There were many critical updates in the reference files after launch, and as a result of the commissioning work. Updating PandExo requires three crucial steps: 


1) Verify pandexo.engine is Current 
````````````````````````````````````

.. code-block:: bash 

    pip install pandexo.engine --upgrade 


2) Verify pandeia.engine version compatible
````````````````````````````````````````````

Currently PandExo requires pandeia.engine==2026.2.

.. code-block:: bash 

    pip install pandeia.engine==2026.2

3) Grab pandeia.engine data 2026.2
````````````````````````````````

The 2026.2 JWST reference data is available here:

- `Pandeia data v2026p2 JWST <https://stsci.box.com/v/pandeia-data-v2026p2-jwst>`_
- `Pandeia PSFs v2026p2 JWST <https://stsci.box.com/v/pandeia-psfs-v2026p2-jwst>`_
