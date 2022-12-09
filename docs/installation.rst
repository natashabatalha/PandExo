.. warning::
    Before reading further, if you do not wish to install PandExo,\
    there is an online version of the code here at \
    `STScI's ExoCTK <https://exoctk.stsci.edu/pandexo/>`_. 


Pre-installation Data Download
==============================

PandExo requires downloading **three folders**: 1) JWST instrument info, 2) stellar SEDs, and 3) normalization bandasses. It also requires setting up **two environment variables**. 

JWST Reference Data
````````````````````
JWST Reference data has been updated to 2.0!

.. warning::
    Reference data for PandExo/Pandeia 1.X is not compatible with PandExo/Pandeia 2.0. Please always ensure that your reference data is up to date with the software (e.g. reference data 2.0 should match software release 2.0)

The new reference data is located `here for v2p0 <https://stsci.app.box.com/v/pandeia-refdata-v2p0-jwst>`_. More information on `pandeia installation can be found here <https://outerspace.stsci.edu/display/PEN/Pandeia+Engine+Installation>`_


After you have downloaded the reference data, create environment variable (`more resources on how to create environment variables are located here <https://natashabatalha.github.io/picaso/installation.html#create-environment-variable>`_). 

.. code-block:: bash 

    echo 'export pandeia_refdata="$USRDIR/pandeia_data"' >>~/.bash_profile

Stellar SEDs  
````````````
PandExo uses Pysynphot's Phoenix stellar atlas, which can be `downloaded here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_pheonix-models_multi_v3_synphot5.tar>`_.

Once untarred, the files will produce a directory tree of `grp/redcat/trds`. The pandeia.engine uses the contents of the `trds` directory.

**Environment variable: $PYSYN_CDBS must point to the trds directory (NOT grp)**

Create your environment variable:

.. code-block:: bash 

    echo 'export PYSYN_CDBS="$USRDIR/grp/redcat/trds"' >>~/.bash_profile

Normalization Files  
````````````````````
New to PandExo 2.0, **users now have to download the master table of all pysynphot throughput tables.** 

`Download the file here <https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar>`_

Once untarred this will also produce a directory tree of `grp/redcat/trds` with two folders `comp` and `mtab`. Place these folders into the folder you created above `$USRDIR/grp/redcat/trds` 

.. code-block:: bash 

    ls $USRDIR/grp/redcat/trds
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
 
There is a `run_test.py` in the `github`. Test that you're code is working: 

.. code-block:: bash 

    python run_test.py
    Starting TEST run
    Running Single Case for: NIRSpec G140H
    Optimization Reqested: Computing Duty Cycle
    Finished Duty Cycle Calc
    Starting Out of Transit Simulation
    End out of Transit
    Starting In Transit Simulation
    End In Transit
    SUCCESS


Troubleshooting-Common Errors
=============================


Multiprocessing
````````````````
Python 2.7 users might need to install multiprocessing

.. code-block:: bash 
    
    pip install multiprocessing

RecursionError: maximum recursion depth exceeded while calling a Python object
````````````````````````````````````````````````````````````````````````````````

There is a known bug with Python 3.6 and Sphinx <1.6. Before updating or installing pandexo do the following:

PIP USERS:

.. code-block:: bash 

    pip install sphinx==1.5.6

CONDA USERS:

.. code-block:: bash 

    conda install sphinx=1.5.6

TypeError: super() argument 1 must be type
````````````````````````````````````````````

This is the same error above with Sphinx, but for Python 2.7 users. The fix is the same: 

PIP USERS:

.. code-block:: bash 

    pip install sphinx==1.5.6

CONDA USERS:

.. code-block:: bash 

    conda install sphinx=1.5.6
    

The Importance of Upgrading PandExo
===================================

It is crucial that your verison of PandExo remain up to date. Especially through commissioning and leading up to launch, there may be crucial changes to the code or the reference data. Updating PandExo requires three crucial steps. 

Verify Reference Data is Current
````````````````````````````````
The link to the reference data is located on `Pandeia's PyPI page <https://pypi.python.org/pypi/pandeia.engine/>`_. Before doing a large batch of calculations, make sure that you have this version. 

Verify pandeia.engine is Current
````````````````````````````````

.. code-block:: bash 

    pip install pandeia.engine --upgrade 

Verify pandexo.engine is Current 
````````````````````````````````

.. code-block:: bash 

    pip install pandexo.engine --upgrade 



