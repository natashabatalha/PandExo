.. warning::
    Before reading further, if you do not wish to install PandExo,\
    there is an online version of the code here at \
    `STScI's ExoCTK <https://exoctk.stsci.edu/pandexo/>`_. 


Pre-installation Data Download
==============================

PandExo requires: JWST instrument info and stellar SEDs. Users must set up these two environment variables before proceeding.

JWST Reference Data
````````````````````
JWST Reference data has been updated to 1.5!

.. warning::
    Reference data for OLD 1.3 is NOT backwards compatible with PandExo/Pandeia 1.4/1.5. The old reference data 
    can be found `here <http://ssb.stsci.edu/pandeia/engine/1.3/pandeia_data-1.3.tar.gz>`_ if you still wish to use PandExo v1.3. It is important to make sure that the version number of `PandExo` matches the version number of your reference data.

The new reference data is located `here for V1p6 <https://stsci.app.box.com/v/pandeia-refdata-v1p6>`_. Hopefully backwards compatibility issues will subside as Pandeia becomes more mature.


After you have downloaded the reference data, create environment variable: 

.. code-block:: bash 

    echo 'export pandeia_refdata="$USRDIR/pandeia_data"' >>~/.bash_profile

Stellar SEDs 
````````````
Likewise, the user may wish to install specific data for use with the PySynPhot package Pandeia uses. We will only be using the Phoenix stellar atlas, which can be downloaded `here <ftp://ftp.stsci.edu/cdbs/tarfiles/synphot5.tar.gz>`_.

More reference files can be downloaded via ftp: 

- ftp archive.stsci.edu
- username "anonymous"
- password is your e-mail address
- cd pub/hst/pysynphot
- download desired files. 

*Note* that the tar.gz files downloaded from the STScI anonymous FTP site will untar into the directory structure "grp/hst/cdbs", with the actual data files in an assortment of directories under "cdbs". pysynphot (and pandeia) expect that the "PYSYN_CDBS" environment variable will point to the "cdbs" directory. As such, you can either move the files out of "grp/hst/" to wherever you would like to store them, or point "PYSYN_CDBS" to "/path/to/data/files/grp/hst/cdbs" in order to allow pysynphot and pandeia to properly detect the reference files.

Finally, create your environment variable:

.. code-block:: bash 

    echo 'export PYSYN_CDBS="$USRDIR/path/to/data/files/grp/hst/cdbs"' >>~/.bash_profile


Installation with Pip or Git
============================

Install STScI specific packages

.. code-block:: bash

    conda config --add channels http://ssb.stsci.edu/astroconda
    conda install pyfftw

Now, install `pandexo`. 

.. code-block:: bash

    pip install pandexo.engine


OR Download PandExo's repository via Github. The Github also has helpful notebooks for getting started!

.. code-block:: bash

    git clone --recursive https://github.com/natashabatalha/pandexo
    cd pandexo
    python setup.py install



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

PyFFTW
````````
PyFFTW is needed to run PandExo. In order to run PyFFTW you need to also isntall fftw. To do so, it is necessary to do so through Homebrew, if you do not have conda. 

.. code-block:: bash 

    brew install fftw
    pip install pyfftw 

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



