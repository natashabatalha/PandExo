.. warning::
    Before reading further, if you do not wish to install PandExo,\
    there is an online version of the code here here at \
    `PSU Science <http://pandexo.science.psu.edu:1111>`_ or here at NASA GSFC. 

Before we start, the scripts written to install PandExo assume the user has installed `Anaconda <https://www.continuum.io/downloads>`_. For diehard pip users, know that STScI tools will be only distributed on conda so PandExo will not be the only package to use it. 

If you are dead set against it you will need to follow the installation guidelines `here <https://natashabatalha.github.io/PandExo/installation.html>`_. 

Pre-installation Data Download
==============================

ss mentioned, the user must provide their own reference data files in order for the engine to understand the instrument with which the calculations will be performed. First step to getting PandExo running is downloading four files. Download them and put them in the same directory. 

JWST Reference Data
```````````````````
For JWST, the user may download the required reference data `here <http://ssb.stsci.edu/pandeia/engine/1.0/pandeia_data-1.0.tar.gz>`_. Unpack and store the data in a choice location.

Stellar SEDs 
````````````
Likewise, the user may wish to install specific data for use with the PySynPhot package Pandeia uses. We will only be using the Phoenix stellar atlas, which can be downloaded `here <ftp://ftp.stsci.edu/cdbs/tarfiles/synphot5.tar.gz>`_. 

PandExo startup bash script
````````````````````````````
A startup script and a python test script was created to do an automatic install of PandExo. You can download both of them from a google drive located `here <https://drive.google.com/open?id=0ByRvbM3zW0d-d2xoSmNMc1RtdTQ>`_.


Installation for Users With Conda
====================================

Open up the file `pandexo.sh`. There are two lines in the beginning which tell PandExo where your reference data will go. Currently this script is set to unpack and store the reference files in your HOME directory. 

If you are not okay with it, **change** $HOME in the first two lines of pandexo.sh to whatever you want the data to be.

Then run

.. code-block:: bash 

    sh pandexo.sh
    Starting TEST run
    Running Single Case for: NIRSpec G140H
    Optimization Reqested: Computing Duty Cycle
    Finished Duty Cycle Calc
    Starting Out of Transit Simulation
    End out of Transit
    Starting In Transit Simulation
    End In Transit
    SUCCESS

Hopefully you have had success but if not the most likely error you will get is the following:

COORDS.PY: TypeError: 'float' object cannot be interpreted as an index
```````````````````````````````````````````````````````````````````````
This is an error within Pandeia, which has not yet been fixed by STScI folk. You can easily fix it by finding where the file /pandeia/engine/coords.py is and changing line 36:

.. code-block:: python 
   
    ones = np.ones((ny, nx))

To this: 

.. code-block:: python

    ones = np.ones((int(ny), int(nx)))


Installation for Users With PIP
===============================

You should have already downloaded the two data tar files, pandexo.sh and run_test.py. If not, go back to `Pre Installation Data Download`. 

Then, start by installing PandExo, which should automatically isntall Pandeia and other dependencies 

.. code-block:: bash
    
    pip install pandexo.engine

OR Download PandExo's repository via Github: 

.. code-block:: bash

    git clone --recursive https://github.com/natashabatalha/pandexo
    cd pandexo
    python setup.py install

Set Environment
```````````````
Open up the file `pandexo.sh`. The first two lines in the beginning tell PandExo where your reference data will go via the variable `$USRDIR`. Currently `$USRDIR` is equal to your HOME directory, meaning your big data files will be untarred there. 

If you are not okay with it, **change** $HOME in the first two lines of pandexo.sh to whatever you want the data to be. It should look like this: 

.. code-block:: bash 

    USRDIR=/I/Want/My/Data/Here
    echo 'USRDIR=/I/Want/My/Data/Here' >>~/.bash_profile

NEXT, since pandexo.sh is set for conda users, delete everything after line 22. We will have to do these manually since we don't have conda. Once you have done these two things you can go ahead and run: 

.. code-block:: bash 

    sh pandexo.sh

Dependencies
````````````
PyFFTW is needed to run PandExo. In order to run PyFFTW you need to also isntall fftw. To do so, it is necessary to do so through Homebrew, if you do not have conda. 

.. code-block:: bash 

    brew install fftw
    pip install pyfftw  

Lastly, Python 2.7 users will need to install multiprocessing: 

.. code-block:: bash 
    
    pip install multiprocessing

Finally try to run the test file to see if there are any additional problems: 

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

Hopefully you have had success but if not the most likely error you will get is the following:

COORDS.PY: TypeError: 'float' object cannot be interpreted as an index
```````````````````````````````````````````````````````````````````````
This is an error within Pandeia, which has not yet been fixed by STScI folk. You can easily fix it by finding where the file /pandeia/engine/coords.py is and changing line 36:

.. code-block:: python 
   
    ones = np.ones((ny, nx))

To this: 

.. code-block:: python

    ones = np.ones((int(ny), int(nx)))

TODOS
===== 

1. Add error messages to the pandas output page 
