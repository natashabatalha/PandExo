Installation
==============
.. warning::
    Before reading further, if you do not wish to install `PandExo`,\
    there is an online version of the code at \
    PSU Science and NASA GSFC

Installation Procedure 
----------------------

1. Download PandExo's repository: 

.. code-block:: bash

    git clone --recursive https://github.com/natashabatalha/pandexo
    cd pandexo
    python setup.py install

2. Download the `Phoenix Stellar Atlas <ftp://ftp.stsci.edu/cdbs/tarfiles/synphot5.tar.gz>`_
in order to easily pull Stellar SED's from phoenix database. 

.. code-block:: bash

    mkdir pysynphot_data
    mv synphot5.tar.gz pysynphot_data
    tar -xvf synphot5.tar.gz
    export PYSYN_CDBS=USRDIR/pysynphot_data

Troubleshooting Pandeia
-----------------------

PYFFTW
``````
Many users experience issues when downloading Pandeia because of it's dependency \
on `pyfftw`. If you experience this problem try these steps:

- If you do not have a non-LLVM based GCC installation on your system, you can obtain one from here but gcc 5.1 does not produce a usable FFTW installation so make sure you download **gcc 4.9 or below**

- STScI created the following script to successfully install `pyfftw`

.. code-block:: bash

    mv $(which gcc) $(which gcc).orig
    curl -O https://bitbucket.org/api/2.0/snippets/jhunkeler/R7gy5/3265aea27175817087ab4a39c21157d926f8afc3/files/build_fftw.sh
    chmod +x build_fftw.sh
    ./build_fftw

Pandeia Reference Data
``````````````````````
Did you forget to point to Pandeia's reference data?? 

- Make sure PandExo knows where the Pandeia reference data is: 

.. code-block:: bash

    export pandeia_refdata=USRDIR/pandeia_data

