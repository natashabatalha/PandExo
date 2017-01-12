Installation
==============
.. warning::
    Before reading further, if you do not wish to install `PandExo`,\
    there is an online version of the code at \
    PSU Science and NASA GSFC. 

Installation Procedure 
----------------------

1. Easiest way to Install: 

.. code-block:: bash
    
    pip install pandexo.engine

OR Download PandExo's repository: 

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

Once you do this your untarred file should automatically have this structure. If you use wget you might 
get something different. So double check this. 

.. code-block:: bash

    ls pysynphot_data/grid/phoenix/
    .DS_Store     catalog.fits  phoenixm05/   phoenixm15/   phoenixm25/   phoenixm35/   phoenixp03/   
    AA_README     phoenixm00/   phoenixm10/   phoenixm20/   phoenixm30/   phoenixm40/   phoenixp05/

3. Download the `JWST Reference Data <http://ssb.stsci.edu/pandeia/engine/1.0/pandeia_data-1.0.tar.gz>`_ . 
This is a big file (6 gigs) so think carefully about where you want to store it. Don't accidentally download 
it on your Mac Air then wonder why you can't save a 32 Kb doc file. 

Then make sure you untar and point to the file so PandExo knows where it is. 

.. code-block:: bash

    tar xf pandeia_data-1.0.tar.gz 
    export pandeia_refdata=USRDIR/pandeia_data

Your Pandeia data reference file should look like this: 

.. code-block:: bash 

    ls pandeia_data/
    .DS_Store      README.md      background/    extinction/    jwst/          sed/           strategy/      
    .git/          VERSION_PSF    devtools/      hst/           normalization/ source/        wfirst/ 
    
Troubleshooting
---------------

Problems with PYFFTW?
`````````````````````
Many users experience issues when downloading Pandeia because of it's dependency \
on `pyfftw`. If you experience this problem try these steps:

- If you do not have a non-LLVM based GCC installation on your system, you can obtain one from here but gcc 5.1 does not produce a usable FFTW installation so make sure you download **gcc 4.9 or below**

- STScI created the following script to successfully install `pyfftw`

.. code-block:: bash

    mv $(which gcc) $(which gcc).orig
    curl -O https://bitbucket.org/api/2.0/snippets/jhunkeler/R7gy5/3265aea27175817087ab4a39c21157d926f8afc3/files/build_fftw.sh
    chmod +x build_fftw.sh
    ./build_fftw

Can't find Pandeia Reference Data
`````````````````````````````````
This usually looks like NoneType errors. 

- Make sure PandExo knows where the Pandeia reference data is: 

.. code-block:: bash

    export pandeia_refdata=USRDIR/pandeia_data

