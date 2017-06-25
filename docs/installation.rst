Installation
==============
.. warning::
    Before reading further, if you do not wish to install PandExo,\
    there is an online version of the code here here at \
    `PSU Science <http://pandexo.science.psu.edu:1111>`_ or here at NASA GSFC. 

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

2. Download the Phoenix Stellar Atlas `FROM THIS LINK <ftp://ftp.stsci.edu/cdbs/tarfiles/synphot5.tar.gz>`_
in order to easily pull Stellar SED's from the phoenix database. Then type the following commands. 
It might be helpful to add the export command to your ~/.bashrc file. 

.. code-block:: bash

    mkdir pysynphot_data
    mv synphot5.tar.gz pysynphot_data
    tar -xvf synphot5.tar.gz
    export PYSYN_CDBS=USRDIR/pysynphot_data

Once you do this your untarred file should automatically have this structure. STScI changed the structure of 
this file but now how pysynphot calls the data. So you will most likely have to change this file structure to what
is below: 

.. code-block:: bash

    ls pysynphot_data/grid/phoenix/
    .DS_Store     catalog.fits  phoenixm05/   phoenixm15/   phoenixm25/   phoenixm35/   phoenixp03/   
    AA_README     phoenixm00/   phoenixm10/   phoenixm20/   phoenixm30/   phoenixm40/   phoenixp05/

3. Download the JWST Reference Data `FROM HERE <http://ssb.stsci.edu/pandeia/engine/1.0/pandeia_data-1.0.tar.gz>`_ . 
This is a big file (6 gigs) so think carefully about where you want to store it. Don't accidentally download 
it on your Mac Air then wonder why you can't save a 32 Kb doc file. 

Then make sure you untar and point to the file so PandExo knows where it is. Like above, it might 
be helpful to put this in your ~/.bashrc file. 

.. code-block:: bash

    tar xf pandeia_data-1.0.tar.gz 
    export pandeia_refdata=USRDIR/pandeia_data-1.0

Your Pandeia data reference file should look like this: 

.. code-block:: bash 

    ls pandeia_data/
    .DS_Store      README.md      background/    extinction/    jwst/          sed/           strategy/      
    .git/          VERSION_PSF    devtools/      hst/           normalization/ source/        wfirst/ 
    
Troubleshooting
---------------
TypeError: 'float' object cannot be interpreted as an index
```````````````````````````````````````````````````````````
This is an error within Pandeia, which has not yet been fixed by STScI folk. You can easily fix it by finding where the file /pandeia/engine/coords.py is and changing line 36:

.. code-block:: python 
   
    ones = np.ones((ny, nx))
    ones = np.ones((int(ny), int(nx)))

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

If that doesn't work Zach Berta-Thompson pointed out that this worked for him: 

.. code-block:: bash 

    brew install fftw
    pip install pyfftw

There are several different conda distributions of pyfftw however, if 
you install version <0.10 it may downgrade 
your version of numpy. Therefore, a version that is greater than 0.10 is best. Numpy should be 
numpy>=1.12 

.. code-block:: bash 

    conda install -c spectraldns pyfftw=0.10.4 

Can't find Pandeia Reference Data
`````````````````````````````````
This usually looks like NoneType errors. 

- Make sure PandExo knows where the Pandeia reference data is: 

.. code-block:: bash

    export pandeia_refdata=USRDIR/pandeia_data
    
Problems Installing Pysynphot
`````````````````````````````

If you are having problems with this 
you can use the astroconda distribution located `here <http://astroconda.readthedocs.io/en/latest/installation.html#install-astroconda>`_. 

Or if you are using conda: 

.. code-block:: bash

    conda install -c astropy photutils=0.3

Problems with Multiprocessing
`````````````````````````````

Multiprocessing seems to throw errors if you are using Python 3. Jonathan Fraine pointed out that 
multiprocessing is automatically included in Python 3. Therefore, if you delete the multiprocessing from 
setup.py your problem will go away. 

To-Do
-----

Below are a list of task items. Please check below for your request before notifying me. 

1. Add error messages to the pandas output page 
