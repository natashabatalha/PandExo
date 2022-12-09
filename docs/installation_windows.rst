Windows Installing Guide for PandExo
====================================

Created by J. N. van Haastere (Physics and Astronomy Bachelorstudent)
Univesity of Amsterdam and VU University Amsterdam

Step 1 
``````

You may want to uninstall older python versions to prevent troubles along the way. 

Step 2
``````

`Download most recent Anaconda version <https://www.anaconda.com/download/>`_ (tested with Python 3.6.4 64 bit) 
    
Step 3
``````

Install Anaconda as main python and add to path environment

.. note::
    If succesfully installed you can open CMD and type `>>python` and `>>quit()` to check version

Step 4
``````
`Install Build Tools for Visual Studio 2017 <https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017>`_.
    
Step 5
``````

`Follow PandExo Installation guide in getting files <https://natashabatalha.github.io/PandExo/installation.html>`_

- get pandeia_data      
- get pysynphot_data (includes phoenix models and throughputs)


Step 6: Setting Reference Data
```````````````````````````````
Using an achiver (eg. Winrar, 7Zip) extract Pandeia and Pysynphot data in two seperate folders.

.. code-block:: bash

    C:/Users/USERNAME/pandeia_data
    C:/Users/USERNAME/pysynphot_data

Check the variable names

.. code-block:: bash

    echo 'export VARIABLE_NAME = ".."'


In case nothing changes these should be `PYSYN_CDBS` & `pandeia_refdata`

Set these in Windows environment variables / registry using SETX. Open CMD or Anaconda Prompt and run:   

.. code-block:: bash

    SETX PYSYN_CDBS 'C:/Users/USERNAME/pysynphot_data'
    SETX pandeia_refdata 'C:/Users/USERNAME/pysynphot_data'

You can check this step by opening a new CMD and run > SET

Step 7: Installing needed Packages
``````````````````````````````````

Install PandExo Engine

.. code-block:: bash

    pip install pandexo.engine

Check troubleshooting if you run into problems. `Troubleshooting here <https://natashabatalha.github.io/PandExo/installation.html#troubleshooting-common-errors>`_.

Step 8: Run Test
````````````````
Navigate to your pandexo-master and run the run_test. Check in the `installation <https://natashabatalha.github.io/PandExo/installation.html#pandexo-startup-bash-script>`_ guide on the exact expected output.

.. code-block:: bash 

    cd '...'
    python run_test.py

Congrats on being persevering/stubborn enough to get it on to work on Windows!
There are some test Jupyter files to try out in `pandexo-master/notebooks`




