.. raw:: html

   <div align="center">
   <img src="docs/logo.png" width="450px">
   </img>
   <br/>
   </div>
   <br/><br/>


Getting Started
---------------

Tools to help the community with planning exoplanet observations.

PandExo is both an online tool and a python package for generating instrument simulations of JWST’s NIRSpec, NIRCam, NIRISS and NIRCam and HST WFC3. It uses throughput calculations from STScI’s Exposure Time Calculator, Pandeia: Pandeia + Exoplanets = PandExo. This documentation contains information on how to download, install and analyze PandExo output.

Should I install PandExo or use the online interface? 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install if... 

- I will be using PandExo for more than 10 runs 
- I want to submit bash runs 
- I want to use the online GUI but don't want to be subject to any slow downs because of high user frequency 
- I want to be able to use plotting functions to analyze my output 

Install EVEN if... 

- I am scared of Python

Do not install if...

- I will only be using PandExo fewer than 10 times

Requires
~~~~~~~~

- Python >=3.8
- Installation: https://natashabatalha.github.io/PandExo/installation.html 
- PandExo now requires downloading **three folders**: 1) JWST instrument info (from pandeia), 2) stellar SEDs (pysynphot), and 3) normalization bandasses (pysynphot). See installation for furthur instructions: https://natashabatalha.github.io/PandExo/installation.html 


.. image:: https://zenodo.org/badge/67237418.svg
   :target: https://zenodo.org/badge/latestdoi/67237418


.. image:: https://travis-ci.org/natashabatalha/PandExo.svg?branch=master
    :target: https://travis-ci.org/natashabatalha/PandExo
