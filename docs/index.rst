.. engine documentation master file, created by
   sphinx-quickstart on Sun Dec 18 17:56:59 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: logo.png
   :alt: PandExo logo.
   :width: 450px
   :align: center

PandExo Documentation
=====================

Tools to help the community with planning exoplanet observations.

.. note::
    Use the `PandExo web interface on STScI's ExoCTK
    <https://exoctk.stsci.edu/pandexo/>`_ for browser-based calculations.

PandExo is both an online tool and a Python package for generating instrument
simulations of JWST's NIRSpec, NIRCam, NIRISS, and MIRI, and HST's WFC3. It
uses throughput calculations from STScI's Exposure Time Calculator, Pandeia:
Pandeia + Exoplanets = PandExo.

Start with :doc:`installation` to configure the local package and required
reference data. Then follow the :doc:`tutorialjwst` or :doc:`tutorialhst`, use
:doc:`jwstinput` to select a supported instrument configuration, and consult
:doc:`jwstdict` when analyzing a completed calculation.

Citation: `PandExo: A Community Tool for Transiting Exoplanet Science with JWST and HST <https://ui.adsabs.harvard.edu/abs/2017PASP..129f4501B>`_

Contents:

.. toctree::
   :maxdepth: 1
   
   includeme 
   installation
   installation_windows
   tutorialjwst
   jwstinput
   jwstdict
   tutorialhst
   engine


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
