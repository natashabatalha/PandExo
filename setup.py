#!/usr/bin/env python

# This sample setup.py can be used as a template for any project using d2to1.
# Simply copy this file and, if desired, delete all the comments.  Also remove
# the 'namespace_packages' and 'packages' arguments to setup.py if this project
# does not contain any packages beloning to a namespace package.

# This import statement attempts to import the setup() function from setuptools
# (this replaces the setup() one uses from distutils.core when using plain
# distutils).
#
# If the import fails (the user doesn't have the setuptools package) it then
# uses the ez_setup bootstrap script to install setuptools, then retries the
# import.  This is common practice for packages using setuptools.
try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup
import sys 
PY_V = sys.version_info
# The standard setup() call.  Notice, however, that most of the arguments
# normally passed to setup() are absent.  They will instead be read from the
# setup.cfg file using d2to1.
#
# In order for this to work it is necessary to specify setup_requires=['d2to1']
# If the user does not have d2to1, this will boostrap it.  Also require
# stsci.distutils to use any of the common setup_hooks included in
# stsci.distutils (see the setup.cfg for more details).
#
# The next line, which defines namespace_packages and packages is only
# necessary if this projet contains a package belonging to the stsci namespace
# package, such as stsci.distutils or stsci.tools.  This is not necessary for
# projects with their own namespace, such as acstools or pyraf.
#
# d2to1=True is required to enable d2to1 and read the remaning project metadata
# from the setup.cfg file.
#
# use_2to3 and zip_safe are common options support by setuptools; these can
# also be placed in the setup.cfg, as will be demonstrated in a future update
# to this sample package.

if sys.version_info < (3,0): 
    pandas_version = '==0.24.0'
else:
    pandas_version = '>=0.25.0'
setup(

    name='pandexo.engine',
    version='2.0.0',
    summary='pandexo transiting exoplanet simulator',
    description_file='README.rst',
    author='Natasha Batalha at NASA Ames',
    author_email='natasha.e.batalha@nasa.gov',
    home_page='https://natashabatalha.github.io/PandExo',
    license='GPL',
    url='https://github.com/natashabatalha/PandExo',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    packages=[
        'pandexo',
        'pandexo.engine',
        'pandexo.engine.reference',
        'pandexo.engine.static',
        'pandexo.engine.static.css',
        'pandexo.engine.static.fonts',
        'pandexo.engine.static.img',
        'pandexo.engine.static.js',
        'pandexo.engine.temp',
        'pandexo.engine.templates',
        'pandexo.engine.utils'
    ],
    package_dir={
        'engine': 'pandexo/engine'
    },
    include_package_data=True,
    # package_data = {
    #            'engine' : ['templates','*.html'],
    #            'engine': ['reference','*.json'],
    #            'engine': ['static','css','*'],
    #            'engine': ['static','fonts','*'],
    #            'engine': ['static','img','*'],
    #            'engine': ['static','js','*']
    # },

    install_requires=[
          'pandeia.engine==2.0',
          'numpy',
          'bokeh==3.0.2',
          'tornado',
          'pandas'+pandas_version,
          'joblib',
          'photutils',
          'astropy',
          'pysynphot',
          'sqlalchemy',
          'astroquery',
          'scipy',
          'batman-package'
          ],
    entry_points = {
        'console_scripts':
            ['start_pandexo=pandexo.engine.run_online:main']
    },
    zip_safe=False,

)
