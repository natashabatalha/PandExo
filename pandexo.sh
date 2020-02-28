#!/bin/bash

#Another manual installation example


USRDIR=$HOME
echo 'USRDIR=$HOME' >>~/.bash_profile

########## BLOCK TO SET STELLAR DATA################
mkdir $USRDIR/pysynphot_data 
mv synphot5.tar* $USRDIR/pysynphot_data 
tar -xvf $USRDIR/pysynphot_data/synphot5.tar || tar -xvf $USRDIR/pysynphot_data/synphot5.tar.gz
echo 'export PYSYN_CDBS="$USRDIR/grp/hst/cdbs/grid"' >>~/.bash_profile

############ BLOCK TO SET PANDEIA REFERENCE DATA #########################

tar -xvf pandeia_data-1.3.tar.gz || tar -xvf pandeia_data-1.3.tar
mv pandeia_data $USRDIR
echo 'export pandeia_refdata="$USRDIR/pandeia_data"' >>~/.bash_profile

##############Make sure your bash profile is sourced
source ~/.bash_profile

########### Configure Conda to find STScI packages ###########
conda config --add channels http://ssb.stsci.edu/astroconda

########### Instal conda/stsci specific packages
conda install pyfftw
conda install numpy synphot joblib scipy astropy pyfftw pysynphot photutils sphinx
conda install bokeh

######### Install PandExo Engine, which will also install Pandeia##############
pip install pandexo.engine 

######### Try multiprocessing for 2.7 users##########
conda install multiprocessing || echo 'OKAY that Multiprocessing not found. Python 3 user'

########## RUN TEST ###############
python run_test.py
