# Only need to change these two variables
#https://gist.github.com/zshaheen/fe76d1507839ed6fbfbccef6b9c13ed9
PKG_NAME=pandexo.engine
USER=natashabatalha

OS=$TRAVIS_OS_NAME-64
mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`date +%Y.%m.%d`
conda build .
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER -l nightly $CONDA_BLD_PATH/$OS/$PKG_NAME-`date +%Y.%m.%d`-0.tar.bz2 --force