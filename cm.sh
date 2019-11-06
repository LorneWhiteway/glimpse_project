#!/bin/bash -i

rm -rf ./build
mkdir build
cd build

# Ignore what follows...
# Note for nfft we use /.libs, not /libs; this is a feature of the locally-built copy of nfft; if we switch back to a module version
# of nfft then probably we will need to go back to using /libs.

cmake .. -DCFITSIO_ROOT_DIR=${FITSIO_DIR} -DFFTW_ROOT=${FFTW_ROOT} -DNFFT_INCLUDE_DIR=${NFFT_HOME}/include -DNFFT_LIBRARY=${NFFT_HOME}/lib/libnfft3.a -DCCFITS_ROOT_DIR=${CCFITS_ROOT_DIR} -DARMADILLO_INCLUDE_DIR=${ARMADILLO_HOME}/include -DARMADILLO_LIBRARY=${ARMADILLO_HOME}/lib64/libarmadillo.so -DCMAKE_BUILD_TYPE=Release

#-DUSE_FFTW=ON -DCOMPILE_NFFT=ON




