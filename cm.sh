#!/bin/bash -i

rm -rf ./build
mkdir build
cd build

# Note that the splinter module file for CCfits doesn't (yet) define an environment variable for the package location.
cmake .. -DCFITSIO_ROOT_DIR=${FITSIO_DIR} -DFFTW_ROOT=${FFTW_ROOT} -DNFFT_INCLUDE_DIR=${NFFT_HOME}/include -DNFFT_LIBRARY=${NFFT_HOME}/lib/libnfft3.a -DCCFITS_ROOT_DIR=/share/splinter/cosmos/modules/may_2017/install_dir/CCfits-2.5/ -DARMADILLO_INCLUDE_DIR=${ARMADILLO_HOME}/include -DARMADILLO_LIBRARY=${ARMADILLO_HOME}/lib64/libarmadillo.so

#-DUSE_FFTW=ON -DCOMPILE_NFFT=ON




