#!/bin/csh

set COMPILER_MODULE            = dev_tools/apr2017/gcc-5.2.0
set CMAKE_MODULE               = dev_tools/oct2016/cmake-3.7.0-rc2
set BOOST_MODULE               = dev_tools/oct2016/boost-1.62.0
set GSL_MODULE                 = science/nov2014/gsl-1.16
set FFTW_MODULE                = fft/may2017/fftw-3.3.6
set NFFT_MODULE                = fft/may2017/nfft-3.3.2
set CFITSIO_MODULE             = astro/oct2016/cfitsio3390
set CCFITS_MODULE              = astro/may2017/CCfits-2.5
set ARMADILLO_MODULE           = linear_algebra/may2017/armadillo-7.900.1

module purge
module load $COMPILER_MODULE
module load $CMAKE_MODULE
module load $BOOST_MODULE
module load $GSL_MODULE
module load $FFTW_MODULE
module load $NFFT_MODULE
module load $CFITSIO_MODULE
module load $CCFITS_MODULE
module load $ARMADILLO_MODULE

setenv FFTWDIR ${FFTW_ROOT}

