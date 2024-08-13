#!/bin/bash

# For Glimpse (both for building and for running - for running perhaps only a subset of these is needed - haven't checked)
COMPILER_MODULE=dev_tools/sep2019/gcc-7.4.0
CMAKE_MODULE=dev_tools/sep2019/cmake-3.15.3
BOOST_MODULE=dev_tools/may2017/boost-1.64.0
GSL_MODULE=science/sep2019/gsl-2.6
FFTW_MODULE=fft/may2017/fftw-3.3.6
NFFT_MODULE=fft/nov2019/nfft-3.5.1
CFITSIO_MODULE=astro/apr2017/cfitsio3390
CCFITS_MODULE=astro/sep2019/CCfits-2.5
ARMADILLO_MODULE=linear_algebra/may2017/armadillo-7.900.1

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
