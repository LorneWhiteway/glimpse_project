#!/bin/csh

# For Glimpse (both for building and for running - for running perhaps only a subset of these is needed - haven't checked)
set COMPILER_MODULE            = dev_tools/sep2019/gcc-7.4.0
set CMAKE_MODULE               = dev_tools/sep2019/cmake-3.15.3
set BOOST_MODULE               = dev_tools/may2017/boost-1.64.0
set GSL_MODULE                 = science/sep2019/gsl-2.6
set FFTW_MODULE                = fft/may2017/fftw-3.3.6
set NFFT_MODULE                = fft/may2017/nfft-3.3.2
set CFITSIO_MODULE             = astro/apr2017/cfitsio3390
set CCFITS_MODULE              = astro/sep2019/CCfits-2.5
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

# Load Python, if desired
set PYTHON_MODULE              = dev_tools/oct2018/python-Anaconda-3-5.3.0
module load $PYTHON_MODULE

# Now set further environment variables as needed by glimpse's cmake.
# Actually I am not sure that these are needed. See cm.sh, which should be used to call cmake.
#setenv FFTWDIR ${FFTW_ROOT}
#setenv CFITSIO_ROOT_DIR ${CFITSIO}

