# glimpse_project

Information useful for running glimpse on the splinter cluster at UCL.

## To build Glimpse in the splinter environment

1. git clone https://github.com/LorneWhiteway/glimpse_project.git
2. cd ./glimpse_project
3. Follow the instructions below to build a version of nfft that runs on the compute nodes.
4. source ./set_environment
5. git clone https://github.com/CosmoStat/Glimpse.git
6. Edit ./Glimpse/src/field.cpp to fix a bug: Find the line 
lensKernel[ind] = 1. ;
and in the same for loop also add the line
lensKernelTrue[ind] = 1. ;
7. cd ./Glimpse
8. ../cm.sh (This will run cmake with the correct settings)
9. cd ./build
10. make

## To build a version of nfft that will run on the compute nodes
The crucial point here is to force the use of -march=x86-64 when building.

1. Start in glimpse_project
2. git clone https://github.com/NFFT/nfft.git
3. cd ./nfft
4. bash
5. ./bootstrap.sh
6. ./configure --enable-all --enable-openmp --with-gcc-arch=x86-64
7. make

## Things to know about glimpse:
1. Command line arguments should be ini file, input catalogue, and output file name (in that order), The flags '--config', '--data' and '--output' should not be provided.
2. In the ini file, the 'size' argument is a diameter, not a radius. 