# glimpse_project

Information useful for running glimpse on the splinter cluster at UCL.

## To build Glimpse in the splinter environment

1. git clone https://github.com/LorneWhiteway/glimpse_project.git
2. cd ./glimpse_project
3. source ./set_environment
4. git clone https://github.com/CosmoStat/Glimpse.git
5. Edit ./Glimpse/src/field.cpp to fix a bug: Find the line 
lensKernel[ind] = 1. ;
and in the same for loop also add the line
lensKernelTrue[ind] = 1. ;
6. cd ./Glimpse
7. ../cm.sh (This will run cmake with the correct settings)
8. cd ./build
9. make

## Things to know about glimpse:
1. Command line arguments should be ini file, input catalogue, and output file name (in that order), The flags '--config', '--data' and '--output' should not be provided.
2. In the ini file, the 'size' argument is a diameter, not a radius. 