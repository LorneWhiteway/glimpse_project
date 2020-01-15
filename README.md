# glimpse_project

## Introduction

This repository contains code to allow [glimpse](https://github.com/CosmoStat/Glimpse) (a mass-mapping algorithm that works on patches of the sky that are assumed to be small and flat) to be run on a larger patch (such as the DES footprint).

The code does this in three steps:
a. *create_cutouts*: an input weak-lensing catalogue is subdivided into many smaller overlapping sub-catalogues (call each one a 'cutout');
b. *run_glimpse*: glimpse is run (possibly in parallel) on each cutout;
c. *merge*: the separate glimpse results are merged together to create an output convergence map (in Healpix format).

The code is structured so that other mass-mapping algorithms can be handled in the future.

## Interface

The primary interface is the Python script `glimpse_on_curved_sky.py` in the `./scripts` directory. Run this with command line arguments as follows:

| Option | Meaning |
| --- | --- |
| -i | Name of configuration file. See below for details of configuration file format. Required. |
| -t | Task. One of "create_cutouts", "run_glimpse" or "merge". Required. |
| -j | Job control. See below for more information. Optional. |

Example:
```
./scripts/glimpse_on_curved_sky.py -i ./output/my_job.ini -t create_cutouts
```

### Configuration file

(TO DO)

### Job Control

(TO DO)

## Information useful for running glimpse on the splinter cluster at UCL.

### To build Glimpse in the splinter environment

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

### Things to know about glimpse:
1. Command line arguments should be ini file, input catalogue, and output file name (in that order), The flags '--config', '--data' and '--output' should not be provided.
2. In the ini file, the 'size' argument is the length of the edge of a square. 