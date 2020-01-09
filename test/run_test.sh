#!/usr/bin/env bash

 ../scripts/glimpse_on_curved_sky.py -i ./glimpse.ini -t create_cutouts
 sbatch --wait --array=0-16 ./job_script_glimpse_caller
 wait
 ../scripts/glimpse_on_curved_sky.py -i ./glimpse.ini -t merge
 