#!/usr/bin/env bash

../scripts/glimpse_on_curved_sky.py -d . -t create_cutouts
sbatch --wait --array=0-16 ./job_script_glimpse_caller
wait
../scripts/glimpse_on_curved_sky.py -d . -t merge
../scripts/compare_regression_results.py ./expected.glimpse.merged.values.dat ./glimpse.merged.values.dat
../scripts/compare_regression_results.py ./expected.glimpse.merged.weights.dat ./glimpse.merged.weights.dat