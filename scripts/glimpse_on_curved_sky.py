#!/usr/bin/env python
# -*- coding: utf-8 -*-


    

def run(ini_filename, task, job_control):

    from pathlib import Path
    import utility
    
    f = Path(ini_filename)
    assert f.is_file(), "ini file {} not found".format(ini_filename)

    tasks = ["create_cutouts", "run_glimpse", "merge"]
    assert (task in tasks), "task should be one of create_cutouts, run_glimpse or merge"
    
    if task == tasks[0]:
        utility.create_cutouts_caller(ini_filename, job_control)
    elif task == tasks[1]:
        utility.run_glimpse(ini_filename)
    elif task == tasks[2]:
        utility.merge_caller(ini_filename, job_control)
    
    


if __name__ == '__main__':

    import argparse
    import sys
    import traceback
    
    try:

        
        parser = argparse.ArgumentParser(description = "Given a weak-lensing catalogue, creates a set of 'cutout' catalogues given one input weak-lensing catalogue.")

        parser.add_argument('-i', '--ini_file', type = str, required = True, help = "Input ini file name.")
        parser.add_argument('-t', '--task', type = str, required = True, help = "Task to be performed; one of create_cutouts, run_glimpse or merge.")
        parser.add_argument('-j', '--job_control', type = str, required = False, default = "::", help = "Numpy slice for which healpixels to process e.g. 2:10:2; default is :: i.e. all. Enclose in double quotes.")
        
        args = parser.parse_args()
        
        run(args.ini_file, args.task, args.job_control)
        
    except Exception as err:
        print(traceback.format_exc())
        sys.exit(1)
        
