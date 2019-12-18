#!/usr/bin/env python
# -*- coding: utf-8 -*-


    

def run(ini_filename, task):

    from pathlib import Path
    import utility
    
    f = Path(ini_filename)
    assert f.is_file(), "ini file {} not found".format(ini_filename)

    tasks = ["create_cutouts", "run_glimpse", "merge"]
    assert (task in tasks), "task should be one of create_cutouts, run_glimpse or merge"
    
    if task == tasks[0]:
        utility.create_cutouts(ini_filename)
    elif task == tasks[1]:
        utility.run_glimpse(ini_filename)
    elif task == tasks[2]:
        utility.merge(ini_filename)
    
    


if __name__ == '__main__':

    import argparse
    import sys
    
    try:

        
        parser = argparse.ArgumentParser(description = "Given a weak-lensing catalogue, creates a set of 'cutout' catalogues given one input weak-lensing catalogue.")

        parser.add_argument('-i', '--ini_file', type = str, required = True, help = "Input ini file name.")
        parser.add_argument('-t', '--task', type = str, required = True, help = "Task to be performed; one of create_cutouts, run_glimpse or merge.")
        
        args = parser.parse_args()
        
        run(args.ini_file, args.task)
        
    except Exception as err:
        print('Error: {0}'.format(str(err)))
        sys.exit(1)
        
