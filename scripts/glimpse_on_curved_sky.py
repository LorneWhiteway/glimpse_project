#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
    Interface file for for 'glimpse on curved sky' project.
    All of the actual functionality is located in utilities.py.
    Author: Lorne Whiteway.
"""


import utility
import argparse
import sys
import traceback
    

def run(directory, task, job_control):


    tasks = ["create_cutouts", "run_glimpse", "merge", "status"]
    assert (task in tasks), "task should be one of create_cutouts, run_glimpse, merge or status"
    
    if task == tasks[0]:
        utility.create_cutouts_caller(directory, job_control)
    elif task == tasks[1]:
        utility.glimpse_caller(directory, job_control)
    elif task == tasks[2]:
        utility.merge_caller(directory, job_control)
    elif task == tasks[3]:
        utility.status_caller(directory)
    
    


if __name__ == '__main__':

    try:
        
        parser = argparse.ArgumentParser(description = "Top-level routine for running glimpse on curved sky; see README.md for details")

        parser.add_argument('-d', '--directory', type = str, required = True, help = "Directory for input and output files.")
        parser.add_argument('-t', '--task', type = str, required = True, help = "Task to be performed; one of create_cutouts, run_glimpse, merge or status.")
        parser.add_argument('-j', '--job_control', type = str, required = False, default = "", help = "Numpy slice for which healpixels to process e.g. 2:10:2; default is empty string i.e. all. Enclose in double quotes. But for run_glimpse this should just be the job number.")
        
        args = parser.parse_args()
        
        run(args.directory, args.task, args.job_control)
        
    except Exception as err:
        print(traceback.format_exc())
        sys.exit(1)
        
