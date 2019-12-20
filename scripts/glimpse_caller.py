#!/usr/bin/env python


def glimpse_caller(ini_file_name, job_id):

    import glob
    import subprocess
    import configparser
    
    config = configparser.ConfigParser()
    config.read(ini_file_name)
    
    output_directory = config["project"].get("directory")
    project = config["project"].get("project_name")
    exe_file = config["project"].get("glimpse_executable")
    
    glimpse_cat_file_pattern = output_directory + project + ".*.glimpse.cat.fits"
    id_list = [int(f.split(".")[-4]) for f in glob.glob(glimpse_cat_file_pattern)]
    id_list.sort()
    this_healpix_pixel_id = id_list[job_id]
    
    print("Healpixel id = {}...".format(this_healpix_pixel_id))
    
    this_healpix_id_as_string = str(this_healpix_pixel_id).zfill(4)
    
    ini_file = output_directory + project  + ".glimpse.ini"
    cat_file = output_directory + project + "." + this_healpix_id_as_string + ".glimpse.cat.fits"
    out_file = output_directory + project + "." + this_healpix_id_as_string + ".glimpse.out.fits"

    subprocess.run([exe_file, ini_file, cat_file, out_file])

    print("Finished.")
    
    
    


if __name__ == '__main__':
    import sys
    if len(sys.argv) <= 2:
        raise RuntimeError("Usage: glimpse_caller.py ini_file job_id")
    main(sys.argv[1], int(sys.argv[2]))