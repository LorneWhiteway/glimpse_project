#!/usr/bin/env python


def main(job_id):

    import glob
    import subprocess
    
    output_directory = "/share/splinter/ucapwhi/glimpse_project/output/"
    project = "Mcal_0.2_1.3"
    exe_file = "/share/splinter/ucapwhi/glimpse_project/Glimpse/build/glimpse"
    
    glimpse_cat_file_pattern = output_directory + project + ".*.glimpse.cat.fits"
    id_list = [int(f.split(".")[1]) for f in glob.glob(glimpse_cat_file_pattern)]
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
    if len(sys.argv) <= 1:
        raise RuntimeError("Usage: glimpse_caller.py job_id")
    main(int(sys.argv[1]))