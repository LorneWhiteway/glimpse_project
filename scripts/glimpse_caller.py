#!/usr/bin/env python


def main(job_id):

    import glob
    import subprocess
    
    glimpse_cat_file_pattern = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.cat.fits"
    id_list = [int(f.split(".")[1]) for f in glob.glob(glimpse_cat_file_pattern)]
    id_list.sort()
    this_healpix_pixel_id = id_list[job_id]
    
    print("This healpix pixel id is " + str(this_healpix_pixel_id) + "...")
    
    root = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192." + str(this_healpix_pixel_id).zfill(4)
    
    ini_file = root + ".glimpse.ini"
    cat_file = root + ".glimpse.cat.fits"
    out_file = root + ".glimpse.out.fits"
    
    exe_file = "/share/splinter/ucapwhi/glimpse_project/Glimpse/build/glimpse"
    
    print("About to call " + exe_file)
    print("with options")
    print(ini_file)
    print(cat_file)
    print(out_file)
    run_info = subprocess.run([exe_file, ini_file, cat_file, out_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(run_info.stdout)
    print("Finished.")
    
    
    


if __name__ == '__main__':
    import sys
    if len(sys.argv) <= 1:
        raise RuntimeError("Usage: glimpse_caller.py job_id")
    main(int(sys.argv[1]))