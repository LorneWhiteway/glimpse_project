#!/usr/bin/env python

# Expecting source_file_name to be of the form .../Buzzard_192.xxxx.glimpse.cat.fits
# where xxxx is a healpix pixel id.
def create_one_file(source_file_name, nside, nest=False):

    import healpy as hp

    healpix_id = int(source_file_name.split(".")[1])
    (centre_ra, centre_dec) = hp.pix2ang(nside, healpix_id, nest, True)
    ini_filename = source_file_name.split(".cat.fits")[0] + ".ini"
    
    
    with open(ini_filename, 'w') as f_out:
        f_out.write("[survey]\n")
        f_out.write("center_ra=" + str(centre_ra) + "\n")
        f_out.write("center_dec=" + str(centre_dec) + "\n")
        f_out.write("size=16.0\n")
        f_out.write("units=degrees\n")
        f_out.write("hdu=1\n")
        f_out.write("flip_e2=true\n")
        f_out.write("ra=RA\n")
        f_out.write("dec=DEC\n")
        f_out.write("e1=E1\n")
        f_out.write("e2=E2\n")
        f_out.write("z=true_z\n")
        f_out.write("[cosmology]\n")
        f_out.write("Omega_m=0.25\n")
        f_out.write("h=0.70\n")
        f_out.write("[field]\n")
        f_out.write("units=arcmin\n")
        f_out.write("pixel_size=3.5\n")
        f_out.write("padding=28\n")
        f_out.write("include_flexion=false\n")
        f_out.write("zlens=-1\n")
        f_out.write("[parameters]\n")
        f_out.write("nrandom=1000\n")
        f_out.write("niter=500\n")
        f_out.write("nreweights=0\n")
        f_out.write("niter_debias=0\n")
        f_out.write("nscales=7\n")
        f_out.write("lambda=3.0\n")
        f_out.write("battle_lemarie_reg=1.\n")
        f_out.write("last_scale_reg=2.\n")

def run_create_ini_files():
    import glob
    
    list_of_augmented_healpix_files = glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.cat.fits")
    
    for f in list_of_augmented_healpix_files:
        print("Working on " + f + "...")
        create_one_file(f, 16)
        
        



if __name__ == '__main__':
    run_create_ini_files()
