#!/usr/bin/env python

# Copied from view_glimpse, and amended
def get_from_glimpse(glimpse_output_file):
    import numpy as np
    import astropy.io.fits as pyfits
    ra = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=1))
    dec = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=2))
    kappa = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=0))
    return (ra, dec, kappa)



# Output is sorted
def healpixel_ids_for_jobs():
    import glob
    filelist = glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.ini")
    filelist.sort()
    return [int(f.split(".")[1]) for f in filelist]
    


def print_list_of_healpixels():
    ids = healpixel_ids_for_jobs()
    for (i, id) in zip(range(len(ids)), ids):
        print(i, id)



def kappa_values_in_one_fine_pixel():

    import healpy as hp
    import glob
    import numpy as np
    import sys
    
    (r,d) = hp.pix2ang(16, 2218, False, True)
    specific_fine_healpix_id = hp.ang2pix(1024, r, d, False, True)
    
    
    j = 0
    for f in sorted(glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.out.fits")):
        (ra, dec, kappa) = get_from_glimpse(f)
        range_k = np.arange(len(kappa))
        fine_healpix_id = hp.pixelfunc.ang2pix(1024, ra, dec, lonlat = True)
        filter = np.where(fine_healpix_id == specific_fine_healpix_id)
        num_in_fine_healpixel = len(ra[filter])
        if num_in_fine_healpixel > 0:
            coarse_id = int(f.split(".")[1])
            for (r,d,k,index,filter_val) in zip(ra[filter], dec[filter], kappa[filter], range(num_in_fine_healpixel), range_k[filter]):
                print(j, coarse_id, r, d, k, index, filter_val)
                j += 1

        
def tester():
    import healpy as hp
    (ra ,dec) = hp.pix2ang(16, 2218, False, True)
    id = hp.ang2pix(1024, ra, dec, False, True)
    print(ra,dec,id)
    
#"RA", "DEC", "true_z", "E1", "E2", "G1", "G2", "k_orig"
def buzzard_data_file_name():
    return '/share/testde/ucapnje/buzzard_desy3_marco/Buzzard_192.fits'


# 'what_to_get' should be a list of field names (case insensitive)
def get_fits_data(file_name, what_to_get):
    import astropy.io.fits as pyfits
    print("Opening input file " + file_name + "...")
    res = []
    if what_to_get: # i.e. if list is not empty
        x = pyfits.open(file_name)
        for w in what_to_get:
            print("Getting field {}".format(w))
            res.append(x[1].data.field(w))
    print("Closing input file " + file_name + ".")
    return res
    



def save_buzzard_truth():

    import healpy as hp
    import numpy as np
    
    do_nest = False
    
    nside = 1024
    num_healpixels = hp.nside2npix(nside)
    sums = np.zeros(num_healpixels)
    count = np.zeros(num_healpixels)
    
    (ra, dec, kappa) = get_fits_data(buzzard_data_file_name(), ["RA", "DEC", "k_orig"])
    num_buzzard_data = len(ra)
    print("num_buzzard_data = {}".format(num_buzzard_data))
    print("About to call ang2pix...")
    id = hp.ang2pix(nside, ra, dec, do_nest, True)
    print("Finished call to ang2pix.")
    i = 0
    for (this_id, this_kappa) in zip(id, kappa):
        sums[this_id] += this_kappa
        count[this_id] += 1
        i += 1
        if i % 1000000 == 0:
            print("{} of {}".format(i, num_buzzard_data))
        
        
    # The following division code treats 0/0 as 0. From https://stackoverflow.com/questions/26248654.
    average_values = np.divide(sums, count, out=np.zeros_like(count, dtype=float), where=count!=0.0)
    
    # Save
    filename = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.nside1024_truth.dat"
    
    hp.write_map(filename, average_values, do_nest, overwrite=True)
    
    
def compare_two_healpix_maps():

    import healpy as hp
    import matplotlib.pyplot as plt

    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    filenames = ["Buzzard_192.nside1024_merged_values.dat", "Buzzard_192.nside1024_truth.dat"]
    maps = []
    titles = []
    for f in filenames:
        maps.append(hp.read_map(path + f))
        titles.append(f)
        
    smooth_truth = True
    if smooth_truth:
        one_arcmin_in_radians = 0.000290888
        smoothing_scale_in_arcmin = 15.0
        maps[1] = hp.smoothing(maps[1], sigma = smoothing_scale_in_arcmin * one_arcmin_in_radians)
        titles[1] += " smoothed at {} arcmin".format(smoothing_scale_in_arcmin)
    
    maps.append(maps[0]-maps[1])
    titles.append("Diff")
    
    for map, title, i in zip(maps, titles, range(len(maps))):
        #hp.mollview(map, fig=i, title=title)
        hp.gnomview(map, fig=i, title=title)
        hp.graticule()
    
    plt.show()
    

def clean_up_edges():

    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    
    do_nest = False

    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    old_filename = "Buzzard_192.nside1024_merged_values.dat"
    new_filename = "Buzzard_192.nside1024_merged_values_masked.dat"
    new_filename_png = "Buzzard_192.nside1024_merged_values_masked.png"
    
    m = hp.read_map(path + old_filename)
    (ra, dec) = get_fits_data(buzzard_data_file_name(), ["RA", "DEC"])
    
    # Create a mask that is zero in healpixels where there are no source galaxies and one elsewhere.
    nside = 1024
    ids = hp.ang2pix(nside, ra, dec, do_nest, True)
    mask = np.zeros(hp.nside2npix(nside))
    for id in ids:
        mask[id] = 1.0 
    
    # Apply the mask
    m *= mask
    
    hp.write_map(path + new_filename, m, overwrite=True)
    
    # Also save as png
    hp.mollview(m, title=new_filename)
    plt.savefig(path + new_filename_png)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


if __name__ == '__main__':

    #print_list_of_healpixels()
    #kappa_values_in_one_fine_pixel()
    #tester()
    #save_buzzard_truth()
    #compare_two_healpix_maps()
    clean_up_edges()
