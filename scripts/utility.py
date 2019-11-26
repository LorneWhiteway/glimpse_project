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
    
    nside = 512
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
    filename = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.nside" + str(nside) + "_truth.dat"
    
    hp.write_map(filename, average_values, do_nest, overwrite=True)
    
    
def compare_two_healpix_maps():

    import healpy as hp
    import matplotlib.pyplot as plt
    
    nside = 512

    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    filenames = ["Buzzard_192.nside" + str(nside) + "_merged_B_values.dat", "Buzzard_192.nside" + str(nside) + "_truth.dat"]
    
    show_weights = True
    if show_weights:
        filenames.append("Buzzard_192.nside" + str(nside) + "_merged_B_weights.dat")
    
    maps = []
    titles = []
    for f in filenames:
        maps.append(hp.read_map(path + f))
        titles.append(f)
        
    smooth_truth = True
    if smooth_truth:
        one_arcmin_in_radians = 0.000290888
        smoothing_scale_in_arcmin = 7.0
        maps[1] = hp.smoothing(maps[1], sigma = smoothing_scale_in_arcmin * one_arcmin_in_radians)
        titles[1] += " smoothed at {} arcmin".format(smoothing_scale_in_arcmin)
    
    show_diff = False
    if show_diff:
        maps.append(maps[0]-maps[1])
        titles.append("Diff")
    
    show_weights = True
        
    
    
    for map, title, i in zip(maps, titles, range(len(maps))):
        #hp.mollview(map, fig=i, title=title)
        hp.gnomview(map, fig=i, title=title, reso=2.0)
        hp.graticule()
    
    plt.show()
    

def clean_up_edges():

    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    
    do_nest = False
    
    nside = 512

    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    old_filename = "Buzzard_192.nside" + str(nside) + "_merged_B_values.dat"
    new_filename = "Buzzard_192.nside" + str(nside) + "_merged_B_values_masked.dat"
    new_filename_png = "Buzzard_192.nside" + str(nside) + "_merged_B_values_masked.png"
    
    m = hp.read_map(path + old_filename)
    (ra, dec) = get_fits_data(buzzard_data_file_name(), ["RA", "DEC"])
    
    # Create a mask that is zero in healpixels where there are no source galaxies and one elsewhere.
    
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
    



def sphere_to_tangent_plane_mapping(ra_centre, dec_centre, ra, dec):
    import numpy as np
    
    # See my notes p. GL90
    
    ra_centre_rad = np.radians(ra_centre)
    dec_centre_rad = np.radians(dec_centre)
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    
    x = np.sin(ra_rad - ra_centre_rad) / (np.cos(dec_centre_rad) * np.cos(ra_rad - ra_centre_rad) + np.tan(dec_rad) * np.sin(dec_centre_rad))
    
    epsilon = np.sin(dec_rad) * np.cos(dec_centre_rad) - np.cos(dec_rad) * np.sin(dec_centre_rad) * np.cos(ra_rad - ra_centre_rad)
    
    y = epsilon * np.sqrt((1.0 + x**2) / (1.0 - epsilon**2))
    
    return (x, y)
    

    
    
def sphere_to_tangent_plane_mapping_test_harness():

    import numpy as np
    
    pixel_id = 2218
    ra_centre = 56.25
    dec_centre = -27.27961273597809

    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192." + str(pixel_id) + ".glimpse.out.fits"
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    array_size = int(np.sqrt(len(ra)))
    
    
    if True:
    
        (x, y) = sphere_to_tangent_plane_mapping(ra_centre, dec_centre, ra, dec)
    
        x_new = np.reshape(x, [array_size, array_size])
        y_new = np.reshape(y, [array_size, array_size])

        step_size = x_new[0, 1] - x_new[0, 0]
        
        max_error = 0.0
        for i in range(array_size):
            for j in range(array_size):
                error_x = x_new[i][j] - step_size*(j - array_size/2 + 0.5)
                error_y = y_new[i][j] - step_size*(i - array_size/2 + 0.5)
                error = np.sqrt(error_x**2 + error_y**2)
                max_error = max(max_error, error)
                
                
        print("Max error = {}".format(max_error))
        
    if True:
    
        ra_new = np.reshape(ra, [array_size, array_size])
        dec_new = np.reshape(dec, [array_size, array_size])
        
        for i in [0, array_size-1]:
            for j in [0, array_size-1]:
                print("RA = {}; DEC = {}".format(ra_new[i][j], dec_new[i][j]))
        
    
    

    
            
    
    
    


if __name__ == '__main__':

    #print_list_of_healpixels()
    #kappa_values_in_one_fine_pixel()
    #tester()
    #save_buzzard_truth()
    #compare_two_healpix_maps()
    #clean_up_edges()
    sphere_to_tangent_plane_mapping_test_harness()
