#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################## Start of one-off utilities ##########################


def compare_two_cutouts():

    import matplotlib.pyplot as plt

    (ra_1, dec_1, e1_1, e2_1) = get_from_fits_file("/share/splinter/ucapwhi/glimpse_project/output/Mcal_0.2_1.3.2218.glimpse.cat.fits", ["RA", "DEC", "E1", "E2"])
    (ra_2, dec_2, e1_2, e2_2) = get_from_fits_file("/share/splinter/ucapwhi/glimpse_project/output_Mcal_badsign/Mcal_0.2_1.3.2218.glimpse.cat.fits", ["RA", "DEC", "E1", "E2"])
    
    for i in range(3159, 3170):
        print(i, ra_1[i], ra_2[i], dec_1[i], dec_2[i])
        



def correct_one_shear_catalogue(catalogue_filename):
    
    from astropy.table import Table
    import astropy.io.fits as pyfits
    
    print("Processing {}...".format(catalogue_filename))
    
    #list_of_field_names = ["RA", "DEC", "E1", "E2", "true_z"]
    list_of_field_names = ["ra_gal", "dec_gal", "e1_gal", "e2_gal", "z"]
    list_of_data_columns = get_from_fits_file(catalogue_filename, list_of_field_names)
        
    list_of_data_columns[2] *= -1.0 # e1_gal
    list_of_data_columns[3] *= -1.0 # e2_gal
    
    #output_filename = catalogue_filename.replace("/output/", "/output1/")
    output_filename = catalogue_filename.replace(".fits", ".negative.fits")
    
    write_to_fits_file(output_filename, list_of_field_names, list_of_data_columns)
    
    
def glimpse_sign_experiment():

    import matplotlib.pyplot as plt

    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/output.fits"
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    
    glimpse_output_file_negative = "/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/output.negative.fits"
    (ra_negative, dec_negative, kappa_negative) = get_from_glimpse(glimpse_output_file_negative)
    
    plt.scatter(kappa, kappa_negative, s=1)
    plt.show()
    
    

def correct_one_shear_catalogue_caller():
    
    import glob
    filelist = glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.cat.fits")
    for f in filelist: 
        correct_one_shear_catalogue(f)
    

    
    

def num_files_in_directory():
    import glob
    print(len(glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Mcal_0.2_1.3.*.glimpse.cat.fits")))



# Written to produce data for Niall; used on 9 Dec 2019
def redshift_histogram():
    
    import numpy as np
    np.set_printoptions(precision=3)
    
    top_z = 2.5
    
    (ra, dec, true_z) = get_from_fits_file(buzzard_data_file_name(), ["RA", "DEC", "true_z"])
    (hist, bin_edges) = np.histogram(true_z, np.linspace(0.0, top_z, int(100*top_z) + 1))
    for (b0, b1, h) in zip(bin_edges[:-1], bin_edges[1:], hist):
        print("{0:.2f}\t{1:.2f}\t{2:d}".format(b0, b1, h))

    
def kappa_histogram():

    import matplotlib.pyplot as plt
    import healpy as hp
    import numpy as np
    
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    f = "Mcal_0.2_1.3_90_110_2048_downgraded_to_1024_masked.glimpse.merged.values.dat"
    m = hp.read_map(path + f)
    
    plt.hist(m[np.where(m != 0.0)], bins = 50)
    plt.show()
    
def show_glimpse_output_as_image():

    import numpy as np
    import math
    import matplotlib.pyplot as plt

    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/output.negative.fits"
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    
    n = int(math.sqrt(kappa.shape[0]))
    
    kappa_as_2d_array = np.reshape(kappa, [n, n])
    
    plt.imshow(kappa_as_2d_array)
    plt.show()
    
    

    
# Written to produce data for Niall; used on 10 Dec 2019
def shear_stdev():
    import numpy as np
    (e1, e2) = get_from_fits_file(buzzard_data_file_name(), ["E1", "E2"])
    e1_and_e2 = np.concatenate((e1, e2))
    print(np.std(e1_and_e2))



def tester():
    import healpy as hp
    (ra ,dec) = hp.pix2ang(16, 1242, False, True)
    id = hp.ang2pix(1024, ra, dec, False, lonlat=True)
    print(ra,dec,id)

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

def add_dummy_redshift_column(file_name):

    import numpy as np
    
    list_of_field_names = ["RA", "DEC", "E1", "E2", "E1_RANDOM", "E2_RANDOM"]
    list_of_data_columns = get_from_fits_file(file_name, list_of_field_names)
    
    dummy_z = np.ones(list_of_data_columns[0].shape[0])
    list_of_field_names.append("DUMMY_Z")
    list_of_data_columns.append(dummy_z)
    
    write_to_fits_file(file_name.replace(".fits", ".new.fits"), list_of_field_names, list_of_data_columns)
    

def save_buzzard_truth():

    import healpy as hp
    import numpy as np
    
    do_nest = False
    
    nside = 512
    num_healpixels = hp.nside2npix(nside)
    sums = np.zeros(num_healpixels)
    count = np.zeros(num_healpixels)
    
    (ra, dec, kappa) = get_from_fits_file(buzzard_data_file_name(), ["RA", "DEC", "k_orig"])
    num_buzzard_data = len(ra)
    print("num_buzzard_data = {}".format(num_buzzard_data))
    id = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
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
    

# One-time test; see p GL97.
def glimpse_array_order():
    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    
    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.1440.glimpse.out.fits"
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    
    plt.plot(ra[:1000], c='red')
    plt.plot(dec[:1000], c='blue')
    plt.show()

########################## End of one-off utilities ##########################


# Credit: Stephen Rauch at https://stackoverflow.com/questions/48494581
# Interprets a string such as '[2:45:3] as a numpy 'slice' object for array indexing.
def slice_from_string(slice_string):
    slices = slice_string.split(',')
    if len(slices) > 1:
        return [slice_from_string(s.strip()) for s in slices]
    return slice(*[int(x) for x in slice_string.split(':')])



# Copied from view_glimpse, and amended
def get_from_glimpse(glimpse_output_file):
    import numpy as np
    import astropy.io.fits as pyfits
    ra = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=1))
    dec = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=2))
    kappa = np.ndarray.flatten(pyfits.getdata(glimpse_output_file, ext=0))
    return (ra, dec, kappa)






def kappa_values_in_one_fine_pixel():

    import healpy as hp
    import glob
    import numpy as np
    import sys
    
    (r,d) = hp.pix2ang(16, 2218, False, True)
    specific_fine_healpix_id = hp.ang2pix(1024, r, d, False, lonlat=True)
    
    
    j = 0
    for f in sorted(glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.out.fits")):
        (ra, dec, kappa) = get_from_glimpse(f)
        range_k = np.arange(len(kappa))
        fine_healpix_id = hp.pixelfunc.ang2pix(1024, ra, dec, lonlat=True)
        filter = np.where(fine_healpix_id == specific_fine_healpix_id)
        num_in_fine_healpixel = len(ra[filter])
        if num_in_fine_healpixel > 0:
            coarse_id = int(f.split(".")[1])
            for (r,d,k,index,filter_val) in zip(ra[filter], dec[filter], kappa[filter], range(num_in_fine_healpixel), range_k[filter]):
                print(j, coarse_id, r, d, k, index, filter_val)
                j += 1

        
    
#"RA", "DEC", "true_z", "E1", "E2", "G1", "G2", "k_orig"
def buzzard_data_file_name():
    return '/share/testde/ucapnje/buzzard_desy3_marco/Buzzard_192.fits'

#"RA", "DEC", "E1", "E2", "E1_RANDOM", "E2_RANDOM", "DUMMY_Z"
def metacal_data_file_name():
    return '/share/testde/ucapnje/year3_des/Mcal_0.2_1.3.fits'
    
    
# 'what_to_get' should be a list of field names (case insensitive)
def get_from_fits_file(file_name, what_to_get):
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
    



def write_to_fits_file(output_filename, list_of_field_names, list_of_data_columns):
    import astropy.io.fits as pyfits
    print("Writing to {}...".format(output_filename))
    column_info = []
    for (field_name, data_column) in zip(list_of_field_names, list_of_data_columns):
        column_info.append(pyfits.Column(name=field_name, format='D', array=data_column))
    tbhdu = pyfits.BinTableHDU.from_columns(column_info)
    tbhdu.writeto(output_filename, overwrite=True)



def fits_catalog_to_healpix_map(catalog_file_name, nside, do_nest):

    import healpy as hp
    import numpy as np

    (ra, dec) = get_from_fits_file(catalog_file_name, ["ra", "dec"])
    
    num_healpixels = hp.nside2npix(nside)
    res = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for id in ids:
        res[id] += 1.0
    return res    
    
def glimpse_output_to_healpix_map(glimpse_output_file, nside, do_nest):

    import healpy as hp
    import numpy as np
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    
    num_healpixels = hp.nside2npix(nside)
    weighted_values = np.zeros(num_healpixels)
    weights = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for (id, k) in zip(ids, kappa):
        weighted_values[id] += k
        weights[id] += 1.0
    return np.divide(weighted_values, weights, out=np.zeros_like(weights, dtype=float), where=weights!=0.0)
    
def glimpse_lattice_points_to_healpix_map(glimpse_output_file, nside, do_nest):

    import healpy as hp
    import numpy as np
    (ra, dec, kappa) = get_from_glimpse(glimpse_output_file)
    
    num_healpixels = hp.nside2npix(nside)
    ret = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for id in ids:
        ret[id] = 1.0
    return ret

# Cleanup ra_rad so that it lies in [0, 2pi)
def cleanup_ra_rad(ra_rad):
    import numpy as np
    two_pi = 2.0 * np.pi
    ret = ra_rad
    # Note that (unexpectedly) the order of the next two lines is crucial
    # for proper handling of the case ra_rad = -epsilon.
    ret[np.where(ret < 0.0)] += two_pi
    ret[np.where(ret >= two_pi)] -= two_pi
    return ret


################
# Example code for timing
    #from timeit import default_timer as timer
    #start = timer()
    #print("step {}, {}".format(1, timer() - start))




# See LW's book p. 114.
# Direction = 1.0 for clockwise on the sky, -1.0 for anticlockwise on the sky
# Rotates about (0, 0)
def rotate_45_degrees(ra_rad, dec_rad, direction):
    import numpy as np

    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    sin_ra = np.sin(ra_rad)
    cos_ra = np.cos(ra_rad)

    ra_rotated_rad = np.arctan2((cos_dec * sin_ra - direction * sin_dec) * np.sqrt(2.0) / 2.0, cos_dec * cos_ra)
    dec_rotated_rad = np.arcsin((sin_dec + direction * cos_dec * sin_ra) * np.sqrt(2.0) / 2.0)

    return (ra_rotated_rad, dec_rotated_rad)
    

def rotate_shear_45_degrees(e1, e2):
    # See LW's notes p. 107 (and 122)
    return (-e2, e1)



def to_standard_position(ra, dec, ra_centre, dec_centre):
    import numpy as np

    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    ra_centre_rad = np.radians(ra_centre)
    dec_centre_rad = np.radians(dec_centre)

    sin_ra_diff = np.sin(ra_rad - ra_centre_rad)
    cos_ra_diff = np.cos(ra_rad - ra_centre_rad)
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
        
    
    # Move (ra_centre, dec_centre) to (0, 0)
    # See LW's book p. 78
    ra_shifted_rad = np.arctan2(cos_dec * sin_ra_diff, cos_dec * np.cos(dec_centre_rad) * cos_ra_diff + sin_dec * np.sin(dec_centre_rad))
    dec_shifted_rad = np.arcsin(sin_dec * np.cos(dec_centre_rad) - cos_dec * np.sin(dec_centre_rad) * cos_ra_diff)
       
    # Rotate by 45 degrees
    (ra_shifted_rotated_rad, dec_shifted_rotated_rad) = rotate_45_degrees(ra_shifted_rad, dec_shifted_rad, 1.0)
        
    # Move (0, 0) to (pi, 0) and cleanup RA to be in [0, 2pi).
    ra_shifted_rotated_rad += np.pi
    ra_shifted_rotated_rad = cleanup_ra_rad(ra_shifted_rotated_rad)
    
    ra_shifted_rotated = np.degrees(ra_shifted_rotated_rad)
    dec_shifted_rotated = np.degrees(dec_shifted_rotated_rad)
    
    
    return (ra_shifted_rotated, dec_shifted_rotated)
    
  

    
def from_standard_position(ra, dec, ra_centre, dec_centre):
    import numpy as np
    
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    ra_centre_rad = np.radians(ra_centre)
    dec_centre_rad = np.radians(dec_centre)

    # Move (pi, 0) to (0, 0)
    ra_rad += np.pi

    # Rotate by -45 degrees
    (ra_rotated_rad, dec_rotated_rad) = rotate_45_degrees(ra_rad, dec_rad, -1.0)
    
    # Move (0, 0) to (ra_centre, dec_centre). See LW's notes p. 85.
    ra_rotated_shifted_rad = ra_centre_rad + np.arctan2(np.cos(dec_rotated_rad) * np.sin(ra_rotated_rad), np.cos(dec_rotated_rad) * np.cos(dec_centre_rad) * np.cos(ra_rotated_rad) - np.sin(dec_rotated_rad) * np.sin(dec_centre_rad))
    dec_rotated_shifted_rad = np.arcsin(np.sin(dec_rotated_rad) * np.cos(dec_centre_rad) + np.cos(dec_rotated_rad) * np.sin(dec_centre_rad) * np.cos(ra_rotated_rad))
    
    # Cleanup RA into [0, 2pi)
    ra_rotated_shifted_rad = cleanup_ra_rad(ra_rotated_shifted_rad)
    
    ra_rotated_shifted = np.degrees(ra_rotated_shifted_rad)
    dec_rotated_shifted = np.degrees(dec_rotated_shifted_rad)
    
    return (ra_rotated_shifted, dec_rotated_shifted)
    

def to_standard_position_test_harness():

    import numpy as np
    import matplotlib.pyplot as plt
    
    ra_centre = 72.0
    dec_centre = -18.0
    
    if False:
        ra = np.linspace(0.0, 359.0, num = 360)
        dec = ra * 0 + 1.0
        theta = ra
    else:
    
        theta = np.linspace(0.0, 2.0 * np.pi, num = 360, endpoint=False)
        r = 10.0
        ra = ra_centre + r * np.cos(theta)
        dec = dec_centre + r * np.sin(theta)
    
    (ra_new, dec_new) = to_standard_position(ra, dec, ra_centre, dec_centre)
    
    plt.scatter(ra_new, dec_new, c = theta)
    plt.show()

def from_standard_position_test_harness():

    import numpy as np
    import matplotlib.pyplot as plt
    
    ra_centre = 72.0
    dec_centre = -38.0
    
    if False:
        ra = np.linspace(0.0, 359.0, num = 360)
        dec = ra * 0.0 + 1.0
        theta = ra
    else:
        
        theta = np.linspace(0.0, 2.0 * np.pi, num = 360, endpoint=False) + np.pi/4.0
        r = 10.0
        ra = r * np.cos(theta) + 180.0
        dec = r * np.sin(theta)
    
    (ra_new, dec_new) = from_standard_position(ra, dec, ra_centre, dec_centre)
    
    plt.scatter(ra_new, dec_new, c = theta)
    plt.show()
    
    
def to_from_standard_position_test_harness():
    
    import numpy as np
    
    num_tests = 1000000

    ra = np.random.uniform(size=num_tests) * 360.0
    dec = np.random.uniform(size=num_tests) * 180.0 - 90.0
    ra_centre = np.random.uniform() * 360.0
    dec_centre = np.random.uniform() * 180.0 - 90.0
    
    (ra_new, dec_new) = from_standard_position(ra, dec, ra_centre, dec_centre)
    (ra_new_new, dec_new_new) = to_standard_position(ra_new, dec_new, ra_centre, dec_centre)
    error = np.sqrt((ra_new_new - ra)**2 + (dec_new_new - dec)**2)
    print(np.max(error))


# See p. GL 140
def rotation_matrix(ra_centre, dec_centre):

    import numpy as np

    ra = np.array([0.0, 90.0, 0.0])
    dec = np.array([0.0, 0.0, 90.0])
    (new_ra, new_dec) = to_standard_position(ra, dec, ra_centre, dec_centre)
    x = np.cos(np.radians(new_dec)) * np.cos(np.radians(new_ra))
    y = np.cos(np.radians(new_dec)) * np.sin(np.radians(new_ra))
    z = np.sin(np.radians(new_dec))
    
    A = np.column_stack((x, y, z))
    return A
    
    
# Set 'to' to True to go to standard position and to False to go from standard position
def standard_position_fast_core(x, y, z, ra_centre, dec_centre, to):

    import numpy as np
    
    A = rotation_matrix(ra_centre, dec_centre)
    if not to:
        A = np.linalg.inv(A)
    new_xyz = np.dot(np.column_stack((x, y, z)), A)
    rotated_ra = np.degrees(cleanup_ra_rad(np.arctan2(new_xyz[:,1], new_xyz[:,0])))
    rotated_dec = np.degrees(np.arcsin(new_xyz[:,2]))
    return (rotated_ra, rotated_dec)

    
def to_standard_position_fast(x, y, z, ra_centre, dec_centre):
    return standard_position_fast_core(x, y, z, ra_centre, dec_centre, True)
    
def from_standard_position_fast(x, y, z, ra_centre, dec_centre):
    return standard_position_fast_core(x, y, z, ra_centre, dec_centre, False)
    
    
def to_from_standard_position_fast_test_harness():

    import numpy as np
    
    to = False

    ra = np.array([0.0, 90.0, 0.0, 34.0, 23.0, 67.0, 332.0, 175.0])
    dec = np.array([0.0, 00.0, 90.0, -34.0, 74.0, 67.0, -39.0, -88.0])
    
    x = np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y = np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z = np.sin(np.radians(dec))
    
    ra_centre = 180.0
    dec_centre = 0.0
    
    if to:
        res1 = to_standard_position_fast(x, y, z, ra_centre, dec_centre)
        res2 = to_standard_position(ra, dec, ra_centre, dec_centre)
    else:
        res1 = from_standard_position_fast(x, y, z, ra_centre, dec_centre)
        res2 = from_standard_position(ra, dec, ra_centre, dec_centre)
    
    
    print(res1)
    print("==========")
    print(res2)
    
    print(res1[0][0])
    


def one_pixel_healpix_map(nside, ra, dec, do_nest):

    import healpy as hp
    import numpy as np
    
    num_healpixels = hp.nside2npix(nside)
    values = np.zeros(num_healpixels)

    index = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    values[index] = 1.0
    
    return values
    
    
def is_in_standard_cutout(ra, dec, side_in_degrees):
    import numpy as np
    return np.where(np.logical_and((np.abs(ra - 180.0) < 0.5*side_in_degrees), (np.abs(dec - 0.0) < 0.5*side_in_degrees)))
    
    
    
# return 0 if the absolute difference is small, else return percentage difference.    
def comparison_function(x, y, tolerance):
    import numpy as np
    r = np.divide((x-y), y, out=np.zeros_like(y, dtype=float), where=y!=0.0)
    r[np.where(np.abs(x-y)<=tolerance)] = 0.0
    return r
    
    
    
def plot_several_healpix_maps():

    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    import skymapper
    
    nside = 1024

    maps = []
    titles = []

    # 1. maps from files
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    #filenames = ["Buzzard_192.90_110_2048_downgraded_to_1024_masked.glimpse.merged.values.dat", "Buzzard_192.90_110_2048_downgraded_to_1024.glimpse.merged.values.dat", "Buzzard_192.nside" + str(nside) + "_truth.dat"]
    #filenames = ["patch_1242_90_110_2048_values.dat"]
    filenames = ["Mcal_0.2_1.3_90_110_2048_downgraded_to_1024_masked.glimpse.merged.values.dat"]
        
    
    weight_maps = []
    for i in weight_maps:
        # Also show weights
        filenames.append(filenames[i].replace("_values", "_weights"))
        
    masked_maps = []
    for i in masked_maps:
        filenames[i] = filenames[i].replace("_values", "_values_masked")
    
    
    for f in filenames:
        print("Using file {}".format(f))
        maps.append(hp.read_map(path + f))
        titles.append(f)
        
    # 2. other maps
        
    if False:
        # Also plot some glimpse input data
        f = "Mcal_0.2_1.3.1689.glimpse.cat.fits"
        maps.append(fits_catalog_to_healpix_map(path + f, nside, False))
        titles.append(f)
    
    if False:
        # Also plot glimpse output lattice
        f = "Mcal_0.2_1.3.1689.glimpse.out.fits"
        maps.append(glimpse_output_to_healpix_map(path + f, nside*4, False))
        titles.append(f + " values")
        
    if False:
        # Also plot glimpse lattice points
        f = "Mcal_0.2_1.3.1689.glimpse.out.fits"
        maps.append(glimpse_lattice_points_to_healpix_map(path + f, nside*4, False))
        titles.append(f + " lattice points")
        
    if False:
        # Also plot just one healpixel
        ra = 56.25
        dec = -27.2796127
        maps.append(one_pixel_healpix_map(nside, ra, dec, False))
        titles.append("RA={}; DEC={}; pixel_id={}".format(ra, dec, hp.ang2pix(nside, ra, dec, nest=False, lonlat=True)))
        

    # To smooth a map, include in maps_to_smooth the index of the map (i.e. the map's index in the array 'maps') 
    maps_to_smooth = []
    for i in maps_to_smooth:
        one_arcmin_in_radians = 0.000290888
        smoothing_scale_in_arcmin = 7.0
        maps[i] = hp.smoothing(maps[i], sigma = smoothing_scale_in_arcmin * one_arcmin_in_radians)
        titles[i] += " smoothed at {} arcmin".format(smoothing_scale_in_arcmin)
    
    
    if False:
        # Show diff between maps[0] and maps[1]
        tolerance = 0.0009
        percentage_diff = comparison_function(maps[0], maps[1], tolerance)
        maps.append(percentage_diff)
        titles.append("Percentage difference")
    
    
    for map, title, i in zip(maps, titles, range(len(maps))):
        if False:
            hp.mollview(map, fig=i, title=title)
        elif False:
            #rot = (326.25, 12.024699, 0.0)
            #rot = (180.0, 0.0, 0.0)
            rot = (75.0, -55.0, 0.0)
            hp.gnomview(map, fig=i, rot=rot, title=title, reso=3.3, xsize=400) # , max=0.104, min=-0.0264)
        else:
            rot = (75.0, -55.0, 0.0)
            hp.cartview(map, fig=i, rot=rot, title=title, xsize=400) # , max=0.104, min=-0.0264)
        
        hp.graticule(dpar=5.0)
    
    plt.show()
    

def downgrade_map():
    
    import healpy as hp
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    f_in = "Mcal_0.2_1.3_90_110_2048.glimpse.merged.values.dat"
    f_out = "Mcal_0.2_1.3_90_110_2048_downgraded_to_1024.glimpse.merged.values.dat"
    
    
    m = hp.read_map(path + f_in)
    m_new = hp.ud_grade(m, 1024)
    hp.write_map(path + f_out, m_new, overwrite=True)



def clean_up_edges():

    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    
    do_nest = False
    nside = 1024
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    project = "Mcal_0.2_1.3_90_110_2048_downgraded_to_1024"
    old_filename = project + ".glimpse.merged.values.dat"
    new_filename = project + "_masked.glimpse.merged.values.dat"
    new_filename_png = project + "_masked.glimpse.merged.values.png"
    
    original_catalogue_file_name = metacal_data_file_name()
    
    m = hp.read_map(path + old_filename)
    
    (ra, dec) = get_from_fits_file(original_catalogue_file_name, ["RA", "DEC"])
    
    # Create a mask that is zero in healpixels where there are no source galaxies and one elsewhere.
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
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
        
    
    
    
def bound_between(x, lower, upper):
    import numpy as np
    return np.minimum(np.maximum(x, lower), upper)


def index_into_glimpse_array_test_harness():
    import numpy as np
    p_shape = [6, 6]
    
    glimpse_lattice_spacing = 1.2
    
    x = 0.001
    y = 0.001
    print(x,y)
    print(bound_between(p_shape[0]//2 + np.floor(x / glimpse_lattice_spacing), 0, p_shape[0]-1) + \
                p_shape[0] * bound_between(p_shape[1]//2 - np.floor(y / glimpse_lattice_spacing) - 1, 0, p_shape[1]-1))
    
    

def ra_dec_to_healpixel_id_test_harness():
    import healpy as hp
    import numpy as np
    
    ra = 0.4464
    dec = -0.2678
    nest = False
    nside = 1024
    
    print(hp.ang2pix(nside, ra, dec, nest, lonlat=True))
    
    
    

# No error if directory already exists
# See https://stackoverflow.com/questions/273192
def create_directory_if_necessary(directory_name):
    import os, errno
    try:
        os.makedirs(directory_name)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

# percentage should be an integer e.g. 60 for 60%. Enter a negative number to clean-up at the end.
def update_percentage_bar(percentage):
    import sys
    if percentage >= 0:
        sys.stdout.write('\r' + str(percentage).rjust(3) + r'%' + ' ' * percentage + '█' * (100 - percentage) + ' ')
        sys.stdout.flush()
    else:
        sys.stdout.write('\r' + str(100) + r'%' + ' ' * 105)
        sys.stdout.flush()
        print('\n')


# The input coordinates describe two points p1 and p2 on sphere (ra in [0, 360] and dec in [-90, 90]).
# Return value is the angle p1->0->p2 in radians.
# See https://en.wikipedia.org/wiki/Haversine_formula
def angular_separation(ra1, dec1, ra2, dec2):

	import numpy as np

	long1 = np.radians(ra1)
	long2 = np.radians(ra2)
	lat1 = np.radians(dec1)
	lat2 = np.radians(dec2)
	half_dlong = (long2 - long1) / 2.0
	half_dlat = (lat2 - lat1) / 2.0
	a = np.sin(half_dlat)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(half_dlong)**2
	return 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    
def angular_separation_test_harness():
    import numpy as np
    print(np.degrees(angular_separation(0.0, 0.0, 5.0, 35.0)))


# Same rules as for angular_separation
def angular_separation_fast(sin_ra1, cos_ra1, sin_dec1, cos_dec1, ra2, dec2):

    import numpy as np
    
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)
    
    sin_ra2 = np.sin(ra2_rad)
    cos_ra2 = np.cos(ra2_rad)
    sin_dec2 = np.sin(dec2_rad)
    cos_dec2 = np.cos(dec2_rad)
    
    # See my book p. GL 143
    a = (0.5 - 0.5 * cos_dec1 * cos_dec2 - 0.5 * sin_dec1 * sin_dec2) + cos_dec1 * cos_dec2 * (0.5 - 0.5 * cos_ra1 * cos_ra2 - 0.5 * sin_ra1 * sin_ra2)
    a = bound_between(a, 0.0, 1.0) # 'a' might fall slightly outside [0, 1] due to numerical issues; fix this up.
    return 2.0 * np.arcsin(np.sqrt(a))
    
    
def angular_separation_fast_test_harness():
    import numpy as np
    ra1 = np.array([0.0, 40.0, 90.0, 180.0, 271.0, 280.0])
    dec1 = np.array([0.0, 40.0, -90.0, 87.0, 0.0, -35.0])
    
    ra2 = 270.0
    dec2 = 0.0
    
    sin_ra1 = np.sin(np.radians(ra1))
    cos_ra1 = np.cos(np.radians(ra1))
    sin_dec1 = np.sin(np.radians(dec1))
    cos_dec1 = np.cos(np.radians(dec1))
    
    res1 = angular_separation_fast(sin_ra1, cos_ra1, sin_dec1, cos_dec1, ra2, dec2)
    res2 = angular_separation(ra1, dec1, ra2, dec2)
    
    print(res1)
    print("=========")
    print(res2)




# Will process ids in the range [ids_to_process_start, ids_to_process_end). Set ids_to_process_end to -1 to mean "to the end"
def create_cutouts(input_catalogue, catformat, raname, decname, shear_names, other_field_names, nside, cutout_side_in_degrees, ids_to_process_slice, output_directory, output_file_root):

    from astropy.table import Table
    import healpy as hp
    import numpy as np
    import math
    import astropy.io.fits as pyfits
    import sys
    
    
    
    
    create_directory_if_necessary(output_directory)
    assert catformat in ("csv", "fits"), "Catalogue format must be either 'fits' or 'csv'; {0} was given".format(catformat)
    assert hp.isnsideok(nside), "nside must be a valid Healpix nside"
    all_field_names = [raname, decname]
    if other_field_names:
        all_field_names.extend(other_field_names.split(","))
    if shear_names:
        parsed_shear_names = shear_names.split(",")
        assert (len(parsed_shear_names)%2==0), "shear_names should be omitted or be a comma-separated list of (an even number of) field names"
        all_field_names.extend(parsed_shear_names)
    
    
    print('Processing ' + input_catalogue + '...')
    data = Table.read(input_catalogue, format=catformat)
    for field_name in all_field_names:
        assert field_name in data.colnames, "{0} is not a column in the input file {1}".format(field_name, input_catalogue)
    data = data[all_field_names]
        
    ra = data[raname]
    dec = data[decname]
    
    assert (all(ra >= 0.0) and all(ra <= 360.0)), "RA should be decimal degrees between 0 and 360."
    assert (all(dec >= -90.0) and all(dec <= 90.0)), "Dec should be decimal degrees between -90 and 90."
    
    num_healpixels = hp.nside2npix(nside)
    num_digits = 1 + int(math.log10(num_healpixels))
    
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    
    sin_ra = np.sin(ra_rad)
    cos_ra = np.cos(ra_rad)
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    
    x = cos_dec * cos_ra
    y = cos_dec * sin_ra
    z = sin_dec
    
    # Create dictionary to hold information (for certain healpixel_ids only) about
    # the minimum angular separation between the catalogue points and the healpixels centre.
    # Key is heaslpixel_id; value is a triple (ra_centre, dec_centre, minimum_separation (in degrees)).
    exclusion_dict = dict()
    
    # If catalogue points are further than critical_angular_separation away from a healpixel centre, then 
    # they will certainly get excluded from this cutout.
    # See p. GL 127. The prefactor allows for a safety margin and handles some 'tangent plane to sphere' distortion...
    
    critical_angular_separation = 1.1 * (cutout_side_in_degrees * np.sqrt(2.0) / 2.0)

    for healpix_id in range(num_healpixels)[slice_from_string(ids_to_process_slice)]:

        print("Processing {} of {}".format(healpix_id, num_healpixels))
        
        if healpix_id < 0 or healpix_id >= num_healpixels:
            print("Excluding healpix_id {} as it is not in range [0, {})".format(healpix_id, num_healpixels))
        else:
        
            (ra_centre, dec_centre) = hp.pix2ang(nside, healpix_id, nest=False, lonlat=True)
            
            # See if we can immediately exclude this pixel based on previous results.
            can_exclude_immediately = False
            for id_prev in exclusion_dict:
                (ra_centre_prev, dec_centre_prev, ang_sep_prev) = exclusion_dict[id_prev]
                # See picture p. GL 147
                separation_between_centres = np.degrees(angular_separation(ra_centre_prev, dec_centre_prev, ra_centre, dec_centre))
                if separation_between_centres < ang_sep_prev - critical_angular_separation:
                    can_exclude_immediately = True
                    print("\tExcluding healpix_id {} as it is close ({} deg) to pixel {} already excluded (distance from catalogue {} deg)".format(healpix_id, separation_between_centres, id_prev, ang_sep_prev))
            
            if not can_exclude_immediately:
            
                
                ang_sep_data_to_centre = np.degrees(angular_separation_fast(sin_ra, cos_ra, sin_dec, cos_dec, ra_centre, dec_centre))
                min_ang_sep_data_to_centre = np.min(ang_sep_data_to_centre)
                if min_ang_sep_data_to_centre > critical_angular_separation:
                    # Now we can exclude because we see that no points in the catalogue are sufficiently
                    # close to the centre of healpix_id.
                    exclusion_dict[healpix_id] = (ra_centre, dec_centre, min_ang_sep_data_to_centre)
                    print("\tExcluding healpix_id {} as its distance from the catalogue is {} deg".format(healpix_id, np.min(min_ang_sep_data_to_centre)))
                else:
                    # Which data lie in a small disc centred on this healpixel?
                    filter1 = np.where(ang_sep_data_to_centre <= critical_angular_separation)
                    
                    # Put these data in standard position
                    (ra_standardised_filter1, dec_standardised_filter1) = to_standard_position_fast(x[filter1], y[filter1], z[filter1], ra_centre, dec_centre)
                    
                    # Now which of the (filtered) data lie in the standard cutout?
                    filter2 = is_in_standard_cutout(ra_standardised_filter1, dec_standardised_filter1, cutout_side_in_degrees)
                    
                    # Put the two filters together to get a filter on the original data set. See p. GL148 for joint filtering.
                    joint_filter = (filter1[0])[filter2]
                    
                    if joint_filter.size == 0:
                        print("\tExcluding healpix_id {} as it has no galaxies in the cutout area".format(healpix_id))
                    else:
                        # Some galaxies to write
                        output_filename = (output_directory + "/" + output_file_root).format(str(healpix_id).zfill(num_digits))
                        
                        field_names = []
                        data_columns = []
                        
                        field_names.append(raname)
                        data_columns.append(ra_standardised_filter1[filter2])
                                    
                        field_names.append(decname)
                        data_columns.append(dec_standardised_filter1[filter2])
                        
                        
                        if shear_names:
                            parsed_shear_names = shear_names.split(",")
                            for (shear_name_1, shear_name_2) in zip(parsed_shear_names[0::2], parsed_shear_names[1::2]):
                                (shear1, shear2) = rotate_shear_45_degrees(data[shear_name_1], data[shear_name_2])
                                
                                field_names.append(shear_name_1)
                                data_columns.append(shear1[joint_filter])
                                
                                field_names.append(shear_name_2)
                                data_columns.append(shear2[joint_filter])
                        
                        for f in other_field_names.split(","):
                            field_names.append(f)
                            data_columns.append(data[f][joint_filter])
                           
                        write_to_fits_file(output_filename, field_names, data_columns)



def create_cutouts_caller(ini_file_name, ids_to_process_slice):

    import configparser
    
    config = configparser.ConfigParser()
    config.read(ini_file_name)

    input_catalogue = config["create_cutouts"].get("input_catalogue")
    catformat = config["create_cutouts"].get("catformat", "fits")
    ra_name = config["create_cutouts"].get("ra_name", "RA")
    dec_name = config["create_cutouts"].get("dec_name", "DEC")
    shear_names = config["create_cutouts"].get("shear_names")
    other_field_names = config["create_cutouts"].get("other_field_names")
    nside = int(config["create_cutouts"].get("nside", "16"))
    cutout_side_in_degrees = float(config["create_cutouts"].get("cutout_side_in_degrees", "16"))
    output_directory = config["project"].get("directory")
    output_file_root = config["project"].get("project_name") + ".{}.glimpse.cat.fits"
    
    create_cutouts(input_catalogue, catformat, ra_name, dec_name, shear_names, other_field_names, nside, cutout_side_in_degrees, ids_to_process_slice, output_directory, output_file_root)


    

##################################### Start of solver code #####################################
# See p. GL 136. Nice idea, didn't work (numerical solution had too much error...)  

def objective_function(x, centre):

    import math
    import numpy as np
    
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    
    cx = centre[0]
    cy = centre[1]
    cz = centre[2]
    
    y0 = a**2 + b**2 + c**2 + d**2 - 1.0
    
    #y1 = (a**2 + b**2 - c**2 - d**2) * cx + 2.0 * (b * c - a* d) * cy + 2.0 * (b * d + a * c) * cz + 1.0
    y1 = (a**2 + b**2 - c**2 - d**2) * cx + 2.0 * (b * c - a * d) * cy + 2.0 * (b * d + a * c) * cz - 1.0
        
    y2 = 2.0 * (b * c + a * d) * cx + (a**2 + c**2 -b**2 -d**2) * cy + 2.0 * (c * d - a * b) * cz
    
    #y3 = 2.0 * (b * c + a * d) * cy - (a**2 + c**2 -b**2 -d**2) * cx - math.sqrt(2.0) / 2.0
    #y3 = 2.0 * (b * c + a * d) * cy - (a**2 + c**2 -b**2 -d**2) * cx - 1.0
    y3 = d
    
    y = np.array([y0, y1, y2, y3])
    
    print(x, " |-->  ", y)
    
    return y



def objective_function_foo(x, centre):

    import math
    import numpy as np
    
    a = x[0]
    b = x[1]
    c = x[2]
    d = x[3]
    
    
    y0 = a**2 + b**2 + c**2 + d**2 - 1.0
    y1 = (a**2 + b**2 - c**2 - d**2) + 1.0
    y2 = 2.0 * (b * d + a * c)
    y3 = a**2 - b**2 - c**2 + d**2 - 1.0
    
    y = np.array([y0, y1, y2, y3])
    
    Jac = 2.0 * np.array([[a, a, -c, a], [b, b, d, -b], [c, -c, -a, -c], [d, -d, b, d]])
    
    
    print(x, " |-->  ", y)
    
    return (y, Jac)
    #return y



def solver():
    
    import scipy.optimize as so
    import numpy as np
    
    print("========================================================================================")


    cx = 0.999
    cy = np.sqrt(1.0 - cx**2)
    cz = 0.0
    
    first_guess = np.array([1.0, 0.0, 0.0, 0.0])
    centre = np.array([cx, cy, cz])
    print(centre)
    
    options = dict()
    options["xtol"] = 1e-8
    
    
    sol = so.root(objective_function_foo, first_guess, centre, "hybr", True, options=options)
    print("===========")
    print("===========")
    print("===========")
    print(sol.x)
    print(sol.message)
    
    a = sol.x[0]
    b = sol.x[1]
    c = sol.x[2]
    d = sol.x[3]
    
    A = np.array([[a**2+b**2-c**2-d**2, 2.0*(b*c-a*d), 2.0*(b*d+a*c)], [2.0*(b*c+a*d), a**2-b**2+c**2-d**2, 2.0*(c*d-a*b)], [2.0*(b*d-a*c), 2.0*(c*d+a*b), a**2-b**2-c**2+d**2]])
    print(A)
    
##################################### End of solver code #####################################    

# See p. GL148
def joint_filter_example():

    import numpy as np
    
    data = np.array(range(8)) + 3.0
    filter1 = np.where(data > 6.0)
    filter2 = np.where(data[filter1] < 5.0)
    joint_filter = (filter1[0])[filter2]
    if joint_filter.size > 0:
        print(joint_filter)
        print(data[joint_filter])
    
    

if __name__ == '__main__':


    #print_list_of_healpixels()
    #kappa_values_in_one_fine_pixel()
    #tester()
    #save_buzzard_truth()
    #sphere_to_tangent_plane_mapping_test_harness()
    #index_into_glimpse_array_test_harness()
    #ra_dec_to_healpixel_id_test_harness()
    #plot_several_healpix_maps()
    #to_standard_position_test_harness()
    #from_standard_position_test_harness()
    #to_from_standard_position_test_harness()
    #correct_one_shear_catalogue_caller()
    #redshift_histogram()
    #downgrade_map()
    #clean_up_edges()
    #shear_stdev()
    #add_dummy_redshift_column(metacal_data_file_name())
    #angular_separation_test_harness()
    #num_files_in_directory()
    #to_from_standard_position_fast_test_harness()
    #angular_separation_fast_test_harness()
    #correct_one_shear_catalogue("/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/cat_3_0.fits")
    #glimpse_sign_experiment()
    #kappa_histogram()
    #show_glimpse_output_as_image()
    compare_two_cutouts()
    #joint_filter_example()
    
