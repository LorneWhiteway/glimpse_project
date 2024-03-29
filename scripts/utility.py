#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
    Core file for 'glimpse on curved sky' project.
    Author: Lorne Whiteway.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.io.fits as pyfits
import glob
import healpy as hp
import math
import sys
import os
import errno
import configparser
import subprocess
import utility
import configparser
import datetime
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
import scipy.stats as st
import scipy.ndimage as ndi



########################## Start of one-off utilities ##########################


def clean_up_edges_caller():
    for root in ["1", "2", "3", "4", "5"]:
        input_catalogue_filename = "/share/splinter/ucapwhi/glimpse_project/data/data_catalogs_weighted_all_bins.fits"
        ra_name = "ra"
        dec_name = "dec"
        input_map_filename = "/share/splinter/ucapwhi/glimpse_project/runs/20201023_unblinded_lambda_{}/glimpse.merged.values.fits".format(root)
        print("Reading from {}".format(input_map_filename))
        map_to_be_masked = hp.read_map(input_map_filename)
        masked_map = clean_up_edges(input_catalogue_filename, ra_name, dec_name, map_to_be_masked)    
        output_filename = "/share/splinter/ucapwhi/glimpse_project/runs/20201023_unblinded_lambda_{}/map_kappa_1024_unblinded_data_20201105_lambda_{}.dat".format(root, root)
        print("Writing to {}".format(output_filename))
        write_healpix_array_to_file(masked_map, output_filename, False)



# This function written by Macro Gatti. See Slack message on 22 December 2019.
def apply_random_rotation(e1_in, e2_in):
    
    np.random.seed() # CRITICAL in multiple processes !
    rot_angle = np.random.rand(len(e1_in)) * 2.0 * np.pi #no need for 2?
    cos = np.cos(rot_angle)
    sin = np.sin(rot_angle)
    e1_out = + e1_in * cos + e2_in * sin
    e2_out = - e1_in * sin + e2_in * cos
    return e1_out, e2_out
    
    
def append_random_shear_to_Buzzard():
    list_of_field_names = ["RA", "DEC", "E1", "E2", "G1", "G2", "k_orig", "true_z"]
    (ra, dec, e1, e2, g1, g2, k_orig, true_z) = get_from_fits_file(buzzard_data_file_name(), list_of_field_names)
    print("{} galaxies".format(len(ra)))
    print(list_of_field_names)
    print("About to call apply_random_rotation")
    (e1_random, e2_random) = apply_random_rotation(e1, e2)
    list_of_data_columns = (ra, dec, e1, e2, g1, g2, k_orig, true_z, e1_random, e2_random)
    list_of_field_names.append("E1_RANDOM")
    list_of_field_names.append("E2_RANDOM")
    print("Finished calling apply_random_rotation")
    print(list_of_field_names)
    output_filename = "/share/splinter/ucapwhi/glimpse_project/Buzzard_192_with_random.fits"
    write_to_fits_file(output_filename, list_of_field_names, list_of_data_columns, True)
    

# See p. GL148
def joint_filter_example():

    data = np.array(range(8)) + 3.0
    filter1 = np.where(data > 6.0)
    filter2 = np.where(data[filter1] < 5.0)
    joint_filter = (filter1[0])[filter2]
    if joint_filter.size > 0:
        print(joint_filter)
        print(data[joint_filter])
    


def compare_two_cutouts():

    (ra_1, dec_1, e1_1, e2_1) = get_from_fits_file("/share/splinter/ucapwhi/glimpse_project/output/Mcal_0.2_1.3.2218.cat.fits", ["RA", "DEC", "E1", "E2"])
    (ra_2, dec_2, e1_2, e2_2) = get_from_fits_file("/share/splinter/ucapwhi/glimpse_project/output_Mcal_badsign/Mcal_0.2_1.3.2218.cat.fits", ["RA", "DEC", "E1", "E2"])
    
    for i in range(3159, 3170):
        print(i, ra_1[i], ra_2[i], dec_1[i], dec_2[i])
        



def correct_one_shear_catalogue(catalogue_filename):
    
    print("Processing {}...".format(catalogue_filename))
    
    list_of_field_names = ["ra_gal", "dec_gal", "e1_gal", "e2_gal", "z"]
    list_of_data_columns = get_from_fits_file(catalogue_filename, list_of_field_names)
        
    list_of_data_columns[2] *= -1.0 # e1_gal
    list_of_data_columns[3] *= -1.0 # e2_gal
    
    output_filename = catalogue_filename.replace(".fits", ".negative.fits")
    
    write_to_fits_file(output_filename, list_of_field_names, list_of_data_columns, True)
    
    
def glimpse_sign_experiment():

    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/output.fits"
    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_file)
    
    glimpse_output_file_negative = "/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/output.negative.fits"
    (ra_negative, dec_negative, kappa_negative) = get_glimpse_output_data(glimpse_output_file_negative)
    
    plt.scatter(kappa, kappa_negative, s=1)
    plt.show()
    
    

def correct_one_shear_catalogue_caller():
    
    filelist = glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.cat.fits")
    for f in filelist: 
        correct_one_shear_catalogue(f)
   
    

def num_files_in_directory():
    print(len(glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Mcal_0.2_1.3.*.cat.fits")))



# Written to produce data for Niall; used on 9 Dec 2019
def redshift_histogram():
    
    np.set_printoptions(precision=3)
    
    top_z = 2.5
    
    (ra, dec, true_z) = get_from_fits_file(buzzard_data_file_name(), ["RA", "DEC", "true_z"])
    (hist, bin_edges) = np.histogram(true_z, np.linspace(0.0, top_z, int(100*top_z) + 1))
    for (b0, b1, h) in zip(bin_edges[:-1], bin_edges[1:], hist):
        print("{0:.2f}\t{1:.2f}\t{2:d}".format(b0, b1, h))


# Written to produce data for Niall; used on 10 Dec 2019
def shear_stdev():
    (e1, e2) = get_from_fits_file(buzzard_data_file_name(), ["E1", "E2"])
    e1_and_e2 = np.concatenate((e1, e2))
    print(np.std(e1_and_e2))

    
def kappa_histogram():

    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    f = "Mcal_0.2_1.3_90_110_2048_downgraded_to_1024_masked.glimpse.merged.values.fits"
    m = hp.read_map(path + f)
    
    plt.hist(m[np.where(m != 0.0)], bins = 50)
    plt.show()
    

def glimpse_output_as_array(glimpse_output_filename):
    
    
    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_filename)
    n = int(math.sqrt(kappa.shape[0]))
    return np.reshape(kappa, [n, n])
    
    

# See also plot_several_glimpse_outputs, which does something similar...
def show_glimpse_output_as_image():

    glimpse_output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/2277.glimpse.out.fits"
    kappa_as_2d_array = glimpse_output_as_array(glimpse_output_filename)
    plt.imshow(kappa_as_2d_array)
    plt.show()
    
    
def define_Buzzard_ten_percent_by_area_subset():

    centre_ra = 40.0
    centre_dec = -30.0
    

    if False:
        print("Getting data...")
        (g_ra, g_dec) = get_from_fits_file(buzzard_data_file_name(), ["RA", "DEC"])

        nside = 16
        
        print("Getting healpixel ids...")
        g_ids = hp.ang2pix(nside, g_ra, g_dec, False, True)
        print("Sorting and removing duplicates...")
        
        s = set(g_ids)
        n_ids = np.array(list(s))
        n_ids.sort()
        
        print("Processing...")
        
        (n_ra, n_dec) = hp.pix2ang(nside, n_ids, False, True)
        
            
        
        n_dist = np.degrees(angular_separation(n_ra, n_dec, centre_ra, centre_dec)) # In degrees
        
        num_n_ids = len(n_ids)
        
        for r in np.linspace(0.0, 30.0, 31, True): # r in degrees
            m_ids = n_ids[np.where(n_dist<r)]
            num_m_ids = len(m_ids)
            print(r, num_m_ids, num_n_ids, num_m_ids/num_n_ids)
    
    if True:
        
        list_of_field_names = ["RA", "DEC", "E1", "E2", "G1", "G2", "k_orig", "true_z", "E1_RANDOM", "E2_RANDOM"]
        
        print("Getting data...")
        (ra, dec, e1, e2, g1, g2, k, z, e1r, e2r) = get_from_fits_file(buzzard_data_file_name(), list_of_field_names)
        
        r = 10.0 # In degrees
        print("Filtering with r = {} degrees...".format(r))
        n_dist = np.degrees(angular_separation(ra, dec, centre_ra, centre_dec)) # In degrees
        filter = np.where(n_dist < r)
        
        print("Writing...")
        output_filename = "/share/testde/ucapnje/buzzard_desy3_marco/Buzzard_192_reduced.fits"
        write_to_fits_file(output_filename, list_of_field_names, [ra[filter], dec[filter], e1[filter], e2[filter], g1[filter], g2[filter], k[filter], z[filter], e1r[filter], e2r[filter]], False)
        
        


def save_buzzard_truth(buzzard_data_filename, output_filename):

    do_nest = False
    
    nside = 1024
    num_healpixels = hp.nside2npix(nside)
    sums = np.zeros(num_healpixels)
    count = np.zeros(num_healpixels)
    
    (ra, dec, kappa) = get_from_fits_file(buzzard_data_filename, ["RA", "DEC", "k_orig"])
    num_buzzard_data = len(ra)
    print("num_buzzard_data = {}".format(num_buzzard_data))
    id = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for (i, this_id, this_kappa) in zip(range(num_buzzard_data), id, kappa):
        sums[this_id] += this_kappa
        count[this_id] += 1
        if i % 1000000 == 0:
            print("{} of {}".format(i, num_buzzard_data))
        
        
    # The following division code treats 0/0 as 0. From https://stackoverflow.com/questions/26248654.
    average_values = np.divide(sums, count, out=np.zeros_like(count, dtype=float), where=count!=0.0)
    
    # Save
    write_healpix_array_to_file(average_values, output_filename, do_nest)
    
    
def run_save_buzzard_truth():
    save_buzzard_truth(buzzard_reduced_data_file_name(), "/share/splinter/ucapwhi/glimpse_project/experiments/truth/true_kappa.fits")

    

# One-time test; see p GL97.
def glimpse_array_order():
    
    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.1440.glimpse.out.fits"
    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_file)
    
    plt.plot(ra[:1000], c='red')
    plt.plot(dec[:1000], c='blue')
    plt.show()
    
    
def create_test_catalogue():

    main_catalogue_filename = "/share/testde/ucapnje/year3_des/Mcal_0.2_1.3.fits"
    list_of_field_names = ["RA", "DEC", "E1", "E2", "DUMMY_Z"]
    (ra, dec, e1, e2, z) = get_from_fits_file(main_catalogue_filename, list_of_field_names)
    
    a = angular_separation(ra, dec, 50.0, -30.0)
    r = np.random.uniform(size=ra.shape[0])
    filter = np.where(np.logical_and(a < 0.01, r > 0.8))
    ra = ra[filter]
    dec = dec[filter]
    e1 = e1[filter]
    e2 = e2[filter]
    z = z[filter]
    list_of_data_columns = [ra, dec, e1, e2, z]
    test_catalogue_filename = "/share/splinter/ucapwhi/glimpse_project/test/catalogue.fits"
    write_to_fits_file(test_catalogue_filename, list_of_field_names, list_of_data_columns, True)
    
    
def weight_test():
    print(hp.ang2pix(16, 30, -30, False, True))
    
    
def set_weights_in_catalogue_file():
    
    input_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/2277.cat.fits"
    (ra, dec, e1, e2, w) = get_from_fits_file(input_filename, ["RA", "DEC", "E1", "E2", "W"])
    
    ones = np.ones(len(ra)) # To be used for redshift
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/A/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "W", "DUMMY_Z"], (ra, dec, e1, e2, ones, ones), True)
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/B/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "W", "DUMMY_Z"], (ra, dec, e1, e2, ones*10.0, ones), True)
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/C/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "W", "DUMMY_Z"], (ra, dec, e1, e2, ones*10.0, ones), True)
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/D/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "DUMMY_Z"], (ra, dec, e1, e2, ones), True)
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/E/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "W", "DUMMY_Z"], (ra, dec, e1, e2, w, ones), True)
    
    output_filename = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/F/2277.cat.fits"
    write_to_fits_file(output_filename, ["RA", "DEC", "E1", "E2", "W", "DUMMY_Z"], (ra, dec, e1, e2, w / np.average(w), ones), True)
    
    


# We assume that the pickle files that we encounter have the following structure:
# Top level dictionary has keys that are non-negative integers that refer to tomographic bins (one of which 
# might be the total across all bins).
# The value associated with one of these keys is another dictionary; it maps field names to
# numpy data arrays.


# Use this to test the format of a pickle file.
def show_pickle_file_structure(filename):

    print("Pickle file structure for {}".format(filename))
    a = np.load(filename, allow_pickle=True, fix_imports=True, mmap_mode='r')
    print("Top level keys:")
    print(a.keys())
    print("Second level keys:")
    for k in a.keys():
        print(k, a[k].keys())
    print("Statistics of dec field for each subkey:")
    for k in a.keys():
        d = a[k]["dec"]
        print(k, len(d), min(d), max(d))
        

def show_pickle_file_structure_caller():
    show_pickle_file_structure("/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_1_dec2020.pkl")

def pickle_to_fits_caller():

    # Set A
    #pickle_filename = "/share/splinter/ucapwhi/glimpse_project/data/data_catalogs_weighted.pkl"
    #list_of_field_names = ["ra", "dec", "e1", "e2", "w"]
    #pickle_dataset = 4
    #output_fits_filename = "/share/splinter/ucapwhi/glimpse_project/data/data_catalogs_weighted_all_bins.fits"
    
    # Set B 
    #pickle_filename = "/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_full.pkl"
    #list_of_field_names = ["ra", "dec", "e1", "e2", "w", "k_orig"]
    #pickle_dataset = 4
    #output_fits_filename = "/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_all_bins.fits"
    
    # Set C
    #pickle_filename = "/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_full_dec2020.pkl"
    #list_of_field_names = ["ra", "dec", "e1", "e2", "w", "k_orig"]
    #pickle_dataset = 4
    #output_fits_filename = "/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_all_bins_dec2020.fits"
    
    # Sets D, E, F, G
    for bin_number in range(4):
        pickle_filename = "/share/splinter/ucapwhi/glimpse_project/data/data_catalogs_weighted.pkl"
        list_of_field_names = ["ra", "dec", "e1", "e2", "w"]
        pickle_dataset = bin_number
        output_fits_filename = "/share/splinter/ucapwhi/glimpse_project/data/data_catalogs_weighted_bin_{}.fits".format(bin_number)

        pickle_to_fits(pickle_filename, pickle_dataset, list_of_field_names, output_fits_filename)
    



# pickle_dataset specifies which key from the top-level dictionary to use.
def pickle_to_fits(pickle_filename, pickle_dataset, list_of_field_names, output_fits_filename):
    a = np.load(pickle_filename, allow_pickle=True, fix_imports=True, mmap_mode='r')
    output_list_of_data_columns = [a[pickle_dataset][field_name] for field_name in list_of_field_names]
    output_list_of_field_names = [field_name.upper() for field_name in list_of_field_names]
    write_to_fits_file(output_fits_filename, output_list_of_field_names, output_list_of_data_columns, verbose=True)
    
    

########################## End of one-off utilities ##########################


    
# job_control will be a string such as "" or ":-1" or "2:17:3"
# Assumes 'a' is one-dimensional
# I tried https://stackoverflow.com/questions/48494581 but it doesn't handle omitted arguments.
def array_slice_from_job_control_string(a, job_control):
    tokens = job_control.split(':')
    tokens.extend(("", "", ""))
    start = int(tokens[0]) if len(tokens[0]) > 0 else 0
    end = int(tokens[1]) if len(tokens[1]) > 0 else len(a)
    stride = int(tokens[2]) if len(tokens[2]) > 0 else 1
    return a[start:end:stride]



    
def array_slice_from_job_control_string_test_harness():
    
    a = np.arange(20)
    job_control = ""
    
    print(array_slice_from_job_control_string(a, job_control))
    
    
    



# Returns the flattened ra, dec and kappa (in that order).
def get_glimpse_output_data(glimpse_output_filename):
    
    pyfile = pyfits.open(glimpse_output_filename)
    ra = pyfile[1].data
    dec = pyfile[2].data
    kappa = pyfile[0].data
    pyfile.close()
    return (ra.reshape(-1), dec.reshape(-1), kappa.reshape(-1))





def kappa_values_in_one_fine_pixel():

    
    (r,d) = hp.pix2ang(16, 2218, False, True)
    specific_fine_healpix_id = hp.ang2pix(1024, r, d, False, lonlat=True)
    
    
    j = 0
    for f in sorted(glob.glob("/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192.*.glimpse.out.fits")):
        (ra, dec, kappa) = get_glimpse_output_data(f)
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
    
#"RA", "DEC", "true_z", "E1", "E2", "G1", "G2", "k_orig"
def buzzard_reduced_data_file_name():
    return '/share/testde/ucapnje/buzzard_desy3_marco/Buzzard_192_reduced.fits'


#"RA", "DEC", "E1", "E2", "E1_RANDOM", "E2_RANDOM", "DUMMY_Z"
def metacal_data_file_name():
    return '/share/testde/ucapnje/year3_des/Mcal_0.2_1.3.fits'
    
    
# 'what_to_get' should be a list of field names (case insensitive).
# Returns a tuple of arrays.
def get_from_fits_file(file_name, what_to_get):
    res = []
    if what_to_get: # i.e. if list is not empty
        x = pyfits.open(file_name)
        for w in what_to_get:
            res.append(x[1].data.field(w))
    return res
    



def write_to_fits_file(output_filename, list_of_field_names, list_of_data_columns, verbose):
    if verbose:
        print("Writing to {}...".format(output_filename))
    column_info = []
    for (field_name, data_column) in zip(list_of_field_names, list_of_data_columns):
        column_info.append(pyfits.Column(name=field_name, format='D', array=data_column))
    tbhdu = pyfits.BinTableHDU.from_columns(column_info)
    tbhdu.writeto(output_filename, overwrite=True)



def fits_catalog_to_healpix_map(catalog_file_name, nside, do_nest):

    (ra, dec) = get_from_fits_file(catalog_file_name, ["ra", "dec"])
    
    num_healpixels = hp.nside2npix(nside)
    res = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for id in ids:
        res[id] += 1.0
    return res    
    
def glimpse_output_to_healpix_map(glimpse_output_file, nside, do_nest):

    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_file)
    
    num_healpixels = hp.nside2npix(nside)
    weighted_values = np.zeros(num_healpixels)
    weights = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for (id, k) in zip(ids, kappa):
        weighted_values[id] += k
        weights[id] += 1.0
    return np.divide(weighted_values, weights, out=np.zeros_like(weights, dtype=float), where=weights!=0.0)
    
def glimpse_lattice_points_to_healpix_map(glimpse_output_file, nside, do_nest):

    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_file)
    
    num_healpixels = hp.nside2npix(nside)
    ret = np.zeros(num_healpixels)
    ids = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    for id in ids:
        ret[id] = 1.0
    return ret

# Cleanup ra_rad so that it lies in [0, 2pi)
def cleanup_ra_rad(ra_rad):
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

    num_healpixels = hp.nside2npix(nside)
    values = np.zeros(num_healpixels)

    index = hp.ang2pix(nside, ra, dec, do_nest, lonlat=True)
    values[index] = 1.0
    
    return values
    
    
def is_in_standard_cutout(ra, dec, side_in_degrees):
    return np.where(np.logical_and((np.abs(ra - 180.0) < 0.5*side_in_degrees), (np.abs(dec - 0.0) < 0.5*side_in_degrees)))
    
    
    
# return 0 if the absolute difference is small, else return percentage difference.    
def comparison_function(x, y, tolerance):
    r = np.divide((x-y), y, out=np.zeros_like(y, dtype=float), where=y!=0.0)
    r[np.where(np.abs(x-y)<=tolerance)] = 0.0
    return r
    
def write_healpix_array_to_file(a, filename, nest):

    if hp.__version__ == "1.10.3":
        # Suppress deprecation warning message about clobber/overwrite
        warnings.simplefilter('ignore', AstropyDeprecationWarning)

    # print("Writing {}".format(filename))
    hp.write_map(filename, a, nest, overwrite=True)
    
    
def plot_several_healpix_maps():

    nside = 1024 # Needed only in the 'other maps' section...

    maps = []
    titles = []

    # 1. maps from healpix files
    path = "/share/splinter/ucapwhi/glimpse_project/comparison/"
    #filenames = ["20201031_hsc_lambda_{}/".format(i, i) for i in [1,2,3,4,5]]
    filenames = ["glimpse.merged.values.fits", "202005_Mcal.lambda_3.glimpse.merged.values.fits", "Mcal_0.2_1.3.signal_shear.glimpse.merged.values.fits"]
    #filenames = [f + "glimpse.merged.values.fits" for f in filenames]
    
    weight_maps = []
    for i in weight_maps:
        # Also show weights
        filenames.append(filenames[i].replace(".values", ".weights"))
        
    masked_maps = []
    for i in masked_maps:
        filenames[i] = filenames[i].replace("_values", "_values_masked")
    
    
    for f in filenames:
        print("Using file {}".format(f))
        maps.append(hp.read_map(path + f, verbose=False))
        titles.append(f.replace("/glimpse.merged.values.fits","").replace("Buzzard_192_signal/truth.values.fits", "Truth").replace("Buzzard_192_signal/", "Buzzard_192_signal_lambda_3/").replace("202005", "Weights"))
        #titles.append(f)
        
    # 2. other maps
        
    if False:
        # Also plot some glimpse input data
        f = "Mcal_0.2_1.3.1689.cat.fits"
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
    
    if False:
        # Also plot true kappa from a simulation
        f = "../../data/hsc_catalogs_weighted_all_bins.fits"
        print("Using file {}".format(f))
        maps.append(kappa_map_from_simulation_catalogue(path + f, nside))
        titles.append(f + " simulation true kappa")
        

    # To smooth a map, include in maps_to_smooth the index of the map (i.e. the map's index in the array 'maps') 
    maps_to_smooth = []
    for i in maps_to_smooth:
        one_arcmin_in_radians = 0.000290888
        smoothing_scale_in_arcmin = 15.0
        maps[i] = hp.smoothing(maps[i], sigma = smoothing_scale_in_arcmin * one_arcmin_in_radians)
        titles[i] += " smoothed at {} arcmin".format(smoothing_scale_in_arcmin)
    
    if True:
        # Show diff between maps[0] and maps[1]
        if False:
            tolerance = 0.0009
            percentage_diff = comparison_function(maps[0], maps[1], tolerance)
            maps.append(percentage_diff)
            titles.append("Percentage difference")
        else:
            maps.append(maps[1] - maps[0])
            titles.append("Diff 1-0")
            maps.append(maps[2] - maps[1])
            titles.append("Diff 2-1")
        
        
    maps_to_flip_sign = []
    for i in maps_to_flip_sign:
        maps[i] = -1.0 * maps[i]
        titles[i] += " TIMES -1"
    
    plt.figure(figsize=(12,7))
    cmap=plt.get_cmap('inferno')
    for this_map, this_title, i in zip(maps, titles, range(len(maps))):
        
        print(np.amin(this_map), np.amax(this_map))
        
        if False:
            # Turn this on when dealing with old maps that indicated 'in masked area' using a value
            # of zero instead of hp.UNSEEN.
            this_map = np.where(np.abs(this_map)<1e-7,  hp.UNSEEN, this_map)
        
        if False:
            hp.mollview(this_map, title=this_title, sub=(1,len(maps),i+1), cmap=cmap)
        elif False:
            #rot = (-3.46154, -51.25581, 0.0) # (pixel 2759 when nside = 16)
            rot = (70.0, -50.0, 0.0)
            hp.gnomview(this_map, rot=rot, title=this_title, reso=5.0, xsize=400, sub=(1,len(maps),i+1)) #, min=-0.04, max=0.06
        elif True:
            rot = (40.0, -30.0, 0.0)
            hp.orthview(this_map, rot=rot, title=this_title, xsize=400, badcolor="grey", half_sky=True, sub=(1,len(maps),i+1), cmap=cmap, max=0.039, min=-0.016) # max=0.1, min=-0.1, 
        
        hp.graticule(dpar=30.0)
    
    if True:
        plt.savefig("./plot.png")
    plt.show()
    
    
def plot_several_glimpse_outputs():

    cmap = plt.get_cmap('inferno')
    fig = plt.figure(figsize=(12,12))
    num_rows_of_subplots = 1
    
    
    titles = []
    
    experiment_names = ["A", "D", "F"]
    
    if False:
        path = "/share/splinter/ucapwhi/glimpse_project/experiments/weight/"
        pixel_id = 2277
        filenames = ["{}/{}.glimpse.out.fits".format(experiment_name, pixel_id) for experiment_name in experiment_names]
    elif False:
        path = "/share/splinter/ucapwhi/glimpse_project/runs/Mcal_signal/"
        filenames = ["{}.glimpse.out.fits".format(pixel_id) for pixel_id in [2277, 2278, 2279, 2280, 2281, 2282]]
    else:
        path = "/share/splinter/ucapwhi/glimpse_project/runs/"
        filenames = ["{}/2277.glimpse.out.fits".format(sd) for sd in ["Mcal_signal", "202005_Mcal_lambda_3"]]
        
    kappa_arrays = [glimpse_output_as_array(path + f) for f in filenames]
    titles = [f.replace(".glimpse.out.fits", "") for f in filenames]
    
    if False:
        # Show diffs
        for (i, j) in [(0,1), (0,2), (1,2)]:
            kappa_arrays.append(kappa_arrays[i] - kappa_arrays[j])
            titles.append("{} - {}".format(experiment_names[i], experiment_names[j]))
    
    (vmin, vmax) = (10e10, -10e10)
    for kappa_array in kappa_arrays:
        vmin = min(vmin, np.amin(kappa_array))
        vmax = max(vmax, np.amax(kappa_array))
    
    for (kappa_a, t, i) in zip(kappa_arrays, titles, range(len(kappa_arrays))):
        ax = fig.add_subplot(num_rows_of_subplots, len(kappa_arrays)/num_rows_of_subplots, i+1)
        im = ax.imshow(kappa_a, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.title.set_text(t)
        plt.colorbar(im, ax=ax, orientation="horizontal")
    
    plt.show()
        



    
    
    
# 'experiment' will be a string such as "B".

def experiment_label(experiment):
    if experiment == "A":
        return "$\lambda$=3"
    if experiment == "B":
        return "nreweights=5"
    if experiment == "C":
        return "$\lambda$=5"
    if experiment == "D":
        return "nscales=3"
    if experiment == "E":
        return "$\lambda$=2"
    if experiment == "F":
        return "niter=1000"
    if experiment == "G":
        return "nrandom=2000"
    if experiment == "H":
        return "in=100,out=80"
    if experiment == "I":
        return "nscales=5"
    if experiment == "J":
        return "nscales=9"
    if experiment == "K":
        return "$\lambda$=4"
    if experiment == "L":
        return "$\lambda$=1"
    if experiment == "M":
        return "in=80,out=60"
    if experiment == "N":
        return "flip_e2=F"
        



def experiment_results_filename(experiment):
    return "/share/splinter/ucapwhi/glimpse_project/experiments/{}/glimpse.merged.values.fits".format(experiment)


def Buzzard_reduced_truth_map_filename():
    return "/share/splinter/ucapwhi/glimpse_project/experiments/truth/true_kappa.fits"
    
    
def pixel_histograms(list_of_maps, list_of_titles):

    fix_max_to_unity = False
    
    plt.yscale('log')
    for (m, t) in zip(list_of_maps, list_of_titles):
        (h, b) = np.histogram(m, bins=100, range=(-0.04, 0.06))
        plt.step(0.5*(b[1:]+b[:-1]), (h/np.max(h) if fix_max_to_unity else h), where="mid", label=t)
    plt.legend()
    plt.xlabel("$\kappa$")
    plt.ylabel("Normalised Pixel Count" if fix_max_to_unity else "Pixel Count")
    if True:
        plt.savefig("./plot.png")
    plt.show()
        
def show_pixel_histograms():

    list_of_maps = []
    list_of_titles = []
    filter = Buzzard_reduced_truth_filter()
    simulation_truth = hp.read_map(Buzzard_reduced_truth_map_filename(), verbose=False)[filter]
    list_of_maps.append(simulation_truth)
    list_of_titles.append("truth")
    for experiment in ["E", "A", "C"]:
        experimental_result = hp.read_map(experiment_results_filename(experiment), verbose=False)[filter]
        list_of_maps.append(experimental_result)
        list_of_titles.append(experiment_label(experiment))
    pixel_histograms(list_of_maps, list_of_titles)
    
    
def filter_zeros(a):
    filter = np.where(np.abs(a) > 1e-15)
    return a[filter]
    

def compare_two_pixel_histograms():

    directory = '/share/splinter/ucapwhi/glimpse_project/runs/'
    filenames = ['Mcal_signal/glimpse.merged.values.fits', '202005_Mcal_lambda_3/glimpse.merged.values.fits']
    list_of_maps = [filter_zeros(hp.read_map(directory + filename, verbose=False)) for filename in filenames]
    titles = ["Original", "New with weights"]
    
    print([len(m) for m in list_of_maps])
    pixel_histograms(list_of_maps, titles)

    
    

def experiment_tests():
    
    filter = Buzzard_reduced_truth_filter()
    
    print("{}\t\t\t{}\t{}\t{}".format("Run", "Pearson r", "Pixel res", "Var ratio"))

    simulation_truth = hp.read_map(Buzzard_reduced_truth_map_filename(), verbose=False)[filter]
    for experiment in ["L", "E", "A", "K", "C", "N"]:
        experimental_result = hp.read_map(experiment_results_filename(experiment), verbose=False)[filter]
        corr = np.corrcoef(np.column_stack((simulation_truth, experimental_result)), rowvar=False)[0,1] # Equation 32 in 1801.08945
        pixel_residual = np.sqrt(np.mean((experimental_result-simulation_truth)**2)) # Equation 34 in 1801.08945
        var_ratio = np.var(experimental_result)/np.var(simulation_truth) # Equation 36 in 1801.08945
        print("{}\t\t{:10.5f}\t{:10.5f}\t{:10.5f}".format(experiment_label(experiment), corr, pixel_residual, var_ratio))
        
def Buzzard_reduced_truth_filter():

    simulation_truth = hp.read_map(Buzzard_reduced_truth_map_filename(), verbose=False)
    return np.where(np.abs(simulation_truth) > 1e-12)
        

def corr_graph():

    from chainconsumer import ChainConsumer

    filter = Buzzard_reduced_truth_filter()
    
    c = ChainConsumer()
    
    list_of_chains = []
    list_of_names = []
    
    simulation_truth = hp.read_map(Buzzard_reduced_truth_map_filename(), verbose=False)[filter]
    list_of_chains.append(simulation_truth)
    list_of_names.append("truth")
    
    for experiment in ["A", "B"]:
        experiment_result = hp.read_map(experiment_results_filename(experiment), verbose=False)[filter]
        list_of_chains.append(experiment_result)
        list_of_names.append(experiment)
        
    c.add_chain(list_of_chains, parameters=list_of_names)
    c.configure(smooth=False, colors="#673AB7")

    fig = c.plotter.plot(display=True)
    


def clean_up_edges(input_catalogue_filename, ra_name, dec_name, map_to_be_masked):

    nside = hp.get_nside(map_to_be_masked)
    
    (ra, dec) = get_from_fits_file(input_catalogue_filename, [ra_name, dec_name])
    
    # Create a mask that is zero in healpixels where there are no source galaxies and one elsewhere.
    ids = hp.ang2pix(nside, ra, dec, nest=False, lonlat=True)
    mask = np.zeros(hp.nside2npix(nside))
    for id in ids:
        mask[id] = 1.0 
    
    return np.where(mask == 0.0, hp.UNSEEN, map_to_be_masked)
    
    



def sphere_to_tangent_plane_mapping(ra_centre, dec_centre, ra, dec):
    # See LW's notes p. GL90
    ra_centre_rad = np.radians(ra_centre)
    dec_centre_rad = np.radians(dec_centre)
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    x = np.sin(ra_rad - ra_centre_rad) / (np.cos(dec_centre_rad) * np.cos(ra_rad - ra_centre_rad) + np.tan(dec_rad) * np.sin(dec_centre_rad))
    epsilon = np.sin(dec_rad) * np.cos(dec_centre_rad) - np.cos(dec_rad) * np.sin(dec_centre_rad) * np.cos(ra_rad - ra_centre_rad)
    y = epsilon * np.sqrt((1.0 + x**2) / (1.0 - epsilon**2))
    return (x, y)
    

    
    
def sphere_to_tangent_plane_mapping_test_harness():

    pixel_id = 2218
    ra_centre = 56.25
    dec_centre = -27.27961273597809

    glimpse_output_file = "/share/splinter/ucapwhi/glimpse_project/output/Buzzard_192." + str(pixel_id) + ".glimpse.out.fits"
    (ra, dec, kappa) = get_glimpse_output_data(glimpse_output_file)
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
    return np.minimum(np.maximum(x, lower), upper)


def index_into_glimpse_array_test_harness():
    p_shape = [6, 6]
    
    glimpse_lattice_spacing = 1.2
    
    x = 0.001
    y = 0.001
    print(x,y)
    print(bound_between(p_shape[0]//2 + np.floor(x / glimpse_lattice_spacing), 0, p_shape[0]-1) + \
                p_shape[0] * bound_between(p_shape[1]//2 - np.floor(y / glimpse_lattice_spacing) - 1, 0, p_shape[1]-1))
    
    

def ra_dec_to_healpixel_id_test_harness():
    ra = 0.4464
    dec = -0.2678
    nest = False
    nside = 1024
    
    print(hp.ang2pix(nside, ra, dec, nest, lonlat=True))
    
    
    

# No error if directory already exists
# See https://stackoverflow.com/questions/273192
def create_directory_if_necessary(directory_name):
    try:
        os.makedirs(directory_name)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

# percentage should be an integer e.g. 60 for 60%. Enter a negative number to clean-up at the end.
def update_percentage_bar(percentage):
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

	long1 = np.radians(ra1)
	long2 = np.radians(ra2)
	lat1 = np.radians(dec1)
	lat2 = np.radians(dec2)
	half_dlong = (long2 - long1) / 2.0
	half_dlat = (lat2 - lat1) / 2.0
	a = np.sin(half_dlat)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(half_dlong)**2
	return 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))
    
def angular_separation_test_harness():
    print(np.degrees(angular_separation(0.0, 0.0, 5.0, 35.0)))


# Same rules as for angular_separation
def angular_separation_fast(sin_ra1, cos_ra1, sin_dec1, cos_dec1, ra2, dec2):

    
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
    
    

def kappa_map_from_simulation_catalogue(simulation_catalogue, nside):

    list_of_field_names = ("ra", "dec", "k_orig")
    (ra, dec, kappa) = get_from_fits_file(simulation_catalogue, list_of_field_names)
    
    num_pixels = hp.nside2npix(nside)
    k_total = np.zeros(num_pixels)
    k_count = np.zeros(num_pixels)
    for (i, k) in zip(hp.ang2pix(nside, ra, dec, False, True), kappa):
        k_total[i] += k
        k_count[i] += 1
        
    # The following division code treats 0/0 as hp.UNSEEN. From https://stackoverflow.com/questions/26248654.
    kmap = np.divide(k_total, k_count, out = np.full(num_pixels, hp.UNSEEN, dtype=float), where=k_count!=0.0)
    return kmap
    
def kappa_map_from_simulation_catalogue_test_harness():
    nside = 1024
    simulation_catalogue = "/share/splinter/ucapwhi/glimpse_project/data/hsc_catalogs_weighted_all_bins.fits"
    k_map = kappa_map_from_simulation_catalogue(simulation_catalogue, nside)
    print(k_map)
    
    

######################### Start of create_cutouts code #########################

# Helper function - see https://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python
def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# Will process ids in the range [ids_to_process_start, ids_to_process_end). Set ids_to_process_end to -1 to mean "to the end"
def create_cutouts(input_catalogue, raname, decname, shear_names, other_field_names, nside, cutout_side_in_degrees, job_control, output_directory, output_file_root):

    
    
    create_directory_if_necessary(output_directory)
    assert hp.isnsideok(nside), "nside must be a valid Healpix nside"
    all_field_names = [raname, decname]
    if other_field_names:
        all_field_names.extend(other_field_names.split(","))
    if shear_names:
        parsed_shear_names = shear_names.split(",")
        assert (len(parsed_shear_names)%2==0), "shear_names should be omitted or be a comma-separated list of (an even number of) field names"
        all_field_names.extend(parsed_shear_names)
    
    
    print('Processing ' + input_catalogue + '...')
    data = Table.read(input_catalogue, format='fits') # Other formats such as cvs are also supported here...
    
    for field_name in all_field_names:
        # field name might be e.g. W.NORMALIZE or DUMMY.1; need to handle these.
        parsed_field_name = field_name.split(':')
        field_name_base = parsed_field_name[0]
        field_name_suffix = "" if len(parsed_field_name) == 1 else parsed_field_name[1]
        if field_name_base[0:5].upper() == "DUMMY":
            assert isfloat(field_name_suffix), "Cannot parse field name {}".format(field_name)
            data[field_name_base] = data[data.colnames[0]]*0.0 + float(field_name_suffix)
        else:
            assert field_name_base in data.colnames, "{} is not a column in the input file {}".format(field_name, input_catalogue)
            # Add additional transform rules here, if you want.
            if field_name_suffix.upper() == "NORMALIZE" or field_name_suffix.upper() == "NORMALISE":
                this_field_values = data[field_name_base]
                data[field_name_base] = this_field_values / np.average(this_field_values)
            elif field_name_suffix.upper() == "FLIPSIGN":
                data[field_name_base] = -1 * data[field_name_base]
            else:
                assert field_name_suffix == "", "Cannot parse field name {}".format(field_name)
        
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

    for healpix_id in array_slice_from_job_control_string(range(num_healpixels), job_control):

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
                    break
            
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
                        output_filename = os.path.join(output_directory, output_file_root).format(str(healpix_id).zfill(num_digits))
                        
                        field_names = []
                        data_columns = []
                        
                        field_names.append(raname)
                        data_columns.append(ra_standardised_filter1[filter2])
                                    
                        field_names.append(decname)
                        data_columns.append(dec_standardised_filter1[filter2])
                        
                        
                        if shear_names:
                            parsed_shear_names = shear_names.split(",")
                            for (shear_name_1, shear_name_2) in zip(parsed_shear_names[0::2], parsed_shear_names[1::2]):
                                shear_name_1 = shear_name_1.split(':')[0]
                                shear_name_2 = shear_name_2.split(':')[0]
                                (shear1, shear2) = rotate_shear_45_degrees(data[shear_name_1], data[shear_name_2])
                                
                                field_names.append(shear_name_1)
                                data_columns.append(shear1[joint_filter])
                                
                                field_names.append(shear_name_2)
                                data_columns.append(shear2[joint_filter])
                        
                        for f in other_field_names.split(","):
                            f = f.split(':')[0]
                            field_names.append(f)
                            data_columns.append(data[f][joint_filter])
                           
                        write_to_fits_file(output_filename, field_names, data_columns, True)



def create_cutouts_caller(directory, job_control):
    
    config = configparser.ConfigParser()
    ini_file_name = os.path.abspath(os.path.join(directory, "control.ini"))
    config.read(ini_file_name)

    input_catalogue = config["create_cutouts"].get("input_catalogue")
    ra_name = config["create_cutouts"].get("ra_name", "RA")
    dec_name = config["create_cutouts"].get("dec_name", "DEC")
    shear_names = config["create_cutouts"].get("shear_names")
    other_field_names = config["create_cutouts"].get("other_field_names")
    nside = int(config["create_cutouts"].get("nside", "16"))
    cutout_side_in_degrees = float(config["create_cutouts"].get("cutout_side_in_degrees", "16"))
    
    output_directory = os.path.dirname(os.path.abspath(ini_file_name))
    output_file_root = "{}.cat.fits"
    
    create_cutouts(input_catalogue, ra_name, dec_name, shear_names, other_field_names, nside, cutout_side_in_degrees, job_control, output_directory, output_file_root)


######################### End of create_cutouts code #########################






######################### Start of glimpse_caller code #########################

def glimpse_caller(directory, job_id):

    config = configparser.ConfigParser()
    ini_file_name = os.path.abspath(os.path.join(directory, "control.ini"))
    config.read(ini_file_name)
    
    output_directory = os.path.dirname(ini_file_name)
    exe_file_name = config["project"].get("glimpse_executable")
    
    glimpse_cat_file_pattern = os.path.join(output_directory, "*.cat.fits")
    id_list = [int((os.path.basename(f)).split(".")[0]) for f in glob.glob(glimpse_cat_file_pattern)]
    id_list.sort()
    this_healpix_pixel_id = id_list[int(job_id)]
    
    print("Healpixel id = {}...".format(this_healpix_pixel_id), flush=True)
    
    this_healpix_id_as_string = str(this_healpix_pixel_id).zfill(4)
    
    cat_file_name = os.path.join(output_directory, this_healpix_id_as_string + ".cat.fits")
    out_file_name = os.path.join(output_directory, this_healpix_id_as_string + ".glimpse.out.fits")

    subprocess.run([exe_file_name, ini_file_name, cat_file_name, out_file_name])

    print("Finished.")
    
######################### End of glimpse_caller code #########################


######################### Start of merge code #########################


# Satisfies f(x0) = 1, f(x1) = 0, f'(x0) = 0, f'(x1) = 0
def smoothed_step_function(x, x0, x1):
    x_scaled = (float(x) - x0) / (x1 - x0)
    return 2.0 * x_scaled**3 - 3.0 * x_scaled**2 + 1.0


# Returns a 1D array of length 'width'. The first outer_border values are 0, and the next
# inner_border-outer_border values are smoothly intermediate between 0 and 1. Then there
# is a series of 1s. The output is symmetric about the midpoint.
# See LW's notes p. GL105.
def one_axis_weight_function(width, outer_border, inner_border):
    assert (outer_border > 0), "outer_border must be positive in one_axis_weight_function"
    assert (outer_border <= inner_border), "outer_border must not exceed inner_border in one_axis_weight_function"
    assert (width > 2 * inner_border), "inner_border is too large in one_axis_weight_function"
    ret = np.zeros(width)
    for i in range(width):
        if i < outer_border:
            pass # ret[i] = 0.0, but we don't need to do this as zero is the default value.
        elif i < inner_border:
            ret[i] = smoothed_step_function(i, inner_border, outer_border - 1)
        elif i < width - inner_border:
            ret[i] = 1.0
        elif i < width - outer_border:
            ret[i] = smoothed_step_function(i, width - inner_border - 1, width - outer_border)
        else:
            pass # ret[i] = 0.0, but we don't need to do this as zero is the default value.
    return ret
    
def one_axis_weight_function_test_harness():
    oawf = one_axis_weight_function(16, 2, 1)
    for x in oawf:
        print(x)



# Returns a 1D array - this is the flattened version of the 2D array that is the
# geometric average of two one_axis_weight_functions, one on each of the two axes.
# Input 'shape' should be the desired shape of the before-being-flattened 2D array.
def two_axis_weight_function(shape, outer_border, inner_border):
    x = one_axis_weight_function(shape[0], outer_border, inner_border)
    y = one_axis_weight_function(shape[1], outer_border, inner_border)
    return np.sqrt(x[np.newaxis].T * y).reshape(-1)

def two_axis_weight_function_test_harness():
    d1 = 100
    d2 = 35
    outer_border = 7
    inner_border = 14
    a = two_axis_weight_function((d1, d2), outer_border, inner_border)
    a = np.reshape(a, (d1, d2))
    plt.imshow(a)
    plt.show()

   
    

def merge(input_file_spec, outer_border, inner_border, output_file_root, intermediate_nside, output_nside, cutouts_nside, apply_galaxy_mask, input_catalogue, ra_name, dec_name, job_control):

    RA = 0
    DEC = 1
    
    # Currently hard-coded to output in RING format; can change this if desired.
    output_nest = False
    
    output_debugging_info = False
    g_index_special = 6318085 # Dump debugging info for this pixel
    if output_debugging_info:
        print("Debugging point: RA={}; DEC={}; index={}".format(hp.pix2ang(intermediate_nside, g_index_special, output_nest, lonlat=True)[0], hp.pix2ang(intermediate_nside, g_index_special, output_nest, lonlat=True)[1], g_index_special))
    
    assert (outer_border <= inner_border), "inner_border must not be larger than outer_border"

    num_healpixels = hp.nside2npix(intermediate_nside)

    # "g_" means "global" i.e. for the whole sphere. Indices into these arrays must repect the 'output_nest' parameter.
    g_weighted_values = np.zeros(num_healpixels)
    g_weights = np.zeros(num_healpixels)
    g_centres = hp.pix2ang(intermediate_nside, np.arange(num_healpixels), output_nest, True)
    
    g_centres_ra_rad = np.radians(g_centres[RA])
    g_centres_dec_rad = np.radians(g_centres[DEC])
    g_sin_centres_ra = np.sin(g_centres_ra_rad)
    g_cos_centres_ra = np.cos(g_centres_ra_rad)
    g_sin_centres_dec = np.sin(g_centres_dec_rad)
    g_cos_centres_dec = np.cos(g_centres_dec_rad)
    
    
    list_of_input_files = glob.glob(input_file_spec)
    list_of_input_files.sort()
    list_of_input_files = array_slice_from_job_control_string(list_of_input_files, job_control)
    
    num_input_files = len(list_of_input_files)

    assert (num_input_files > 0), "No input files found"
    
    is_first_file = True

    for (filename, file_index) in zip(list_of_input_files, range(num_input_files)):

        if not output_debugging_info:
            print("Processing {} of {}: {}...".format(file_index, num_input_files, filename))
            
        # "p_" means "for this one picture"
        (p_ra, p_dec, p_val) = get_glimpse_output_data(filename)
        num_glimpse_pixels_on_side = int(np.sqrt(p_ra.shape[0]))
        p_shape = (num_glimpse_pixels_on_side, num_glimpse_pixels_on_side)

        
        if is_first_file:
            # For speed these data can be cached
            # (as they will be the same for all cutouts (as all cutouts are in standard position)).
            p_ra_rad = np.radians(p_ra)
            p_dec_rad = np.radians(p_dec)
            p_cos_dec = np.cos(p_dec_rad)
            p_x = p_cos_dec * np.cos(p_ra_rad)
            p_y = p_cos_dec * np.sin(p_ra_rad)
            p_z = np.sin(p_dec_rad)
            is_first_file = False
       
        
        # Here we are assuming the standard filename convention, so we can get
        # the coarse healpixel id from a certain location in the file name.
        # File name will be of the form (various}/XXXX.glimpse.out.fits
        #where XXXX is a healpix id (nside = cutouts_nside)
        assert (filename[-17:] == ".glimpse.out.fits"), "Filename {} was not of the expected format (should end with .glimpse.out.fits)".format(filename)
        cutout_healpix_id = int(os.path.basename(filename).split('.')[0])
        (p_centre_ra, p_centre_dec) = hp.pix2ang(cutouts_nside, cutout_healpix_id, nest=False, lonlat=True)
        
        (p_ra, p_dec) = from_standard_position_fast(p_x, p_y, p_z, p_centre_ra, p_centre_dec)

        if any(np.isnan(p_val)):
            print("\t\tRejected as it contains NaNs")
        else:
            # The radius of the glimpse lattice (on the sphere).
            radius_rad = angular_separation(p_centre_ra, p_centre_dec, p_ra[0], p_dec[0])

            # The angular separation between each healpixel and the centre of the current picture.
            #g_ang_sep_rad = angular_separation(g_centres[RA], g_centres[DEC], p_centre_ra, p_centre_dec)
            g_ang_sep_rad = angular_separation_fast(g_sin_centres_ra, g_cos_centres_ra, g_sin_centres_dec, g_cos_centres_dec, p_centre_ra, p_centre_dec)
            
            
            # Only include those healpixels that are sensibly close to the centre of the lattice.
            # This is a speedup, and also prevents problems with the 'project sphere to tangent plane'
            # routine (which would fail e.g. if the point on the sphere were on the wrong hemisphere).
            g_include = np.where(g_ang_sep_rad <= radius_rad)
            
            g_centres_ra_filtered = g_centres[RA][g_include]
            g_centres_dec_filtered = g_centres[DEC][g_include]
            g_indices_filtered = np.arange(num_healpixels)[g_include]
            
            # Map the (relevant) healpixel centres to the tangent plane for this picture.
            (g_x_filtered, g_y_filtered) = sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, g_centres_ra_filtered, g_centres_dec_filtered)
            
            # We know that the glimpse lattice points will map to a perfectly regular lattice on the tangent plane,
            # so no need to actually do this computation. However, we do need to know the spacing of the glimpse lattice
            # points in the tangent plane. So do a simple experiment to discover this.
            glimpse_lattice_spacing = sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra[1], p_dec[1])[0] - sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra[0], p_dec[0])[0]
            
            if False:
                (p_lattice_x, p_lattice_y) = sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra, p_dec)
                plt.scatter(p_lattice_x, p_lattice_y)
                plt.scatter(p_lattice_x[0], p_lattice_y[0], c="red", s=10)
                plt.scatter(p_lattice_x[30], p_lattice_y[30], c="blue", s=10)
                plt.show()
            
                        
            # For each relevant healpixel, find the index of the nearest lattice point. See LW's notes p. GL94.
            # Note that this function hard-codes specific information about how the points in the glimpse lattice
            # are ordered in the glimpse output. If you wanted to remove the dependency on this information,
            # could replace this code with a routine to 'find the closest glimpse lattice point'. This would be
            # slower but would be more future-proof against possible changes in the glimpse layout.
            
            # Ha! These changes came sooner than you expected. See p. GL121
            g_x_filtered_rotated = (g_x_filtered - g_y_filtered) * np.sqrt(2.0) / 2.0
            g_y_filtered_rotated = (g_x_filtered + g_y_filtered) * np.sqrt(2.0) / 2.0
            glimpse_lattice_spacing *= np.sqrt(2.0)
            
            g_index_into_glimpse_array = (bound_between(p_shape[0]//2 + np.floor(g_x_filtered_rotated / glimpse_lattice_spacing), 0, p_shape[0]-1) + \
                p_shape[0] * bound_between(p_shape[1]//2 + np.floor(g_y_filtered_rotated / glimpse_lattice_spacing) - 1, 0, p_shape[1]-1)).astype(int)
                
            p_weight = two_axis_weight_function(p_shape, outer_border, inner_border)
            
            # See 'Assigning values to indexed arrays' at https://www.numpy.org/devdocs/user/basics.indexing.html
            # to see why the naive approach to what follows fails:
            for (p_index, g_index) in zip(g_index_into_glimpse_array, g_indices_filtered):
                v = p_val[p_index]
                w = p_weight[p_index]
                
                g_weighted_values[g_index] += (w * v)
                g_weights[g_index] += w
                
                if output_debugging_info:
                    if g_index == g_index_special:
                        print(g_index, p_index, v, w, filename)
       

    # The following division code treats 0/0 as 0. From https://stackoverflow.com/questions/26248654.
    g_average_values = np.divide(g_weighted_values, g_weights, out=np.zeros_like(g_weights, dtype=float), where=g_weights!=0.0)
    
    if output_debugging_info:
        print("g_weighted_values[{}]={}".format(g_index_special, g_weighted_values[g_index_special]))
        print("g_weights[{}]={}".format(g_index_special, g_weights[g_index_special]))
        print("g_average_values[{}]={}".format(g_index_special, g_average_values[g_index_special]))
        
    # Downgrade to final output resolution
    if output_nside != intermediate_nside:
        g_average_values = hp.ud_grade(g_average_values, output_nside)
        g_weights = hp.ud_grade(g_weights, output_nside)
        
    # Mask out pixels where there were no input galaxies
    if apply_galaxy_mask:
        print("Applying galaxy mask...")
        g_average_values = clean_up_edges(input_catalogue, ra_name, dec_name, g_average_values)
        
    # Save data to healpix format file
    for (filename_part, array) in zip(["values", "weights"], [g_average_values, g_weights]):
        dat_filename = output_file_root + "glimpse.merged." + filename_part + ".fits"
        print("Writing {}...".format(dat_filename))
        write_healpix_array_to_file(array, dat_filename, output_nest)



def merge_caller(directory, job_control):
    
    
    config = configparser.ConfigParser()
    ini_file_name = os.path.abspath(os.path.join(directory, "control.ini"))
    config.read(ini_file_name)
    
    outer_border = int(config["merge"].get("outer_border"))
    inner_border = int(config["merge"].get("inner_border"))
    intermediate_nside = int(config["merge"].get("intermediate_nside"))
    output_nside = int(config["merge"].get("output_nside"))
    
    input_file_spec = os.path.join(directory, "*.glimpse.out.fits")
    output_file_root = os.path.join(directory, "")
    
    # Info needed to parse the output glimpse file names
    cutouts_nside = int(config["create_cutouts"].get("nside"))
    
    # Info needed for masking
    apply_galaxy_mask = config['merge'].getboolean('apply_galaxy_mask?', fallback=False)
    input_catalogue = config["create_cutouts"].get("input_catalogue")
    ra_name = config["create_cutouts"].get("ra_name", "RA")
    dec_name = config["create_cutouts"].get("dec_name", "DEC")

    
    merge(input_file_spec, outer_border, inner_border, output_file_root, intermediate_nside, output_nside, cutouts_nside, apply_galaxy_mask, input_catalogue, ra_name, dec_name, job_control)

######################### End of merge code #########################







######################### Start of status code #########################


# Helper function for 'status' routine
def report_whether_file_exists(file_description, file_name):
    print("File {} {} '{}'".format(("exists:" if os.path.isfile(file_name) else "DOES NOT exist: no"), file_description, file_name))

def plural_suffix(count):
    return ("" if count==1 else "s")


def status(directory):
    
    # date and time
    print("Status as of {}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    
    # inifile present?
    ini_file_name = os.path.abspath(os.path.join(directory, "control.ini"))
    report_whether_file_exists("ini file", ini_file_name)
    
    # cutouts present?
    cutouts_filespec = os.path.join(directory, "*.cat.fits")
    num_cutouts = len(glob.glob(cutouts_filespec))
    if num_cutouts > 0:
        print("Files exist: {} cutout file{}".format(num_cutouts, plural_suffix(num_cutouts)))
    else:
        print("Files DO NOT exist: no cutout files")
    
    # glimpse outputs present?
    glimpse_filespec = os.path.join(directory, "*.glimpse.out.fits")
    glimpse_filelist = glob.glob(glimpse_filespec)
    if len(glimpse_filelist) > 0:
        file_size_dir = {} # key is file size; value is number of files with that size.
        for f in glimpse_filelist:
            f_size = os.path.getsize(f)
            if f_size not in file_size_dir:
                file_size_dir[f_size] = 0
            file_size_dir[f_size] += 1
        file_sizes_string = ""
        for (i, f_size) in zip(range(len(file_size_dir)), file_size_dir):
            file_sizes_string += ("" if i==0 else "; ") + "{} file{} of size {}".format(file_size_dir[f_size], plural_suffix(file_size_dir[f_size]), f_size)
        print("Files exist: {} glimpse output file{} ({})".format(len(glimpse_filelist), plural_suffix(len(glimpse_filelist)), file_sizes_string))
    else:
        print("Files DO NOT exist: no glimpse output files")
        
    # Merge files present?
    report_whether_file_exists("glimpse output values file", os.path.abspath(os.path.join(directory, "glimpse.merged.values.fits")))
    report_whether_file_exists("glimpse output weights file", os.path.abspath(os.path.join(directory, "glimpse.merged.weights.fits")))
    

def status_caller(directory):
    #directory = os.path.dirname(os.path.abspath(ini_file_name))
    status(directory)
    

######################### End of status code #########################


######################### Start of bins comparison code #########################


def comparison():

    files = ["202005_Mcal_lambda_3", "202005_Mcal_lambda_3_bin_0", "202005_Mcal_lambda_3_bin_1", "202005_Mcal_lambda_3_bin_2", "202005_Mcal_lambda_3_bin_3"]
    files = ["/share/splinter/ucapwhi/glimpse_project/runs/" + f + "/glimpse.merged.values.fits" for f in files]
    maps = [hp.read_map(f, verbose=False) for f in files]
    masked_maps = [np.where(m == hp.UNSEEN, 0, m) for m in maps]
    #joint_masked_map = 0.25* (masked_maps[1] + masked_maps[2] + masked_maps[3] + masked_maps[4])
    joint_masked_map = masked_maps[1]
    #plt.scatter(masked_maps[0], joint_masked_map)
    plt.hist2d(masked_maps[0], joint_masked_map, bins=500, cmax=10000000, range=[[-0.005, 0.005], [-0.005, 0.005]])
    plt.show()
    
    



######################### End of bins comparison code #########################








    

if __name__ == '__main__':

    # Test harnesses
    #sphere_to_tangent_plane_mapping_test_harness()
    #index_into_glimpse_array_test_harness()
    #ra_dec_to_healpixel_id_test_harness()
    #to_standard_position_test_harness()
    #from_standard_position_test_harness()
    #to_from_standard_position_test_harness()
    #array_slice_from_job_control_string_test_harness()
    #angular_separation_test_harness()
    #one_axis_weight_function_test_harness()
    #define_Buzzard_ten_percent_by_area_subset()
    #kappa_values_in_one_fine_pixel()
    #run_save_buzzard_truth()
    #plot_several_glimpse_outputs()
    #experiment_tests()
    #show_pixel_histograms()
    #corr_graph()
    #create_test_catalogue()
    #correct_one_shear_catalogue_caller()
    #redshift_histogram()
    #clean_up_edges()
    #shear_stdev()
    #num_files_in_directory()
    #to_from_standard_position_fast_test_harness()
    #angular_separation_fast_test_harness()
    #correct_one_shear_catalogue("/share/splinter/ucapwhi/glimpse_project/shear_sign_experiment/cat_3_0.fits")
    #glimpse_sign_experiment()
    #kappa_histogram()
    #show_glimpse_output_as_image()
    #compare_two_cutouts()
    #joint_filter_example()
    #append_random_shear_to_Buzzard()
    #set_weights_in_catalogue_file()
    #compare_two_pixel_histograms()
    #pickle_to_fits_caller()
    #pickle_test()
    #show_pickle_file_structure_caller()
    #clean_up_edges_caller()
    #plot_several_healpix_maps()
    comparison()
    pass
    
    
