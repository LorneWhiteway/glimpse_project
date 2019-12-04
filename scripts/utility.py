#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

        
def tester():
    import healpy as hp
    (ra ,dec) = hp.pix2ang(16, 2218, False, True)
    id = hp.ang2pix(1024, ra, dec, False, lonlat=True)
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
    
def fits_catalog_to_healpix_map(catalog_file_name, nside, do_nest):

    import healpy as hp
    import numpy as np

    (ra, dec) = get_fits_data(catalog_file_name, ["ra", "dec"])
    
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

# Cleanup ra_rad so that it luies in [0, 2pi)
def cleanup_ra(ra_rad):
    import numpy as np
    two_pi = 2.0 * np.pi
    ret = ra_rad
    ret[np.where(ret >= two_pi)] -= two_pi
    ret[np.where(ret < 0.0)] += two_pi
    return ret
    
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
    # See LW's notes p. 107.
    return (e2, -e1)


def to_standard_position(ra, dec, ra_centre, dec_centre):
    import numpy as np
    
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    ra_centre_rad = np.radians(ra_centre)
    dec_centre_rad = np.radians(dec_centre)
    
    # Move (ra_centre, dec_centre) to (0, 0)
    # See LW's book p. 78
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    sin_ra_diff = np.sin(ra_rad - ra_centre_rad)
    cos_ra_diff = np.cos(ra_rad - ra_centre_rad)
    ra_shifted_rad = np.arctan2(cos_dec * sin_ra_diff, cos_dec * np.cos(dec_centre_rad) * cos_ra_diff + sin_dec * np.sin(dec_centre_rad))
    dec_shifted_rad = np.arcsin(sin_dec * np.cos(dec_centre_rad) - cos_dec * np.sin(dec_centre_rad) * cos_ra_diff)
    
    # Rotate by 45 degrees
    (ra_shifted_rotated_rad, dec_shifted_rotated_rad) = rotate_45_degrees(ra_shifted_rad, dec_shifted_rad, 1.0)
    
    # Move (0, 0) to (pi, 0) and cleanup RA to be in [0, 2pi).
    ra_shifted_rotated_rad += np.pi
    ra_shifted_rotated_rad = cleanup_ra(ra_shifted_rotated_rad)
    
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
    ra_rotated_shifted_rad = cleanup_ra(ra_rotated_shifted_rad)
    
    ra_rotated_shifted = np.degrees(ra_rotated_shifted_rad)
    dec_rotated_shifted = np.degrees(dec_rotated_shifted_rad)
    
    return (ra_rotated_shifted, dec_rotated_shifted)
    

def to_standard_position_test_harness():

    import numpy as np
    import matplotlib.pyplot as plt
    
    ra_centre = 72.0
    dec_centre = -88.0
    
    #ra = np.linspace(0.0, 359.0, num = 360)
    #dec = ra * 0 + 1.0
    
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
    
    #ra = np.linspace(0.0, 359.0, num = 360)
    #dec = ra * 0 + 1.0
    
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
    
    nside = 1024

    maps = []
    titles = []

    # 1. maps from files
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    #filenames = ["Buzzard_192.nside" + str(nside) + "_merged_B_values.dat", "Buzzard_192.nside" + str(nside) + "_truth.dat"]
    #filenames = ["foo2_values.dat", "Buzzard_192.nside" + str(nside) + "_truth.dat"]
    #filenames = ["foo1_values.dat", "foo2_values.dat"]
    #filenames = ["foo2218_values.dat", "foo1440_values.dat"]
    filenames = ["foo2218_values.dat",]
    
    
    weight_maps = [0]
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
        
    if True:
        # Also plot some glimpse input data
        f = "Buzzard_192.2218.glimpse.cat.fits"
        maps.append(fits_catalog_to_healpix_map(path + f, nside, False))
        titles.append(f)
    
    if True:
        # Also plot glimpse output lattice
        f = "Buzzard_192.2218.glimpse.out.fits"
        maps.append(glimpse_output_to_healpix_map(path + f, nside*4, False))
        titles.append(f + " values")
        
    if True:
        # Also plot glimpse lattice points
        f = "Buzzard_192.2218.glimpse.out.fits"
        maps.append(glimpse_lattice_points_to_healpix_map(path + f, nside*4, False))
        titles.append(f + " lattice points")
        
    if True:
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
        tolerance = 0.0009 # 0.0012 is the best so far
        percentage_diff = comparison_function(maps[0], maps[1], tolerance)
        maps.append(percentage_diff)
        titles.append("Percentage difference")
    
    
    for map, title, i in zip(maps, titles, range(len(maps))):
        if False:
            #hp.mollview(map, fig=i, title=title)
            pass
        else:
            rot = (56.25, -27.2796127, 0.0)
            hp.gnomview(map, fig=i, rot=rot, title=title, reso=5.0, xsize=400)
        hp.graticule()
    
    plt.show()
    

def clean_up_edges():

    import healpy as hp
    import numpy as np
    import matplotlib.pyplot as plt
    
    do_nest = False
    
    path = "/share/splinter/ucapwhi/glimpse_project/output/"
    
    nside = 1024

    if False:
        old_filename = "Buzzard_192.nside" + str(nside) + "_merged_B_values.dat"
        new_filename = "Buzzard_192.nside" + str(nside) + "_merged_B_values_masked.dat"
        new_filename_png = "Buzzard_192.nside" + str(nside) + "_merged_B_values_masked.png"
    else:
        old_filename = "foo2_values.dat"
        new_filename = "foo2_values_masked.dat"
        new_filename_png = "foo2_values_masked.png"
    
    
    m = hp.read_map(path + old_filename)
    (ra, dec) = get_fits_data(buzzard_data_file_name(), ["RA", "DEC"])
    
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





def create_cutouts(input_catalogue, catformat, raname, decname, shear_names, other_field_names, nside, nest, cutout_side_in_degrees, output_directory, output_file_root):

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
        
    print('Creating output files...')
    old_percentage = -1
    loop_length = num_healpixels
    i = 0
    for healpix_id in np.arange(num_healpixels):
        print("Processing {} of {}".format(healpix_id, num_healpixels))
        
        # Progress bar
        percentage = int(100*(i+1)/float(loop_length))
        if percentage > old_percentage:
            #update_percentage_bar(percentage)
            old_percentage = percentage
        
        
        
        (ra_centre, dec_centre) = hp.pix2ang(nside, healpix_id, nest, lonlat=True)
        
        # Put the data in standard position relative to this centre point
        print("Standardising...")
        (ra_standardised, dec_standardised) = to_standard_position(ra, dec, ra_centre, dec_centre)
        
        # Which data lie in the standard cutout?
        print("Filtering...")
        filter = is_in_standard_cutout(ra_standardised, dec_standardised, cutout_side_in_degrees)
        
        if filter[0].size > 0:
            # Some galaxies to write
            print("Writing...")
            output_file_name = (output_directory + "/" + output_file_root).format(str(healpix_id).zfill(num_digits))
            
            column_info = []
            column_info.append(pyfits.Column(name=raname, format='D', array=ra_standardised[filter]))
            column_info.append(pyfits.Column(name=decname, format='D', array=dec_standardised[filter]))
            
            if shear_names:
                parsed_shear_names = shear_names.split(",")
                for (shear_name_1, shear_name_2) in zip(parsed_shear_names[0::2], parsed_shear_names[1::2]):
                    print(shear_name_1, shear_name_2)
                    (shear1, shear2) = rotate_shear_45_degrees(data[shear_name_1], data[shear_name_2])
                    column_info.append(pyfits.Column(name=shear_name_1, format='D', array=shear1[filter]))
                    column_info.append(pyfits.Column(name=shear_name_2, format='D', array=shear2[filter]))
            
            for f in other_field_names.split(","):
                column_info.append(pyfits.Column(name=f, format='D', array=data[f][filter]))

            tbhdu = pyfits.BinTableHDU.from_columns(column_info)
            tbhdu.writeto(output_file_name, overwrite=True)
            if False:
                update_percentage_bar(-1)
                sys.exit()
            
        
        i += 1
       
    update_percentage_bar(-1)


def create_cutouts_test_harness():

    input_catalogue = buzzard_data_file_name()
    catformat = "fits"
    raname = "RA"
    decname = "DEC"
    shear_names = "E1,E2"
    other_field_names = "true_z"
    nside = 16
    nest = False
    cutout_side_in_degrees = 16
    output_directory = "/share/splinter/ucapwhi/glimpse_project/output1"
    output_file_root = "Buzzard_192.{}.glimpse.cat.fits"

    create_cutouts(input_catalogue, catformat, raname, decname, shear_names, other_field_names, nside, nest, cutout_side_in_degrees, output_directory, output_file_root)

if __name__ == '__main__':

    #print_list_of_healpixels()
    #kappa_values_in_one_fine_pixel()
    #tester()
    #save_buzzard_truth()
    #clean_up_edges()
    #sphere_to_tangent_plane_mapping_test_harness()
    #index_into_glimpse_array_test_harness()
    #ra_dec_to_healpixel_id_test_harness()
    
    #plot_several_healpix_maps()
    #to_standard_position_test_harness()
    #from_standard_position_test_harness()
    #to_from_standard_position_test_harness()
    create_cutouts_test_harness()
