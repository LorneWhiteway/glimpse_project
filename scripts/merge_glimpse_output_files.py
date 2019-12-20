#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
    Merges output files from Glimpse into one healpix file.
    Author: Lorne Whiteway.
"""


#RA and Dec assumed to be in decimal degrees [0, 360] and [-90, 90].
def ra_dec_to_healpixelid(nside, nest, ra, dec):
    import healpy as hp
    return hp.ang2pix(nside, ra, dec, nest, lonlat=True)


# Satisfies f(x0) = 1, f(x1) = 0, f'(x0) = 0, f'(x1) = 0
def smoothed_step_function(x, x0, x1):
    x_scaled = (float(x) - x0) / (x1 - x0)
    return 2.0 * x_scaled**3 - 3.0 * x_scaled**2 + 1.0


# Returns a 1D array of length 'width'. The first outer_border values are 0, and the next
# inner_border-outer_border values are smoothly intermediate between 0 and 1. Then there
# is a series of 1s. The output is symmetric about the midpoint.
# See LW's notes p. GL105.
def one_axis_weight_function(width, outer_border, inner_border):
    import numpy as np
    assert (outer_border >= 0), "outer_border must be non-negative in one_axis_weight_function"
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
    oawf = one_axis_weight_function(16, 2, 5)
    for x in oawf:
        print(x)



# Returns a 1D array - this is the flattened version of the 2D array that is the
# geometric average of two one_axis_weight_functions, one on each of the two axes.
# Input 'shape' should be the desired shape of the before-being-flattened 2D array.
def two_axis_weight_function(shape, outer_border, inner_border):
    import numpy as np
    x = one_axis_weight_function(shape[0], outer_border, inner_border)
    y = one_axis_weight_function(shape[1], outer_border, inner_border)
    return np.sqrt(x[np.newaxis].T * y).reshape(-1)

def two_axis_weight_function_test_harness():
    import matplotlib.pyplot as plt
    import numpy as np
    d1 = 100
    d2 = 35
    outer_border = 7
    inner_border = 14
    a = two_axis_weight_function((d1, d2), outer_border, inner_border)
    a = np.reshape(a, (d1, d2))
    plt.imshow(a)
    plt.show()


    
# Returns the original shape, then the flattened ra, dec and values (in that order).
def get_glimpse_output_data(filename):
    import astropy.io.fits as pyfits
    pyfile = pyfits.open(filename)
    ra = pyfile[1].data
    dec = pyfile[2].data
    val = pyfile[0].data
    pyfile.close()
    return (val.shape, ra.reshape(-1), dec.reshape(-1), val.reshape(-1))

  
    

def merge_glimpse_output_files(input_file_spec, outer_border, inner_border, output_file_root, output_nside, cutouts_nside, job_control):
    import math
    import astropy.io.fits as pyfits
    import healpy as hp
    import numpy as np
    import glob
    import sys
    import matplotlib.pyplot as plt
    import utility
    

    RA = 0
    DEC = 1
    
    # Currently hard-coded to output in RING format; can change this if desired.
    output_nest = False
    
    output_debugging_info = False
    g_index_special = 6318085 # Dump debugging info for this pixel
    if output_debugging_info:
        print("Debugging point: RA={}; DEC={}; index={}".format(hp.pix2ang(output_nside, g_index_special, output_nest, lonlat=True)[0], hp.pix2ang(output_nside, g_index_special, output_nest, lonlat=True)[1], g_index_special))
    
    assert (outer_border <= inner_border), "inner_border must not be larger than outer_border"

    num_healpixels = hp.nside2npix(output_nside)

    # "g_" means "global" i.e. for the whole sphere. Indices into these arrays must repect the 'output_nest' parameter.
    g_weighted_values = np.zeros(num_healpixels)
    g_weights = np.zeros(num_healpixels)
    g_centres = hp.pix2ang(output_nside, np.arange(num_healpixels), output_nest, True)
    
    g_centres_ra_rad = np.radians(g_centres[RA])
    g_centres_dec_rad = np.radians(g_centres[DEC])
    g_sin_centres_ra = np.sin(g_centres_ra_rad)
    g_cos_centres_ra = np.cos(g_centres_ra_rad)
    g_sin_centres_dec = np.sin(g_centres_dec_rad)
    g_cos_centres_dec = np.cos(g_centres_dec_rad)
    
    

    list_of_input_files = glob.glob(input_file_spec)
    list_of_input_files.sort()
    list_of_input_files = list_of_input_files[utility.slice_from_string(job_control)]
    num_input_files = len(list_of_input_files)

    assert (num_input_files > 0), "No input files found"
    
    is_first_file = True

    for (filename, file_index) in zip(list_of_input_files, range(num_input_files)):

        if not output_debugging_info:
            print("Processing {} of {}: {}...".format(file_index+1, num_input_files, filename))
            
        # "p_" means "for this one picture"
        (p_shape, p_ra, p_dec, p_val) = get_glimpse_output_data(filename)
        
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
        # File name will be of the form (various}.XXXX.glimpse.out.fits
        #where XXXX is a healpix id (nside = cutouts_nside)
        assert (filename[-17:] == ".glimpse.out.fits", "Filename {} was not of the expected format (should end with .glimpse.out.fits)".format(filename)
        
        index_of_coarse_healpix_id_in_file_name = -4
        cutout_healpix_id = int(filename.split('.')[index_of_coarse_healpix_id_in_file_name])
        (p_centre_ra, p_centre_dec) = hp.pix2ang(cutouts_nside, cutout_healpix_id, nest=False, lonlat=True)
        
        (p_ra, p_dec) = utility.from_standard_position_fast(p_x, p_y, p_z, p_centre_ra, p_centre_dec)

        if any(np.isnan(p_val)):
            print("\t\tRejected as it contains NaNs")
        else:
            # The radius of the glimpse lattice (on the sphere).
            radius_rad = utility.angular_separation(p_centre_ra, p_centre_dec, p_ra[0], p_dec[0])

            # The angular separation between each healpixel and the centre of the current picture.
            #g_ang_sep_rad = utility.angular_separation(g_centres[RA], g_centres[DEC], p_centre_ra, p_centre_dec)
            g_ang_sep_rad = utility.angular_separation_fast(g_sin_centres_ra, g_cos_centres_ra, g_sin_centres_dec, g_cos_centres_dec, p_centre_ra, p_centre_dec)
            
            
            # Only include those healpixels that are sensibly close to the centre of the lattice.
            # This is a speedup, and also prevents problems with the 'project sphere to tangent plane'
            # routine (which would fail e.g. if the point on the sphere were on the wrong hemisphere).
            g_include = np.where(g_ang_sep_rad <= radius_rad)
            
            g_centres_ra_filtered = g_centres[RA][g_include]
            g_centres_dec_filtered = g_centres[DEC][g_include]
            g_indices_filtered = np.arange(num_healpixels)[g_include]
            
            # Map the (relevant) healpixel centres to the tangent plane for this picture.
            (g_x_filtered, g_y_filtered) = utility.sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, g_centres_ra_filtered, g_centres_dec_filtered)
            
            # We know that the glimpse lattice points will map to a perfectly regular lattice on the tangent plane,
            # so no need to actually do this computation. However, we do need to know the spacing of the glimpse lattice
            # points in the tangent plane. So do a simple experiment to discover this.
            glimpse_lattice_spacing = utility.sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra[1], p_dec[1])[0] - utility.sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra[0], p_dec[0])[0]
            
            if False:
                (p_lattice_x, p_lattice_y) = utility.sphere_to_tangent_plane_mapping(p_centre_ra, p_centre_dec, p_ra, p_dec)
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
            
            g_index_into_glimpse_array = (utility.bound_between(p_shape[0]//2 + np.floor(g_x_filtered_rotated / glimpse_lattice_spacing), 0, p_shape[0]-1) + \
                p_shape[0] * utility.bound_between(p_shape[1]//2 + np.floor(g_y_filtered_rotated / glimpse_lattice_spacing) - 1, 0, p_shape[1]-1)).astype(int)
                
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
   
    for (filename_part, array) in zip(["_values", "_weights"], [g_average_values, g_weights]):
        dat_filename = output_file_root + filename_part + ".dat"
        # Save data to healpix format file
        utility.write_healpix_array_to_file(array, output_nest, dat_filename)



def merge_glimpse_output_files_caller(ini_file_name, job_control):
    import configparser
    
    config = configparser.ConfigParser()
    config.read(ini_file_name)
    
    outer_border = int(config["merge"].get("outer_border"))
    inner_border = int(config["merge"].get("inner_border"))
    output_nside = int(config["merge"].get("output_nside"))
    input_file_spec = config["project"].get("directory") + config["project"].get("project_name") + ".{}.glimpse.out.fits"
    output_file_root = config["project"].get("directory") + config["project"].get("project_name")
    cutouts_nside = int(config["create_cutouts"].get("nside"))

    merge_glimpse_output_files(input_file_spec, outer_border, inner_border, output_file_root, output_nside, cutouts_nside, job_control)

# Start of checker functions for argparse
# See https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse    
def get_input_bool(arg):
    ua = str(arg).upper()
    if 'TRUE'.startswith(ua):
       return True
    elif 'FALSE'.startswith(ua):
       return False
    else:
        msg = "Input should be True or False but {0} was given".format(arg)
        raise argparse.ArgumentTypeError(msg)


def get_input_power_of_two(arg):
    import math
    n = int(arg)
    if (n <= 0) or ((2 ** int(math.log(n, 2))) != n):
        msg = "Input should be a power of two but {0} was given".format(arg)
        raise argparse.ArgumentTypeError(msg)
    return n


if __name__ == '__main__':

    import sys
    import argparse

    try:

        parser = argparse.ArgumentParser(description = "Merges output files from Glimpse into one healpix file. Set 'input_file_spec' to the file specification for the Glimpse files. These files are assumed to be in FITS format, with two-dimensional blocks for data, ra and dec (in that order). The data block therefore looks like a picture and we assume that pixels near the edge of this picture may be bad; downweight them to zero (using 'outer_border') or to less than one (using 'inner_border'). Use 'output_file_root' to determine where the output files will be written; two files will be created here, with name suffixes '_values.dat' and '_weights.dat'. The former contains the weighted average output values, while latter contains the weights that were used. Use 'nside' and 'nest' to determine the parameters of the output files. Use 'show_output_plots' to have the routine show the output maps that were created. Use 'save_output_plots' to have the routine save the output maps that were created in png format with name suffixes '_values.png' and '_weights.png'.")
        parser.add_argument('-i', '--input_file_spec', type = str, required = True, help = "File specification (may include wildcards) for Glimpse maps to be used as input.")
        parser.add_argument('-b', '--outer_border', type = int, required = True, help = "Outer border (in pixels). Pixels this close (or closer) to the edge will have zero weight.")
        parser.add_argument('-c', '--inner_border', type = int, required = True, help = "Inner border (in pixels). Pixels this close (or closer) to the edge will have weight less than one.")
        parser.add_argument('-t', '--output_file_root', type = str, required = True, help = "File name prefix for output files.")
        parser.add_argument('-n', '--nside', type = get_input_power_of_two, required = True, help = "Nside (power of two) to use in the output healpix files.")
        parser.add_argument('-e', '--nest', type = get_input_bool, required = False, default = False, help = "(Optional) style of the output healpix files - True for nested and False (this is the default) for ring.")
        parser.add_argument('-p', '--show_output_plots', type = get_input_bool, required = False, default = False, help = "(Optional) True to show plots of output data; False (this is the default) to not show plots.")
        parser.add_argument('-s', '--save_output_plots', type = get_input_bool, required = False, default = False, help = "(Optional) True to save plots of output data in png format; False (this is the default) to not save plots.")

        args = parser.parse_args()
        merge_glimpse_output_files(args.input_file_spec, args.outer_border, args.inner_border, args.output_file_root, args.nside, args.nest, args.show_output_plots, args.save_output_plots)

    except Exception as err:
        print('Error: {0}'.format(str(err)))
        sys.exit(1)
