#!/usr/bin/env python
# -*- coding: utf-8 -*-



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

def get_input_catalogue_format(arg):
    if not arg in ("csv", "fits"):
        msg = "Input format type should be csv or fits but {0} was given".format(arg)
        raise argparse.ArgumentTypeError(msg)
    return arg
    
def get_input_power_of_two(arg):
    import math
    n = int(arg)
    if (n <= 0) or ((2 ** int(math.log(n, 2))) != n):
        msg = "Input should be a power of two but {0} was given".format(arg)
        raise argparse.ArgumentTypeError(msg)
    return n
    
# End of checker functions for argparse

### STILL TO BE DONE: description text, and help for four arguments.   
    
if __name__ == '__main__':


    import argparse
    import sys
    import utility
    
    try:

        
        parser = argparse.ArgumentParser(description = "Given a weak-lensing catalogue, creates a set of 'cutout' catalogues given one input weak-lensing catalogue.")

        parser.add_argument('-i', '--input_catalogue', type = str, required = True, help = "Input catalogue file name. Allowed formats are fits and csv. At a minimum the file should have columns containing RA (decimal degrees [0, 360]) and Dec (decimal degrees [-90, 90]).")
        parser.add_argument('-c', '--cat_format', type = get_input_catalogue_format, required = True, help = "Input catalogue file format (fits or csv).") 
        parser.add_argument('-r', '--ra_name', type = str, required = True, help = "Field name for RA.")
        parser.add_argument('-d', '--dec_name', type = str, required = True, help = "Field name for Dec.")
        
        parser.add_argument('-s', '--shear_names', type = str, required = True, help = "##########TO DO########")
        
        parser.add_argument('-o', '--other_field_names', type = str, required = False, default = "", help = "(Optional) list of field names (in addition to RA and Dec) to be included in output file; comma-delimited.")
        parser.add_argument('-n', '--nside', type = get_input_power_of_two, required = True, help = "Coarse nside (power of two). This determines the number of augmented pixels.")
        parser.add_argument('-u', '--cutout_side_in_degrees', type = float, required = True, help = "##########TO DO########")
        
        parser.add_argument('-p', '--ids_to_process_begin', type = int, required = False, default = 0, help = "##########TO DO########")
        parser.add_argument('-q', '--ids_to_process_end', type = int, required = False, default = -1, help = "##########TO DO########")
        
        parser.add_argument('-y', '--output_directory', type = str, required = True, help = "Directory for output files. Will be created if it does not already exist.")
        parser.add_argument('-t', '--output_file_root', type = str, required = True, help = "File name prefix for output files.")
        
        args = parser.parse_args()
        
        utility.create_cutouts(args.input_catalogue, args.cat_format, args.ra_name, args.dec_name, args.shear_names, args.other_field_names, args.nside, args.cutout_side_in_degrees,
            args.ids_to_process_begin, args.ids_to_process_end, args.output_directory, args.output_file_root)
        
    except Exception as err:
        print('Error: {0}'.format(str(err)))
        sys.exit(1)
        
