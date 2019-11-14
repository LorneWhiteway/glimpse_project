#!/usr/bin/env python



def run_extract():
    import h5py
    import numpy as np
    import astropy.io.fits as pyfits
    import sys
    
    sys.path.insert(0, "/share/splinter/ucapwhi/cosmic_web/hong_dey")
    import cosmic_web_utilities as cwu

    # Open and read in from fits file
    
    path_fits = '/share/testde/ucapnje/buzzard_desy3_marco/Buzzard_192.fits'
    print("Open input file " + path_fits + "...")
    x = pyfits.open(path_fits)
    
    ra = x[1].data.field("RA")
    dec = x[1].data.field("DEC")
    true_z = x[1].data.field("true_z")
    e1 = x[1].data.field("E1")
    e2 = x[1].data.field("E2")
    g1 = x[1].data.field("G1")
    g2 = x[1].data.field("G2")
    kappa = x[1].data.field("k_orig")
    
    # Filter to a subset
    
    
    #centre = (55.0, -30.0) # RA, DEC (decimal degrees)
    centre = (110.45454545, -57.3995236) # RA, DEC (decimal degrees)
    if True:
        # Filter method 1: circle on the sky
        radius_in_degrees = 8.0
        separation_in_degrees = np.degrees(cwu.angular_separation(ra, dec, centre[0], centre[1]))
        filter_condition = (0.0, radius_in_degrees, 0.0, 360.0, -90.0, 90.0)
        (f_separation_in_degrees, f_ra, f_dec, f_true_z, f_e1, f_e2, f_g1, f_g2, f_kappa) = cwu.triple_filter(separation_in_degrees, ra, dec, filter_condition, true_z, e1, e2, g1, g2, kappa)
    else:
        # Filter method 2: square on the sky
        side_in_degrees = 4.0
        distortion_factor = 1.0 / np.cos(np.radians(centre[1])) # See p. GL40
        filter_condition = (centre[0] - distortion_factor*side_in_degrees/2.0, centre[0] + distortion_factor*side_in_degrees/2.0, centre[1] - side_in_degrees/2.0, centre[1] + side_in_degrees/2.0, 0.0, 10000.0)
        (f_ra, f_dec, f_true_z, f_e1, f_e2, f_g1, f_g2, f_kappa) = cwu.triple_filter(ra, dec, true_z, filter_condition, e1, e2, g1, g2, kappa)
    
    print(len(ra), len(f_ra))
    
    # Save the filtered results.
    output_catalogue_name = '/share/splinter/ucapwhi/glimpse_project/scripts/sample_buzzard_disk_2821.glimpse.cat.fits'
    print("Open output file " + output_catalogue_name + "...")

    y = pyfits.BinTableHDU.from_columns([
        pyfits.Column(name='ra', format='D', array=f_ra),
        pyfits.Column(name='dec', format='D', array=f_dec),
        pyfits.Column(name='z', format='D', array=f_true_z),
        pyfits.Column(name='e1', format='D', array=f_e1),
        pyfits.Column(name='e2', format='D', array=f_e2),
        pyfits.Column(name='g1', format='D', array=f_g1),
        pyfits.Column(name='g2', format='D', array=f_g2),
        pyfits.Column(name='kappa', format='D', array=f_kappa)])

    y.writeto(output_catalogue_name, overwrite=True)
    
    x.close()



if __name__ == '__main__':
    run_extract()