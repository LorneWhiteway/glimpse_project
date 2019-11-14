#!/usr/bin/env python

def run():

    import numpy as np
    import astropy.io.fits as pyfits
    
    num_lensed_galaxies = 100
    
    if True:
        # One lens
        r_max = 0.16667 # degrees
        centre = (0.0, 0.0) #degrees
        
        
    else:
        # Two lenses
        r_max = 0.16667/2.0 # degrees
        indicator = np.random.randint(0, 2, num_lensed_galaxies) * 2 - 1 # Will be +1 or -1
        centre = (indicator*r_max, -indicator*r_max)
    
    
    k = 1.0 # scaling constant for shear
    shape_noise_factor = 0.05 # Set to zero for no shape noise

    alpha = np.random.uniform(0.0, 2.0 * np.pi, num_lensed_galaxies) # radians
    r = np.sqrt(np.random.uniform(0.8, 1.0, num_lensed_galaxies)) * r_max # degrees

    ra = centre[0] + r * np.cos(alpha) # degrees
    dec = centre[1] + r * np.sin(alpha) # degrees    

    e1 = -k * ((r/r_max)**(-3.0)) * np.cos(2.0 * alpha) + np.random.uniform(-shape_noise_factor*0.5, shape_noise_factor*0.5, num_lensed_galaxies)
    e2 = -k * ((r/r_max)**(-3.0)) * np.sin(2.0 * alpha) + np.random.uniform(-shape_noise_factor*0.5, shape_noise_factor*0.5, num_lensed_galaxies)
    
    z = np.ones(num_lensed_galaxies)

    # OK, write this out in FITS format

    # Save to FITS file
    output_catalogue_name = '/share/splinter/ucapwhi/glimpse_project/scripts/toy_data.glimpse.cat.fits'
    print("Writing to " + output_catalogue_name + "...")

    tbhdu = pyfits.BinTableHDU.from_columns([
        pyfits.Column(name='ra', format='D', array=ra),
        pyfits.Column(name='dec', format='D', array=dec),
        pyfits.Column(name='alpha', format='D', array=alpha),
        pyfits.Column(name='r', format='D', array=r),
        pyfits.Column(name='z', format='D', array=z),
        pyfits.Column(name='e1', format='D', array=e1),
        pyfits.Column(name='e2', format='D', array=e2)])

    tbhdu.writeto(output_catalogue_name, overwrite=True)




if __name__ == '__main__':
    run()