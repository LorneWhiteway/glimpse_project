#!/usr/bin/env python

# This scriot probably won't be useful as the data file /share/testde/ucapnje/y3_sims/Y3_mastercat_w_y3rmg_b3_v1.9.2.h5
# is incomplete (in that it contains internal pointers to other files that don't exist).




def run_extract():
    import h5py
    import numpy as np
    import astropy.io.fits as pyfits

    # Open and read in from Buzzard file

    path_Buzzard = '/share/testde/ucapnje/y3_sims/Y3_mastercat_w_y3rmg_b3_v1.9.2.h5'
    print("Open input file " + path_Buzzard + "...")

    mute = h5py.File(path_Buzzard, "r")
    buz_ra = np.array(mute['catalog']['gold']['ra'])
    buz_dec = np.array(mute['catalog']['gold']['dec'])
    buz_e1 = np.array(mute['catalog']['metacal'][u'unsheared']['e1'])
    buz_e2 = np.array(mute['catalog']['metacal'][u'unsheared']['e2'])
    buz_g1 = np.array(mute['catalog']['metacal'][u'unsheared']['g1'])
    buz_g2 = np.array(mute['catalog']['metacal'][u'unsheared']['g2'])
    buz_kappa = np.array(mute['catalog']['metacal'][u'unsheared']['kappa'])
    buz_id= np.array(mute['catalog']['metacal'][u'unsheared']['coadd_object_id'])
    buz_id2=np.array(mute['catalog']['bpz'][u'unsheared'][u'coadd_object_id'])
    buz_zmc=np.array(mute['catalog']['bpz'][u'unsheared']['zmc_sof'])
    buz_zmean=np.array(mute['catalog']['bpz'][u'unsheared']['zmean_sof'])
    buz_z=np.array(mute['catalog']['bpz'][u'unsheared']['z'])

    print("Num records = " + str(len(buz_ra)))

    # Save to FITS file
    output_catalogue_name = '/share/splinter/ucapwhi/glimpse_project/scripts/out.fits'
    print("Open output file " + output_catalogue_name + "...")

    tbhdu = pyfits.BinTableHDU.from_columns([
        pyfits.Column(name='ra', format='D', array=buz_ra),
        pyfits.Column(name='dec', format='D', array=buz_dec),
        pyfits.Column(name='z', format='D', array=buz_z)])

    tbhdu.writeto(output_catalogue_name, clobber=True)



if __name__ == '__main__':
    run_extract()