#!/usr/bin/env python

RA = 0
DEC = 1
KAPPA = 2

SIMULATION = 0
GLIMPSE = 1


def run_view(file_root):
    import numpy as np
    import astropy.io.fits as pyfits
    import matplotlib.pyplot as plt
    #from astropy.visualization import astropy_mpl_style
    #plt.style.use(astropy_mpl_style)
    
    
    # TODO: get these from the ini file
    if file_root == "their_example":
        ra_name = "ra_gal"
        dec_name = "dec_gal"
        e1_name = "e1_gal"
        e2_name = "e2_gal"
        factor = 0.5
    else:
        ra_name = "ra"
        dec_name = "dec"
        e1_name = "e1"
        e2_name = "e2"
        factor = 0.005
    

    

    fig = plt.figure()
    
    # Input catalogue
    path_cat = './' + file_root + '.glimpse.cat.fits'
    x = pyfits.open(path_cat)
    ra = x[1].data.field(ra_name)
    dec = x[1].data.field(dec_name)
    e1 = x[1].data.field(e1_name)
    e2 = x[1].data.field(e2_name)
    x.close()
    
    w = e1 + 1j*e2
    w_prime = np.abs(w) * np.exp(1j * np.angle(w) / 2.0)
    e1_prime = np.real(w_prime)
    e2_prime = np.imag(w_prime)
    
    stretch = 1.0/np.cos(np.radians(np.mean(dec)))
    ax1 = fig.add_subplot(1,2,1, adjustable='box', aspect=stretch)
    
    if not True:
        ax1.scatter(ra, dec)
    else:
        ax1.plot([ra-e1_prime*factor, ra+e1_prime*factor], [dec-e2_prime*factor*stretch, dec+e2_prime*factor*stretch], c='black')

    
    # Output glimpse results
    # From https://docs.astropy.org/en/stable/generated/examples/io/plot_fits-image.html
    path_out = './' + file_root + '.glimpse.out.fits'
    image_data = pyfits.getdata(path_out, ext=0)
    ax2 = fig.add_subplot(1,2,2)
    ax2.imshow(image_data, cmap='gray', origin='lower') # https://stackoverflow.com/questions/14320159/matplotlib-imshow-data-rotated

    plt.show()
    
    



    
    
def hist2d_test_harness():
        
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.stats as st
    
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    bins = 2
    x = np.array([1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 2.0])
    y = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0])
    w = np.array([5.0, 3.0, 1.0, 2.0, 10.0, 4.0, -3.0])
    
    g = st.binned_statistic_2d(x, y, w, 'mean', bins=2)
    h = ax1.imshow(g[0], origin='lower') # https://stackoverflow.com/questions/14320159/matplotlib-imshow-data-rotated
    plt.colorbar(h, ax=ax1)
    
    
    
    if False:
        h_w = plt.hist2d(x, y, bins, weights=w)
        h = plt.hist2d(x, y, bins)
        d = {}
        d['x'] = x
        d['y'] = y
        d['weights'] = h_w[0]/h[0]
        h_r = ax1.hist2d('x', 'y', bins, weights='weights', data=d)
        plt.colorbar(h_r[3], ax=ax1)
        print(h_r)
    plt.show()




if __name__ == '__main__':
    import sys
    if len(sys.argv) <= 1:
        raise RuntimeError("Usage: view_glimpse file_root")
    if True:
        run_view(sys.argv[1])
    else:
        hist2d_test_harness()