#!/usr/bin/env python


######### START OF SECTION COPIED FROM https://raw.githubusercontent.com/NiallJeffrey/DeepMass/master/deepmass/lens_data.py #########
def ks_fourier_matrix(size):
    """

    :param size: size of square image
    :return: diagonal of forward (kappa to shear) Fourier operator
    """
    
    from scipy import fftpack
    import numpy as np

    k_modes = fftpack.fftfreq(size)
    k1_grid = np.dstack(np.meshgrid(k_modes, k_modes))[:, :, 0]
    k1_vector = np.reshape(k1_grid, -1)
    k2_vector = np.reshape(k1_grid.T, -1)

    k_squared = k1_vector * k1_vector + k2_vector * k2_vector
    k_squared2 = np.where(k_squared == 0, 1e-18, k_squared)

    A_ft_diagonal = (k1_vector ** 2 - k2_vector ** 2 + 1j * 2.0 * k1_vector * k2_vector) / k_squared2
    return np.where(k_squared != 0.0, A_ft_diagonal, 1.0)



def ks(shear_map, fourier_forward_matrix=None):
    """
    Kaiser squires 1993
    :param shear_map: complex square shear map (e1 + i e2)
    :param fourier_forward_matrix: if matrix precalculated use it for speed
    :return: kappa map
    """
    
    import numpy as np

    if fourier_forward_matrix is None:
        fourier_forward_matrix = ks_fourier_matrix(shear_map.shape[0])

    fourier_shear_vector = np.fft.fft2(shear_map).flatten()

    return np.fft.ifft2(np.reshape(fourier_shear_vector/fourier_forward_matrix,shear_map.shape))


######### END OF SECTION COPIED FROM https://raw.githubusercontent.com/NiallJeffrey/DeepMass/master/deepmass/lens_data.py #########

def ks_test_harness():
    import numpy as np
    np.set_printoptions(precision=2)
    n_rows_half = 2
    shear_map = np.zeros((2*n_rows_half+1, 2*n_rows_half+1),dtype=np.complex_)
    shear_map[n_rows_half, n_rows_half] = 10.0 + 0.0j
    kappa = ks(shear_map)
    print(shear_map)
    print(kappa)



    
def buzzard_ks_analysis(bins):

    import astropy.io.fits as pyfits
    import scipy.stats as st
    import scipy.ndimage as ndi
    import numpy as np
    
    path_fits = '/share/splinter/ucapwhi/glimpse_project/scripts/sample_buzzard_rectangle.glimpse.cat.fits'
    print("Open input file " + path_fits + "...")
    x = pyfits.open(path_fits)
    ra = x[1].data.field("ra")
    dec = x[1].data.field("dec")
    true_z = x[1].data.field("z")
    e1 = x[1].data.field("e1")
    e2 = x[1].data.field("e2")
    g1 = x[1].data.field("g1")
    g2 = x[1].data.field("g2")
    kappa = x[1].data.field("kappa")
    (binned_g1, x_edge_g1, y_edge_g1, binnumber_g1) = st.binned_statistic_2d(ra, dec, g1, 'mean', bins)
    (binned_g2, x_edge_g2, y_edge_g2, binnumber_g2) = st.binned_statistic_2d(ra, dec, g2, 'mean', bins)
    (binned_kappa, x_edge_kappa, y_edge_kappa, binnumber_kappa) = st.binned_statistic_2d(ra, dec, kappa, 'mean', bins)
    
    binned_g1 = np.nan_to_num(binned_g1)
    binned_g2 = np.nan_to_num(binned_g2)
    binned_kappa = np.nan_to_num(binned_kappa)
    
    shear_map = -binned_g1 + binned_g2 * (0.0 + 1.0j) # Note the sign conventions used here.
    ks_kappa = ks(shear_map)
    real_ks_kappa = np.real(ks_kappa)
    
    if False:
        smoothing_scale_in_pixels = 1
        real_ks_kappa = ndi.gaussian_filter(real_ks_kappa, smoothing_scale_in_pixels)
        binned_kappa = ndi.gaussian_filter(binned_kappa, smoothing_scale_in_pixels)
        
    return (binned_kappa, real_ks_kappa)
    
# Input should be a two element list: truth, then reconstruction
def display_buzzard_ks_analysis(results):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    this_min = min(np.amin(results[0]), np.amin(results[1]))
    this_max = min(np.amax(results[0]), np.amax(results[1]))
    
    plot_title = ["Buzzard true kappa", "Kaiser-Squires kappa using noiseless shear"]
    
    
    fig = plt.figure()
    for i in range(2):
    
        ax = fig.add_subplot(1,2,(i+1))
        h = ax.imshow(results[i], cmap='viridis', origin='lower', vmin=this_min, vmax=this_max)
        plt.title(plot_title[i])
        ax.set_adjustable('box-forced')
        #plt.colorbar(h, ax=ax)
    
    
    plt.show()



    

if __name__ in '__main__':
    bins = 400
    display_buzzard_ks_analysis(buzzard_ks_analysis(bins))

    

