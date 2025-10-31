#Script para comparar los resultados de z_phot-z_espec frente a los bines de photo_z+
#corregido z_spec y añadido barras de error con bootstrap
##nuevos resultados con la definición de zphot y zspec bien
#ver los resultados del match wise y vhs y comparar dnf con anf
#Ejecución: 
#python plots_comparar_vhs_bootstrap_zspec.py dnf_vhs_bandas_optico_corregido.fits dnf_vhs_corregido.fits
#python incertidumbres_wise.py anf_wise_bandas_optico.fits anf_wise_W1_W2.fits anf_wise_cony.fits

#Mirar los filtros como los pongo y en funcion del resultado los _1 etc

import numpy as np
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
import pylab
import weighted_kde as kde
from scipy.stats import bootstrap
from matplotlib.colors import LogNorm

# Read files
file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]

plt.rcParams.update({
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.titlesize': 16
})


def process_data(file, filters):
    mag_name = []
    magerr_name = []
    nfilters = len(filters)

    for f in filters:
        mag_name.append(f'BDF_MAG_{f}_CORRECTED')
        magerr_name.append(f'BDF_MAG_ERR_{f}')

    t = Table.read(file)
    t['Z_MEAN'] = t['DNF_Z']
    t['Z_MC'] = t['DNF_ZN']
    t['Z_SIGMA'] = t['DNF_ZSIGMA']
    
    # Cleaning
    sel = t['Z_MEAN'] != -99.0
    t = t[sel]
    sel = t['Z_MEAN'] > 0.0
    t = t[sel]
    sel = ~np.isnan(t['Z_MEAN'])
    t = t[sel]
    Nvalid = len(t)
    print(Nvalid)
    M = np.zeros((Nvalid, nfilters), dtype='double')
    Merr = np.zeros((Nvalid, nfilters), dtype='double')

    for k in range(nfilters):
        M[:, k] = np.array([t[mag_name[k]]])
        Merr[:, k] = np.array([t[magerr_name[k]]])

    sel = np.ones(Nvalid, dtype='bool')
    for k in range(nfilters):
        selaux = np.logical_and(8 < M[:, k], M[:, k] < 26.0)
        sel = np.logical_and(sel, selaux)

    t = t[sel]
    Nvalid = len(t)
    print(Nvalid)
    # Calculate bias and sigma, lo comento para el TFM bueno
    #bias = np.mean(t['Z'] - t['Z_MEAN'])
    #sigma = np.std(t['Z'] - t['Z_MEAN'])
    #zerrabs = np.abs(t['Z'] - t['Z_MEAN'])
    #zerrSort = np.sort(zerrabs)
    #zerrsigma68 = zerrSort[int(Nvalid * 68 / 100)]
    bias = np.mean(np.abs(t['Z_MEAN'] - t['Z']))
    sigma = np.std(t['Z_MEAN'] - t['Z'])
    zerrabs = np.abs(t['Z_MEAN'] - t['Z'])
    zerrSort = np.sort(zerrabs)
    zerrsigma68 = zerrSort[int(Nvalid * 68 / 100)]

    return t, bias, sigma, zerrsigma68

filters1 = ['G', 'R', 'I', 'Z']  # Example filters for file1
filters2 = ['G', 'R', 'I', 'Z', 'W1', 'W2']  # Example filters for file2
filters3 = ['G', 'R', 'I', 'Z', 'Y','W1', 'W2']  # Example filters for file3




t1, bias1, sigma1, zerrsigma68_1 = process_data(file1, filters1)
t2, bias2, sigma2, zerrsigma68_2 = process_data(file2, filters2)
t3, bias3, sigma3, zerrsigma68_3 = process_data(file3, filters3)

print('longitudes:' , len(t1), len(t2), len(t3))
print('bias:' , bias1, bias2, bias3)
print('sigma:' , sigma1, sigma2, sigma3)
print('sigma68:' , zerrsigma68_1, zerrsigma68_2, zerrsigma68_3)

def calculate_means_and_errorsMAD(t, z_min=0.2, z_max=0.8):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)  
    delta_z = np.ma.filled((t['Z_MEAN'] - t['Z']), np.nan)  

    mask = (photo_z >= z_min) & (photo_z < z_max)
    dz_bin = delta_z[mask]

    if len(dz_bin) > 0:
        mean_val = np.mean(np.abs(dz_bin))
        res = bootstrap((np.abs(dz_bin),), np.mean, n_resamples=100, confidence_level=0.95)
        error_val = res.standard_error
    else:
        mean_val = np.nan
        error_val = np.nan

    return mean_val, error_val



meansMAD1, errorsMAD1 = calculate_means_and_errorsMAD(t1,0.2,0.8)
meansMAD2, errorsMAD2 = calculate_means_and_errorsMAD(t2,0.2,0.8)
meansMAD3, errorsMAD3 = calculate_means_and_errorsMAD(t3,0.2,0.8)


print(f'Catalog 1 (filters {filters1}): MAD = {meansMAD1:.4f} ± {errorsMAD1:.4f}')
print(f'Catalog 2 (filters {filters2}): MAD = {meansMAD2:.4f} ± {errorsMAD2:.4f}')
print(f'Catalog 3 (filters {filters3}): MAD = {meansMAD3:.4f} ± {errorsMAD3:.4f}')


meansMAD1, errorsMAD1 = calculate_means_and_errorsMAD(t1,0.8,1.4)
meansMAD2, errorsMAD2 = calculate_means_and_errorsMAD(t2,0.8,1.4)
meansMAD3, errorsMAD3 = calculate_means_and_errorsMAD(t3,0.8,1.4)

print(f'Catalog 1 (filters {filters1}): MAD = {meansMAD1:.4f} ± {errorsMAD1:.4f}')
print(f'Catalog 2 (filters {filters2}): MAD = {meansMAD2:.4f} ± {errorsMAD2:.4f}')
print(f'Catalog 3 (filters {filters3}): MAD = {meansMAD3:.4f} ± {errorsMAD3:.4f}')

meansMAD1, errorsMAD1 = calculate_means_and_errorsMAD(t1,0.2,1.4)
meansMAD2, errorsMAD2 = calculate_means_and_errorsMAD(t2,0.2,1.4)
meansMAD3, errorsMAD3 = calculate_means_and_errorsMAD(t3,0.2,1.4)

print(f'Catalog 1 (filters {filters1}): MAD = {meansMAD1:.4f} ± {errorsMAD1:.4f}')
print(f'Catalog 2 (filters {filters2}): MAD = {meansMAD2:.4f} ± {errorsMAD2:.4f}')
print(f'Catalog 3 (filters {filters3}): MAD = {meansMAD3:.4f} ± {errorsMAD3:.4f}')
def calculate_sigma_global(t, zmin=0.2, zmax=0.8):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma_val = []
    sigma_err = []
    mask = (photo_z >= zmin) & (photo_z <= zmax)
    dz_range = delta_z[mask]

    if len(dz_range) > 0:
        sigma_val = np.std(dz_range)
        res = bootstrap((dz_range,), np.std, n_resamples=100, confidence_level=0.95)
        sigma_err = res.standard_error
    else:
        sigma_val = np.nan
        sigma_err = np.nan

    return sigma_val, sigma_err


# Ejemplo de uso
sigma1, sigma_err1 = calculate_sigma_global(t1)
sigma2, sigma_err2 = calculate_sigma_global(t2)
sigma3, sigma_err3 = calculate_sigma_global(t3)

print(f"Catalog 1 sigma (0.2-0.8): {sigma1:.4f} ± {sigma_err1:.4f}")
print(f"Catalog 2 sigma (0.2-0.8): {sigma2:.4f} ± {sigma_err2:.4f}")
print(f"Catalog 3 sigma (0.2-0.8): {sigma3:.4f} ± {sigma_err3:.4f}")
sigma1, sigma_err1 = calculate_sigma_global(t1,0.8,1.4)
sigma2, sigma_err2 = calculate_sigma_global(t2,0.8,1.4)
sigma3, sigma_err3 = calculate_sigma_global(t3,0.8,1.4)

print(f"Catalog 1 sigma (0.2-0.8): {sigma1:.4f} ± {sigma_err1:.4f}")
print(f"Catalog 2 sigma (0.2-0.8): {sigma2:.4f} ± {sigma_err2:.4f}")
print(f"Catalog 3 sigma (0.2-0.8): {sigma3:.4f} ± {sigma_err3:.4f}")
sigma1, sigma_err1 = calculate_sigma_global(t1,0.2,1.4)
sigma2, sigma_err2 = calculate_sigma_global(t2,0.2,1.4)
sigma3, sigma_err3 = calculate_sigma_global(t3,0.2,1.4)

print(f"Catalog 1 sigma (0.2-0.8): {sigma1:.4f} ± {sigma_err1:.4f}")
print(f"Catalog 2 sigma (0.2-0.8): {sigma2:.4f} ± {sigma_err2:.4f}")
print(f"Catalog 3 sigma (0.2-0.8): {sigma3:.4f} ± {sigma_err3:.4f}")

def calculate_sigma68_global_with_error(t, zmin=0.2, zmax=0.8):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)

    mask = (photo_z >= zmin) & (photo_z <= zmax)
    dz_range = delta_z[mask]

    if len(dz_range) > 10:
        lower = np.percentile(dz_range, 16)
        upper = np.percentile(dz_range, 84)
        sigma68_val = (upper - lower) / 2

        # Bootstrap para el error
        res = bootstrap((dz_range,), lambda x: (np.percentile(x, 84) - np.percentile(x, 16)) / 2, 
                        n_resamples=100, confidence_level=0.95)
        sigma68_err = res.standard_error
    else:
        sigma68_val = np.nan
        sigma68_err = np.nan

    return sigma68_val, sigma68_err


# Uso para cada catálogo
sigma68_1, sigma68_err1 = calculate_sigma68_global_with_error(t1)
sigma68_2, sigma68_err2 = calculate_sigma68_global_with_error(t2)
sigma68_3, sigma68_err3 = calculate_sigma68_global_with_error(t3)

print(f"Catalog 1 sigma68 (0.2-0.8): {sigma68_1:.4f} ± {sigma68_err1:.4f}")
print(f"Catalog 2 sigma68 (0.2-0.8): {sigma68_2:.4f} ± {sigma68_err2:.4f}")
print(f"Catalog 3 sigma68 (0.2-0.8): {sigma68_3:.4f} ± {sigma68_err3:.4f}")
sigma68_1, sigma68_err1 = calculate_sigma68_global_with_error(t1,0.8,1.4)
sigma68_2, sigma68_err2 = calculate_sigma68_global_with_error(t2,0.8,1.4)
sigma68_3, sigma68_err3 = calculate_sigma68_global_with_error(t3,0.8,1.4)

print(f"Catalog 1 sigma68 (0.2-0.8): {sigma68_1:.4f} ± {sigma68_err1:.4f}")
print(f"Catalog 2 sigma68 (0.2-0.8): {sigma68_2:.4f} ± {sigma68_err2:.4f}")
print(f"Catalog 3 sigma68 (0.2-0.8): {sigma68_3:.4f} ± {sigma68_err3:.4f}")
sigma68_1, sigma68_err1 = calculate_sigma68_global_with_error(t1,0.2,1.4)
sigma68_2, sigma68_err2 = calculate_sigma68_global_with_error(t2,0.2,1.4)
sigma68_3, sigma68_err3 = calculate_sigma68_global_with_error(t3,0.2,1.4)

print(f"Catalog 1 sigma68 (0.2-0.8): {sigma68_1:.4f} ± {sigma68_err1:.4f}")
print(f"Catalog 2 sigma68 (0.2-0.8): {sigma68_2:.4f} ± {sigma68_err2:.4f}")
print(f"Catalog 3 sigma68 (0.2-0.8): {sigma68_3:.4f} ± {sigma68_err3:.4f}")
