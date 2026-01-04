
#Script para comparar los resultados de z_phot-z_espec frente a los bines de photo_z+
#corregido z_spec y añadido barras de error con bootstrap
##nuevos resultados con la definición de zphot y zspec bien
#ver los resultados del match wise y vhs y comparar dnf con anf
#Ejecución: 
#python plots_comparar_vhs_bootstrap_zspec.py dnf_vhs_bandas_optico_corregido.fits dnf_vhs_corregido.fits
#python plots_comparar_vhs_bootstrap_zspec_JKS.py anf_vhs_bandas_optico_corregido.fits anf_vhs_corregido.fits anf_vhs_JKs.fits

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


filters1 = ['G', 'R', 'I', 'Z']  
filters2 = ['G', 'R', 'I', 'Z', 'J', 'H', 'KS'] 
filters3 = ['G', 'R', 'I', 'Z', 'J', 'KS'] 




t1, bias1, sigma1, zerrsigma68_1 = process_data(file1, filters1)
t2, bias2, sigma2, zerrsigma68_2 = process_data(file2, filters2)
t3, bias3, sigma3, zerrsigma68_3 = process_data(file3, filters3)

print('longitudes:' , len(t1), len(t2), len(t3))
print('bias:' , bias1, bias2, bias3)
print('sigma:' , sigma1, sigma2, sigma3)
print('sigma68:' , zerrsigma68_1, zerrsigma68_2, zerrsigma68_3)

def calc_stats_halves(t, bins):
    photo_z = t['Z_MEAN']
    delta_z = t['Z_MEAN'] - t['Z']

    total_bias = np.mean(np.abs(delta_z))
    total_sigma = np.std(delta_z)
    zerrabs = np.abs(delta_z)
    total_sigma68 = np.percentile(zerrabs, 68)  # 68% percentil

    mid = (bins[0] + bins[-1]) / 2

    mask_half1 = (photo_z >= bins[0]) & (photo_z < mid)
    mask_half2 = (photo_z >= mid) & (photo_z < bins[-1])

    bias_half1 = np.mean(np.abs(delta_z[mask_half1]))
    bias_half2 = np.mean(np.abs(delta_z[mask_half2]))

    sigma_half1 = np.std(delta_z[mask_half1])
    sigma_half2 = np.std(delta_z[mask_half2])

    sigma68_half1 = np.percentile(np.abs(delta_z[mask_half1]), 68)
    sigma68_half2 = np.percentile(np.abs(delta_z[mask_half2]), 68)

    print("Total:")
    print(f" Bias = {total_bias:.4f}, Sigma = {total_sigma:.4f}, Sigma68 = {total_sigma68:.4f}")
    print("Primera mitad:")
    print(f" Bias = {bias_half1:.4f}, Sigma = {sigma_half1:.4f}, Sigma68 = {sigma68_half1:.4f}")
    print("Segunda mitad:")
    print(f" Bias = {bias_half2:.4f}, Sigma = {sigma_half2:.4f}, Sigma68 = {sigma68_half2:.4f}\n")

    return (total_bias, total_sigma, total_sigma68), \
           (bias_half1, sigma_half1, sigma68_half1), \
           (bias_half2, sigma_half2, sigma68_half2)

bins = np.linspace(0.2, 1.4, 10)

print("Stats para t1:")
calc_stats_halves(t1, bins)
print("Stats para t2:")
calc_stats_halves(t2, bins)
print("Stats para t3:")
calc_stats_halves(t3, bins)

"""
import numpy as np
from scipy.stats import bootstrap

def bootstrap_error(data, func, n_resamples=100):
    if len(data) > 10:
        res = bootstrap((data,), func, n_resamples=n_resamples, confidence_level=0.95)
        return res.standard_error
    else:
        return np.nan

def calculate_bias(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    
    def bias_func(x):
        return np.mean(np.abs(x))

    def calc(data):
        data = data[~np.isnan(data)]
        if len(data) == 0:
            return np.nan, np.nan
        val = bias_func(data)
        err = bootstrap_error(data, bias_func)
        return val, err

    total = calc(delta_z)
    mid = (bins[0] + bins[-1]) / 2
    half1 = calc(delta_z[(photo_z >= bins[0]) & (photo_z < mid)])
    half2 = calc(delta_z[(photo_z >= mid) & (photo_z < bins[-1])])

    print("=== Bias ===")
    print(f"Total: {total[0]:.4f} ± {total[1]:.4f}")
    print(f"1ª mitad: {half1[0]:.4f} ± {half1[1]:.4f}")
    print(f"2ª mitad: {half2[0]:.4f} ± {half2[1]:.4f}\n")

def calculate_sigma(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    
    def sigma_func(x):
        return np.std(x)

    def calc(data):
        data = data[~np.isnan(data)]
        if len(data) == 0:
            return np.nan, np.nan
        val = sigma_func(data)
        err = bootstrap_error(data, sigma_func)
        return val, err

    total = calc(delta_z)
    mid = (bins[0] + bins[-1]) / 2
    half1 = calc(delta_z[(photo_z >= bins[0]) & (photo_z < mid)])
    half2 = calc(delta_z[(photo_z >= mid) & (photo_z < bins[-1])])

    print("=== Sigma ===")
    print(f"Total: {total[0]:.4f} ± {total[1]:.4f}")
    print(f"1ª mitad: {half1[0]:.4f} ± {half1[1]:.4f}")
    print(f"2ª mitad: {half2[0]:.4f} ± {half2[1]:.4f}\n")

def calculate_sigma68(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    
    def sigma68_func(x):
        return (np.percentile(x, 84) - np.percentile(x, 16)) / 2

    def calc(data):
        data = data[~np.isnan(data)]
        if len(data) == 0:
            return np.nan, np.nan
        val = sigma68_func(data)
        err = bootstrap_error(data, sigma68_func)
        return val, err

    total = calc(delta_z)
    mid = (bins[0] + bins[-1]) / 2
    half1 = calc(delta_z[(photo_z >= bins[0]) & (photo_z < mid)])
    half2 = calc(delta_z[(photo_z >= mid) & (photo_z < bins[-1])])

    print("=== Sigma 68 ===")
    print(f"Total: {total[0]:.4f} ± {total[1]:.4f}")
    print(f"1ª mitad: {half1[0]:.4f} ± {half1[1]:.4f}")
    print(f"2ª mitad: {half2[0]:.4f} ± {half2[1]:.4f}\n")

def calculate_outliers(t, bins, threshold=0.15):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    
    def outlier_func(x):
        return np.sum(np.abs(x) > threshold) / len(x)

    def calc(data):
        data = data[~np.isnan(data)]
        if len(data) == 0:
            return np.nan, np.nan
        val = outlier_func(data)
        err = bootstrap_error(data, outlier_func)
        return val, err

    total = calc(delta_z)
    mid = (bins[0] + bins[-1]) / 2
    half1 = calc(delta_z[(photo_z >= bins[0]) & (photo_z < mid)])
    half2 = calc(delta_z[(photo_z >= mid) & (photo_z < bins[-1])])

    print("=== Outliers (±{:.2f}) ===".format(threshold))
    print(f"Total: {total[0]:.4f} ± {total[1]:.4f}")
    print(f"1ª mitad: {half1[0]:.4f} ± {half1[1]:.4f}")
    print(f"2ª mitad: {half2[0]:.4f} ± {half2[1]:.4f}\n")

# Uso ejemplo para t1, t2, t3

bins = np.linspace(0.2, 1.4, 10)

print("Resultados para t1:")
calculate_bias(t1, bins)
calculate_sigma(t1, bins)
calculate_sigma68(t1, bins)
calculate_outliers(t1, bins)

print("Resultados para t2:")
calculate_bias(t2, bins)
calculate_sigma(t2, bins)
calculate_sigma68(t2, bins)
calculate_outliers(t2, bins)

print("Resultados para t3:")
calculate_bias(t3, bins)
calculate_sigma(t3, bins)
calculate_sigma68(t3, bins)
calculate_outliers(t3, bins)



bins = np.linspace(0.2, 1.4, 10)  # Bins completos
bin_centers = (bins[:-1] + bins[1:]) / 2


def calculate_means_and_errorsMAD(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)  
    delta_z = np.ma.filled((t['Z_MEAN'] - t['Z']), np.nan)  
    means = []
    errors = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 0:
            mean_val = np.mean(np.abs(dz_bin))
            res = bootstrap((np.abs(dz_bin),), np.mean, n_resamples=100, confidence_level=0.95)
            error_val = res.standard_error
        else:
            mean_val = np.nan
            error_val = np.nan

        means.append(mean_val)
        errors.append(error_val)

    return means, errors


meansMAD1, errorsMAD1 = calculate_means_and_errorsMAD(t1, bins)
meansMAD2, errorsMAD2 = calculate_means_and_errorsMAD(t2, bins)
meansMAD3, errorsMAD3 = calculate_means_and_errorsMAD(t3, bins)

bias1_total = np.nanmean(meansMAD1)
err1_total = np.nanmean(errorsMAD1)
bias2_total = np.nanmean(meansMAD2)
err2_total = np.nanmean(errorsMAD2)
bias3_total = np.nanmean(meansMAD3)
err3_total = np.nanmean(errorsMAD3)

mid = len(bins) // 2

# Primera mitad bins y segunda mitad bins
bins_h1 = bins[:mid+1]
bins_h2 = bins[mid:]

# t1 mitades
mask1_h1 = (t1['Z_MEAN'] >= bins_h1[0]) & (t1['Z_MEAN'] < bins_h1[-1])
mask1_h2 = (t1['Z_MEAN'] >= bins_h2[0]) & (t1['Z_MEAN'] < bins_h2[-1])
meansMAD1_h1, errorsMAD1_h1 = calculate_means_and_errorsMAD(t1[mask1_h1], bins_h1)
meansMAD1_h2, errorsMAD1_h2 = calculate_means_and_errorsMAD(t1[mask1_h2], bins_h2)
bias1_h1 = np.nanmean(meansMAD1_h1)
err1_h1 = np.nanmean(errorsMAD1_h1)
bias1_h2 = np.nanmean(meansMAD1_h2)
err1_h2 = np.nanmean(errorsMAD1_h2)

# t2 mitades
mask2_h1 = (t2['Z_MEAN'] >= bins_h1[0]) & (t2['Z_MEAN'] < bins_h1[-1])
mask2_h2 = (t2['Z_MEAN'] >= bins_h2[0]) & (t2['Z_MEAN'] < bins_h2[-1])
meansMAD2_h1, errorsMAD2_h1 = calculate_means_and_errorsMAD(t2[mask2_h1], bins_h1)
meansMAD2_h2, errorsMAD2_h2 = calculate_means_and_errorsMAD(t2[mask2_h2], bins_h2)
bias2_h1 = np.nanmean(meansMAD2_h1)
err2_h1 = np.nanmean(errorsMAD2_h1)
bias2_h2 = np.nanmean(meansMAD2_h2)
err2_h2 = np.nanmean(errorsMAD2_h2)

# t3 mitades
mask3_h1 = (t3['Z_MEAN'] >= bins_h1[0]) & (t3['Z_MEAN'] < bins_h1[-1])
mask3_h2 = (t3['Z_MEAN'] >= bins_h2[0]) & (t3['Z_MEAN'] < bins_h2[-1])
meansMAD3_h1, errorsMAD3_h1 = calculate_means_and_errorsMAD(t3[mask3_h1], bins_h1)
meansMAD3_h2, errorsMAD3_h2 = calculate_means_and_errorsMAD(t3[mask3_h2], bins_h2)
bias3_h1 = np.nanmean(meansMAD3_h1)
err3_h1 = np.nanmean(errorsMAD3_h1)
bias3_h2 = np.nanmean(meansMAD3_h2)
err3_h2 = np.nanmean(errorsMAD3_h2)


# Imprimir resultados
print(f"t1 - Total: {bias1_total:.4f} ± {err1_total:.4f}")
print(f"t1 - 1ª mitad: {bias1_h1:.4f} ± {err1_h1:.4f}")
print(f"t1 - 2ª mitad: {bias1_h2:.4f} ± {err1_h2:.4f}")
print(f"t2 - Total: {bias2_total:.4f} ± {err2_total:.4f}")
print(f"t2 - 1ª mitad: {bias2_h1:.4f} ± {err2_h1:.4f}")
print(f"t2 - 2ª mitad: {bias2_h2:.4f} ± {err2_h2:.4f}")
print(f"t3 - Total: {bias3_total:.4f} ± {err3_total:.4f}")
print(f"t3 - 1ª mitad: {bias3_h1:.4f} ± {err3_h1:.4f}")
print(f"t3 - 2ª mitad: {bias3_h2:.4f} ± {err3_h2:.4f}\n")



def calculate_sigma(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma_values = []
    sigma_err = []
   

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 0:
            sigma_values.append(np.std(dz_bin))
            res = bootstrap((dz_bin,), np.std, n_resamples=100, confidence_level=0.95)
            sigma_err = res.standard_error
        else:
            sigma_values.append(np.nan)  
            sigma_err.append(np.nan)  

    return bin_centers, sigma_values, sigma_err


# Calcular sigma 
bin_centers1, sigma1, sigma_err1 = calculate_sigma(t1, bins)
bin_centers2, sigma2, sigma_err2 = calculate_sigma(t2, bins)
bin_centers3, sigma3, sigma_err3 = calculate_sigma(t3, bins)

mid = len(sigma1) // 2 

sigma1_total, sigma1_err_total = weighted_mean_and_error(sigma1, sigma_err1)
sigma2_total, sigma2_err_total = weighted_mean_and_error(sigma2, sigma_err2)
sigma3_total, sigma3_err_total = weighted_mean_and_error(sigma3, sigma_err3)
# t1
sigma1_h1, sigma1_err_h1 = weighted_mean_and_error(sigma1[:mid], sigma_err1[:mid])
sigma1_h2, sigma1_err_h2 = weighted_mean_and_error(sigma1[mid:], sigma_err1[mid:])
# t2
sigma2_h1, sigma2_err_h1 = weighted_mean_and_error(sigma2[:mid], sigma_err2[:mid])
sigma2_h2, sigma2_err_h2 = weighted_mean_and_error(sigma2[mid:], sigma_err2[mid:])
# t3
sigma3_h1, sigma3_err_h1 = weighted_mean_and_error(sigma3[:mid], sigma_err3[:mid])
sigma3_h2, sigma3_err_h2 = weighted_mean_and_error(sigma3[mid:], sigma_err3[mid:])

print("=== Sigma ===")
print(f"t1 - Total: {sigma1_total:.4f} ± {sigma1_err_total:.4f}")
print(f"t1 - 1ª mitad: {sigma1_h1:.4f} ± {sigma1_err_h1:.4f}")
print(f"t1 - 2ª mitad: {sigma1_h2:.4f} ± {sigma1_err_h2:.4f}")
print(f"t2 - Total: {sigma2_total:.4f} ± {sigma2_err_total:.4f}")
print(f"t2 - 1ª mitad: {sigma2_h1:.4f} ± {sigma2_err_h1:.4f}")
print(f"t2 - 2ª mitad: {sigma2_h2:.4f} ± {sigma2_err_h2:.4f}")
print(f"t3 - Total: {sigma3_total:.4f} ± {sigma3_err_total:.4f}")
print(f"t3 - 1ª mitad: {sigma3_h1:.4f} ± {sigma3_err_h1:.4f}")
print(f"t3 - 2ª mitad: {sigma3_h2:.4f} ± {sigma3_err_h2:.4f}\n")


def calculate_sigma68_with_error(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma68 = []
    sigma68_err = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 10:  
            lower = np.percentile(dz_bin, 16)
            upper = np.percentile(dz_bin, 84)
            sigma = (upper - lower) / 2
            sigma68.append(sigma)

            # Bootstrap para estimar el error
            res = bootstrap((dz_bin,), lambda x: (np.percentile(x, 84) - np.percentile(x, 16)) / 2, 
                            n_resamples=100, confidence_level=0.95)
            sigma68_err.append(res.standard_error)

        else:
            sigma68.append(np.nan)
            sigma68_err.append(np.nan)

    return bin_centers, np.array(sigma68), np.array(sigma68_err)


# Calcular sigma68 y errores para cada archivo
bin_centers1, sigma68_1, sigma68_err1 = calculate_sigma68_with_error(t1, bins)
bin_centers2, sigma68_2, sigma68_err2 = calculate_sigma68_with_error(t2, bins)
bin_centers3, sigma68_3, sigma68_err3 = calculate_sigma68_with_error(t3, bins)

sigma681_total, sigma681_err_total = weighted_mean_and_error(sigma68_1, sigma68_err1)
sigma682_total, sigma682_err_total = weighted_mean_and_error(sigma68_2, sigma68_err2)
sigma683_total, sigma683_err_total = weighted_mean_and_error(sigma68_3, sigma68_err3)

# t1
sigma68_1_h1, sigma68_err1_h1 = weighted_mean_and_error(sigma68_1[:mid], sigma68_err1[:mid])
sigma68_1_h2, sigma68_err1_h2 = weighted_mean_and_error(sigma68_1[mid:], sigma68_err1[mid:])
# t2
sigma68_2_h1, sigma68_err2_h1 = weighted_mean_and_error(sigma68_2[:mid], sigma68_err2[:mid])
sigma68_2_h2, sigma68_err2_h2 = weighted_mean_and_error(sigma68_2[mid:], sigma68_err2[mid:])
# t3
sigma68_3_h1, sigma68_err3_h1 = weighted_mean_and_error(sigma68_3[:mid], sigma68_err3[:mid])
sigma68_3_h2, sigma68_err3_h2 = weighted_mean_and_error(sigma68_3[mid:], sigma68_err3[mid:])

print("=== Sigma₆₈ ===")
print(f"t1 - Total: {sigma681_total:.4f} ± {sigma681_err_total:.4f}")
print(f"t1 - 1ª mitad: {sigma681_h1:.4f} ± {sigma681_err_h1:.4f}")
print(f"t1 - 2ª mitad: {sigma681_h2:.4f} ± {sigma681_err_h2:.4f}")
print(f"t2 - Total: {sigma682_total:.4f} ± {sigma682_err_total:.4f}")
print(f"t2 - 1ª mitad: {sigma682_h1:.4f} ± {sigma682_err_h1:.4f}")
print(f"t2 - 2ª mitad: {sigma682_h2:.4f} ± {sigma682_err_h2:.4f}")
print(f"t3 - Total: {sigma683_total:.4f} ± {sigma683_err_total:.4f}")
print(f"t3 - 1ª mitad: {sigma683_h1:.4f} ± {sigma683_err_h1:.4f}")
print(f"t3 - 2ª mitad: {sigma683_h2:.4f} ± {sigma683_err_h2:.4f}\n")


###outliers de banerji
def calculate_outlier_fraction_Banerji(t, bins, sigma68_values):
    photo_z = t['Z_MEAN']
    delta_z = t['Z_MEAN'] - t['Z']
    outlier_fractions_Banerji = []
    bin_centers_Banerji = (bins[:-1] + bins[1:]) / 2
    outlier_errors_Banerji = []
    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z'][mask]

        if len(dz_bin) > 0:
            threshold = 0.15
            outlier_mask = (np.abs(dz_bin) / (1 + z_mean_bin)) > threshold
            outlier_fraction_Banerji = np.mean(outlier_mask)
            outlier_fractions_Banerji.append(outlier_fraction_Banerji)

            # Calcular error con bootstrap
            res = bootstrap((outlier_mask,), np.mean, n_resamples=100, confidence_level=0.95)
            outlier_errors_Banerji.append(res.standard_error)  # Simétrico arriba y abajo
        else:
            outlier_fractions_Banerji.append(np.nan)
            outlier_errors_Banerji.append(np.nan)  # Evitar errores en bins vacíos


    return bin_centers_Banerji, outlier_fractions_Banerji, outlier_errors_Banerji



bin_centers1_Banerji, outlier_fraction1_Banerji, outlier_errors1_Banerji = calculate_outlier_fraction_Banerji(t1, bins, sigma68_1)
bin_centers2_Banerji, outlier_fraction2_Banerji, outlier_errors2_Banerji = calculate_outlier_fraction_Banerji(t2, bins, sigma68_2)
bin_centers3_Banerji, outlier_fraction3_Banerji, outlier_errors3_Banerji = calculate_outlier_fraction_Banerji(t3, bins, sigma68_3)

outlier1_total, outlier1_err_total = weighted_mean_and_error(outlier_fraction1_Banerji, outlier_errors1_Banerji)
outlier2_total, outlier2_err_total = weighted_mean_and_error(outlier_fraction2_Banerji, outlier_errors2_Banerji)
outlier3_total, outlier3_err_total = weighted_mean_and_error(outlier_fraction3_Banerji, outlier_errors3_Banerji)

# t1
outlier1_h1, outlier_err1_h1 = weighted_mean_and_error(outlier_fraction1_Banerji[:mid], outlier_errors1_Banerji[:mid])
outlier1_h2, outlier_err1_h2 = weighted_mean_and_error(outlier_fraction1_Banerji[mid:], outlier_errors1_Banerji[mid:])
# t2
outlier2_h1, outlier_err2_h1 = weighted_mean_and_error(outlier_fraction2_Banerji[:mid], outlier_errors2_Banerji[:mid])
outlier2_h2, outlier_err2_h2 = weighted_mean_and_error(outlier_fraction2_Banerji[mid:], outlier_errors2_Banerji[mid:])
# t3
outlier3_h1, outlier_err3_h1 = weighted_mean_and_error(outlier_fraction3_Banerji[:mid], outlier_errors3_Banerji[:mid])
outlier3_h2, outlier_err3_h2 = weighted_mean_and_error(outlier_fraction3_Banerji[mid:], outlier_errors3_Banerji[mid:])

print("=== Outlier Fraction (Banerji) ===")
print(f"t1 - Total: {outlier1_total:.4f} ± {outlier1_err_total:.4f}")
print(f"t1 - 1ª mitad: {outlier1_h1:.4f} ± {outlier1_err_h1:.4f}")
print(f"t1 - 2ª mitad: {outlier1_h2:.4f} ± {outlier1_err_h2:.4f}")
print(f"t2 - Total: {outlier2_total:.4f} ± {outlier2_err_total:.4f}")
print(f"t2 - 1ª mitad: {outlier2_h1:.4f} ± {outlier2_err_h1:.4f}")
print(f"t2 - 2ª mitad: {outlier2_h2:.4f} ± {outlier2_err_h2:.4f}")
print(f"t3 - Total: {outlier3_total:.4f} ± {outlier3_err_total:.4f}")
print(f"t3 - 1ª mitad: {outlier3_h1:.4f} ± {outlier3_err_h1:.4f}")
print(f"t3 - 2ª mitad: {outlier3_h2:.4f} ± {outlier3_err_h2:.4f}")





otras métricas


def calculate_means_and_errors(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)  
    delta_z = np.ma.filled((t['Z_MEAN'] - t['Z']), np.nan)  
    means = []
    errors = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 0:
            mean_val = np.mean(dz_bin)
            res = bootstrap((dz_bin,), np.mean, n_resamples=100, confidence_level=0.95)
            error_val = res.standard_error
        else:
            mean_val = np.nan
            error_val = np.nan

        means.append(mean_val)
        errors.append(error_val)

    return means, errors


means1, errors1 = calculate_means_and_errors(t1, bins)
means2, errors2 = calculate_means_and_errors(t2, bins)
means3, errors3 = calculate_means_and_errors(t3, bins)

bin_centers_offset_1=bin_centers-0.01
bin_centers_offset_2=bin_centers+0.01
bin_centers_offset_3=bin_centers+0.03

plt.figure(figsize=(12, 8))
plt.errorbar(bin_centers_offset_1, means1, yerr=errors1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, means2, yerr=errors2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, means3, yerr=errors3, fmt='o', label='vhs J,Ks', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\langle z_\mathrm{photo} - z_\mathrm{spec} \rangle$')
plt.legend(loc='best')
plt.grid()
plt.savefig('3catalogos/F_2_b_vhs.png', dpi=300)
#plt.show()




def calculate_sigma68_with_error(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma68 = []
    sigma68_err = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 10:  
            lower = np.percentile(dz_bin, 16)
            upper = np.percentile(dz_bin, 84)
            sigma = (upper - lower) / 2
            sigma68.append(sigma)

            # Bootstrap para estimar el error
            res = bootstrap((dz_bin,), lambda x: (np.percentile(x, 84) - np.percentile(x, 16)) / 2, 
                            n_resamples=100, confidence_level=0.95)
            sigma68_err.append(res.standard_error)

        else:
            sigma68.append(np.nan)
            sigma68_err.append(np.nan)

    return bin_centers, np.array(sigma68), np.array(sigma68_err)


# Calcular sigma68 y errores para cada archivo
bin_centers1, sigma68_1, sigma68_err1 = calculate_sigma68_with_error(t1, bins)
bin_centers2, sigma68_2, sigma68_err2 = calculate_sigma68_with_error(t2, bins)
bin_centers3, sigma68_3, sigma68_err3 = calculate_sigma68_with_error(t3, bins)
bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
# Graficar sigma_68 vs. z_phot con barras de error
plt.figure(figsize=(12, 8))
plt.errorbar(bin_centers_offset_1, sigma68_1, yerr=sigma68_err1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma68_2, yerr=sigma68_err2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma68_3, yerr=sigma68_err3, fmt='o', label='vhs J,Ks', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma_{68}$')
plt.legend()
plt.grid()
plt.savefig("3catalogos/F_3_b_vhs.png", dpi=300)
#plt.show()


# Función para calcular sigma en cada bin
def calculate_sigma(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma_values = []
    sigma_err = []
   

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 0:
            sigma_values.append(np.std(dz_bin))
            res = bootstrap((dz_bin,), np.std, n_resamples=100, confidence_level=0.95)
            sigma_err = res.standard_error
        else:
            sigma_values.append(np.nan)  
            sigma_err.append(np.nan)  

    return bin_centers, sigma_values, sigma_err


# Calcular sigma 
bin_centers1, sigma1, sigma_err1 = calculate_sigma(t1, bins)
bin_centers2, sigma2, sigma_err2 = calculate_sigma(t2, bins)
bin_centers3, sigma3, sigma_err3 = calculate_sigma(t3, bins)
bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
# Graficar sigma vs. z_phot 
plt.figure(figsize=(12, 8))
plt.errorbar(bin_centers_offset_1, sigma1, yerr=sigma_err1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma2, yerr=sigma_err2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma3, yerr=sigma_err3, fmt='o', label='vhs J,Ks', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma$')
plt.legend()
plt.grid()
#plt.title(r'$\sigma$ en función de $z_{\text{phot}}$')
plt.savefig('3catalogos/F_5_b_vhs.png', dpi=300)
#plt.show()

###prueba de outliers según sigma
def calculate_outlier_fraction(t, bins, sigma68_values):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    outlier_fractions = []
    outlier_errors = []


    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z'][mask]

        if len(dz_bin) > 0:
            threshold = 4 * sigma68_values[i]
            outlier_mask = (np.abs(dz_bin) / (1 + z_mean_bin)) > threshold
            outlier_fraction = np.mean(outlier_mask)
            outlier_fractions.append(outlier_fraction)

           
            res = bootstrap((outlier_mask,), np.mean, n_resamples=100, confidence_level=0.95)
            outlier_errors.append(res.standard_error) 
        else:
            outlier_fractions.append(np.nan)
            outlier_errors.append(np.nan)

    return bin_centers, outlier_fractions, outlier_errors


 
bin_centers1, outlier_fraction1, outlier_errors1 = calculate_outlier_fraction(t1, bins, sigma68_1)
bin_centers2, outlier_fraction2, outlier_errors2 = calculate_outlier_fraction(t2, bins, sigma68_2)
bin_centers3, outlier_fraction3, outlier_errors3 = calculate_outlier_fraction(t3, bins, sigma68_3)
bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
# Graficar outliers vs. z_spec
plt.figure(figsize=(12, 8))
plt.errorbar(bin_centers_offset_1, outlier_fraction1, yerr=outlier_errors1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, outlier_fraction2, yerr=outlier_errors2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, outlier_fraction3, yerr=outlier_errors3, fmt='o', label='vhs J,Ks', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction')
plt.legend()
plt.grid()
#plt.title('Fracción de outliers vs. $z_{\text{spec}}$')
#plt.show()
plt.savefig('3catalogos/F_6_b_vhs.png', dpi=300)





#Otros outliers de prueba
#Banerji
def calculate_outlier_fraction(t, bins, sigma68_values, n_resamples=1000):
    photo_z = np.ma.filled(t['Z_1'], np.nan)
    delta_z = np.ma.filled(t['Z_1'] - t['Z_MEAN'], np.nan)
    outlier_fractions = []
    outlier_errors = []
    bin_centers = (bins[:-1] + bins[1:]) / 2  # Calcular centros de bin

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = np.ma.filled(t['Z_MEAN'][mask], np.nan)  # Evitar problemas con NaN

        if len(dz_bin) > 0:
            threshold = 4 * sigma68_values[i]
            outlier_mask = (np.abs(dz_bin) / (1 + z_mean_bin)) > threshold
            outlier_fraction = np.mean(outlier_mask)
            outlier_fractions.append(outlier_fraction)

            # Calcular error con bootstrap
            res = bootstrap((outlier_mask,), np.mean, n_resamples=n_resamples, confidence_level=0.68)
            outlier_errors.append(res.standard_error)  # Simétrico arriba y abajo
        else:
            outlier_fractions.append(np.nan)
            outlier_errors.append(np.nan)  # Evitar errores en bins vacíos

    return bin_centers, outlier_fractions, outlier_errors

# Calcular fracción de outliers con error
bin_centers1, outlier_fraction1, outlier_errors1 = calculate_outlier_fraction(t1, bins, sigma68_1)
bin_centers2, outlier_fraction2, outlier_errors2 = calculate_outlier_fraction(t2, bins, sigma68_2)

# Graficar con barras de error
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers1, outlier_fraction1, yerr=outlier_errors1, fmt='o', capsize=3, label='dnf', color='blue')
plt.errorbar(bin_centers2, outlier_fraction2, yerr=outlier_errors2, fmt='o', capsize=3, label='anf', color='red')

plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction')
plt.legend()
plt.grid()
plt.show()
def calculate_outlier_fraction(t, bins, sigma68_values):
    photo_z = t['Z_1']
    delta_z = np.abs(t['Z_1'] - t['Z_MEAN'])/(1+photo_z)
    outlier_fractions = []
    bin_centers = (bins[:-1] + bins[1:]) / 2

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z_1'][mask]

        if len(dz_bin) > 0:
            threshold = 2.5 * sigma68_values[i]
            outlier_fraction = np.sum((dz_bin) > threshold) / len(dz_bin)
            outlier_fractions.append(outlier_fraction)
        else:
            outlier_fractions.append(np.nan)  # Evitar errores si el bin está vacío

    return bin_centers, outlier_fractions


# Calcular fracción de outliers para cada conjunto de datos
bin_centers1, outlier_fraction1 = calculate_outlier_fraction(t1, bins, sigma68_1)
bin_centers2, outlier_fraction2 = calculate_outlier_fraction(t2, bins, sigma68_2)

# Graficar fracción de outliers vs. z_spec
plt.figure(figsize=(8, 5))
plt.scatter(bin_centers1, outlier_fraction1, label='dnf', color='blue')
plt.scatter(bin_centers2, outlier_fraction2, label='anf', color='red')
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction')
plt.legend()
plt.grid()
#plt.title('Fracción de outliers vs. $z_{\text{spec}}$')
plt.show()


def calculate_outlier_fraction_Banerji(t, bins, sigma68_values):
    photo_z = t['Z_1']
    delta_z = t['Z_1'] - t['Z_MEAN']
    outlier_fractions_Banerji = []
    bin_centers_Banerji = (bins[:-1] + bins[1:]) / 2

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z_1'][mask]

        if len(dz_bin) > 0:
            threshold = 0.15
            outlier_fraction_Banerji = np.sum((np.abs(dz_bin) / (1 + z_mean_bin)) > threshold) / len(dz_bin)
            outlier_fractions_Banerji.append(outlier_fraction_Banerji)
        else:
            outlier_fractions_Banerji.append(np.nan)  # Evitar errores si el bin está vacío

    return bin_centers_Banerji, outlier_fractions_Banerji


# Calcular fracción de outliers para cada conjunto de datos
bin_centers1_Banerji, outlier_fraction1_Banerji = calculate_outlier_fraction_Banerji(t1, bins, sigma68_1)
bin_centers2_Banerji, outlier_fraction2_Banerji = calculate_outlier_fraction_Banerji(t2, bins, sigma68_2)

# Graficar fracción de outliers vs. z_spec
plt.figure(figsize=(8, 5))
plt.scatter(bin_centers1_Banerji, outlier_fraction1_Banerji, label='dnf', color='blue')
plt.scatter(bin_centers2_Banerji, outlier_fraction2_Banerji, label='anf', color='red')
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction Banerji')
plt.legend()
plt.grid()
#plt.title('outliers vs. $z_{\text{spec}}$')
plt.show()

"""
