#Script para comparar los resultados de z_phot-z_espec frente a los bines de photo_z+
#corregido z_spec y añadido barras de error con bootstrap
#ver los resultados del match wise y vhs y comparar dnf con anf
#Ejecución: 
#python anf_plots_comparar_wise_bootstrap_zspec_y.py anf_wise_bandas_optico.fits anf_wise_W1_W2.fits anf_wise_cony.fits



import numpy as np
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
import pylab
import weighted_kde as kde
from scipy.stats import bootstrap
from matplotlib.colors import LogNorm

# Read three files
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

# Function to process data
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

# Specify filters for each file
filters1 = ['G', 'R', 'I', 'Z']  # Example filters for file1
filters2 = ['G', 'R', 'I', 'Z', 'W1', 'W2']  # Example filters for file2
filters3 = ['G', 'R', 'I', 'Z', 'Y','W1', 'W2']  # Example filters for file3




# Process all three files
t1, bias1, sigma1, zerrsigma68_1 = process_data(file1, filters1)
t2, bias2, sigma2, zerrsigma68_2 = process_data(file2, filters2)
t3, bias3, sigma3, zerrsigma68_3 = process_data(file3, filters3)
print('longitudes:' , len(t1), len(t2))
print('bias:' , bias1, bias2)
print('sigma:' , sigma1, sigma2)
print('sigma68:' , zerrsigma68_1, zerrsigma68_2)

# Create the first plot: z vs DNF_Z for all three files
plt.figure(figsize=(8, 5))
plt.hist2d(t1['Z'], t1['Z_MEAN'], bins=1000, norm=LogNorm(), label='optical', alpha=1)
plt.hist2d(t2['Z'], t2['Z_MEAN'], bins=1000, norm=LogNorm(), label='wise', alpha=1)
plt.hist2d(t3['Z'], t3['Z_MEAN'], bins=1000, norm=LogNorm(), label='y', alpha=1)
plt.colorbar()
pylab.xlim(0.0, 1.5)
pylab.ylim(0.0, 1.5)
pylab.xlabel('Z')
pylab.ylabel('DNF_Z')
plt.legend(loc='best')
#plt.title(f"\ndnf: bias={bias1:.4f}, sigma={sigma1:.4f}\nanf: bias={bias2:.4f}, sigma={sigma2:.4f}", fontsize=8)
plt.savefig('cony/F_1_b_wise.png', dpi=300)
#pylab.show()




# Create a scatter plot comparing the mean delta_z for all three files in bins
bins = np.linspace(0.2, 1.4, 10)  # Define the bins
bin_centers = (bins[:-1] + bins[1:]) / 2  # Bin centers


# Calculate mean delta_z for each bin in all three files
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

bin_centers_offset_1_MAD=bin_centers-0.01
bin_centers_offset_2_MAD=bin_centers+0.01
bin_centers_offset_3_MAD=bin_centers+0.03
# Plot comparison
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1_MAD, meansMAD1, yerr=errorsMAD1, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2_MAD, meansMAD2, yerr=errorsMAD2, fmt='o', label='DES+WISE', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3_MAD, meansMAD3, yerr=errorsMAD3, fmt='o', label='DES+WISE (with Y)', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\langle |z_\mathrm{photo} - z_\mathrm{spec}| \rangle$')
plt.legend(loc='best')
plt.grid()
plt.savefig('cony/F_8_b_wise.png', dpi=300)
#plt.show()


def calculate_sigma68_with_error(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma68 = []
    sigma68_err = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]

        if len(dz_bin) > 10:  # Evita bins con pocos datos
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


# Función para calcular sigma_68 / (1 + Z_MEAN)
def calculate_sigma68_ratio(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z'], np.nan)
    sigma68_ratio = []
    sigma68_ratio_err = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z'][mask]

        if len(dz_bin) > 0:
            lower = np.percentile(dz_bin, 16)
            upper = np.percentile(dz_bin, 84)
            sigma_68 = (upper - lower) / 2
            sigma68_ratio.append(sigma_68 / np.mean(1 + z_mean_bin))
            # Aplicamos el bootstrap usando lambda para calcular sigma_68 y dividir por (1 + Z_MEAN)
            res = bootstrap(
                (dz_bin,), 
                lambda x: (np.percentile(x, 84) - np.percentile(x, 16)) / 2 / np.mean(1 + z_mean_bin),
                n_resamples=100, 
                confidence_level=0.95
            )


            # Almacenar la relación sigma68 / (1 + Z_MEAN) para el bin
            sigma68_ratio_err.append(res.standard_error)

        else:
            sigma68_ratio.append(np.nan)  # Para evitar errores si no hay datos en el bin
            sigma68_ratio_err.append(np.nan)  # Si no hay datos, el error también es NaN

    return bin_centers, sigma68_ratio, sigma68_ratio_err



# Calcular sigma68 / (1 + Z_MEAN) para cada archivo
bin_centers1, sigma68_ratio1, sigma68_ratio_err1 = calculate_sigma68_ratio(t1, bins)
bin_centers2, sigma68_ratio2, sigma68_ratio_err2 = calculate_sigma68_ratio(t2, bins)
bin_centers3, sigma68_ratio3, sigma68_ratio_err3 = calculate_sigma68_ratio(t3, bins)

bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
# Graficar sigma_68 / (1 + Z_MEAN) vs. z_phot para los tres archivos
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, sigma68_ratio1, yerr=sigma68_ratio_err1, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma68_ratio2, yerr=sigma68_ratio_err2, fmt='o', label='DES+WISE', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma68_ratio3, yerr=sigma68_ratio_err3, fmt='o', label='DES+WISE (with Y)', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma_{68} / (1 + z_\mathrm{spec})$')
plt.legend()
plt.grid()
plt.savefig('cony/F_4_b_wise.png', dpi=300)
#plt.title(r'Relación $\sigma_{68} / (1 + Z_{\text{MEAN}})$ en función de $z_{\text{phot}}$')


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


# Calcular fracción de outliers para cada conjunto de datos
bin_centers1_Banerji, outlier_fraction1_Banerji, outlier_errors1_Banerji = calculate_outlier_fraction_Banerji(t1, bins, sigma68_1)
bin_centers2_Banerji, outlier_fraction2_Banerji, outlier_errors2_Banerji = calculate_outlier_fraction_Banerji(t2, bins, sigma68_2)
bin_centers3_Banerji, outlier_fraction3_Banerji, outlier_errors3_Banerji = calculate_outlier_fraction_Banerji(t3, bins, sigma68_3)
bin_centers_offset_1=bin_centers1_Banerji-0.01
bin_centers_offset_2=bin_centers2_Banerji+0.01
bin_centers_offset_3=bin_centers3_Banerji+0.03
# Graficar fracción de outliers vs. z_spec
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, outlier_fraction1_Banerji, yerr=outlier_errors1_Banerji, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, outlier_fraction2_Banerji, yerr=outlier_errors2_Banerji, fmt='o', label='DES+WISE', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, outlier_fraction3_Banerji, yerr=outlier_errors3_Banerji, fmt='o', label='DES+WISE (with Y)', color='green', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction (Banerji)')
plt.legend()
plt.grid()
#plt.title('outliers banerji vs. $z_{\text{spec}}$')
#plt.show()
plt.savefig('cony/F_7_banerji_b_wise.png', dpi=300)
