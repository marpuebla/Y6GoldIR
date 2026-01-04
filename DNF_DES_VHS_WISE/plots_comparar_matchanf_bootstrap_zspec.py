
#Script para comparar los resultados de z_phot-z_espec frente a los bines de photo_z+
#corregido z_spec y añadido barras de error con bootstrap
#ver los resultados del match wise y vhs y comparar dnf con anf
#Ejecución: 
#python plots_comparar_matchanf_bootstrap_zspec.py anf_match_optico.fits anf_match_vhs.fits anf_match_wise.fits anf_match.fits  

#Mirar los filtros como los pongo y en funcion del resultado los _1 etc

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
file4 = sys.argv[4]
plt.rcParams.update({
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'figure.titlesize': 16
})

# Function to process data
def process_data(file, filters, e_filters):
    mag_name = []
    magerr_name = []
    nfilters = len(filters)
    n_efilters = len(e_filters)

    for f in filters:
        mag_name.append(f'BDF_MAG_{f}')
    for f_e in e_filters:
        magerr_name.append(f'BDF_MAG_ERR_{f_e}')


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
    bias = np.mean(np.abs(t['Z_MEAN'] - t['Z_1']))
    sigma = np.std(t['Z_MEAN'] - t['Z_1'])
    zerrabs = np.abs(t['Z_MEAN'] - t['Z_1'])
    zerrSort = np.sort(zerrabs)
    zerrsigma68 = zerrSort[int(Nvalid * 68 / 100)]

    return t, bias, sigma, zerrsigma68

# Specify filters for each file
filters1 = ['G_CORRECTED_1', 'R_CORRECTED_1', 'I_CORRECTED_1', 'Z_CORRECTED_1']  # NO PQ AŚÍ LOS ERRORES??
filters2 = ['G_CORRECTED_1', 'R_CORRECTED_1', 'I_CORRECTED_1', 'Z_CORRECTED_1', 'J_CORRECTED', 'H_CORRECTED', 'KS_CORRECTED']  # NO PQ AŚÍ LOS ERRORES??
filters3 = ['G_CORRECTED_1', 'R_CORRECTED_1', 'I_CORRECTED_1', 'Z_CORRECTED_1', 'W1_CORRECTED','W2_CORRECTED']  # Example filters for file2
filters4 = ['G_CORRECTED_1', 'R_CORRECTED_1', 'I_CORRECTED_1', 'Z_CORRECTED_1', 'J_CORRECTED', 'H_CORRECTED', 'KS_CORRECTED', 'W1_CORRECTED','W2_CORRECTED']  # Example filters for file2
e_filters1 = ['G_1', 'R_1', 'I_1', 'Z_1']  # NO PQ AŚÍ LOS ERRORES??
e_filters2 = ['G_1', 'R_1', 'I_1', 'Z_1', 'J', 'H', 'KS']  # NO PQ AŚÍ LOS ERRORES??
e_filters3 = ['G_1', 'R_1', 'I_1', 'Z_1','W1','W2']  # Example filters for file2
e_filters4 = ['G_1', 'R_1', 'I_1', 'Z_1', 'J', 'H', 'KS', 'W1','W2']  # Example filters for file2


# Process all three files
t1, bias1, sigma1, zerrsigma68_1 = process_data(file1, filters1, e_filters1)
t2, bias2, sigma2, zerrsigma68_2 = process_data(file2, filters2, e_filters2)
t3, bias3, sigma3, zerrsigma68_3 = process_data(file3, filters3, e_filters3)
t4, bias4, sigma4, zerrsigma68_4 = process_data(file4, filters4, e_filters4)

print('longitudes:' , len(t1), len(t2))
print('bias:' , bias1, bias2, bias3, bias4)
print('sigma:' , sigma1, sigma2,sigma3, sigma4)
print('sigma68:' , zerrsigma68_1, zerrsigma68_2,  zerrsigma68_3, zerrsigma68_4)
# Create the first plot: z vs DNF_Z for all three files
plt.figure(figsize=(8, 5))
plt.hist2d(t1['Z_1'], t1['Z_MEAN'], bins=1000, norm=LogNorm(), label='optical', alpha=1)
plt.hist2d(t2['Z_1'], t2['Z_MEAN'], bins=1000, norm=LogNorm(), label='vhs', alpha=1)
plt.hist2d(t3['Z_1'], t3['Z_MEAN'], bins=1000, norm=LogNorm(), label='wise', alpha=1)
plt.hist2d(t4['Z_1'], t4['Z_MEAN'], bins=1000, norm=LogNorm(), label='vhs+wise', alpha=1)
plt.colorbar()
pylab.xlim(0.0, 1.5)
pylab.ylim(0.0, 1.5)
pylab.xlabel('Z_1')
pylab.ylabel('DNF_Z')
plt.legend(loc='best')
#plt.title(f"\ndnf: bias={bias1:.4f}, sigma={sigma1:.4f}\nanf: bias={bias2:.4f}, sigma={sigma2:.4f}", fontsize=8)
plt.savefig('PRUEBAS/F_1_b_match_completo.png', dpi=300)
#pylab.show()


# Create a scatter plot comparing the mean delta_z for all three files in bins
bins = np.linspace(0.2, 1.4, 10)  # Define the bins
bin_centers = (bins[:-1] + bins[1:]) / 2  # Bin centers

def calculate_means_and_errorsMAD(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)  
    delta_z = np.ma.filled((t['Z_MEAN'] - t['Z_1']), np.nan)  
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
meansMAD4, errorsMAD4 = calculate_means_and_errorsMAD(t4, bins)

bin_centers_offset_1_MAD=bin_centers-0.01
bin_centers_offset_2_MAD=bin_centers+0.01
bin_centers_offset_3_MAD=bin_centers+0.03
bin_centers_offset_4_MAD=bin_centers+0.05
# Plot comparison
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1_MAD, meansMAD1, yerr=errorsMAD1, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2_MAD, meansMAD2, yerr=errorsMAD2, fmt='D', label='DES+VHS', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3_MAD, meansMAD3, yerr=errorsMAD3, fmt='^', label='DES+WISE', color='mediumpurple', capsize=4)
plt.errorbar(bin_centers_offset_4_MAD, meansMAD4, yerr=errorsMAD4, fmt='*', label='DES+VHS+WISE', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\langle |z_\mathrm{photo} - z_\mathrm{spec}| \rangle$')
plt.legend(loc='best')
#plt.grid()
plt.savefig('PRUEBAS/F_8_b_match.pdf', dpi=300)
#plt.show()
def calculate_sigma68_with_error(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z_1'], np.nan)
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
bin_centers4, sigma68_4, sigma68_err4 = calculate_sigma68_with_error(t4, bins)

# Función para calcular sigma_68 / (1 + Z_MEAN)
def calculate_sigma68_ratio(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z_1'], np.nan)
    sigma68_ratio = []
    sigma68_ratio_err = []

    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z_1'][mask]

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
bin_centers4, sigma68_ratio4, sigma68_ratio_err4 = calculate_sigma68_ratio(t4, bins)

bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
bin_centers_offset_4=bin_centers4+0.05
# Graficar sigma_68 / (1 + Z_MEAN) vs. z_phot para los tres archivos
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, sigma68_ratio1, yerr=sigma68_ratio_err1, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma68_ratio2, yerr=sigma68_ratio_err2, fmt='D', label='DES+VHS', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma68_ratio3, yerr=sigma68_ratio_err3, fmt='^', label='DES+WISE', color='mediumpurple', capsize=4)
plt.errorbar(bin_centers_offset_4, sigma68_ratio4, yerr=sigma68_ratio_err4, fmt='*', label='DES+VHS+WISE', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma_{68} / (1 + z_\mathrm{spec})$')
plt.legend()
#plt.grid()
plt.savefig('PRUEBAS/F_4_b.png', dpi=300)
#plt.title(r'Relación $\sigma_{68} / (1 + Z_{\text{MEAN}})$ en función de $z_{\text{phot}}$')


def calculate_outlier_fraction_Banerji(t, bins, sigma68_values):
    photo_z = t['Z_MEAN']
    delta_z = t['Z_MEAN'] - t['Z_1']
    outlier_fractions_Banerji = []
    bin_centers_Banerji = (bins[:-1] + bins[1:]) / 2
    outlier_errors_Banerji = []
    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z_1'][mask]

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
bin_centers4_Banerji, outlier_fraction4_Banerji, outlier_errors4_Banerji = calculate_outlier_fraction_Banerji(t4, bins, sigma68_4)
bin_centers_offset_1=bin_centers1_Banerji-0.01
bin_centers_offset_2=bin_centers2_Banerji+0.01
bin_centers_offset_3=bin_centers3_Banerji+0.03
bin_centers_offset_4=bin_centers4_Banerji+0.05
# Graficar fracción de outliers vs. z_spec
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, outlier_fraction1_Banerji, yerr=outlier_errors1_Banerji, fmt='o', label='DES', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, outlier_fraction2_Banerji, yerr=outlier_errors2_Banerji, fmt='D', label='DES+VHS', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, outlier_fraction3_Banerji, yerr=outlier_errors3_Banerji, fmt='^', label='DES+WISE', color='mediumpurple', capsize=4)
plt.errorbar(bin_centers_offset_4, outlier_fraction4_Banerji, yerr=outlier_errors4_Banerji, fmt='*', label='DES+VHS+WISE', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction (Banerji)')
plt.legend()
#plt.grid()
#plt.title('outliers banerji vs. $z_{\text{spec}}$')
#plt.show()
plt.savefig('PRUEBAS/F_7_banerji_b_match.png', dpi=300)
"""

# Calculate mean delta_z for each bin in all three files
def calculate_means_and_errors(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)  
    delta_z = np.ma.filled((t['Z_MEAN'] - t['Z_1']), np.nan)  
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
means4, errors4 = calculate_means_and_errors(t4, bins)
bin_centers_offset_1=bin_centers-0.01
bin_centers_offset_2=bin_centers+0.01
bin_centers_offset_3=bin_centers+0.03
bin_centers_offset_4=bin_centers+0.05
# Plot comparison
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, means1, yerr=errors1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, means2, yerr=errors2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, means3, yerr=errors3, fmt='o', label='wise', color='green', capsize=4)
plt.errorbar(bin_centers_offset_4, means4, yerr=errors4, fmt='o', label='vhs+wise', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\langle z_\mathrm{photo} - z_\mathrm{spec} \rangle$')
plt.legend(loc='best')
plt.grid()
plt.savefig('PRUEBAS/F_2_b.png', dpi=300)
#plt.show()

bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
bin_centers_offset_4=bin_centers4+0.05
# Graficar sigma_68 vs. z_phot con barras de error
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, sigma68_1, yerr=sigma68_err1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma68_2, yerr=sigma68_err2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma68_3, yerr=sigma68_err3, fmt='o', label='wise', color='green', capsize=4)
plt.errorbar(bin_centers_offset_4, sigma68_4, yerr=sigma68_err4, fmt='o', label='vhs+wise', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma_{68}$')
plt.legend()
plt.grid()
plt.savefig("PRUEBAS/F_3_b.png", dpi=300)
#plt.show()


#plt.show()


# Función para calcular sigma en cada bin
def calculate_sigma(t, bins):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z_1'], np.nan)
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
            sigma_values.append(np.nan)  # Para evitar errores si no hay datos en el bin
            sigma_err.append(np.nan)  # Para evitar errores si el bin está vacío

    return bin_centers, sigma_values, sigma_err


# Calcular sigma para cada archivo
bin_centers1, sigma1, sigma_err1 = calculate_sigma(t1, bins)
bin_centers2, sigma2, sigma_err2 = calculate_sigma(t2, bins)
bin_centers3, sigma3, sigma_err3 = calculate_sigma(t3, bins)
bin_centers4, sigma4, sigma_err4 = calculate_sigma(t4, bins)

bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
bin_centers_offset_4=bin_centers4+0.05

# Graficar sigma vs. z_phot para los tres archivos
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, sigma1, yerr=sigma_err1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, sigma2, yerr=sigma_err2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, sigma3, yerr=sigma_err3, fmt='o', label='wise', color='green', capsize=4)
plt.errorbar(bin_centers_offset_4, sigma4, yerr=sigma_err4, fmt='o', label='vhs+wise', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel(r'$\sigma$')
plt.legend()
plt.grid()
#plt.title(r'$\sigma$ en función de $z_{\text{phot}}$')
plt.savefig('PRUEBAS/F_5_b.png', dpi=300)
#plt.show()

def calculate_outlier_fraction(t, bins, sigma68_values):
    photo_z = np.ma.filled(t['Z_MEAN'], np.nan)
    delta_z = np.ma.filled(t['Z_MEAN'] - t['Z_1'], np.nan)
    outlier_fractions = []
    outlier_errors = []


    for i in range(len(bins) - 1):
        mask = (photo_z >= bins[i]) & (photo_z < bins[i + 1])
        dz_bin = delta_z[mask]
        z_mean_bin = t['Z_1'][mask]

        if len(dz_bin) > 0:
            threshold = 4 * sigma68_values[i]
            outlier_mask = (np.abs(dz_bin) / (1 + z_mean_bin)) > threshold
            outlier_fraction = np.mean(outlier_mask)
            outlier_fractions.append(outlier_fraction)

            # Calcular error con bootstrap
            res = bootstrap((outlier_mask,), np.mean, n_resamples=100, confidence_level=0.95)
            outlier_errors.append(res.standard_error)  # Simétrico arriba y abajo
        else:
            outlier_fractions.append(np.nan)
            outlier_errors.append(np.nan)  # Evitar errores en bins vacíos

    return bin_centers, outlier_fractions, outlier_errors


# Calcular fracción de outliers para cada conjunto de datos
bin_centers1, outlier_fraction1, outlier_errors1 = calculate_outlier_fraction(t1, bins, sigma68_1)
bin_centers2, outlier_fraction2, outlier_errors2 = calculate_outlier_fraction(t2, bins, sigma68_2)
bin_centers3, outlier_fraction3, outlier_errors3 = calculate_outlier_fraction(t3, bins, sigma68_3)
bin_centers4, outlier_fraction4, outlier_errors4 = calculate_outlier_fraction(t4, bins, sigma68_4)
bin_centers_offset_1=bin_centers1-0.01
bin_centers_offset_2=bin_centers2+0.01
bin_centers_offset_3=bin_centers3+0.03
bin_centers_offset_4=bin_centers4+0.05
# Graficar fracción de outliers vs. z_spec
plt.figure(figsize=(8, 5))
plt.errorbar(bin_centers_offset_1, outlier_fraction1, yerr=outlier_errors1, fmt='o', label='optical', color='blue', capsize=4)
plt.errorbar(bin_centers_offset_2, outlier_fraction2, yerr=outlier_errors2, fmt='o', label='vhs', color='red', capsize=4)
plt.errorbar(bin_centers_offset_3, outlier_fraction3, yerr=outlier_errors3, fmt='o', label='wise', color='green', capsize=4)
plt.errorbar(bin_centers_offset_4, outlier_fraction4, yerr=outlier_errors4, fmt='o', label='vhs+wise', color='orange', capsize=4)
plt.xlabel('Photo-z bin')
plt.ylabel('Outlier Fraction')
plt.legend()
plt.grid()
#plt.title('Fracción de outliers vs. $z_{\text{spec}}$')
#plt.show()
plt.savefig('PRUEBAS/F_6_b.png', dpi=300)


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
