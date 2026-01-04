#### 
#EJECUCIoN:: python grafico_imag.py abril_vhs/match_20nov_vhs_columnasbuenas_DES.fits abril_wise/20novmatchwise_cambio_nombres_mag.fits abril_match_wise_vhs/match_wise_vhs.fits 



import numpy as np
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

plt.rcParams.update({
    'axes.labelsize': 26,
    'xtick.labelsize': 22,
    'ytick.labelsize': 22,
    'legend.fontsize': 24,
    'figure.titlesize': 26
})

def plot_mag_histogram(file1, file2, file3):
    data1 = Table.read(file1)
    data2 = Table.read(file2)
    data3 = Table.read(file3)

    mag1 = data1['BDF_MAG_I_CORRECTED']
    mag2 = data2['BDF_MAG_I_CORRECTED']
    mag3 = data3['BDF_MAG_I_CORRECTED_2']

    mag1 = mag1[np.isfinite(mag1)]
    mag2 = mag2[np.isfinite(mag2)]
    mag3 = mag3[np.isfinite(mag3)]

    mag1 = mag1[(mag1 > 9) & (mag1 < 25)]
    mag2 = mag2[(mag2 > 9) & (mag2 < 25)]
    mag3 = mag3[(mag3 > 9) & (mag3 < 25)]

    bins = np.linspace(9, 25, 38)  

    plt.figure(figsize=(12,8))

    plt.hist(mag2, bins=bins, color='mediumpurple', alpha=1, label='DES+WISE', edgecolor='black')
    plt.hist(mag1, bins=bins, color='green', alpha=1, label='DES+VHS', edgecolor='black')
    plt.hist(mag3, bins=bins, color='darkorange', alpha=0.9, label='DES+VHS+WISE', edgecolor='black')
    
    plt.xlabel('i-band magnitude / mag')
    plt.ylabel('Number of sources')


    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(4, 4)) 

    plt.legend()
    plt.tight_layout()
    plt.savefig('imag_2.pdf', dpi=500)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py file1.fits file2.fits file3.fits")
    else:
        plot_mag_histogram(sys.argv[1], sys.argv[2], sys.argv[3])


