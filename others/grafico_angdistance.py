
##EJECUCIÖN:: python grafico_angdistance.py abril_vhs/match_20nov_vhs_columnasbuenas_DES.fits abril_wise/20novmatchwise_cambio_nombres_mag.fits
import numpy as np
from astropy.table import Table
import sys
import matplotlib.pyplot as plt
plt.rcParams.update({
    'axes.labelsize': 26,
    'xtick.labelsize': 22,
    'ytick.labelsize': 22,
    'legend.fontsize': 24,
    'figure.titlesize': 26
})

def plot_angdist_histogram(file1, file2):
    # Leer los archivos
    data1 = Table.read(file1)
    data2 = Table.read(file2)

    # Extraer la columna 'angdist'
    angdist1 = data1['angDist']
    angdist2 = data2['angDist']

    # Crear histograma
    plt.figure(figsize=(12, 8))
    plt.hist(angdist2, bins=50, color='darkorange', alpha=0.6, label='DES+WISE', edgecolor='black')
    plt.hist(angdist1, bins=50, color='mediumpurple', alpha=0.6, label='DES+VHS', edgecolor='black')
    

    # Etiquetas y leyenda
    plt.xlabel('Angular Distance / arcsec')
    plt.ylabel('Number of matches')
    plt.legend()
    #plt.title('Histogram of Angular Distances')
    #plt.grid(True)
    plt.tight_layout()
    plt.savefig('angDist.png',dpi=500)
    plt.show()
# Ejecutar si el script se llama directamente
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py file1.fits file2.fits")
    else:
        plot_angdist_histogram(sys.argv[1], sys.argv[2])

