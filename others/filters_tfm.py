###### Ejecucion::: python filters_tfm.py 

import matplotlib.pyplot as plt
import speclite as speclite
from speclite import filters
import astropy
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from speclite import filters

plt.rcParams.update({
    'axes.labelsize': 26,
    'xtick.labelsize': 22,
    'ytick.labelsize': 22,
    'legend.fontsize': 24,
    'figure.titlesize': 26
})
wise = filters.load_filters('wise2010-*')
decam = filters.load_filters('decamDR1-*')

plt.figure(figsize=(14, 8))

colors_decam = ['yellow', 'orange', 'magenta', 'steelblue', 'blue']
labels_decam = ['g', 'r', 'i', 'z', 'Y']
for decam_filter, color, label in zip(decam[:5], colors_decam, labels_decam):
    wl_micron = decam_filter.wavelength / 10000.0
    plt.fill_between(
        wl_micron,
        decam_filter.response,
        facecolor=color,
        edgecolor=color,
        linewidth=3,
        alpha=0.3,
        label=label
    )
    plt.plot(wl_micron, decam_filter.response, color=color)

colors_wise = ['limegreen', 'cyan', 'blueviolet', 'deeppink']
labels_wise = ['W1', 'W2', 'W3', 'W4']
for wise_filter, color, label in zip(wise[:4], colors_wise, labels_wise):
    wl_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(
        wl_micron,
        wise_filter.response,
        facecolor=color,
        edgecolor=color,
        linewidth=3,
        alpha=0.3,
        label=label
    )
    plt.plot(wl_micron, wise_filter.response, color=color)

plt.xscale('log')
plt.xlim(0.3, 30)
plt.ylim(0, 1.05)
plt.xlabel('Wavelength / μm')
plt.ylabel('Filter Response')
#plt.title('Filtros DES (DECam) + WISE combinados', fontsize=18)
plt.legend(loc='upper right')
plt.grid(False)
plt.tight_layout()
plt.savefig('wisedecam.pdf',dpi=500)
plt.show()


wise = filters.load_filters('wise2010-*')
decam = filters.load_filters('decamDR1-*')
data_Y = np.loadtxt('Paranal_VISTA.Y.dat')
data_J = np.loadtxt('Paranal_VISTA.J.dat')
data_H = np.loadtxt('Paranal_VISTA.H.dat')
data_Ks = np.loadtxt('Paranal_VISTA.Ks.dat')

plt.figure(figsize=(14, 8))

colors_decam = ['yellow', 'orange', 'magenta', 'steelblue', 'blue']
labels_decam = ['g', 'r', 'i', 'z', 'Y_DES']
for decam_filter, color, label in zip(decam[:5], colors_decam, labels_decam):
    wl_micron = decam_filter.wavelength / 10000.0
    plt.fill_between(wl_micron, decam_filter.response, facecolor=color, edgecolor=color, linewidth=3, alpha=0.3, label=label)
    plt.plot(wl_micron, decam_filter.response, color=color)

colors_vista = ['red', 'green', 'blue', 'purple']
labels_vista = ['Y', 'J', 'H', 'Ks']
for data, color, label in zip([data_Y, data_J, data_H, data_Ks], colors_vista, labels_vista):
    wl_micron = data[:, 0] / 10000.0
    plt.fill_between(wl_micron, data[:, 1], facecolor=color, edgecolor=color, linewidth=3, alpha=0.3, label=label)
    plt.plot(wl_micron, data[:, 1], color=color)

plt.xscale('log')
plt.xlim(0.3, 3)
plt.ylim(0, 1.05)
plt.xlabel('Wavelength / μm')
plt.ylabel('Filter Response')
#plt.title('Filtros DES (DECam) + VHS (VISTA) combinados', fontsize=18)
plt.legend(loc='upper right')
plt.grid(False)
plt.tight_layout()
plt.savefig('todos.pdf',dpi=500)
plt.show()




wise = filters.load_filters('wise2010-*')
decam = filters.load_filters('decamDR1-*')
data_Y = np.loadtxt('Paranal_VISTA.Y.dat')
data_J = np.loadtxt('Paranal_VISTA.J.dat')
data_H = np.loadtxt('Paranal_VISTA.H.dat')
data_Ks = np.loadtxt('Paranal_VISTA.Ks.dat')

plt.figure(figsize=(14, 8))

colors_decam = ['yellow', 'orange', 'magenta', 'steelblue', 'blue']
labels_decam = ['g', 'r', 'i', 'z', 'Y_DES']
for decam_filter, color, label in zip(decam[:5], colors_decam, labels_decam):
    wl_micron = decam_filter.wavelength / 10000.0
    plt.fill_between(wl_micron, decam_filter.response, facecolor=color, edgecolor=color, linewidth=3, alpha=0.3, label=label)
    plt.plot(wl_micron, decam_filter.response, color=color)

colors_vista = ['red', 'green', 'blue', 'purple']
labels_vista = ['Y', 'J', 'H', 'Ks']
for data, color, label in zip([data_Y, data_J, data_H, data_Ks], colors_vista, labels_vista):
    wl_micron = data[:, 0] / 10000.0
    plt.fill_between(wl_micron, data[:, 1], facecolor=color, edgecolor=color, linewidth=3, alpha=0.3, label=label)
    plt.plot(wl_micron, data[:, 1], color=color)

colors_wise = ['limegreen', 'cyan', 'blueviolet', 'deeppink']
labels_wise = ['W1', 'W2', 'W3', 'W4']
for wise_filter, color, label in zip(wise[:4], colors_wise, labels_wise):
    wl_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(wl_micron, wise_filter.response, facecolor=color, edgecolor=color, linewidth=3, alpha=0.3, label=label)
    plt.plot(wl_micron, wise_filter.response, color=color)

plt.xscale('log')
plt.xlim(0.3, 30)
plt.ylim(0, 1.05)
plt.xlabel('Wavelength / μm')
plt.ylabel('Filter Response')
#plt.title('Filtros DECam + VISTA + WISE combinados', fontsize=18)
plt.legend(loc='upper right')
plt.grid(False)
plt.tight_layout()
plt.savefig('desvhs.pdf',dpi=500)
plt.show()




wise = filters.load_filters('wise2010-*')
decam = filters.load_filters('decamDR1-*')

plt.figure(figsize=(14, 8))

colors_wise = ['limegreen', 'cyan', 'blueviolet', 'deeppink']
colors_decam = ['yellow', 'orange', 'magenta', 'steelblue', 'lightgreen']

for wise_filter, color in zip(wise, colors_wise):
    wl_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(wl_micron, wise_filter.response, color=color, alpha=0.4, label=wise_filter.name)
    plt.plot(wl_micron, wise_filter.response, color=color)

for decam_filter, color in zip(decam, colors_decam):
    wl_micron = decam_filter.wavelength / 10000.0
    plt.fill_between(wl_micron, decam_filter.response, color=color, alpha=0.4, label=decam_filter.name)
    plt.plot(wl_micron, decam_filter.response, color=color)

plt.xlim(0.3, 30)
plt.ylim(0, 1.05)
plt.xlabel('Wavelength / μm', fontsize=16)
plt.ylabel('Filter Response', fontsize=16)
plt.title('Filtros DES (DECam) + WISE combinados', fontsize=18)
plt.legend(loc='upper right', fontsize=10)
plt.grid(False)
plt.tight_layout()
plt.show()




wise = filters.load_filters('wise2010-*')
colors = ['limegreen', 'cyan', 'blueviolet' , 'deeppink']
plt.figure(figsize=(12, 8))
for wise_filter, color in zip(wise, colors):
    # Convertir longitud de onda de Angstroms a micrones
    wavelength_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(wavelength_micron, wise_filter.response, color=color, alpha=0.4, label=wise_filter.name)
    plt.plot(wavelength_micron, wise_filter.response, color=color)
plt.xscale('log')
plt.gca().set_xticks([2.5, 5, 10, 20, 30])
plt.gca().set_xticklabels([2.5, 5, 10, 20, 30])
plt.xlim(2.5, 30)
plt.ylim(0,1)
plt.xlabel('Wavelength / $\mu m$', fontsize=26, labelpad=10)
plt.ylabel('Filter Response', fontsize=26, labelpad=15)
plt.title('WISE', fontsize=26)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::1], ['W1', 'W2', 'W3', 'W4'], loc='upper right')
plt.grid(False)
plt.subplots_adjust(bottom=0.17)
plt.savefig('wisefilters.pdf',dpi=500)
plt.show()
plt.close()




decam = filters.load_filters('decamDR1-*')
colors = ['yellow', 'limegreen', 'cyan', 'blueviolet' , 'deeppink']
plt.figure(figsize=(12, 8))
for decam_filter, color in zip(decam, colors):
    # Convertir longitud de onda de Angstroms a micrones
    #wavelength_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(decam_filter.wavelength, decam_filter.response, color=color, alpha=0.4, label=wise_filter.name)
    plt.plot(decam_filter.wavelength, decam_filter.response, color=color)
#plt.xscale('log')
#plt.gca().set_xticks([2.5, 5, 10, 20, 30])
#plt.gca().set_xticklabels([2.5, 5, 10, 20, 30])
#plt.xlim(2.5, 30)
plt.ylim(0,0.6)
plt.xlabel('Wavelength / Angstrom', fontsize=26, labelpad=10)
plt.ylabel('Filter Response', fontsize=26, labelpad=15)
plt.title('DECAM', fontsize=26)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::1], ['g', 'r', 'i', 'z','Y'], loc='upper left')
plt.grid(False)
plt.subplots_adjust(bottom=0.17)
plt.savefig('decamfilters.pdf',dpi=500)
plt.show()
plt.close()


#### VHS directamente leyendo los archivos

# Leer los archivos .dat
data_Y = np.loadtxt('Paranal_VISTA.Y.dat')
data_J = np.loadtxt('Paranal_VISTA.J.dat')
data_H = np.loadtxt('Paranal_VISTA.H.dat')
data_Ks = np.loadtxt('Paranal_VISTA.Ks.dat')

# Crear el gráfico
plt.figure(figsize=(12, 8))

# Colores para cada filtro
colors = ['limegreen', 'cyan', 'blueviolet', 'deeppink']

# Graficar los tres filtros sin añadir ceros
plt.fill_between(data_Y[:, 0], data_Y[:, 1], color=colors[0], alpha=0.4, label='Y')
plt.plot(data_Y[:, 0], data_Y[:, 1], color=colors[0])

plt.fill_between(data_J[:, 0], data_J[:, 1], color=colors[1], alpha=0.4, label='J')
plt.plot(data_J[:, 0], data_J[:, 1], color=colors[1])

plt.fill_between(data_H[:, 0], data_H[:, 1], color=colors[2], alpha=0.4, label='H')
plt.plot(data_H[:, 0], data_H[:, 1], color=colors[2])

plt.fill_between(data_Ks[:, 0], data_Ks[:, 1], color=colors[3], alpha=0.4, label='Ks')
plt.plot(data_Ks[:, 0], data_Ks[:, 1], color=colors[3])

# Configuración del gráfico
plt.ylim(0, 1)
plt.xlabel('Wavelength / Angstrom', fontsize=26, labelpad=10)
plt.ylabel('Filter Response', fontsize=26, labelpad=15)
plt.title('VISTA', fontsize=26)

# Mostrar leyenda
plt.legend(loc='upper right')
plt.grid(False)
plt.subplots_adjust(bottom=0.17)
# Guardar el gráfico
plt.savefig('vistafilters_no_modification.pdf',dpi=500)
plt.show()
plt.close()

"""
wise = speclite.filters.load_filters('wise2010-*')
speclite.filters.plot_filters(wise, wavelength_limits=(2, 30),
    wavelength_unit=astropy.units.micron, wavelength_scale='log')
plt.gca().set_xticks([2, 5, 10, 20, 30])
plt.gca().set_xticklabels([2, 5, 10, 20, 30])

###DES con extinción
with_atm = speclite.filters.load_filters('decamDR1-*')
speclite.filters.plot_filters(with_atm)

###Para nacho:

# Estilo


# Cargar filtros
euclid = filters.load_filters('Euclid-VIS', 'Euclid-Y', 'Euclid-J', 'Euclid-H')
decam = filters.load_filters('decamDR1-*')
all_filter_names = euclid.names + decam.names
all_filters = filters.load_filters(*all_filter_names)

# Crear figura más grande
plt.figure(figsize=(12, 6))

# Graficar filtros
filters.plot_filters(all_filters)

# Obtener el objeto Axes actual
ax = plt.gca()

# Nombres personalizados para la leyenda
custom_labels = [
    'Euclid VIS', 'Euclid Y', 'Euclid J', 'Euclid H',
    'DECam g', 'DECam r', 'DECam i', 'DECam z', 'DECam Y'
]

# Obtener las líneas graficadas (en el mismo orden en que se grafican)
lines = ax.get_lines()

# Crear handles y leyenda con los mismos colores de las curvas
handles = [Line2D([0], [0], color=line.get_color(), lw=2) for line in lines[:len(custom_labels)]]

# Añadir la leyenda con colores correctos
ax.legend(handles, custom_labels, loc='upper right')



# Mostrar gráfico
plt.show()

####cosas mal
import matplotlib.pyplot as plt
import speclite as speclite
from speclite import filters
import astropy
import numpy as np
from astropy import units as u
###wise
plt.rcParams.update({
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 18,
    'figure.titlesize': 18
    
})


import matplotlib.pyplot as plt
from speclite import filters

# Estilo
plt.rcParams.update({
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14,
    'figure.titlesize': 18
})

import matplotlib.pyplot as plt
from speclite import filters
from matplotlib.lines import Line2D  # Para crear handles de la leyenda

# Estilo
plt.rcParams.update({
    'axes.labelsize': 18,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14,
    'figure.titlesize': 18
})

# Cargar filtros
euclid = filters.load_filters('Euclid-VIS', 'Euclid-Y', 'Euclid-H', 'Euclid-J')
decam = filters.load_filters('decamDR1-*')
all_filter_names = euclid.names + decam.names
all_filters = filters.load_filters(*all_filter_names)

# Crear figura más grande
plt.figure(figsize=(12, 6))

# Graficar filtros (esto dibuja las curvas, pero no devuelve los handles)
filters.plot_filters(all_filters)

# Obtener el objeto Axes actual
ax = plt.gca()

# Nombres personalizados para la leyenda
custom_labels = [
    'Euclid VIS', 'Euclid Y', 'Euclid H', 'Euclid J',
    'DECam g', 'DECam r', 'DECam i', 'DECam z', 'DECam Y'
]

# Crear handles manuales para la leyenda con colores que usa speclite
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
handles = [
    Line2D([0], [0], color=color, lw=2) for color in colors[:len(custom_labels)]
]

# Añadir la leyenda
ax.legend(handles, custom_labels, loc='upper right')

# Título
plt.grid()
# Mostrar
plt.show()









import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import speclite.filters

# Leer los archivos .dat
data_J = np.loadtxt('Paranal_VISTA.J.dat')
data_H = np.loadtxt('Paranal_VISTA.H.dat')
data_Ks = np.loadtxt('Paranal_VISTA.Ks.dat')

# Crear los filtros con los datos leídos
# Asegurarse de que la respuesta empiece y termine en 0
wavelength_J = data_J[:, 0] * u.Angstrom
response_J = np.concatenate(([0], data_J[:, 1], [0]))  # Añadir 0 al principio y al final

wavelength_H = data_H[:, 0] * u.Angstrom
response_H = np.concatenate(([0], data_H[:, 1], [0]))  # Añadir 0 al principio y al final

wavelength_Ks = data_Ks[:, 0] * u.Angstrom
response_Ks = np.concatenate(([0], data_Ks[:, 1], [0]))  # Añadir 0 al principio y al final

# Crear los filtros con las longitudes de onda y respuestas
vhs_J = speclite.filters.FilterResponse(
    wavelength=wavelength_J,
    response=response_J,
    meta=dict(group_name='vhs', band_name='J')
)

vhs_H = speclite.filters.FilterResponse(
    wavelength=wavelength_H,
    response=response_H,
    meta=dict(group_name='vhs', band_name='H')
)

vhs_Ks = speclite.filters.FilterResponse(
    wavelength=wavelength_Ks,
    response=response_Ks,
    meta=dict(group_name='vhs', band_name='Ks')
)

# Cargar los filtros
vhs = [vhs_J, vhs_H, vhs_Ks]
colors = ['limegreen', 'cyan', 'blueviolet']

# Crear el gráfico
plt.figure(figsize=(10, 6))
for vhs_filter, color in zip(vhs, colors):
    plt.fill_between(vhs_filter.wavelength, vhs_filter.response, color=color, alpha=0.4, label=vhs_filter.meta['band_name'])
    plt.plot(vhs_filter.wavelength, vhs_filter.response, color=color)

# Configurar el gráfico
plt.ylim(0, 0.6)
plt.xlabel('Wavelength / Angstrom', fontsize=16)
plt.ylabel('Filter Response', fontsize=16)
plt.title('VISTA Filters', fontsize=16)

# Mostrar leyenda
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::1], ['J', 'H', 'Ks'])
plt.grid(False)

# Guardar el gráfico
plt.savefig('vistafilters.png')
plt.show()



#### VHS
data_J = np.loadtxt('Paranal_VISTA.J.dat') 
data_H = np.loadtxt('Paranal_VISTA.H.dat')
data_Ks = np.loadtxt('Paranal_VISTA.Ks.dat')
vhs_J=speclite.filters.FilterResponse(
    wavelength = data_J[:, 0] * u.Angstrom,
    response = data_J[:, 1], meta=dict(group_name='vhs', band_name='J'))
vhs_H=speclite.filters.FilterResponse(
    wavelength = data_H[:, 0] * u.Angstrom,
    response = data_H[:, 1], meta=dict(group_name='vhs', band_name='H'))
vhs_Ks=speclite.filters.FilterResponse(
    wavelength = data_Ks[:, 0] * u.Angstrom,
    response = data_Ks[:, 1], meta=dict(group_name='vhs', band_name='Ks'))
vhs=speclite.filters.load_filters('vhs_J','vhs_H','vhs_Ks')
colors = ['limegreen', 'cyan', 'blueviolet']
plt.figure(figsize=(10, 6))
for vhs_filter, color in zip(vhs, colors):
    # Convertir longitud de onda de Angstroms a micrones
    #wavelength_micron = wise_filter.wavelength / 10000.0
    plt.fill_between(vhs_filter.wavelength, vhs_filter.response, color=color, alpha=0.4, label=wise_filter.name)
    plt.plot(vhs_filter.wavelength, vhs_filter.response, color=color)
#plt.xscale('log')
#plt.gca().set_xticks([2.5, 5, 10, 20, 30])
#plt.gca().set_xticklabels([2.5, 5, 10, 20, 30])
#plt.xlim(2.5, 30)
plt.ylim(0,0.6)
plt.xlabel('Wavelength / Angstrom', fontsize=16)
plt.ylabel('Filter Response', fontsize=16)
plt.title('VISTA', fontsize=16)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::1], ['J','H','Ks'])
plt.grid(False)
plt.savefig('vistafilters.png')
plt.show()

"""
