import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress,pearsonr
from scipy.optimize import curve_fit



def convert_float(value):
    # Convert bytes to string if necessary
    if isinstance(value, bytes):
        value = value.decode('utf-8')
    return float(value.replace(',', '.'))

def load_data_from_file(filename):
    # Define converters for each column to handle commas as decimal separators
    converters = {0: convert_float, 1: convert_float, 2: convert_float, 3: convert_float}
    
    # Load the data from the file using the converters and specifying the tab delimiter
    data = np.loadtxt(filename, delimiter='\t', converters=converters)
    
    # Extract each column from the data
    x = data[:, 0]
    y = data[:, 1]
    errx = data[:, 2]
    erry = data[:, 3]
    
    return x, y, errx, erry

# Cambiar nombre del archivo de datos para cambiar los datos recogidos
x, y ,errx,erry= load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/datosPC_h=0.001.txt")
x1, y1,errx1,erry1 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/datosJOEL_-lm-O3_h=0.001.txt")
x2, y2,errx2,erry2 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/datosJOEL_-lm-O2_h=0.001.txt")

#slope1, intercept1 = linregress(x1, y1)
#slope2, intercept2 = linregress(x2, y2)

#def line_fit(x, slope, intercept):
  #  return slope * x + intercept

#Valores de x para la línea ajustada
line_x = np.array([min(x), max(x)])
line_x1 = np.array([min(x1), max(x1)])
line_x2 = np.array([min(x2), max(x2)])

# Valores de y usando la función de ajuste
line_y = np.array([min(y), max(y)])
line_y1 = np.array([min(y1), max(y1)])
line_y2 = np.array([min(y2), max(y2)])



#Cálculo del coef pearson
pearson_coef, _ = pearsonr(x, y)
pearson_coef1, _ = pearsonr(x1, y1)
pearson_coef2, _ = pearsonr(x2, y2)
print (pearson_coef)
print (pearson_coef1)
print (pearson_coef2)

num_bins = 50


plt.figure(figsize=(8, 6))

plt.errorbar(x, y, xerr=errx, yerr=erry, fmt='o', color='blue', capsize=5)  # fmt='o' para puntos
plt.errorbar(x1, y1, xerr=errx1, yerr=erry1, fmt='o', color='red', capsize=5)  # fmt='o' para puntos
plt.errorbar(x2, y2, xerr=errx2, yerr=erry2, fmt='o', color='yellow', capsize=5)  # fmt='o' para puntos

plt.plot(x, y, color='blue', label='-lm -O3 13th Gen Intel(R) Core(TM) i7-1360P 2.20 GHz  ', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x1, y1, color='red', label='-lm -O3 Joel', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x2, y2, color='yellow', label='-lm -O2 Joel', alpha=0.5)  # Dibuja los puntos de datos con transparencia

#plt.hist(y, bins=num_bins, alpha=0.5, label='Módulo de la velocidad (v)', color='blue', edgecolor='black')
#plt.hist(y1, bins=num_bins, alpha=0.5, label='Componente x de la velocidad (vx)', color='red', edgecolor='black')
#plt.hist(y2, bins=num_bins, alpha=0.5, label='Componente y de la velocidad (vy)', color='yellow', edgecolor='black')

#plt.plot(line_x, line_y, color='blue', label='Energía cinética')  # Dibuja la línea de ajuste
#plt.plot(line_x1, line_y1, color='red', label='Energía potencial') #para hacer mas de un ajuste en la misma grafica
#plt.plot(line_x2, line_y2, color='yellow', label='Energía total')  # Dibuja la línea de ajuste

# Configurar las etiquetas del eje x
plt.xscale('log')
plt.xticks([10, 100, 1000, 10000, 100000], ['10', '100', '1000', '10000', '100000'])

plt.autoscale(enable=True, axis='both', tight=False)


plt.xlabel('$t_f$ (s)')
plt.ylabel(' t (s)')
plt.title('Rendimiento de 13th Gen Intel(R) Core(TM) i7-1360P 2.20 GHz vs Joel, para distintas formas de optimización y para un $h=0.001$ ')

#plt.xlabel('$\phi$ (rad)')
#plt.ylabel('$\\psi$ (rad)')
#plt.title('Mapa de Poincaré para los parámetros $\phi$ y $\dot{\phi}$, y para distintas energías y condiciones iniciales $\phi=0.1rad$ y $\psi=0.2rad$')



plt.legend()
plt.grid(True)  #Poner cuadrilla
#plt.savefig("H:/Escritorio/UGR/Año 3/Electromagnetismo/Prácticas/Práctica 3/Fotos/3")  #Guardar plot en ruta
plt.show()

#se_intercept = std_err * np.sqrt(np.sum(x**2) / len(x))
#print(f"La ecuación de la línea es y = ({slope:.9f} ± {std_err:.9f})x + ({intercept:.9f} ± {se_intercept:.9f})")