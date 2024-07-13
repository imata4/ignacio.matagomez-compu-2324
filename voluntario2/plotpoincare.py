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
    converters = {0: convert_float, 1: convert_float}#, 2: convert_float, 3: convert_float}
    
    # Load the data from the file using the converters and specifying the tab delimiter
    data = np.loadtxt(filename, delimiter=',', converters=converters)
    
    # Extract each column from the data
    x = data[:, 0]
    y = data[:, 1]
    #errx = data[:, 2]
    #erry = data[:, 3]
    
    return x, y#, errx, erry

# Cambiar nombre del archivo de datos para cambiar los datos recogidos
x, y = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/poincare_phi_psi_E=1.txt")
x1, y1 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/poincare_phi_psi_E=3.txt")
x2, y2 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/poincare_phi_psi_E=5.txt")
x3, y3 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/poincare_phi_psi_E=10.txt")
x4, y4 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario2/poincare_phi_psi_E=15.txt")
#slope, intercept = linregress(x, y)
#slope1, intercept1 = linregress(x1, y1)
#slope2, intercept2 = linregress(x2, y2)

#def line_fit(x, slope, intercept):
  #  return slope * x + intercept

#Valores de x para la línea ajustada
line_x = np.array([min(x), max(x)])
line_x1 = np.array([min(x1), max(x1)])
line_x2 = np.array([min(x2), max(x2)])
line_x3 = np.array([min(x3), max(x3)])
line_x4 = np.array([min(x4), max(x4)])
# Valores de y usando la función de ajuste
line_y = np.array([min(y), max(y)])
line_y1 = np.array([min(y1), max(y1)])
line_y2 = np.array([min(y2), max(y2)])
line_y3 = np.array([min(y3), max(y3)])
line_y4 = np.array([min(y4), max(y4)])



#Cálculo del coef pearson
pearson_coef, _ = pearsonr(x, y)
pearson_coef1, _ = pearsonr(x1, y1)
pearson_coef2, _ = pearsonr(x2, y2)
print (pearson_coef)
print (pearson_coef1)
print (pearson_coef2)

num_bins = 50

plt.figure(figsize=(8, 6))

#plt.errorbar(x, y, xerr=errx, yerr=erry, fmt='o', color='blue', capsize=5)  # fmt='o' para puntos
#plt.errorbar(x1, y1, xerr=errx1, yerr=erry1, fmt='o', color='red', capsize=5)  # fmt='o' para puntos
#plt.errorbar(x2, y2, xerr=errx2, yerr=erry2, fmt='o', color='yellow', capsize=5)  # fmt='o' para puntos

plt.plot(x, y, color='blue', label='E=1', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x1, y1, color='red', label='E=3', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x2, y2, color='purple', label='E=5', alpha=0.5)  # Dibuja los puntos de datos con transparencia
plt.plot(x3, y3, color='orange', label='E=10', alpha=0.5)
plt.plot(x4, y4, color='green', label='E=15', alpha=0.5)

#plt.hist(y, bins=num_bins, alpha=0.5, label='Módulo de la velocidad (v)', color='blue', edgecolor='black')
#plt.hist(y1, bins=num_bins, alpha=0.5, label='Componente x de la velocidad (vx)', color='red', edgecolor='black')
#plt.hist(y2, bins=num_bins, alpha=0.5, label='Componente y de la velocidad (vy)', color='yellow', edgecolor='black')

#plt.plot(line_x, line_y, color='blue', label='Energía cinética')  # Dibuja la línea de ajuste
#plt.plot(line_x1, line_y1, color='red', label='Energía potencial') #para hacer mas de un ajuste en la misma grafica
#plt.plot(line_x2, line_y2, color='yellow', label='Energía total')  # Dibuja la línea de ajuste

plt.xlabel('$\phi$ (rad)')
plt.ylabel('$\psi$ (rad)')
plt.title('Mapa de Poincaré para los parámetros $\phi$ y $\psi$, y para distintas energías y condiciones iniciales $\phi=0.1rad$ y $\psi=0.2rad$')

#plt.xlabel('$\phi$ (rad)')
#plt.ylabel('$\\psi$ (rad)')
#plt.title('Mapa de Poincaré para los parámetros $\phi$ y $\dot{\phi}$, y para distintas energías y condiciones iniciales $\phi=0.1rad$ y $\psi=0.2rad$')



plt.legend()
plt.grid(True)  #Poner cuadrilla
#plt.savefig("H:/Escritorio/UGR/Año 3/Electromagnetismo/Prácticas/Práctica 3/Fotos/3")  #Guardar plot en ruta
plt.show()

#se_intercept = std_err * np.sqrt(np.sum(x**2) / len(x))
#print(f"La ecuación de la línea es y = ({slope:.9f} ± {std_err:.9f})x + ({intercept:.9f} ± {se_intercept:.9f})")