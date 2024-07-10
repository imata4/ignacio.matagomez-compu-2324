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

# Usage example, replace the path with your actual data file path
x, y = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario1/sepcuadmed.txt")
#slope, intercept = linregress(x, y)
#slope1, intercept1 = linregress(x1, y1)
#slope2, intercept2 = linregress(x2, y2)

#def line_fit(x, slope, intercept):
  #  return slope * x + intercept

#Valores de x para la línea ajustada
line_x = np.array([min(x), max(x)])
# Valores de y usando la función de ajuste
line_y = np.array([min(y), max(y)])


#Cálculo del coef pearson
pearson_coef, _ = pearsonr(x, y)
print (pearson_coef)


plt.figure(figsize=(8, 6))

#plt.errorbar(x, y, xerr=errx, yerr=erry, fmt='o', color='blue', capsize=5)  # fmt='o' para puntos
#plt.errorbar(x1, y1, xerr=errx1, yerr=erry1, fmt='o', color='red', capsize=5)  # fmt='o' para puntos
#plt.errorbar(x2, y2, xerr=errx2, yerr=erry2, fmt='o', color='yellow', capsize=5)  # fmt='o' para puntos

plt.plot(x, y, color='blue')  # Dibuja los puntos de datos si pongo scatter, si pongo plot es una línea

#plt.plot(line_x, line_y, color='blue', label='Energía cinética')  # Dibuja la línea de ajuste
#plt.plot(line_x1, line_y1, color='red', label='Energía potencial') #para hacer mas de un ajuste en la misma grafica
#plt.plot(line_x2, line_y2, color='yellow', label='Energía total')  # Dibuja la línea de ajuste

plt.xlabel('Tiempo')
plt.ylabel('<$(\\Delta r_{ij}(t))^2$>')
plt.title('Separación cuadrática media')

plt.legend()
plt.grid(True)  #Poner cuadrilla
#plt.savefig("H:/Escritorio/UGR/Año 3/Electromagnetismo/Prácticas/Práctica 3/Fotos/3")  #Guardar plot en ruta
plt.show()

#se_intercept = std_err * np.sqrt(np.sum(x**2) / len(x))
#print(f"La ecuación de la línea es y = ({slope:.9f} ± {std_err:.9f})x + ({intercept:.9f} ± {se_intercept:.9f})")