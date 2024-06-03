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

# Usage example, replace the path with your actual data file path
x, y, errx, erry = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/obligatorio1/tiempopc.txt")
x1, y1, errx1, erry1 = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/obligatorio1/tiempojoel.txt")
slope, intercept, r_value, p_value, std_err = linregress(x, y)
slope1, intercept1, r_value1, p_value1, std_err1 = linregress(x1, y1)

def line_fit(x, slope, intercept):
    return slope * x + intercept

#Valores de x para la línea ajustada
line_x = np.array([min(x), max(x)])
line_x1 = np.array([min(x1), max(x1)])
# Valores de y usando la función de ajuste
line_y = line_fit(line_x, slope, intercept)
line_y1 = line_fit(line_x1, slope1, intercept1)



#Cálculo del coef pearson
pearson_coef, _ = pearsonr(x, y)
pearson_coef1, _ = pearsonr(x1, y1)
print (pearson_coef)
print (pearson_coef1)

plt.figure(figsize=(8, 6))

plt.errorbar(x, y, xerr=errx, yerr=erry, fmt='o', color='blue', capsize=5)  # fmt='o' para puntos
plt.errorbar(x1, y1, xerr=errx1, yerr=erry1, fmt='o', color='red', capsize=5)  # fmt='o' para puntos

plt.scatter(x, y, color='blue', label='PC')  # Dibuja los puntos de datos
plt.scatter(x1, y1, color='red', label='JOEL')  # Dibuja los puntos de datos

plt.plot(line_x, line_y, color='blue', label='Ajuste PC')  # Dibuja la línea de ajuste
plt.plot(line_x1, line_y1, color='red', label='Ajuste JOEL') #para hacer mas de un ajuste en la misma grafica


plt.xlabel('Número de planetas')
plt.ylabel('t(s)')
plt.title('gcc planetasnuevos.c –o planetasnuevos.exe -lm –O3')

plt.legend()
plt.grid(True)  #Poner cuadrilla
#plt.savefig("H:/Escritorio/UGR/Año 3/Electromagnetismo/Prácticas/Práctica 3/Fotos/3")  #Guardar plot en ruta
plt.show()

se_intercept = std_err * np.sqrt(np.sum(x**2) / len(x))
print(f"La ecuación de la línea es y = ({slope:.9f} ± {std_err:.9f})x + ({intercept:.9f} ± {se_intercept:.9f})")