import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.stats import linregress,pearsonr

m=1
k=1
N=20


def convert_float(value):
    # Convert bytes to string if necessary
    if isinstance(value, bytes):
        value = value.decode('utf-8')
    return float(value.replace(',', '.'))

def load_data_from_file(filename):
    # Define converters for each column to handle commas as decimal separators
    converters = {0: convert_float, 1: convert_float, 2: convert_float}

    data = np.loadtxt(filename, delimiter=',', converters=converters)
    
    # Extract each column from the data
    v = data[:, 0]
    vx = data[:, 1]
    vy = data[:, 2]
    
    
    return v, vx, vy


# Usage example, replace the path with your actual data file path
v, vx, vy = load_data_from_file("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario1/velocidades.txt")
print("v vale:")
print(v)
print("vx vale:")
print(vx)
print("vy vale:")
print(vy)



Temp=0
for i in range(100):
    Temp= Temp+0.5*(v[i]**2) 
Temp=Temp/100




def distribution_module (v):
    return (m/(k*Temp))*v*np.exp(-m*v*v/(2*k*Temp))

def distribution_velocity (v):
    return np.sqrt(2*m/(np.pi*k*Temp))*np.exp(-m*v*v/(2*k*Temp))

fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
v_uni = np.linspace(0, np.max(v), 100)
vx_uni = np.linspace(0, np.max(vx), 100)
vy_uni = np.linspace(0, np.max(vy), 100) 
plt.grid(zorder=0)
plt.hist(v, bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma de velocidades')
plt.plot(v_uni, distribution_module(v_uni), label=f'Distribución de velocidades ', zorder = 4)
plt.ylabel('nº de cuenta ')
plt.xlabel('velocidad')
plt.legend()
plt.title(f'Histograma del módulo de la velocidad, con temperatura T = {Temp:.2f}')

fig=plt.figure(figsize=(10,7)) 
ax=fig.add_subplot(111)
plt.grid(zorder=0)
plt.hist(v, bins=100 , rwidth=0.7, zorder = 3, density=True, label='Histograma para el modulo de la velocidad', alpha=0.7, color='tab:blue')
plt.hist(np.abs(vx), bins=100 , rwidth=0.7, zorder = 4, density=True, label='Histograma para el eje x', alpha=0.7, color='tab:purple')
plt.hist(np.abs(vy), bins=100 , rwidth=0.7, zorder = 5, density=True, label='Histograma para el eje y', alpha=0.7, color='tab:red')

plt.plot(v_uni, distribution_module(v_uni), label=f'Distribución del módulo de velocidades ', zorder = 7, color='tab:blue')
plt.plot(v_uni, distribution_velocity(vx_uni), label=f'Distribución de la componente x', zorder = 8, color='tab:purple')
plt.plot(v_uni, distribution_velocity(vy_uni), label=f'Distribución de la componente y', zorder = 9, color='tab:red')


plt.axvline(x=1, color='g', linestyle='--', label='|v| = 1')

# Definir los límites del eje y para asegurar que la línea vertical se muestre con altura 1
plt.ylim(0, 1)  # Ajusta el límite superior del eje y a 1 para que la línea vertical tenga la altura deseada
plt.ylabel('Número de cuentas')
plt.xlabel('Velocidad')
plt.legend()
plt.title(f'Histograma del valor absoluto de la velocidad, a temperatura T = {Temp:.2f}')

plt.show()