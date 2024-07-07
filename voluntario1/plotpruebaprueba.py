import numpy as np
import matplotlib.pyplot as plt

# Función para leer los datos del archivo
def leer_datos(archivo):
    datos = np.loadtxt(archivo, delimiter=',')
    x = datos[:, 0]  # Primera columna
    y1 = datos[:, 1]  # Segunda columna
    y2 = datos[:, 2]  # Tercera columna
    y3 = datos[:, 3]  # Cuarta columna
    return x, y1, y2, y3

# Función para graficar los datos
def graficar_datos(x, y1, y2, y3):
    plt.figure(figsize=(10, 6))
    plt.hist(y1, bins=100, color='b', alpha=0.5, label='Módulo de Velocidad')
    plt.hist(y2, bins=100, color='g', alpha=0.5, label='vx')
    plt.hist(y3, bins=100, color='r', alpha=0.5, label='vy')
    plt.xlabel('Valores de X')
    plt.ylabel('Valores de Y')
    plt.title('Gráfico de Dispersión de Y1, Y2 y Y3 en función de X')
    plt.legend()
    plt.show()
# Ruta al archivo de datos
archivo_datos = "C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/voluntario1/v1.txt"

# Leer los datos
x, y1, y2, y3 = leer_datos(archivo_datos)

# Graficar los datos
graficar_datos(x, y1, y2, y3)