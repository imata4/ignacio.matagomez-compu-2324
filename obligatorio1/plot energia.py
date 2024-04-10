from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
import numpy as np

# Listas para almacenar los datos de las dos variables
t = []
E = []

# Leer el archivo y extraer los datos
with open('C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/obligatorio1/ejemplo.txt', 'r') as file:
    for line in file:
        # Dividir la línea en dos valores
        values = line.split()
        # Validar que hay dos valores en la línea
        if len(values) == 2:
            try:
                # Convertir los valores a números y agregarlos a las listas
                t.append(float(values[0]))
                E.append(float(values[1]))
            except ValueError:
                print("Error: No se pueden convertir los valores a números flotantes en la línea:", line)
        else:
            print("Error: La línea no contiene dos valores:", line)

# Crear el gráfico de dispersión
plt.scatter(t, E)


# Agregar etiquetas a los ejes y título al gráfico
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (E)')
plt.title('Gráfico de Dispersión')

# Mostrar el gráfico
plt.show()