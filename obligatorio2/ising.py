import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


#Asigno una T entre 0 y 5
T=3.0

#Establezco N como el tamaño de la matriz

N=10

#Creo una matriz NxN y la relleno de 1 y -1

s=np.random.choice([1,-1],(N,N))

#Pongo un rango de tiempo máximo para mi bucle e inicializo el contador de tiempo a 0
tf=10000
t=0

#Defino una variable de fichero donde guardar los datos
archivo=open("C:/Users/ignac/OneDrive/Escritorio/Fisica23-24/FisicaComputacional/Git/ignacio.matagomez-compu-2324/obligatorio2/isingresultados.txt","w")

while t<tf:

    for i in range (N*N):
    #Genero un punto al azar de la red, para ello, elijo un numero de fila al azar n y un numero de columna aleatorio m

        n=np.random.randint(0,N)
        m=np.random.randint(0,N)

        #Introducimos las condiciones de contorno 

        if (n==N-1): #Mi array va hasta N-1
            aux1=0
        else:
            aux1=n+1

        if (n==0): 
            aux2=n
        else:
            aux2=n-1


        if (m==N-1): 
            aux3=0
        else:
            aux3=m+1

        if (m==0): 
            aux4=m
        else:
            aux4=m-1


        E=2.0*s[n][m]*(s[aux1][m]+s[aux2][m]+s[n][aux3]+s[n][aux4])

        p=min(1.0,np.exp(-E/T))

        if (p>np.random.randint(0,1)):
            s[n][m]=-s[n][m]
    
    for i in range(N): #Columnas
        for j in range(N-1): #Filas
            archivo.write(f'{s[i][j]},')  
        archivo.write(f'{s[i][j]}')
        archivo.write('\n')  #Salto de línea
    archivo.write('\n')


    t=t+1
