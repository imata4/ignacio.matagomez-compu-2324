#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

void aleatorios(int s[][3],int N);


void main(void)
{
     int i,j;
     int N; //Tamaño de la red cuadrada
     int s[N][N]; //Spin--Valores 1 y -1
     int T;

     //Asignamos una T entre 0 y 5
     T=1;

     //Le doy un valor a N
     N=3;

    //Asignamos a cada elemento de la matriz de spins, (que será un spin) un valor aleatorio 1 o -1 con una probabilidad de 1/2

    
    // Inicializar la semilla usando el tiempo actual
    srand(time(NULL));
    
   // Generar un número aleatorio entre 0 y 2*RAND_MAX, es decir entre 0 y 2
  int aleatorio = rand() % 2 == 0 ? 1 : -1; // Generar aleatoriamente -1 o 1

    // Restar 1 para obtener un número entre -1 y 1
    aleatorio -= 1;



    aleatorios(s,N);

    return;

}

void aleatorios(int s[][3],int N)
{
    int i,j;
    double probabilidad;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            probabilidad = (double)rand() / RAND_MAX; // Generar número aleatorio entre 0 y 1
            s[i][j] = (probabilidad <= 0.5) ? 1 : -1; // Asignar 1 o -1 dependiendo de la probabilidad
        }
    }

    return;
}