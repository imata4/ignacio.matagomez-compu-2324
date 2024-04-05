#include <stdio.h>
#include <math.h>

//Declaración de funciones a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

//void reescalamiento(double *r,double *m );

void main(void)
{
    int i;
    int size;
    double t; //Tiempo
    double x; 
    double h; 
    double m[5]; //masa
    double r[5]; //posiciones
    double v[5]; // velocidades 
    int flag;

    FILE *fe, *fs; // Punteros para ficheros de entrada y salida


    fe=fopen("masas.txt","r");
    

     // Verificamos si el archivo se abrió correctamente
    if (fe == NULL) {
        printf("No se pudo abrir el archivo.\n");
        return 1;
    }


   
    i=0;
    while (fscanf(fe,"%lf",&x)!=EOF)
    {
        m[i]=x;
        printf("%lf",x);
        printf("\n");
        i=i+1;
    } 

    fclose(fe);

    for (i=0;i<5;i++)
    {
         printf("%.3e\n", m[i]);
    }

    //fe=fopen("posiciones_planetas.txt","r");


//Hago una función que permita el reescalamiento 
//void reescalamiento(double r[],double m[], int size) //Necesito pasarle las posiciones de los planetas y sus masas

    return ;
}