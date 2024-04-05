#include <stdio.h>
#include <math.h>

//Declaración de funciones a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

//void reescalamiento(double *r,double *m );
void leerArchivo(const char *nombre_archivo, double vector[], int *num_elementos);

//Se trata de una función a la que le pasaremos el nombre del archivo, un vector y el número de elementos de dicho vector. 
//Su objetivo es leer datos de un archivo y almacenarlos en un vector.
void leerArchivo(const char *nombre_archivo, double vector[], int *num_elementos)
{
    FILE *archivo;
    int num_leidos=0; //Indica el número de datos que se han leído hasta el momento
    double valor; //

    archivo=fopen(nombre_archivo,"r");
    

     // Verificamos si el archivo se abrió correctamente
    if (archivo == NULL) {
        printf("No se pudo abrir el archivo.\n");
    }
   
    num_leidos=0;
    while (fscanf(archivo,"%lf",&valor)!=EOF) 
    {
        vector[num_leidos]=valor; //Asingamos el valor de cada dato que va leyendo, al vector en la posición equivalente al número de datos leídos del vector
        printf("%lf",valor); //Vemos en pantalla que se ha leído bien el archivo y guardado en el vector.
        printf("\n");
        num_leidos=num_leidos+1; //Incrementamos el número de datos leídos en uno para actualizar el siguiente dato que ha de leerse.
    } 

    fclose(archivo);

    *num_elementos=num_leidos; //Hago que el puntero de num_elementos valga num_leidos, es decir que me devuelva el número de elementos que ha leido

                                //Esto quiere decir, que al introducir el número de elementos a leer, la función los recogerá de la variable num_leidos


    return;

}

//Hago una función que permita el reescalamiento 
//void reescalamiento(double r[],double m[], int size) //Necesito pasarle las posiciones de los planetas y sus masas



void main(void)
{
    double h; 
    int i,j; //Variables auxiliares para bucles
    int size;
    double t; //Tiempo

    double m[5]; //masa

    double r_x[5]; //Vector para obtener la componente x del vector de posición en las condiciones iniciales
    double r[5][5]; //Matriz de posiciones, primera columna tiene componente x, y la segunda la componente y

    double v[5][5]; //Matriz de velocidades, primera columna tiene componente x, y la segunda la componente y
    double v_y[5]; // Vector para obtener la componente y del vector velocidad en las condiciones iniciales



    int num_elementos=0; //Me sirve para leer los datos de los archivos y obtener cuantos elementos hay





//Consideraremos que el radio orbital es la posición inicial de todos los planetas
//Entonces todos estarán en la posición x=radio orbital y=0. Con velocidad v_y=velocidad orbital v_x=0
//La posición tendrá componentes x e y, por lo tanto tendré un vector de dos columnas (matriz) donde la primera es la posición en x y la segunda la posición en y

    leerArchivo("masas.txt",m,&num_elementos);

    printf("Se han leído %d elementos del archivo:\n", num_elementos); //%d hace referencia a num_elementos
  
        for (i = 0; i <num_elementos; i++) 
        {
            printf("%lf ", m[i]);
        }
        printf("\n");

//En este fichero se encuentran en la primera fila todos los valores del radio orbital en m
//Leo y actualizo el vector de posiciones con los datos iniciales, en este caso solo tendrán componente x, así que trabajo con el vector x
leerArchivo("posiciones_iniciales.txt",r_x,&num_elementos);

printf("Se han leído %d elementos del archivo:\n", num_elementos); //%d hace referencia a num_elementos
  
        for (i = 0; i < num_elementos; i++) 
        {
            printf("%lf ", r_x[i]);
        }
        printf("\n");
   
   
//En este fichero se encuentran en la primera fila todos los valores de la velocidad en m/s
leerArchivo("velocidades_iniciales.txt",v_y,&num_elementos);

printf("Se han leído %d elementos del archivo:\n", num_elementos); //%d hace referencia a num_elementos
  
        for (i = 0; i < num_elementos; i++) 
        {
            printf("%lf ", v_y[i]);
        }
        printf("\n");

//Voy a trabajar con una matriz de posiciones y otra de velocidades, no obstante me he valido de vectores para obtener los datos iniciales.
//Tanto en posición como en velocidad, he supuesto una componente nula, ahora reescribo mi matriz a partir de los vectores obtenidos.
//Hago esto despues de haber machacado la variable num_elementos, pues siempre voy a tener a ser igual en posiciones y en velocidades, e igual al número de cuerpos


for(i=0;i<num_elementos;i++)
   {
        r[i][0]=r_x[i];
        v[i][0]=v_y[i];
   }

//Compruebo que este bien la matriz de posiciones (Nota: Al no haber asignado valor a la componente y, adquiere el valor 0 por defecto, lo que queríamos)

 for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j <2; j++) {
            printf("%lf\t", r[i][j]); // Imprimir cada elemento de la matriz
        }
        printf("\n"); // Salto de línea al final de cada fila
    }

//Compruebo que también haya quedado bien la matriz de velocidades

for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j <2; j++) {
            printf("%lf\t", v[i][j]); // Imprimir cada elemento de la matriz
        }
        printf("\n"); // Salto de línea al final de cada fila
    }
    


}