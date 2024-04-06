#include <stdio.h>
#include <math.h>

//Declaración de funciones a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

void reescalamiento(double r[][2],double m[],double t, int size);
void leerArchivo(const char *nombre_archivo, double vector[], int *num_elementos);
void modulo_distancia(double r[][2], double mod_dist[][5], int size);
void aceleracion(double r[][2],double m[],double a[][2],int size);
void resultados(const char *nombre_archivo, double r[][2], int num_elementos);

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
void reescalamiento(double r[][2],double m[],double t, int size) //Necesito pasarle las posiciones de los planetas y sus masas
{
    int i,j;

    for(i=0;i<size;i++)
    {
        m[i]=m[i]/Ms;   //Reescalamiento de la masa
        for(j=0;j<2;j++)
        {
            r[i][j]=r[i][j]/c;  //Reescalamiento de la posición
        }
    }

    t=t*sqrt(G*Ms/(pow(c,3)));  //Reescalamiento del tiempo

    return;
}

void modulo_distancia(double r[][2], double mod_dist[][5], int size)
{
    int i,j;

    //La matriz mod_dist guardará el módulo de la distancia entre la posición de la masa i-ésima y la de la masa j-ésima, cumpliéndose que i!=j.
    
    for ( i = 0; i < size; i++ )
    {
        
         for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                mod_dist[i][j]=sqrt(pow(r[i][0]-r[j][0],2)+pow(r[i][1]-r[j][1],2));
                  printf("%lf",mod_dist[i][j]);
            }
            
         }
          printf("\n");
    

            
    }
    
    return;
}


void aceleracion(double r[][2],double m[],double a[][2],int size)
{
    int i,j;
    double mod_dist[5][5];

    //Calculo los módulos de las distancias que se me guardarán en la matriz mod_dist, que usaré para calcular la aceleración
    modulo_distancia(r,mod_dist,size);

    for ( i = 0; i < size; i++)
    {
        for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                a[i][0]-=(m[j]*(r[i][0]-r[j][0]))/pow(mod_dist[i][j],3);
                a[i][1]-=(m[j]*(r[i][1]-r[j][1]))/pow(mod_dist[i][j],3);
            }
            
        }
        
    }
    return;
}


void resultados(const char *nombre_archivo, double r[][2], int num_elementos)
{
    FILE *archivo = fopen(nombre_archivo, "w"); // Abre el archivo en modo escritura (sobreescribe si ya existe)

    if (archivo == NULL) {
        printf("No se pudo abrir el archivo.\n");
        return;
    }

    // Escribe la matriz en el archivo
    for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j < 2; j++) {
            fprintf(archivo, "%lf ", r[i][j]); // Escribe cada elemento de la matriz
        }
        fprintf(archivo, "\n"); // Nueva línea al final de cada fila
    }

    fclose(archivo); // Cierra el archivo

    return;
}



void main(void)
{
    double h; 
    int i,j; //Variables auxiliares para bucles
    int size;
    double t; //Tiempo

    double m[5]; //masa

    double r_x[5]; //Vector para obtener la componente x del vector de posición en las condiciones iniciales
    double r[5][2]; //Matriz de posiciones, primera columna tiene componente x, y la segunda la componente y

    double v[5][2]; //Matriz de velocidades, primera columna tiene componente x, y la segunda la componente y
    double v_y[5]; // Vector para obtener la componente y del vector velocidad en las condiciones iniciales

    double mod_dist[5][5]; //Matriz que almacena el modulo de la distancia

    double a[5][2]; //Matriz que guarda las aceleraciones de los 5 planetas, componente x en columna 1 y a_y en columna 2

    double w[5][2];

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
    

//Hago el reescalamiento
//Tengo que meter un tiempo distinto de 0 porque sino, el reescalamiento no va a funcionar

printf("Reescalamiento\n");
t=0;

reescalamiento(r,m,t,num_elementos);

for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j <2; j++) {
            printf("%e\t", r[i][j]); // Imprimir cada elemento de la matriz
        }
        printf("\n"); // Salto de línea al final de cada fila
    }

    for (int i = 0; i < num_elementos; i++) {
            printf("%e\t", m[i]); // Imprimir cada elemento de la matriz
        printf("\n"); // Salto de línea al final de cada fila
    }

    printf("%e",t); //Revisar luego el reescalamiento temporal

    printf("\n");

    /*##############################
    ########### PASO 1 #############
    ##############################*/
    
    //Calculamos el módulo de la distancia entre dos cuerpos y lo guardo en una matriz
    //Para ello utilizo una función que me calcula el modulo de la distancia entre los 2 cuerpos.

 modulo_distancia(r,mod_dist,num_elementos); //La matriz mod dist

//Compruebo que quedan bien almacenados los datos en la matriz mod_dist, la matriz tiene 0 cuando no entra al bucle for, es decir cuando i!=j, pero solo imprimo los valores no nulos
 printf("\n");

   for ( i = 0; i < num_elementos; i++ )
    {
       
            for ( j = 0; j < num_elementos; j++)
            {   if (i!=j)
                printf("%lf",mod_dist[i][j]);
            }
            printf("\n");
       
            
    }
    
    printf("\n");
    //Ya puedo calcular la aceleración
    aceleracion(r,m,a,num_elementos);

    //Veo como queda la matriz de aceleraciones

   for ( i = 0; i < num_elementos; i++)
   {
        for ( j = 0; j < 2; j++)
            {
                printf("%lf",a[i][j]);
            }
    printf("\n");
   }
    

      /*##############################
    ########### PASO 2 #############
    ##############################*/


     //Evalúo las posiciones a tiempo t y también w
    
    for ( i = 0; i < num_elementos; i++)
   {
        for ( j = 0; j < 2; j++)
            {
               r[i][j]=r[i][j]+h*v[i][j]+0.5*h*h*a[i][j]; //Actualizo la posición del vector de posiciones
               w[i][j]=v[i][j]+0.5*h*a[i][j];
            }
    
   }


   /*##############################
    ########### PASO 3 #############
    ##############################*/

    //Evalúo las aceleraciones con mis nuevas posiciones 

    aceleracion(r,m,a,num_elementos); //La única diferencia es que ahora mi r está modificado, así que obtendré otra aceleración

    /*##############################
    ########### PASO 4 #############
    ##############################*/

    //Evalúo las velocidades con mis nuevas aceleraciones

      for ( i = 0; i < num_elementos; i++)
   {
        for ( j = 0; j < 2; j++)
            {
               v[i][j]=v[i][j]+w[i][j]+0.5*h*a[i][j];
            }
    
   }
t=t+h;
   //METO LOS RESULTADOS EN UN ARCHIVO

   resultados("resultados.txt",r, num_elementos);




//PROGRAMA ORIGINAL



/*
    //Evalúo las posiciones a tiempo t y también w
    
    for ( i = 0; i < num_elementos; i++)
   {
        for ( j = 0; j < 2; j++)
            {
               r[i][j]=r[i][j]+h*v[i][j]+0.5*h*h*a[i][j]; //Actualizo la posición del vector de posiciones
               w[i][j]=v[i][j]+0.5*h*a[i][j];
            }
    
   }
*/

   /*##############################
    ########### PASO 3 #############
    ##############################*/

    //Evalúo las aceleraciones con mis nuevas posiciones 

   // aceleracion(r,m,a,num_elementos); //La única diferencia es que ahora mi r está modificado, así que obtendré otra aceleración

    /*##############################
    ########### PASO 4 #############
    ##############################*/

    //Evalúo las velocidades con mis nuevas aceleraciones
/*
      for ( i = 0; i < num_elementos; i++)
   {
        for ( j = 0; j < 2; j++)
            {
               
               v[i][j]=v[i][j]+w[i][j]+0.5*h*a[i][j];
            }
    
   }
*/
     /*##############################
    ########### PASO 5 #############
    ##############################*/

    //Elevo un paso en el tiempo
   // t=t+h;




}