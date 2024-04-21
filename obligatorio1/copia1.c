#include <stdio.h>
#include <math.h>

//Declaración de funciones a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

void reescalamiento(double r[][2],double m[],double v[][2],double t, int size);
void leerArchivo(const char *nombre_archivo, double vector[], int *num_elementos);
void modulo_distancia(double r[][2], double mod_dist[][6], int size);
void aceleracion(double r[][2],double m[],double a[][2],int size);
void escribeMatriz(const char *nombre_archivo, double r[][2], int num_elementos);
void escribeVector(const char *nombre_archivo, double v[], int num_elementos);
double calculoT(double v[][2],double m[], int size);
double calculoV(double r[][2],double m[],int size);
double calculoE(double V,double T);
double calculoL(double r[][2],double m[],double v[][2],int size);
void calculoPeriodo(double r[][2],double aux[],double tiempo,double T[],int size);




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
        //printf("%lf",valor); //Vemos en pantalla que se ha leído bien el archivo y guardado en el vector.
       // printf("\n");
        num_leidos=num_leidos+1; //Incrementamos el número de datos leídos en uno para actualizar el siguiente dato que ha de leerse.
    } 

    fclose(archivo);

    *num_elementos=num_leidos; //Hago que el puntero de num_elementos valga num_leidos, es decir que me devuelva el número de elementos que ha leido

                                //Esto quiere decir, que al introducir el número de elementos a leer, la función los recogerá de la variable num_leidos


    return;

}

//Hago una función que permita el reescalamiento 
void reescalamiento(double r[][2],double m[],double v[][2],double t, int size) //Necesito pasarle las posiciones de los planetas y sus masas
{
    int i,j;

       t=t*sqrt((G*Ms)/(pow(c,3)));  //Reescalamiento del tiempo

    for(i=0;i<size;i++)
    {
        m[i]=m[i]/Ms;   //Reescalamiento de la masa
        for(j=0;j<2;j++)
        {
            r[i][j]=r[i][j]/c;//Reescalamiento de la posición
            v[i][j]=v[i][j]*sqrt(c/(G*Ms));
        }
    }

    

    return;
}

void modulo_distancia(double r[][2], double mod_dist[][6], int size)
{
    int i,j;

    //La matriz mod_dist guardará el módulo de la distancia entre la posición de la masa i-ésima y la de la masa j-ésima, cumpliéndose que i!=j.
    
    for ( i = 0; i < size; i++ )
    {
        
         for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                mod_dist[i][j]=sqrt(pow((r[i][0]-r[j][0]),2)+pow((r[i][1]-r[j][1]),2));
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
    double mod_dist[size][size];
    double mod;

    //Calculo los módulos de las distancias que se me guardarán en la matriz mod_dist, que usaré para calcular la aceleración
    modulo_distancia(r,mod_dist,size);


    for ( i = 0; i < size; i++)
    {
        a[i][0]=0;
        a[i][1]=0;
        for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                a[i][0]-=(m[j]*(r[i][0]-r[j][0]))/(pow(mod_dist[i][j],3));
                a[i][1]-=(m[j]*(r[i][1]-r[j][1]))/(pow(mod_dist[i][j],3));
               //mod=sqrt(pow(r[i][0]-r[j][0],2)+pow(r[i][1]-r[j][1],2));

              // a[i][0]-=(m[j]*(r[i][0]-r[j][0]))/(pow(mod,3));
               //a[i][1]-=(m[j]*(r[i][1]-r[j][1]))/(pow(mod,3));
            }
            
        }
        
    }
    return;
}


void escribeMatriz(const char *nombre_archivo, double r[][2], int num_elementos)
{
    FILE *archivo = fopen(nombre_archivo, "a"); // Abre el archivo en modo añadir,de forma que tras cada iteración añada, nuevos datos al archivo
                                                //Si lo abriese en mode w "write", lo que ocurriría es que los datos se sobreescribirián, y solo quedarían los de la última iteración.

    if (archivo == NULL) {
        printf("No se pudo abrir el archivo.\n");
        return;
    }

    // Escribe la matriz en el archivo
    for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j < 2; j++)
        {
            fprintf(archivo, "%lf ", r[i][j]); // Escribe cada elemento de la matriz
            if (j==0)
            {
                fprintf(archivo,",");
            }
            
            
        }
        fprintf(archivo, "\n"); // Nueva línea al final de cada fila
    }

    fprintf(archivo, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones

    fclose(archivo); // Cierra el archivo

    return;
}


void escribeVector(const char *nombre_archivo, double v[], int num_elementos)

{
     FILE *archivo = fopen(nombre_archivo, "a"); 
                                              
    if (archivo == NULL) {
        printf("No se pudo abrir el archivo.\n");
        return;
    }

    // Escribe la matriz en el archivo
    for (int i = 0; i < num_elementos; i++) 
    {
        fprintf(archivo, "%lf ", v[i]); // Escribe cada elemento del vector

        //fprintf(archivo, "\n"); // Nueva línea al final de cada fila
    }

    fprintf(archivo, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones

    fclose(archivo); // Cierra el archivo

    return;
}

double calculoT(double v[][2],double m[],int size)
{
    int i,j;
    double T=0;
    

    for ( i = 0; i < size; i++)
    {
            T=0.5*m[i]*sqrt(pow(v[i][0],2)+pow(v[i][1],2)); //El módulo lo calculo como la raíz cuadrada de la suma de las componentes al cuadrado
        
    }

    return T;
}


double calculoV(double r[][2],double m[],int size)
{
     int i,j;
    double V=0;

    double mod_dist[size][size];
     //Calculo los módulos de las distancias que se me guardarán en la matriz mod_dist, que usaré para calcular la aceleración
    modulo_distancia(r,mod_dist,size);

    for ( i = 0; i < size; i++)
    {
        for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                V-=m[i]*m[j]/mod_dist[i][j];
            }
            
        }
        return V;
    }
  
}



double calculoE(double V,double T)
{   
    double E=V+T;

    return E;
}

double calculoL(double r[][2],double m[],double v[][2],int size)
{
    int i,j;

    double L=0; //Inicializo el momento

    for ( i = 0; i < size; i++)
    {
           L+=m[i]*(r[i][0]*v[i][1]-r[i][1]*v[i][0]); //Calculo el momento a partir del producto vectorial, será un vector con dirección +z
    }
    

    return L;
}



void calculoPeriodo(double r[][2],double aux[],double tiempo,double T[],int size)
{
    int i;

    for ( i = 0; i < size; i++)
    {
        if (aux[i]<0 && r[i][1]>0 && T[i]==0)
        {
            T[i]=tiempo;
        }
        
    }
    
    return;

}

void main(void)
{
    double h=0.006;  //Inicializo h a 0.1 que será el paso temporal cda vez que se ejecute el algoritmo
    int i,j; //Variables auxiliares para bucles
    int size;
    double t; //Tiempo

    int n=6;    //Número de planetas

    double m[n]; //masa

    double r_x[n]; //Vector para obtener la componente x del vector de posición en las condiciones iniciales
    double r[n][2]; //Matriz de posiciones, primera columna tiene componente x, y la segunda la componente y

    double v[n][2]; //Matriz de velocidades, primera columna tiene componente x, y la segunda la componente y
    double v_y[n]; // Vector para obtener la componente y del vector velocidad en las condiciones iniciales

    double mod_dist[n][n]; //Matriz que almacena el modulo de la distancia

    double a[n][2]; //Matriz que guarda las aceleraciones de los 5 planetas, componente x en columna 1 y a_y en columna 2

    double w[n][2];

    int num_elementos=0; //Me sirve para leer los datos de los archivos y obtener cuantos elementos hay
    double tf;
    double E,T,V=0; //Variables para energía mecánica, potencial y cinética, respectivamente, e inicializadas a 0.

    double L=0; //Momento inicializado a 0.

    double periodo[n];

    double aux[n];

    FILE *f1,*f2;
    f1= fopen("energia.txt", "a"); 

    f2= fopen("momento.txt", "a"); 
   



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
        v[i][1]=v_y[i];
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
tf=1*365*24*3600;

reescalamiento(r,m,v,tf,num_elementos);

for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j <2; j++) {
            printf("%lf\t", r[i][j]); // Imprimir cada elemento de la matriz
        }
        printf("\n"); // Salto de línea al final de cada fila
    }

    for (int i = 0; i < num_elementos; i++) {
            printf("%e\t", m[i]); // Imprimir cada elemento de la matriz
        printf("\n"); // Salto de línea al final de cada fila
    }

    printf("%lf",tf); //Revisar luego el reescalamiento temporal

    for (int i = 0; i < num_elementos; i++) {
        for (int j = 0; j <2; j++) {
            printf("%lf\t", v[i][j]); // Imprimir cada elemento de la matriz
        }
        printf("\n"); // Salto de línea al final de cada fila
    }

    printf("\n");

    /*##############################
    ########### PASO 1 #############
    ##############################*/
    
    //Calculamos el módulo de la distancia entre dos cuerpos y lo guardo en una matriz
    //Para ello utilizo una función que me calcula el modulo de la distancia entre los 2 cuerpos.

 modulo_distancia(r,mod_dist,num_elementos); //La matriz mod dist

//Compruebo que quedan bien almacenados los datos en la matriz mod_dist, la matriz tiene 0 cuando no entra al bucle for, es decir cuando i!=j, pero solo imprimo los valores no nulos
 printf("\n");
printf("Módulo de la distancia");
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
    
for ( t = 0; t < tf; t+=h)
{
   //Almaceno las posiciones en la variable auxiliar. Útil para sacar el período.
   for ( i = 0; i < num_elementos; i++)
   {
        aux[i]=r[i][1];     
   }
   


    //Calculo el momento
    L=calculoL(r,m,v,num_elementos);

    //Calculo la energía potencial
    V=calculoV(r,m,num_elementos);

     //Calculo la energía cinética.
    T=calculoT(v,m,num_elementos);

     //Calculo la energía mecánica.
    E=calculoE(V,T);


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

   calculoPeriodo(r,aux,t,periodo,num_elementos); //Calculo los períodos si es que ya han recorrido un periodo

   

   //METO LOS RESULTADOS EN UN ARCHIVO
 
   escribeMatriz("resultados.txt",r, num_elementos);

   // Escribo la energía mecánica en energia.txt
    
            fprintf(f1, "%lf, %lf ", t, E); 
   
        fprintf(f1, "\n"); // Nueva línea al final de cada fila
 
// Escribo el momento angular en momento.txt
    fprintf(f2, "%lf, %lf ", t, L); 
   
        fprintf(f2, "\n"); // Nueva línea al final de cada fila

    

   
}

fclose(f1);//Cierro el archivo donde guardo las energias
fclose(f2);//Cierro el archivo donde guardo los momentos
//Imprimo en un archivo los periodos
escribeVector("periodos.txt",periodo,num_elementos);



}