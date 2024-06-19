#include <stdio.h>
#include <math.h>
#include <time.h>

//Declaración de constantes a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

void leerArchivo(const char *nombre_archivo, double vector[]);
void reescalamiento(double r[][2],double m[],double v[][2], int size);
void aceleracion(double r[][2],double m[],double a[][2],int size);
void geocentrico(double r[][2],double rgeo[][2],int size);
double calculoT(double v[][2],double m[],int size);
double calculoV(double r[][2],double m[],int size);
double calculoE(double V,double T);
double calculoL(double r[][2],double m[],double v[][2],int size);
void calculoPeriodo(double r[][2],double aux[],double tiempo,double T[],int size);

void main(void)
{ 
     //Pongo esto para ver el tiempo de compilación
    clock_t begin= clock();

    int i,j; //Variables auxiliares para bucles
    double h; //Paso temporal
    h=0.01;

    int n; //Número de planetas
    n=10; 

    double m[n];    //Vector de masas

    //Defino un vector de posicion rx y otro vy para leer las condiciones iniciales
    double rx[n];
    double vy[n];

    double r[n][2]; //Matriz de posiciones
    double v[n][2]; //Matriz de velocidades
    double a[n][2]; //Matriz de aceleraciones
    double w[n][2]; //Matriz auxiliar para calcular velocidad

    double t,tf; //Tiempo para el bucle temporal t, y tf tiempo final del bucle temporal

    double E,T,V=0; //Variables para energía mecánica, potencial y cinética, respectivamente, e inicializadas a 0.

    double L=0; //Momento inicializado a 0.

    double periodo[n];

    double aux[n];

    double rgeo[n][2];


    FILE *f1,*f2,*f3,*f4,*f5;
    f1= fopen("resultados.txt", "w"); 
    f2= fopen("geocentrico.txt", "w");
    f3= fopen("energia.txt", "w"); 
    f4= fopen("momento.txt", "w");
    f5= fopen("periodos.txt", "w");  
    

     //Inicializo el vector de periodos
    for ( i = 0; i < n; i++)
        periodo[i]=0;


    //Leo las masas, la posición en x y la velocidad en y de los archivos
    leerArchivo("masas.txt",m);
    leerArchivo("posiciones_iniciales.txt",rx);
    leerArchivo("velocidades_iniciales.txt",vy);
   

    //Completo la matriz de posiciones y velocidades, suponiendo que inicialmente los planetas están situados en el eje x y con velocidad en y

    for ( i = 0; i < n; i++)
    {
        r[i][0]=rx[i];
        r[i][1]=0;
        v[i][0]=0;
        v[i][1]=vy[i];
    }

    //Hago el reescalamiento 

    reescalamiento(r,m,v,n);

    //Establezco el tiempo final de mi simulación en años

    tf=300;
    //Paso el tiempo a segundos y lo reescalo
    tf=tf*365*24*3600;
    tf=tf*sqrt((G*Ms)/(pow(c,3)));
    
    //Calculo las aceleraciones
    aceleracion(r,m,a,n);

    

    for ( t = 0; t < tf; t+=h)
    {   

        for ( i = 0; i < n; i++)
            aux[i]=r[i][1];     
   

         //Calculo el momento
        L=calculoL(r,m,v,n);

        //Calculo la energía potencial
        V=calculoV(r,m,n);

        //Calculo la energía cinética.
        T=calculoT(v,m,n);

        //Calculo la energía mecánica.
        E=calculoE(V,T);

        //Utilizo esta función para calcular la posición geocéntrica
        geocentrico(r,rgeo,n); 

        //METO LOS RESULTADOS GEOCÉNTRICOS Y HELIOCÉNTRICOS ARCHIVOS

        // Escribe la matriz en el archivo
        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < 2; j++)
            {
                fprintf(f1, "%lf ", r[i][j]); // Escribe cada elemento de la matriz
                fprintf(f2, "%lf ", rgeo[i][j]); 
                if (j==0)
                {
                    fprintf(f1,",");
                    fprintf(f2,",");
                }
            
            
            }
            fprintf(f1, "\n"); // Nueva línea al final de cada fila
            fprintf(f2, "\n");
        }
        fprintf(f1, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones
        fprintf(f2, "\n"); 

        //Escribo la energía y el momento para cada instante tiempo en archivos

         // Escribo la energía mecánica en energia.txt
    
        fprintf(f3, "%lf, %lf ", t, E); 
   
        fprintf(f3, "\n"); // Nueva línea al final de cada fila
 
        // Escribo el momento angular en momento.txt
        fprintf(f4, "%lf, %lf ", t, L); 
   
        fprintf(f4, "\n"); // Nueva línea al final de cada fila


        //Actualizo las posiciones

        for ( i = 0; i < n; i++)
            for ( j = 0; j < 2; j++)
                {
                    r[i][j]=r[i][j]+h*v[i][j]+(h*h*a[i][j])/2.0; //Actualizo la posición del vector de posiciones
                    w[i][j]=v[i][j]+(h*a[i][j])/2.0;
                }
    

        //Evalúo las aceleraciones con mis nuevas posiciones 

        aceleracion(r,m,a,n); //La única diferencia es que ahora mi r está modificado, así que obtendré otra aceleración

        //Evalúo las velocidades con mis nuevas aceleraciones

        for ( i = 0; i < n; i++)
            for ( j = 0; j < 2; j++)
                {
                    v[i][j]= w[i][j]+(h*a[i][j])/2.0;
                }

        calculoPeriodo(r,aux,t,periodo,n); //Calculo los períodos si es que ya han recorrido un periodo

        
    }

    //Cambio las unidades del periodo de segundos reescalados a años
    for ( i = 0; i < n; i++)
    {
        periodo[i]=periodo[i]/(sqrt((G*Ms)/(pow(c,3)))*(3600*24*365));
    }
    
    //Imprimo los periodos en un archivo
    for ( i = 0; i < n; i++)
        {
            if(periodo[i]!=0)
            {
                fprintf(f5, "%lf", periodo[i]);
                fprintf(f5, "\n");              
            }
                
        }
    //Cierro todos los archivos
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    //Calculo el tiempo de ejecución
    clock_t end=clock();
    double tiempo= (double) (end-begin)/ CLOCKS_PER_SEC;
    printf("Tiempo de compilacion=%lf \n", tiempo);

    return;
}


void leerArchivo(const char *nombre_archivo, double vector[])
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
        vector[num_leidos]=valor; 
        num_leidos=num_leidos+1; 
    } 

    fclose(archivo);


    return;

}

void reescalamiento(double r[][2],double m[],double v[][2], int size) //Necesito pasarle las posiciones, velocidades de los planetas y sus masas
{
    int i,j;

    for(i=0;i<size;i++)
    {
        m[i]=m[i]/Ms;   //Reescalamiento de la masa

        for(j=0;j<2;j++)
        {
            r[i][j]=r[i][j]/c;//Reescalamiento de la posición
            v[i][j]=v[i][j]*sqrt(c/(G*Ms)); //Reescalamiento de la velocidad
        }
    }

    return;
}

void aceleracion(double r[][2],double m[],double a[][2],int size)
{
    int i,j;
    double mod;

    for ( i = 0; i < size; i++)
        for ( j = 0; j < 2; j++)
        {
            a[i][j]=0;
        }
    


    for ( i = 0; i < size; i++)
    {
       
        for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
               mod=sqrt(pow(r[i][0]-r[j][0],2)+pow(r[i][1]-r[j][1],2));

               a[i][0]-=(m[j]*(r[i][0]-r[j][0]))/(pow(mod,3));
               a[i][1]-=(m[j]*(r[i][1]-r[j][1]))/(pow(mod,3));
            }
            
        }
        
    }
    return;
}

void geocentrico(double r[][2],double rgeo[][2],int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < 2; j++)
        {
            rgeo[i][j]=r[i][j]-r[3][j];
        }
        
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
    double mod;

    for ( i = 0; i < size; i++)
        for ( j = 0; j < size; j++)
        {
            if (i!=j)
            {
                mod=sqrt(pow(r[i][0]-r[j][0],2)+pow(r[i][1]-r[j][1],2));
                V-=m[i]*m[j]/mod;
            }
            
        }
        
    return V;
  
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
