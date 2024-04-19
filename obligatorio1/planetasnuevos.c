#include <stdio.h>
#include <math.h>

//Declaración de constantes a utilizar
#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11

void escribirEnArchivo(FILE *archivo, int dato);
void leerArchivo(const char *nombre_archivo, double vector[]);
void reescalamiento(double r[][2],double m[],double v[][2], int size);
void aceleracion(double r[][2],double m[],double a[][2],int size);
void geocentrico(double r[][2],double rgeo[][2],int size);
void escribeMatriz(const char *nombre_archivo, double r[][2], int num_elementos);

void main(void)
{   
    int i,j; //Variables auxiliares para bucles
    double h; //Paso temporal
    h=0.01;

    int p=0;

    int n; //Número de planetas
    n=6; 

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


    FILE *f1,*f2;
    f1= fopen("resultados.txt", "w"); 
    f2= fopen("geocentrico.txt", "w"); 

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

    tf=1;
    //Paso el tiempo a segundos y lo reescalo
    tf=tf*365*24*3600;
    tf=tf*sqrt((G*Ms)/(pow(c,3)));
    
    //Calculo las aceleraciones
    aceleracion(r,m,a,n);


    for ( t = 0; t < tf; t+=h)
    {
        //Esta función la utilizo para calcular la posición geocéntrica
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

        //calculoPeriodo(r,aux,t,periodo,n); //Calculo los períodos si es que ya han recorrido un periodo

   

        //METO LOS RESULTADOS EN UN ARCHIVO
 
        

   

        

        
 
        


    }
    //Cierro todos los archivos
    fclose(f1);
    fclose(f2);

    return;
}

void escribirEnArchivo(FILE *archivo, int dato) 
{
    fprintf(archivo, "%d\n", dato);
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
void escribeMatriz(const char *nombre_archivo, double r[][2], int num_elementos)
{
    FILE *archivo = fopen(nombre_archivo, "w"); // Abre el archivo en modo añadir,de forma que tras cada iteración añada, nuevos datos al archivo
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