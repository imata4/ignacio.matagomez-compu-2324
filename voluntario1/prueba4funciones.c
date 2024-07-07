#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define pi 3.141592

const int N=20; //Número de partículas
const float L=10; //Tamaño del sistema
const double dt=0.002; //Paso temporal
const double dr=L*0.5; //Desviación máxima de la posición
const double v0=1; //Velocidad inicial máxima 

void distribucion_aleatoria(double x[],double y[]);
void distribucion_cuadrada(double x[],double y[]);
void distribucion_hexagonal(double x[],double y[]);
void generar_velocidades(double vx[],double vy[],double v[]);
double actualiza_posiciones(double x[], double y[],double vx[],double vy[],double wx[], double wy[],double momento);
void aceleracion(double x[],double y[],double ax[],double ay[]);
double Ecinetica(double vx[],double vy[],int size);
double Epotencial(double x[],double y[],int size,float L);
void promedios (double vx[],double vy[],double v[],double *vxm,double *vym,double *vm);
double distancia(double *rx, double *ry);
double fluctuaciones(double x[], double y[],double x0[],double y0[]);
double sepcuadmed(double x[],double y[]);

void print_v_stats(double distribucion_v[100][4]);
void estadistica_velocidades(double vx[],double vy[],double distribucion_v[100][4]);


int main (void)
{
    int i,j;

    double t,tf=50; //Contador de tiempo y tiempo final

    double x[N],y[N]; //Componentes x e y del vector posición
    double x0[N],y0[N]; //Componentes x e y del vector posición en el instante inicial   
    double vx[N],vy[N]; //Componentes x e y del vector velocidad
    double v[N]; //Módulo de la velocidad
    double vxm,vym,vm; //Promedio para la componente x,y de la velocidad y su módulo. Promedio a todas las partículas
    double ax[N],ay[N]; //Componentes x e y del vector aceleración
    double momento=0;
    double modulo;
    double theta[N]; //Ángulo
    double aux1,aux2;

    double Ec,Ep,E; //Energía cinética, potencial y energía total
    double wx[N],wy[N]; //Vectores auxiliares para algoritmo de Verlet
    double T;//Temperatura
    double P; //Presión
    double fluctuacion; 

    double distribucion_v[100][4];

    FILE *f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9,*f10,*f11,*f12;
    f1= fopen("trayectorias.txt", "w"); 
    f2= fopen("energiacinetica.txt", "w");
    f3= fopen("energiapotencial.txt", "w"); 
    f4= fopen("energiatotal.txt", "w");
    f5= fopen("temperatura.txt", "w");
    f6= fopen("v.txt", "w");
    f7= fopen("vx.txt", "w");
    f8= fopen("vy.txt", "w");
    f9= fopen("presion.txt", "w");
    f10= fopen("ecuacionestado.txt", "w");
    f11= fopen("ecuacionestadomedia.txt", "a");
    f12=fopen("fluctuaciones.txt", "w");

    srand(time(NULL));


    distribucion_aleatoria(x,y);
    //distribucion_cuadrada(x,y);
    //distribucion_hexagonal(x,y);

    generar_velocidades( vx, vy, v); 

    //Asigno a los vectores x0 y0 las posiciones iniciales de las partículas, para después estudiar las fluctuaciones
    for ( i = 0; i < N; i++)
    {
        x0[i]=x[i];
        y0[i]=y[i];
    }
    //Inicializo el array de distribucionv
    for(int k=0; k<100; k++)
    {
        distribucion_v[k][1] = 0;
        distribucion_v[k][2] = 0;
        distribucion_v[k][3] = 0;
    }

    //Calculo las aceleraciones
    aceleracion(x,y,ax,ay);

    for ( t = 0; t <= tf; t+=dt)
    {   
        //Calculo las energías:
        Ec=Ecinetica(vx,vy,N)/N;
        Ep=Epotencial(x,y,N,L)/N;
        E=Ec+Ep;

        for (int i = 0; i < N; i++) 
        {  
            fprintf(f1, "%lf", x[i]); 
            fprintf(f1,"," ); 
            fprintf(f1, "%lf", y[i]); 

            fprintf(f1, "\n"); // Nueva línea al final de cada fila
        }
        fprintf(f1, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones
     
        fprintf(f2, "%lf, %lf ", t, Ec); 
        fprintf(f2, "\n"); 

        fprintf(f3, "%lf, %lf ", t, Ep); 
        fprintf(f3, "\n");

        fprintf(f4, "%lf, %lf ", t, E); 
        fprintf(f4, "\n");

        //Calculo la temperatura para un determinado intervalo de tiempo
        if(t>=20 && t<=50)
        {
            T=Ecinetica(vx,vy,N)/N;
            aux1+=T; //Para calcular la temperatura media en un intervalo de tiempo, quedaría dividir entre el intervalo de tiempo en este caso 30s
            fprintf(f5, "%lf, %lf ", t, T); 
            fprintf(f5, "\n");
        }

        //Imprimo los promedios de las velocidades (v,vx y vy) para un tiempo concreto
        if (t>=20 && t<=50)
        {
            promedios(vx,vy,v,&vxm,&vym,&vm);
            estadistica_velocidades(vx,vy,distribucion_v);
            
            //fprintf(f6, "%lf, %lf,%lf",vm,vxm,vym); 
            fprintf(f6, "%lf, %lf ", t, vm); 
            fprintf(f6, "\n");
            fprintf(f7, "%lf, %lf ", t, vxm); 
            fprintf(f7, "\n");
            fprintf(f8, "%lf, %lf ", t, vym); 
            fprintf(f8, "\n");
        }

    

        //Imprimo la presión en función del tiempo
        if (t>=20 && t<=50)
        {
            
            P=momento/(t*L*L);
            aux2+=P; //Para calcular la presion media en un intervalo de tiempo, quedaría dividir entre el intervalo de tiempo en este caso 30s
            fprintf(f9, "%lf, %lf", t, P); 
            fprintf(f9, "\n");

            //Imprimo la ecuacion de estado
            fprintf(f10, "%lf, %lf", T, P); 
            fprintf(f10, "\n");
        }

        //Para el apartado 6 imprimo las fluctuaciones en función de la temperatura, para ver como cambian
        if(t>=20 && t<=50)
        {
            fluctuacion=fluctuaciones(x,y,x0,y0);
            fprintf(f12, "%lf,%lf", T, fluctuacion); 
            fprintf(f12, "\n");
        }
        
        //Para el apartado 7 imprimo la separación cuadrática media

        for ( i = 0; i < N; i++)
        {
            wx[i]=vx[i]+(dt*ax[i])/2.0;
            wy[i]=vy[i]+(dt*ay[i])/2.0;
        } 
        //Establezco las condiciones de contorno periódicas para actualizar las posiciones y calculo del momento
    
        momento=actualiza_posiciones(x,y,vx,vy,wx,wy,momento);

        //Evalúo las aceleraciones con mis nuevas posiciones 

        aceleracion(x,y,ax,ay);
    
        //Evalúo las velocidades con mis nuevas aceleraciones

        for ( i = 0; i < N; i++)
        {
            vx[i]= wx[i]+(dt*ax[i])/2.0;
            vy[i]= wy[i]+(dt*ay[i])/2.0;
            v[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]);//Calculo también el módulo de la velocidad
        }

        //Para poner las condiciones para un t concreto hay que tener en cuenta la naturaleza 
        //de los números en coma flotante, por tanto es posible que t=20 no se dé nunca
        //incluiré un pequeño intervalo donde es válida la condición
        //Apartado 6    
        /*if( ((t >= 20) && (t < 20+dt)) || ((t >= 30) && (t < 30+dt)) || ((t >= 35) && (t < 35+dt)) || ((t >= 45) && (t < 45+dt)))
        {
            for ( i = 0; i < N; i++)
            {
                vx[i]*=1.5;
                vy[i]*=1.5;
            }
        }    
        */
        //Apartado 7, reescalamiento de velocidades
        //Fijo la tolerancia en la mitad del paso temporal
        /*if (fabs(fmod(t, 60) < 0.5*dt))
        {
            for ( i = 0; i < N; i++)
            {
                vx[i]*=1.1;
                vy[i]*=1.1;
            }
        }*/

    }
    fprintf(f10, "%lf, %lf", dt*aux1/30, dt*aux2/30); 
        fprintf(f10, "\n");

    fprintf(f11, "%lf, %lf", dt*aux1/30, dt*aux2/30); 
            fprintf(f11, "\n");
    print_v_stats(distribucion_v);

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    fclose(f6);
    fclose(f7);
    fclose(f8);
    fclose(f9);
    fclose(f10);
    fclose(f11);
    fclose(f12);
    
    return 0;
}

double distancia(double *rx, double *ry)
{
    if (fabs(*rx) > 0.5 * L)
        *rx = (L - fabs(*rx)) * (-*rx) / fabs(*rx);

    if (fabs(*ry) > 0.5 * L)
        *ry = (L - fabs(*ry)) * (-*ry) / fabs(*ry);

    return sqrt((*rx) * (*rx) + (*ry) * (*ry));
}

void distribucion_aleatoria(double x[],double y[])
{
    double min_distance = 0.7;
    double rx, ry;
    bool condition = false;
    srand(time(NULL));

    double rnd1 = (double)rand() / (double)(RAND_MAX + 1);
    double rnd2 = (double)rand() / (double)(RAND_MAX + 1);
    x[0] = 0.5*L+2*(rnd1-0.5)*dr;
    y[0] = 0.5*L+2*(rnd2-0.5)*dr;
    for (int i = 1; i < N; i++)
    {
        condition = true;
        while (condition == true)
        {
            double rnd1 = (double)rand() / (double)(RAND_MAX + 1);
            double rnd2 = (double)rand() / (double)(RAND_MAX + 1);

            x[i] = 0.5*L+2*(rnd1-0.5)*dr;
            y[i] = 0.5*L+2*(rnd2-0.5)*dr;

            condition = false;

            for (int j = i - 1; j >= 0; j--)
            {
                rx = x[i] - x[j];
                ry = y[i] - y[j];

                if (distancia(&rx,&ry) < min_distance)
                {
                    condition = true;
                    break;
                }
            }
        }
    }
}

void distribucion_cuadrada(double x[],double y[])
{
    int j=0;
    for (int i = 0; i < N; i++)
    {
        if ((i % 4 == 0) && (i != 0))
            j++;
        x[i] = L / 8.0 + (L / 4.0) * (i % 4);
        y[i] = L / 8.0 + (L / 4.0) * j;
    }
}
void distribucion_hexagonal(double x[],double y[])
{
    int j=0;
    for (int i = 0; i < N; i++)
    {
        if ((i % 4 == 0) && (i != 0))
            j++;
        x[i] = L / 8.0 + (L / 4.0) * (i % 4) - (L / 8.0) * (j % 2);
        y[i] = L / 8.0 + (L / 4.0) * j;
    }
}

void generar_velocidades(double vx[],double vy[],double v[])
{   
    srand(time(NULL));
    double theta[N];
    for (int i = 0; i < N; i++)
    {
        double rnd3 = (double)rand() / (double)(RAND_MAX + 1);
        //theta[i]=rnd3*2*pi;
        //theta[i]=rnd3*pi+3*0.5*pi; //Si quiero que vx sea positiva pondría esto

        //Pongo velocidad aleatoria entre 0 y v0
        vx[i]=v0*cos(rnd3*2*pi);
        vy[i]=v0*sin(rnd3*2*pi);
        
        //Para apartado 4
        //vx[i]=vy[i]=0;

        v[i]=sqrt(vx[i]*vx[i] + vy[i]*vy[i]); 
        
    }
    return;
    
}




void aceleracion(double x[],double y[],double ax[],double ay[])
{
    int i,j;
    double r,rx,ry;
    i=j=0;
    r=rx=ry=0;
    for ( i = 0; i < N; i++){
        ax[i]=ay[i]=0;
    }
    
    for ( i = 0; i < (N-1); i++) //De esta forma solo calculo una vez la interaccion entre i y j
    {
        for ( j = i+1; j < N; j++)
        {  
                rx=x[i]-x[j];

                ry=y[i]-y[j];
                  
                r=distancia(&rx,&ry);
                
                if (r<=3)
                {
                    ax[i]+=24*((2.0/pow(r,13)-1.0/pow(r,7))-(2.0/pow(3,13)-1.0/pow(3,7)))*(rx/r);
                    ay[i]+=24*((2.0/pow(r,13)-1.0/pow(r,7))-(2.0/pow(3,13)-1.0/pow(3,7)))*(ry/r);
                    ax[j]-=24*((2.0/pow(r,13)-1.0/pow(r,7))-(2.0/pow(3,13)-1.0/pow(3,7)))*(rx/r);
                    ay[j]-=24*((2.0/pow(r,13)-1.0/pow(r,7))-(2.0/pow(3,13)-1.0/pow(3,7)))*(ry/r);  
                }
            
        }   
    }
    return;
}

double actualiza_posiciones(double x[], double y[],double vx[],double vy[],double wx[], double wy[],double momento)
{   

    for ( int i = 0; i < N; i++)
        {
            x[i]+=dt*wx[i];
            y[i]+=dt*wy[i];
            if (x[i]>L)
            {
                x[i]=fmod(x[i], L);
                momento+=fabs(2*vx[i]);
            }              

            if (x[i]<0)
            {
                x[i]=L-fmod(fabs(x[i]),L);
                momento+=fabs(2*vx[i]);
            }    
                
            if (y[i]>L)
            {
                y[i]=fmod(y[i],L);
                momento+=fabs(2*vy[i]);
            }

            if (y[i]<0)
            {
                y[i]=L-fmod(fabs(y[i]),L); 
                momento+=fabs(2*vy[i]);
            }

        }  
    return momento;  
}




double Ecinetica(double vx[],double vy[],int size)
{
    int i;
    double Ecinetica=0;

    for ( i = 0; i < size; i++)
        Ecinetica+=0.5*(vx[i]*vx[i]+vy[i]*vy[i]); //m=1 en nuestras unidades

    return Ecinetica;
}

double Epotencial(double x[],double y[],int size,float L)
{
    int i,j;
    double V=0;
    double r,rx,ry;

    for ( i = 0; i < size-1; i++)
        for ( j = i+1; j < size; j++)
        {
            if (i!=j)
            {
                rx=fabs(x[i]-x[j]);

                ry=fabs(y[i]-y[j]);
                  
                r=distancia(&rx,&ry);

                if(r<=3)
                    V+=(4/pow(r,12)-4/pow(r,6))-(4/pow(3,12)-4/pow(3,6));
                    //Resto el potencial a partir de la distancia de corte para mejorar la eficiencia computacional y evitar cálculos para partículas muy lejanas
            }
            
        }
        
    return V;
  
}


void promedios (double vx[],double vy[],double v[],double *vxm,double *vym,double *vm)
{
    int i;
    *vxm=*vym=*vm=0;

    for ( i = 0; i < N; i++)
    {
        *vxm+=vx[i];
        *vym+=vy[i];
        *vm+=v[i];
    }
    return;
}

void estadistica_velocidades(double vx[],double vy[],double distribucion_v[100][4])
{
    double v_module;
    double size = 0.2; // Size of the intervals

    for(int k=0; k<100; k++)
    {
        distribucion_v[k][0] = (k-50)*size; //Creo 100 intervalos en el histograma
    }

    // Registramos la velocidad de cada partícula (|V|, Vx and Vy)
    for(int i=0; i<N; i++)
    {
        v_module = sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
        
        for(int j=0; j<100; j++)
        {
            if( (v_module >= size*(j-50)) && (v_module < size*(j-49)))
            {
                distribucion_v[j][1] += 1;
            }
            if( (vx[i] >= size*(j-50)) && (vx[i] < size*(j-49)))
            {
                distribucion_v[j][2] += 1;
            }
            if( (vy[i] >= size*(j-50)) && (vy[i] < size*(j-49)))
            {
                distribucion_v[j][3] += 1;
            }
        }
    }
}

void print_v_stats(double distribucion_v[100][4])
{
    FILE *f1,*f2,*f3; 
    double sum_v = 0, sum_vx = 0, sum_vy = 0;

    for(int i=0; i<100; i++)
    {
        sum_v += distribucion_v[i][1];
        sum_vx += distribucion_v[i][2];
        sum_vy += distribucion_v[i][3];
    }
    
    // We normalize the total data collected by taken into account the total number of values and the size of the intervals (0.2)
    for(int j=0; j<100; j++)
    {
        distribucion_v[j][1] *= 1.0/(sum_v * 0.2) ;
        distribucion_v[j][2] *= 1.0/(sum_vx * 0.2);
        distribucion_v[j][3] *= 1.0/(sum_vy * 0.2);
    }

    //cout << "Sum_v = " << sum_v << " | Sum_vx = " << sum_vx << " | Sum_vy = " << sum_vy << endl;

    // We print the data for the Python script to read
    f1=fopen("v1.txt","w");
    f2=fopen("vx1.txt","w");
    f3=fopen("vy1.txt","w");
    for(int k=0; k<100; k++)
    {
        //v_stats << v_distribution[k][0] << " " << v_distribution[k][1] << " " << v_distribution[k][2] << " " << v_distribution[k][3] << endl;
        /*fprintf(f1,"%lf,%lf",fabs(distribucion_v[k][0]),fabs(distribucion_v[k][1]));
        fprintf(f1, "\n");
        fprintf(f2,"%lf,%lf",fabs(distribucion_v[k][0]),fabs(distribucion_v[k][2]));
        fprintf(f2, "\n");
        fprintf(f3,"%lf,%lf",fabs(distribucion_v[k][0]),fabs(distribucion_v[k][3]));
        fprintf(f3, "\n");*/
        fprintf(f1,"%lf,%lf,%lf,%lf",fabs(distribucion_v[k][0]),fabs(distribucion_v[k][1]),fabs(distribucion_v[k][2]),fabs(distribucion_v[k][3]));
        fprintf(f1, "\n");
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);

}

void reescalar_velocidad(double vx[],double vy[],double t)
{
    double factor;
    factor=1.5; //Apartado 6
    //factor=1.1; //Apartado 7

    if( ((t >= 20) && (t < 20+dt)) || ((t >= 30) && (t < 30+dt)) || ((t >= 35) && (t < 35+dt)) || ((t >= 45) && (t < 45+dt))) //Apartado 6
    //if (fabs(fmod(t, 60) < 0.5*dt)) //Apartado 7
    {   
        for(int i=0; i<N; i++)
        {
            vx[i] *= factor;
            vy[i] *= factor;
        }

    }
    return;    
    
}

double fluctuaciones(double x[], double y[],double x0[],double y0[])
{
    double fluctuacion; //De la primera partícula (Cambiar índice de x[],y[] para considerar otra)
    double rx,ry;
    //fluctuacion=pow(sqrt(x[0]*x[0]+y[0]*y[0])-sqrt(x0[0]*x0[0]+y0[0]*y0[0]),2);
    rx=fabs(x[0]-x0[0]);
    ry=fabs(y[0]-y0[0]);
    fluctuacion=distancia(&rx,&ry);

    return fluctuacion;
}

double sepcuadmed(double x[],double y[])
{
    double sepcuadmed;
    //Separación cuadrática media entre las partículas 1 y 2, cambiar índices si se quieren considerar otras partículas   
    sepcuadmed+=pow(sqrt(x[0]*x[0]+y[0]*y[0])-sqrt(x[1]*x[1]+y[1]*y[1]),2);
        
    return sepcuadmed;
    
}

