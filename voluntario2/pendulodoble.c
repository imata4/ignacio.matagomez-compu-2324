#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define g 9.80665
#define pi 3.14156


double f_vphi(double vphi);
double f_vpsi(double vpsi);
double f_aphi(double phi, double psi, double vphi, double vpsi);
double f_apsi(double phi, double psi, double vphi, double vpsi);

int main(void)
{
    clock_t begin= clock(); //Empiezo a contar el tiempo de simulación

    int i, j; //Contadores
    double h; //Paso temporal
    double t; //Contador de tiempo
    double tf; //Tiempo final de la simulación

    double y[4]; //Vector que contiene las variables phi, psi, vphi, vpsi respectivamente
    double k[4][4]; //Matriz que almacena los 4 k de las 4 variables correspondientes al algoritmo Runge-Kutta de cuarto orden
    double E; //Energía del sistema
    
    double x1, y1; //Posición en x e y del primer péndulo
    double x2, y2; //Posición en x e y del segundo péndulo

    //Variables para calculas los coeficientes de Lyapunov

    double ypert[4]; //Vector y perturbado
    double perturbacion[4] = {-0.05, -0.05, -0.05, -0.05}; //Vector que introduce las perturbaciones
    double sumadivergente = 0.0;
    double norma_perturbacion; //Lo utilizaremos para normalizar
    double exponente_lyapunov;
    int n; //Contador de iteraciones al calcular el exponente de Lyapunov

    FILE *f1,*f2,*f3,*f4;
    f1=fopen("pendulodoble.txt", "w"); 
    f2=fopen("poincare_phi_psi_E=15.txt", "w");
    f3=fopen("poincare_phi_phipunto_E=15.txt", "w");
    f4=fopen("lyapunov_E=15.txt", "w");

    //Establecemos las variables que regulan el tiempo de simulación
    h=0.01;
    tf=100;
    t=0;
    //Condiciones iniciales
    E=15;
    y[0]=0.10; //Phi
    y[1]=0.2; //Psi

    //Para calcular pphi fijamos la velocidad en psi a 0
    y[2]=sqrt(E+2*g*cos(y[0])+g*cos(y[1])); //Velocidad de Phi  //Hay que introducir valores iniciales que permitan que la siguiente raíz no sea nula
    y[3]=0; //Velocidad de Psi

    //Calculo e imprimo las posiciones iniciales
    x1=sin(y[0]);
    y1=-cos(y[0]);
    x2=x1+sin(y[1]);
    y2=y1-cos(y[1]);

    fprintf(f1, "%lf, %lf \n", x1, y1);
    fprintf(f1, "%lf, %lf \n", x2, y2);
    fprintf(f1, "\n");

    //Imprimo los valores iniciales de phi y psi para hacer el mapa de Poincaré
    fprintf(f2, "%lf, %lf", y[0], y[1]);
    fprintf(f2, "\n");

    //Imprimo los valores iniciales de phi y phi punto para hacer el mapa de Poincaré
    fprintf(f3, "%lf, %lf", y[0], y[2]);
    fprintf(f3, "\n");

  

    for (i = 0; i < 4; i++) 
        ypert[i] = y[i] + perturbacion[i];
    

    for ( tf = 1000; tf <= 10000; tf+=1000)
    {
        t=0;
        sumadivergente=0;
        n=0;
        while(t<tf)
        {
            //Calculo de k1:
            k[0][0]=h*f_vphi(y[2]);
            k[0][1]=h*f_vpsi(y[3]);
            k[0][2]=h*f_aphi(y[0], y[1], y[2], y[3]);
            k[0][3]=h*f_apsi(y[0], y[1], y[2], y[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_vphi(y[2]+0.5*k[0][2]);
            k[1][1]=h*f_vpsi(y[3]+0.5*k[0][3]);
            k[1][2]=h*f_aphi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            k[1][3]=h*f_apsi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            
            //Calculo de k3:
            k[2][0]=h*f_vphi(y[2]+0.5*k[1][2]);
            k[2][1]=h*f_vpsi(y[3]+0.5*k[1][3]);
            k[2][2]=h*f_aphi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            k[2][3]=h*f_apsi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);

            //Calculo de k4
            k[3][0]=h*f_vphi(y[2]+k[2][2]);
            k[3][1]=h*f_vpsi(y[3]+k[2][3]);
            k[3][2]=h*f_aphi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            k[3][3]=h*f_apsi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            
            //Calculo de la nueva y
            for(i=0;i<4;i++)
                y[i]=y[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;
            
            
            x1=sin(y[0]);
            y1=-cos(y[0]);
            x2=x1+sin(y[1]);
            y2=y1-cos(y[1]);

            

            //Imprimo las trayectorias en un fichero y los mapas de Poincaré en otro
            fprintf(f1, "%lf, %lf \n", x1, y1);
            fprintf(f1, "%lf, %lf \n", x2, y2);
            fprintf(f1, "\n");

            fprintf(f2, "%lf, %lf", y[0], y[1]);
            fprintf(f2, "\n");

            fprintf(f3, "%lf, %lf", y[0], y[2]);
            fprintf(f3, "\n");

            
            //Procedo a computar las trayectorias perturbadas

            //Calculo de k1:
            k[0][0]=h*f_vphi(ypert[2]);
            k[0][1]=h*f_vpsi(ypert[3]);
            k[0][2]=h*f_aphi(ypert[0], ypert[1], ypert[2], ypert[3]);
            k[0][3]=h*f_apsi(ypert[0], ypert[1], ypert[2], ypert[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_vphi(ypert[2]+0.5*k[0][2]);
            k[1][1]=h*f_vpsi(ypert[3]+0.5*k[0][3]);
            k[1][2]=h*f_aphi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            k[1][3]=h*f_apsi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            
            //Calculo de k3:
            k[2][0]=h*f_vphi(ypert[2]+0.5*k[1][2]);
            k[2][1]=h*f_vpsi(ypert[3]+0.5*k[1][3]);
            k[2][2]=h*f_aphi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);
            k[2][3]=h*f_apsi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);

            //Calculo de k4
            k[3][0]=h*f_vphi(ypert[2]+k[2][2]);
            k[3][1]=h*f_vpsi(ypert[3]+k[2][3]);
            k[3][2]=h*f_aphi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            k[3][3]=h*f_apsi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            
            //Calculo de la nueva ypert
            for(i=0;i<4;i++)
                ypert[i]=ypert[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;



            for (i = 0; i < 4; i++) 
            {
                perturbacion[i] = ypert[i] - y[i]; //Diferencia entre las trayectorias
            }

            norma_perturbacion = sqrt(perturbacion[0] * perturbacion[0] + perturbacion[1] * perturbacion[1] +
                                    perturbacion[2] * perturbacion[2] + perturbacion[3] * perturbacion[3]);
            sumadivergente += log(fabs(norma_perturbacion)/ 0.05);

            for (i = 0; i < 4; i++) 
            {
                perturbacion[i] *= 0.05/ norma_perturbacion;
                ypert[i] = y[i] + perturbacion[i];
            }

            n++;
      
           //Añado un paso temporal
            t+=h;
        }

        exponente_lyapunov=sumadivergente/(n*h) ;
        fprintf(f4,"%lf, %lf",tf, exponente_lyapunov);
        fprintf(f4, "\n");
    }
            
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);

     //Calculo el tiempo de ejecución
    clock_t end=clock();
    double tiempo= (double) (end-begin)/ CLOCKS_PER_SEC;
    printf("Tiempo de compilacion=%lf \n", tiempo);

    return 0;
}

double f_vphi(double vphi)
{
    return vphi;
}

double f_vpsi(double vpsi)
{
    return vpsi;
}


double f_aphi(double phi, double psi, double vphi, double vpsi) //Función para calcular la aceleración en phi (phi doble punto)
{
    double aphi= (g*sin(psi)*cos(phi-psi) - 2*g*sin(phi) - pow(vphi, 2)*cos(phi-psi)*sin(phi-psi) - pow(vpsi, 2)*sin(phi - psi))/(2 - pow(cos(phi-psi), 2));
    
    return aphi;
}

double f_apsi(double phi, double psi, double vphi, double vpsi) //Función para calcular la aceleración en psi (psi doble punto)
{
    double apsi = (g*sin(phi)*cos(phi-psi) - g*sin(psi) + 0.5*pow(vpsi, 2)*cos(phi-psi)*sin(phi-psi) + pow(vphi, 2)*sin(phi-psi))/(1 - 0.5*pow(cos(phi-psi), 2));
  
    return apsi;
}
