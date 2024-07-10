#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define g 9.80665
#define pi 3.14156


double f_phi(double phi, double psi, double pphi, double ppsi);
double f_psi(double phi, double psi, double pphi, double ppsi);
double f_pphi(double phi, double psi, double pphi, double ppsi);
double f_ppsi(double phi, double psi, double pphi, double ppsi);

int main(void)
{
    clock_t begin= clock(); //Empiezo a contar el tiempo de simulación

    int i, j; //Contadores
    double h; //Paso temporal
    double t; //Contador de tiempo
    double tf; //Tiempo final de la simulación

    double y[4]; //Vector que contiene las variables phi, psi, pphi, ppsi respectivamente
    double k[4][4]; //Matriz que almacena los 4 k de las 4 variables correspondientes al algoritmo Runge-Kutta de cuarto orden
    double E; //Energía del sistema
    
    double x1, y1; //Posición en x e y del primer péndulo
    double x2, y2; //Posición en x e y del segundo péndulo

    FILE *f1,*f2,*f3,*f4;
    f1=fopen("pendulodoble1.txt", "w"); 
    f2=fopen("poincare2_phi_psi_E=5.txt", "w");
    f3=fopen("poincare2_phi_phipunto_E=5.txt", "w");
    f4=fopen("lyapunov E=1.txt", "w");

    //Establecemos las variables que regulan el tiempo de simulación
    h=0.0001;
    tf=100;
    t=0;
    //Condiciones iniciales
    E=15;
    y[0]=0.1; //Phi
    y[1]=0.3; //Psi

    //Para calcular pphi fijamos la velocidad en psi a 0
    y[2]=pow((E-2*g*(1-cos(y[1]))-g*(1-cos(y[0]))),0.5); //Momento de Phi  //Hay que introducir valores iniciales que permitan que la siguiente raíz no sea nula
    y[3]=y[2]*cos(y[1]-y[0]); //Momento de Psi

    y[2]=2*y[2];

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

    //LYAPUNOV

    int num_steps = (int)(tf / h);
    double ypert[4];
    double perturbacion[4] = {0, -0.05, 0, 0};
    double sumadivergente = 0.0;
    double norm_perturbation;
    double exponente_lyapunov;

    for (i = 0; i < 4; i++) 
        ypert[i] = y[i] + perturbacion[i];
    

   /* for ( tf = 1000; tf <= 10000; tf+=100)
    {
        t=0;
        sumadivergente=0;*/
        while(t<tf)
        {
            //Calculo de k1:
            k[0][0]=h*f_phi(y[0], y[1], y[2], y[3]);
            k[0][1]=h*f_psi(y[0], y[1], y[2], y[3]);
            k[0][2]=h*f_pphi(y[0], y[1], y[2], y[3]);
            k[0][3]=h*f_ppsi(y[0], y[1], y[2], y[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_phi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            k[1][1]=h*f_psi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            k[1][2]=h*f_pphi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            k[1][3]=h*f_ppsi(y[0]+0.5*k[0][0], y[1]+0.5*k[0][1], y[2]+0.5*k[0][2], y[3]+0.5*k[0][3]);
            
            //Calculo de k3:
            k[2][0]=h*f_phi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            k[2][1]=h*f_psi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            k[2][2]=h*f_pphi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);
            k[2][3]=h*f_ppsi(y[0]+0.5*k[1][0], y[1]+0.5*k[1][1], y[2]+0.5*k[1][2], y[3]+0.5*k[1][3]);

            //Calculo de k4
            k[3][0]=h*f_phi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            k[3][1]=h*f_psi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            k[3][2]=h*f_pphi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            k[3][3]=h*f_ppsi(y[0]+k[2][0], y[1]+k[2][1], y[2]+k[2][2], y[3]+k[2][3]);
            
            //Calculo de la nueva y
            for(i=0;i<4;i++)
                y[i]=y[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;
            
            
            x1=sin(y[1]);
            y1=-cos(y[1]);
            x2=x1+sin(y[0]);
            y2=y1-cos(y[0]);

            //Aqui iba paso temporal

            //Imprimo las trayectorias en un fichero y los mapas de Poincaré en otro
            fprintf(f1, "%lf, %lf \n", x1, y1);
            fprintf(f1, "%lf, %lf \n", x2, y2);
            fprintf(f1, "\n");

            fprintf(f2, "%lf, %lf", y[0], y[1]);
            fprintf(f2, "\n");

            fprintf(f3, "%lf, %lf", y[0], y[2]);
            fprintf(f3, "\n");

            //Procedo a computar las trayectorias perturbadas

            /*//Calculo de k1:
            k[0][0]=h*f_phi(ypert[0], ypert[1], ypert[2], ypert[3]);
            k[0][1]=h*f_psi(ypert[0], ypert[1], ypert[2], ypert[3]);
            k[0][2]=h*f_pphi(ypert[0], ypert[1], ypert[2], ypert[3]);
            k[0][3]=h*f_ppsi(ypert[0], ypert[1], ypert[2], ypert[3]);
        
            //Calculo de k2:
            k[1][0]=h*f_phi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            k[1][1]=h*f_psi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            k[1][2]=h*f_pphi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            k[1][3]=h*f_ppsi(ypert[0]+0.5*k[0][0], ypert[1]+0.5*k[0][1], ypert[2]+0.5*k[0][2], ypert[3]+0.5*k[0][3]);
            
            //Calculo de k3:
            k[2][0]=h*f_phi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);
            k[2][1]=h*f_psi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);
            k[2][2]=h*f_pphi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);
            k[2][3]=h*f_ppsi(ypert[0]+0.5*k[1][0], ypert[1]+0.5*k[1][1], ypert[2]+0.5*k[1][2], ypert[3]+0.5*k[1][3]);

            //Calculo de k4
            k[3][0]=h*f_phi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            k[3][1]=h*f_psi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            k[3][2]=h*f_pphi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            k[3][3]=h*f_ppsi(ypert[0]+k[2][0], ypert[1]+k[2][1], ypert[2]+k[2][2], ypert[3]+k[2][3]);
            
            //Calculo de la nueva ypert
            for(i=0;i<4;i++)
                ypert[i]=ypert[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;

*/
        /*  double dist = 0.0;
        for (i = 0; i < 4; i++)
        {
            dist += pow(ypert[i] - y[i], 2);
        }
        dist = sqrt(dist);
        
        if (dist > 0)
        {
            // Normalización de la perturbación
            for (i = 0; i < 4; i++)
            {
                ypert[i] = y[i] + (perturbacion[i] / dist) * 0.05;
            }

            // Acumular el logaritmo natural de la distancia
            sumadivergente += log(dist);
        }

        // Calcular el exponente de Lyapunov en este punto
        exponente_lyapunov = sumadivergente / (t + h);

        // Guardar el valor del exponente de Lyapunov en el archivo
        fprintf(f4, "%lf, %lf\n", t, exponente_lyapunov);

        t += h;    */

            /*for (i = 0; i < 4; i++) 
            {
                perturbacion[i] = ypert[i] - y[i]; //Diferencia entre las trayectorias
            }

            norm_perturbation = sqrt(perturbacion[0] * perturbacion[0] + perturbacion[1] * perturbacion[1] +
                                    perturbacion[2] * perturbacion[2] + perturbacion[3] * perturbacion[3]);
            sumadivergente += log(fabs(norm_perturbation)/ 0.05);

            for (i = 0; i < 4; i++) 
            {
                perturbacion[i] *= 0.05/ norm_perturbation;
                ypert[i] = y[i] + perturbacion[i];
            }

            if((int)t%100==0)
            {
                exponente_lyapunov=sumadivergente/t ;
                fprintf(f4,"%lf, %lf",t, exponente_lyapunov);
                fprintf(f4, "\n");
            }*/

            //Añado un paso temporal
            t+=h;

        }

      /*  exponente_lyapunov=sumadivergente/tf ;
        fprintf(f4,"%lf, %lf",tf, exponente_lyapunov);
        fprintf(f4, "\n");
    }*/
            
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

double f_phi(double phi, double psi, double pphi, double ppsi)
{
    double fphi=(pphi-ppsi*cos((psi-phi)))/(2-pow(cos((phi-psi)),2));

    return fphi;
}

double f_psi(double phi, double psi, double pphi, double ppsi)
{
    double fpsi=(2*ppsi-pphi*cos((psi-phi)))/(2-pow(cos((phi-psi)),2));

    return fpsi;
}



double f_pphi(double phi, double psi, double pphi, double ppsi)
{
    double fpphi=-2*g*sin(phi)+((ppsi*pphi*pow(cos(psi-phi),2)-(2*pow(ppsi,2)+pow(pphi,2)*cos(psi-phi))+2*pphi*ppsi)/(pow((2-pow(cos(psi-phi),2)),2)))*2*sin(phi-psi);

    return fpphi;
}

double f_ppsi(double phi, double psi, double pphi, double ppsi)
{
    double fppsi=-g*sin(psi)+((ppsi*pphi*pow(cos(psi-phi),2)-(2*pow(ppsi,2)+pow(pphi,2)*cos(psi-phi))+2*pphi*ppsi)/(pow((2-pow(cos(psi-phi),2)),2)))*2*sin(psi-phi);

    return fppsi;
}
