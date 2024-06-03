#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#define pi 3.14159

void main(void)
{   
    int n,i,j;
    double h=0.1; //Espaciado en la discretización espacial
    double s; //Paso temporal
    int N=1000; //Números de puntos del retículo
    int tmax=10000;

    double nciclos=250;//Toma valores de 0 a N/4
    double k_0;
    double lambda=0.2;
    double norma;
    double aux;

    double complex phi[N+1]; //Va d 0 a N
    double complex alpha[N];
    double complex beta[N];
    double complex chi[N+1];

    double V[N+1];
    double modulo[N+1];
    double complex b[N];

    double complex aux2;

    FILE *f1,*f2,*f3;
	f1 = fopen("sch.txt","w");
	f2 = fopen("pot.txt","w");
    f3 = fopen("animacion.txt","w");

    k_0=(2*pi*nciclos)/(N*1.0);

    //Cálculo del potencial
    
    for (i=0;i<N;i++)
	{
		if(i>=0.4*N && i<=0.6*N)
		    V[i] = lambda*k_0*k_0; 
		else V[i] =0;
	}

    //Defino s tilde
    s=1.0/(4*k_0*k_0);


    //Condiciones de contorno para phi
    phi[0]=0;
    phi[N]=0;

    norma = 0.0;
	for (i=1;i<N;i++) //Iniciar la función de onda según la fórmula
	{
        aux2=cexp(I*k_0*i);// Asigna valores a aux2
		aux=exp(-8.0*pow((4.0*i-N), 2.0)/(N*N));
		
		phi[i]=aux*aux2; 
		norma=norma+pow(cabs(phi[i]),2.0); //Esta función se usa para calcular el módulo de un número complejo
		
	}	
    /*
    //Normalizo la función de onda
    for (i=1;i<N;i++)
		phi[i]= (phi[i])/(sqrt(norma));
    */
		
	//Impongo condición para calcular alpha por recurrencia
	alpha[N-1]= 0;

    //Calculamos las alpha iniciales
	for (i=N-1;i>=0;i--)
	{
        alpha[i-1]=(-1.0)/((-2+I*(2/s)-V[i])+alpha[i]); //Es V[i+1] porque i va de 1 a N-1
	}

      //Imprimo los modulos (lo hago al principio para imprimir también los iniciales)
		
		for (i=0;i<=N;i++)
		{
			fprintf(f2,"%i\t%lf\n",i,V[i]);	
            fprintf(f3,"%f,%lf,%lf",i*h,pow(cabs(phi[i]), 2),V[i]);		
            //fprintf(f3,"%f,%lf",i*h,pow(cabs(phi[i]), 2));		
            fprintf(f3,"\n");
		}

        fprintf(f1,"%i,%lf",n,norma);
		fprintf(f1,"\n");
		fprintf(f2,"\n");
		fprintf(f2,"\n");
        fprintf(f3,"\n");

    //Comienzo el bucle temporal
    for ( n = 0; n <tmax ; n++)
    {
       

        //Calculo los b_j

        for ( i = 0; i < N; i++) 
            b[i]=I*4*phi[i]/s;

        //Calculo los coeficientes beta

        //Condición inicial
        beta[N-1]=0;

        for(i=N-1;i>=0;i--)	 
			beta[i-1]=(b[i]-beta[i])/((-2+I*(2/s)-V[i])+alpha[i]);
		
        //Calculo chi
        chi[N]=0;
        chi[0]=0;

        for ( i = 0; i < N; i++)
            chi[i+1]=alpha[i]*chi[i]+beta[i];
        
        //Calculo la función de onda para un instante posterior

        for ( i = 0; i <= N; i++)
            phi[i]=chi[i]-phi[i];

         //Calculo la norma de phi para cada distancia---> comprobar conservación
        norma=0.0;

        for (i=0;i<=N;i++)
		{
            norma = norma+ pow(cabs(phi[i]),2.0);
			
		}

          //Imprimo los modulos (lo hago al principio para imprimir también los iniciales)
		
		for (i=0;i<=N;i++)
		{
            fprintf(f3,"%f,%lf,%lf",i*h,pow(cabs(phi[i]), 2),V[i]);		
            fprintf(f3,"\n");
		}

        fprintf(f1,"%i,%lf",n,norma);
		fprintf(f1,"\n");
		fprintf(f2,"\n");
		fprintf(f2,"\n");
        fprintf(f3,"\n");
                
              
    }
    
    fclose(f1);
    fclose(f2);
    fclose(f3);

    return;
}