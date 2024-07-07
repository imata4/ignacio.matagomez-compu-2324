#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "FuncionesPenduloDoble.h"
#define pi 3.14

int main(){
    FILE *fichero_out;
    fichero_out=fopen("PenduloDoble.txt", "w");
    double g=9.81, h, t, Tmax, m;
    double y[4];
    double k[4][4];
    double E;
    int i, j;
    double x1, y1, x2, y2;
    x1=0;
    y1=0;
    x2=0;
    y2=0;

    Tmax=1000;
    t=0;
    //Condiciones iniciales
    E=1;
    y[0]=0.3; //Psi
    y[1]=0.1; //Phi
    y[3]=sqrt(E-2*g*(1-cos(y[1]))-g*(1-cos(y[0]))); //Momento de Phi
    y[2]=y[3]*cos(y[0]-y[1]); //Momento de Psi

    y[3]=2*y[3];

    h=0.001;  
    while(t<Tmax){
        //Calculo de k1:
        k[0][0]=h*fPsi(y[1], y[3], y[0], y[2]);
        k[0][1]=h*fPhi(y[1], y[3], y[0], y[2]);
        k[0][2]=h*fMomentoPsi(y[1], y[3], y[0], y[2]);
        k[0][3]=h*fMomentoPhi(y[1], y[3], y[0], y[2]);
    
        //Calculo de k2:
        k[1][0]=h*fPsi(y[1]+k[0][1]/2, y[3]+k[0][3]/2, y[0]+k[0][0]/2, y[2]+k[0][2]/2);
        k[1][1]=h*fPhi(y[1]+k[0][1]/2, y[3]+k[0][3]/2, y[0]+k[0][0]/2, y[2]+k[0][2]/2);
        k[1][2]=h*fMomentoPsi(y[1]+k[0][1]/2, y[3]+k[0][3]/2, y[0]+k[0][0]/2, y[2]+k[0][2]/2);
        k[1][3]=h*fMomentoPhi(y[1]+k[0][1]/2, y[3]+k[0][3]/2, y[0]+k[0][0]/2, y[2]+k[0][2]/2);
        
        //Calculo de k3:
        k[2][0]=h*fPsi(y[1]+k[1][1]/2, y[3]+k[1][3]/2, y[0]+k[1][0]/2, y[2]+k[1][2]/2);
        k[2][1]=h*fPhi(y[1]+k[1][1]/2, y[3]+k[1][3]/2, y[0]+k[1][0]/2, y[2]+k[1][2]/2);
        k[2][2]=h*fMomentoPsi(y[1]+k[1][1]/2, y[3]+k[1][3]/2, y[0]+k[1][0]/2, y[2]+k[1][2]/2);
        k[2][3]=h*fMomentoPhi(y[1]+k[1][1]/2, y[3]+k[1][3]/2, y[0]+k[1][0]/2, y[2]+k[1][2]/2);   

        //Calculo de k4
        k[3][0]=h*fPsi(y[1]+k[2][1], y[3]+k[2][3], y[0]+k[2][0], y[2]+k[2][2]);
        k[3][1]=h*fPhi(y[1]+k[2][1], y[3]+k[2][3], y[0]+k[2][0], y[2]+k[2][2]);
        k[3][2]=h*fMomentoPsi(y[1]+k[2][1], y[3]+k[2][3], y[0]+k[2][0], y[2]+k[2][2]);
        k[3][3]=h*fMomentoPhi(y[1]+k[2][1], y[3]+k[2][3], y[0]+k[2][0], y[2]+k[2][2]);    
        
        //Calculo de la nueva y
        for(i=0;i<4;i++){
            y[i]=y[i]+(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;
        }
        
        x1=sin(y[1]);
        y1=-cos(y[1]);
        x2=x1+sin(y[0]);
        y2=y1-cos(y[0]);

        t=t+h;

        fprintf(fichero_out, "%lf, %lf \n", x1, y1);
        fprintf(fichero_out, "%lf, %lf \n", x2, y2);
        fprintf(fichero_out, "\n");
        }
    fclose(fichero_out);


    return 0;
}