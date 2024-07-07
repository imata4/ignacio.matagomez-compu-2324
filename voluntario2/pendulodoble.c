#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define g 9.80665
#define pi 3.14156

void main(void)
{
    double h=0.001;
    double tf=1000; //Tiempo máximo de la simulación

    double E; //Energía total
	double phi,vphi,psi,vpsi; //Definición de ángulos del primer péndulo (phi) y del segundo péndulo (psi) junto con sus velocidades angulares
	double pphi,ppsi; //Momentos generalizados
	double k1phi,k1pphi,k1psi,k1ppsi;
	double k2phi,k2pphi,k2psi,k2ppsi;
	double k3phi,k3pphi,k3psi,k3ppsi;
	double k4phi,k4pphi,k4psi,k4ppsi; //ki para todas las coordenadas
	double aux,aux1;
	int i,j;
	double t=0; //tiempo
	
	double x1,y1; //Componentes x e y de la posición del primer péndulo  
	double x2,y2;//Componentes x e y de la posición del segundo péndulo  
	double m1,m2; //Masa del primer y el segundo péndulo
	double l1,l2; //Longitudes del primer y segundo péndulo

	FILE *f1;
	
	f1 = fopen("pendulo.txt","w");

    //Establecemos los parámetros a 1
    l1=l2=m1=m2=1;
	x1=x2=y1=y2=0;

    //Condiciones iniciales
	vpsi=0;
    phi=0.3;
    psi=0.1;
	E=2;
	//Hay que introducir valores iniciales que permitan que la siguiente raíz no sea nula
    pphi=sqrt(E-2*g*(1-cos(phi))-g*(1-cos(psi)));
	ppsi=pphi*cos(psi-phi);

	pphi=pphi*2;

	
	//pphi=2*vphi+vpsi*cos(psi-phi);//Salen de derivar el lagrangiano respecto a las velocidades angulares
	//ppsi=vpsi+vphi*cos(psi-phi);

    x1 = l1*sin(phi);
    y1 = -l1*cos(phi);
    x2 = x1 + l2*sin(psi);
    y2 = y1 - l2*cos(psi);

	fprintf (f1,"%lf\t",x1);
    fprintf(f1,",");
    fprintf (f1,"%lf\t",y1);
    fprintf(f1,"\n");
    fprintf (f1,"%lf\t",x2);
    fprintf(f1,",");
    fprintf (f1,"%lf\n",y2);
    fprintf(f1,"\n"); 

	j=0;
    while (t<tf)
    {
        //Calculamos las k1
		k1phi = h*(pphi-ppsi*cos(psi-phi))/(2-pow(cos(psi-phi),2));
		k1psi = h*(2*ppsi-pphi*cos(psi-phi))/(2-pow(cos(psi-phi),2));
	
		k1pphi = h*(-2*g*sin(phi)+2*sin(phi-psi)*((pphi*ppsi*pow(cos(psi-phi),2)-cos(psi-phi)*(2*pow(ppsi,2)+pow(pphi,2))+2*pphi*ppsi)/pow(2-pow(cos(psi-phi),2),2)));
	
		k1ppsi = h*(-g*sin(psi)+2*sin(psi-phi)*((pphi*ppsi*pow(cos(psi-phi),2)-cos(psi-phi)*(2*pow(ppsi,2)+pow(pphi,2))+2*pphi*ppsi)/pow(2-pow(cos(psi-phi),2),2)));
		
		//Calculamos las k2 usando las k1:
		
		k2phi = h*((pphi+k1pphi*0.5)-(ppsi+k1ppsi*0.5)*cos((ppsi+k1ppsi*0.5)-(pphi+k1pphi*0.5)))/(2-pow(cos((psi+k1psi*0.5)-(phi+k1phi*0.5)),2));
		k2psi = h*(-(pphi+0.5*k1pphi)*cos((psi+0.5*k1psi)-(phi+0.5*k1phi))+2*(ppsi+0.5*k1ppsi))/(2-pow(cos((psi+0.5*k1psi)-(phi+0.5*k1phi)),2));

		k2pphi = h*(-2*g*sin(phi+0.5*k1phi)+2*sin((phi+0.5*k1phi)-(psi+0.5*k1psi))*(((pphi+0.5*k1pphi)*(ppsi+0.5*k1ppsi)*pow(cos((psi+0.5*k1psi)-(phi+0.5*k1phi)),2)-cos((psi+0.5*k1psi)-(phi+0.5*k1phi))*(2*pow(ppsi+0.5*k1ppsi,2)+pow(pphi+0.5*k1pphi,2))+2*(pphi+0.5*k1pphi)*(ppsi+0.5*k1ppsi))/pow(2-pow(cos((psi+0.5*k1psi)-(phi+0.5*k1phi)),2),2)));
		k2ppsi = h*(-g*sin(psi+0.5*k1psi)+2*sin((psi+0.5*k1psi)-(phi+0.5*k1phi))*(((pphi+0.5*k1pphi)*(ppsi+0.5*k1ppsi)*pow(cos((psi+0.5*k1psi)-(phi+0.5*k1phi)),2)-cos((psi+0.5*k1psi)-(phi+0.5*k1phi))*(2*pow(ppsi+0.5*k1ppsi,2)+pow(pphi+0.5*k1pphi,2))+2*(pphi+0.5*k1pphi)*(ppsi+0.5*k1ppsi))/pow(2-pow(cos((psi+0.5*k1psi)-(phi+0.5*k1phi)),2),2)));
			
		//Calculamos las k3 usando las k2:
			
		k3phi = h*((pphi+k2pphi*0.5)-(ppsi+k2ppsi*0.5)*cos((ppsi+k2ppsi*0.5)-(pphi+k2pphi*0.5)))/(2-pow(cos((psi+k2psi*0.5)-(phi+k2phi*0.5)),2));
		k3psi = h*(-(pphi+0.5*k2pphi)*cos((psi+0.5*k2psi)-(phi+0.5*k2phi))+2*(ppsi+0.5*k2ppsi))/(2-pow(cos((psi+0.5*k2psi)-(phi+0.5*k2phi)),2));

		k3pphi = h*(-2*g*sin(phi+0.5*k2phi)+2*sin((phi+0.5*k2phi)-(psi+0.5*k2psi))*(((pphi+0.5*k2pphi)*(ppsi+0.5*k2ppsi)*pow(cos((psi+0.5*k2psi)-(phi+0.5*k2phi)),2)-cos((psi+0.5*k2psi)-(phi+0.5*k2phi))*(2*pow(ppsi+0.5*k2ppsi,2)+pow(pphi+0.5*k2pphi,2))+2*(pphi+0.5*k2pphi)*(ppsi+0.5*k2ppsi))/pow(2-pow(cos((psi+0.5*k2psi)-(phi+0.5*k2phi)),2),2)));
		k3ppsi = h*(-g*sin(psi+0.5*k2psi)+2*sin((psi+0.5*k2psi)-(phi+0.5*k2phi))*(((pphi+0.5*k2pphi)*(ppsi+0.5*k2ppsi)*pow(cos((psi+0.5*k2psi)-(phi+0.5*k2phi)),2)-cos((psi+0.5*k2psi)-(phi+0.5*k2phi))*(2*pow(ppsi+0.5*k2ppsi,2)+pow(pphi+0.5*k2pphi,2))+2*(pphi+0.5*k2pphi)*(ppsi+0.5*k2ppsi))/pow(2-pow(cos((psi+0.5*k2psi)-(phi+0.5*k2phi)),2),2)));
			
				
		//Calculamos las k4 usando las k3:
			
		k4phi = h*((pphi+k3pphi)-(ppsi+k3ppsi)*cos((ppsi+k3ppsi)-(pphi+k3pphi)))/(2-pow(cos((psi+k3psi)-(phi+k3phi)),2));
		k4psi = h*(-(pphi+k3pphi)*cos((psi+k3psi)-(phi+k3phi))+2*(ppsi+k3ppsi))/(2-pow(cos((psi+k3psi)-(phi+k3phi)),2));

		k4pphi = h*(-2*g*sin(phi+k3phi)+2*sin((phi+k3phi)-(psi+k3psi))*(((pphi+k3pphi)*(ppsi+k3ppsi)*pow(cos((psi+k3psi)-(phi+k3phi)),2)-cos((psi+k3psi)-(phi+k3phi))*(2*pow(ppsi+k3ppsi,2)+pow(pphi+k3pphi,2))+2*(pphi+k3pphi)*(ppsi+k3ppsi))/pow(2-pow(cos((psi+k3psi)-(phi+k3phi)),2),2)));
		k4ppsi = h*(-g*sin(psi+k3psi)+2*sin((psi+k3psi)-(phi+k3phi))*(((pphi+k3pphi)*(ppsi+k3ppsi)*pow(cos((psi+k3psi)-(phi+k3phi)),2)-cos((psi+k3psi)-(phi+k3phi))*(2*pow(ppsi+k3ppsi,2)+pow(pphi+k3pphi,2))+2*(pphi+k3pphi)*(ppsi+k3ppsi))/pow(2-pow(cos((psi+k3psi)-(phi+k3phi)),2),2)));
			

        //Calculamos las nuevas coordenadas:

		phi=phi+(k1phi+2*k2phi+2*k3phi+k4phi)/6.0; 
        psi=psi+(k1psi+2*k2psi+2*k3psi+k4psi)/6.0; 
		pphi=pphi+(k1pphi+2*k2pphi+2*k3pphi+k4pphi)/6.0; 
        ppsi=ppsi+(k1ppsi+2*k2ppsi+2*k3ppsi+k4ppsi)/6.0; 

        //Aumento el tiempo en un paso h
        t+=h;

		//Computamos las ecuaciones del movimiento:
		x1 = l1*sin(phi);
    	y1 = -l1*cos(phi);
    	x2 = x1 + l2*sin(psi);
    	y2 = y1 - l2*cos(psi);

        //if(j==30)
        //{
            fprintf (f1,"%lf\t",x1);
            fprintf(f1,",");
            fprintf (f1,"%lf\t",y1);
            fprintf(f1,"\n");
            fprintf (f1,"%lf\t",x2);
            fprintf(f1,",");
            fprintf (f1,"%lf\n",y2);
            fprintf(f1,"\n");  
            //j=0;             
        //}

        //else j=j+1;

    }
    fclose(f1);

}
    

