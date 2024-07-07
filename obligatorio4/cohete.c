#include <stdio.h>
#include <math.h>

#define PI 3.14159
#define w 2.6617*pow(10,-6)

int main (void)
{
    double h;
	h = 30;

	double r,phi,pr,pphi; //Coordenadas
	double k1r,k1phi,k1pr,k1pphi;
	double k2r,k2phi,k2pr,k2pphi;
	double k3r,k3phi,k3pr,k3pphi;
	double k4r,k4phi,k4pr,k4pphi; //ki para todas las coordenadas
	double aux,aux1;
	int i,j;
	double delta,mu,rprima; // constantes
	double t; //tiempo
	
	double xTierra,yTierra; //Posciones x e y de la Tierra,cohete y Luna.  
	double xcohete,ycohete;
	double xLuna, yLuna;
	double v0,theta;

    FILE *f1;
	
	f1 = fopen("cohete.txt","w");

    //Valores iniciales
	t =0;
	r = 0.0166; //posición inicial del cohete rescalada. Lo situamos a una distancia de 1R_T
	phi =PI/2.0; //lanzamos el cohete desde el polo norte

	v0= 0.000029; //velocidad de escape del cohete reescalada
	theta = 56.3*PI/180.0; // Lo lanzo desde granada
	
	pr= v0*cos(theta-phi);
	pphi = r*v0*sin(theta-phi);
	
    //Coloco inicialmente a la Luna y el cohete en el eje y
    //La luna se encuentra a una distancia tierra luna y el cohete sale desde la superficie de la Tierra
	xLuna = 1;
	yLuna = 0;
	xcohete = r;
	ycohete = 0;
	
	fprintf (f1,"%lf\t",xTierra);
    fprintf(f1,",");
	fprintf (f1,"%lf\t",yTierra);
    fprintf(f1,"\n");
	fprintf (f1,"%lf\t",xcohete);
    fprintf(f1,",");
	fprintf (f1,"%lf\t",ycohete);
    fprintf(f1,"\n");
	fprintf (f1,"%lf\t",xLuna);
    fprintf(f1,",");
	fprintf (f1,"%lf\n",yLuna);
    fprintf(f1,"\n");

    //Cálculo de Constantes
	
	delta = 7.014744145*pow(10,-12); 
	mu	= 0.0123; //cociente entre la masa de la Luna y la Tierra

    j=0; //Esta variable la uso para no imprimir en cada iteración sino cada las que yo quiera, porque sino tarda mucho la animación
    for ( i = 0; i < 70000; i++)
    {
       //Calculamos las k1
		k1r = h*pr;
		k1phi = h*pphi/(r*r);
	
		rprima = pow(1+r*r-2*r*cos(phi-2.6617*pow(10,-6)*t),0.5);
		k1pr = h*((pphi*pphi)/pow(r,3)-delta*(1/(r*r)+(mu/pow(rprima,3))*(r-cos(phi-2.6617*pow(10,-6)*t))));
	
		k1pphi = h*(-1.0)*delta*mu*(r/pow(rprima,3))*sin(phi-2.6617*pow(10,-6)*t);
		
		//Calculamos las k2 usando las k1:
		
		k2r = h*(pr+k1pr*0.5);
		k2phi = h*(pphi+k1pphi*0.5)/((r+k1r*0.5)*(r+k1r*0.5));

		rprima = pow(1+(r+k1r*0.5)*(r+k1r*0.5)-2*(r+k1r*0.5)*cos((phi+k1phi*0.5)-2.6617*pow(10,-6)*(t+h*0.5)),0.5);
			
		aux = (pphi+k1pphi*0.5)*(pphi+k1pphi*0.5)/pow((r+k1r*0.5),3) ;
		aux1 = 1.0/((r+k1r*0.5)*(r+k1r*0.5));
		
		k2pr = h*(aux-delta*(aux1+(mu/pow(rprima,3))*((r+k1r*0.5)-cos((phi+k1phi*0.5)-2.6617*pow(10,-6)*(t+h*0.5)))));
	
		k2pphi = h*(-1.0)*delta*mu*((r+k1r*0.5)/pow(rprima,3))*sin((phi+k1phi*0.5)-2.6617*pow(10,-6)*(t+h*0.5));
			
		//Calculamos las k3 usando las k2:
			
		k3r = h*(pr+k2pr*0.5);
		k3phi = h*(pphi+k2pphi*0.5)/((r+k2r*0.5)*(r+k2r*0.5));
	
		rprima = pow(1+(r+k2r*0.5)*(r+k2r*0.5)-2*(r+k2r*0.5)*cos((phi+k2phi*0.5)-w*(t+h*0.5)),0.5);
		
		aux = (pphi+k2pphi*0.5)*(pphi+k2pphi*0.5)/pow((r+k2r*0.5),3) ;
		aux1 = 1/((r+k2r*0.5)*(r+k2r*0.5));
		k3pr = h*(aux-delta*(aux1+(mu/pow(rprima,3))*((r+k2r*0.5)-cos((phi+k2phi*0.5)-w*(t+h*0.5)))));
	
		k3pphi = h*(-1.0)*delta*mu*((r+k2r*0.5)/pow(rprima,3))*sin((phi+k2phi*0.5)-w*(t+h*0.5));
			
		//Calculamos las k4 usando las k3:
			
		k4r = h*(pr+k3pr);
		k4phi = h*(pphi+k3pphi)/((r+k3r)*(r+k3r));
	
		rprima = pow(1+(r+k3r)*(r+k3r)-2*(r+k3r)*cos((phi+k3phi)-w*(t+h)),0.5);
			
		aux = (pphi+k3pphi)*(pphi+k3pphi)/pow((r+k3r),3) ;
		aux1 = 1/((r+k3r)*(r+k3r));
		k4pr = h*(aux-delta*(aux1+(mu/pow(rprima,3))*((r+k3r)-cos((phi+k3phi)-w*(t+h)))));
	
		k4pphi = h*(-1.0)*delta*mu*((r+k3r)/pow(rprima,3))*sin((phi+k3phi)-w*(t+h));
			

        //Calculamos las nuevas coordenadas:

        r=r+(k1r+2*k2r+2*k3r+k4r)/6.0;
        phi=phi+(k1phi+2*k2phi+2*k3phi+k4phi)/6.0; 
        pr=pr+(k1pr+2*k2pr+2*k3pr+k4pr)/6.0;
        pphi=pphi+(k1pphi+2*k2pphi+2*k3pphi+k4pphi)/6.0; 

        //Aumento el tiempo en un paso h
        t+=h;

        //Computamos las ecuaciones del movimiento:
		xTierra = 0;
		yTierra = 0;
		xcohete = r*cos(phi);
		ycohete = r*sin(phi);
		xLuna = cos(w*t);
		yLuna = sin(w*t);

        if(j==30)
        {
            fprintf (f1,"%lf\t",xTierra);
            fprintf(f1,",");
            fprintf (f1,"%lf\t",yTierra);
            fprintf(f1,"\n");
            fprintf (f1,"%lf\t",xcohete);
            fprintf(f1,",");
            fprintf (f1,"%lf\t",ycohete);
            fprintf(f1,"\n");
            fprintf (f1,"%lf\t",xLuna);
            fprintf(f1,",");
            fprintf (f1,"%lf\n",yLuna);
            fprintf(f1,"\n");  
            j=0;             
        }

        else j=j+1;

    }
    fclose(f1);
	return 0;
}