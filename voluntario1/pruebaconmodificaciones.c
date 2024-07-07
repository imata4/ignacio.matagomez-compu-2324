#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define pi 3.141592
// --- Constant Parameters -----------------------------------------

const int PartN = 16; // Number of particles
const double sigma = 1; // Parameter sigma
const double ε = 1; // Energy factor
const double h = 0.002; // Time accuracy
const int L = 4; // Length of the simulation box's side

const double rc = 3; // Cutoff distance

// --- Introduction -------------------------------------------------

// Adrián Marín Boyero, Física UGR

/*

For this code, we will store our data in a 2-dim array, data[PartN][6],
in such a way that we'll have, for a particle i:

data[i][0] = x_position
data[i][1] = y_position
data[i][2] = x_velocity
data[i][3] = y_velocity
data[i][4] = x_acceleration
data[i][5] = y_acceleration

Also, we are assuming that beyond a distance rc between two particles, their interaction is zero. That way, we'll have
a modified force, Fm, which will be the calculated force, F, minus the force of the rc distance, Fc. That way we'll have a continuous function. 

Also, we are using periodic boundary conditions. 

*/

// --- Functions and main sequence ----------------------------------

// Initialization of the data matrix with empty values
void initialize_matrix(double data[PartN][6])
{
    for (int i = 0; i < PartN; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            data[i][j] = 0;
        }
    }
}

// Function that calculates the closest distance between 2 particles by making use of periodic boundary conditions (LxL box)
long double distance(double dx, double dy)
{
    if(abs(dx) > L/2)
    {
        dx = (L - abs(dx)) * (-dx)/abs(dx);
    }
    if(abs(dy) > L/2)
    {
        dy = (L - abs(dy)) * (-dy)/abs(dy);
    }

    return sqrt(dx*dx + dy*dy);
}

// Function that generates a random distribution of particles, separated by an arbitrary minimum distance
void random_initial_distribution(double data[PartN][6])
{
    double min_distance = 0.7;
    long double dx, dy;
    bool condition = false;
    double dr=0.5*L;

    double rnd1 = (double)rand() / (double)(RAND_MAX + 1);
    double rnd2 = (double)rand() / (double)(RAND_MAX + 1);
    data[0][0] = 0.5*L+2*(rnd1-0.5)*dr;
    data[0][1] = 0.5*L+2*(rnd2-0.5)*dr;

    for (int i = 1; i < PartN; i++)
    {
        condition = true;
        while (condition == true)
        {
            rnd1 = (double)rand() / (double)(RAND_MAX + 1);
            rnd2 = (double)rand() / (double)(RAND_MAX + 1);
            data[i][0] = 0.5*L+2*(rnd1-0.5)*dr;
            data[i][1] = 0.5*L+2*(rnd2-0.5)*dr;
            condition = false;

            for (int j = i - 1; j >= 0; j--)
            {
                dx = data[i][0] - data[j][0];
                dy = data[i][1] - data[j][1];

                if (distance(dx, dy) < min_distance)
                {
                    condition = true;
                    break;
                }
            }
        }
    }
    printf("Successful random distribution\n"); 
}

// Function that generates a aquare-net initial distribution 
void square_initial_distribution(double data[PartN][6])
{
    int j = 0;
    for (int i = 0; i < PartN; i++)
    {
        if ((i % 4 == 0) && (i != 0))
        {
            j++;
        }
        data[i][0] = L / 8.0 + L / 4.0 * (i % 4);
        data[i][1] = L / 8.0 + L / 4.0 * j;
    }
}

// Function that generates an hexagonal initial distribution 
void honeycomb_initial_distribution(double data[PartN][6])
{
    int j = 0;
    for (int i = 0; i < PartN; i++)
    {
        if ((i % 4 == 0) && (i != 0))
        {
            j++;
        }
        data[i][0] = L / 8.0 + L / 4.0 * (i % 4) - L / 8.0 * (j % 2);
        data[i][1] = L / 8.0 + L / 4.0 * j;
    }
}

// Function that gives each particle a velocity of a certain module, but in random directions
void random_initial_velocities_module(double data[PartN][6], double v_module)
{
    double theta;
    for (int i = 0; i < PartN; i++)
    {
        double rnd3 = (double)rand() / (double)(RAND_MAX + 1);
        theta=rnd3*2*pi;
        data[i][2] = cos(theta) * v_module;
        data[i][3] = sin(theta) * v_module; 

       /*Para cuando quiera vx positiva
       theta[i]=rnd3*pi+3*0.5*pi;
       data[i][2] = cos(theta) * v_module;;
       data[i][3] = 0;
       */

       /*Para apartado 4 ambas velocidades son nulas
       data[i][2] = 0;
       data[i][3] = 0;
       */
    }
}

// Function to print data to be read by the Python script for visualization
void print_data_for_python(double data[PartN][6], const char *nombre_archivo)
{
    FILE *f1;
    f1=fopen(nombre_archivo,"w");
    for ( int i = 0; i < PartN; i++)
    {
        for (int j = 0; i < 2; i++) 
        {  
            fprintf(f1, "%lf", data[i][j]); 
            if (j!=1)
                fprintf(f1,"," ); 
        }
        fprintf(f1, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones
    }
    fclose(f1);
    
     
}

// Function that calculates the forces between particles due to the Lennard-Jones potential
void force_calculations(double data[PartN][6],double Potential_Energy)
{
    Potential_Energy = 0;
    for(int k=0; k<PartN; k++)
    {
        for(int l=4; l<6; l++)
        {
            data[k][l] = 0;
        }
    }
    
    long double dx, dy, r, F, Fc, Fm, Vm, V, Vc;
    Fc=(48/pow(rc, 13)) - (24/pow(rc, 7));
    Vc=(4/pow(rc, 12)) - (4/pow(rc, 6));
    for(int i=0; i< (PartN - 1); i++)  // This way, we only need to calculate interactions between i and j only once
    {
        for(int j=i+1; j<PartN; j++)
        {
            dx = data[i][0] - data[j][0];
            dy = data[i][1] - data[j][1];
            r = distance(dx, dy);

            if( r <= rc)
            {
                F = (48/pow(r, 13)) - (24/pow(r, 7));
                Fm = F - Fc;

                data[i][4] += Fm * (dx/r);
                data[i][5] += Fm * (dy/r);
                data[j][4] -= Fm * (dx/r);
                data[j][5] -= Fm * (dy/r);
                V = (4/pow(r, 12)) - (4/pow(r, 6));
                Vm = V - Vc + r*Fc - rc*Fc;
                Potential_Energy += Vm;
            }
        }
    }
}

// Function resposible for updating the position of particles by making use of periodic boundary conditions
void update_positions(double data[PartN][6],double omega[PartN][2], double *momentum)
{
    momentum = 0;
    for(int i=0; i<PartN; i++)
    {
        for(int j=0; j<2; j++)
        {
            data[i][j] += h*omega[i][j];

            if(data[i][j] > L)
            {
                momentum += abs(data[i][j + 2]);
                data[i][j] = fmod(data[i][j], L);
            }
            if(data[i][j] < 0)
            {
                momentum += abs(data[i][j + 2]); //DEBERÍA DE PONER FABS????
                data[i][j] = L - fmod(fabs(data[i][j]), L);
            }
        }
    }
}

// Function that updates the velocity of each particle
void update_velocities(double data[PartN][6],double omega[PartN][2],double Kinetic_Energy)
{
    Kinetic_Energy = 0;
    for(int i=0; i<PartN; i++)
    {
        for(int j=0; j<2; j++)
        {
            data[i][j + 2] = omega[i][j] + (h/2.0) * data[i][j + 4];
        }
        Kinetic_Energy +=(data[i][2]*data[i][2] + data[i][3]*data[i][3]);
    }
    Kinetic_Energy=Kinetic_Energy*0.5;
}

// Function that makes use of the velvet method to update de positions, velocities and accelerations of each particle in the system
void velvet_integration(double data[PartN][6],double Potential_Energy,double Kinetic_Energy, double *momentum)
{
    // We define omega[PartN][2] as an auxiliary 2-dim array
    double omega[PartN][2];

    // We evaluate omega using omega = vel + (h/2)*accel
    for (int i = 0; i < PartN; i++) 
    {
        for (int j = 0; j < 2; j++)
        {
            omega[i][j] = data[i][j + 2] + (h / 2) * data[i][j + 4];
        }
    }

    update_positions(data, omega, momentum);
    force_calculations(data, Potential_Energy);
    update_velocities(data, omega, Kinetic_Energy);
}

// This function serves the purpose of treating data related to some parameters of the simulation, as well as printing them onto a .txt
void print_other_parameters(double data[PartN][6], double t, double Potential_E, double Kinetic_E, const char *nombre_archivo, double *momentum, double initial_r_0[2])
{
    double Average_PE, Average_KE, Total_E, sum_of_v = 0, fluctuation, dx, dy, distance_1p, distance_2p, distance_3p, dp, p = 0, sum_p = 0;
    double Fc=(48/pow(rc, 13)) - (24/pow(rc, 7));

    // We calculate the average PE, KE and Total Energy
    Average_PE = Potential_E/PartN;
    Average_KE = Kinetic_E/PartN;
    Total_E = Average_PE + Average_KE;

    // We calculate the distance between particle 1 and its initial position
    dx = data[0][0] - initial_r_0[0];
    dy = data[0][1] - initial_r_0[1];
    fluctuation = pow(distance(dx, dy), 2);

    // We calculate the distance between particles 1 and 2
    dx = data[0][0] - data[1][0];
    dy = data[0][1] - data[1][1];
    distance_1p = pow(distance(dx, dy), 2);

    dx = data[0][0] - data[5][0];
    dy = data[0][1] - data[5][1];
    distance_2p = pow(distance(dx, dy), 2);

    dx = data[0][0] - data[7][0];
    dy = data[0][1] - data[7][1];
    distance_3p = pow(distance(dx, dy), 2);

    // For an alternative definition of pressure, we need to calculate the average sum of distances between one particle and the rest
    for(int j=0; j<PartN; j++)
    {
        for(int k=0; k<PartN; k++)
        {
            if(j!=k)
            {
            dx = data[j][0] - data[k][0];
            dy = data[j][1] - data[k][1];
            dp = pow(distance(dx, dy), 2);
            
            p += dp * ((48/pow(dp, 13)) - (24/pow(dp, 7)) - Fc);
            }
        }
        sum_p = p/(1.0*PartN);
        p = 0;
    }

    // We calculate the mean velocity of the particles
    for(int i=0; i<PartN; i++)
    {
        sum_of_v += pow(data[i][2], 2) + pow(data[i][3], 2);
    }
    sum_of_v *= 1.0/PartN;

    // We print the data for the Python script to read
    //output << t << " " << Average_PE << " " << Average_KE << " " << Total_E << " " << sum_of_v << " " << 2.0 * momentum  << " " << fluctuation << " " << distance_1p <<  " " << distance_2p << " " << distance_3p << " " << sum_p << endl;
    print_data_for_python(data,"resultados.txt");
}

// This function registers the |V|, Vx, and Vy of the particles in a series of intervals, so as to be represented in bars by the Python script
void v_statistics(double data[PartN][6], double v_distribution[PartN][4])
{
    double v_module;
    double size = 0.2; // Size of the intervals

    for(int k=0; k<PartN; k++)
    {
        v_distribution[k][0] = (k-50)*size;
    }

    // We register the velocities of each particle (|V|, Vx and Vy)
    for(int i=0; i<PartN; i++)
    {
        v_module = sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
        
        for(int j=0; j<100; j++)
        {
            if( (v_module >= size*(j-50)) && (v_module < size*(j-49)))
            {
                v_distribution[j][1] += 1;
            }
            if( (data[i][2] >= size*(j-50)) && (data[i][2] < size*(j-49)))
            {
                v_distribution[j][2] += 1;
            }
            if( (data[i][3] >= size*(j-50)) && (data[i][3] < size*(j-49)))
            {
                v_distribution[j][3] += 1;
            }
        }
    }
}

// This function serves the purpose of treating data related to particles' velocities, as well as printing it onto a .txt
void print_v_stats(double v_distribution[PartN][4], const char *nombre_archivo)
{
    FILE *f1;
    f1=fopen(nombre_archivo,"w");
    double sum_v = 0, sum_vx = 0, sum_vy = 0;

    for(int i=0; i<PartN; i++)
    {
        sum_v += v_distribution[i][1];
        sum_vx += v_distribution[i][2];
        sum_vy += v_distribution[i][3];
    }
    
    // We normalize the total data collected by taken into account the total number of values and the size of the intervals (0.2)
    for(int j=0; j<PartN; j++)
    {
        v_distribution[j][1] *= 1.0/(sum_v * 0.2) ;
        v_distribution[j][2] *= 1.0/(sum_vx * 0.2);
        v_distribution[j][3] *= 1.0/(sum_vy * 0.2);
    }

    //cout << "Sum_v = " << sum_v << " | Sum_vx = " << sum_vx << " | Sum_vy = " << sum_vy << endl;

    // We print the data for the Python script to read
    /*for(int k=0; k<100; k++)
    {
        v_stats << v_distribution[k][0] << " " << v_distribution[k][1] << " " << v_distribution[k][2] << " " << v_distribution[k][3] << endl;
    }
    v_stats.close();*/

    for ( int i = 0; i < PartN; i++)
    {
        for (int j = 0; i < 4; i++) 
        {  
            fprintf(f1, "%lf", v_distribution[i][j]); 
            if (j<4)
                fprintf(f1,"," ); 
        }
        fprintf(f1, "\n"); //Añade nueva línea en blanco para separar datos entre iteraciones
    }
    fclose(f1);
    

}

// This function increases the velocity of the particles by a certain factor, at specific timestamps 
void increase_velocity(double data[PartN][6], double t)
{
    double factor = 1.1;

    bool condition = false;
    for(int j=0; j<20; j++)
    {
        if(condition == false)
        {
            if( abs(t - (20 + j*5)) < h/2.0 )
            //if( ((t >= 20) && (t < 20.002)) || ((t >= 30) && (t < 30.002)) || ((t >= 35) && (t < 35.002)) || ((t >= 45) && (t < 45.002)))
            {
                condition = true;
                for(int i=0; i<PartN; i++)
                {
                    data[i][2] *= factor;
                    data[i][3] *= factor;
                }

                //cout << "Velocity changed at t: " << t << " by a factor of " << factor << endl;
            }
        }
    }
}


// General algorithm for the iterations carried out for the simulations
void general_algorithm(double data[PartN][6], const char *archivo_parametros, const char *archivo_salida, int total_iterations, double initial_r_0[2])
{
    //FILE *parametros,*salida;

    double t = 0.0, v_distribution[PartN][4], *momentum; int show = 0;

    // We initilalize the v_distribution array (for future velocity analysis)
    for(int k=0; k<PartN; k++)
    {
        v_distribution[k][1] = 0;
        v_distribution[k][2] = 0;
        v_distribution[k][3] = 0;
    }
   

    double Potential_Energy, Kinetic_Energy;

    //salida=fopen(archivo_salida,"w");
    //parametros=fopen(archivo_parametros,"w");
    //output.open(output_name);
    //parameters_output.open("parameters.txt");

    // We print our initial data into the .txt
    print_data_for_python(data,archivo_salida);
    force_calculations(data, Potential_Energy);

    // We carry out the following functions a number of arbitrary iterations
    for (int steps = 0; steps < total_iterations; steps++)
    {
        if((t > 20) && (t < 50))
        {
            v_statistics(data, v_distribution);
        }
        velvet_integration(data, Potential_Energy, Kinetic_Energy, momentum);
        print_other_parameters(data, t, Potential_Energy, Kinetic_Energy,archivo_parametros, momentum, initial_r_0);

        // Interval time at which we want to study the velocities of the particles
        //if((t > 20) && (t < 50))

        // We increase the time value
        t += h;

        // We print the state of the system each x iterations (so as to avoid an excesively long .txt)
        if (show == 20)
        {
            print_data_for_python(data,archivo_salida);
            show = 0;
        }
        show++;

        //We update the velocities of the particles at certain timestamps to study some of the system's properties
        //increase_velocity(data, t);
    }

    //fclose(salida);
    //fclose(parametros);

 

    // We print the data for particles' velocites and the pair correlation functions
    
}

int main()
{
    
    srand(time(NULL));
    double data[PartN][6];
    initialize_matrix(data);

    // --- Possible initial configurations -----------------------

    random_initial_distribution(data);
    //square_initial_distribution(data);
    //honeycomb_initial_distribution(data);

    double v_module = 1; // Module of the velocities of each particle
    random_initial_velocities_module(data, v_module);

    double initial_r_0[2];
    initial_r_0[0] = data[0][0]; initial_r_0[1] = data[0][1];
    
    // --- Iterations and general algorithm -----------------------

    int total_iterations = 25000;
    general_algorithm(data,"parametros.txt","lennardjonesoutput.txt", total_iterations, initial_r_0); 

    return 0;
}