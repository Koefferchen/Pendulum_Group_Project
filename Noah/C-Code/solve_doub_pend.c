# include "header.h"




    // implementation of the differential equation
void derhs_doub_pend( int nDifEqu, double t, double y[], double y_dot[], double params[] )
{
        // nDifEqu = 4   
        // "y[0]"                   is  theta1
        // "y_dot[0]" and "y[1]"    are theta1_dot
        // "y_dot[1]"               is  theta1_dot_dot
        
        // "y[2]"                   is  theta2
        // "y_dot[2]" and "y[3]"    are theta2_dot
        // "y_dot[3]"               is  theta2_dot_dot  
    
    double g_grav       = params[2];
    double mass_1       = params[3];
    double mass_2       = params[4];
    double length_1     = params[5];
    double length_2     = params[6];    
    double mass_ratio   = mass_2/(mass_1+mass_2);
    double pre_factor   = 1.0/( 1 - cos(y[0]-y[2])*cos(y[0]-y[2]) * mass_ratio );
    double delta_theta  = y[0] - y[2];

    y_dot[0] = y[1];

    y_dot[1] =  pre_factor *
                (
                    g_grav/length_1 * (cos(delta_theta) *sin(y[2]) *mass_ratio - sin(y[0]))
                    - mass_ratio *sin(delta_theta) * ( length_2/length_1 * y[3]*y[3] + cos(delta_theta) * y[1]*y[1] )
                );

    y_dot[2] = y[3];                       
    
    y_dot[3] =  pre_factor * 
                ( 
                    g_grav/length_2 * (cos(delta_theta) *sin(y[0]) - sin(y[2])) 
                    + sin(delta_theta) * ( length_1/length_2 * y[1]*y[1] + mass_ratio * cos(delta_theta) * y[3]*y[3]  ) 
                ); 

}

int solve_doub_pend( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[] ) 
{
    int     nDifEqu = 4;                    // here: 2 DEs
    double  t_end = params[0];
    double  h     = params[1];
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];                     // Initialisation 
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];
    t_values[0]         = steps;            // first entry = length of array
    E_values[0]         = steps;
    theta1_sol[0]       = steps;      
    theta1_dot_sol[0]   = steps;
    theta2_sol[0]       = steps;
    theta2_dot_sol[0]   = steps;

    y[0] = params[7];                         // Initial conditions at t_0
    y[1] = params[8];                     
    y[2] = params[9];
    y[3] = params[10];

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_doub_pend, params);
        
        t_values[j+1]       = t;            // save solution for each time step in arrays
        theta1_sol[j+1]      = y[0];
        theta1_dot_sol[j+1]  = y[1];
        theta2_sol[j+1]      = y[2];
        theta2_dot_sol[j+1]  = y[3];
        calc_doub_energy( y, params, &E_values[j+1] );
    
        t = t + h; 
    } 
    
    return 0;
}


int calc_doub_energy( double y[], double params[], double *E_value )
{
        // unpacking parameters
    double g_grav       = params[2];
    double mass_1       = params[3];
    double mass_2       = params[4];
    double length_1     = params[5];
    double length_2     = params[6];    
    double theta1        = y[0];
    double theta1_dot    = y[1];
    double theta2        = y[2];
    double theta2_dot    = y[3];

    *E_value =  0.5 * (mass_1+mass_2) * length_1*length_1 * theta1_dot*theta1_dot + 
                0.5 * mass_2 * length_2*length_2 * theta2_dot*theta2_dot + 
                mass_2 * length_1*length_2 * theta1_dot*theta2_dot * cos(theta1-theta2) +
                (mass_1+mass_2) * g_grav * length_1 * (1 - cos(theta1)) +
                mass_2 * g_grav * length_2 * (1 - cos(theta2));

    return 0;
}