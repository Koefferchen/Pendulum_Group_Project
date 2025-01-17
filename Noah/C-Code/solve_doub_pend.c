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
    
    double g_grav2      = params[2];
    double mass_up      = params[3];
    double mass_down    = params[4];
    double length_up    = params[5];
    double length_down  = params[6];    
    double mass_ratio = mass_down/(mass_up+mass_down);
    double pre_factor = 1.0/( 1 - cos(y[0]-y[2])*cos(y[0]-y[2]) * mass_ratio );
    double delta_theta = y[0] - y[2];

    y_dot[0] = y[1];

    y_dot[1] =  pre_factor *
                (
                    g_grav2/length_up * (cos(delta_theta) *sin(y[2]) *mass_ratio - sin(y[0]))
                    - mass_ratio *sin(delta_theta) * ( length_down/length_up * y[3]*y[3] + cos(delta_theta) * y[1]*y[1] )
                );

    y_dot[2] = y[3];                       
    
    y_dot[3] =  pre_factor * 
                ( 
                    g_grav2/length_down * (cos(delta_theta) *sin(y[0]) - sin(y[2])) 
                    + sin(delta_theta) * ( length_up/length_down * y[1]*y[1] + mass_ratio * cos(delta_theta) * y[3]*y[3]  ) 
                ); 

}

int solve_doub_pend( double params[], double t_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[] ) 
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
    
        t = t + h; 
    } 
    
    return 0;
}

