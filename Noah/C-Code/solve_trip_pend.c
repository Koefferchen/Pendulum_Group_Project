# include "header.h"



    // implementation of the differential equation
void derhs_trip_pend( int nDifEqu, double t, double y[], double y_dot[], double params[] )
{
        // nDifEqu = 6   
        
        // "y[0]"                   is  theta1
        // "y_dot[0]" and "y[1]"    are theta1_dot
        // "y_dot[1]"               is  theta1_dot_dot
        
        // "y[2]"                   is  theta2
        // "y_dot[2]" and "y[3]"    are theta2_dot
        // "y_dot[3]"               is  theta2_dot_dot 

        // "y[4]"                   is  theta3
        // "y_dot[4]" and "y[5]"    are theta3_dot
        // "y_dot[5]"               is  theta3_dot_dot  
    double g_grav        = params[2];
    double mass_1        = params[3];
    double mass_2        = params[4];
    double mass_3        = params[5];
    double length_1      = params[6];
    double length_2      = params[7];
    double length_3      = params[8];

    double theta_12 = y[0] - y[2];
    double theta_13 = y[0] - y[4];
    double theta_23 = y[2] - y[4];
    double mass_123 = mass_1 + mass_2 + mass_3;
    double mass_23  = mass_2 + mass_3; 
    double a_1 = mass_123   * length_1;
    double a_2 = mass_23    * length_2 * cos(theta_12);
    double a_3 = mass_3     * length_3 * cos(theta_13);
    double b_1 = mass_23    * length_1 * cos(theta_12);
    double b_2 = mass_23    * length_2;
    double b_3 = mass_3     * length_3 * cos(theta_23);
    double c_1 = mass_3     * length_1 * cos(theta_13);
    double c_2 = mass_3     * length_2 * cos(theta_23);
    double c_3 = mass_3     * length_3;
    double g_1 = mass_23 * length_2 * y[3]*y[3] * sin(theta_12) + mass_3 * length_3 * y[5]*y[5] * sin(theta_13) + mass_123 * g_grav * sin(y[0]);
    double g_2 = - mass_23 * length_1 * y[1]*y[1] * sin(theta_12) + mass_3 * length_3 * y[5]*y[5] * sin(theta_23) + mass_23 * g_grav * sin(y[2]);
    double g_3 = - mass_3 * length_1 * y[1]*y[1] * sin(theta_13) - mass_3 * length_2 * y[3]*y[3] * sin(theta_23) + mass_3 * g_grav * sin(y[4]);
    double e_2 = 1.0 - b_1 * a_2 / (a_1 * b_2);
    double e_3 = 1.0 - c_1 * a_3 / (c_3 * a_1);
    double k_12 = g_1 * b_1 / (a_1 * b_2) - g_2/b_2;
    double k_13 = g_1 * c_1 / (a_1 * c_3) - g_3/c_3;
    double v_2 = a_2 * c_1 / (a_1 * c_3) - c_2/c_3;
    double v_3 = a_3 * b_1 / (a_1 * b_2) - b_3/b_2;

    y_dot[5] = ( k_13/e_3 + v_2 * k_12 / (e_2 * e_3) )*( 1.0 / (1.0 - v_2*v_3 / (e_2*e_3)) );

    y_dot[3] =  k_12/e_2 + v_3/e_2 * y_dot[5]; 

    y_dot[1] =  -1.0/a_1 * ( a_2*y_dot[3] + a_3*y_dot[5] );

    y_dot[0] = y[1];

    y_dot[2] = y[3];                       

    y_dot[4] = y[5];


}

int solve_trip_pend( double params[], double t_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[], double theta3_sol[], double theta3_dot_sol[] ) 
{
    int     nDifEqu = 6;                    // here: 6 DEs
    double  t_end   = params[0];
    double  h       = params[1];
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
    theta3_sol[0]       = steps;
    theta3_dot_sol[0]   = steps;

    y[0] = params[9];                         // Initial conditions at t_0
    y[1] = params[10];                     
    y[2] = params[11];
    y[3] = params[12];
    y[4] = params[13];
    y[5] = params[14];

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_trip_pend, params);
        
        t_values[j+1]       = t;            // save solution for each time step in arrays
        theta1_sol[j+1]      = y[0];
        theta1_dot_sol[j+1]  = y[1];
        theta2_sol[j+1]      = y[2];
        theta2_dot_sol[j+1]  = y[3];
        theta3_sol[j+1]      = y[4];
        theta3_dot_sol[j+1]  = y[5];
    
        t = t + h; 
    } 
    
    return 0;
}