# include "header.h"



    // implementation of the differential equation
void derhs_trip_pend( int n_ODE, double t, double y[], double y_dot[], double params[] )
{
        // n_ODE = 6   
        
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

    y_dot[5] = ( k_13/e_3 + v_2 * k_12 / (e_2 * e_3) ) / (1.0 - v_2*v_3 / (e_2*e_3)) ;

    y_dot[3] =  k_12/e_2 + v_3/e_2 * y_dot[5]; 

    y_dot[1] =  -1.0/a_1 * ( a_2*y_dot[3] + a_3*y_dot[5] + g_1 );

    y_dot[0] = y[1];

    y_dot[2] = y[3];                       

    y_dot[4] = y[5];


}

int solve_trip_pend( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[], double theta3_sol[], double theta3_dot_sol[],
                        void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ) 
{
    int     n_ODE   = 6;                    
    double  t_end   = params[0];
    double  h       = params[1];
    int     steps   = (int)(t_end/h);
    double  t       = 0;
    double  y[n_ODE];                     

         // first entry = length of array
    t_values[0]         = steps;           
    E_values[0]         = steps;
    theta1_sol[0]       = steps;      
    theta1_dot_sol[0]   = steps;
    theta2_sol[0]       = steps;
    theta2_dot_sol[0]   = steps;
    theta3_sol[0]       = steps;
    theta3_dot_sol[0]   = steps;

        // Initial conditions at t_0
    y[0] = params[9];                         
    y[1] = params[10];                     
    y[2] = params[11];
    y[3] = params[12];
    y[4] = params[13];
    y[5] = params[14];

    for(int j = 0; j < steps; j++)
    {
            // solve set of ODE's for one time step 
        (num_solver)( n_ODE, h, t, y, &derhs_trip_pend, params);
            // save solution for each time step in arrays
        t_values[j+1]       = t;           
        theta1_sol[j+1]      = y[0];
        theta1_dot_sol[j+1]  = y[1];
        theta2_sol[j+1]      = y[2];
        theta2_dot_sol[j+1]  = y[3];
        theta3_sol[j+1]      = y[4];
        theta3_dot_sol[j+1]  = y[5];
        calc_trip_energy( y, params, &E_values[j+1] );
        t = t + h; 
    } 
    
    return 0;
}

int calc_trip_energy( double y[], double params[], double *E_value )
{
        // unpacking parameters
    double g_grav        = params[2];
    double mass_1        = params[3];
    double mass_2        = params[4];
    double mass_3        = params[5];
    double length_1      = params[6];
    double length_2      = params[7];
    double length_3      = params[8];
    double theta1        = y[0];
    double theta1_dot    = y[1];
    double theta2        = y[2];
    double theta2_dot    = y[3];
    double theta3        = y[4];
    double theta3_dot    = y[5];
    double mass_123 = mass_1 + mass_2 + mass_3;
    double mass_23  = mass_2 + mass_3; 

    *E_value =  0.5 * mass_123 * length_1*length_1 * theta1_dot*theta1_dot +
                0.5 * mass_23 * length_2*length_2 * theta2_dot*theta2_dot + 
                0.5 * mass_3 * length_3*length_3 * theta3_dot*theta3_dot + 
                mass_23 * length_1*length_2 * theta1_dot*theta2_dot * cos(theta1-theta2) +
                mass_3 * length_1*length_3 * theta1_dot*theta3_dot * cos(theta1-theta3) + 
                mass_3 * length_2*length_3 * theta2_dot*theta3_dot * cos(theta2-theta3) +
                mass_123 * g_grav * length_1 * (1 - cos(theta1)) +
                mass_23 * g_grav * length_2 * (1 - cos(theta2)) + 
                mass_3 * g_grav * length_3 * (1 - cos(theta3));
 
    return 0;
}