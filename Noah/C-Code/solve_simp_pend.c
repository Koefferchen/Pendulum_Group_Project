
# include "header.h"



    // implementation of the differential equation for simple pendulum
void derhs_simp_pend( int n_ODE, double t, double y[], double y_dot[], double params[] )
{        
        // n_ODE = 2
        // "y[0]" is theta
        // "y_dot[0]" and "y[1]" are theta_dot
        // "y_dot[2]" is theta_dot_dot             
    double g_grav = params[2];
    double length = params[3];
    double omega  = pow( g_grav/length, 0.5);

    y_dot[0] = y[1] ;                           // 1. DE: theta(t)_dot      = theta_dot(t)
    y_dot[1] = - pow(omega, 2)  * sin(y[0]) ;   // 2. DE: theta_dot(t)_dot  = -omega**2 * theta(t)    
}

    // implementation of the approximated equation for simple pendulum: sin(x) = x
void derhs_analyt_pend( int n_ODE, double t, double y[], double y_dot[], double params[] )
{        
        // n_ODE = 2
        // "y[0]" is theta
        // "y_dot[0]" and "y[1]" are theta_dot
        // "y_dot[2]" is theta_dot_dot             
    double g_grav = params[2];
    double length = params[3];
    double omega  = pow( g_grav/length, 0.5);

    y_dot[0] = y[1] ;                           // 1. DE: theta(t)_dot      = theta_dot(t)
    y_dot[1] = - pow(omega, 2)  * y[0] ;        // 2. DE: theta_dot(t)_dot  = -omega**2 * theta(t)    
}
   
    // compute initial angular velocity to make pendulum stand upright
double stand_up_theta_dot( double length, double g_grav, double theta_0 )
{
    return sqrt(g_grav / length) * pow( 2 * (1+cos(theta_0)), 0.5 );
}




    // fully solve the DE of a simple pendulum numerically
int solve_simp_pend( double params[], double t_values[], double theta_sol[], double theta_dot_sol[], 
                        void (*derhs)( int, double, double[], double[], double[] ), 
                        void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ) 
{
    int     n_ODE = 2;                   // here: 2 DEs
    double  t_end   = params[0];
    double  h       = params[1];
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[n_ODE];                     // Initialisation 

    t_values[0]         = steps;            // first entry = length of array
    theta_sol[0]        = steps;      
    theta_dot_sol[0]    = steps;

    y[0] = params[4];                       // Initial conditions at t_0
    y[1] = params[5];                       // stand_up_theta_dot(omega, theta_0);
    
    for(int j = 0; j < steps; j++)
    {
        (num_solver)( n_ODE, h, t, y, derhs, params);
        
        t_values[j+1]       = t;            // save solution for each time step in arrays
        theta_sol[j+1]      = y[0];
        theta_dot_sol[j+1]  = y[1];

        t = t + h; 
    } 
    
    return 0;
}

int solve_analyt_pend( double params[], double theta_analyt_sol[] )
{
    double t_end        = params[0];
    double h            = params[1];
    double g_grav       = params[2];
    double length       = params[3];
    double omega        = sqrt( g_grav/length );
    double theta_0      = params[4];
    double theta_dot_0  = params[5];
    int    steps        = (int)(t_end/h);
    theta_analyt_sol[0] = (double)steps;

    for(int j=0; j < steps; j++ )
    {
        theta_analyt_sol[j+1] = theta_0 * cos( omega * (j*h) ) + theta_dot_0/omega * sin(omega *(j*h) );
    }

    return 0;
}

    // calculates the maximum deviation between elements of 2 arrays
double num_max_deviation( double theta2_num[], double theta2_ana[] )
{
    int length      = (int)theta2_ana[0];
    double max_dev  = 0.0;
    double temp_dev;

    for(int i = 1; i<length+1; i++) 
    {
        temp_dev = fmin( fabs(theta2_ana[i] - theta2_num[i]), fabs(2*M_PI - fabs(theta2_ana[i] - theta2_num[i])) );
        if( temp_dev > max_dev )
        {
            max_dev = temp_dev;
        }
    }

    return max_dev;
}