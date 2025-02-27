# include "header.h"


    // implementation of the differential equation
void derhs_doub_pend( int n_ODE, double t, double y[], double y_dot[], double params[] )
{
        // n_ODE = 4   
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

    // solves double pendulum for given (initial cond., numeric solver, given ODE) and returns solution arrays
int solve_doub_pend( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[],
                        void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ) 
{
    int     n_ODE = 4;                    // here: 4 DEs
    double  t_end = params[0];
    double  h     = params[1];
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[n_ODE];                     
        
        // first entry = length of array
    t_values[0]         = steps;            
    E_values[0]         = steps;
    theta1_sol[0]       = steps;      
    theta1_dot_sol[0]   = steps;
    theta2_sol[0]       = steps;
    theta2_dot_sol[0]   = steps;

        // Initial conditions at t_0
    y[0] = params[7];       // theta1                    
    y[1] = params[8];       // theta1_dot      
    y[2] = params[9];       // theta2
    y[3] = params[10];      // theta2_dot

    for(int j = 0; j < steps; j++)
    {
            // solve for one time step
        (num_solver)( n_ODE, h, t, y, &derhs_doub_pend, params);
        
            // save solution for each time step in arrays
        t_values[j+1]       = t;            
        theta1_sol[j+1]      = y[0];
        theta1_dot_sol[j+1]  = y[1];
        theta2_sol[j+1]      = y[2];
        theta2_dot_sol[j+1]  = y[3];
        calc_doub_energy( y, params, &E_values[j+1] );
    
        t = t + h; 
    } 
    
    return 0;
}

    // calculates current energy from current position in phase space for the double pendulum
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

    // find the root of a the energy function (E - E_0) with respect to theta_2 
double find_doub_theta2_max ( double theta2_min, double theta2_max, const double tol, double params[] )  
{
        // use altered bisection method 
    double intersection = (theta2_max + theta2_min)/2 ;         
    double theta1       = params[7];
    double theta1_dot   = params[8];
    double theta2_dot   = 0.0;
    double E_0          = params[11];
    double state_min[]  = {theta1, theta1_dot, theta2_min, theta2_dot};
    double state_int[]  = {theta1, theta1_dot, intersection, theta2_dot};
    double E_min;
    double E_int;
    calc_doub_energy(state_min, params, &E_min);
    calc_doub_energy(state_int, params, &E_int);

    if ( (theta2_max - theta2_min)/2 < tol )                   
    {
        return intersection;
    }
    if ( (E_min - E_0) * (E_int - E_0) <= 0 )   
    {
        return find_doub_theta2_max( theta2_min, intersection, tol, params );
    }
    else                                       
    {
        return find_doub_theta2_max( intersection, theta2_max, tol, params );
    }

}

    // calculates the theta2_dot_0 from given (theta1_0, theta1_dot_0, theta2_0, E) from energy conservation
int calc_doub_thetadot2_0( double params[], double parity )
{
        // unpack parameters
    double g_grav       = params[2];
    double mass_1       = params[3];
    double mass_2       = params[4];
    double length_1     = params[5];
    double length_2     = params[6]; 
    double theta1_0     = params[7];
    double theta1_dot_0 = params[8];
    double theta2_0     = params[9];
    double* theta2_dot_0= &params[10];
    double E_value      = params[11];

    double alpha    = 0.5 * mass_2 * length_2*length_2;
    double p        = mass_2 * theta1_dot_0 * length_1*length_2 * cos(theta1_0-theta2_0);
    double q        = 0.5 * (mass_1+mass_2) * length_1*length_1 * theta1_dot_0*theta1_dot_0 
                        +(mass_1+mass_2) * g_grav * length_1 * (1.0 - cos(theta1_0)) 
                        +mass_2 * g_grav * length_2 * (1.0 - cos(theta2_0)) 
                        -E_value;
    if( p == 0.0)
    {
        *theta2_dot_0 = parity* sqrt(-q/alpha);
    } else {
        *theta2_dot_0   = p/(2.0*alpha) * ( parity * sqrt(1.0 - 4.0*alpha*q /p /p) - 1.0 );
    }

    return 0;
}

    // simulates a poincare section of the double pendulum for given initial conditions
int solve_doub_poincare( double params[], double t_values[], double theta2_sol[], double theta2_dot_sol[],
                        void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ) 
{
    int     n_ODE       = 4;                    // here: 4 DEs
    double  t_end       = params[0];
    double  h           = params[1];
    int     steps       = (int)(t_end/h);
    int     max_points  = (int)params[13];
    double  t           = 0;
    
    double  y[n_ODE];                    
    double last_theta1;
    double last_theta1_dot;
    double last_theta2;      
    double last_theta2_dot;

    double new_theta1;
    double new_theta1_dot;
    double new_theta2;      
    double new_theta2_dot;

    y[0] = params[7];                       
    y[1] = params[8];                     
    y[2] = params[9];
    y[3] = params[10];

    int i = 0;
    for(int j = 0; j < steps; j++)
    {
        last_theta1     = modulus_s(y[0]+M_PI, 2*M_PI) - M_PI;
        last_theta1_dot = y[1];
        last_theta2     = modulus_s(y[2]+M_PI, 2*M_PI) - M_PI;
        last_theta2_dot = y[3];

        (num_solver)( n_ODE, h, t, y, &derhs_doub_pend, params);

        new_theta1     = modulus_s(y[0]+M_PI, 2*M_PI) - M_PI;
        new_theta1_dot = y[1];
        new_theta2     = modulus_s(y[2]+M_PI, 2*M_PI) - M_PI;
        new_theta2_dot = y[3];

        if( (last_theta1 < 0.0) && (new_theta1 > 0.0) && ( (last_theta1_dot+new_theta1_dot)/2.0 > 0.0) && ( i < max_points) )
        {
            t_values[i+1]       = t + h/2.0;           
            theta2_sol[i+1]     = (last_theta2+new_theta2)/2.0;
            theta2_dot_sol[i+1] = (last_theta2_dot+new_theta2_dot)/2.0;
            i++;
        }
    
        t = t + h; 
        if( i > max_points){ break; printf("--- Double Poincare: data full (%d/%d) \n", max_points, max_points); }
    } 
        // first entry = length of array
    t_values[0]         = i;            
    theta2_sol[0]       = i;
    theta2_dot_sol[0]   = i;

    return 0;
}

    // finds "t" when theta1/2 flips over for the 1st time for each combination of initial conditions (theta1_0, theta2_0)
int solve_doub_flip( double params[], double theta1_0_list[], double theta2_0_list[], double **theta1_sol_matrix, double **theta2_sol_matrix,
                        void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ) 
{

        // unpacking params
    int     n_ODE       = 4;                    
    double  t_max       = params[0];
    double  h           = params[1];
    int     t_steps     = (int)(t_max/h);
    
        // initialisation
    double t;
    double y[n_ODE];                    
    double last_theta1;
    double last_theta1_dot;
    double last_theta2;      
    double last_theta2_dot;

    double new_theta1;
    double new_theta1_dot;
    double new_theta2;      
    double new_theta2_dot;

    double now_theta1_dot;
    double now_theta2_dot;
    int    continue_loop[2];

        // repeat for each set of initial conditions
    for( int x = 0; x < theta1_0_list[0]; x++)
    {
        for( int z = 0; z < theta2_0_list[0]; z++)
        {
            y[0] = theta1_0_list[x+1];                       
            y[1] = params[9];                     
            y[2] = theta2_0_list[z+1];
            y[3] = params[12];

                // reset initial values for next simulation
            continue_loop[0] = 1;   continue_loop[1] = 1;    t = 0.0;

            for(int j = 0; j < t_steps; j++)
            {
                    // state at t_{n}
                last_theta1     = modulus_s(y[0]+M_PI, 2*M_PI) - M_PI;
                last_theta1_dot = y[1];
                last_theta2     = modulus_s(y[2]+M_PI, 2*M_PI) - M_PI;
                last_theta2_dot = y[3];

                    // calculate state at t_{n+1}
                (num_solver)( n_ODE, h, t, y, &derhs_doub_pend, params);

                    // state at t_{n+1}
                new_theta1     = modulus_s(y[0]+M_PI, 2*M_PI) - M_PI;
                new_theta1_dot = y[1];
                new_theta2     = modulus_s(y[2]+M_PI, 2*M_PI) - M_PI;
                new_theta2_dot = y[3];

                    // average velocity at flip over
                now_theta1_dot = (new_theta1_dot+last_theta1_dot)/2.0;
                now_theta2_dot = (new_theta2_dot+last_theta2_dot)/2.0;

                    // if theta1 flips over for the 1st time, save that time 
                if( (now_theta1_dot*last_theta1 > 0.0) && (now_theta1_dot*new_theta1 < 0.0) && (continue_loop[0])  )
                {
                    theta1_sol_matrix[z+1][x] = t + h/2.0;
                    continue_loop[0] = 0;  
                }
                    // if theta2 flips over for the 1st time, save that time 
                if( (now_theta2_dot*last_theta2 > 0.0) && (now_theta2_dot*new_theta2 < 0.0) && (continue_loop[1]) )
                {
                    theta2_sol_matrix[z+1][x] = t + h/2.0;
                    continue_loop[1] = 0;
                }
            
                t = t + h;
                
            } 
        }
    }
    merge_arrays(14, params, theta1_0_list[0], theta1_sol_matrix[0] );
    merge_arrays(14, params, theta1_0_list[0], theta2_sol_matrix[0] );

    return 0;
}