
#include "header.h"

int main(void)
{
    simple_pendulum();

    double_pendulum();
    double_poincare();

    triple_pendulum();
    triple_chaos();
    triple_compare();

    test_num_solvers();

    erase_last_line();
    printf("Simulation Data generated\n");
    return 0;
}

// -----------------------------------------------------------------------------

int simple_pendulum(void)
{
    erase_last_line();
    printf("--- Simple Pendulum \n");

    double t_end        = 10.0;             // simulation time [seconds]            
    double h            = 0.0001;            // step size [steps per second]
    double g_grav       = 9.81;
    double length       = 0.7;       
    double theta_0      = 0.3 * M_PI; 
    double theta_dot_0  = stand_up_theta_dot(length, g_grav, theta_0); 

    int    steps = (int)(t_end/h);      // Initialisation
    double params[] = {t_end, h, g_grav, length, theta_0, theta_dot_0};
    double params_extended  [steps+1];
    double t_values         [steps+1];
    double theta            [steps+1];              
    double theta_dot        [steps+1];          
    double theta_analy      [steps+1];
    double null             [steps+1]; zeros(null, steps+1);

    solve_simp_pend( params, t_values, theta, theta_dot, &derhs_simp_pend, &RuKu_6 );       
    solve_analyt_pend( params, theta_analy );                                  
    
    merge_arrays(6, params, steps, params_extended);
    save_numb_list7(t_values, theta, theta_dot, theta_analy, params_extended, null, null, "../data/data_simp_pend.txt" );   // save data in .txt

    return 0;
}


// -----------------------------------------------------------------------------


int double_pendulum(void)
{
    erase_last_line();
    printf("--- Double Pendulum \n");
    
    double t_end        = 10.0;             // simulation time [seconds]            
    double h            = 0.001;            // step size [steps per second]
    double g_grav       = 9.80665;
    double mass_1       = 1.0;
    double mass_2       = 2.0;
    double length_1     = 1.0;
    double length_2     = 0.7;
    double theta1_0     = 1.0*M_PI;
    double theta1_dot_0 = 0.0;
    double theta2_0     = 0.5*M_PI;
    double theta2_dot_0 = 0.0;
    
    int    steps = (int)(t_end/h);      // Initialisation
    double params[] = {t_end, h, g_grav, mass_1, mass_2, length_1, length_2, theta1_0, theta1_dot_0, theta2_0, theta2_dot_0};
    double params_extended  [steps+1];
    double t_values         [steps+1];
    double E_values         [steps+1];
    double theta1           [steps+1];
    double theta1_dot       [steps+1];
    double theta2           [steps+1];
    double theta2_dot       [steps+1];
    double null             [steps+1]; zeros(null, steps+1);
    merge_arrays(11, params, steps, params_extended);
  
    solve_doub_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, &RuKu_4 );
    save_numb_list7( t_values, theta1, theta1_dot, theta2, theta2_dot, params_extended, E_values, "../data/data_doub_pend.txt" );

    return 0;
}


int double_poincare(void)
{    
    double t_end        = 1000.0;            // simulation time [seconds]            
    double h            = 0.01;              // step size [steps per second]
    double g_grav       = 9.81;
    double mass_1       = 1.0;
    double mass_2       = 1.0;
    double length_1     = 1.0;
    double length_2     = 1.0;
    double theta1_0     = 0.0 *M_PI;        // fixed for given Poincare-Section
    double theta1_dot_0 = 0.0 *M_PI;        // fixed for given Poincare-Section
    double theta2_0;                        // iterates through [0, theta2_max]
    double theta2_dot_0;                    // calculated to make the energy stay the same
    
    double E_value      = 1.0;                          // energy for whole simulation
    int    repitions    = 5;                            // number of iterations for diff. (theta2_0, theta2_dot_0) & same energy
    int    max_points   = 10000;                        // max. number of data_points per iteration
    int    steps        = (int)(t_end/h);      
    int    n_points     = (int)fmin(steps,max_points);

    double params[] = {t_end, h, g_grav, mass_1, mass_2, length_1, length_2, theta1_0, theta1_dot_0, theta2_0, theta2_dot_0, E_value, repitions, max_points};
    double theta2_max   = find_doub_theta2_max( 0.0, M_PI, 0.0001, params );
    double **data = create_2d_matrix( 4*repitions, n_points+1, 0.0);

    for(int j = 0; j < repitions; j++)
    {
        erase_last_line();
        printf("--- Double Poincare (%d/%d) \n", j+1, repitions);
            
            // choose location for saving solution
        double *t_values    = data[4*j +0];
        double *theta2      = data[4*j +1];
        double *theta2_dot  = data[4*j +2];
        double *params_ext  = data[4*j +3];

            // set (theta2_0, theta2_dot_0) to new values
        params[9] = theta2_max * j/(double)(repitions); 
        calc_doub_theta2_0(params, E_value);
        
            // hunt for points in the poincare section and safe them
        merge_arrays(14, params, n_points, params_ext);
        solve_doub_poincare( params, t_values, theta2, theta2_dot, &RuKu_6 );
    }
        
    save_matrix( data, 4*repitions, n_points+1, "../data/data_doub_poincare.txt");
    free_2d_matrix(data);

    return 0;
}


// -----------------------------------------------------------------------------


int triple_pendulum(void)
{
    erase_last_line();
    printf("--- Triple Pendulum \n");
    
    double t_end         = 10.0;             // simulation time [seconds]            
    double h             = 0.001;            // step size [steps per second]
    double g_grav        = 9.81;
    double mass_1        = 1.0;
    double mass_2        = 100.0;
    double mass_3        = 1.0;
    double length_1      = 1.0;
    double length_2      = 1.0;
    double length_3      = 1.0;
    double theta_1_0     = 0.0 *M_PI;
    double theta_1_dot_0 = 0.0 *M_PI;
    double theta_2_0     = 0.5 *M_PI;
    double theta_2_dot_0 = 0.0 *M_PI;
    double theta_3_0     = 0.0 *M_PI;
    double theta_3_dot_0 = 0.0 *M_PI;

    int    steps = (int)(t_end/h); 
    double params[] = { t_end, h, g_grav, mass_1, mass_2, mass_3, length_1, length_2, length_3, theta_1_0, theta_1_dot_0, theta_2_0, theta_2_dot_0, theta_3_0, theta_3_dot_0 };
    double params_extended  [steps+1];
    double t_values         [steps+1];
    double E_values         [steps+1];   
    double theta1           [steps+1];
    double theta1_dot       [steps+1];
    double theta2           [steps+1];
    double theta2_dot       [steps+1];
    double theta3           [steps+1];
    double theta3_dot       [steps+1];
    double null             [steps+1]; zeros(null, steps+1);

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    
    merge_arrays(15, params, steps, params_extended);
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_pend.txt" );

    return 0;
}

int triple_chaos(void)
{
    erase_last_line();
    printf("--- Triple Pendulum Chaos \n");
    
    double t_end         = 10.0;             // simulation time [seconds]            
    double h             = 0.001;            // step size [steps per second]
    double g_grav        = 9.81;
    double mass_1        = 1.0;
    double mass_2        = 1.0;
    double mass_3        = 1.0;
    double length_1      = 1.0;
    double length_2      = 1.0;
    double length_3      = 1.0;
    double theta_1_0     = 1.0 *M_PI;
    double theta_1_dot_0 = 0.0 *M_PI;
    double theta_2_0     = 0.0 *M_PI;
    double theta_2_dot_0 = 0.0 *M_PI;
    double theta_3_0     = 0.1 *M_PI;
    double theta_3_dot_0 = 0.3 *M_PI;

    double theta_1_eps   = 0.0001 *M_PI;

    int    steps = (int)(t_end/h); 
    double params[] = { t_end, h, g_grav, mass_1, mass_2, mass_3, length_1, length_2, length_3, theta_1_0, theta_1_dot_0, theta_2_0, theta_2_dot_0, theta_3_0, theta_3_dot_0 };
    double params_extended  [steps+1];
    double t_values         [steps+1];
    double E_values         [steps+1];   
    double theta1           [steps+1];
    double theta1_dot       [steps+1];
    double theta2           [steps+1];
    double theta2_dot       [steps+1];
    double theta3           [steps+1];
    double theta3_dot       [steps+1];
    double null             [steps+1]; zeros(null, steps);
    merge_arrays(15, params, steps, params_extended);

        // testing sensitivity on initial condition
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_a.txt" );

    params[9] = theta_1_0 + theta_1_eps; 
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_b.txt" );

    params[9] = theta_1_0 - theta_1_eps;
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_c.txt" );

        // testing sensitivity on numerical solver
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &Euler );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_Euler.txt" );

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_4 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_RuKu4.txt" );

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_RuKu6.txt" );

    return 0;
}


int triple_compare(void)
{
    erase_last_line();
    printf("--- Triple Pendulum Compared \n");

        // shared parameters of simple and triple pendulum
    double t_end         = 10.0;                       
    double h             = 0.001;            
    double g_grav        = 9.81;
    double theta_1_0     = 0.2 *M_PI;
    double theta_1_dot_0 = 0.1 *M_PI;
    double theta_2_0     = 0.3 *M_PI;
    double theta_2_dot_0 = 0.0 *M_PI;
    double theta_3_0     = 0.0 *M_PI;
    double theta_3_dot_0 = 0.0 *M_PI;
    double mass_1       = 100.0;
    double mass_2       = 1.0;
    double mass_3       = 1.0;
    double length_1     = 1.0;
    double length_2     = 0.1;
    double length_3     = 0.1;
    int    steps        = (int)(t_end/h); 
    
        // initialisation of arrays
    double t_values         [steps+1];
    double E_values         [steps+1];   
    double trip_theta1      [steps+1];
    double trip_theta1_dot  [steps+1];
    double trip_theta2      [steps+1];
    double trip_theta2_dot  [steps+1];
    double trip_theta3      [steps+1];
    double trip_theta3_dot  [steps+1];
    double doub_theta1      [steps+1];
    double doub_theta1_dot  [steps+1];
    double doub_theta2      [steps+1];
    double doub_theta2_dot  [steps+1];
    double theta            [steps+1];              
    double theta_dot        [steps+1]; 
    double trip_params_ex   [steps+1];
    double null             [steps+1];  zeros(null, steps+1);


    double params_trip[] = {t_end, h, g_grav, mass_1, mass_2, mass_3, length_1, length_2, length_3, theta_1_0, theta_1_dot_0, theta_2_0, theta_2_dot_0, theta_3_0, theta_3_dot_0 };
    double params_simp[] = {t_end, h, g_grav, length_1, theta_1_0, theta_1_dot_0};
    merge_arrays( 14, params_trip, steps, trip_params_ex);

        // solving triple and simple pendulum
    solve_trip_pend( params_trip, t_values, E_values, trip_theta1, trip_theta1_dot, trip_theta2, trip_theta2_dot, trip_theta3, trip_theta3_dot, &RuKu_6 );
    solve_simp_pend( params_simp, t_values, theta, theta_dot, &derhs_simp_pend, &RuKu_6 );       
    save_numb_list7(t_values, theta, trip_theta1, trip_params_ex, null, null, null, "../data/data_trip_compare_3-1.txt" );  

    mass_1 = 100;
    mass_2 = 100;
    mass_3 = 0.1;
    length_1 = 1.0;
    length_2 = 1.0;
    length_3 = 0.1;
    
    double params_doub[] = {t_end, h, g_grav, mass_1, mass_2, length_1, length_2, theta_1_0, theta_1_dot_0, theta_2_0, theta_2_dot_0};
    double params_trip_2[] = { t_end, h, g_grav, mass_1, mass_2, mass_3, length_1, length_2, length_3, theta_1_0, theta_1_dot_0, theta_2_0, theta_2_dot_0, theta_3_0, theta_3_dot_0 };
    merge_arrays( 14, params_trip_2, steps, trip_params_ex);

        // solving triple and double pendulum
    solve_trip_pend( params_trip_2, t_values, E_values, trip_theta1, trip_theta1_dot, trip_theta2, trip_theta2_dot, trip_theta3, trip_theta3_dot, &RuKu_6 );
    solve_doub_pend( params_doub, t_values, E_values, doub_theta1, doub_theta1_dot, doub_theta2, doub_theta2_dot, &RuKu_6 );
    save_numb_list7(t_values, doub_theta1, doub_theta2, trip_theta1, trip_theta2, trip_params_ex, null, "../data/data_trip_compare_3-2.txt" );  

}


// -----------------------------------------------------------------------------


int test_num_solvers(void)
{
    erase_last_line();
    printf("--- Test Numerical Solvers \n");
    
    double t_end        = 20.0;         
    double g_grav       = 9.81;
    double length       = 1.0;       
    double theta_0      = 0.1 * M_PI; 
    double theta_dot_0  = 0.0 * M_PI; 
    int    steps;
    double h;
    int    h_steps      = 200;
    double h_low        = 0.005;
    double h_max        = 0.1;
    double h_delta      = h_max / (double)h_steps;
    double h_array      [h_steps+1];
    double dev_array_Eul[h_steps+1];    
    double dev_array_RK4[h_steps+1];
    double dev_array_RK6[h_steps+1];
    double null         [h_steps+1];    zeros(null, h_steps);
    double params_ext   [h_steps+1];
    double params[]     = {t_end, h, g_grav, length, theta_0, theta_dot_0};

    h_array[0]       = h_steps;
    dev_array_Eul[0] = h_steps;
    dev_array_RK4[0] = h_steps;
    dev_array_RK6[0] = h_steps;

    for( int i = 0; i < h_steps; i++ )
    {
            // choose diff stepsizes h (logarithmic equally distributed)
        h = h_max * exp( - log(h_max/h_low) * i/(double)h_steps );              
        params[1]   = h;
        steps       = (int)(t_end/h);      
        double *t_values      = (double *)malloc( (steps+1) * sizeof(double) );
        double *theta_Eul     = (double *)malloc( (steps+1) * sizeof(double) ); 
        double *theta_RK4     = (double *)malloc( (steps+1) * sizeof(double) ); 
        double *theta_RK6     = (double *)malloc( (steps+1) * sizeof(double) );         
        double *theta_dot     = (double *)malloc( (steps+1) * sizeof(double) );    
        double *theta_analyt  = (double *)malloc( (steps+1) * sizeof(double) );            
        
            // solve analyt pendel for diff. h
        solve_analyt_pend( params, theta_analyt );
        modulus_array(theta_analyt, -M_PI, +M_PI);
        h_array[i+1]    = h;

            // for each numerical solver: solve the approx. (analyt) pendulum for different h
        solve_simp_pend( params, t_values, theta_Eul, theta_dot, derhs_analyt_pend, &Euler );       
        modulus_array(theta_Eul, -M_PI, +M_PI);
        dev_array_Eul[i+1]  = num_max_deviation(theta_Eul, theta_analyt);        

        solve_simp_pend( params, t_values, theta_RK4, theta_dot, derhs_analyt_pend, &RuKu_4 );       
        modulus_array(theta_RK4, -M_PI, +M_PI);
        dev_array_RK4[i+1]  = num_max_deviation(theta_RK4, theta_analyt);        

        solve_simp_pend( params, t_values, theta_RK6, theta_dot, derhs_analyt_pend, &RuKu_6 );       
        modulus_array(theta_RK6, -M_PI, +M_PI);
        dev_array_RK6[i+1]  = num_max_deviation(theta_RK6, theta_analyt);        


        free(theta_analyt); free(t_values); free(theta_Eul); free(theta_RK4); free(theta_RK6); free(theta_dot); 
    }

    merge_arrays(6, params, h_steps, params_ext);                           
    save_numb_list7(h_array, params_ext, dev_array_Eul, dev_array_RK4, dev_array_RK6, null, null, "../data/data_test_num_solver.txt" );   // save data in .txt
    
    return 0;
}