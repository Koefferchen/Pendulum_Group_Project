
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

    printf("Simulation Data generated\n");
    return 0;
}

int simple_pendulum(void)
{

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
    double* null = create_null(steps);  
    merge_arrays(6, params, steps, params_extended);

    solve_simp_pend( params, t_values, theta, theta_dot, &derhs_simp_pend, &RuKu_6 );       
    solve_analyt_pend( params, theta_analy );                                  
    save_numb_list7(t_values, theta, theta_dot, theta_analy, params_extended, null, null, "../data/data_simp_pend.txt" );   // save data in .txt

    free(null);
    return 0;
}

int double_pendulum(void)
{
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
    double* null = create_null(steps);
    merge_arrays(11, params, steps, params_extended);
  
    solve_doub_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, &RuKu_4 );
    save_numb_list7( t_values, theta1, theta1_dot, theta2, theta2_dot, params_extended, E_values, "../data/data_doub_pend.txt" );

    free(null);
    return 0;
}

int triple_pendulum(void)
{
    double t_end         = 20.0;             // simulation time [seconds]            
    double h             = 0.001;            // step size [steps per second]
    double g_grav        = 9.81;
    double mass_1        = 100.0;
    double mass_2        = 10.0;
    double mass_3        = 1.0;
    double length_1      = 1.0;
    double length_2      = 1.0;
    double length_3      = 1.0;
    double theta_1_0     = 0.5 *M_PI;
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
    double* null = create_null(steps);
    merge_arrays(15, params, steps, params_extended);

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, &RuKu_6 );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_pend.txt" );

    free(null);
    return 0;
}


int triple_chaos(void)
{
    double t_end         = 20.0;             // simulation time [seconds]            
    double h             = 0.001;            // step size [steps per second]
    double g_grav        = 9.81;
    double mass_1        = 1.0;
    double mass_2        = 1.0;
    double mass_3        = 1.0;
    double length_1      = 1.0;
    double length_2      = 1.0;
    double length_3      = 1.0;
    double theta_1_0     = 0.5 *M_PI;
    double theta_1_dot_0 = 0.0 *M_PI;
    double theta_2_0     = 0.3 *M_PI;
    double theta_2_dot_0 = 0.2 *M_PI;
    double theta_3_0     = 0.0 *M_PI;
    double theta_3_dot_0 = 0.0 *M_PI;

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
    double* null = create_null(steps);
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

    free(null);
    return 0;
}

int double_poincare(void)
{    
    
    double t_end        = 50000.0;            // simulation time [seconds]            
    double h            = 0.1;              // step size [steps per second]
    double g_grav2      = 9.81;
    double mass_1       = 1.0;
    double mass_2       = 1.0;
    double length_1     = 1.0;
    double length_2     = 1.0;
    
    double theta1_0     = 0.0 *M_PI;       // fixed for given Poincare-Section
    double theta1_dot_0 = 0.0 *M_PI;              // fixed for given Poincare-Section
    double theta2_0;                        // iterates through [0, Pi]
    double theta2_dot_0;                    // calculated to make the energy stay the same
    double E_value      = 0.01;
    int    repitions    = 1;
    
    int    steps = (int)(t_end/h);      // Initialisation
    double params[] = {t_end, h, g_grav2, mass_1, mass_2, length_1, length_2, theta1_0, theta2_0, theta1_dot_0, theta2_dot_0, E_value, repitions};
    double **data = create_2d_matrix( 4*repitions, steps+1, 0.0);

    double condition_eq_1[] = {0.0, 0.0, M_PI, 0.0};
    double condition_eq_2[] = {M_PI, 0.0, 0.0, 0.0};
    double condition_eq_3[] = {M_PI, 0.0, M_PI, 0.0};
    //calc_doub_energy(condition_eq_1, params, &params[11]);


    for(int j = 0; j < repitions; j++)
    {
        double *t_values    = data[4*j +0];
        double *theta2      = data[4*j +1];
        double *theta2_dot  = data[4*j +2];
        double *params_ext  = data[4*j +3];
        params[8] = M_PI/20.0 + j * 2*M_PI / (double)(repitions); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        calc_doub_theta2_0(params, E_value);
        merge_arrays(13, params, steps, params_ext);
        solve_doub_poincare( params, t_values, theta2, theta2_dot, &RuKu_6 );
    }

    save_matrix( data, 4*repitions, steps+1, "../data/data_doub_poincare.txt");
    return 0;
}


int test_num_solvers(void)
{
    double t_end        = 100.0;         
    double g_grav       = 9.81;
    double length       = 1.0;       
    double theta_0      = 0.5 * M_PI; 
    double theta_dot_0  = 0.2 * M_PI; 
    int    steps;
    double h;
    int    h_steps      = 200;
    double h_end        = 1.0;
    double h_delta      = h_end / (double)h_steps;
    double h_array      [h_steps+1];
    double dev_array_Eul[h_steps+1];    
    double dev_array_RK4[h_steps+1];
    double dev_array_RK6[h_steps+1];
    double null         [h_steps+1];    zeros(null, h_steps);
    double params_ext   [h_steps+1];
    double params[] = {t_end, h, g_grav, length, theta_0, theta_dot_0};

    h_array[0]       = h_steps;
    dev_array_Eul[0] = h_steps;
    dev_array_RK4[0] = h_steps;
    dev_array_RK6[0] = h_steps;

    for( int i = 0; i < h_steps; i++ )
    {
        h = h_end - i * h_delta;
        params[1]   = h;
        steps       = (int)(t_end/h);      
        double *t_values      = (double *)malloc( (steps+1) * sizeof(double) );
        double *theta         = (double *)malloc( (steps+1) * sizeof(double) );         
        double *theta_dot     = (double *)malloc( (steps+1) * sizeof(double) );    
        double *theta_analyt  = (double *)malloc( (steps+1) * sizeof(double) );            
        
        solve_analyt_pend( params, theta_analyt );
        modulus(theta_analyt, +2*M_PI);
        h_array[i+1]    = h;

        solve_simp_pend( params, t_values, theta, theta_dot, derhs_analyt_pend, &Euler );       
        modulus(theta, +2*M_PI);
        dev_array_Eul[i+1]  = fmin( fabs(theta[steps] - theta_analyt[steps]), fabs( 2*M_PI - fabs(theta[steps] - theta_analyt[steps]) ) );          // average_diff( theta, theta_analyt);

        solve_simp_pend( params, t_values, theta, theta_dot, derhs_analyt_pend, &RuKu_4 );       
        modulus(theta, +2*M_PI);
        dev_array_RK4[i+1]  = fmin( fabs(theta[steps] - theta_analyt[steps]), fabs( 2*M_PI - fabs(theta[steps] - theta_analyt[steps]) ) );          // average_diff( theta, theta_analyt);

        solve_simp_pend( params, t_values, theta, theta_dot, derhs_analyt_pend, &RuKu_6 );       
        modulus(theta, +2*M_PI);
        dev_array_RK6[i+1]  = fmin( fabs(theta[steps] - theta_analyt[steps]), fabs( 2*M_PI - fabs(theta[steps] - theta_analyt[steps]) ) );          // average_diff( theta, theta_analyt);

        free(theta_analyt); free(t_values); free(theta); free(theta_dot); 
    }
    

    merge_arrays(6, params, h_steps, params_ext);                           
    save_numb_list7(h_array, params_ext, dev_array_Eul, dev_array_RK4, dev_array_RK6, null, null, "../data/data_test_num_solver.txt" );   // save data in .txt
    
    return 0;
}


int triple_compare(void)
{
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
    double null[steps+1];   zeros(null, steps+1);
    
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