
#include "header.h"

int main(void)
{
    simple_pendulum();
    double_pendulum();
    triple_pendulum();

    triple_chaos();

    double_poincare();
    
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

    solve_simp_pend( params, t_values, theta, theta_dot );       
    solve_analyt_pend( params, theta_analy );                                  
    save_numb_list7(t_values, theta, theta_dot, theta_analy, params_extended, null, null, "../data/data_simp_pend.txt" );   // save data in .txt

    free(null);
    return 0;
}

int double_pendulum(void)
{
    double t_end        = 5.0;             // simulation time [seconds]            
    double h            = 0.001;            // step size [steps per second]
    double g_grav2      = 9.81;
    double mass_1       = 1.0;
    double mass_2       = 1.0;
    double length_1     = 1.0;
    double length_2     = 1.0;
    double theta1_0     = 3.0;
    double theta1_dot_0 = 0.0;
    double theta2_0     = 3.0;
    double theta2_dot_0 = 0.0;
    
    int    steps = (int)(t_end/h);      // Initialisation
    double params[] = {t_end, h, g_grav2, mass_1, mass_2, length_1, length_2, theta1_0, theta1_dot_0, theta2_0, theta2_dot_0};
    double params_extended  [steps+1];
    double t_values         [steps+1];
    double E_values         [steps+1];
    double theta1           [steps+1];
    double theta1_dot       [steps+1];
    double theta2           [steps+1];
    double theta2_dot       [steps+1];
    double* null = create_null(steps);
    merge_arrays(11, params, steps, params_extended);
  
    solve_doub_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot );
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

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot );
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
    double mass_2        = 20.0;
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

    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_a.txt" );

    params[9] = theta_1_0 + theta_1_eps;
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_b.txt" );

    params[9] = theta_1_0 - theta_1_eps;
    solve_trip_pend( params, t_values, E_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot );
    save_numb_list9( t_values, theta1, theta1_dot, theta2, theta2_dot, theta3, theta3_dot, params_extended, E_values, "../data/data_trip_chaos_c.txt" );


    free(null);
    return 0;
}

int double_poincare(void)
{    
    
    double t_end        = 400.0;            // simulation time [seconds]            
    double h            = 0.01;             // step size [steps per second]
    double g_grav2      = 9.81;
    double mass_1       = 1.0;
    double mass_2       = 1.0;
    double length_1     = 1.0;
    double length_2     = 1.0;

    double theta1_0     = -0.5 *M_PI;       // fixed for given Poincare-Section
    double theta1_dot_0 = 0.0;              // fixed for given Poincare-Section
    double theta2_0;                        // iterates through [0, Pi]
    double theta2_dot_0;                    // calculated to make the energy stay the same
    double E_value      = 10.0;
    int    repitions    = 40;
    
    int    steps = (int)(t_end/h);      // Initialisation
    double params[] = {t_end, h, g_grav2, mass_1, mass_2, length_1, length_2, theta1_0, theta2_0, theta1_dot_0, theta2_dot_0, E_value, repitions};
    double **data = create_2d_matrix( 4*repitions, steps+1, 0.0);

    for(int j = 0; j < repitions; j++)
    {
        double *t_values    = data[4*j +0];
        double *theta2      = data[4*j +1];
        double *theta2_dot  = data[4*j +2];
        double *params_ext  = data[4*j +3];
        params[8] = j * M_PI / (double)(repitions);
        params[10] = calc_doub_theta2_0(params, E_value);
        merge_arrays(13, params, steps, params_ext);
        solve_doub_poincare( params, t_values, theta2, theta2_dot );

    }

    save_matrix( data, 4*repitions, steps+1, "../data/data_doub_poincare.txt");
    return 0;
}