
#ifndef __My_super_Header__
#define __My_super_Header__

    #include <math.h>           // for math operations "sin", "exp"
    #define M_PI 3.14159265358979323846

    #include <stdio.h>          // for using "printf"
    #include <stdlib.h>         // for using "maloc"
    #include <string.h>         // for using "memcopy"

        // numerical solver
    void    RuKu_2              ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    void    RuKu_4              ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    void    RuKu_6              ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    int     test_num_solvers    ( void );
    int     test_numeric_solver ( double params[], double h_array[], double deviation[],
                                    void (*analyt_sol)( double[], double[]),
                                    void (*derhs)( int, double, double[], double[], double[] ), 
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    double  num_max_deviation   ( double theta2_num[], double theta2_ana[] );


        // simple pendulum
    int     simple_pendulum     (void);
    double  stand_up_theta_dot  ( double length, double g_grav, double theta_0 );
    void    solve_analyt_pend   ( double params[], double* y_analytic );
    void    derhs_simp_pend     ( int nDifEqu, double t, double y[], double k[], double params[] );
    void    derhs_analyt_pend   ( int nDifEqu, double t, double y[], double k[], double params[] );
    int     solve_simp_pend     ( double params[], double t_values[], double theta_sol[], double theta_dot_sol[], 
                                    void (*derhs)( int, double, double[], double[], double[] ), 
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    
        // double pendulum
    int     double_pendulum     (void);
    int     double_flip_over    (void);
    int     double_poincare     (void);
    void    derhs_doub_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_doub_energy    ( double y[], double params[], double *E_value );
    int     calc_doub_thetadot2_0( double params[], double parity );
    double  find_doub_theta2_max( double theta2_min, double theta2_max, const double tol, double params[] );  
    int     solve_doub_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    int     solve_doub_poincare ( double params[], double t_values[], double theta2_sol[], double theta2_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    int     solve_doub_flip     ( double params[], double theta1_0_list[], double theta2_0_list[], double **theta1_sol_matrix, double **theta2_sol_matrix,
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) );


        // triple pendulum
    int     triple_compare      (void);
    int     triple_chaos        (void);    
    int     triple_pendulum     (void);
    void    derhs_trip_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_trip_energy    ( double y[], double params[], double *E_value );
    int     solve_trip_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[], double theta3_sol[], double theta3_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) );


        // helper functions
    double  *create_null        ( int length);
    int     display_numb_list   ( double* numb_list);
    int     save_numb_list7     ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, char* save_as );
    int     save_numb_list9     ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, double* numb_list8, double* numb_list9, char* save_as );
    int     merge_arrays        ( int old_length, double old_array[], int new_length, double new_array[] );
    double  **create_2d_matrix  (int x_dim, int y_dim, double initial_value);
    int     save_matrix         ( double **matrix, int x_dim, int y_dim, char* save_as );
    int     free_2d_matrix      ( double **matrix );
    double  average_diff        ( double array1[], double array2[] );
    int     modulus_array       ( double array[], double limit_low, double limit_up );
    int     modulus             ( double array[], double limit_up );
    double  modulus_s           ( double value, double limit_up );
    double  *add_IP             ( double *array1, double *array2, double *result, int length);
    double  *scale_IP           ( double *array1, double scalar, double *result, int length);
    double  *linear_comb_arrays ( double** arrays, double* coeffs, double *result, int array_length, int array_numb );
    int     zeros               ( double array[], int length );
    int     copy_array          ( double array[], double copy[], int length);
    int     erase_last_line     (void);
    int     print_array         ( double array[] );
    int     fill_linspace       ( double array[], double min_range, double max_range, int length );


#endif
