
#ifndef __My_super_Header__
#define __My_super_Header__

    #include <math.h>           // for math operations "sin", "exp"
    #define M_PI 3.14159265358979323846

    #include <stdio.h>          // for using "printf"
    #include <stdlib.h>         // for using "maloc"
    #include <string.h>         // for using "memcopy"

        // numerical solver
    void    Euler               ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    void    RuKu_4              ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    void    RuKu_6              ( int n_ODE, double h, double t, double y[], 
                                    void (*derhs) ( int, double, double[], double[], double[] ), double params[] );
    int     test_num_solvers    ( void );

        // simple pendulum
    int     simple_pendulum     (void);
    double  stand_up_theta_dot  ( double length, double g_grav, double theta_0 );
    int     solve_simp_pend     ( double params[], double t_values[], double theta_sol[], double theta_dot_sol[], 
                                    void (*derhs)( int, double, double[], double[], double[] ), 
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    int     solve_analyt_pend   ( double params[], double* y_analytic );
    void    derhs_simp_pend     ( int nDifEqu, double t, double y[], double k[], double params[] );
    void    derhs_analyt_pend   ( int nDifEqu, double t, double y[], double k[], double params[] );
        
        // double pendulum
    int     double_pendulum     (void);
    int     solve_doub_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
    void    derhs_doub_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_doub_energy    ( double y[], double params[], double *E_value );
    int     double_poincare     (void);
    int     calc_doub_theta2_0  ( double params[], double E_value );
    int     solve_doub_poincare ( double params[], double t_values[], double theta2_sol[], double theta2_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) ); 
<<<<<<< HEAD
=======
    double  find_doub_theta2_max( double theta2_min, double theta2_max, const double tol, double params[] );  
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e


        // triple pendulum
    int     triple_compare      (void);
    int     triple_chaos        (void);    
    int     triple_pendulum     (void);
    int     solve_trip_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[], double theta3_sol[], double theta3_dot_sol[],
                                    void (*num_solver)( int, double, double, double[], void (*derhs)(int,double,double[],double[],double[]), double[] ) );
    void    derhs_trip_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_trip_energy    ( double y[], double params[], double *E_value );
    

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
    int     modulus             ( double array[], double limit_up );
<<<<<<< HEAD
=======
    double  modulus_s           ( double value, double limit_up );
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e
    double  *add_IP             ( double *array1, double *array2, double *result, int length);
    double  *scale_IP           ( double *array1, double scalar, double *result, int length);
    double  *linear_comb_arrays ( double** arrays, double* coeffs, double *result, int array_length, int array_numb );
    int     zeros               ( double array[], int length );
    int     copy_array          ( double array[], double copy[], int length);
<<<<<<< HEAD
=======
    int     erase_last_line     (void);
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e

#endif
