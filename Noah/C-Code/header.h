
#ifndef __My_super_Header__
#define __My_super_Header__

    #include <math.h>           // for math operations "sin", "exp"
    #define M_PI 3.14159265358979323846

    #include <stdio.h>          // for using "printf"
    #include <stdlib.h>         // for using "maloc"
    #include <string.h>         // for using "memcopy"

    void RuKu_4 ( int nDifEqu,  // number of Differential Equations (DE)
              double h,         // step size 
              double t,         // curve parameter (time) 
              double y[],       // trajectory [nDifEqu] with input: y(t) and output; y(t+h) 
              double yh[],      // help array [nDifEqu] 
              double k1[], double k2[], double k3[], double k4[],           // help parameters [nDifEqu] 
              void (*derhs) ( int, double, double[], double[], double[] ),  // function that calculates the right-hand-side of the DE 
              double params[]               
    );  

        // simple pendulum
    int     simple_pendulum     (void);
    double  stand_up_theta_dot  ( double omega, double theta_0 );
    int     solve_simp_pend     ( double params[], double t_values[], double y1_sol[], double y2_sol[] );
    int     solve_analyt_pend   ( double params[], double* y_analytic );
    void    derhs_pend          ( int nDifEqu, double t, double y[], double k[], double params[] );
        
        // double pendulum
    int     double_pendulum     (void);
    int     solve_doub_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[] ); 
    void    derhs_doub_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_doub_energy    ( double y[], double params[], double *E_value );
    
        // triple pendulum
    int     triple_chaos        (void);    
    int     triple_pendulum     (void);
    int     solve_trip_pend     ( double params[], double t_values[], double E_values[], double theta1_sol[], double theta1_dot_sol[], double theta2_sol[], double theta2_dot_sol[], double theta3_sol[], double theta3_dot_sol[] ); 
    void    derhs_trip_pend     ( int nDifEqu, double t, double y[], double y_dot[], double params[] );
    int     calc_trip_energy    ( double y[], double params[], double *E_value );


        // helper functions
    double  *create_null        ( int length);
    int     display_numb_list   ( double* numb_list);
    int     save_numb_list7     ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, char* save_as );
    int     save_numb_list9     ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, double* numb_list8, double* numb_list9, char* save_as );
    int     merge_arrays        ( int old_length, double old_array[], int new_length, double new_array[] );
#endif
