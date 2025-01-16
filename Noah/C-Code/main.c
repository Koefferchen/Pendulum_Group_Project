
#include "header.h"

int main(void)
{

    double t_end    = 50.0;                 // simulation time [seconds]            
    double h        = 0.001;                 // step size [steps per second]
    
    int    steps = (int)(t_end/h);          // Initialisation
    double t_values[steps+1];
    double theta[steps+1];                     // theta(t)
    double theta_dot[steps+1];                     // theta_dot(t)
    double theta_analy[steps+1];
    double* null = create_null(steps);      // null

    solve_DE_pend( h, t_end, t_values, theta, theta_dot );       
    solve_analyt_pend( h, t_end, theta_analy);                                  // solve DE using RK4
    save_numb_list7(t_values, theta, theta_dot, theta_analy, null, null, null, "../data/data_pend.txt" );   // save data in .txt

    return 0;
}