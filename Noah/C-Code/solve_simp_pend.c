
# include "header.h"

double omega = 0.4;        
double theta_0 = 0.3 * M_PI; 
double theta_dot_0 = -0.2;

void derhs_pend( int nDifEqu, double t, double y[], double k[] )
{
    // here:    nDifEqu = 2    
    k[0] = y[1] ;                           // 1. DE: x(t)_dot = v(t)
    k[1] = - pow(omega, 2)  * sin(y[0]) ;   // 2. DE: v(t)_dot = -omega**2 * x(t)    

}

double stand_up_theta_dot( double omega, double theta_0 )
{
    return omega * pow( 2 * (1+cos(theta_0)), 0.5 );
}

int solve_DE_pend( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[] ) 
{
    int     nDifEqu = 2;            // here: 2 DEs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];             // Initialisation 
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];
    t_values[0] = steps;            // first entry = length of array
    y1_sol[0] = steps;      
    y2_sol[0] = steps;

    y[0] = theta_0;                     // Initial conditions at t_0
    y[1] = theta_dot_0; //stand_up_theta_dot(omega, theta_0);
    
    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_pend);
        
        t_values[j+1] = t;          // save solution for each time step in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];

        t = t + h; 
    } 
    
    return 0;
}

int solve_analyt_pend( double h, double t_end, double* y_analytic )
{
    int steps = (int)(t_end/h);
    y_analytic[0] = steps;

    for(int j=0; j < steps; j++ )
    {
        y_analytic[j+1] = theta_0 * cos( omega * (j*h) ) + theta_dot_0/omega * sin(omega *(j*h) );
    }

    return 0;
}