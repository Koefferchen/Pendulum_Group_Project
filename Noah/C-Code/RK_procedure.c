
#include "header.h"

      // Calculates 1 timestep of Euler methode (2. Order) to solve system of ODE's:      Y' = f(Y, t)  --> derhs
void Euler ( int n_ODE,          // number of Ordinary Differential Equations
              double h,          // step size 
              double t,          // curve parameter (time) 
              double y[],        // trajectory [n_ODE] with input: y(t) and output; y(t+h) 
              void (*derhs) ( int, double, double[], double[], double[] ),    // implementation of the ODE
              double params[] )  // all constants & initial conditions
{

      double result[n_ODE];    zeros(result, n_ODE);
      double k1[n_ODE];        zeros(k1, n_ODE);
      double k2[n_ODE];        zeros(k2, n_ODE);        
      double *k_array[] = {y, k1, k2};

      double coeffs_1[] = {1.0,  0.0,     0.0};
      double coeffs_2[] = {1.0,  h/2,     0.0};
     
      double coeffs_y[] = {1.0,  0.0,     h};

      double coeffs_h[] = {0.0, 1/2.0};

      (derhs)( n_ODE, t + h *coeffs_h[0], linear_comb_arrays( k_array, coeffs_1, result, n_ODE, 3), k1, params );
      (derhs)( n_ODE, t + h *coeffs_h[1], linear_comb_arrays( k_array, coeffs_2, result, n_ODE, 3), k2, params );

      linear_comb_arrays(k_array, coeffs_y, result, n_ODE, 3);
      copy_array( result, y, n_ODE);
} 



      // Calculates 1 timestep of Runge-Kutta (4. Order) to solve system of ODE's:      Y' = f(Y, t)  --> derhs
void RuKu_4 ( int n_ODE,         // number of Ordinary Differential Equations
              double h,          // step size 
              double t,          // curve parameter (time) 
              double y[],        // trajectory [n_ODE] with input: y(t) and output; y(t+h) 
              void (*derhs) ( int, double, double[], double[], double[] ),    // implementation of the ODE
              double params[] )  // all constants & initial conditions
{

      double result[n_ODE];    zeros(result, n_ODE);
      double k1[n_ODE];        zeros(k1, n_ODE);
      double k2[n_ODE];        zeros(k2, n_ODE);        
      double k3[n_ODE];        zeros(k3, n_ODE);
      double k4[n_ODE];        zeros(k4, n_ODE);
      double *k_array[] = {y, k1, k2, k3, k4};

      double coeffs_1[] = {1.0,  0.0,     0.0,     0.0,     0.0};
      double coeffs_2[] = {1.0,  h/2,     0.0,     0.0,     0.0};
      double coeffs_3[] = {1.0,  0.0,     h/2,     0.0,     0.0};
      double coeffs_4[] = {1.0,  0.0,     0.0,     h,       0.0};
      
      double coeffs_y[] = {1.0,  h/3.0,   h/6.0,   h/6.0,   h/3.0};

      double coeffs_h[] = {0.0,  1/2.0,   1/2.0,   1.0};

      (derhs)( n_ODE, t + h *coeffs_h[0], linear_comb_arrays( k_array, coeffs_1, result, n_ODE, 5), k1, params );
      (derhs)( n_ODE, t + h *coeffs_h[1], linear_comb_arrays( k_array, coeffs_2, result, n_ODE, 5), k2, params );
      (derhs)( n_ODE, t + h *coeffs_h[2], linear_comb_arrays( k_array, coeffs_3, result, n_ODE, 5), k3, params );
      (derhs)( n_ODE, t + h *coeffs_h[3], linear_comb_arrays( k_array, coeffs_4, result, n_ODE, 5), k4, params );
      
      linear_comb_arrays(k_array, coeffs_y, result, n_ODE, 5);
      copy_array( result, y, n_ODE);
} 



      // Calculates 1 timestep of Runge-Kutta (6. Order) to solve system of ODE's:      Y' = f(Y, t)  --> derhs
void RuKu_6 ( int n_ODE,         // number of Ordinary Differential Equations
              double h,          // step size 
              double t,          // curve parameter (time) 
              double y[],        // trajectory [n_ODE] with input: y(t) and output; y(t+h) 
              void (*derhs) ( int, double, double[], double[], double[] ),    // implementation of the ODE
              double params[] )  // all constants & initial conditions
{

      double result[n_ODE];    zeros(result, n_ODE);
      double k1[n_ODE];        zeros(k1, n_ODE);
      double k2[n_ODE];        zeros(k2, n_ODE);        
      double k3[n_ODE];        zeros(k3, n_ODE);
      double k4[n_ODE];        zeros(k4, n_ODE);
      double k5[n_ODE];        zeros(k5, n_ODE);
      double k6[n_ODE];        zeros(k6, n_ODE);
      double k7[n_ODE];        zeros(k7, n_ODE);
      double *k_array[] = {y, k1, k2, k3, k4, k5, k6, k7};

      double coeffs_1[] = {1.0,  0.0,           0.0,        0.0,        0.0,        0.0,     0.0,     0.0};
      double coeffs_2[] = {1.0,  h/3.0,         0.0,        0.0,        0.0,        0.0,     0.0,     0.0};
      double coeffs_3[] = {1.0,  0.0,           h*2/3.0,    0.0,        0.0,        0.0,     0.0,     0.0};
      double coeffs_4[] = {1.0,  h/12.0,        h/3.0,      -h/12.0,    0.0,        0.0,     0.0,     0.0};
      double coeffs_5[] = {1.0,  -h/16.0,       h*9/8.0,    -h*3/16.0,  -h*3/8.0,   0.0,     0.0,     0.0};
      double coeffs_6[] = {1.0,  0.0,           h*9/8.0,    -h*3/8.0,   -h*3/4.0,   h/2.0,   0.0,     0.0};
      double coeffs_7[] = {1.0,  h*9/44.0,      -h*9/11.0,  h*63/44.0,  h*18/11.0,  0.0,     -h*16/11.0,0.0};

      double coeffs_y[] = {1.0,  h*11/120.0,0.0,     h*27/40.0, h*27/40.0, -h*4/15.0, -h*4/15.0, h*11/120.0};
      
      double coeffs_h[] = {0.0,  1/3.0,   2/3.0,   1/3.0,   1/2.0,   1/2.0,   1.0};

      (derhs)( n_ODE, t + h *coeffs_h[0], linear_comb_arrays( k_array, coeffs_1, result, n_ODE, 8), k1, params );
      (derhs)( n_ODE, t + h *coeffs_h[1], linear_comb_arrays( k_array, coeffs_2, result, n_ODE, 8), k2, params );
      (derhs)( n_ODE, t + h *coeffs_h[2], linear_comb_arrays( k_array, coeffs_3, result, n_ODE, 8), k3, params );
      (derhs)( n_ODE, t + h *coeffs_h[3], linear_comb_arrays( k_array, coeffs_4, result, n_ODE, 8), k4, params );
      (derhs)( n_ODE, t + h *coeffs_h[4], linear_comb_arrays( k_array, coeffs_5, result, n_ODE, 8), k5, params );
      (derhs)( n_ODE, t + h *coeffs_h[5], linear_comb_arrays( k_array, coeffs_6, result, n_ODE, 8), k6, params );
      (derhs)( n_ODE, t + h *coeffs_h[6], linear_comb_arrays( k_array, coeffs_7, result, n_ODE, 8), k7, params );
      
      linear_comb_arrays(k_array, coeffs_y, result, n_ODE, 8);
      copy_array( result, y, n_ODE);      
}

