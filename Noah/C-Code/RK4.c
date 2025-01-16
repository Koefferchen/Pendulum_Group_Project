
#include "header.h"

// Calculates 1 timestep of Runge-Kutta (4. Order) to solve system of DE's:      Y' = f(Y, t)  --> derhs

void RuKu_4 ( int nDifEqu,      // number of Differential Equations (DE)
              double h,         // step size 
              double t,         // curve parameter (time) 
              double y[],       // trajectory [nDifEqu] with input: y(t) and output; y(t+h) 
              double yh[],      // help array [nDifEqu] 
              double k1[], double k2[], double k3[], double k4[],       // help parameters [nDifEqu] 
              void (*derhs) ( int, double, double[], double[] )         // function that calculates the right-hand-side of the DE 
){

      double h2 = 0.5 * h;
      int i;                       

      // Tangent at initial point 
      (derhs)( nDifEqu, t, y, k1 );        
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k1[i];        // y(t + h/2) 
      }
      // Tangent at the help point      
      (derhs)( nDifEqu, t+h2, yh, k2 );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k2[i];        // y(t + h/2) improved
      }
      // improved Tangent at initial point
      (derhs)( nDifEqu, t+h2, yh, k3 );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h  * k3[i];        // y(t + h) 
      }
      // improved Tangent at the help point
      (derhs)( nDifEqu, t+h , yh, k4 );
      for (i=0; i<nDifEqu; i++) { 
         y[i] += ( h2 * (k1[i]+k4[i]) + h * (k2[i]+k3[i]) ) / 3;       // y(t + h) improved
      } 

}

