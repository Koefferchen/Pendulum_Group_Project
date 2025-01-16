// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


#include "header.h"

// Berechnet 1 Runge-Kutta Schritt (4. Ordnung) zur Loesung eines DGL-Systems:      Y' = f(Y, t) = derhs

void RuKu_4 ( int nDifEqu,      // # der Differentialgleichungen 
              double h,         // Schrittweite 
              double t,         // Kurvenparameter 
              double y[],       // Bahnkurve [nDifEqu] mit Eingabe: y(t) und Rueckgabe; y(t+h) 
              double yh[],      // Hilfsfeld [nDifEqu] 
              double k1[], double k2[], double k3[], double k4[],       // Hilfssteigungen [nDifEqu] 
              void (*derhs) ( int, double, double[], double[] )         // Funktion zur Berechnung der rechten Seite der DGL f(Y, t) 
){

      double h2 = 0.5 * h;
      int i;                       

      // Tangente an den Anfangspunkt
      (derhs)( nDifEqu, t, y, k1 );        
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k1[i];        // y(t + h/2) 
      }
      // Tangente an das Hilfsfeld
      (derhs)( nDifEqu, t+h2, yh, k2 );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k2[i];        // y(t + h/2) verbessert
      }
      // verbesserte Tangente an den Anfangspunkt
      (derhs)( nDifEqu, t+h2, yh, k3 );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h  * k3[i];        // y(t + h) 
      }
      // verbessserte Tangente an das Hilfsfeld
      (derhs)( nDifEqu, t+h , yh, k4 );
      for (i=0; i<nDifEqu; i++) { 
         y[i] += ( h2 * (k1[i]+k4[i]) + h * (k2[i]+k3[i]) ) / 3;       // y(t + h) verbessert
      } 

}




void RuKu_4_adapted ( int nDifEqu,      // # der Differentialgleichungen 
              double h,         // Schrittweite 
              double t,         // Kurvenparameter 
              double y[],       // Bahnkurve [nDifEqu] mit Eingabe: y(t) und Rueckgabe; y(t+h) 
              double yh[],      // Hilfsfeld [nDifEqu] 
              double k1[], double k2[], double k3[], double k4[],       // Hilfssteigungen [nDifEqu] 
              void (*derhs) ( int, double, double[], double[], double[], double ),        // Funktion zur Berechnung der rechten Seite der DGL f(Y, t) 
              double vaccinated[]        // zur Verwendung vergangenen Impfstandes
){

      double h2 = 0.5 * h;
      int i;                       

      // Tangente an den Anfangspunkt
      (derhs)( nDifEqu, t, y, k1, vaccinated, h );        
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k1[i];        // y(t + h/2) 
      }
      // Tangente an das Hilfsfeld
      (derhs)( nDifEqu, t+h2, yh, k2, vaccinated, h );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h2 * k2[i];        // y(t + h/2) verbessert
      }
      // verbesserte Tangente an den Anfangspunkt
      (derhs)( nDifEqu, t+h2, yh, k3, vaccinated, h );
      for (i=0; i<nDifEqu; i++) { 
         yh[i] = y[i] + h  * k3[i];        // y(t + h) 
      }
      // verbessserte Tangente an das Hilfsfeld
      (derhs)( nDifEqu, t+h , yh, k4, vaccinated, h );
      for (i=0; i<nDifEqu; i++) { 
         y[i] += ( h2 * (k1[i]+k4[i]) + h * (k2[i]+k3[i]) ) / 3;       // y(t + h) verbessert
      } 

}
