// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


#ifndef __My_super_Header__
#define __My_super_Header__

    #include <math.h>           // notwendig für trig, exp & mehr
    #define M_PI 3.14159265358979323846

    #include <stdio.h>          // notwendig für printf
    #include <stdlib.h>         // notwendig für malloc
    #include <string.h>         // notwendig für memcopy

    void RuKu_4 (   int nDifEqu,      // # der Differentialgleichungen 
                    double h,         // Schrittweite 
                    double t,         // Kurvenparameter 
                    double y[],       // Bahnkurve [nDifEqu] mit Eingabe: y(t) und Rueckgabe; y(t+h) 
                    double yh[],      // Hilfsfeld [nDifEqu] 
                    double k1[], double k2[], double k3[], double k4[],       // Hilfsfelder [nDifEqu] 
                    void (*derhs) ( int, double, double[], double[] )         // Funktion zur Berechnung der rechten Seite: 
        );

    void RuKu_4_adapted (   int nDifEqu,              // # der Differentialgleichungen 
                            double h,                 // Schrittweite 
                            double t,                 // Kurvenparameter 
                            double y[],               // Bahnkurve [nDifEqu] mit Eingabe: y(t) und Rueckgabe; y(t+h) 
                            double yh[],              // Hilfsfeld [nDifEqu] 
                            double k1[], double k2[], double k3[], double k4[],                         // Hilfssteigungen [nDifEqu] 
                            void (*derhs) ( int, double, double[], double[], double[], double ),        // Funktion zur Berechnung der rechten Seite der DGL f(Y, t) 
                            double vaccinated[]       // zur Verwendung vergangenen Impfstandes
        );

        // Definitionen der DGLs 1. Ordnung 
    void derhs_bsp      ( int nDifEqu, double t, double y[], double k[] );
    void derhs_SIR_1    ( int nDifEqu, double t, double y[], double k[] );
    void derhs_SIR_2    ( int nDifEqu, double t, double y[], double k[] );
    void derhs_SIRV     ( int nDifEqu, double t, double y[], double k[] );
    void derhs_SIRVD    ( int nDifEqu, double t, double y[], double k[] );
    void derhs_SIRVD_2  ( int nDifEqu, double t, double y[], double k[], double vaccinated[], double h);
    void derhs_SIRVDZ   ( int nDifEqu, double t, double y[], double k[], double vaccinated[], double h);

        // Verwendung der DGLs und des RK_4 zur numerischen Lösung
    int solve_bsp       ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[] );
    int solve_SIR_1     ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[] ); 
    int solve_SIR_2     ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[] ); 
    int solve_SIRV      ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[] ); 
    int solve_SIRVD     ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[] ); 
    int solve_SIRVD_2   ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[] ); 
    int solve_SIRVDZ    ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[], double y6_sol[] );  



    double beta_2   ( double t );   // infection rate (t)
    double v        ( double t );   // vaccination rate (t)

        // Hilfsfunktionen zum Speichern und Interpolieren der berechneten Lösungen
    double  *create_null        ( int length);
    double  interpolation       ( double array[], double index );
    int     display_numb_list   ( double* numb_list);
    int     save_numb_list7     ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, char* save_as );

#endif
