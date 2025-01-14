// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


# include "header.h"


// ------------------------- Aufgabe 1.1 ---------------------------

int solve_SIR_1( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[] ) 
{
    int     nDifEqu = 3;    // hier: 3 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;
    y[2] = 0.0;

    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIR_1);
        
        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];

        t = t + h; 
    } 

    return 0;
}


// ------------------------- Aufgabe 1.2 ---------------------------

int solve_SIR_2( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[] ) 
{
    int     nDifEqu = 3;    // hier: 3 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;
    y[2] = 0.0;                 // keine Genesenen zu Beginn

    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIR_2);
        
        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];

        t = t + h; 
    } 

    return 0;
}


// ------------------------- Aufgabe 1.3 ---------------------------

int solve_SIRV( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[] ) 
{
    int     nDifEqu = 4;    // hier: 4 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;               // 0.1% Infizierte zu Beginn
    y[2] = 0.0;                 // keine Genenesen zu Beginn
    y[3] = 0.0;                 // keine Geimpfte zu Beginn

    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;
    y4_sol[0] = steps;

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRV);
        
        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];
        y4_sol[j+1] = y[3];
        t = t + h; 
    } 

    return 0;
}


// ------------------------- Aufgabe 1.4 ---------------------------

int solve_SIRVD( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[] ) 
{
    int     nDifEqu = 5;    // hier: 5 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;               // 0.1% Infizierte zu Beginn
    y[2] = 0.0;                 // keine Genenesen zu Beginn
    y[3] = 0.0;                 // keine Geimpfte zu Beginn
    y[4] = 0.0;                 // keine Tote zu Beginn


    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;
    y4_sol[0] = steps;
    y5_sol[0] = steps;

    for(int j = 0; j < steps; j++)
    {
        RuKu_4( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRVD);
        
        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];
        y4_sol[j+1] = y[3];
        y5_sol[j+1] = y[4];
        t = t + h; 
    } 

    return 0;
}


// ------------------------- Aufgabe 2 -----------------------------

int solve_SIRVD_2( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[] ) 
{
    int     nDifEqu = 5;    // hier: 5 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;               // 0.1% Infizierte zu Beginn
    y[2] = 0.0;                 // keine Genenesen zu Beginn
    y[3] = 0.0;                 // keine Geimpfte zu Beginn
    y[4] = 0.0;                 // keine Tote zu Beginn


    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;
    y4_sol[0] = steps;
    y5_sol[0] = steps;
    RuKu_4_adapted( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRVD_2, y4_sol);

    for(int j = 0; j < steps; j++)
    {
        t = t + h; 
        RuKu_4_adapted( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRVD_2, y4_sol);
        
        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];
        y4_sol[j+1] = y[3];
        y5_sol[j+1] = y[4];

    } 

    return 0;
}


// ------------------------- Aufgabe 3 -----------------------------

int solve_SIRVDZ( double h, double t_end, double t_values[], double y1_sol[], double y2_sol[], double y3_sol[], double y4_sol[], double y5_sol[], double y6_sol[] ) 
{
    int     nDifEqu = 6;    // hier: 6 DGLs
    int     steps = (int)(t_end/h);
    double  t = 0;
    double  y[nDifEqu];     // Initialisierung
    double  yh[nDifEqu];
    double  k1[nDifEqu];
    double  k2[nDifEqu];
    double  k3[nDifEqu];
    double  k4[nDifEqu];

    y[0] = 0.999;               // Anfangsbedingungen in t_0
    y[1] = 0.001;               // 0.1% Infizierte zu Beginn
    y[2] = 0.0;                 // keine Genenesen zu Beginn
    y[3] = 0.0;                 // keine Geimpfte zu Beginn
    y[4] = 0.0;                 // keine Tote zu Beginn
    y[5] = 0.0;                 // keine Zombies zu Beginn


    t_values[0] = steps;
    y1_sol[0] = steps;      // erster Eintrag = Länge des arrays
    y2_sol[0] = steps;
    y3_sol[0] = steps;
    y4_sol[0] = steps;
    y5_sol[0] = steps;
    y6_sol[0] = steps;
    
    RuKu_4_adapted( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRVDZ, y4_sol);

    for(int j = 0; j < steps; j++)
    {
        t = t + h; 
        RuKu_4_adapted( nDifEqu, h, t, y, yh, k1, k2, k3, k4, &derhs_SIRVDZ, y4_sol);

        t_values[j+1] = t;      // speichere Lösungen der DGL in arrays
        y1_sol[j+1] = y[0];
        y2_sol[j+1] = y[1];
        y3_sol[j+1] = y[2];
        y4_sol[j+1] = y[3];
        y5_sol[j+1] = y[4];
        y6_sol[j+1] = y[5];

    } 

    return 0;
}