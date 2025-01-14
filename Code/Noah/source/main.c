// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


#include "header.h"

int main(void)
{

    double t_end_365 = 365.0;                   // betrachte 365 Tage             
    double h_365     = 0.05  ;                    // berechne 20 Werte pro Tag
    
    int    steps_365 = (int)(t_end_365/h_365);  // Initialisierung der Parameter
    double t_values_365[steps_365+1];
    double y1_sol_365[steps_365+1];             // Susceptible (S)
    double y2_sol_365[steps_365+1];             // Infected (I)
    double y3_sol_365[steps_365+1];             // Recovered (R)
    double y4_sol_365[steps_365+1];             // Vaccinated (V)
    double y5_sol_365[steps_365+1];             // Dead (D)
    double y6_sol_365[steps_365+1];             // Zombified (Z)
    double* null = create_null(steps_365);      // (Platzhalter)

        // Aufgabe 1.1
    solve_SIR_1( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365 );                  // Lösung der DGL mithilfe RK_4
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, null, null, null, "./data/SIR_1_data365.txt" );     // Speichern der Datenpunkte in txt
     
        // Aufgabe 1.2
    solve_SIR_2( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365 ); 
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, null, null, null, "./data/SIR_2_data365.txt" );

        // Aufgabe 1.3
    solve_SIRV( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365 ); 
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, null, null, "./data/SIRV_data365.txt" );

        // Aufgabe 1.4
    solve_SIRVD( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365 ); 
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365, null, "./data/SIRVD_data365.txt" );

        // Aufgabe 2
    solve_SIRVD_2( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365); 
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365, null, "./data/SIRVD_2_data365.txt" );

        // Aufgabe 3
    solve_SIRVDZ( h_365, t_end_365, t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365, y6_sol_365 ); 
    save_numb_list7(t_values_365, y1_sol_365, y2_sol_365, y3_sol_365, y4_sol_365, y5_sol_365, y6_sol_365, "./data/SIRVDZ_data365.txt" );



    return 0;
}