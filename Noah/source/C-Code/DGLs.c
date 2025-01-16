// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


#include "header.h"


// ------------------------- Aufgabe 1.1 ---------------------------

double beta = 0.4;          // Infection rate (constant)
double gammas = 0.05;       // healing rate (constant)

void derhs_SIR_1( int nDifEqu, double t, double y[], double k[] )
{
    // hier:    nDifEqu = 3    
    k[0] = - beta * y[1] * y[0] ;                    // 1. DGL: Susceptible (S)
    k[1] = + beta * y[1] * y[0] - gammas * y[1] ;    // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] ;                         // 3. DGL: Recovered (R)
}


// ------------------------- Aufgabe 1.2 ---------------------------

double beta_2( double t ){ return 0.045 + 0.355 * pow( cos(3.0 * M_PI * t /200.0), 2.0 ) ; }     // Infectionrate

void derhs_SIR_2( int nDifEqu, double t, double y[], double k[] )
{
    // hier:    nDifEqu = 3    
    k[0] = - beta_2(t) * y[1] * y[0] ;                    // 1. DGL: Susceptible (S)
    k[1] = + beta_2(t) * y[1] * y[0] - gammas * y[1] ;    // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] ;                              // 3. DGL: Recovered (R)
}


// ------------------------- Aufgabe 1.3 ---------------------------

double v( double t ){ return 0.05 * pow( sin(4.0 * M_PI * t /200.0), 2.0 ) ; }   // vaccination rate

void derhs_SIRV( int nDifEqu, double t, double y[], double k[] )
{
    // hier:    nDifEqu = 4    
    k[0] = - beta_2(t) * y[1] * y[0] - v(t) * y[0];       // 1. DGL: Susceptible (S)
    k[1] = + beta_2(t) * y[1] * y[0] - gammas * y[1] ;    // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] ;                              // 3. DGL: Recovered (R)
    k[3] = v(t) * y[0];                                   // 4. DGL: Vaccinated (V)
}


// ------------------------- Aufgabe 1.4 ---------------------------

double psi = 0.005;         // deathrate /day /Infected

void derhs_SIRVD( int nDifEqu, double t, double y[], double k[] )
{
    // hier:    nDifEqu = 5    
    k[0] = - beta_2(t) * y[1] * y[0] - v(t) * y[0];                     // 1. DGL: Susceptible (S)
    k[1] = + beta_2(t) * y[1] * y[0] - gammas * y[1] - psi * y[1];      // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] ;                                            // 3. DGL: Recovered (R)
    k[3] = v(t) * y[0];                                                 // 4. DGL: Vaccinated (V)
    k[4] = psi * y[1];                                                  // 5. DGL: Dead (D)
}


// ------------------------- Aufgabe 2 ----------------------------

double epsilon = 0.01;      // Vaccinated getting Susceptible /day /Vaccinated
double delta_t = 15.0;      // charactistic duration of vaccine-protection

void derhs_SIRVD_2( int nDifEqu, double t, double y[], double k[], double vaccinated[], double h)
{
    
    double prior_vaccinated;                                // berechne % der Geimpften zum Zeitpunkt t - delta_t
    if ( t - delta_t <= 0 ) 
    { 
        prior_vaccinated = 0.0;
    }                  
    else 
    { 
        prior_vaccinated = interpolation(vaccinated, (t-delta_t)/h);
    }

    // hier:    nDifEqu = 5    
    k[0] = - beta_2(t) * y[1] * y[0] - v(t) * y[0] + epsilon * prior_vaccinated;            // 1. DGL: Susceptible (S)
    k[1] = + beta_2(t) * y[1] * y[0] - gammas * y[1] - psi * y[1];                          // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] ;                                                                // 3. DGL: Recovered (R)
    k[3] = v(t) * y[0] - epsilon * prior_vaccinated;                                        // 4. DGL: Vaccinated (V)
    k[4] = psi * y[1];                                                                      // 5. DGL: Dead (D)
}


// ------------------------- Aufgabe 3 ----------------------------

double z = 10e-12;      // Zombification rate by Vaccine  /day /Vaccinated
double zz = 0.1;        // Zombification rate by Zombies  /day /Zombified


void derhs_SIRVDZ( int nDifEqu, double t, double y[], double k[], double vaccinated[], double h)
{
    
    double prior_vaccinated;                                // berechne % der Geimpften zum Zeitpunkt t - delta_t
    if ( t - delta_t <= 0 ) 
    { 
        prior_vaccinated = 0.0;
    }                  
    else 
    { 
        prior_vaccinated = interpolation(vaccinated, (t-delta_t)/h);
    }

    // hier:    nDifEqu = 6    
    k[0] = - beta_2(t) * y[1] * y[0] - v(t) * y[0] + epsilon * prior_vaccinated - zz * y[5] * y[0];     // 1. DGL: Susceptible (S)
    k[1] = + beta_2(t) * y[1] * y[0] - gammas * y[1] - psi * y[1] - zz * y[5] * y[1];                   // 2. DGL: Infected (I)
    k[2] = + gammas * y[1] - zz * y[5] * y[2];                                                          // 3. DGL: Recovered (R)
    k[3] = v(t) * y[0] - epsilon * prior_vaccinated - zz * y[5] * y[3];                                 // 4. DGL: Vaccinated (V)
    k[4] = psi * y[1];                                                                                  // 5. DGL: Dead (D)
    k[5] = z * y[3] +  zz * y[5] * ( y[0] + y[1] + y[2] + y[3] );                                       // 6. DGL: Zombified (Z)                                                                    
}                                   
