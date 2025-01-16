// "gcc -o ./out/made main.c my_RK4.c solve_DGLs.c DGLs.c tools.c header.h -lm"           (kompilieren)
// "./out/made"                                                                           (ausfuehren)
// Tobias Neuhoff, Noah Reinhardt                                                         (Autoren)

// alternativ: "make"                                                                     (kompilieren & ausfuehren)


#include "header.h"

    // gibt Inhalt einer double-liste aus
int display_numb_list( double* numb_list)               // should read out all elements of the list
{
    int length = (int)numb_list[0];                     // again, the first element of the list is chosen to be its length
    printf("Die Laenge der Liste ist %d. \n", length ); 
    
    for ( int i = 1 ; i < length+1 ; i++ )              // print all elements bit by bit
        {
            printf("%lf \n", numb_list[i] );
        }
    
    return 0;
}

    // speichert 7 double-Listen spaltenweise in einer txt
int save_numb_list7 ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, char* save_as )
{
    int length = (int)numb_list1[0];                 // extract length from first element of list
   
    FILE *file_pointer;                             // open file in write-mode and save its address
    file_pointer = fopen(save_as, "w"); 
    
    for( int i = 0 ; i < length+1 ; i++ )
    {
        fprintf( file_pointer, "%lf    %lf    %lf    %lf    %lf    %lf    %lf\n" , numb_list1[i], numb_list2[i], numb_list3[i], numb_list4[i], numb_list5[i], numb_list6[i], numb_list7[i] );
    }

    fclose(file_pointer);
    return 0;
}

    // Interpoliert linear zwischen den Werten eines Arrays --> Aufruf mit nicht-ganzzaligen Indizes möglich
double interpolation( double array[], double index )
{
    int length = (int)array[0];
    for(int i = 1; i <= length ; i++ )
    {
        if( ((double)i <= index) && (index < (double)(i+1)) )
        {
            return ( array[i] + ( array[i+1] - array[i] )*(index - i) );
        }
    }
    if ((0.0 <= index) && (index < 1.0))
    {
        return array[1];
    } 
    else
    {
        printf("Etwas ist schief gelaufen bei der Interpolation einer arrays!\n");
        return 0.0;
    }
}


    // erstellt einen double-array voller Nullen gewünschter Länge
double *create_null( int length )
{
    double *array = (double *)malloc( (length+1) * sizeof(double) );
    array[0] = (double)length;
    for(int i ; i < length ; i++)
    {
        array[i+1] = 0.0;
    }
    return array;
}