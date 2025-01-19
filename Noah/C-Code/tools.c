
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

    // saves 7 double-lists by columns in a .txt
int save_numb_list7 ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, char* save_as )
{
    int length = (int)numb_list1[0];                 // extract length from first element of list
   
    FILE *file_pointer;                             // open file in write-mode and save its address
    file_pointer = fopen(save_as, "w"); 
    
    for( int i = 0 ; i < length+1 ; i++ )
    {
        fprintf( file_pointer, "%+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e\n" , numb_list1[i], numb_list2[i], numb_list3[i], numb_list4[i], numb_list5[i], numb_list6[i], numb_list7[i] );
    }

    fclose(file_pointer);
    return 0;
}

    // saves 7 double-lists by columns in a .txt
int save_numb_list9 ( double* numb_list1, double* numb_list2, double* numb_list3, double* numb_list4, double* numb_list5, double* numb_list6, double* numb_list7, double* numb_list8, double* numb_list9, char* save_as )
{
    int length = (int)numb_list1[0];                 // extract length from first element of list
   
    FILE *file_pointer;                             // open file in write-mode and save its address
    file_pointer = fopen(save_as, "w"); 
    
    for( int i = 0 ; i < length+1 ; i++ )
    {
        fprintf( file_pointer, "%+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e    %+.8e\n" , numb_list1[i], numb_list2[i], numb_list3[i], numb_list4[i], numb_list5[i], numb_list6[i], numb_list7[i], numb_list8[i], numb_list9[i] );
    }

    fclose(file_pointer);
    return 0;
}

    // creates an array full of "0" of given length
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

    // copies the entries of the shorter array "old_array" into the empty longer array "new_array" (from left to right)
int merge_arrays( int old_length, double old_array[], int new_length, double new_array[] )
{
    new_array[0] = new_length;

    for( int i = 0; i < old_length; i++ )
    {
        new_array[i+1] = old_array[i];
    }
    for( int i = old_length; i < new_length; i++ )
    {
        new_array[i+1] = 0.0;
    }

    return 0;
}