
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
        fprintf( file_pointer, "%+.10e    %+.10e    %+.10e    %+.10e    %+.10e    %+.10e    %+.10e\n" , numb_list1[i], numb_list2[i], numb_list3[i], numb_list4[i], numb_list5[i], numb_list6[i], numb_list7[i] );
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

    // returns pointer-pointer to matrix of  dimension [x_dim, y_dim], initialised with "initial_value"
double **create_2d_matrix (int x_dim, int y_dim, double initial_value)
{
    double *linear_matrix      = (double *)malloc(x_dim * y_dim * sizeof(double));       
    double **quadratic_matrix  = (double **)malloc(x_dim * sizeof(double *));  

    if( linear_matrix == NULL || quadratic_matrix == NULL ){ printf("Matrix konnte nicht angelegt werden!\n"); }    

        // set initial value
    for( int j = 0 ; j < x_dim * y_dim ; j++ )
    {
        linear_matrix[j] = initial_value;       
    }
        // introduce [x][y]-referencing
    for( int i = 0 ; i < x_dim ; i++ )
    {
        quadratic_matrix[i] = linear_matrix + i * y_dim;        
    }

    return quadratic_matrix;
}

    // frees the allocated memory of a 2D matrix
int free_2d_matrix( double **matrix )
{
    free(matrix[0]);        // free "linear_matrix"
    free(matrix);           // free "quadratic_matrix"
    return 0;
}

    // speichert Matrix in txt
int save_matrix( double **matrix, int x_dim, int y_dim, char* save_as )
{
    FILE *file_pointer;                             // open file in write-mode and save its address
    file_pointer = fopen(save_as, "w"); 

    for( int y = 0 ; y < y_dim ; y++ )
    {
        for( int x = 0 ; x < x_dim ; x++ )
        {
            fprintf(file_pointer, "%+.8e ", matrix[x][y]);
        }
        fprintf(file_pointer, "\n");
    }

    fclose(file_pointer);
    return 0;
}

    // calculate the average of the absolute difference between array1 and array2 (elementwise)
double average_diff( double array1[], double array2[] )
{
        // assume: array1[0] = array2[0] = length(array1) = length(array2)
    double diff = 0.0;  

    for(int i = 0; i < array1[0]; i++ )     //
    {
        diff = diff + fmin( fabs( array1[i+1] - array2[i+1] ), fabs( array1[i+1] + array2[i+1] )  );
    }

    return diff/array1[0];
}

    // translates all values of an array into the range [0, limit_up] the way modulo should
int modulus( double array[], double limit_up )
{
    for(int i = 0; i < array[0]; i++)
    {
        if( array[i+1] > 0)
        {
            array[i+1] = fmod(array[i+1], limit_up ) ;
        } else {
            array[i+1] = fmod(array[i+1] , limit_up ) + limit_up;
        }
    }
    return 0;
}

    // translates all values of an array into the range [0, limit_up] the way modulo should
int modulus_array( double array[], double limit_low, double limit_up )
{
    for(int i = 0; i < array[0]; i++)
    {
        if( array[i+1] >= 0)
        {
            array[i+1] = fmod(array[i+1] - limit_low, limit_up-limit_low ) + limit_low ;
        } else {
            array[i+1] = fmod(array[i+1] + limit_low , limit_up-limit_low ) - limit_low;
        }
    }
    return 0;
}

    // translates a value into the range [0, limit_up] the way modulo should
double modulus_s( double value, double limit_up )
{
    double result;
    if( value > 0)
    {
        result = fmod(value, limit_up ) ;
    } else {
        result = fmod(value , limit_up ) + limit_up;
    }

    return result;
}

    // calculates the result array of the linearcombination of (# = "array_numb") given arrays of size "array_length" and their scalar coefficients. Returns result array
double *linear_comb_arrays( double** arrays, double* coeffs, double *result, int array_length, int array_numb )
{
    for( int i = 0; i < array_numb; i++ )
    {
        for( int j = 0; j < array_length; j++ )
        {
            if( i == 0 ){ result[j] = 0.0; }
            result[j] = result[j] + coeffs[i] * arrays[i][j];
        }
    }
    return result;
}

    // fill the given array of length "length" with "0"
int zeros( double array[], int length )
{
    for( int i = 0; i < length; i++ )
    {
        array[i] = 0.0;
    }
    return 0;
}

    // copy the array "array" into the array "copy" of same length "length"
int copy_array( double array[], double copy[], int length)
{
    for( int i = 0; i < length; i++ )
    {
        copy[i] = array[i];
    }
    return 0;
}

    // prints out an array 
int print_array( double array[] )
{
    printf("[ ");
    for( int i = 0; i < array[0]+1; i++)
    {
        printf("%+.2lf, ", array[i]);
    }
    printf("] \n");

    return 0;
}

    // erase last line printed in shell)
int erase_last_line(void)
{        
    printf("\r");
    printf("\x1b[1F"); // Move to beginning of previous line
    printf("\x1b[2K"); // Clear entire line
    return 0;
}

    // fills up array with lineary spaced values between "min_range" and "max_range"
int fill_linspace( double array[], double min_range, double max_range, int length )
{
    array[0] = length;
    for( int j = 0; j < length; j++)
    {
        array[j+1] = min_range + (max_range-min_range) * j / (double)(length-1);
    }
    return 0;
}