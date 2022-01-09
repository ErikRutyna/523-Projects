#include <stdio.h>
#include "headers.h"

void PrintMatrix(double* A, int N, int M){
    // Prints out NxM matrix A in rows -> columns

    //printf("Number of rows is: %i \n", N);
    //printf("Number of cols is: %i \n", M);

    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            printf("%f ", A[i * M + j]);
        }
        printf("\n");
    }
}

void SaveMatrix(double* A, int N, int M, FILE* fname){
    // Saves the matrix of N rows and M columns in CSV format to the filename
    int Gi;

    for(int i = 0; i < N; i ++) {
        for(int j = 0; j < M; j++) {
            Gi = i * M + j;

            fprintf(fname, "%1.6lf,", A[Gi]);
        }
        fprintf(fname, "\n");
    }
}