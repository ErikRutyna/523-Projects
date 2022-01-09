#include <stdlib.h>
#include <stdio.h>
#include "headers.h"

double* Zeros(int N, int M){
    // Initializes a zero matrix for N rows and M columns
    int sizeM = N * M;

    // Calloc automatically sets values to zero
    double* A = (double*) calloc(sizeM, sizeof(double));
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            A[(i * M + j)] = 0.0;
        }
    }

    return A;
}

double* UGhostUpdate(double* U, int uW, int N, int M) {
    // Updates the ghost values for the x-velocity matrix

    int Ti;
    int Bi;
    int Li;
    int Ri;

    // Top & bottom walls
    for(int j = 1; j < (M-1); j++){
        Ti = (N-1) * M + j;
        Bi = j;
        U[Ti] = 2 * uW - U[Ti - M]; // Top
        U[Bi] = -U[Bi + M]; // Bottom
    }

    // Left & right walls
    for(int i = 1; i < (N-1); i++){
        Li = i * M;
        Ri = (i+1) * M - 1;
        U[Li] = U[Li + 2]; // Left
        U[Ri] = U[Ri - 2]; // Right
    }

    return U;
}

double* VGhostUpdate(double* V, int uW, int N, int M) {
    // Updates the ghost values for the y-velocity matrix

    int Ti;
    int Bi;
    int Li;
    int Ri;

    // Top & bottom walls
    for(int j = 1; j < (M-1); j++) {
        Ti = (N-1) * M + j;
        Bi = j;
        V[Ti] = V[Ti - 2 * M]; // Top
        V[Bi] = V[Bi + 2 * M]; // Bottom
    }

    // Left & right walls
    for(int i = 2; i < (N-2); i++) {
        Li = i * M;
        Ri = (i + 1) * M - 1;
        V[Li] = 2 * uW -V[Li + 1]; // Left
        V[Ri] = 2 * uW -V[Ri - 1]; // Right
    }

    return V;
}

double* PGhostUpdate(double* P, int N, int M){
    // Updates the ghost values for the pressure matrix
    //printf("Dimensions - Rows: %i\tCols: %i\n", N, M);
    // Top & bottom walls
    for(int j = 1; j < M; j++){
        P[(N-1) * M + j] = P[(N-2) * M + j];// Top
        P[j] = P[M + j]; // Bottom
    }    

    // Left & right walls
    for(int i = 1; i < N; i++){
        P[i * M] = P[i * M + 1];
        P[(i+1) * M - 1] = P[(i+1) * M - 2];
    }

    return P;
}