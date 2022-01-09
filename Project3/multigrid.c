#include "headers.h"
#include <math.h>
#include <stdio.h>

double* PPESolve3(double* E, double* R, double* AU, double N){
    // Recursive Multigrid
    /*
    if( (int) N == 16) {
        E = Smoother(E, R, N, 100);
        return E;
    }
    */
    // Pre-smooth
    E = Smoother(E, R, N, 20);

    int Gi = 0;
    // Restriction
    for(int i = 1; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++) {
            Gi = (int) i * (N+2) + j;

            R[Gi] -= AU[Gi]; 
        }
    }
    double* rR = Restriction(R, N, N);

    // Create new matrices for next level
    AU = Zeros((N/2)+2, (N/2)+2);
    double* E2 = Zeros((N/2)+2, (N/2)+2);

    E2 = Smoother(E2, rR, (N/2), 100);

    // Enter recursion
    //E2 = PPESolve3(E2, R, AU, N/2);

    // Prolongate back out the error
    double* pE2 = DirectProlongation(E2, N/2, N/2);

    // Add it to the current level's error
    /*Gi = 0;
    for(int i = 1; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++) {
            Gi = (int) i * (N+2) + j;

            E[Gi] += E2[Gi]; 
        }
    }*/

    // Post-smoothing
    E = Smoother(pE2, R, N, 20);

    return E;
}

double* Restriction(double* A, int N, int M){ 
    // Restricts down from NxM to (N/2)x(M/2)

    //printf("\nUnrestricted Size: %i x %i\n", N+2, M+2);

    int rCol = M / 2 + 2;
    int rRow = N / 2 + 2;
    double* Anew = Zeros(rRow, rCol); // Restricted size + ghost

    //printf("\nError matrix after initialization\n");
    //PrintMatrix(Anew, rCol, rRow);
    //printf("Restricted Size: %i x %i\n", rCol, rRow);
    //PrintMatrix(Anew, rCol, rRow);
    
    int iBig = 1;
    int jBig = 1;

    int GiSmall;
    int GiBig;

    for(int iSmall = 1; iSmall < (rRow-1); iSmall++){
        for(int jSmall = 1; jSmall < (rCol-1); jSmall++){
            GiSmall = iSmall * (rCol) + jSmall;
            GiBig = iBig * (M+2) + jBig;

            //printf("Small Index: %i\t Big Index: %i\n", GiSmall, GiBig);
            //printf("Pts: %f\t %f\t %f\t %f\n", A[GiBig], A[GiBig+1], A[GiBig + M + 2], A[GiBig + (M+2) + 1]);
            Anew[GiSmall] = (A[GiBig] + A[GiBig + 1] + A[GiBig + (M+2)] + A[GiBig + (M+2) + 1]) / 4;

            jBig += 2;
        }
        iBig += 2;
        jBig = 1;
    }
    //printf("After restriction: \n");
    //PrintMatrix(Anew, rCol, rRow);
    return Anew;
}

double* DirectProlongation(double* A, int N, int M){ 
    // Prolongates from NxM to (2N)x(2M) via direct injection

    //printf("Restricted Size: %i x %i\n", N, M);
    int pCol = 2 * M + 2;
    int pRow = 2 * N + 2;
    double* Anew = Zeros(pRow, pCol);
    //printf("\nUnrestricted Size: %i x %i\n", pRow-2, pCol-2);

    int iBig = 1;
    int jBig = 1;

    int GiSmall;
    int GiBig;

    for(int iSmall = 1; iSmall < (N+1); iSmall++){
        for(int jSmall = 1; jSmall < (M+1); jSmall++){
            GiSmall = iSmall * (M+2) + jSmall;
            GiBig = iBig * (pCol) + jBig;

            //printf("Small Index: %i\t Big Indices: %i\t %i\t %i\t %i\n", GiSmall, GiBig, GiBig+1, GiBig+pCol, GiBig+pCol+1);
            
            Anew[GiBig] = A[GiSmall];
            Anew[GiBig + 1] = A[GiSmall];
            Anew[GiBig + pCol] = A[GiSmall];
            Anew[GiBig + pCol + 1] = A[GiSmall];

            jBig += 2;
        }
        iBig += 2;
        jBig = 1;
    }

    return Anew;
}

double* Smoother(double* P, double* divV, double N, int Z) {
    // RB-GS Smoothering for Z-number of iterations

    int Gi;
    int N2 = (int) N;
    double w = 1.6;
    double h = 1.0 / N;
    double Pterm;
    double Vterm;

    for(int Zi = 0; Zi < Z; Zi++) {
        // Red Update
        for(int i = 1; i < (N+1); i++){
            for(int j = 1; j < (N+1); j++){
                if (((i + j) % 2) == 0) {
                    Gi = (int) i * (N+2) + j;

                    Pterm = (P[Gi + 1] + P[Gi - 1] + P[Gi + (N2+2)] + P[Gi - (N2+2)]);

                    Vterm = h * h * divV[Gi];

                    P[Gi] = w / 4 * (Pterm - Vterm) + (1 - w) * P[Gi];
                }
            }
        }
        // Black Update
        for(int i = 1; i < (N+1); i++){
            for(int j = 1; j < (N+1); j++){
                if (((i + j) % 2) == 1) {
                    Gi = (int) i * (N+2) + j;
                    
                    Pterm = (P[Gi + 1] + P[Gi - 1] + P[Gi + (N2+2)] + P[Gi - (N2+2)]);

                    Vterm = h * h * divV[Gi];

                    P[Gi] = w / 4 * (Pterm - Vterm) + (1 - w) * P[Gi];
                }
            }
        }
        // Ghost Update
        P = PGhostUpdate(P, (N2+2), (N2+2));
    }
    return P;
}

double* RHSCalc(double* U, double* V, double dt, int N){
    // Solves for the RHS of the PPE at each pressure cell

    int Gi;
    double* divV = Zeros(N+2, N+2);
    double h = 1.0 / N;

    for(int i = 1; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++) {
            Gi = (int) i * (N+2) + j;
            
            divV[Gi] = (U[Gi + i + 1] - U[Gi + i] + V[Gi + (N+2)] - V[Gi]) / (h * dt);
        }
    }
    return divV;
}

double* DivError(double* P, double* divV, double N){
    // Calculates residual error "r" in r = div(V) - Laplace(P)

    double* r = Zeros(N+2, N+2);
    double h = 1.0 / N;
    int Gi;

    for(int i = 1; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++) {
            Gi = (int) i * (N+2) + j;

            r[Gi] = divV[Gi] - P[Gi];
        }
    }
    return r;
}

double* AUCalc(double* A, int N){

    double* AU = Zeros(N+2, N+2);
    int Gi;
    double h = 1.0 / N;

    for(int i = 1; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++) {
            Gi = (int) i * (N+2) + j;

            AU[Gi] = A[Gi + 1] + A[Gi - 1] + A[Gi + (N+2)] + A[Gi - (N+2)] - 4 * A[Gi] / (h * h);
        }
    }
    return AU;
}
