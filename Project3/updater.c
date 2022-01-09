#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "headers.h"

double* UProjection(double* U, double* F, double* HX, double dt, double N){
    // Updates the x-velocity to n+1/2
    int Gi;
    int N2 = N;
    double Fx;
    double HXy;
    double h = 1 / (double) N;

    for(int i = 1; i < (N+1); i++){
        for(int j = 2; j < (N+1); j++){
            Gi = (int) i * (N+3) + j;

            Fx = (F[Gi - i] - F[Gi -i - 1]) / h;

            HXy = (HX[Gi + (N2+3)] - HX[Gi]) / h;

            U[Gi] -= dt * (Fx + HXy);
        }
    }
    return U;
}

double* VProjection(double* V, double* G, double* HY, double dt, double N){
    // Updates the y-velocity to n+1/2
    int Gi;
    int N2 = N;
    double h = 1 / N;
    double Gy;
    double HYx; 

    for(int i = 2; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++){
            Gi = (int) i * (N+2) + j;

            Gy = (G[Gi] - G[Gi - (N2+2)]) / h;

            HYx = (HY[Gi + i + 1] - HY[Gi + i]) / h;

            V[Gi] -= dt * (Gy + HYx);
        }
    }
    return V;
}

double* PPESolve(double* P, double* U, double* V, double dt, double N){
    // Corrects the pressure field to make our velocity field divergence free
    // by using SOR RB-GS w/ recursive multigrid

    int Gi;
    int N2 = (int) N;
    double w = 1.6;
    double h = 1 / N;
    double Pterm;
    double Vterm;
    double divError = 1.0;
    double tolerance = pow(10, -4);

    //for(int Z = 0; Z < 100; Z++) {
    while(divError > tolerance) {
        for(int i = 1; i < (N+1); i++){
            for(int j = 1; j < (N+1); j++){
                Gi = (int) i * (N+2) + j;

                Pterm = (P[Gi + 1] + P[Gi - 1] + P[Gi + (N2+2)] + P[Gi - (N2+2)]);

                Vterm = h / dt * (U[Gi + i + 1] - U[Gi + i] + V[Gi + (N2+2)] - V[Gi]);

                P[Gi] = w / 4 * (Pterm - Vterm) + (1 - w) * P[Gi];
            }
        }
        divError = GSError(P, U, V, dt, N);
        P = PGhostUpdate(P, (N2+2), (N2+2));
    }
    //printf("GS Error: %1.15f\n", GSError(P, U, V, dt, N));
    return P;
}

double* PPESolve2(double* P, double* U, double* V, double dt, double N){
    // Does SOR RB-GS w/ an exit condition to solve for the PPE instead of SOR-GS
    // and a fixed number of iterations

    int Gi;
    int N2 = (int) N;
    double w = 1.6;
    double h = 1 / N;
    double Pterm;
    double Vterm;
    double divError = 1.0;
    double tolerance = pow(10, -6);

    while (divError > tolerance) {
    //for(int Z = 0; Z < 100; Z++) {
        // Red Update
        for(int i = 1; i < (N+1); i++){
            for(int j = 1; j < (N+1); j++){
                if (((i + j) % 2) == 0) {
                    Gi = (int) i * (N+2) + j;
                    //printf("Gi: %i\n", Gi);
                    Pterm = (P[Gi + 1] + P[Gi - 1] + P[Gi + (N2+2)] + P[Gi - (N2+2)]);

                    Vterm = h / dt * (U[Gi + i + 1] - U[Gi + i] + V[Gi + (N2+2)] - V[Gi]);

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

                    Vterm = h / dt * (U[Gi + i + 1] - U[Gi + i] + V[Gi + (N2+2)] - V[Gi]);

                    P[Gi] = w / 4 * (Pterm - Vterm) + (1 - w) * P[Gi];
                }
            }
        }
        // Ghost Update
        P = PGhostUpdate(P, (N2+2), (N2+2));

        // Calculate RB-GS Error
        divError = GSError(P, U, V, dt, N);
    }
    //printf("RB-GS Error: %1.15f\n", divError);
    return P;

}

double* UCorrection(double* U, double* P, double dt, double N){
    // Uses the pressure field to correct the x-velocity

    int Gi;
    int N2 = (int) N;
    double h = 1 / N;

    for(int i = 1; i < (N+1); i++) {
        for(int j = 2; j < (N+1); j++) {
            Gi = i * (N2+3) + j;

            U[Gi] -= dt / h * (P[Gi - i] - P[Gi - i - 1]);
        }
    }
    return U;
}

double* VCorrection(double* V, double* P, double dt, double N){
    // Uses the pressure field to correct the y-velocity

    int Gi;
    int N2 = (int) N;
    double h = 1 / N;

    for(int i = 2; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++){
            Gi = i * (N2+2) + j;

            V[Gi] -= dt / h * (P[Gi] - P[Gi - (N2+2)]);
        }
    }
    return V;
}
