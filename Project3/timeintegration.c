#include "headers.h"
#include <stdio.h>
#include <math.h>

void CavityFlow(double N, double Re, int uW, double tol, double beta){

    // Initialize various flow parameters
    int i = 0;                  // Iteration number
    double t = 0;               // Time elapsed in simulation
    double h = 1 / N;           // Cell length in both x & y directions
    double dt;                  // Timestep - constant for each iteration
    FILE *vFile, *uFile, *rFile;// Pointers to output files
    double PPEE = 1; 
    double PPEtol = pow(10, -4);
    int Gi;

    // Check which timestep to use
    double dt1 = pow(h, 2) / (4 * 1 / Re);
    double dt2 = (4 * 1 / Re) / pow(uW, 2);

    if (dt1 <= dt2) {
        dt = beta * dt1;
    }
    else {
        dt = beta * dt2;
    }

    printf("\nInitialization of non-field values passed!\n");

    // Initialize the fields based on N
    double* U = Zeros(N+2, N+3);
    double* V = Zeros(N+3, N+2);
    double* P = Zeros(N+2, N+2);
    double* RHS;
    double* BigR;
    double* AU;
    double* E;
    printf("Initialization of state-fields passed!\n");

    // Initialize the fluxes
    double* F = Zeros(N+2, N+2);
    double* G = Zeros(N+2, N+2);
    double* HX = Zeros(N+3, N+3);
    double* HY = Zeros(N+3, N+3);
    printf("Initialization of other matrices passed!\n");

    // Generate output files
    uFile = fopen("uVelocity.csv", "w");
    vFile = fopen("vVelocity.csv", "w");
    rFile = fopen("residuals.csv", "w");
    printf("Initialization of output files passed!\n");

    // Residual term
    double error = 1;

    printf("\nGrid Size w/ ghost: %i x %i\n\n", (int) N+2, (int) N+2);
    while (error > tol){
        // Print a nice statement to terminal to track simulation progress
        if (i % 25 == 0){
            printf("Iteration: %i at time t = %f s - Residual: %f\n", i, t, error);
        }
        // Exit early if an instability causes residuals to increase
        if (error > 100){
            printf("Instability detected! The solution is no longer converging to steady state.\n");
            printf("Exiting simulation on iteration %i and time t = %f s.\n", i, t);
            break;
        }

        // Enforce ghost cell updates based on the BC's
        U = UGhostUpdate(U, uW, N+2, N+3);
        V = VGhostUpdate(V, 0, N+3, N+2);
        P = PGhostUpdate(P, N+2, N+2);

        // Calculate the fluxes for the projection step
        F = FCalc(F, U, 1/Re, N);
        G = GCalc(G, V, 1/Re, N);
        HX = HXCalc(HX, U, V, 1/Re, N);
        HY = HYCalc(HY, U, V, 1/Re, N);

        // Projects the velocity using the fluxes
        U = UProjection(U, F, HX, dt, N);
        V = VProjection(V, G, HY, dt, N);

        // Find RHS of PPE
        RHS = RHSCalc(U, V, dt, N);        
        P = Smoother(P, RHS, N, 50);

        // ================================
        // Failed Multigrid experimentation
        // ================================
        /*

        // Calculate the AU matrix
        AU = AUCalc(P, N);

        // Calculate error at each cell
        BigR = DivError(P, RHS, N);


        E = Zeros((N)+2, (N)+2);
        // Find the divergence-free pressure field using projected velocity
        //while(PPEE > PPEtol) {
        //E = PPESolve3(E, Restriction(BigR, N, N), Restriction(AU, N, N), N/2);
        E = Smoother(E, BigR, N, 1000);

        // Prolongate the error
        //E = DirectProlongation(E, N/2, N/2);
        
        // Add the error back into the pressure
        for(int i = 1; i < (N+1); i++) {
            for(int j = 1; j < (N+1); j++) {
                Gi = (int) i * (N+2) + j;

                P[Gi] += E[Gi];
            }
        }            
        
        // Some final smoothing
        //P = Smoother(P, RHS, N, 5);

        // Check if we have "solved" the PPE to our tolerance
        //PPEE = GSError(P, U, V, dt, N);
        //}
        // ================================
        // END OF MULTIGRID EXPERIMENTATION
        // ================================
        */

        // Correct the velocity fields to n+1
        U = UCorrection(U, P, dt, N);
        V = VCorrection(V, P, dt, N);

        // Calculate the residual
        error = ResidualCalc(P, F, G, HX, HY, N);
        
        /*
        printf("\nU Matrix \n");
        PrintMatrix(U, N+2, N+3);
        printf("\nV Matrix \n");
        PrintMatrix(V, N+3, N+2);
        printf("\nPressure Matrix \n");
        PrintMatrix(P, N+2, N+2);
        printf("\nPressure2 Matrix \n");
        PrintMatrix(P2, N+2, N+2);
        printf("\nF Flux Matrix\n");
        PrintMatrix(F, N+2, N+2);
        printf("\nG Flux Matrix\n");
        PrintMatrix(G, N+2, N+2);
        printf("\nHX Flux Matrix\n");
        PrintMatrix(HX, N+3, N+3);
        printf("\nHY Flux Matrix\n");
        PrintMatrix(HY, N+3, N+3);*/

        // Update residual, simulation time, iteration number and save residual
        //error = 0;
        i++;
        t += dt;
        fprintf(rFile, "%lf", error);
        fprintf(rFile, "\n");
    }
    // Save the velocity fields to outputs files
    SaveMatrix(U, (int) N+2, (int) N+3, uFile);
    SaveMatrix(V, (int) N+3, (int) N+2, vFile);

}