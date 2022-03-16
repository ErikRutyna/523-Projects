#include "headers.h"
#include <stdio.h>
#include <math.h>

void CavityFlow(double N, double Re, int uW, double tol, double beta){

    // Initialize various flow parameters
    int i = 0;                          // Iteration number
    double t = 0;                       // Time elapsed in simulation
    double h = 1 / N;                   // Cell length in both x & y directions
    double dt;                          // Timestep - constant for each iteration
    FILE *vFile, *uFile, *rFile, *pfile;// Pointers to output files
    int Gi;
    FILE *pfile2;

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
    //double* P2 = Zeros(N+2, N+2);
    double* RHS;
    double* error;
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
    pfile = fopen("pressure.csv", "w");
    pfile2 = fopen("pressure2.csv", "w");
    printf("Initialization of output files passed!\n");

    // Residual term
    double residual = 1;

    printf("\nGrid Size w/ ghost: %i x %i\n\n", (int) N+2, (int) N+2);
    while (residual > tol){
        // Print a nice statement to terminal to track simulation progress
        if (i % 25 == 0){
            printf("Iteration: %i at time t = %f s - Residual: %f\n", i, t, residual);
        }
        // Exit early if an instability causes residuals to increase
        if (residual > 100){
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

        // Start w/ 2x Multigrid, then make it recursive
        // =========================
        // Multigrid experimentation
        // =========================

        //printf("\nPressure Matrix \n");
        //PrintMatrix(P, N+2, N+2);
        //printf("\nPressure (Multigrid) Matrix \n");
        //PrintMatrix(P2, N+2, N+2);

        // Find RHS of PPE
        RHS = RHSCalc(U, V, dt, N);

        // Pre-smooth a few iterations
        P = Smoother(P, RHS, N, 10);

        //printf("\nPressure Matrix \n");
        //PrintMatrix(P, N+2, N+2);
        //printf("\nPressure (Multigrid) Matrix \n");
        //PrintMatrix(P2, N+2, N+2);

        
        // Compute the initial residual for Au=b, where r = b - Au
        error = PPEResidual(P, RHS, N);

        // Restrict the error to be on the corase grid
        error = Restriction(error, N, N);
        
        // Make the "e" matrix
        E = Zeros((N/2)+2, (N/2)+2);

        
        // Smooth on the restricted domain
        E = Smoother(E, error, N/2, 10);
        
        // Prolongate E
        E = DirectProlongation(E, N/2, N/2);

        // Add the error E back into P
        P = AddMatrix(P, E, N);
        
        // Post Smoothing Operation
        P = Smoother(P, RHS, N, 5);

        //printf("\nPressure Matrix \n");
        //PrintMatrix(P, N+2, N+2);
        //printf("\nPressure (Multigrid) Matrix \n");
        //PrintMatrix(P2, N+2, N+2);

        // ================
        // END OF MULTIGRID 
        // ================

        // Correct the velocity fields to n+1
        U = UCorrection(U, P, dt, N);
        V = VCorrection(V, P, dt, N);

        // Calculate the residual
        residual = ResidualCalc(P, F, G, HX, HY, N);
        
        /*
        printf("\nU Matrix \n");
        PrintMatrix(U, N+2, N+3);
        printf("\nV Matrix \n");
        PrintMatrix(V, N+3, N+2);
        printf("\nPressure Matrix \n");
        PrintMatrix(P, N+2, N+2);
        printf("\nPressure (Multigrid) Matrix \n");
        PrintMatrix(P2, N+2, N+2);
        printf("\nF Flux Matrix\n");
        PrintMatrix(F, N+2, N+2);
        printf("\nG Flux Matrix\n");
        PrintMatrix(G, N+2, N+2);
        printf("\nHX Flux Matrix\n");
        PrintMatrix(HX, N+3, N+3);
        printf("\nHY Flux Matrix\n");
        PrintMatrix(HY, N+3, N+3);
        */

        // Update residual, simulation time, iteration number and save residual
        //residual = 0;
        i++;
        t += dt;
        fprintf(rFile, "%lf", residual);
        fprintf(rFile, "\n");
    }
    // Save the velocity fields to outputs files
    SaveMatrix(U, (int) N+2, (int) N+3, uFile);
    SaveMatrix(V, (int) N+3, (int) N+2, vFile);
    SaveMatrix(P, (int) N, (int) N, pfile);
    //SaveMatrix(P2, (int) N, (int) N, pfile);

}