#pragma once
#include <stdio.h>

// Time Integration
void CavityFlow(double N, double Re, int uW, double tol, double beta);

/*======================
===== BC Functions =====
======================*/

// Initialized a zero matrix for rows/cols
double* Zeros(int N, int M);

// Ghost cell enforcement
double* UGhostUpdate(double* U, int uW, int N, int M);
double* VGhostUpdate(double* V, int uW, int N, int M);
double* PGhostUpdate(double* P, int N, int M);

/* =======================
===== Flux Functions =====
========================*/

// Flux caluclator function for entire field
double* FCalc(double* F, double* U, double nu, int N);
double* GCalc(double* G, double* V, double nu, int N);
double* HXCalc(double* HX, double* U, double* V, double nu, int N);
double* HYCalc(double* HY, double* U, double* V, double nu, int N);

// Limiter
double SMART(double phiL, double phi, double phiR);


/*=============================
===== Velocity Projection =====
=============================*/

// Projects the velocity field to n+1/2 using the fluxes
double* UProjection(double* U, double* F, double* HX, double dt, double N);
double* VProjection(double* V, double* G, double* HY, double dt, double N);


// Solves for the pressure field that makes the velocity field divergence free
double* PPESolve(double* P, double* U, double* V, double dt, double N); // SOR-GS
double* PPESolve2(double* P, double* U, double* V, double dt, double N); // SOR-RBGS
double  GSError(double* P, double* U, double* V, double dt, double N); 

// Corrects the velocity field to n+1 using the divergence free pressure field
double* UCorrection(double* U, double* P, double dt, double N);
double* VCorrection(double* V, double* P, double dt, double N);

/*===================
===== Multigrid =====
===================*/

// Main multigrid driver
double* PPESolve3(double* E, double* R, double* AU, double N);

// Restriction
double* Restriction(double* A, int N, int M); // Restricts down from NxM to (N/2)x(M/2)

// Prolongation/Injection
double* DirectProlongation(double* A, int N, int M); // Prolongates from NxM to (2N)x(2M) via direct injection

// Divergence calc
double* RHSCalc(double* U, double* V, double dt, int N); // div(V)/dt for the fields

// SOR RB-GS Smoother
double* Smoother(double* P, double* divV, double N, int Z);

// Residual vector calc
double* PPEResidual(double* P, double* divV, int N);

// The mock-laplacian thing
double* AUCalc(double* A, int N); 


/*========================
===== Misc Debugging =====
========================*/

// Add two matrices
double* AddMatrix(double* K, double* L, int N);

// Calculates the residual at the timestep
double ResidualCalc(double* P, double* F, double* G, double* HX, double* HY, double N);

// Prints a matrix of N rows and M cols to terminal
void PrintMatrix(double* A, int N, int M);

// Saves a matrix to the current directory
void SaveMatrix(double* A, int N, int M, FILE* fname);