#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "headers.h"

// Set of primary flux calculating functions for velocity projection
double* FCalc(double* F, double* U, double nu, int N){
    // Computes x-flux for all cells
    
    double q;
    double phi;
    double ux;
    double h = 1 / (double) N;
    int N2 = (int) N;
    int Gi;

    for(int i = 1; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++){
            Gi = i * (N+2) + j;

            // Find the velocity
            q = (U[Gi + i] + U[Gi + i + 1]) / 2;
            
            // Apply SMART
            if (q > 0) {
                phi = SMART(U[Gi + i - 1], U[Gi + i], U[Gi + i + 1]);
            }
            if (q < 0) {
                phi = SMART(U[Gi + i + 2], U[Gi + i + 1], U[Gi + i]);
            }
            if (q == 0) {
                phi = 0;
            }            
            // Grab our local derivative
            ux = (U[Gi + i + 1] - U[Gi + i]) / h;

            // Update flux matrix
            F[Gi] = q * phi - nu * ux;
        }
    }
    return F;
}

double* GCalc(double* G, double* V, double nu, int N){
    // Computes y-flux for all cells

    double q;
    double phi;
    double vy;
    double h = 1 / (double) N;
    int N2 = (int) N;
    int Gi;

    for(int i = 1; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++){
            Gi = i * (N2 + 2)+ j;

            // Find the velocity
            q = (V[Gi + (N2+2)] + V[Gi]) / 2;

            // Apply SMART
            if (q > 0) {
                phi = SMART(V[Gi - (N2+2)], V[Gi], V[Gi + (N+2)]);
            }
            else if (q < 0) {
                phi = SMART(V[Gi + 2 * (N2+2)], V[Gi + (N2+2)], V[Gi]);
            }
            else {
                phi = 0;
            }
            // Compute the local derivative
            vy = (V[Gi + (N2+2)] - V[Gi]) / h;

            // Update flux matrix
            G[Gi] = q * phi - nu * vy;
        }
    }
    return G;
}

double* HXCalc(double* HX, double* U, double* V, double nu, int N){
    // Computes cross-directional HX fluxes for all nodes
    
    double q;
    double phi;
    double uy;
    double h = 1 / (double) N;
    int M = N + 2;
    int G;

    for(int i = 1; i < M; i++){
        for(int j = 1; j < M; j++){
            G =  i * (N+3) + j;

            // Find the velocity
            q = (V[G - i - 1] + V[G - i]) / 2;
            
            // Apply SMART
            if (q > 0) {
                phi = SMART(U[G - 2 * (N+3)], U[G - (N+3)], U[G]);
            }
            if (q < 0) {
                phi = SMART(U[G + (N+3)], U[G], U[G - (N+3)]);
            }
            if (q == 0) {
                phi = 0;
            }

            // Grab our local derivative
            uy = (U[G] - U[G - (N+3)]) / h;

            // Update flux matrix
            HX[G] = (double) q * phi - nu * uy;
        }
    }
    return HX;
}

double* HYCalc(double* HY, double* U, double* V, double nu, int N){
    // Computes cross-directional HY fluxes for all nodes

    double q;
    double phi;
    double vx;
    double h = 1 / (double) N;
    int M = N + 2;
    int G;

    for(int i = 1; i < M; i++){
        for(int j = 1; j < M; j++){
            G = i * (N+3) + j;
            
            // Find the velocity
            q = (U[G] + U[G - (N+3)]) / 2;
            
            // Apply SMART
            if (q > 0) {
                phi = SMART(V[G - i - 2], V[G - i - 1], V[G - i]);
            }
            if (q < 0) {
                phi = SMART(V[G - i + 1], V[G - i], V[G - i - 1]);
            }
            if (q == 0) {
                phi = 0;
            }

            // Grab our local derivative
            vx = (V[G - i] - V[G - i - 1]) / h;

            // Update flux matrix
            HY[G] = q * phi - nu * vx;
        }
    }
    return HY;
}

// Limiter
double SMART(double phiL, double phi, double phiR){
    // Computes SMART limiter for the set of values

    double phiHalf;
    double phiHat;

    // Edge case check
    if (fabs(phiL - phiR) < 0.00000001){
        return phi;
    }

    phiHat = (phi - phiL) / (phiR - phiL);
    //printf("phiL: %f\t, Phi: %f\t phiR: %f\t PhiHat: %f\n", phiL, phi, phiR, phiHat);

    if ((phiHat <= 0.0) || (phiHat >= 1.0)){
        phiHalf = phiHat;
        //printf("Check passed. \t phiHalf: %f\t phiHat: %f \n", phiHalf, phiHat);
    }
    else if ((phiHat > 0.0) && (phiHat <= (1.0/6.0))){
        phiHalf = 3.0 * phiHat;
    }
    else if ((phiHat > (1.0/6.0)) && (phiHat <= (5.0/6.0))){
        phiHalf = 3.0 / 8.0 * (2.0 * phiHat + 1.0);
        //printf("Phihat: %f\t Phihalf: %f\t PhiHat Supposed: %f\n",phiHat, phiHalf, (3.0 / 8.0 * (2.0 * phiHat + 1.0)));
    }
    else {
        phiHalf = 1.0;
    }

    //printf("PhiHalf: %f\t, phiL %f\t phi %f\t phiR %f\n", phiHalf, phiL, phi, phiR);
    phiHalf = phiL + (phiR - phiL) * phiHalf;

    //printf("PhiHalf Return: %f\n", phiHalf);
    return phiHalf;
}

// Residual Calculation
double ResidualCalc(double* P, double* F, double* G, double* HX, double* HY, double N){
    // Calculates the residual for the system at the current timestep

    int Gi;
    int N2 = (int) N;
    double h = 1 / N;
    double L1Normx = 0;
    double L1Normy = 0;

    for(int i = 1; i < (N+1); ++i) {
        for(int j = 2; j < (N+1); ++j) {
            Gi = i * (N+3) + j;
            L1Normx += fabs(h * (F[Gi - i] - F[Gi - i - 1]) +
             h * (P[Gi - i] - P[Gi - i - 1]) +
              h * (HX[Gi + (N2+3)] - HX[Gi]));
        }
    }

    for(int i = 2; i < (N+1); i++) {
        for(int j = 1; j < (N+1); j++) {
            Gi = i * (N+2) + j;
            L1Normy += fabs(h * (G[Gi] - G[Gi - (N2+2)]) +
             h * (P[Gi] - P[Gi - (N2+2)]) +
              h * (HY[Gi + i + 1] - HY[Gi + i]));
        }
    }
    return L1Normx + L1Normy;
}

// RB_GS Error Calculation
double GSError(double* P, double* U, double* V, double dt, double N){
    // Calculates the L1Error norm via the error equation for RB-GS

    int Gi;
    int N2 = (int) N;
    double h = 1 / N;
    double laplaceP;
    double divV;
    double errorNorm = 0.0;

    for(int i = 1; i < (N+1); i++){
        for(int j = 1; j < (N+1); j++){
            Gi = (int) i * (N+2) + j;

            laplaceP = (P[Gi + 1] + P[Gi - 1] + P[Gi + (N2+2)] + P[Gi - (N2+2)] - 4 * P[Gi]) / h / h;

            divV = (U[Gi + i + 1] - U[Gi + i] + V[Gi + (N2+2)] - V[Gi]) / h / dt;

            errorNorm += fabs(divV - laplaceP);
        }
    }
    return errorNorm;
}

// Add two matrices
double* AddMatrix(double* K, double* L, int N) {
    // Adds matrix L into matrix K assuming they're both square of size N

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) {
            K[i * N + j] += L[i * N + j];
        }
    }

    return K;
}