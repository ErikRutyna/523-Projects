#include <stdio.h>
#include <time.h>
#include <math.h>
#include "headers.h"

int main(){

	// Givens for the set problem like grid size and Reynolds Number
	double N;
	double Re;
	int uW = 1;
	double tol = pow(10, -5);
	double beta = 0.1;

	// Getting simulation to run and printing inputs to check
	printf("\nPlease enter the grid size (powers of 2 only): ");
	scanf("%lf", &N);

	printf("\nPlease enter the Reynolds Number: ");
	scanf("%lf", &Re);

	printf("\nSimulation w/ grid size N = %i for Re = %i.\n", (int) N, (int) Re);
	
	// Want to test times - so need total real-world clock time
	time_t startTime = time(NULL);

	CavityFlow(N, Re, uW, tol, beta);

	time_t endTime = time(NULL);

	// Difference in time -> total run-time
	printf("The simulation converged in %ld seconds.\n", (endTime -  startTime));

	return 0;
}
