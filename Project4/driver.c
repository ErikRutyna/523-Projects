#include "header.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Primary driver file - runs the CFD program for the given .gri file
// and initial conditions (Minf, alpha, atmosphere properties)
int main(){

    // TODO POST FINISHED: Replace input process with JSON-based config file
    char mesh;
    double tolerance = pow(10, -6);
    double AoA;
    double Minf;
    double altitude;

    // Get the mesh file
    mesh = *(getString("\nPlease enter the mesh (.gri) filename: "));

    // Get mach number
    Minf = getDouble("\nPlease enter the Mach number: ");

    // Get AoA
    AoA = getDouble("\nPlease enter the angle of attack (in degrees): ");

    // Calculate atmospheric properties based on height
    // via formulas: http://www.braeunig.us/space/atmmodel.htm
    altitude = getDouble("\nPlease enter the altitude (in meters): ");

    // Simulation start time
    time_t startTime = time(NULL);

    // Simulation function
    SupersonicCFD(*mesh, Minf, AoA, altitude, tolerance);

    // Simulation end time
    time_t endTime = time(NULL);

    printf("Time til convergence: %ld", startTime - endTime);

    return 0;
}



