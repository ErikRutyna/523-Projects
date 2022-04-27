#include "header.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Primary driver file - runs the CFD program for the given .gri file
// and initial conditions (Minf, alpha, atmosphere properties)
int main(){

    // TODO POST FINISHED: Replace input process with JSON-based config file
    char mesh[50];
    char* meshname;
    double AoA;
    double Minf;
    double altitude;
    double tolerance = pow(10, -6);

    // Get the mesh file
    meshname = (getString(&mesh[0], "\nPlease enter the mesh (.gri) filename: "));

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
    SupersonicCFD(meshname, Minf, AoA, altitude, tolerance);

    // Simulation end time
    time_t endTime = time(NULL);

    printf("\nTime til convergence: %ld\n", startTime - endTime);

    return 0;
}



