#include "header.h"
#include <stdio.h>
#include <math.h>

void SupersonicCFD(char* meshFilename, double freesteamMach, double angleOfAttack, double altitude, double tolerance){

    printf("\nInside the CFD function.\n");

    printf("\nMesh filename: %s", meshFilename);
    printf("\nFreestream Mach number: %lf", freesteamMach);
    printf("\nAngle of Attack (degrees): %d", angleOfAttack);
    printf("\nAltitude (m): %lf", altitude);
    printf("\nConvergance Tolerance: %lf", tolerance);
    

    return;
}