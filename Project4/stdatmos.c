#include "header.h"
#include <math.h>

//TODO: Convert altitude to geopotential altitude

double localPressure(double altitude){
    // Returns the local ambient pressure in Pascals (Pa) from given altitude (in meters)

    double altitude_km = altitude / 1000;
    double P; // Pressure in Pa

    if (altitude_km >= 0 && altitude_km <= 11) {

    }
    else if (altitude_km > 11 && altitude_km <= 20) {

    }
    else if (altitude_km > 20 && altitude_km <= 32) {

    }
    else if (altitude_km > 32 && altitude_km <= 47) {

    }
    else if (altitude_km > 47 && altitude_km <= 51) {

    }
    else if (altitude_km > 51 && altitude_km <= 71) {
        
    }
    else if (altitude_km > 71) {

    }

    return P;
}

double localTemperature(double altitude){
    // Returns the local ambient temperature in Kelvin (K) from given altitude (in meters)

    double altitude_km = altitude / 1000;
    double T; // Pressure in Pa

    if (altitude_km >= 0 && altitude_km <= 11) {

    }
    else if (altitude_km > 11 && altitude_km <= 20) {

    }
    else if (altitude_km > 20 && altitude_km <= 32) {

    }
    else if (altitude_km > 32 && altitude_km <= 47) {

    }
    else if (altitude_km > 47 && altitude_km <= 51) {

    }
    else if (altitude_km > 51 && altitude_km <= 71) {
        
    }
    else if (altitude_km > 71) {
        
    }

    return T;
}
