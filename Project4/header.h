#pragma once
#include <string.h>

/*
=====================================
Driving function that run the program
=====================================
*/
int main(); // Driver function to execute the program
void SupersonicCFD(char* meshFilename, double freesteamMach, double angleOfAttack, double altitude, double tolerance); // Runs compressible Euler on the given mesh and input parameters

/* 
=============================================
Functions to prompt user for input parameters
=============================================
*/

char* getString(char* external, char* prompt); // Prompts user for an input string
int getInt(char* prompt); // Prompts user for an input integer
double getDouble(char* prompt); // Prompts user for an input double

/*
====================================
Atmospheric & Gas Dynamics Functions
====================================
*/

double localPressure(double altitude); // Returns the local ambient pressure in Pascals (Pa) from given altitude (in meters)
double localTemperature(double altitude); // Returns the local ambient temperature in Kelvin (K) from given altitude (in meters)




/*
===============
Data Structures
===============
*/
typedef struct ambientConditions {
    double altitude;
    double pressure;
    double temperature;
    double density;
} ambCond;