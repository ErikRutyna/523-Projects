#include "header.h"
#include <string.h>
#include <stdio.h>

char* getString(char* prompt){
    // Prompts user for an input string
    char response[99];

    printf(prompt);
    scanf("%s", &response);

    return &response;
} 

int getInt(char* prompt){
    // Prompts user for an input integer
    int response;

    printf(prompt);
    scanf("%d", &response);

    return response;
}

double getDouble(char* prompt){
    // Prompts user for an input double
    double response;

    printf(prompt);
    scanf("%lf", &response);

    return response;
}