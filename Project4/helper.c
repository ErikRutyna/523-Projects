#include "header.h"
#include <string.h>
#include <stdio.h>

char* getString(char* external, char* prompt){
    // Prompts user for an input string

    printf("%s", prompt);
    scanf("%s", external);

    return external;
} 

int getInt(char* prompt){
    // Prompts user for an input integer
    int response;

    printf("%s", prompt);
    scanf("%d", &response);

    return response;
}

double getDouble(char* prompt){
    // Prompts user for an input double
    double response;

    printf("%s", prompt);
    scanf("%lf", &response);

    return response;
}