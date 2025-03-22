// Arrays are used to store different values for different variables 

#include <stdlib.h>
#include <stdio.h>

int main() {

    // 1 Dimension Array (Type Intiger)
    int num1[10] = {};
    num1[0] = 5;
    num1[3] = 2;
    num1 [10] = 55;

    printf("Here's your 2 element fo the array: %d\n", num1[1]);

    // 2 Dimension Array (Type Double)
    double table[2][3] = {              // means rows numbers 0 and 1 and colums numbers 0,1,2
        {2,3,4},                // First row
        {50,12,0}               // Second row
    };                    
    table[0][2] = 4.56;
    table[1][0] = 0.3333;

    printf("Some Elements:\n%.1lf\n%.2lf\n%lf",table[0][2],table[1][0],table[0][0]);
}