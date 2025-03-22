// How to iterate a number 

#include <stdio.h>                                  // Opens the stdio library (ins and outs)
#include <stdlib.h>                                 // Opens the stdlib library (library with standart functions)

int main() {

    int x = 10;                                     // Creates the variable x and assigns it's value to 10
    int y;                                          // Creates the variable y without assigning a value to it

    // Iteration cycle 
    for (int i = 0; i < 10; i++) {                  // Creates a for cycle that runs from i = 0 to i = x-1 and that after the cycle ir runned iterates i
        x++;                                        // Iterates the x value
        y = x--;                                    // Assigns the value x to y and iterates negativly x
        printf("Your number x is: %d\n", x);        // Prints x to the terminal
        printf("Your number y is: %d\n", y);        // Prints y to the terminal
    }

    printf("Final value for x is: %d\n", x);        // Prints x to the terminal
    printf("Final value for y is: %d\n", y);        // Prints y to the terminal

}

// The iteration cycle will be runned 10 times
// For each cycle the value of x is increased by one unit and decreased by one unit, meaning that the x final value will be the initial, x = 10
// The y value will be the final value of x increased by one unit, because the last statment assigns x to y and decreases x by one, y = 11