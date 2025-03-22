// There are several data types that we can use when writing code

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

// To store values or characters or other types of data on a program we need to use variables
// To create a variable, first we need to especify the varaiable type

int main() {
    // Intiger
    int num1 = 3;
    printf("num1 = %d\n",num1);

    // Double
    double num2 = 3.57;
    printf("num2 = %lf\n",num2);

    // Float
    float num3 = 1/3;
    printf("num3 = %f\n",num3);

    // Character
    char c = 'A';
    printf("c = %c\n",c);

    // Boolian  (to use this, we need the library 'stdbool.h')
    _Bool b = false;
    printf("b = %d\n",b);

    // String
    char name[] = "Diogo MM";
    printf("name = %s\n",name);

    // How much can a variable hold?
    int x = sizeof(int);                                            // Get the number of bytes in an intiger an store it inside x
    printf("The number of bytes of an intiger is: %d\n",x);           // Prints the number of bytes of an intiger
    printf("The number of bits in a byte is: 4x8=32 bits\n");         // Prints the number of bits in a byte
    printf("How much we can store on a variable type int: 2^32\n");   // 4294967296
}

// To see the address of a certain variable in the memory use the symbol '&' before the variable name 

// To limite the number of decimals in a variable, before specifying the variable type in 'printf' write '.x'
// This will limit the number of decimals of the variable to x decimals ('%.4lf' - get 4 decimals)