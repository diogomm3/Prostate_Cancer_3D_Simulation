// How to use Strings?
// Strings are used between double quotes while characters are used between single quotes

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

int main() {

    char FirstName[] = "Diogo";     // Size of the string = 6 ; "D" "i" "o" "g" "o" "\0"
    char SecondName[] = "Martins";  // Size of the string = 8 ; "M" "a" "r" "t" "i" "n" "s" "\0"

    char Garbage[5] = {};           // Creating a string with 5 characters '\0'

    Garbage[0] = 'f';
    Garbage[1] = 'g';
    Garbage[2] = 'h';
    Garbage[3] = 'A';

    printf("Here's garbage: %s\n", Garbage);
    printf("Here's it's length: %d\n", strlen(Garbage));
    printf("Here's it's bytes: %d\n",sizeof(Garbage));

    char CopyFirstName[5] = {};
    strcpy(CopyFirstName,FirstName);
    printf("FirstNameCopy: %s\n", CopyFirstName);

    char FullName [strlen(FirstName) + strlen(SecondName) + 1];
    strcpy(FullName, FirstName);
    strcat(FullName, SecondName);
    printf("Here's your full name: %s\n", FullName);

    printf("Compare 2 strings: %d\n", strcmp(FirstName,SecondName));
    printf("Compare 2 strings: %d\n", strcmp(FirstName,CopyFirstName));

    // Returns 0 if the strings are equal, -1 if the first string e less than the second and 1 if the first string is greater than the second

    char str1[13] = {};
    gets(str1);
    printf("String: %s\n", str1);
    printf("String size: %d\n", strlen(str1));
    printf("String bytes: %d\n", sizeof(str1));

}