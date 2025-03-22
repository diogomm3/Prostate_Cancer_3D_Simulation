// To get the user input for a certain variable we use the 'scanf' function

#include <stdlib.h>
#include <stdio.h>

int main() {

    int x;
    char c[] = "";
    printf("What's your number? ");
    scanf("%d", &x);
    fflush(stdin);          // This function gets rid of everything that was not used for the x input
    printf("Whats your name? ");
    scanf("%s", c);

    printf("Your number is: %d\n", x);
    printf("Your name is: %s\n", c);
}