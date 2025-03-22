// A logical operator is used when teh logic can not be satisfied with only one condition
// The operator '&&' is used as an and operator
// The operator '||' is used as an or operator
// The operator '!' is used as a not operator
// These operators return a bool value as a result, meaning it's either true or false

#include <iostream>

int main() {
    
    // Operator '&&'
    int hunger = true;
    int anger = true;
    // Write the code below:
    if (hunger && anger) {
        std::cout << "Hangry\n";
    }
    
    // Operator '||'
    int day = 6;
    // Write the code below:
    if (day==6 || day==7) {
        std::cout << "Weekend\n";
    }

    // Operator '!'
    bool logged_in = false;
    // Write the code below:
    if (!logged_in) {
        std::cout << "Try again.\n";
    }
    
    return 0;
}