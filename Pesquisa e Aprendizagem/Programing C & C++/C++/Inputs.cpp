// Inputs

// A lot of times we want the user to give us the value for a certain variable
// For that we need to ask him for the variable value through the  input 'cin' instead of 'cout'

#include <iostream>

int main() {

    double price;
    
    std::cout <<"Enter the product price: ";
    std::cin >> price;

    std::cout << "You paid " << price << " dollars for the product.";

    return 0;
}

// We use the variable type 'double' for floating-point numbers (with decimals)
// For the 'cin' comand we need to use the characters '>>' (inverted when compared to cout) 