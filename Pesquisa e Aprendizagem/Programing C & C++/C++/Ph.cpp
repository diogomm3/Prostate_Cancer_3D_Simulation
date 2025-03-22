// This code is a ph classifier according to the value in the chemistry scale

#include <iostream>

int main() {
  
    double ph;

    std::cout << "What is the ph value? ";
    std::cin >> ph;
  
    // Write the if, else if, else here:

    if (ph > 7) {
        std::cout << "Basic";
    } else if (ph < 7) {
        std::cout << "Acidic";
    } else {
        std::cout << "Neutral";
    }

    return 0;
}