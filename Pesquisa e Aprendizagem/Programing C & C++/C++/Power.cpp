// Program that returns a number raised to the power of ten

#include <iostream>
#include <cmath>

int tenth_power(int num) {
    return pow(num, 10);
}

int main() {
  
    std::cout << tenth_power(0) << "\n";
    std::cout << tenth_power(1) << "\n";
    std::cout << tenth_power(2) << "\n";
  
}