// Program to ask the user the temperature in Fahrenheit and return the temperature in Degrees Celsius

#include <iostream>

int main() {

    double tempf;
    double tempc;

    std::cout << "The temperature (Fahrenheit) in NY is: ";
    std::cin >> tempf;

    tempc = (tempf-32)/1.8;

    std::cout << "The temperature in Degrees Celsius is: " << tempc << " C.\n";

    return 0;
}