// Program that gives access to an account if the given pin is correct
// If it's not correct it will ask for the pin again and again

#include <iostream>

int main() {
  
    int pin = 0;
    int tries = 0;
  
    std::cout << "BANK OF CODECADEMY\n";
  
    std::cout << "Enter your PIN: ";
    std::cin >> pin;

    tries++;

    while (pin != 1234 && tries < 3) {
    
        std::cout << "Enter your PIN: ";
        std::cin >> pin;
        tries++;
    } 
  
    if (pin == 1234) {
    
        std::cout << "PIN accepted!\n";
        std::cout << "You now have access.\n"; 
    }
  
    return 0;
}