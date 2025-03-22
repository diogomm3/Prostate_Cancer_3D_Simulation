// This program determines if an year given by the user is a leap year or not

#include <iostream>

int main() {
  
    int year;

    std::cout << "What Year? ";
    std::cin >> year;

    if (year>=1000 && year<=9999) {
        if (year%4==0 && year%100!=0 || year%400==0) {
            std::cout << "Falls on a leap year.\n";
        } else {
            std::cout << "Is not a leap year.\n";
        }
    
    } else {
        std::cout << "Four Digit Number Is Required.\n";
    }
  
    return 0;
}