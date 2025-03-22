// Code to change the Earth weight to the weight from another planet on the Solar System

#include <iostream>

int main() {
  
    double weight;
    int number;

    std::cout << "Earth Weight? ";
    std::cin >> weight;

    std::cout << "Planet Fight Number? ";
    std::cin >> number;
  
    switch(number){

        case 1:
            std::cout << "Mercury Weight is: " << weight*.38 << "kg.\n";
            break;
        case 2:
            std::cout << "Venus Weight is: " << weight*.91 << "kg.\n";
            break;
        case 3:
            std::cout << "Mars Weight is: " << weight*.38 << "kg.\n";
            break;
        case 4:
            std::cout << "Jupiter Weight is: " << weight*2.38 << "kg.\n";
            break;
        case 5:
            std::cout << "Saturn Weight is: " << weight*1.06 << "kg.\n";
            break;
        case 6:
            std::cout << "Uranus Weight is: " << weight*.92 << "kg.\n";
            break;
        case 7:
            std::cout << "Neptune Weight is: " << weight*1.19 << "kg.\n";
            break;
        default:
            std::cout << "Unknown Planet.\n";
            break;
  }
  
  return 0;
}