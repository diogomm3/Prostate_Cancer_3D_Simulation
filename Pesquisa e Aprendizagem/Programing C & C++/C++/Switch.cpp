// The 'switch' command is used in C/C++ to write multiple outputs
// This command can only be used with integer numbers
// The argumment in the switch function is going to be analised for every case described in the function
// The break command makes the code jump out of the switch function after the case is observed 
// If none of the cases is observed there should be a default statment (like an else statment in an if cicle)

#include <iostream>

int main() {
  
  int number;

  std::cout << "Species Number? ";
  std::cin >> number;
  
  switch(number) {
    
    case 1 :
      std::cout << "Bulbusaur\n";
      break;
    case 2 :
      std::cout << "Ivysaur\n";
      break;
    case 3 :
      std::cout << "Venusaur\n";
      break;
    case 4 :
      std::cout << "Charmander\n";
      break;
    case 5 :
      std::cout << "Charmeleon\n";
      break;
    case 6 :
      std::cout << "Charizard\n";
      break;
    case 7 :
      std::cout << "Squirtle\n";
      break;
    case 8 :
      std::cout << "Wartortle\n";
      break;
    case 9 :
      std::cout << "Blastoise\n";
      break;
    default :
      std::cout << "Unknown\n";
      break;
  }

    return 0;
}