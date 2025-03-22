// Basic while loop to return the square of the numbers from 0 to 9

#include <iostream>

int main() {
  
    int i = 0;
    int square = 0;
  
    while (i < 10) {
        square = i*i;
        std::cout << i << "  " << square << "\n";
        i++;
    }
  
    return 0;
}