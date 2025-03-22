#include <iostream>
#include <string>

int main() {
    std::string name;

    while (name.empty()) {
        std::cout << "Name: ";
        std::getline(std::cin, name);
    } 
    std::cout << name;
} 