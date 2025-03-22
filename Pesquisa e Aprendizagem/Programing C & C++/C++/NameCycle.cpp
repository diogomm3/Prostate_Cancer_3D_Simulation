// Create a program that repeats a given name x times

#include <iostream>

void name_cycle(std::string name , int x) {

    while (x > 0) {
        std::cout << name << "\n";
        x--;
    }

}

int main() {

    std::string name;
    int number;

    std::cout << "What's your name? ";
    std::cin >> name;
    std::cout << "How many times? ";
    std::cin >> number;

    std::cout << "\nHere's your cycle: \n";
    name_cycle(name , number);

}