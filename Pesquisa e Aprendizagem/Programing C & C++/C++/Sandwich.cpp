// Code to create a complex sandwich

#include <iostream>

std::string sandwich_func(int protein) {

    std::string sandwich = "";
    sandwich += "Bread\n";

    if (protein == 0) {
        sandwich += "Eggs\n";
    } 
    else {
        sandwich += "Meat\n";
    }

    sandwich += "Salad\n";

    return sandwich;
}

int main() {
    int protein;
    std::cout << "Eggs - 0 ; Meat - 1\n";
    std::cin >> protein;
    std::string my_sandwich = sandwich_func(protein);
    std::cout << "Here's your sandwich:\n" << my_sandwich;
}