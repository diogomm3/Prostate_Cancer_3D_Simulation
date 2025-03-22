// Vectors are used to store data so that it can be utilized or analised further on the code
// In order to create a vector we need the following command: 'std::vector<type> VectorName;'
// The types of vectors are the same as variables, integer, floating numbers, characters
// We can also assign values to the vector when we create it: 'std::vector<type> VectorName = {10, 20, 30, 40};'
// When the values we want to assign the vector are unknown but we know the vector length we can write 'std::vector<type> VectorName(n);'
// For the 'VectorName(n)' command, it creates a vector with n cells with the number 0 in every singe one of them
// To use vectors we will now need a new library called 'vector'
// To add a new element to the end or the 'back' of the vector we use the '.push_back' function
// To remove an element from the end of a vector we use the '.pop_back' function

#include <iostream>
#include <vector>

int main() {
   
    std::vector<double> subway_adult = {800, 1200, 1500};
  
    std::vector<double> subway_child = {400, 600, 750};

    std::vector<double> subway_baby(3);
  
    std::cout << subway_baby[1] << "\n";
    std::cout << subway_adult[0] << "\n";
    std::cout << subway_child[1] << "\n";
    return 0;
}