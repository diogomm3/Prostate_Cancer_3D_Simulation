// The 'if' statment is used to test an expression for truth and execute the code in it
// if (condition) {code we want executed if the condition is true}
// If the condition is false, the code within it is skipped and the program continues to run

#include <iostream>

int main() {

    int grade;

    std::cout << "Student Grade: ";
    std::cin >> grade;

    if (grade > 60) {
        std::cout <<"Approved";
    }
    else {
        std::cout <<"Disapproved";
    }

    return 0;
}

// This code takes in the grade from a student introduced by the user 
// If the grade is greater than 60 the student is approved
// If the grade is smaller than 60 the student is disapproved

// If we want to add more than two conditions, simply use the command 'else if {code}' before the 'else' command