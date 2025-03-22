// The void function is used when we want to have a function but no input or output
// This functions can be used inside other functions when need 
// It is important to follow the DRY (Don't Repeat Yourself) rule

// Remembering the basic format of a function

/*  
return_type function_name(parameters) {

    Code Block

    return Output
}
*/

// Defining a void function and calling it in a main function after

#include <iostream>

void oscar_wilde_quote() {

    std::cout << "The highest, as the lowest, form of criticism is a mode of autobiography.";

}

 
int main() {
    oscar_wilde_quote();
}