// The iostream contains functions for basic input and output operations
// Using the comand include we are including iostream in our code file (is used to include libraries)
#include<iostream>

// The main function is were the program begins, it is followed by parenthesis and curly brackets    
int main() {
    // The command 'std' means standard and is used to return an output
    // It is followed by '::' and 'cout', wich means 'c' for character and 'out' for output 
    // Followed by the previous instructions we write between '<<' the characters we want to be returned
    // After a statment, we write a ';' to let the compiler know that the statement is compelted
    std::cout << "Hello World!\n";
    std::cout << "My name is Diogo!";
    // If the command return 0 is reached, that means that the code was completely correct
    // If the number 1 is returned that means that the code has a issue that needs to be fixed
    return 0;
}

// We compile the program so that the human-language written in the 'cpp' file can be read by the CPU
// In order to compile the 'ccp' file we will need to write in the terminal the following command: g++ FileName.cpp -o WantedName
// This will create a new file with the 'exe' extension and the name we choose: WantedName
// The extension file created when we compile de 'cpp' file is an 'exe' (Executable File)
// In order to execute the 'exe' file and receive the output from the original 'cpp' file we use de command: ./WantedName