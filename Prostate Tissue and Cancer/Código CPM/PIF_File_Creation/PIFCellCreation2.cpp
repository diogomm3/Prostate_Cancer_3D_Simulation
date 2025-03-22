// Program to create the PIF file with the dimensions of the various cells

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main() {

    std::ofstream MyFile("PIF_File_Creation\\Prostate3DCells2.txt");                // Create and open the described text file
    std::ofstream NewFile("Prostate3DCells2.txt");
    
    /*
    Lumen - Cell Type 1
    Luminal Cells - Cell Type 2
    Basal Cells - Cell Type 3
    Estroma - Cell Type 4

    Text file structure:
    Cell Number; Cell Type; p_ini; p_fin; theta_ini; theta_fin; z_ini; z_fin 
    */

    int cell = 1;

    // Define the Lumen
    MyFile << cell << " 1 0 91 0 2 0 120\n";
    NewFile << cell << " 1 0 91 0 2 0 120\n";
    cell++;
    
    // Define the Luminal Cells
    for (float i = 0; i < 120; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Basal Cells
    for (float i = 0; i < 120; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Stroma
    MyFile << cell << " 4 109 200 0 2 0 120";
    NewFile << cell << " 4 109 200 0 2 0 120";
    cell++; 

    MyFile.close();                                             // Closes the file open in the beginning of the program
    NewFile.close();

    return 0;
}