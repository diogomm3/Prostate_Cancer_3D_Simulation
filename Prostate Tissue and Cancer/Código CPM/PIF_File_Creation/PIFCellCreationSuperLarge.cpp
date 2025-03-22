// Program to create the PIF file with the dimensions of the various cells

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main() {

    std::ofstream MyFile("PIF_File_Creation\\Prostate3DCellsSuperLarge.txt");                // Create and open the described text file
    std::ofstream NewFile("Prostate3DCellsSuperLarge.txt");
    
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
    MyFile << cell << " 1 0 140 0 2 0 120\n";
    NewFile << cell << " 1 0 140 0 2 0 120\n";
    cell++;
    
    // Define the Luminal Cells
    for (int i = 0; i < 120; i += 5) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 2 140 160 " << j << " " << (j+0.05) << " " << i << " " << (i+5) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 2 140 160 " << j << " " << (j+0.05) << " " << i << " " << (i+5) << "\n";
            cell++;
        }
    }

    // Define the Basal Cells
    for (int i = 0; i < 120; i += 5) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 3 160 168 " << j << " " << (j+0.05) << " " << i << " " << (i+5) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 3 160 168 " << j << " " << (j+0.05) << " " << i << " " << (i+5) << "\n";
            cell++;
        }
    }

    // Define the Stroma
    MyFile << cell << " 4 168 285 0 2 0 120";
    NewFile << cell << " 4 168 285 0 2 0 120";
    cell++; 

    MyFile.close();                                             // Closes the file open in the beginning of the program
    NewFile.close();

    return 0;
}