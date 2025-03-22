// Program to create the PIF file with the dimensions of the various cells

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main() {

    std::ofstream MyFile("PIF_File_Creation\\Prostate3DCellsFuture.txt");                // Create and open the described text file
    std::ofstream NewFile("Prostate3DCellsFuture.txt");
    
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
    MyFile << cell << " 1 0 91 0 2 0 40\n";
    NewFile << cell << " 1 0 91 0 2 0 40\n";
    cell++;
    
    // Define the Luminal Cells
    for (float i = 0; i < 40; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Basal Cells
    for (float i = 0; i < 40; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Stroma
    MyFile << cell << " 4 109 200 0 2 0 40\n";
    NewFile << cell << " 4 109 200 0 2 0 40\n";
    cell++; 

    // Define the Lumen
    MyFile << cell << " 1 0 91 0 2 0 40\n";
    NewFile << cell << " 1 0 91 0 2 0 40\n";
    cell++;
    
    // Define the Luminal Cells
    for (float i = 0; i < 40; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 2 91 104 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Basal Cells
    for (float i = 0; i < 40; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 3 104 109 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Stroma
    MyFile << cell << " 4 109 200 0 2 0 40\n";
    NewFile << cell << " 4 109 200 0 2 0 40\n";
    cell++; 

    // Cartesian Coordiantes (cell_type xi xf yi yf zi zf)
    // Define Connection Lumen
    MyFile << cell << " 6 0 10 0 2 0 100\n";
    NewFile << cell << " 6 0 10 0 2 0 100\n";
    cell++;

    // Define the Luminal Cells
    for (float i = 0; i < 100; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 2 10 23 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 2 10 23 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    // Define the Basal Cells
    for (float i = 0; i < 100; i += 4) {
        for (float j = 0; j < 1.99; j += 0.05) {
            MyFile << std::fixed << std::setprecision(2) << cell << " 3 23 28 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            NewFile << std::fixed << std::setprecision(2) << cell << " 3 23 28 " << j << " " << (j+0.05) << " " << i << " " << (i+4) << "\n";
            cell++;
        }
    }

    MyFile.close();                                             // Closes the file open in the beginning of the program
    NewFile.close();

    return 0;
}