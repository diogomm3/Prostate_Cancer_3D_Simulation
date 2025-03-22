// Determine the volume of each cell to bu used in the model

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

int main() {
    std::ifstream MyFile("Cell_Volume_Analysis/Prostate_Cell_Characteristics.txt");

    std::string Line;

    int type;
    std::string vol;
    int pixel;
    std::string size;
    int siz;

    std::vector<int> cell2;
    std::vector<int> cell3;

    while (getline(MyFile, Line)) {
        std::stringstream sLine;
        sLine << Line;

        sLine >> type;
        sLine >> vol;
        sLine >> pixel;
        sLine >> size;
        sLine >> siz;

        if (pixel < 3000 && pixel > 2000) {
            cell2.push_back(pixel);
        }
        else if (pixel < 2000 && pixel > 1000) {
            cell3.push_back(pixel);
        }
        //std::cout << pixel << "\n";
    }
    

    int pixel2 = 0;
    for (int i = 0; i < cell2.size(); i++) {
        pixel2 += cell2[i];
    }
    std::cout << pixel2/cell2.size() << "\n";

    int pixel3 = 0;
    for (int i = 0; i < cell3.size(); i++) {
        pixel3 += cell3[i];
    }
    std::cout << pixel3/cell3.size() << "\n";

    return 0;
}