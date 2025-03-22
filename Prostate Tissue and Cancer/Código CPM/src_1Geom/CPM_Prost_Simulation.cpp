#include "CPM3D.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <chrono>


#define PI 3.141592653589793238462643383279

int main(int argc, char* argv[]) {
	// The argc variable will have a value of 6
	// argv[1] - directory name
	// argv[2] - 7826345
	// argv[3] - 31
	// argv[4] - 2.55
	// argv[5] - 2
	// argv[6] - 1

	CPM cpm;									// Creates the cpm object from the CPM class


	// Get the input arguments when running the code
	std::string dir(argv[1]);					// First argument - directory name
	int seed = atoi(argv[2]);   				// Second argument - random number seed and converts it to int using the atoi function
	int ifile = atoi(argv[3]);      			// Third argument - created files name
	double J11 = atof(argv[4]); 				// Fourth argument - adhesion energy
	int tumor_layer = atoi(argv[5]); 			// Fifth argument - tummor beggining layer
	int id = atoi(argv[6]);     				// Sixth argument - run index (allows to make different runs to compare results)


	// Domain dimensions
	int Nx, Ny, Nz;								// Creates the integer variables to store the domain dimensions
	int* n_cell;
  	std::ifstream fs;							// Declares fs as the standart library to read inputs form a file 

	// Small 3D cells geometry
	// Nx = 400;  Ny = 400;  Nz = 15;  			// Small: 400x400x15 domain dimension
	// double vtarget1 = 922755;				// Lumen - 9922755 pixels
	// double vtarget2 = 2356;  				// Luminal cells - 2356 pixels
	// double vtarget3 = 1032;					// Basal cells - 1032 pixels
	// double vtarget4 = 1070565;				// Stroma - 1070565 pixels
	// double vtarget5 = 2356;					// Tumoral cells - 2356 pixels
	// fs.open("Prostate3DCellsSmall.txt"); 	// Open the file with the cell geometry (small)
	// std::string  geom = "Small";

	// Large 3D cells geometry
	// Nx = 400;  Ny = 400;  Nz = 45;  			// Large: 400x400x45 domain dimension
	// double vtarget1 = 2768265;				// Lumen - 9922755 pixels
	// double vtarget2 = 2356;  				// Luminal cells - 2356 pixels
	// double vtarget3 = 1032;					// Basal cells - 1032 pixels
	// double vtarget4 = 3211695;				// Stroma - 1070565 pixels
	// double vtarget5 = 2356;					// Tumoral cells - 2356 pixels
	// fs.open("Prostate3DCellsLarge.txt");	// Open the file with the cell geometry (large)
	// std::string geom = "Large";


	// SuperLarge 3D cells geometry	
	// Nx = 400; Ny = 400; Nz = 120;				// SuperLarge : 400x400x120 domain dimension
	// double vtarget1 = 7382040;  				// Lumen - 922755 pixels
	// double vtarget2 = 2356;  				// Luminal cells - 2356 pixels
	// double vtarget3 = 1032;					// Basal cells - 1032 pixels
	// double vtarget4 = 8564520;					// Stroma - 1070565 pixels
	// double vtarget5 = 2356;					// Tumoral cells - 2356 pixels
	// fs.open("Prostate3DCellsSuperLarge.txt");	// Open the file with the cell geometry (superlarge)
	// std::string geom = "SuperLarge";

	// Real 3D cells geometry
	// Nx = 400; Ny = 400; Nz = 60;				// Real: 400x400x60 domain dimension
	// double vtarget1 = 7382040;  				// Lumen - 922755 pixels
	// double vtarget2 = 2356;  				// Luminal cells - 2356 pixels
	// double vtarget3 = 1032;					// Basal cells - 1032 pixels
	// double vtarget4 = 8564520;					// Stroma - 1070565 pixels
	// double vtarget5 = 2356;					// Tumoral cells - 2356 pixels
	// fs.open("Prostate3DCellsReal.txt");		// Open the file with the cell geometry (real)
	// std::string geom = "Real";

	// Second type 3D cells geometry
	Nx = 250; Ny = 250; Nz = 40;
	double vtarget1 = 1039400;
	double vtarget2 = 795; 
	double vtarget3 = 340; 
	double vtarget4 = 1008600; 
	double vtarget5 = 795;					
	fs.open("Prostate3DCells2.txt");
	std::string geom = "SecondGeom";
	
	//Adhesion matrix creation for all the cells (including the extracelular matrix)
	// First parameter - adhesion to extracellular matrix
	// Second parameter - adhesion to cells of the same type (lower value meaning higher adhesion)
	// All the ohter parameters - adhesion to other cells (higher value meaning lower adhesion)
	double Jxx = J11*0.90;													// Adhesion for cells of the same type
	double J01 [6] = {1.8000,    Jxx,    J11,    J11,    J11,    J11}; 		// Lumen adhesion matrix
	double J02 [6] = {J01[0], J01[2],    Jxx,    J11,    J11,    J11}; 		// Luminal cells adhesion matrix
	double J03 [6] = {J01[0], J01[3], J02[3],    Jxx,    J11,    J11};		// Basel cells adhseion matrix
	double J04 [6] = {J01[0], J01[4], J02[4], J03[4],    Jxx,    J11}; 		// Stroma adhesion matrix
	double J05 [6] = {J01[0], J01[5], J02[5], J03[5], J04[5], 	 Jxx}; 		// Tumoral cells adhesion matrix
	

	// Hamiltonian Parameters - Cell characteristics and Parameters Weight
	// Target volume deviation parameters
	double lambdav1 = 50;  		// Lumen
	double lambdav2 = 7500;		// Luminal cells
	double lambdav3 = 10000;	// Basal cells
	double lambdav4 = 5e7;		// Stroma
	double lambdav5 = 20000;	// Tumoral cells
	
	// Target length deviation parameter (set to 0 because it is not being used)
	double lambdal = 0;  	// All the cells	
	
	// Target length for the different cells (the value does not have any implication because the deviation parameter is set to 0)
	double ltarget1 = 1;  	// Lumen
	double ltarget2 = 1;  	// Luminal cells
	double ltarget3 = 1;	// Basal cells
	double ltarget4 = 1;	// Stroma
	double ltarget5 = 1;	// Tumoral cells
	
	// Parameters that are not currently being used
	double lambdap1 = 10;  	// Amplitude of polarization condition
	double Chi = 20;  		// Chemotaxis parameter


	// Pariodic conditions
	bool x_periodic = false;  				// x axis periodic condition
	bool y_periodic = false;				// y axis periodic condition
	bool z_periodic = true;					// z axis periodic condition 

	// System dynamics parameters
	int n_MCS = 1000;   					// Number of monte carlo steps
	double L = 0.0000025;  					// Pixel lateral size - 
	double acc = 1E-5;						// Not being used
	int class_range = 1;					// Not being used
	double class_tol = 0.35;				// Not being used
	int mesh_back_diff = 1;					// Not being used
	int mesh_range = 0;						// Not being used
	double vol_mitosis = 51;  				// Not being used
	bool flag_change = true;				// Not being used
	bool flag_first = true;					// Not being used


	char buffer [6], num [10];

	sprintf(buffer, "%d", id);
	std::string id_s(buffer);

	std::ofstream NewFile("Simulation_Data\\"+geom+"_Cell_Number_EVO_"+id_s+".txt");
	
	// Creates a file called with the specifications of the code running
	std::ofstream log_file;  																			
	log_file.open((dir+"/Settings_File_"+id_s+".log").c_str());
	log_file << "*** Simulation Configuration ***\n\n - GENERAL - "+geom+"\n";
	log_file << "Seed: " << seed << "\n";
	log_file << "Number of Monte Carlo Steps: " << n_MCS << "\n";
	log_file << "Grid size: " << Nx << "X" << Ny << "X" << Nz << "\n";
	log_file << "Grid element side length: " << L << "\n";
	log_file << "Periodic Condition: x - " << x_periodic << "\n";
	log_file << "Periodic Condition: y - " << y_periodic << "\n";
	log_file << "Periodic Condition: z - " << z_periodic << "\n";
	log_file << "\n - CELLULAR POTTS MODEL\n";
	log_file << "Luminal cells target volume: " << vtarget2 << "\n";
	log_file << "Luminal cells inelasticity constant (lambdav2): " << lambdav2 << "\n";
	log_file << "Basal cells target voluem: " << vtarget3 << "\n";
	log_file << "Basal cells inelasticity constant (lambdav3): " << lambdav3 << "\n";
	log_file << "Lumen inelasticity constant (lambdav1): " << lambdav1 << "\n";
	log_file << "Stroma inelasticity constante (lambdav4): " << lambdav4 << "\n";
	log_file << "Tumoral cells target volume: " << vtarget5 << "\n";
	log_file << "Tumoral cells inelasticity constant (lambdav5): " << lambdav5 << "\n";
	log_file << "Tumoral cells initial layer: " << tumor_layer << "\n";
	log_file << "Cell1-Cell1: " << J01[1] << "\n";
	log_file << "Cell1-Cell2 adhesion: " << J01[2] << "\n";
	log_file.close();


	sprintf(buffer, "%d", ifile);
	std::string id_s2(buffer);
	std::ofstream file_xCells (("Simulation_data/"+geom+"_Cells_TumorDynamics"+id_s+".dat").c_str());

    
	cpm.setup(Nx, Ny, Nz, seed); // setup of CPM
	cpm.set_periodicity(x_periodic, y_periodic, z_periodic);

	cpm.init();   // add all cell types with their respective parameters


	// Cell creation with the various parameters
	cpm.add_type(lambdav1, vtarget1, lambdal, ltarget1, 0, &J01[0], 0);		// Lumen
	cpm.add_type(lambdav2, vtarget2, lambdal, ltarget2, 0, &J02[0], 0);  	// Luminal cells
	cpm.add_type(lambdav3, vtarget3, lambdal, ltarget3, 0, &J03[0], 0); 	// Basal cells
	cpm.add_type(lambdav4, vtarget4, lambdal, ltarget4, 0, &J04[0], 0); 	// Stroma
	cpm.add_type(lambdav5, vtarget5, lambdal, ltarget5, 0, &J05[0], 0); 	// Tumoral cells
	

	srand(seed);

  	
  	int xcell, xtype;		// Create the integer variables to store the cell characteristics (cylindrical coordinates)
	float ti, tf, pi, pf, zi, zf;							// Create the float variables to store the cell angle variations 
  	int index;								// Create an integer variable to store the cell index
  	
  	std::string Linha;						// Declares nhe as the standart string library

  	while (getline(fs, Linha)) {			// Reads the file with the cell geometry line by line and stores each line in the Linha Variable
		std::stringstream sLinha;			// Creates the sLinha object from the std::stringstream class 
		sLinha << Linha;					// Extracts the values from Linha to sLinha
	
		// Extracts the values from each line separated by spaces
		sLinha >> xcell;					// Cell number			
		sLinha >> xtype; 					// Cell type
		sLinha >> pi; 						// Cell rho start
		sLinha >> pf; 						// Cell rho end
		sLinha >> ti;	 					// Cell theta start
		sLinha >> tf; 						// Cell theta end
		sLinha >> zi; 						// Cell z start
		sLinha >> zf; 						// Cell z end

		// Set tumoral cell based on the starting membrane
		// Small Cells Dimensions
		if (geom.compare("Small") == 0) {
			if ( (tumor_layer==2) && (xcell==60) ) xtype = 5;  // Luminal Membrane
			if ( (tumor_layer==3) && (xcell==180) ) xtype = 5; // Basal Membrane
		}

		//Large Cells Dimensions
		if (geom.compare("Large") == 0) {
			if ( (tumor_layer==2) && (xcell==180)) xtype = 5;	// Luminal Membrane
			if ( (tumor_layer==3) && (xcell==540)) xtype = 5;	// Basal Membrane
		}

		//SuperLarge Cells Dimensions
		if (geom.compare("SuperLarge") == 0) {
			if ( (tumor_layer==2) && (xcell==500)) xtype = 5;	// Luminal Membrane
			if ( (tumor_layer==3) && (xcell==1460)) xtype = 5;	// Basal Membrane
		}

		//Real Cells Dimensions
		// if (geom.compare("Real")) {
			// if ( (tumor_layer==2) && (xcell==)) xtype = 5;		// Luminal Membrane
			// if ( (tumor_layer==3) && (xcell==)) xtype = 5;		// Basal Membrane
		//}

		//Second Geometry Cells Dimensions 
		if (geom.compare("SecondGeom") == 0) {
			if ( (tumor_layer==2) && (xcell==220)) xtype = 5; 
			if ( (tumor_layer==3) && (xcell==620)) xtype = 5;
		}


		if (xtype==5) std::cout << "Tumor_layer (" << tumor_layer << "): CellNum = " << xcell << "; Coordinates: " << pi << "  " << pf << "  " << ti << "  " << tf << "  " << zi << "  " << zf << std::endl; 
		cpm.add_cell_vol(pi, pf, ti, tf, zi, zf, xtype, xcell, Ny);
		
  	}
  	fs.close();								// Closes the opened file with the cells geometry

	
	int tumor_volume;						// Variável para guardar o volume do tumor

	
	void** data = new void* [5];
	data[0] = &cpm;
	
	cpm.set_dataref(&data[0]);

	int i, nchanged = 0;
	//std::cout << " ndead1 " << ndead1 << std::endl;
	
	sprintf(buffer, "%d", ifile);
	std::string id_si(buffer);
	std::string term("Simulation_Data/"+geom+"_Test_TumorDynamics"+id_s+".dat");
	cpm.print_types(term);  // types


	std::cout << "\nCPM BEGIN :: MCS\n" << std::endl; 
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int cell_deaths;	

	for (i=0; i<=n_MCS; ++i) { 				// Cycle to run all Monte Carlo Steps
		
	auto start1 = std::chrono::high_resolution_clock::now();

	cpm.MCS(); 								// Pixel copy atempts

	auto start2 = std::chrono::high_resolution_clock::now();

	tumor_volume = cpm.proliferation(vol_mitosis);  // try cell proliferation

	auto start3 = std::chrono::high_resolution_clock::now();

	cell_deaths += cpm.death_pressure();

	auto start4 = std::chrono::high_resolution_clock::now();

	if( nchanged>0 ) std::cout << " nchanged = " << nchanged << std::endl;
		
		if (i%1==0) {  // printout
			cpm.print_xCells(i, file_xCells);  // count number of cells
			n_cell = cpm.n_cell();
			std::cout << "MCS = " << i << "; Tumor_Cells_Num = " << n_cell[0] << "; Tumor_Volume = " << tumor_volume << std::endl;
			NewFile << i << " " << n_cell[0] << " " << n_cell[1] << " " << n_cell[2] << std::endl;
		}

		if (i%100==0) {
			sprintf(buffer, "%d", ifile);
			std::string id_si(buffer);
			std::string term("Simulation_Data/"+geom+"_Test_TumorDynamics"+id_s+".dat");
			cpm.print_types(term);  // types
		}

        	std::chrono::duration<double>elapsed1 = start2 - start1;
        	std::chrono::duration<double>elapsed2 = start3 - start2;
        	std::chrono::duration<double>elapsed3 = start4 - start3;
        	
        	std::cout << "Pixel_Copys = " << elapsed1.count() << "; Cell_Proliferation = " << elapsed2.count() << "; Cells_Death = " << elapsed3.count() << ";\n" << std::endl;

	}  // end MAIN i cycle (MCS)


	std::cout << "Number of cells murdered by pressure: " << cell_deaths << std::endl;

	file_xCells.close();

	delete [] data;
	cpm.free_mem();

	NewFile.close();
}