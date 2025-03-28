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
//#define GNUPLOT "gnuplot"

#define PI 3.141592653589793238462643383279

int main(int argc, char* argv[]) {

        //FILE *gp;
        //gp = popen(GNUPLOT, "w"); // open gnuplots for final cells' distribution

	CPM cpm;

	//for (int ir=0; ir<10; ir++) std::cout << " seed = " << ((int)(genrand_int31()  % 10000000)) << std::endl;

	std::string dir(argv[1]);
	int seed = atoi(argv[2]);   // random number seed
	int ifile = atoi(argv[3]);      // not being used

	ifile = ifile + 9250;		// Parameter TESTS

	double J11 = atof(argv[4]); // adhesion energy
	int tumor_layer = atoi(argv[5]); // 1: Basal, 2-6: Intermediate, 7: Umbrella
	///tumor_layer = 0; // no tumor!
	int id = atoi(argv[6]);     // run
	//int ifile = atoi(argv[9]);  // file
	int Nx, Ny, Nz, ncell;
	// pifUrot3D_CPM.txt
	//Nx = 400;  Ny = 300;  Nz = 20;  // 400x300x20 domain dimension
	// pifUrot3D_CPM_large.txt
	Nx = 400;  Ny = 300;  Nz = 80;  // 400x300x80 domain dimension
	double Jxx = J11*0.90;
	double Jyy = J11/0.90;

//	adhesion matrix (for all cell types)
	double J01 [10] = {1.8000,    Jxx,    J11,    J11,    J11,    J11,    J11,    J11,    J11, J11}; // 1.25 1.50  stroma
	double J02 [10] = {J01[0], J01[2],    Jxx,    J11,    J11,    J11,    J11,    J11,    J11, J11}; // 1.25 1.50 membrane
	double J03 [10] = {J01[0], J01[3], J02[3],    Jxx,    J11,    J11,    J11,    J11,    J11, J11}; // 1.25 1.50 basal
	double J04 [10] = {J01[0], J01[4], J02[4], J03[4],    Jxx,    J11,    J11,    J11,    J11, J11}; // 1.25 1.50 stem
//	double J02 [10] = {1.00, J11, Jxx, J11, 0.1, J11, J11, J11, J11, J11}; // 1.25 1.50 membrane
//	double J03 [10] = {1.00, J11, J11, Jxx, 0.1, J11, J11, J11, J11, J11}; // 1.25 1.50 basal
//	double J04 [10] = {1.00, J11, 0.1, 0.1, Jxx, J11, J11, J11, J11, J11}; // 1.25 1.50 stem
//	double J05 [10] = {1.00, J11, J11, J11, J11, Jxx, J11, J11, J11, J11}; // 1.25 1.50 intermediate
//	double J06 [10] = {1.00, J11, J11, J11, J11, J11, Jxx, J11, J11, J11}; // 1.25 1.50 umbrella
	double J05 [10] = {J01[0], J01[5], J02[5], J03[5], J04[5],    Jxx,    J11,    J11,    J11, J11}; // 1.25 1.50 intermediate
	double J06 [10] = {J01[0], J01[6], J02[6], J03[6], J04[6], J05[6],    Jxx,    J11,    J11, J11}; // 1.25 1.50 umbrella
	//double J07 [10] = {1.00, J11, J11, J11, J11, J11, J11, Jxx, J11, J11}; // 1.25 1.50 urine
	double J07 [10] = {J01[0], J01[7], J02[7], J03[7], J04[7], J05[7], J06[7],    Jxx,    J11, J11}; // 1.25 1.50 urine
	//double J08 [10] = {1.00, J11, J11, J11, J11, J11, J11, J11, Jxx, J11}; // 1.25 1.50 tumor
	double J08 [10] = {J01[0], J01[8], J02[8], J03[8], J04[8], J05[8], J06[8], J07[8],    Jxx, J11}; // 1.25 1.50 tumor
	double J09 [10] = {J01[0], J01[9], J02[9], J03[9], J04[9], J05[9], J06[9], J07[9], J08[9], Jxx}; // 1.25 1.50 capillary

/*
	Jyy = 2.55*1.1;  // 2.55*2  2.55*1.5
	for (int icell=1; icell<10; icell++) J08 [icell] = Jyy;	// tumor cell with smaller adhesion
	J01[8] = Jyy; J02[8] = Jyy; J03[8] = Jyy; J04[8] = Jyy; J05[8] = Jyy;
	J06[8] = Jyy; J07[8] = Jyy; J08[8] = Jxx; J09[8] = Jyy;
*/
/*
	Jyy = 2.55*1.1;  // 2.55*2  2.55*1.5
	for (int icell=3; icell<=6; icell++) J08 [icell] = Jyy;	// tumor cell with smaller adhesion
	J03[8] = Jyy;  J04[8] = Jyy;  J05[8] = Jyy;  J06[8] = Jyy;
*/
	double lambdav1 = 250;  // 500 250 amplitude of volume condition
	double lambdav2 = 120;  // 500 2000 amplitude of volume condition
	// old pif
	double vtarget01 = 25000, vtarget02 = 500, vtarget03 = 25, vtarget04 = 25,
	       vtarget05 = 100, vtarget06 = 250, vtarget07 = 97000, vtarget08 = 25, vtarget09 = 100;  // target cell volumes
	// new pif
	//double vtarget01 = 25000, vtarget02 = 504, vtarget03 = 144, vtarget04 = 144, vtarget05 = 192, vtarget06 = 216, vtarget07 =89208, vtarget08 = 192, vtarget09 = 100;  // target cell volumes
	double lambdal1 = 50;  // 500  amplitude of length condition
	double lambdal2 = 50;  // 500 amplitude of length condition
	// old pif
	double ltarget01 = 502, ltarget02 = 500, ltarget03 = 7, ltarget04 = 7,
	       ltarget05 = 14, ltarget06 = 50, ltarget07 = 535, ltarget08 = 5, ltarget09 = 10;  // target cell lengths
	// new pif
	//double ltarget01 = 506, ltarget02 = 504, ltarget03 = 17, ltarget04 = 17, ltarget05 = 20, ltarget06 = 22, ltarget07 = 535, ltarget08 = 20, ltarget09 = 10;  // target cell lengths
	double lambdap1 = 10;  // 50  amplitude of polarization condition
	double Chi = 20;  // 20 chemotaxis parameter

	bool v_periodic = false;  // true periodicity of border condition
	bool h_periodic = true;	// x
	bool d_periodic = true;	// z
	int n_MCS = 1002;   // 1002 5000 501 851  number of simulation step
	double L = 0.0000025;  // 2.5 um pixel size
	double acc = 1E-5;


	int class_range = 1;
	double class_tol = 0.35;
	int mesh_back_diff = 1;
	int mesh_range = 0;
	double vol_mitosis = 51;  // 51
	bool flag_change, flag_first = true;

	char buffer [6], num [10];

	sprintf(buffer, "%d", id);
	std::string id_s(buffer);

	std::ofstream log_file;  // information for log file
	log_file.open((dir+"config_"+id_s+".log").c_str());
	log_file << "*** Simulation Configuration ***\n\n - GENERAL\n";
	log_file << "seed: " << seed << "\n";
	log_file << "number of Monte Carlo steps: " << n_MCS << "\n";
	log_file << "grid size: " << Nx << "X" << Ny << "X" << Nz << "\n";
	log_file << "grid element side length: " << L << "\n";
	log_file << "periodic: h, v - " << h_periodic << "," << v_periodic <<"," << d_periodic << "\n";
	log_file << "\n - CELLULAR POTTS MODEL\n";
	log_file << "target volume: " << vtarget01 << "\n";
	log_file << "inelasticity constant (lambdav1): " << lambdav1 << "\n";
	log_file << "inelasticity constant (lambdav2): " << lambdav2 << "\n";
	log_file << "chemotaxis constant (Chi): " << Chi << "\n";
	log_file << "cell1-ECM adhesion: " << J01[0] << "\n";
	log_file << "cell1-cell1 adhesion: " << J01[1] << "\n";
	log_file << "cell1-cell2 adhesion: " << J01[2] << "\n";
	log_file << "cell2-ECM adhesion: " << J02[0] << "\n";
	log_file << "cell2-cell1 adhesion: " << J02[1] << "\n";
	log_file << "cell2-cell2 adhesion: " << J02[2] << "\n";
	log_file.close();


	sprintf(buffer, "%d", ifile);
	std::string id_s2(buffer);
        std::ofstream file_xCells (("tests_ur3D/cells"+id_s2+".dat").c_str());

        //std::ofstream file_xCells ("/tmp/jcarlos/tests_ur/cells.dat");

	cpm.setup(Nx, Ny, Nz, seed); // setup of CPM
	cpm.set_periodicity(h_periodic, v_periodic, d_periodic);

	cpm.init();   // add all cell types with their respective parameters
/*
	// pifUrot3D_CPM_large.txt
	//cpm.add_type(120, 1536e3, 0, ltarget01, 0, &J01[0], 0);		// stroma 25 lambdal1
	cpm.add_type(120, 1504e3, 0, ltarget01, 0, &J01[0], 0);		// stroma 25 lambdal1
	//cpm.add_type(120, 32e3, 200000, 32e3, 0, &J02[0], 0);  	// basal membrane 2500 20000 || 120 200k
	cpm.add_type(120, 2*32e3, 200000, 2*32e3, 0, &J02[0], 0);  	// basal membrane 2500 20000 || 120 200k
	cpm.add_type(120, 125, 0, ltarget03, 0, &J03[0], Chi); 	// basal lambdal2
	cpm.add_type(120, 125, 0, ltarget04, 0, &J04[0], Chi); 		// stem 500 100
	cpm.add_type(120, 1e3, 0, ltarget05, 0, &J05[0], Chi); 	// intermediate lambdal2
	cpm.add_type(120, 10e3, 100, 4900, 0, &J06[0], Chi); 	// umbrella lambdal2=1000
	cpm.add_type(5, 6048e3, 0, ltarget07, 0, &J07[0], 0); 		// urine 120 lambdal2
	//cpm.add_type(1, 6048e3, 0, ltarget07, 0, &J07[0], 0); 		// urine 120 lambdal2
	cpm.add_type(2000, 125, 0, ltarget08, 0, &J08[0], Chi); 	// tumor, lambdav=1800  500 | 150 0
	cpm.add_type(250, 8e3, 0, ltarget09, 0, &J09[0], 0); 	//  capillary lambdal2 | 250 0
*/

	// pifUrot3D_CPM_large.txt
	cpm.add_type(5e7, 1504e3, 0, ltarget01, 0, &J01[0], 0);		// stroma 25 lambdal1
	cpm.add_type(2400, 2*32e3, 1.5e9, 2*32e3, 0, &J02[0], 0);  	// basal membrane 2500 20000 || 120 200k
	cpm.add_type(750, 125, 0, ltarget03, 0, &J03[0], Chi); 		// basal lambdal2
	cpm.add_type(750, 125, 0, ltarget04, 0, &J04[0], Chi); 		// stem 500 100
	cpm.add_type(2400, 1e3, 0, ltarget05, 0, &J05[0], Chi); 	// intermediate lambdal2
	cpm.add_type(2400, 1e4, 100, 4900, 0, &J06[0], Chi); 		// umbrella lambdal2=1000
	cpm.add_type(2400, 6048e3, 0, ltarget07, 0, &J07[0], 0); 	// urine 120 lambdal2
	cpm.add_type(3000, 125, 0, ltarget08, 0, &J08[0], Chi); 	// tumor, lambdav=1800  500 | 150 0
	cpm.add_type(1e5, 8e3, 0, ltarget09, 0, &J09[0], 0); 		//  capillary lambdal2 | 250 0
/*
	// pifUrot3D_CPM_large.txt: DIFFERENT dH_vol2
	cpm.add_type(2.2e-5, 1504e3, 0, ltarget01, 0, &J01[0], 0);	// stroma 25 lambdal1
	cpm.add_type(5.9e-7, 2*32e3, 2e6, 2*32e3, 0, &J02[0], 0);  	// basal membrane 2500 20000 || 120 200k
	cpm.add_type(0.1536, 125, 0, ltarget03, 0, &J03[0], Chi); 	// basal lambdal2
	cpm.add_type(0.1536, 125, 0, ltarget04, 0, &J04[0], Chi); 	// stem 500 100
	cpm.add_type(2.4e-3, 1e3, 0, ltarget05, 0, &J05[0], Chi); 	// intermediate lambdal2
	cpm.add_type(2.4e-5, 10e3, 100, 4900, 0, &J06[0], Chi); 	// umbrella lambdal2=1000
	cpm.add_type(6.6e-11, 6048e3, 0, ltarget07, 0, &J07[0], 0); 	// urine 120 lambdal2
	cpm.add_type(0.192, 125, 0, ltarget08, 0, &J08[0], Chi); 	// tumor, lambdav=1800  500 | 150 0
	cpm.add_type(1.56e-3, 8e3, 0, ltarget09, 0, &J09[0], 0); 	//  capillary lambdal2 | 250 0
*/
	srand(seed);

  	std::ifstream fs;
  	//fs.open("pifUrot3D_CPM.txt");  // file with the cell geometry
  	fs.open("pifUrot3D_CPM_large.txt");  // file with the cell geometry
  	//fs.open("pifUrotVol2DVascularCPM2.txt");  // file with the cell geometry (new geometry)
  	int xcell, xtype, xi, xf, yi, yf, zi, zf;
  	int index;
  	std::string nhe;
  	// read the text file with data line by line
  	std::string Linha;

  	while (getline(fs, Linha)) {
    		std::stringstream sLinha;
    		sLinha << Linha;

    		sLinha >> xcell;
    		sLinha >> xtype; // cell type
    		sLinha >> xi; // cell y start
    		sLinha >> xf; // cell x end
    		sLinha >> yi; // cell y start
    		sLinha >> yf; // cell y end
    		sLinha >> zi; // cell z start
    		sLinha >> zf; // cell z end
/*		// pifUrot3D_CPM.txt
		if ( (tumor_layer==1) && (xcell==125) ) xtype = 8; // Basal
		if ( (tumor_layer==2) && (xcell==343) ) xtype = 8; // Intermediate 1
		if ( (tumor_layer==3) && (xcell==423) ) xtype = 8;
		if ( (tumor_layer==4) && (xcell==503) ) xtype = 8;
		if ( (tumor_layer==5) && (xcell==583) ) xtype = 8;
		if ( (tumor_layer==6) && (xcell==663) ) xtype = 8;
		if ( (tumor_layer==7) && (xcell==727) ) xtype = 8; // Umbrella
*/
		// pifUrot3D_CPM_large.txt
		if ( (tumor_layer==1) && (xcell==603) ) xtype = 8; // Basal
		if ( (tumor_layer==2) && (xcell==1423) ) xtype = 8; // Intermediate 1
		if ( (tumor_layer==3) && (xcell==1743) ) xtype = 8;
		if ( (tumor_layer==4) && (xcell==2063) ) xtype = 8;
		if ( (tumor_layer==5) && (xcell==2383) ) xtype = 8;
		if ( (tumor_layer==6) && (xcell==2707) ) xtype = 8;
		if ( (tumor_layer==7) && (xcell==2887) ) xtype = 8; // Umbrella

		if (xtype==8) std::cout << " tumor_layer " << xcell << "  " << xtype << "  " << xi << "  " << xf << "  " << yi << "  " << yf << "  " << zi << "  " << zf << std::endl;
    		///cpm.add_cell(index/Ny, index%Ny, types[i]);
    		cpm.add_cell_vol(xi, xf, yi, yf, zi, zf, xtype, Ny);
  	}
  	fs.close();


	int N2 = 0, tumor_volume; // número de células do tipo 1


	void** data = new void* [5];
	data[0] = &cpm;

	cpm.set_dataref(&data[0]);

	int i, ndead1 = 0, ndead2 = 0, nchanged = 0, n_new = 0, n_to_replace = 0;
	//std::cout << " ndead1 " << ndead1 << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (i=0; i<n_MCS; ++i) { // MAIN cycle (on the number of simulation steps)

		auto start1 = std::chrono::high_resolution_clock::now();
/*
               	if (i%10==0) {  // printout
                        sprintf(buffer, "%d", ifile);
                        std::string id_si(buffer);
                        std::string term("tests_ur3D/test"+id_si+".dat");
                        cpm.print_types(term);  // types
  		}
*/
		auto start2 = std::chrono::high_resolution_clock::now();
		///cpm.MCS(dH_durotaxis); // pixel copy atempts
		cpm.MCS(); // pixel copy atempts
		//std::cout << " i " << i << std::endl;
		auto start3 = std::chrono::high_resolution_clock::now();

		tumor_volume = cpm.proliferation(vol_mitosis);  // try cell proliferation
		auto start4 = std::chrono::high_resolution_clock::now();
		//std::cout << " colony_area " << colony_area << std::endl;

		ncell = cpm.n_cell();
		auto start5 = std::chrono::high_resolution_clock::now();
		//std::cout << " ncell " << ncell << std::endl;
		///std::cout << " MCS = " << i << " , n_cell = " << ncell << " , colony area = " << colony_area << std::endl;

/* no oxygen run
		cpm.diffusion(Nx, Ny, Nz);

		if (i>1) {
			ndead1 = cpm.checkO2(0.0001);  // death from O2 starvation
			if( ndead1>0 ) std::cout << " ndead1 = " << ndead1 << std::endl;
		}
*/
		//std::cout << " debug 1 " << std::endl;

		ndead2 = cpm.random_death(); // random cell death
		auto start6 = std::chrono::high_resolution_clock::now();
		n_to_replace = n_to_replace + ndead1 + ndead2;  // number of cells needed for replacement
		if( ndead2>0 ) std::cout << " random ndead2 = " << ndead2 << " n_to_replace " << n_to_replace << std::endl;

		cpm.death_pressure(); // cell death due to size...

/*
		if (i>700) {  // 500
			ndead2 = cpm.chemo_death(); // random tumor cell death due to chemotherapy
			if( ndead2>0 ) std::cout << " chemo_death ndead2 = " << ndead2 << std::endl;
		}
*/
/*
		if (i>500) {  // 500
			ndead2 = cpm.bcg_death(); // random tumor cell death due to immunotherapy (BCG), after ablation
			if( ndead2>0 ) std::cout << " bcg_death ndead2 = " << ndead2 << std::endl;
		}
*/
		//std::cout << " debug 2 " << std::endl;
		n_new = cpm.proliferation_stasis(n_to_replace); // to keep homeostasis
		auto start7 = std::chrono::high_resolution_clock::now();
		//std::cout << " debug 3 " << std::endl;
		n_to_replace = n_to_replace - n_new; // decrease the ones already duplicated

		nchanged = cpm.contact(); // change type due to contact interactions (differentiation)
		auto start8 = std::chrono::high_resolution_clock::now();
		if( nchanged>0 ) std::cout << " nchanged = " << nchanged << std::endl;

/*
		if (i%100==0) {  // printout
		//if (i%10==0) {  // printout
			std::cout << " MCS = " << i << std::endl;

			sprintf(buffer, "%d", i);
			std::string id_si(buffer);
			std::string term("_"+id_si+".dat");

			cpm.print_types(dir+"types18x"+term);  // types
		}
*/

		if (i%1==0) {  // printout
	                cpm.print_xCells(i, file_xCells);  // count number of cells
			ncell = cpm.n_cell();
			if(i>0) std::cout << " MCS = " << i << " , n_cell = " << ncell << " , tumor volume = " << tumor_volume <<
				     std::endl;
			//	     " , BM length = " << cpm.BM_length() << std::endl;
		}
		auto start9 = std::chrono::high_resolution_clock::now();

		if (i%100==0) {
                        sprintf(buffer, "%d", ifile);
                        std::string id_si(buffer);
                        std::string term("tests_ur3D/test"+id_si+".dat");
                        cpm.print_types(term);  // types

                        std::string term3("tests_ur3D/test3_"+id_si+".dat");
			cpm.print_tags(term3);  // tags
//		}

//		if (i%10==0) {  // printout
			/*
			sprintf(buffer, "%d", ifile);
			std::string id_si(buffer);
			std::string term("tests_ur3D/test"+id_si+".dat");
			cpm.print_types(term);  // types
			*/
			//cpm.print_outlinesX(term);  // types
			//std::cout << " ifile " << ifile << " term " << term << std::endl;

			///cpm.print_types("tests_ur/test3.dat");  // cell types
			//cpm.print_types("/tmp/jcarlos/tests_ur/test3.dat");  // cell types

			//cpm.print_tags("tests_ur3D/test3.dat");  // tags

			//cpm.print_O2("tests_ur/test3.dat");  // oxygen

	                //cpm.print_nCells(i,"tests_ur/cells.dat");  // count number of cells
	                ///cpm.print_xCells(i, file_xCells);  // count number of cells
/*
	                fprintf(gp, "set term png\n");
        	        fprintf(gp, "set output ('tests_ur3D/framey%d_0%d.png')\n", ifile, i);
                	fprintf(gp, "unset colorbox\n");  // new
                	fprintf(gp, "set size ratio 0.6\n");  // new
                	fprintf(gp, "unset xtics\n");  // new
                	fprintf(gp, "unset ytics\n");  // new
			fprintf(gp, "set cbrange [0:9]\n");
			fprintf(gp, "set palette defined ( 0 'black', 1 'white', 2 'cyan', 3 'blue', 4 'magenta', 5 'gray', 6 'green', 7 'yellow', 8 'orange', 9 'red' )\n");
	                //fprintf(gp, "splot 'tests_ur3D/test%d.dat' matrix u ($1+0.5):($2+0.5):($3+0.5):4 notitle w image\n", ifile);
			fprintf(gp, "splot 'tests_ur3D/test%d.dat' u 1:2:3:4 w pm3d notitle\n", ifile);
        	        fprintf(gp, "replot\n");
                	fflush(gp);
*/

/* no oxygen run
			cpm.print_O2("tests_ur/test2.dat");  // oxygen

	                fprintf(gp, "set term png\n");
        	        fprintf(gp, "set output ('tests_ur/frameo38_0%d.png')\n", i);
                	fprintf(gp, "unset colorbox\n");  // new
                	//fprintf(gp, "set size square 1,1\n");  // new
                	fprintf(gp, "set size ratio 0.6\n");  // new
                	fprintf(gp, "unset xtics\n");  // new
                	fprintf(gp, "unset ytics\n");  // new
	                fprintf(gp, "plot 'tests_ur/test2.dat' matrix u ($1+0.5):($2+0.5):3 notitle w image\n");
        	        fprintf(gp, "replot\n");
                	fflush(gp);
*/
		}
		auto start10 = std::chrono::high_resolution_clock::now();
/*
                if (i==500) {  // chirurgical ablation of papillary tumor 500
                        int ncells_ablation = cpm.ablation();
                        ///std::cout << " radical_ablation ! " << std::endl;
                        ///int ncells_ablation = cpm.radical_ablation(); // remove tumor and cells below
                        std::cout << " ncells_ablation = " << ncells_ablation << std::endl;
                }
*/

        	std::chrono::duration<double>elapsed1 = start2 - start1;
        	std::chrono::duration<double>elapsed2 = start3 - start2;
        	std::chrono::duration<double>elapsed3 = start4 - start3;
        	std::chrono::duration<double>elapsed4 = start5 - start4;
        	std::chrono::duration<double>elapsed5 = start6 - start5;
        	std::chrono::duration<double>elapsed6 = start7 - start6;
        	std::chrono::duration<double>elapsed7 = start8 - start7;
        	std::chrono::duration<double>elapsed8 = start9 - start8;
        	std::chrono::duration<double>elapsed9 = start10 - start9;
        	std::cout << " chrono t1 = " << elapsed1.count() << " t2 = " << elapsed2.count() << " t3 = " << elapsed3.count() <<
        			    " t4 = " << elapsed4.count() << " t5 = " << elapsed5.count() << " t6 = " << elapsed6.count() <<
        			    " t7 = " << elapsed7.count() << " t8 = " << elapsed8.count() << " t9 = " << elapsed9.count() << std::endl;

	}  // end MAIN i cycle (MCS)


	ncell = cpm.n_cell();
	std::cout << " n_cell = " << ncell << " , tumor volume = " << tumor_volume << " at " << i << std::endl;

        file_xCells.close();

	delete [] data;
	cpm.free_mem();
}
