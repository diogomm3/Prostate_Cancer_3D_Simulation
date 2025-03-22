// Used to prevent recursive inclusions
#ifndef CPM_H							// Checks if the header file has been previously declared (if not generates true)
#define CPM_H							// If the previous result is true (header file is not declared) this statment will declare the header file
// In the end of the code is an #endif statment that closes the scope of #ifndef

#include <string>						// Includes the <string> library used to manipulate the string data type 
#include <cstdlib>						// Includes the <cstdlib> library wich defines several general purpose functions (memory management p.e.)
#include <fstream>						// Includes the <fstream> library that allows to create files, write information on files and read files p.e.
#include <iostream>						// Includes the <iostream> library that allows for input and output statments
#include <vector>						// Includes the <vector> library that allows for the user to manipulate the data type vectors
#include "mt19937ar.h"					// Includes the header file "mt19937ar.h" file
#include <set>							// Includes the <set> library that can take any data type and manipulate it

struct type;							// User defined data type to combine different standart data types 

struct cell {							
	type* tau;															
	int vol;															
	double len;															
	////(double W0;
	int timeProlif;														
	std::vector<int> voxels;											
	std::set<int> cell_content;											

	cell() : tau(NULL), vol(0), len(0), timeProlif(), voxels(0) {}		 
};				

struct pixel {							

	cell* tag;																	
	double O2;																	
	////pixel* nbr [26];								
	pixel* nbr [6];																

	////pixel() : tag(NULL), O2()  {int i; for (i=0; i<26; ++i) nbr[i] = NULL;}
	pixel() : tag(NULL), O2()  {int i; for (i=0; i<6; ++i) nbr[i] = NULL;}		
};


struct type {							

	int index;																	

	double lambdaV;																
	double V_target;															
	double lambdaL;																
	double L_target;															
	double Chi;																	
	double lambdaP;																
	double* J;																	

	type() : lambdaV(0), V_target(0), lambdaL(0), L_target(0), Chi(0), lambdaP(0), J(NULL) {}	 
};

struct l_cell {							

	l_cell* next;																
	cell CELL;																	

	l_cell() : next(NULL) {}													
};

struct l_type {							

	l_type* next;																
	type TYPE;																	

	l_type() : next(NULL) {}													
};

class CPM {								
private:													

	double T;												
	int Nx, Ny, Nz, NE;										
	pixel *PIXELS;											
	l_type *t_list;											
	int n_types;											
	l_cell *c_list;											

	bool h_periodic, v_periodic, d_periodic;				

	void** dataref;											

        //int C[18][5] = {{2,5,7,0,0},  {2,5,9,0,0}, {2,5,11,0,0},  {2,5,13,0,0}, {4,1,2,3,4},  {2,7,13,0,0},{4,1,6,8,14},{2,7,9,0,0},
        //                {4,2,8,10,15},{2,9,11,0,0},{4,3,10,12,16},{2,11,13,0,0},{4,4,6,12,17},{2,7,18,0,0},{2,9,18,0,0},{2,11,18,0,0},
        //                {2,13,18,0,0},{4,14,15,16,17}};
        int C[18][5] = { {2,4,6,0,0},  {2,4,8,0,0}, {2,4,10,0,0}, {2,4,12,0,0}, {4,0,1,2,3},  {2,6,12,0,0},{4,0,5,7,13},{2,6,8,0,0},	
                         {4,1,7,9,14}, {2,8,10,0,0},{4,2,9,11,15},{2,10,12,0,0},{4,3,5,11,16},{2,6,17,0,0},{2,8,17,0,0},{2,10,17,0,0},
                         {2,12,17,0,0},{4,13,14,15,16} };

public:

	cell* nullcell;															

	void setup(int, int, int, int);											
	
	void setup(int nx, int ny, int nz){ setup(nx, ny, nz, 123); }			
	
	void setup(int seed){setup(Nx, Ny, Nz, seed);}							
	
	void setup(){setup(Nx, Ny, Nz, 123);}									
		
	void free_mem();														
		
	void init();															
	
	void set_dataref(void** value){ dataref = value; }						
	
	void set_periodicity (bool valh, bool valv, bool vald) { h_periodic = valh; v_periodic = valv;  d_periodic = vald; }	
	
	void set_T(double value) { T = value; }
	
	void add_type(double, double, double, double, double, double*, double);
	
	bool add_cell(int, int);
	
	bool add_cell_vol(float, float, float, float, float, float, int, int, int);
	
	int proliferation(double);
	
	int proliferation_stasis(int);
  	
	bool mitosis(cell*);
	
	int checkO2(double);
	
	double O2average(cell*);
	
	void cell_death (cell*);
	
	int random_death();
	
	int chemo_death();
	
	int bcg_death();
	
	void remove_cell(cell*);
	
	void remove_cell(pixel* p){remove_cell(p->tag);}
	
	void remove_cell(int index){if (index<NE) remove_cell(&(PIXELS[index]));}
	
	int contact();

	int check_contact(cell*);
	
	void change_type(cell*, int);
	
	bool add_cell(int i, int j, int k, int type) { return add_cell(i*Ny*Nz+j*Nz+k, type); }
	
	bool add_rand_cell(int);
	
	void add_multiple_rand_cells(int, int);
	
	void step(int, int);

	//static double null_f(void** data,int a,int b){ return 0.0; }
	//void step(int index,int nbr) { step(index,nbr,null_f); }

	void rand_step();

	//void rand_step() { rand_step(null_f); }
	//void take_step(pixel*, cell*, double, double);
	void take_step(pixel*, cell*, cell* ,int, double, double);
		
	void MCS() { 
		for (int i=0; i<NE; ++i) {
			rand_step(); 
			// std::cout << i << "\n";
		}
	}
	////void MCS() { for (int i=0; i<index_border.size(); ++i) rand_step(); }
	//void MCS() { MCS(null_f); }

	bool check_step(pixel*, cell*, int) const;
	bool check_stem(cell*, int);
	bool CCA_check(pixel*) const;
	bool CCA_check3D(pixel*, int) const;

	double dH_vol(cell*, cell*) const;
	double dH_vol2(cell*, cell*) const;
	double dH_len(cell*, cell*, int, double*, double*) const;
	double dH_adh(pixel*, cell*) const;
	double dH_area(cell*, cell*, int) const;
	double dH_chemotaxis(int, pixel*, cell*) const;
	double dH_polar(cell*, pixel*, int) const;

	void print_tags(std::string) const;
	void print_types(std::string) const;
	void print_O2(std::string) const;
	//void print_strains(std::string) const;
	void print_outlines(std::string) const;
	void print_outlinesX(std::string) const;
        int* n_cell();
        double BM_length();
        int ablation();
        int radical_ablation();
        int death_pressure();

	pixel* access_pixel(int index) const{ return &(PIXELS[index]); }
	pixel* access_pixels() const{ return PIXELS; }
	l_cell* access_cell_list() const { return c_list; }
	int get_size(int coordinate) const { if(coordinate==0) return Nx; else if(coordinate==1) return Ny; else return Nz;}

        void get_content(cell* current);
        void get_list_content ();
        void print_nCells(int mcs, std::string filename);
        void print_xCells(int, std::ofstream&);
        void diffusion(int, int, int);

	double calc_length (cell*, int) const;
};

#endif
