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

struct cell {							// Initiates a struct called cell 

	type* tau;															// Creates a pointer variable called tau that stores the address in the memory of type   
	int vol;															// Defines a integer variable called vol
	double len;															// Defines a double variable called len
	////(double W0;
	int timeProlif;														// Defines a integer variable called timeProlif
	std::vector<int> voxels;											// Defines a std::vector data type that accepts integers and is called voxels
        std::set<int> cell_content;										// Defines a std::set data type that accepts integers and is called cell_content

	cell() : tau(NULL), vol(0), len(0), timeProlif(), voxels(0) {}		// Constructor for the cell struct that initializes the values to default; timeProlif() and cell_content is not initialized 
};				

struct pixel {							// Initiates a struct called pixel

	cell* tag;																	// Creates a pointer named tag that stores the cell' memory adress
	double O2;																	// Creates a double variable called O2
	////pixel* nbr [26];								
	pixel* nbr [6];																// Stores the pixel's memory adress in the nbr vector on the 7th position 

	////pixel() : tag(NULL), O2()  {int i; for (i=0; i<26; ++i) nbr[i] = NULL;}
	pixel() : tag(NULL), O2()  {int i; for (i=0; i<6; ++i) nbr[i] = NULL;}		// Constructor for the pixel struct that initializes the values to default
};


struct type {							// Initiates a struct called type

	int index;																	// Creates a integer variable called index

	double lambdaV;																// Creates a double variable called lambdaV
	double V_target;															// Creates a double variable called V_target 
	double lambdaL;																// Creates a double variable called lambdaL
	double L_target;															// Creates a double variable called L_target
	double Chi;																	// Creates a double variable called Chi
	double lambdaP;																// Creates a double variable called lambdaP
	double* J;																	// Pointer to a double 

	type() : lambdaV(0), V_target(0), lambdaL(0), L_target(0), Chi(0), lambdaP(0), J(NULL) {}	// Initiates the default values for the type struct 
};

struct l_cell {							// Initiates a struct called l_cell

	l_cell* next;																// Creates a pointer for l_cell called next
	cell CELL;																	// Creates a variable CELL of type cell

	l_cell() : next(NULL) {}													// Constructor for the l_cell struct that assigns default values
};

struct l_type {							// Initiates a struct called l_type

	l_type* next;																// Creates a pointer for l_type called next
	type TYPE;																	// Creates a variable TYPE of type type

	l_type() : next(NULL) {}													// Constructor for the l_type struct that assigns default values
};

class CPM {								// Creates a class named CPM
private:													// Starts a private section (not accessed to the user)

	double T;												// Creates a double variable called T
	int Nx, Ny, Nz, NE;										// Creates four integer variables called Nx, Ny, Nz, NE 
	pixel *PIXELS;											// Creates a PIXELS variable that stores the memory adress of the pixel struct
	l_type *t_list;											// Creates a t_list variable that stores the memory adress of the l_type struct
	int n_types;											// Creates a integer variable called n_types
	l_cell *c_list;											// Creates a c_list variable that stores the memory adress of the l_cell struct
	////std::vector<int> index_border;

	bool h_periodic, v_periodic, d_periodic;				// Creates three boolean variables called h_periodic, v_periodic and d_periodic 

	void** dataref;											// Declares a variable named dataref of type pointer to a pointer to void

        //int C[18][5] = {{2,5,7,0,0},  {2,5,9,0,0}, {2,5,11,0,0},  {2,5,13,0,0}, {4,1,2,3,4},  {2,7,13,0,0},{4,1,6,8,14},{2,7,9,0,0},
        //                {4,2,8,10,15},{2,9,11,0,0},{4,3,10,12,16},{2,11,13,0,0},{4,4,6,12,17},{2,7,18,0,0},{2,9,18,0,0},{2,11,18,0,0},
        //                {2,13,18,0,0},{4,14,15,16,17}};
        int C[18][5] = { {2,4,6,0,0},  {2,4,8,0,0}, {2,4,10,0,0}, {2,4,12,0,0}, {4,0,1,2,3},  {2,6,12,0,0},{4,0,5,7,13},{2,6,8,0,0},	// Defines a C vector of 18 rows and 5 columns and initializes it
                         {4,1,7,9,14}, {2,8,10,0,0},{4,2,9,11,15},{2,10,12,0,0},{4,3,5,11,16},{2,6,17,0,0},{2,8,17,0,0},{2,10,17,0,0},
                         {2,12,17,0,0},{4,13,14,15,16} };

public:

	cell* nullcell;															// Declares a variable nullcell that stores the memory adress of the cell variable

	void setup(int, int, int, int);											
	// Void function with 4 integers as inputs

	void setup(int nx, int ny, int nz){ setup(nx, ny, nz, 123); }			
	// Void function that accepts 3 integers as inputs and calls the previous setup function with the last value as 123

	void setup(int seed){setup(Nx, Ny, Nz, seed);}							
	// Void function that accepts 1 integer as input and calls the previous setup function 

	void setup(){setup(Nx, Ny, Nz, 123);}									
	// Void function that accepts no input and calls the previous setup function (Why so many setup functions?)
	
	void free_mem();														
	// Void function called free_mem to free the memory of the program
	
	void init();															
	// Void function called init (Does what??)
	
	void set_dataref(void** value){ dataref = value; }						
	// Void function called set_dataref that accpets the input pointer to a pointer to void and sets the variable dataref to the input value  

	void set_periodicity (bool valh, bool valv, bool vald) { h_periodic = valh; v_periodic = valv;  d_periodic = vald; }	
	// Void function called set_periodicity that accepts 3 boolean tye inputs and sets the variables defined in the private part of the class to the input values
	
	void set_T(double value) { T = value; }
	// Void function called set_T that takes a double variable as input and assigns the value of T (defined in the private) to the input 

	void add_type(double, double, double, double, double, double*, double);
	// Void function called add_type that takes the following inputs: 5 double type values, 1 pointer to double type and 1 double type value

	bool add_cell(int, int);
	// Void function called add_cell that takes as inputs two integer variables

	bool add_cell_vol(int, int, int, int, int, int, int, int);
	// Void function called add_cell_vol that takes as inputs 8 integer type values

	int proliferation(double);
	// Function called proliferation that takes an input of type double and returns an integer type value

	int proliferation_stasis(int);
  	// Function called proliferation_stasis that takes an input of type integer and also returns an integer type value

	bool mitosis(cell*);
	// Function called mitosis that takes as input the cell pointer and returns a bollean value

	int checkO2(double);
	// Function called checkO2 that takes as input a variable of type double and returns an integer value

	double O2average(cell*);
	// Function called O2average that takes as input the cell pointer and returns a double value

	void cell_death (cell*);
	// Void function called cell_death that takes as input the cell pointer

	int random_death();
	// Function called random_death that takes no inputs and returns an integer type value

	int chemo_death();
	// Function called chemo_death that takes no inputs and returns an integer type value

	int bcg_death();
	// Function called bcg_death that takes no inputs and returns an integer type value

	void remove_cell(cell*);
	// Void function called remove_cell that takes as input the cell pointer

	void remove_cell(pixel* p){remove_cell(p->tag);}
	// Void function called remove_cell that uses recursion to call it self ??????

	void remove_cell(int index){if (index<NE) remove_cell(&(PIXELS[index]));}
	// ????????????????????????????? 

	int contact();
	// Function called contact that takes no input and retutns an integer value

	int check_contact(cell*);
	// Function called check_contact that takes as input a pointer to the cell struct and returns an integer

	void change_type(cell*, int);
	// Void function called change_type that takes as input a pointer to the cell struct and an integer

	bool add_cell(int i, int j, int k, int type) { return add_cell(i*Ny*Nz+j*Nz+k, type); }
	// Boolean type function called add_cell that takes 4 integer inputs and ????????????????????????

	bool add_rand_cell(int);
	// Voolean type function called add_rand_cell that takes an integer input returning a boolean value

	void add_multiple_rand_cells(int, int);
	// Void function called add_multiple_rand_cells that takes two integer values as input

	void step(int, int);
	// Void function called step that takes two integers as inputs 

	//static double null_f(void** data,int a,int b){ return 0.0; }
	//void step(int index,int nbr) { step(index,nbr,null_f); }

	void rand_step();
	// Void function called rand_step that takes no inputs 

	//void rand_step() { rand_step(null_f); }
	//void take_step(pixel*, cell*, double, double);
        void take_step(pixel*, cell*, cell* ,int, double, double);
		
	void MCS() { for (int i=0; i<NE; ++i) rand_step(); }
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
        int n_cell();
        double BM_length();
        int ablation();
        int radical_ablation();
        void death_pressure();

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
