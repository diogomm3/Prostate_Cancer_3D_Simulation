#include "CPM3D.h"												// Includes the header file "CPM3D.h"
#include <cmath>												// Includes <cmath> library that allows for basic math operations
#include "disjoint_set.h"										// Includes the header file "disjoint_set.h" (Does what?)
#include <fstream>												// Includes the <fstream> library that allows to create files, write information on files and read files p.e.
#include <cstdlib>												// Includes the <cstdlib> library wich defines several general purpose functions (memory management p.e.)
#include <iostream>												// Includes the <iostream> library that allows for input and output statments
#include <vector>												// Includes the <vector> library that allows for the user to manipulate the data type vectors
#include <chrono>												// Includes the <chrono> library which is used for time counting
#include "mt19937ar.h"											// Includes the header file "mt19937ar.h" that is used to generate random numbers 

#define PI 3.141592653589793238462643383279						// Defines the variable PI as the following value

void CPM::setup(int nx, int ny, int nz, int seed) { 			// Function to setup the CPM
	// Accepts as inputs the grid dimensions and a seed

	T = 1;						// Temperature of the model?
	Nx = nx;
	Ny = ny;
	Nz = nz;
	NE = Nx*Ny*Nz;				// Number of voxels is set to the dimensions limits product
	n_types = 0;				// Sets the number of type cells to 0
	c_list = NULL;				// Initializes the list of cells
	t_list = NULL;				// Initializes the list of cells types
	PIXELS = NULL;				// Initializes the list of pixels and sets it to e empty
	nullcell = NULL;			// Pointer that is currently not pointing to any cell
	init_genrand(seed);			// Initializes a random number generator using the seed

	// Sets the periodicity
	h_periodic = true;  		
	v_periodic = false;			
	d_periodic = true;			

	dataref = NULL;				// Pointer dataref is not poiting to any data
}

void CPM::free_mem(void) {  	// Function to free all memory

	l_cell* aux;				// Declares a aux pointer to 

	// Free the memory of all the cells in c_list
	while (c_list != NULL) {

		aux = c_list;			// Assigns the current c_list to hte aux pointer
		c_list = c_list->next;	// Moves to the next value of c_list
		delete aux;				// Deletes the aux variable

	}

	// Free the memory of all cell types in the t_list
	l_type *aux2, *curr;		// Declares two pointers to the cell type called aux2 and curr	
	double* Jaux;
	curr = t_list;				// Sets curr to the type list

	while (curr != NULL) {

		Jaux = curr->TYPE.J;	// Gets the address of the array for the current cell type
		curr = curr->next;		// Moves to the next adress
		delete [] Jaux;			// Deletes the address of the array for the current cell type

	}

	// Free the memory of all cell types in the t_list
	while (t_list != NULL) {

		aux2 = t_list;			// Set aux2 to the current type in the list
		t_list = t_list->next;	// Move to the next type in the list 
		delete aux2;			// Frees the memory occupied by the current type

	}

	// Free the memory occupied by the PIXELS array
	if (PIXELS!=NULL) delete [] PIXELS;
	// Free the memory occupied by the nullcell
	if (nullcell!=NULL) delete nullcell;
}

void CPM::init() { 				// Initialization of CPM

	int i, j, k, i1, j1, k1, n, ix1, ix2, iy1, iy2, iz1, iz2;
	pixel *curr, *nbr;			// Pointers to a pixel 
	++n_types;

	t_list = new l_type; 		// Creates the new list for the parameters of each cell type
	t_list->next = NULL;
	t_list->TYPE.index = 0;
	t_list->TYPE.lambdaV = 0;
	t_list->TYPE.lambdaL = 0;
	t_list->TYPE.lambdaP = 0;
	t_list->TYPE.V_target = 0;
	t_list->TYPE.L_target = 0;
	t_list->TYPE.J = new double [1];
	t_list->TYPE.J[0] = 0;

	nullcell = new cell;		// Creates an empty cell
	nullcell->vol = 0;
	nullcell->len = 0;
	nullcell->timeProlif = 0;
	nullcell->tau = &(t_list->TYPE);
	nullcell->voxels.clear();

	////index_border.clear();

	PIXELS = new pixel [NE];
	for(i=0; i<NE; ++i) {
		PIXELS[i].tag = nullcell;
		////for (j=0; j<26; ++j) PIXELS[i].nbr[j] = NULL;  // each cell can have many pixels
		for (j=0; j<6; ++j) PIXELS[i].nbr[j] = NULL;  // each cell can have many pixels
	}


	for(i=0; i<Nx; ++i) {
		ix1 = i + 1;	ix2 = i - 1;
		if (ix1==Nx) ix1 = 0;
		if (ix2==-1) ix2 = Nx - 1;
		for (j=0; j<Ny; ++j) {
			iy1 = j + 1;	iy2 = j - 1;
			if (iy1==Ny) iy1 = Ny - 1;
			if (iy2==-1) iy2 = 0;
			for (k=0; k<Nz; ++k) {
				iz1 = k + 1;	iz2 = k - 1;
				if (iz1==Nz) iz1 = 0;
				if (iz2==-1) iz2 = Nz - 1;
				curr = &PIXELS[i*Ny*Nz + j*Nz + k];
				/*
				n = 0;

				for (i1=i-1; i1<i+2; ++i1) {
					for (j1=j-1; j1<j+2; ++j1) {
						for (k1=k-1; k1<k+2; ++k1) {

							///nbr = &PIXELS[((Nx+i1)%Nx)*Ny + ((Ny+j1)%Ny)];
							nbr = &PIXELS[((Nx+i1)%Nx)*Ny*Nz + ((Ny+j1)%Ny)*Nz + ((Nz+k1)%Nz)];
							//std::cout << " n " << n << " k1 " << k1 << std::endl;
							if (nbr!=curr) {
								curr->nbr[n] = nbr;
								//nbr->nbr[25-n] = curr;
								++n;
							}
						} // k1
					} // j1
				} // i1
				*/

				//nbr = &PIXELS[((Nx+i-1)%Nx)*Ny*Nz + j*Nz + k];
				nbr = &PIXELS[ix2*Ny*Nz + j*Nz + k];
				curr->nbr[0] = nbr;	//nbr->nbr[5] = curr;
				//nbr = &PIXELS[((Nx+i+1)%Nx)*Ny*Nz + j*Nz + k];
				nbr = &PIXELS[ix1*Ny*Nz + j*Nz + k];
				curr->nbr[1] = nbr;     //nbr->nbr[4] = curr;
				//nbr = &PIXELS[i*Ny*Nz + ((Ny+j-1)%Ny)*Nz + k];
				nbr = &PIXELS[i*Ny*Nz + iy2*Nz + k];
				curr->nbr[2] = nbr;     //nbr->nbr[3] = curr;
				//nbr = &PIXELS[i*Ny*Nz + ((Ny+j+1)%Ny)*Nz + k];
				nbr = &PIXELS[i*Ny*Nz + iy1*Nz + k];
				curr->nbr[3] = nbr;     //nbr->nbr[2] = curr;
				//nbr = &PIXELS[i*Ny*Nz + j*Nz + ((Nz+k-1)%Nz)];
				nbr = &PIXELS[i*Ny*Nz + j*Nz + iz2];
				curr->nbr[4] = nbr;     //nbr->nbr[1] = curr;
				//nbr = &PIXELS[i*Ny*Nz + j*Nz + ((Nz+k+1)%Nz)];
				nbr = &PIXELS[i*Ny*Nz + j*Nz + iz1];
				curr->nbr[5] = nbr;     //nbr->nbr[0] = curr;

			} // k
		} // j
	} // i

}

void CPM::add_type(double lambdaV, double V_target, double lambdaL, double L_target, double lambdaP, double* J, double Chi) {
// add new cell type, with all the necessary parameters

	l_type* newtype = new l_type;
	newtype->next = NULL;
	newtype->TYPE.index = n_types;
	newtype->TYPE.lambdaV = lambdaV;
	newtype->TYPE.V_target = V_target;
	newtype->TYPE.lambdaL = lambdaL;
	newtype->TYPE.L_target = L_target;
	newtype->TYPE.lambdaP = lambdaP;
	newtype->TYPE.Chi = Chi;
	newtype->TYPE.J = new double [n_types + 1];

	int i;
	for (i=0; i<n_types+1; ++i) newtype->TYPE.J[i] = J[i];

	int n;
	double* Jaux = new double [n_types];
	double* aux;
	l_type* curr = t_list;
	l_type* prev = NULL;
	n = 0;

	while(curr) {

		for (i=0; i<n_types; ++i) Jaux[i] = curr->TYPE.J[i];
		aux = (curr->TYPE.J);
		delete [] aux;
		curr->TYPE.J = new double [n_types + 1];
		for (i=0; i<n_types; ++i) curr->TYPE.J[i] = Jaux[i];
		curr->TYPE.J[n_types] = J[n];

		++n;
		prev = curr;
		curr = curr->next;
	}

	prev->next = newtype;
	delete [] Jaux;
	std::cout << "\nAdd_type: n_types = " << n_types << "; V_target = " << V_target << std::endl;
	for (i=0; i<n_types+1; ++i) std::cout << J[i] << "  ";
	std::cout << "  " << std::endl;
	++n_types;

}


bool CPM::add_cell(int index, int newtau) { 							// Adds a new cell 

	bool a, b;
	int i;
	l_cell* aux;

	type* ptr;
	l_type* curr;

	if (newtau>0 && newtau < n_types) {

		curr = t_list;
		for (i=0; i<newtau; ++i) curr = curr->next;
		ptr = &(curr->TYPE);

		a = PIXELS[index].tag == nullcell;
		b = PIXELS[index].tag->vol > 1;

		if (a||b) {
			aux = new l_cell;

			aux->CELL.tau = ptr;
			aux->CELL.vol = 0;
			aux->CELL.len = 0;
			////aux->CELL.W0 = 0;
			aux->CELL.timeProlif = 0;
			///aux->CELL.cell_content.insert(index);

			aux->next = c_list;
			c_list = aux;

			PIXELS[index].tag = &(c_list->CELL);

			return true;
		}
	}

	return false;
}

bool CPM::add_cell_vol(float pi, float pf, float ti, float tf, float zi, float zf, int newtau, int cellnum, int Ny) { 		// Function to add the cells within the given area
	// The function accepts the radius, angle and height variations as well as 
	// the cell type and the number of voxels in the y axis 

	bool a, b;											// Booleans to check if a pixel is a nullcell or a cell with a volume greater than 1
	int i;

	l_cell* aux;										// aux is a pointer to a new cell
	type* ptr;											// ptr is a pointer to a specific cell type
	l_type* curr;										// curr is a pointer to a list of cells

	if (newtau > 0 && newtau < n_types) {											// Checks if the cell type is one of the given in the model

		// If the the cell type is valid assigns the cell type to ptr
		curr = t_list;
		for (i=0; i<newtau; ++i) curr = curr->next;									
		ptr = &(curr->TYPE);

		int indexy;										// Variable used to access a specific pixel

		a = PIXELS[indexy].tag == nullcell;				// True if the pixel is a nullcell
		b = PIXELS[indexy].tag->vol > 1;				// True if the pixel belongs to a cell with volume greater than 1

		if (a||b) {																	// Checks if the 
			aux = new l_cell;							// Adds a new cell

			// Initializes the values for the new cell
			aux->CELL.tau = ptr;				
			aux->CELL.vol = 0;
			aux->CELL.len = 0;
			aux->CELL.timeProlif = 0;
			aux->CELL.voxels.clear();

			
			aux->next = c_list;
			c_list = aux;

			// Calculates the center of the cell grid
      		double xmed = Nx/2;
			double ymed = Ny/2;
			double rho, phi;

			// Converts the angles to multiples of PI
			ti = ti*PI;  tf = tf*PI; 	

			// Cicle 
			for (int ix=0; ix<Nx; ix++) {
				for (int iy=0; iy<Ny; iy++) {
					for (int iz=0; iz<Nz; iz++) {
						rho = sqrt((ix-xmed)*(ix-xmed) + (iy-ymed)*(iy-ymed));
						phi = atan2(iy-ymed, ix-xmed);
   						if (phi < 0) {phi += 2*PI;}
						if (rho>=pi && rho<pf && phi>=ti && phi<tf && iz>=zi && iz<zf) {
							indexy = ix*Ny*Nz + iy*Nz + iz;
							PIXELS[indexy].tag = &(c_list->CELL);
							aux->CELL.vol++;
							aux->CELL.voxels.push_back(indexy);
						}
					}
				}
			}
			
			std::cout << "CellNum = " << cellnum << "; CellType = " <<aux->CELL.tau->index << "; CellVolume = " << aux->CELL.vol << "; CellVoxels = " << aux->CELL.voxels.size() << std::endl;

			// Creates the cartesian variables for the cell length and adds 200 (domain is only with positive values)
			int xi = 100+round(pi*cos(ti*PI));
			int xf = 100+round(pf*cos(tf*PI));
			int yi = 100+round(pi*sin(ti*PI));
			int yf = 100+round(pf*sin(tf*PI));

			// Determines the length of the cell
			double len = sqrt((xf-xi+1)*(xf-xi+1) + (yf-yi+1)*(yf-yi+1) + (zf-zi+1)*(zf-zi+1));
			aux->CELL.len = len;

			return true;
		}
	}

	return false;
}

int CPM::proliferation(double vol_mitosis) { // cell proliferation

	l_cell* current_node;
	cell* current_cell;

	int type, volume, prob, VolumeT = 0, flag_ur;
	double lambdaV, V_target, probC;
	double* J;
	bool flag;

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);
		type = current_cell->tau->index;

		if ( type == 5 ) {  // tumor cell
			volume = current_cell->vol;
			current_cell->timeProlif += 1;
			VolumeT += volume;\
			//lambdaV = current_cell->tau->lambdaV;
			V_target = current_cell->tau->V_target;
			//J = current_cell->tau->J;

			double u1 = genrand_real3(), u2 = genrand_real3();
			double timeLimit = 25 + 5*sqrt(-2*log(u1))*cos(2*PI*u2);  // gaussian mu=24, sigma=3
			///double timeLimit = 48 + 6*sqrt(-2*log(u1))*cos(2*PI*u2);  // gaussian mu=24, sigma=3
			///double timeLimit = 36 + 4.5*sqrt(-2*log(u1))*cos(2*PI*u2);  // gaussian mu=24, sigma=3

			//if (type == 8) std::cout << " V " << volume << " VT "<< V_target << " time " << current_cell->timeProlif << std::endl;

			//if ( (type == 8) && (volume > 0.90*V_target) && ((current_cell->timeProlif) > 24) ) {  // 90  tumor cell
			///if ( (type == 8) && (volume > 0.90*V_target) && ((current_cell->timeProlif) > timeLimit) ) {  // 90  tumor cell
			std::cout << "TimeLimit: " << timeLimit << ";\tCellTime: " << current_cell->timeProlif << ";\tCellVol: " << current_cell->vol << std::endl;
			if ( (volume > 1.5*V_target) && ((current_cell->timeProlif) > timeLimit) ) {  //
/*
				int ytumor = -999;
				for (int j=0; j<Ny; ++j) {
                			for (int i=0; i<Nx; ++i) {
						if (PIXELS[i*Ny+j].tag == current_cell) { // find pixel y-position of tumor cell
							ytumor = j;
							goto endloop; // jump to out of all cycles...
						}
					}
				}
endloop:
				if (ytumor<34+40) { //probability of proliferation depends on the distance to the capillary (@y=34 pixels)
					probC = 100; // close enough to capillary
				}
				else if (ytumor>34+80) { // too far way...
					probC = 0;
				}
				else {
					probC = 100 - 2.5*(ytumor-34-40); // probability decreases with distance to capillary
				}

				if (probC > ((int)(genrand_int31() % 100))) // proliferate with some probability
					flag = mitosis (current_cell);	//  5
*/
				flag = mitosis (current_cell);	//  5
			}
		}

		//std::cout << "Cell with address " << current_cell << " has volume " << volume << " and is of type " << type << ".\n";
		//std::cout << "Its type has volume deformation cost " << lambdaV << ", target volume " << V_target << " and the address where the adhesion costs are stored is at " << J << "\n";
		//std::cout << "It has J" << type <<"0 = "<< J[0] << " contact cost with the ECM and J" << type << type << " = " << J[type] << " with itself.\n";
		//std::cout << "\n";

		//std::cout << " proliferation Cell volume " << volume << ", target volume " << V_target << " and prob " << prob << "\n";

		current_node = current_node->next;
	}
	return VolumeT;

}

int CPM::proliferation_stasis(int n_dead) {  // homeostasis from stem cells

	l_cell* current_node;
	cell* current_cell;

	int type, volume, prob, n_new = 0, Time, n_possible = 0;
	bool flag;

	//std::cout << " proliferation_stasis n_new = " << n_new << " n_dead " << n_dead << std::endl;
	while (n_new<n_dead) {

		current_node = c_list;
		n_possible = 0;

		while (current_node) {

			current_cell = &(current_node->CELL);
			type = current_cell->tau->index;
			volume = current_cell->vol;

			if ( (type==4) && (volume>=4) ) {  // stem cell, minimum volume
				flag = false;
				Time = current_cell->timeProlif;

				prob = ((int)(genrand_int31() % (7000)));
				//if (type==4) std::cout << " proliferation_stasis prob = " << prob << " volume " << volume << " Time " << Time << " n_new " << n_new << " n_dead " << n_dead << "  n_possible " << n_possible << std::endl;

				if ( (type==4) && (n_new<n_dead) && (prob<volume*Time) ) {  //  stem cell, minimum time
					flag = mitosis (current_cell);	//  5
					//std::cout << " proliferation_stasis type = " << type << std::endl;
					if (flag) n_new++;
				}
				if ( (volume*Time>0) && (flag) ) n_possible++;
			}

			//std::cout << " proliferation Cell volume " << volume << ", target volume " << V_target << " and prob " << prob << "\n";

			current_node = current_node->next;
		}
		if (n_possible==0) return n_new;
	}

	return n_new;
}


bool CPM::mitosis (cell* current) { // cell division

	double x=0, y=0, z=0, n=0, xmin=Nx, xmax=0, ymin=Ny, ymax=0, zmin=Nz, zmax=0;
	int i, j, k, index, index_erase[5000];

	int type = current->tau->index;
/*
        for (int k=0; k<Nz; ++k) {

	        for (int j=0; j<Ny; ++j) {

	                for (int i=0; i<Nx; ++i) {
				if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current) {
*/
	for (int iv=0; iv < (current->voxels.size()); iv++) {
		index = current->voxels[iv];
		i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;

		x = x + i;   y = y + j;   z = z + k;   n++;
		if (i<xmin) xmin = i;   if (i>xmax) xmax = i;
		if (j<ymin) ymin = j;   if (j>ymax) ymax = j;
		if (k<zmin) zmin = k;   if (k>zmax) zmax = k;
		//std::cout << "i i " << i << " j " << j << " current " << current << std::endl;
	}

	double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
	double xmed = x/n, ymed = y/n, zmed = z/n, dx1 = dx, dz1 = dz;
	int inew, knew;

	if ( (dy==0) && (type==4) ) return false;
	if ( (dy==0) && (dx==0) && (dz==0)) return false;

	if ( (dx1>Nx/2) || (dz1>Nz/2) ) { // when the cell, due to periodic boundary conditions, cross the boundary...
		xmin = 2*Nx;  xmax = 0;  zmin = 2*Nz;  zmax = 0;
		x = 0;	z = 0;
		for (int iv=0; iv < (current->voxels.size()); iv++) {
			index = current->voxels[iv];
			i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;
			if (dx1>Nx/2) {
				inew = i;
				if (i<Nx/2) inew = i + Nx; // shift
				x = x + inew;
				if (inew<xmin) xmin = inew;   if (inew>xmax) xmax = inew;
			}
			if (dz1>Nz/2) {
				knew = k;
				if (k<Nz/2) knew = k + Nz; // shift
				z = z + knew;
				if (knew<zmin) zmin = knew;   if (knew>zmax) zmax = knew;
			}
		}
		if (dx1>Nx/2) {
			dx = xmax - xmin;  xmed = x/n;
		}
		if (dz1>Nz/2) {
			dz = zmax - zmin;  zmed = z/n;
		}
	}

	l_cell* aux;
        aux = new l_cell;
        aux->next = c_list;
        c_list = aux;

	cell* new_cell = &(aux->CELL);
	aux->CELL.voxels.clear();

	double len1 = 0, len2 = 0, lT;
	int i1, j1, i2, j2, n1 = 0;

/*
        for (int k=0; k<Nz; ++k) {
	        for (int j=0; j<Ny; ++j) {
        	        for (int i=0; i<Nx; ++i) {
				if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current) {
*/
	for (int iv=0; iv < (current->voxels.size()); iv++) {
		index = current->voxels[iv];
		i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;

		if (type==4) {
			if ( (dy>0) && (j<ymed) ) {
				PIXELS[index].tag = &(c_list->CELL);
				aux->CELL.voxels.push_back(index);
				//current->voxels.erase(current->voxels.begin()+iv);
				index_erase[n1] = iv;
				n1++;
			}
		}
		else {
			inew = i;
			if ( (dx1>Nx/2) && (i<Nx/2) ) inew = inew + Nx;
			knew = k;
			if ( (dz1>Nz/2) && (k<Nz/2) ) knew = knew + Nz;

			if ( (dx>=dy && dx>=dz && inew>xmed) || (dy>dx && dy>=dz && j>ymed) || (dz>dx && dz>dy && knew>zmed) ) {
				PIXELS[index].tag = &(c_list->CELL);
				aux->CELL.voxels.push_back(index);
				//current->voxels.erase(current->voxels.begin()+iv);
				index_erase[n1] = iv;
				n1++;
			}
		}

	}

	n = n1;

	std::cout << "Mitosis: vol_c = " <<  current->vol << "; vol_n = " << n << "; size = " << current->voxels.size() << std::endl;
        aux->CELL.tau = current->tau;
        aux->CELL.vol = n;
        ////aux->CELL.len = calc_length(new_cell, 0);
        aux->CELL.timeProlif = 0;

	current->vol = current->vol - n;
	////current->len = calc_length(current, 0);
	current->timeProlif = 0;
	//for (int iv=0; iv<n1; iv++) current->voxels.erase(current->voxels.begin()+index_erase[iv]);
	for (int iv=n1-1; iv>-1; iv--) current->voxels.erase(current->voxels.begin() + index_erase[iv]);  // erase voxels from mother cell
/*
	for (int iv=0; iv<n1-1; iv++) {
		if (index_erase[iv+1]<=index_erase[iv])
		std::cout << " ERROR mitosis index_erase = " << index_erase[iv] << "  " << index_erase[iv+1] << std::endl;  // check voxels ordering
	}
*/
	//std::cout << " mitosis len_c " <<  current->len << " len_n " << aux->CELL.len
	//          << " mitosis vol_c " <<  current->vol << " vol_n " << aux->CELL.vol << std::endl;
	//std::cout << " mitosis type " << type << " vol_c " <<  current->vol << " size " << current->voxels.size()
	//          << " new size " << aux->CELL.voxels.size() << std::endl;


	if (type==4) { //asymmetric division
		change_type(current, 3);
		//std::cout << " 4 mitosis n " << n << " ymin " << ymin << " ymax " << ymax << " ymed " << ymed << " xmin " << xmin << " xmax " << xmax << " xmed " << xmed << 
		//	" current->vol " << current->vol << " current->tau->index " << current->tau->index << std::endl;
	}

	return true;
}


double CPM::O2average (cell* current) {  // average of Oxygen concentration

	double totO2=0;
	int n = 0;

        for (int k=0; k<Nz; ++k) {

	        for (int j=0; j<Ny; ++j) {

        	        for (int i=0; i<Nx; ++i) {
				if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current) {
					totO2 = totO2 + PIXELS[i*Ny*Nz+j*Nz+k].O2;
					n++;
				}
			}
		}
	}

	if (n>0) totO2 = totO2/n;
	return totO2;
}

int CPM::checkO2(double O2limit) { // check if oxygen concentration is above limit

	l_cell* current_node;
	cell* current_cell;

	int ncell_dead = 0;
	double cellO2;
	cell* to_be_removed[1000];

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);

		if ( (current_cell->tau->index>=3) && (current_cell->tau->index<=6) ) {  // only normal cells can die
			cellO2 = O2average (current_cell);
			//std::cout << " checkO2 cellO2 = " << cellO2 << std::endl;

			if (cellO2<O2limit) {
				to_be_removed[ncell_dead] = current_cell;
				ncell_dead++;
			}
		}

		current_node = current_node->next;
	}

	for (int ir=ncell_dead-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

	return ncell_dead;

}


int CPM::random_death() { // cell random death

	l_cell* current_node;
	cell* current_cell;
	cell* to_be_removed[1000];

	int ncell_dead = 0, perc, nbasal=0;

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);

		if (current_cell->tau->index==3) {  // basal cells can die

			perc = ((int)(genrand_int31() % (4*500*90/10))); // 500*90/10
			if (perc<1) {  // 2 %
			//if (perc<0) {  // 2 %
				to_be_removed[ncell_dead] = current_cell;
				//for (int ir=0; ir<ncell_dead; ir++) {
				//	if (to_be_removed[ir] == current_cell) std::cout << " random_death ERROR repeated remove " << std::endl;
				//}
				ncell_dead++;
			}
			//nbasal++; std::cout << " perc " << perc << " nbasal " << nbasal << std::endl;
		}
		else if (current_cell->tau->index==5) {  // intermediate cells can die

			perc = ((int)(genrand_int31() % (4*500*250/10))); // 500*250/10
			if (perc<1) {  // 1 %
			//if (perc<0) {  // 1 %
				to_be_removed[ncell_dead] = current_cell;
				//for (int ir=0; ir<ncell_dead; ir++) {
				//	if (to_be_removed[ir] == current_cell) std::cout << " random_death ERROR repeated remove " << std::endl;
				//}
				ncell_dead++;
			}
		}
		else if (current_cell->tau->index==6) {  // umbrella cells can die

			perc = ((int)(genrand_int31() % (4*500*10/10))); // 500*10/10
			if (perc<1) {  // 5 %
			//if (perc<0) {  // 5 %
				to_be_removed[ncell_dead] = current_cell;
				//for (int ir=0; ir<ncell_dead; ir++) {
				//	if (to_be_removed[ir] == current_cell) std::cout << " random_death ERROR repeated remove " << std::endl;
				//}
				ncell_dead++;
			}
		}
/*		else if (current_cell->tau->index==8) {  // tumor cells can die

			perc = ((int)(genrand_int31() % (500*250/10))); // small probability
			if (perc<1) {  // 5 %
				to_be_removed[ncell_dead] = current_cell;
				ncell_dead++;
			}
		}
*/
		current_node = current_node->next;
	}

	for (int ir=ncell_dead-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

	return ncell_dead;

 }

int CPM::chemo_death() { // tumor cell random death due to chemotherapy

	l_cell* current_node;
	cell* current_cell;
	cell* to_be_removed[1000];

	int ncell_dead = 0, perc, nbasal=0;

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);

		if (current_cell->tau->index==8) {  // tumor cells can die

			perc = ((int)(genrand_int31() % (108))); // small probability 108
			if (perc<1) {  // 5 %
				to_be_removed[ncell_dead] = current_cell;
				ncell_dead++;
			}
		}

		current_node = current_node->next;
	}

	for (int ir=ncell_dead-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

	return ncell_dead;

 }


int CPM::bcg_death() { // tumor cell random death due to BCG immunotherapy (bladder surface treatment)

	l_cell* current_node;
	cell* current_cell;
	cell* to_be_removed[1000];

	int ncell_dead = 0, perc, nbasal=0, flag, ymax = 100;

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);

		if (current_cell->tau->index==8) {  // tumor cells can die

			flag = 0;
                        for (int k=0; k<Nz; ++k) {
    	                    for (int j=0; j<Ny; ++j) {
        	                        for (int i=0; i<Nx; ++i) {
                	                        if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
                        	                        if ( (j>ymax) && (flag==0) ) {  // is it on the top surface?
                                	                        flag = 1;
	                                                }
                                                }
                                        }
                                }
                        }

			perc = ((int)(genrand_int31() % (48))); // medium probability 48
			if ( (perc<1) && (flag==1) ) {  // 5 % only top tumor cells
				to_be_removed[ncell_dead] = current_cell;
				ncell_dead++;
			}
		}

		current_node = current_node->next;
	}

	for (int ir=ncell_dead-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

	return ncell_dead;

 }


void CPM::cell_death (cell* current) { // kill cell

	int npixels=0;
        for (int k=0; k<Nz; ++k) {
	        for (int j=0; j<Ny; ++j) {
        	        for (int i=0; i<Nx; ++i) {
				if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current) {
					PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index = 0;
					npixels++;
				}
			}
		}
	}
	std::cout << " cell_death current " << current << " vol " << current->vol << " npixels " << npixels << std::endl;

}

void CPM::remove_cell(cell* dying_cell) { // remove cell

	// should actually make no difference since its not in c_list.
	// and will not change anything in grid cleanup.
	// Better safe than sorry, though
	if (dying_cell==nullcell) return;

	dying_cell->vol = 0;
	// Since the 2nd part is optional,
	// This makes sure we can track dead cells,
	// if we don't remove them from c_list.

	/* Pixel grid cleanup */
	for (int i=0; i<NE; ++i) {
		if (PIXELS[i].tag==dying_cell)
			PIXELS[i].tag = nullcell;
	}

	dying_cell->voxels.clear();

	/* Remove from c_list */
	// Although it is not mandatory, this ensures the tag numbers
	// have no gaps when printing.
	// Easier to understand how many cells are there, by just counting
	// the number of items in the list. (or highest printed tag)

	if (c_list==NULL) return; // guard

	l_cell *aux;
	if (&(c_list->CELL) == dying_cell) {

		aux = c_list;
		c_list = aux->next;
		delete aux;

		return;
	}
	// If the program makes it here, then it is not the first item.
	// We can do,

	l_cell *curr, *prev;
	prev = c_list;
	curr = prev->next;
	while (curr!=NULL) {

		if (&(curr->CELL) == dying_cell) {

			prev->next = curr->next; // remove from list
			delete curr;
			return;
		}

		prev = curr;
		curr = prev->next;
	}

	return;
}

int CPM::contact() {  // contact interation (differentiation of cells by contact)

	l_cell* current_node;
	cell* current_cell;

	int nchanged = 0, jup, jdo, ile, iri, kfr, kba, flag_ur, flag_ba, type, i, j, k, index;

	current_node = c_list;

	while (current_node) {

		current_cell = &(current_node->CELL);
		type = current_cell->tau->index;

		if ( (type==3) || (type==5) || (type==6) ) {  // only basal or intermediate or umbrella
			flag_ur = 0;   flag_ba = 0; // flag urine and basal membrane
/*
		        for (int k=0; k<Nz; ++k) {
			        for (int j=0; j<Ny; ++j) {
        	        		for (int i=0; i<Nx; ++i) {
						if (PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
*/
			for (int iv=0; iv<current_cell->voxels.size(); iv++) {
				index = current_cell->voxels[iv];
				i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;

				jup = j+1; if (jup==Ny) jup = Ny-1;
				jdo = j-1; if (jdo==-1) jdo = 0;
				iri = i+1; if (iri==Nx) iri = 0;
				ile = i-1; if (ile==-1) ile = Nx-1;  // periodicity
				kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
				kba = k-1; if (kba==-1) kba = Nz-1;
				/* int index26[26] = {ile*Ny*Nz+jup*Nz+kfr, ile*Ny*Nz+jup*Nz+k,   ile*Ny*Nz+jup*Nz+kba, ile*Ny*Nz+j*Nz+kfr, ile*Ny*Nz+j*Nz+k,
						   ile*Ny*Nz+j*Nz+kba,   ile*Ny*Nz+jdo*Nz+kfr, ile*Ny*Nz+jdo*Nz+k,   ile*Ny*Nz+jdo*Nz+kba,
				                   i*Ny*Nz+jup*Nz+kfr,   i*Ny*Nz+jup*Nz+k,     i*Ny*Nz+jup*Nz+kba,   i*Ny*Nz+j*Nz+kfr,
						   i*Ny*Nz+j*Nz+kba,     i*Ny*Nz+jdo*Nz+kfr,   i*Ny*Nz+jdo*Nz+k,     i*Ny*Nz+jdo*Nz+kba,
				                   iri*Ny*Nz+jup*Nz+kfr, iri*Ny*Nz+jup*Nz+k,   iri*Ny*Nz+jup*Nz+kba, iri*Ny*Nz+j*Nz+kfr, iri*Ny*Nz+j*Nz+k,
						   iri*Ny*Nz+j*Nz+kba,   iri*Ny*Nz+jdo*Nz+kfr, iri*Ny*Nz+jdo*Nz+k,   iri*Ny*Nz+jdo*Nz+kba};  // 26 neighbors
				if ( (type==5) || (type==6) || (type==3) ) {
					for (int ii=0; ii<26; ii++) { if(PIXELS[index26[ii]].tag->tau->index == 7) flag_ur = 1; }  // contact urine
					for (int ii=0; ii<26; ii++) { if(PIXELS[index26[ii]].tag->tau->index == 2) flag_ba = 1; }  // contact membrane
				}*/
				int index6[6] = {ile*Ny*Nz+j*Nz+k,
				                   i*Ny*Nz+jup*Nz+k, i*Ny*Nz+j*Nz+kfr,
						   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k,
				                 iri*Ny*Nz+j*Nz+k };  // 6 neighbors
				if ( (type==5) || (type==6) || (type==3) ) {
					for (int ii=0; ii<6; ii++) { if(PIXELS[index6[ii]].tag->tau->index == 7) flag_ur = 1; }  // contact urine
					for (int ii=0; ii<6; ii++) { if(PIXELS[index6[ii]].tag->tau->index == 2) flag_ba = 1; }  // contact membrane
				}
			}

			if ( (flag_ur==1)  && (type!=6) ) { // not umbrella
				change_type(current_cell, 6); // change to umbrella
				std::cout << " contact change to umbrella, from " << type << " to " << current_cell->tau->index << std::endl;
				nchanged++;
				/*
				std::cout << " umbrella ";
        			for (k=0; k<Nz; ++k) { for (j=0; j<Ny; ++j) { for (i=0; i<Nx; ++i) {
					if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
						std::cout << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
					}
				} } }
				std::cout << " " << std::endl; */
			}
			if ( (flag_ur==0)  && (type==6) ) { // umbrella
				change_type(current_cell, 5); // change to intermediate
				std::cout << " U contact change to intermediate, from " << type << " to " << current_cell->tau->index << std::endl;
				nchanged++;
				/*
				std::cout << " intermediate ";
        			for (k=0; k<Nz; ++k) { for (j=0; j<Ny; ++j) { for (i=0; i<Nx; ++i) {
					if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
						std::cout << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
					}
				} } }
				std::cout << " " << std::endl; */
			}
			if ( (flag_ba==1) && (type!=3) ) { // not basal
				change_type(current_cell, 3);  // change to basal
				std::cout << " contact change to basal, from " << type << " to " << current_cell->tau->index << std::endl;
				nchanged++;
				/*
				std::cout << " basal ";
        			for (k=0; k<Nz; ++k) { for (j=0; j<Ny; ++j) { for (i=0; i<Nx; ++i) {
					if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
						std::cout << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
					}
				} } }
				std::cout << " " << std::endl; */
			}
			if ( (flag_ba==0) && (type==3) ) { // basal
				change_type(current_cell, 5);  // change to intermediate
				std::cout << " B contact change to intermediate, from " << type << " to " << current_cell->tau->index << std::endl;
				nchanged++;
				/*
				std::cout << " intermediate ";
        			for (k=0; k<Nz; ++k) { for (j=0; j<Ny; ++j) { for (i=0; i<Nx; ++i) {
					if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
						std::cout << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
					}
				} } }
				std::cout << " " << std::endl; */
			}
		}

		current_node = current_node->next;
	}
	return nchanged;

}


int CPM::check_contact(cell* current_cell) {  // check contact interation

	int jup, jdo, ile, iri, kfr, kba, flag_ur, type;

	//type = current_cell->tau->index;

	//if ( type==8 ) {  // only tumor cells
	flag_ur = 0;
	for (int k=0; k<Nz; ++k) {
	        for (int j=0; j<Ny; ++j) {
               		for (int i=0; i<Nx; ++i) {
				if (PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
					/* jup = j+1; if (jup==Ny) jup = Ny-1;
					jdo = j-1; if (jdo==-1) jdo = 0;
					ile = i+1; if (ile==Nx) ile = 0;  // periodicity
					iri = i-1; if (iri==-1) iri = Nx-1;
					kfr = i+1; if (kfr==Nz) kfr = 0;  // periodicity
					kba = i-1; if (kba==-1) kba = Nz-1;
					int index26[26] = {ile*Ny*Nz+jup*Nz+kfr, ile*Ny*Nz+jup*Nz+k,   ile*Ny*Nz+jup*Nz+kba, ile*Ny*Nz+j*Nz+kfr, ile*Ny*Nz+j*Nz+k,
							   ile*Ny*Nz+j*Nz+kba,   ile*Ny*Nz+jdo*Nz+kfr, ile*Ny*Nz+jdo*Nz+k,   ile*Ny*Nz+jdo*Nz+kba,
					                   i*Ny*Nz+jup*Nz+kfr,   i*Ny*Nz+jup*Nz+k,     i*Ny*Nz+jup*Nz+kba,   i*Ny*Nz+j*Nz+kfr,
							   i*Ny*Nz+j*Nz+kba,     i*Ny*Nz+jdo*Nz+kfr,   i*Ny*Nz+jdo*Nz+k,     i*Ny*Nz+jdo*Nz+kba,
					                   iri*Ny*Nz+jup*Nz+kfr, iri*Ny*Nz+jup*Nz+k,   iri*Ny*Nz+jup*Nz+kba, iri*Ny*Nz+j*Nz+kfr, iri*Ny*Nz+j*Nz+k,
							   iri*Ny*Nz+j*Nz+kba,   iri*Ny*Nz+jdo*Nz+kfr, iri*Ny*Nz+jdo*Nz+k,   iri*Ny*Nz+jdo*Nz+kba};  // 26 neighbors

					for (int ii=0; ii<26; ii++) { if(PIXELS[index26[ii]].tag->tau->index == 7) flag_ur = 1; }  // contact urine
					*/
					int index6[6] = {ile*Ny*Nz+j*Nz+k,
					                   i*Ny*Nz+jup*Nz+k, i*Ny*Nz+j*Nz+kfr,
							   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k,
					                 iri*Ny*Nz+j*Nz+k };  // 6 neighbors
					for (int ii=0; ii<6; ii++) { if(PIXELS[index6[ii]].tag->tau->index == 7) flag_ur = 1; }  // contact urine
				}
			}
		}
	}
	//}

	return flag_ur;

}


bool CPM::check_stem(cell* current_cell, int index_t) { // check stem cell state/position, keep in contact with BM

	l_cell* current_node;
	//cell* current_cell;

	int jup, jdo, ile, iri, kfr, kba, flag_mem = 0, type, index, i, j, k;

	current_node = c_list;

	//current_cell = &(current_node->CELL);
	type = current_cell->tau->index;
	//std::cout << " check_stem type " << type << " index_t " << index_t << std::endl;

	if ( type == 4 ) {  // only stem cell (niche in contact with basal membrane)
/*
	        for (int k=0; k<Nz; ++k) {
		        for (int j=0; j<Ny; ++j) {
        	       		for (int i=0; i<Nx; ++i) {
					index = i*Ny*Nz + j*Nz + k;
					if( (PIXELS[index].tag == current_cell) && (index!=index_t) ) {  // if exclude the target pixel...
*/
						//std::cout << " check_stem index " << index << " index_t " << index_t << std::endl;
						/*int index26[26] = {ile*Ny*Nz+jup*Nz+kfr, ile*Ny*Nz+jup*Nz+k,   ile*Ny*Nz+jup*Nz+kba, ile*Ny*Nz+j*Nz+kfr, ile*Ny*Nz+j*Nz+k,
								   ile*Ny*Nz+j*Nz+kba,   ile*Ny*Nz+jdo*Nz+kfr, ile*Ny*Nz+jdo*Nz+k,   ile*Ny*Nz+jdo*Nz+kba,
						                   i*Ny*Nz+jup*Nz+kfr,   i*Ny*Nz+jup*Nz+k,     i*Ny*Nz+jup*Nz+kba,   i*Ny*Nz+j*Nz+kfr,
								   i*Ny*Nz+j*Nz+kba,     i*Ny*Nz+jdo*Nz+kfr,   i*Ny*Nz+jdo*Nz+k,     i*Ny*Nz+jdo*Nz+kba,
						                   iri*Ny*Nz+jup*Nz+kfr, iri*Ny*Nz+jup*Nz+k,   iri*Ny*Nz+jup*Nz+kba, iri*Ny*Nz+j*Nz+kfr, iri*Ny*Nz+j*Nz+k,
								   iri*Ny*Nz+j*Nz+kba,   iri*Ny*Nz+jdo*Nz+kfr, iri*Ny*Nz+jdo*Nz+k,   iri*Ny*Nz+jdo*Nz+kba};  // 26 neighbors

						for (int ii=0; ii<26; ii++) { if(PIXELS[index26[ii]].tag->tau->index == 2) flag_mem = 1; }  // contact membrane
						*/

		for (int iv=0; iv < (current_cell->voxels.size()); iv++) {
			index = current_cell->voxels[iv];
			if ( index != index_t ) {
				i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;

				jup = j+1; if (jup==Ny) jup = Ny-1;
				jdo = j-1; if (jdo==-1) jdo = 0;
				iri = i+1; if (iri==Nx) iri = 0;
				ile = i-1; if (ile==-1) ile = Nx-1;  // periodicity
				kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
				kba = k-1; if (kba==-1) kba = Nz-1;

				int index6[6] = {ile*Ny*Nz+j*Nz+k,
				       	           i*Ny*Nz+jup*Nz+k, i*Ny*Nz+j*Nz+kfr,
						   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k,
				       	  	 iri*Ny*Nz+j*Nz+k };  // 6 neighbors
				for (int ii=0; ii<6; ii++) { if(PIXELS[index6[ii]].tag->tau->index == 2) flag_mem = 1; }  // contact membrane
			}
		}
	}

	return flag_mem;

}

void CPM::change_type(cell* current_cell, int new_type) { // cell differentiation

        if (!(new_type>0 && new_type<n_types) ) // if not valid
                return;

        l_type* current_type = t_list;
        while (new_type--)
		current_type = current_type->next;

        current_cell->tau = &(current_type->TYPE);
}


bool CPM::add_rand_cell(int type) { // add a random cell

	int a = 1, b = 1, c = 1;
	int index;

	index = ((int)(genrand_int31() % NE));
	int i = index/(Ny*Nz);
	int j = (index - i*Ny*Nz)/Nz;
	int k = index - i*Ny*Nz - j*Nz;

	///if (!h_periodic) a = (index/(Ny*Nz))%(Nx-1);
	///if (!v_periodic) b = (index%Ny)%(Ny-1);
	if (!h_periodic) a = i%(Nx-1);
	if (!v_periodic) b = j%(Ny-1);
	if (!d_periodic) c = k%(Nz-1);
	if (a*b*c) {

		return add_cell(index, type);
	}

	return false;
}

void CPM::add_multiple_rand_cells(int type, int N) { // add N random cell

	bool done;
	while(N) {

		done = add_rand_cell(type);
		if (done) --N;
	}
}

void CPM::step(int index_t, int nbr) { // pixel copy attempt

	//for (int i=0;i<Nx;i++){for (int j=0;j<Ny;j++){for (int k=0;k<Nz;k++){
	//	pixel* test = &PIXELS[i*Ny*Nz+j*Nz+k];
	//	int test1 = test->tag->tau->index;
	//	if( (test1<0)&&(test1>9))std::cout<<"test1 "<<test1<<" i "<<i<<" j "<<j<<" k "<<k<<std::endl;
	//}}}
	pixel* target = &PIXELS[index_t];
	pixel* source = target->nbr[nbr];
	pixel* curr;
	//std::cout << " index_t " << index_t << " nbr " << nbr << std::endl;

	double dH, p, lens, lent;
	int typet, types;
	int is, js, ks, it, jt, kt, index_s;
	cell *tcell, *scell;
	tcell = target->tag;
	scell = source->tag;

	typet = tcell->tau->index;
	types = scell->tau->index;

	// if ( (types==1) || (types==2) ) return; // source cannot be stroma or BM (passive fibers...)

	//if(types!=1&&types!=7){
	//	for (int iv; iv<scell->voxels.size(); iv++) std::cout<<" types "<<types<<" iv "<<iv<<" voxel "<<scell->voxels[iv]<<std::endl;
	//}
	//if ( (typet==4) || (types==4) ) return;
	/////if ( typet==4 ) return;
	//std::cout << " types " << types << std::endl;

	bool check = check_step(target, scell, index_t); // check if the step is possible
	bool check2;
	//if (((typet==1)||(typet==7))&&(types==8)) std::cout << " typet " << typet << " check " << check << std::endl;


	if (check) {

		//std::cout << " typet " << typet << " types " << types << " nbr " << nbr << std::endl;
		//std::cout << " dH_vol 1 " << dH_vol(tcell, scell) << std::endl;
		//std::cout << " dH_adh 1 " << dH_adh(target, scell) << std::endl;
		//std::cout << " dH_len 1 " << dH_len(tcell, scell, index_t, &lent, &lens) << std::endl;
		//if (((typet==1)||(typet==7))&&(types==8)) std::cout << " typet " << typet << " dH_vol " << dH_vol(tcell, scell) <<
		//	" dH_adh " << dH_adh(target, scell) << " dH_area " << dH_area(scell, tcell, index_t) << std::endl;

		////dH = dH_vol(tcell, scell) + dH_adh(target, scell) + dH_len(tcell, scell, index_t, &lent, &lens);  // change in the Hamiltonian (with different contributions)
		//dH = dH_vol(tcell, scell) + dH_adh(target, scell);  // change in the Hamiltonian (with different contributions)
		dH = dH_vol(tcell, scell) + dH_adh(target, scell);  // change in the Hamiltonian (with different contributions)
		//dH = dH_vol2(tcell, scell) + dH_adh(target, scell) + dH_area(scell, tcell, index_t);  // change in the Hamiltonian (with different contributions)
		//std::cout << " dH 2 " << dH << " types " << types << std::endl;
		
		if (dH>0) { // if it increases, can still have a chance...

			p = exp(-dH/T);
			if (p<genrand_real3()) check = false;
		}

		if (check) {
			typet = tcell->tau->index;
			//std::cout << " type " << type << " index_t " << index_t << std::endl;
			check2 = true;
			// if (typet==4) check2 = check_stem(tcell, index_t);

			//take_step(target, scell, lent, lens);
			 take_step(target, scell, tcell, index_t, lent, lens);
						// if (tcell->tau->index == 5) {std::cout << dH_vol(tcell, scell) << "\t" << dH_adh(target, scell) << "\t" << dH_area(scell, tcell, index_t) << "\n" << std::endl;}
				//if ((typet==8)||(types==8)) std::cout << " types " << types << " typet " << typet <<
				//	" dH " << dH << std::endl;
				/* it = index_t/(Ny*Nz); jt = (index_t - it*Ny*Nz)/Nz; kt = index_t - it*Ny*Nz - jt*Nz;
				for (int i=0; i<NE; ++i) {
					curr = &PIXELS[i];
					if (curr == source) {
						index_s = i;
						break;
					}
				}
				is = index_s/(Ny*Nz); js = (index_s - is*Ny*Nz)/Nz; ks = index_s - is*Ny*Nz - js*Nz;
				std::cout << " i " << is << " " << it << " j " << js << " " << jt <<" k " << ks << " " << kt <<
					     " tt " << typet << " ts " << types << std::endl;
				*/
			
		}
	}
}

void CPM::rand_step() {  // take a random pixel copy attempt
	int index_t;
	double r;
	
	index_t = ((int)(genrand_int31() % NE));	
	
	int i = index_t/(Ny*Nz);
	int j = (index_t - i*Ny*Nz)/Nz;
	int k = index_t - i*Ny*Nz - j*Nz;

	// Radius limit condition for voxels copy (Optimize the running time)
	r = sqrt(pow(Nx/2-i,2)+pow(Ny/2-j,2));
	// std::cout << r << std::endl;
	// if (r < 100 || r > 200) {return;}

	// Coordinate xx condition for voxels copy (Optimize the running time)
	// if (i > 100) return;
	// if (i > 200) return;
	
	int nbr;
	int a = 1, b = 1, c = 1;

	if (!h_periodic) a = i%(Nx-1);
	if (!v_periodic) b = j%(Ny-1);
	if (!d_periodic) c = k%(Nz-1);


	// If all the target pixel is at the edge the value of a, b, c is equal to 0
	if (a*b*c) {
		
		////nbr = ((int)(genrand_int31() % 26));
		nbr = ((int)(genrand_int31() % 6));
		//std::cout << " nbr " << nbr << std::endl;
		pixel* target = &PIXELS[index_t];
		pixel* source = target->nbr[nbr];

		if ((target->tag)!=(source->tag)) { // different cells
			step(index_t, nbr);
		}
	}
}

//void CPM::take_step (pixel* target, cell* scell, double lent, double lens) {
void CPM::take_step (pixel* target, cell* scell, cell* tcell, int index_t,  double lent, double lens) { // make the pixel copy

	if ( (target->tag->vol > 1) || (target->tag == nullcell) ) {
		////int jup, jdo, ile, iri, kfr, kba, i, j, k, flag, flag2, indexy, indexy1;

		if (target->tag != nullcell) {
			--(target->tag->vol); // target looses a pixel
			target->tag->len = lent;
			//if (target->tag->tau->index==2) std::cout << " take_step lent " << lent << std::endl;
		        ///tcell->cell_content.erase(index);
		}
		if (scell != nullcell) {
			++(scell->vol); // source gains a pixel
			scell->len = lens;
			//if (target->tag->tau->index==2) std::cout << " take_step lens " << lens << std::endl;
	        	///scell->cell_content.insert(index);
		}

		target->tag = scell;

		scell->voxels.push_back(index_t); // add pixel to source cell

		for (int iv=0; iv<tcell->voxels.size(); iv++) {
			if (tcell->voxels[iv]==index_t) {
				tcell->voxels.erase(tcell->voxels.begin() + iv); // remove target pixel from target cell
				break;
			}
		}
		//std::cout << " s vol " << scell->vol << "  " << scell->voxels.size() <<
		//	     " t vol " << tcell->vol << "  " << tcell->voxels.size() << std::endl;

		/*
		for (int iv=0; iv<tcell->voxels.size(); iv++) {
			indexy = tcell->voxels[iv];
			i = indexy/(Ny*Nz);   j = (indexy - i*Ny*Nz)/Nz;   k = indexy - i*Ny*Nz - j*Nz;

			jup = j+1; if (jup==Ny) jup = Ny-1;
			jdo = j-1; if (jdo==-1) jdo = 0;
			ile = i+1; if (ile==Nx) ile = 0;  // periodicity
			iri = i-1; if (iri==-1) iri = Nx-1;
			kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
			kba = k-1; if (kba==-1) kba = Nz-1;
			int index6[6] = {ile*Ny*Nz+j*Nz+k,
			                   i*Ny*Nz+jup*Nz+k, i*Ny*Nz+j*Nz+kfr,
					   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k,
			                 iri*Ny*Nz+j*Nz+k };  // 6 neighbors
			flag = 0;
			for (int ii = 0; ii<6; ii++) {
				for (int iv1=0; iv1<tcell->voxels.size(); iv1++) {
					indexy1 = tcell->voxels[iv1];
					if (index6[ii]==indexy1) flag++; // is it part of the same cell?
				}
			}
			flag2 = 0;
			for (int iv2=0; iv2<index_border.size(); iv2++) {
				if (indexy == index_border[iv2]) {
					flag2 = 1; // it already is on the list
					break;
				}
			}

			if ( (flag<6) && (flag2==0) ) index_border.push_back(indexy);
			if ( (flag==6) && (flag2==1) ) tcell->voxels.erase(tcell->voxels.begin()+iv);
		}


		for (int iv=0; iv<scell->voxels.size(); iv++) {
			indexy = scell->voxels[iv];
			i = indexy/(Ny*Nz);   j = (indexy - i*Ny*Nz)/Nz;   k = indexy - i*Ny*Nz - j*Nz;

			jup = j+1; if (jup==Ny) jup = Ny-1;
			jdo = j-1; if (jdo==-1) jdo = 0;
			ile = i+1; if (ile==Nx) ile = 0;  // periodicity
			iri = i-1; if (iri==-1) iri = Nx-1;
			kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
			kba = k-1; if (kba==-1) kba = Nz-1;
			int index6[6] = {ile*Ny*Nz+j*Nz+k,
			                   i*Ny*Nz+jup*Nz+k, i*Ny*Nz+j*Nz+kfr,
					   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k,
			                 iri*Ny*Nz+j*Nz+k };  // 6 neighbors
			flag = 0;
			for (int ii = 0; ii<6; ii++) {
				for (int iv1=0; iv1<scell->voxels.size(); iv1++) {
					indexy1 = scell->voxels[iv1];
					if (index6[ii]==indexy1) flag++;
				}
			}
			flag2 = 0;
			for (int iv2=0; iv2<index_border.size(); iv2++) {
				if (indexy == index_border[iv2]) {
					flag2 = 1; // it already is on the list
					break;
				}
			}

			if ( (flag<6) && (flag2==0) ) index_border.push_back(indexy);
			if ( (flag==6) && (flag2==1) ) scell->voxels.erase(tcell->voxels.begin()+iv);
		}
		*/


	}
}

bool CPM::check_step (pixel* target, cell* scell, int index_t) const { // check if pixel copy is possible

	cell* tcell;
	cell *curr, *prev;
	int n, i, typet, types;
	bool check;
	////cell* nbrs [26];
	cell* nbrs [6];

	tcell = target->tag;

	typet = tcell->tau->index;
	//types = scell->tau->index;
	//if ( (typet==4) || (types==4) ) return false;
	/////if ( typet==4 ) return false;
        //if (typet!=types) std::cout<< " typet " << typet << " types " << types << std::endl;

	//std::cout << " .types " << types << std::endl;

	if (tcell!=scell) { // source and target are from different cells

		if (tcell != nullcell) { // and not the medium

			if ((tcell->vol) > 1) { // and not disapear

				n = 0;

				///nbrs[0] = target->nbr[0]->tag;	nbrs[1] = target->nbr[1]->tag;	nbrs[2] = target->nbr[2]->tag;
				///nbrs[3] = target->nbr[4]->tag;	nbrs[4] = target->nbr[7]->tag;	nbrs[5] = target->nbr[6]->tag;
				///nbrs[6] = target->nbr[5]->tag;	nbrs[7] = target->nbr[3]->tag;

				////for (int in=0; in<26; in++) nbrs[in] = target->nbr[in]->tag; // ??
				for (int in=0; in<6; in++) nbrs[in] = target->nbr[in]->tag; // ??
				//std::cout << " .nbrs " << nbrs[25] << std::endl;

				////prev = nbrs[25]; // ??
				//prev = nbrs[5]; // ???

				////for (i=0; i<26; ++i) {
				for (i=0; i<6; ++i) {
					curr = nbrs[i];  // neighbor cells of target pixel
					//if (curr != prev && curr == tcell) {
					if (curr == tcell) {
						++n;  // number of target cells pixels on the 6 neighbors
					}
					//prev = curr; //???
				}
				//std::cout << " .n " << n << std::endl;

				////return true; // ??
				if ( (n==0) || (n==6) ) { return false; }  // can't do it, isolate pixel
				if ( (n==1) || (typet==1) || (typet==4) ) { return true; } // no connectivity in urine or stroma if only a pixel don't need to test connectivity!
				//if ( (typet==1) || (typet==7) ) return true; // no connectivity in urine or stroma
				else {
					check = CCA_check3D(target, index_t);	// connectivity
					return check;
				}

			} // tcell->vol > 1
			else {
				return false;
			}
		}
		else {	// no problem when target is ECM...
			return true;
		} // tcell != nullcell
	} // tcell!=scell

	return false;
}


bool CPM::CCA_check(pixel* target) const {	// cell connectivity check

	cell* this_cell = target->tag;
	pixel* curr;


        int type = this_cell->tau->index;
	//if (type==2) return false;   // dirty trick to avoid rupture of membrane. not nice (JC)
	/////if (type==4) return false;   // dirty trick to avoid rupture of membrane. not nice (JC)


	int i, j;
	int NE = Nx*Ny*Nz;
	int* label = new int [NE];
	for (i=0; i<NE; ++i) {

		curr = &PIXELS[i];

		if (curr == target) {

			label[i] = 0;
		}
		else if (curr->tag == this_cell) {

			label[i] = 1;
		}
		else {

			label[i] = 0;
		}
	}

	int N = 0;
	int index, index2;
	int k;

	if (label[0]) { ++N; label[0] = N; }
	for (j=1; j<Ny; ++j) {

		if (label[j]) {

			if (label[j-1]) label[j] = label[j-1];
			else { ++N; label[j] = N; }
		}
	}

	for (i=1; i<Nx; ++i) {

		if (label[i*Ny]) {

			if (label[(i-1)*Ny]) label[i*Ny] = label[(i-1)*Ny];
			else { ++N; label[i*Ny] = N; }
		}
	}

	i = Nx-1;
	for (j=1; j<Ny-1; ++j) {

		if (label[i*Ny+j]) {

			if (label[i*Ny+j-1]) label[i*Ny+j] = label[i*Ny+j-1];
			else { ++N; label[i*Ny+j] = N; }
		}
	}

	j = Ny-1;
	for (i=1; i<Nx-1; ++i) {

		if (label[i*Ny+j]) {

			if (label[(i-1)*Ny+j]) label[i*Ny+j] = label[(i-1)*Ny+j];
			else { ++N; label[i*Ny+j] = N; }
		}
	}

	i = Nx-1;
	j = Ny-1;
	if (label[i*Ny+j]) {

		if (label[i*Ny+j-1]) label[i*Ny+j] = label[i*Ny+j-1];
		else if (label[(i-1)*Ny+j]) label[i*Ny+j] = label[(i-1)*Ny+j];
		else { ++N; label[i*Ny+j] = N; }
	}

	for (i=1; i<Nx-1; ++i) {

		for (j=1; j<Ny-1; ++j) {

			index = i*Ny+j;
			if (label[index]) {

				for (k=0; k<4; ++k) {

					index2 = (i-1+(k/3))*Ny+(j-1+(k%3));

					if (label[index2]) {

						label[index] = label[index2];
						k = 10;
					}
				}

				if (k<10) {

					++N;
					label[index] = N;
				}
			}
		}
	}

	set* set_list = new set [N];
	set *a, *b;
	int ii, jj;

	for (i=0; i<Nx; ++i) {

		for (j=0; j<Ny; ++j){

			index = i*Ny + j;
			if (label[index]) {

				a = &set_list[label[index]-1];

				for (k=0; k<9; ++k) {

					if (k==4) ++k;
					ii = i - 1 + (k/3);
					jj = j - 1 + (k%3);
					index2 = ((ii+Nx)%Nx)*Ny + (Ny+jj)%Ny;
					if (label[index2]) {

						b = &set_list[label[index2]-1];
						a->unite(b);
					}
				}
			}
		}
	}

	delete [] label;

	int n = 0;
	set* curr_set;
	for (i=0; i<N; ++i) {

		curr_set = &set_list[i];
		if (curr_set->find() == curr_set) ++n;
	}

	delete [] set_list;

	//if (type==2) std::cout << " CCA_check n " << n << std::endl;

	if (n==1) return true;
	else return false;
}


bool CPM::CCA_check3D(pixel* target, int index_t) const {	// cell connectivity checkc in 3D

	cell* this_cell = target->tag;
	pixel* curr;

	int i, j, k, ile, iri, jup, jdo, kfr, kba;

	//vector<vector<vector<int> > > A (3,vector<vector<int> >(3,vector <int>(3,0)));
	//vector <int> x(18);
	//vector<vector<int> > C (18,vector <int>(5,0));

	////int C[18][5] = {{2,5,7,0,0},  {2,5,9,0,0}, {2,5,11,0,0},  {2,5,13,0,0}, {4,1,2,3,4},  {2,7,13,0,0},{4,1,6,8,14},{2,7,9,0,0},
	////		    {4,2,8,10,15},{2,9,11,0,0},{4,3,10,12,16},{2,11,13,0,0},{4,4,6,12,17},{2,7,18,0,0},{2,9,18,0,0},{2,11,18,0,0},
	////		    {2,13,18,0,0},{4,14,15,16,17}};

	//auto start1 = std::chrono::high_resolution_clock::now();
/*
	for (int ip=0; ip<NE; ip++) {
		curr = &PIXELS[ip];
		if (target==curr) {
			index_t = ip;
			break;
		}
	}
*/
	//auto start2 = std::chrono::high_resolution_clock::now();
	i = index_t/(Ny*Nz);
	j = (index_t - i*Ny*Nz)/Nz;
	k = index_t - i*Ny*Nz - j*Nz;

	jup = j+1; if (jup==Ny) jup = Ny-1;
	jdo = j-1; if (jdo==-1) jdo = 0;
	iri = i+1; if (iri==Nx) iri = 0;  // periodicity
	ile = i-1; if (ile==-1) ile = Nx-1;
	kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
	kba = k-1; if (kba==-1) kba = Nz-1;
/*
	int A[3][3][3] = {{{ile*Ny*Nz+jup*Nz+kfr, ile*Ny*Nz+j*Nz+kfr, ile*Ny*Nz+jdo*Nz+kfr},
	      		   {i*Ny*Nz+jup*Nz+kfr,   i*Ny*Nz+j*Nz+kfr,   i*Ny*Nz+jdo*Nz+kfr},
   	      		   {iri*Ny*Nz+jup*Nz+kfr, iri*Ny*Nz+j*Nz+kfr, iri*Ny*Nz+jdo*Nz+kfr}},
	     		  {{ile*Ny*Nz+jup*Nz+k,   ile*Ny*Nz+j*Nz+k,   ile*Ny*Nz+jdo*Nz+k},
	      		   {i*Ny*Nz+jup*Nz+k,     i*Ny*Nz+j*Nz+k,     i*Ny*Nz+jdo*Nz+k},
   	      		   {iri*Ny*Nz+jup*Nz+k,   iri*Ny*Nz+j*Nz+k,   iri*Ny*Nz+jdo*Nz+k}},
	     		  {{ile*Ny*Nz+jup*Nz+kba, ile*Ny*Nz+j*Nz+kba, ile*Ny*Nz+jdo*Nz+kba},
	      		   {i*Ny*Nz+jup*Nz+kba,   i*Ny*Nz+j*Nz+kba,   i*Ny*Nz+jdo*Nz+kba},
   	      		   {iri*Ny*Nz+jup*Nz+kba, iri*Ny*Nz+j*Nz+kba, iri*Ny*Nz+jdo*Nz+kba}}};  // 27 voxels
*/

        int A[3][3][3] = { { {-1,                   ile*Ny*Nz+jdo*Nz+k, -1},
                             {ile*Ny*Nz+j*Nz+kba,   ile*Ny*Nz+j*Nz+k,   ile*Ny*Nz+j*Nz+kfr},
                             {-1,                   ile*Ny*Nz+jup*Nz+k, -1} },
                           { {i*Ny*Nz+jdo*Nz+kba,   i*Ny*Nz+jdo*Nz+k,   i*Ny*Nz+jdo*Nz+kfr},
                             {i*Ny*Nz+j*Nz+kba,                 -1,     i*Ny*Nz+j*Nz+kfr},
                             {i*Ny*Nz+jup*Nz+kba,   i*Ny*Nz+jup*Nz+k,   i*Ny*Nz+jup*Nz+kfr} },
                           { {-1,                   iri*Ny*Nz+jdo*Nz+k, -1},
                             {iri*Ny*Nz+j*Nz+kba,   iri*Ny*Nz+j*Nz+k,   iri*Ny*Nz+j*Nz+kfr},
                             {-1,                   iri*Ny*Nz+jup*Nz+k, -1} } };  // 27 voxels

        for (int ii=0; ii<3; ii++) { for (int jj=0; jj<3; jj++) { for (int kk=0; kk<3; kk++) {
                if (A[ii][jj][kk]>-1) {
                        curr = &PIXELS[A[ii][jj][kk]];
			A[ii][jj][kk] = 0;
			if (curr->tag != nullcell) {
                        	if (curr->tag == this_cell) A[ii][jj][kk] = 1; // this pixel belongs to the cell
                        	else A[ii][jj][kk] = 0;
	                }
                }
        } } }


    	// pixels ordering
    	int x[18] = { A[0][1][0], A[1][2][0], A[2][1][0], A[1][0][0], A[1][1][0],
	     	      A[0][0][1], A[0][1][1], A[0][2][1], A[1][2][1], A[2][2][1],
	     	      A[2][1][1], A[2][0][1], A[1][0][1], A[0][1][2], A[1][2][2],
	     	      A[2][1][2], A[1][0][2], A[1][1][2] };


	int vpixel[18] = {0};
    	int np = 0;
    	for (int ip=0; ip<18; ip++) { // count number of on pixels
        	if (x[ip]==1) {
            		vpixel[np] = ip;
            		np++;
		}
	}

	if (np<2) std::cout << " CCA_check3D np = " << np << std::endl;

        int vflag[18] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        //std::fill(vflag,vflag+(end(vflag) - begin(vflag)),1);

        //std::cout << " vflag[0] " << vflag[0] << " vflag[1] " << vflag[1] << " vflag[2] " << vflag[2] << std::endl;
        //std::cout << " C[0][0] " << C[0][0] << " C[0][1] " << C[0][1] << " C[0][2] " << C[0][2] <<
        //             " C[1][0] " << C[1][0] << " C[1][1] " << C[1][1] << " C[1][2] " << C[1][2] << std::endl;

	//auto start3 = std::chrono::high_resolution_clock::now();

	// the algorithm! 2
	vflag[0] = 0;
	int change = 1, nbr;
	while (change>0) {
    		change = 0;
    		for (int ip=0; ip<np; ip++) {
        		if (vflag[ip]==0) {
            			for (int inbr=0; inbr<C[vpixel[ip]][0]; inbr++) {
                			nbr = C[vpixel[ip]][inbr+1];
                    			for (int ipn=0; ipn<np; ipn++) {
                        			if ( (ip!=ipn) && (vpixel[ipn]==nbr) && (vflag[ipn]!=0) ) {
                            				vflag[ipn] = 0;
                            				change = 1;
                        			}
                    			}
            			}
        		}
    		}
	}
/*
	auto start4 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double>elapsed1 = start2 - start1;
	std::chrono::duration<double>elapsed2 = start3 - start2;
	std::chrono::duration<double>elapsed3 = start4 - start3;
	std::cout << " chrono t1 = " << elapsed1.count() << " t2 = " << elapsed2.count() << " t3 = " << elapsed3.count() << std::endl;
*/
    	// check decision
    	int flag_con = 0;
    	for (int ip=0; ip<np; ip++) {
        	if (vflag[ip]!=0) {
            		flag_con = 1;
        	}
    	}

     	if (flag_con==0) return true;
	else return false;

}


double CPM::dH_vol(cell* tcell, cell* scell) const {	// change of energy due to change of cell volume

	double dh = 0;
	int vol, type, types = 0;
	double V_T, Vmax = 50, k = 0.2;
	double l;

	if (scell != nullcell) {

		vol = scell->vol;
		V_T = scell->tau->V_target;
		l = scell->tau->lambdaV;

		types = scell->tau->index;
		if (types==5) { // tumor cell with dynamic target volume
			int VT = scell->tau->V_target;
			V_T = VT + 20000*(scell->timeProlif); // linear growth
			if (V_T>VT*2) V_T = VT*2; // volume limit	
			// std::cout << "V_T: " << V_T << std::endl;
		}

		dh += l*(1+2*(vol-V_T))/(V_T*V_T); // Hamiltonian change
	}

	if (tcell != nullcell) {

		vol = tcell->vol;
		V_T = tcell->tau->V_target;
		l = tcell->tau->lambdaV;

		type = tcell->tau->index;
		if (type==5) { // tumor cell with dynamic target volume
			int VT = tcell->tau->V_target;
			V_T = VT + 20000*(tcell->timeProlif); // linear growth
			if (V_T>VT*2) V_T = VT*2; // volume limit
			// std::cout << "V_T: " << V_T << std::endl;
		}

		dh += l*(1-2*(vol-V_T))/(V_T*V_T); // Hamiltonian 
	}
	return dh;
}

double CPM::dH_len(cell* tcell, cell* scell, int index_t, double* lent, double* lens) const {   // cell length. needs to be optimized! and checked...

	double dh = 0;
	double L_T;
	double l, lT, len, len1;
	//int it, jt, index, index1, i1, j1, i2, j2;
	//std::cout << " len 1 " << len << std::endl;
	if (scell != nullcell) { // not medium

		len = scell->len;
		L_T = scell->tau->L_target;
		l = scell->tau->lambdaL;
		//len1 = 0;
		//std::cout << " len 2 " << len << std::endl;

		if (l>0) {
			len1 = calc_length(scell, index_t);

			//std::cout << " len1 1 " << len1 << std::endl;
			//if (len1<len) {
			//	len1 = len;
				dh += l*( len1*len1 - len*len - 2*L_T*(len1-len) )/(L_T*L_T);
				*lens = len1;
			//}
		} // if l>0
		//if (scell->tau->index==2) std::cout << " dH_len s len1 = " << len1 << " l " << l << std::endl;
	}

	if (tcell != nullcell) { // not medium

		len = tcell->len;
		L_T = tcell->tau->L_target;
		l = tcell->tau->lambdaL;
		//len1 = 0;
		//std::cout << " len 3 " << len << std::endl;

		if (l>0) {
			len1 = calc_length(tcell, -index_t);
			//std::cout << " len1 2 " << len1 << std::endl;

			dh += l*( len1*len1 - len*len - 2*L_T*(len1-len) )/(L_T*L_T);
			*lent = len1;
		} // if l>0
		//if (tcell->tau->index==2) std::cout << " dH_len t len1 = " << len1 << " l " << l << std::endl;
	}

	//if (((scell->tau->index==2)||(tcell->tau->index==2))&&(dh!=0)) std::cout << " dH_len dh = " << dh << std::endl;
	return dh;
}


double CPM::dH_chemotaxis(int index_t, pixel* source, cell* scell) const {  // chemotaxis

	int NE = Nx*Ny*Nz,  index_s, i;
	double l = scell->tau->Chi;

	pixel* curr;

	for (i=0; i<NE; ++i) {

		curr = &PIXELS[i];

		if (curr == source) {
			index_s = i;
			break;
		}
	}


	double dh = -l*(PIXELS[index_t].O2 - PIXELS[index_s].O2);
	//std::cout << " dH_chemotaxis index_t " << index_t << " index_s " << index_s << " abs " << abs(index_t-index_s) << std::endl;

	return dh;
}


double CPM::dH_adh(pixel* target, cell* scell) const { // cell adhesion condition

	double dh = 0;
	cell* tcell = target->tag;
	cell* curr;
	int i, j;
	double *Js, *Jt;

	Js = scell->tau->J;
	Jt = tcell->tau->J;
	//std::cout << " scell->tau->index " << scell->tau->index << " [0] " << Js[0] << " [1] " << Js[1]<< " [2] " << Js[2] << " [3] " << Js[3] << " [4] " << Js[4] << " [5] " << Js[5] << std::endl;

	////for (i=0; i<26; ++i) {
	for (i=0; i<6; ++i) {

		curr = target->nbr[i]->tag;
		j = curr->tau->index;
		//std::cout << " j " << j << std::endl;
		//std::cout << " j " << j << " scell " << scell->tau->index << std::endl;
		//std::cout << " j " << j << " scell " << scell->tau->index << " Js[j] " << Js[j] << std::endl;

		if (curr!=scell) dh += Js[j];
		if (curr!=tcell) dh -= Jt[j];
	}

	return dh;
}

double CPM::dH_polar(cell* scell, pixel* source, int index_t) const {   // cell polarization

	double px = 0, py = 1, pz = 0, dh = 0;  // polarization versor

	double l = scell->tau->lambdaP; // amplitude of polarization energy
	if (l == 0) return dh;

	int it, jt, kt, is, js, ks, dx, dy, dz;

	int NE = Nx*Ny*Nz, index_s, i;

	pixel* curr;

	for (i=0; i<NE; ++i) {
		curr = &PIXELS[i];
		if (curr == source) {
			index_s = i;
			break;
		}
	}


	is = index_s/(Ny*Nz);
	js = (index_s - is*Ny*Nz)/Nz; // x, y and z of source pixel
	ks =  index_s - is*Ny*Nz - js*Nz;

	it = index_t/(Ny*Nz);
	jt = (index_t - it*Ny*Nz)/Nz; // x, y and z of target pixel
	kt =  index_t - it*Ny*Nz - jt*Nz;

	dx = it - is;
	dy = jt - js;
	dz = kt - ks;

	if ( (fabs(dx)>1) || (fabs(dy)>1) || (fabs(dz)>1) ) std::cout << " dH_polar ERROR is " << is << " js " << js << " it " << it << " jt " << jt << std::endl;

	dh = -l*(dx*px + dy*py + dz*pz);

	return dh;
}

void CPM::print_tags(std::string filename) const {  // cell ID, sigma

	int i, j, k, n;
	l_cell* curr;
	std::ofstream file;
	file.open(filename.c_str());
/*
	for (k=0; k<Nz; ++k) {
		for (j=0; j<Ny; ++j) {
			for (i=0; i<Nx; ++i) {
				if (PIXELS[i*Ny*Nz+j*Nz+k].tag != nullcell) {
					n = 0;
					curr = c_list;

					while (curr!=NULL) {
						++n;
						if (PIXELS[i*Ny*Nz+j*Nz+k].tag == &(curr->CELL)) curr = NULL;
						if (curr) { curr = curr->next; }
					}
					file << n << " ";
				}
				else {
					file << 0 << " ";
				}
			}
		}
		file << "\n";
	}
*/
	for (i=0; i<Nx; ++i) {
		for (j=0; j<Ny; ++j) {
			for (k=0; k<Nz; ++k) {
				n = 0;
				if (PIXELS[i*Ny*Nz+j*Nz+k].tag != nullcell) {
					curr = c_list;

					while (curr!=NULL) {
						++n; // type
						if (PIXELS[i*Ny*Nz+j*Nz+k].tag == &(curr->CELL)) curr = NULL;
						if (curr) { curr = curr->next; }
					}
				}
				file << i << " " << j << " " << k << " " << n << "\n";
			}
		}
		file << "\n";
	}

	file.close();
}

void CPM::print_types(std::string filename) const {   // cell type, tau

	int i, j, k;
	std::ofstream file;
	file.open(filename.c_str());

	for (i=0; i<Nx; ++i) {
		for (j=0; j<Ny; ++j) {
			for (k=0; k<Nz; ++k) {
				///file << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
				file << i << " " << j << " " << k << " " << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << "\n";
			}
		}
		file << "\n";
	}
	file.close();
}

void CPM::print_O2(std::string filename) const {   // O2 concentration

	int i, j, k;
	std::ofstream file;
	file.open(filename.c_str());

	for (k=0; k<Nz; ++k) {

		for (j=0; j<Ny; ++j) {

			for (i=0; i<Nx; ++i) {

				file << PIXELS[i*Ny*Nz+j*Nz+k].O2 << " ";
			}
		}

		file << "\n";
	}

	file.close();
}


void CPM::print_outlines(std::string filename) const {

	int i, j, k;
	std::ofstream file;
	file.open(filename.c_str());

	cell *curr, *beside, *below;

	for (k=0; k<Nz; ++k) {

	for (j=0; j<Ny; ++j) {

		for (i=0; i<Nx; ++i) {

			curr = PIXELS[i*Ny*Nz+j*Nz+k].tag;

			if (j != Ny-1) {

				beside = PIXELS[i*Ny*Nz+(j+1)*Nz+k].tag;
				if (curr != beside) {

					file << i+1 << " " << j+1 << " " << -1 << " " << 0 << "\n";
				}
			}

			if (i != Nx-1) {

				below = PIXELS[(i+1)*Ny*Nz+j*Nz+k].tag;
				if (curr != below) {

					file << i+1 << " " << j+1 << " " << 0 << " " << -1 << "\n";
				}
			}
		}

	}
	}
	file.close();
}


void CPM::print_outlinesX(std::string filename) const {

	int i, j, k, flag;
	std::ofstream file;
	file.open(filename.c_str());

	cell *curr, *beside, *below, *front;

	for (k=0; k<Nz; ++k) {

		for (j=0; j<Ny; ++j) {

			for (i=0; i<Nx; ++i) {

				curr = PIXELS[i*Ny*Nz+j*Nz+k].tag;
				flag = 0;

				if (k != Nz-1) {
					front = PIXELS[i*Ny*Nz+j*Nz+k+1].tag;
					if (curr != front) flag = 1;
				}

				if (j != Ny-1) {
					beside = PIXELS[i*Ny*Nz+(j+1)*Nz+k].tag;
					if (curr != beside) flag = 1;
				}

				if (i != Nx-1) {
					below = PIXELS[(i+1)*Ny*Nz+j*Nz+k].tag;
					if (curr != below) flag = 1;
				}

				if (flag==0) {
					file << PIXELS[i*Ny*Nz+j*Nz+k].tag->tau->index << " ";
				}
				else {
					file << " 0 ";
				}
			} // cycle i
		file << "\n";

		} // cycle j

	} // cycle k

	file.close();
}


int* CPM::n_cell() { // number of cells
		int * cellsNum = new int[3];
        l_cell* current_node;
        cell* current_cell;

        int type, volume, vol_max = 0;
		int tumor = 0;
		int luminal = 0;
		int basal = 0;

        current_node = c_list;

        while (current_node) {

                current_cell = &(current_node->CELL);

                volume = current_cell->vol;
                type = current_cell->tau->index;

                current_node = current_node->next;
				if (type == 5) tumor++;
				if (type == 2) luminal++;
				if (type == 3) basal++;
		///std::cout << " CPM::n_cell n = " << n << " W0 = " << current_cell->W0 << " reProg = " << current_cell->reProg << std::endl;
        }
	//std::cout << " vol_max = " << vol_max << std::endl;
	cellsNum[0] = tumor;
	cellsNum[1] = luminal;
	cellsNum[2] = basal;
		
	return cellsNum;

}


void CPM::get_list_content () {

        l_cell* cell_list = access_cell_list();
        l_cell* current_node = cell_list;
        cell* current_cell;

        while (current_node) {

                current_cell = &(current_node->CELL);
                std::cout << "VOLUME = " << (current_cell->vol) << std::endl;
                get_content(current_cell);
                current_node = current_node->next;
        }

}


void CPM::get_content(cell* current) {          //corre cada elemento de cell_content e imprime o indice

        std::set<int> myset;
        myset = current->cell_content;
        std::set<int>::iterator it;

       	std::cout << " This cell contains the indexes: ";
        for (it=myset.begin(); it!=myset.end(); ++it)
        	std::cout << ' ' << *it;
       	std::cout << '\n';

}


void CPM::print_nCells(int mcs, std::string filename)  { // count number of cells of different types

        int nBasal = 0; int nStem = 0; int nInter = 0; int nUmbr = 0; int nTumor = 0;
	double x = 0, y = 0, n = 0, ymin = Ny, ymax = 0, ymed = 0;
        l_cell* current_node;
	cell* current_cell;
        current_node = c_list;

        while (current_node) {
           	current_cell = &(current_node->CELL);

           	if(current_cell->tau->index == 3) ++nBasal;
                if(current_cell->tau->index == 4) ++nStem;
                if(current_cell->tau->index == 5) ++nInter;
                if(current_cell->tau->index == 6) ++nUmbr;
                //if(current_cell->tau->index == 8) ++nTumor;
                if(current_cell->tau->index == 8) {
			++nTumor;
		        for (int k=0; k<Nz; ++k) {
			        for (int j=0; j<Ny; ++j) {
        	        		for (int i=0; i<Nx; ++i) {
						if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
							y = y + j;   n++;
							if (j<ymin) ymin = j;
							if (j>ymax) ymax = j;
						}
					}
				}
			}

		}

		if(current_cell->vol < 1) std::cout << " print_nCells vol " << current_cell->vol << " index " << current_cell->tau->index << std::endl;
        	current_node = current_node->next;
    	}

	if (n>0) ymed = y/n;

        std::ofstream file;
        file.open(filename.c_str(), std::ios::app);
        file << mcs << " " << nBasal << " " << nStem << " " << nInter << " " << nUmbr << " " << nTumor << " " << ymax << " " << ymin << " " << ymed << " " << n << "\n";
        file.close();

}


int CPM::ablation()  { // chirurgical removal of papillary tumor

        int ncell_ablation = 0, flag = 0;
        double ymax = 110;

        l_cell* current_node;
        cell* current_cell;
        cell* to_be_removed[4000];
        current_node = c_list;

        while (current_node) {
                current_cell = &(current_node->CELL);
                flag = 0;

                if(current_cell->tau->index == 8) {
                        for (int k=0; k<Nz; ++k) {
	                        for (int j=0; j<Ny; ++j) {
        	                        for (int i=0; i<Nx; ++i) {
                	                        if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
                        	                        if ( (j>ymax) && (flag==0) ) {
                                	                        to_be_removed[ncell_ablation] = current_cell;
                                        	                ncell_ablation++;
                                                	        flag = 1;
	                                                }
                                                }
                                        }
                                }
                        }

                }

                current_node = current_node->next;
        }


	for (int ir=ncell_ablation-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

        return ncell_ablation;

}


int CPM::radical_ablation()  { // chirurgical radical removal of papillary tumor and neighbor cells

        int ncell_ablation = 0, flag = 0;
        double ymax = 110, ymin = 50, xmin = Nx, xmax = 0, zmin = Nz, zmax = 0;

        l_cell* current_node;
        cell* current_cell;
        cell* to_be_removed[4000];
        current_node = c_list;

        while (current_node) { // look for tumor cells above the umbrella layer (visible tumor)
                current_cell = &(current_node->CELL);
                flag = 0;

                if(current_cell->tau->index == 8) {  // tumor cell
                        for (int k=0; k<Nz; ++k) {
	                        for (int j=0; j<Ny; ++j) {
        	                        for (int i=0; i<Nx; ++i) {
                	                        if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
                        	                        if ( (j>ymax) && (flag==0) ) {
                                	                        to_be_removed[ncell_ablation] = current_cell;
                                        	                ncell_ablation++;
                                                	        flag = 1;
	                                                	if ( i>xmax ) xmax = i;
	                                        	        if ( i<xmin ) xmin = i;
	                                                	if ( k>zmax ) zmax = k;
	                                        	        if ( k<zmin ) zmin = k;
	                                                }
                                                }
                                        }
                                }
                        }

                }

                current_node = current_node->next;
        }

	std::cout << " radical_ablation xmin = " << xmin << " xmax " << xmax << " zmin = " << zmin << " zmax " << zmax << std::endl;
	xmin = xmin + 10;   xmax = xmax - 10; // not so radical, border kept
	zmin = zmin + 10;   zmax = zmax - 10; // not so radical, border kept

        current_node = c_list;
        while (current_node) { // look for cells under the tumor (above the basal mabrane)
                current_cell = &(current_node->CELL);
                flag = 0;
		int type = current_cell->tau->index;

                for (int ir=0; ir<ncell_ablation; ++ir) { if(to_be_removed[ir]==current_cell) flag == 1; } //check if already removed

                if( ( (type==3) || (type==4) || (type==5) || (type==6) || (type==8) ) && (flag==0) ) {  // normal, stem and tumor cells
                        for (int k=0; k<Nz; ++k) {
	                        for (int j=0; j<Ny; ++j) {
        	                        for (int i=0; i<Nx; ++i) {
                	                        if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
                        	                        if ( (j>ymin) && (i>=xmin) && (i<=xmax) && (k>=zmin) && (k<=zmax) && (flag==0) ) {
                                	                        to_be_removed[ncell_ablation] = current_cell;
                                        	                ncell_ablation++;
                                                	        flag = 1;
	                                                }
                                                }
                                        }
                                }
                        }

                }

                current_node = current_node->next;
        }


	for (int ir=ncell_ablation-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

        return ncell_ablation;

}


void CPM::print_xCells(int mcs, std::ofstream &file)  { // count number of cells of different types

        int nBasal = 0, nStem = 0, nInter = 0, nUmbr = 0, nTumor = 0, nNull = 0, nStroma = 0, nUrine = 0,
		flag_inv = 0, ileft, iright, jdown, jup, kback, kfront;
        double x = 0, y = 0, n = 0, ymin = Ny, ymax = 0, ymed = 0;
        l_cell* current_node;
	cell* current_cell;
        current_node = c_list;

        while (current_node) {
           	current_cell = &(current_node->CELL);

           	if(current_cell->tau->index == 0) ++nNull;
           	if(current_cell->tau->index == 1) ++nStroma;
           	if(current_cell->tau->index == 3) ++nBasal;
                if(current_cell->tau->index == 4) ++nStem;
                if(current_cell->tau->index == 5) ++nInter;
                if(current_cell->tau->index == 6) ++nUmbr;
                if(current_cell->tau->index == 7) ++nUrine;
                if(current_cell->tau->index == 8) {  // tumor cell
                        ++nTumor;
/*                        for (int k=0; k<Nz; ++k) {
	                        for (int j=0; j<Ny; ++j) {
        	                        for (int i=0; i<Nx; ++i) {
                	                        if(PIXELS[i*Ny*Nz+j*Nz+k].tag == current_cell) {
                        	                        y = y + j;   n++;
                                	                if (j<ymin) ymin = j;
                                        	        if (j>ymax) ymax = j;

							ileft = i - 1;  iright = i + 1; // check for stroma invasion
							jdown = j - 1;  jup = j + 1;    // (contact of tumor cell with stroma)
							kback = k - 1;  kfront = k + 1;    // (contact of tumor cell with stroma)
							if (ileft==-1) ileft = Nx-1; // periodic in x
							if (iright==Nx) iright = 0;
							if (kback==-1) kback = Nz-1; // periodic in z
							if (kfront==Nz) kfront = 0;
							if (jdown==-1) jdown = j;   // non-periodic in y
							if (jup==Ny) jup = j;
							if ( (PIXELS[ileft*Ny*Nz+j*Nz+k].tag->tau->index==1) || (PIXELS[iright*Ny*Nz+j*Nz+k].tag->tau->index==1) ||
							     (PIXELS[i*Ny*Nz+jdown*Nz+k].tag->tau->index==1) || (PIXELS[i*Ny*Nz+jup*Nz+k].tag->tau->index==1) ||
							     (PIXELS[i*Ny*Nz+j*Nz+kback].tag->tau->index==1) || (PIXELS[i*Ny*Nz+j*Nz+kfront].tag->tau->index==1) ) {
								flag_inv = 1; // invasion! Neumann neighborhood
	                                        	}
                                        	}
					}
                                }
                        }
*/
                }

		if(current_cell->vol < 1) std::cout << " print_nCells vol " << current_cell->vol <<
						       " index " << current_cell->tau->index << std::endl;
        	current_node = current_node->next;
    	}

	if (n>0) ymed = y/n;

        //std::ofstream file;
        //file.open(filename.c_str(), std::ios::app);
        //file << mcs << " " << nBasal << " " << nStem << " " << nInter << " " << nUmbr << " " << nTumor << "\n" ;
        file << mcs << " " << nNull << " " << nStroma << " " << nBasal << " " << nStem << " " << nInter << " " << nUmbr <<
		       " " << nUrine << " " << nTumor << "\n";
        //     << ymax << " " << ymin << " " << ymed << " " << n << " " << flag_inv << "\n";
        //file.close();
        //std::cout << mcs << " " << nNull << " " << nStroma << " " << nBasal << " " << nStem << " " << nInter << " " << nUmbr <<
	//		    " " << nUrine << " " << nTumor << std::endl;

}


double CPM::BM_length()  { // check Basal Membrane length (debug)

        l_cell* current_node;
	cell* current_cell;
        current_node = c_list;
	double len;

        while (current_node) {
           	current_cell = &(current_node->CELL);

           	if(current_cell->tau->index == 2) { // basal membrane
			len = current_cell->len; // just look into what value is registered
			std::cout << " BM_vol " << current_cell->voxels.size() << std::endl;
			return len;
                }

		current_node = current_node->next; // next node
    	}
	return len;

}



void CPM::diffusion(int Nx, int Ny, int Nz) {          // diffusion of oxygen

    // Initial Conditions

    // 0-Medium, 1-Fat, 2-Membrane, 3-Basal, 4-Stem, 5-Intermediate, 6-Umbrella, 7-Urine, 8-Tumour, 9-Artery
    //ox_consumption = {0: 0, 1: 0, 2: 0, 3: 100, 4: 125, 5:100, 6: 50, 7: 0, 8:0, 10: 150};

    double tmp, diff, valor_artery, decay_const, diff_const;
    double val_additional;
    int current_cell_type, index, xlef, xrig, ybot, ytop, zbac, zfro;
    double tol = 0.01, err = 1.0;
    valor_artery = 1;
    decay_const = 0.39;
    diff_const = 2500;
    double ox_consumption[] = {0, 0, 0, 0.1, 0.12, 0.10, 0.05, 0, 0.14, 0};

    double* diffusion_matrix = new double [Nx*Ny*Nz];

    while (err > tol) {
	for (int z=0; z<Nz; z++) {
		for (int y=0; y<Ny; y++) {
        		for (int x=0; x<Nx; x++) {
				diffusion_matrix[x*Ny*Nz+y*Nz+z] = PIXELS[x*Ny*Nz+y*Nz+z].O2;
			}
		}
	}
	err = 0;
        for (int z=0; z<Nz; z++) {
	    zbac = z - 1;
	    if(zbac<0) zbac = Nz-1;
	    zfro = z + 1;
	    if(zfro==Nz) zfro = 0;

            for (int y=0; y<Ny; y++) {
	      ybot = y - 1;
	      if(ybot<0) ybot = Ny-1;
	      ytop = y + 1;
	      if(ytop==Ny) ytop = 0;

              for (int x=0; x<Nx; x++) {
		index = x*Ny*Nz + y*Nz + z;
                current_cell_type = PIXELS[x*Ny*Nz+y*Nz+z].tag->tau->index;
                if (current_cell_type == 9) {
                    diffusion_matrix[index] = valor_artery;
                }
		else if (current_cell_type == 7 || current_cell_type == 0) {
                    diffusion_matrix[index] = 0;
                }
		else {
                    val_additional = ox_consumption[current_cell_type];
		    xlef = x - 1;
		    if(xlef<0) xlef = Nx-1;
		    xrig = x + 1;
		    if(xrig==Nx) xrig = 0;

                    tmp = diffusion_matrix[index];
                    diffusion_matrix[index] = ((diffusion_matrix[xlef*Ny*Nz+y*Nz+z] + diffusion_matrix[xrig*Ny*Nz+y*Nz+z] +
						diffusion_matrix[x*Ny*Nz+y*Nz+zbac] + diffusion_matrix[x*Ny*Nz+y*Nz+zfro] +
						diffusion_matrix[x*Ny*Nz+ybot*Nz+z] + diffusion_matrix[x*Ny*Nz+ytop*Nz+z])*diff_const - 
						val_additional*diffusion_matrix[x*Ny*Nz+y*Nz+z])/(decay_const + 6*diff_const);
                    diff = diffusion_matrix[index] - tmp;  // defference new - old
                    //err += diff*diff/(diffusion_matrix[index]*diffusion_matrix[index]);
                    err += diff*diff;
                  }
                }
            }
        }
        err = sqrt(err);
        ///std::cout << " diffusion err = " << err << std::endl;

	for (int z=0; z<Nz; z++) {
  		for (int y=0; y<Ny; y++) {
        		for (int x=0; x<Nx; x++) {
				PIXELS[x*Ny*Nz+y*Nz+z].O2 = diffusion_matrix[x*Ny*Nz+y*Nz+z];
			}
		}
	}

    }

    //return_val = err;

}


double CPM::calc_length (cell* current, int it) const {   // calculate cell length...

	//std::cout << " calc_length it " << it << std::endl;
	//std::vector<int> xx(100000), yy(100000);
	std::vector<int> xx(10000), yy(10000), zz(10000);
	int npixel = 0, flagx0 = 0, flagNx = 0, flagz0 = 0, flagNz = 0, flag, type;
	double lT, len = 0, l;

	type = current->tau->index;
	l = current->tau->lambdaL;
	//std::cout << " calc_length l " << l << " type " << type << std::endl;

	if (l>0) { // only relevant cells

          for (int k=0; k<Nz; ++k) {
          	for (int j=0; j<Ny; ++j) {
                	for (int i=0; i<Nx; ++i) { // less one pixel for target cell
				if((PIXELS[i*Ny*Nz+j*Nz+k].tag == current) && ( (it>=0) || (it!=-(i*Ny*Nz+j*Nz+k)) ) ) {
					xx[npixel] = i;  yy[npixel] = j;  zz[npixel] = k;
					if (i==0) flagx0 = 1; // pixel on border?
					if (i==Nx-1) flagNx = 1; // pixel on the other border?
					if (k==0) flagz0 = 1; // pixel on border?
					if (k==Nz-1) flagNz = 1; // pixel on the other border?
					npixel++; // number of cell pixels
				}
			}
		}
	  }
	  if (it>0) { // new pixel (for source cell)
		xx[npixel] = it/(Ny*Nz);  yy[npixel] = (it-xx[npixel]*Ny*Nz)/Nz;  zz[npixel] = it-xx[npixel]*Ny*Nz-yy[npixel]*Nz;
		npixel++;
	  }
	  //std::cout << " calc_length npixel " << npixel << std::endl;

	  // periodic boundary conditions in x and (take into account cell split)
	  if ( ((type==3)||(type==4)||(type==5)||(type==6)||(type==8)||(type==9)) && (flagx0==1) && (flagNx==1) ) {
        	for (int j=0; j<npixel; ++j) { // cell spreads through the boundary...
			if (xx[j]<Nx/2) xx[j] = xx[j] + Nx;
		}
	  }

	  // periodic boundary conditions in x and (take into account cell split)
	  if ( ((type==3)||(type==4)||(type==5)||(type==6)||(type==8)||(type==9)) && (flagz0==1) && (flagNz==1) ) {
        	for (int j=0; j<npixel; ++j) { // cell spreads through the boundary...
			if (zz[j]<Nz/2) zz[j] = zz[j] + Nz;
		}
	  }

	  ///if ( (type==2) && (flag0==1) && (flagNx==1) ) { // basal membrane (special case)
	  if ( type==2 ) { // basal membrane (special case)
	        for (int k=0; k<Nz; ++k) { // cycle all x points
		        for (int j=0; j<Nx; ++j) { // cycle all x points
				flag = 0;
         			for (int i=0; i<npixel; ++i) {
					if ( (xx[i] == j) || (zz[i] == j) ) flag = 1; // bin on ??
				}
				if (flag==1) len++; // length of all segments (even disconnected ones)
			}
		}
		//std::cout << " calc_length 2 len = " << len << std::endl;
	  }
	  else {
            for (int j=0; j<npixel; ++j) {
                for (int i=j+1; i<npixel; ++i) {
			//std::cout << " calc_length i " << i << " j " << j << std::endl;
			lT = sqrt( (fabs(xx[i]-xx[j])+1)*(fabs(xx[i]-xx[j])+1) + (fabs(yy[i]-yy[j])+1)*(fabs(yy[i]-yy[j])+1) +
			           (fabs(zz[i]-zz[j])+1)*(fabs(zz[i]-zz[j])+1) );
			if (lT>len) len = lT; // highest distance between 2 pixels on cell
	        }
            }
	  } // if type==2

	} // if l>0

	//if (type==2) std::cout << " calc_length len " << len << " type " << type << " npixel " << npixel << std::endl;
	return len;
}


double CPM::dH_area (cell* scell, cell* tcell, int index_t) const {   // calculate cell area...

	double dH = 0;
	int types = scell->tau->index;	int typet = tcell->tau->index;

	if ( (types==2) || (types==6) || (typet==2) || (typet==6) ) {  // only BM or umbrella

		///std::cout << " types " << types << " typet " << typet << std::endl;
		std::vector<int> xx, yy, zz;
		xx.clear();   yy.clear();   zz.clear();

		int npixel, flaga = 0, flagb = 0, area_a, area_b, index;
		double area_t, lambda;

		int jup, jdo, ile, iri, kfr, kba, i, j, k;

/*
        	for (k=0; k<Nz; ++k) {
	        	for (j=0; j<Ny; ++j) {
               			for (i=0; i<Nx; ++i) {
					if ( (types==2) || (types==6) ) {  // only BM or umbrella
						if( PIXELS[i*Ny*Nz+j*Nz+k].tag == scell ) {
							xxs[npixels] = i;  yys[npixels] = j;  zzs[npixels] = k;
							npixels++; // number of source cell pixels
						}
						//if (npixels>9990) std::cout << " npixels " << npixels << std::endl;
					}

					if ( (typet==2) || (typet==6) ) {  // only BM or umbrella
						if( PIXELS[i*Ny*Nz+j*Nz+k].tag == tcell ) {
							xxt[npixelt] = i;  yyt[npixelt] = j;  zzt[npixelt] = k;
							npixelt++; // number of target cell pixels
						}
						//if (npixelt>9990) std::cout << " npixelt " << npixelt << std::endl;
					}
				} // i
			} // j
		} // k
*/
		if ( (types==2) || (types==6) ) {  // only BM or umbrella
			npixel = 0;
			for (int iv=0; iv<scell->voxels.size(); iv++) {
				index = scell->voxels[iv];
				i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;
				xx.push_back(i);  yy.push_back(j);  zz.push_back(k);
				npixel++; // number of source cell pixels
			}

			area_b = 0;  area_a = 0;
			i = index_t/(Ny*Nz);   j = (index_t - i*Ny*Nz)/Nz;   k = index_t - i*Ny*Nz - j*Nz;
			xx.push_back(i);   yy.push_back(j);   zz.push_back(k); // add new pixel (from target cell)
			npixel++; // number of cell pixels
			for (int ip=0; ip<xx.size(); ip++) {
				i = xx[ip];  j = yy[ip];  k = zz[ip];
				index = i*Ny*Nz + j*Nz + k;

				jup = j+1; if (jup==Ny) jup = Ny-1;
				jdo = j-1; if (jdo==-1) jdo = 0;
				iri = i+1; if (iri==Nx) iri = 0;
				ile = i-1; if (ile==-1) ile = Nx-1;  // periodicity
				kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
				kba = k-1; if (kba==-1) kba = Nz-1;
				int index6[6] = {ile*Ny*Nz+j*Nz+k,   i*Ny*Nz+jup*Nz+k,   i*Ny*Nz+j*Nz+kfr,
						   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k, iri*Ny*Nz+j*Nz+k };  // 6 neighbors
				flagb = 0;  flaga = 0;

				for (int ii=0; ii<6; ii++) {
					if( (PIXELS[index6[ii]].tag != scell) && (ip<npixel-1) ) flagb = 1;
					if( (PIXELS[index6[ii]].tag != scell) && (index6[ii]!=index_t) ) flaga = 1;
				}
				if (flagb==1) area_b++;
				if (flaga==1) area_a++;
				/*
					for (int ii=0; ii<6; ii++) {
						if( (PIXELS[index6[ii]].tag != scell) && (ip<npixel-1) ) flagb++;
						if( (PIXELS[index6[ii]].tag != scell) && (index6[ii]!=index_t) ) flaga++;
					}
					area_b = area_b + flagb;
					area_a = area_a + flaga;
				*/
			} // ip

			lambda = scell->tau->lambdaL;	area_t = scell->tau->L_target;
			dH = dH + lambda*(area_a-area_b)*(area_a+area_b-2*area_t)/(area_t*area_t);

		} // (types==2) || (types==6)

		if ( (typet==2) || (typet==6) ) {  // only BM or umbrella
			npixel = 0;
			for (int iv=0; iv<tcell->voxels.size(); iv++) {
				index = tcell->voxels[iv];
				i = index/(Ny*Nz);   j = (index - i*Ny*Nz)/Nz;   k = index - i*Ny*Nz - j*Nz;
				//xx[npixel] = i;  yy[npixel] = j;  zz[npixel] = k;
				xx.push_back(i);  yy.push_back(j);  zz.push_back(k);
				npixel++; // number of source cell pixels
			}

			area_b = 0;  area_a = 0;
			for (int ip=0; ip<xx.size(); ip++) {
				i = xx[ip];  j = yy[ip];  k = zz[ip];
				index = i*Ny*Nz + j*Nz + k;

				jup = j+1; if (jup==Ny) jup = Ny-1;
				jdo = j-1; if (jdo==-1) jdo = 0;
				iri = i+1; if (iri==Nx) iri = 0;
				ile = i-1; if (ile==-1) ile = Nx-1;  // periodicity
				kfr = k+1; if (kfr==Nz) kfr = 0;  // periodicity
				kba = k-1; if (kba==-1) kba = Nz-1;
				int index6[6] = {ile*Ny*Nz+j*Nz+k,   i*Ny*Nz+jup*Nz+k,   i*Ny*Nz+j*Nz+kfr,
						   i*Ny*Nz+j*Nz+kba, i*Ny*Nz+jdo*Nz+k, iri*Ny*Nz+j*Nz+k };  // 6 neighbors
				flagb = 0;  flaga = 0;
				for (int ii=0; ii<6; ii++) {
					if( PIXELS[index6[ii]].tag != tcell ) flagb = 1; // before
					if( ( (PIXELS[index6[ii]].tag != tcell) && (index!=index_t) ) ||
					    ( (index6[ii]==index_t) ) ) flaga = 1;
				}
				if (flagb==1) area_b++;
				if (flaga==1) area_a++;
				/*
					for (int ii=0; ii<6; ii++) {
						if( PIXELS[index6[ii]].tag != tcell ) flagb++; // before
						if( ( (PIXELS[index6[ii]].tag != tcell) && (index!=index_t) ) ||
						    ( (index6[ii]==index_t) ) ) flaga++;
					}
					area_b = area_b + flagb;
					area_a = area_a + flaga;
				*/
			} // ip

			lambda = tcell->tau->lambdaL;	area_t = tcell->tau->L_target;
			dH = dH + lambda*(area_a-area_b)*(area_a+area_b-2*area_t)/(area_t*area_t);

		} // (typet==2) || (typet==6)

		//std::cout << " npixels " << npixels << " npixelt " << npixelt << std::endl;
	} // (types==2) || (types==6) || (typet==2) || (typet==6)

	return dH;

}


int CPM::death_pressure()  { // removal of small cells (due to tumor pressure)

	int ncell_ablation = 0, type;
	double volume, volumeT;

	l_cell* current_node;
	cell* current_cell;
	cell* to_be_removed[4000];
	current_node = c_list;

	while (current_node) {
	current_cell = &(current_node->CELL);
	type = current_cell->tau->index;

		 if (type != 1 && type != 4 && type != 5) { 		// If the cell is stroma or lumen skip to the next cell 
			volume = current_cell->vol;
			volumeT = current_cell->tau->V_target;
			if (volume<0.2*volumeT) {  // too small?
				std::cout << "DEATH - Cell Type: " << type << "; Volume: " << volume << std::endl;
				to_be_removed[ncell_ablation] = current_cell;
				ncell_ablation++;
			}
		}

		current_node = current_node->next;
	}


	for (int ir=ncell_ablation-1; ir>=0; ir--) remove_cell(to_be_removed[ir]); // remove cells

	return ncell_ablation;
}


