// Program to calculate overlap matrix between two wave functions produced from Gaussian09 software package in WFX format.  The program is intended for use in calculating the Kohn-Sham non-adiabatic coupling terms that are input to a subsequent TDDFT NAC code.  The progarm reads .wfx files and outputs a nbnd by nbnd overlap matrix, where nbnd is the total number of occupied+virtual orbitals included in the input_ovlap file.

// Required for input:
// 1) wfx files - two wfx files with the same number of orbitals - must be named 'wfx1' and 'wfx2'
// 2) input_ovlap - file containing the subset of bands for which overlap integral will be calculated on first line, MD time step used to generate configurations on second line (delta_t)

// Output:
// 1) output_coord  - list of nuclear coordinates for first and second configuration
// 2) output_ovlap - overlap matrix of dimensions determined by band_range input
// 3) output_coupling - coupling matrix of dimensions determined by band_range input  - (<psi_a(t)|psi_b(t+delta_t)> - <psi_b(t)|psi_a(t+delta_t)>)/(2*delta_t)
// ** Use TDDFT NAC code to calculate TDDFT NAC (Tavernelli et al 2009) for each pair of adjacent time steps using output_ovlap data.**

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <mpi.h>

using namespace std;

//==============================================
//= PRE-MAIN DELCARATIONS AND PROTOTYPES =======
//==============================================

// FUNCTION PROTOTYPES:
int func_fact(int n); // Compute factorial of input integer
int func_dfact(int n); // Compute double factorial of input integer
int func_bicoeff(int m, int n); // Compute binomial coefficient for two input integers
double func_dist(double x1, double y1, double z1, double x2, double y2, double z2); // Compute 3D distance

// GLOBAL CONSTANTS: Input constants that will be not be changed
const double pi = 3.1415926535;
const double temp_ev = .0257; // room temperature in eV
const double hbar_evfs = 0.6582119; // hbar planck's constant in eV*fs
const double hbar_Js = 1.0546E-34; // hbar planck's constant in J*s
const double mass_au_to_kg = 1.6726E-27; // convert mass from atomic units to kg
const double energy_ev_to_J = 1.619E-19; // convert eV to Joules
const double length_ang_to_m = 1E-10; // convert angstroms to meters
const double ev_to_cm = 8065.73; // convert ev to cm by multiplying by this number

// GLOBAL VARIABLES:
string temp; // temporary file used while looking for string using getline
int elevels_min; // minimum and maximum bands read from band_range
int elevels_max;
int num_elevels; // total number of bands used
double delta_t; // MD time step

//=============================
//= MAIN PROGRAM ------ =======
//=============================

int main()
{

  // Input KS orbital range to be used for overlap matrix
  // File must contain one line - first number is lowest MO, second is highest MO you wish to use, based on Gaussian index
  ifstream input_elevels;
  input_elevels.open("input_ovlap");
  // Input first wave function file (.wfx)
  ifstream input_wfx1;
  input_wfx1.open("wfx1");
  // Input second wave function file (.wfx)
  ifstream input_wfx2;
  input_wfx2.open("wfx2");

  // Output nuclear coordinates in Bohrs
  ofstream out_coord;
  out_coord.open("output_coord");

  // Output overlap matrix
  ofstream out_ovlap_matrix;
  out_ovlap_matrix.open("output_ovlap");

  // Output coupling matrix
  ofstream out_coupling_matrix;
  out_coupling_matrix.open("output_coupling");

  // Testing output
  ofstream out_testing;
  out_testing.open("output_testing");

  // Read min and max bands to use for overlap calculation
  input_elevels >> elevels_min >> elevels_max;
  input_elevels >> delta_t;

  num_elevels = elevels_max - elevels_min + 1;

  cout << "MO Min: " << elevels_min << "  MO Max: " << elevels_max << endl;
  cout << "Time Step: " << delta_t << " fs" << endl;

  // ---READ IN FIRST WAVE FUNCTION DATA---:
  string discard;
  string search_coord = "<Number of Nuclei>"; // search for lines that contain this string, then input levels after this

  // Number of atoms:
  int natoms1;
  int k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  input_wfx1 >> natoms1;
	  k = k+1;
	}
    }

  cout << "----- First Wave Function -----" << endl;
  cout << "Number of atoms: " << natoms1 << endl;

  // Number of total orbitals in wfx file:
  int nmol1;
  search_coord = "<Number of Occupied Molecular Orbitals>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          input_wfx1 >> nmol1;
          k = k+1;
        }
    }

  cout << "Total number of molecular orbitals (MO): " << nmol1 << endl;

  // Nuclear cartesian coordinates:
  double cart_coord1[natoms1][3];
  search_coord = "<Nuclear Cartesian Coordinates>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
	  for(int i =0; i < natoms1; i++)
	    {
	      input_wfx1 >> cart_coord1[i][0] >> cart_coord1[i][1] >> cart_coord1[i][2];
	    }
          k = k+1;
        }
    }

  out_coord << "Wave Function 1 Nuclear Coordinates (Bohrs):" << endl;
  for(int i = 0; i < natoms1; i++)
    {
      out_coord << cart_coord1[i][0] << " " << cart_coord1[i][1] << " " << cart_coord1[i][2] << endl;
    }

  // Number of primitives
  int nprim1;
  search_coord = "<Number of Primitives>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          input_wfx1 >> nprim1;
          k = k+1;
        }
    }

  cout << "Total number of primitive basis functions: " << nprim1 << endl;

  // Primitive centers (atom index for center of each basis functions):
  int prim_cent1[nprim1];
  search_coord = "<Primitive Centers>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
	  for(int i = 0; i < nprim1; i++)
	    {
	      input_wfx1 >> prim_cent1[i];
	    }
          k = k+1;
        }
    }

  // Primitive types (code for orbital symmetry for each primitive function):
  int prim_type1[nprim1];
  search_coord = "<Primitive Types>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          for(int i = 0; i < nprim1; i++)
            {
              input_wfx1 >> prim_type1[i];
            }
          k = k+1;
        }
    }

  // Primitive types: convert code to angular momentum numbers to use for overlap integral later
  // NOTE: ONLY DEFINED THROUGH d-ORBITALS AT THIS POINTS - need to add more for higher angular momenta
  int a1_x[nprim1], a1_y[nprim1], a1_z[nprim1];
  for(int i = 0; i < nprim1; i++)
    {
      if(prim_type1[i] == 1)
	{
	  a1_x[i] = a1_y[i] = a1_z[i] = 0;
	}
      else if(prim_type1[i] == 2)
	{
	  a1_x[i] = 1;
	  a1_y[i] = a1_z[i] = 0;
	}
      else if(prim_type1[i] == 3)
        {
          a1_y[i] = 1;
          a1_x[i] = a1_z[i] = 0;
        }
      else if(prim_type1[i] == 4)
        {
          a1_z[i] = 1;
          a1_x[i] = a1_y[i] = 0;
        }
      else if(prim_type1[i] == 5)
        {
          a1_x[i] = 2;
          a1_y[i] = a1_z[i] = 0;
        }
      else if(prim_type1[i] == 6)
        {
          a1_y[i] = 2;
          a1_x[i] = a1_z[i] = 0;
        }
      else if(prim_type1[i] == 7)
        {
          a1_z[i] = 2;
          a1_x[i] = a1_y[i] = 0;
        }
      else if(prim_type1[i] == 8)
        {
          a1_x[i] = a1_y[i] = 1;
          a1_z[i] = 0;
        }
      else if(prim_type1[i] == 9)
        {
          a1_x[i] = a1_z[i] = 1;
          a1_y[i] = 0;
        }
      else if(prim_type1[i] == 10)
        {
          a1_y[i] = a1_z[i] = 1;
          a1_x[i] = 0;
        }
      else if(prim_type1[i] == 11)
        {
          a1_x[i] = 3;
          a1_y[i] = a1_z[i] = 0;
        }
      else if(prim_type1[i] == 12)
        {
          a1_y[i] = 3;
          a1_x[i] = a1_z[i] = 0;
        }
      else if(prim_type1[i] == 13)
        {
          a1_z[i] = 3;
          a1_y[i] = a1_x[i] = 0;
        }
      else if(prim_type1[i] == 14)
        {
          a1_x[i] = 2;
          a1_y[i] = 1;
          a1_z[i] = 0;
        }
      else if(prim_type1[i] == 15)
        {
          a1_x[i] = 2;
          a1_y[i] = 0;
          a1_z[i] = 1;
        }
      else if(prim_type1[i] == 16)
        {
          a1_x[i] = 0;
          a1_y[i] = 2;
          a1_z[i] = 1;
        }
      else if(prim_type1[i] == 17)
        {
          a1_x[i] = 1;
          a1_y[i] = 2;
          a1_z[i] = 0;
        }
      else if(prim_type1[i] == 18)
        {
          a1_x[i] = 1;
          a1_y[i] = 0;
          a1_z[i] = 2;
        }
      else if(prim_type1[i] == 19)
        {
          a1_x[i] = 0;
          a1_y[i] = 1;
          a1_z[i] = 2;
        }
      else if(prim_type1[i] == 20)
        {
          a1_x[i] = 1;
          a1_y[i] = 1;
          a1_z[i] = 1;
        }
    }

  // Primitive exponents:
  double prim_exp1[nprim1];
  search_coord = "<Primitive Exponents>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          for(int i = 0; i < nprim1; i++)
            {
              input_wfx1 >> prim_exp1[i];
            }
          k = k+1;
        }
    }

  // Molecular Orbital Coefficients (contraction coefficients):
  double mo_coeff1[nmol1][nprim1];
  search_coord = "</MO Number>";
  k = 0;
  while(k < nmol1)
    {
      getline(input_wfx1,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          for(int i = 0; i < nprim1; i++)
            {
              input_wfx1 >> mo_coeff1[k][i];
            }
          k = k+1;
        }
    }
  // ---END READING OF FIRST WAVE FUNCTION---

  // ---READ IN SECOND WAVE FUNCTION DATA---:
  search_coord = "<Number of Nuclei>"; // search for lines that contain this string, then input levels after this

  // Number of atoms:
  int natoms2;
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          input_wfx2 >> natoms2;
          k = k+1;
        }
    }

  cout << "----- Second Wave Function -----" << endl;
  cout << "Number of atoms: " << natoms2 << endl;

  // Number of total orbitals in wfx file:
  int nmol2;
  search_coord = "<Number of Occupied Molecular Orbitals>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          input_wfx2 >> nmol2;
          k = k+1;
        }
    }

  cout << "Total number of molecular orbitals (MO): " << nmol2 << endl;

  // Nuclear cartesian coordinates:
  double cart_coord2[natoms2][3];
  search_coord = "<Nuclear Cartesian Coordinates>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
        {
          for(int i =0; i < natoms2; i++)
            {
              input_wfx2 >> cart_coord2[i][0] >> cart_coord2[i][1] >> cart_coord2[i][2];
            }
          k = k+1;
        }
    }
     
  out_coord << " " << endl;
  out_coord << "Wave Function 2 Nuclear Coordinates (Bohrs):" << endl;
  for(int i = 0; i < natoms2; i++)
    {
      out_coord << cart_coord2[i][0] << " " << cart_coord2[i][1] << " " << cart_coord2[i][2] << endl;
    }

  // Number of primitives
  int nprim2;
  search_coord = "<Number of Primitives>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  input_wfx2 >> nprim2;
	  k = k+1;
	}
    }

  cout << "Total number of primitive basis functions: " << nprim2 << endl;

  // Primitive centers (atom index for center of each basis functions):
  int prim_cent2[nprim2];
  search_coord = "<Primitive Centers>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  for(int i = 0; i < nprim2; i++)
	    {
	      input_wfx2 >> prim_cent2[i];
	    }
	  k = k+1;
	}
    }

  // Primitive types (code for orbital symmetry for each primitive function):
  int prim_type2[nprim2];
  search_coord = "<Primitive Types>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  for(int i = 0; i < nprim2; i++)
	    {
	      input_wfx2 >> prim_type2[i];
	    }
	  k = k+1;
	}
    }

  // Primitive types: convert code to angular momentum numbers to use for overlap integral later
  // NOTE: ONLY DEFINED THROUGH d-ORBITALS AT THIS POINTS - need to add more for higher angular momenta
  int a2_x[nprim2], a2_y[nprim2], a2_z[nprim2];
  for(int i = 0; i < nprim2; i++)
    {
      if(prim_type2[i] == 1)
        {
          a2_x[i] = a2_y[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 2)
        {
          a2_x[i] = 1;
          a2_y[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 3)
        {
          a2_y[i] = 1;
          a2_x[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 4)
        {
          a2_z[i] = 1;
          a2_x[i] = a2_y[i] = 0;
        }
      else if(prim_type2[i] == 5)
        {
          a2_x[i] = 2;
          a2_y[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 6)
        {
          a2_y[i] = 2;
          a2_x[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 7)
        {
          a2_z[i] = 2;
          a2_x[i] = a2_y[i] = 0;
        }
      else if(prim_type2[i] == 8)
        {
          a2_x[i] = a2_y[i] = 1;
          a2_z[i] = 0;
        }
      else if(prim_type2[i] == 9)
        {
          a2_x[i] = a2_z[i] = 1;
          a2_y[i] = 0;
        }
      else if(prim_type2[i] == 10)
        {
          a2_y[i] = a2_z[i] = 1;
          a2_x[i] = 0;
        }
      else if(prim_type2[i] == 11)
        {
          a2_x[i] = 3;
          a2_y[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 12)
        {
          a2_y[i] = 3;
          a2_x[i] = a2_z[i] = 0;
        }
      else if(prim_type2[i] == 13)
        {
          a2_z[i] = 3;
          a2_y[i] = a2_x[i] = 0;
        }
      else if(prim_type2[i] == 14)
        {
          a2_x[i] = 2;
          a2_y[i] = 1;
          a2_z[i] = 0;
        }
      else if(prim_type2[i] == 15)
        {
          a2_x[i] = 2;
          a2_y[i] = 0;
          a2_z[i] = 1;
        }
      else if(prim_type2[i] == 16)
        {
          a2_x[i] = 0;
          a2_y[i] = 2;
          a2_z[i] = 1;
        }
      else if(prim_type2[i] == 17)
        {
          a2_x[i] = 1;
          a2_y[i] = 2;
          a2_z[i] = 0;
        }
      else if(prim_type2[i] == 18)
        {
          a2_x[i] = 1;
          a2_y[i] = 0;
          a2_z[i] = 2;
        }
      else if(prim_type2[i] == 19)
        {
          a2_x[i] = 0;
          a2_y[i] = 1;
          a2_z[i] = 2;
        }
      else if(prim_type2[i] == 20)
        {
          a2_x[i] = 1;
          a2_y[i] = 1;
          a2_z[i] = 1;
        }
    }

  // Primitive exponents:
  double prim_exp2[nprim2];
  search_coord = "<Primitive Exponents>";
  k = 0;
  while(k < 1)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  for(int i = 0; i < nprim2; i++)
	    {
	      input_wfx2 >> prim_exp2[i];
	    }
	  k = k+1;
	}
    }

  // Molecular Orbital Coefficients (contraction coefficients):
  double mo_coeff2[nmol2][nprim2];
  search_coord = "</MO Number>";
  k = 0;
  while(k < nmol2)
    {
      getline(input_wfx2,temp); // read each line in wfx file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  for(int i = 0; i < nprim2; i++)
	    {
	      input_wfx2 >> mo_coeff2[k][i];
	    }
	  k = k+1;
	}
    }

  // ---END READING OF SECOND WAVE FUNCTION---

  // -----------------------------------------
  // ----COUPLING MATRIX CALCULATION----------
  // -----------------------------------------

  // Calculate (<psi_a(t)|psi_b(t+delta_t)> - <psi_a(t+delta_t)|psi_b(t)>)/(2*delta_t)
  // delta_t is time step used in MD simulation

  // ---CALCULATE OVERLAP INTEGRAL: <psi_a(t)|psi_b(t+delta_t)>---
  // ---Parallelized version over outer loop of bands---

  double ovlap_mat[num_elevels][num_elevels];  // overlap matrix calculated in each parallel process
  double ovlap_mat_tot[num_elevels][num_elevels]; // final overlap matrix summed over parallel processes
  int final_size = num_elevels*num_elevels;
  double coupling_mat[num_elevels][num_elevels];

  // Intermediate variables to calculate overlap between each MO that are overwritten each time through loop
  double ov_gauss; // Gaussian component of overlap integral
  double dist; // Distance between nuclei in Gaussian component of overlap
  double ov_x;  // x-component of overlap integral
  double ov_y;  // y-component of overlap integral
  double ov_z;  // z-component of overlap integral
  double ov_x_sum;
  double ov_y_sum;
  double ov_z_sum;
  double px;   // cartesian componenets of subsection of overlap integral
  double py;
  double pz;

  for (int i = 0; i < num_elevels; i++)
    {
      for (int j = 0; j < num_elevels; j++)
        {
	  ovlap_mat[i][j] = 0.0;
	  ovlap_mat_tot[i][j] = 0.0;
	}
    }

  //Parallelization initialization
  int rank, size, start, end;

  MPI_Init(NULL, NULL);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  start = rank*(num_elevels/size);
  if(rank == (size - 1))
    {
      end = num_elevels;
    }
  else
    {
      end = start + (num_elevels/size);
    }
  
  //Overlap Integral - Analytical Calculation
  for (int i = start; i < end; i++)
    {
      for (int j = 0; j < num_elevels; j++)
	{
	  for (int l = 0; l < nprim1; l++)
	    {
	      for (int m = 0; m < nprim2; m++)
		{
		  ov_gauss = 0.0;
		  dist = 0.0;
		  ov_x = 0.0;
		  ov_y = 0.0;
		  ov_z = 0.0;

		  dist = func_dist(cart_coord1[prim_cent1[l]-1][0],cart_coord1[prim_cent1[l]-1][1],cart_coord1[prim_cent1[l]-1][2],cart_coord2[prim_cent2[m]-1][0],cart_coord2[prim_cent2[m]-1][1],cart_coord2[prim_cent2[m]-1][2]);
		  ov_gauss = exp(-((prim_exp1[l]*prim_exp2[m])/(prim_exp1[l] + prim_exp2[m]))*dist*dist);
		  	  
		  // x-component sum
		  ov_x_sum = 0.0;
		  px = 0.0;
		  px = (prim_exp1[l]*cart_coord1[prim_cent1[l]-1][0] + prim_exp2[m]*cart_coord2[prim_cent2[m]-1][0])/(prim_exp1[l] + prim_exp2[m]);
		  for (int p = 0; p <= a1_x[l]; p++)
		    {
		      for (int q = 0; q <= a2_x[m]; q++)
			{
			  if ((p+q) % 2 == 0)
			    {
			      ov_x_sum = ov_x_sum + func_bicoeff(a1_x[l],p)*func_bicoeff(a2_x[m],q)*(func_dfact(p+q-1)/(pow(2*(prim_exp1[l]+prim_exp2[m]),0.5*(p+q))))*(pow((px - cart_coord1[prim_cent1[l]-1][0]),(a1_x[l]-p)))*(pow((px - cart_coord2[prim_cent2[m]-1][0]),(a2_x[m]-q)));
			    }
			  else
			    {
			      ov_x_sum = ov_x_sum + 0.0;
			    }
			}
		    }      
		  ov_x = ov_x_sum;

		  // y-component sum
		  ov_y_sum = 0.0;
		  py = 0.0;
		  py = (prim_exp1[l]*cart_coord1[prim_cent1[l]-1][1] + prim_exp2[m]*cart_coord2[prim_cent2[m]-1][1])/(prim_exp1[l] + prim_exp2[m]);
		  for (int p =0; p <= a1_y[l]; p++)
		    {
		      for (int q = 0; q <= a2_y[m]; q++)
			{
			  if ((p+q) % 2 == 0)
			    {
			      ov_y_sum = ov_y_sum + func_bicoeff(a1_y[l],p)*func_bicoeff(a2_y[m],q)*(func_dfact(p+q-1)/(pow(2*(prim_exp1[l]+prim_exp2[m]),0.5*(p+q))))*(pow((py - cart_coord1[prim_cent1[l]-1][1]),(a1_y[l]-p)))*(pow((py - cart_coord2[prim_cent2[m]-1][1]),(a2_y[m]-q)));
			    }
			  else
			    {
			      ov_y_sum = ov_y_sum + 0.0;
			    }
			}
		    }
		  ov_y = ov_y_sum;

		  // z-component sum
		  ov_z_sum = 0.0;
		  pz = 0.0;
		  pz = (prim_exp1[l]*cart_coord1[prim_cent1[l]-1][2] + prim_exp2[m]*cart_coord2[prim_cent2[m]-1][2])/(prim_exp1[l] + prim_exp2[m]);
		  for (int p = 0; p <= a1_z[l]; p++)
		    {
		      for (int q = 0; q <= a2_z[m]; q++)
			{
			  if ((p+q) % 2 == 0)
			    {
			      ov_z_sum = ov_z_sum + func_bicoeff(a1_z[l],p)*func_bicoeff(a2_z[m],q)*(func_dfact(p+q-1)/(pow(2*(prim_exp1[l]+prim_exp2[m]),0.5*(p+q))))*(pow((pz - cart_coord1[prim_cent1[l]-1][2]),(a1_z[l]-p)))*(pow((pz - cart_coord2[prim_cent2[m]-1][2]),(a2_z[m]-q)));
			    }
			  else
			    {
			      ov_z_sum = ov_z_sum + 0.0;
			    }
			}
		    }
		  ov_z = ov_z_sum;
		  
		  // Sum overlap matrix as product of gaussian, cartesian components, and MO coefficients:
		  ovlap_mat[i][j] = ovlap_mat[i][j] + pow(sqrt(pi/(prim_exp1[l]+prim_exp2[m])),3)*ov_gauss*ov_x*ov_y*ov_z*mo_coeff1[i+elevels_min-1][l]*mo_coeff2[j+elevels_min-1][m];
		}
	    }
	}
    }

  // Combine all ovlap_mat into final ovlap_mat_tot containing overlap integrals of all band combinations
  MPI_Reduce(&ovlap_mat, &ovlap_mat_tot, final_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // Output overlap matrix

  if(rank == 0)
    {
      for (int i = 0; i < num_elevels; i++)
	{
	  for (int j = 0; j < num_elevels; j++)
	    {
	      out_ovlap_matrix << scientific << setprecision(6) << ovlap_mat_tot[i][j] << " ";
	    }
	  out_ovlap_matrix << endl;
	}
    
      // CALCULATE COUPLING MATRIX - since wave functions are real, don't need to worry about complex conjugate
      // Multiplied by hbar to give units of energy - wfx wave function is normalized, so hbar*(1/delta_t) == energy

      for (int i = 0; i < num_elevels; i++)
	{
	  for (int j = 0; j < num_elevels; j++)
	    {
	      coupling_mat[i][j] = hbar_evfs*(ovlap_mat_tot[i][j] - ovlap_mat_tot[j][i])/(2*delta_t);
	    }
	}

      // Output coupling matrix
      for (int i = 0; i < num_elevels; i++)
	{
	  for (int j = 0; j < num_elevels; j++)
	    {
	      out_coupling_matrix << scientific << setprecision(6) << coupling_mat[i][j] << " ";
	    }
	  out_coupling_matrix << endl;
	}
    }
  // End MPI section of code:
  MPI_Finalize();

  return 0;
}

//=============================
// FUNCTIONS:
//=============================

// Factorial function
int func_fact(int n)
{
  int func_fact = 1;
  int i;
  for (i = 1; i <=n; i++)
    {
      func_fact*=i;
    }
  return func_fact;
}

// Double factorial function
int func_dfact(int n)
{
  int func_dfact = 1;
  int i;
  for (i = n; i >= 1 ; i -=2)
    {
      func_dfact*=i;
    }
  return func_dfact;
}

// Binomial coefficient function
int func_bicoeff(int m,int n)
{
  int func_factM = 1;
  int i;
  for (i = 1; i <=m; i++)
    {
      func_factM*=i;
    }
  int func_factN = 1;
  int j;
  for (j = 1; j <=n; j++)
    {
      func_factN*=j;
    }
  int func_factM_N = 1;
  int k;
  for (k = 1; k <=(m-n); k++)
    {
      func_factM_N*=k;
    }
  int func_bicoeff;
  func_bicoeff = func_factM/(func_factN*func_factM_N);

  return func_bicoeff;
}

// Distance function
double func_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
  double dist;
  dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
  
  return dist;
}
