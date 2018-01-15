
// The primary function of this code is to calculate the Huang Rhys factors for each pair of excited states.  These factors are then used as input for the Franck-Condon-weighted density of states code (FCDOS) that calculates the lineshape function giving the phonon density states that enters the equation for relaxation rates.

// Program takes as input the concatenated Gaussian outputs of scf calculations for all excited states (and ground state) relating to states involved in electronic transitions and calculates Huang-Rhys factors related to each mode for each transition.  The code first calculates the coordination change (DeltaR) corresponding to the difference in each excited state's equilibrium geometries.  The code then calculates the inner product with the normal mode eigenvectors of the ground state system find the Huang-Rhys factors corresponding to each mode and state combination.  The code outputs the Huang-Rhys factor for each normal mode for each band combination - convert to unitless value by multiplying by (mw/hbar) (output_HRF_matrix).  The Huang-Rhys factor output is used as input to calculate the Franck-Condon-weighted density of states using a separate Matlab code. 

// Required for input:
// 1) energy_pop: header with min and max excited state index, followed by list of excitation energies from Gaussiasn TDDFT output
// 2) gauss_out_coord: Concatenated Gaussian output files from scf calculation that includes list of atomic coordinates
// 3) gauss_out_freq: outptu from gaussian frequency calculation using ground state
// 4) modes_cm: single column of modes in cm-1 listed in same order as given in gauss_out_freq
// 5) input_fcf: input file giving the number of atoms (natoms), number of atom types (ntypes) and the atomic symbol, number, mass, and number of atoms for each type.  Number of total atoms and number of types must be given in first two lines, followed by a line for each atom.  For example:
// natoms=38
// ntypes=3
// C 6 12.01 20
// H 1 1.00  14
// N 7 14.01 4

// Output:
// 1) output_HRF_info: organized by band combination, then lists the energy, projected displacement, HRF, and reorganization energy for each mode pertaining to this transition.  This file is helpful for plotting the HRF vs mode energy to identify dominant modes
// 2) output_HR_matrix: matrix output used as input for FCDOS calculation in separate code.  Each column is one band combination (beginning with 0-1, 0-2, ..., 1-2, 1-3,... and so on). Each row is the Huang Rhys factor for one mode.
// 2) output_deltaR: lists the displacement between each excited state combination
// 3) output_phonon: lists the phonon energies and their reduced mass (which enter into HRF calculation)

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

//==============================================
//= PRE-MAIN DELCARATIONS AND PROTOTYPES =======
//==============================================

// FUNCTION PROTOTYPES:
int func_fact(int n); // Compute factorial of input integer

// Convergence parameters:

// GLOBAL CONSTANTS: Input constants that will be not be changed
const double pi = 3.14159;
const double temp_ev = .0257; // room temperature in eV
const double hbar_evfs = 0.6582119; // hbar planck's constant in eV*fs
const double hbar_Js = 1.0546E-34; // hbar planck's constant in J*s
const double mass_au_to_kg = 1.6726E-27; // convert mass from atomic units to kg
const double energy_ev_to_J = 1.619E-19; // convert eV to Joules
const double length_ang_to_m = 1E-10; // convert angstroms to meters
const double ev_to_cm = 8065.73; // convert ev to cm by multiplying by this number

//=============================
//= MAIN PROGRAM ------ =======
//=============================

int main()
{
  //----INPUT FILES----//
  // Input data of energy levels taken from Gaussiasn TDDFT energy levels
  ifstream input_elevels;
  input_elevels.open("energy_pop");
  // Input data of coordinate positions from Gaussian for ground state and all excited states from lowest to highest energy
  ifstream coord_file;
  coord_file.open("gauss_out_coord");
  // Input data for phonon modes to use for inner product
  ifstream ph_mod_file; // OUTCAR file containing normal modes
  ph_mod_file.open("gauss_out_freq");
  // Input data from modes_cm file to get gaussian phonon energies - in cm-1 units
  ifstream modes_cm_file;
  modes_cm_file.open("modes_cm");
  // Input data used to calculate reduced mass
  ifstream input_fcf_file;
  input_fcf_file.open("input_fcf");

  //----OUTPUT FILES----//
  // Output transition energies
  ofstream out_trans_e_file;
  out_trans_e_file.open("output_trans_energies");

  // Output data of all coordinates and forces input
  ofstream out_phonon_file;
  out_phonon_file.open("output_phonon");
  // Output displacement shift in potential energy surface minima between states and HR factor for all band combos
  ofstream out_disp_norm;
  out_disp_norm.open("output_HRF_info");
  // Output HR matrix (phonon mode by transition energy) for input to matlab code
  ofstream out_HR_matlab;
  out_HR_matlab.open("output_HR_matrix");
  // Output deltaR for each band combination
  ofstream out_deltaR;
  out_deltaR.open("output_deltaR");
  // Output unweighted DOS times BE+1 for each transition energy
  ofstream out_dos_total;
  out_dos_total.open("output_dos_unweighted");
  // Output for code testing
  ofstream out_testing;
  out_testing.open("output_testing");
  
  // -------------------- //
  // ---- READ INPUT ---- //
  // -------------------- //

  // Read first two lines of input_fcf file for number of atoms and number of atom types:
  string input_line1;
  string input_line2;
  int natoms; // number of total atoms in system
  int ntypes; // number of atom types in system
  int input_counter = 0;

  getline(input_fcf_file,input_line1); // Read in first two lines of input for number of atoms and atom types
  getline(input_fcf_file,input_line2);
  // Check if first input line has total number of atoms or atom types and assign to variables:
  if(input_line1.find("natoms") != string::npos)
    {
      size_t eq_pos1 = input_line1.find("=");
      string str_natoms = input_line1.substr(eq_pos1+1);
      stringstream(str_natoms) >> natoms;
      input_counter = input_counter + 1;
    } else if(input_line1.find("ntypes") != string::npos)
    {
      size_t eq_pos1 = input_line1.find("=");
      string str_ntypes = input_line1.substr(eq_pos1+1);
      stringstream(str_ntypes) >> ntypes;
      input_counter = input_counter + 1;
    }
  // Check if second input line has total number of atoms or atom types and assign to variables:  
  if(input_line2.find("natoms") != string::npos)
    {
      size_t eq_pos2 = input_line2.find("=");
      string str_natoms = input_line2.substr(eq_pos2+1);
      stringstream(str_natoms) >> natoms;
      input_counter = input_counter + 1;
      cout << natoms << endl;
    } else if(input_line2.find("ntypes") != string::npos)
    {
      size_t eq_pos2 = input_line2.find("=");
      string str_ntypes = input_line2.substr(eq_pos2+1);
      stringstream(str_ntypes) >> ntypes;
      input_counter = input_counter + 1;
      cout << ntypes << endl;
    }
  // Quit program if both natoms and ntypes not assigned:
  if(input_counter != 2)
    {
      cout << "Both total atoms (natoms) and atom types (ntypes) have not been provided in first two lines of input!" << endl;
      exit(1);
    }

  cout << "Number of atoms: " << natoms << endl;
  cout << "Number of atom types: " << ntypes << endl;

  // Read in information about each atom (symbol, number, and mass) and calculate total mass
  string symb_atom[ntypes]; // array of atomic symbols
  int numb_atom[ntypes]; // array of atomic numbers
  double mass_atom[ntypes]; // array of atomic masses
  int natom_type[ntypes]; // array of number of atoms of each type
  double mass_tot = 0.0; // total mass of system

  for(int i = 0; i < ntypes; i++)
    {
      input_fcf_file >> symb_atom[i] >> numb_atom[i] >> mass_atom[i] >> natom_type[i];
      cout << "Symbol: " << symb_atom[i] << " Atomic Number: " << numb_atom[i] << " Mass: " << mass_atom[i] << " Number: " << natom_type[i] << endl;
      mass_tot = mass_tot + mass_atom[i]*natom_type[i];
    }

  cout << "Total Mass: " << mass_tot << endl;

  // Read in energy levels from energy_pop:
  int elevels_min; // minimum and maximum bands read from energy_pop
  int elevels_max;
  int num_elevels; // total number of bands used

  input_elevels >> elevels_min >> elevels_max;
  num_elevels = elevels_max - elevels_min + 1;
  int band_comb = (num_elevels*(num_elevels-1))/2; // number of unique transition energies
  double elevels[num_elevels][3]; // 3d array to input band number, energy, and occupation
  for(int i = 0; i < num_elevels; i++)
    {
      input_elevels >> elevels[i][0] >> elevels[i][1] >> elevels[i][2];
    }

  // Create matrix of transition energies of size num_elevels by num_elevels
  double band_ediff[num_elevels][num_elevels];
  double band_ediff_onerow[band_comb];
  int q = 0;
  for (int m = 0; m < num_elevels; m++)
    {
      for (int n = m+1; n < num_elevels; n++)
	{
	  band_ediff_onerow[q] = elevels[n][1] - elevels[m][1];
	  band_ediff[m][n] = elevels[n][1] - elevels[m][1];
	  band_ediff[n][m] = -band_ediff[m][n];
	  q = q +1;
	}
      band_ediff[m][m] = 0.0;
    }
  
  for (int m =0; m < num_elevels; m++)
    {
      for (int n = 0; n < num_elevels; n++)
	{
	  out_trans_e_file << band_ediff[m][n] << " ";
	}
      out_trans_e_file << endl;
    }

  // Read in phonon energies and conver to meV:
  int nmodes = 3*natoms-6; // number of modes of system, gaussian removes 3 translation and 3 rotational
  double modes_mev[nmodes]; // array of all normal mode energies in mev
  for (int k = 0; k < nmodes; k++)
    {
      modes_cm_file >> modes_mev[k];
      modes_mev[k] = modes_mev[k]/ev_to_cm*1000; // convert gaussian cm-1 input to mev
    }

  // Read in cartesian positions of atoms for lower and higher energy states and store associated atomic mass
  string temp; // temporary file used while looking for string using getline
  double coord_cart[natoms][4][num_elevels];
  int coord_info[natoms][3][num_elevels];
  string discard;
  string search_coord = "Coordinates (Angstroms)"; // search for lines that contain this string, then input levels after this

  int k = 0;
  while(k < (num_elevels))
    {
      getline(coord_file,temp); // read each line in MD file
      if(temp.find(search_coord) != string::npos) // if line contains the search term, begin loop
	{
	  coord_file >> discard >> discard >> discard >> discard >> discard >> discard; // discard next two lines
	  coord_file >> discard;
	  for(int i = 0; i < natoms; i++)
	    {
	      // input cartesian positions given excited state k
	      coord_file >> coord_info[i][0][k] >> coord_info[i][1][k] >> coord_info[i][2][k] >> coord_cart[i][0][k] >> coord_cart[i][1][k] >> coord_cart[i][2][k];
	      // Input atomic mass corresponding to each cartesian coordinate into coord_cart[i][3][k]:
	      for(int j = 0; j < ntypes; j++)
		{
		  if(coord_info[i][1][k] == numb_atom[j])
		    {
		      coord_cart[i][3][k] = mass_atom[j];
		    }
		}
	    }
	  k = k + 1;
	}
    }

  // ----------------- //
  // ----END INPUT---- //
  // ----------------- //

  // Calculate center of mass along each coordinate axis for each excited state
  double COM_x[num_elevels];
  double COM_y[num_elevels];
  double COM_z[num_elevels];

  for(int k = 0; k < num_elevels; k++)
    {
      COM_x[k] = 0.0;
      COM_y[k] = 0.0;
      COM_z[k] = 0.0;
      for(int i = 0; i < natoms; i++)
	{
	  COM_x[k] = COM_x[k] + coord_cart[i][3][k]*coord_cart[i][0][k];
	  COM_y[k] = COM_y[k] + coord_cart[i][3][k]*coord_cart[i][1][k];
	  COM_z[k] = COM_z[k] + coord_cart[i][3][k]*coord_cart[i][2][k];
	}
      COM_x[k]= COM_x[k]/mass_tot;
      COM_y[k]= COM_y[k]/mass_tot;
      COM_z[k]= COM_z[k]/mass_tot;
    }

  // Redefine cartesian coordinates as shifted by center of mass

  for(int k = 0; k < num_elevels; k++)
    {
      for(int i =0; i < natoms; i++)
	{
	  coord_cart[i][0][k] = coord_cart[i][0][k] - COM_x[k];
	  coord_cart[i][1][k] = coord_cart[i][1][k] - COM_y[k];
	  coord_cart[i][2][k] = coord_cart[i][2][k] - COM_z[k];
	}
    }

  // Input eigenvector for each normal mode from OUTCAR_ph

  double cart_modes[natoms][3][nmodes]; // cartesian coordinates of each mode
  double eig_modes[natoms][3][nmodes]; // eigenvector of each mode
  string discard1, discard2, discard3, discard4, discard5, discard6;

  int l = 0;
  string search_ph = "Atom  AN"; // search for lines that contain this string, then input levels after this
  while (l < nmodes)
    {
      getline(ph_mod_file,temp); // read each line in OUTCAR_ph
      if(temp.find(search_ph) != string::npos) // if line contains the search term, begin loop
	{
	  for(int i = 0; i < natoms; i++)
	    {
	      // input eigenvectors of each mode - gaussian output has three modes per section
	      ph_mod_file >> discard1 >> discard2 >> eig_modes[i][0][l] >> eig_modes[i][1][l] >> eig_modes[i][2][l] >> eig_modes[i][0][l+1] >> eig_modes[i][1][l+1] >> eig_modes[i][2][l+1] >> eig_modes[i][0][l+2] >> eig_modes[i][1][l+2] >> eig_modes[i][2][l+2];
	    }
	  l = l+3;
	}
    }

  // Renormalize eigenmodes so that the dot product preserves the length of coordinate difference between excited states
  // We must multiply each eigenvector component from Gaussian output by sqrt of mass and renormalize to do this
  // Also multiply each cartesian displacement by square root of mass of given atom, then don't multiply by reduced mass for HRF calc

  double eig_modes_renorm[natoms][3][nmodes];
  double coord_cart_sqrtm[natoms][3][num_elevels];

  if (ntypes == 0)
    {
      cout << "Need number of atom types to be greater than or equal to 1 for program to run correctly!" << endl;
      return 1;
    }
  else{
      for (int k = 0; k < nmodes; k++)
	{
	  for (int i = 0; i < natoms; i++)
	    {
	      eig_modes_renorm[i][0][k] = eig_modes[i][0][k]*(sqrt(coord_cart[i][3][0]));
              eig_modes_renorm[i][1][k] = eig_modes[i][1][k]*(sqrt(coord_cart[i][3][0]));
              eig_modes_renorm[i][2][k] = eig_modes[i][2][k]*(sqrt(coord_cart[i][3][0]));
	    }
	}
      for (int k = 0; k < num_elevels; k++)
        {
          for (int i = 0; i < natoms; i++)
            {
              coord_cart_sqrtm[i][0][k] = coord_cart[i][0][k]*(sqrt(coord_cart[i][3][k]));
              coord_cart_sqrtm[i][1][k] = coord_cart[i][1][k]*(sqrt(coord_cart[i][3][k]));
              coord_cart_sqrtm[i][2][k] = coord_cart[i][2][k]*(sqrt(coord_cart[i][3][k]));
            }
        }
    }

  // Calculate norm of each normal mode using renormalized values
  double norm_modes[nmodes]; // norm of each normal mode
  double mode_sum;

  for (int k = 0; k < nmodes; k++)
    {
      mode_sum = 0.0;
      for (int i = 0; i < natoms; i++)
        {
          mode_sum = mode_sum + eig_modes_renorm[i][0][k]*eig_modes_renorm[i][0][k] + eig_modes_renorm[i][1][k]*eig_modes_renorm[i][1][k] + eig_modes_renorm[i][2][k]*eig_modes_renorm[i][2][k];
        }
      norm_modes[k] = sqrt(mode_sum);
    }

  // Calculate reduced mass of each normal mode: mu = (Sum_i[(c_i)^2/m_i])^-1
  double red_mass[nmodes]; // reduced mass of normal mode k
  double red_mass_sum[nmodes]; // sum of product of coefficients and inverse mass that is inverted to calculate reduced mass

  // Note separate loop for increasing number of n_type - only coded for up to n_type = 4!
  if (ntypes == 0)
    {
      cout << "Need number of atom types to be greater than or equal to 1 for program to run correctly!" << endl;
      return 1;
    }
  else{
      for (int k = 0; k < nmodes; k++)
	{
	  red_mass[k] = 0.0;
	  red_mass_sum[k] = 0.0;
	  for (int i = 0; i < natoms; i++)
	    {
	      for (int j = 0; j < 3; j++)
		{
		  red_mass_sum[k] = red_mass_sum[k] + (eig_modes[i][j][k]*eig_modes[i][j][k])*coord_cart[i][3][0];
		}
	    }
	  red_mass[k] = red_mass_sum[k];
	}
    }

  // Output phonon information:
  out_phonon_file << "Mode Number:      " << "Mode energy:      " << "Mode Reduced Mass:    " << endl;
  out_phonon_file << "-------------------------------------------------------------" << endl;
  for (int k = 0; k < nmodes; k++)
    {
      out_phonon_file << k+1 << "                 " << modes_mev[k] << "           " << red_mass[k] << endl;
    }

  // Calculate coordinate difference (deltaR) for each excited state combination
  double deltaR[natoms][3][band_comb];
  int band_count = 0;
      for(int m = 0; m < num_elevels; m++)
	{  
	  for(int n = m+1; n < num_elevels; n++)
	    {
	      for(int i = 0; i < natoms; i++)
		{
		  deltaR[i][0][band_count] = coord_cart[i][0][m] - coord_cart[i][0][n];
		  deltaR[i][1][band_count] = coord_cart[i][1][m] - coord_cart[i][1][n];
		  deltaR[i][2][band_count] = coord_cart[i][2][m] - coord_cart[i][2][n];
		}
	      band_count = band_count + 1;
	    }
	}

  // Calculate length of each deltaR vector - deltaR_length
  double deltaR_length[band_comb];
  for(int p = 0; p < band_comb; p++)
    {
      double deltaR_sum = 0.0;
      for(int i = 0; i < natoms; i++)
	{
	  deltaR_sum = deltaR_sum + deltaR[i][0][p]*deltaR[i][0][p] + deltaR[i][1][p]*deltaR[i][1][p] + deltaR[i][2][p]*deltaR[i][2][p];
        }
      deltaR_length[p] = sqrt(deltaR_sum);
    }

  for(int p = 0; p < band_comb; p++)
    {
      out_deltaR << "Band Combination " << p + 1 << endl;
      out_deltaR << "--------------------------------------------------------------------------------" << endl;
      out_deltaR << "DeltaR Length: " << deltaR_length[p] << endl; 
      out_deltaR << "--------------------------------------------------------------------------------" << endl;
      for(int i = 0; i < natoms; i++)
	{
	  out_deltaR << deltaR[i][0][p] << "  " << deltaR[i][1][p] << "  " << deltaR[i][2][p] << endl;
	}
    }

  // Calculate deltaR with mass-weighted cartesian coordinates
  double deltaR_sqrtm[natoms][3][band_comb];

  band_count = 0;
  for(int m = 0; m < num_elevels; m++)
    {
      for(int n = m+1; n < num_elevels; n++)
	{
	  for(int i = 0; i < natoms; i++)
	    {
	      deltaR_sqrtm[i][0][band_count] = coord_cart_sqrtm[i][0][m] - coord_cart_sqrtm[i][0][n];
	      deltaR_sqrtm[i][1][band_count] = coord_cart_sqrtm[i][1][m] - coord_cart_sqrtm[i][1][n];
	      deltaR_sqrtm[i][2][band_count] = coord_cart_sqrtm[i][2][m] - coord_cart_sqrtm[i][2][n];
	    }
	  band_count = band_count + 1;
	}
    }


  // Calculate dot product of deltaR vector with each mode to find displacement that enters Huang-Rhys formula
  double deltaR_mode_DP[natoms][nmodes][band_comb];
  double R_mode_dp_sum[nmodes][band_comb];
  double R_mode_dp_sum_norm[nmodes][band_comb];

  double deltaR_mode_DP_sqrtm[natoms][nmodes][band_comb];
  double R_mode_dp_sum_sqrtm[nmodes][band_comb];
  double R_mode_dp_sum_sqrtm_norm[nmodes][band_comb];

  for (int p = 0; p < band_comb; p++)
    {
      for (int k = 0; k < nmodes; k++)
	{
	  R_mode_dp_sum[k][p] = 0.0;
	  for (int i = 0; i < natoms; i++)
	    {
	      deltaR_mode_DP[i][k][p] = deltaR[i][0][p]*eig_modes_renorm[i][0][k] + deltaR[i][1][p]*eig_modes_renorm[i][1][k] + deltaR[i][2][p]*eig_modes_renorm[i][2][k];
	      R_mode_dp_sum[k][p] = R_mode_dp_sum[k][p] + deltaR_mode_DP[i][k][p];
	      deltaR_mode_DP_sqrtm[i][k][p] = deltaR_sqrtm[i][0][p]*eig_modes_renorm[i][0][k] + deltaR_sqrtm[i][1][p]*eig_modes_renorm[i][1][k] + deltaR_sqrtm[i][2][p]*eig_modes_renorm[i][2][k];
	      R_mode_dp_sum_sqrtm[k][p] = R_mode_dp_sum_sqrtm[k][p] + deltaR_mode_DP_sqrtm[i][k][p];
	    }
	  R_mode_dp_sum_norm[k][p] = R_mode_dp_sum[k][p]/norm_modes[k];
	  R_mode_dp_sum_sqrtm_norm[k][p] = R_mode_dp_sum_sqrtm[k][p]/norm_modes[k];
	}
    }

  // Calculate Huang Rhys Factor: S = (mode_mev/(2*hbar^2))*R_mode_dp_sum_norm^2 - using mass-weighted cartesian coordinates, so don't need to multiply by reduced mass!
  // Convert units - red_mass to kg, modes_mev to Joules, then use hbar in J*s, R_mode_dep_sum_norm from Ang to meters
  double HRfactor[nmodes][band_comb];
  double reorg_energy[band_comb];

  for (int p = 0; p < band_comb; p++)
    {
      reorg_energy[p] = 0.0;
      for (int k = 0; k < nmodes; k++)
	{
	  HRfactor[k][p] = 0.5*(mass_au_to_kg*(modes_mev[k]/1000)*energy_ev_to_J)/(hbar_Js*hbar_Js)*(R_mode_dp_sum_sqrtm_norm[k][p]*R_mode_dp_sum_sqrtm_norm[k][p]*length_ang_to_m*length_ang_to_m);
	  reorg_energy[p] = reorg_energy[p] + HRfactor[k][p]*(modes_mev[k]/1000);
	}
    }

  // Output all HR factors to 2D matrix (rows: transition energies, columns: phonon modes)
  for (int k = 0; k < nmodes; k++)
    {
      for (int p = 0; p < band_comb; p++)
        {
	  out_HR_matlab << HRfactor[k][p] << "  ";
        }
      out_HR_matlab << endl;
    }

  int p = 0;
  for (int m = 0; m < num_elevels; m++)
    {
      for (int n = m+1; n < num_elevels; n++)
	{
	  out_disp_norm << "Bands: " << elevels_min + m << " and " << elevels_min + n << endl;
	  out_disp_norm << "----------------------------------------------------------------------------------------------------------------" << endl;
	  out_disp_norm << "Normal mode:     Mode Energy (meV):     Mode Energy (cm-1):      Abs. Displacement (Ang):    Huang-Rhys Factor:      Reorg. Energy (eV): " << fixed << setprecision(6) << reorg_energy[p] << endl; 
	  out_disp_norm << "----------------------------------------------------------------------------------------------------------------" << endl;
	  for (int k = 0; k < nmodes; k++)
	    {
	      out_disp_norm << fixed << setprecision(6) << k + 1 << setprecision(6) << "                " << modes_mev[k] << "              " << (modes_mev[k]/1000)*ev_to_cm << scientific << setprecision(6) << "                  " << R_mode_dp_sum_norm[k][p] << "               " << HRfactor[k][p] << "             " << HRfactor[k][p]*(modes_mev[k]/1000) << endl;
	    }
	  p = p + 1;
	}
    }

  // Calculate unweighted density of states times BE+1 to compare to results without Franck Condon Factor
  double BE_dist_ab[nmodes];
  double pow_gauss_unweighted[nmodes][band_comb];
  double dos_unweighted_sum[band_comb];
  double smearing;
  smearing = 2*temp_ev*temp_ev; // smearing for gaussian

  band_count = 0;
  for (int m = 0; m < num_elevels; m++)
    {
      for (int n = m+1; n < num_elevels; n++)
	{
	  dos_unweighted_sum[band_count] = 0.0;
	  for (int k = 0 ; k < nmodes; k++)
	    {
	      pow_gauss_unweighted[k][band_count] = band_ediff_onerow[band_count]-(modes_mev[k]/1000);
	      BE_dist_ab[k] = 1/abs(exp(-(modes_mev[k]/1000)/temp_ev)-1);
	      dos_unweighted_sum[band_count] = dos_unweighted_sum[band_count] + BE_dist_ab[k]*(1/(sqrt(2*pi)*temp_ev))*exp(-(pow(pow_gauss_unweighted[k][band_count],2))/smearing);
	    }
	  band_count = band_count + 1;
	}
    }

  // Output lineshape matrix to file
  double dos_unweighted_matrix[num_elevels][num_elevels];
  band_count = 0;
  for (int m = 0; m < num_elevels; m++)
    {
      for (int n = m+1; n < num_elevels; n++)
	{
	  dos_unweighted_matrix[m][n] = dos_unweighted_sum[band_count];
	  dos_unweighted_matrix[n][m] = dos_unweighted_matrix[m][n];
	  band_count = band_count + 1;
	}
      dos_unweighted_matrix[m][m] = 0.0;
    }

  for (int m = 0; m < num_elevels; m++)
    {
      for (int n = 0; n < num_elevels; n++)
	{
	  out_dos_total << dos_unweighted_matrix[m][n] << " ";
	}
      out_dos_total << endl;
    }


  return 0;
}

// FUNCTIONS:

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
