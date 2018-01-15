// Program calculates the nonadiabatic coupling between excited state described by linear response TDDFT using Casida's description of the excited state wave function.  The code calculates electron-phonon, nonadiabatic coupling based on Tavernelli et al 2009, 124107.  The Kohn-Sham single-particle orbital overlaps are calculated in a previous step and used as input here along with LR-TDDFT excitation information.  The DFT wave functions must be in WFX format, and the TDDFT excitations and eigenvectors must come from Gaussian09's implementation of Casida LR-TDDFT.
// ----INPUT----
// 1) input file: file called from command line with code execution that must include the following lines: (must only have one command per line, otherwise code will not read correctly)
      //dft_dir_name = scf          / Name of directory prefix where scf calculations and K-S overlap matrix is located (e.g., scf1, scf2, etc.)
      //td_dir_name = tddft         / Name of directory within dft directory containing LR-TDDFT calculation output
      //overlap_file_name = output_ovlap       / Name of K-S overlap matrix file within dft_dir_name
      //tddft_file_name = gauss_lrtddft.out    / Name of LR-TDDFT output file

      //nsteps = 300                / Total number of MD timesteps to use to calculate coupling
      //gauss_homo = 81             / Gaussian KS index of HOMO
      //gauss_emin = 40             / Minimumum Gaussian KS index included in overlap matrix (output_ovlap)
      //gauss_emax = 120            / Maximum Gaussian KS index included in overlap matrix
      //num_exc = 6                 / Total number of excitations included in LR-TDDFT calculation
      //time_step = 0.25            / Time step of MD trajectory

      //dft_ordering = yes          / Turn on ordering of Kohn-Sham states from step to step by comparing overlaps
      //td_ordering = yes           / Turn on ordering of LR-TDDFT states from step to step by comparing overlaps
      //same_sign = yes             / Turn on option to ensure the same sign of each coefficient at each time step
      //use_threshold = no          / Turn on threshold to not use any coefficient in coupling that changes by 'threshold'
      //threshold = 0.05            / Threshold value for thresholding out large changes in coefficients

// 2) OUTPUT: all output for each time step is placed in an scf($i)/coupling folder created by the program
// Primary output: In base directory in which code was run: 
//                 final_NAC_matrix_elec, final_NAC_matrix_hole, and final_NAC_matrix_total are the final nonadiabatic coupling matrices (already squared) for total, electron, and hole components that should be used as the square of the electronic coupling matrix in relaxation dynamics.
//                 NAC_timeseries_elec, NAC_timeseries_hole, NAC_timeseries_total give the cumulative, averaged nonadiabatic coupling between each excitation combination as a function of time-step.  Each row corresponds to excitation combination (starting with Excitation 1-1, 1-2, 1-3, etc.), and each column is a time step.  Use this data to check for convergence of coupling values with respect to total number of time steps.
// Secondary output: within each scf($i)/coupling folder, additional output files contain information about coupling at that time step:
//                   output_coeff_exc/output_coeff_de_exc: matrices of LR-TDDFT coefficients read in from Gaussian LR-TDDFT output file for (i) and (i+1) time steps.
//                   output_elec_details/output_hole_details: List of nonadiabatic coupling values between each K-S state making up each LR-TDDFT excitation.  Analyze for more information about physics behind coupling values between excited states.
//                   output_NAC_total/output_NAC_elec/output_NAC_hole: unsquared, nonadiabatic coupling matrix for this time step only.
//                   output_KS_flip/output_TD_flip/output_TD_overlap/output_TD_overlap_ordered: optional outputs given if 'dft_ordering' and 'td_ordering' turned on.  Provide information about overlap matrices and list which states have been flipped in ordering at this time step.  

// Assume that the code is run in base directory from which all DFT calculations are in subdirectory down (scf1, scf2, etc.) and TDDFT output is in subdirectory within each of these scf folders.

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <sys/stat.h>

using namespace std;

//==============================================
//= PRE-MAIN DELCARATIONS AND PROTOTYPES =======
//==============================================

// FUNCTION PROTOTYPES:
void readInput(const char * input); // Read main input file
double** alloc2Darray(int x, int y); // Allocate 2D array filled with zeroes
int** alloc2DarrayInt(int x, int y); // Allocate 2D array filled with ints
double*** alloc3Darray(int x, int y, int z); // Allocate 3D array filled with zeroes
void destroy2Darray(double** array, int x, int y); // Delete 2D array
void destroy2DarrayInt(int** array, int x, int y);
void destroy3Darray(double*** array, int x, int y, int z); // Delete 3D array
void readKSoverlap(double ** oMat, string dirName, int index, string oName, int eLevels); // Read KS overlap file
void readTDcoeff(double ** TDexcEnergy, double *** TDexcmat, double *** TDdeMat, string dirName, int index, string TDdirName, string TDoutputName, int exc, int eMin, int eLevels); // Read in LR-TDDFT output file
void switchTDsign(double *** TDmat1, double *** TDmat2, int exc, int eLevels); // Switch TD coefficient signs
void outputTDcoeff(double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, string output_path); // Output LR-TDDFT coefficients
void calculateNACelec(double ** NACmat, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path); // Electronic nonadiabatic coupling calculation
void calculateNAChole(double ** NACmat, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path); // Hole nonadiabatic coupling calculation
void calculateNACground(double * NACvec, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path); // Excited-ground state nonadiabatic coupling calculation
void constructNACmatrix(double *** NACmatTotalTimeStep, double *** NACmatElecTimeStep, double *** NACmatHoleTimeStep, double ** NACmatElec, double ** NACmatHole, double * NACvecGr, int tStep, int exc, string output_path);
void thermalAvgNACmatrix(double ** avgNACtotal, double ** avgNACelec, double ** avgNAChole, double *** NACmatTotalTimeStep, double *** NACmatElecTimeStep, double *** NACmatHoleTimeStep, int exc, int totalSteps); // Calculate averaged NAC matrix over MD time steps and output final matrix files in root directory where code was run
void orderDFTstates(double ** oMat, int ** switchMat, double *** TDmat1, double *** TDdeMat1, double *** TDmat2, double *** TDdeMat2, int eLevels, int exc, int tStep, int eMin, string output_path); // Function that orders KS states from timestep to timestep using KS overlaps
void orderTDDFTstates(double ** TDexcEn1, double ** TDexcEn2, double ** TDexcEnOrd1, double ** TDexcEnOrd2, double ** TDoverlap, int ** switchMat, double *** TDmat1, double *** TDdeMat1, double *** TDmat2, double *** TDdeMat2, int eLevels, int exc, int tStep, int eMin, string output_path); // Function that orders LR-TDDFT states from timestep to timestep using LR-TDDFT overlaps
void swapDouble(double * first, double * second); // swap function
// GLOBAL CONSTANTS: Input constants that will be not be changed
const double pi = 3.14159;
const double tempEv = .0257; // room temperature in eV
const double hbarEvFs = 0.6582119; // hbar planck's constant in eV*fs
const double hbarJs = 1.0546E-34; // hbar planck's constant in J*s
const double massAuToKg = 1.6726E-27; // convert mass from atomic units to kg
const double energyEvToJ = 1.619E-19; // convert eV to Joules
const double lengthAngToM = 1E-10; // convert angstroms to meters
const double evToCm = 8065.73; // convert ev to cm by multiplying by this number
const double hartToEv = 27.2107; // convert Hartree to eV

// GLOBAL VARIABLES:
// Input variables: set as global variables so they can be set in readInput function and used globally
string dftDirName; // directory base name containing SCF output and KS overlap matrix
string tdDirName; // name of subdirectory within dftDirName that contains LR-TDDFT data
string overlapName; // name of file containing KS overlap matrix in each dftDirName
string tddftName; // name of LR-TDDFT output file in tdDirName
string str_nSteps; // total number of MD timesteps to be used in LR-TDDFT coupling calculation
int nSteps; // integer conversion of string nsteps
string str_gaussHomo; // Gaussian09 index label for HOMO
int gaussHomo;
int gaussLumo;
string str_gaussEmin; // Gaussian09 index label for minimum energy eigenvalue
int gaussEmin;
string str_gaussEmax; // Gaussian09 index label for maximum energy eigenvalue
int gaussEmax;
int nElevels; // number of Kohn-Sham energy levels
string str_numExc; // Number of LR-TDDFT excited states to include in coupling calculation
int numExc;
string str_timeStep; // MD time step in femtoseconds (fs)
double timeStep; // Double conversion of time step
string dftOrdering; // Turns on or off ordering of DFT Kohn-Sham states from timestep to timestep
string tdOrdering; // Turns on or off ordering of TDDFT states to preserve order
string sameSign; // Turns on or off keeping the sign of all TDDFT coefficients the same every step
string useThreshold; // Turns on or off using a threshold that cuts off any coefficients that change above threshold
string str_thresh; // Value of threshold for cutting off TDDFT coefficients if turned on
double threshold; 

//=============================
//= MAIN PROGRAM ------ =======
//=============================

int main(int argc, char *argv[])
{

  // Read input file and output info to screen
  readInput(argv[1]);

  // Output general information
  cout << "*** LR-TDDFT Coupling Input Parameters ***" << endl; 
  cout << "--Filenames--" << endl;
  cout << "DFT Directory Prefix: " << dftDirName << endl;
  cout << "TDDFT Subdirectory Name: " << tdDirName << endl;
  cout << "Kohn-Sham Overlap Matrix Filename: " << overlapName << endl;
  cout << "LR-TDDFT Ouput Filename: " << tddftName << endl;
  cout << "Total MD Timesteps: " << nSteps << endl;
  cout << endl;
  cout << "--Parameters--" << endl;
  cout << "HOMO Index: " << gaussHomo << endl;
  cout << "Minimum Index: " << gaussEmin << endl;
  cout << "Maximum Index: " << gaussEmax << endl;
  cout << "Number of K-S Energy Levels: " << nElevels << endl;
  cout << "Number of LR-TDDFT excitations: " << numExc << endl;
  cout << "MD timestep (fs): " << timeStep << endl;
  cout << endl;
  cout << "--Other Options--" << endl;
  cout << "DFT State Ordering: " << dftOrdering << endl;
  cout << "LR-TDDFT State Ordering: " << tdOrdering << endl;
  cout << "Keep Same Sign of LR-TDDFT Coefficients: " << sameSign << endl;
  cout << "Use Threshold for Coefficient Changes: " << useThreshold << endl;
  cout << "Threshold Value: " << threshold << endl;
  cout << "*** End Input Parameters ***" << endl;
  cout << endl;

  // ----ARRAY ALLOCATION---- //
  // Create pointers to 1st and 2nd KS overlap matrces and allocate 2D arrays
  // Note: The values pointed to will be written over each loop iteration for each pair of timesteps
  // First time step (t)
  double** KSoMat = alloc2Darray(nElevels, nElevels);

  // Create 2D array for LR-TDDFT excitation energies (in nm) and oscillator strengths
  double** TDenergies1 = alloc2Darray(numExc,2);
  double** TDenergies2 = alloc2Darray(numExc,2);
  double** TDenergiesOrd1 = alloc2Darray(numExc,2);
  double** TDenergiesOrd2 = alloc2Darray(numExc,2);

  // Create pointers to 1st and 2nd LR-TDDFT excitation coefficient matrices and allocate arrays with zeroes
  // Note: The values pointed to will be written over each loop iteration for each pair of timesteps
  double ***TDexcCoeff1 = alloc3Darray(numExc, nElevels, nElevels);
  double ***TDexcCoeff2 = alloc3Darray(numExc, nElevels, nElevels);

  // Same for de-excitation coefficients:
  double ***TDdeCoeff1 = alloc3Darray(numExc, nElevels, nElevels);
  double ***TDdeCoeff2 = alloc3Darray(numExc, nElevels, nElevels);

  // Pointer to 2D array of electronic components of nonadiabatic coupling (see Tavernelli Eq 58) calculated at each time step
  double** tsNACelec = alloc2Darray(numExc, numExc);
  double** tsNAChole = alloc2Darray(numExc, numExc);

  // Pointer to 1D array of nonadiabatic coupling between ground state and each excited state
  double * tsNACground;
  tsNACground = new double[numExc];

  // Pointer to 3D array containing LR-TDDFT nonadiabatic coupling matrices for each time step
  // --dimensions of numExc+1 because first row and column are ground state couplings
  // --these are then squared and averaged across the time dimension for the thermally averaged coupling matrix
  double ***totalNACall = alloc3Darray(nSteps, numExc+1, numExc+1);
  double ***totalNACelec = alloc3Darray(nSteps, numExc+1, numExc+1);
  double ***totalNAChole = alloc3Darray(nSteps, numExc+1, numExc+1);
  double **avgSqNACall = alloc2Darray(numExc+1, numExc+1);
  double **avgSqNACelec = alloc2Darray(numExc+1, numExc+1);
  double **avgSqNAChole = alloc2Darray(numExc+1, numExc+1);

  // Pointer to 2D array containing binary elements where 1 indicates these columns have been switched
  int **KSswitchMat = alloc2DarrayInt(nSteps, nElevels);

  // Point to 2D array containing LR-TDDFT overlap elements
  double ** TDoMat = alloc2Darray(numExc, numExc); 
  int **TDswitchMat = alloc2DarrayInt(nSteps, numExc);

  // ----LR-TDDFT COUPLING CALCULATION---- //

  cout << "-----LR-TDDFT COUPLING CALCULATION-----" << endl;

  // Loop over all MD timesteps to calculate coupling for each and store in matrix
  for (int i = 0; i < nSteps; i++)
    {
      cout << "--Time Step " << i+1 << "--" << endl;

      // Create folder in 'i'th scf directory called 'coupling' where all output files will go:
      // Convert directory index to string
      stringstream ss_mainIndex;
      ss_mainIndex << (i+1);
      string strMainIndex = ss_mainIndex.str();
      string outputPath = dftDirName + strMainIndex + "/coupling/";
      // Make directory
      mkdir(outputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH| S_IXOTH);

      // Read in KS overlap matrix of ith and (i+1)th time steps and store in KSoMat
      readKSoverlap(KSoMat, dftDirName, i+1, overlapName, nElevels);

      // Read in 1st and 2nd LR-TDDFT excitation coefficients of ith/(i+1)th time steps
      cout << "Reading in LR-TDDFT energies for time step t (" << i+1 << "):" << endl;
      readTDcoeff(TDenergies1, TDexcCoeff1, TDdeCoeff1, dftDirName, i+1, tdDirName, tddftName, numExc, gaussEmin, nElevels);
      cout << "Reading in LR-TDDFT energies for time step t + delta (" << i+2 << "):" << endl;
      readTDcoeff(TDenergies2, TDexcCoeff2, TDdeCoeff2, dftDirName, i+2, tdDirName, tddftName, numExc, gaussEmin, nElevels);

      // Output LR-TDDFT excitation and de-excitation coefficients to output files folder within LR-TDDFT output file
      outputTDcoeff(TDexcCoeff1, TDexcCoeff2, TDdeCoeff1, TDdeCoeff2, outputPath);

      // If dft_ordering selected, this step checks that the diagonal elements of each Kohn-Sham overlap matrix are the largest, otherwise it is assumed that states have changed ordering.  Columns are switched to ensure the same ordering and a matrix is filled with 1's to indicate which states have switched so the same flipping scheme is applied to each subsequent time step.
      if (dftOrdering == "yes")
	{
	  cout << "You have turned on DFT state ordering across all time steps!" << endl;
	  orderDFTstates(KSoMat, KSswitchMat, TDexcCoeff1, TDdeCoeff1, TDexcCoeff2, TDdeCoeff2, nElevels, numExc, i, gaussEmin, outputPath);
	}

      // Check that coefficients have not switched sign:
      // For an unknown reason, sometimes all the signs of the TDDFT coefficients for a given excited state will change all signs: this will lead to unphysically large coupling values.  Here, we check that signs have not changed, and if they have, we change them back.  This could affect states with small coefficients that do physically switch from positive to negative values, but these will contribute little to the coupling anyway and will introduce an acceptable level of error.
      if (sameSign == "yes")
	{
	  cout << "You have turned on the fix to force coefficients from adjacent time steps to have the same sign!" << endl;
	  switchTDsign(TDexcCoeff1, TDexcCoeff2, numExc, nElevels);
	  switchTDsign(TDdeCoeff1, TDdeCoeff2, numExc, nElevels);
	}

      // If td_ordering is selected, this step checks the overlap between all LR-TDDFT excitations between adjacent timesteps and switches state ordering to ensure that largest overlaps are along the diagonal of the overlap matrix.  Columns are switched to ensure the same ordering and a matrix is filled with 1's to indicate which states have been flipped so the cumulative switching can be applied to future time steps.
      if (tdOrdering == "yes")
	{
	  cout << "You have turned on LR-TDDFT state ordering across all time steps!" << endl;
	  orderTDDFTstates(TDenergies1, TDenergies2, TDenergiesOrd1, TDenergiesOrd2, TDoMat, TDswitchMat, TDexcCoeff1, TDdeCoeff1, TDexcCoeff2, TDdeCoeff2, nElevels, numExc, i, gaussEmin, outputPath);
	}

      // Calculate nonadiabatic coupling using LR-TDDFT excitations based on Tavernelli Eq 58
      // First term - electronic component of NAC
      calculateNACelec(tsNACelec, TDexcCoeff1, TDexcCoeff2, TDdeCoeff1, TDdeCoeff2, KSoMat, numExc, gaussHomo, gaussLumo, gaussEmin, gaussEmax, timeStep, outputPath);
      // Second term - hole component of NAC
      calculateNAChole(tsNAChole, TDexcCoeff1, TDexcCoeff2, TDdeCoeff1, TDdeCoeff2, KSoMat, numExc, gaussHomo, gaussLumo, gaussEmin, gaussEmax, timeStep, outputPath);
      // Excited state - ground state coupling
      calculateNACground(tsNACground, TDexcCoeff1, TDexcCoeff2, TDdeCoeff1, TDdeCoeff2, KSoMat, numExc, gaussHomo, gaussLumo, gaussEmin, gaussEmax, timeStep, outputPath);

      // Create nonadiabatic coupling matrices for this time step (combining excited and ground state coupling)
      // Stored in 3D array with time as one dimension so the 2D arrays can be reused for each step of loop
      constructNACmatrix(totalNACall, totalNACelec, totalNAChole, tsNACelec, tsNAChole, tsNACground, i, numExc, outputPath);
      

    }

  // Thermally averaged nonadiabatic coupling
  // After loop, square the nonadiabatic coupling matrix elements at each time step and then take average over total number of time steps
  thermalAvgNACmatrix(avgSqNACall, avgSqNACelec, avgSqNAChole, totalNACall, totalNACelec, totalNAChole, numExc, nSteps);
  
  return 0;
}

//=============================
// FUNCTIONS:
//=============================

//----Read input file----//
void readInput(const char * input)
{
  // Input data file
  ifstream inputFile;
  inputFile.open(input);

  string inputLine;
  // Read through all lines, read in relevant variable if keyword found
  while (getline(inputFile,inputLine))
    {
      if(inputLine.find("dft_dir_name") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  dftDirName = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  dftDirName.erase(remove(dftDirName.begin(),dftDirName.end(),' '),dftDirName.end());
	}
      if(inputLine.find("td_dir_name") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  tdDirName = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  tdDirName.erase(remove(tdDirName.begin(),tdDirName.end(),' '),tdDirName.end());
	}
      if(inputLine.find("overlap_file_name") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  overlapName = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  overlapName.erase(remove(overlapName.begin(),overlapName.end(),' '),overlapName.end());
	}
      if(inputLine.find("tddft_file_name") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  tddftName = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  tddftName.erase(remove(tddftName.begin(),tddftName.end(),' '),tddftName.end());
	}
      if(inputLine.find("nsteps") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_nSteps = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_nSteps.erase(remove(str_nSteps.begin(),str_nSteps.end(),' '),str_nSteps.end());
	  stringstream(str_nSteps) >> nSteps;
	}
      if(inputLine.find("gauss_homo") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_gaussHomo = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_gaussHomo.erase(remove(str_gaussHomo.begin(),str_gaussHomo.end(),' '),str_gaussHomo.end());
	  stringstream(str_gaussHomo) >> gaussHomo;
	  gaussLumo = gaussHomo + 1;
	}
      if(inputLine.find("gauss_emin") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_gaussEmin = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_gaussEmin.erase(remove(str_gaussEmin.begin(),str_gaussEmin.end(),' '),str_gaussEmin.end());
	  stringstream(str_gaussEmin) >> gaussEmin;
	}
      if(inputLine.find("gauss_emax") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_gaussEmax = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_gaussEmax.erase(remove(str_gaussEmax.begin(),str_gaussEmax.end(),' '),str_gaussEmax.end());
	  stringstream(str_gaussEmax) >> gaussEmax;
	}
      if(inputLine.find("num_exc") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_numExc = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_numExc.erase(remove(str_numExc.begin(),str_numExc.end(),' '),str_numExc.end());
	  stringstream(str_numExc) >> numExc;
	}
      if(inputLine.find("time_step") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_timeStep = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_timeStep.erase(remove(str_timeStep.begin(),str_timeStep.end(),' '),str_timeStep.end());
	  stringstream(str_timeStep) >> timeStep;
	}
      if(inputLine.find("dft_ordering") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  dftOrdering = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  dftOrdering.erase(remove(dftOrdering.begin(),dftOrdering.end(),' '),dftOrdering.end());
	}
      if(inputLine.find("td_ordering") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  tdOrdering = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  tdOrdering.erase(remove(tdOrdering.begin(),tdOrdering.end(),' '),tdOrdering.end());
	}
      if(inputLine.find("same_sign") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  sameSign = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  sameSign.erase(remove(sameSign.begin(),sameSign.end(),' '),sameSign.end());
	}
      if(inputLine.find("use_threshold") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  useThreshold = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  useThreshold.erase(remove(useThreshold.begin(),useThreshold.end(),' '),useThreshold.end());
	}
      if(inputLine.find("threshold") != string::npos)
	{
	  size_t eqPos1 = inputLine.find("=");
	  str_thresh = inputLine.substr(eqPos1+1);
	  // Delete white space around string
	  str_thresh.erase(remove(str_thresh.begin(),str_thresh.end(),' '),str_thresh.end());
	  stringstream(str_thresh) >> threshold;
	}
    }
  nElevels = gaussEmax - gaussEmin + 1; // set total number of energy levels
  inputFile.close(); // close input file
}

//----Create dynamic 2D array filled with zeroes----//
double** alloc2Darray(int x, int y)
{
  double** array = new double*[x];
  for (int i = 0; i < x; i++)
    {
      array[i] = new double[y];
      for (int j = 0; j < y; j++)
	{
	  array[i][j] = 0.0;
	}
    } 

  return array;
}

//----Destroy dynamic 2D array----//
void destroy2Darray(double** array, int x, int y)
{
  for (int i = 0; i < x; i++)
    {
      delete [] array[i];
    }
  delete [] array;
}

//----Create dynamic 2D array filled with integer zeroes----//
int** alloc2DarrayInt(int x, int y)
{
  int** array = new int*[x];
  for (int i = 0; i < x; i++)
    {
      array[i] = new int[y];
      for (int j = 0; j < y; j++)
	{
	  array[i][j] = 0;
	}
    } 

  return array;
}

//----Destroy dynamic 2D array of ints----//
void destroy2DarrayInt(int** array, int x, int y)
{
  for (int i = 0; i < x; i++)
    {
      delete [] array[i];
    }
  delete [] array;
}


//----Create dynamic 3D array filled with zeroes----//
double*** alloc3Darray(int x, int y, int z)
{
  double ***array = new double**[x];
  for (int i = 0; i < x; i++)
    {
      array[i] = new double*[y];
      for (int j = 0; j < y; j++)
	{
	  array[i][j] = new double[z];
	  for (int k = 0; k < z; k++)
	    {
	      array[i][j][k] = 0.0;
	    }
	}
    }

  return array;
}

//----Destroy dynamic 3D array----//
void destroy3Darray(double*** array, int x, int y, int z)
{
  for (int i = 0; i < x; i++)
    {
      for (int j = 0; j < y; j++)
	{
	  delete [] array[i][j];
	}
      delete [] array[i];
    }
  delete [] array;
}

//----Read KS orbital overlap file----//
void readKSoverlap(double ** oMat, string dirName, int index, string oName, int eLevels)
{
  // Convert directory index to string
  stringstream ss_index;
  ss_index << index;
  string str_index = ss_index.str();

  // Create path name to overlap matrix
  string path = dirName + str_index + "/" + oName;

  // Open and read in overlap matrix
  ifstream KSoverlapFile;
  KSoverlapFile.open(path.c_str());
  for (int i = 0; i < eLevels; i++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  oMat[i][j] = 0.0;
	  KSoverlapFile >> oMat[i][j];
	}
    }
  KSoverlapFile.close();
}

//----Read excitation and de-excitation coefficients from G09 LR-TDDFT output file----//
void readTDcoeff(double ** TDexcEnergy, double *** TDexcMat, double *** TDdeMat, string dirName, int index, string TDdirName, string TDoutputName, int exc, int eMin, int eLevels)
{
  // Convert directory index to string
  stringstream ss_index;
  ss_index << index;
  string str_index = ss_index.str();

  // Create path name to G09 LR-TDDFT output file
  string path = dirName + str_index + "/" + TDdirName + "/" + TDoutputName;

  // Open and read in output file
  ifstream TDoutputFile;
  TDoutputFile.open(path.c_str());

  // Initialize variables:
  string searchCoord = "Excited State ";
  string temp;
  size_t posArrow;
  string arrow;
  string strTmpOcc;
  int tmpOcc;
  string strTmpUnoccCoeff;
  string strTmpUnocc;
  int tmpUnocc;
  string strTmpCoeff;
  double tmpCoeff;
  int tmpOccShift;
  int tmpUnoccShift;
  int lineBreak;
  string strAllLine;
  string strTmpEnergy;
  string strTmpEnergy2;
  double tmpEnergy;
  string strTmpOS;
  string strTmpOS2;
  double tmpOS;

  // Since this matrix is replaced in each loop of time steps, set all elements to 0.0 to start
  for (int i = 0; i < exc; i++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  for (int k = 0; k < eLevels; k++)
	    {
	      TDexcMat[i][j][k] = 0.0;
	      TDdeMat[i][j][k] = 0.0;
	    }
	}
    }

  // Read in indices and coefficients for first time step (t)
  // First, input coefficients from first excitation - this is different from rest because it has a string right after the last coefficient rather than a space, so we must treat differently to tell program when to end
  int firstK = 0;
    while (firstK < 1)
    {
      getline(TDoutputFile,temp);
      if(temp.find(searchCoord) != string::npos)
	{
	  // Store energy (in nm) and oscillator strength to 1D array to be used in LR-TDDFT state ordering
	  strTmpEnergy = temp.substr(temp.find("eV"));
	  strTmpEnergy2 = strTmpEnergy.substr(2,strTmpEnergy.find("nm")-2);
	  stringstream(strTmpEnergy2) >> tmpEnergy;
	  strTmpOS = temp.substr(temp.find("f="),temp.find("<"));
	  strTmpOS2 = strTmpOS.substr(2,strTmpOS.find(" "));
	  stringstream(strTmpOS2) >> tmpOS;
	  cout << temp << endl;
	  TDexcEnergy[0][0] = tmpEnergy;
	  TDexcEnergy[0][1] = tmpOS;
	  lineBreak = 1;
	  while (lineBreak > 0)
	    {
	      getline(TDoutputFile,strAllLine);
	      if(strAllLine.find("This") != string::npos) // End loop after running through all coefficients and hitting the word "This" in line directly below
		{
		  lineBreak = 0;
		}
	      else if(strAllLine.find("->") != string::npos) // Excitations: break string into variables needed
		{
		  posArrow = strAllLine.find("->");
		  strTmpOcc = strAllLine.substr(0,posArrow); // occupied KS orbital index
		  strTmpUnoccCoeff = strAllLine.substr(posArrow,strAllLine.length()); // unoccupied KS orbital index and coefficient
		  strTmpUnocc = strTmpUnoccCoeff.substr(2,strTmpUnoccCoeff.find("  "));  // Break unocc and coeff string
		  strTmpCoeff = strTmpUnoccCoeff.substr(strTmpUnoccCoeff.find("  "));

		  stringstream(strTmpOcc) >> tmpOcc; // Convert to int
		  stringstream(strTmpUnocc) >> tmpUnocc; // Convert to int
		  stringstream(strTmpCoeff) >> tmpCoeff; // Convert to double
		  tmpCoeff = tmpCoeff*sqrt(2.0); // Normalize coefficients to 1 (Gaussian has them normalized to 0.5)
		  tmpOccShift = tmpOcc - eMin + 1; // shift indices with respect to lowest Gaussiasn index as 1
		  tmpUnoccShift = tmpUnocc - eMin + 1;
		  TDexcMat[firstK][tmpOccShift-1][tmpUnoccShift-1] = tmpCoeff; // Fill excitation coefficient matrix for first excitation
		}
	      else if(strAllLine.find("<-") != string::npos) // De-excitations: break string into variables needed
		{
		  posArrow = strAllLine.find("<-");
		  strTmpOcc = strAllLine.substr(0,posArrow); // occupied KS orbital index
		  strTmpUnoccCoeff = strAllLine.substr(posArrow,strAllLine.length()); // unoccupied KS orbital index and coefficient
		  strTmpUnocc = strTmpUnoccCoeff.substr(2, strTmpUnoccCoeff.find("  ")); // Break unocc and coeff string
		  strTmpCoeff = strTmpUnoccCoeff.substr(strTmpUnoccCoeff.find("  "));

		  stringstream(strTmpOcc) >> tmpOcc; // Convert to int
		  stringstream(strTmpUnocc) >> tmpUnocc; // Convert to int
		  stringstream(strTmpCoeff) >> tmpCoeff; // Convert to double
		  tmpCoeff = tmpCoeff*sqrt(2.0); // Normalize coefficients to 1 (Gaussian has them normalized to 0.5)
		  tmpOccShift = tmpOcc - eMin + 1; // shift indices with respect to lowest Gaussiasn index as 1
		  tmpUnoccShift = tmpUnocc - eMin + 1;
		  TDdeMat[firstK][tmpOccShift-1][tmpUnoccShift-1] = tmpCoeff; // Fill de-excitation coefficient matrix
		}
	    }
	  firstK = 1; // Break loop after final de-excitation coefficient
	}
    }

  // Next, go through all excitations and fill in coefficient matrices (excitations and de-excitations)
  int secondK = 0;

  while (secondK < exc)
    {
      getline(TDoutputFile,strAllLine);
      if(strAllLine.find("Excited State ") != string::npos)  // If line contains "Excited State ", then signals next excited state, increment k
	{
	  secondK = secondK+1;
	  // Store energy (in nm) and oscillator strength to 1D array to be used in LR-TDDFT state ordering
	  strTmpEnergy = strAllLine.substr(strAllLine.find("eV"));
	  strTmpEnergy2 = strTmpEnergy.substr(2,strTmpEnergy.find("nm")-2);
	  stringstream(strTmpEnergy2) >> tmpEnergy;
	  strTmpOS = strAllLine.substr(strAllLine.find("f="),strAllLine.find("<"));
	  strTmpOS2 = strTmpOS.substr(2,strTmpOS.find(" "));
	  stringstream(strTmpOS2) >> tmpOS;
	  cout << strAllLine << endl;
	  TDexcEnergy[secondK][0] = tmpEnergy;
	  TDexcEnergy[secondK][1] = tmpOS;
	}
      else if(strAllLine.find("SavETr:") != string::npos) // Signals end  of excited states, incremen k by one to end loop
	{
	  secondK = secondK + 1;
	}
      else if(strAllLine.find("->") != string::npos) // Excitations: break strings into variables needed
	{
	  posArrow = strAllLine.find("->");
	  strTmpOcc = strAllLine.substr(0,posArrow); // occupied KS orbital index
	  strTmpUnoccCoeff = strAllLine.substr(posArrow,strAllLine.length()); // unoccupied KS orbital index and coefficient
	  strTmpUnocc = strTmpUnoccCoeff.substr(2,strTmpUnoccCoeff.find("  "));  // Break unocc and coeff string
	  strTmpCoeff = strTmpUnoccCoeff.substr(strTmpUnoccCoeff.find("  "));

	  stringstream(strTmpOcc) >> tmpOcc; // Convert to int
	  stringstream(strTmpUnocc) >> tmpUnocc; // Convert to int
	  stringstream(strTmpCoeff) >> tmpCoeff; // Convert to double
	  tmpCoeff = tmpCoeff*sqrt(2.0); // Normalize coefficients to 1 (Gaussian has them normalized to 0.5)
	  tmpOccShift = tmpOcc - eMin + 1; // shift indices with respect to lowest Gaussiasn index as 1
	  tmpUnoccShift = tmpUnocc - eMin + 1;
	  TDexcMat[secondK][tmpOccShift-1][tmpUnoccShift-1] = tmpCoeff; // Fill coefficient matrix
	}
      else if(strAllLine.find("<-") != string::npos) // De-excitations: break strings into variables needed
	{
	  posArrow = strAllLine.find("<-");
	  strTmpOcc = strAllLine.substr(0,posArrow); // occupied KS orbital index
	  strTmpUnoccCoeff = strAllLine.substr(posArrow,strAllLine.length()); // unoccupied KS orbital index and coefficient
	  strTmpUnocc = strTmpUnoccCoeff.substr(2, strTmpUnoccCoeff.find("  ")); // Break unocc and coeff string
	  strTmpCoeff = strTmpUnoccCoeff.substr(strTmpUnoccCoeff.find("  "));

	  stringstream(strTmpOcc) >> tmpOcc; // Convert to int
	  stringstream(strTmpUnocc) >> tmpUnocc; // Convert to int
	  stringstream(strTmpCoeff) >> tmpCoeff; // Convert to double
	  tmpCoeff = tmpCoeff*sqrt(2.0); // Normalize coefficients to 1 (Gaussian has them normalized to 0.5)
	  tmpOccShift = tmpOcc - eMin + 1; // shift indices with respect to lowest Gaussiasn index as 1
	  tmpUnoccShift = tmpUnocc - eMin + 1;
	  TDdeMat[secondK][tmpOccShift-1][tmpUnoccShift-1] = tmpCoeff; // Fill coefficient matrix
	}
    }
  TDoutputFile.close();
}

//----Switch signs of all LR-TDDFT excitation coefficients to match those of previous time step----//
void switchTDsign(double *** TDmat1, double *** TDmat2, int exc, int eLevels)
{
  for (int i = 0; i < exc; i ++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  for (int k = 0; k < eLevels; k++)
	    {
	      if((TDmat2[i][j][k] < 0.0 && TDmat1[i][j][k] > 0.0) || (TDmat2[i][j][k] > 0.0 && TDmat1[i][j][k] < 0.0))
		{
		  TDmat2[i][j][k] = -1*TDmat2[i][j][k];
		}
	      else
		{
		  TDmat2[i][j][k] = TDmat2[i][j][k];
		}
	    }
	}
    }
}

//----Output LR-TDDFT coefficients from G09 output file to file in matrix format----//
void outputTDcoeff(double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, string output_path)
{
  // Open and read in output files for excitation and de-excitation coefficients for t and t+Delta time steps
  // Excitation coefficients of first time step (t)
  ofstream TDexcOutputFile1;
  TDexcOutputFile1.open((output_path+"output_coeff_exc_1").c_str());
  // Excitation coefficients of second time step (t + delta)
  ofstream TDexcOutputFile2;
  TDexcOutputFile2.open((output_path+"output_coeff_exc_2").c_str());
  // De-excitation coefficients of second time step (t)
  ofstream TDdeOutputFile1;
  TDdeOutputFile1.open((output_path+"output_coeff_de_exc_1").c_str());
  // De-excitation coefficients of second time step (t + delta)
  ofstream TDdeOutputFile2;
  TDdeOutputFile2.open((output_path+"output_coeff_de_exc_2").c_str());
  
  for (int i = 0; i < numExc; i++)
    {
      TDexcOutputFile1 << "Excitation " << i+1 << ":" << endl;
      TDexcOutputFile2 << "Excitation " << i+1 << ":" << endl;
      TDdeOutputFile1 << "Excitation " << i+1 << ":" << endl;
      TDdeOutputFile2 << "Excitation " << i+1 << ":" << endl;

      for (int j = 0; j < nElevels; j++)
	{
	  for (int k = 0; k < nElevels; k++)
	    {
	      TDexcOutputFile1 << scientific << setprecision(6) << TDmat1[i][j][k] << '\t';
	      TDexcOutputFile2 << scientific << setprecision(6) << TDmat2[i][j][k] << '\t';
	      TDdeOutputFile1 << scientific << setprecision(6) << TDdeMat1[i][j][k] << '\t';
	      TDdeOutputFile2 << scientific << setprecision(6) << TDdeMat2[i][j][k] << '\t';
	    }
	  TDexcOutputFile1 << endl;
	  TDexcOutputFile2 << endl;
	  TDdeOutputFile1 << endl;
	  TDdeOutputFile2 << endl;
	}
    }
  TDexcOutputFile1.close();
  TDexcOutputFile2.close();
  TDdeOutputFile1.close();
  TDdeOutputFile2.close();
}

//----Calculate the electronic component of the nonadiabatic coupling (based on Tavernelli Eq 58)----//
void calculateNACelec(double ** NACmat, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path)
{
  // Output file that lists all transitions contributing to electronic component of NAC
  ofstream outputElecCompFile;
  outputElecCompFile.open((output_path+"output_elec_details").c_str());

  // Shift energy levels to loop over correct number of orbitals, starting form zero
  int HOMOshift = homo - eMin + 1;
  int LUMOshift = lumo - eMin + 1;
  int eMaxShift = eMax - eMin + 1;

  // Create temporary double to calculate coupling to add to sum
  double coupling;

  for (int k = 0; k < exc; k++)
    {
      for (int m = 0; m < exc; m++)
        {
          // Output header for current excitation combination
          outputElecCompFile << "Excitation " << k+1 << "-" << m+1 << endl;
	  outputElecCompFile << setw(8) << "KS Occ" << setw(16) << "KS Unocc 1" << setw(16) << "KS Unocc 2" << setw(16) << "Coeff 1 (t)" << setw(16) << "Coeff 2 (t+Dt)" << setw(16) << "Overlap" << setw(16) << "Coeff 1 (t+Dt)" << setw(16) << "Coeff 2 (t)" << setw(16) << "Overlap" << setw(16) << "LR-TDDFT NAC" << endl;

	  NACmat[k][m] = 0.0;
	  // Calculate first elec NAC term in Tavernelli Eq 58 - XOX
          for (int i = 0; i < HOMOshift; i++)
            {
              for (int j = (LUMOshift-1); j < eMaxShift; j++)
                {
                  for (int p = (LUMOshift-1); p < eMaxShift; p++)
                    {
		      // elec component of TDDFT NAC, multipled by -h/(2step) for units and finite difference
		      coupling = (-1*hbarEvFs/(2*step))*(TDmat1[k][i][j]*TDmat2[m][i][p]*oMat[j][p] - TDmat2[k][i][j]*TDmat1[m][i][p]*oMat[p][j]);
                      NACmat[k][m] = NACmat[k][m] + coupling; 
		      // only output NAC values great than 1E-10
		      if (abs(coupling) > 1.0e-10)
                        {
                          outputElecCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDmat1[k][i][j] << setw(16) << TDmat2[m][i][p] << setw(16) << oMat[j][p] << setw(16) << TDmat2[k][i][j] << setw(16) << TDmat1[m][i][p] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }
	  
	  // Calculate second elec NAC term in Tavernelli Eq 58 YOX
	  for (int i = 0; i < HOMOshift; i++)
            {
              for (int j = (LUMOshift-1); j < eMaxShift; j++)
                {
                  for (int p = (LUMOshift-1); p < eMaxShift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDdeMat1[k][i][j]*TDmat2[m][i][p]*oMat[j][p] - TDdeMat2[k][i][j]*TDmat1[m][i][p]*oMat[p][j]);
		      NACmat[k][m] = NACmat[k][m] + coupling; 
                      if (abs(coupling) > 1.0e-10)
                        {
                          outputElecCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDdeMat1[k][i][j] << setw(16) << TDmat2[m][i][p] << setw(16) << oMat[j][p] << setw(16) << TDdeMat2[k][i][j] << setw(16) << TDmat1[m][i][p] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }

	  // Calculate third elec NAC term in Tavernelli Eq 58 - XOY
          for (int i = 0; i < HOMOshift; i++)
            {
              for (int j = (LUMOshift-1); j < eMaxShift; j++)
                {
                  for (int p = (LUMOshift-1); p < eMaxShift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDmat1[k][i][j]*TDdeMat2[m][i][p]*oMat[j][p] - TDmat2[k][i][j]*TDdeMat1[m][i][p]*oMat[p][j]);
                      NACmat[k][m] = NACmat[k][m] + coupling; 
                      if (abs(coupling) > 1.0e-10) 
                        {
                          outputElecCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDmat1[k][i][j] << setw(16) << TDdeMat2[m][i][p] << setw(16) << oMat[j][p] << setw(16) << TDmat2[k][i][j] << setw(16) << TDdeMat1[m][i][p] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }
	  // Calculate fourth elec NAC term in Tavernelli Eq 58 - YOY
          for (int i = 0; i < HOMOshift; i++)
            {
              for (int j = (LUMOshift-1); j < eMaxShift; j++)
                {
                  for (int p = (LUMOshift-1); p < eMaxShift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDdeMat1[k][i][j]*TDdeMat2[m][i][p]*oMat[j][p] - TDdeMat2[k][i][j]*TDdeMat1[m][i][p]*oMat[p][j]);
                      NACmat[k][m] = NACmat[k][m] + coupling;
                      if (abs(coupling) > 1.0e-10)
                        {
                          outputElecCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDdeMat1[k][i][j] << setw(16) << TDdeMat2[m][i][p] << setw(16) << oMat[j][p] << setw(16) << TDdeMat2[k][i][j] << setw(16) << TDdeMat1[m][i][p] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }
	}
    } 
  outputElecCompFile.close();
}

//----Calculate the hole component of the nonadiabatic coupling (based on Tavernelli Eq 58)----//
void calculateNAChole(double ** NACmat, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path)
{
  // Output file that lists all transitions contributing to hole component of NAC
  ofstream outputHoleCompFile;
  outputHoleCompFile.open((output_path+"output_hole_details").c_str());

  // Shift energy levels to loop over correct number of orbitals, starting form zero
  int HOMOshift = homo - eMin + 1;
  int LUMOshift = lumo - eMin + 1;
  int eMaxShift = eMax - eMin + 1;

  // Create temporary double to calculate coupling to add to sum
  double coupling;

  for (int k = 0; k < exc; k++)
    {
      for (int m = 0; m < exc; m++)
        {
          // Output header for current excitation combination
          outputHoleCompFile << "Excitation " << k+1 << "-" << m+1 << endl;
          outputHoleCompFile << setw(8) << "KS Unocc" << setw(16) << "KS Occ 1" << setw(16) << "KS Occ 2" << setw(16) << "Coeff 1 (t)" << setw(16) << "Coeff 2 (t+Dt)" << setw(16) << "Overlap" << setw(16) << "Coeff 1 (t+Dt)" << setw(16) << "Coeff 2 (t)" << setw(16) << "Overlap" << setw(16) << "LR-TDDFT NAC" << endl;

	  NACmat[k][m] = 0.0;
	  // Calculate first hole NAC term of Taverneli Eq 58 - XOX
          for (int i = (LUMOshift-1); i < eMaxShift; i++)
            {
              for (int j = 0; j < HOMOshift; j++)
                {
                  for (int p = 0; p < HOMOshift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDmat1[k][j][i]*TDmat2[m][p][i]*oMat[j][p] - TDmat2[k][j][i]*TDmat1[m][p][i]*oMat[p][j]);
                      NACmat[k][m] = NACmat[k][m] + coupling;
                      if (abs(coupling) > 1.0e-10)
                        {
                          outputHoleCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDmat1[k][j][i] << setw(16) << TDmat2[m][p][i] << setw(16) << oMat[j][p] << setw(16) << TDmat2[k][j][i] << setw(16) << TDmat1[m][p][i] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }

	  // Calculate second hole NAC term of Taverneli Eq 58 - YOX
          for (int i = (LUMOshift-1); i < eMaxShift; i++) 
            {
              for (int j = 0; j < HOMOshift; j++)
                {
                  for (int p = 0; p < HOMOshift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDdeMat1[k][j][i]*TDmat2[m][p][i]*oMat[j][p] - TDdeMat2[k][j][i]*TDmat1[m][p][i]*oMat[p][j]);
		      NACmat[k][m] = NACmat[k][m] + coupling;
                      if (abs(coupling) >  1.0e-10) // only output NAC values > 1E-10
                        {
                          outputHoleCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDdeMat1[k][j][i] << setw(16) << TDmat2[m][p][i] << setw(16) << oMat[j][p] << setw(16) << TDdeMat2[k][j][i] << setw(16) << TDmat1[m][p][i] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }

	  // Calculate third hole NAC term of Taverneli Eq 58 - XOY
          for (int i = (LUMOshift-1); i < eMaxShift; i++) 
            {
              for (int j = 0; j < HOMOshift; j++)
                {
                  for (int p = 0; p < HOMOshift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDmat1[k][j][i]*TDdeMat2[m][p][i]*oMat[j][p] - TDmat2[k][j][i]*TDdeMat1[m][p][i]*oMat[p][j]);
		      NACmat[k][m] = NACmat[k][m] + coupling;
                      if (abs(coupling) > 1.0e-10)
                        {
                          outputHoleCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDmat1[k][j][i] << setw(16) << TDdeMat2[m][p][i] << setw(16) << oMat[j][p] << setw(16) << TDmat2[k][j][i] << setw(16) << TDdeMat1[m][p][i] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }

	  // Calculate fourth hole NAC term of Taverneli Eq 58 - YOY
          for (int i = (LUMOshift-1); i < eMaxShift; i++) 
            {
              for (int j = 0; j < HOMOshift; j++)
                {
                  for (int p = 0; p < HOMOshift; p++)
                    {
		      coupling = (-1*hbarEvFs/(2*step))*(TDdeMat1[k][j][i]*TDdeMat2[m][p][i]*oMat[j][p] - TDdeMat2[k][j][i]*TDdeMat1[m][p][i]*oMat[p][j]);
		      NACmat[k][m] = NACmat[k][m] + coupling;
                      if (abs(coupling) > 1.0e-10)
                        {
                          outputHoleCompFile << setw(8) << (i+eMin) << setw(16) << (j+eMin) << setw(16) << (p+eMin) << scientific << setprecision(4) << setw(16) << TDdeMat1[k][j][i] << setw(16) << TDdeMat2[m][p][i] << setw(16) << oMat[j][p] << setw(16) << TDdeMat2[k][j][i] << setw(16) << TDdeMat1[m][p][i] << setw(16) << oMat[p][j] << setw(16) << coupling << endl;
                        }
                    }
                }
            }
	}
    }
  outputHoleCompFile.close();
}

//----Calculate the nonadiabatic coupling between the ground state and all excited states (based on Tavernelli Eq 47)----//
void calculateNACground(double * NACvec, double *** TDmat1, double *** TDmat2, double *** TDdeMat1, double *** TDdeMat2, double ** oMat, int exc, int homo, int lumo, int eMin, int eMax, double step, string output_path)
{
  // Output file that lists all transitions contributing to ground-excited NAC
  ofstream outputGroundCompFile;
  outputGroundCompFile.open((output_path+"output_exc_gs_details").c_str());

  // Shift energy levels to loop over correct number of orbitals, starting form zero
  int HOMOshift = homo - eMin + 1;
  int LUMOshift = lumo - eMin + 1;
  int eMaxShift = eMax - eMin + 1;

  // Create temporary double to calculate coupling to add to sum
  double coupling;

  for (int i = 0; i < exc; i++)
    {
      outputGroundCompFile << "Excitation " << i+1 << "-Ground State" << endl;
      outputGroundCompFile << setw(8) << "KS Unocc" << setw(16) << "KS Occ" << setw(16) << "Coeff (t+Dt)" << setw(16) << "Overlap" << setw(16) << "Coeff (t)" << setw(16) << "Overlap" << setw(16) << "LR-TDDFT NAC" << endl;
      NACvec[i] = 0.0;
      for (int j = (LUMOshift-1); j < eMaxShift; j++)
        {
          for (int k = 0; k < HOMOshift; k++)
            {
	      coupling = (-1*hbarEvFs/(2*step))*(TDmat2[i][k][j]*oMat[k][j] - TDmat1[i][k][j]*oMat[j][k]);
              NACvec[i] = NACvec[i] + coupling;
              if (abs(coupling) > 1.0e-10)
                {
                  outputGroundCompFile << setw(8) << (j+eMin) << setw(16) << (k+eMin) << setw(16) << scientific << setprecision(6) << setw(16) << TDmat2[i][k][j] << setw(16) << oMat[k][j] << setw(16) << TDmat1[i][k][j] << setw(16) << oMat[j][k] << setw(16) << coupling  << endl;;
                }
            }
        }
      for (int j = (LUMOshift-1); j < eMaxShift; j++)
        {
          for (int k = 0; k < HOMOshift; k++)
            {
	      coupling = (-1*hbarEvFs/(2*step))*(TDdeMat2[i][k][j]*oMat[k][j] - TDdeMat1[i][k][j]*oMat[j][k]);
              NACvec[i] = NACvec[i] + coupling;
              if (abs(coupling) > 1.0e-10)
                {
                  outputGroundCompFile << setw(8) << (j+eMin) << setw(16) << (k+eMin) << setw(16) << scientific << setprecision(6) << setw(16) << TDdeMat2[i][k][j] << setw(16) << oMat[k][j] << setw(16) << TDdeMat1[i][k][j] << setw(16) << oMat[j][k] << setw(16) << coupling << endl;;
                }
            }
        }
    }
  outputGroundCompFile.close();
}

//----Construct 1) total 2) electronic and 3) hole nonadiabatic coupling matrix for one time step pair----//
void constructNACmatrix(double *** NACmatTotalTimeStep, double *** NACmatElecTimeStep, double *** NACmatHoleTimeStep, double ** NACmatElec, double ** NACmatHole, double * NACvecGr, int tStep, int exc, string output_path)
{
  // Output file with total nonadiabatic coupling matrix for each time step
  ofstream outputNACallFile;
  outputNACallFile.open((output_path+"output_NAC_total").c_str());
  ofstream outputNACelecFile;
  outputNACelecFile.open((output_path+"output_NAC_elec").c_str());
  ofstream outputNACholeFile;
  outputNACholeFile.open((output_path+"output_NAC_hole").c_str());
  

  // Place ground state coupling in first row and column of each matrix
  NACmatTotalTimeStep[tStep][0][0] = 0.0;
  NACmatElecTimeStep[tStep][0][0] = 0.0;
  NACmatHoleTimeStep[tStep][0][0] = 0.0;
  for (int i = 0; i < exc; i++)
    {
      NACmatTotalTimeStep[tStep][0][i+1] = NACvecGr[i];
      NACmatElecTimeStep[tStep][0][i+1] = NACvecGr[i];
      NACmatHoleTimeStep[tStep][0][i+1] = NACvecGr[i];
      NACmatTotalTimeStep[tStep][i+1][0] = -1*NACvecGr[i];
      NACmatElecTimeStep[tStep][i+1][0] = -1*NACvecGr[i];
      NACmatHoleTimeStep[tStep][i+1][0] = -1*NACvecGr[i];
    }
  // Fill total, electronic, and hole nonadiabatic coupling matrices
  for (int i = 0; i < exc; i++)
    {
      for (int j = 0; j < exc; j++)
	{
	  NACmatTotalTimeStep[tStep][i+1][j+1] = NACmatElec[i][j] - NACmatHole[i][j];
	  NACmatElecTimeStep[tStep][i+1][j+1] = NACmatElec[i][j];
	  NACmatHoleTimeStep[tStep][i+1][j+1] = NACmatHole[i][j];
	}
    }

  //Output to file
  for (int i = 0; i < (exc+1); i++)
    {
      for (int j = 0; j < (exc+1); j++)
	{
	  //Output to file
	  outputNACallFile << scientific << setprecision(6) << NACmatTotalTimeStep[tStep][i][j] <<"\t";
	  outputNACelecFile << scientific << setprecision(6) << NACmatElecTimeStep[tStep][i][j] << "\t";
	  outputNACholeFile << scientific << setprecision(6) << NACmatHoleTimeStep[tStep][i][j] << "\t";
	}
      outputNACallFile << endl;
      outputNACelecFile << endl;
      outputNACholeFile << endl;
    }
  outputNACallFile.close();
  outputNACelecFile.close();
  outputNACholeFile.close();
}

//----Takes thermal averaging of nonadiabatic coupling matrix over all MD time steps----//
void thermalAvgNACmatrix(double ** avgNACtotal, double ** avgNACelec, double ** avgNAChole, double *** NACmatTotalTimeStep, double *** NACmatElecTimeStep, double *** NACmatHoleTimeStep, int exc, int totalSteps)
{
  // Output files for final, averaged, squared nonadiabatic couplings for total, electron, and hole
  ofstream outputAvgNACallFile;
  outputAvgNACallFile.open("final_NAC_matrix_total");
  ofstream outputAvgNACelecFile;
  outputAvgNACelecFile.open("final_NAC_matrix_elec");
  ofstream outputAvgNACholeFile;
  outputAvgNACholeFile.open("final_NAC_matrix_hole");

  ofstream NACtimeseriesAllFile;
  NACtimeseriesAllFile.open("NAC_timeseries_total");
  ofstream NACtimeseriesElecFile;
  NACtimeseriesElecFile.open("NAC_timeseries_elec");
  ofstream NACtimeseriesHoleFile;
  NACtimeseriesHoleFile.open("NAC_timeseries_hole");


  
  // Loop through each excitation combination, square NAC value at given time step, then average over all time steps
  for (int i = 0; i < (exc+1); i++)
    {
      for (int j = 0; j < (exc+1); j++)
	{
	  avgNACtotal[i][j] = 0.0;
	  avgNACelec[i][j] = 0.0;
	  avgNAChole[i][j] = 0.0;
	  for (int k = 0; k < totalSteps; k++)
	    {
	      avgNACtotal[i][j] = avgNACtotal[i][j] + NACmatTotalTimeStep[k][i][j]*NACmatTotalTimeStep[k][i][j];
	      avgNACelec[i][j] = avgNACelec[i][j] + NACmatElecTimeStep[k][i][j]* NACmatElecTimeStep[k][i][j];
	      avgNAChole[i][j] = avgNAChole[i][j] + NACmatHoleTimeStep[k][i][j]*NACmatHoleTimeStep[k][i][j];

	      // Output current squared average to time series file to check for convergence w/r/t MD time steps
	      double dk = double(k);
	      NACtimeseriesAllFile << scientific << setprecision(6) << setw(16) << avgNACtotal[i][j]/(dk+1);
	    }
	  avgNACtotal[i][j] = avgNACtotal[i][j]/totalSteps;
	  avgNACelec[i][j] = avgNACelec[i][j]/totalSteps;
	  avgNAChole[i][j] = avgNAChole[i][j]/totalSteps;
	  
	  // Output final, thermally averaged nonadiabatic coupling matrix to input into relaxation dynamics
	  outputAvgNACallFile << scientific << setprecision(6) << avgNACtotal[i][j] << "\t";
	  outputAvgNACelecFile << scientific << setprecision(6) << avgNACelec[i][j] << "\t";
	  outputAvgNACholeFile << scientific << setprecision(6) << avgNAChole[i][j] << "\t";

	  NACtimeseriesAllFile << endl;
	}
      outputAvgNACallFile << endl;
      outputAvgNACelecFile << endl;
      outputAvgNACholeFile << endl;
    }  
}


//----Function to order DFT states from time step to time step based on KS overlap matrix
void orderDFTstates(double ** oMat, int ** switchMat, double *** TDmat1, double *** TDdeMat1, double *** TDmat2, double *** TDdeMat2, int eLevels, int exc, int tStep, int eMin, string output_path)
{
  ofstream outputKSflip;
  ofstream outputOvlapOrdered;
  ofstream outputTDcoeffOrdered1;
  ofstream outputTDdeCoeffOrdered1;
  ofstream outputTDcoeffOrdered2;
  ofstream outputTDdeCoeffOrdered2;

  outputKSflip.open((output_path+"output_KS_flip").c_str());
  outputOvlapOrdered.open((output_path+"output_ovlap_ordered").c_str());
  outputTDcoeffOrdered1.open((output_path+"output_coeff_exc_1_KSordered").c_str());
  outputTDdeCoeffOrdered1.open((output_path+"output_coeff_de_exc_1_KSordered").c_str());
  outputTDcoeffOrdered2.open((output_path+"output_coeff_exc_2_KSordered").c_str());
  outputTDdeCoeffOrdered2.open((output_path+"output_coeff_de_exc_2_KSordered").c_str());

  // First, implement all swaps from previous time steps to keep KS overlap matrix up to date by swapping rows (ith time step)
  for (int i = 0; i < tStep; i++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  if (switchMat[i][j] != 0 )
	    {
	      // Flip rows based on previous switches
	      for (int k = 0; k < eLevels; k++)
		{
		  swapDouble(&oMat[switchMat[i][j]-1][k],&oMat[j][k]);
		  // Flip TD coefficients based on previous wtiches
		  for (int l = 0; l < exc; l++)
		    {
		      swapDouble(&TDmat1[l][switchMat[i][j]-1][k],&TDmat1[l][j][k]);
		      swapDouble(&TDdeMat1[l][switchMat[i][j]-1][k],&TDdeMat1[l][j][k]);
		      swapDouble(&TDmat1[l][k][switchMat[i][j]-1],&TDmat1[l][k][j]);
		      swapDouble(&TDdeMat1[l][k][switchMat[i][j]-1],&TDdeMat1[l][k][j]);
		    }
		}
	      // Flip columns based on previous witches
	      for (int k = 0; k < eLevels; k++)
		{
		  swapDouble(&oMat[k][j], &oMat[k][switchMat[i][j]-1]);
		}
	    }
	}
    }

  // Next, check diagonal elements of current KS overlap matrix to determine switching for current time step
  // Loop through all column elements for each row to determine maximum element
  for (int i = 0; i < eLevels; i++)
    {
      int maxIndex = 0;
      double max = abs(oMat[i][0]);
      for (int j = 1; j < eLevels; j++)
	{
	  if (abs(oMat[i][j]) > max)
	    {
	      max = abs(oMat[i][j]);
	      maxIndex = j;
	    }
	}

      // If maximum value not on diagonal, then swap column with column containing max value to make diagonal value maximum
      double temp = 0.0;
      if ((abs(oMat[i][i]) != max) && (max > 0.3))
	{
	  // Save index of column to be switched in switching matrix to apply to future time steps
	  // Note that index saved with (+1) so that 0 can indicate no switching occurred!
	  switchMat[tStep][i] = maxIndex + 1;
	  for (int k = 0; k < eLevels; k++)
	    {
	      swapDouble(&oMat[k][i], &oMat[k][maxIndex]);
	    }
	  // Output columns flipped
	  outputKSflip << "Flipped " << eMin+i << " and " << eMin+maxIndex << " columns " << endl;
      	}
    }
  
  // Output new overlap matrix to output file in coupling subfolder
  for (int i = 0; i < eLevels; i++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  outputOvlapOrdered << scientific << setprecision(6) << oMat[i][j] << " "; 
	}
      outputOvlapOrdered << endl;
    }


  // Update second LR-TDDFT coefficients (i+1th time step) to taken into account changes in ordering
  // For first coefficient matrix, implement swaps from all previous time steps to keep up to date
  for (int i = 0; i < (tStep+1); i++)
    {
      for (int j = 0; j < eLevels; j++)
	{
	  if (switchMat[i][j] != 0 )
	    {
	      for (int k = 0; k < eLevels; k++)
		{
		  for (int l = 0; l < exc; l++)
		    {
		      swapDouble(&TDmat2[l][switchMat[i][j]-1][k],&TDmat2[l][j][k]);
		      swapDouble(&TDdeMat2[l][switchMat[i][j]-1][k],&TDdeMat2[l][j][k]);
		      swapDouble(&TDmat2[l][k][switchMat[i][j]-1],&TDmat2[l][k][j]);
		      swapDouble(&TDdeMat2[l][k][switchMat[i][j]-1],&TDdeMat2[l][k][j]);
		    }
		}
	    }
	}
    }  

  // Output ordered LR-TDDFT excitation coefficient matrices to each coupling folder
  for (int i = 0; i < exc; i++)
    {
      outputTDcoeffOrdered1 << "Excitation " << i+1 << ":" << endl;
      outputTDdeCoeffOrdered1 << "Excitation " << i+1 << ":" << endl;
      outputTDcoeffOrdered2 << "Excitation " << i+1 << ":" << endl;
      outputTDdeCoeffOrdered2 << "Excitation " << i+1 << ":" << endl;

      for (int j = 0; j < eLevels; j++)
	{
	  for (int k = 0; k < eLevels; k++)
	    {
	      outputTDcoeffOrdered1 << scientific << setprecision(6) << TDmat1[i][j][k] << '\t';
	      outputTDcoeffOrdered2 << scientific << setprecision(6) << TDmat2[i][j][k] << '\t';
	      outputTDdeCoeffOrdered1 << scientific << setprecision(6) << TDdeMat1[i][j][k] << '\t';
	      outputTDdeCoeffOrdered2 << scientific << setprecision(6) << TDdeMat2[i][j][k] << '\t';
	    }
	  outputTDcoeffOrdered1 << endl;
	  outputTDcoeffOrdered2 << endl;
	  outputTDdeCoeffOrdered1 << endl;
	  outputTDdeCoeffOrdered2 << endl;
	}
    }

  outputKSflip.close();
  outputOvlapOrdered.close();
  outputTDcoeffOrdered1.close();
  outputTDcoeffOrdered2.close();
  outputTDdeCoeffOrdered1.close();
  outputTDdeCoeffOrdered2.close();

}

//----Function to order LR-TDDFT states from time step to time step based on LR-TDDFT overlap matrix
void orderTDDFTstates(double ** TDexcEn1, double ** TDexcEn2, double ** TDexcEnOrd1, double ** TDexcEnOrd2, double ** TDoverlap, int ** switchMat, double *** TDmat1, double *** TDdeMat1, double *** TDmat2, double *** TDdeMat2, int eLevels, int exc, int tStep, int eMin, string output_path)
{
  ofstream outputTDoverlap;
  ofstream outputTDoverlapOrdered;
  ofstream outputTDflip;
  ofstream outputTDcoeffTDOrdered1;
  ofstream outputTDdeCoeffTDOrdered1;
  ofstream outputTDcoeffTDOrdered2;
  ofstream outputTDdeCoeffTDOrdered2;

  outputTDoverlap.open((output_path+"output_TD_overlap").c_str());
  outputTDflip.open((output_path+"output_TD_flip").c_str());
  outputTDoverlapOrdered.open((output_path+"output_TD_overlap_ordered").c_str());
  outputTDcoeffTDOrdered1.open((output_path+"output_coeff_exc_1_TDordered").c_str());
  outputTDdeCoeffTDOrdered1.open((output_path+"output_coeff_de_exc_1_TDordered").c_str());
  outputTDcoeffTDOrdered2.open((output_path+"output_coeff_exc_2_TDordered").c_str());
  outputTDdeCoeffTDOrdered2.open((output_path+"output_coeff_de_exc_2_TDordered").c_str());

  // Fill ordered TD excitation arrays with original to start:
  for (int j = 0; j < exc; j++)
    {
      TDexcEnOrd1[j][0] = TDexcEn1[j][0];
      TDexcEnOrd1[j][1] = TDexcEn1[j][1];
      TDexcEnOrd2[j][0] = TDexcEn2[j][0];
      TDexcEnOrd2[j][1] = TDexcEn2[j][1];
    }

  // Calculate LR-TDDFT overlap matrix between current, adjacent time steps (i and i+1)
  for (int j = 0; j < exc; j++)
    {
      for (int k = 0; k < exc; k++)
	{
	  TDoverlap[j][k] = 0.0;
	  for (int l = 0; l < eLevels; l++)
	    {
	      for (int m = 0; m < eLevels; m++)
		{
		  TDoverlap[j][k] = TDoverlap[j][k] + TDmat1[j][l][m]*TDmat2[k][l][m] - TDdeMat1[j][l][m]*TDdeMat2[k][l][m];
		}
	    }
	  outputTDoverlap << scientific << setprecision(6) << TDoverlap[j][k] << '\t';	  
	}
      outputTDoverlap << endl;
    }

  // Then, implement all swaps from previous time steps to keep LR-TDDFT overlap matrix up to date by swapping rows (ith time step)
  for (int i = 0; i < tStep; i++)
    {
      for (int j = 0; j < exc; j++)
	{
	  if (switchMat[i][j] != 0 )
	    {
	      // Swap LR-TDDFT excitatione nergy ordering
	      swapDouble(&TDexcEnOrd1[switchMat[i][j]-1][0],&TDexcEnOrd1[j][0]);
	      swapDouble(&TDexcEnOrd1[switchMat[i][j]-1][1],&TDexcEnOrd1[j][1]);
	      // Swap LR-TDDFT coefficient ordering
	      for (int k = 0; k < eLevels; k++)
		{
		  for (int l = 0; l < eLevels; l++)
		    {
		      swapDouble(&TDmat1[switchMat[i][j]-1][k][l],&TDmat1[j][k][l]);
		      swapDouble(&TDdeMat1[switchMat[i][j]-1][k][l],&TDdeMat1[j][k][l]);
		    }
		}
	    }
	}
    }

  // Next, check diagonal elements of current TD overlap matrix to determine switching for current time step
  // Loop through all column elements for each row to determine maximum element
  for (int i = 0; i < exc; i++)
    {
      int maxIndex = 0;
      double max = abs(TDoverlap[i][0]);
      for (int j = 1; j < exc; j++)
	{
	  if (abs(TDoverlap[i][j]) > max)
	    {
	      max = abs(TDoverlap[i][j]);
	      maxIndex = j;
	    }
	}

      // If maximum value not on diagonal, then swap column witih column containing max value to make diagonal value maximum
      double temp = 0.0;
      if ((abs(TDoverlap[i][i]) != max) && (max > 0.3))
	{
	  // Save index of column to be switched in switching matrix to apply to future time steps
	  // Note that index saved with (+1) so that 0 can indicate no switching occurred!
	  switchMat[tStep][i] = maxIndex + 1;
	  for (int k = 0; k < exc; k++)
	    {
	      swapDouble(&TDoverlap[k][i], &TDoverlap[k][maxIndex]);
	    }
	  // Output columns flipped
	  outputTDflip << "Flipped " << i+1 << " and " << maxIndex+1 << " columns " << endl;
      	}
    }

  // Update second LR-TDDFT coefficients (i+1th time step) to taken into account changes in ordering
  // For first coefficient matrix, implement swaps from all previous time steps to keep up to date
  for (int i = 0; i < (tStep+1); i++)
    {
      for (int j = 0; j < exc; j++)
	{
	  if (switchMat[i][j] != 0 )
	    {
	      // Swap LR-TDDFT excitatione nergy ordering
	      swapDouble(&TDexcEnOrd2[switchMat[i][j]-1][0],&TDexcEnOrd2[j][0]);
	      swapDouble(&TDexcEnOrd2[switchMat[i][j]-1][1],&TDexcEnOrd2[j][1]);
	      // Swap LR-TDDFT coefficient ordering
	      for (int k = 0; k < eLevels; k++)
		{
		  for (int l = 0; l < eLevels; l++)
		    {
		      swapDouble(&TDmat2[switchMat[i][j]-1][k][l],&TDmat2[j][k][l]);
		      swapDouble(&TDdeMat2[switchMat[i][j]-1][k][l],&TDdeMat2[j][k][l]);
		    }
		}
	    }
	}
    }  

  // Recalculate LR-TDDFT overlap matrix and output to new file
  for (int j = 0; j < exc; j++)
    {
      for (int k = 0; k < exc; k++)
	{
	  TDoverlap[j][k] = 0.0;
	  for (int l = 0; l < eLevels; l++)
	    {
	      for (int m = 0; m < eLevels; m++)
		{
		  TDoverlap[j][k] = TDoverlap[j][k] + TDmat1[j][l][m]*TDmat2[k][l][m] - TDdeMat1[j][l][m]*TDdeMat2[k][l][m];
		}
	    }
	  outputTDoverlapOrdered << scientific << setprecision(6) << TDoverlap[j][k] << '\t';	  
	}
      outputTDoverlapOrdered << endl;
    }

  // Output ordered LR-TDDFT excitation coefficient matrices to each coupling folder
  for (int i = 0; i < exc; i++)
    {
      outputTDcoeffTDOrdered1 << "Excitation " << i+1 << ":" << endl;
      outputTDdeCoeffTDOrdered1 << "Excitation " << i+1 << ":" << endl;
      outputTDcoeffTDOrdered2 << "Excitation " << i+1 << ":" << endl;
      outputTDdeCoeffTDOrdered2 << "Excitation " << i+1 << ":" << endl;

      for (int j = 0; j < eLevels; j++)
	{
	  for (int k = 0; k < eLevels; k++)
	    {
	      outputTDcoeffTDOrdered1 << scientific << setprecision(6) << TDmat1[i][j][k] << '\t';
	      outputTDcoeffTDOrdered2 << scientific << setprecision(6) << TDmat2[i][j][k] << '\t';
	      outputTDdeCoeffTDOrdered1 << scientific << setprecision(6) << TDdeMat1[i][j][k] << '\t';
	      outputTDdeCoeffTDOrdered2 << scientific << setprecision(6) << TDdeMat2[i][j][k] << '\t';
	    }
	  outputTDcoeffTDOrdered1 << endl;
	  outputTDcoeffTDOrdered2 << endl;
	  outputTDdeCoeffTDOrdered1 << endl;
	  outputTDdeCoeffTDOrdered2 << endl;
	}
    }

  // Close files
  outputTDoverlap.close();
  outputTDoverlapOrdered.close();
  outputTDflip.close();
  outputTDcoeffTDOrdered1.close();
  outputTDcoeffTDOrdered2.close();
  outputTDdeCoeffTDOrdered1.close();
  outputTDdeCoeffTDOrdered2.close();

  // Output raw and ordered excitation energies to screen for comparison
  cout << "Pre- and Post-Ordering of LR-TDDFT Excited States:" << endl;
  cout << "Time Step " << tStep+1 << ":" << endl;
  for (int j = 0; j < exc; j++)
    {
      cout << "Excited State  " << j+1 << ":  " << fixed << setprecision(4) << TDexcEn1[j][0] << " nm " << TDexcEn1[j][1] << " Ordered: " << TDexcEnOrd1[j][0] << " nm " << TDexcEnOrd1[j][1] << endl;
    }
  cout << "Time Step " << tStep+2 << ":" << endl;
  for (int j = 0; j < exc; j++)
    {
      cout << "Excited State  " << j+1 << ":  " << fixed << setprecision(4) << TDexcEn2[j][0] << " nm " << TDexcEn2[j][1] << " Ordered: " << TDexcEnOrd2[j][0] << " nm " << TDexcEnOrd2[j][1] << endl;
    }

}

//----Function to swap two doubles----//
void swapDouble(double * first, double * second)
{
  double temp = *first;
  *first = *second;
  *second = temp;
}

