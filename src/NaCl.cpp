#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "tsin.h"
#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"
#include "GetInput.h"

#include "erf.h"

//#include "potfit.h"

// excess electron includes
#include "Potential.h"
#include "DVR.h"

using namespace std;

//////////////////////////////
//
//  IonTypes: 
//
const int Na_ion = 1;
const int Cl_ion = 2;


// Functions defined here
int GetIonCoordinates(int nLines, char **coorlines, double *IonPos, int *IonType);
void PrintIonCoordinates(int nIons, double *IonPos, int *IonType, double scale);



///////////////////////////////////////////////
//
//  single point for NaCl
//  has so far just the VEA
//
void NaCl_sp(TSIN Input)
{

  printf("-----------------\nNaCl single point\n-----------------\n");

  // read the number and coordinates of the Na and Cl ions from Input
  char **coorlines;
  int nLines = Input.GetGroup("NaCl", &coorlines);
  double *IonPos = new double[3*nLines];  // coordinates for Na and Cl
  int *IonType = new int[nLines];
  int nIons = GetIonCoordinates(nLines, coorlines, IonPos, IonType);
  if (nIons < nLines)
    printf("Found %i ions reading %i NaCl-lines.\n", nIons, nLines);
  printf("Cluster Geometry [Angs]\n");
  PrintIonCoordinates(nIons, IonPos, IonType, Bohr2Angs);

  //
  //  compute the energy of the neutral cluster
  //
  
  // 
  //  Now, dealing with the excess electron is more involved
  //

  // set up the Potential for the excess electron
  // for NaCl, Velectron is independent of the model for the ion/ion interactions
  double *Charges = new double[nIons];
  int *iqs = new int[nIons];
  for (int i = 0; i < nIons; ++i) {
    iqs[i] = i;
    switch (IonType[i])
      {
      case Na_ion:
	Charges[i] = 1.0;
	break;
      case Cl_ion:
	Charges[i] = -1.0;
	break;
      default:
	printf("This should never happen\n");
	exit(42);
      }
  }
  class Potential Vel;
  double PotPara[100];
  int PotFlag = Input.GetInt("ElectronPotential", "Potential", 101);
  int PotVerb = Input.GetInt("ElectronPotential", "Verbose", 1);
  Vel.SetVerbose(PotVerb);
  switch (PotFlag)
    {
    case 101:
      PotPara[0] = Input.GetDouble("ElectronPotential", "rNa", 2.9);
      PotPara[1] = Input.GetDouble("ElectronPotential", "rCl", 1.0);
      Vel.SetupBloomfield(nIons, IonPos, nIons, Charges, iqs, PotPara);
      break;
    default:
      cout << "NaCl Single Point: unknown potential for the excess electron: " << PotFlag << "\n";
      exit(1);
    }




  // set up a grid (a DVR) for the wavefunction of the electron
  int ngrid[3] = {0,0,0}; 
  Input.GetIntArray("GridDef", "NoOfGridPoints", ngrid, 3);
  int gridverbose = Input.GetInt("GridDef", "Verbose", 0);
  int dvrtype = Input.GetInt("GridDef", "DVRType", 3);      // 1 HO,  other = Sine 
  double gpara[3] = {0,0,0};  // this is for the grid length or the oscillator freq
  switch (dvrtype)
    {
    case 1:
      if (gridverbose > 0)
	printf("Harmonic oscillator DVR\n");
      Input.GetDoubleArray("GridDef", "Omega0", gpara, 3);
      break;
    default:
      if (gridverbose > 0)
	printf("Particle-in-a-box (sine) DVR\n");
      Input.GetDoubleArray("GridDef", "Length", gpara, 3);
      break;
    }

  class DVR Hel;
  Hel.SetupDVR(ngrid, dvrtype, 1, gpara, gridverbose);
  Hel.ComputePotential(Vel);




  // now get the input for the diagonalizer
  int nStates = Input.GetInt("Diag", "NEigenValues", 1);
  int DiagFlag = Input.GetInt("Diag", "Method", 1);     // 1=Lanczos  2=Davidson  3= 0th-order Jacobi-Davidson
  if (nStates > 1 && DiagFlag != 1) {
    cout << "So far only Lanczos can compute multiple eigen values\n";
    DiagFlag = 1;
  }
  int maxSub = Input.GetInt("Diag", "maxSubspace", 20);   // max no of micro-iterations
  int maxIter = Input.GetInt("Diag", "MaxIter", 10);      // max no of macro-iterations 
  int ptol = Input.GetDouble("Diag", "pTol", 5);          // tolerance = 10^-pTol 
  int istartvec = Input.GetInt("Diag", "StartVector", 1); // this is different for Lanczos + Davidson and needs work
  int diagverbose = Input.GetInt("Diag", "Verbose", 1);
  Hel.SetVerbose(diagverbose);  // resets verbose from gridverbose to diagverbose

  double *energies = new double[nStates];
  Hel.DiagonalizeSetup(nStates, DiagFlag, maxSub, maxIter, ptol);
  Hel.Diagonalize(istartvec, energies);
  cout << "Eel = " << energies[0] << " Hartree = " << energies[0]*AU2MEV << " meV\n";


  // plot the wavefunction
  int *NucCharge = new int[nIons];
  for (int k = 0; k < nIons; ++k) {
    switch (IonType[k])
      {
      case Na_ion:
	NucCharge[k] = 11;
	break;
      case Cl_ion:
	NucCharge[k] = 17;
	break;
      default:
	printf("This should never happen\n");
	exit(42);
      }
  }
  Hel.WriteCubeFile(1, "WaveFn.cub", nIons, NucCharge, IonPos, 1);
  Hel.WriteCuts();



  cout << "All done.\n"; cout.flush();

}





////////////////////////////////////////////////////////////////////////////////////
//
//   takes nIons lines and reads the coordinates and the ion type
//
//
//
//   ignores empty lines and lines starting with #
//   returns the number ions actually found
//
int GetIonCoordinates(int nLines, char **coorlines, double *IonPos, int *IonType)
{
  int nIons = 0;
  char atomsymbol[5];
  double x = 0, y = 0, z = 0;
  char *inputline = 0;

  for (int iLine = 0; iLine < nLines; ++iLine) {
    inputline = coorlines[iLine];
    if (*inputline != '#' && strlen(inputline) > 1) {
      if (sscanf(inputline," %s %lf %lf %lf", atomsymbol, &x, &y, &z) != 4) {
	cout << "Error reading input group NaCl, line" << iLine << ":" << inputline << "\n";
      exit(1);
      }
      IonPos[3*nIons+0] = x * Angs2Bohr;
      IonPos[3*nIons+1] = y * Angs2Bohr;
      IonPos[3*nIons+2] = z * Angs2Bohr;
      switch (atomsymbol[0])
	{
	case 'N':
	  IonType[nIons] = Na_ion;
	  break;
	case 'C':
	  IonType[nIons] = Cl_ion;
	  break;
	default:
	  cout << "Error reading input group NaCl, line" << iLine << ":" << inputline << "\n";
	  exit(1);
	}
      nIons ++;
    }
  }
  return nIons;
}


void PrintIonCoordinates(int nIons, double *IonPos, int *IonType, double scale)
{
  for (int i = 0; i < nIons; ++i) {
    switch (IonType[i])
      {
      case Na_ion:
	printf("Na ");
	break;
      case Cl_ion:
	printf("Cl ");
	break;
      default:
	cout << "Error in PrintIonCoordinates: unknown atomtype for atom " << i << ": " << IonType[i] << "\n";
	exit(1);
      }
    printf("  %12.8f  %12.8f  %12.8f\n", IonPos[3*i]*scale, IonPos[3*i+1]*scale, IonPos[3*i+2]*scale);
  }
}






/////////////////////////////////////////////////////////////////////////////////////////////
//
//  Bloomfield's potential (JCP 102, 4965 (1995)
//
//  charges damped with an error function
//  parameter[0] rNa  (2.9 in the paper)
//  parameter[1] rCl  (1.0 in the paper)
//
void Potential::SetupBloomfield(int nr, const double *r, int nq, const double *q, const int *iq, const double *PotPara)
{ 
  cout << "Enter SetupBl\n";

  PotFlags[0] = 101;
  CationDamping = 1.0 / PotPara[0];
  AnionDamping = 1.0 / PotPara[1];

  if (verbose > 0) {
    cout << "------------------------\n";
    cout << "Bloomfield's potential\n";
    cout << "PotPara 0: " << PotPara[0] << " (rNa) [Bohr]\n";
    cout << "PotPara 1: " << PotPara[1] << " (rCl) [Bohr]\n";
    cout << "------------------------\n";
  }

  SetSites(nr, r);
  SetCharges(nq, q, iq);

  SetupMinMax();
  cout << "Leaving SetupBl\n";

}


double Potential::EvaluateBloomfield(const double *x)
{
  // compute a list of distances of all sites to the point r=(x,y,z)
  //for (int i = 0; i < nSites; ++i) {
  //  Rx[i] = x[0] - Site[3*i+0];
  //  Ry[i] = x[1] - Site[3*i+1];
  //  Rz[i] = x[2] - Site[3*i+2];
  //  R2[i] = Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i] + 1e-15; 
  //  R[i] = sqrt(R2[i]);
  //  Rminus3[i] = 1.0 / (R2[i] * R[i]);
  //}

  //  potential due to damped point charges
  double Vel = 0;
  double damp = 0;
  for (int i = 0; i < nCharges; ++i) {
    int iqs = ChargeSite[i];
    if (Charge[i] > 0)
      damp = CationDamping;
    else
      damp = AnionDamping;
    Vel += -Charge[i] * erf(damp*R[iqs]) / R[iqs];
  }


  if (verbose > 0)
    CheckMinMax(&Vel);

  if (verbose > 15)
    printf("%10.6f %10.6f %10.6f %14.6e\n", 
	   x[0], x[1], x[2], Vel);

  return Vel;
}
