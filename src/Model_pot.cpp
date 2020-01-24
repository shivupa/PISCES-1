#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "tsin.h"
#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"
#include "GetInput.h"
#include "lapackblas.h"
//#include "potfit.h"

// excess electron includes
#include "Potential.h"
#include "DVR.h"




using namespace std;

//////////////////////////////
//
//  IonTypes: 
//







///////////////////////////////////////////////
///
///  single point for C60-electron
///  has so far just the VEA
///
void SphericalPotential(TSIN Input)
{

  printf("-----------------\nSpherical well in bohrs-----------------\n");

  class Potential Vel;
  double PotPara[100];
  int PotVerb = Input.GetInt("SphericalPotential", "Verbose", 1);
  int PotFlag = Input.GetInt("SphericalPotential", "Potential", 1);
  Vel.SetVerbose(PotVerb);
  switch (PotFlag)
    {
    case 1023:
      PotPara[0] = Input.GetDouble("SphericalPotential", "Depth", 1.0) ;
      PotPara[1] = Input.GetDouble("SphericalPotential", "Width", 1.0) ;
      PotPara[2] = Input.GetDouble("SphericalPotential", "ForceConstant", 1.0) ;
      Vel.SetupSphericalPotential(PotPara);
      break;
    default:
      cout << " SetupSphericalSingle Point: unknown potential for the excess electron: " << PotFlag << "\n";
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
  int nStates = Input.GetInt("Diag", "nStates", 1);
  int DiagFlag = Input.GetInt("Diag", "Method", 1);     // 1=Lanczos  2=Davidson  3= 0th-order Jacobi-Davidson
  int maxSub = Input.GetInt("Diag", "maxSubspace", 20);   // max no of micro-iterations
  int MaxIter = Input.GetInt("Diag", "MaxIter", 10);      // max no of macro-iterations 
  int ptol = Input.GetDouble("Diag", "pTol", 5);          // tolerance = 10^-pTol 
  int istartvec = Input.GetInt("Diag", "StartVector", 1); // this is different for Lanczos + Davidson and needs work
  int diagverbose = Input.GetInt("Diag", "Verbose", 1);
  Hel.SetVerbose(diagverbose);  // resets verbose from gridverbose to diagverbose

  double *energies = new double[nStates];
  Hel.DiagonalizeSetup(nStates, DiagFlag, maxSub, MaxIter, ptol);
    cout << "Diagverbose = " << diagverbose << " \n";
  //Hel.DiagonalizeSetup(nStates, DiagFlag, maxSub, 50, ptol);
  Hel.Diagonalize(istartvec, energies);
  for (int istate = 0; istate < nStates ; ++istate){
    cout << "Eel = " << energies[istate] << " Hartree = " << energies[istate]*AU2MEV << " meV\n";
  }


  Hel.WriteOneDCuts() ; 

  cout << "All done.\n"; cout.flush();

}





void Potential::SetupSphericalPotential( const double *PotPara)
{
  PotFlags.resize(1) ;
  PotFlags[0] =  1023 ;
  WellDepth = PotPara[0] ;
  WellWidth =  PotPara[1] ; 
  ForceConstant =  PotPara[2] ; 
  if (verbose > 0) {
    cout << "------------------------\n";
    cout << "FullerElec potential\n";
    cout << "PotPara 0: " << PotPara[0] << " (Depth) [Hartree]\n";
    cout << "PotPara 1: " << PotPara[1] << " (Width) [Bohrs]\n";
    cout << "PotPara 2: " << PotPara[2] << " (ForceConstant, k) [Bohrs]\n";
    cout << "------------------------\n";
  }
  //  SetupMinMax();
  cout << "Leaving SetupSphericalPot\n";
}

double Potential::EvaluateSphericalPotential(const double *x)
{
  double Vtotal = 0.0 ; 
  double dist = sqrt(x[0]*x[0]+ x[1]*x[1] + x[2]*x[2]) ;
  //Vtotal = WellDepth/(exp(WellWidth*dist)+exp(-WellWidth*dist)) ;
  Vtotal = 0.5*ForceConstant*dist*dist ; 
  //if (dist < WellWidth)
    //Vtotal = WellDepth ;
  //else
     //Vtotal = 0.0  ;
  //printf("%10.6f %10.6f %10.6f %14.6e %14.6e %14.6e %14.6e \n", x[0], x[1], x[2], Vel, VRep, VPol, Vtotal);
  return Vtotal ;
}


