#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>


#include "tsin.h"
#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"
#include "GetInput.h"
#include "ReadCubeFile.h"

// DPP incudes
#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"
#include "DPP.h"

// excess electron includes
#include "Potential.h"
#include "DVR.h"


// for potfit
#include "lapackblas.h"
#include "WriteCubeFile.h"


using namespace std;
//
//  all the variables and data structures chisquared needs to compute Eel
//
namespace{
  iVec PotFlag(32);
  int nSites = 0;
  int nCharges = 0;
  int nDipoles = 0;
  int nPntPols = 0;
  dVec Sites;
  dVec Charges;
  dVec Dipoles;
  dVec Alphas;
  dVec Epc;
  dVec DmuByDR;
  iVec iqs;
  iVec ids;
  iVec ips;
  double PotPara[32];
  
  dVec energies;
  
  class Potential Velfit;
  class DVR Helfit;
  
  int nCubePts = 0;
  double *eomcube = 0;
  double *weights = 0;
  double *dvrcube = 0;
  
  int nParaOpt = 0;
  int *mapping = 0;
  double *ParaOpt = 0;
  
  int ChiOutput = 1;
}


//
//   Function defined in this file:
//
double chisquared(double *fitpara);
double ComputeChiSquared(int npts, double *dvrcube, double *eomcube, double *weights);

void ComputeWeights(const int *ngrid, const double *x0, const double *steplen, 
		    int nWaters, const double *WaterCoor, 
		    double *weights, const double *eomcube,
		    double LowCutOff, double HighCutOff, double exponent);

double rand01(void);
double randm11(void);


//////////////////////////////////////////////////////////////////////////////////
///
///  PotFit:
///
///
///
void potfit(const dVec& WaterCoor0, const Parameters& InP)
{
  cout << "=======================================================================================\n";
  cout << "This is potfit\n";
  cout << "Make sure the EOM cube file and the cube file produced by pisces look exactly identical\n";
  cout << "==================================================================================== ==\n";

  int nWater = InP.nWater;
  const double *WaterPos = &WaterCoor0[0];

  //
  //  compute the energy of the neutral cluster
  //
  WaterCluster WaterN;
  WaterN.SetStructure(nWater, WaterPos, 1, InP.CenterFlag, InP.KTFlag, InP.WMverbose);   // 1=Angs
  double E0 = WaterN.CalcEnergy(InP.WMverbose);
  cout << "Energy of the neutral cluster E0 = " << E0 << "\n";

  //
  // get charges and dipoles the electron sees from the WaterCluster
  //
  WaterN.ReportNoOfSites(nSites, nCharges, nDipoles);
  if (InP.WMverbose > 0)
    cout << "The electron sees " << nCharges << " charges, " << nDipoles << " dipoles, and there are " 
	 <<  nSites << " sites in the gradient calculations.\n"; 
  Sites.resize(3*nSites);
  Charges.resize(nCharges);
  Dipoles.resize(3*nDipoles);
  DmuByDR.resize(3*nSites*3*nDipoles);
  iqs.resize(nCharges);
  ids.resize(nDipoles);
  WaterN.GetLists(nSites, &Sites[0], nCharges, &Charges[0], &iqs[0], nDipoles, &Dipoles[0], &ids[0], &DmuByDR[0], 1);

  // needed for writing cube files
  iVec NucCharge(3*nWater);
  for (int k = 0; k < nWater; ++k) {
    NucCharge[3*k+0] = 8;
    NucCharge[3*k+1] = 1;
    NucCharge[3*k+2] = 1;
  }




  //
  //  read the cube file with the data to be fitted against
  //  all points are weighted according to their distance to the nearest O atom
  //  too close must have weight 0, as the EOM NO will oscillate in this region
  //  too close to the edge of the grid will suffer from non-identical boundary conditions
  //  in the DVR and EOM calculations
  //  then the number of points should increase quadratically 
  //  (so inverse square weights will weigh all distances equally) 
  //
  cout << "\nReading the reference wave function from the cube file StartVector.cub\n";
  int ncube[3] = {0,0,0};
  double x0cube[3] = {0,0,0};
  double dxcube[3] = {0,0,0};
  ReadCubeFile(ncube, x0cube, dxcube, &eomcube);
  nCubePts = ncube[0] * ncube[1] * ncube[2];
  dvrcube = new double[nCubePts];
  weights = new double[nCubePts];
  int ione = 1;
  double eomnrm = dnrm2(&nCubePts, eomcube, &ione);
  cout << "DVR-norm of the EOM Wf = " << eomnrm << "\n";
  eomnrm = 1.0 / eomnrm;
  dscal(&nCubePts, &eomnrm, eomcube, &ione);
  eomnrm = dnrm2(&nCubePts, eomcube, &ione);
  cout << "DVR-norm of the EOM Wf = " << eomnrm << "\n";
  WriteCubeFile("EOMNO_renorm.cub", 3*nWater, &NucCharge[0], WaterPos, ncube, x0cube, dxcube, eomcube, 0); 
  // compute weights
  ComputeWeights(ncube, x0cube, dxcube, nWater, WaterPos, weights, eomcube, InP.LowCutOff, InP.HighCutOff, InP.wexp);
  WriteCubeFile("Weights.cub", 3*nWater, &NucCharge[0], WaterPos, ncube, x0cube, dxcube, weights, 0);



  // set up a grid (a DVR) for the wavefunction of the electron
  cout << "\nSetting up the DVR grid\n";
  Helfit.SetupDVR(InP.ngrid, InP.DVRType, InP.Sampling, InP.gpara, InP.gridverbose);
  //
  //  output for checking whether the grids in eomcube and dvrcube are compatible
  //  this is just a bare bones check
  //  make sure the atoms are at identical positions
  //
  cout << "c.f. Grid from cub-file:\n";
  for (int idim = 0; idim < 3; ++idim) {
    cout << idim << ": " << ncube[idim] << " gridpoints, from " << x0cube[idim] << " to " 
       << x0cube[idim]+(ncube[idim]-1)*dxcube[idim] << " dx = " << dxcube[idim] << "\n";
  }
  Velfit.SetVerbose(InP.PotVerbose);
  for (int p=0; p < 32; ++p)
    PotPara[p] = InP.PotPara[p];
  PotFlag = InP.PotFlag;
  Velfit.Setup(PotFlag, nSites, &Sites[0], 
	       nCharges, &Charges[0], &iqs[0], 
	       nDipoles, &Dipoles[0], &ids[0], 
	       nPntPols, &Alphas[0], &ips[0], &Epc[0],
	       &DmuByDR[0], PotPara);


  // determine which parameters will be optimized
  nParaOpt =  InP.nParaOpt;
  mapping = new int[nParaOpt];
  ParaOpt = new double[nParaOpt];
  cout << "\nPotFit optimizes " << nParaOpt << " parameters:\n";
  for (int ip = 0; ip < nParaOpt; ++ip) {
    mapping[ip] = InP.mapping[ip];
    ParaOpt[ip] = InP.PotPara[mapping[ip]];
    cout << " no. " << mapping[ip] << " = " << ParaOpt[ip] << "\n";
  }



  // compute the potential on the grid
  cout << "\nEvaluate Vel on grid\n"; cout.flush();
  Helfit.ComputePotential(Velfit);
  cout << "\nEvaluate Vel on grid\n"; cout.flush();

  // compute the DVR wavefunction
  Helfit.SetVerbose(InP.diagverbose);  // resets verbose from gridverbose to diagverbose
  energies.resize(InP.nStates);
  Helfit.DiagonalizeSetup(InP.nStates, InP.DiagMethod, InP.maxSub, InP.maxIter, InP.ptol);
  Helfit.Diagonalize(InP.istartvec, &energies[0]);
  cout << "Eel = " << energies[0] << " Hartree = " << energies[0]*AU2MEV << " meV\n";
  Helfit.GetWaveFnCube(1, dvrcube);

  // evaluate weighted difference with EOM-NO
  double chi2 = ComputeChiSquared(nCubePts, dvrcube, eomcube, weights);
  cout << "Initial Chi**2 = " << chi2 << "\n";

  // plot the initial wavefunction
  Helfit.WriteCubeFile(1, "InitialWaveFn.cub", 3*nWater, &NucCharge[0], WaterPos, Angs2Bohr);

  // settings for subsequent calls of chisquared
  Helfit.SetVerbose(0);
  Velfit.SetVerbose(0);

  
  if (InP.minimizer == 0) {
    ChiOutput = 1;

    // MC optimization
    double maxstep = 0.05;     // 5% 
    double tolerance = 0.002;  // 0.2%
    int nSteps = 50 * nParaOpt;

    double *TrialPara = new double[nParaOpt];

    while (maxstep > tolerance) {
      int iStep = 0;
      while (iStep < nSteps) {
	cout << "iStep=" << iStep << "\n";
	for (int ip = 0; ip < nParaOpt; ++ip) {
	  cout << "ip=" << ip << "\n";
	  TrialPara[ip] = ParaOpt[ip] * (1.0 + maxstep*randm11());
	}
	double chi2Trial = chisquared(TrialPara);
	if (chi2Trial < chi2) {
	  for (int ip = 0; ip < nParaOpt; ++ip)
	    ParaOpt[ip] = TrialPara[ip];
	  chi2 = chi2Trial;
	  iStep = 0;
	}
	else
	  iStep += 1;
      }
      maxstep *= 0.5;
    }
    cout << "Optimized Parameters:\n";
    for (int ip = 0; ip < nParaOpt; ++ip)
      cout << ip << "   " << ParaOpt[ip] << "\n";

    delete[] TrialPara;

  }

// else if (InP.minimizer == 1)
// {
//   // call amoeba minimizer
//   ChiOutput = 1;
//   int np = nParaOpt;
//   int mp = np + 1;
//   int j, k;
//   double ftol = 2e-5;
//   int iter;
//   double *p = new double[np*mp];          // array of vertices
//   double *y = new double[mp];             // chisquared at vertices
//   double *para = new double[np];          // passing parameters to chisquared
//   
//   
//   // set all vertices to the start parameters
//   for (j = 0; j < mp; ++j)
//     for (k = 0; k < np; ++k)
//       p[j*np + k] = ParaOpt[k];
//   
//   // add 5% + eps to the (k-1)-th element of vertex k=1,np
//   const double eps = 1e-3;
//   double RelStep = 0.05;
//   for (k = 1; k < mp; ++k)
//     p[k*np+(k-1)] += RelStep * p[k*np+(k-1)] + eps;
//   
//   // evaluate chi_square at the vertices
//   for (j = 0; j < mp; ++j) {
//     for (k = 0; k < np; ++k)
//       para[k] = p[j*np + k];
//     y[j] = chisquared(para);
//   }
//   printf("\nInitial vertices:\n");
//   printf("    chi^2\n");
//   for (j = 0; j < mp; ++j) {
//     printf(" %10.3e :", y[j]);
//     if (np < 8 || InP.optverbose > 0)
//       for (k = 0; k < np; ++k)
//         printf("%10.3e ", p[j*np + k]);
//     printf("\n");
//   }
//
//
//   printf("BeginAmoeba\n");    
//   printf("# Chi**2,  Parameters:\n#------------------------------\n");
//   ChiOutput = 1;
//   int info = amoeba(p, y, np, ftol, chisquared, &iter);
//   if (info != 0) {
//     cout << "error in amoeba : " << info << endl;
//   }
//   printf("EndAmoeba\n");
//   
//   printf("\nFinal vertices after %d iterations:\n", iter);
//   printf("    chi^2\n");
//   int imin = 0;
//   double chimin = y[0];
//   for (j = 0; j < mp; ++j) {
//     printf(" %10.6e :", y[j]);
//     if (y[j] < chimin) {
//       chimin = y[j];
//       imin = j;
//     }
//     if (np < 8 || InP.optverbose > 0)
//       for (k = 0; k < np; ++k)
//         printf("%10.4e ", p[np*j + k]);
//     printf("\n");
//   }
//   
//   cout << "\nThe best vertex is no. " << imin << endl;
//   cout << "Chi**2 = " << y[imin] << endl;
//   for (k = 0; k < np; ++k) {
//     ParaOpt[k] =  p[np*imin + k];
//     cout << "Parameter no " << k << " = " << ParaOpt[k] << endl;
//   }
//
//   delete[] para; delete[] p; delete[] y;
//
// }
//
// else if (InP.minimizer == 2) {
//   cout << "call Powell minimizer\n"; cout.flush();
//   int niter = 0;
//   int itmax = 50;
//   double finalchi2 = 0;
//   double ftol = 1e-6;
//   finalchi2 = Powell(nParaOpt, ParaOpt, ftol, &niter, itmax, chisquared, 1);
// }
// 
  else {
    cout << "Unknown miniization method\n";
    exit(42);
  }
 
 
 
  // at the optimized point do the calc again, and print Eel
  ChiOutput = 0;
  cout << "Final Chi**2 = " << chisquared(ParaOpt) << "\n";  
  cout << "Eel = " << energies[0] << " Hartree = " << energies[0]*AU2MEV << " meV\n";
  // plot the final wavefunction
  Helfit.WriteCubeFile(1, "OptWaveFn.cub", 3*nWater, &NucCharge[0], WaterPos, Angs2Bohr);



}




////////////////////////////////////////////////////////////////////////////
///
///  double chisquared(double *fitpara) given to simplex or Powell minimizers
///
///  ChiOutput = 0   no output
///              1   output for table
///              3   debug output
///
double chisquared(double *fitpara)
{

  // copy fitpara into PotPara according to mapping
  if (ChiOutput >= 3)
    cout << "Evaluate Chi2\n";
  for (int ip = 0; ip < nParaOpt; ++ip) {
    PotPara[mapping[ip]] = fitpara[ip];
    if (ChiOutput >= 3)
      cout << "  Parameter " << ip << " = " << fitpara[ip] << " is PotPara[" << mapping[ip] << "]\n";
  }

  // compute wavefunction of the electron with fitpara setup potential
  if (ChiOutput >= 3) cout << "Call Vel.Setup\n";
  Velfit.Setup(PotFlag, nSites, &Sites[0], 
	       nCharges, &Charges[0], &iqs[0], 
	       nDipoles, &Dipoles[0], &ids[0], 
	       nPntPols, &Alphas[0], &ips[0], &Epc[0],	       
	       &DmuByDR[0], PotPara);
  if (ChiOutput >= 3) cout << "Call Hel.ComputePotential\n";
  Helfit.ComputePotential(Velfit);
  if (ChiOutput >= 3) cout << "Call Hel.Diagonalize\n";
  Helfit.Diagonalize(2, &energies[0]);
  if (ChiOutput >= 3) cout << "Call Hel.GetWaveFnCube\n";
  Helfit.GetWaveFnCube(1, dvrcube);

  // evaluate weighted difference with EOM-NO
  double chi2 = ComputeChiSquared(nCubePts, dvrcube, eomcube, weights);

  // debug output
  if (ChiOutput == 1) {
    printf("%12.6e", chi2);
    for (int ip = 0; ip < nParaOpt; ++ip)
      printf("  %12.6f", fitpara[ip]);
    printf("\n");
  }
  else if (ChiOutput >= 3) {
    for (int ip = 0; ip < nParaOpt; ++ip)
      cout << "P" << ip << "= " << fitpara[ip] << "  ";
    cout << "  Chi2= " << chi2 << "\n";
  }


  return chi2;

}















////////////////////////////////////
///
///  chi**2 of a DVR and a EOM or ADC cube with npts points
///  both files are orbitals; 
///  the difference of the densities is computed, i.e., the orbitals are squared
///  therefore the difference itself is not again squared, but just fabsed 
///
double ComputeChiSquared(int npts, double *dvrcube, double *eomcube, double *weights)
{
  double chisq = 0;
  double delta = 0;
  for (int i = 0; i < npts; ++i) {
    if (weights[i] > 0) {
      delta = dvrcube[i]*dvrcube[i] - eomcube[i]*eomcube[i];
      chisq += fabs(delta) * weights[i];
    }
  }
  return chisq;
}



//////////////////////////////////////////////////////////
///
///  weights of the grid points for fitting
///  - anything too close to an O, i.e., < 2 Bohr must have weight=0 (EOM NO oscillates)
///  - anything too far cannot have a large weight owing to the different boundary
///    conditions in DVR and GTO wavefunctions
///
///
///  for dO between Low and High cutoffs  w = dO^exponent
///
void ComputeWeights(const int *ngrid, const double *x0, const double *steplen, 
		    int nWaters, const double *WaterCoor, 
		    double *weights, const double *eomcube,
		    double LowCutOff, double HighCutOff, double exponent)
{
  int nx = ngrid[0];
  int ny = ngrid[1];
  int nz = ngrid[2];
  double *weightptr = weights;
  const double *eomptr = eomcube;


  int nWeights = 0;

  cout << "\nComputing weights; use points between " << LowCutOff << " and " << HighCutOff << "\n"; 


  for (int ix = 0; ix < nx; ++ix) {
    double xp = x0[0] + ix * steplen[0];
    for (int iy = 0; iy < ny; ++iy) {
      double yp = x0[1] + iy * steplen[1];
      for (int iz = 0; iz < nz; ++iz) {
	double zp = x0[2] + iz * steplen[2];
	// find smallest distance to any O
	double dOmin = 10000000.0;
	for (int iw = 0; iw < nWaters; iw++) {
	  double dx = xp - WaterCoor[9*iw];
	  double dy = yp - WaterCoor[9*iw+1];
	  double dz = zp - WaterCoor[9*iw+2];
	  double dO = dx*dx + dy*dy + dz*dz;
	  if (dO < dOmin)
	    dOmin = dO;
	}
	dOmin = sqrt(dOmin);
	// formula for the weights
	if (dOmin > LowCutOff && dOmin < HighCutOff) {
	  *weightptr = pow(dOmin,exponent) * sqrt(fabs(*eomptr));
	  nWeights ++;
	}
	else
	  *weightptr = 0;
	weightptr ++;
	eomptr ++;
      }
    }
  }
  cout << "Number of grid points with non-zero weights: " << nWeights 
       << " (" << 100* nWeights / double(nx*ny*nz) << "%)\n"; 


  // normalize weights
  double sum = 0;
  for (int i = 0; i < nx*ny*nz; ++i)
    sum += weights[i];
  sum = 1.0 / sum;
  for (int i = 0; i < nx*ny*nz; ++i)  
    weights[i] *= sum;
}




//////////////////////////////////////////////////////////////////
///
/// get random number in [0,1)
///
double rand01(void)
{
  return double(rand()) / (RAND_MAX + 1.0);
}

/////////////////////////////////////////
///
/// get random number in (-1,1)
///
double randm11(void)
{
  if (rand() > RAND_MAX/2)
    return double(rand()) / (RAND_MAX + 1.0);
  else
    return double(-rand()) / (RAND_MAX + 1.0);
}

