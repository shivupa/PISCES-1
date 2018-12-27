#include <iomanip>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cassert>

#ifdef PISCES_USE_MPI
#include <mpi.h>
#endif

#include "timer.hpp"
#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"
#include "optimize.h"

// DPP incudes
#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"
#include "DPP.h"
#include "Molecule.h"

// excess electron includes
#include "Potential.h"
#include "DVR.h"
#include "ClusterAnion.h"

#include <boost/fortran/prototype.hpp>

using namespace std;


// lbfgs function
typedef double (*LBFGS_FUN)(double *);
typedef void (*LBFGS_FFUN)(double *, double *);
typedef void   (*LBFGS_DFUN)(const double *, double *);

BOOST_FORTRAN_SUBROUTINE( LBFGS, lbfgs,  LBFGS, 
                         (int*)//n
                         (int*)//m
                         (double*) //p
                         (LBFGS_FUN) 
                         (LBFGS_FFUN) 
                         (LBFGS_DFUN) 
                         (double *)//
                         );



namespace {
  WaterCluster WaterNN;
  ClusterAnion Wn ;
}





void Optimizer(dVec WaterPos, const Parameters InP)
{

  if (InP.nStates > 1) {
    cout << "Optimizer called with nStates=" << InP.nStates << "\nOptimizing only the lowest state so far\n";
    exit(1);
  }

  if (InP.PotFlag[3] != 3 && InP.PotFlag[3] != 6 ) {
    cout << "*** Error :: Only Polarization=3 gradients are implemented.\n Cannot use Polarization = " << InP.PotFlag[3]<< endl;
    exit(1);
  }
  WaterNN.SetStructure(InP.nWater, &WaterPos[0], 1, InP.CenterFlag, InP.KTFlag, InP.WMverbose);  
  Wn.SetUpClusterAnion(WaterPos, InP);
  int nWater = InP.nWater;
  dVec 	WaterConf(nWater*6);
  WaterNN.GetConfiguration(nWater, &WaterConf[0], 1.0, 0);
  double gtol = InP.gtol;
  int ndim = nWater*6;
  int mdim = 5; // no of corrections for BFGS (recommended 3 to 7)


  LBFGS(&ndim, &mdim, &WaterConf[0], GetClusterEnergy, GetGeometry, GetAnalGrad, &gtol);
  //  LBFGS(&ndim, &mdim, &WaterConf[0], GetClusterEnergy, GetNumGrad, &gtol);
  

//  GetClusterEnergy( &WaterConf[0]);
//   dVec numgrad(nWater*6);
//   Wn.GetNumGrad( &WaterConf[0], &numgrad[0]);
//  GetClusterEnergy( &WaterConf[0]);
//   dVec analgrad(nWater*6);
//   GetAnalGrad( &WaterConf[0], &analgrad[0]);

}




///////////////////////////////////////////////
//
//  this computes the total energy of the cluster defined in WaterN, Vel, and Hel
//  it assumes that all parameters have been set in an earlier single-point like evaluation
//  input is a configuration, that is center-of-mass coordinates and Euler angles for each monomer
//  maybe there is merit in building a super-class, anion or so
//

double GetClusterEnergy(double *configuration)
{
   return Wn.EnergyFromConfiguration(configuration);;
}

void GetAnalGrad(const double *configuration, double *analgrad)
{
    Wn.GetAnalGrad(configuration, analgrad);
}





void GetGeometry(double *configuration, double *energy)

{
  int nWater = WaterNN.ReportNoOfWaters();
  WaterNN.SetConfiguration(nWater, configuration, 1, 0);
  int nSites = 0, nCharges = 0, nDipoles = 0;
  WaterNN.ReportNoOfSites(nSites, nCharges, nDipoles);
  int nPntPol = WaterNN.ReportNoOfPolSites();

  static dVec Sites; Sites.resize(3*nSites);
  static dVec Charges; Charges.resize(nCharges);
  static dVec Dipoles; Dipoles.resize(3*nDipoles);
  static dVec DmuByDR; DmuByDR.resize(3*nSites*3*nDipoles);
  static iVec iqs; iqs.resize(nCharges);
  static iVec ids; ids.resize(nDipoles);
  static dVec Alphas; Alphas.resize(nPntPol);
  static iVec ips; ips.resize(nPntPol);
  static dVec Epc; Epc.resize(3*nPntPol);
  WaterNN.GetLists(nSites, &Sites[0], nCharges, &Charges[0], &iqs[0], nDipoles, &Dipoles[0], &ids[0],
                  &DmuByDR[0], 0);

  DumpGeometry(nWater, energy[0], &Sites[0] );
}

void DumpGeometry(int nWater, double energy, double *coor)
{
  static std::ofstream  gdump;
  if ( !gdump.is_open() ) 
    gdump.open( "Optimizing.xyz" );
  gdump
    << 3*nWater << "\n"
    << "  E = " <<  energy * AU2MEV << "e meV\n";
  for (int i = 0; i < nWater; ++i) {
    int add = 12*i;
    gdump 
      << "O " << coor[add+0]/Angs2Bohr << " " << coor[add+1]/Angs2Bohr << " " << coor[add+2]/Angs2Bohr << "\n"
      << "H " << coor[add+3]/Angs2Bohr << " " << coor[add+4]/Angs2Bohr << " " << coor[add+5]/Angs2Bohr << "\n"
      << "H " << coor[add+6]/Angs2Bohr << " " << coor[add+7]/Angs2Bohr << " " << coor[add+8]/Angs2Bohr << "\n"
      ;
  }
  gdump << std::flush;
}
