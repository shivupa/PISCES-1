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

#include <fstream>
#include <boost/fortran/prototype.hpp>

using namespace std;



//////////////////////////////////////////////////////////
///
///  common-block-like structure used by
///  optimize Optimizer(), single_point(), GetClusterEnergy(), GetAnalGrad(), and, GetNumGrad()
///  the one reason to use this construct is the lbfgs fortran function employed
///
namespace {
  int nMolecules;
  Molecule *Molecules;
  WaterCluster WaterN;
  Potential Vel;
  DVR Hel;
  int NoOfGradients;
  int optverbose;
}






///////////////////////////////////////////////////////////////////
///
///  Compute the energy at a single point specified in WaterPos
///
void polplot(const dVec WaterPos, const Parameters InP)
{

  // the plan is to get away from WaterCluster, and to replace it slowly with Molecule
  nMolecules = InP.nWater;
  Molecules = new Molecule[nMolecules];
  for (int im = 0; im < nMolecules; ++im) {
    cout << im << "\n"; cout.flush();
    Molecules[im].SetupDPPWater(&WaterPos[9*im]);
  }

  //  compute the energy of the neutral cluster
  WaterN.SetStructure(InP.nWater, &WaterPos[0], 1, InP.CenterFlag, InP.KTFlag, InP.WMverbose);   // 1=Angs
  double E0 = WaterN.CalcEnergy(InP.WMverbose);
  cout << "Energy of the neutral cluster E0 = " << E0 << "\n";
  
  //
  //  Now, dealing with the excess electron is more involved
  //
  
  // get the charges and dipoles the electron sees from the WaterCluster
  int nSites = 0, nCharges = 0, nDipoles = 0;
  WaterN.ReportNoOfSites(nSites, nCharges, nDipoles);
  if (InP.PotVerbose > 0) 
    cout << "The electron sees " << nCharges << " charges, " << nDipoles << " dipoles, and there are "
	 <<  nSites << " sites in the gradient calculations.\n";
  
  // 
  int nPntPol = WaterN.ReportNoOfPolSites();

  static dVec Sites; Sites.resize(3*nSites);
  static dVec Charges; Charges.resize(nCharges);
  static dVec Dipoles; Dipoles.resize(3*nDipoles);
  static dVec DmuByDR; DmuByDR.resize(3*nSites*3*nDipoles);
  static iVec iqs; iqs.resize(nCharges);
  static iVec ids; ids.resize(nDipoles);
  static dVec Alphas; Alphas.resize(nPntPol);
  static iVec ips; ips.resize(nPntPol);
  static dVec Epc; Epc.resize(3*nPntPol);
  WaterN.GetLists(nSites, &Sites[0], nCharges, &Charges[0], &iqs[0], nDipoles, &Dipoles[0], &ids[0], 
		  &DmuByDR[0], 1);

  // Get all the dipoles
  nDipoles *= 3;
  Dipoles.resize(3*nDipoles);
  ids.resize(nDipoles);
  WaterN.ReportInducedDipoles(nDipoles, &Dipoles[0], &ids[0], 1);

  // get the polarizable sites and the point-charge field (only used for Polarization=3)
  WaterN.ReportPolSitesAndField(nPntPol, &Alphas[0], &ips[0], &Epc[0], 1); 


  // set up the Potential for the excess electron
  Vel.SetVerbose(InP.PotVerbose);
  Vel.Setup(InP.PotFlag, nSites, &Sites[0], 
	    nCharges, &Charges[0], &iqs[0], 
	    nDipoles, &Dipoles[0], &ids[0], 
	    nPntPol, &Alphas[0], &ips[0], &Epc[0], 
	    &DmuByDR[0], InP.PotPara);


//vkv

// declaring required variables

   double Unit_vec[3]; 
   double R_dist = 0.0 ; 
   double Pos_vec[3];
   

// evaluating the unit vector 
   for (int i = 0; i < 3; ++i) {
	R_dist += (InP.direction[i] - InP.origin[i])*(InP.direction[i] - InP.origin[i]) ;
   }
   R_dist = sqrt(R_dist); 
   for (int i = 0; i < 3; ++i) {
	Unit_vec[i] = (InP.direction[i] - InP.origin[i])/R_dist ;
   }


// calculation of Vpol along the Unit_vec directions at the specified number of points 
// with the specifed step_size
        FILE  * pfile ;
        pfile = fopen("gnuplot.txt","w") ;
   	fprintf(pfile ,"# Vpc;\t\tVind;\t\tVrep;\t\tVpol;\t\tVtotal(eV);| Dist(Ang) Mag		x		y		z \n");


   for (int i = InP.polplot_pts*-1; i < InP.polplot_pts*1+1 ; ++i) {
	double MagPos_vec = 0.0 ;	
	double MagPos_vec2  = 0.0 ;	
	double energies[5] ;

	for (int j =0 ; j < 3 ; ++j) {
	   Pos_vec[j] = InP.origin[j] + Unit_vec[j]*InP.polplot_stepsize*i ;
	   MagPos_vec2 += (Pos_vec[j] - InP.origin[j]) * (Pos_vec[j] - InP.origin[j]) ;
	   Pos_vec[j] = Pos_vec[j]*Angs2Bohr 	;
	}

	MagPos_vec = sqrt(MagPos_vec2) ; 	
	if (i < 0)
	   MagPos_vec = -1.0*MagPos_vec;
	double Vtotal = Vel.Evaluate(Pos_vec) ; 
   	
	
	Vel.ReportEnergies(5, energies);


	fprintf(pfile," %f \t  %f \t  %f \t  %f \t  %f \t  %f \t    %f \t   %f \t   %f\n", energies[0]*AU2EV, energies[1]*AU2EV , energies[2]*AU2EV,  energies[3]*AU2EV, energies[4]*AU2EV  , MagPos_vec , Pos_vec[0]*Bohr2Angs, Pos_vec[1]*Bohr2Angs, Pos_vec[2]*Bohr2Angs) ;
   }

	fclose (pfile);
//vkv
}



