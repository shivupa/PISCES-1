#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>

#include <mpi.h>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>

#include "constants.h"
#include "ttm_data.h"
#include "Potential.h"
#include "lapackblas.h"
#include "timer.hpp"

// DPP incudes
#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"
#include "DPP.h"
#include "Molecule.h"

#include "erf.h"

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

//
//  verbose > 0  : print min and max of all contributions
//  verbose > 15 : print potential at every grid point
//
void Potential::SetVerbose(int v)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   verbose = v;
   if(rank==0)cout << "Potential::SetVerbose called with " << v << "\n";

}

////////////////////////////////////////////////////////////////////////////
//
//  initialize the potential of the excess electron
//  set Potential.h for comments on the variables
//  
void Potential::Setup(const iVec& potflaginp, int nSites, const double *Sites, 
		      int nCharges, const double *Charges, const int *iqs, 
		      int nDipoles, const double *Dipoles, const int *ids,
		      int nPntPols, const double *Alphas, const int *ips, const double *Epc,
		      const double *DmuByDR, const double *PotPara)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  PotFlags = potflaginp;

  switch (PotFlags[0])
    {
    case 1:
    case 2:
    case 3:
    case 4:
      SetupDPPGTOP(nSites, Sites, 
		   nCharges, Charges, iqs, 
		   nDipoles, Dipoles, ids, 
		   nPntPols, Alphas, ips, Epc, 
		   DmuByDR, PotPara);
      break;
    case 11:
      SetupDPPTB(nSites, Sites, 
		 nCharges, Charges, iqs, 
		 nDipoles, Dipoles, ids,
		 DmuByDR, PotPara);
      break;
    default:
      // this should not happen (tested in GetInput)
      if(rank==0)cout << "Potential::Setup: unknown potential for the excess electron: " << PotFlags[0] << "\n";
      exit(1);
    }

  SetupMinMax(); // min and max over the grid can be printed
}



double Potential::MinDistCheck(const double *relectron)
{
 int intflag;
 double Rx, Ry, Rz, R2, R;
 double MinR=50.0;
  for (int i = 0; i < nSites; ++i) {
    Rx = relectron[0] - Site[3*i+0];
    Ry = relectron[1] - Site[3*i+1];
    Rz = relectron[2] - Site[3*i+2];
    R2 = Rx*Rx + Ry*Ry + Rz*Rz;
    R = sqrt(R2);
   // if(rank==0)cout<<"R = "<<R<<endl;
    if (MinR > R) MinR=R;
  }

 return MinR;
}





///////////////////////////////////////////////////////////////////////////////
//
//  this function evaluates the potential that has been setup most recently
//  it is called by the DVR class to compute the potential 
//
//  for Evaluation there is one generic function for DPP electrostatics
//  with a GTO repulsive core and polarizable sites
//
//  for other functional forms, other functions (and setup-functions) are needed 
//
double Potential::Evaluate(const double *relectron)
{
  
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // compute a list of distances of all sites to the point r=(x,y,z)
  // this is identical for all potentials
  for (int i = 0; i < nSites; ++i) {
    Rx[i] = relectron[0] - Site[3*i+0];
    Ry[i] = relectron[1] - Site[3*i+1];
    Rz[i] = relectron[2] - Site[3*i+2];
    R2[i] = Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i];
    R[i] = sqrt(R2[i]);
  //  cout<<"relectron= "<<relectron[0]<<" "<<relectron[1]<<" "<<relectron[2]<<endl;
  //  cout<<"Site ="<<Site[3*i+0]<<" "<<Site[3*i+1]<<" "<<Site[3*i+2]<<endl;
  //  cout<<"R["<<i<<"]="<<Rx[i]<<" "<<Ry[i]<<" "<<Rz[i]<<" "<<R[i]<<endl;
    Rminus3[i] = 1.0 / (R2[i] * R[i]);
  }
  
 // cout<<"--------------------------------------"<<endl;
//  cout<<"relectron : "<<relectron[0]<<" "<< relectron[1]<<" "<<relectron[2]<< endl;


  switch (PotFlags[0])
    {
    case 0:
      if(rank==0)cout << "Error in EvaluateVel: PotFlag=0, so, no potential has been set up so far.\n";
      exit(1);
      break;
    case 1:
    case 2:
    case 3:
    case 4:
      return EvaluateDPPGTOP(relectron);
      break;
    case 11:
      return EvaluateDPPTB(relectron);
      break;
    case 101 :
      return EvaluateBloomfield(relectron);
      break;
    case 1021 :
      return EvaluateCsixty(relectron);
      break;
    case 1023 :
      return EvaluateSphericalPotential(relectron);
      break ;  
    default:
      if(rank==0)cout << "Error in EvaluateVel: unknown PotFlag = " << PotFlags[0] << "; this should not happen\n";
      exit(1);
    } 
  // never ever get here
  return 0;
}

/*

void Potential:: EvaluateGradient(
                                  const double *relectron, 
                                  dVec& Grad, dVec& tmu, dVec& mu_cross_mu, double wavefn,
                                  class WaterCluster &WaterN) 
{

   // compute a list of distances of all sites to the point r=(x,y,z)
   for (int i = 0; i < nSites; ++i) {
      Rx[i] = relectron[0] - Site[3*i+0];
      Ry[i] = relectron[1] - Site[3*i+1];
      Rz[i] = relectron[2] - Site[3*i+2];
      R2[i] = Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i];
      R[i] = sqrt(R2[i]);
      Rminus3[i] = 1.0 / (R2[i] * R[i]);
   }




   switch (PotFlags[0])
   {
   case 0:
      if(rank==0)cout << "Error in EvaluateVel: so far no potential has been set up.\n";
      exit(1);
      break;
   case 1: // this should work for 2 as well
   case 2:
   case 3:
      EvaluateDPPGTOPGradient(relectron,Grad, tmu, mu_cross_mu, wavefn,  WaterN);
      break;
   default:
      if(rank==0)cout << "Error in EvaluateGradient: unknown PotFlags[0]; this should not happen\n";
      exit(1);
   }

   // never ever get here
}

*/
////////////////////////////////////////////////////////////////////////////////////
//
//  A DPP based potential, that is, the Sites are O1, H1a, H1b, M1, O2, H2a, H2b, M2, ....
//
//  the charges and induced dipoles are provided by the DPP module
//  for the neutral cluster
//  polarizable sites with alpha = PotPara[0] and damping PotPara[1] are put on the M-sites
//  the 6S repulsive core is put on the Hs and the scaling factor is 1/PotPara[2]
//
void Potential::SetupDPPGTOP(
                            int nr, 
                            const double *r, 
                            int nq, 
                            const double *q, 
                            const int *iq, 
                            int nd, 
                            const double *d, 
                            const int *id,
			    int npp, 
			    const double *am, 
			    const int *ip, 
			    const double *epc,
                            const double *DmuByDR, 
                            const double *PotPara
                            )
{ 

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  nReturnEnergies = 5;
  ReturnEnergies.resize(nReturnEnergies); // Vpc, Vind, Vrep, Vpol, and Vtotal

  //progress_timer tmr("SetupDPPGTOP:", verbose);
  

  // Flags
  DampType        = PotFlags[1];  // GTO or effective-R
  SigmaOFlag      = PotFlags[2];  // Separate scaling for GTOs on O  
  PolType         = PotFlags[3];  // none, single, distributed, or fully self-consistent polarizable sites
  AdiabaticPolPot = PotFlags[4];  // adiabatic potential or its 1st-order expression:  -0.5*a/R**4 


  //Set the distance cut-off for interpolation, Rtol
   Rtol = PotPara[10];
  

  // set parameters as called for by flags
  Alpha = PotPara[0];
  VRepScale = PotPara[2];
  if (SigmaOFlag == 1)
    VRepScaleO = PotPara[5];
  switch (DampType)
    {
    case 1:
      ChargeDamping = 1.0 / (PotPara[3]*PotPara[3]);   // Gaussian damping width
      DipoleDamping = 1.0 / (PotPara[4]*PotPara[4]);   // the damping factor contains terms of the form
      PolDamping    = 1.0 / (PotPara[1]*PotPara[1]);   // (1-exp(-r*r/damp/damp)
      break;
    case 2:
      ChargeDamping = fabs(PotPara[3]);   // distance r0 used for distance damping
      DipoleDamping = fabs(PotPara[4]);   // if r < r0, r is replaced by 
      PolDamping = fabs(PotPara[1]);      // -0.5*r^4/r0^3 + r^3/r0^2 + 0.5*r0   (do not ponder, plot)
      break;
    default: 
      if(rank==0)cout << " Potential::SetupDPPGTOP DampFlag this should never happen.\n"; exit(1);
    }
  if (AdiabaticPolPot == 1)
    EpsDrude = PotPara[6];


  if (verbose > 0) {
    switch (PotFlags[0])
      {
      case 1: if(rank==0)cout << "DPP6SP potential:  6s GTO Vrep plus Polarization\n"; break;
      case 2: if(rank==0)cout << "DPP12SP potential:  12s GTO Vrep plus Polarization\n"; break;
      case 3: if(rank==0)cout << "DPP12SP potential:  4s STO Vrep plus Polarization\n"; break;
      case 4: if(rank==0)cout << "DPP12SP potential:  no Vrep, but Polarization\n"; break;
      default: if(rank==0)cout << "Potential::SetupDPPGTOP PotFlag[0] this should never happen.\n"; exit(1);
      } 
    if(rank==0)cout << "PotPara 0: " << PotPara[0] << " (alpha)\n";
    if(rank==0)cout << "PotPara 1: " << PotPara[1] << " (damping length for polarizable sites [Bohr])\n";
    if(rank==0)cout << "PotPara 2: " << PotPara[2] << " (sigma for Vrep scaling)\n";
    if(rank==0)cout << "PotPara 3: " << PotPara[3] << " (Coulomb damping length [Bohr])\n";
    if(rank==0)cout << "PotPara 4: " << PotPara[4] << " (Dipole damping length [Bohr])\n";
    if (SigmaOFlag == 1)
	if(rank==0)cout << "PotPara 4: " << PotPara[5] << " (separate sigma for O-sites scaling)\n";
    if(rank==0)cout << "DampType = " << DampType; 
    switch (DampType)
      {
      case 1: if(rank==0)cout << " (GTO damping as in the Drude code)\n"; break;
      case 2: if(rank==0)cout << " (Distance Damping)\n"; break;
      default: if(rank==0)cout << " Potential::SetupDPPGTOP DampFlag this should never happen.\n"; exit(1);
      } 
  }

  
  int nWater = nr / 4;  // DPP uses 4 Sites per Water in the order O, H, H, M
  if (nr % 4 != 0) 
    {if(rank==0)cout << "SetupDPPGTOP: no of sites " << nr << " is not a multiple of 4\n"; exit(1);}
  SetSites(nr, r);
  SetCharges(nq, q, iq);
  if (PolType == 3 ) {
    if (verbose > 0) {
    if(rank==0)cout << "Warning: " << nd << " point dipoles from the water model are ignored.\n"
	 << "  This assumes that these dipoles come from water-water induction.\n"
	 << "  If this is not true, changes in the code are needed.\n";
    }
  }
  else
    SetDipoles(nd, d, id); // for self-consistent polarization ignore dipoles from water-water induction (no double counting)


  // polarization potential
  switch (PolType)
    {
    case 1:
      if (AdiabaticPolPot == 1)
	if(rank==0)cout << " Full adiabatic expression for the polarization potential with EpsDrude = " << EpsDrude << " Hartree\n";
      else 
	if(rank==0)cout << " Standard polarization potential of non-interacting sites\n";
      // single site per water isotropic polarisabilities
      nPPS = nWater;
      PPS.resize(nPPS);
      PPSSite.resize(nPPS);
      for (int i = 0; i < nPPS; ++i) {
	PPSSite[i] = 4*i+3;  // this is the i-th M-site
	PPS[i] = Alpha;
      }
      break;
    case 2:
      if(rank==0)cout << " Interacting atomic polarizability on each water, but non-interacting waters\n";

      // very dirty, just for testing, relies on O, H, H, M order of Sites in DPP
      if (nMolPol != nWater) {
	if (nMolPol > 0)
	  delete[] MolPol;
	nMolPol = nWater;
	MolPol = new DistributedPolarizabilities[nMolPol];
      }
      
      for (int im = 0; im < nMolPol; ++im) {
	MolPol[im].nAtoms = 3;

	MolPol[im].SiteList.resize(3);
	MolPol[im].SiteList[0] = 4*im;
	MolPol[im].SiteList[1] = 4*im + 1;
	MolPol[im].SiteList[2] = 4*im + 2;

	MolPol[im].Alpha.resize(3);
	MolPol[im].Alpha[0] = TTM_POL_O;
	MolPol[im].Alpha[1] = TTM_POL_H;
	MolPol[im].Alpha[2] = TTM_POL_H;

	MolPol[im].InvA.resize(81);  // 3*nAtoms x 3*nAtoms	
 	void ComputeInvA(int nSites, const double *R, const double *alpha, double aThole, double *InvA, int verbose); // defined in MolPolAux.cpp

	// so far aThole=0.3  this should be different for DPP and DPP2
	
	ComputeInvA(3, &r[12*im], &MolPol[im].Alpha[0], 0.3, &MolPol[im].InvA[0], 1);
      }
      break;
    case 3:   // self-consistent-polarizability 
    case 33:  // self-consistent-polaribility-33
    case 5:  // self-consistent-polaribility-5
    case 6:  // self-consistent-polaribility-5
     
      {
	if (verbose > 0) if(rank==0)cout << "\nSelf-consistent water-water and water-electron polarizabilities\n";
	
	// treat as one big MolPol
	if (nMolPol != 1) {
	  if (nMolPol > 0)
	    delete[] MolPol;
	  nMolPol = 1;
	  MolPol = new DistributedPolarizabilities[nMolPol];
	}

	// copy data from arguments into MolPol[0]
	MolPol[0].nAtoms = npp;
	MolPol[0].SiteList.resize(npp+1);
	MolPol[0].Alpha.resize(npp);
	MolPol[0].InvA.resize(3*npp*3*npp);
	MolPol[0].Epc.resize(3*npp);
//	MolPol[0].dTij.resize(npp*npp*27);
	MolPol[0].dT_x.resize(npp*npp*9);
	MolPol[0].dT_y.resize(npp*npp*9);
	MolPol[0].dT_z.resize(npp*npp*9);
	MolPol[0].Tensor.resize(npp*npp*9);
 
	dVec Rpps; Rpps.resize(3*npp);  // needed for calling ComputeInvA
	for (int i = 0; i < npp; ++i) {
	  MolPol[0].SiteList[i] = ip[i];
	  MolPol[0].Alpha[i] = am[i];
	  Rpps[3*i+0] = r[3*ip[i]+0];
	  Rpps[3*i+1] = r[3*ip[i]+1];
	  Rpps[3*i+2] = r[3*ip[i]+2];
	  MolPol[0].Epc[3*i+0] = epc[3*i+0];
	  MolPol[0].Epc[3*i+1] = epc[3*i+1];
	  MolPol[0].Epc[3*i+2] = epc[3*i+2];
	}
	MolPol[0].SiteList.back() = -1; // invalid value to skip M
	
	void ComputeInvA(int nSites, const double *R, const double *alpha, double aThole, double *InvA, int verbose); // defined in PolAux.cpp

	// so far aThole=0.3  this should be different for DPP and DPP2
	ComputeInvA(npp, &Rpps[0], &MolPol[0].Alpha[0], 0.3, &MolPol[0].InvA[0], 1);
        
        // Calculate Tensor Tae Hoon Choi
        //void BuildTensor(int nSites, const double *R, const double *alpha, double aThole, double *Tensor);
	//BuildTensor(npp, &Rpps[0], &MolPol[0].Alpha[0], 0.3, &MolPol[0].Tensor[0]);

        //vkv 3rd march 2012
	// building the supermatrix for the derivatives of Tij, which is defined in DerivDDTensor
        //
        void BuildSuperdT(int nPolSites, const double *R, const double *alpha, double aThole, double *dT_x, double *dT_y, double *dT_z) ;	
        BuildSuperdT(npp, &Rpps[0], &MolPol[0].Alpha[0], 0.3, &MolPol[0].dT_x[0] , &MolPol[0].dT_y[0] , &MolPol[0].dT_z[0]) ;  


	// compute the Water-Water induced dipoles in Rpps
	// and then the Water-Water polarization energy as an energy offset 
	// compute dipoles induced by Efield
	double dzero = 0.0;
	double done = 1.0;
	int one = 1;
	int n = 3*npp;
      //dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &epc[0], &one, &dzero, &Rpps[0], &one);
      //VpolWaterWater = -0.5 * ddot(&n, &epc[0], &one, &Rpps[0], &one);
	if(rank==0)cout << "Water-Water polarization energy = " << VpolWaterWater << "\n";

/*testing the supermatrix
        for (int i = 0; i < npp ; ++i)
            { for (int j = 0 ; j < npp ; ++j)
                  { for (int k = 0; k < 27 ; ++k)
                        {if(rank==0)cout << "A[k] = "  << MolPol[0].dTij[i*npp*27 + j*27 + k] << endl; 
                        } 
                   
                  }
            }
        exit(1) ;
*/
	if (verbose > 1) {
	  if(rank==0)cout << "List of Polarizable sites, and the external (point-charge) field\n"; 
	  for (int i = 0; i < MolPol[0].nAtoms; ++i)
	    if(rank==0)printf("%3i at Site %i alpha=%9.4f  R=(%f, %f, %f)  Epc=(%f, %f, %f)\n", 
		   i, MolPol[0].SiteList[i], MolPol[0].Alpha[i],
		   Rpps[3*i+0],  Rpps[3*i+1], Rpps[3*i+2],
		   MolPol[0].Epc[3*i+0],  MolPol[0].Epc[3*i+1],  MolPol[0].Epc[3*i+2]);
	}
      }
      break;

     case 4:
     {

	if(rank==0)cout << "\nSelf-consistent water-water and water-electron polarizabilities\n"
	     << "with one polarizable site per water\n";
	
	npp = nWater;
	// treat as one big MolPol
	if (nMolPol != 1) {
	  if (nMolPol > 0)
	    delete[] MolPol;
	  nMolPol = 1;
	  MolPol = new DistributedPolarizabilities[nMolPol];
	}
	// copy data from arguments into MolPol[0]
	MolPol[0].nAtoms = npp;
	MolPol[0].SiteList.resize(npp);
	MolPol[0].Alpha.resize(npp);
	MolPol[0].InvA.resize(3*npp*3*npp);
	MolPol[0].Epc.resize(3*npp);
	dVec Rpps; Rpps.resize(3*npp);  // needed for calling ComputeInvA: list with coors of pol-sites only
	for (int i = 0; i < npp; ++i) {
	  int msite = 4*i+3;
	  MolPol[0].SiteList[i] = msite;  // M sites 
	  MolPol[0].Alpha[i] = Alpha;
	  Rpps[3*i+0] = r[3*msite+0];
	  Rpps[3*i+1] = r[3*msite+1];
	  Rpps[3*i+2] = r[3*msite+2];
	  //if(rank==0)cout << "M-site???? " << i << " (" << Rpps[3*i+0] * Bohr2Angs << ", " << Rpps[3*i+1] * Bohr2Angs << ", "  << Rpps[3*i+2] * Bohr2Angs << ")\n"; 
	}

	// compute Fields of point charges at M-sites
	// positions of the M-sites are in Rpps
//	double GroupCutoff = 2.0; // cutoff for point charges polarizing point polarizabilities  (M-site H distance is 1.56 Bohr)
//        double tholeDC_damp = 0.23 ; 
//	if(rank==0)cout << "GroupCutOff = " << GroupCutoff << "\n";
//	if(rank==0)cout << "Alpha = " << Alpha << "\n";
//	for (int i = 0; i < npp; ++i) {
//	  if(rank==0)cout << "Field at polarizable Site " << i << ", which is site " << MolPol[0].SiteList[i] << "\n";
//	  MolPol[0].Epc[3*i+0] = 0.0;
//	  MolPol[0].Epc[3*i+1] = 0.0;
//	  MolPol[0].Epc[3*i+2] = 0.0;
//	  for (int jq = 0; jq < nq; ++jq) {
//	    int iqsite = iq[jq];
//	    double qj = q[jq];
//	    if(rank==0)cout << "  Charge " << jq << " (q=" << qj << ") at Site " << iqsite;
//	    double Rx = Rpps[3*i+0] - r[3*iqsite+0];  // vector from point charge to polarizable site
//	    double Ry = Rpps[3*i+1] - r[3*iqsite+1];
//	    double Rz = Rpps[3*i+2] - r[3*iqsite+2];
//	    double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
//	    if(rank==0)cout << " is at distance of " << distance << "\n";
//	    if (distance > GroupCutoff) {
//	      void FieldOfCharge(const double q, const double *rq, const double *r, double *e, 
//                                const double damp, double alpha_Mol);
//	      double field[3];
//	      FieldOfCharge(qj, &r[3*iqsite], &Rpps[3*i], field, tholeDC_damp , Alpha);
//	      MolPol[0].Epc[3*i+0] += field[0];
//	      MolPol[0].Epc[3*i+1] += field[1];
//	      MolPol[0].Epc[3*i+2] += field[2];
//	      if(rank==0)cout << "        making a field of (" << field[0] << ", " << field[1] << ", " << field[2] << ")\n";
//	    }	
//	  }
//	}
	
	void ComputeInvA(int nSites, const double *R, const double *alpha, double aThole, double *InvA, int verbose); // defined in PolAux.cpp
	// so far aThole=0.3  this is taken from DPP
	ComputeInvA(npp, &Rpps[0], &MolPol[0].Alpha[0], 0.3, &MolPol[0].InvA[0], 1);       

	// compute the Water-Water induced dipoles in Rpps
	// and then the Water-Water polarization energy as an energy offset 
	// compute dipoles induced by Efield
	double dzero = 0.0;
	double done = 1.0;
	int one = 1;
	int n = 3*npp;
//	dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &epc[0], &one, &dzero, &Rpps[0], &one);
//	VpolWaterWater = -0.5 * ddot(&n, &epc[0], &one, &Rpps[0], &one);
	if(rank==0)cout << "Water-Water polarization energy = " << VpolWaterWater << "\n";

	if (verbose > 1) {
	  if(rank==0)cout << "List of Polarizable sites, and the external (point-charge) field\n"; 
	  for (int i = 0; i < MolPol[0].nAtoms; ++i)
	    if(rank==0)printf("%3i at Site %i alpha=%9.4f  Mu=(%f, %f, %f)  Epc=(%f, %f, %f)\n", 
		   i, MolPol[0].SiteList[i], MolPol[0].Alpha[i],
		   Rpps[3*i+0],  Rpps[3*i+1], Rpps[3*i+2],
		   MolPol[0].Epc[3*i+0],  MolPol[0].Epc[3*i+1],  MolPol[0].Epc[3*i+2]);
	}
      }
      break;
    default:
      if(rank==0)cout << " (No polarization potential used.)\n";
    }
    
  // repulsive potential
  switch (PotFlags[0])
    {
    case 1:
      Set6GTORepCore(nWater);
      RepCoreType = 1;
      break;
    case 2:
      Set12GTORepCore(nWater);
      RepCoreType = 1;
      break;
    case 3:
      Set4STORepCore(nWater);
      RepCoreType = 2;
      break;
    case 4:
      if(rank==0)cout << "No explicit repulsive potential to set up.\n";
      RepCoreType = 0;
      break;
    default:
      if(rank==0)cout << " Potential::SetupDPPGTOP PotFlags[0] this should never happen.\n"; exit(1);
    }

}



// Here for historical reasons (should not be used any more)
// 6S repulsive potential, 3S from STO3G on each H
// the coefficients come supposedly from a Schnitker-Rossky-like procedure
void Potential::Set6GTORepCore(int nWater)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  nGauss = 6 * nWater;
  Gauss.resize(2*nGauss);
  GaussSite.resize(nGauss);
  for (int i = 0; i < nWater; ++i) {
    int i6 = 6*i;
    int i12 = 12*i;
    GaussSite[i6+0] = 4*i + 1;  // Ha site
    GaussSite[i6+1] = 4*i + 1;   
    GaussSite[i6+2] = 4*i + 1; 
    GaussSite[i6+3] = 4*i + 2;  // Hb site
    GaussSite[i6+4] = 4*i + 2;  
    GaussSite[i6+5] = 4*i + 2; 
    Gauss[i12+0] = 33.87;            // exponent  
    Gauss[i12+1] = 2.84230876640242; // Seymour's coefficient
    Gauss[i12+2] = 5.095;            // exponent
    Gauss[i12+3] = 4.55401933663180; // Seymour's coefficient
    Gauss[i12+4] = 1.159;            // exponent
    Gauss[i12+5] = 6.91243747981459; // Seymour's coefficient
    Gauss[i12+6] = 33.87;            // exponent  
    Gauss[i12+7] = 2.84230876640242; // Seymour's coefficient
    Gauss[i12+8] = 5.095;            // exponent
    Gauss[i12+9] = 4.55401933663180; // Seymour's coefficient
    Gauss[i12+10]= 1.159;            // exponent
    Gauss[i12+11]= 6.91243747981459; // Seymour's coefficient
    // scale all by common RepCore scaling factor
    double sigma = 1.0 / VRepScale;
    for (int j = 1; j < 12; j += 2){
      Gauss[i12+j] *= sigma;
      //   if(rank==0)cout<<" Gauss["<< i12 + j<<" ] = "<< Gauss[i12+j]<<endl; 
    }
  }
}

// 12S repulsive potential
// 6S for O and 3S for H (exponents from 3-21G basis)
// the coefficients come from the Schnitker-Rossky procedure (coefficients for p orbitals are ignored)
void Potential::Set12GTORepCore(int nWater)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  const int nGpW = 12;
  nGauss = nGpW * nWater;
  GaussSite.resize(nGauss);
  Gauss.resize(2*nGauss);  // two parameters per Gaussian (exp and coefficient)        
  double sigmaH = 1.0 / VRepScale;  // scaling factor for the repulsive core, or just for H
  double sigmaO = sigmaH;
  if (SigmaOFlag)
    sigmaO = 1.0 / VRepScaleO; // separate scaling factor for O
  for (int i = 0; i < nWater; ++i) {
    int isite = nGpW * i;
    int ipara = 2 * nGpW * i;
    GaussSite[isite+0] = 4*i + 0;  // O site
    GaussSite[isite+1] = 4*i + 0;   
    GaussSite[isite+2] = 4*i + 0; 
    GaussSite[isite+3] = 4*i + 0;
    GaussSite[isite+4] = 4*i + 0;  
    GaussSite[isite+5] = 4*i + 0; 
    GaussSite[isite+6] = 4*i + 1;  // Ha site
    GaussSite[isite+7] = 4*i + 1;   
    GaussSite[isite+8] = 4*i + 1; 
    GaussSite[isite+9] = 4*i + 2;  // Hb site
    GaussSite[isite+10]= 4*i + 2;  
    GaussSite[isite+11]= 4*i + 2; 
    Gauss[ipara+0]  = 322.037;              // exponent  
    Gauss[ipara+1]  = 41.69216064 * sigmaO; // coefficient
    Gauss[ipara+2]  = 48.4308;              // exponent
    Gauss[ipara+3]  = 58.35946922 * sigmaO; // coefficient
    Gauss[ipara+4]  = 10.4206;              // exponent
    Gauss[ipara+5]  = 35.19878892 * sigmaO; // coefficient
    Gauss[ipara+6]  = 7.40294;              // exponent  
    Gauss[ipara+7]  = -3.27944113 * sigmaO; // coefficient
    Gauss[ipara+8]  = 1.5762;               // exponent
    Gauss[ipara+9]  = 4.44843633 * sigmaO;  // coefficient
    Gauss[ipara+10] = 0.373684;             // exponent  (last on O)
    Gauss[ipara+11] = 2.36213092 * sigmaO;  // coefficient
    Gauss[ipara+12] = 5.447178;             // exponent  (first on Ha)
    Gauss[ipara+13] = 0.69233452 * sigmaH;  // coefficient
    Gauss[ipara+14] = 0.824547;             // exponent
    Gauss[ipara+15] = 0.66686241 * sigmaH;  // coefficient
    Gauss[ipara+16] = 0.183192;             // exponent
    Gauss[ipara+17] = 0.07193908 * sigmaH;  // coefficient
    Gauss[ipara+18] = 5.447178;             // exponent  (first on Hb)
    Gauss[ipara+19] = 0.69233452 * sigmaH;  // coefficient
    Gauss[ipara+20] = 0.824547;             // exponent
    Gauss[ipara+21] = 0.66686241 * sigmaH;  // coefficient
    Gauss[ipara+22] = 0.183192;             // exponent
    Gauss[ipara+23] = 0.07193908 * sigmaH;  // coefficient
  }
}

// 4 STOs made from fitting the 12S repulsive potential
// 2S for O and 1S for H 
// the coefficients of the original GTOs come from the Schnitker-Rossky procedure
void Potential::Set4STORepCore(int nWater)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  const int nGpW = 4;   // 4 s-type Slaters per water
  nGauss = nGpW * nWater;
  GaussSite.resize(nGauss);
  Gauss.resize(2*nGauss);  // two parameters per Slater (exp and coefficient)        
  double sigmaH = 1.0 / VRepScale;  // scaling factor for the repulsive core, or just for H
  double sigmaO = sigmaH;
  if (SigmaOFlag)
    sigmaO = 1.0 / VRepScaleO; // separate scaling factor for O
  for (int i = 0; i < nWater; ++i) {
    int isite = nGpW * i;
    int ipara = 2 * nGpW * i;
    GaussSite[isite+0] = 4*i + 0;  // O site
    GaussSite[isite+1] = 4*i + 0;   
    GaussSite[isite+2] = 4*i + 1;  // Ha site
    GaussSite[isite+3] = 4*i + 2;  // Hb site
    Gauss[ipara+0]  = 9.98526;           // exponent     O1 
    Gauss[ipara+1]  = 192.608 * sigmaO;  // coefficient
    Gauss[ipara+2]  = 1.67075;           // exponent     O2
    Gauss[ipara+3]  = 13.8937 * sigmaO;  // coefficient
    Gauss[ipara+4]  = 1.64121;           // exponent     H
    Gauss[ipara+5]  = 1.7643 * sigmaH;   // coefficient
    Gauss[ipara+6]  = 1.64121;           // exponent     H
    Gauss[ipara+7]  = 1.7643 * sigmaH;   // coefficient
  }
}



////////////////////
//
//  this works so far only for DPPnSP 
//
void Potential::UpdateParameters(const double *PotPara)
{ 
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  if (PotFlags[0] == 1 || PotFlags[0] == 2) 
    {
      if (verbose > 0) {
	if(rank==0)cout << "New parameters\n";
	if(rank==0)cout << "PotPara 0: " << PotPara[0] << " (alpha)\n";
	if(rank==0)cout << "PotPara 1: " << PotPara[1] << " (b: damping length for polarizable sites [Bohr])\n";
	if(rank==0)cout << "PotPara 2: " << PotPara[2] << " (sigma for Vrep scaling)\n";
	if(rank==0)cout << "PotPara 3: " << PotPara[3] << " (Coulomb damping length [Bohr])\n";
	if(rank==0)cout << "PotPara 4: " << PotPara[4] << " (Dipole damping length [Bohr])\n";
      }

      Alpha = PotPara[0];
      PolDamping = fabs(PotPara[1]);  // use fabs to keep Powell from doing stupid stuff
      VRepScale = fabs(PotPara[2]);
      ChargeDamping = 1.0 / (PotPara[3]*PotPara[3]);   // Gaussian damping width (default may still be 0.2 Bohr)
      DipoleDamping = 1.0 / (PotPara[4]*PotPara[4]);   // Gaussian damping width (default may still be 0.4 Bohr)
      
      // update polarization potential
      for (int i = 0; i < nPPS; ++i)
	PPS[i] = Alpha;

      // update repulsive cores:
      if (PotFlags[0] == 1) {
	int nWater = nGauss / 6;
	double sigma = 1.0 / VRepScale;
	for (int i = 0; i < nWater; ++i) {
	  int i12 = 12*i;
	  Gauss[i12+1] = 2.84230876640242; // Seymour's coefficient
	  Gauss[i12+3] = 4.55401933663180; // Seymour's coefficient
	  Gauss[i12+5] = 6.91243747981459; // Seymour's coefficient
	  Gauss[i12+7] = 2.84230876640242; // Seymour's coefficient
	  Gauss[i12+9] = 4.55401933663180; // Seymour's coefficient
	  Gauss[i12+11]= 6.91243747981459; // Seymour's coefficient
	  // scale all by common RepCore scaling factor
	  for (int j = 1; j < 12; j += 2){
	    Gauss[i12+j] *= sigma;
	  }
	}
      }

      if (PotFlags[0] == 2) {
	const int nGpW = 12;   // 12 s-type Gaussians per water
	int nWater = nGauss / nGpW;
	double sigma = 1.0 / VRepScale; // overall scaling factor for Vrep
	for (int i = 0; i < nWater; ++i) {
	  int ipara = 2 * nGpW * i;
	  Gauss[ipara+1]  = 41.69216064; // coefficient
	  Gauss[ipara+3]  = 58.35946922; // coefficient
	  Gauss[ipara+5]  = 35.19878892; // coefficient
	  Gauss[ipara+7]  = -3.27944113; // coefficient
	  Gauss[ipara+9]  = 4.44843633;  // coefficient
	  Gauss[ipara+11] = 2.36213092;  // coefficient
	  Gauss[ipara+13] = 0.69233452;  // coefficient
	  Gauss[ipara+15] = 0.66686241;  // coefficient
	  Gauss[ipara+17] = 0.07193908;  // coefficient
	  Gauss[ipara+19] = 0.69233452;  // coefficient
	  Gauss[ipara+21] = 0.66686241;  // coefficient
	  Gauss[ipara+23] = 0.07193908;  // coefficient
	  // scale all coefficients by common RepCore scaling factor
	  for (int j = 1; j < 2*nGpW; j += 2){
	    Gauss[ipara+j] *= sigma;
	  }
	}
      }
    }
  else {
    if(rank==0)cout << "UpdateParameters: So far Potential " << PotFlags[0] << " does not update\n";
    exit(1);
  }
}


void Potential::UpdateChargeDipole(int nr, double *r, int nq, double *q, int *iq, int nd, double *d, int *id)
{ 

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if (nSites != nr || nCharges != nq || nDipoles != nd) {
    if(rank==0)cout << "You should not be calling UpdateChargeDipole. Call SetupXXX instead.\n";
    exit(1);
  }

   SetSites(nr, r);
   SetCharges(nq, q, iq);
   SetDipoles(nd, d, id);

}




//


/////////////////////////////////////////////////////////////////////////////////////
//
//  this evalues the DPP + Schnitker-Rossky motivated potential for the electron.
//  Sites can have charges and point dipoles (both damped), polarizabilities, and 
//  a repulsive potential that is expressed as a sum over GTOs (s-type only)  
//  it can also be used to evaluate other DPP6SP and DPP12SP as only the number of GTOs changes
//
double Potential::EvaluateDPPGTOP(const double *x)
{

  // progress_timer tmr("EvaluateDPPGTOP:", verbose);

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


 //  cout<<"x ="<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;

  double Vpc, Vind,  Vrep, Vpol; 
 if (PolType != 5 ) {
  //  potential due to point charges
   Vpc =0;
    Vpc = EvaluateChargePotential(x, DampType, ChargeDamping);
  //   cout<<"Vpc = "<<Vpc<<endl;
  
  //  potential due to point dipoles
   Vind=0;
  Vind = EvaluateDipolePotential(x, DampType, DipoleDamping);
 //    cout<<"Vind = "<<Vind<<endl;

  //  repulsive potential expressed as a sum over s-type GTOs or STOs
   Vrep = 0;
  if (RepCoreType > 0)
    Vrep = EvaluateRepulsivePotential(x, RepCoreType);
  //   cout<<"Vrep = "<<Vrep<<endl;
 }
  //  polarization potential
   Vpol = 0; 
  switch (PolType)
    {
    case 0:
      break; // no polarization potential in the model
    case 1:  // non-interacting point-polarizable sites
      Vpol = EvaluatePolPot(x, DampType, PolDamping);
      if (AdiabaticPolPot)
	Vpol = EpsDrude - sqrt(EpsDrude*(EpsDrude + Vpol));  // adiabatic potential of Drude oscillators
      else
	Vpol *= -0.5;  // this is the usual -0.5a/R**4, the 1st-order expansion of the adiabatic potential for large R
      break;
    case 2:  // non-interacting molecules with interacting atomic point-polarizabilities 
      Vpol = EvaluateMolecularPolarizableSites(x, DampType, PolDamping);
      if (AdiabaticPolPot)
	Vpol = EpsDrude - sqrt(EpsDrude*(EpsDrude + Vpol));  // adiabatic potential of Drude oscillators
      else
	Vpol *= -0.5;  // this is the usual -0.5a/R**4, the 1st-order expansion of the adiabatic potential for large R
      break;
    case 3:  // fully interacting atomic polarizabilities
      Vpol = SelfConsistentPolarizability(x, DampType, PolDamping);
    //   cout<<"Vpol = "<<Vpol<<endl;
      break;
    case 4:  // case 3: distributed, case 4: single sites
//      Vpol = SelfConsistentPolarizability(x, DampType, PolDamping);
//      break;
    case 5:  // fully interacting atomic polarizabilities
      Vpol = SelfConsistentPolarizability(x, DampType, PolDamping);
      Vind = 0.0;
      break;
    case 6:  // fully interacting atomic polarizabilities
      Vpol = SelfConsistentPolarizability(x, DampType, PolDamping);
      Vind = 0.0;
      break;
    case 33:
      Vpol = SelfConsistentPolarizability_33(x, DampType, PolDamping) ;
      break; 
    default:
      if(rank==0)cout << "EvaluateDPPGTOP: PolType = " << PolType << ", this should not happen.\n";
      exit(1);
    }


 if (x[0] == 0 && x[2] == 0.0) {
  cout<<" Y: Vpc, Vrep, Vind, Vpol = "<<x[1]<<" "<<Vpc<<" "<<Vrep<<" "<<Vind<<" "<<Vpol<<" "<<endl;
 }

/*
 if (x[0] == 0 && x[1] == 0.0) {
  cout<<" Z: Vpc, Vrep, Vind, Vpol = "<<x[2]<<" "<<Vpc<<" "<<Vrep<<" "<<Vind<<" "<<Vpol<<" "<<endl;
 }
*/


// TESTING E-WATER ELECTROSTATICS   
//     Vpc = 0.00000000 ;    //
//     Vrep = 0.00000000 ;     //
//     Vind = 0.00000000 ;    //
//     Vpol = 0.00000000 ;    //

  double Vtotal = Vpc + Vind + Vrep + Vpol;
// if(rank==0)cout<<"Vpc ="<<Vpc<<endl;
//  if(rank==0)cout<<"Vind ="<<Vind<<endl;
// if(rank==0)cout<<"Vrep ="<<Vrep<<endl;
// if(rank==0)cout<<"Vpol ="<<Vpol<<endl;

    ReturnEnergies[0] = Vpc;
    ReturnEnergies[1] = Vind;
    ReturnEnergies[2] = Vrep;
    ReturnEnergies[3] = Vpol;
    ReturnEnergies[4] = Vtotal;

  if (verbose > 0) {
    double v[5] = {Vpc, Vind, Vrep, Vpol, Vtotal};
    CheckMinMax(v);
  }

  if (verbose > 1)
    if(rank==0)printf("%10.6f %10.6f %10.6f %14.6e %14.6e %14.6e %14.6e %14.6e\n", 
	   x[0], x[1], x[2], Vpc, Vind, Vrep, Vpol, Vtotal);
  
  return Vtotal;
}








/////////////////////////////////////////////////////////////////////////////////////////////
//
//  A DPP based potential, that is, the Sites are O1, H1a, H1b, M1, O2, H2a, H2b, M2, ....
//
//  the charges and induced dipoles are provided by the DPP module
//  for the neutral cluster
//  polarizable sites with alpha = PotPara[0] and damping PotPara[1] are put on the M-sites
//  The repulsive potential has the form suggested by Turi and Borgis (JCP 117, 6186 (2002))
//  PotPara[2] : A1O
//  PotPara[3] : A1H  (not used)
//  PotPara[4] : B1O
//  PotPara[5] : B1H
//  PotPara[6] : B2O
//  PotPara[7] : B2H
//  PotPara[8] : B3O
//  PotPara[9] : B3H
//
void Potential::SetupDPPTB(int nr, const double *r,                         // site coordinates
			   int nq, const double *q, const int *iq,          // charges and their sites
			   int nd, const double *d, const int *id,          // dipoles and their sites
			   const double *DmuByDR,                           // devivatives 
			   const double *PotPara)                           // Parameters 
{ 
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  PotFlags.resize(1);
  PotFlags[0] = 11;

  if (verbose > 0) {
    if(rank==0)cout << "DPPTB potential: DPP charges plus Polarization plus Coulomb damping and Vrep from Turi \n";
    if(rank==0)cout << "PotPara 0: " << PotPara[0] << " (alpha)\n";
    if(rank==0)cout << "PotPara 1: " << PotPara[1] << " (b for polarization damping)\n";
    if(rank==0)cout << "PotPara 2: " << PotPara[2] << " (A1 parameter for Coulomb damping for O)\n";
    if(rank==0)cout << "PotPara 3: " << PotPara[3] << " (A1 parameter for Coulomb damping for H (not used))\n";
    if(rank==0)cout << "PotPara 4: " << PotPara[4] << " (B1O parameter of Vrep)\n";
    if(rank==0)cout << "PotPara 5: " << PotPara[5] << " (B1H parameter of Vrep)\n";
    if(rank==0)cout << "PotPara 6: " << PotPara[6] << " (B2O parameter of Vrep)\n";
    if(rank==0)cout << "PotPara 7: " << PotPara[7] << " (B2H parameter of Vrep)\n";
    if(rank==0)cout << "PotPara 8: " << PotPara[8] << " (B3O parameter of Vrep)\n";
    if(rank==0)cout << "PotPara 9: " << PotPara[9] << " (B3H parameter of Vrep)\n";
  }

  ChargeDamping = PotPara[2];  // use only one A for both H and O
  DipoleDamping = 0.5 * PotPara[2];
  PolDamping = PotPara[1];

  int nWater = nr / 4;
  if (nr % 4 != 0) 
    {if(rank==0)cout << "SetupDPPTB: no of sites " << nr << " is not a multiple of 4\n"; exit(1);}
  
  SetSites(nr, r);
  SetCharges(nq, q, iq);
  SetDipoles(nd, d, id);
  
  // polarization potential
  nPPS = nWater;
  PPS.resize(nPPS);
  PPSSite.resize(nPPS);
  for (int i = 0; i < nWater; ++i) {
    PPSSite[i] = 4*i+3;  // this is the i-th M-site
    PPS[i] = PotPara[0]; // alpha
  }

  // this is the Turi-Borgis repulsion, a combination of two error functions
  // at each center we need three parameters B1, B2, and B3
  // Vrep = B1/r * (erf(B2*r) - erf(B3*r))
  // this is put on H and O
  
  nGauss = 3 * nWater;
  Gauss.resize(3*nGauss);
  GaussSite.resize(nGauss);

  double B1O = PotPara[4];
  double B1H = PotPara[5];
  double B2O = PotPara[6];
  double B2H = PotPara[7];
  double B3O = PotPara[8];
  double B3H = PotPara[9];

  for (int i = 0; i < nWater; ++i) {
    int i3 = 3*i;
    int i9 = 9*i;
    GaussSite[i3+0] = 4*i + 0;  // O site
    GaussSite[i3+1] = 4*i + 1;  // Ha site
    GaussSite[i3+2] = 4*i + 2;  // Hb site
    Gauss[i9+0] = B1O;
    Gauss[i9+1] = B2O;
    Gauss[i9+2] = B3O;
    Gauss[i9+3] = B1H;
    Gauss[i9+4] = B2H;
    Gauss[i9+5] = B3H;
    Gauss[i9+6] = B1H;
    Gauss[i9+7] = B2H;
    Gauss[i9+8] = B3H;
  }
}


//
//  this evalues the DPP + Turi-Borgis motivated poential for the electron
//  Sites can have charges and point dipoles (both damped), polarizabilities, and 
//  a repulsive potential that is expressed as a sum over GTOs (s-type only so far)  
//
double Potential::EvaluateDPPTB(const double *x)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // this list has been computed in Evaluate()
  // compute a list of distances of all sites to the point r=(x,y,z)
  // for (int i = 0; i < nSites; ++i) {
  //   Rx[i] = x[0] - Site[3*i+0];
  //   Ry[i] = x[1] - Site[3*i+1];
  //   Rz[i] = x[2] - Site[3*i+2];
  //   R2[i] = Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i];
  //   R[i] = sqrt(R2[i]);
  //   Rminus3[i] = 1.0 / (R2[i] * R[i]);
  // }

  //
  //  potential due to point charges (damped with an error function)
  //  ChargeDamping is Turi's parameter A1, but here we have only one parameter
  //  instead of two for O and H separately
  //
  double Vpc = 0;
  for (int i = 0; i < nCharges; ++i) {
    int iqs = ChargeSite[i];
    Vpc += -Charge[i]/R[iqs] * erf(ChargeDamping * R2[iqs]);
  }


  //
  //  potential due to point dipoles (damped with an error function)
  //
  double Vind = 0;
  for (int i = 0; i < nDipoles; ++i) {
    int ids = DipoleSite[i];
    double *mu = &Dipole[3*i];
    double damp = erf(DipoleDamping * R2[ids]);
    Vind += -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2]) * Rminus3[ids] * damp * damp;
  }

  //
  //  repulsive potential from Turi
  //
  double Vrep = 0;
  for (int i = 0; i < nGauss; ++i) {
    int igs = GaussSite[i];
    double B1 = Gauss[3*i];
    double B2 = Gauss[3*i+1];
    double B3 = Gauss[3*i+2];
    Vrep +=  B1 / R[igs] * (erf(B2*R[igs]) - erf(B3*R[igs]));
  }


  //
  //  polarization potential: local, simple, Drude-like damping
  //
  //              alpha *e*e
  //  V(r,R) = - ------------ ( 1 - exp(-b(|r-R|^2) )^2
  //              2 |r-R|^4
  //
  double Vpol = 0;
  for (int i = 0; i < nPPS; ++i) {
    int ips = PPSSite[i];
    double alpha = PPS[i];
    double bdamping = PolDamping;
    double damp = 1.0 - exp(-bdamping*R2[ips]);
    Vpol += alpha * damp * damp / (R2[ips] * R2[ips]);
  }
  Vpol *= -0.5;


  //
  // done
  //

  double Vtotal = Vpc + Vind + Vrep + Vpol;
  if (verbose > 0) {
    double v[5] = {Vpc, Vind, Vrep, Vpol, Vtotal};
    CheckMinMax(v);
  }

  if (verbose > 15)
    if(rank==0)printf("%10.6f %10.6f %10.6f %14.6e %14.6e %14.6e %14.6e %14.6e\n", 
	   x[0], x[1], x[2], Vpc, Vind, Vrep, Vpol, Vtotal);

  return Vtotal;
}





/////////////////////////////////////////////////////////
//
//  evaluate potential due to point charges in the lists:
//  Charge[i], ChargeSite
//
//  DampFlag = 1 :  damping with 1 - GTO; DampParameter is Gaussian exponent
//           = 2 :  effective distance damping; DampParameter is cut-off distance  
//
double Potential::EvaluateChargePotential(const double *x, int DampFlag, double DampParameter)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

 //  cout<<"x ="<<x[0]<<" "<<x[1]<<" "<<x[2]<<endl;
   //progress_timer tmr("EvaluateChargePotential:", verbose); 
   double Vpc = 0;
   switch (DampFlag)
     {
     case 1:
       {
	 for (int i = 0; i < nCharges; ++i) {
	   int iqs = ChargeSite[i];
	   Vpc += -Charge[i]/R[iqs] * (1.0 - exp(-DampParameter * R2[iqs]));
	 }
       }
       break;
     case 2:
       {
	 for (int i = 0; i < nCharges; ++i) {
	   double Reff = R[ChargeSite[i]];
	   if (Reff < DampParameter) {
	     double ror0 = Reff / DampParameter;
	     Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
	     // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
	   }
	   Vpc += -Charge[i]/Reff;
      //    if (x[0] == 0 && x[1] == 0.0 && x[2]==10) { 
      //     cout<<"Charge["<<i<<"] = "<<Charge[i]<<", Reff = "<<Reff<<"Vpc = "<<Vpc<<endl;
      //   }
          // cout<<"Vpc = "<<Vpc<<endl;
	 }
       }
       break;
     default:
       if(rank==0)cout << "Potential::EvaluateChargePotential, DampFlag=" << DampFlag << ", this should not happen.\n"; 
       exit(1); 
     }
   return Vpc;
}



/////////////////////////////////////////////////////////
//
//  evaluate potential due to the point dipoles in the lists:
//  Dipole[i], DipoleSite
//
//  DampFlag = 1 :  damping with (1 - GTO)^2; DampParameter is Gaussian exponent
//           = 2 :  effective distance damping; DampParameter is cut-off distance  
//
double Potential::EvaluateDipolePotential(const double *x, int DampFlag, double DampParameter)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  //progress_timer tmr("EvaluateDipolePotential:", verbose);

  double Vpd = 0;

 // if(rank==0)cout <<"DampFlag: "<<DampFlag<<endl;

   switch (DampFlag)
     {
     case 1:
       {
	 for (int i = 0; i < nDipoles; ++i) {
	   int ids = DipoleSite[i];
	   double *mu = &Dipole[3*i];
	   double damp = 1.0 - exp(-DampParameter * R2[ids]);
	   Vpd += -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2]) * Rminus3[ids] * damp * damp;
	 }
       }
       break;
     case 2:
       {
	 for (int i = 0; i < nDipoles; ++i) {
	   int ids = DipoleSite[i];
	   double Rm3 = Rminus3[ids];
	   double Reff = R[ids];
	   if (Reff < DampParameter) {
	     double ror0 = Reff / DampParameter;
	     Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
	     // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
	     Rm3 = 1.0 / (Reff * Reff * Reff);
	   }
	   double *mu = &Dipole[3*i];
	   Vpd += -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2]) * Rm3;
	 }

       }
       break;
     default:
       if(rank==0)cout << "Potential::EvaluateDipolePotential, DampFlag=" << DampFlag << ", this should not happen.\n"; 
       exit(1); 
     }
   return Vpd;

}

/////////////////////////////////////////////////////////
//
//  evaluate repulsive potential consisting of GTOs or STOs 
//  the parameters are in the lists Gauss[i], GaussSite[i]
//  Gauss[2*i]   is the i-th exponent
//  Gauss[2*i+1] is the i-th coefficient
//
//  RepCoreFlag = 1 :  use Gaussians
//              = 2 :  use Slaters
//
double Potential::EvaluateRepulsivePotential(const double *x, int RepCoreFlag)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  //progress_timer tmr("EvaluateRepulsivePotential:", verbose);

   double Vrep = 0;
   switch (RepCoreFlag)
     {
     case 1:
       {
	 for (int i = 0; i < nGauss; ++i) {
	   int igs = GaussSite[i];
	   Vrep +=  Gauss[2*i+1] * exp(-Gauss[2*i] * R2[igs]);
	 }
       }
       break;
     case 2:
       {
	 for (int i = 0; i < nGauss; ++i) {
	   int igs = GaussSite[i];
          // cout<<" Gauss[2*i+1] , Gauss[2*i] "<<Gauss[2*i+1]<<" "<<Gauss[2*i]<<endl;
	   Vrep +=  Gauss[2*i+1] * exp(-Gauss[2*i] * R[igs]);
	 }
       }
       break;
     default:
       if(rank==0)cout << "Potential::EvaluateRepulsivePotential: RepCoreFlag = " << RepCoreFlag << ", this should not happen\n"; 
       exit(1);
     }
   return Vrep;
}

////////////////////////////////////////////////////////////////////////////////
//
//  for evaluating a potential consisting of non-interacting point-polarizabilities
//
//  what is evaluated is Spol = Sum alpha/R^4 with appropriate damping functions
//  Spol can then be used to evaluate the adiabatic potential Vpol = eps - sqrt(eps^2 + eps*Spol)
//  or the more usual 1st-order expression Vpol = -0.5*Spol 
//
//  the parameters are in the lists PPS[i], PPSSite[i]
//  PPS[i] is the polarizability alpha
//
//  DampFlag = 1 :  damping with (1 - GTO)^2; DampParameter is Gaussian exponent
//           = 2 :  effective distance damping; DampParameter is cut-off distance  
//
double Potential::EvaluatePolPot(const double *x, int DampFlag, double DampParameter)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   double Spol = 0;
   switch (DampFlag)
     {
     case 1:
       {
	 for (int i = 0; i < nPPS; ++i) {
	   int ips = PPSSite[i];
	   double alpha = PPS[i];
	   double Rsq = R2[ips];
	   double damp = 1.0 - exp(-DampParameter*Rsq);
	   Spol += alpha * damp * damp / (Rsq * Rsq);
	 }
       }
       break;
     case 2:
       {       
	 for (int i = 0; i < nPPS; ++i) {
	   int ips = PPSSite[i];
	   double alpha = PPS[i];
	   double Reff = R[ips];
	   double Rsq = R2[ips];
	   if (Reff < DampParameter) {
	     double ror0 = Reff / DampParameter;
	     Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
	     // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
	     double Rsqeff = Reff * Reff;
             Spol += alpha * Rsq / (Rsqeff * Rsqeff * Rsqeff);
	   }
           else
	   Spol += alpha / (Rsq * Rsq);
	 }	 
       }
       break;
     default:
       if(rank==0)cout << "Potential::EvaluatePolPot, DampFlag=" << DampFlag << ", this should not happen.\n"; 
       exit(1); 
     }
   return Spol;
}




//////////////////////////////////////////////////////////////////////////////////////
//
//  Evaluate polarization of molecules with distribted atomic polarizabilities
//  this is still one-molecule-at-a-time and decoupled from any water-water polarization
//  but the polarizability has the correct anisotropy
//
//  the no of Molecules is nMolPol
//  each Molecules has its data in a struct DistributedPolarizabilities MolPol
//  
//
double Potential::EvaluateMolecularPolarizableSites(const double *x, int DampFlag, double DampParameter)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  double Spol = 0;

  static dVec Efield;
  static dVec mu;


  for (int imol = 0; imol < nMolPol; ++imol) 
    {
      // compute Efield at atomic sites
      int nAtoms = MolPol[imol].nAtoms;
      int n = 3*nAtoms;  // dimension of InvA and Efield, and mu
      Efield.resize(n);
      mu.resize(n);
      double gij;
      for (int i = 0; i < nAtoms; ++i) {
	int SiteIndex = MolPol[imol].SiteList[i];
	// Rij[k] = RSites[3*i+k] - Rel[k];
	// rij2 = Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2];
	// rij  = sqrt(rij2);
	// rij3 = R[SiteInex] * R2[SiteIndex];
	// fij = 1 - exp(-DampParameter * rij3);
	// gij = -fij/rij3;
	//
	// here -gij is computed, since (Rx, Ry, Rz) = -Rij
	gij = Rminus3[SiteIndex];
	switch (DampFlag)
	  {
	  case 1:
	    // Thole, Gauss damping as in the Drude code is not enough
	    gij *= (1.0 - exp(-DampParameter * R[SiteIndex] * R2[SiteIndex]));
	    if(rank==0)printf("%20.12f", gij);
	    break;
	  case 2:
	    // effective-r damping
	    {
	      double Reff = R[SiteIndex];
	      if (Reff < DampParameter) {
		double ror0 = Reff / DampParameter;
		Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
		// old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
		gij = 1.0 / (Reff * Reff * Reff);
	      }
	    }
	    break;
	  default:
	    if(rank==0)cout << "EvaluateMolecularPolarizableSites: DampFlag=" << DampFlag << ", this should never happen\n";
	    exit(1);
	  }
	Efield[3*i+0] = gij * Rx[SiteIndex];
	Efield[3*i+1] = gij * Ry[SiteIndex];
	Efield[3*i+2] = gij * Rz[SiteIndex];
      }
      // compute dipoles induced by Efield
      double dzero = 0.0;
      double done = 1.0;
      int one = 1;
      double *InvA = &(MolPol[imol].InvA[0]);

      dgemv("N", &n, &n, &done, InvA, &n, &Efield[0], &one, &dzero, &mu[0], &one);
      // compute Efield * dipoles
      Spol += ddot(&n, &Efield[0], &one, &mu[0], &one);
            
      int debug = 0;
      if (debug > 0) {
	if(rank==0)printf("Electron is at (%f, %f, %f)\n", x[0], x[1], x[2]);
	for (int i = 0; i < nAtoms; ++i) {
	  int SiteIndex = MolPol[imol].SiteList[i];
	  if(rank==0)printf("Dipole induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
		 i, imol, SiteIndex, mu[3*i], mu[3*i+1], mu[3*i+2]);
	}
	double e = ddot(&n, &Efield[0], &one, &mu[0], &one);
	if(rank==0)printf("Epol = %f\n", e);
	if(rank==0)printf("Reversing field\n");
	double mo = -1.0;
	dscal(&n, &mo, &Efield[0], &one);


	dgemv("N", &n, &n, &done, InvA, &n, &Efield[0], &one, &dzero, &mu[0], &one);
	for (int i = 0; i < nAtoms; ++i) {
	  int SiteIndex = MolPol[imol].SiteList[i];
	  if(rank==0)printf("Dipole induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
		 i, imol, SiteIndex, mu[3*i], mu[3*i+1], mu[3*i+2]);
	}
	double e2 = ddot(&n, &Efield[0], &one, &mu[0], &one);
	if(rank==0)printf("Epol = %f   RelError = %e\n", e2, 2.0*fabs(e-e2)/(e+e2));
      
	if(rank==0)cout.flush();

      }
    } // endfor (int imol = 0; imol < nMolPol; ++imol)

  return Spol;
}
///////////////////////////////////////////////////////////////////////
//
//  full self-consistent polarizable sites
//  this is just like one molecular polarizability plus extra field
//
//
//
double Potential::SelfConsistentPolarizability(const double *x, int DampFlag, double DampParameter)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  dVec Efield;
  dVec mu;


   //progress_timer tmr("SelfConsistentPolarizability:", verbose);

  // compute electron's Efield at atomic sites
  int nAtoms = MolPol[0].nAtoms;
  int n = 3*nAtoms;  // dimension of InvA and Efield, and mu



  // if(rank==0)cout<< "nAtoms*3:" << n << endl;

  Efield.resize(n);
  mu.resize(n);
  double gij;
  for (int i = 0; i < nAtoms; ++i) {
    int SiteIndex = MolPol[0].SiteList[i];
    gij = Rminus3[SiteIndex];
    switch (DampFlag)
      {
      case 1:
	// cubic Thole, Gauss damping as in the Drude code is too weak
	gij *= 1.0 - exp(-DampParameter * R[SiteIndex] * R2[SiteIndex]);
	break;
      case 2:
	// effective-r damping
	{
	  double Reff = R[SiteIndex];
	  if (Reff < DampParameter) {
	    double ror0 = Reff / DampParameter;
	    Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
	    // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
	    gij = 1.0 / (Reff * Reff * Reff);
	  }
	}
	break;
      default:
	if(rank==0)cout << "SelfConsistentPolarizability: DampFlag=" << DampFlag << ", this should never happen\n";
	exit(1);
      }

    Efield[3*i+0] = gij * Rx[SiteIndex];
    Efield[3*i+1] = gij * Ry[SiteIndex];
    Efield[3*i+2] = gij * Rz[SiteIndex];

  }

  // add the point charges' field here
  // but use a blas call instead of this:

  for (int k = 0; k < n; ++k) {
    Efield[k] += MolPol[0].Epc[k];
   // Efield[k] = MolPol[0].Epc[k];
    //if(rank==0)cout<<"Efield["<<k<<"]="<<Efield[k]<<endl;
  }

// Self-consistent calculation for dipole moment --- Tae Hoon Choi

// if(rank==0)cout<<" Start sef-consistent calculation for dipole moment --- Tae Hoon Choi"<<endl;
//void SCF_inducedDipole(int nSites, const double *Tensor, const double *alpha, double aThole, double *Efield, double *mu);
//   SCF_inducedDipole(nAtoms, &MolPol[0].Tensor[0], &MolPol[0].Alpha[0], 0.3, &Efield[0], &mu[0]);

 // for (int i =0; i < n; ++i)
 //  if(rank==0)cout<<"first mu["<<i<<"i] ="<<mu[i]<<endl;

  //
  //
  // compute dipoles induced by Efield
  double dzero = 0.0;
  double done = 1.0;
  int one = 1;
  
 // if(rank==0)cout << "start dgemv in SCP --------"<<endl;

//  for (int i = 0; i < n; ++i)
//    for (int j = 0; j < n; ++j)
//    if(rank==0)cout<<"MolPol[0].InvA["<<n*i+j<<"] ="<<MolPol[0].InvA[n*i+j]<<endl;


//  dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
  dsymv("L", &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);

//   for (int i =0; i < n; ++i)
//    if(rank==0)cout<<"Second mu["<<i<<"i] ="<<mu[i]<<endl;

  // compute Efield * dipoles
  double Vpol = -0.5 * ddot(&n, &Efield[0], &one, &mu[0], &one) - VpolWaterWater;


//  double Vpol = -0.5 * ddot(&n, &Efield[0], &one, &mu[0], &one) ;
//  if(rank==0)cout<< "VpolWaterWater ="<<VpolWaterWater<<"\n";  



//  for (int k = 0; k < n; ++k)
//   if(rank==0)cout<<"in energy mu["<<k<<"]="<<mu[k]<<"\n";
  
  int debug = 0;
  if (debug > 0) {
    if(rank==0)printf("Electron is at (%f, %f, %f) Angs\n", x[0], x[1], x[2]);
    for (int i = 0; i < nAtoms; ++i) {
      int SiteIndex = MolPol[0].SiteList[i];
      if(rank==0)printf("Dipole, in Debye, induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
	     i, 0, SiteIndex, mu[3*i]*AU2Debye, mu[3*i+1]*AU2Debye, mu[3*i+2]*AU2Debye);
    }
    if(rank==0)printf("Epol = %f Hartree = %f eV\n", Vpol, Vpol*AU2EV);
    //if(rank==0)printf("Reversing field\n");
    //double mo = -1.0;
    //dscal(&n, &mo, &Efield[0], &one);
    //dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
    //for (int i = 0; i < nAtoms; ++i) {
    //  int SiteIndex = MolPol[0].SiteList[i];
    //  if(rank==0)printf("Dipole, in Debye, induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
    //	     i, 0, SiteIndex, mu[3*i]*AU2Debye, mu[3*i+1]*AU2Debye, mu[3*i+2]*AU2Debye);
    //}
    //double e2 = ddot(&n, &Efield[0], &one, &mu[0], &one);
    //if(rank==0)printf("Epol = %f   RelError = %e\n", e2, 2.0*fabs(e-e2)/(e+e2));
    
    //if(rank==0)cout.flush();
  }


  return Vpol;


}

///////////////////////////////////////////////////////////////////////
//
//  full self-consistent polarizable sites
//  this is just like one molecular polarizability plus extra field
//
//
//
double Potential::SelfConsistentPolarizability_33(const double *x, int DampFlag, double DampParameter)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  dVec Efield;
  dVec mu;



  // compute electron's Efield at atomic sites
  int nAtoms = MolPol[0].nAtoms;
  int n = 3*nAtoms;  // dimension of InvA and Efield, and mu
  Efield.resize(n);
  mu.resize(n);
  double gij;
  for (int i = 0; i < nAtoms; ++i) {
    int SiteIndex = MolPol[0].SiteList[i];
    gij = Rminus3[SiteIndex];
    switch (DampFlag)
      {
      case 1:
	// cubic Thole, Gauss damping as in the Drude code is too weak
	gij *= 1.0 - exp(-DampParameter * R[SiteIndex] * R2[SiteIndex]);
	break;
      case 2:
	// effective-r damping
	{
	  double Reff = R[SiteIndex];
	  if (Reff < DampParameter) {
	    double ror0 = Reff / DampParameter;
	    Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));  
	    // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
	    gij = 1.0 / (Reff * Reff * Reff);
	  }
	}
	break;
      default:
	if(rank==0)cout << "SelfConsistentPolarizability: DampFlag=" << DampFlag << ", this should never happen\n";
	exit(1);
      }
    Efield[3*i+0] = gij * Rx[SiteIndex];
    Efield[3*i+1] = gij * Ry[SiteIndex];
    Efield[3*i+2] = gij * Rz[SiteIndex];
  }

  // add the point charges' field here
  // but use a blas call instead of this:

//  for (int k = 0; k < n; ++k)
//    Efield[k] += MolPol[0].Epc[k];

  // compute dipoles induced by Efield
  double dzero = 0.0;
  double done = 1.0;
  int one = 1;
  dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
  // compute Efield * dipoles
// double Vpol = -0.5 * ddot(&n, &Efield[0], &one, &mu[0], &one) - VpolWaterWater;
  double Vpol = -0.5 * ddot(&n, &Efield[0], &one, &mu[0], &one) ;
//  if(rank==0)cout<< "VpolWaterWater ="<<VpolWaterWater<<"\n";



//  for (int k = 0; k < n; ++k)
//   if(rank==0)cout<<"in energy mu["<<k<<"]="<<mu[k]<<"\n";
  
  int debug = 0;
  if (debug > 0) {
    if(rank==0)printf("Electron is at (%f, %f, %f) Angs\n", x[0], x[1], x[2]);
    for (int i = 0; i < nAtoms; ++i) {
      int SiteIndex = MolPol[0].SiteList[i];
      if(rank==0)printf("Dipole, in Debye, induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
	     i, 0, SiteIndex, mu[3*i]*AU2Debye, mu[3*i+1]*AU2Debye, mu[3*i+2]*AU2Debye);
    }
    if(rank==0)printf("Epol = %f Hartree = %f eV\n", Vpol, Vpol*AU2EV);
    //if(rank==0)printf("Reversing field\n");
    //double mo = -1.0;
    //dscal(&n, &mo, &Efield[0], &one);
    //dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
    //for (int i = 0; i < nAtoms; ++i) {
    //  int SiteIndex = MolPol[0].SiteList[i];
    //  if(rank==0)printf("Dipole, in Debye, induced on Site %i of Molecule %i (SiteIndex = %i) is (%f, %f, %f)\n", 
    //	     i, 0, SiteIndex, mu[3*i]*AU2Debye, mu[3*i+1]*AU2Debye, mu[3*i+2]*AU2Debye);
    //}
    //double e2 = ddot(&n, &Efield[0], &one, &mu[0], &one);
    //if(rank==0)printf("Epol = %f   RelError = %e\n", e2, 2.0*fabs(e-e2)/(e+e2));
    
    //if(rank==0)cout.flush();
  }


  return Vpol;


}


void Potential::SetupMinMax()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   switch (PotFlags[0]) 
   {
   case 1:
   case 2:
   case 3:
   case 4:
   case 11:
      nContributions = 5;  // Vpc, Vind, Vrep, Vpol, Vtotal 
      break;
   case 101:
     nContributions = 1; // Vps only
     break;
   default:
      if(rank==0)cout << "SetupMinMax: unknown PotFlag " << PotFlags[0] << "\n"; 
      exit(1);
   }

   MaxV.resize(nContributions, -std::numeric_limits<double>::max());
   MinV.resize(nContributions, std::numeric_limits<double>::max());


}

void Potential::CheckMinMax(const double *v)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


   for (int i = 0; i < nContributions; ++i) {
      if (v[i] > MaxV[i]) 
         MaxV[i] = v[i]; 
      if (v[i] < MinV[i]) 
         MinV[i] = v[i]; 
   }
}

void Potential::PrintMinMax()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   if (nContributions > 0) {  // that means MinMax has been initialized at some point
      switch (PotFlags[0])
      {
      case 1:
      case 2:
      case 3:
      case 4:
         if (verbose > 0) {    
           if(rank==0)printf("Range of the potential energy at the grid points:\n");
           if(rank==0)printf("  Total potential          :  %10.3e  -  %10.3e\n", MinV[4], MaxV[4]);
           if(rank==0)printf("  Point charge potential   :  %10.3e  -  %10.3e\n", MinV[0], MaxV[0]);
           if(rank==0)printf("  Induced dipole potential :  %10.3e  -  %10.3e\n", MinV[1], MaxV[1]);
           if(rank==0)printf("  Repulsive potential      :  %10.3e  -  %10.3e\n", MinV[2], MaxV[2]);
           if(rank==0)printf("  Polarization potential   :  %10.3e  -  %10.3e\n", MinV[3], MaxV[3]);
         }
         break;
      default:
         if(rank==0)cout << "PrintMinMax: unknown PotFlag " << PotFlags[0] << "\n"; 
         exit(1);
      }
   }
}







//
// provide a list of sites (cartesian coordinates)
// other objects have no coordinates, but a reference index to the site list
//
void Potential::SetSites(int n, const double *sites)
{

   nSites = n;
   Site.resize(3*n);
   std::copy(sites, sites+3*n, Site.begin());

   Rx.resize(n);
   Ry.resize(n);
   Rz.resize(n);
   R.resize(n);
   R2.resize(n);
   Rminus3.resize(n);

   RVec.resize(n/4*3);
   X.resize(n);
   Y.resize(n);
   Z.resize(n);
}

//
//  list of charges, q, with reference to the Site list
//
void Potential::SetCharges(int n, const double *q, const int *iq)
{

   nCharges = n;
   Charge.resize(n);
   ChargeSite.resize(n);
   std::copy(q, q+n, Charge.begin());
   std::copy(iq, iq+n, ChargeSite.begin());
}

//
//  list of dipoles, mx, my, mz, with reference to the Site list
//
void Potential::SetDipoles(int n, const double *d, const int *id)
{

   nDipoles = n;
   Dipole.resize(3*n);
   DipoleSite.resize(n);
   std::copy(d, d+3*n, Dipole.begin());
   std::copy(id, id+n, DipoleSite.begin());
}


//
//  list of polarizabilities, alpha, damping, with reference to the Site list
//

void Potential::SetPPS(int n, const double *a, const int *ia)
{
   nPPS = n;
   PPS.resize(2*n);
   PPSSite.resize(n);
   std::copy( a, a+2*n, PPS.begin() );
   std::copy( ia, ia+n, PPSSite.begin());
}



//
//  list of s-type gaussians, i,.e., exponents and coefficient, with reference to the Site list
//
void Potential::SetGauss(int n, const double *expcoeff, const int *ig)
{
   nGauss = n;
   Gauss.resize(2*n);
   GaussSite.resize(n);
   std::copy(expcoeff, expcoeff+2*n, Gauss.begin());
   std::copy(ig, ig+n, GaussSite.begin());
}

void Potential::ReportEnergies(int n, double *e)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if (n != nReturnEnergies) {
    if(rank==0)cout << "Nonono in ReportEnergies.";
    exit(1);
  }
  for (int i = 0; i < nReturnEnergies; ++i) 
    e[i] = ReturnEnergies[i];

}

double Potential::getRtol()
{

  return Rtol; 

}

int Potential::getPolType()
{

  if ( PotFlags[0] == 1021) return PolFlagCsixty;
  else  return PolType; 

}


void FieldOfCharge(const double q, const double *rq, const double *r, double *e, double a, double alpha)
{
  double Rx = r[0] - rq[0]; // vector from point charge to polarizable site
  double Ry = r[1] - rq[1]; 
  double Rz = r[2] - rq[2];
  double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
  double distance3 = distance *distance*distance ;
  double tensor = q / distance3 ;
  double damp = 1.0 - exp(-a * distance3 / alpha ) ;
  
  e[0] = Rx * tensor * damp; 
  e[1] = Ry * tensor * damp; 
  e[2] = Rz * tensor * damp;
}
