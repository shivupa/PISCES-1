////////////////////////////////////////////////////////////////////////////////////
//
//  Here are functions of WaterCluster:: that compute the intermolecular potential
//  and those that deal with the input
//  we use the TTM/Albert water model and this will be pretty much hardwired here
//  nevertheless, should you want a different model, all you need to change is here
//
//  this replaces Seymour's function for the DC model that were mostly in Cluster::
//
/////////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <time.h>

#include <iomanip>
#include <fstream>
#include <iostream>


#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"


#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"
#include "DPP.h"

using namespace std;

//
//  functions defined and used (so far) only in this file
//

void watergeo(double ox, double oy, double oz, double h1x, double h1y, double h1z,
              double h2x, double h2y, double h2z);

void swap(double &a, double &b);




void FieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az,  
                   double &bx, double &by, double &bz, double &ex, double &ey, double &ez,double alphaA, double alphaB, double dampDD);
void dFieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az,
                    double &bx, double &by, double &bz, double &dxex, double &dxey, double &dxez,
                    double &dyex, double &dyey, double &dyez, double &dzex, double &dzey, double &dzez,
                    double alphaA, double alphaB, double dampDD);
void dFieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az,
                    double &bx, double &by, double &bz, double &Txx, double &Tyy, double &Tzz,
                    double &Txy, double &Txz, double &Tyz, double alphaA, double alphaB, double dampDD);
void InverseMatrix(double* Input, double* Output, int order, int type );

double WaterWaterElectrostaticEnergy(Water &W1, Water &W2, int verbose);
double WaterWaterElectrostaticGradient(Water &W1, Water &W2, double *Grad1, double *Grad2, int verbose);
double ElectrostaticDCDamp(double distance, double alphai, double alphaj, double a);
double ElectrostaticDDDamp(double distance, double alphai, double alphaj, double a);
double WaterWaterPolarizationEnergy(Water &W1, Water &W2, int verbose);
double WaterWaterPolarizationGradient(Water &W1, Water &W2, double *Grad1, double *Grad2, int firstCheck, int verbose);
double WaterWaterDispersionEnergy(Water &W1, Water &W2, int verbose);
double WaterWaterDispersionGradient(Water &W1, Water &W2, double *Grad1,double *Grad2, int verbose);


int WaterCluster::ReportNoOfWaters(void)
{
   return nwaters;
}

// all Sites (4*nwaters for DPP) 
int WaterCluster::ReportNoOfSitesPerMonomer(void)
{
   return nSitesPerMonomer;
}

///////////////////////////////////////////////////////
//
//  this returns the numbers appropriate for the DPP model
//  the electron sees only the vector sum of the induced dipoles, i.e., one dipole per monomer
//  
void WaterCluster::ReportNoOfSites(int &nSites, int &nCharges, int &nDipoles)
{
   nSites = nwaters * nSitesPerMonomer;
   nCharges = nwaters * nChargesPerMonomer;
   nDipoles = nwaters; 
}

///////////////////////////////////////////////////////
//
//  return no of point-polarizable sites
//  
int WaterCluster::ReportNoOfPolSites(void)
{
   return nwaters * 3;
}


double WaterCluster::ReportPolarization(void)
{
   return PolarizationEnergy; 
}

//////////////////////////////////////////////////////////////////////////
//
//  Initialize the WaterCluster with nW cartesian coordinates
//  unit == 1 coordinates are in Angs
//  unit == 2 coordinates are in Bohr 
//  CenterFlag == 1 changes to center of mass origin
//  KTCharges == 1 use DPP-like model with KT charges
//
void WaterCluster::SetStructure(int nW, const double WCoorIn[], int unit, int CenterFlag, int KTCharges, int verbose)
{

  KTFlag = KTCharges;

  if (verbose > 7) {
    cout << "Input coordinates of "<< nW << " waters in WaterCluster::SetStructure:\n";
    for (int i = 0; i < nW; i++) {
      int add = i*9;
      cout << "Water #" << i << endl;
      if (unit == 1)
	watergeo(WCoorIn[add+0], WCoorIn[add+1], WCoorIn[add+2],
		 WCoorIn[add+3], WCoorIn[add+4], WCoorIn[add+5],
		 WCoorIn[add+6], WCoorIn[add+7], WCoorIn[add+8]);
      else if (unit == 2)
	watergeo(WCoorIn[add+0] * Bohr2Angs, WCoorIn[add+1] * Bohr2Angs, WCoorIn[add+2] * Bohr2Angs,
		 WCoorIn[add+3] * Bohr2Angs, WCoorIn[add+4] * Bohr2Angs, WCoorIn[add+5] * Bohr2Angs,
		 WCoorIn[add+6] * Bohr2Angs, WCoorIn[add+7] * Bohr2Angs, WCoorIn[add+8] * Bohr2Angs);	
      else
	{cout << "Unknown unit = " << unit << "\n"; exit(1);}
    }
  }
  
  static dVec WCoor; 
  WCoor.resize(9*nW);
  std::copy(WCoorIn, &WCoorIn[9*nW], &WCoor[0]);
  if (CenterFlag == 1)
    PutIntoCenterOfMassSystem(nW, WCoor);
  
  //
  //  if the no of waters, nwaters, changes delete old arrays
  //  on first cal (and after change) allocate space for waters
  //
  if (nW != nwaters) {
    if (nwaters > 0) {
      for (int i = 0; i < nwaters; ++i)
	delete[] waters[i];
      delete[] waters;
    }
    nwaters = nW;
    waters = new Water*[nwaters];
    for (int i = 0; i < nwaters; ++i) {
      if (KTFlag == 1)
	waters[i] = new Water(0, 0, 0, 12, 1, 1.0); // KT charges, at bit of TTM and DPP 
      else
	waters[i] = new Water(0, 0, 0,  2, 1, 1.0); // standard DPP 
    }
  }
  
  //  the internal function to go from cartesian input to
  //  Configuration, and adjusted internal xyz coordinates is Water::Reposition
  //  I guess Seymour called it reposition as all Waters are initally created 
  //  at the origin in standard orientation
  //
  if (unit == 2) { 
    for (int i = 0; i < nwaters; i++) {
      int add = i*9;
      waters[i]->Reposition(WCoor[add+0], WCoor[add+1], WCoor[add+2],  // Bohrs
			    WCoor[add+3], WCoor[add+4], WCoor[add+5],
			    WCoor[add+6], WCoor[add+7], WCoor[add+8]);
    }
  }
  else if (unit == 1) {
    for (int i = 0; i < nwaters; i++) {
      int add = i*9;
      waters[i]->Reposition_A(WCoor[add+0], WCoor[add+1], WCoor[add+2], // Angstroms
			      WCoor[add+3], WCoor[add+4], WCoor[add+5],
			      WCoor[add+6], WCoor[add+7], WCoor[add+8]);
    }
  }
  else 
    {cout << "WaterCluster::SetStructure: Unit for coordinates " << unit << " unknown\n"; exit(1);}
  
  
  if (verbose > 7) {
    double *conf = new double[6*nwaters];
    Water *w;
    this->GetRigidConf(conf);
    cout << "Projected Water coordinates in WaterCluster::SetStructure:\n";
    for (int i = 0; i < nwaters; i++) {
      w = waters[i];
      cout << "Water #" << i << endl;
      watergeo(w->AtomCenters[0].x*Bohr2Angs, w->AtomCenters[0].y*Bohr2Angs, w->AtomCenters[0].z*Bohr2Angs,
	       w->AtomCenters[1].x*Bohr2Angs, w->AtomCenters[1].y*Bohr2Angs, w->AtomCenters[1].z*Bohr2Angs,
	       w->AtomCenters[2].x*Bohr2Angs, w->AtomCenters[2].y*Bohr2Angs, w->AtomCenters[2].z*Bohr2Angs);
      printf("  Center = %f, %f, %f\n", conf[6*i+0]*Bohr2Angs, conf[6*i+1]*Bohr2Angs, conf[6*i+2]*Bohr2Angs); 
      printf("  Orient = %f, %f, %f\n", conf[6*i+3], conf[6*i+4], conf[6*i+5]); 
    }
    delete[] conf;
  }
  
  // flag whether the induced dipoles have been computed
  DipoleSCFFlag = 0;
}


//////////////////////////////////////////////////////////////////////////
//
//  Report internal cartesian coordinates (Bohr) scaled with scale
//  nW is the maximal no of waters to return
//  the actual no of water molecules returned is returned 
//
int WaterCluster::GetStructure(int nW, double *WCoor, double scale, int verbose)
{
   int nwret = 0;

   if (nW >= nwaters)
      nwret = nwaters;
   else {
      nwret = nW; 
      cout << "Warning: not all waters are returned by GetStructure\n";
   }

   Water *w;
   double *coor = WCoor;
   for (int i = 0; i < nwret; i++) {
      w = waters[i];
      *coor = w->AtomCenters[0].x * scale; coor++;
      *coor = w->AtomCenters[0].y * scale; coor++;
      *coor = w->AtomCenters[0].z * scale; coor++;
      *coor = w->AtomCenters[1].x * scale; coor++;
      *coor = w->AtomCenters[1].y * scale; coor++;
      *coor = w->AtomCenters[1].z * scale; coor++;
      *coor = w->AtomCenters[2].x * scale; coor++;
      *coor = w->AtomCenters[2].y * scale; coor++;
      *coor = w->AtomCenters[2].z * scale; coor++;
   }
   return nwret;
}

/////////////////////////////////////////////////////////////////
//
//  Set the structure of a rigid water model by specifying CoM coordinates
//  of each monomer and Euler Angles
//  Conf is at least nW*6 doubles
//  scale should be set such that Conf[0]*scale is in Bohr 
//
void WaterCluster::SetConfiguration(int nW, const double *WConf, int scale, int verbose)
{

   //
   //  if the no of waters, nwaters, changes delete old arrays
   //  on first cal (and after change) allocate space for waters
   //
   if (nW != nwaters) {
      if (nwaters > 0) {
         for (int i = 0; i < nwaters; ++i)
            delete[] waters[i];
         delete[] waters;
      }
      nwaters = nW;
      waters = new Water*[nwaters];
      for (int i = 0; i < nwaters; ++i) {
	if (KTFlag == 1)
	  waters[i] = new Water(0, 0, 0, 12, 1, 1.0); // KT charges, at bit of TTM and DPP 
	else
         waters[i] = new Water(0, 0, 0, 2, 1, 1.0); // standard DPP 
      }
   }


   Water *w;
   const double *pconf = WConf;
   for (int i = 0; i < nwaters; i++, pconf+=6) {
      w = waters[i];
      w->x_origin = pconf[0] * scale;
      w->y_origin = pconf[1] * scale;
      w->z_origin = pconf[2] * scale;
      w->phi   = pconf[3];
      w->theta = pconf[4];
      w->psi   = pconf[5];
      w->UpdatePositions();
   }

   if (verbose > 7) {
      double *conf = new double[6*nwaters];
      this->GetRigidConf(conf);
      cout << "Projected Water coordinates in WaterCluster::SetStructure:\n";
      for (int i = 0; i < nwaters; i++) {
         w = waters[i];
         cout << "Water #" << i << endl;
         watergeo(w->AtomCenters[0].x*Bohr2Angs, w->AtomCenters[0].y*Bohr2Angs, w->AtomCenters[0].z*Bohr2Angs,
            w->AtomCenters[1].x*Bohr2Angs, w->AtomCenters[1].y*Bohr2Angs, w->AtomCenters[1].z*Bohr2Angs,
            w->AtomCenters[2].x*Bohr2Angs, w->AtomCenters[2].y*Bohr2Angs, w->AtomCenters[2].z*Bohr2Angs);
         printf("  Center = %f, %f, %f\n", conf[6*i+0]*Bohr2Angs, conf[6*i+1]*Bohr2Angs, conf[6*i+2]*Bohr2Angs); 
         printf("  Orient = %f, %f, %f\n", conf[6*i+3], conf[6*i+4], conf[6*i+5]); 
      }
      delete[] conf;
   }

   // induced dipoles are unknown
   DipoleSCFFlag = 0;
}




//////////////////////////////////////////////////////////////////////////
//
//  Report internal configuration (CoM position and Euler Angles
//  where the position (internal in Bohr) is scaled with scale
//  nW is the maximal no of waters to return
//  the actual no of water molecules returned is returned 
//  conf has at least the size of nW*6 doubles
//
int WaterCluster::GetConfiguration(int nW, double *WConf, double scale, int verbose)
{
   int nwret = 0;

   if (nW >= nwaters)
      nwret = nwaters;
   else {
      nwret = nW; 
      cout << "Warning: not all waters are returned by GetConfiguration\n";
   }

   Water *w;
   double *conf = WConf;
   for (int i = 0; i < nwret; i++) {
      w = waters[i];
      *conf = w->x_origin * scale; conf++;
      *conf = w->y_origin * scale; conf++;
      *conf = w->z_origin * scale; conf++;
      *conf = w->phi; conf++;
      *conf = w->theta; conf++;
      *conf = w->psi; conf++;
   }
   return nwret;
}

double WaterCluster::CalcEnergy(int verbose)
{
   if (nwaters > 0) {
      if (DipoleSCFFlag == 0)
         CalcInducedDipoles(verbose);
      return CalcIntermolecularPotential(verbose);
   }
   else
      cout << "WaterCluster::CalcEnergy: Warning: there are no waters\n";
   return 0;
}

//////////////////////////////////////////
//
//  this is the cartesian gradient wrt to all sites (4 for DPP)
//  for a rigid model it should really be the gradient for the configuration, 
//
double WaterCluster::CalcGradient(int n, double* Grad, double* PolGrad, int verbose)
{
   if (n < 3*nSitesPerMonomer*nwaters) {
      cout << "WaterCluster::CalcGradient: not enough memory for gradients of " 
         << nSitesPerMonomer*nwaters << " sites\n"; 
      exit(1);
   }

   if (nwaters > 0) {
      if (DipoleSCFFlag == 0)
         CalcInducedDipoles(verbose);
      return CalcIntermolecularGradient(Grad, PolGrad, verbose);
   }
   else
      cout << "WaterCluster::CalcEnergy: Warning: there are no waters\n";
   return 0;
}


//////////////////////////////////////////////////////////////////////////////
//
//  return lists with the objects needed for the excess electron Hamiltonian
//
//  Sites[nSites]  : 3*nSites doubles (this is the list of sites the gradients pertain to)
//                   xS1, yS1, zS1, .....  
//
//  q[nCharges]   : q1, q2, q3,
//  iqs[nCharges] : Site index of q1, site index of q2, etc
//
//  d[nDipoles]:   3*nDipoles: d1x, d1y, d1z, in au 
//  ids[nDipoles]: Site index of d1, site index of d2, etc
//
//  DmuByDR : Derivates of the dipole moments wrt the sites 
//            this is a matrix 3*nDipoles x 3*nSites 
//            the fast running index is that of the dipoles
//            the slow running index is that of the Site, i.e.:
//            Dm1x/DxS1, Dm1y/DxS1, Dm1z/DxS1, ...  Dm1x/DyS1, Dm1y/DyS1, ...
//            
//
//
//
void WaterCluster::GetLists(int nSites, double *Sites, int nCharges, double *q, int *iqs, 
                            int nDipoles, double *d, int *ids, double *DmuByDR, int verbose)
{
   if (verbose > 0)
      cout << "Lists for the electron:\n";

   class Water *W ;
   // check whether nSites, nCharges, and nDipoles match if you want to be careful

   // list of sites, charges and dipoles
   int iSite = 0;
   int iCharge = 0;

   for (int iw = 0; iw < nwaters; ++iw) {
      W = waters[iw];
      for (int is = 0; is < nSitesPerMonomer; ++is) {
         Sites[3*iSite+0] = W->AtomCenters[is].x;
         Sites[3*iSite+1] = W->AtomCenters[is].y;
         Sites[3*iSite+2] = W->AtomCenters[is].z;
         double qhere = W->AtomCenters[is].charge;
         if (qhere != 0) {
            q[iCharge] = qhere;
            iqs[iCharge] = iSite;
            iCharge ++;
         }
         iSite ++;
      }
   }

   // things that should not happen
   if (iSite != nSitesPerMonomer * nwaters) 
   {cout << "WaterCluster::GetLists: This should not happen in the site loop\n"; exit(1);} 
   if (iCharge != nChargesPerMonomer * nwaters) 
   {cout << "WaterCluster::GetLists: This should not happen in the charge loop\n"; exit(1);} 

   if (verbose > 0) {
     cout << "List of " << iSite << " Sites [Angs]:\n";
     for (int i = 0; i < iSite; ++ i)
       printf ("%-3i %10.6f  %10.6f  %10.6f\n", i,                  
	       Sites[3*i+0]*Bohr2Angs, Sites[3*i+1]*Bohr2Angs, Sites[3*i+2]*Bohr2Angs);
     cout << "List of " << iCharge << " Charges [au]:\n";
     for (int i = 0; i < iCharge; ++ i)
       printf ("%-3i %8.4f   at site %i\n", i, q[i], iqs[i]);
   }


   // for the dipole list get the induced dipoles and collapse into  
   // a single dipole for each monomer, as a center use the M site:
   // mu1x, mu1y, mu1z, mu2x, ... 
   int iDipole = 0;
   for (int iw = 0; iw < nwaters; ++iw) {
      W = waters[iw];
      double mux = 0;
      double muy = 0;
      double muz = 0;
      for (int is = 0; is < nSitesPerMonomer; ++is) {
         //  contribution from point charges
         //  if you ever need the total dipoles, here is how:
         //      double q = W->AtomCenters[is].charge;
         // if (q != 0) {
         //	mux += q * W->AtomCenters[is].x;
         //	muy += q * W->AtomCenters[is].y;
         //	muz += q * W->AtomCenters[is].z;
         //}
         // contribution from induced dipoles
         double a = W->AtomCenters[is].alpha;
         if (a != 0) {
            mux += (W->AtomCenters[is]).mx;
            muy += (W->AtomCenters[is]).my;
            muz += (W->AtomCenters[is]).mz;
         }
      }
      d[3*iDipole+0] = mux;
      d[3*iDipole+1] = muy;
      d[3*iDipole+2] = muz;
      ids[iDipole] = nSitesPerMonomer*iDipole + 3;  // this should be the M-site
      iDipole ++;
   }
   if (iDipole != nwaters) 
   {cout << "WaterCluster::GetLists: This should not happen in the dipole loop\n"; exit(1);} 
   if (verbose > 0) {
      cout << "List of " << iDipole << " Dipoles [Debye]\n";
      for (int i = 0; i < iDipole; ++ i)
         printf ("%-3i  %10.6f  %10.6f  %10.6f  at  site %i\n", i, 
         d[3*i]*AU2Debye, d[3*i+1]*AU2Debye, d[3*i+2]*AU2Debye, ids[i]);  

   }

   // put the d mu / d R matrix into DmuByDR array
   // this is needed to compute the gradient of the electronic energy

   //    Calc_dInducedDipoles(0);




}

////////////////////////////////////////////////////////////////////////////////
//
//  report all induced dipoles (self-consistent wrt water-water)
//  ("all" as opposed to collapsing the atomic dipoles to a molecular dipole)
//
void WaterCluster::ReportInducedDipoles(int nDipoles, double *d, int *ids, int verbose)
{
   class Water *W;
   int iDipole = 0;
   for (int iw = 0; iw < nwaters; ++iw) {
      W = waters[iw];
      for (int is = 0; is < nSitesPerMonomer; ++is) {
         if (W->AtomCenters[is].alpha != 0) {
            d[3*iDipole+0] = (W->AtomCenters[is]).mx;
            d[3*iDipole+1] = (W->AtomCenters[is]).my;
            d[3*iDipole+2] = (W->AtomCenters[is]).mz;
	    ids[iDipole] = nSitesPerMonomer*iw + is;
	    iDipole ++;
         }
      }
   }
   if (iDipole != 3*nwaters) 
     {cout << "WaterCluster::ReportDipoles: This should not happen in the dipole loop\n"; exit(1);} 
   if (verbose > 0) {
      cout << "List of all dipoles for the electron [Debye]:\n";      
      for (int i = 0; i < iDipole; ++ i)
         printf ("%-3i  %10.6f  %10.6f  %10.6f  at  site %i\n", 
		 i, d[3*i]*AU2Debye, d[3*i+1]*AU2Debye, d[3*i+2]*AU2Debye, ids[i]);  
   }
}


////////////////////////////////////////////////////////////////////////////////
//
//  report the polarizable sites, and the E-fields due to permanent changes at these sites
//  this is essentially the first part of CalcInducedDipoles
// 
//  for each point-polarizable site report: alpha, index in the site list, and the field 
//  this is needed for self-consistent water-water/water-electron induction
//
void WaterCluster::ReportPolSitesAndField(int nSites, double *alpha, int *ipp, double *Epc, int verbose)
{
  
  int PolSites = waters[0]->PolSites;
  if (nSites != nwaters * PolSites) {
    cout << "WaterCluster::ReportPolSitesAndField: nSites = " << nSites 
	 << " does not match nwaters * PolSites = " << nwaters * PolSites << "\n";
    exit(1);
  }
  
  // this is the relevant data structures used in CalcInducedDipoles
  //double *alpha = new double[nsites];      //  polarizabilities (isotropic)
  //double *xsite = new double[nsites];      //  positions of all sites
  //double *ysite = new double[nsites];
  //double *zsite = new double[nsites];
  //double *field1x = new double[nsites];    //  field due to permanent charges
  //double *field1y = new double[nsites];
  //double *field1z = new double[nsites];
  
  //double *Coeff = new double[3*nsites*3*nsites];
  //double *InvCoeff = new double[3*nsites*3*nsites];
  //double *dipole1 = new double[3*nsites];
  //double *dipole2 = new double[3*nsites];
  
  double dampCD = waters[0]->Damping[1];
  //double dampDD = waters[0]->Damping[2];
  
  dVec R; R.resize(3*nSites);
  
  // make lists of polarizable sites, initialize the E-field
  int SiteIndex = 0;
  int isite = 0;
  for (int j = 0; j < nwaters; ++j) {
    for (int k = 0; k < waters[j]->nAtomCenters; ++k) {
      double a = (waters[j]->AtomCenters[k]).alpha;
      if (a != 0) {
	if (verbose > 15)
	  cout << "Init polarizable center#" << isite << " at water#" << j << " AtomCenter#" << k << endl;
	alpha[isite] = a;
	R[3*isite+0] = (waters[j]->AtomCenters[k]).x; 
	R[3*isite+1] = (waters[j]->AtomCenters[k]).y; 
	R[3*isite+2] = (waters[j]->AtomCenters[k]).z;
	Epc[3*isite+0] = 0;
	Epc[3*isite+1] = 0;
	Epc[3*isite+2] = 0;
	ipp[isite] = SiteIndex;
	++isite;
      }
      SiteIndex ++;
    }
  }
  if (isite != nSites) { 
    cout << "WaterCluster::ReportPolSitesAndField: nSites = << " << nSites 
	 << " and isite = " << isite << ";  this should not happen.\n"; 
    exit(1);
  }
  
  // compute the fields due to the permanent charges:   Epc = Sum_i Q_i*R_i / |R_i|^3
  if (verbose > 18) cout << "Computing fields due to permanent charges\n";
  for (int j = 0; j < nwaters; ++j) {
    if (verbose > 18) cout << "Water #" << j << endl;
    for (int k = 0; k < waters[j]->nAtomCenters; ++k) {
      AtomCenter &PointCharge = waters[j]->AtomCenters[k];
      double q = PointCharge.charge;
      double alphak = PointCharge.sigma;
//      cout<< " water["<<j<<"] Atomcenter["<<k<<"] -------\n";
      if (verbose > 18) 
	cout << "  AtomCenter#" << k << " with q=" << q 
	     << " R = " << PointCharge.x << " " << PointCharge.y << " " << PointCharge.z << endl;
      if (q != 0) {
	for (int i = 0; i < nSites; ++i) {
	  if (i / PolSites == j)  // no interactions between sites on same water
	    continue;
	  double Rx = R[3*i+0] - PointCharge.x;
	  double Ry = R[3*i+1] - PointCharge.y;
	  double Rz = R[3*i+2] - PointCharge.z;
	  double tensor = q * pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5);         //  q = charge/R^3
	  
	  if (verbose > 18) 
	    cout << "    acting on site#" << i << " at " << R[3*i+0] << " " << R[3*i+1] << " " << R[3*i+2] 
		 << " tensor = " << tensor << endl
		 << "      Rik = " << Rx << " " << Ry << " " << Rz<<" alpha[i]"<<alpha[i] << endl;
	  
	  
	  double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
	  double damp = ElectrostaticDCDamp(distance, alpha[i], alphak, dampCD);
	  
	  Epc[3*i+0] += tensor * Rx * damp;
	  Epc[3*i+1] += tensor * Ry * damp;
	  Epc[3*i+2] += tensor * Rz * damp;
	}
      }
    }
  }
  
//  cout << "Field due to point charges\n";
//  for (int i = 0; i < nSites; ++i) {
//      cout << i << "  " << Epc[3*i+0] << endl;
//      cout << i << "  " << Epc[3*i+1] << endl;
//      cout << i << "  " << Epc[3*i+2] << endl;
//   }

}



void WaterCluster::GetDInducedDipoles(double *DmuByDR, double *dEfield, double *dTensor)
{

   Calc_dInducedDipoles(0, DmuByDR, dEfield, dTensor);

}


//////////////////////////////////////////
//
//  GetRigidConf
//
//  gets the vector with x, y, z, phi, theta, psi for all water monomers
//
void WaterCluster::GetRigidConf(double *configuration)
{
   for (int k=0; k < nwaters; k++) {
      configuration[k*6]  =waters[k]->x_origin;
      configuration[k*6+1]=waters[k]->y_origin;
      configuration[k*6+2]=waters[k]->z_origin;
      configuration[k*6+3]=waters[k]->phi;
      configuration[k*6+4]=waters[k]->theta;
      configuration[k*6+5]=waters[k]->psi;
   }
}






/////////////////////////////////////////////////////////////////////
//
//  ComputeInducedDipolesAndIntermolecularEnergy
//
//  this function deals with the neutral cluster
//   - compute the intermolecular induced dipoles self consistently
//   - compute the total intermolecular energy (electrostatic, polarization, and van der Waals energies)
//   - condense the three dipoles of the TTM model into a single for the electron
//     
//
//
//  Seymour's old code for the Dang Chang water model (still in Cluster.cpp)
//
//  this->CalculateInterMoleculeInteractions(2);
//  CalculateNuclearRepulsionEnergy();
//  CalculateLJEnergy();  
//
//
double WaterCluster::ComputeInducedDipolesAndIntermolecularEnergy(int verbose) 
{
   CalcInducedDipoles(verbose);
   double ime = CalcIntermolecularPotential(verbose);
   AddUpTTMDipoles(verbose);
   return ime;
}  

/////////////////////////////////////////////
//
//  AddUpTTMDipoles
// 
//  the electron sees only a single polarizable site on each water monomer
//  centered at the M-site
//  the associated dipole is set to the vector sum of the three induced TTM dipoles
//
void WaterCluster::AddUpTTMDipoles(int verbose) 
{
   double mux, muy, muz;
   int j, k, tester;
   if (verbose > 1) 
      cout << "Adding up the induced TTM dipoles\n";
   for (j = 0; j < nwaters; ++j) {
      if (verbose > 1) 
         cout << " Water#" << j;
      mux = muy = muz = 0;
      tester = 0;
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            tester ++;
            mux += (waters[j]->AtomCenters[k]).mx;
            muy += (waters[j]->AtomCenters[k]).my;
            muz += (waters[j]->AtomCenters[k]).mz;
         }
      }
      if (tester != 3) {cout << "This is not good\n"; exit(1);}
      if (verbose > 1) {
         double dpm = sqrt(mux*mux + muy*muy + muz*muz);
         cout << " " << dpm << "  (" << mux << ", " << muy << ", " << muz << ")\n";
      }
      //oscs[j]->SetDipole(mux, muy, muz);
   }
}





/////////////////////////////////////////////////////////////////////////
//
//  WaterCluster::CalcInducedDipoles
//
//  step #1 of the water model: compute the induced dipoles on the water molecules self consistently
//
//  in the TTM each water has 3 polarizable sites
//  the polarizability of each site is initialized in Water::SetupCharge,
//  and this is alos were the order of sites: O, H1, H2, M is defined
//  
//  this functions starts by initializing the positions of the polarizable sites by reading
//  the coordinates of H1, H1, and M
//  then the E-field due to the permanent charges is computed --> field1 
//  and the E-field due to the induced dipoles (field2) is iterated until self consistent  
//
//
void WaterCluster::CalcInducedDipoles(int verbose)
{
   int i, j, k;

   // establish lists with all polarizable sites and the fields at these positions
   // after convergence these data go back into waters[j], and can be used in the energy expressions
   int nsites = nwaters * waters[0]->PolSites;
   int PolSites = waters[0]->PolSites;

   double *alpha = new double[nsites];      //  polarizabilities (isotropic)
   double *xsite = new double[nsites];      //  positions of all sites
   double *ysite = new double[nsites];
   double *zsite = new double[nsites];
   double *dipolex = new double[nsites];    //  current dipoles of all sites
   double *dipoley = new double[nsites];
   double *dipolez = new double[nsites];
   double *olddplx = new double[nsites];    //  old dipoles for convergence checking
   double *olddply = new double[nsites];
   double *olddplz = new double[nsites];
   double *field1x = new double[nsites];    //  field due to permanent charges
   double *field1y = new double[nsites];
   double *field1z = new double[nsites];
   double *field2x = new double[nsites];    //  field due to induced dipoles
   double *field2y = new double[nsites];
   double *field2z = new double[nsites];


   double *Coeff = new double[3*nsites*3*nsites];
   double *InvCoeff = new double[3*nsites*3*nsites];
   double *dipole1 = new double[3*nsites];
   double *dipole2 = new double[3*nsites];

   bool Inversion_check=false;

   //double newfac = 0.925;			//from xantheas
   //double oldfac = 0.075;			//from xantheas			these two factors are needed to add dipoles
   double newfac = 0.5;
   double oldfac = 1 - newfac;
   //                                    newdip = oldfac*olddip + newfac*alpha*field

   double dampCD = waters[0]->Damping[1];
   double dampDD = waters[0]->Damping[2];

  // cout<<" dampCD and dampDD" << dampCD <<" "<<dampDD<<endl;
   
   // initialize the polarizable sites
   i = 0;
   for (j = 0; j < nwaters; ++j) {
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            if (verbose > 15)
               cout << "Init polarizable center#" << i << " at water#" << j << " AtomCenter#" << k << endl;
            alpha[i] = a;
            xsite[i] = (waters[j]->AtomCenters[k]).x; 
            ysite[i] = (waters[j]->AtomCenters[k]).y; 
            zsite[i] = (waters[j]->AtomCenters[k]).z;
            field1x[i] = 0; field1y[i] = 0; field1z[i] = 0;
            field2x[i] = 0; field2y[i] = 0; field2z[i] = 0;
            dipolex[i] = 0; dipoley[i] = 0; dipolez[i] = 0;
            ++i;
         }
      }
   }
   if (i != nsites) { cout << "WaterCluster::CalcInducedDipoles:  This should not happen.\n"; exit(1);}



   // compute the fields due to the permanent charges --> field1   E = Sum_i Q_i*R_i / |R_i|^3
   if (verbose > 18) cout << "Computing fields due to permanent charges\n";
   for (j = 0; j < nwaters; ++j) {
      if (verbose > 18) cout << "Water #" << j << endl;
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         AtomCenter &PointCharge = waters[j]->AtomCenters[k];
         double q = PointCharge.charge;
         double alphak = PointCharge.sigma;
         if (verbose > 18) 
            cout << "  AtomCenter#" << k << " with q=" << q 
            << " R = " << PointCharge.x << " " << PointCharge.y << " " << PointCharge.z << endl;
         if (q != 0) {
            for (i = 0; i < nsites; ++i) {
               if (i / PolSites == j)  // no interactions between sites on same water
                  continue;
               double Rx = xsite[i] - PointCharge.x;
               double Ry = ysite[i] - PointCharge.y;
               double Rz = zsite[i] - PointCharge.z;
               double tensor = q * pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5);         //  q = charge/R^3

               if (verbose > 18) 
                  cout << "    acting on site#" << i << " at " << xsite[i] << " " << ysite[i] << " " << zsite[i] 
               << " tensor = " << tensor << endl
                  << "      Rik = " << Rx << " " << Ry << " " << Rz<<" alpha[i]"<<alpha[i] << endl;


               double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
               double damp = ElectrostaticDCDamp(distance, alpha[i], alphak, dampCD);
               //	  cout<<" alpha[i] and alphak = "<<alpha[i]<<" "<<alphak<<"\n";
               //double damp = 1.0;

               field1x[i] += tensor * Rx * damp;
               field1y[i] += tensor * Ry * damp;
               field1z[i] += tensor * Rz * damp;
            }
         }
      }
   }

   // compute induced dipoles due to the permanent charges to start the self consitent iterations
   //cout << "Fields:\n";
   for (i = 0; i < nsites; ++i) {
      dipolex[i] = newfac * alpha[i] * field1x[i];
      dipoley[i] = newfac * alpha[i] * field1y[i];
      dipolez[i] = newfac * alpha[i] * field1z[i];
      dipole1[i*3] = -alpha[i] * field1x[i];
      dipole1[i*3+1] = -alpha[i] * field1y[i];
      dipole1[i*3+2] = -alpha[i] * field1z[i];
      //    cout<<" newfac = "<<newfac<<" dipolex["<<i<<"]= "<< dipolex[i]<<"\n";
      //    cout<<" newfac = "<<newfac<<" dipoley["<<i<<"]= "<< dipoley[i]<<"\n";
      //    cout<<" newfac = "<<newfac<<" dipolez["<<i<<"]= "<< dipolez[i]<<"\n";
      //cout << i << "  " << field1x[i] << endl;
     // cout << i << "  " << field1y[i] << endl;
     // cout << i << "  " << field1z[i] << endl;

   }




   if (verbose > 18) {
      printf("Permanent charge induced Dipoles [Debye]:\n");
      for (i = 0; i < nsites; ++i) {
         double absdp = sqrt(dipolex[i]*dipolex[i] + dipoley[i]*dipoley[i] + dipolez[i]*dipolez[i]);
         printf("   PolSite #%d  mu = %10.4f   at (%10.4f %10.4f %10.4f)  mu = (%10.8f %10.8f %10.8f)\n",
            i, absdp, xsite[i]*Bohr2Angs, ysite[i]*Bohr2Angs, zsite[i]*Bohr2Angs, 
            dipolex[i]*AU2Debye, dipoley[i]*AU2Debye, dipolez[i]*AU2Debye);
      }
   }


   if (Inversion_check){
      for (i = 0;  i < nsites; ++i) {             // at site i
         for (j = 0;  j < nsites; ++j) {           // dipole j creates field E
            if (j == i){ 
               Coeff[i*3*3*nsites+j*3]=-1.;
               Coeff[i*3*3*nsites+j*3+1]=0.;
               Coeff[i*3*3*nsites+j*3+2]=0.;
               Coeff[(i*3+1)*3*nsites+j*3]=0.;
               Coeff[(i*3+1)*3*nsites+j*3+1]=-1.;
               Coeff[(i*3+1)*3*nsites+j*3+2]=0.;
               Coeff[(i*3+2)*3*nsites+j*3]=0.;
               Coeff[(i*3+2)*3*nsites+j*3+1]=0.;
               Coeff[(i*3+2)*3*nsites+j*3+2]=-1.;
               continue;
            }

            double Txx,Tyy,Tzz;
            double Txy,Txz,Tyz;

            dFieldOfDipole(dipolex[j], dipoley[j], dipolez[j], xsite[j], ysite[j], zsite[j],
               xsite[i], ysite[i], zsite[i], Txx, Tyy, Tzz, Txy, Txz, Tyz, alpha[i],alpha[j], dampDD);

            double Tyx = Txy;
            double Tzx = Txz;
            double Tzy = Tyz;
            Coeff[i*3*3*nsites+j*3]=alpha[i]*Txx;
            Coeff[i*3*3*nsites+j*3+1]=alpha[i]*Txy;
            Coeff[i*3*3*nsites+j*3+2]=alpha[i]*Txz;
            Coeff[(i*3+1)*3*nsites+j*3]=alpha[i]*Tyx;
            Coeff[(i*3+1)*3*nsites+j*3+1]=alpha[i]*Tyy;
            Coeff[(i*3+1)*3*nsites+j*3+2]=alpha[i]*Tyz;
            Coeff[(i*3+2)*3*nsites+j*3]=alpha[i]*Tzx;
            Coeff[(i*3+2)*3*nsites+j*3+1]=alpha[i]*Tzy;
            Coeff[(i*3+2)*3*nsites+j*3+2]=alpha[i]*Tzz;

         }
      }

      InverseMatrix(Coeff, InvCoeff, nsites*3, 1 );


      for (i = 0;  i < 3*nsites; ++i) {
         dipole2[i]=0.0;
         for (j = 0;  j < 3*nsites; ++j) {
            dipole2[i] += InvCoeff[i*3*nsites+j]*dipole1[j];
         }
      } 
   }

   if(!Inversion_check) {

      //
      // self consistent induced dipoles iteration 
      //
      int converged = 0;
      int max_iterations = 5*(nsites+10);  // I just made this up
      double thresh = 1e-12;
      double tester = 0;
      for (int iteration = 0; iteration < max_iterations; ++ iteration) 
      {
         // compute the field due to the induced dipoles; result is field2
         for (i = 0;  i < nsites; ++i)
            field2x[i] = field2y[i] = field2z[i] = 0;
         for (i = 0;  i < nsites; ++i) {             // at site i
            for (j = 0;  j < nsites; ++j) {           // dipole j creates field E
               if (j == i) 
                  continue;

               double Ex, Ey, Ez;
               FieldOfDipole(dipolex[j], dipoley[j], dipolez[j], xsite[j], ysite[j], zsite[j],
                  xsite[i], ysite[i], zsite[i], Ex, Ey, Ez,alpha[i],alpha[j], dampDD);
               field2x[i] += Ex; 
               field2y[i] += Ey; 
               field2z[i] += Ez; 
            }
         }
         // compute the new dipoles
         for (i = 0;  i < nsites; ++i) {
            olddplx[i] = dipolex[i];   dipolex[i] = oldfac*olddplx[i] + newfac * alpha[i] * (field1x[i] + field2x[i]);
            olddply[i] = dipoley[i];   dipoley[i] = oldfac*olddply[i] + newfac * alpha[i] * (field1y[i] + field2y[i]);
            olddplz[i] = dipolez[i];   dipolez[i] = oldfac*olddplz[i] + newfac * alpha[i] * (field1z[i] + field2z[i]);
               //  cout<<" new  dipole["<<i<<"]= "<< dipolex[i]<<" "<<dipoley[i]<<" "<<dipolez[i]<<"\n";
               //  cout<<" old  dipole["<<i<<"]= "<< olddplx[i]<<" new : "<<dipolex[i]<<"\n";
         }
             // cout<<" new  dipole[1]= "<< dipolex[1]<<" "<<dipoley[1]<<" "<<dipolez[1]<<"\n";
         // check convergence
         tester = 0;
         for (i = 0;  i < nsites; ++i) {
            double tx = olddplx[i] - dipolex[i];
            double ty = olddply[i] - dipoley[i];
            double tz = olddplz[i] - dipolez[i];
            tester += tx*tx + ty*ty + tz*tz;
         }
         tester = sqrt(tester)/3/nsites;  // this is the average change
         if (verbose > 15) {
            printf("Iteration #%d  tester = %9.2e:\n", iteration, tester);
            if (verbose > 18) 
               for (i = 0; i < nsites; ++i) {
                  double absdp = sqrt(dipolex[i]*dipolex[i] + dipoley[i]*dipoley[i] + dipolez[i]*dipolez[i]);
                  printf("   Site #%d  mu = %10.4f   mu = (%10.8f %10.8f %10.8f)\n",
                     i, absdp*AU2Debye, dipolex[i]*AU2Debye, dipoley[i]*AU2Debye, dipolez[i]*AU2Debye);
               }
         } 
         if (tester < thresh) {
            converged = 1;
            break;
         }
      }

      if (converged == 0) {
         cout << "**** Warning in WaterCluster::CalcInducedDipoles():";
         cout << " after " << max_iterations << " iterations the induced dipoles converged only to " 
            << tester << endl;
      }
      else
         if (verbose > 1) cout << "Induced dipoles converged to " << thresh << endl;


      // put the results into the waters
      i = 0;
      for (j = 0; j < nwaters; ++j)
         for (k = 0; k < waters[j]->nAtomCenters; ++k) {
            double a = (waters[j]->AtomCenters[k]).alpha;
            if (a != 0) {
               (waters[j]->AtomCenters[k]).mx = dipolex[i];
               (waters[j]->AtomCenters[k]).my = dipoley[i];
               (waters[j]->AtomCenters[k]).mz = dipolez[i];

               //      cout<<" AtomCenters["<<k<<"].mx= "<<waters[j]->AtomCenters[k].mx<<"  dipolex["<<i<<"]="<<dipolex[i]<<endl;
               //      cout<<"  dipole2["<<i<<"]="<<dipole2[i]<<endl;
               ++i;
            }
         }

   }

   if (Inversion_check){
      i=0;
      for (j = 0; j < nwaters; ++j)
         for (k = 0; k < waters[j]->nAtomCenters; ++k) {
            double a = (waters[j]->AtomCenters[k]).alpha;
            if (a != 0) {
               (waters[j]->AtomCenters[k]).mx = dipole2[i*3];
               (waters[j]->AtomCenters[k]).my = dipole2[i*3+1];
               (waters[j]->AtomCenters[k]).mz = dipole2[i*3+2];
               //      cout<<" AtomCenters["<<k<<"].mx= "<<waters[j]->AtomCenters[k].mx<<endl;
               ++i;
            }
         }
   }

   if (verbose > 1) 
      ReportInducedDipoles();

   delete[] field2z; delete[] field2y; delete[] field2x;
   delete[] field1z; delete[] field1y; delete[] field1x;
   delete[] dipolez; delete[] dipoley; delete[] dipolex;
   delete[] olddplz; delete[] olddply; delete[] olddplx;
   delete[] zsite;   delete[] ysite;   delete[] xsite; delete[] alpha;


   delete[] Coeff; delete[] InvCoeff;
   delete[] dipole1; delete[] dipole2; 

   DipoleSCFFlag = 1;

}


//
//
//  Find the derivative of Induced Dipole moment analytically
//                                          by Tae Hoon Choi
//
///////////////////////////////////////////////////////////////
//void WaterCluster::Calc_dInducedDipoles(int verbose, double *DmuByDx, double *DmuByDy,double *DmuByDz)
void WaterCluster::Calc_dInducedDipoles(int verbose, double *DmuByDR, double *dEfield, double *dTensor)
{
   int i, j, k;

   // establish lists with all polarizable sites and the fields at these positions
   // after convergence these data go back into waters[j], and can be used in the energy expressions
   int nsites = nwaters * waters[0]->PolSites;
   int PolSites = waters[0]->PolSites;

   double *alpha = new double[nsites];      //  polarizabilities (isotropic)
   double *xsite = new double[nsites];      //  positions of all sites
   double *ysite = new double[nsites];
   double *zsite = new double[nsites];
   double *dipolex = new double[nsites];    //  current dipoles of all sites
   double *dipoley = new double[nsites];
   double *dipolez = new double[nsites];

   double *dxfield2x = new double[nsites];    //  field due to induced dipoles
   double *dxfield2y = new double[nsites];
   double *dxfield2z = new double[nsites];
   double *dyfield2x = new double[nsites];    
   double *dyfield2y = new double[nsites];
   double *dyfield2z = new double[nsites];
   double *dzfield2x = new double[nsites];    
   double *dzfield2y = new double[nsites];
   double *dzfield2z = new double[nsites];

   double *dxfield1x = new double[nsites];    //  field due to permanent charges
   double *dxfield1y = new double[nsites];
   double *dxfield1z = new double[nsites];
   double *dyfield1x = new double[nsites];    
   double *dyfield1y = new double[nsites];
   double *dyfield1z = new double[nsites];
   double *dzfield1x = new double[nsites];    
   double *dzfield1y = new double[nsites];
   double *dzfield1z = new double[nsites];

   double *Coeff = new double[3*nsites*3*nsites];
   double *dxCoeff = new double[3*nsites];
   double *dyCoeff = new double[3*nsites];
   double *dzCoeff = new double[3*nsites];



   double dampCD = waters[0]->Damping[1];
   double dampDD = waters[0]->Damping[2];


   for(int i=0; i<nwaters; i++){
      for (int j=0;j<waters[i]->nAtomCenters;j++){
         if (waters[i]->AtomCenters[j].dxmx == NULL)
            waters[i]->AtomCenters[j].dxmx = new double[nwaters];
         if (waters[i]->AtomCenters[j].dxmy == NULL) 
            waters[i]->AtomCenters[j].dxmy = new double[nwaters];
         if (waters[i]->AtomCenters[j].dxmz == NULL) 
            waters[i]->AtomCenters[j].dxmz = new double[nwaters];
         if (waters[i]->AtomCenters[j].dymx == NULL) 
            waters[i]->AtomCenters[j].dymx = new double[nwaters];
         if (waters[i]->AtomCenters[j].dymy == NULL) 
            waters[i]->AtomCenters[j].dymy = new double[nwaters];
         if (waters[i]->AtomCenters[j].dymz == NULL) 
            waters[i]->AtomCenters[j].dymz = new double[nwaters];
         if (waters[i]->AtomCenters[j].dzmx == NULL) 
            waters[i]->AtomCenters[j].dzmx = new double[nwaters];
         if (waters[i]->AtomCenters[j].dzmy == NULL) 
            waters[i]->AtomCenters[j].dzmy = new double[nwaters];
         if (waters[i]->AtomCenters[j].dzmz == NULL) 
            waters[i]->AtomCenters[j].dzmz = new double[nwaters];
      }
   }



   // initialize the polarizable sites
   i = 0;
   for (j = 0; j < nwaters; ++j) {
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            if (verbose > 15)
               cout << "Init polarizable center#" << i << " at water#" << j << " AtomCenter#" << k << endl;
            alpha[i] = a;
            xsite[i] = (waters[j]->AtomCenters[k]).x; 
            ysite[i] = (waters[j]->AtomCenters[k]).y; 
            zsite[i] = (waters[j]->AtomCenters[k]).z;
            dipolex[i] = 0; dipoley[i] = 0; dipolez[i] = 0;
            ++i;
         }
      }
   }
   if (i != nsites) { cout << "WaterCluster::CalcInducedDipoles:  This should not happen.\n"; exit(1);}



   // compute the fields due to the permanent charges --> field1   E = Sum_i Q_i*R_i / |R_i|^3

   i = 0;
   for (j = 0; j < nwaters; ++j)
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            dipolex[i]=(waters[j]->AtomCenters[k]).mx ;
            dipoley[i]=(waters[j]->AtomCenters[k]).my ;
            dipolez[i]=(waters[j]->AtomCenters[k]).mz ;
            ++i;
         }
      }


      //
      // self consistent induced dipoles iteration 
      //
      // compute the field due to the induced dipoles; result is field2
      for (i = 0;  i < nsites; ++i) {             // at site i
         for (j = 0;  j < nsites; ++j) {           // dipole j creates field E
            if (j == i){ 
               Coeff[i*3*3*nsites+j*3]=-1.;
               Coeff[i*3*3*nsites+j*3+1]=0.;
               Coeff[i*3*3*nsites+j*3+2]=0.;
               Coeff[(i*3+1)*3*nsites+j*3]=0.;
               Coeff[(i*3+1)*3*nsites+j*3+1]=-1.;
               Coeff[(i*3+1)*3*nsites+j*3+2]=0.;
               Coeff[(i*3+2)*3*nsites+j*3]=0.;
               Coeff[(i*3+2)*3*nsites+j*3+1]=0.;
               Coeff[(i*3+2)*3*nsites+j*3+2]=-1.;
            }
            else{
               double Txx,Tyy,Tzz;
               double Txy,Txz,Tyz;

               dFieldOfDipole(dipolex[j], dipoley[j], dipolez[j], xsite[j], ysite[j], zsite[j],
                  xsite[i], ysite[i], zsite[i], Txx, Tyy, Tzz, Txy, Txz, Tyz, alpha[i],alpha[j], dampDD);

               double Tyx = Txy;
               double Tzx = Txz;
               double Tzy = Tyz;

               //cout<<"alpha[i]="<<alpha[i]<<"\n";
               Coeff[i*3*3*nsites+j*3]=alpha[i]*Txx;
               Coeff[i*3*3*nsites+j*3+1]=alpha[i]*Txy;
               Coeff[i*3*3*nsites+j*3+2]=alpha[i]*Txz;
               Coeff[(i*3+1)*3*nsites+j*3]=alpha[i]*Tyx;
               Coeff[(i*3+1)*3*nsites+j*3+1]=alpha[i]*Tyy;
               Coeff[(i*3+1)*3*nsites+j*3+2]=alpha[i]*Tyz;
               Coeff[(i*3+2)*3*nsites+j*3]=alpha[i]*Tzx;
               Coeff[(i*3+2)*3*nsites+j*3+1]=alpha[i]*Tzy;
               Coeff[(i*3+2)*3*nsites+j*3+2]=alpha[i]*Tzz;
            }
         }
      } 


      double *InvCoeff = new double[3*nsites*3*nsites];
      InverseMatrix(Coeff, InvCoeff, nsites*3, 1 );


      double *dxdipole = new double[3*nsites];
      double *dydipole = new double[3*nsites];
      double *dzdipole = new double[3*nsites];

      int BaseR;
      int BaseR2;
      int TotalSites=nsites*9;
      int Tindex=0;


      int index=0;
      for (int jj = 0; jj < nwaters; ++jj) {
         for (int kk = 0; kk < waters[jj]->nAtomCenters; ++kk) {
            if (kk ==3) BaseR=-1; 
            else BaseR = jj*3+kk;
            BaseR2 = jj*waters[jj]->nAtomCenters+kk;

            for (int ii = 0; ii < nsites; ++ii) {
               dxfield1x[ii] = 0.; dxfield1y[ii] = 0.; dxfield1z[ii] = 0.;
               dyfield1x[ii] = 0.; dyfield1y[ii] = 0.; dyfield1z[ii] = 0.;
               dzfield1x[ii] = 0.; dzfield1y[ii] = 0.; dzfield1z[ii] = 0.;
               dxfield2x[ii] = 0.; dxfield2y[ii] = 0.; dxfield2z[ii] = 0.;
               dyfield2x[ii] = 0.; dyfield2y[ii] = 0.; dyfield2z[ii] = 0.;
               dzfield2x[ii] = 0.; dzfield2y[ii] = 0.; dzfield2z[ii] = 0.;
            }

            for (j = 0; j < nwaters; ++j) {
               for (k = 0; k < waters[j]->nAtomCenters; ++k) {
                  AtomCenter &PointCharge = waters[j]->AtomCenters[k];
                  double q = PointCharge.charge;
                  double alphak = PointCharge.sigma;
                  if (q != 0) {
                     for (i = 0; i < nsites; ++i) {
                        if (i / PolSites == j)  // no interactions between sites on same water
                           continue;
                        if (BaseR == i || BaseR2 == j*waters[j]->nAtomCenters+k){
                           double Rx = xsite[i] - PointCharge.x;
                           double Ry = ysite[i] - PointCharge.y;
                           double Rz = zsite[i] - PointCharge.z;
                           double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
                           double damp = ElectrostaticDCDamp(distance, alpha[i], alphak, dampCD);
                           double tensor = q *damp* pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5);         //  q = charge/R^3
                           double de6 = q * 3.0*dampCD*pow(alpha[i]*alphak,-0.5)*(1.0-damp)/(Rx*Rx + Ry*Ry + Rz*Rz);
                           double de7 = -q * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0;	  

                           if (BaseR == i){
                              dxfield1x[i] += tensor+(de6+de7)*Rx*Rx ;
                              dxfield1y[i] +=        (de6+de7)*Ry*Rx ;
                              dxfield1z[i] +=        (de6+de7)*Rz*Rx ;
                              dyfield1x[i] +=        (de6+de7)*Rx*Ry ;
                              dyfield1y[i] += tensor+(de6+de7)*Ry*Ry ;
                              dyfield1z[i] +=        (de6+de7)*Rz*Ry ;
                              dzfield1x[i] +=        (de6+de7)*Rx*Rz ;
                              dzfield1y[i] +=        (de6+de7)*Ry*Rz ;
                              dzfield1z[i] += tensor+(de6+de7)*Rz*Rz ;
                           }
                           else if (BaseR2 == j*waters[j]->nAtomCenters+k){
                              dxfield1x[i] -= tensor+(de6+de7)*Rx*Rx ;
                              dxfield1y[i] -= (de6+de7)* Ry*Rx ;
                              dxfield1z[i] -= (de6+de7)* Rz*Rx ;
                              dyfield1x[i] -= (de6+de7)*Rx*Ry ;
                              dyfield1y[i] -= tensor+(de6+de7)* Ry*Ry ;
                              dyfield1z[i] -= (de6+de7)* Rz*Ry ;
                              dzfield1x[i] -= (de6+de7)*Rx*Rz ;
                              dzfield1y[i] -= (de6+de7)*Ry*Rz ;
                              dzfield1z[i] -= tensor+(de6+de7)* Rz*Rz ;
                           }
                        }
                     }
                  }
               }
            }


            for (i = 0;  i < nsites; ++i) {             // at site i
               for (j = 0;  j < nsites; ++j) {           // dipole j creates field E
                  if (j == i) 
                     continue;

                  if(BaseR == i || BaseR==j){
                     double dxEx, dxEy, dxEz;
                     double dyEx, dyEy, dyEz;
                     double dzEx, dzEy, dzEz;

                     dFieldOfDipole(dipolex[j], dipoley[j], dipolez[j], xsite[j], ysite[j], zsite[j],
                        xsite[i], ysite[i], zsite[i], dxEx, dxEy, dxEz,dyEx,dyEy,dyEz,
                        dzEx,dzEy,dzEz, alpha[i],alpha[j], dampDD);
                     if (BaseR == i){
                        dxfield2x[i] += dxEx; 
                        dxfield2y[i] += dxEy; 
                        dxfield2z[i] += dxEz; 
                        dyfield2x[i] += dyEx; 
                        dyfield2y[i] += dyEy; 
                        dyfield2z[i] += dyEz; 
                        dzfield2x[i] += dzEx; 
                        dzfield2y[i] += dzEy; 
                        dzfield2z[i] += dzEz; 
                     }
                     else if (BaseR == j){
                        dxfield2x[i] -= dxEx; 
                        dxfield2y[i] -= dxEy; 
                        dxfield2z[i] -= dxEz; 
                        dyfield2x[i] -= dyEx; 
                        dyfield2y[i] -= dyEy; 
                        dyfield2z[i] -= dyEz; 
                        dzfield2x[i] -= dzEx; 
                        dzfield2y[i] -= dzEy; 
                        dzfield2z[i] -= dzEz; 
                     }
                  }
               }
            }


            // compute the new dipoles
            //      for (i = 0;  i < 3*nsites*3*nsites; ++i) {
            //	      cout<<"Coeff["<<i<<"]= "<<Coeff[i]<<"\n";
            //      }
            for (i = 0;  i < nsites; ++i) {
               dxCoeff[i*3]   = -alpha[i] *(dxfield1x[i] + dxfield2x[i]);
               dxCoeff[i*3+1] = -alpha[i] *(dxfield1y[i] + dxfield2y[i]);
               dxCoeff[i*3+2] = -alpha[i] *(dxfield1z[i] + dxfield2z[i]);
               dyCoeff[i*3]   = -alpha[i] *(dyfield1x[i] + dyfield2x[i]);
               dyCoeff[i*3+1] = -alpha[i] *(dyfield1y[i] + dyfield2y[i]);
               dyCoeff[i*3+2] = -alpha[i] *(dyfield1z[i] + dyfield2z[i]);
               dzCoeff[i*3]   = -alpha[i] *(dzfield1x[i] + dzfield2x[i]);
               dzCoeff[i*3+1] = -alpha[i] *(dzfield1y[i] + dzfield2y[i]);
               //	  dzCoeff[i*3+1] = 0;
               dzCoeff[i*3+2] = -alpha[i] *(dzfield1z[i] + dzfield2z[i]);


               dEfield[TotalSites*Tindex + i*3]   = dxfield1x[i] ;
               dEfield[TotalSites*Tindex + i*3+1] = dxfield1y[i] ;
               dEfield[TotalSites*Tindex + i*3+2] = dxfield1z[i] ;
               dEfield[TotalSites*Tindex + 3*nsites + i*3]   =dyfield1x[i] ;
               dEfield[TotalSites*Tindex + 3*nsites + i*3+1] =dyfield1y[i] ;
               dEfield[TotalSites*Tindex + 3*nsites + i*3+2] =dyfield1z[i] ;
               dEfield[TotalSites*Tindex + 6*nsites + i*3]   =dzfield1x[i] ;
               dEfield[TotalSites*Tindex + 6*nsites + i*3+1] =dzfield1y[i] ;
               dEfield[TotalSites*Tindex + 6*nsites + i*3+2] =dzfield1z[i] ;

               dTensor[TotalSites*Tindex + i*3]   =  dxfield2x[i];
               dTensor[TotalSites*Tindex + i*3+1] =  dxfield2y[i];
               dTensor[TotalSites*Tindex + i*3+2] =  dxfield2z[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3]   = dyfield2x[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3+1] = dyfield2y[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3+2] = dyfield2z[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3]   = dzfield2x[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3+1] = dzfield2y[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3+2] = dzfield2z[i];


//     cout<<" dEfield["<<TotalSites*Tindex + i*3<<"] :" << dEfield[TotalSites*Tindex + i*3]<<"\n";
//    cout<<" dEfield[6*nsites + 3*i+0] :" << dEfield[6*nsites + 3*i+0]<<"\n";
//    exit(0);

            
            }

            Tindex = Tindex+1;

            for (i = 0;  i < 3*nsites; ++i) {
               dxdipole[i]=0.0;
               dydipole[i]=0.0;
               dzdipole[i]=0.0;
               for (j = 0;  j < 3*nsites; ++j) {
                  dxdipole[i] += InvCoeff[i*3*nsites+j]*dxCoeff[j];
                  dydipole[i] += InvCoeff[i*3*nsites+j]*dyCoeff[j];
                  dzdipole[i] += InvCoeff[i*3*nsites+j]*dzCoeff[j];
               }
            } 


            // put the results into the waters


            // sum results of dmu/dx for all polarizable sites of each water into one
            // common polarizable site at the M site (for the electron all induced dipoles
            // are summed into one dipole at the M site)

            for (i = 0;  i < nwaters; ++i) {

               waters[jj]->AtomCenters[kk].dxmx[i]= dxdipole[i*3*3]+dxdipole[i*3*3+3]+dxdipole[i*3*3+6];
               waters[jj]->AtomCenters[kk].dxmy[i]= dxdipole[i*3*3+1]+dxdipole[i*3*3+1+3]+dxdipole[i*3*3+1+6];
               waters[jj]->AtomCenters[kk].dxmz[i]= dxdipole[i*3*3+2]+dxdipole[i*3*3+2+3]+dxdipole[i*3*3+2+6];
               waters[jj]->AtomCenters[kk].dymx[i]= dydipole[i*3*3]+dydipole[i*3*3+3]+dydipole[i*3*3+6];
               waters[jj]->AtomCenters[kk].dymy[i]= dydipole[i*3*3+1]+dydipole[i*3*3+1+3]+dydipole[i*3*3+1+6];
               waters[jj]->AtomCenters[kk].dymz[i]= dydipole[i*3*3+2]+dydipole[i*3*3+2+3]+dydipole[i*3*3+2+6];
               waters[jj]->AtomCenters[kk].dzmx[i]= dzdipole[i*3*3]+dzdipole[i*3*3+3]+dzdipole[i*3*3+6];
               waters[jj]->AtomCenters[kk].dzmy[i]= dzdipole[i*3*3+1]+dzdipole[i*3*3+1+3]+dzdipole[i*3*3+1+6];
               waters[jj]->AtomCenters[kk].dzmz[i]= dzdipole[i*3*3+2]+dzdipole[i*3*3+2+3]+dzdipole[i*3*3+2+6];

/*
               DmuByDx[index]   = dxdipole[i*3*3]+dxdipole[i*3*3+3]+dxdipole[i*3*3+6];
               DmuByDx[index+1] = dxdipole[i*3*3+1]+dxdipole[i*3*3+1+3]+dxdipole[i*3*3+1+6];
               DmuByDx[index+2] = dxdipole[i*3*3+2]+dxdipole[i*3*3+2+3]+dxdipole[i*3*3+2+6];
               DmuByDy[index]   = dydipole[i*3*3]+dydipole[i*3*3+3]+dydipole[i*3*3+6];
               DmuByDy[index+1] = dydipole[i*3*3+1]+dydipole[i*3*3+1+3]+dydipole[i*3*3+1+6];
               DmuByDy[index+2] = dydipole[i*3*3+2]+dydipole[i*3*3+2+3]+dydipole[i*3*3+2+6];
               DmuByDz[index]   = dzdipole[i*3*3]+dzdipole[i*3*3+3]+dzdipole[i*3*3+6];
               DmuByDz[index+1] = dzdipole[i*3*3+1]+dzdipole[i*3*3+1+3]+dzdipole[i*3*3+1+6];
               DmuByDz[index+2] = dzdipole[i*3*3+2]+dzdipole[i*3*3+2+3]+dzdipole[i*3*3+2+6];
               index = index+3;

*/
               DmuByDR[index]   = dxdipole[i*3*3]+dxdipole[i*3*3+3]+dxdipole[i*3*3+6];
               DmuByDR[index+1] = dxdipole[i*3*3+1]+dxdipole[i*3*3+1+3]+dxdipole[i*3*3+1+6];
               DmuByDR[index+2] = dxdipole[i*3*3+2]+dxdipole[i*3*3+2+3]+dxdipole[i*3*3+2+6];
               DmuByDR[3*nwaters + index]   = dydipole[i*3*3]+dydipole[i*3*3+3]+dydipole[i*3*3+6];
               DmuByDR[3*nwaters + index+1] = dydipole[i*3*3+1]+dydipole[i*3*3+1+3]+dydipole[i*3*3+1+6];
               DmuByDR[3*nwaters + index+2] = dydipole[i*3*3+2]+dydipole[i*3*3+2+3]+dydipole[i*3*3+2+6];
               DmuByDR[6*nwaters + index]   = dzdipole[i*3*3]+dzdipole[i*3*3+3]+dzdipole[i*3*3+6];
               DmuByDR[6*nwaters + index+1] = dzdipole[i*3*3+1]+dzdipole[i*3*3+1+3]+dzdipole[i*3*3+1+6];
               DmuByDR[6*nwaters + index+2] = dzdipole[i*3*3+2]+dzdipole[i*3*3+2+3]+dzdipole[i*3*3+2+6];
    
               index = index+3;


            }
//         cout<<" DmuByDR["<<index<<"] :" << DmuByDR[index] <<"\n";

         }
      }





      if (verbose > 1) 
         ReportInducedDipoles();

      delete[] dxfield2z; delete[] dxfield2y; delete[] dxfield2x;
      delete[] dyfield2z; delete[] dyfield2y; delete[] dyfield2x;
      delete[] dzfield2z; delete[] dzfield2y; delete[] dzfield2x;
      delete[] dxfield1z; delete[] dxfield1y; delete[] dxfield1x;
      delete[] dyfield1z; delete[] dyfield1y; delete[] dyfield1x;
      delete[] dzfield1z; delete[] dzfield1y; delete[] dzfield1x;
      delete[] Coeff; delete[] InvCoeff;
      delete[] dxCoeff; delete[] dyCoeff; delete[] dzCoeff;
      delete[] dipolez; delete[] dipoley; delete[] dipolex;
      delete[] zsite;   delete[] ysite;   delete[] xsite; delete[] alpha;
      delete[] dxdipole; delete[] dydipole; delete[] dzdipole;

}


void WaterCluster::Calc_dTenosrDotDipoles(int verbose, double *mu, double *dTensor)
{
   int i, j, k;


   // establish lists with all polarizable sites and the fields at these positions
   // after convergence these data go back into waters[j], and can be used in the energy expressions
   int nsites = nwaters * waters[0]->PolSites;
   int PolSites = waters[0]->PolSites;

   double *alpha = new double[nsites];      //  polarizabilities (isotropic)
   double *xsite = new double[nsites];      //  positions of all sites
   double *ysite = new double[nsites];
   double *zsite = new double[nsites];
   double *dipolex = new double[nsites];    //  current dipoles of all sites
   double *dipoley = new double[nsites];
   double *dipolez = new double[nsites];

   double *dxfield2x = new double[nsites];    //  field due to induced dipoles
   double *dxfield2y = new double[nsites];
   double *dxfield2z = new double[nsites];
   double *dyfield2x = new double[nsites];    
   double *dyfield2y = new double[nsites];
   double *dyfield2z = new double[nsites];
   double *dzfield2x = new double[nsites];    
   double *dzfield2y = new double[nsites];
   double *dzfield2z = new double[nsites];



   double dampCD = waters[0]->Damping[1];
   double dampDD = waters[0]->Damping[2];


   // initialize the polarizable sites
   i = 0;
   for (j = 0; j < nwaters; ++j) {
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            if (verbose > 15)
               cout << "Init polarizable center#" << i << " at water#" << j << " AtomCenter#" << k << endl;
            alpha[i] = a;
            xsite[i] = (waters[j]->AtomCenters[k]).x; 
            ysite[i] = (waters[j]->AtomCenters[k]).y; 
            zsite[i] = (waters[j]->AtomCenters[k]).z;
            dipolex[i] = 0; dipoley[i] = 0; dipolez[i] = 0;
            ++i;
         }
      }
   }
   if (i != nsites) { cout << "WaterCluster::CalcInducedDipoles:  This should not happen.\n"; exit(1);}



   // compute the fields due to the permanent charges --> field1   E = Sum_i Q_i*R_i / |R_i|^3

    
   for (i = 0; i < nwaters*3; ++i) {
      dipolex[i]=mu[i*3+0] ;
      dipoley[i]=mu[i*3+1] ;
      dipolez[i]=mu[i*3+2] ;
   }        


      //
      // Dot product of derivative of Tensor and Dipolemoments 
      //

      double *dxdipole = new double[3*nsites];
      double *dydipole = new double[3*nsites];
      double *dzdipole = new double[3*nsites];

      int BaseR;
      int BaseR2;
      int TotalSites=nsites*9;
      int Tindex=0;


      int index=0;
      for (int jj = 0; jj < nwaters; ++jj) {
         for (int kk = 0; kk < waters[jj]->nAtomCenters; ++kk) {
            if (kk ==3) BaseR=-1; 
            else BaseR = jj*3+kk;
            BaseR2 = jj*waters[jj]->nAtomCenters+kk;

            for (int ii = 0; ii < nsites; ++ii) {
               dxfield2x[ii] = 0.; dxfield2y[ii] = 0.; dxfield2z[ii] = 0.;
               dyfield2x[ii] = 0.; dyfield2y[ii] = 0.; dyfield2z[ii] = 0.;
               dzfield2x[ii] = 0.; dzfield2y[ii] = 0.; dzfield2z[ii] = 0.;
            }


            for (i = 0;  i < nsites; ++i) {             // at site i
               for (j = 0;  j < nsites; ++j) {           // dipole j creates field E
                  if (j == i) 
                     continue;

                  if(BaseR == i || BaseR==j){
                     double dxEx, dxEy, dxEz;
                     double dyEx, dyEy, dyEz;
                     double dzEx, dzEy, dzEz;

                     dFieldOfDipole(dipolex[j], dipoley[j], dipolez[j], xsite[j], ysite[j], zsite[j],
                        xsite[i], ysite[i], zsite[i], dxEx, dxEy, dxEz,dyEx,dyEy,dyEz,
                        dzEx,dzEy,dzEz, alpha[i],alpha[j], dampDD);
                     if (BaseR == i){
                        dxfield2x[i] += dxEx; 
                        dxfield2y[i] += dxEy; 
                        dxfield2z[i] += dxEz; 
                        dyfield2x[i] += dyEx; 
                        dyfield2y[i] += dyEy; 
                        dyfield2z[i] += dyEz; 
                        dzfield2x[i] += dzEx; 
                        dzfield2y[i] += dzEy; 
                        dzfield2z[i] += dzEz; 
                     }
                     else if (BaseR == j){
                        dxfield2x[i] -= dxEx; 
                        dxfield2y[i] -= dxEy; 
                        dxfield2z[i] -= dxEz; 
                        dyfield2x[i] -= dyEx; 
                        dyfield2y[i] -= dyEy; 
                        dyfield2z[i] -= dyEz; 
                        dzfield2x[i] -= dzEx; 
                        dzfield2y[i] -= dzEy; 
                        dzfield2z[i] -= dzEz; 
                     }
                  }
               }
            }


            for (i = 0;  i < nsites; ++i) {
               dTensor[TotalSites*Tindex + i*3]   =  dxfield2x[i];
               dTensor[TotalSites*Tindex + i*3+1] =  dxfield2y[i];
               dTensor[TotalSites*Tindex + i*3+2] =  dxfield2z[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3]   = dyfield2x[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3+1] = dyfield2y[i];
               dTensor[TotalSites*Tindex + 3*nsites + i*3+2] = dyfield2z[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3]   = dzfield2x[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3+1] = dzfield2y[i];
               dTensor[TotalSites*Tindex + 6*nsites + i*3+2] = dzfield2z[i];

            }
            Tindex = Tindex+1;
         }
      }


      delete[] dxfield2z; delete[] dxfield2y; delete[] dxfield2x;
      delete[] dyfield2z; delete[] dyfield2y; delete[] dyfield2x;
      delete[] dzfield2z; delete[] dzfield2y; delete[] dzfield2x;
      delete[] dipolez; delete[] dipoley; delete[] dipolex;
      delete[] zsite;   delete[] ysite;   delete[] xsite; delete[] alpha;
      delete[] dxdipole; delete[] dydipole; delete[] dzdipole;

}



/////////////////////////////////////////////////////////////////////
//
//  CalcIntermolecularPotential step#2 of the water potential
//
//  compute the intermolecular energy, i.e, sum all
//  - point-charge point-charge interactions  (electrostatics)
//  - point-charge induced-dipole interactions (polarization)
//  - van-der-Waals terms (dispersion)
//
//
double WaterCluster::CalcIntermolecularPotential(int verbose)
{
   if (verbose > 1)
      cout << "\nIntermolecular interactions in eV and kcal/mol:\n";
   NuclearRepulsionEnergy = 0; // these three are public member variables of Cluster with a somewhat misleading names
   LJEnergy = 0;               // the 1st is electrostatics+polarization, the 2nd is the van-der-Waals interaction
   ElectrostaticEnergy = 0;    // the 3rd is only the electrostatic contribution
   PolarizationEnergy = 0;    // water-water polarization energy

   double elec=0;		//used to output the total elec. only
   double pol=0;		//output the total polarizaiton only
   double vdw=0;		//output the total dispersion only


   for (int i = 0; i < nwaters - 1; ++i) {
      for (int j = i+1; j < nwaters; ++j) {
         double electrostatic = WaterWaterElectrostaticEnergy(*waters[i], *waters[j], verbose);
         double polarization = WaterWaterPolarizationEnergy(*waters[i], *waters[j], verbose);
         double dispersion = WaterWaterDispersionEnergy(*waters[i], *waters[j], verbose);
         ElectrostaticEnergy += electrostatic;
         NuclearRepulsionEnergy += electrostatic + polarization;
         LJEnergy += dispersion;

         elec += electrostatic;
         pol  += polarization;
         vdw += dispersion;

         if (verbose > 2)
            printf("  #%-2d -- #%-2d  Ves = %12.8f(%11.7f)  Vpol = %12.8f(%11.7f)  Vdis = %12.8f(%11.7f)\n", 
            i, j, electrostatic*AU2EV, electrostatic*AU2KCAL, polarization*AU2EV, polarization*AU2KCAL,
            dispersion*AU2EV, dispersion*AU2KCAL); 
      }
   }
   if (verbose > 1) {
      double bind=elec+pol+vdw;
      printf("  -----------------------------\n");
      printf("  Total:  Ves = %12.8f(%11.7f)  Vpol = %12.8f(%11.7f)  Vdis = %12.8f(%11.7f)\n", 
         elec*AU2EV, elec*AU2KCAL, pol*AU2EV, pol*AU2KCAL, vdw*AU2EV, vdw*AU2KCAL); 
      printf("  Total energy of the neutral cluster:  %12.8f eV = %12.8f kcal\n", bind*AU2EV, bind*AU2KCAL);
   }
 //   cout<<"total energy "<<NuclearRepulsionEnergy+LJEnergy<<endl;
 //   cout<<"polarization energy "<<pol<<endl;
    cout<<"ElectrostaticEnergy = "<<ElectrostaticEnergy<<endl;
    PolarizationEnergy=pol;

   return NuclearRepulsionEnergy + LJEnergy;
}


/////////////////////////////////////////////////////////////////////
//
//  CalcIntermolecularGradient ------ March 06, 2007 by Tae Hoon Choi
//
//  compute the intermolecular gradient, i.e, sum all
//  - point-charge point-charge interactions  (electrostatics)
//  - point-charge induced-dipole interactions (polarization)
//  - induced-dipole induced-dipole interactions (polarization)
//  - van-der-Waals terms (dispersion)
//

double WaterCluster::CalcIntermolecularGradient(double *Gradient, double *PolGrad, int verbose)
{

   double Grad1[12];
   double Grad2[12];

   if (verbose > 1)
      cout << "\nIntermolecular interactions in eV and kcal/mol:\n";
   NuclearRepulsionEnergy = 0; // these three are public member variables of Cluster with a somewhat misleading names
   LJEnergy = 0;               // the 1st is electrostatics+polarization, the 2nd is the van-der-Waals interaction
   ElectrostaticEnergy = 0;    // the 3rd is only the electrostatic contribution

   double elec=0;		//used to output the total elec. only
   double pol=0;		//output the total polarizaiton only
   double vdw=0;		//output the total dispersion only
   int firstCheck=0; 

   for (int k = 0; k < 4 ; ++k) {
      Grad1[k*3]=0.00; Grad1[k*3+1]=0.00; Grad1[k*3+2]=0.00;
      Grad2[k*3]=0.00; Grad2[k*3+1]=0.00; Grad2[k*3+2]=0.00;
   }

   for (int i = 0; i < nwaters - 1; ++i) {
      for (int j = i+1; j < nwaters; ++j) {
         double polarization = WaterWaterPolarizationGradient(*waters[i], *waters[j],Grad1,Grad2, firstCheck, verbose);
         for (int k=0; k < 4 ; ++k){
            PolGrad[(i*4+k)*3] += Grad1[k*3];
            PolGrad[(i*4+k)*3+1] += Grad1[k*3+1];
            PolGrad[(i*4+k)*3+2] += Grad1[k*3+2];
            PolGrad[(j*4+k)*3] += Grad2[k*3];
            PolGrad[(j*4+k)*3+1] += Grad2[k*3+1];
            PolGrad[(j*4+k)*3+2] += Grad2[k*3+2];
         }

         double electrostatic = WaterWaterElectrostaticGradient(*waters[i], *waters[j], Grad1,Grad2, verbose);
         double dispersion = WaterWaterDispersionGradient(*waters[i], *waters[j], Grad1,Grad2, verbose);
         //      double electrostatic = WaterWaterElectrostaticEnergy(*waters[i], *waters[j], verbose);
         //      double dispersion = WaterWaterDispersionEnergy(*waters[i], *waters[j], verbose);
         ElectrostaticEnergy += electrostatic;
         NuclearRepulsionEnergy += electrostatic + polarization;
         LJEnergy += dispersion;

         elec += electrostatic;
         pol  += polarization;
         vdw += dispersion;

         if (verbose > 2)
            printf("  #%-2d -- #%-2d  Ves = %12.8f(%11.7f)  Vpol = %12.8f(%11.7f)  Vdis = %12.8f(%11.7f)\n", 
            i, j, electrostatic*AU2EV, electrostatic*AU2KCAL, polarization*AU2EV, polarization*AU2KCAL,
            dispersion*AU2EV, dispersion*AU2KCAL); 
         // 4 is nAtomCenters (should be the same for every water member)
         for (int k=0; k < 4 ; ++k){
            Gradient[(i*4+k)*3] += Grad1[k*3]; 
            Gradient[(i*4+k)*3+1] += Grad1[k*3+1]; 
            Gradient[(i*4+k)*3+2] += Grad1[k*3+2];
            Gradient[(j*4+k)*3] += Grad2[k*3]; 
            Gradient[(j*4+k)*3+1] += Grad2[k*3+1]; 
            Gradient[(j*4+k)*3+2] += Grad2[k*3+2];
            Grad1[k*3]=0.00; Grad1[k*3+1]=0.00; Grad1[k*3+2]=0.00;
            Grad2[k*3]=0.00; Grad2[k*3+1]=0.00; Grad2[k*3+2]=0.00;
         }
         if ( i==0 && j == i+1) firstCheck =1;
      }
      if ( i==0) firstCheck = 2;
   }

   //  for (int i = 0; i < 8 ; ++i) 
   //    cout<<"Gradient="<<Gradient[i*3]<<" "<<Gradient[i*3+1]<<" "<<Gradient[i*3+2]<<endl;



   if (verbose > 1) {
      double bind=elec+pol+vdw;
      printf("  -7777----------------------------\n");
      printf("  Total:  Ves = %12.8f(%11.7f)  Vpol = %12.8f(%11.7f)  Vdis = %12.8f(%11.7f)\n", 
         elec*AU2EV, elec*AU2KCAL, pol*AU2EV, pol*AU2KCAL, vdw*AU2EV, vdw*AU2KCAL); 
      printf("  Total energy of the neutral cluster:  %12.8f eV = %12.8f kcal\n", bind*AU2EV, bind*AU2KCAL);
   }

   return NuclearRepulsionEnergy + LJEnergy;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   AAD 4-11-06
//      ttm2-u6 no longer damps charge-charge
//


double ElectrostaticDCDamp(double distance, double alphai, double alphaj, double a)
{
   double A=pow(alphai*alphaj,1.0/6.0);
   double x=a*pow(distance/A,3);
      return 1.0-exp(-x);
   //return 1.0;
}

double ElectrostaticDDDamp(double distance, double alphai, double alphaj, double a)
{
   double A=pow(alphai*alphaj,1.0/6.0);
   double x=a*pow(distance/A,3);
   //return -3.0*x*exp(-x)+2.0*(1.0-exp(-x));  //how it looks in thole, not quite correct for us
   return 1.0-(1+x)*exp(-x);
   //return 1.0;
}


/////////////////////////////////////////////////////////////////////
//
//  WaterWaterElectrostaticEnergy
//
//  this is the point-charge/point-charge interaction, i.e., the modeled electrostatic interaction 
//
double WaterWaterElectrostaticEnergy(Water &W1, Water &W2, int verbose)
{
   double distance = 0;
   double energy = 0;
   int count = 0;
   for (int i = 0; i < W1.nAtomCenters; ++i) {
      double x1 = W1.AtomCenters[i].x;
      double y1 = W1.AtomCenters[i].y;
      double z1 = W1.AtomCenters[i].z;
      double q1 = W1.AtomCenters[i].charge;
      //double alpha1 = W1.AtomCenters[i].sigma;
      if (q1 != 0) 
      {
         for (int j = 0; j < W2.nAtomCenters; ++j) {
            double dx = x1 - W2.AtomCenters[j].x;
            double dy = y1 - W2.AtomCenters[j].y;
            double dz = z1 - W2.AtomCenters[j].z;
            double q2 = W2.AtomCenters[j].charge;
            //double alpha2 = W2.AtomCenters[j].sigma;
            if (q2 != 0)
            {
               //			cout << "-------------hey!!!!!!!---------------" << endl;
               //			cout << "q1=" << q1 << endl << "q2=" << q2 << endl;
               //			cout << "---------------done-------------------" << endl << endl;
               distance=sqrt(dx*dx + dy*dy + dz*dz);
               // cout << "############## hey!!!!! interaction "<< count+1 << " ###########" << endl;
               // cout << "# site 1:  " << i << "  and  " << "site 2:  " << j << endl;
               //  cout << "# -----------------------------------" << endl;

               //double damp=1.0;
               //damp=1.0;
               // energy += damp*(q1*q2/distance);
               energy += (q1*q2 / distance);


               //  cout << "# site 1 charge: " << W1.AtomCenters[i].charge << "	" << " site 2 charge: " << W2.AtomCenters[j].charge << endl;
               //     printf(" distance   %20.14f\n",distance*Bohr2Angs/10.);
               //  cout << "------------" << endl << " Energy:	" << (q1*q2)/distance << endl;
               //  cout << "# Damped Energy: " << damp*(q1*q2)/ distance << endl;
               //  cout << "############### done #################" << endl << endl << endl;
               count++;
            }
         }
      }
   }
   //cout << "----------hey!!!-----------" << endl;
   //cout << "I ran " << count << " electrostatic interactions" << endl;
   //cout << "----------done-------------" << endl;
   return energy;
}



/////////////////////////////////////////////////////////////////////
//
//  WaterWaterElectrostaticGradient --------- March 06, 2007 by Tae Hoon Choi
//
//  this is the point-charge/point-charge interaction, i.e., the modeled electrostatic interaction 
//
double WaterWaterElectrostaticGradient(Water &W1, Water &W2, double *Grad1, double *Grad2,int verbose)
{
   double distance = 0;
   double energy = 0;
   double de = 0.00;
   int count = 0;
   for (int i = 0; i < W1.nAtomCenters; ++i) {
      double x1 = W1.AtomCenters[i].x;
      double y1 = W1.AtomCenters[i].y;
      double z1 = W1.AtomCenters[i].z;
      double q1 = W1.AtomCenters[i].charge;
      //double alpha1 = W1.AtomCenters[i].sigma;
      if (q1 != 0) 
      {
         for (int j = 0; j < W2.nAtomCenters; ++j) {
            double dx = x1 - W2.AtomCenters[j].x;
            double dy = y1 - W2.AtomCenters[j].y;
            double dz = z1 - W2.AtomCenters[j].z;
            double q2 = W2.AtomCenters[j].charge;
            //double alpha2 = W2.AtomCenters[j].sigma;
            if (q2 != 0)
            {
               distance=sqrt(dx*dx + dy*dy + dz*dz);

               energy += (q1*q2 / distance);
               de = -q1*q2 / pow(distance,3);
               Grad1[i*3]   = Grad1[i*3]   + de*dx;
               Grad1[i*3+1] = Grad1[i*3+1] + de*dy;
               Grad1[i*3+2] = Grad1[i*3+2] + de*dz;
               Grad2[j*3]   = Grad2[j*3]   - de*dx;
               Grad2[j*3+1] = Grad2[j*3+1] - de*dy;
               Grad2[j*3+2] = Grad2[j*3+2] - de*dz;

               count++;
            }
         }
      }
   }
   return energy;
}




/////////////////////////////////////////////////////////////////////
//
//  WaterWaterPolariztionEnergy
//
//  this is the point-charge/induced dipole interaction, i.e., the modeled polarization interaction
//
//  all terms have the form Mu*R / |R|^3 
//
//
//
double WaterWaterPolarizationEnergy(Water &W1, Water &W2, int verbose)
{
   double energy = 0;
   int j;
   int counter=0;
   //  cout << " water 1 has " << W1.nAtomCenters << " and water 2 has " << W2.nAtomCenters << " centers" << endl;
   for (int k = 0; k < W1.nAtomCenters; ++k) {
      double x1, y1, z1, x2, y2, z2;        // positions of the polarizable sites for both waters
      double m1x, m1y, m1z, m2x, m2y, m2z;  // induced dipoles at these sites
      double a = W1.AtomCenters[k].alpha;
      double b = W2.AtomCenters[k].alpha;
      double damping = W1.Damping[1];
      if (a != 0) {
         x1 = W1.AtomCenters[k].x;
         y1 = W1.AtomCenters[k].y;
         z1 = W1.AtomCenters[k].z;
         x2 = W2.AtomCenters[k].x;
         y2 = W2.AtomCenters[k].y;
         z2 = W2.AtomCenters[k].z;
         m1x = W1.AtomCenters[k].mx;
         m1y = W1.AtomCenters[k].my;
         m1z = W1.AtomCenters[k].mz;
         m2x = W2.AtomCenters[k].mx;
         m2y = W2.AtomCenters[k].my;
         m2z = W2.AtomCenters[k].mz;
         // W2 polarizing W1
         for (j = 0; j < W2.nAtomCenters; ++j) {
            double Rx = x1 - W2.AtomCenters[j].x;
            double Ry = y1 - W2.AtomCenters[j].y;
            double Rz = z1 - W2.AtomCenters[j].z;
            double q2 = W2.AtomCenters[j].charge;
            double a2 = W2.AtomCenters[j].sigma; //polarizability used in damping
            if (q2 != 0) {
               double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
               double damp = ElectrostaticDCDamp(distance,a,a2,damping);
               energy -= q2 * (m1x*Rx + m1y*Ry + m1z*Rz) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;

               //	  cout<<" first energy = "<<energy<<" m1x = "<<m1x<<" m1y = "<<m1y<<" m1z = "<<m1z<<"\n";
               counter++;
               //double eh = -q2 * (m1x*Rx + m1y*Ry + m1z*Rz) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5);
               //cout << " ---- Calculating Polarization Energy ---- " << endl;
               //cout << "   # first loop interaction " << counter << "   # " << endl;
               //cout << "     Water 1, site " << k << "  at " << x1*Bohr2Angs << " " << y1*Bohr2Angs << " " << z1*Bohr2Angs << endl;
               //cout << "	dipole: " << m1x*AU2Debye << " " << m1y*AU2Debye << " " << m1z*AU2Debye << endl;
               //cout << "     Water 2, site " << j << "  at " << W2.AtomCenters[j].x*Bohr2Angs << " " << W2.AtomCenters[j].y*Bohr2Angs << " " << W2.AtomCenters[j].z*Bohr2Angs << endl;
               //cout << "	charge: " << q2 << endl;
               //cout << endl;
               //cout << "  dist:   " << distance << endl;
               //cout << "  damp:   " << damp << endl;
               //cout << "  Energy: " << eh*AU2KCAL << endl;
               //cout << " ----------------------------------------- " << endl;
            }
         }
         // W1 polarizing W2
         for (j = 0; j < W1.nAtomCenters; ++j) {
            double Rx = x2 - W1.AtomCenters[j].x;
            double Ry = y2 - W1.AtomCenters[j].y;
            double Rz = z2 - W1.AtomCenters[j].z;
            double q1 = W1.AtomCenters[j].charge;
            double a1 = W1.AtomCenters[j].sigma; //polarizability used in damping
            if (q1 != 0) {
               double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
               double damp = ElectrostaticDCDamp(distance,b,a1,damping);
               energy -= q1 * (m2x*Rx + m2y*Ry + m2z*Rz) * damp * pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;

               //	  cout<<" second energy = "<<energy<<" m2x = "<<m2x<<" m2y = "<<m2y<<" m2z = "<<m2z<<"\n";
               counter++;
               //double eh = -q1 * (m2x*Rx + m2y*Ry + m2z*Rz) * damp * pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5);
               //cout << " ---- Calculating Polarization Energy ---- " << endl;
               //cout << "   # second loop interaction " << counter << "   # " << endl;
               //cout << "     Water 1, site " << j << "  at " << W1.AtomCenters[j].x*Bohr2Angs << " " << W1.AtomCenters[j].y*Bohr2Angs << " " << W1.AtomCenters[j].z*Bohr2Angs << endl;
               //cout << "       charge: " << q1 << endl;
               //cout << "     Water 2, site " << k << "  at " << x2*Bohr2Angs << " " << y2*Bohr2Angs << " " << z2*Bohr2Angs << endl;
               //cout << "       dipole: " << m2x*AU2Debye << " " << m2y*AU2Debye << " " << m2z*AU2Debye << endl;
               //cout << endl;
               //cout << "  dist:   " << distance << endl;
               //cout << " dist.dip:" << (m2x*Rx + m2y*Ry + m2z*Rz) << endl;
               //cout << "  damp:   " << damp << endl;
               //cout << "  Energy: " << eh*AU2KCAL << endl;
               //cout << " ----------------------------------------- " << endl;

            }
         }
      }
   }
   //cout << "####### There are " << counter << " polarization interactions #######" << endl;
   return energy;
}

/////////////////////////////////////////////////////////////////////
//
//  WaterWaterPolariztionGradient  ---- March 6, 2007 by Tae Hoon Choi
//
//  this is the point-charge/induced dipole interaction gradient
//
//  
//
double WaterWaterPolarizationGradient(Water &W1, Water &W2, double *Grad1, double *Grad2, int firstCheck, int verbose)
{
   double energy = 0;
   int j;
   int counter=0;
   double de0=0.00;
   double de1=0.00;
   double de2=0.00;
   double de3=0.00;
   double de4=0.00;
   double de5=0.00;
   double de6=0.00;
   double de7=0.00;
   double de8=0.00;
   //cout << " water 1 has " << W1.nAtomCenters << " and water 2 has " << W2.nAtomCenters << " centers" << endl;
   for (int k = 0; k < W1.nAtomCenters; ++k) {
      double x1, y1, z1, x2, y2, z2;        // positions of the polarizable sites for both waters
      double m1x, m1y, m1z, m2x, m2y, m2z;  // induced dipoles at these sites
      double m1xa, m1ya, m1za, m2xa, m2ya, m2za;  // induced dipoles at these sites
      double a = W1.AtomCenters[k].alpha;
      double b = W2.AtomCenters[k].alpha;
      //cout<<" k = "<<k<<" avalue= "<<a<<" b = "<<b<<"\n";

      double damping = W1.Damping[1];
      double damping2= W1.Damping[2];
      if (a != 0) {
         x1 = W1.AtomCenters[k].x;
         y1 = W1.AtomCenters[k].y;
         z1 = W1.AtomCenters[k].z;
         x2 = W2.AtomCenters[k].x;
         y2 = W2.AtomCenters[k].y;
         z2 = W2.AtomCenters[k].z;
         m1x = W1.AtomCenters[k].mx;
         m1y = W1.AtomCenters[k].my;
         m1z = W1.AtomCenters[k].mz;
         m2x = W2.AtomCenters[k].mx;
         m2y = W2.AtomCenters[k].my;
         m2z = W2.AtomCenters[k].mz;
         // 
         //    W1 and W1 interaction part 
         //    This part is calculated just one time
         //
         if (firstCheck == 0 ){

            for (j = k+1; j < W1.nAtomCenters; ++j) {
               double Rx = x1 - W1.AtomCenters[j].x;
               double Ry = y1 - W1.AtomCenters[j].y;
               double Rz = z1 - W1.AtomCenters[j].z;
               double a2 = W1.AtomCenters[j].sigma; //polarizability used in damping
               double a3 = W1.AtomCenters[j].alpha; //polarizability used in damping
               if (a2 == 0) 
                  a2 = W1.AtomCenters[j].alpha; 
               m1xa = W1.AtomCenters[j].mx;
               m1ya = W1.AtomCenters[j].my;
               m1za = W1.AtomCenters[j].mz;
               double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
               double damp = ElectrostaticDCDamp(distance,a,a2,damping2);
               //
               //   induced dipole and induced dipole interaction part 
               //
               if (a3 != 0){
                  de0 = -2.0* (m1x*m1xa+m1y*m1ya+m1z*m1za) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0/2.0;
                  de4 = -2.0* (m1x*Rx+m1y*Ry+m1z*Rz)*(m1xa*Rx+m1ya*Ry+m1za*Rz)*damping2*damping2*3.0*3.0*(1-damp)/distance/(a*a2)/2.0;
                  de5 = 2.0* (m1x*m1xa+m1y*m1ya+m1z*m1za) *damping2*3.0*pow(a*a2,-0.5)*(1-damp)/(Rx*Rx + Ry*Ry + Rz*Rz)/2.0;
                  de6 = 2.0*(m1x*Rx+m1y*Ry+m1z*Rz)*(m1xa*Rx+m1ya*Ry+m1za*Rz)*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp))*15.0*pow((Rx*Rx + Ry*Ry + Rz*Rz),-3.5)/2.0;
                  de7=-2.0*(m1xa*Rx+m1ya*Ry+m1za*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp));
                  de8=-2.0*(m1x*Rx+m1y*Ry+m1z*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp));
               }
               else {de0= 0.00; de4=0.00; de5=0.00; de6=0.00; de7=0.00; de8=0.00;}
               Grad1[k*3]   = Grad1[k*3]   + (de0+de4+de5+de6)*Rx + de7*m1x+de8*m1xa;
               Grad1[k*3+1] = Grad1[k*3+1] + (de0+de4+de5+de6)*Ry + de7*m1y+de8*m1ya; 
               Grad1[k*3+2] = Grad1[k*3+2] + (de0+de4+de5+de6)*Rz + de7*m1z+de8*m1za;
               Grad1[j*3]   = Grad1[j*3]   - (de0+de4+de5+de6)*Rx - de7*m1x-de8*m1xa;
               Grad1[j*3+1] = Grad1[j*3+1] - (de0+de4+de5+de6)*Ry - de7*m1y-de8*m1ya; 
               Grad1[j*3+2] = Grad1[j*3+2] - (de0+de4+de5+de6)*Rz - de7*m1z-de8*m1za;

            }
         }


         // W2 polarizing W1
         for (j = 0; j < W2.nAtomCenters; ++j) {
            double Rx = x1 - W2.AtomCenters[j].x;
            double Ry = y1 - W2.AtomCenters[j].y;
            double Rz = z1 - W2.AtomCenters[j].z;
            double q2 = W2.AtomCenters[j].charge;
            double a2 = W2.AtomCenters[j].sigma; //polarizability used in damping
            double a3 = W2.AtomCenters[j].alpha; //polarizability used in damping
            if (a2 == 0) 
               a2 = W2.AtomCenters[j].alpha; 
            m1xa = W2.AtomCenters[j].mx;
            m1ya = W2.AtomCenters[j].my;
            m1za = W2.AtomCenters[j].mz;
            double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
            double damp = ElectrostaticDCDamp(distance,a,a2,damping2);
            //	
            //     
            //   induced dipole and induced dipole interaction part 
            //
            if (a3 != 0){
               de0 = -2.0* (m1x*m1xa+m1y*m1ya+m1z*m1za) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0/2.0;
               de4 = -2.0* (m1x*Rx+m1y*Ry+m1z*Rz)*(m1xa*Rx+m1ya*Ry+m1za*Rz)*damping2*damping2*3.0*3.0*(1-damp)/distance/(a*a2)/2.0;
               de5 = 2.0* (m1x*m1xa+m1y*m1ya+m1z*m1za) *damping2*3.0*pow(a*a2,-0.5)*(1-damp)/(Rx*Rx + Ry*Ry + Rz*Rz)/2.0;
               de6 = 2.0*(m1x*Rx+m1y*Ry+m1z*Rz)*(m1xa*Rx+m1ya*Ry+m1za*Rz)*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp))*15.0*pow((Rx*Rx + Ry*Ry + Rz*Rz),-3.5)/2.0;
               de7=-2.0*(m1xa*Rx+m1ya*Ry+m1za*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp));
               de8=-2.0*(m1x*Rx+m1y*Ry+m1z*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(a*a2,-0.5)*(1-damp));
            }
            else {de0= 0.00; de4=0.00; de5=0.00; de6=0.00; de7=0.00; de8=0.00;}
            //
            //    Charge - induced dipole interaction term
            //
            damp = ElectrostaticDCDamp(distance,a,a2,damping);
            if (q2 != 0) {
               energy -= q2 * (m1x*Rx + m1y*Ry + m1z*Rz) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;
               de1 = -q2 * 2.0*damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;
               de2 = -q2 * 2.0*(m1x*Rx + m1y*Ry + m1z*Rz) * 3.0*damping*pow(a*a2,-0.5)*(1.0-damp)/(Rx*Rx + Ry*Ry + Rz*Rz)/2.0;
               de3 =  q2 * 2.0*(m1x*Rx + m1y*Ry + m1z*Rz) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0/2.0;
               counter++;
            }
            else {de1= 0.00; de2=0.00; de3=0.00;}

            Grad1[k*3]   = Grad1[k*3]   + de1*m1x+(de2+de3+de0+de4+de5+de6)*Rx + de7*m1x+de8*m1xa;
            Grad1[k*3+1] = Grad1[k*3+1] + de1*m1y+(de2+de3+de0+de4+de5+de6)*Ry + de7*m1y+de8*m1ya; 
            Grad1[k*3+2] = Grad1[k*3+2] + de1*m1z+(de2+de3+de0+de4+de5+de6)*Rz + de7*m1z+de8*m1za;
            Grad2[j*3]   = Grad2[j*3]   - de1*m1x-(de2+de3+de0+de4+de5+de6)*Rx - de7*m1x-de8*m1xa;
            Grad2[j*3+1] = Grad2[j*3+1] - de1*m1y-(de2+de3+de0+de4+de5+de6)*Ry - de7*m1y-de8*m1ya; 
            Grad2[j*3+2] = Grad2[j*3+2] - de1*m1z-(de2+de3+de0+de4+de5+de6)*Rz - de7*m1z-de8*m1za;
            //cout<<" Grad1[k*3]= "<<Grad1[k*3]<<" m1x = "<<m1x<<" m1y = "<<m1y<<" m1z="<<m1z<<"\n"; 
            //cout<<" Rx = "<<Rx<<" Ry = "<<Ry<<" Rz="<<Rz<<" tempA = "<<tempA<<" de1 = "<<de1<<"\n"; 

         }


         // W2 and W2 interaction term
         if (firstCheck == 0 || firstCheck == 1){

            for (j = k+1; j < W2.nAtomCenters; ++j) {
               double Rx = x2 - W2.AtomCenters[j].x;
               double Ry = y2 - W2.AtomCenters[j].y;
               double Rz = z2 - W2.AtomCenters[j].z;
               double a1 = W2.AtomCenters[j].sigma; //polarizability used in damping
               double a3 = W2.AtomCenters[j].alpha; //polarizability used in damping
               if (a1 == 0) 
                  a1 = W2.AtomCenters[j].alpha; 
               m2xa = W2.AtomCenters[j].mx;
               m2ya = W2.AtomCenters[j].my;
               m2za = W2.AtomCenters[j].mz;
               double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
               double damp = ElectrostaticDCDamp(distance,b,a1,damping2);
               if (a3 != 0){
                  de0 = -2.0* (m2x*m2xa+m2y*m2ya+m2z*m2za) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0/2.0;
                  de4 = -2.0* (m2x*Rx+m2y*Ry+m2z*Rz)*(m2xa*Rx+m2ya*Ry+m2za*Rz)*damping2*damping2*3.0*3.0*(1-damp)/distance/(b*a1)/2.0;
                  de5 = 2.0* (m2x*m2xa+m2y*m2ya+m2z*m2za) *damping2*3.0*pow(b*a1,-0.5)*(1-damp)/(Rx*Rx + Ry*Ry + Rz*Rz)/2.0;
                  de6 = 2.0*(m2x*Rx+m2y*Ry+m2z*Rz)*(m2xa*Rx+m2ya*Ry+m2za*Rz)*(damp-damping2*pow(distance,3)*pow(b*a1,-0.5)*(1-damp))*15.0*pow((Rx*Rx + Ry*Ry + Rz*Rz),-3.5)/2.0;
                  de7=-2.*(m2xa*Rx+m2ya*Ry+m2za*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(b*a1,-0.5)*(1-damp));
                  de8=-2.*(m2x*Rx+m2y*Ry+m2z*Rz)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5)*3.0/2.0*(damp-damping2*pow(distance,3)*pow(b*a1,-0.5)*(1-damp));
               }
               else {de0= 0.00; de4=0.00; de5=0.00; de6=0.00; de7=0.00; de8=0.00;}
               Grad2[k*3]   = Grad2[k*3]   + (de0+de4+de5+de6)*Rx + de7*m2x+de8*m2xa;
               Grad2[k*3+1] = Grad2[k*3+1] + (de0+de4+de5+de6)*Ry + de7*m2y+de8*m2ya; 
               Grad2[k*3+2] = Grad2[k*3+2] + (de0+de4+de5+de6)*Rz + de7*m2z+de8*m2za;
               Grad2[j*3]   = Grad2[j*3]   - (de0+de4+de5+de6)*Rx - de7*m2x-de8*m2xa;
               Grad2[j*3+1] = Grad2[j*3+1] - (de0+de4+de5+de6)*Ry - de7*m2y-de8*m2ya; 
               Grad2[j*3+2] = Grad2[j*3+2] - (de0+de4+de5+de6)*Rz - de7*m2z-de8*m2za;
               //cout<<" second Grad1[k*3]= "<<Grad1[k*3]<<"\n"; 

            }
         }

         // W1 polarizing W2
         for (j = 0; j < W1.nAtomCenters; ++j) {
            double Rx = x2 - W1.AtomCenters[j].x;
            double Ry = y2 - W1.AtomCenters[j].y;
            double Rz = z2 - W1.AtomCenters[j].z;
            double q1 = W1.AtomCenters[j].charge;
            double a1 = W1.AtomCenters[j].sigma; //polarizability used in damping
            if (a1 == 0) 
               a1 = W1.AtomCenters[j].alpha; 
            m2xa = W1.AtomCenters[j].mx;
            m2ya = W1.AtomCenters[j].my;
            m2za = W1.AtomCenters[j].mz;
            double distance = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
            double damp = ElectrostaticDCDamp(distance,b,a1,damping);
            de0= 0.00; de4=0.00; de5=0.00; de6=0.00; de7=0.00; de8=0.00;

            if (q1 != 0) {
               energy -= q1 * (m2x*Rx + m2y*Ry + m2z*Rz) * damp * pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;
               de1 = -q1 * 2.0*damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -1.5)/2.0;
               de2 = -q1 * 2.0*(m2x*Rx + m2y*Ry + m2z*Rz) * 3.0*damping*pow(b*a1,-0.5)*(1.0-damp)/(Rx*Rx + Ry*Ry + Rz*Rz)/2.0;
               de3 =  q1 * 2.0*(m2x*Rx + m2y*Ry + m2z*Rz) * damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*3.0/2.0;
               counter++;
            }
            else {de1= 0.00; de2=0.00; de3=0.00;}
            Grad2[k*3]   = Grad2[k*3]   + de1*m2x+(de2+de3+de0+de4+de5+de6)*Rx + de7*m2x+de8*m2xa;
            Grad2[k*3+1] = Grad2[k*3+1] + de1*m2y+(de2+de3+de0+de4+de5+de6)*Ry + de7*m2y+de8*m2ya; 
            Grad2[k*3+2] = Grad2[k*3+2] + de1*m2z+(de2+de3+de0+de4+de5+de6)*Rz + de7*m2z+de8*m2za;
            Grad1[j*3]   = Grad1[j*3]   - de1*m2x-(de2+de3+de0+de4+de5+de6)*Rx - de7*m2x-de8*m2xa;
            Grad1[j*3+1] = Grad1[j*3+1] - de1*m2y-(de2+de3+de0+de4+de5+de6)*Ry - de7*m2y-de8*m2ya; 
            Grad1[j*3+2] = Grad1[j*3+2] - de1*m2z-(de2+de3+de0+de4+de5+de6)*Rz - de7*m2z-de8*m2za;


         }
      }
   }
   //cout << "####### There are " << counter << " polarization interactions #######" << endl;
   return energy;
}




/////////////////////////////////////////////////////////////////////
//
//  WaterWaterDispersionEnergy
//
//  this is long-range R^6 (really dispersion) plus a short-range repulsion
//  the functional form and the parameters come from Albert's improved TTM 
//
//   AAD 4-11-06
//      ttm2-u6 has damped dispersion and all-atom repulsion
//         tang-toennies damping function with delta=2.23 
//
//
//
double WaterWaterDispersionEnergy(Water &W1, Water &W2, int verbose)
{


   double energy = 0.0;

   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         double dx = W1.AtomCenters[i].x - W2.AtomCenters[j].x; 
         double dy = W1.AtomCenters[i].y - W2.AtomCenters[j].y;
         double dz = W1.AtomCenters[i].z - W2.AtomCenters[j].z;
         double R2 = dx*dx + dy*dy + dz*dz;
         double R = sqrt(R2);
         double delta = W1.DispersionOO[1];

         if ( (i == 0) && (j ==0) )
         {
            double R6 = R2*R2*R2;

            double damp = 1-(1.0 + delta*R + pow(delta*R,2)/2.0 + pow(delta*R,3)/6.0 + pow(delta*R,4)/24.0 + pow(delta*R,5)/120.0 + pow(delta*R,6)/720.0) * exp(-delta*R);
            //cout << "DAMP: " << damp;

            //cout << "DISP: " << damp*W1.DispersionOO[0]/R6 * AU2KCAL << endl;
            energy = energy + W1.RepulsionOO[0] * exp(-W1.RepulsionOO[1]*R) + damp*W1.DispersionOO[0]/R6;
            //				cout<<"W1.RepulsionOO[0]="<<W1.RepulsionOO[0]<<"  W1.RepulsionOO[1]="<<W1.RepulsionOO[1]<<" W1.DispersionOO[0]="<<W1.DispersionOO[0]<<" damp="<<damp<<endl;
            //cout << "rep-disp: " << energy << endl;
         }

         else if ( ((i == 1) && (j == 1)) ||  ((i == 2) && (j == 2)) ||  ((i == 1) && (j == 2)) || ((i == 2) && (j == 1)))
         {

            energy = energy +W1.RepulsionHH[0] * exp(-W1.RepulsionHH[1]*R);
         }

         else
         {

            energy = energy + W1.RepulsionOH[0] * exp(-W1.RepulsionOH[1]*R);

         }
      }


      //double energy = TTM_A * exp(-TTM_B*R)  +  TTM_C8 / R8  +  TTM_C6 / R6;   
      return energy;
}

/////////////////////////////////////////////////////////////////////
//
//  WaterWaterDispersionGradient -- March 6, 2007 by Tae Hoon Choi
// 
//  
//  

double WaterWaterDispersionGradient(Water &W1, Water &W2,double *Grad1, double *Grad2, int verbose)
{


   double energy = 0.0;
   double de = 0.0;

   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         double dx = W1.AtomCenters[i].x - W2.AtomCenters[j].x; 
         double dy = W1.AtomCenters[i].y - W2.AtomCenters[j].y;
         double dz = W1.AtomCenters[i].z - W2.AtomCenters[j].z;
         double R2 = dx*dx + dy*dy + dz*dz;
         double R = sqrt(R2);
         double delta = W1.DispersionOO[1];
         if ( (i == 0) && (j ==0) )
         {
            double R6 = R2*R2*R2;

            double temp = (1.0 + delta*R + pow(delta*R,2)/2.0 + pow(delta*R,3)/6.0 + pow(delta*R,4)/24.0 + pow(delta*R,5)/120.0 + pow(delta*R,6)/720.0) * exp(-delta*R);
            double damp = 1-temp;
            double damp1 =  -exp(-delta*R)*(delta + pow(delta,2)*R + delta*pow(delta*R,2)/2.0 + delta*pow(delta*R,3)/6.0 + delta*pow(delta*R,4)/24.0 + delta*pow(delta*R,5)/120.0) + temp*delta;

            energy = energy + W1.RepulsionOO[0] * exp(-W1.RepulsionOO[1]*R) + damp*W1.DispersionOO[0]/R6;
            de =  -W1.RepulsionOO[0]*W1.RepulsionOO[1]*exp(-W1.RepulsionOO[1]*R) - 6.0*damp*W1.DispersionOO[0]/(R6*R) + W1.DispersionOO[0]*damp1/R6;
         }

         else if ( ((i == 1) && (j == 1)) ||  ((i == 2) && (j == 2)) ||  ((i == 1) && (j == 2)) || ((i == 2) && (j == 1)))
         {

            energy = energy + W1.RepulsionHH[0] * exp(-W1.RepulsionHH[1]*R);
            de = (-W1.RepulsionHH[1])*W1.RepulsionHH[0] * exp(-W1.RepulsionHH[1]*R);
         }

         else
         {

            energy = energy + W1.RepulsionOH[0] * exp(-W1.RepulsionOH[1]*R);
            de = (-W1.RepulsionOH[1])*W1.RepulsionOH[0] * exp(-W1.RepulsionOH[1]*R);

         }
         Grad1[i*3]   = Grad1[i*3]   + de*dx/R;
         Grad1[i*3+1] = Grad1[i*3+1] + de*dy/R; 
         Grad1[i*3+2] = Grad1[i*3+2] + de*dz/R;
         Grad2[j*3]   = Grad2[j*3]   - de*dx/R;
         Grad2[j*3+1] = Grad2[j*3+1] - de*dy/R; 
         Grad2[j*3+2] = Grad2[j*3+2] - de*dz/R;

      }

      return energy;
}





//////////////////////////////////////////////////////////////
//
//  we have a dipole (mx, my, mz) at (ax, ay, az)
//  and want to know its field (ex, ey, ez) at (bx, by, bz) 
//
//  this uses Stone's tensor T to be able to do the TTM damping 
//
//
void FieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az, 
                   double &bx, double &by, double &bz, double &ex, double &ey, double &ez,double alphaA, double alphaB, double dampDD)
{
   double Rx = bx - ax;
   double Ry = by - ay;
   double Rz = bz - az;
   double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
   double OneOverR5 = pow(R2, -2.5); 

   double damp1 = ElectrostaticDDDamp(sqrt(R2),alphaA,alphaB,dampDD);
   double damp2 = ElectrostaticDCDamp(sqrt(R2),alphaA,alphaB,dampDD);
   //double damp = 1.0;


   double Txx = (damp1*3.*Rx*Rx - damp2*R2) * OneOverR5;
   double Tyy = (damp1*3.*Ry*Ry - damp2*R2) * OneOverR5;
   double Tzz = (damp1*3.*Rz*Rz - damp2*R2) * OneOverR5;

   double Txy = damp1*3.*Rx*Ry * OneOverR5;    double Tyx = Txy;
   double Txz = damp1*3.*Rx*Rz * OneOverR5;    double Tzx = Txz;
   double Tyz = damp1*3.*Ry*Rz * OneOverR5;    double Tzy = Tyz;

  //  printf("Tij(%14.7f)= %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f\n", sqrt(R2), Txx, Tyy, Tzz, Txy, Txz, Tyz);

   ex = Txx*mx + Txy*my + Txz*mz;
   ey = Tyx*mx + Tyy*my + Tyz*mz;
   ez = Tzx*mx + Tzy*my + Tzz*mz;

   //printf("FieldOfDipole: A=(%10.8f,%10.8f,%10.8f)  B=(%10.8f,%10.8f,%10.8f)  R2=%10.8f  E=(%10.8f,%10.8f,%10.8f)\n",
   //        ax, ay, az, bx, by, bz, R2, ex, ey, ez);


}


void dFieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az, 
                    double &bx, double &by, double &bz, double &Txx, double &Tyy, double &Tzz,
                    double &Txy, double &Txz, double &Tyz, double alphaA, double alphaB, double dampDD)
{
   double Rx = bx - ax;
   double Ry = by - ay;
   double Rz = bz - az;
   double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
   double OneOverR5 = pow(R2, -2.5); 

   double damp1 = ElectrostaticDDDamp(sqrt(R2),alphaA,alphaB,dampDD);
   double damp2 = ElectrostaticDCDamp(sqrt(R2),alphaA,alphaB,dampDD);
   //double damp = 1.0;


   Txx = (damp1*3.*Rx*Rx - damp2*R2) * OneOverR5;
   Tyy = (damp1*3.*Ry*Ry - damp2*R2) * OneOverR5;
   Tzz = (damp1*3.*Rz*Rz - damp2*R2) * OneOverR5;

   Txy = damp1*3.*Rx*Ry * OneOverR5;   
   Txz = damp1*3.*Rx*Rz * OneOverR5;   
   Tyz = damp1*3.*Ry*Rz * OneOverR5;   


   //printf("FieldOfDipole: A=(%10.8f,%10.8f,%10.8f)  B=(%10.8f,%10.8f,%10.8f)  R2=%10.8f  E=(%10.8f,%10.8f,%10.8f)\n",
   //        ax, ay, az, bx, by, bz, R2, ex, ey, ez);


}



void dFieldOfDipole(double &mx, double &my, double &mz, double &ax, double &ay, double &az, 
                    double &bx, double &by, double &bz, double &dxex, double &dxey, double &dxez,
                    double &dyex, double &dyey, double &dyez, double &dzex, double &dzey, double &dzez,
                    double alphaA, double alphaB, double dampDD)
{
   double Rx = bx - ax;
   double Ry = by - ay;
   double Rz = bz - az;
   double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
   double distance =sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

   double damp = ElectrostaticDCDamp(sqrt(R2),alphaA,alphaB,dampDD);
   //double damp = 1.0;

   double de0 = -dampDD*3.0*pow(alphaA*alphaB,-0.5)*(1-damp)/(Rx*Rx + Ry*Ry + Rz*Rz);
   double de1 = -damp *pow(Rx*Rx + Ry*Ry + Rz*Rz, -2.5)*(-3.0);

   double de2 = dampDD*dampDD*3.0*3.0*(1-damp)/distance/(alphaA*alphaB);
   double de3 = (3*damp-3.0*dampDD*pow(distance,3)*pow(alphaA*alphaB,-0.5)*(1-damp))*pow((Rx*Rx + Ry*Ry + Rz*Rz),-2.5);
   double de4 = (3*damp-3.0*dampDD*pow(distance,3)*pow(alphaA*alphaB,-0.5)*(1-damp))*(-5.0)*pow((Rx*Rx + Ry*Ry + Rz*Rz),-3.5);


   double dTxx = (de0+de1)*Rx + de3*2.*Rx + (de2+de4)*Rx*Rx*Rx;
   double dTyy = (de0+de1)*Rx +             (de2+de4)*Rx*Ry*Ry;
   double dTzz = (de0+de1)*Rx +             (de2+de4)*Rx*Rz*Rz;
   double dTxy =                de3*Ry +    (de2+de4)*Rx*Rx*Ry; double dTyx = dTxy;
   double dTxz =                de3*Rz +    (de2+de4)*Rx*Rx*Rz; double dTzx = dTxz;
   double dTyz =                            (de2+de4)*Rx*Ry*Rz; double dTzy = dTyz;

   dxex = dTxx*mx + dTxy*my + dTxz*mz;
   dxey = dTyx*mx + dTyy*my + dTyz*mz;
   dxez = dTzx*mx + dTzy*my + dTzz*mz;


   dTxx = (de0+de1)*Ry +             (de2+de4)*Ry*Rx*Rx;
   dTyy = (de0+de1)*Ry + de3*2.*Ry + (de2+de4)*Ry*Ry*Ry;
   dTzz = (de0+de1)*Ry +             (de2+de4)*Ry*Rz*Rz;
   dTxy =                de3*Rx +    (de2+de4)*Ry*Rx*Ry;  dTyx = dTxy;
   dTxz =                            (de2+de4)*Ry*Rx*Rz;  dTzx = dTxz;
   dTyz =                de3*Rz +    (de2+de4)*Ry*Ry*Rz;  dTzy = dTyz;

   dyex = dTxx*mx + dTxy*my + dTxz*mz;
   dyey = dTyx*mx + dTyy*my + dTyz*mz;
   dyez = dTzx*mx + dTzy*my + dTzz*mz;


   dTxx = (de0+de1)*Rz +             (de2+de4)*Rz*Rx*Rx;
   dTyy = (de0+de1)*Rz +             (de2+de4)*Rz*Ry*Ry;
   dTzz = (de0+de1)*Rz + de3*2.*Rz + (de2+de4)*Rz*Rz*Rz;
   dTxy =                            (de2+de4)*Rz*Rx*Ry;  dTyx = dTxy;
   dTxz =                de3*Rx +    (de2+de4)*Rz*Rx*Rz;  dTzx = dTxz;
   dTyz =                de3*Ry +    (de2+de4)*Rz*Ry*Rz;  dTzy = dTyz;

   dzex = dTxx*mx + dTxy*my + dTxz*mz;
   dzey = dTyx*mx + dTyy*my + dTyz*mz;
   dzez = dTzx*mx + dTzy*my + dTzz*mz;





   //printf("FieldOfDipole: A=(%10.8f,%10.8f,%10.8f)  B=(%10.8f,%10.8f,%10.8f)  R2=%10.8f  E=(%10.8f,%10.8f,%10.8f)\n",
   //        ax, ay, az, bx, by, bz, R2, ex, ey, ez);


}



//////////////////////////////////////////////////////////
//
//  Print the induced dipole of all waters in the cluster
//
//
void WaterCluster::ReportInducedDipoles()
{
   printf("Intermolecular induced dipoles in Debye (coordinates of sites in Angs)\n");
   int j, k;
   double x, y, z, mx, my, mz, absdp;
   for (j = 0; j < nwaters; ++j) {
      printf(" Water #%i\n", j);
      for (k = 0; k < waters[j]->nAtomCenters; ++k) {
         double a = (waters[j]->AtomCenters[k]).alpha;
         if (a != 0) {
            mx = (waters[j]->AtomCenters[k]).mx*AU2Debye;
            my = (waters[j]->AtomCenters[k]).my*AU2Debye;
            mz = (waters[j]->AtomCenters[k]).mz*AU2Debye;
            x = (waters[j]->AtomCenters[k]).x*Bohr2Angs;
            y = (waters[j]->AtomCenters[k]).y*Bohr2Angs;
            z = (waters[j]->AtomCenters[k]).z*Bohr2Angs;
            absdp = sqrt(mx*mx + my*my + mz*mz);
            printf("   PolSite #%d  mu = %10.4f   at (%10.4f %10.4f %10.4f)  mu = (%10.8f %10.8f %10.8f)\n",
               k, absdp, x, y, z, mx, my, mz);
         }
      }
   }

}



void InverseMatrix(double* Input, double* Output, int order, int type )
{

   int n=order;
   int m=n;
   double **a = new double*[n];
   double **b = new double*[n];
   for (int ii=0; ii<n; ++ii){
      a[ii]= new double[n];
      b[ii]= new double[m];
   }

   int index;

   if (type==0){
      index=0;
      for (int ii=0; ii< n; ++ii){
         for (int ij=ii; ij< n; ++ij){
            a[ii][ij] = Input[index++];
            if (ii != ij) a[ij][ii]=a[ii][ij];
         }
      }
   }
   else if (type ==1){
      index=0;
      for (int ii=0; ii< n; ++ii){
         for (int ij=0; ij< n; ++ij){
            a[ij][ii] = Input[index++];
         }
      }
   } 




   /* 
   I'm going to find the Inverse Matrix using Gauss-Jordan Elimination method
   this source code came from Numerical recipes: the art of scientific computing
   by William H. Press
   */

   // cout<<"start finding inverse matrix ------------"<<"\n";
   int i, icol=0, irow=0, j, k, l, ll;
   double big,dum,pivinv;
   int *indxc, *indxr, *ipiv;
   indxc = new int[n];
   indxr = new int[n];
   ipiv = new int[n];
   for (j=0;j<n;j++) ipiv[j]=0;
   for (i=0;i<n;i++){
      big=0.0;
      for (j=0;j<n;j++)
         if (ipiv[j] !=1)
            for (k=0;k<n;k++){
               if (ipiv[k]==0){
                  if ( abs(a[j][k]) >= big){
                     big=abs(a[j][k]);
                     irow=j;
                     icol=k;
                  }
               }
            }
            ++(ipiv[icol]);
            if( irow != icol){
               for(l=0;l<n;l++) swap(a[irow][l],a[icol][l]);
               for(l=0;l<m;l++) swap(b[irow][l],b[icol][l]);
            }
            indxr[i]=irow;
            indxc[i]=icol;
            if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
            pivinv=1.0/a[icol][icol];
            a[icol][icol]=1.0;
            for (l=0;l<n;l++) a[icol][l] *= pivinv;
            for (l=0;l<m;l++) b[icol][l] *= pivinv;
            for (ll=0;ll<n;ll++)
               if(ll != icol){
                  dum=a[ll][icol];
                  a[ll][icol]=0.0;
                  for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                  for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
               }
               //	  cout<<"-- first checking  -- \n";
   }
   for(l=n-1; l>=0;l--){
      if (indxr[l] != indxc[l])
         for (k=0;k<n;k++)
            swap(a[k][indxr[l]],a[k][indxc[l]]);
      //	  cout<<"-- second checking "<<l<<" -- \n";
   }

   if (type==0){
      index=0;
      for (int ii=0; ii< n; ++ii){
         for (int ij=ii; ij< n; ++ij){
            Output[index++]=a[ii][ij];
         }
      }
   }
   else if (type==1){
      index=0;
      for (int ii=0; ii< n; ++ii){
         for (int ij=0; ij< n; ++ij){
            Output[index++]=a[ij][ii];
         }
      }
   }

   for (int ii=0; ii<n; ++ii){
      delete [] a[ii];
      delete [] b[ii];
   }

   delete [] a;
   delete [] b;

}


inline void swap(double &a, double &b) {
   double swap = a; a=b; b=swap;
}


void watergeo(double ox, double oy, double oz, double h1x, double h1y, double h1z,
              double h2x, double h2y, double h2z)
{
   printf("  Cartesians\n");
   printf("    O  %12.8f  %12.8f  %12.8f\n", ox, oy, oz);
   printf("    H  %12.8f  %12.8f  %12.8f\n", h1x, h1y, h1z);
   printf("    H  %12.8f  %12.8f  %12.8f\n", h2x, h2y, h2z);

   double r1x = ox - h1x; double r1y = oy - h1y; double r1z = oz - h1z;
   double r2x = ox - h2x; double r2y = oy - h2y; double r2z = oz - h2z;
   double r1 = sqrt(r1x*r1x + r1y*r1y + r1z*r1z);
   double r2 = sqrt(r2x*r2x + r2y*r2y + r2z*r2z);
   double phi = acos((r1x*r2x + r1y*r2y + r1z*r2z)/r1/r2);
   printf("    R1 = %11.8f,  R2 = %11.8f,  phi = %10.7f\n", r1, r2, phi);
}


void WaterCluster::ConvertF (double *Gradient, double *gradEular)
{


   double phi,theta;
   double c1,c2,s1,s2;
   double xterm,yterm,zterm;
   double tau[3],ephi[3],etheta[3],echi[3];
   int i, k ;

   double *Cartesian = new double[nwaters*3*4];
   double *configuration = new double[nwaters*6];
   for(int k=0; k < nwaters*4 ;++k){
      Cartesian[k*3]=0.00; Cartesian[k*3+1]=0.00; Cartesian[k*3+2]=0.00;
      if ( k < nwaters*2) {
         gradEular[k*3]=0.00; gradEular[k*3+1]=0.00; gradEular[k*3+2]=0.00;
      }
   }

   GetXYZM(Cartesian);
   GetRigidConf(configuration);

   for( k=0; k<nwaters;++k){
      phi=configuration[(2*k+1)*3];
      theta=configuration[(2*k+1)*3+1];
      c1=cos(phi), c2=cos(theta);
      s1=sin(phi), s2=sin(theta);
      //
      //  construct the unit vector that will convert the torgue from Cartesian to Eler angle
      //
      ephi[0] = 0.0;
      ephi[1] = 0.0;
      ephi[2] = 1.0;
      etheta[0] = -s1;
      etheta[1]= c1;
      etheta[2]= 0.0;
      echi[0] = s2 * c1;
      echi[1] = s2 * s1;
      echi[2] = c2;
      for (i=0; i<4;++i){
         gradEular[2*k*3] = gradEular[2*k*3]+Gradient[(4*k+i)*3];
         gradEular[2*k*3+1] = gradEular[2*k*3+1]+Gradient[(4*k+i)*3+1];
         gradEular[2*k*3+2] = gradEular[2*k*3+2]+Gradient[(4*k+i)*3+2];
      }
      for (i=0; i<3;++i)
         tau[i]=0.00;


      for (i=0; i<4;++i){
         xterm = Cartesian[(4*k+i)*3]-configuration[2*k*3];
         yterm = Cartesian[(4*k+i)*3+1]-configuration[2*k*3+1];
         zterm = Cartesian[(4*k+i)*3+2]-configuration[2*k*3+2];
         tau[0] = tau[0] + yterm*Gradient[(4*k+i)*3+2] - zterm*Gradient[(4*k+i)*3+1];
         tau[1] = tau[1] + zterm*Gradient[(4*k+i)*3]   - xterm*Gradient[(4*k+i)*3+2];
         tau[2] = tau[2] + xterm*Gradient[(4*k+i)*3+1] - yterm*Gradient[(4*k+i)*3];
      }
      for (i=0; i<3;++i){
         gradEular[(2*k+1)*3] = gradEular[(2*k+1)*3] + tau[i]*ephi[i];
         gradEular[(2*k+1)*3+1] = gradEular[(2*k+1)*3+1] + tau[i]*etheta[i];
         gradEular[(2*k+1)*3+2] = gradEular[(2*k+1)*3+2] + tau[i]*echi[i];
      }
   }

   delete [] Cartesian;

}




void WaterCluster::GetXYZM(double *coordinates)
{       
   for (int k=0; k < nwaters; k++) {
      coordinates[k*12]  =waters[k]->AtomCenters[0].x;
      coordinates[k*12+1]=waters[k]->AtomCenters[0].y;
      coordinates[k*12+2]=waters[k]->AtomCenters[0].z;
      coordinates[k*12+3]=waters[k]->AtomCenters[1].x;
      coordinates[k*12+4]=waters[k]->AtomCenters[1].y;
      coordinates[k*12+5]=waters[k]->AtomCenters[1].z;
      coordinates[k*12+6]=waters[k]->AtomCenters[2].x;
      coordinates[k*12+7]=waters[k]->AtomCenters[2].y;
      coordinates[k*12+8]=waters[k]->AtomCenters[2].z;
      coordinates[k*12+9]=waters[k]->AtomCenters[3].x;
      coordinates[k*12+10]=waters[k]->AtomCenters[3].y;
      coordinates[k*12+11]=waters[k]->AtomCenters[3].z;
   }  
}

void WaterCluster::PutIntoCenterOfMassSystem(int nw, dVec& pos)
{
  const double mH = 1.0;
  const double mO = 16.0;
  double sumx = 0, sumy = 0, sumz = 0, summ = 0;
  for (int iw = 0; iw < nw; ++iw) {
    int iadd = 9*iw;
    sumx += mO * pos[iadd+0] + mH * (pos[iadd+3] + pos[iadd+6]);
    sumy += mO * pos[iadd+1] + mH * (pos[iadd+4] + pos[iadd+7]);
    sumz += mO * pos[iadd+2] + mH * (pos[iadd+5] + pos[iadd+8]);
    summ += mO + 2 * mH;
    }
  double comx = sumx / summ;
  double comy = sumy / summ;
  double comz = sumz / summ;
  for (int iatom = 0; iatom < 3*nw; ++iatom) {
    int iadd = 3*iatom;
    pos[iadd+0] -= comx;
    pos[iadd+1] -= comy;
    pos[iadd+2] -= comz;
  }
}
