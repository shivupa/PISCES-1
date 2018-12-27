#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <iomanip>
#include <fstream>
#include <iostream>

#include "constants.h"
#include "ttm_data.h"

#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"

#define max(x,y) ((x) < (y)? (y): (x))
#define min(x,y) ((x) > (y)? (y): (x))

using namespace std;


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//TS**************************************
//
//  Water-Constructor
//
//  each Water object has a RepulsiveCore object
//  the point charges and LJ-terms are also stored in the RepulsiveCore object
//
Water::Water(double x, double y, double z, int chargeselector, int RepcoreSelector, double divider)
{

  nSites = 4; // DPP


  //TS
  //    this is where each rigid water monomer keeps its center
  //
  x_origin = x;
  y_origin = y;
  z_origin = z;

  phi =0;
  theta = 0;
  psi = 0;

  nAtomCenters = 0;
  AtomCenters = 0;

  //TS
  //    this initializes the water model for the neutral cluster
  //    it adds "BasisFunctionCenters" however, these have no basis functions
  //    it merely uses the same data structures to store charges and polarizabilities
  //    associated with each site
  //    these BasisFunctionCenters are AtomCenters belonging to RepCore     
  //
  SetupCharge(chargeselector);

  //TS
  //    here the basis functions representing the repulsive core are added
  //    the MOs needed for Rossky's potential are added
  //    and Rossky's C and K are computed
  //    the RepulsiveCore.AtomCenters are directly changed, the order of sites from SetupCharge is used
  //
  //  this part is not needed for stand-alone DPP
  //
  // AddRepulsiveCore(RepcoreSelector);
  // this->divider = divider;
  // begin calculating;
  //	Calculate_K();
  //	Calculate_C();
}

Water::~Water()
{
  // each atom center takes care of destoying itself
  if (AtomCenters > 0) 
    delete[] AtomCenters;
}




//TS*************************************************************
//
//  Water::SetupCharge
//
//  calls RepulsiveCore::Setup(nAtoms, nMOs, x, y, z)  which initializes its repulsive core
//  calls RepulsiveCore::AddBasisfynctionCenter
//        these centers have the DC charges and the DC sigma/epsilons
//        they do not have any GTOs on them
//
void Water::SetupCharge(int chargeselector)
{
  switch(chargeselector)
    {
    case 1: 
      // TTM2-u6 models
      Setup(nSites, x_origin, y_origin, z_origin); //TS  sets up its repulsive core
      // cout << "four-site TTM/A water monomer are used.\n";
      RefPosOx = 0;
      RefPosOy = 0;
      RefPosOz = TTM_Oz*Angs2Bohr;
      RefPosHx = 0;
      RefPosHy = TTM_Hy*Angs2Bohr;
      RefPosHz = TTM_Hz*Angs2Bohr;
      RefPosMx = 0;
      RefPosMy = 0;
      RefPosMz = TTM_Mz*Angs2Bohr;
      //TS then 4 basis function centers are added to the repulsive core
      //   where x, y, and z are relative to x_origin from above
      //   the order is important since the Drudes are put on Atom3 (the forth atom)
      //
      //   NOTE:  the sigma for each atom center is used as polarizability for electrostatic damping!!!!
      //                                                                      
      //                  index  x, y,               z,       nGTOs, charge, sigma,  eps,   alpha                
      AddBasisfunctionCenter(0, RefPosOx,  RefPosOy, RefPosOz, 0,   0,    0.0,       0.0, TTM_POL_O);  // O
      AddBasisfunctionCenter(1, RefPosHx,  RefPosHy, RefPosHz, 0, TTM_QH, TTM_POL_H, 0.0, TTM_POL_H);  // H1
      AddBasisfunctionCenter(2, RefPosHx, -RefPosHy, RefPosHz, 0, TTM_QH, TTM_POL_H, 0.0, TTM_POL_H);  // H2
      AddBasisfunctionCenter(3, RefPosMx,  RefPosMy, RefPosMz, 0, TTM_QM, TTM_POL_O, 0.0, 0.0);        // M
      // the positions of the polarizabilities are set later, after reorientation
      // again the order of the centers is used: POL1 = Pos[3]  POL2 = Pos[1]  POL3 = Pos[2]
	
		
      //Setting the neutral model parameters
      DispersionOO = new double[2];
      DispersionOO[0] = u6_c6; DispersionOO[1] = u6_delta;
      //DispersionOH = NULL;
      //DispersionHH = NULL;
      
      RepulsionOO = new double[2];
      RepulsionOO[0] = u6_A1; RepulsionOO[1] = u6_b1;
      RepulsionOH = new double[2];
      RepulsionOH[0] = u6_A3; RepulsionOH[1] = u6_b3;
      RepulsionHH = new double[2];
      RepulsionHH[0] = u6_A2; RepulsionHH[1] = u6_b2;
      
      Damping = new double[3];
      Damping[0] = 0.; Damping[1] = 0.3; Damping[2] = 0.3;
      
      PolSites = TTM_NoOfPolSites;
      
      
      break;
    case 11: 
      // KT charges, TTM polarizability and neutral parameters
      Setup(nSites, x_origin, y_origin, z_origin);
      RefPosOx = 0;
      RefPosOy = 0;
      RefPosOz = TTM_Oz*Angs2Bohr;
      RefPosHx = 0;
      RefPosHy = TTM_Hy*Angs2Bohr;
      RefPosHz = TTM_Hz*Angs2Bohr;
      RefPosMx = 0;
      RefPosMy = 0;
      RefPosMz = DC_Mz*Angs2Bohr;
      cout << "four-site TTM2-u6 water monomer with KT charges and DC M-site.\n";
      {
      double KT_pol_O = TTM_POL_O / W_ALPHA_EXP * W_ALPHA_HF;  // 1st guess
      double KT_pol_H = TTM_POL_H / W_ALPHA_EXP * W_ALPHA_HF;
      //
      //                  index  x, y,               z,      nGTOs, charge,  sigma, eps                
      AddBasisfunctionCenter(0, RefPosOx,  RefPosOy, RefPosOz, 0,   0,    0.0,     0.0, KT_pol_O);  // O
      AddBasisfunctionCenter(1, RefPosHx,  RefPosHy, RefPosHz, 0, KT_QH, KT_pol_H, 0.0, KT_pol_H);  // H1
      AddBasisfunctionCenter(2, RefPosHx, -RefPosHy, RefPosHz, 0, KT_QH, KT_pol_H, 0.0, KT_pol_H);  // H2
      AddBasisfunctionCenter(3, RefPosMx,  RefPosMy, RefPosMz, 0, KT_QM, KT_pol_O, 0.0, 0.0);       // M
      }
		  
		  
		  
      //Setting the neutral model parameters
      DispersionOO = new double[2];
      DispersionOO[0] = u6_c6; DispersionOO[1] = u6_delta;
      //DispersionOH = NULL;
      //DispersionHH = NULL;
      
      RepulsionOO = new double[2];
      RepulsionOO[0] = u6_A1; RepulsionOO[1] = u6_b1;
      RepulsionOH = new double[2];
      RepulsionOH[0] = u6_A3; RepulsionOH[1] = u6_b3;
      RepulsionHH = new double[2];
      RepulsionHH[0] = u6_A2; RepulsionHH[1] = u6_b2;
      
      Damping = new double[3];
      Damping[0] = 0.; Damping[1] = 0.3; Damping[2] = 0.3;
      
      PolSites = TTM_NoOfPolSites;
      
      
      break;
      
    case 2: 
      // DPP model
      Setup(nSites, x_origin, y_origin, z_origin); //TS  sets up its repulsive core
      
      RefPosOx = 0;
      RefPosOy = 0;
      RefPosOz = TTM_Oz*Angs2Bohr;
      RefPosHx = 0;
      RefPosHy = TTM_Hy*Angs2Bohr;
      RefPosHz = TTM_Hz*Angs2Bohr;
      RefPosMx = 0;
      RefPosMy = 0;
      RefPosMz = TTM_Mz*Angs2Bohr;
      //TS then 4 basis function centers are added to the repulsive core
      //   where x, y, and z are relative to x_origin from above
      //   the order is important since the Drudes are put on Atom3 (the forth atom)
      //
      //   NOTE:  the sigma for each atom center is used as polarizability for electrostatic damping!!!!
      //                                                                      
      //                  index  x, y,               z,       nGTOs, charge, sigma,  eps,   alpha                
      AddBasisfunctionCenter(0, RefPosOx,  RefPosOy, RefPosOz, 0,   0,    0.0,       0.0, TTM_POL_O);  // O
      AddBasisfunctionCenter(1, RefPosHx,  RefPosHy, RefPosHz, 0, TTM_QH, TTM_POL_H, 0.0, TTM_POL_H);  // H1
      AddBasisfunctionCenter(2, RefPosHx, -RefPosHy, RefPosHz, 0, TTM_QH, TTM_POL_H, 0.0, TTM_POL_H);  // H2
      AddBasisfunctionCenter(3, RefPosMx,  RefPosMy, RefPosMz, 0, TTM_QM, TTM_POL_O, 0.0, 0.0);        // M
      // the positions of the polarizabilities are set later, after reorientation
      // again the order of the centers is used: POL1 = Pos[3]  POL2 = Pos[1]  POL3 = Pos[2]
      
      
      //Setting the neutral model parameters
      DispersionOO = new double[2];
      DispersionOO[0] = u6_c6; DispersionOO[1] = u6_delta;
      //DispersionOH = NULL;
      //DispersionHH = NULL;
      
      RepulsionOO = new double[2];
      RepulsionOO[0] = DPP_A1; RepulsionOO[1] = DPP_b1;
      RepulsionOH = new double[2];
      RepulsionOH[0] = DPP_A3; RepulsionOH[1] = DPP_b3;
      RepulsionHH = new double[2];
      RepulsionHH[0] = DPP_A2; RepulsionHH[1] = DPP_b2;
      
      Damping = new double[3];
   //   Damping[0] = 0.; Damping[1] = 0.23; Damping[2] = 0.3;
      Damping[0] = 0.; Damping[1] = 0.21; Damping[2] = 0.3;
      
      PolSites = TTM_NoOfPolSites;
      break;
      
      
    case 12: 
      // KT charges, TTM polarizability, DPP neutral parameters
      Setup(nSites, x_origin, y_origin, z_origin);
      RefPosOx = 0;
      RefPosOy = 0;
      RefPosOz = TTM_Oz*Angs2Bohr;
      RefPosHx = 0;
      RefPosHy = TTM_Hy*Angs2Bohr;
      RefPosHz = TTM_Hz*Angs2Bohr;
      RefPosMx = 0;
      RefPosMy = 0;
      RefPosMz = DC_Mz*Angs2Bohr;
      cout << "four-site DPP water monomer with KT charges and DC M-site.\n";
      {
	double KT_pol_O = TTM_POL_O / W_ALPHA_EXP * W_ALPHA_HF;  // 1st guess
	double KT_pol_H = TTM_POL_H / W_ALPHA_EXP * W_ALPHA_HF;
	//
	//                  index  x, y,               z,      nGTOs, charge,  sigma, eps                
	AddBasisfunctionCenter(0, RefPosOx,  RefPosOy, RefPosOz, 0,   0,    0.0,     0.0, KT_pol_O);  // O
	AddBasisfunctionCenter(1, RefPosHx,  RefPosHy, RefPosHz, 0, KT_QH, KT_pol_H, 0.0, KT_pol_H);  // H1
	AddBasisfunctionCenter(2, RefPosHx, -RefPosHy, RefPosHz, 0, KT_QH, KT_pol_H, 0.0, KT_pol_H);  // H2
	AddBasisfunctionCenter(3, RefPosMx,  RefPosMy, RefPosMz, 0, KT_QM, KT_pol_O, 0.0, 0.0);       // M
      }
      
      //Setting the neutral model parameters
      DispersionOO = new double[2];
      DispersionOO[0] = u6_c6; DispersionOO[1] = u6_delta;
      //DispersionOH = NULL;
      //DispersionHH = NULL;
      
      RepulsionOO = new double[2];
      RepulsionOO[0] = DPP_A1; RepulsionOO[1] = DPP_b1;
      RepulsionOH = new double[2];
      RepulsionOH[0] = DPP_A3; RepulsionOH[1] = DPP_b3;
      RepulsionHH = new double[2];
      RepulsionHH[0] = DPP_A2; RepulsionHH[1] = DPP_b2;
      
      Damping = new double[3];
   //   Damping[0] = 0.; Damping[1] = 0.23; Damping[2] = 0.3;
       Damping[0] = 0.; Damping[1] = 0.21; Damping[2] = 0.3;
      
      PolSites = TTM_NoOfPolSites;
      break;
      
      
    default:
      cout << "Error in Water::SetupCharge(): unknown chargeselector";
      exit(1);
    }
}



void Water::Setup(int nAtoms, double x, double y, double z)
{
  x_origin = x;
  y_origin = y;
  z_origin = z;
  nAtomCenters = nAtoms;
  AtomCenters = new AtomCenter[nAtomCenters];
}

///////////////////
//
//  nGTOs is ignored
//
//
//
void Water::AddBasisfunctionCenter(int AtomIndex, double x, double y, double z, int nGTOs,
				   double charge, double sigma, double epsilon, double alpha)
{
  if(AtomIndex > this->nAtomCenters)
    {printf("Added more basis function center than reserved"); exit(1);}
  AtomCenters[AtomIndex].x = x + x_origin;
  AtomCenters[AtomIndex].y = y + y_origin;
  AtomCenters[AtomIndex].z = z + z_origin;
  AtomCenters[AtomIndex].charge = charge;
  AtomCenters[AtomIndex].identifier = AtomIndex;
  AtomCenters[AtomIndex].sigma = sigma;  AtomCenters[AtomIndex].epsilon = epsilon;
  AtomCenters[AtomIndex].alpha   = alpha;
}

////void Water::AddRepulsiveCore(int RepcoreSelector)
////{
////  switch(RepcoreSelector)
////    {
////    case 1:
////      //cout << "-CC-PVTZ -F core is selected for A symmetry occupied orbitals\n";
////      AddRepulsiveCCPVTZ();
////      break;
////    case 2:
////      //cout << "A fast 6S core is selected\n";
////      AddRepulsive_6S();
////      break;
////    case 3:
////      //cout << "CC-PVTZ -F core is selected for all of the occupied orbitals\n";
////      AddRepulsiveCCPVTZ_all();
////      break;
////    default:
////      printf(" Water::AddRepulsiveCore: Unrecognized core selection: %i\n", RepcoreSelector);
////      exit(42);
////    }
////}
////
////void Water::AddRepulsive_6S()
////{
////
////	if (AtomCenters[1].GTOs!=NULL)	 delete[] AtomCenters[1].GTOs;
////	AtomCenters[1].nGTOs = 3;
////	AtomCenters[1].GTOs = new GTO[3];
////	AtomCenters[1].dxRGTOs = new GTO[3];
////	AtomCenters[1].dxRGTOs_1 = new GTO[3];
////	AtomCenters[1].dyRGTOs = new GTO[3];
////	AtomCenters[1].dyRGTOs_1 = new GTO[3];
////	AtomCenters[1].dzRGTOs = new GTO[3];
////	AtomCenters[1].dzRGTOs_1 = new GTO[3];
////
////	if (AtomCenters[2].GTOs!=NULL)	 delete[] AtomCenters[2].GTOs;
////	AtomCenters[2].nGTOs = 3;
////	AtomCenters[2].GTOs = new GTO[3];
////	AtomCenters[2].dxRGTOs = new GTO[3];
////	AtomCenters[2].dxRGTOs_1 = new GTO[3];
////	AtomCenters[2].dyRGTOs = new GTO[3];
////	AtomCenters[2].dyRGTOs_1 = new GTO[3];
////	AtomCenters[2].dzRGTOs = new GTO[3];
////	AtomCenters[2].dzRGTOs_1 = new GTO[3];
////
////	this->nGTOs = 3+3;
////
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(1,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(1,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(1,2,'s',1.159);
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(2,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(2,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(2,2,'s',1.159);
////
////	// finish setting up
////	FinishAddingAtoms();
////	double coefficients[3][6]=	//enhance 10 times
////	{{ -2.91026974E-05, 8.09061918E-04, 1.36181551E-03,
////	   -2.91026974E-05, 8.09061918E-04,  1.36181551E-03 } 
////
////	 ,{ 2.59087634E-02, 1.69713119E-01, 7.77840383E-01,
////	    2.59087634E-02, 1.69713119E-01, 7.77840383E-01}
////
////	 ,{ -2.33292227E-02, -1.59884747E-01,-7.47223050E-01,
////	    -2.33292227E-02,-1.59884747E-01,-7.47223050E-01}};
////
////
////	//TS just copy the coefficients into the RepulsiveCore object
////	BuildMO(0,-20.55554,coefficients[0]);
////	BuildMO(1,-1.34646,coefficients[1]);
////	BuildMO(2,-0.57847,coefficients[2]);
////
////	Calculate_K();
////	Calculate_C();
////}
////
////
////void Water::AddRepulsiveCCPVTZ()
////{
////	if (AtomCenters[0].GTOs!=NULL)		delete[] AtomCenters[0].GTOs;
////	if (AtomCenters[0].dxRGTOs!=NULL)	delete[] AtomCenters[0].dxRGTOs;
////	if (AtomCenters[0].dxRGTOs_1!=NULL)	delete[] AtomCenters[0].dxRGTOs_1;
////	if (AtomCenters[0].dyRGTOs!=NULL)	delete[] AtomCenters[0].dyRGTOs;
////	if (AtomCenters[0].dyRGTOs_1!=NULL)	delete[] AtomCenters[0].dyRGTOs_1;
////	if (AtomCenters[0].dzRGTOs!=NULL)	delete[] AtomCenters[0].dzRGTOs;
////	if (AtomCenters[0].dzRGTOs_1!=NULL)	delete[] AtomCenters[0].dzRGTOs_1;
////	AtomCenters[0].nGTOs = 37;
////	AtomCenters[0].GTOs = new GTO[37];
////	AtomCenters[0].dxRGTOs = new GTO[37];
////	AtomCenters[0].dxRGTOs_1 = new GTO[37];
////	AtomCenters[0].dyRGTOs = new GTO[37];
////	AtomCenters[0].dyRGTOs_1 = new GTO[37];
////	AtomCenters[0].dzRGTOs = new GTO[37];
////	AtomCenters[0].dzRGTOs_1 = new GTO[37];
////
////	if (AtomCenters[1].GTOs!=NULL)		delete[] AtomCenters[1].GTOs;
////	if (AtomCenters[1].dxRGTOs!=NULL)	delete[] AtomCenters[1].dxRGTOs;
////	if (AtomCenters[1].dxRGTOs_1!=NULL)	delete[] AtomCenters[1].dxRGTOs_1;
////	if (AtomCenters[1].dyRGTOs!=NULL)	delete[] AtomCenters[1].dyRGTOs;
////	if (AtomCenters[1].dyRGTOs_1!=NULL)	delete[] AtomCenters[1].dyRGTOs_1;
////	if (AtomCenters[1].dzRGTOs!=NULL)	delete[] AtomCenters[1].dzRGTOs;
////	if (AtomCenters[1].dzRGTOs_1!=NULL)	delete[] AtomCenters[1].dzRGTOs_1;
////	AtomCenters[0].nGTOs = 37;
////	AtomCenters[1].nGTOs = 11;
////	AtomCenters[1].GTOs = new GTO[11];
////	AtomCenters[1].dxRGTOs = new GTO[11];
////	AtomCenters[1].dxRGTOs_1 = new GTO[11];
////	AtomCenters[1].dyRGTOs = new GTO[11];
////	AtomCenters[1].dyRGTOs_1 = new GTO[11];
////	AtomCenters[1].dzRGTOs = new GTO[11];
////	AtomCenters[1].dzRGTOs_1 = new GTO[11];
////
////	if (AtomCenters[2].GTOs!=NULL)		delete[] AtomCenters[2].GTOs;
////	if (AtomCenters[2].dxRGTOs!=NULL)	delete[] AtomCenters[2].dxRGTOs;
////	if (AtomCenters[2].dxRGTOs_1!=NULL)	delete[] AtomCenters[2].dxRGTOs_1;
////	if (AtomCenters[2].dyRGTOs!=NULL)	delete[] AtomCenters[2].dyRGTOs;
////	if (AtomCenters[2].dyRGTOs_1!=NULL)	delete[] AtomCenters[2].dyRGTOs_1;
////	if (AtomCenters[2].dzRGTOs!=NULL)	delete[] AtomCenters[2].dzRGTOs;
////	if (AtomCenters[2].dzRGTOs_1!=NULL)	delete[] AtomCenters[2].dzRGTOs_1;
////	AtomCenters[2].nGTOs = 11;
////	AtomCenters[2].GTOs = new GTO[11];
////	AtomCenters[2].dxRGTOs = new GTO[11];
////	AtomCenters[2].dxRGTOs_1 = new GTO[11];
////	AtomCenters[2].dyRGTOs = new GTO[11];
////	AtomCenters[2].dyRGTOs_1 = new GTO[11];
////	AtomCenters[2].dzRGTOs = new GTO[11];
////	AtomCenters[2].dzRGTOs_1 = new GTO[11];
////
////	this->nGTOs = 37+11+11;
////
////		// add base function to O
//////	cout << "RHF/CC-PVTZ core was positioned at "<<x<<','<<y<<','<<z<<" and scaled by "<< divider<<"\n";
////
////	AddBasisfuntionsToBasisfunctionCenter(0,0,'s',15330.0);
////	AddBasisfuntionsToBasisfunctionCenter(0,1,'s',2299.0);
////	AddBasisfuntionsToBasisfunctionCenter(0,2,'s',522.4);
////	AddBasisfuntionsToBasisfunctionCenter(0,3,'s',147.3);
////	AddBasisfuntionsToBasisfunctionCenter(0,4,'s',47.55);
////	AddBasisfuntionsToBasisfunctionCenter(0,5,'s',16.76);
////	AddBasisfuntionsToBasisfunctionCenter(0,6,'s',6.207);
////	AddBasisfuntionsToBasisfunctionCenter(0,7,'s',1.752);
////	AddBasisfuntionsToBasisfunctionCenter(0,8,'s',0.6882);
////	AddBasisfuntionsToBasisfunctionCenter(0,9,'s',0.2384);
////	AddBasisfuntionsToBasisfunctionCenter(0,10,'p',34.46);
////	AddBasisfuntionsToBasisfunctionCenter(0,13,'p',7.749);
////	AddBasisfuntionsToBasisfunctionCenter(0,16,'p',2.280);
////	AddBasisfuntionsToBasisfunctionCenter(0,19,'p',0.7156);
////	AddBasisfuntionsToBasisfunctionCenter(0,22,'p',0.214);
////	AddBasisfuntionsToBasisfunctionCenter(0,25,'d',2.314);
////	AddBasisfuntionsToBasisfunctionCenter(0,31,'d',0.645);
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(1,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(1,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(1,2,'s',1.159);
////	AddBasisfuntionsToBasisfunctionCenter(1,3,'s',0.3258);
////	AddBasisfuntionsToBasisfunctionCenter(1,4,'s',0.1027);
////	AddBasisfuntionsToBasisfunctionCenter(1,5,'p',1.407);
////	AddBasisfuntionsToBasisfunctionCenter(1,8,'p',0.388);
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(2,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(2,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(2,2,'s',1.159);
////	AddBasisfuntionsToBasisfunctionCenter(2,3,'s',0.3258);
////	AddBasisfuntionsToBasisfunctionCenter(2,4,'s',0.1027);
////	AddBasisfuntionsToBasisfunctionCenter(2,5,'p',1.407);
////	AddBasisfuntionsToBasisfunctionCenter(2,8,'p',0.388);
////
////	// finish setting up
////	FinishAddingAtoms();
////	double coefficients[3][59]=
////	{{ 5.08010194E-04, 3.92865633E-03, 2.02382748E-02, 7.91746700E-02, 2.30647254E-01,
////  4.32886147E-01, 3.51313657E-01, 4.74025537E-02,-1.90929741E-02,-4.36816501E-03,
////  0.00000000E+00, 0.00000000E+00,-2.01635343E-04, 0.00000000E+00, 0.00000000E+00,
//// -1.15176683E-03, 0.00000000E+00, 0.00000000E+00,-3.10631500E-04, 0.00000000E+00,
////  0.00000000E+00, 9.30256433E-06, 0.00000000E+00, 0.00000000E+00, 2.76626079E-04,
////  2.12616028E-03, 2.08177312E-03, 2.12069260E-03, 0.00000000E+00, 0.00000000E+00,
////  0.00000000E+00, 4.86430330E-03, 4.53265037E-03, 4.63276657E-03, 0.00000000E+00,
////  0.00000000E+00, 0.00000000E+00,-2.91026974E-06, 8.09061918E-05, 1.36181551E-04,
////  1.49758188E-04, 3.12749475E-04, 0.00000000E+00,-1.56669577E-04, 1.23255784E-04,
////  0.00000000E+00, 4.35544023E-05,-7.26832654E-05,-2.91026974E-06, 8.09061918E-05,
////  1.36181551E-04, 1.49758188E-04, 3.12749475E-04, 0.00000000E+00, 1.56669577E-04,
////  1.23255784E-04, 0.00000000E+00,-4.35544023E-05,-7.26832654E-05} 
////
////	,{-1.03302931E-04,-8.05537730E-04,-4.14714298E-03,-1.68839982E-02,-5.20574306E-02,
////	-1.23361302E-01, -1.49727314E-01, 1.30616813E-01, 5.11453976E-01, 2.88648129E-01,
////	0.00000000E+00, 0.00000000E+00,-1.93820617E-03, 0.00000000E+00, 0.00000000E+00,
////	-1.09939727E-02, 0.00000000E+00, 0.00000000E+00,-3.66971125E-02, 0.00000000E+00,
////	0.00000000E+00, -6.55733902E-02, 0.00000000E+00, 0.00000000E+00,-1.45687845E-02,
////	6.03155225E-04, 8.08209221E-03, 5.94873719E-03, 0.00000000E+00, 0.00000000E+00,
////	0.00000000E+00, -3.82242440E-03, 8.22967004E-03, 7.01680112E-03, 0.00000000E+00,
////	0.00000000E+00, 0.00000000E+00, 2.59087634E-03, 1.69713119E-02, 7.77840383E-02,
////	6.73694026E-02, 5.46307802E-03, 0.00000000E+00,-1.16264283E-02, 8.53441289E-03,
////	0.00000000E+00, -5.96866122E-03, 2.05934811E-03, 2.59087634E-03, 1.69713119E-02,
////	7.77840383E-02, 6.73694026E-02, 5.46307802E-03, 0.00000000E+00, 1.16264283E-02,
////	8.53441289E-03, 0.00000000E+00, 5.96866122E-03, 2.05934811E-03}
////
////	,{ -3.41452189E-05,-2.64616296E-04,-1.37563258E-03,
//// -5.53724398E-03,-1.74541078E-02,-4.05484590E-02,-5.55960357E-02, 4.48196141E-02,
////  1.92494902E-01, 2.98804327E-01, 0.00000000E+00, 0.00000000E+00, 1.30480468E-02,
////  0.00000000E+00, 0.00000000E+00, 8.07162087E-02, 0.00000000E+00, 0.00000000E+00,
////  2.50056552E-01, 0.00000000E+00, 0.00000000E+00, 3.90593708E-01, 0.00000000E+00,
////  0.00000000E+00, 2.80136401E-01, 2.44927119E-03,-3.75890957E-03,-8.60496162E-03,
////  0.00000000E+00, 0.00000000E+00, 0.00000000E+00, 3.63554726E-04,-6.82699669E-03,
//// -3.71690916E-02, 0.00000000E+00, 0.00000000E+00, 0.00000000E+00,-2.33292227E-03,
//// -1.59884747E-02,-7.47223050E-02,-1.29267307E-01,-2.97368131E-02, 0.00000000E+00,
////  9.52655473E-03, 1.55291787E-05, 0.00000000E+00, 1.43355100E-02, 1.51810310E-02,
//// -2.33292227E-03,-1.59884747E-02,-7.47223050E-02,-1.29267307E-01,-2.97368131E-02,
////  0.00000000E+00,-9.52655473E-03, 1.55291787E-05, 0.00000000E+00,-1.43355100E-02,
////  1.51810310E-02}};
////
////	BuildMO(0,-20.55554,coefficients[0]);
////	BuildMO(1,-1.34646,coefficients[1]);
////	BuildMO(2,-0.57847,coefficients[2]);
////
////	Calculate_K();
////	Calculate_C();
////}
////
////void Water::AddRepulsiveCCPVTZ_all()
////{
////
////
////
////	if (AtomCenters[0].GTOs!=NULL)		delete[] AtomCenters[0].GTOs;
////	AtomCenters[0].nGTOs = 37;
////	AtomCenters[0].GTOs = new GTO[37];
////	AtomCenters[0].dxRGTOs = new GTO[37];
////	AtomCenters[0].dxRGTOs_1 = new GTO[37];
////	AtomCenters[0].dyRGTOs = new GTO[37];
////	AtomCenters[0].dyRGTOs_1 = new GTO[37];
////	AtomCenters[0].dzRGTOs = new GTO[37];
////	AtomCenters[0].dzRGTOs_1 = new GTO[37];
////
////	if (AtomCenters[1].GTOs!=NULL)		delete[] AtomCenters[1].GTOs;
////	AtomCenters[1].nGTOs = 11;
////	AtomCenters[1].GTOs = new GTO[11];
////	AtomCenters[1].dxRGTOs = new GTO[11];
////	AtomCenters[1].dxRGTOs_1 = new GTO[11];
////	AtomCenters[1].dyRGTOs = new GTO[11];
////	AtomCenters[1].dyRGTOs_1 = new GTO[11];
////	AtomCenters[1].dzRGTOs = new GTO[11];
////	AtomCenters[1].dzRGTOs_1 = new GTO[11];
////
////	if (AtomCenters[2].GTOs!=NULL)		delete[] AtomCenters[2].GTOs;
////	AtomCenters[2].nGTOs = 11;
////	AtomCenters[2].GTOs = new GTO[11];
////	AtomCenters[2].dxRGTOs = new GTO[11];
////	AtomCenters[2].dxRGTOs_1 = new GTO[11];
////	AtomCenters[2].dyRGTOs = new GTO[11];
////	AtomCenters[2].dyRGTOs_1 = new GTO[11];
////	AtomCenters[2].dzRGTOs = new GTO[11];
////	AtomCenters[2].dzRGTOs_1 = new GTO[11];
////
////	this->nGTOs = 37+11+11;
////
////		// add base function to O
//////	cout << "RHF/CC-PVTZ core was positioned at "<<x<<','<<y<<','<<z<<" and scaled by "<< divider<<"\n";
////
////	AddBasisfuntionsToBasisfunctionCenter(0,0,'s',15330.0);
////	AddBasisfuntionsToBasisfunctionCenter(0,1,'s',2299.0);
////	AddBasisfuntionsToBasisfunctionCenter(0,2,'s',522.4);
////	AddBasisfuntionsToBasisfunctionCenter(0,3,'s',147.3);
////	AddBasisfuntionsToBasisfunctionCenter(0,4,'s',47.55);
////	AddBasisfuntionsToBasisfunctionCenter(0,5,'s',16.76);
////	AddBasisfuntionsToBasisfunctionCenter(0,6,'s',6.207);
////	AddBasisfuntionsToBasisfunctionCenter(0,7,'s',1.752);
////	AddBasisfuntionsToBasisfunctionCenter(0,8,'s',0.6882);
////	AddBasisfuntionsToBasisfunctionCenter(0,9,'s',0.2384);
////	AddBasisfuntionsToBasisfunctionCenter(0,10,'p',34.46);
////	AddBasisfuntionsToBasisfunctionCenter(0,13,'p',7.749);
////	AddBasisfuntionsToBasisfunctionCenter(0,16,'p',2.280);
////	AddBasisfuntionsToBasisfunctionCenter(0,19,'p',0.7156);
////	AddBasisfuntionsToBasisfunctionCenter(0,22,'p',0.214);
////	AddBasisfuntionsToBasisfunctionCenter(0,25,'d',2.314);
////	AddBasisfuntionsToBasisfunctionCenter(0,31,'d',0.645);
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(1,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(1,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(1,2,'s',1.159);
////	AddBasisfuntionsToBasisfunctionCenter(1,3,'s',0.3258);
////	AddBasisfuntionsToBasisfunctionCenter(1,4,'s',0.1027);
////	AddBasisfuntionsToBasisfunctionCenter(1,5,'p',1.407);
////	AddBasisfuntionsToBasisfunctionCenter(1,8,'p',0.388);
////	// add base function to H
////	AddBasisfuntionsToBasisfunctionCenter(2,0,'s',33.87);
////	AddBasisfuntionsToBasisfunctionCenter(2,1,'s',5.095);
////	AddBasisfuntionsToBasisfunctionCenter(2,2,'s',1.159);
////	AddBasisfuntionsToBasisfunctionCenter(2,3,'s',0.3258);
////	AddBasisfuntionsToBasisfunctionCenter(2,4,'s',0.1027);
////	AddBasisfuntionsToBasisfunctionCenter(2,5,'p',1.407);
////	AddBasisfuntionsToBasisfunctionCenter(2,8,'p',0.388);
////
////	// finish setting up
////	FinishAddingAtoms();
////	double coefficients[3][59]=
////  
////   {{5.08010194E-04,  3.92865633E-03,  2.02382747E-02,  7.91746700E-02,  2.30647254E-01,
////    4.32886147E-01,  3.51313657E-01,  4.74025552E-02,  1.90929750E-02, -4.36816663E-03,
////    0.00000000E+00,  0.00000000E+00, -2.01635343E-04,  0.00000000E+00,  0.00000000E+00,
////   -1.15176682E-03,  0.00000000E+00,  0.00000000E+00, -3.10631954E-04,  0.00000000E+00,
////    0.00000000E+00,  9.30289139E-06,  0.00000000E+00,  0.00000000E+00,  2.76626613E-04,
////    2.12616021E-03,  2.08177323E-03,  2.12069279E-03,  0.00000000E+00,  0.00000000E+00,
////    0.00000000E+00,  4.86430378E-03,  4.53265069E-03,  4.63276662E-03,  0.00000000E+00,
////    0.00000000E+00,  0.00000000E+00, -2.91025907E-06,  8.09063978E-05,  1.36181702E-04,
////    1.49758785E-04,  3.12749562E-04,  0.00000000E+00,  1.56669625E-04,  1.23255903E-04,
////    0.00000000E+00, -4.35542402E-05, -7.26831619E-05, -2.91025907E-06,  8.09063978E-05,
////    1.36181702E-04,  1.49758785E-04,  3.12749562E-04,  0.00000000E+00, -1.56669625E-04,
////    1.23255903E-04, 0.00000000E+00,  4.35542402E-05,-7.26831619E-05}
////
////
////   ,{ -1.03302923E-04,
////     -8.05537668E-04,  -4.14714265E-03,  -1.68839969E-02,  -5.20574262E-02,  -1.23361293E-01,
////     -1.49727296E-01, 1.30616794E-01,  -5.11453896E-01, 2.88648015E-01,  -0.00000000E+00,
////     0.00000000E+00,  -1.93820967E-03,  0.00000000E+00,  0.00000000E+00,  -1.09939928E-02,
////     0.00000000E+00,  0.00000000E+00,  -3.66971789E-02,  0.00000000E+00,  0.00000000E+00,
////     -6.55734739E-02,  0.00000000E+00,  0.00000000E+00,  -1.45687930E-02, 6.03154543E-04,
////    8.08209194E-03, 5.94874861E-03,  0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
////     -3.82242858E-03, 8.22964875E-03, 7.01681403E-03,  0.00000000E+00,  0.00000000E+00,
////     0.00000000E+00, 2.59087841E-03, 1.69713282E-02, 7.77840912E-02, 6.73694331E-02,
////    5.46306557E-03,  0.00000000E+00, 1.16264323E-02, 8.53442324E-03,  0.00000000E+00,
////    5.96866117E-03, 2.05935108E-03, 2.59087841E-03, 1.69713282E-02, 7.77840912E-02,
////    6.73694331E-02, 5.46306557E-03,  0.00000000E+00,  -1.16264323E-02, 8.53442324E-03,
////     0.00000000E+00,  -5.96866117E-03, 2.05935108E-03}
////
////  
////   ,{ -3.41452433E-05, -2.64616485E-04, -1.37563356E-03,
////    -5.53724793E-03, -1.74541202E-02, -4.05484875E-02, -5.55960742E-02,  4.48196372E-02,
////    -1.92495049E-01,  2.98804472E-01,  0.00000000E+00,  0.00000000E+00,  1.30480462E-02,
////     0.00000000E+00,  0.00000000E+00,  8.07162043E-02,  0.00000000E+00,  0.00000000E+00,
////     2.50056540E-01,  0.00000000E+00,  0.00000000E+00,  3.90593676E-01,  0.00000000E+00,
////     0.00000000E+00,  2.80136275E-01,  2.44927190E-03, -3.75890283E-03, -8.60497391E-03,
////     0.00000000E+00,  0.00000000E+00,  0.00000000E+00,  3.63554651E-04, -6.82699166E-03,
////     -3.71691082E-02,  0.00000000E+00,  0.00000000E+00,  0.00000000E+00, -2.33292289E-03,
////    -1.59884829E-02, -7.47223160E-02, -1.29267311E-01, -2.97367232E-02,  0.00000000E+00,
////    -9.52655917E-03,  1.55256296E-05,  0.00000000E+00, -1.43355105E-02,  1.51810370E-02,
////    -2.33292289E-03, -1.59884829E-02, -7.47223160E-02, -1.29267311E-01, -2.97367232E-02,
////     0.00000000E+00,  9.52655917E-03,  1.55256296E-05,  0.00000000E+00,  1.43355105E-02,
////     1.51810370E-02}};
////
////
////	BuildMO(0,-20.55554,coefficients[0]);
////	BuildMO(1,-1.34646,coefficients[1]);
////	BuildMO(2,-0.57847,coefficients[2]);
////
////	Calculate_K();
////	Calculate_C();
////}
////

void Water::Reposition_A(double o_x, double o_y, double o_z, 
			 double h1_x, double h1_y, double h1_z, 
			 double h2_x, double h2_y, double h2_z) 
{
  Reposition(o_x*Angs2Bohr, o_y*Angs2Bohr, o_z*Angs2Bohr, 
	     h1_x*Angs2Bohr, h1_y*Angs2Bohr, h1_z*Angs2Bohr, 
	     h2_x*Angs2Bohr, h2_y*Angs2Bohr, h2_z*Angs2Bohr);
}



//TS////////////////////////////////////////////////////
//
//  takes 9 cartesian coordinates and produces x_origin, y, z, and phi, theta, psi (Euler angles)
//  the internal cartesians are then computed from the internal configuration
//
void Water::Reposition(double o_x, double o_y, double o_z, double h1_x, double h1_y, double h1_z, double h2_x, double h2_y, double h2_z)
{
	// convert water from cartitian coordinate to euler coordinate system.
	// based on fer's code.
	// water is in yz plan while z axis start from mass center and pointing toward 
	// oxygen.
	// x,y,z form a right handed system.
	// This should be the same coordinate system as the one used in gaussian98

	// variable's name are based on fer's code.

	// the 2 Hs are treated as different.

	double Rp1x, Rp1y, Rp1z, Rp3x, Rp3y, Rp3z;
	Rp1x = h1_x - o_x, Rp1y = h1_y - o_y, Rp1z = h1_z - o_z;
	Rp3x = h2_x - o_x, Rp3y = h2_y - o_y, Rp3z = h2_z - o_z;

	double zaxis_m, zaxis_x, zaxis_y, zaxis_z;
	zaxis_x = Rp1x + Rp3x;
	zaxis_y = Rp1y + Rp3y;
	zaxis_z = Rp1z + Rp3z;
	zaxis_m = sqrt(zaxis_x * zaxis_x + zaxis_y * zaxis_y + zaxis_z * zaxis_z);
	zaxis_x/=zaxis_m; zaxis_y/=zaxis_m; zaxis_z/=zaxis_m;
	zaxis_x=-zaxis_x; zaxis_y=-zaxis_y; zaxis_z=-zaxis_z;
//	Got normalized z axis.
//	printf("%f  %f  %f\n", zaxis_x, zaxis_y, zaxis_z); 

	double Rp1_p; // p stands for projection
	Rp1_p = Rp1x*zaxis_x + Rp1y*zaxis_y + Rp1z*zaxis_z;
	double yaxis_m, yaxis_x, yaxis_y, yaxis_z;
	yaxis_x = Rp1x - Rp1_p*zaxis_x;
	yaxis_y = Rp1y - Rp1_p*zaxis_y;
	yaxis_z = Rp1z - Rp1_p*zaxis_z;

	yaxis_m = sqrt(yaxis_x * yaxis_x + yaxis_y * yaxis_y + yaxis_z * yaxis_z);
	yaxis_x/=yaxis_m; yaxis_y/=yaxis_m; yaxis_z/=yaxis_m;
//	Got normalized y axis.
//	printf("%f  %f  %f\n", yaxis_x, yaxis_y, yaxis_z); exit(1);

	double xaxis_y, xaxis_z; // 
	double xaxis_x = yaxis_y*zaxis_z-yaxis_z*zaxis_y; // It seems to me for water it is never used. 
	xaxis_y = yaxis_z*zaxis_x-yaxis_x*zaxis_z;
	xaxis_z = yaxis_x*zaxis_y-yaxis_y*zaxis_x;
//	Got normalized x axis.
	//printf("%f  %f  %f\n", xaxis_x, xaxis_y, xaxis_z); exit(1);

//	Now determine euler angle based on x,y,z axises.
//	the three euler angle is phi, theta, psi.
//	it use the same defination as the rotation functions.
	
	double phi, theta, psi;  // shadow the internal parameters, O dear,
	double sintheta;
	if (fabs(zaxis_z) > 1.0) {
	  cout << "Water::Reposition  zaxis_z = " << zaxis_z << endl; cout.flush();
	}

	theta = acos(zaxis_z);	// 0<= theta < pi
	sintheta = sin(theta);
//	if (large(fabs(theta)))
	if (fabs(sintheta) > 1e-5)
	{
	  // search for nans:
	  // phi = atan2( zaxis_y, zaxis_x);  
	  // psi = atan2( yaxis_z,-xaxis_z);
	  //		sintheta = sin(theta);
	  double tempcos;
	  tempcos = -xaxis_z/sintheta; 
	  
	 // 	cout << "tempcos = " << tempcos << endl; // tempcos = -1 is a problem --> psi = nan
	  if (tempcos > 1.0) {
	    cout << "Water::Reposition  tempcos-1 = " << tempcos-1.0 << endl; cout.flush();
	  }
          double aCosTemp = acos(tempcos);
	  
	  if (tempcos <= -1) aCosTemp=PI;
	  if (yaxis_z >= 0) psi = aCosTemp;	// sintheta always greater than 0
	  else psi = 2*PI - aCosTemp;
	  //printf("%f\n", psi); exit(1);
	  
	  tempcos = zaxis_x/sintheta;
	  tempcos = min(tempcos,1);	//this is useful under very rare circumstance (zaxis_z almost 1)
	  tempcos = max(tempcos,-1);

	  if (fabs(tempcos) > 1.0) {
	    cout << "Water::Reposition  second tempcos = " << tempcos << endl; cout.flush();
	  }


	  if (zaxis_y >= 0) 
	    phi = acos(tempcos);
	  else 
	    phi = 2*PI - acos(tempcos);
	}
	else  // that means we are close to north or south pole, i.e. theta = 0 or theta = PI
	{
	  
	  if (fabs(yaxis_y) > 1.0) {
	    cout << "Water::Reposition  yaxis_y = " << yaxis_y << endl; cout.flush();
	  }

	  phi = 0.0;
	  if (xaxis_y >= 0) psi = acos(yaxis_y);	// for water it does not matter
	  else psi = 2*PI - acos(yaxis_y);			// we code like this for generalness
	}


	//
	//  choice of reference coordiante system
	//  this is defined in ttm_data and needs to be consistent with it
	//
	// use the mass center as the origin
	this->x_origin = (h1_x+h2_x+16.0*o_x) / 18.0;
	this->y_origin = (h1_y+h2_y+16.0*o_y) / 18.0;
	this->z_origin = (h1_z+h2_z+16.0*o_z) / 18.0;
	// this is the ab initio origin: center of charge
	//this->x_origin = (h1_x+h2_x+8.0*o_x) / 10.0;
	//this->y_origin = (h1_y+h2_y+8.0*o_y) / 10.0;
	//this->z_origin = (h1_z+h2_z+8.0*o_z) / 10.0;
	// use O atom as origin
	//this->x_origin = o_x;
	//this->y_origin = o_y;
	//this->z_origin = o_z;

		
	this->phi = phi; 
	this->theta = theta;
	this->psi = psi;  
	this->UpdatePositions();	// this is a better approach. 
								// it will fix the geo.
}


void Water::PrintPositions_InA(ostream& astream)
{
  astream 
    << setprecision(8)
    << "O  " << AtomCenters[0].x*Bohr2Angs <<"   "<< AtomCenters[0].y*Bohr2Angs <<"   "<< AtomCenters[0].z*Bohr2Angs<< "\n"
    << "H  " << AtomCenters[1].x*Bohr2Angs <<"   "<< AtomCenters[1].y*Bohr2Angs <<"   "<< AtomCenters[1].z*Bohr2Angs<< "\n"
    << "H  " << AtomCenters[2].x*Bohr2Angs <<"   "<< AtomCenters[2].y*Bohr2Angs <<"   "<< AtomCenters[2].z*Bohr2Angs<< "\n";
  //	astream << 'M' << "  " << AtomCenters[3].x <<"   "<< AtomCenters[3].y <<"   "<< AtomCenters[3].z<< "\n";


  //  printf("O %10.6f %10.6f %10.6f\n", AtomCenters[0].x*Bohr2Angs, AtomCenters[0].y*Bohr2Angs, AtomCenters[0].z*Bohr2Angs);
  //  printf("H %10.6f %10.6f %10.6f\n", AtomCenters[1].x*Bohr2Angs, AtomCenters[1].y*Bohr2Angs, AtomCenters[1].z*Bohr2Angs);
  //  printf("H %10.6f %10.6f %10.6f\n", AtomCenters[2].x*Bohr2Angs, AtomCenters[2].y*Bohr2Angs, AtomCenters[2].z*Bohr2Angs);

}


////////////////////////////////////////////////////
//
//  this is a Conf2Cart function
//  the internal cartesian coordinates are computed from the interal
//  center-of-mass position (xyz-origin) and the internal Euler angles (phi, theta, psi)
//
void Water::UpdatePositions()
{
  double phi, theta, psi; // why, Oh, why Seymour?
  phi = this->phi, theta = this->theta, psi = this->psi;
  ResetPositions();
  // this->Reorient(phi, theta, psi);
  // ApplyEuler only rotates the atoms, not the basis functions of the repulvive core
  ApplyEuler(phi, theta, psi);
}

////////////////////////////////////////////
//
//  put standard orientation, translated to the current center-of-mass coordinate
//  into internal xyz (storred in AtomCenter)
//  so you need Waters[i].AtomCenters[j].x in WaterCluster 
//
void Water::ResetPositions()
{
  this->AtomCenters[0].x = x_origin + RefPosOx;
  this->AtomCenters[0].y = y_origin + RefPosOy;
  this->AtomCenters[0].z = z_origin + RefPosOz;
  
  this->AtomCenters[1].x = x_origin + RefPosHx;
  this->AtomCenters[1].y = y_origin + RefPosHy;
  this->AtomCenters[1].z = z_origin + RefPosHz;
  
  this->AtomCenters[2].x = x_origin + RefPosHx;
  this->AtomCenters[2].y = y_origin - RefPosHy;
  this->AtomCenters[2].z = z_origin + RefPosHz;
  
  this->AtomCenters[3].x = x_origin + RefPosMx;
  this->AtomCenters[3].y = y_origin + RefPosMy;
  this->AtomCenters[3].z = z_origin + RefPosMz;

  // TS not needed for stand-alone DPP
  //this->UpdateBasisfunctions();
  //this->Calculate_C();
  
  this->phi   = 0.0;
  this->theta = 0.0;
  this->psi   = 0.0;		
}



void Water::ApplyEuler(double angle1, double angle2, double angle3)
{
  // 
  //  this is used to rotate from the reference orientation (Euler: 0,0,0)
  //  to the Euler angles 1, 2, and 3
  //
  // Seymour:
  // Eular angle defination in Mathematical methods for physists. George B. Arfken
  // We are interested in what is the old coordinate after we rotate the coordinate system
  // We want to convert from new to old!.
  // The matrix in Arfken's book convert old to new.
  //
  //TS
  //     another place where this rotation matrix can be found is Wilson, Bright, Decius, Cross; Appendix 1
  //
  //   in WBDC angle1 = phi,   angle2 = theta,   angle3 = chi
  //
  //   since, the calculation is done in space-fixed coordinates, I assume that 
  //   old is molecule-fixed and new is space-fixed (with identical origin)

  
  this->phi += angle1, this->theta += angle2, this->psi += angle3;
  double NewtoOld[9] ;
  double c1=cos(angle1), c2=cos(angle2), c3=cos(angle3);
  double s1=sin(angle1), s2=sin(angle2), s3=sin(angle3);
  

  NewtoOld[0]=c1*c2*c3-s1*s3, NewtoOld[1]=-c1*c2*s3-s1*c3, NewtoOld[2]=c1*s2;
  NewtoOld[3]=s1*c2*c3+c1*s3, NewtoOld[4]=-s1*c2*s3+c1*c3, NewtoOld[5]=s1*s2;
  NewtoOld[6]=-s2*c3, NewtoOld[7]=s2*s3, NewtoOld[8]=c2;
 
  double newx, newy, newz;
  double oldx, oldy, oldz;
  int icount;
  for(icount=0; icount<nAtomCenters; icount++)
    {
      oldx = AtomCenters[icount].x - x_origin;  // origin is the center-of-mass position of the molecule
      oldy = AtomCenters[icount].y - y_origin;  // so oldx is the position in the molecule system
      oldz = AtomCenters[icount].z - z_origin;

      newx = NewtoOld[0]*oldx+NewtoOld[1]*oldy+NewtoOld[2]*oldz;
      newy = NewtoOld[3]*oldx+NewtoOld[4]*oldy+NewtoOld[5]*oldz;
      newz = NewtoOld[6]*oldx+NewtoOld[7]*oldy+NewtoOld[8]*oldz;

      AtomCenters[icount].x = newx + x_origin;
      AtomCenters[icount].y = newy + y_origin;
      AtomCenters[icount].z = newz + z_origin;
    }
}
