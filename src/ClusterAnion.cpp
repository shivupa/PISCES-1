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
                                 
#include <mpi.h>                 

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

using namespace std;

//namespace {
//DVR Hel ;} ;



//ClusterAnion::~ClusterAnion()
//{
//   int rank;
//  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//  if(rank==0) std::cout << "Hello from ClusterAnion\n";
//}


void ClusterAnion::SetUpClusterAnion(const dVec WaterPos, const Parameters InP)
{
   int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  //  Sets-up WaterN, Vel, and Hel
  Para = InP ;
  int nMolecules = Para.nWater;
  Molecules = new Molecule[nMolecules];
  for (int im = 0; im < nMolecules; ++im) {
    if(rank==0)cout << im << "\n"; cout.flush();
    Molecules[im].SetupDPPWater(&WaterPos[9*im]);
  }
  WaterN.SetStructure(Para.nWater, &WaterPos[0], 1, Para.CenterFlag, Para.KTFlag, Para.WMverbose);
  Hel.SetupDVR(Para.ngrid, Para.DVRType, Para.Sampling, Para.gpara, Para.gridverbose);
  Hel.SetVerbose(Para.gridverbose);
  Hel.DiagonalizeSetup(Para.nStates, Para.DiagMethod, Para.maxSub, Para.maxIter, Para.ptol);
  Vel.SetVerbose(Para.PotVerbose);
  //delete[] Molecules ; 
} 
 
double ClusterAnion::SinglePoint(const dVec WaterPos)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );



  int nWater = Para.nWater ; 
  dVec  WaterConf(nWater*6);
  WaterN.GetConfiguration(nWater, &WaterConf[0], 1.0, 0);
  double Etotal = EnergyFromConfiguration(&WaterConf[0]) ; 
 

  if(rank==0)cout << "\nElectron binding energies:\n";
  for (int istate = 0; istate < Para.nStates; ++istate)
    if(rank==0) printf("State %-3i  %12.8f Hartree =  %10.2f meV\n", istate, ebes[istate], ebes[istate]*AU2MEV);
  if (nconverged < Para.nStates)
    {if(rank==0)cout << "single-point:: Only " << nconverged << " States of " << Para.nStates << " converged.\n"; exit(1);}

  if (Para.nStates > 0)
    Hel.ExpectationValues(1);

  //Hel.EnergyPartitioning(Vel) ;
  // prepare cubefiles showing the wavefunctions
  if (Para.CubeFile > 0) {
    dVec WaterCoor(9*Para.nWater);
    WaterN.GetStructure(Para.nWater, &WaterCoor[0], Bohr2Angs, 0);
    iVec NucCharge(3*Para.nWater);
    for (int k = 0; k < Para.nWater; ++k) {
      NucCharge[3*k+0] = 8;
      NucCharge[3*k+1] = 1;
      NucCharge[3*k+2] = 1;
    }
    if(rank==0)cout << "\nWriting the following gaussian-like cube files:\n";
    for (int istate = 1; istate <= Para.nStates; ++istate) {
      char fname[30];
      if(rank==0) sprintf(fname, "WaveFn%02i.gcube", istate);
      if(rank==0)cout << "  " << fname << "\n";
      Hel.WriteCubeFile(istate, fname, 3*Para.nWater, &NucCharge[0], &WaterCoor[0]);
    }
  }

  // prepare cuts (for gnuplot splot)
  if (Para.WfCuts > 0)
    Hel.WriteCuts();



  return Etotal ;
  
}


double ClusterAnion::PrintConf(int nWater, double *WaterConf)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  if(rank==0)cout << "number of waters in  ClusterAnion::PrintConf =  " << nWater << "\n";
  for(int k=0; k < nWater ;++k){
      if(rank==0)cout<<" ClusterAnion_PrintConf "<< WaterConf[k*6]<<" "<<WaterConf[k*6+1]<<" "<< WaterConf[k*6+2]<<endl;
      if(rank==0)cout<<" ClusterAnion_PrintConf "<< WaterConf[k*6+3]<<" "<< WaterConf[k*6+4]<<" "<< WaterConf[k*6+5]<<endl;
   } ;
 return nWater*2.0 ;
}



/////////////////////////////////////////////////
//
// Computes the energy from a configuration (CoM coordinates and Euler angles)
// It assumes that a SinglePoint has been run to setup all parameters
// of the potential and the diagonalizer
// printing is minimized by setting all verbose parameters explicitly to zero
//
double ClusterAnion::EnergyFromConfiguration(double *WaterConf)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  int nWater = WaterN.ReportNoOfWaters();
//  if(rank==0)cout << "number of waters in  ClusterAnion::EnergyFromConfiguration =  " << nWater << "\n";
//  for(int k=0; k < nWater ;++k){
//      if(rank==0)cout<<" ClusterAnion_EnergyFromConfiguration "<< WaterConf[k*6]<<" "<<WaterConf[k*6+1]<<" "<< WaterConf[k*6+2]<<endl;
//      if(rank==0)cout<<" ClusterAnion_EnergyFromConfiguration "<< WaterConf[k*6+3]<<" "<< WaterConf[k*6+4]<<" "<< WaterConf[k*6+5]<<endl;
//   } ;
  WaterN.SetConfiguration(nWater, WaterConf, 1, 0);
  double E0 = WaterN.CalcEnergy(Para.WMverbose);

// Save WaterWater Polarization ----- Tae Hoon Choi
  Vel.VpolWaterWater=WaterN.ReportPolarization();

  if(rank==0)cout << "Energy of the neutral cluster E0 = " << E0 << "\n";
  int nSites = 0 , nCharges = 0, nDipoles = 0;
  WaterN.ReportNoOfSites(nSites, nCharges, nDipoles);
  int nPntPol = WaterN.ReportNoOfPolSites() ;

  if (Para.PotVerbose > 0)
    if(rank==0)cout << "The electron sees " << nCharges << " charges, " << nDipoles << " dipoles, and there are "
         <<  nSites << " sites in the gradient calculations.\n";

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
                  &DmuByDR[0], 0);

  nDipoles *= 3;
  Dipoles.resize(3*nDipoles);
  ids.resize(nDipoles);
  int Tngrid[3];
  double Tgpara[3];

  WaterN.ReportInducedDipoles(nDipoles, &Dipoles[0], &ids[0], Para.PotVerbose);
  WaterN.ReportPolSitesAndField(nPntPol, &Alphas[0], &ips[0], &Epc[0], 1);

// Now we are checking it is for dual-grid or full-gird  -- Tae Hoon Choi
// 5 and 6 is for dual-grid, but we are not using 5 which is using non-SC as the initial guess
  int Icombine = Para.PotFlag[3];
  if (Icombine == 5) Para.PotFlag[3]= 1;

   int OriginTol=Para.ptol;

// 6 is for dual-grid method, if it is even grid then triple spacing 
// if it is odd number gird, then it will be double spacing  -- Tae Hoon Choi
  if (Icombine == 6) {
    Para.PotFlag[3]= 3;
    double spacing = Para.gpara[0]/(Para.ngrid[0]+1);
    if(rank==0)cout<<"spacing : "<<spacing<<endl;
    if (Para.ngrid[0]%2 == 1 ) {
      Tngrid[0]=(Para.ngrid[0]+1)/2;
      Tngrid[1]=(Para.ngrid[1]+1)/2;
      Tngrid[2]=(Para.ngrid[2]+1)/2;
      Tgpara[0]=Tngrid[0]*2*spacing+2*spacing;
      Tgpara[1]=Tngrid[1]*2*spacing+2*spacing;
      Tgpara[2]=Tngrid[2]*2*spacing+2*spacing;
    }
    else {
      Tngrid[0]=(Para.ngrid[0]/2+1)/3*2;
      Tngrid[1]=(Para.ngrid[1]/2+1)/3*2;
      Tngrid[2]=(Para.ngrid[2]/2+1)/3*2;
      Tgpara[0]=Tngrid[0]*3*spacing+3*spacing;
      Tgpara[1]=Tngrid[1]*3*spacing+3*spacing;
      Tgpara[2]=Tngrid[2]*3*spacing+3*spacing;
    }
     
   //Para.ptol=6;

   Hel.SetupDVR2(Tngrid, Para.DVRType, Para.Sampling, Tgpara, Para.gridverbose);
   Hel.DiagonalizeSetup(Para.nStates, Para.DiagMethod, Para.maxSub, Para.maxIter, Para.ptol);

  }

  Vel.Setup(Para.PotFlag, nSites, &Sites[0],
            nCharges, &Charges[0], &iqs[0],
            nDipoles, &Dipoles[0], &ids[0],
            nPntPol, &Alphas[0], &ips[0], &Epc[0],
           &DmuByDR[0],Para.PotPara);
   if(rank==0)cout << "  ************  BEGIN ENERGY ***************" << endl;
   Hel.ComputePotential(Vel);

  ebes.resize(Para.nStates);
  if(rank==0)cout<< "Para.istartvec = "<<Para.istartvec<<endl;
  nconverged = Hel.Diagonalize(Para.istartvec, &ebes[0]);
  
  // Now we can save the wavefn - Tae Hoon Choi
  Para.istartvec = 0;
   // nconverged = Hel.Diagonalize(Para.istartvec, &ebes[0]);

// This is for non-selfconsistent polarization as a initial guess 
// but it not effecctive, so not using it  -Tae Hoon Choi
  if (Icombine == 5) {
    Para.PotFlag[3]= 5;
    Vel.Setup(Para.PotFlag, nSites, &Sites[0],
            nCharges, &Charges[0], &iqs[0],
            nDipoles, &Dipoles[0], &ids[0],
            nPntPol, &Alphas[0], &ips[0], &Epc[0],
           &DmuByDR[0],Para.PotPara);
     if(rank==0)cout << "  ************  refining ENERGY ***************" << endl;
     Hel.ComputePotential(Vel);

    ebes.resize(Para.nStates);
    nconverged = Hel.Diagonalize(Para.istartvec, &ebes[0]);
  }

// this if for dual-grid method from Tae Hoon Choi
  if (Icombine == 6) {

     if(rank==0)cout<<"  ***  Energies For pre-calculations    ***  "<<endl;
     if(rank==0)cout << "       EBE      = "<<  ebes[0]*AU2MEV        << " meV"<<endl;

    Para.PotFlag[3]= 6;
    Para.ptol=OriginTol;


   Hel.SetupDVR2(Para.ngrid, Para.DVRType, Para.Sampling, Para.gpara, Para.gridverbose);
   Hel.DiagonalizeSetup(Para.nStates, Para.DiagMethod, Para.maxSub, Para.maxIter, Para.ptol);

    Vel.Setup(Para.PotFlag, nSites, &Sites[0],
            nCharges, &Charges[0], &iqs[0],
            nDipoles, &Dipoles[0], &ids[0],
            nPntPol, &Alphas[0], &ips[0], &Epc[0],
           &DmuByDR[0],Para.PotPara);
     if(rank==0)cout << "  ************  refining ENERGY ***************" << endl;
     Hel.ComputePotential(Vel);
    ebes.resize(Para.nStates);
    Para.istartvec = 3; // using the converged wavefn with coarse grid + interpolations
    nconverged = Hel.Diagonalize(Para.istartvec, &ebes[0]);
    Para.istartvec = 4; // using the converged coarse grid wavefn for optimization 
   // nconverged = Hel.Diagonalize(Para.istartvec, &ebes[0]);
   }
  double EBE = ebes[0];

   if (nconverged < 1)
    {if(rank==0)cout << "GetClusterEnergy:: Convergence failure.\n"; exit(1);}


     if(rank==0)cout<<"                                              "<<endl;
     if(rank==0)cout<<"  ------------------------------------------  "<<endl;
     if(rank==0)cout<<"  ***  Energies For The Current Geometry ***  "<<endl;
     if(rank==0)cout<<"  ------------------------------------------  "<<endl;
     if(rank==0)cout << "       Eneutral = "<<  E0*AU2MEV         << " meV"<<endl;
     if(rank==0)cout << "       EBE      = "<<  EBE*AU2MEV        << " meV"<<endl;
     if(rank==0)cout << "       Etotal   = "<<  (EBE+E0)*AU2MEV   << " meV"<<endl;
     if(rank==0)cout<<"  ------------------------------------------  "<<endl;
     if(rank==0)cout << "  **************  END ENERGY ***************" << endl;
     if(rank==0)cout<<"                                              "<<endl;

//   DumpGeometry(nWater, E0+EBE, &Sites[0] );

   return E0+EBE;
}

void ClusterAnion::GetAnalGrad( const double *WaterConf, double *analgrad) 
{


  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  int nWater = WaterN.ReportNoOfWaters();
  int nSites = 0, nCharges = 0, nDipoles = 0;
  WaterN.ReportNoOfSites(nSites, nCharges, nDipoles);
  int nPntPol = WaterN.ReportNoOfPolSites();
  static dVec Sites; Sites.resize(3*nSites);
  static dVec Charges; Charges.resize(nCharges);
  static dVec Dipoles; Dipoles.resize(3*nDipoles);
  static dVec DmuByDR; DmuByDR.resize(3*nSites*3*nDipoles);
  static dVec dEfield; dEfield.resize(3*nSites*nCharges*3);
  static dVec dxWEfield; dxWEfield.resize(3*nCharges*nSites); //vkv
  static dVec dyWEfield; dyWEfield.resize(3*nCharges*nSites); //vkv
  static dVec dzWEfield; dzWEfield.resize(3*nCharges*nSites); //vkv
  static dVec dT_x; dT_x.resize(9*nPntPol*nPntPol);
  static dVec dT_y; dT_y.resize(9*nPntPol*nPntPol);
  static dVec dT_z; dT_z.resize(9*nPntPol*nPntPol);


  static dVec dTensor; dTensor.resize(3*nSites*3*nDipoles*3);
  static iVec iqs; iqs.resize(nCharges);
  static iVec ids; ids.resize(nDipoles);
  static dVec Alphas; Alphas.resize(nPntPol);
  static iVec ips; ips.resize(nPntPol);
  static dVec Epc; Epc.resize(3*nPntPol);
  WaterN.GetLists(nSites, &Sites[0], nCharges, &Charges[0], &iqs[0], nDipoles, &Dipoles[0], &ids[0],
                  &DmuByDR[0], 0);

  // Get all the dipoles
  nDipoles *= 3;
  Dipoles.resize(3*nDipoles);
  ids.resize(nDipoles);
  int verbose = 0 ;
  WaterN.ReportInducedDipoles(nDipoles, &Dipoles[0], &ids[0], verbose);

  static dVec Gradient; Gradient.resize(nSites*3);
  static dVec PolGrad; PolGrad.resize(nSites*3);
  for (int i=0; i < nSites*3; i++){
     Gradient[i]=0.0 ;
     PolGrad[i]=0.0 ;
 }


  WaterN.SetConfiguration(nWater, WaterConf , 1, 0);


//  WaterN.GetDInducedDipoles( &DmuByDR[0], &dEfield[0], &dTensor[0] );
   WaterN.Calc_dInducedDipoles(0, &DmuByDR[0], &dEfield[0], &dTensor[0] );

//  if(rank==0)cout << "ncharges, nSites, nPntPol, nDipoles" << nCharges << " " << nSites <<" "<< nPntPol << " " << nDipoles << endl ; 

  for (int j=0; j < nSites; j++){
    for (int i=0; i < nCharges; i++){
    dxWEfield[j*nCharges*3 + i*3+0] = dEfield[nCharges*9*j + 3*i+0];
    dxWEfield[j*nCharges*3 + i*3+1] = dEfield[nCharges*9*j + 3*i+1];
    dxWEfield[j*nCharges*3 + i*3+2] = dEfield[nCharges*9*j + 3*i+2];
    dyWEfield[j*nCharges*3 + i*3+0] = dEfield[nCharges*9*j + 3*nCharges + 3*i+0];
    dyWEfield[j*nCharges*3 + i*3+1] = dEfield[nCharges*9*j + 3*nCharges + 3*i+1];
    dyWEfield[j*nCharges*3 + i*3+2] = dEfield[nCharges*9*j + 3*nCharges + 3*i+2];
    dzWEfield[j*nCharges*3 + i*3+0] = dEfield[nCharges*9*j + 6*nCharges + 3*i+0];
    dzWEfield[j*nCharges*3 + i*3+1] = dEfield[nCharges*9*j + 6*nCharges + 3*i+1];
    dzWEfield[j*nCharges*3 + i*3+2] = dEfield[nCharges*9*j + 6*nCharges + 3*i+2];
    }
  }

  double E0 = WaterN.CalcGradient(nSites*3, &Gradient[0], &PolGrad[0], 0);

// get the polarizable sites and the point-charge field (only used for Polarization=3)
  WaterN.ReportPolSitesAndField(nPntPol, &Alphas[0], &ips[0], &Epc[0], 1);

  dVec Rpps; Rpps.resize(3*nPntPol);  // needed for calling ComputeInvA
  for (int i = 0; i < nPntPol; ++i) {
    Rpps[3*i+0] = Sites[3*ips[i]+0];
    Rpps[3*i+1] = Sites[3*ips[i]+1];
    Rpps[3*i+2] = Sites[3*ips[i]+2];
  }

void BuildSuperdT(int nPolSites, const double *R, const double *alpha, double aThole, double *dT_x, double *dT_y, double *dT_z) ;
 BuildSuperdT(nPntPol, &Rpps[0], &Alphas[0], 0.3, &dT_x[0] , &dT_y[0] , &dT_z[0]) ;


   if(rank==0)cout << "  **********  BEGIN GRADIENT  ************" << endl;
   Hel.ComputeGradient(Vel, nSites, &Gradient[0], &dT_x[0], &dT_y[0], &dT_z[0], &PolGrad[0], WaterN , &dEfield[0]);


   WaterN.ConvertF(&Gradient[0], analgrad);

   // printing gradients
  if(rank==0)cout<<"                                              "<<endl;
  if(rank==0)cout<<"  ------------------------------------------  "<<endl;
  if(rank==0)cout<<"  *** Gradients For The Current Geometry ***  "<<endl;
  if(rank==0)cout<<"  ------------------------------------------  "<<endl;
   for(int k=0; k < nWater*2 ;++k){
  if(rank==0)cout<<"  AnalGrad "<<analgrad[k*3]<<" "<<analgrad[k*3+1]<<" "<<analgrad[k*3+2]<<endl;
   }
  if(rank==0)cout<<"  ------------------------------------------  "<<endl;
  if(rank==0)cout << "  **********  END GRADIENT  **************" << endl;
  if(rank==0)cout<<"                                              "<<endl;
//  exit(1);



}


void ClusterAnion::GetNumGrad(const double *configuration, double *numgrad)
{


  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  int nWater = WaterN.ReportNoOfWaters();
//  double numgrad[nWater*2*3] ;
  for(int k=0; k < nWater*2 ;++k){
     if ( k < nWater*2) {
           numgrad[k*3]=0.00; numgrad[k*3+1]=0.00; numgrad[k*3+2]=0.00; // initializing gradients
     }
  }

  double *Conf = new double[nWater*6];          //copying the configuration to Conf          
   for (int i = 0; i < nWater*6; ++i)           // which would be used to modify
              Conf[i] = configuration[i];

  double transstep = 0.002;                     // Bohr  0.005
  double rotstep   = 0.02 / 180.0 * PI;         // deg   0.05
  double eplus, eminus, h2;
  int igrad = 0;

  double E0;


   for (int i = 0; i < nWater; ++i) {                    // each rigid-water has 6 grads
     for (int k = 0; k < 3; ++k) {                       // 3 due to trans_steps
       Conf[igrad] = configuration[igrad] + transstep;
       eplus = EnergyFromConfiguration(Conf);
       h2 = Conf[igrad];

       Conf[igrad] = configuration[igrad] - transstep;
       eminus = EnergyFromConfiguration(Conf);
       h2 -= Conf[igrad];

       Conf[igrad] = configuration[igrad];
       numgrad[igrad] = (eplus-eminus) / h2;
         igrad ++;
     }
     for (int k = 0; k < 3; ++k) {                      // 3 due to rot_steps
       Conf[igrad] = configuration[igrad] + rotstep;
       eplus = EnergyFromConfiguration(Conf);
       h2 = Conf[igrad];

       Conf[igrad] = configuration[igrad] - rotstep;
       eminus = EnergyFromConfiguration(Conf);
       h2 -= Conf[igrad];
    
       Conf[igrad] = configuration[igrad];
       numgrad[igrad] = (eplus-eminus) / h2;
         igrad ++;
     }
   }

  for(int k=0; k < nWater*2 ;++k){
     if(rank==0)cout<<" Numgrad "<<numgrad[k*3]<<" "<<numgrad[k*3+1]<<" "<<numgrad[k*3+2]<<endl;
  }

  delete[] Conf;


}

