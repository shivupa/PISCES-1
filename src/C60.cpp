#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
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
#include "C60.h"
#include "ChargeDipPol.h"
#include "timer.hpp"




using namespace std;

//////////////////////////////
//
//  IonTypes: 
//
const int C_atom = 1 ;
const int H_atom = 2 ;
const int CH_atom = 3 ;







///////////////////////////////////////////////
//
//  single point for C60-electron
//  has so far just the VEA
//
void C60_sp(TSIN Input)
{

  printf("-----------------\nC60 single point\n-----------------\n");
  printf("-----------------\nInput Geom should be in Bohrs\n-----------------\n");

  // read the number and coordinates of the C atoms from Input
  char **coorlines;
  int nLines     = Input.GetGroup("Fuller", &coorlines);
  double *IonPos = new double[3*nLines];  // coordinates for Catoms
  int *IonType   = new int[nLines];
  int nIons      = GetAtomCoordinates(nLines, coorlines, IonPos, IonType);
  cout << "nIons = " << nIons << endl;
  if (nIons < nLines)
  printf("Found %i atoms reading %i Fuller-lines.\n", nIons, nLines);
  printf("Cluster Geometry [bohrs]\n");
  PrintAtomCoordinates(nIons, IonPos, IonType, 1.00);

  // compute the energy of the neutral cluster, presently nothing for 
  // neutral cluster
  //
  //  Now, dealing with the excess electron is more involved
  //
  // set up the Potential for the excess electron

  double *Charges = new double[nIons];
  int *iqs = new int[nIons];
  for (int i = 0; i < nIons; ++i) 
  {
    iqs[i] = i;
    switch (IonType[i]) 
    {
      case C_atom:
        Charges[i] = 0.0;
        break;
      case H_atom:
        Charges[i] = 0.0;
        break;
      case CH_atom:
        Charges[i] = 0.0;
        break;
      default:
      printf("This should never happen\n");
      exit(42);
    }
  }
  class Potential Vel;
  double PotPara[100];
  int CubeFile   = Input.GetInt("Job", "CubeFile", 0);
  int PotFlag = Input.GetInt("CsixtyPotential", "Potential", 1);
  int PotVerb = Input.GetInt("CsixtyPotential", "Verbose", 1);
  Vel.SetVerbose(PotVerb);
  int nMolecules  = Input.GetInt("CsixtyPotential", "nMolecules", 1);
  int nAtomsArray[nMolecules]; 
  double DipoleArray[nMolecules];
  double CenterArray[nMolecules*3];
  int sumatoms = 0 ; 
  switch (PotFlag)
    {
    case 1021:
      PotPara[0] = Input.GetDouble("CsixtyPotential", "EsDamp", 1.7724) ; // I  think this is for radius, and should be in bohr
      PotPara[1] = Input.GetDouble("CsixtyPotential", "RepScale", 13.8) ; //repulsive core scaling paramter
      PotPara[2] = Input.GetDouble("CsixtyPotential", "GaussExp", 1.62) ; //repulsive core exponent
      PotPara[3] = Input.GetDouble("CsixtyPotential", "PolDamp", 1.77) ; //polarization damping paramter
      PotPara[4] = Input.GetInt("CsixtyPotential", "Polarization", 1) ; //polarization flag
      PotPara[5] = Input.GetInt("CsixtyPotential", "nMolecules", 1) ; //number of molecules
      PotPara[6] = Input.GetInt("CsixtyPotential", "Verbose", 0) ; //number of molecules
      PotPara[7] = Input.GetDouble("CsixtyPotential", "ExtraNucQ", 0.00) ; //extra nuclear charge
      PotPara[8] = Input.GetDouble("CsixtyPotential", "RepFlag", 1) ; //extra nuclear charge
      PotPara[10] = Input.GetDouble("CsixtyPotential", "Rtol", 10.0) ; //extra nuclear charge
      PotPara[11] = Input.GetDouble("CsixtyPotential", "CarbonType", 1) ; //extra nuclear charge

      PotPara[12] = Input.GetDouble("CsixtyPotential", "AlphaComp", 4) ; //A switch that determines which component of the Polarizability will be evaluated, 1 - charge flow, 2 - cross-term, 3 - Induced Dipole, 4 - Sum.
      
      Input.GetIntArray("CsixtyPotential", "nAtomsArray", nAtomsArray, nMolecules);
      Input.GetDoubleArray("CsixtyPotential", "DipoleArray", DipoleArray, nMolecules);
      Input.GetDoubleArray("CsixtyPotential", "CenterArray", CenterArray, nMolecules*3);
      for (int i = 0; i < nMolecules ; ++i) {
        sumatoms += nAtomsArray[i] ;
        cout << "number of atoms in Molecule  " << i+1 << " = " << nAtomsArray[i] << "\n" ; 
        cout << "Dipole moment of atoms in Molecule  " << i+1 << " = " << DipoleArray[i] << "\n" ; 
        } 
      if (sumatoms == nIons){
         cout << "total number of atoms = " << sumatoms << "\n";  
        }
      else { 
        cout << "Error: something wrong with number of atoms " << "\n" ;
        exit(1);
        }

      Vel.SetupFullerElec(nIons, IonPos, PotPara, nAtomsArray, DipoleArray, CenterArray, IonType);
      break;
    default:
      cout << " Fullerene Single Point: unknown potential for the excess electron: " << PotFlag << "\n";
      exit(1);
    }



  // set up a grid (a DVR) for the wavefunction of the electron
  int ngrid[3] = {0,0,0}; 
  Input.GetIntArray("GridDef", "NoOfGridPoints", ngrid, 3);
  int gridverbose = Input.GetInt("GridDef", "Verbose", 0);
  int DVRType = Input.GetInt("GridDef", "DVRType", 0);      // 1 HO,  other = Sine 
  double gpara[3] = {0,0,0};  // this is for the grid length or the oscillator freq
  switch (DVRType)
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

  int Tngrid[3];
  double Tgpara[3];

  Hel.SetupDVR(ngrid, DVRType, 1, gpara, gridverbose);

// 6 is for dual-grid method, if it is even grid then triple spacing 
// if it is odd number gird, then it will be double spacing  -- Tae Hoon Choi

  int Icombine = PotPara[4];
  if (Icombine == 6) {
    PotPara[4]= 3;
    double spacing = gpara[0]/(ngrid[0]+1);

   Vel.SetupFullerElec(nIons, IonPos, PotPara, nAtomsArray, DipoleArray, CenterArray, IonType);
    if (ngrid[0]%2 == 1 ) {
      Tngrid[0]=(ngrid[0]+1)/2;
      Tngrid[1]=(ngrid[1]+1)/2;
      Tngrid[2]=(ngrid[2]+1)/2;
      Tgpara[0]=Tngrid[0]*2*spacing+spacing*2;
      Tgpara[1]=Tngrid[1]*2*spacing+spacing*2;
      Tgpara[2]=Tngrid[2]*2*spacing+spacing*2;
    }
    else {
      Tngrid[0]=(ngrid[0]/2+1)/3*2;
      Tngrid[1]=(ngrid[1]/2+1)/3*2;
      Tngrid[2]=(ngrid[2]/2+1)/3*2;
      Tgpara[0]=Tngrid[0]*3*spacing+3*spacing;
      Tgpara[1]=Tngrid[1]*3*spacing+3*spacing;
      Tgpara[2]=Tngrid[2]*3*spacing+3*spacing;
    }

   Hel.SetupDVR2(Tngrid, DVRType, 1, Tgpara, gridverbose);
   cout<<" finish SEtupDVR2 "<<endl;
  }
  Hel.ComputePotential(Vel);




  // now get the input for the diagonalizer
  int nStates = Input.GetInt("Diag", "nStates", 1);
  int DiagFlag = Input.GetInt("Diag", "Method", 1);     // 1=Lanczos  2=Davidson  3= 0th-order Jacobi-Davidson
  int maxSub = Input.GetInt("Diag", "maxSubspace", 20);   // max no of micro-iterations
  int maxIter = Input.GetInt("Diag", "MaxIter", 10);      // max no of macro-iterations 
  int ptol = Input.GetDouble("Diag", "pTol", 5);          // tolerance = 10^-pTol 
  int istartvec = Input.GetInt("Diag", "StartVector", 1); // this is different for Lanczos + Davidson and needs work
  int diagverbose = Input.GetInt("Diag", "Verbose", 1);
  Hel.SetVerbose(diagverbose);  // resets verbose from gridverbose to diagverbose

  double *energies = new double[nStates];
  //  Hel.DiagonalizeSetup(nStates, DiagFlag, maxSub, maxIter, ptol);
  Hel.DiagonalizeSetup(nStates, DiagFlag, maxSub, 50, ptol);
  Hel.Diagonalize(istartvec, energies);
  for (int istate = 0; istate < nStates ; ++istate){
    cout << "Eel = " << energies[istate] << " Hartree = " << energies[istate]*AU2MEV << " meV\n";
  }

  istartvec = 5;

  if (Icombine == 6) {
    PotPara[4]= 6;

    Vel.SetupFullerElec(nIons, IonPos, PotPara, nAtomsArray, DipoleArray, CenterArray, IonType);

    Hel.SetupDVR2(ngrid, DVRType, 1, gpara, gridverbose);
    Hel.ComputePotential(Vel);

    istartvec = 3;   // using the converged wavefn with coarse grid + interpolation
    Hel.Diagonalize(istartvec, energies);
    istartvec = 0;   // for saving converged wavefn

    for (int istate = 0; istate < nStates ; ++istate){
      cout << "Eel = " << energies[istate] << " Hartree = " << energies[istate]*AU2MEV << " meV\n";
  }
 }
//  if (nStates > 0)
//    Hel.ExpectationValues(1);
//  Hel.EnergyPartitioning(Vel) ;
  // plot the wavefunction
  //
  int *NucCharge = new int[nIons];
  for (int k = 0; k < nIons; ++k) {
    switch (IonType[k])
      {
      case C_atom:
	NucCharge[k] = 6;
	break;
      case H_atom:
	NucCharge[k] = 1;
	break;
      case CH_atom:
        NucCharge[k] = 6;
        break;
      default:
	printf("This should never happen\n");
	exit(42);
      }
  }
  for (int i = 0 ; i < 3*nIons; ++i) {
    IonPos[i] = IonPos[i]*Bohr2Angs ;
  }

  cout<<"CubeFile : "<<CubeFile<<endl;

  if (CubeFile == 1 ) {
   cout<<"writing cude -------------------"<<endl; 
   for (int istate = 1; istate <= nStates ; ++istate) {
    char fname[15];
      sprintf(fname, "WaveFn%02i.cube", istate);
      Hel.WriteCubeFile(istate, fname, nIons , &NucCharge[0], &IonPos[0]);
    }
    Hel.WriteOneDCuts() ; 
   cout << "All done.\n"; cout.flush();
 }

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
int GetAtomCoordinates(int nLines, char **coorlines, double *IonPos, int *IonType)
{
  int nIons = 0;
  char atomsymbol[5];
  double x = 0, y = 0, z = 0;
  char *inputline = 0;

  for (int iLine = 0; iLine < nLines; ++iLine) {
    inputline = coorlines[iLine];
    if (*inputline != '#' && strlen(inputline) > 1) {
      if (sscanf(inputline," %s %lf %lf %lf", atomsymbol, &x, &y, &z) != 4) {
	cout << "Error reading input group Fuller, line" << iLine << ":" << inputline << "\n";
      exit(1);
      }
      //      IonPos[3*nIons+0] = x * Angs2Bohr;
      //      IonPos[3*nIons+1] = y * Angs2Bohr;
      //      IonPos[3*nIons+2] = z * Angs2Bohr;
      IonPos[3*nIons+0] = x ;
      IonPos[3*nIons+1] = y ;
      IonPos[3*nIons+2] = z ;
      switch (atomsymbol[0])
	{
	case 'C':
	  IonType[nIons] = C_atom;
          switch (atomsymbol[1])
           {
	   case 'H':
	    IonType[nIons] = CH_atom;
            break;
           }
	  break;
	case 'H':
	  IonType[nIons] = H_atom;
	  break;
	default:
	  cout << "Error reading input group Fuller, line" << iLine << ":" << inputline << "\n";
	  exit(1);
	}
      nIons ++;
    }
  }
  return nIons;
}


void PrintAtomCoordinates(int nIons, double *IonPos, int *IonType, double scale)
{
  for (int i = 0; i < nIons; ++i) {
    switch (IonType[i])
      {
      case C_atom :
	printf("C ");
	break;
      case H_atom :
	printf("H ");
	break;
      case CH_atom :
	printf("CH");
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
// vkv3: new potential for C60-electron
// Vpol is using fluctuating charges and dipoles
// 
/////////////////////////////////////////////////////////////////////////////////////////////
void Potential::SetupFullerElec(int nr, const double *rSites, const double *PotPara, const int *NoAtomsArray, const double *DipoleArray, const double *CenterArray, int *IonType)
{ 
  cout << "Enter SetupFullerElec\n";
  PotFlags.resize(1) ;
  PotFlags[0] =  1021 ;
  EsDampCsixty = PotPara[0];
  RepScaleCsixty = PotPara[1];
  GaussExpCsixty = PotPara[2];
  PolDampCsixty = PotPara[3];
  PolFlagCsixty = PotPara[4];
  nFullerenes =  PotPara[5] ;
  verbose = PotPara[6] ;
  ExtraNucQ = PotPara[7] ;
  RepFlagCsixty =  PotPara[8] ; 
  Rtol =  PotPara[10] ; 
  CarbonType = PotPara[11] ;
  AlphaComp = PotPara[12] ;

  AtomType  = new int[nr];

  for (int i; i < nr; i++)
  AtomType[i]=IonType[i];

  nAtomsArray.resize(nFullerenes) ;
  nReturnEnergies = 5;
  ReturnEnergies.resize(nReturnEnergies) ;
  for (int i = 0 ; i < nFullerenes ; ++i) {
    nAtomsArray[i] = NoAtomsArray[i] ; 
  }
  //  if (verbose > 0) {
    cout << "------------------------\n";
    cout << "FullerElec potential\n";
    cout << "PotPara 0: " << PotPara[0] << " (EsDamp) [Bohr]\n";
    cout << "PotPara 1: " << PotPara[1] << " (RepScale) [Hartree]\n";
    cout << "PotPara 2: " << PotPara[2] << " (GaussExp) [Bohr]\n";
    cout << "PotPara 3: " << PotPara[3] << " (PolDamp) [Bohr]\n";
    cout << "PotPara 4: " << PotPara[4] << " (PolFlag)=0 (no polarization) = 1 (ind+pol) =2 (only ind) =3 (induction with charge flow between fullerenes)\n";
    cout << "PotPara 5: " << PotPara[5] << " (Number of fullerene molecules) \n";
    cout << "PotPara 6: " << PotPara[6] << " (Verbosity) \n";
    cout << "PotPara 7: " << PotPara[7] << " (Extra Nuclear Charge for KT stabilization) \n";
    cout << "PotPara 8: " << PotPara[8] << " (Repulsive flag) = 0 (no repulsion) =1 (gaussian repulsion) \n";
    cout << "------------------------\n";
  //}
  SetSites(nr, rSites);
  SetPtDipoles( NoAtomsArray, DipoleArray, CenterArray) ; 
  SetAngleArray( NoAtomsArray,  CenterArray) ; 
  cout<<"start SetMoleculeFields "<<endl;
  SetMoleculeFields( NoAtomsArray) ;
  MolPol = new DistributedPolarizabilities[1];
// for the coronene
 if (CarbonType != 1) SetCoroneneGDMA() ;

  if (PolFlagCsixty == 1 || PolFlagCsixty == 2 || PolFlagCsixty == 8) {
    MolPol[0].InvA.resize((4*nr + nFullerenes)*(4*nr + nFullerenes));
    InvChargeDipolePol(nSites, &rSites[0], &MolPol[0].InvA[0], nFullerenes, &NoAtomsArray[0], CarbonType, IonType, AlphaComp);
   // for (int i=0;  i < (4*nr + nFullerenes)*(4*nr + nFullerenes); i++) 
   //    cout<<" MolPol[0].InvA[0] = "<< MolPol[0].InvA[i]<<endl;
    InterFullerPol = InterFullerPolPot(&EfieldOnAtoms[0]) ; 
  //  cout << "Intermolecular polarization is" << InterFullerPol << "au \n" ; 
    }
  if (PolFlagCsixty == 3 || PolFlagCsixty == 6 || PolFlagCsixty == 7) {
    cout << "Using new Polarization =3 that has intermolecular charge flow "<< PolFlagCsixty << endl ; 
    iVec NoAtomsArray_pol3 ; NoAtomsArray_pol3.resize(1) ; 
    NoAtomsArray_pol3[0] = nr ; 

    cout<<"(4*nr + 1)*(4*nr + 1) ="<<(4*nr + 1)*(4*nr + 1)<<endl;
    MolPol[0].InvA.resize((4*nr + 1)*(4*nr + 1)); 

    cout<<"(nr*(nr-1.0)/2.0)*(nr*(nr-1.0)/2.0) ="<<(nr*(nr-1.0)/2.0)*(nr*(nr-1.0)/2.0)<<endl;
 //   MolPol[0].InvA.resize((nr*(nr-1.0)/2.0)*(nr*(nr-1.0)/2.0)); 

   cout << "Start InvChargeDipolePol\n";

    InvChargeDipolePol(nSites, &rSites[0], &MolPol[0].InvA[0], 1, &NoAtomsArray_pol3[0],CarbonType, IonType, AlphaComp); 
  //  for (int i=0 ; i < (4*nr + 1)*(4*nr + 1) ; i++) 
   //    cout<<" MolPol[0].InvA[0] = "<< MolPol[0].InvA[i]<<endl;
    }

  // cout << "nr" << nr << endl;
  // int n4plus1 = nr*4+1 ; 
  // for (int i = 0; i < n4plus1*n4plus1; ++i){
  //   cout << "InvA" << i << "=" << MolPol[0].InvA[i] << endl;
  // }
  //nReturnEnergies = 5;
  //ReturnEnergies.resize(nReturnEnergies); 
  //  SetupMinMax();
  cout << "Leaving SetupFullerElec\n";
}

  // for now coding only distance damping method 
  // and electrostatics, charge-flow within this
void Potential::SetPtDipoles(const int *NoAtomsArray, const double *DipoleArray, const double *CenterArray)
{ 
  // do not damp these
  PtDip.resize(nSites*3) ;  
  //dVec origin ; origin.resize(3) ; 
  //origin[0] = 0.0 ;  origin[1] = 0.0 ;  origin[2] = 0.0 ;
  int current_atom = 0;
  int old_atom = 0;
  for (int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom = current_atom ; 
    current_atom += NoAtomsArray[iFuller] ;
    //cout<<"old_atom = "<<old_atom<<" current_atom = "<<current_atom<<endl;
    for (int i =old_atom ; i < current_atom; ++i){
      double RSiteCenter = sqrt((Site[i*3+0]-CenterArray[iFuller*3+0])*(Site[i*3+0]-CenterArray[iFuller*3+0]) +
                                (Site[i*3+1]-CenterArray[iFuller*3+1])*(Site[i*3+1]-CenterArray[iFuller*3+1]) +
                                (Site[i*3+2]-CenterArray[iFuller*3+2])*(Site[i*3+2]-CenterArray[iFuller*3+2])) ;
      PtDip[3*i+0] = DipoleArray[iFuller]*(Site[i*3+0] - CenterArray[iFuller*3+0])/RSiteCenter ; 
      PtDip[3*i+1] = DipoleArray[iFuller]*(Site[i*3+1] - CenterArray[iFuller*3+1])/RSiteCenter ;
      PtDip[3*i+2] = DipoleArray[iFuller]*(Site[i*3+2] - CenterArray[iFuller*3+2])/RSiteCenter ;
     // cout<<"PtDip ["<<3*i+0<<","<<3*i+1<<","<<3*i+2<<"]="<<PtDip[3*i+0]<<", "<<PtDip[3*i+1]<<", "<<PtDip[3*i+2]<<endl;
     // cout<<"CenterArray ["<<iFuller*3+0<<","<<iFuller*3+1<<","<<iFuller*3+2<<"]="<<CenterArray[iFuller*3+0]<<", "<<CenterArray[iFuller*3+1]<<", "<<CenterArray[iFuller*3+2]<<endl;

    }
  }
}

void Potential::SetAngleArray(const int *NoAtomsArray, const double *CenterArray)
{ 
  // do not damp these
  AngleArray.resize(nSites*3) ;  
  //dVec origin ; origin.resize(3) ; 
  //origin[0] = 0.0 ;  origin[1] = 0.0 ;  origin[2] = 0.0 ;
  int current_atom = 0;
  int old_atom = 0;
  for (int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom = current_atom ; 
    current_atom += NoAtomsArray[iFuller] ;
    //cout<<"old_atom = "<<old_atom<<" current_atom = "<<current_atom<<endl;
    for (int i =old_atom ; i < current_atom; ++i){
/*
      double RSiteCenter = sqrt((Site[i*3+0]-CenterArray[iFuller*3+0])*(Site[i*3+0]-CenterArray[iFuller*3+0]) +
                                (Site[i*3+1]-CenterArray[iFuller*3+1])*(Site[i*3+1]-CenterArray[iFuller*3+1]) +
                                (Site[i*3+2]-CenterArray[iFuller*3+2])*(Site[i*3+2]-CenterArray[iFuller*3+2])) ;
      AngleArray[3*i+0] = (Site[i*3+0] - CenterArray[iFuller*3+0])/RSiteCenter ; 
      AngleArray[3*i+1] = (Site[i*3+1] - CenterArray[iFuller*3+1])/RSiteCenter ;
      AngleArray[3*i+2] = (Site[i*3+2] - CenterArray[iFuller*3+2])/RSiteCenter ;
*/
      double RSiteCenter = sqrt((Site[i*3+2])*(Site[i*3+2])) ;
      AngleArray[3*i+0] = 0.0 ; 
      AngleArray[3*i+1] = 0.0 ;
      AngleArray[3*i+2] = (Site[i*3+2])/RSiteCenter ;



     // cout<<"AngleArray ["<<3*i+0<<","<<3*i+1<<","<<3*i+2<<"]="<<AngleArray[3*i+0]<<", "<<AngleArray[3*i+1]<<", "<<AngleArray[3*i+2]<<endl;
     // cout<<"CenterArray ["<<iFuller*3+0<<","<<iFuller*3+1<<","<<iFuller*3+2<<"]="<<CenterArray[iFuller*3+0]<<", "<<CenterArray[iFuller*3+1]<<", "<<CenterArray[iFuller*3+2]<<endl;

    }
  }
}

void Potential::SetMoleculeFields(const int *NoAtomsArray)
{ 
  // calculating the fields on the atoms of a given molecule from all other molecules
  int current_atom = 0;
  int old_atom = 0;
  double qpart[8];
  qpart[0]=-0.0235;
  qpart[1]= 0.0235;
  qpart[2]=-0.0235;
  qpart[3]= 0.0235;
  EfieldOnAtoms.resize(nSites*4+nFullerenes) ;
  for (int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom = current_atom ; 
    current_atom += NoAtomsArray[iFuller] ;
    cout<<"nfullerenes old_atom current_atom "<<nFullerenes<<" "<<old_atom<<" "<<current_atom<<endl;
    for (int i =old_atom ; i < current_atom; ++i){
      EfieldOnAtoms[nSites+i*3+0] = 0.0 ;
      EfieldOnAtoms[nSites+i*3+1] = 0.0 ; 
      EfieldOnAtoms[nSites+i*3+2] = 0.0 ;
      EfieldOnAtoms[i] = 0.0 ; 

      double ChargePen=-1.75;
      double NuCharge=4.0;
      double electrostatic = 0.0;
      double R_const = 0.68620399*1.889725989 ;
      double Rqq = R_const*sqrt(2.0)  ;
 //

// for the tsst of the electrostatic potential

      double space_x[60];
      double space_y[60];
      double space_z[60];

/*

char str[999];
double myvariable=0.0;
double temp=0.0;
int myint;
FILE * myfile;
myfile = fopen( "electro.data" , "r");
for(int kk =0; kk< 60; kk++) {
  fscanf(myfile,"%lf %lf %lf %lf \n",&temp, &space_x[kk],&space_y[kk],&space_z[kk]);
  cout<<"space="<<space_x[kk]<<" "<<space_y[kk]<<" "<<space_z[kk]<<endl;
}
*/


     for (int kk=0; kk< 60; ++kk){
     space_x[kk]	=  0.0;
     space_y[kk]	=  0.0;
     space_z[kk]	=  double(kk)*0.25;
    }

/*


     for (int ii = 0 ; ii < 60; ++ii){
      electrostatic = 0.0;

        for (int j = 0  ; j < nSites ; ++j){
          double Rx[i] = space_x[ii] - Site[j*3+0] ;
          double Ry[i] = space_y[ii] - Site[j*3+1] ;
          double Rz[i] = space_z[ii] - Site[j*3+2] ;
          double rSite   = sqrt(Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i]) ;
          double rSite2  = rSite*rSite ;
          double rSite3  = rSite2*rSite ;
          double MuDotR  = PtDip[3*j+0]*Rx[i] + PtDip[3*j+1]*Ry[i] + PtDip[3*j+2]*Rz[i] ;  
          //potential on atoms  due to dipoles
          electrostatic  += MuDotR/rSite3 ;
         // electrostatic  += MuDotR/rSite3 - NuCharge*exp(ChargePen*rSite)*erf(rSite/Rqq)/rSite;
        //  cout<<" MuDotR/rSite3 [ "<<j<<" ]= "<< MuDotR/rSite3<<endl;
        }
          cout<<" electrostatic potentail[ "<<ii<<" ]= "<<electrostatic<<endl;
     }

     exit(0);
*/

      for (int j = 0  ; j < nSites ; ++j){
        // avoids fields on a molecules from itself
       // cout << "i atom is "<< i << " j atom is " << j << "\n";
      //  if (i != j ) {
        if (j < old_atom | j >= current_atom) {
        //  cout << "i atom is "<< i << " j atom is " << j << "\n";
           Rx[i] = Site[i*3+0] - Site[j*3+0] ;
           Ry[i] = Site[i*3+1] - Site[j*3+1] ;
           Rz[i] = Site[i*3+2] - Site[j*3+2] ;
          double rSite   = sqrt(Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i]) ;
          double rSite2  = rSite*rSite ;
          double rSite3  = rSite2*rSite ;
          double MuDotR  = PtDip[3*j+0]*Rx[i] + PtDip[3*j+1]*Ry[i] + PtDip[3*j+2]*Rz[i] ;  
         //  cout<<"PtDip[3*j+2] = "<<PtDip[3*j+2]<<endl;
         //  cout<<"Rz[i] = "<<Rz[i]<<endl;
         //  cout<<"MuDotR = "<<MuDotR<<endl;
          // potential on atoms  due to dipoles
          EfieldOnAtoms[i]            += MuDotR/rSite3;
        // EfieldOnAtoms[i]            +=  -  NuCharge*exp(ChargePen*rSite)*erf(rSite/Rqq)/rSite;
         // EfieldOnAtoms[i]            -= (1.0-exp(ChargePen*rSite))*erf(rSite/Rqq)/rSite*NuCharge*exp(ChargePen*rSite);
        // EfieldOnAtoms[i]            -=  erf(rSite/Rqq)/rSite*NuCharge*exp(ChargePen*rSite);

          
          //fields on atoms due to dipoles
          EfieldOnAtoms[nSites+i*3+0] += (3.0*MuDotR*Rx[i]/rSite2 - PtDip[3*j+0])/rSite3 ; 
          EfieldOnAtoms[nSites+i*3+1] += (3.0*MuDotR*Ry[i]/rSite2 - PtDip[3*j+1])/rSite3 ; 
          EfieldOnAtoms[nSites+i*3+2] += (3.0*MuDotR*Rz[i]/rSite2 - PtDip[3*j+2])/rSite3 ; 
          //cout<<" EfieldOnAtoms[ "<<nSites+i*3+0<<" ]= "<<EfieldOnAtoms[nSites+i*3+0]<<endl;
          //cout<<" EfieldOnAtoms[ "<<nSites+i*3+1<<" ]= "<<EfieldOnAtoms[nSites+i*3+1]<<endl;
          //cout<<" EfieldOnAtoms[ "<<nSites+i*3+2<<" ]= "<<EfieldOnAtoms[nSites+i*3+2]<<endl;

        double fac, fac1, fac2, fac_pq ;
        fac = erf(rSite/Rqq)-2.0*rSite*exp(-1.0*(rSite/Rqq)*(rSite/Rqq))/(sqrt(PI)*Rqq) ;
        fac_pq = fac/(rSite*rSite*rSite) ;



         // EfieldOnAtoms[nSites+i*3+0] += Rx[i]*fac_pq*NuCharge*exp(ChargePen*rSite);
        //  EfieldOnAtoms[nSites+i*3+1] += Ry[i]*fac_pq*NuCharge*exp(ChargePen*rSite);
        //  EfieldOnAtoms[nSites+i*3+2] += Rz[i]*fac_pq*NuCharge*exp(ChargePen*rSite);


         // double Rz[i]_p = Site[i*3+2]+0.5 - Site[j*3+2] ;
        //  double Rz[i]_m = Site[i*3+2]-0.5 - Site[j*3+2] ;
          double rSite_p   = sqrt((Rx[i]+0.5)*(Rx[i]+0.5) + (Ry[i]+0.5)*(Ry[i]+0.5) + (Rz[i]+0.5)*(Rz[i]+0.5)) ;
          double rSite_m   = sqrt((Rx[i]-0.5)*(Rx[i]-0.5) + (Ry[i]-0.5)*(Ry[i]-0.5) + (Rz[i]-0.5)*(Rz[i]-0.5)) ;
          double rSite3_p  = rSite_p*rSite_p*rSite_p ;
          double rSite3_m  = rSite_m*rSite_m*rSite_m ;


        //  EfieldOnAtoms[i]            += qpart[j]/rSite_p - qpart[j]/rSite_m;
        //  EfieldOnAtoms[nSites+i*3+2] +=  qpart[j]*Rz[i]_p/rSite3_p -  qpart[j]*Rz[i]_m/rSite3_m; 
        //    cout<<"qpart["<<j<<"]="<<qpart[j]<<endl;
         // cout<<" EfieldOnAtoms[ "<< i <<" ]= "<<EfieldOnAtoms[i]<<endl;
        //  cout<<" EfieldOnAtoms[ "<<nSites+i*3+2<<" ]= "<<EfieldOnAtoms[nSites+i*3+2]<<endl;
        }

      }
          //cout<<" EfieldOnAtoms[ "<< i <<" ]= "<<EfieldOnAtoms[i]<<endl;


    }
  }
  if (verbose > 0){
    for (int j = 0  ; j < nSites ; ++j)
      {cout << "Potential and Efield on atom "<< j << " are  " << EfieldOnAtoms[j] << "  "
                                                               << EfieldOnAtoms[nSites+j*3+0] << " " 
                                                               << EfieldOnAtoms[nSites+j*3+1] << " " 
                                                               << EfieldOnAtoms[nSites+j*3+2] << "\n";} 
  }
}
   
void Potential::SetCoroneneGDMA()
{

   Q00.resize(nSites);
   Q20.resize(nSites);
   Q21c.resize(nSites);
   Q21s.resize(nSites);
   Q22c.resize(nSites);
   Q22s.resize(nSites);

   Q30.resize(nSites);
   Q31c.resize(nSites);
   Q31s.resize(nSites);
   Q32c.resize(nSites);
   Q32s.resize(nSites);
   Q33c.resize(nSites);
   Q33s.resize(nSites);

   Q40.resize(nSites);
   Q41c.resize(nSites);
   Q41s.resize(nSites);
   Q42c.resize(nSites);
   Q42s.resize(nSites);
   Q43c.resize(nSites);
   Q43s.resize(nSites);
   Q44c.resize(nSites);
   Q44s.resize(nSites);

  std::string word,word2;
  int rank ;
  ifstream GDMA;
  GDMA.open( "DBC.punch");

  for(int kk =0; kk< 4; ++kk) {
    std::getline(GDMA,word);
  }

  for(int kk =0; kk< nSites; ++kk) {
    std::getline(GDMA,word);

    GDMA >> Q00[kk];
    GDMA >> PtDip[3*kk+2] >>  PtDip[3*kk+0] >>  PtDip[3*kk+1];
    GDMA >> Q20[kk] >> Q21c[kk] >> Q21s[kk] >> Q22c[kk] >> Q22s[kk];

    //printf("%.10lf %.10lf %.10lf \n",coordX, coordY, coordZ);
    //printf("Q00 : %.10lf \n",Q00[kk]);
    //printf(" PtDip %.10lf %.10lf %.10lf \n", PtDip[3*kk+0], PtDip[3*kk+1], PtDip[3*kk+2]);
    //printf(" Q20 ~~ %.10lf %.10lf %.10lf %.10lf %.10lf\n", Q20[kk],Q21c[kk],Q21s[kk],Q22c[kk],Q22s[kk]);

    GDMA >> Q30[kk] >> Q31c[kk] >> Q31s[kk] >> Q32c[kk] >> Q32s[kk];
    GDMA >> Q33c[kk] >> Q33s[kk];
    GDMA >> Q40[kk] >> Q41c[kk] >> Q41s[kk] >> Q42c[kk] >> Q42s[kk];
    GDMA >> Q43c[kk] >> Q43s[kk] >> Q44c[kk] >> Q44s[kk];

/*
    printf(" Q30 ~~ %.10lf %.10lf %.10lf %.10lf %.10lf\n", Q30[kk],Q31c[kk],Q31s[kk],Q32c[kk],Q32s[kk]);
    printf(" Q33 ~~ %.10lf %.10lf \n", Q33c[kk],Q33s[kk]);
    printf(" Q40 ~~ %.10lf %.10lf %.10lf %.10lf %.10lf\n", Q40[kk],Q41c[kk],Q41s[kk],Q42c[kk],Q42s[kk]);
    printf(" Q43 ~~ %.10lf %.10lf %.10lf %.10lf\n", Q43c[kk],Q43s[kk], Q44c[kk],Q44s[kk]);
*/

    for(int i =0; i< 5; ++i) {
      std::getline(GDMA,word);
    } 

  }
}



double Potential::EvaluateCsixty(const double *x)
{

//  progress_timer tmr("EvaluateCsixty:", verbose);
  double Vtotal = 0.0 ; 
  double VRep = EvaluateCsixtyRepPot(x, GaussExpCsixty, RepScaleCsixty) ; 
  double Vel;

  if (CarbonType == 1 ) Vel = EvaluateCsixtyElecPot(x);
  else  Vel = EvaluateColoreneElecPot(x);

  double VPol =  EvaluateCsixtyPolPot(x) ;
 // cout<<"out VPol colorene"<<endl;
  ReturnEnergies[0] = Vel;
  ReturnEnergies[2] = VRep;
  ReturnEnergies[3] = VPol;
  if (verbose > 0)
    CheckMinMax(&Vel);
  if (verbose > 15)
    printf("%10.6f %10.6f %10.6f %14.6e\n", x[0], x[1], x[2], Vel);

  if (PolFlagCsixty == 7 || PolFlagCsixty == 8) Vtotal = VPol;
  else Vtotal = Vel + VPol + VRep ; 
 // else Vtotal = VPol ; 
  


  //Vtotal = Vel + VRep ; 
  if (verbose > 0)
    printf("%10.6f %10.6f %10.6f %14.6e %14.6e %14.6e %14.6e \n", x[0], x[1], x[2], Vel, VRep, VPol, Vtotal);
  return Vtotal ;
}

double Potential::EvaluateColoreneElecPot(const double *x)
{ 

// progress_timer tmr("EvaluateColoreneElecPot:", verbose);

  double Velec = 0.0 ;
  double ElsPot= 0.0 ;
  double coordX, coordY, coordZ;


  for (int i =0 ; i < nSites; ++i){
    // effective-r damping
    double DampParameter = EsDampCsixty ;
    double Reff = R[i];
    double gij = Rminus3[i] ; 
    double g5 = 1/(Reff*Reff*Reff*Reff*Reff) ; 
    double g7 = 1/(Reff*Reff*Reff*Reff*Reff*Reff*Reff) ; 


    double Qxx, Qyy,Qzz;
    double Qxy,Qyx,Qxz,Qzx,Qyz,Qzy;

    double Qxxx,Qxxy,Qyxx,Qxyx,Qxyy,Qyxy,Qyyx,Qyyy,Qxxz,Qzxx,Qxzx,Qxyz,Qyzx,Qzxy;
    double Qyxz,Qxzy,Qzyx,Qyyz,Qyzy,Qzyy,Qxzz,Qzxz,Qzzx,Qyzz,Qzyz,Qzzy,Qzzz;



    Qxx=-0.5*Q20[i]+0.5*sqrt(3.0)*Q22c[i];
    Qyy=-0.5*Q20[i]-0.5*sqrt(3.0)*Q22c[i];
    Qzz=Q20[i];
    Qxy=0.5*sqrt(3.0)*Q22s[i];
    Qyx=Qxy;
    Qxz=0.5*sqrt(3.0)*Q21c[i];
    Qzx=Qzx;
    Qyz=0.5*sqrt(3.0)*Q21s[i];
    Qzy=Qyz;

    Qxxx=sqrt(5.0/8.0)*Q33c[i]-sqrt(3.0/8.0)*Q31c[i];
    Qxxy=sqrt(5.0/8.0)*Q33s[i]-sqrt(1.0/24.0)*Q31s[i];
    Qyxx=Qxxy;
    Qxyx=Qxxy;
    Qxyy=sqrt(5.0/8.0)*Q33c[i]-sqrt(1.0/24.0)*Q31c[i];
    Qyxy=Qxyy;
    Qyyx=Qxyy;
    Qyyy=sqrt(5.0/8.0)*Q33s[i]-sqrt(1.0/24.0)*Q31s[i];
    Qxxz=sqrt(5.0/12.0)*Q32c[i]-sqrt(1.0/2.0)*Q30[i];
    Qzxx=Qxxz;
    Qxzx=Qxxz;
    Qxyz=sqrt(5.0/12.0)*Q32s[i];
    Qyzx=Qxyz;
    Qzxy=Qxyz;
    Qyxz=Qxyz;
    Qxzy=Qxyz;
    Qzyx=Qxyz;
    Qyyz=-sqrt(5.0/12.0)*Q32c[i]-sqrt(1.0/2.0)*Q30[i];
    Qyzy=Qyyz;
    Qzyy=Qyyz;
    Qxzz=sqrt(2.0/3.0)*Q31c[i];
    Qzxz=Qxzz;
    Qzzx=Qxzz;
    Qyzz=sqrt(2.0/3.0)*Q31s[i];
    Qzyz=Qyzz;
    Qzzy=Qyzz;
    Qzzz=Q30[i];

   double MuDotR=0.0;
   double  rQr=0.0;
   double  rrQr=0.0;

 //   cout<<"rrQr = "<<rrQr<<endl;

     MuDotR  = PtDip[3*i+0]*Rx[i] + PtDip[3*i+1]*Ry[i] +PtDip[3*i+2]*Rz[i];


     rQr=      Rx[i]*Qxx*Rx[i] + Rx[i]*Qxy*Ry[i] +Rx[i]*Qxz*Rz[i] \
             + Ry[i]*Qyx*Rx[i] + Ry[i]*Qyy*Ry[i] +Ry[i]*Qyz*Rz[i] \
             + Rz[i]*Qzx*Rx[i] + Rz[i]*Qzy*Ry[i] +Rz[i]*Qzz*Rz[i];


     rrQr=   Rx[i]*Rx[i]*Qxxx*Rx[i] + Rx[i]*Rx[i]*Qxxy*Ry[i] +Rx[i]*Rx[i]*Qxxz*Rz[i]\
           + Rx[i]*Ry[i]*Qxyx*Rx[i] + Rx[i]*Ry[i]*Qxyy*Ry[i] +Rx[i]*Ry[i]*Qxyz*Rz[i]\
           + Rx[i]*Rz[i]*Qxzx*Rx[i] + Rx[i]*Rz[i]*Qxzy*Ry[i] +Rx[i]*Rz[i]*Qxzz*Rz[i]\
           + Ry[i]*Rx[i]*Qyxx*Rx[i] + Ry[i]*Rx[i]*Qyxy*Ry[i] +Ry[i]*Rx[i]*Qyxz*Rz[i]\
           + Ry[i]*Ry[i]*Qyyx*Rx[i] + Ry[i]*Ry[i]*Qyyy*Ry[i] +Ry[i]*Ry[i]*Qyyz*Rz[i]\
           + Ry[i]*Rz[i]*Qyzx*Rx[i] + Ry[i]*Rz[i]*Qyzy*Ry[i] +Ry[i]*Rz[i]*Qyzz*Rz[i]\
           + Rz[i]*Rx[i]*Qzxx*Rx[i] + Rz[i]*Rx[i]*Qzxy*Ry[i] +Rz[i]*Rx[i]*Qzxz*Rz[i]\
           + Rz[i]*Ry[i]*Qzyx*Rx[i] + Rz[i]*Ry[i]*Qzyy*Ry[i] +Rz[i]*Ry[i]*Qzyz*Rz[i]\
           + Rz[i]*Rz[i]*Qzzx*Rx[i] + Rz[i]*Rz[i]*Qzzy*Ry[i] +Rz[i]*Rz[i]*Qzzz*Rz[i];


    if (Reff < DampParameter) {
      double ror0 = Reff / DampParameter;
      Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
      // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
 
      gij = 1.0 / (Reff * Reff * Reff);
      g5 = 1.0 / (Reff * Reff * Reff*Reff*Reff);
      g7 = 1.0 / (Reff * Reff * Reff*Reff*Reff*Reff*Reff);

    }

    // v_dipole_charge_interaction = dot(dipole[i],R[i])/R^3
    //Velec += (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*
    //          gij ;
 //   cout<<" Q00["<<i<<"] and Reff = "<<Q00[i]<<" "<<Reff<<endl; 
 //   cout<<"dipole: "<<PtDip[3*i+0]<<" "<< PtDip[3*i+1] <<" "<< PtDip[3*i+2]<<endl;  
 //   cout<<"Reff: "<<Rx[i]<<" "<< Ry[i] <<" "<< Rz[i]<<endl;  
   
   Velec+=-1.0*Q00[i]/Reff ;
   Velec+= -(Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*gij  ; 
   Velec+= -rQr*g5 ; 
   Velec+= - rrQr*g7 ; 

  //  Velec += -1.0*ExtraNucQ/Reff + (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*gij ;
    
  }
 // cout<<" :  Velec = "<<Velec*27.2114<<endl;
  return Velec ;
   
}


double Potential::EvaluateCsixtyElecPot(const double *x)
{ 

// progress_timer tmr("EvaluateCsixtyElecPot:", verbose);

  double Velec = 0.0 ;
  for (int i =0 ; i < nSites; ++i){
    // effective-r damping
    double DampParameter = EsDampCsixty ;
    double Reff = R[i];
    double gij = Rminus3[i] ; 
    if (Reff < DampParameter) {
      double ror0 = Reff / DampParameter;
      Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
      // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
      gij = 1.0 / (Reff * Reff * Reff);
    }

    // v_dipole_charge_interaction = dot(dipole[i],R[i])/R^3
    //Velec += (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*
    //          gij ;
   // cout<<" ExtranucQ = "<<ExtraNucQ<<endl; 
    Velec += -1.0*ExtraNucQ/Reff + (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*gij ;
    
  }
  return Velec ;
   
}

double Potential::EvaluateCsixtyRepPot(const double *x, double GaussExp, double RepScale)
{
// progress_timer tmr("EvaluateCsixtyRepPot", verbose);
  double Vrep = 0.0 ;
 // cout<<"GaussExp = "<<GaussExp<<"  RepScale = "<<RepScale<<endl;
  switch (RepFlagCsixty)
    {
      case 0:
        Vrep = 0.0 ; //no repulsion
      break;
      case 1:
        for (int i = 0; i < nSites; ++i) {
           if (AtomType[i] == 2){
            // RepScale = 0.284244;  
           //  GaussExp = 1.64121;
             Vrep +=  RepScale*exp(-GaussExp * R2[i]); //gaussians
           }
           else 
             Vrep +=  RepScale*exp(-GaussExp * R2[i]); //gaussians
        }
      break ;
      default:
      printf("EvaluateCsixtyRep: This should never happen\n");
      exit(42);

    }
  return Vrep ; 
}


double Potential::EvaluateCsixtyPolPot(const double *x )
{ 

// progress_timer tmr("EvaluateCsixtyPolPot:", verbose);
  dVec Efield; Efield.resize(nSites*4+nFullerenes);
  dVec EfieldNet_new ; EfieldNet_new.resize(nSites*4+nFullerenes);
  dVec mu;     mu.resize(nSites*4+nFullerenes);
  double Vpol = 0.0 ;
  double dzero = 0.0;
  double done = 1.0;
  int one = 1;
  int n = nSites*4+nFullerenes ;


  dVec Efield_pol3 ; Efield_pol3.resize(nSites*4+1);
  dVec EfieldNet_pol3 ; EfieldNet_pol3.resize(nSites*4+1);
  dVec mu_pol3;     mu_pol3.resize(nSites*4+1);
  int n_pol3 = nSites*4+1 ;

  int ncomb=nSites*(nSites-1.0)/2.0;
  dVec Efield_comb ; Efield_comb.resize(n);
  dVec SijTotal ; SijTotal.resize(nSites*4+1 );
  dVec SijEach ; SijEach.resize(nSites*4+1 );
  dVec mu_comb;     mu_comb.resize(ncomb);



//  cout << "PolFlagCsixty " <<  PolFlagCsixty << "\n" ;
//  cout << "nSites " <<  nSites << "\n" ;
//  cout << "nFullerenes " <<  nFullerenes << "\n" ;
  
  switch (PolFlagCsixty)
    {
      case 0:
        Vpol = 0.0;
      break;
      case 1:
        // computing polarization and induction
        for (int i =0 ; i < nSites; ++i){
          double Reff = R[i];
          double Reffminus3 = Rminus3[i] ;
          if (Reff < PolDampCsixty) {
            double ror0 = Reff / PolDampCsixty;
            Reff = PolDampCsixty * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
            // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
            Reffminus3 = 1.0 / (Reff * Reff * Reff);
          }
          if (AtomType[i] == 3 ) Efield[i] = -1.0/Reff - 0.2451537;
          else if(AtomType[i] == 2 ) Efield[i] = -1.0/Reff - 0.21150515;
          else if(AtomType[i] == 1 ) Efield[i] = -1.0/Reff;

          Efield[nSites+3*i+0] = Reffminus3 * Rx[i];
          Efield[nSites+3*i+1] = Reffminus3 * Ry[i];
          Efield[nSites+3*i+2] = Reffminus3 * Rz[i];


        }
        // adding fields from other fullerenes
        // make sure fields are ADDED

     //   Efield[nSites] = 0.0000;
     //   n=nSites+1;

        for (int i =0 ; i < n; ++i){
         // cout<<"Efield["<<i<<"]="<<Efield[i]<<" and  EfieldOnAtoms[i]: "<< EfieldOnAtoms[i]<<endl; 
          //cout<<"EfieldOnAtoms[i]: "<< EfieldOnAtoms[i]<<endl; 
          Efield[i] -= EfieldOnAtoms[i]; 
          // EfieldNet_pol3[i] = 0.5*Efield[i]- EfieldOnAtoms[i]; 
          // cout<<"Efield["<<i<<"]="<<EfieldNet_pol3[i]<<endl; 
        }
        //cout<<"InterFullerPol = "<<InterFullerPol<<endl;
        //exit(0);
        // compute dipoles induced by Efield
       // dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
        dsymv("U", &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
        Vpol = -0.5*ddot(&n, &Efield[0], &one, &mu[0], &one) - InterFullerPol;
        //Vpol = -ddot(&n, &Efield[0], &one, &mu[0], &one) ;
       // for (int i =0 ; i < n ; ++i){
      //     cout<<"mu[ "<<i<<" ] = "<<mu[i]<<endl;
      //  }
      // exit(0);
       // cout<<"Vpol = "<<Vpol<<endl;
        if (verbose > 1) {
          int current_atom = 0;
          int old_atom = 0;
          for (int i = 0 ; i < nFullerenes ; ++i){
            old_atom = current_atom ;
            current_atom += nAtomsArray[i] ;
            double net_charge = 0.0 ;
            for (int j =old_atom ; j < current_atom; ++j){
              net_charge += mu[j] ;
            }
            cout << "net charge on fullerene " <<  i << " is " << net_charge << "\n" ;
          }    
        }
        // Vpol = -0.5*ddot(&n, &Efield[0], &one, &mu[0], &one) ;
      break;
      case 2:
        // computing only induction
        for (int i =0 ; i < nSites; ++i){
          double Reff = R[i];
          double Reffminus3 = Rminus3[i] ;
          if (Reff < PolDampCsixty) {
            double ror0 = Reff / PolDampCsixty;
            Reff = PolDampCsixty * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
            // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
            Reffminus3 = 1.0 / (Reff * Reff * Reff);
          }
          Efield[i] = -1.0/Reff ;
          Efield[nSites+3*i+0] = Reffminus3 * Rx[i];
          Efield[nSites+3*i+1] = Reffminus3 * Ry[i];
          Efield[nSites+3*i+2] = Reffminus3 * Rz[i];


        }
        //dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldOnAtoms[0], &one, &dzero, &mu[0], &one);
        dsymv("L", &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldOnAtoms[0], &one, &dzero, &mu[0], &one);
        Vpol = -ddot(&n, &Efield[0], &one, &mu[0], &one) ; // not using 0.5 as nothing is double counted
      break ;
      case 3:
      case 6:
        // induction with intermolecular charge flow
        // computing only induction
        for (int i =0 ; i < nSites; ++i){
          double Reff = R[i];
          double Reffminus3 = Rminus3[i] ;
          if (Reff < PolDampCsixty) {
            double ror0 = Reff / PolDampCsixty;
            Reff = PolDampCsixty * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
            // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
            Reffminus3 = 1.0 / (Reff * Reff * Reff);
          }
          Efield_pol3[i] = -1.0/Reff   ;
          Efield_pol3[nSites+3*i+0] = Reffminus3 * Rx[i];
          Efield_pol3[nSites+3*i+1] = Reffminus3 * Ry[i];
          Efield_pol3[nSites+3*i+2] = Reffminus3 * Rz[i];
        }
        for (int i =0 ; i < n_pol3 ; ++i){
          //EfieldNet_pol3[i] = 0.5*Efield_pol3[i] + EfieldOnAtoms[i];
        //  cout<<"Efield_pol3["<<i<<"]="<<Efield_pol3[i]<<" and  EfieldOnAtoms[i]: "<< EfieldOnAtoms[i]<<endl; 
          EfieldNet_pol3[i] =  0.5*Efield_pol3[i] - EfieldOnAtoms[i] ;
          //cout<<"EfieldOnAtoms[i]: "<< EfieldOnAtoms[i]<<endl; 
        //  cout<<"Efield["<<i<<"]="<<EfieldNet_pol3[i]<<endl; 
        }
       // exit(0);
       // cout<<" n_pol3 = "<<n_pol3<<endl;
        //dgemv("N", &n_pol3, &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &EfieldNet_pol3[0], &one, &dzero, &mu_pol3[0], &one);
        dsymv("L", &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &EfieldNet_pol3[0], &one, &dzero, &mu_pol3[0], &one);
       //  for (int i =0 ; i < n_pol3 ; ++i){
      //     cout<<"mu[ "<<i<<" ] = "<<2.0*mu_pol3[i]<<endl;
      //  }

        /*
        for (int i =0 ; i < n_pol3 ; ++i){
           cout<<"mu_pol3 = "<<mu_pol3[i]<<endl;
        }
        */
         
        Vpol = -ddot(&n_pol3, &Efield_pol3[0], &one, &mu_pol3[0], &one) ; // not using 0.5 as nothing is double counted
       // cout<<"Vpol = "<<Vpol<<endl;
      break;
      case 7: {
      // cout<<"here is the elecric filed ----- "<<endl;
 


      double Ex=0.000;
      double Ey=0.000;
      double Ez=0.001;

        for (int i =0 ; i < n_pol3 ; ++i) {
          Efield_comb[i] =  0.00;
          SijTotal[i] =  0.00;
          SijEach[i] =  0.00;
        }
/*
      int icount=0;
      for (int i =0 ; i < nSites; ++i){
         for (int j = i+1 ; j < nSites; ++j){
           Efield_comb[icount] = 0.5*(Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez) ; 
           Efield_comb[icount] -= 0.5*(Rx[j]*Ex+Ry[j]*Ey+ Rz[j]*Ez) ; 
          // Efield_comb[icount] = 0.5*(Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez) - EfieldOnAtoms[i] ; 
         //  Efield_comb[icount] -= 0.5*(Rx[j]*Ex+Ry[j]*Ey+ Rz[j]*Ez) - EfieldOnAtoms[j]; 
           icount++;
         }
      }


*/


        for (int i =0 ; i < nSites; ++i){
          Efield_pol3[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez; 
         // cout<<"Rz["<<i<<"]="<<Rz[i]<<endl;
          //Efield_pol3[i] = Site[3*i+0]*Ex+Ry[3*i+1]*Ey+ Rz[3*i+2]*Ez; 
          //
          //
     //   if (AtomType[i] == 3 ) Efield_pol3[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 0.2451537;
     //   else if(AtomType[i] == 2 ) Efield_pol3[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 0.21150515;
     //    else if(AtomType[i] == 1 ) Efield_pol3[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez ;

          Efield_pol3[nSites+3*i+0] = Ex;
          Efield_pol3[nSites+3*i+1] = Ey;
          Efield_pol3[nSites+3*i+2] = Ez;

          Efield_comb[nSites+3*i+0] = Ex;
          Efield_comb[nSites+3*i+1] = Ey;
          Efield_comb[nSites+3*i+2] = Ez;
          cout<<"Efield_pol3["<<i<<"]="<<Efield_pol3[i]<<endl; 

         }



     double Sij_p;

     double c1, a1, x1, y1, z1;
     double c2, a2, x2, y2, z2;
     double cosx, cosy, cosz;
     int l1, m1, n1, l2, m2, n2 ;

     c1=1;
     c2=1;
     a1=0.9554;
     a2=0.9554;

//     a1=1.46;
//     a2=1.46;

     x1= 0.0;
     y1= 0.0;
     z1= 0.0;

     x2= 0.0;
     y2= 0.0;
     z2= 0.0;

     l1=0;
     m1=0;
     n1=1;

     l2=0;
     m2=0;
     n2=1;

   for ( int i=0; i< 100 ; ++i){
      x1 = double(i*0.1);
    // cout<<"x1="<<x1<<endl;
     Sij_p=CalculateOverlapIntegral_S(c1, a1, l1, m1, n1, x1, y1, z1, c2, a2, l2, m2, n2, x2, y2, z2);
    cout<<"Slator Sij_p ["<<x1<<"]="<<Sij_p<<endl;
   }


     int Min1=0;
     int Min2=0;
     int Min3=0;
     int CurrentI=61;
     for (int i = 0; i< nSites ; i++){
        find_closed_aoms(i, Min1, Min2, Min3);
     //   cout<<i<<" closed to "<<Min1<<" "<<Min2<<" "<<Min3<<endl;
     }

  int halfnSites=nSites/2;

  halfnSites=halfnSites-6; //this is for only butadiene

  char jobz = 'V';
  char uplo = 'U';
  int lwork = -1;
  double wkopt;


// this is for vibrational conversion from h2o to hod : Tae Hoon Choi
   halfnSites=9;
  double* eigendata= new double[45];


  double* eigenvalues= new double[halfnSites*halfnSites];
  double* eigenvectors= new double[halfnSites*halfnSites];



  double* OccOrbital= new double[30*nSites];
  double* Overlapping= new double[nSites];
  double* VirtOrbital= new double[30*nSites];
  double* Sign= new double[30*30*nSites];

  double* wr= new double[halfnSites];
  double* wi= new double[halfnSites];
  double* vl= new double[halfnSites*2];
  double* vr= new double[halfnSites*2];







  //double alpha=-0.00389;
  //double beta=-0.3319;
 // double alpha=-0.09;
//  double beta=-0.318;

// double alpha=-0.137;
// double beta=-0.233;

  double alpha=0;
  double beta=-1;


  double ratio=2.64762/2.75508;

  for (int i =0 ; i < halfnSites; ++i){
    find_closed_aoms(i, Min1, Min2, Min3);
    eigenvectors[i*halfnSites +i]= alpha;
    eigenvectors[i*halfnSites+Min1]= beta;
   if (i ==0 || i==3)
    eigenvectors[i*halfnSites +Min2]= 0.0;
   else
    eigenvectors[i*halfnSites +Min2]= beta;
 //   eigenvectors[i*halfnSites +Min3]= beta*ratio;
  } 
 
eigendata[0]=	0.64681005	;
eigendata[1]=	-0.075984283	;
eigendata[2]=	0.47753783	;
eigendata[3]=	-1.22606E-16	;
eigendata[4]=	-4.54301E-17	;
eigendata[5]=	3.08801E-32	;
eigendata[6]=	-0.4773874	;
eigendata[7]=	-0.1034917	;
eigendata[8]=	1.10866E-16	;
eigendata[9]=	0.50146701	;
eigendata[10]=	-0.16354941	;
eigendata[11]=	-0.084786534	;
eigendata[12]=	4.42804E-17	;
eigendata[13]=	0.11981771	;
eigendata[14]=	0.078341811	;
eigendata[15]=	1.18534E-16	;
eigendata[16]=	3.19954E-17	;
eigendata[17]=	-2.83318E-32	;
eigendata[18]=	-1.17876E-16	;
eigendata[19]=	-3.4512E-17	;
eigendata[20]=	2.85185E-32	;
eigendata[21]=	-0.16942264	;
eigendata[22]=	0.17947599	;
eigendata[23]=	1.17402E-17	;
eigendata[24]=	-0.024079608	;
eigendata[25]=	0.043731695	;
eigendata[26]=	-6.58169E-19	;
eigendata[27]=	0.19350225	;
eigendata[28]=	0.23953369	;
eigendata[29]=	-0.39275129	;
eigendata[30]=	1.14969E-18	;
eigendata[31]=	-0.01632601	;
eigendata[32]=	0.006444723	;
eigendata[33]=	2.51668E-18	;
eigendata[34]=	-0.22320768	;
eigendata[35]=	0.38630657	;
eigendata[36]=	4.07195E-18	;
eigendata[37]=	1.34347E-17	;
eigendata[38]=	-2.54831E-33	;
eigendata[39]=	7.01012E-18	;
eigendata[40]=	-9.76837E-18	;
eigendata[41]=	-1.86703E-34	;
eigendata[42]=	-1.10821E-17	;
eigendata[43]=	-3.66636E-18	;
eigendata[44]=	2.73501E-33	;

/*
eigendata[0]=	3.47987E-32	;
eigendata[1]=	-7.5594E-17	;
eigendata[2]=	0.16382575	;
eigendata[3]=	-4.54121E-18	;
eigendata[4]=	1.42815E-14	;
eigendata[5]=	-0.24951615	;
eigendata[6]=	-1.8028E-32	;
eigendata[7]=	3.89462E-17	;
eigendata[8]=	-3.1342E-18	;
eigendata[9]=	1.96017E-32	;
eigendata[10]=	3.80102E-17	;
eigendata[11]=	-0.081912875	;
eigendata[12]=	0.011713178	;
eigendata[13]=	-4.18264E-17	;
eigendata[14]=	0.089168978	;
eigendata[15]=	-2.68661E-17	;
eigendata[16]=	0.063144455	;
eigendata[17]=	0.12475808	;
eigendata[18]=	1.65784E-17	;
eigendata[19]=	-0.037428817	;
eigendata[20]=	-0.038040848	;
eigendata[21]=	-1.67707E-32	;
eigendata[22]=	3.66478E-17	;
eigendata[23]=	7.67541E-18	;
eigendata[24]=	-1.57367E-33	;
eigendata[25]=	3.81621E-18	;
eigendata[26]=	1.02877E-17	;
eigendata[27]=	1.83443E-32	;
eigendata[28]=	3.75838E-17	;
eigendata[29]=	-0.081912875	;
eigendata[30]=	-0.011713178	;
eigendata[31]=	2.88015E-18	;
eigendata[32]=	-0.007256104	;
eigendata[33]=	-0.025715638	;
eigendata[34]=	-4.0464E-17	;
eigendata[35]=	0.089168978	;
eigendata[36]=	3.14073E-17	;
eigendata[37]=	-0.063144455	;
eigendata[38]=	0.12475808	;
eigendata[39]=	-1.34442E-17	;
eigendata[40]=	0.025715638	;
eigendata[41]=	-0.086717227	;
eigendata[42]=	-1.79631E-17	;
eigendata[43]=	0.037428817	;
eigendata[44]=	-0.038040848	;
*/

/*
 
eigendata[0]=	1.44062E-31	;
eigendata[1]=	-3.11887E-16	;
eigendata[2]=	0.67591423	;
eigendata[3]=	8.16151E-18	;
eigendata[4]=	2.13573E-15	;
eigendata[5]=	0.44843273	;
eigendata[6]=	-7.59025E-32	;
eigendata[7]=	1.60685E-16	;
eigendata[8]=	-9.65811E-17	;
eigendata[9]=	8.21432E-32	;
eigendata[10]=	1.59592E-16	;
eigendata[11]=	-0.33795712	;
eigendata[12]=	0.20046479	;
eigendata[13]=	-1.73671E-16	;
eigendata[14]=	0.36728548	;
eigendata[15]=	-1.24293E-16	;
eigendata[16]=	0.26052215	;
eigendata[17]=	-0.22421637	;
eigendata[18]=	1.10224E-16	;
eigendata[19]=	-0.23049347	;
eigendata[20]=	0.21252298	;
eigendata[21]=	-6.81597E-32	;
eigendata[22]=	1.51202E-16	;
eigendata[23]=	8.84196E-17	;
eigendata[24]=	-6.24064E-33	;
eigendata[25]=	1.40795E-17	;
eigendata[26]=	1.40689E-17	;
eigendata[27]=	7.44004E-32	;
eigendata[28]=	1.52295E-16	;
eigendata[29]=	-0.33795712	;
eigendata[30]=	-0.20046479	;
eigendata[31]=	1.29864E-17	;
eigendata[32]=	-0.029328367	;
eigendata[33]=	-0.030028678	;
eigendata[34]=	-1.65281E-16	;
eigendata[35]=	0.36728548	;
eigendata[36]=	1.16132E-16	;
eigendata[37]=	-0.26052215	;
eigendata[38]=	-0.22421637	;
eigendata[39]=	-1.36433E-17	;
eigendata[40]=	0.030028678	;
eigendata[41]=	0.01169339	;
eigendata[42]=	-1.02489E-16	;
eigendata[43]=	0.23049347	;
eigendata[44]=	0.21252298	;
*/



  double mass1=16.0;
  double mass2=2.0;
  double mass3=1.0;
  double massi;
  double massj;
  int ictest=0;
  

  for (int i =0 ; i < 9; ++i){
     for (int j = 0 ; j < i+1; ++j){
       if (i < 3) massi=mass1;
       else if ( i >=3 && i < 6) massi=mass2;
       else massi=mass3;
       if (j < 3) massj=mass1;
       else if ( j >=3 && j < 6) massj=mass2;
       else massj=mass3;

   
       eigenvectors[i*9 +j]=eigendata[ictest]/pow(massi,0.5)/pow(massj,0.5) ;
       eigenvectors[j*9 +i]=eigendata[ictest]/pow(massi,0.5)/pow(massj,0.5) ;
      // cout<<"eigendata[ictest]="<<eigendata[ictest]<<endl;
      // cout<<"eigenvectors[i*9 +j]="<<eigenvectors[i*9 +j]<<endl;
       cout<<i<<" "<<j<<" massi massj= "<<massi<<" "<<massj<<endl;
       ictest++;
     }
  }



//  for (int i =0 ; i < halfnSites; ++i){
//     for (int j =0 ; j < halfnSites; ++j){
//    cout<<"symmetric check "<<i<<" "<<j<<"  ="<< eigenvectors[i*halfnSites+j]<<" "<<eigenvectors[i*halfnSites+j] <<endl;
//    }
//  }

   for (int i =0 ; i < halfnSites; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          cout<<eigenvectors[i*halfnSites+j]<<" ";
       }
       cout<<endl;
   }

  

  int info=0;
  dsyev(&jobz, &uplo, &halfnSites, eigenvectors, &halfnSites, eigenvalues, &wkopt, &lwork, &info);

 cout<<"wkopt="<<wkopt<<endl;
  lwork=(int)wkopt;



//  lwork=12*nSites;

  double* work;
  work = (double*)malloc( lwork*sizeof(double) );
  dsyev(&jobz, &uplo, &halfnSites, eigenvectors, &halfnSites, eigenvalues, work, &lwork, &info);

  if (info != 0) {
    cout << "DMDCBase::CallDSYEV: DSYEV returned with info = " << info << endl;
    exit(1);
  }


   for (int i =0 ; i < halfnSites; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          cout<<eigenvectors[i*halfnSites+j]<<" ";
       }
       cout<<endl;
   }

   for (int i =0 ; i < halfnSites*halfnSites; ++i)
        cout<<"eigenvalues = "<<eigenvalues[i]<<endl;

     
   exit(1);
      
//   for (int i =0 ; i < halfnSites; ++i)
//      cout<<"eigenvalues = "<<eigenvalues[i]<<endl;
    //  cout<<"eigenvalues = "<<wr[i]<<" "<<wi[i]<<endl;

  

/*
   for (int i =30 ; i < 60; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          VirtOrbital[(i-30)*nSites+j]=eigenvectors[i*halfnSites+j];
       }
   }
   for (int i =0 ; i < 30; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          OccOrbital[(i-0)*nSites+j]=eigenvectors[i*halfnSites+j];
       }
   }
*/


   for (int i =2 ; i < 4; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          VirtOrbital[(i-2)*nSites+j]=eigenvectors[i*halfnSites+j];
       }
   }
   for (int i =0 ; i < 2; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          OccOrbital[(i-0)*nSites+j]=eigenvectors[i*halfnSites+j];
       }
   }



   for (int i =0 ; i < halfnSites; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          cout<<"eigenvector["<<i*halfnSites+j<<"]="<<eigenvectors[i*halfnSites+j]<<endl;
       }
       cout<<endl;
   }

   for (int i =0 ; i < halfnSites*halfnSites; ++i)
      eigenvectors[i]=0.0;

 //  halfnSites=nSites/2;  // this is for only butadiene

  for (int i =0 ; i < halfnSites; ++i){
    find_closed_aoms(i+nSites/2, Min1, Min2, Min3);
//    cout<<"min: "<<Min1<<" "<<Min2<<" "<<Min3<<endl;


    eigenvectors[i*halfnSites +i]= alpha;
    eigenvectors[i*halfnSites+Min1-nSites/2]= beta;
//    cout<<"i*halfnSites+Min1-nSites/2 ="<<i*halfnSites+Min1-nSites/2<<endl;
//    cout<<"i*halfnSites+Min2-nSites/2 ="<<i*halfnSites+Min2-nSites/2<<endl;
    if (i ==0 || i==3)
     eigenvectors[i*halfnSites +Min2-nSites/2]= 0.0;
    else
     eigenvectors[i*halfnSites +Min2-nSites/2]= beta;


   // eigenvectors[i*halfnSites +Min2-halfnSites]= beta*ratio;
  //  eigenvectors[i*halfnSites +Min3-halfnSites]= beta*ratio;
  }

//   halfnSites=halfnSites-6;  // this is for only butadiene

   for (int i =0 ; i < halfnSites; ++i){
       for (int j =0 ; j < halfnSites; ++j){
          cout<<eigenvectors[i*halfnSites+j];
         // cout<<eigenvectors[i*halfnSites+j];
       }
       cout<<endl;
   }

  dsyev(&jobz, &uplo, &halfnSites, eigenvectors, &halfnSites, eigenvalues, work, &lwork, &info);

  if (info != 0) {
    cout << "DMDCBase::CallDSYEV: DSYEV returned with info = " << info << endl;
    exit(1);
  }

//   for (int i =0 ; i < halfnSites; ++i)
//      cout<<"eigenvalues = "<<eigenvalues[i]<<endl;
    
/*
 
   for (int i =30 ; i < 60; ++i)
       for (int j =0 ; j < halfnSites; ++j)
          VirtOrbital[(i-30)*nSites+j+halfnSites]=eigenvectors[i*halfnSites+j];
   for (int i =0 ; i < 30; ++i)
       for (int j =0 ; j < halfnSites; ++j)
          OccOrbital[(i-0)*nSites+j+halfnSites]=eigenvectors[i*halfnSites+j];
*/


 
   for (int i =2 ; i < 4; ++i)
       for (int j =0 ; j < halfnSites; ++j)
          VirtOrbital[(i-2)*nSites+j+nSites/2]=eigenvectors[i*halfnSites+j];
   for (int i =0 ; i < 2; ++i)
       for (int j =0 ; j < halfnSites; ++j)
          OccOrbital[(i-0)*nSites+j+nSites/2]=eigenvectors[i*halfnSites+j];




//   for (int j =0 ; j < 2*nSites; ++j)
//      cout<<"VirtOrbital["<<j<<"]="<<VirtOrbital[j]<<endl;
   
//   for (int j =0 ; j < 2*nSites; ++j)
//      cout<<"OccOrbital["<<j<<"]="<<OccOrbital[j]<<endl;
  
   for (int i =25 ; i < 30; ++i)
     for (int j =30 ; j < 33; ++j){
       cout<<" i and j "<<i<<" "<<j<<endl;
       for (int ii =0 ; ii < halfnSites; ++ii){
          for (int jj =0 ; jj < halfnSites; ++jj){
           Overlapping[ii] += eigenvectors[i*halfnSites+ii]*eigenvectors[j*halfnSites+jj]; 
          }
       }
      }

//   for (int j =0 ; j < halfnSites; ++j)
//      cout<<"Overlapping["<<j<<"]="<<Overlapping[j]<<endl;





   //   double A = 0.5;
  //    double B = 0.05 ;
      double Be = 1 ;
      double Sijsum = 0.0 ;
      int icount=0;
      int Nocc = 0;
      int Nvir = 0;
      for (Nocc = 1; Nocc < 2 ; ++Nocc)
       for (Nvir=0; Nvir < 1 ; ++Nvir) 
        for (int i =0 ; i < nSites; ++i){
          for (int j =0 ; j < nSites; ++j){
             if ( (i < nAtomsArray[0] && j < nAtomsArray[0]) || (i >= nAtomsArray[0] && j >= nAtomsArray[0]) )  {
            //    Efield_comb[i] += Efield_pol3[i]- Efield_pol3[j];
            }
            else {
                  
           //  cout<<"x AngleArray: "<<i<<" and "<<j <<"= "<<AngleArray[i*3+0]<<" "<<AngleArray[j*3+0]<<endl;
           //  cout<<"y AngleArray: "<<i<<" and "<<j <<"= "<<AngleArray[i*3+1]<<" "<<AngleArray[j*3+1]<<endl;
           //  cout<<"z AngleArray: "<<i<<" and "<<j <<"= "<<AngleArray[i*3+2]<<" "<<AngleArray[j*3+2]<<endl;


             Sij_p =c1*AngleArray[i*3+0]*AngleArray[j*3+0]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 1, 0, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 1, 0, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+1]*AngleArray[j*3+1]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 1, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 1, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+2]*AngleArray[j*3+2]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 0, 1, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 0, 1, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
/*
             Sij_p =c1*AngleArray[i*3+0]*AngleArray[j*3+0]*
                    CalculateOverlapIntegral_S(c1, a1, 1, 0, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 1, 0, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+1]*AngleArray[j*3+1]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 1, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 1, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+2]*AngleArray[j*3+2]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 0, 1, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 0, 1, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
*/
             SijEach[i] += Sij_p;
           //  Efield_comb[i] += Sij_p*(Efield_pol3[i]- Efield_pol3[j]);
          //   cout<<"Sij_p="<<Sij_p<<endl;
            }
          }
          if (SijEach[i] < 0.0 ) Sign[icount]=-1;
          else Sign[icount]=1;

      //     cout<<" 1 SijEach["<<i<<"] = "<<SijEach[i]<<" "<<Sign[icount]<<" "<<icount<<endl;
           icount++;

           SijTotal[i] += abs(SijEach[i]);
        //   SijTotal[i] += abs(SijEach[i])/(eigenvalues[Nvir+30]-eigenvalues[Nocc])*0.780456;
         // SijTotal[i] += abs(SijEach[i])*exp(-Be*(eigenvalues[Nvir+30]-eigenvalues[Nocc]));
           SijEach[i] = 0.0;
        //   cout<<" 1 SijTotal["<<i<<"]= "<<SijTotal[i]<<endl; 
        }
      
      for (int i =0 ; i < nSites; ++i) {SijTotal[i]=0.0; SijEach[i]=0.0;}

      icount=0;
      double Const=1;
      Nocc = 0;
      Nvir = 0;

      double A = 0.5 ;
      double B = 0.04 ;

      for (Nocc = 1; Nocc < 2 ; ++Nocc)
       for (Nvir=0; Nvir < 1 ; ++Nvir)
        for (int i =0 ; i < nSites; ++i){
          for (int j =0 ; j < nSites; ++j){
             if ( (i < nAtomsArray[0] && j < nAtomsArray[0]) || (i >= nAtomsArray[0] && j >= nAtomsArray[0]) )  {
         //      if (Nocc == 1 && Nvir == 1) {
                Efield_comb[i] += Efield_pol3[i]- Efield_pol3[j];
                SijTotal[i] += 1;
          //     }
            }
            else {

             double Rij=sqrt((Rx[i]-Rx[j])*(Rx[i]-Rx[j])+(Ry[i]-Ry[j])*(Ry[i]-Ry[j])+(Rz[i]-Rz[j])*(Rz[i]-Rz[j]));
             double Sij = A*(1.0+ B*Rij +B*B*Rij*Rij/3.0)*exp(-B*Rij);


             Sij_p =c1*AngleArray[i*3+0]*AngleArray[j*3+0]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 1, 0, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 1, 0, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+1]*AngleArray[j*3+1]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 1, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 1, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+2]*AngleArray[j*3+2]*OccOrbital[nSites*Nocc+i]*VirtOrbital[nSites*Nvir+j]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 0, 1, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 0, 1, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
/*
             Sij_p =c1*AngleArray[i*3+0]*AngleArray[j*3+0]*
                    CalculateOverlapIntegral_S(c1, a1, 1, 0, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 1, 0, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+1]*AngleArray[j*3+1]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 1, 0, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 1, 0, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
             Sij_p +=c1*AngleArray[i*3+2]*AngleArray[j*3+2]*
                    CalculateOverlapIntegral_S(c1, a1, 0, 0, 1, Site[i*3+0], Site[i*3+1], Site[i*3+2], c2, a2, 0, 0, 1, Site[j*3+0], Site[j*3+1], Site[j*3+2]);
*/

           //  Sij_p=1;
           // SijEach[i] += Sij_p;
           // Efield_comb[i] += Sij_p*(Efield_pol3[i]- Efield_pol3[j]);
           SijEach[i] += Const*Sign[icount]*Sij_p;
             Efield_comb[i] += Const*Sign[icount]*Sij_p*(Efield_pol3[i]- Efield_pol3[j]);

         //  cout<<"Const*Sign[icount]*Sij_p = " <<Const*Sign[icount]*Sij_p<<endl;
         //    SijEach[i] += Const*Sij_p;
          //   Efield_comb[i] += Const*Sij_p*(Efield_pol3[i]- Efield_pol3[j]);

          //   SijEach[i] += Sij;
          //   Efield_comb[i] += Sij*(Efield_pol3[i]- Efield_pol3[j]);
            }
          }
          //  cout<<" 2 SijEach[i] = "<<SijEach[i]<<" "<<Sign[icount]<<endl;
           icount++;
           SijTotal[i] += SijEach[i];
        //  SijTotal[i] += SijEach[i]/(eigenvalues[Nvir+30]-eigenvalues[Nocc])*0.780456;
           SijEach[i] = 0.0;
         //  cout<<" 2 SijTotal["<<i<<"]= "<<SijTotal[i]<<endl; 
        }




        for (int i =0 ; i < nSites ; ++i) {
               Efield_comb[i] = Efield_comb[i]/SijTotal[i];
             //  Efield_comb[i] = Efield_comb[i]/60;
        }



       //   cout<<"Rx["<<i<<"]= "<<Rx[i]<<endl;
       //   cout<<"Ry["<<i<<"]= "<<Ry[i]<<endl;
      //    cout<<"Rz["<<i<<"]= "<<Rz[i]<<endl;
      //    cout<<"Site["<<i*3<<"]= "<<Site[i*3]<<endl;
      //    cout<<"Site["<<i*3+1<<"]= "<<Site[i*3+1]<<endl;
      //    cout<<"Site["<<i*3+2<<"]= "<<Site[i*3+2]<<endl;

     //  Efield[nSites] = 0.0000;
     //  n_pol3=nSites+1;


        for (int i =0 ; i < n_pol3 ; ++i){
          EfieldNet_pol3[i] =  0.5*Efield_pol3[i] ;
         // EfieldNet_pol3[i] =  - EfieldOnAtoms[i] ;
        //  EfieldNet_pol3[i] =  0.5*Efield_pol3[i] ;
        //  cout<<"EfieldOnAtoms[i]: "<< EfieldOnAtoms[i]<<endl; 
        // cout<<"Efield_pol3["<<i<<"]="<<Efield_pol3[i]<<endl; 
       //  cout<<"Efield_comb["<<i<<"]="<<Efield_comb[i]<<endl; 
       //  EfieldNet_pol3[i] =  0.5*Efield_comb[i] - EfieldOnAtoms[i] ;
       //   Efield_comb[i] =  0.5*Efield_comb[i] ;
        } 

      
//        for (int i =0 ; i < n_pol3 ; ++i)
//          cout<<"Efield_comb[i]: "<< Efield_comb[i]<<endl; 

 //       for (int i =0 ; i < n_pol3 ; ++i)
 //         cout<<"MolPol[0].InvA[0]: "<< MolPol[0].InvA[i]<<endl;


        dgemv("N", &n_pol3, &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &Efield_comb[0], &one, &dzero, &mu_pol3[0], &one);
     //   dgemv("N", &n_pol3, &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &EfieldNet_pol3[0], &one, &dzero, &mu_pol3[0], &one);

        //dsymv("L", &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &EfieldNet_pol3[0], &one, &dzero, &mu_pol3[0], &one);
        
         double summu=0.0;


        for (int i =0 ; i < nSites; ++i){
           cout<<"mu[ "<<i<<" ] = "<<mu_pol3[i]<<endl;
           summu += mu_pol3[i];
        }



/*
        for (int i =0 ; i < ncomb; ++i){
           cout<<"mu[ "<<i<<" ] = "<<mu_comb[i]<<endl;
           summu += mu_comb[i];
        }
*/

    //    cout<<"summu="<<summu<<endl;


/* 
// I'm not sure about that
double Q_xx = 0.0;
double Q_yy = 0.0;
double Q_zz = 0.0;
double Q_xy = 0.0;
double Q_xz = 0.0;
double Q_yz = 0.0;

  for (int i = 0; i < nSites; ++i){
    double Ri2 = Site[3*i+0]*Site[3*i+0] + Site[3*i+1]*Site[3*i+1] + Site[3*i+2]*Site[3*i+2];
    Q_xx += 3.0*2.0*mu_pol3[i]*Site[3*i+0]*Site[3*i+0] - Ri2;
    Q_yy += 3.0*2.0*mu_pol3[i]*Site[3*i+1]*Site[3*i+1] - Ri2;
    Q_zz += 3.0*2.0*mu_pol3[i]*Site[3*i+2]*Site[3*i+2] - Ri2;
    Q_xy += 3.0*2.0*mu_pol3[i]*Site[3*i+0]*Site[3*i+1];
    Q_xz += 3.0*2.0*mu_pol3[i]*Site[3*i+0]*Site[3*i+2];
    Q_yz += 3.0*2.0*mu_pol3[i]*Site[3*i+1]*Site[3*i+2];
    //cout<<"mu_pol3["<<i<<"] = "<<mu_pol3[i]<<endl;
    //cout<<"Site["<<3*i<<"] = "<<Site[3*i]<<endl;
  }

 cout<<"---Quadrupole moment  -------"<<endl;
 cout<<"Q_xx = "<<Q_xx<<endl;
 cout<<"Q_yy = "<<Q_yy<<endl;
 cout<<"Q_zz = "<<Q_zz<<endl;
 cout<<"Q_xy = "<<Q_xy<<endl;
 cout<<"Q_xz = "<<Q_xz<<endl;
 cout<<"Q_yz = "<<Q_yz<<endl;
*/






      /*
        for (int i =0 ; i < n_pol3 ; ++i){
           if (i%3 == 0) cout<<"mu_x["<<i<<"] = "<<mu_pol3[i]<<endl;
           if (i%3 == 1) cout<<"mu_y["<<i<<"] = "<<mu_pol3[i]<<endl;
           if (i%3 == 2) cout<<"mu_z["<<i<<"] = "<<mu_pol3[i]<<endl;
        }
        for (int i =0 ; i < n_pol3 ; ++i){
            cout<<"Efield_pol3["<<i<<"] = "<<Efield_pol3[i]<<endl;
        }
       */
       // Vpol = -ddot(&n_pol3, &Efield_pol3[0], &one, &mu_pol3[0], &one) ; // not using 0.5 as nothing is double counted
        Vpol = -ddot(&n_pol3, &Efield_comb[0], &one, &mu_pol3[0], &one) ; // not using 0.5 as nothing is double counted
    //    cout<<"Vpol = "<<Vpol<<endl;
        exit(0);

      break;
   }
      case 8:{
        // computing polarization and induction
        //

        for (int i =0 ; i < n ; ++i) {
          Efield_comb[i] =  0.00;
          SijTotal[i] =  0.00;
        }

       
      double Ex=0.00;
      double Ey=0.00;
      double Ez=0.001;

  //      cout<<"CarbonType = "<<CarbonType<<endl;
        for (int i =0 ; i < nSites; ++i){
          Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez; 
          Efield_comb[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez; 
          cout<<"Rz["<<i<<"]="<<Rz[i]<<endl;
        //cout<<"AtomType["<<i<<"] = "<<AtomType[i]<<endl;
/*
         if (AtomType[i] == 3 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 0.2451537;
        else if(AtomType[i] == 2 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 0.21150515; 
         else if(AtomType[i] == 1 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez ; 
*/
          
      //   if (AtomType[i] == 3 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 24.51537 ;
      //   else if(AtomType[i] == 2 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez - 0.150515; 
      //   else if(AtomType[i] == 1 ) Efield[i] = Rx[i]*Ex+Ry[i]*Ey+ Rz[i]*Ez  ;
              
                   
          Efield[nSites+3*i+0] = Ex;
          Efield[nSites+3*i+1] = Ey;
          Efield[nSites+3*i+2] = Ez;
         // cout<<"Ry="<<Ry[i]<<endl;

          Efield_comb[nSites+3*i+0] = Ex;
          Efield_comb[nSites+3*i+1] = Ey;
          Efield_comb[nSites+3*i+2] = Ez;

         }
 

        for (int i =0 ; i < n; ++i){
          // Efield[i]  += EfieldOnAtoms[i]; 
         // EfieldNet_new[i] = 0.5*Efield[i] - EfieldOnAtoms[i]; 
           cout<<"Efield["<<i<<"]= "<<Efield[i]<<endl; 
        }

   //    Efield[nSites] = 0.0000;
   //   n=nSites+nFullerenes;
    cout<<" n = "<<n<<endl;
 

double MSij[20][20];


        for (int i =0 ; i < 20; ++i)
         for (int j =0 ; j < 20; ++j)
             MSij[i][j]=0.0;
           

MSij[0][10]=0.0562701;
MSij[0][11]=-0.0291688;
MSij[0][12]=-0.0185307;
MSij[0][13]=0.0142048;
MSij[1][10]=0.0291688;
MSij[1][11]=-0.0214933;
MSij[1][12]=-0.0171377;
MSij[1][13]=0.0185307;
MSij[2][10]=0.0185307;
MSij[2][11]=-0.0171377;
MSij[2][12]=-0.0214933;
MSij[2][13]=0.0291688;
MSij[3][10]=0.0142048;
MSij[3][11]=-0.0185307;
MSij[3][12]=-0.0291688;
MSij[3][13]=0.0562701;


MSij[10][0]=0.0562701;
MSij[10][1]=-0.0291688;
MSij[10][2]=-0.0185307;
MSij[10][3]=0.0142048;
MSij[11][0]=0.0291688;
MSij[11][1]=-0.0214933;
MSij[11][2]=-0.0171377;
MSij[11][3]=0.0185307;
MSij[12][0]=0.0185307;
MSij[12][1]=-0.0171377;
MSij[12][2]=-0.0214933;
MSij[12][3]=0.0291688;
MSij[13][0]=0.0142048;
MSij[13][1]=-0.0185307;
MSij[13][2]=-0.0291688;
MSij[13][3]=0.0562701;


/*
       A = 1.0;
       B = 0.25 ;
        for (int i =0 ; i < nSites; ++i){
        double temp=0.0;
          for (int j =0 ; j < nSites; ++j){
            if ( (i < nAtomsArray[0] && j < nAtomsArray[0]) || (i >= nAtomsArray[0] && j >= nAtomsArray[0]) )  {
            // cout<<" i and j "<<i<<" "<<j<<endl;
             Efield_comb[i] += Efield[i]- Efield[j];
              SijTotal[i] += 1.0;
           //  cout<<" Efield[i] - Efield[j]"<<Efield[i]<<" - "<<Efield[j]<<endl;
          //   cout<<" Efield_comb[i] = "<<Efield_comb[i]<<endl;
            }
            else {
           //  cout<<" S i and j "<<i<<" "<<j<<endl;
             double Rij=sqrt((Rx[i]-Rx[j])*(Rx[i]-Rx[j])+(Ry[i]-Ry[j])*(Ry[i]-Ry[j])+(Rz[i]-Rz[j])*(Rz[i]-Rz[j]));
           //  double Sij = A*(1.0+ B*Rij +B*B*Rij*Rij/3.0)*exp(-B*Rij);
             double Sij = 0.0;
           //   SijTotal[i] += MSij[i][j];
               cout<<"MSij "<<MSij[i][j]<<endl;
              SijTotal[i] += Sij;
           //    cout<<"Rij and Sij "<<Rij<<" "<<Sij<<endl;
               Efield_comb[i] += Sij*(Efield[i]- Efield[j]);
            //   Efield_comb[i] += MSij[i][j]*(Efield[i]- Efield[j]);
            }
          }
        //  cout<<"SijTotal["<<i<<"]="<<SijTotal[i]<<endl;
         // cout<<"Efield_comb["<<i<<"]="<<Efield_comb[i]<<endl;
         }

*/


/*

        for (int i =0 ; i < nSites ; ++i) {
               Efield_comb[i] = Efield_comb[i]/SijTotal[i];
        }
*/

        for (int i =0 ; i < n ; ++i) {
              cout<<"Efield["<<i<<"]="<<Efield[i]<<endl;
        }
        for (int i =0 ; i < n ; ++i) {
              cout<<"Efield_comb["<<i<<"]="<<Efield_comb[i]<<endl;
        }

/*
        for (int i =0 ; i < n*n ; ++i)
         cout<<"MolPol[0].InvA["<<i<<"]="<<MolPol[0].InvA[i]<<endl;  
*/


        dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield_comb[0], &one, &dzero, &mu[0], &one);
     //   dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
       // dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldNet_new[0], &one, &dzero, &mu[0], &one);

         double summu=0.0;
        for (int i =0 ; i < 60; ++i){
           cout<<"mu[ "<<i<<" ] = "<<mu[i]<<endl;
           summu += mu[i];        
        }
        cout<<"summu="<<summu<<endl;


        Vpol = -0.5*ddot(&n, &Efield[0], &one, &mu[0], &one) ;
        cout<<"1 Vpol = "<<Vpol<<endl;
        Vpol = Vpol  - InterFullerPol;
        cout<<"2 Vpol = "<<Vpol<<endl;

       // Vpol = -ddot(&n, &Efield[0], &one, &mu[0], &one) ;
      //  cout<<"3 Vpol = "<<Vpol<<endl;
        exit(0);

      break;
      printf("EvaluateCsixtyPol: This should never happen\n");
      }
      default:
      exit(42);
    }
  
  //cout << "electron position =  " << x[0] << x[1]  << x[2] <<endl ; 
  //for (int i =0 ; i < 4*nSites+1; ++i){
  //cout << "Efield" << i << " = " << Efield[i] << endl;  } 
  //cout << "mu" << i << " = " << mu[i] << endl;  } 
  //cout << "Vpol= " << Vpol << endl;   
  return Vpol  ;
}

double Potential::InterFullerPolPot(const double *EfieldVec )
{
  dVec Efield; Efield.resize(nSites*4+nFullerenes);
  dVec mu;     mu.resize(nSites*4+nFullerenes);

  double Vpol = 0.0 ;
  double dzero = 0.0;
  double done = 1.0;
  int one = 1;
  int n = nSites*4+nFullerenes ;
//  for (int ii=0; ii<n; ii++) 
//      cout<<"EfieldVec ["<<ii<<"]="<<EfieldVec[ii]<<endl;


  dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldVec[0], &one, &dzero, &mu[0], &one);
  //dsymv("L", &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldVec[0], &one, &dzero, &mu[0], &one);
  if (verbose > 5){
    for (int i = 0 ; i < nSites ; ++i)
    {
      cout << "this is the intermolecular induced dipoles on atom " << i << " = " << mu[i] << "\n"; 
      cout << "this is the intermolecular induced dipoles on atom " << i << " = " << mu[nSites+3*1+0] << "\n"; 
      cout << "this is the intermolecular induced dipoles on atom " << i << " = " << mu[nSites+3*1+1] << "\n"; 
      cout << "this is the intermolecular induced dipoles on atom " << i << " = " << mu[nSites+3*1+2] << "\n"; 
    }
  }
  Vpol = -0.5*ddot(&n, &EfieldVec[0], &one, &mu[0], &one);
  return Vpol  ;
}


/*

double Potential::CalculateOverlapIntegral_S(double c1, double a1, int l1, int m1, int n1,
                                            double x1, double y1, double z1,
                                            double c2, double a2, int l2, int m2, int n2,
                                            double x2, double y2, double z2)
{
        double AB2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        double K = exp(-a1*a2/(a1+a2)*AB2) * pow(PI/(a1+a2),1.5);
        double Ix = 0.0;        double Iy = 0.0;        double Iz = 0.0;

        double PAx = a2/(a1+a2) * (x2-x1), PBx = a1/(a1+a2) * (x1-x2);
        double PAy = a2/(a1+a2) * (y2-y1), PBy = a1/(a1+a2) * (y1-y2);
        double PAz = a2/(a1+a2) * (z2-z1), PBz = a1/(a1+a2) * (z1-z2);


        int icount, jcount ;
        double tempvariable;
        for ( icount = 0; icount <= (l1 + l2)/2 ; icount++ )
        {
                tempvariable = 1;
                for ( jcount = 3; jcount <= (2*icount - 1) ; jcount += 2)
                        tempvariable *= jcount;
                Ix +=
                tempvariable / pow((2*(a1+a2)),icount) * fil1l2PAPB(2*icount,l1,l2,PAx,PBx);
        }
        for ( icount = 0; icount <= (m1 + m2)/2 ; icount++)
        {
                tempvariable = 1;
                for ( jcount = 3; jcount <= (2*icount - 1) ; jcount += 2)
                        tempvariable *= jcount;
                Iy +=
                tempvariable / pow((2*(a1+a2)),icount) * fil1l2PAPB(2*icount,m1,m2,PAy,PBy);
        }
        for ( icount = 0; icount <= (n1 + n2)/2 ; icount++)
        {
                tempvariable = 1;
                for ( jcount = 3; jcount <= (2*icount - 1) ; jcount += 2)
                        tempvariable *= jcount;
                Iz +=
                tempvariable / pow((2*(a1+a2)),icount) * fil1l2PAPB(2*icount,n1,n2,PAz,PBz);
        }
        return c1*c2*Ix*Iy*Iz*K;
}
*/


double Potential::CalculateOverlapIntegral_S(double c1, double a1, int l1, int m1, int n1,
                                            double x1, double y1, double z1,
                                            double c2, double a2, int l2, int m2, int n2,
                                            double x2, double y2, double z2)
{
        double R = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
        double Overlap;
        double omega=R*a1;
        double omega2=omega*omega;
        double omega3=omega2*omega;
        double omega4=omega3*omega;

        if (n1 == 1 && n2 == 1 ) {
           Overlap=exp(-omega)*(1+omega+omega2/5.0-2*omega3/15.0-omega4/15.0);
        }
        else {
           Overlap=exp(-omega)*(1+omega+2*omega2/15.0+omega3/15.0);
        }

        return c1*c2*Overlap;
}




double Potential::fil1l2PAPB(int i, int l1, int l2, double PA, double PB)
{
        double result  = 0.0;
        if (l2 < l1)
        {
                int temp = l1;
                l1 =l2;
                l2 = temp;
                double tempdouble = PA;
                PA = PB;
                PB = tempdouble;
        }

        int jcount;
        double factoriall1=1, factoriall2=1;
        for ( jcount = 1; jcount < l1;factoriall1 *= ++jcount);
        for ( jcount = 1; jcount < l2;factoriall2 *= ++jcount);

        for (int lamda = 0 ; lamda <= l1 ; lamda++)
        {
                int u = i - lamda;
                if ((u >= 0) && (u <= l2))
                {
                        double factoriallamda=1, factorialu=1;
                        for ( jcount = 1; jcount < lamda;factoriallamda *= ++jcount);
                        for ( jcount = 1; jcount < u;factorialu *= ++jcount);
                        int lamda2 = l1 - lamda;
                        int u2 = l2 - u;
                        double factoriallamda2=1, factorialu2=1;
                        for ( jcount = 1; jcount < lamda2;factoriallamda2 *= ++jcount);
                        for ( jcount = 1; jcount < u2;factorialu2 *= ++jcount);

                        result +=  factoriall1 * factoriall2 / (factoriallamda * factoriallamda2 * factorialu * factorialu2)
                                * pow(PA,lamda2) * pow(PB,u2);
                }
        }
        return result;
}




void Potential::find_closed_aoms(int Index, int &Min1, int &Min2, int &Min3)
{
       
  dVec distcheck; distcheck.resize(nSites);
  iVec distIndex; distIndex.resize(nSites);
  int current_atom = 0;
  int old_atom = 0;
  double swap;
  int iswap;
  for (int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom = current_atom ;
    current_atom += nAtomsArray[iFuller] ;


    if ( (old_atom <= Index) && (Index  < current_atom) ) {
       current_atom=current_atom-6; //this is for only butadiene
      // cout<<" Index ="<<Index<<" old_atom ="<<old_atom<<"current_atom="<<current_atom<<endl;
       for (int j =old_atom ; j < current_atom; ++j){
           distcheck[j]=sqrt( (Site[Index*3+0]-Site[j*3+0])* (Site[Index*3+0]-Site[j*3+0])+
                          (Site[Index*3+1]-Site[j*3+1])* (Site[Index*3+1]-Site[j*3+1])+
                          (Site[Index*3+2]-Site[j*3+2])* (Site[Index*3+2]-Site[j*3+2]));
           distIndex[j]=j;
       }
      // cout<<" old_atom ="<<old_atom<<"current_atom="<<current_atom<<endl;
       for (int c = old_atom ; c < ( current_atom - 1 ); c++) {
          // cout<< " c  "<<c<<endl;
         for (int d = old_atom ; d < current_atom - (c - old_atom) - 1; d++) {
           if (distcheck[d] > distcheck[d+1]) /* For decreasing order use < */
           {
            swap       = distcheck[d];
            distcheck[d]   = distcheck[d+1];
            distcheck[d+1] = swap;
         //  cout<<"distcheck["<<d<<"]="<<distcheck[d]<<" distIndex ="<<distIndex[d]<<endl;
            iswap       = distIndex[d];
            distIndex[d]   = distIndex[d+1];
            distIndex[d+1] = iswap;
           }
         }
       }
       Min1=distIndex[old_atom+1];
       Min2=distIndex[old_atom+2];
       Min3=distIndex[old_atom+3];
       
    //   for (int j =old_atom ; j < current_atom; ++j){
    //       cout<<"distcheck["<<j<<"]="<<distcheck[j]<<" distIndex ="<<distIndex[j]<<endl;
    //   }



    }
   }

 //  cout<<"Min1= "<<Min1<<endl;
 }





