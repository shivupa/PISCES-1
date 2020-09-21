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
#include "C60.h"
#include "ChargeDipPol.h"




using namespace std;

//////////////////////////////
//
//  IonTypes: 
//
const int C_atom = 1 ;







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
      default:
      printf("This should never happen\n");
      exit(42);
    }
  }
  class Potential Vel;
  double PotPara[100];
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

      Vel.SetupFullerElec(nIons, IonPos, PotPara, nAtomsArray, DipoleArray, CenterArray);
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
  Hel.SetupDVR(ngrid, DVRType, 1, gpara, gridverbose);
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

//  if (nStates > 0)
//    Hel.ExpectationValues(1);
  Hel.EnergyPartitioning(Vel) ;
  // plot the wavefunction
  int *NucCharge = new int[nIons];
  for (int k = 0; k < nIons; ++k) {
    switch (IonType[k])
      {
      case C_atom:
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

    for (int istate = 1; istate <= nStates ; ++istate) {
      char fname[15];
      sprintf(fname, "WaveFn%02i.cube", istate);
      Hel.WriteCubeFile(istate, fname, nIons , &NucCharge[0], &IonPos[0]);
    }
    Hel.WriteOneDCuts() ; 

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
void Potential::SetupFullerElec(int nr, const double *rSites, const double *PotPara, const int *NoAtomsArray, const double *DipoleArray, const double *CenterArray)
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
  SetMoleculeFields( NoAtomsArray) ;
  MolPol = new DistributedPolarizabilities[1];
  
  if (PolFlagCsixty == 1 || PolFlagCsixty == 2) {
    MolPol[0].InvA.resize((4*nr + nFullerenes)*(4*nr + nFullerenes));
    InvChargeDipolePol(nSites, &rSites[0], &MolPol[0].InvA[0], nFullerenes, &NoAtomsArray[0]);
    InterFullerPol = InterFullerPolPot(&EfieldOnAtoms[0]) ; 
    cout << "Intermolecular polarization is" << InterFullerPol << "au \n" ; 
    }
  if (PolFlagCsixty == 3) {
    cout << "Using new Polarization =3 that has intermolecular charge flow "<< PolFlagCsixty << endl ; 
    iVec NoAtomsArray_pol3 ; NoAtomsArray_pol3.resize(1) ; 
    NoAtomsArray_pol3[0] = nr ; 
    MolPol[0].InvA.resize((4*nr + 1)*(4*nr + 1)); 
    InvChargeDipolePol(nSites, &rSites[0], &MolPol[0].InvA[0], 1, &NoAtomsArray_pol3[0]); 
    }

  //  int n4plus1 = nr*4+1 ; 
  //  for (int i = 0; i < n4plus1*n4plus1; ++i){
  //    cout << "InvA" << i << "=" << MolPol[0].InvA[i] << endl;
  //  }
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
    for (int i =old_atom ; i < current_atom; ++i){
      double RSiteCenter = sqrt((Site[i*3+0]-CenterArray[iFuller*3+0])*(Site[i*3+0]-CenterArray[iFuller*3+0]) +
                                (Site[i*3+1]-CenterArray[iFuller*3+1])*(Site[i*3+1]-CenterArray[iFuller*3+1]) +
                                (Site[i*3+2]-CenterArray[iFuller*3+2])*(Site[i*3+2]-CenterArray[iFuller*3+2])) ;
      PtDip[3*i+0] = DipoleArray[iFuller]*(Site[i*3+0] - CenterArray[iFuller*3+0])/RSiteCenter ; 
      PtDip[3*i+1] = DipoleArray[iFuller]*(Site[i*3+1] - CenterArray[iFuller*3+1])/RSiteCenter ;
      PtDip[3*i+2] = DipoleArray[iFuller]*(Site[i*3+2] - CenterArray[iFuller*3+2])/RSiteCenter ;
    }
  }
}

void Potential::SetMoleculeFields(const int *NoAtomsArray)
{ 
  // calculating the fields on the atoms of a given molecule from all other molecules
  int current_atom = 0;
  int old_atom = 0;
  EfieldOnAtoms.resize(nSites*4+nFullerenes) ;
  for (int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom = current_atom ; 
    current_atom += NoAtomsArray[iFuller] ;
    for (int i =old_atom ; i < current_atom; ++i){
      EfieldOnAtoms[nSites+i*3+0] = 0.0 ;
      EfieldOnAtoms[nSites+i*3+1] = 0.0 ; 
      EfieldOnAtoms[nSites+i*3+2] = 0.0 ;
      EfieldOnAtoms[i] = 0.0 ; 
      for (int j = 0  ; j < nSites ; ++j){
        // avoids fields on a molecules from itself
        if (j < old_atom | j >= current_atom) {
          //{cout << "i atom is "<< i << " j atom is " << j << "\n";}
          double rSite_x = Site[i*3+0] - Site[j*3+0] ;
          double rSite_y = Site[i*3+1] - Site[j*3+1] ;
          double rSite_z = Site[i*3+2] - Site[j*3+2] ;
          double rSite   = sqrt(rSite_x*rSite_x + rSite_y*rSite_y + rSite_z*rSite_z) ;
          double rSite2  = rSite*rSite ;
          double rSite3  = rSite2*rSite ;
          double MuDotR  = PtDip[3*j+0]*rSite_x + PtDip[3*j+1]*rSite_y + PtDip[3*j+2]*rSite_z ;  
          // potential on atoms  due to dipoles
          EfieldOnAtoms[i]            += MuDotR/rSite3;
          //fields on atoms due to dipoles
          EfieldOnAtoms[nSites+i*3+0] += (3.0*MuDotR*rSite_x/rSite2 - PtDip[3*j+0])/rSite3 ; 
          EfieldOnAtoms[nSites+i*3+1] += (3.0*MuDotR*rSite_y/rSite2 - PtDip[3*j+1])/rSite3 ; 
          EfieldOnAtoms[nSites+i*3+2] += (3.0*MuDotR*rSite_z/rSite2 - PtDip[3*j+2])/rSite3 ; 
        }
      }
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
   
double Potential::EvaluateCsixty(const double *x)
{
  double Vtotal = 0.0 ; 
  double VRep = EvaluateCsixtyRepPot(x, GaussExpCsixty, RepScaleCsixty) ; 
  double Vel = EvaluateCsixtyElecPot(x);
  double VPol =  EvaluateCsixtyPolPot(x) ;
  ReturnEnergies[0] = Vel;
  ReturnEnergies[2] = VRep;
  ReturnEnergies[3] = VPol;
  if (verbose > 0)
    CheckMinMax(&Vel);
  if (verbose > 15)
    printf("%10.6f %10.6f %10.6f %14.6e\n", x[0], x[1], x[2], Vel);
  Vtotal = Vel + VPol + VRep ; 
  //Vtotal = Vel + VRep ; 
  if (verbose > 0)
    printf("%10.6f %10.6f %10.6f %14.6e %14.6e %14.6e %14.6e \n", x[0], x[1], x[2], Vel, VRep, VPol, Vtotal);
  return Vtotal ;
}


double Potential::EvaluateCsixtyElecPot(const double *x)
{ 
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
    Velec += -1.0*ExtraNucQ/Reff + (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*gij ;
    
  }
  return Velec ;
   
}

double Potential::EvaluateCsixtyRepPot(const double *x, double GaussExp, double RepScale)
{
  double Vrep = 0.0 ;
  switch (RepFlagCsixty)
    {
      case 0:
        Vrep = 0.0 ; //no repulsion
      break;
      case 1:
        for (int i = 0; i < nSites; ++i) {
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
  dVec Efield; Efield.resize(nSites*4+nFullerenes);
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
          Efield[i] = -1.0/Reff   ; 
          Efield[nSites+3*i+0] = Reffminus3 * Rx[i];
          Efield[nSites+3*i+1] = Reffminus3 * Ry[i];
          Efield[nSites+3*i+2] = Reffminus3 * Rz[i];
        }
        // adding fields from other fullerenes
        // make sure fields are ADDED
        for (int i =0 ; i < n; ++i){
          Efield[i] -= EfieldOnAtoms[i]; 
        }
        // compute dipoles induced by Efield
        dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);
        Vpol = -0.5*ddot(&n, &Efield[0], &one, &mu[0], &one) - InterFullerPol;
        //Vpol = -0.5*ddot(&n, &Efield[0], &one, &mu[0], &one) ;
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
          Efield[i] = -1.0/Reff   ;
          Efield[nSites+3*i+0] = Reffminus3 * Rx[i];
          Efield[nSites+3*i+1] = Reffminus3 * Ry[i];
          Efield[nSites+3*i+2] = Reffminus3 * Rz[i];
        }
        dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldOnAtoms[0], &one, &dzero, &mu[0], &one);
        Vpol = -ddot(&n, &Efield[0], &one, &mu[0], &one) ; // not using 0.5 as nothing is double counted
      break ;
      case 3:
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
          EfieldNet_pol3[i] =  0.5*Efield_pol3[i] - EfieldOnAtoms[i] ;
        }
        dgemv("N", &n_pol3, &n_pol3, &done, &(MolPol[0].InvA[0]), &n_pol3, &EfieldNet_pol3[0], &one, &dzero, &mu_pol3[0], &one);
        Vpol = -ddot(&n_pol3, &Efield_pol3[0], &one, &mu_pol3[0], &one) ; // not using 0.5 as nothing is double counted
      break;
      default:
      printf("EvaluateCsixtyPol: This should never happen\n");
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
  dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &EfieldVec[0], &one, &dzero, &mu[0], &one);
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

