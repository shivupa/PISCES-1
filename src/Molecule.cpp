#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>

#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>

#include "constants.h"
#include "ttm_data.h"
#include "GeoAux.h"
#include "Molecule.h"

using namespace std;

/////////////////////////
//
// input:  WaterCoor: xO, yO, zO, xH1, yH1, zH1, xH2, yH2, zH2  in Angstrom
//
//
void Molecule::SetupDPPWater(const double *WaterCoor)
{

  Type = 1;
  nSites = 4;
  Tag = new char[nSites];
  Tag[0] = 'O';
  Tag[1] = 'H';
  Tag[2] = 'H';
  Tag[3] = 'M';

  Coor.resize(3*nSites); 

  nCharges = 3;
  ChargeList.resize(nCharges);
  ChargeSite.resize(nCharges);
  ChargeList[0] = TTM_QH;  ChargeSite[0] = 1;
  ChargeList[1] = TTM_QH;  ChargeSite[1] = 2;
  ChargeList[2] = TTM_QM;  ChargeSite[2] = 3;

  
  // these dipoles must be determined in a DPP calculation, or totally self-consistently with the electron
  // as yet, they are not needed
  //nDipoles = 3;
  //DipoleList.resize(3*nDipoles);
  //DipoleSite.resize(nDipoles);
  //DipoleSite[0] = 0; 
  //DipoleSite[1] = 1;
  //DipoleSite[2] = 2;
  //for (int i = 0; i < 3*nDipoles; ++i)
  //  DipoleList[i] = 0.0;

  nPntPols = 3;
  PntPolList.resize(nPntPols);
  PntPolSite.resize(nPntPols);
  PntPolSite[0] = 0; PntPolList[0] = 0.837 * Angs2Bohr*Angs2Bohr*Angs2Bohr;
  PntPolSite[1] = 1; PntPolList[1] = 0.496 * Angs2Bohr*Angs2Bohr*Angs2Bohr;
  PntPolSite[2] = 2; PntPolList[2] = 0.496 * Angs2Bohr*Angs2Bohr*Angs2Bohr;

  int debug = 1;

  if (debug) {
    printf("Input structure [Angs]:\n");
    for (int j = 0; j < 3; ++j)
      printf("  %c  %12.7f  %12.7f  %12.7f\n", Tag[j], WaterCoor[3*j],  WaterCoor[3*j+1], WaterCoor[3*j+2]);
  }
  cout.flush();

  // make rigid water: go to Euler angles, and rotate reference structure
  // if I ever needed the Euler angles, or the rotation matrix in Molecule
  // RigidWater returns them, so far they are discarded

  double Conf[6];
  double MSiteDisplacement = DPP_DeltaM;
  RigidWater(WaterCoor, Angs2Bohr, Conf, &Coor[0], MSiteDisplacement); 
  
  if (debug) {
    printf("Rigid structure [Angs]:\n");
    for (int j = 0; j < 4; ++j)
      printf("  %c  %12.7f  %12.7f  %12.7f\n", Tag[j], Coor[3*j]*Bohr2Angs,  Coor[3*j+1]*Bohr2Angs, Coor[3*j+2]*Bohr2Angs);
    printf("  center = %f, %f, %f\n", Conf[0]*Bohr2Angs, Conf[1]*Bohr2Angs, Conf[2]*Bohr2Angs);
    printf("  Euler = %f, %f, %f\n", Conf[3], Conf[4], Conf[5]);
    //exit(1);
  }
  cout.flush();

}
