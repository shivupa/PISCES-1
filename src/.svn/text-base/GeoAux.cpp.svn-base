#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 
#include <iostream>

#include "constants.h"
#include "ttm_data.h"

////////////////////////////////////////////////////////////////////////////////
//
// make water rigid: go to Euler angles, and rotate back from reference structure
//
// input:  WaterCoor: xO, yO, zO, xH1, yH1, zH1, xH2, yH2, zH2 
//         scale:  factor to convert input coordinates to Bohr
//         MdZ  displacement of the M site from oxygen
//
//
// output: Conf[6]: Center-of-mass coordiantes and Euler angles
//         NewCoor: Rigid coordinates O, H, H, and M
//
void RigidWater(const double *WaterCoor, const double scale, double *Conf, double *NewCoor, const double MdZ) 
{

  int debug = 0;

  if (debug) std::cout << "\n-----------------------------RigidWater-------------------------\n";


  // center of mass
  for (int k = 0; k <3; ++k)
    Conf[k] = (16.0 * WaterCoor[k]*scale + WaterCoor[3+k]*scale + WaterCoor[6+k]*scale) / 18.0;


  // to compute the Euler angles establish a molecule fixed corrdinate system

  double OH1[3], OH2[3];      // vectors: O->H
  double z[3];                // z-axis is along the bisector
  double nrm = 0;
  for (int k = 0; k <3; ++k) {
    OH1[k] = (WaterCoor[3+k] - WaterCoor[k]) * scale;
    OH2[k] = (WaterCoor[6+k] - WaterCoor[k]) * scale;
    z[k] = OH1[k] + OH2[k];
    nrm += z[k]*z[k];
  }
  nrm = 1.0 / sqrt(nrm);
  for (int k = 0; k <3; ++k)
    z[k] *= -nrm;

  if (debug) printf("molecule-fixed z = (%f  %f  %f)\n", z[0], z[1], z[2]); 

  double y[3]; // y-axis from y = OH1 - z*(OH1*z)
  double ovl = OH1[0]*z[0] + OH1[1]*z[1] + OH1[2]*z[2];
  nrm = 0;
  for (int k = 0; k <3; ++k) {
    y[k] = OH1[k] - ovl*z[k];
    nrm += y[k]*y[k];
  }
  nrm = 1.0 / sqrt(nrm);
  for (int k = 0; k <3; ++k)
    y[k] *= nrm;

  if (debug) printf("molecule-fixed y = (%f  %f  %f)\n", y[0], y[1], y[2]);


  double x[3]; // x = y cross-product z   (should come out normalized)
  x[0] = y[1]*z[2] - y[2]*z[1];
  x[1] = y[2]*z[0] - y[0]*z[2];
  x[2] = y[0]*z[1] - y[1]*z[0];

  if (debug) printf("molecule-fixed x = (%f  %f  %f)\n", x[0], x[1], x[2]); 

  double theta = acos(z[2]);
  double sintheta = sin(theta);
  double psi, phi, tempcos, aCosTemp;

  if (debug) {
    std::cout << "theta = " << theta << "\n";
    std::cout << "sin(theta) = " << sintheta << "\n";
  }

  if (fabs(sintheta) > 1e-5)  {
    // psi
    tempcos = -x[2]/sintheta; 
    if (tempcos >= 1.0)
      aCosTemp = 0.0;
    else if (tempcos <= -1.0) 
      aCosTemp = PI;
    else
      aCosTemp = acos(tempcos);
    if (y[2] >= 0) 
      psi = aCosTemp;
    else 
      psi = 2*PI - aCosTemp;
    // phi
    tempcos = z[0]/sintheta;
    if (tempcos > 1.0) tempcos = 1.0;
    if (tempcos < -1.0) tempcos = -1.0;
    if (z[1] >= 0) 
      phi = acos(tempcos);
    else 
      phi = 2*PI - acos(tempcos);
  }
  else { // that means we are close to north or south pole, i.e. theta = 0 or theta = PI
    phi = 0.0;
    if (x[1] >= 0) 
      psi = acos(y[1]);
    else 
      psi = 2*PI - acos(y[1]);
  }
  
  if (debug) {
    std::cout << "psi = " << psi << "\n";
    std::cout << "phi = " << phi << "\n";
  }
  
  Conf[3] = phi;
  Conf[4] = theta;
  Conf[5] = psi;


  // Now start with rigid reference position, rotate with Euler angles, and then shift
  double RefCoor[12];
  RefCoor[0]  = 0.0;
  RefCoor[1]  = 0.0;
  RefCoor[2]  = TTM_Oz * Angs2Bohr;
  RefCoor[3]  = 0.0;
  RefCoor[4]  = TTM_Hy * Angs2Bohr;
  RefCoor[5]  = TTM_Hz * Angs2Bohr;
  RefCoor[6]  = 0.0;
  RefCoor[7]  = -TTM_Hy * Angs2Bohr;
  RefCoor[8]  = TTM_Hz * Angs2Bohr;
  RefCoor[9]  = 0.0;
  RefCoor[10] = 0.0;
  RefCoor[11] = (TTM_Oz - MdZ) * Angs2Bohr;

  // this is the rotation matrix in terms of sines and cosines of the three Euler angles
  double c1 = cos(phi), c2 = cos(theta), c3 = cos(psi);
  double s1 = sin(phi), s2 = sin(theta), s3 = sin(psi);
  double NewtoOld[9];
  NewtoOld[0] = c1*c2*c3-s1*s3, NewtoOld[1] = -c1*c2*s3-s1*c3, NewtoOld[2] = c1*s2;
  NewtoOld[3] = s1*c2*c3+c1*s3, NewtoOld[4] = -s1*c2*s3+c1*c3, NewtoOld[5] = s1*s2;
  NewtoOld[6] = -s2*c3,         NewtoOld[7] = s2*s3,           NewtoOld[8] = c2;

  for(int i = 0; i < 4; ++i) {
    NewCoor[3*i+0] = NewtoOld[0] * RefCoor[3*i+0] + NewtoOld[1] * RefCoor[3*i+1] + NewtoOld[2] * RefCoor[3*i+2];
    NewCoor[3*i+1] = NewtoOld[3] * RefCoor[3*i+0] + NewtoOld[4] * RefCoor[3*i+1] + NewtoOld[5] * RefCoor[3*i+2];
    NewCoor[3*i+2] = NewtoOld[6] * RefCoor[3*i+0] + NewtoOld[7] * RefCoor[3*i+1] + NewtoOld[8] * RefCoor[3*i+2];  
  }

  // shift to center of mass
  NewCoor[0]  += Conf[0];
  NewCoor[1]  += Conf[1];
  NewCoor[2]  += Conf[2];
  NewCoor[3]  += Conf[0];
  NewCoor[4]  += Conf[1];
  NewCoor[5]  += Conf[2];
  NewCoor[6]  += Conf[0];
  NewCoor[7]  += Conf[1];
  NewCoor[8]  += Conf[2];
  NewCoor[9]  += Conf[0];
  NewCoor[10] += Conf[1];
  NewCoor[11] += Conf[2];


}
