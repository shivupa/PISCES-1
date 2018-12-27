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
      cout << "Error in EvaluateVel: so far no potential has been set up.\n";
      exit(1);
      break;
   case 1: // this should work for 2 as well
   case 2:
   case 3:
      EvaluateDPPGTOPGradient(relectron,Grad, tmu, mu_cross_mu, wavefn,  WaterN);
      break;
   default:
      cout << "Error in EvaluateGradient: unknown PotFlags[0]; this should not happen\n";
      exit(1);
   }

   // never ever get here
}



void Potential::EvaluateDPPGTOPGradient(
                                       const double *x,
                                       dVec& Grad, dVec& tmu,dVec& mu_cross_mu, double wavefn,
                                       class WaterCluster &WaterN) 
{

// gradients  due to point charge electrostatic potential



//  double DampParameter = 1.7724;

        for (int i = 0; i < nCharges; ++i) {
            double VpcGrad = 0.0 ; 
            int iqs = ChargeSite[i];
            double Reff = R[iqs];
            if (Reff < ChargeDamping) {
              double ror0 = Reff / ChargeDamping;
              Reff = ChargeDamping * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
              VpcGrad = Charge[i] * (3.0*ror0*ror0 - 2.0*ror0*ror0*ror0) / (Reff*Reff);
            }
            else{
              VpcGrad = Charge[i] / (R[iqs]*R[iqs]);
            }
    
            Grad[iqs*3]   += -VpcGrad*Rx[iqs] / R[iqs] ;
            Grad[iqs*3+1] += -VpcGrad*Ry[iqs] / R[iqs] ;
            Grad[iqs*3+2] += -VpcGrad*Rz[iqs] / R[iqs];
      }

// gradients  for exchange repulsion potential
  for (int i = 0; i < nGauss; ++i) {
     int igs = GaussSite[i];
     double de1 = Gauss[2*i+1]*(-Gauss[2*i]) * exp(-Gauss[2*i] * R[igs]);
     Grad[igs*3] += -de1*Rx[igs]/R[igs];
     Grad[igs*3+1] += -de1*Ry[igs]/R[igs];
     Grad[igs*3+2] += -de1*Rz[igs]/R[igs];
  }

 
 


// we need this part again due to dipole calculation
// compute electron's Efield at atomic sites 

//  DampParameter = 4.895;
//  double DampParameter = 5.005;

  int nAtoms = MolPol[0].nAtoms;
  int n = 3*nAtoms;  // dimension of InvA and Efield, and mu

  dVec Efield; Efield.resize(n);
  dVec mu;     mu.resize(n);

  double gij;
  double de1;
  double de2;
  for (int i = 0; i < nAtoms; ++i) {
    int SiteIndex = MolPol[0].SiteList[i];
    gij = Rminus3[SiteIndex];
    double Reff = R[SiteIndex];
    if (Reff < PolDamping ) {
       double ror0 = Reff / PolDamping;
       Reff = PolDamping * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
       gij = 1.0 / (Reff * Reff * Reff);
    }

    Efield[3*i+0] = gij * Rx[SiteIndex];
    Efield[3*i+1] = gij * Ry[SiteIndex];
    Efield[3*i+2] = gij * Rz[SiteIndex];

  }

  for (int k = 0; k < n; ++k)
    Efield[k] += MolPol[0].Epc[k];
  
  // compute dipoles induced by Efield
  double dzero = 0.0;
  double done = 1.0;
  int one = 1;
  int dim_x = 3 ;
  dgemv("N", &n, &n, &done, &(MolPol[0].InvA[0]), &n, &Efield[0], &one, &dzero, &mu[0], &one);


//////////////////////////////////////////////////////////////////////
//
//start gradient calculations  written by Tae Hoon Choi 2011, 12,14 
//
//
//modified by vkv;
//the gradient formula is
//dVpol/dR = mu * (2*dEfield/dR + dT/dR *mu )
//////////////////////////////////////////////////////////////////////

   nSites = nAtoms/3*4;
   // Calculate derivative of Ee which is from interaction between  electron and atoms
   // dE(elec)/dR
   DerivElecField ( nSites, PolDamping, &Rminus3[0] , &R[0] , &Grad[0] , &mu[0]) ;

   // this just calculates a kind of cross product of mu^T x mu and accumulates density wighted mu
   // but is the most time consuming

   MuCrossMu( nAtoms, &mu[0] , tmu , mu_cross_mu, wavefn) ;  


}

void Potential::DerivElecField ( int nSites, double DampParameter, double *Rminus3 , double *R , double *Grad , double *mu)
{
   double gij;
   int n = nSites*9/4 ;
   int one = 1 ;
   dVec dxEfield; dxEfield.resize(n);
   dVec dyEfield; dyEfield.resize(n);
   dVec dzEfield; dzEfield.resize(n);

   int index =0;
   for (int j = 0; j < nSites ; ++j) {
  
      int SiteIndex = MolPol[0].SiteList[index];
 
      if (SiteIndex == j ) {
      gij = Rminus3[SiteIndex];
      double Reff = R[SiteIndex];
      double de_gij = -3.0/(Reff*Reff*Reff*Reff);
      if (Reff < DampParameter) {
         double ror0 = Reff / DampParameter;
         Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
         de_gij = -3.0/(Reff*Reff*Reff*Reff);
         double de_Reff = 3.0*ror0 * ror0  - 0.5 * 4.0*ror0*ror0*ror0;
         gij = 1.0 / (Reff * Reff * Reff);
         de_gij = de_Reff*de_gij;
      }


      dxEfield[3*index+0] = -de_gij*Rx[SiteIndex]/R[SiteIndex]*Rx[SiteIndex] - gij ;
      dxEfield[3*index+1] = -de_gij*Rx[SiteIndex]/R[SiteIndex]*Ry[SiteIndex]       ;
      dxEfield[3*index+2] = -de_gij*Rx[SiteIndex]/R[SiteIndex]*Rz[SiteIndex]       ;

      dyEfield[3*index+0] = -de_gij*Ry[SiteIndex]/R[SiteIndex]*Rx[SiteIndex]       ;
      dyEfield[3*index+1] = -de_gij*Ry[SiteIndex]/R[SiteIndex]*Ry[SiteIndex] - gij ;
      dyEfield[3*index+2] = -de_gij*Ry[SiteIndex]/R[SiteIndex]*Rz[SiteIndex]       ;

      dzEfield[3*index+0] = -de_gij*Rz[SiteIndex]/R[SiteIndex]*Rx[SiteIndex]       ;
      dzEfield[3*index+1] = -de_gij*Rz[SiteIndex]/R[SiteIndex]*Ry[SiteIndex]       ;
      dzEfield[3*index+2] = -de_gij*Rz[SiteIndex]/R[SiteIndex]*Rz[SiteIndex] - gij ;

      int three=3;

      Grad[(j)*3]   += -0.5 * 2.0* ddot(&three, &dxEfield[index*3], &one, &mu[index*3], &one) ;
      Grad[(j)*3+1] += -0.5 * 2.0* ddot(&three, &dyEfield[index*3], &one, &mu[index*3], &one) ;
      Grad[(j)*3+2] += -0.5 * 2.0* ddot(&three, &dzEfield[index*3], &one, &mu[index*3], &one) ;

       index++; 
      } 
   }
}

void Potential::MuCrossMu (int nAtoms, double *mu , dVec& tmu , dVec& mu_cross_mu, double wavefn)
{

 
 for (int i=0; i< nAtoms*3; ++i) {
    double  t0 = wavefn*wavefn*mu[i];
    tmu[i] += t0;
    for (int j=0; j< nAtoms*3; ++j)
      mu_cross_mu[i*nAtoms*3+j] += t0*mu[j];
 }

}

void Potential::FinalGradient( int nAtoms, double *Gradient, double *mu,  double  *mu_cross_mu, 
                               double *dT_x, double *dT_y, double *dT_z, double *dEfield )
{

   int nSites = nAtoms/3*4;
   int nAtoms3 = nAtoms*3;
   int nAtoms9 = nAtoms*9;
   int one = 1;

   for (int j = 0; j < nSites ; ++j) {
      for (int dim = 0; dim < 3 ; ++dim) {
          Gradient[(j)*3 + dim]   += -0.5 * 2.0* ddot(&nAtoms3, &dEfield[nAtoms*3* ((j)*3+dim) ], &one, &mu[0], &one) ;
      }

      int atom_j = 3*int(j/4) + j%4;
      if ((j+1)%4 != 0 ) {
          Gradient[j*3]   += -0.5 * 2.0* ddot(&nAtoms9, &(dT_x[atom_j*nAtoms*9]), &one, &mu_cross_mu[atom_j*nAtoms*9], &one) ;
          Gradient[j*3+1] += -0.5 * 2.0* ddot(&nAtoms9, &(dT_y[atom_j*nAtoms*9]), &one, &mu_cross_mu[atom_j*nAtoms*9], &one) ;
          Gradient[j*3+2] += -0.5 * 2.0* ddot(&nAtoms9, &(dT_z[atom_j*nAtoms*9]), &one, &mu_cross_mu[atom_j*nAtoms*9], &one) ;
      }
   

   }
   

}

void Potential::SubtractWWGradient (int nSites, double *PolGrad, double *Gradient)
{
//    printf(">>Before Water-Water Gradient subtracted:\n");
//    for (int j = 0; j < nSites*3 ; j+=3) 
//    printf("  % 20.10E   % 20.10E   % 20.10E\n", Gradient[j+0], Gradient[j+1], Gradient[j+2]);
//    printf(">>END\n");

   for (int j = 0; j < nSites ; ++j) {
      Gradient[j*3]   -= PolGrad[j*3] ;
      Gradient[j*3+1]   -= PolGrad[j*3+1] ;
      Gradient[j*3+2]   -= PolGrad[j*3+2] ;
   }

}


/*

void Potential::EvaluateDPP6SPGradient(
                                       const double *x,
                                       dVec& Grad
                                       ) 
{


   //
   //  potential due to point charges (damped with 0.2)
   //
   double Vpc = 0;
   const double charge_damping_with = 0.2;
   const double chargeexp = 1.0 / (charge_damping_with * charge_damping_with);
   for (int i = 0; i < nCharges; ++i) {
      int iqs = ChargeSite[i];
      Vpc += -Charge[i]/R[iqs] * (1.0 - exp(-chargeexp * R2[iqs]));
      double de1 = Charge[i]/R2[iqs]*(1.0 - exp(-chargeexp * R2[iqs]));
      double de2 = -Charge[i]/R[iqs] *(2.0*R[iqs]*chargeexp * exp(-chargeexp * R2[iqs]));
      double de3 = (de1+de2)/R[iqs];
      Grad[iqs*3] += -de3*Rx[iqs];
      Grad[iqs*3+1] += -de3*Ry[iqs];
      Grad[iqs*3+2] += -de3*Rz[iqs];
   }

   // cout<<"============================= "<<endl;

   //
   //  potential due to point dipoles (damped with 0.4)
   //
   //



   // don't use static! this function must be thread-safe
   // an optimization might be possible using an OpenMP call, but the performance penalty is not that bad yet

   double Vind = 0;
   const double dipole_damping_with = 0.4;
   const double dipoleexp = 1.0 / (dipole_damping_with * dipole_damping_with);
   for (int i = 0; i < nDipoles; ++i) {
      int ids = DipoleSite[i];
      const double *mu = &Dipole[3*i];
      double damp = 1.0 - exp(-dipoleexp * R2[ids]);
      Vind += -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2]) * Rminus3[ids] * damp * damp;

      double de1 = -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2])*(-3.0)/(R2[ids]*R2[ids])*damp*damp;
      double de2 = -(Rx[ids]*mu[0] + Ry[ids]*mu[1] + Rz[ids]*mu[2])* Rminus3[ids] * 2.0 * damp * (2.0 * R[ids]*dipoleexp * (1-damp));
      double de3 = -damp*damp * Rminus3[ids];
      double de4 = (de1 + de2) / R[ids];

      Grad[ids*3]  +=   - de4*Rx[ids] - de3*mu[0];
      Grad[ids*3+1] +=  - de4*Ry[ids] - de3*mu[1];
      Grad[ids*3+2] +=  - de4*Rz[ids] - de3*mu[2];

      RVec[i*3] = de3*Rx[ids];
      RVec[i*3+1] = de3*Ry[ids];
      RVec[i*3+2] = de3*Rz[ids];
   }


   char transa='T';
   double alpha = 1.0, beta = 0.0;
   int incx=1, incy =1;
   int order = nSites/4*3;
   int order2 = nSites;


   //   cout<<"aDmuByDx = "<<DmuByDx[0]<<endl;

   dgemv(&transa, &order, &order2, &alpha, &DmuByDx[0], &order, &RVec[0], &incx, &beta, &X[0], &incy);
   dgemv(&transa, &order, &order2, &alpha, &DmuByDy[0], &order, &RVec[0], &incx, &beta, &Y[0], &incy);
   dgemv(&transa, &order, &order2, &alpha, &DmuByDz[0], &order, &RVec[0], &incx, &beta, &Z[0], &incy);

   for (int i = 0; i < nSites; ++i) {
      Grad[i*3] += X[i];
      Grad[i*3+1] += Y[i];
      Grad[i*3+2] += Z[i];
   }





   //
   //  repulsive potential expressed as a sum over s-type GFs
   //
   double Vrep = 0;
   for (int i = 0; i < nGauss; ++i) {
      int igs = GaussSite[i];
      Vrep +=  Gauss[2*i+1] * exp(-Gauss[2*i] * R2[igs]);
      double de1 = Gauss[2*i+1]*(-2.0*Gauss[2*i]*R[igs])*exp(-Gauss[2*i] * R2[igs]);
      Grad[igs*3] += -de1*Rx[igs]/R[igs];
      Grad[igs*3+1] += -de1*Ry[igs]/R[igs];
      Grad[igs*3+2] += -de1*Rz[igs]/R[igs];
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

      double de1 = 2.0 * damp*(2.0*bdamping*exp(-bdamping*R2[ips]))/(R2[ips] * R[ips]);
      double de2 = -4.0*damp*damp/ (R2[ips] * R2[ips]*R[ips]);
      double de3= (de1+de2)/R[ips];

      Grad[ips*3] +=  0.5* alpha*de3*Rx[ips];
      Grad[ips*3+1] += 0.5* alpha*de3*Ry[ips];
      Grad[ips*3+2] += 0.5* alpha*de3*Rz[ips];


      Vpol += alpha * damp * damp / (R2[ips] * R2[ips]);
   }
   Vpol *= -0.5;


   //
   // done
   //


}
*/


