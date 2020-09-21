#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 

#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;

#include "vecdefs.h"
#include "lapackblas.h"

//
//  this is for setting up the derivative of the dipole_dipole interaction
//  tensor. Essentially, it will be a supermatrix of dimension 
//  nPolSites x nPolSites. Each element "ij" would correspond to the derivative
//  of the tensor T(ij) and hence of dimensions 3x3x3. Due to the symmetry of 
//  the T(ij) element, only 6x3 matrix would suffice in principle
//
void BuildSuperdT(int nPolSites, const double *R, const double *alpha, double aThole, double *A_x ,double *A_y , double *A_z );
void CalcdT(double *TDerivXX , double *TDerivYY , double *TDerivZZ ,
                                double *TDerivXY , double *TDerivYZ , double *TDerivXZ ,
                                double *Rij, double alpha_i, double alpha_j, double aThole) ; 


void BuildSuperdT(int nPolSites, const double *R, const double *alpha, double aThole, double *A_x ,double *A_y , double *A_z )
{

  int n = nPolSites;
// dimension of super-matrix A would be n x n x 27 

  for (int k = 0; k < n*n*9; ++k)
{      A_x[k] = 0.0;
      A_y[k] = 0.0;
      A_z[k] = 0.0; }

  for (int j = 0; j < n; ++j) {
      for (int i =0 ; i < n ; ++i) {
        if (i != j){
          double Rij[3]; 
          double TDerivXX[3] ;
          double TDerivXY[3] ;
          double TDerivXZ[3] ;
          double TDerivYY[3] ;
          double TDerivYZ[3] ;
          double TDerivZZ[3] ; 
          Rij[0] = R[3*j+0] - R[3*i+0];
          Rij[1] = R[3*j+1] - R[3*i+1];
          Rij[2] = R[3*j+2] - R[3*i+2];

          CalcdT( TDerivXX , TDerivYY , TDerivZZ ,
                  TDerivXY , TDerivYZ , TDerivXZ ,
                  Rij, alpha[i], alpha[j], aThole) ; 

/*

          A_x[j*n*9+i*9+0]  = TDerivXX[0] ;
          A_x[j*n*9+i*9+1]  = TDerivXY[0] ;
          A_x[j*n*9+i*9+2]  = TDerivXZ[0] ;
          A_x[j*n*9+i*9+3]  = TDerivXY[0] ;
          A_x[j*n*9+i*9+4]  = TDerivYY[0] ;
          A_x[j*n*9+i*9+5]  = TDerivYZ[0] ;
          A_x[j*n*9+i*9+6]  = TDerivXZ[0] ;
          A_x[j*n*9+i*9+7]  = TDerivYZ[0] ;
          A_x[j*n*9+i*9+8]  = TDerivZZ[0] ;
          

          A_y[j*n*9+i*9+0] = TDerivXX[1] ;
          A_y[j*n*9+i*9+1] = TDerivXY[1] ;
          A_y[j*n*9+i*9+2] = TDerivXZ[1] ;
          A_y[j*n*9+i*9+3] = TDerivXY[1] ;
          A_y[j*n*9+i*9+4] = TDerivYY[1] ;
          A_y[j*n*9+i*9+5] = TDerivYZ[1] ;
          A_y[j*n*9+i*9+6] = TDerivXZ[1] ;
          A_y[j*n*9+i*9+7] = TDerivYZ[1] ;
          A_y[j*n*9+i*9+8] = TDerivZZ[1] ;

      
          A_z[j*n*9+i*9+0] = TDerivXX[2] ;
          A_z[j*n*9+i*9+1] = TDerivXY[2] ;
          A_z[j*n*9+i*9+2] = TDerivXZ[2] ;
          A_z[j*n*9+i*9+3] = TDerivXY[2] ;
          A_z[j*n*9+i*9+4] = TDerivYY[2] ;
          A_z[j*n*9+i*9+5] = TDerivYZ[2] ;
          A_z[j*n*9+i*9+6] = TDerivXZ[2] ;
          A_z[j*n*9+i*9+7] = TDerivYZ[2] ;
          A_z[j*n*9+i*9+8] = TDerivZZ[2] ;

*/
          A_x[j*n*9+       i*3+0]  = TDerivXX[0] ;
          A_x[j*n*9+       i*3+1]  = TDerivXY[0] ;
          A_x[j*n*9+       i*3+2]  = TDerivXZ[0] ;
          A_x[j*n*9+ n*3+  i*3+0]  = TDerivXY[0] ;
          A_x[j*n*9+ n*3+  i*3+1]  = TDerivYY[0] ;
          A_x[j*n*9+ n*3+  i*3+2]  = TDerivYZ[0] ;
          A_x[j*n*9+ n*3*2+i*3+0]  = TDerivXZ[0] ;
          A_x[j*n*9+ n*3*2+i*3+1]  = TDerivYZ[0] ;
          A_x[j*n*9+ n*3*2+i*3+2]  = TDerivZZ[0] ;

          A_y[j*n*9+       i*3+0]  = TDerivXX[1] ;
          A_y[j*n*9+       i*3+1]  = TDerivXY[1] ;
          A_y[j*n*9+       i*3+2]  = TDerivXZ[1] ;
          A_y[j*n*9+ n*3+  i*3+0]  = TDerivXY[1] ;
          A_y[j*n*9+ n*3+  i*3+1]  = TDerivYY[1] ;
          A_y[j*n*9+ n*3+  i*3+2]  = TDerivYZ[1] ;
          A_y[j*n*9+ n*3*2+i*3+0]  = TDerivXZ[1] ;
          A_y[j*n*9+ n*3*2+i*3+1]  = TDerivYZ[1] ;
          A_y[j*n*9+ n*3*2+i*3+2]  = TDerivZZ[1] ;

          A_z[j*n*9+       i*3+0]  = TDerivXX[2] ;
          A_z[j*n*9+       i*3+1]  = TDerivXY[2] ;
          A_z[j*n*9+       i*3+2]  = TDerivXZ[2] ;
          A_z[j*n*9+ n*3+  i*3+0]  = TDerivXY[2] ;
          A_z[j*n*9+ n*3+  i*3+1]  = TDerivYY[2] ;
          A_z[j*n*9+ n*3+  i*3+2]  = TDerivYZ[2] ;
          A_z[j*n*9+ n*3*2+i*3+0]  = TDerivXZ[2] ;
          A_z[j*n*9+ n*3*2+i*3+1]  = TDerivYZ[2] ;
          A_z[j*n*9+ n*3*2+i*3+2]  = TDerivZZ[2] ;





/*
          A_x[j*n*9+       i*3+0]  = TDerivXX[0] ;
          A_x[j*n*9+ n*3+  i*3+0]  = TDerivXY[0] ;
          A_x[j*n*9+ n*3*2+i*3+0]  = TDerivXZ[0] ;
          A_x[j*n*9+       i*3+1]  = TDerivXY[0] ;
          A_x[j*n*9+ n*3+  i*3+1]  = TDerivYY[0] ;
          A_x[j*n*9+ n*3*2+i*3+1]  = TDerivYZ[0] ;
          A_x[j*n*9+       i*3+2]  = TDerivXZ[0] ;
          A_x[j*n*9+ n*3+  i*3+2]  = TDerivYZ[0] ;
          A_x[j*n*9+ n*3*2+i*3+2]  = TDerivZZ[0] ;
          

          A_y[j*n*9+       i*3+0] = TDerivXX[1] ;
          A_y[j*n*9+ n*3+  i*3+0] = TDerivXY[1] ;
          A_y[j*n*9+ n*3*2+i*3+0] = TDerivXZ[1] ;
          A_y[j*n*9+       i*3+1] = TDerivXY[1] ;
          A_y[j*n*9+ n*3+  i*3+1] = TDerivYY[1] ;
          A_y[j*n*9+ n*3*2+i*3+1] = TDerivYZ[1] ;
          A_y[j*n*9+       i*3+2] = TDerivXZ[1] ;
          A_y[j*n*9+ n*3+  i*3+2] = TDerivYZ[1] ;
          A_y[j*n*9+ n*3*2+i*3+2] = TDerivZZ[1] ;

      
          A_z[j*n*9+       i*3+0] = TDerivXX[2] ;
          A_z[j*n*9+ n*3+  i*3+0] = TDerivXY[2] ;
          A_z[j*n*9+ n*3*2+i*3+0] = TDerivXZ[2] ;
          A_z[j*n*9+       i*3+1] = TDerivXY[2] ;
          A_z[j*n*9+ n*3+  i*3+1] = TDerivYY[2] ;
          A_z[j*n*9+ n*3*2+i*3+1] = TDerivYZ[2] ;
          A_z[j*n*9+       i*3+2] = TDerivXZ[2] ;
          A_z[j*n*9+ n*3+  i*3+2] = TDerivYZ[2] ;
          A_z[j*n*9+ n*3*2+i*3+2] = TDerivZZ[2] ;

*/ 

/*
          for (int m=0 ; m < 3 ; ++m){
          cout << TDerivXX[m] << endl ;
          cout << TDerivXY[m] << endl ;
          cout << TDerivXZ[m] << endl ;
          cout << TDerivXY[m] << endl ;
          cout << TDerivYY[m] << endl ;
          cout << TDerivYZ[m] << endl ;
          cout << TDerivXZ[m] << endl ;
          cout << TDerivYZ[m] << endl ;
          cout << TDerivZZ[m] << endl ;         
          }
*/
 
        }
      }
   }

//          exit(1) ; 
}

    
void CalcdT(double *TDerivXX , double *TDerivYY , double *TDerivZZ ,
                                double *TDerivXY , double *TDerivYZ , double *TDerivXZ ,
                                double *Rij, double alpha_i, double alpha_j, double aThole)
{
  double rij2 = Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2];
  double rij  = sqrt(rij2);
  double rij3 = rij * rij2;
  double rij5 = rij2 * rij3;
  double rij7 = rij2 * rij5;

  double TholeExp = aThole * rij3 / sqrt(alpha_i * alpha_j);
  double TholeFac = exp(-TholeExp);
  double DerivTholeExp = aThole * 3.0*rij / sqrt(alpha_i * alpha_j);
  double DerivTholeFac = TholeFac *(-1.0)*DerivTholeExp;

  double f3 = 1 - TholeFac;
  double f5 = 1 - (1+TholeExp) * TholeFac;
  double f3Deriv = -DerivTholeFac ;
  double f5Deriv = - (DerivTholeExp*TholeFac + (1+TholeExp)*DerivTholeFac) ;


  double fac1 = f5 * 3.0 / rij5;
  double fac2 = f3 / rij3;

  double fac1Deriv = f5Deriv*3.0/rij5 + f5*3.0*(-5.0)/rij7 ;
  double fac1DerivComp[3] ;
         fac1DerivComp[0] =  fac1Deriv*Rij[0] ;
         fac1DerivComp[1] =  fac1Deriv*Rij[1] ;
         fac1DerivComp[2] =  fac1Deriv*Rij[2] ;

  double fac2Deriv = f3Deriv/rij3 + f3*(-3.0)/rij5 ;
  double fac2DerivComp[3] ;
         fac2DerivComp[0] =  fac2Deriv*Rij[0] ;
         fac2DerivComp[1] =  fac2Deriv*Rij[1] ;
         fac2DerivComp[2] =  fac2Deriv*Rij[2] ;



  TDerivXX[0] = fac1DerivComp[0]*Rij[0]*Rij[0] + fac1*2.0*Rij[0] - fac2DerivComp[0] ;
  TDerivXX[1] = fac1DerivComp[1]*Rij[0]*Rij[0]                   - fac2DerivComp[1] ;
  TDerivXX[2] = fac1DerivComp[2]*Rij[0]*Rij[0]                   - fac2DerivComp[2] ;


  TDerivYY[0] = fac1DerivComp[0]*Rij[1]*Rij[1]                   - fac2DerivComp[0] ;
  TDerivYY[1] = fac1DerivComp[1]*Rij[1]*Rij[1] + fac1*2.0*Rij[1] - fac2DerivComp[1] ;
  TDerivYY[2] = fac1DerivComp[2]*Rij[1]*Rij[1]                   - fac2DerivComp[2] ;

  TDerivZZ[0] = fac1DerivComp[0]*Rij[2]*Rij[2]                   - fac2DerivComp[0] ;
  TDerivZZ[1] = fac1DerivComp[1]*Rij[2]*Rij[2]                   - fac2DerivComp[1] ;
  TDerivZZ[2] = fac1DerivComp[2]*Rij[2]*Rij[2] + fac1*2.0*Rij[2] - fac2DerivComp[2] ;


  TDerivXY[0] = fac1DerivComp[0]*Rij[0]*Rij[1] + fac1*Rij[1] ;
  TDerivXY[1] = fac1DerivComp[1]*Rij[0]*Rij[1] + fac1*Rij[0] ;
  TDerivXY[2] = fac1DerivComp[2]*Rij[0]*Rij[1]               ;

  TDerivYZ[0] = fac1DerivComp[0]*Rij[1]*Rij[2]               ;
  TDerivYZ[1] = fac1DerivComp[1]*Rij[1]*Rij[2] + fac1*Rij[2] ;
  TDerivYZ[2] = fac1DerivComp[2]*Rij[1]*Rij[2] + fac1*Rij[1] ;

  TDerivXZ[0] = fac1DerivComp[0]*Rij[0]*Rij[2] + fac1*Rij[2] ;
  TDerivXZ[1] = fac1DerivComp[1]*Rij[0]*Rij[2]               ;
  TDerivXZ[2] = fac1DerivComp[2]*Rij[0]*Rij[2] + fac1*Rij[0] ;
}





