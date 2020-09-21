#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 

#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;
#include "ChargeDipPol.h"
#include "vecdefs.h"
#include "lapackblas.h"
#include "constants.h"

//
//  this is for setting up a molecule with interacting atomic polarizability
//  InvA is the matrix needed to compute the induced dipoles from the "external" E field
//

//////////////////////////////////////////////////////////////
//
//   InvA = (alpha^-1 - T)^-1 
//
//   this is the matrix that can be multiplied with the electric field at
//   the polarizable sites to yield the induced dipoles:
//
//   mu = InvA * E
//
//   dim(A) = 3*nSiyes x 3*nSites
//
void InvChargeDipolePol(int nSites, const double *R, double *InvA, int nFullerenes, const int *NoAtomsArray)
{

  BuildPolarizationMatrix(nSites, R, InvA, nFullerenes, NoAtomsArray);  
  InvertCDMatrix(4*nSites+nFullerenes, InvA);

}


//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray)
{

  int n = nSites;
  int n4 = n*4 ;
  int n4plus1 = 4*n + nFullerenes ; // dimension of super-matrix A and leading dimension of any block B

  for (int k = 0; k < n4plus1*n4plus1 ; ++k)
    A[k] = 0.0;
  
//  double R_const = 0.62*1.889725989 ; 
  double R_const = 0.68620399*1.889725989 ; 
  double chi = sqrt(2.0/PI)/R_const  ;  
  double Rqq = R_const*sqrt(2.0)  ;
  //cout << "R_const , R_qq" << R_const << " "<<Rqq << endl ; 
  //Tqq
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (i == j) 
        {A[n4plus1*i+j] = chi; }
        //cout << "A" << n4plus1*i+j << "= " << A[n4plus1*i+j] << endl; }
      else
        {double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ; }
        //cout << "A" << n4plus1*i+j << "= " << A[n4plus1*i+j] << endl; }
    }
  }
/*
  //Tlagrangei
  cout <<"  warning, no restricted charge flow " << endl ; 
  for (int i = 0; i < n; ++i){
    A[i*n4plus1+n4plus1-1] = 1.0 ;
    A[i+n4plus1*(n4plus1-1)] = 1.0 ;
  //    cout << "A" << i*n4plus1+n4plus1-1 << "= " << A[i*n4plus1+n4plus1-1] << endl;
  //    cout << "A" << i+n4plus1*(n4plus1-1) << "= " << A[i+n4plus1*(n4plus1-1)] << endl; 
  }
  A[n4plus1*n4plus1 -1]= 0.0 ;
*/     

  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  for ( int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom      = current_atom ;
    current_atom += NoAtomsArray[iFuller] ;
    for (int i = old_atom; i < current_atom ; ++i) {
      A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
    }
  }
  
  double alpha = 3.0*sqrt(PI/2)*R_const*R_const*R_const ;
  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){
      if(i==j){
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = 1.0/alpha ;
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = 1.0/alpha ;
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = 1.0/alpha ;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq ; 
        double Txy, Txz, Tyz ;
        fac = erf(Rij/Rqq)-2.0*Rij*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq) ;
        fac_pq = fac/(Rij*Rij*Rij) ; 
        fac1 = fac_pq/(Rij*Rij) ; 
        fac2 = 4.0*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq*Rqq*Rqq*Rij*Rij) ;
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+1] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+2] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2] = -1.0*Tyz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1] = -1.0*Tyz ;
         
        //Tpq
        A[i*n4plus1+n+3*j+0] = rij[0]*fac_pq*-1.0;
        A[i*n4plus1+n+3*j+1] = rij[1]*fac_pq*-1.0;
        A[i*n4plus1+n+3*j+2] = rij[2]*fac_pq*-1.0;
        A[n*n4plus1+i*3*n4plus1                  +j] = rij[0]*fac_pq*-1.0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1         +j] = rij[1]*fac_pq*-1.0;
        A[n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j] = rij[2]*fac_pq*-1.0;
 
      }
     
    }
  }
  //for (int i = 0; i < n4plus1*n4plus1; ++i){
  //  cout << "A" << i << "=" << A[i] << endl;
  //}

}

/////////////////////////////////////////////////
//
//  symmetric matrix A is inverted in place
//
void InvertCDMatrix(int n, double *A)
{

  int N = n;
  char uplo = 'U';
  static iVec ipiv; ipiv.resize(n);    // this is called repeatedly to invert 9x9 matrices, but only during setup, so
  static dVec work; work.resize(n*n);  // static vectors should be fine here as long as nobody calls this in mutiple threads
  int lwork = n*n;
  int info;

  dgetrf(&N, &N, A, &n, &ipiv[0], &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrf = " << info << "\n";
    exit(1);
  }

  dgetri(&N, A, &n, &ipiv[0], &work[0],  &lwork, &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrfi = " << info << "\n";
    exit(1);
  }
}


