#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 

#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;

#include "timer.hpp"

#include "vecdefs.h"
#include "lapackblas.h"

//
//  this is for setting up a molecule with interacting atomic polarizability
//  InvA is the matrix needed to compute the induced dipoles from the "external" E field
//
void ComputeInvA(int nSites, const double *R, const double *alpha, double aThole, double *InvA, int verbose);
void BuildPolarizationMatrix(int nSites, const double *R, const double *alpha, double aThole, double *A);
void BuildTensor(int nSites, const double *R, const double *alpha, double aThole, double *Tensor);
void SCF_inducedDipole(int nSites, const double *Tensor, const double *alpha, double aThole, double *Efield, double *mu);
void CalcDDTensor(double *T, double *Rij, double alpha_i, double alpha_j, double aThole);
void InvertMatrix(int n, double *A);

//////////////////////////////////////////////////////////////
///
///   \f$A^{-1} = (\alpha^{-1} - T)^{-1}\f$ 
///
///   this is the matrix that can be multiplied with the electric field at
///   the polarizable sites to yield the induced dipoles:
///
///   \f$\mu = A^{-1} * E\f
///
///   dim(A) = 3*nSites x 3*nSites
///
void ComputeInvA(int nSites, const double *R, const double *alpha, double aThole, double *InvA, int verbose)
{


  BuildPolarizationMatrix(nSites, R, alpha, aThole, InvA);  
  InvertMatrix(3*nSites, InvA);

  // compute molecular polarizability 
  if (verbose > 1) {
    int n = nSites;
    int n3 = 3*n; // dimension of super-matrix A and leading dimension of any block B
    double am[9];
    cout << "\nMolecular polarizability tensor in Bohr^3\n";
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
	am[3*a+b] = 0;
	for (int i = 0; i < n; ++i) {
	  for (int j = 0; j < n; ++j) {
	    am[3*a+b] += InvA[n3*(3*i+a)+(3*j+b)];
	  }
	}
	cout << "  " << am[3*a+b];
      }
      cout << "\n";
    }
    cout << "Average = " << (am[0]+am[4]+am[8])/3.0 << "\n";

    cout << "\nMolecular polarizability tensor in Angs^3\n";
    const double Bohr2Angs = 0.52917720859;
    const double fac = Bohr2Angs * Bohr2Angs * Bohr2Angs;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
	cout << "  " << am[3*a+b]*fac ;
      }
      cout << "\n";
    }
    cout << "Average = " << fac*(am[0]+am[4]+am[8])/3.0 << "\n";

  }
}


//////////////////////////////////////////////////////////////
///
///  build matrix \f$A = (\alpha^{-1} - T)\f$
///
///  this is a super matrix
///  each diagonal 3x3 block has \f$\alpha^{-1}\f$ on its diagonal
///  each off-diagonal 3x3 block is the dipole-dipole tensor \f$T_{ij}\f$
///  as each \f$T_{ij}\f$ block is symmetric, the whole matrix is symmetric 
///
///   dim(A) = 3*nSites x 3*nSites
///
///   input:  nSites        : no of point polarizable atoms
///           R[3*nSites]   : position of sites
///           alpha[nSites] : atomic polarizabilities
///           aThole        : Thole's damping parameter
///   output: A       
///
///   calls: CalcDDTensor() 
///
void BuildPolarizationMatrix(int nSites, const double *R, const double *alpha, double aThole, double *A)
{

  int n = nSites;
  int n3 = 3*n; // dimension of super-matrix A and leading dimension of any block B


  for (int k = 0; k < n3*n3; ++k)
    A[k] = 0.0;

  // super-diagonal: (alpha)^-1 
  for (int i = 0; i < n; ++i) {
    double *B = A + n3*3*i + 3*i;  // pointer to current 3x3 block 
    B[0] = B[n3+1] = B[2*n3+2] = 1.0 / alpha[i]; 
    B[1] = B[2] = B[n3+0] = B[n3+2] = B[2*n3+0] = B[2*n3+1] = 0.0;
  }

  // diagonal blocks: Tij
  double rij[3];
  double T[9];

  cout<<"nSites: "<<nSites<<endl;

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      double *B  = A + n3*3*j + 3*i;  // pointer to current 3x3 block (ldb = n)
      double *Bt = A + n3*3*i + 3*j;  // and to the symmetry-related block
      rij[0] = R[3*i+0] - R[3*j+0];
      rij[1] = R[3*i+1] - R[3*j+1];
      rij[2] = R[3*i+2] - R[3*j+2];
      CalcDDTensor(T, rij, alpha[i], alpha[j], aThole);
      for (int a = 0; a < 3; ++a)
	for (int b = 0; b < 3; ++b) 
	  B[a*n3+b] = Bt[a*n3+b] = -T[3*a+b];
    }
  }
}


void BuildTensor(int nSites, const double *R, const double *alpha, double aThole, double *Tensor)
{

  int n = nSites;
  int n3 = 3*n; // dimension of super-matrix A and leading dimension of any block B

  double rij[3];
  double T[9];
  int iter=0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
       if (j == i)
         continue;

      rij[0] = R[3*i+0] - R[3*j+0];
      rij[1] = R[3*i+1] - R[3*j+1];
      rij[2] = R[3*i+2] - R[3*j+2];
      CalcDDTensor(T, rij, alpha[i], alpha[j], aThole);
      for (int a = 0; a < 3; ++a)
	for (int b = 0; b < 3; ++b) {
            Tensor[9*iter+3*a+b]=T[3*a+b];
        }
      iter++;
    }
  }
}


void SCF_inducedDipole(int nSites, const double *Tensor, const double *alpha, double aThole, double *Efield, double *mu)
{
  int n = nSites;
  int n3 = 3*n; // dimension of super-matrix A and leading dimension of any block B

  // diagonal blocks: Tij
  double rij[3];
  double T[9];
  double *field2 = new double[n3];
  double *dipole = new double[n3];
  double *olddipole = new double[n3];
  double newfac = 0.5;
  double oldfac = 1 - newfac;
  double tester = 0;
  int converged = 0;
  int iter ;
  int iteration;
  int max_iterations = 100;  // I just made this up
   double thresh = 1e-5;


  for (int k = 0; k < n3; ++k) {
    dipole[k] = 0;
    dipole[k] = newfac * alpha[k/3] * Efield[k];
   // cout<<"alpha[k/3]:"<<alpha[k/3]<<"  "<<Efield[k]<<endl;
 }



  for (iteration = 0; iteration < max_iterations; ++ iteration)
 {

  iter = 0;
  for (int i = 0; i < n3; ++i) 
    field2[i] = 0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
       if (j == i)
         continue;

      for (int a = 0; a < 3; ++a)
	for (int b = 0; b < 3; ++b) {
            field2[i*3+a] += Tensor[9*iter+3*a+b]*dipole[j*3+b];
           // cout<<i<<" "<<j<<" : T="<<Tensor[9*iter+3*a+b]<<" d="<<dipole[j*3+b]<<endl;          
	   //cout<<" i : field2["<<i*3+a<<"]: "<<"dipole["<<j*3+b<<"]"<<endl;
	   //cout<<" j : field2["<<j*3+a<<"]: "<<"dipole["<<j*3+b<<"]"<<endl;
	   //cout<<" Tensor["<<9*iter+ 3*a+b<<"]: "<<Tensor[9*iter+3*a+b]<<endl;
	}
       iter++;
    }
  }

   tester = 0;

  for (int i = 0; i < n3; ++i) {
     olddipole[i] = dipole[i];
     dipole[i] = oldfac*olddipole[i] + newfac * alpha[i/3] * (Efield[i] + field2[i]);
     //cout<<" field2["<<i<<"]= "<<field2[i]<<endl;
     //cout<<" olddipole:"<<olddipole[i]<<" newdipole "<< dipole[i]<<endl;
     double t = olddipole[i] - dipole[i];
      tester += t*t;
  }
  // cout<<" ------------------------------------------------------ "<<endl;
  
  //check convergence
    tester = sqrt(tester)/3/nSites;  // this is the average change
    // cout<<"tester ="<<tester<<" thresh:"<<thresh<<endl;
    if (tester < thresh) {
        converged = 1;
        break;
    }

  }

  for (int i = 0; i < n3; ++i) 
     mu[i] = dipole[i];

   if (converged == 0) {
       cout << "**** Warning in WaterCluster::CalcInducedDipoles():";
       cout << " after " << max_iterations << " iterations the induced dipoles converged only to "
          << tester << endl;
    }
    else
       cout << "Induced dipoles converged to " << thresh << "after "<<iteration<<endl;


  delete[] field2; delete[] dipole; delete[] olddipole;

}


/////////////////////////////////////////////////
///
///  symmetric matrix A is inverted in place
///
void InvertMatrix(int n, double *A)
{

  int verbose=1;

  int N = n;
  char uplo = 'U';
  static iVec ipiv; ipiv.resize(n);    // this is called repeatedly to invert 9x9 matrices, but only during setup, so
  static dVec work; work.resize(n*n);  // static vectors should be fine here as long as nobody calls this in mutiple threads
  int lwork = n*n;
  int info;

  dsytrf("U", &N, A, &n, &ipiv[0], &work[0], &lwork, &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dsytrf = " << info << "\n";
    exit(1);
  }

  dsytri("U", &N, A, &n, &ipiv[0], &work[0], &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dsytri = " << info << "\n";

    exit(1);
  }

  for (int i = 0; i < n; ++i)
    for (int j = i+1; j < n; ++j)
      A[n*i+j] = A[i+n*j];

}


/////////////////////////////////////////////////////////////////////////
/// 
///  compute tensor for dipole-dipole interaction
/// 
///  T is stored in fortran convention, that is, by column, T[0]=T11, T[1]=T21
///  0 3 6
///  1 4 7
///  2 5 8
///
///  Rij[3] is the vector from dipole-j to site-i (where it acts on dipole-i)  
///  OORij3 = 1/r_ij^3
///  OORij5 = 1/r_ij^5
///  f3 and f5 are Thole-type damping functions
///
///  
void CalcDDTensor(double *T, double *Rij, double alpha_i, double alpha_j, double aThole)
{

  double rij2 = Rij[0]*Rij[0] + Rij[1]*Rij[1] + Rij[2]*Rij[2];
  double rij  = sqrt(rij2);
  double rij3 = rij * rij2;
  double rij5 = rij2 * rij3;

  double TholeExp = aThole * rij3 / sqrt(alpha_i * alpha_j);
  double TholeFac = exp(-TholeExp);
  double f3 = 1 - TholeFac;
  double f5 = 1 - (1+TholeExp) * TholeFac;
  
  double fac1 = f5 * 3.0 / rij5;
  double fac2 = f3 / rij3;

  T[0] = fac1 * Rij[0]*Rij[0] - fac2;
  T[4] = fac1 * Rij[1]*Rij[1] - fac2;
  T[8] = fac1 * Rij[2]*Rij[2] - fac2;

  T[1] = T[3] = fac1 * Rij[0]*Rij[1];
  T[2] = T[6] = fac1 * Rij[0]*Rij[2];
  T[5] = T[7] = fac1 * Rij[1]*Rij[2];

  //printf("Second: Tij(%14.7f)= %14.7f %14.7f %14.7f %14.7f %14.7f %14.7f\n", rij, T[0], T[4], T[8], T[1], T[2], T[5]);

}







