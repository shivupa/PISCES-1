//    
//     Driver for Lanczos-Arnoldi using the ARPACK  
//     
//     many variable names are consistent with the examples provided with the ARPACK library
//
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>


#include "constants.h"
#include "Potential.h"
#include "DVR.h"
#include "lapackblas.h"

//#define HAVE_LAPACK

//#ifndef HAVE_LAPACK
#define HAVE_ARPACK
#ifndef HAVE_ARPACK

int DVR::larnoldi(int ng, int nev, int maxsub, int maxiter, int ptol, double *ev)
{
  std::cout << "Error in DVR::larnoldi(): the ARPACK library is not available.\n"
       << "If you need this method, install it, and recompile pisces.\n"; 
  exit(0);
}

#else

#include <boost/fortran/prototype.hpp>

BOOST_FORTRAN_SUBROUTINE(DSAUPD, dsaupd, DSAUPD, 
                         (int*) // ido
                         (char*) // bmat
                         (int*) // n
                         (char*) // which
                         (int*) // nev
                         (double*) //tol
                         (double*) //resid
                         (int*) //ncv
                         (double*) //v
                         (int*) //ldv
                         (int*) //iparam
                         (int*) //ipntr
                         (double*) //workd
                         (double*) //workl
                         (int*) //lworkl
                         (int*) // info
                         );
BOOST_FORTRAN_SUBROUTINE(DSEUPD, dseupd, DSEUPD, 
                         (int*) //rvec
                         (char*) //howmny
                         (int*) //select
                         (double*) //ev
                         (double*) //z
                         (int*) //ldz
                         (double*) //sigma
                         (char*) // bmat
                         (int*)//n
                         (char*) //which
                         (int*)//nev
                         (double*)//tol
                         (double*)//resid
                         (int*)//ncv
                         (double*)//v
                         (int*)//ldv
                         (int*)//iparam
                         (int*)//ipntr
                         (double*)//workd
                         (double*)//workl
                         (int*)//lworkl
                         (int*)//info
                         );

int DVR::larnoldi(int ng, int nev, int maxsub, int maxiter, int ptol, double *ev)
{

  int ncv = maxsub; 
  int info = 1; // start vector is in wavefn
  int lmode = 1;
  double shift = 0; 
  double tol = std::pow(10.,-ptol); 
  // parameters for the ARPACK
  char *which;
  char sa[] = "SA";
  char lm[] = "LM";
  char bmat[] = "I";
  char all[] = "A";
  
   //  if (verbose > 0) {
  if (verbose > 0) {
    printf("Diagonalizing using the ARPACK\n");
   }
    switch (lmode)
      {
      case 1:
      if(verbose >0)
        printf(" (normal Lanczos mode)\n");
        which = sa;
        break;
      case 3:
      if(verbose >0){
        printf(" (shift-inverted mode)\n");
        printf(" Shift = %f\n", shift);
      }
        which = lm;
        printf(" (sorry, not yet implemented)\n");
        exit(1);
	break;
      default:
        printf(" unknown lmode %d; this should not happen", lmode);
        exit(1);
      }
   if (verbose > 0) {
    printf(" requested eigenvalues %d\n", nev);
    printf(" maximal subspace size %d\n", ncv);
    printf(" max no. of iterations %d\n", maxiter);
  }
   //  }

  // arrays for the ARPACK
  int lworkl = ncv * (ncv + 8);
  static iVec select; select.resize(ncv);       // the bloody, allegedly not referenced array
  static dVec v; v.resize(ncv*ng);              // Lanczos basis
  static dVec workd; workd.resize(3*ng);        // Lanczos vectors
  static dVec workl; workl.resize(lworkl);   // work space
  static dVec resid; resid.resize(ng);       // residual vector
  
  int iparam[11], ipntr[14]; 
  iparam[0] = 1;        //  exact shifts 
  iparam[2] = maxiter;
  iparam[3] = 1;        //  blocksize (only 1 works)
  iparam[6] = lmode;


  int nconv = 0;

  // copy startvector 
  std::copy(&wavefn[0], &wavefn[ng-1], resid.begin());

  
  // Start of Lanczos Iteration

  // only for shift-inverted   int n_solve = 0;      // no of time solve has been called
  int ido = 0;          // for the reverse communication interface
  int n_mtx = 0;
  while (ido != 99) {
    fflush(stdout);
    DSAUPD(&ido, bmat, &ng, which, &nev, &tol, &resid[0], &ncv, &v[0], &ng,
	   iparam, ipntr, &workd[0], &workl[0], &lworkl, &info);
    switch (ido)
      {
      case  1:
      case -1:
	switch (lmode)
	  {
          case 1:
            //TV: Calling the FFT
            if(dvrtype == 3)
                 VectorFFT(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
	    else
               MatrixTimesVector(&workd[ipntr[0]-1], &workd[ipntr[1]-1]); 
            n_mtx ++;
            if (verbose > 0) {
              if (n_mtx %  10 == 0) printf(".");
              if (n_mtx % 500 == 0) printf("\n");
            }
            break;
          case 3:
            printf("Case 3. This should not happen\n");
	    exit(31);
	    /* solve(ng, &workd[ipntr[1]-1], &workd[ipntr[0]-1], shift);
	    if (++n_solve > 49)
	      {printf(".\n"); n_solve = 0;}
	    else
	      printf(".");
	      break; */
	  }
	break;
      case 99:
	break;
      default:
	printf("IDO = %d, INFO = %d; what is happening?\n", ido, info);
        exit(1);
      }
  }
  // std::cout << "verbose: " << n_mtx << std::endl;

  if (info < 0) {
    printf("INFO = %d; see dsaupd.f\n", info);
    exit(1);
  }
  
  nconv = iparam[4];
  if (verbose > 0) {
    printf("\n");
    printf(" Arnoldi updates       %d\n", iparam[2]);
    printf(" y = OP * x operations %d\n", iparam[8]);
    printf(" converged eigenvalues %d\n", iparam[4]);
  }
  //if (lmode == 3)
  //  means();

  int one;
  DSEUPD(&one, all, &select[0], ev, &v[0], &ng, &shift, 
	 bmat, &ng, which, &nev, &tol, &resid[0], &ncv, &v[0], &ng,
	 iparam, ipntr, &workd[0], &workl[0], &lworkl, &info);
  if (info < 0) {
    printf("INFO = %d; see documentation in the ARPACK subroutine dseupd.f\n", info);
    exit(1);
  }

  /*   write to stdout  */
   printf("\n");
   if (verbose > 0) {
    printf("\nConverged Eigenvalues in Hartree and meV:\n");
    for (int i = 0; i < nconv; i++)
      printf("%3d  %14.9f  %12.4f\n", i+1, ev[i], ev[i]*AU2MEV);
  }

  if (verbose > 1) {
    printf("\nPseudo spectrum:\n       eigenvalue              error\n");
    if (nconv > 0)
      printf("Careful; this information may be rubbish if nconverged=0:\n");
    for (int i = 0; i < ncv; i++) {
      printf("%3d  %21.14e  %15.8e\n", i+1, workl[ipntr[7] + i - 1], workl[ipntr[8] + i - 1]);
    }
  }
  std::copy(&v[0], &v[ng*nconv], wavefn.begin());

  return nconv;
}

#endif
