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

#include <mpi.h>

#include "constants.h"
#include "Potential.h"
#include "DVR.h"
#include "lapackblas.h"
#include "timer.hpp"


//#define HAVE_LAPACK

//#ifndef HAVE_LAPACK
#define HAVE_ARPACK
#ifndef HAVE_ARPACK

int DVR::larnoldi(int ng, int nev, int maxsub, int maxiter, int ptol, double *ev)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if(rank==0)std::cout << "Error in DVR::larnoldi(): the ARPACK library is not available.\n"
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

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

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
  
 if(rank==0) progress_timer t("larnoldi", verbose);

  DiagCount++;

   //  if (verbose > 0) {
  if (verbose > 0) {
    if(rank==0)printf("Diagonalizing using the ARPACK\n");
   }
    switch (lmode)
      {
      case 1:
      if(verbose >0)
        if(rank==0)printf(" (normal Lanczos mode)\n");
        which = sa;
        break;
      case 3:
      if(verbose >0){
        if(rank==0)printf(" (shift-inverted mode)\n");
        if(rank==0)printf(" Shift = %f\n", shift);
      }
        which = lm;
        if(rank==0)printf(" (sorry, not yet implemented)\n");
        exit(1);
	break;
      default:
        if(rank==0)printf(" unknown lmode %d; this should not happen", lmode);
        exit(1);
      }
   if (verbose > 0) {
    if(rank==0)printf(" requested eigenvalues %d\n", nev);
    if(rank==0)printf(" maximal subspace size %d\n", ncv);
    if(rank==0)printf(" max no. of iterations %d\n", maxiter);
  }
   //  }

//if (DiagCount != 1) ncv=20;

  // arrays for the ARPACK
  int lworkl = ncv * (ncv + 8);
  static iVec select; select.resize(ncv);       // the bloody, allegedly not referenced array
  static dVec v; v.resize(ncv*ng*2);              // Lanczos basis
  static dVec workd; workd.resize(3*ng*2);        // Lanczos vectors
  static dVec workl; workl.resize(lworkl*2);   // work space
  static dVec resid; resid.resize(ng*2);       // residual vector

//  int iparam[11], ipntr[14];
  iparam[0] = 1;        //  exact shifts 
  iparam[2] = maxiter;
  iparam[3] = 1;        //  blocksize (only 1 works)
  iparam[6] = lmode;

  int nconv = 0;
  int nthreads, tid;
  int ido = 0;          // for the reverse communication interface



if (DiagCount == 1) {

  SaveV.resize(ncv*ng);              // Lanczos basis
  SaveWorkd.resize(3*ng);        // Lanczos vectors
  SaveWorkl.resize(lworkl);   // work space
  SaveResid.resize(ng);
  std::copy(&wavefn[0], &wavefn[ng-1], resid.begin());
} 
else {
  info=1;
  ido =0;
  std::copy(&wavefn[0], &wavefn[ng-1], resid.begin());

}

   if(rank==0)cout<<"DiagCount = "<<DiagCount<<endl;

  if(rank==0)std::cout << "dvrtype: " << dvrtype << std::endl;

  // only for shift-inverted   int n_solve = 0;      // no of time solve has been called
  int n_mtx = 0;
  int icount = 0;

  while (ido != 99) {
    fflush(stdout);

    icount++;
  
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
            if(dvrtype == 3) {
               VectorFFT(&workd[ipntr[0]-1], &workd[ipntr[1]-1]);
            }
	    else
               MatrixTimesVector(&workd[ipntr[0]-1], &workd[ipntr[1]-1]); 
            n_mtx ++;
            if (verbose > 0) {
              if (n_mtx %  10 == 0 && rank==0)printf(".");
              if (n_mtx % 500 == 0 && rank==0)printf("\n");
            }
            break;
          case 3:
            if(rank==0)printf("Case 3. This should not happen\n");
	    exit(31);
	  }
	break;
      case 99:
      //  convN++;
	break;
      case 199:
	break;
      default:
	if(rank==0)printf("IDO = %d, INFO = %d; what is happening?\n", ido, info);
        exit(1);
      }

if (DiagCount == 1 && icount%ncv ==2 ) {
  std::copy(&resid[0], &resid[ng-1], SaveResid.begin());
  std::copy(&v[0], &v[ncv*ng-1], SaveV.begin());
  std::copy(&workd[0], &workd[3*ng-1], SaveWorkd.begin());
  std::copy(&workl[0], &workl[lworkl-1], SaveWorkl.begin());
  for (int i=0; i < 11; i++) SaveIparam[i]=iparam[i];
//  for (int i=0; i < 14; i++) SaveIpntr[i]=ipntr[i];
}

 }


 if (info < 0) {
    if(rank==0)printf("INFO = %d; see dsaupd.f\n", info);
    exit(1);
  }
  

  nconv = iparam[4];
//  nconv= 0;
//  if(rank==0)cout<<" info= "<<info<<endl;

  if (verbose > 0) {
    if(rank==0)printf("\n");
    if(rank==0)printf(" Arnoldi updates       %d\n", iparam[2]);
    if(rank==0)printf(" y = OP * x operations %d\n", iparam[8]);
    if(rank==0)printf(" converged eigenvalues %d\n", iparam[4]);
  }


  int one;
  DSEUPD(&one, all, &select[0], ev, &v[0], &ng, &shift, 
	 bmat, &ng, which, &nev, &tol, &resid[0], &ncv, &v[0], &ng,
	 iparam, ipntr, &workd[0], &workl[0], &lworkl, &info);


  if (info < 0) {
    if(rank==0)printf("INFO = %d; see documentation in the ARPACK subroutine dseupd.f\n", info);
    exit(1);
  }

  /*   write to stdout  */
   if(rank==0)printf("\n");
   if (verbose > 0) {
    if(rank==0)printf("\nConverged Eigenvalues in Hartree and meV:\n");
    for (int i = 0; i < nconv; i++)
      if(rank==0)printf("%3d  %14.9f  %12.4f\n", i+1, ev[i], ev[i]*AU2MEV);
  }


  if (verbose > 1) {
    if(rank==0)printf("\nPseudo spectrum:\n       eigenvalue              error\n");
    if (nconv > 0)
      if(rank==0)printf("Careful; this information may be rubbish if nconverged=0:\n");
    for (int i = 0; i < ncv; i++) {
      if(rank==0)printf("%3d  %21.14e  %15.8e\n", i+1, workl[ipntr[7] + i - 1], workl[ipntr[8] + i - 1]);
    }
  }
  if(rank==0)cout<< "number of larnold iterations : "<<icount<<endl;
  std::copy(&v[0], &v[ng*nconv], wavefn.begin());

  return nconv;
}

#endif
