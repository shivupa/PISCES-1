///   
///    Driver for reverse-interface Davidson  
///
///    int ng      : number of grid points = matrix dimension
///    int maxsub  : maximal subspace size in subspace iterations
///    int maxiter : maximal no of macro-iterations (up to maxsub*maxiter mtx operations)
///    int ptol    : convergence threshold = 10^-ptol
///    int corrflag: which corrections vectors: 0, 1, or 2 see doc in Davidson.cpp
///    double *ev  : computed energy
///
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "Davidson.h"
#include "Potential.h"
#include "DVR.h"
#include "VectorFFT.hpp"



using namespace std;

/// keep a local copy of the class variables so this function may be useful outside of class DVR at some point 

int DVR::davdriver(int ng, int nstates, int maxsub, int maxiter, int ptol, int corrflag, double *ev)
{

  if (verbose > 0)
    cout << "Diagonalizing using the reverse-interface Davidson\n";

  static dVec B; B.resize(ng * (maxsub)); // basis for the subspace plus the residual vector (next correction vector)
  static dVec Z; Z.resize(ng * maxsub); // Hamilton matrix times basis vectors
  static dVec diag; diag.resize(ng);    // diagonal of H
  static dVec davwork;
  int workmem = DavidsonWorkSize(ng, maxsub, nstates, corrflag);
  davwork.resize(workmem);

  // copy start vectors to B 
  std::copy(&wavefn[0], &wavefn[nstates*ng], B.begin());

  //  compute diagonal
  ComputeDiagonal(&diag[0]);

  cout<<"define fft_engine"<<endl;
  class VectorFFT fft_engine(&n_1dbas[0]);


  //
  //  reverse-interface iteration
  // 
  int ido = 1;
  int n_mtx = 0;
  int inout[3] = {0, 0, 0};  // vector-block for mtx: "pointers" to B and Z, and no of vectors 
  int nConv = 0;
  while (ido != 0) {
    fflush(stdout);
    ido = Davidson(ng, maxsub, nstates, maxiter, ptol, corrflag, verbose, ev, nConv, 
		   &B[0], &Z[0], &diag[0], &davwork[0], inout);
    switch (ido)
      {
      case  1:
	// inefficient, for the block-Davidson mtx should call dgemm instead of dspmv
	// but then the block-Davidson is only for single-points, so what's the point
	for (int ivec = 0; ivec < inout[2]; ++ivec) {
	  int iB = inout[0]+ivec;
	  int iZ = inout[1]+ivec;
	  
          //TV: Calling the FFT
            if(dvrtype == 3)
              //  VectorFFT(&B[iB*ng], &Z[iZ*ng]);
               fft_engine.apply(&B[iB*ng], &Z[iZ*ng], &v_diag[0], &KE_diag[0]);
	    else
                MatrixTimesVector(&B[iB*ng], &Z[iZ*ng]);
	  n_mtx ++;
	}
	break;
      case 0:
	break;
      default:
	printf("ido = %d, this should not happen.\n", ido);
	exit(1);
      }
  }

  std::copy(&B[0], &B[nstates*ng], wavefn.begin());

  if (verbose > 0) {
    printf("-----------------------------------------------\n");
    printf("Davidson finishes after %i matrix-times-vector operations\n", n_mtx);
    printf("%i states have been converged.\n", nConv);
  }
  return nConv;
}
