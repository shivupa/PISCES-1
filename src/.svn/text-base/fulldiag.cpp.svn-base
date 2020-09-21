//
//  this is only for testing
//  
//  build the Hamilton matrix in hmat[ng*ng]
//  and diagonalise it with a LAPACK routine
//
#include <cstdio> 
#include <stdexcept>
#include <sstream>
#include "constants.h"
#include "lapackblas.h"

#include "DVR.h"


using std::printf;

void addt(int n1, int inc1, int nod, int ndod[], int incr[], double t[], double hop[], int ng);

////////////////////////////////////////////////////////////////////////
//
//  take v_diag and e_kin(i) to build the DVR Hamilitonian explicitly
//
//  this follows very much the function mtx()
//
void DVR::build_h(double *hmat)
{
   if (tformat != 0) throw std::runtime_error( 
      "So far fulldiag works only with triangles (tformat=0).\n" 
      );
   
   int k, id;

   // initialise with zero
   for (k = 0; k < ngp*ngp; ++k)
      hmat[k] = 0.0;
   
   // put the potential on the diagonal
   for (k = 0; k < ngp; ++k)
      hmat[k + ngp*k] = v_diag[k];
   
   // add each T(i) n times where n is the product of the no of grid points in the other dimensions
   for (id = 0; id < no_dim; id++) {
    k = n_1dbas[id]; n_1dbas[id] = n_1dbas[0]; n_1dbas[0] = k;
    k = incv[id]; incv[id] = incv[0]; incv[0] = k;
    addt (n_1dbas[0], incv[0], no_dim-1, &n_1dbas[1], &incv[1], e_kin[id], hmat, ngp);
    k = n_1dbas[id]; n_1dbas[id] = n_1dbas[0]; n_1dbas[0] = k;
    k = incv[id]; incv[id] = incv[0]; incv[0] = k;
  }
}

////////////////////////////////////////////////////////////////////
//
//  look into mtx.cpp titxbl() documentation
//  T(1) is added to H:    H = H + T(1)
//
void addt(int n1, int inc1, int nod, int ndod[], int incr[], double t[], double hop[], int ng)
{
  int iod, irc, nrc, incrc, k, i1, j1;
  int ir[20]; // this is MAXDIM
  double *tt;

  nrc = 1;
  for (iod = 0; iod < nod; iod++)
    nrc *= ndod[iod];

  ir[0] = -1;
  for (iod = 1; iod < nod; iod++) 
    ir[iod] = 0;

  for (irc = 0; irc < nrc; irc++) {
    for (k = 0; k < nod; k++) {
      ir[k] ++;
      if (ir[k] == ndod[k])
	ir[k] = 0;
      else
	break;
    }
    incrc = 0;
    for (k = 0; k < nod; k++)
      incrc += ir[k] * incr[k];
    {
      for (i1 = 0; i1 < n1; i1 ++) {
	int ig = incrc + i1*inc1;
	tt = &t[i1*(i1+1)/2];
	hop[ig + ng*ig] += tt[i1];
	for (j1 = 0; j1 < i1; j1 ++) {
	  int jg = incrc + j1*inc1;
	  hop[ig + ng * jg] = tt[j1];
	  hop[jg + ng * ig] = tt[j1];
	}
      }
    }
  }
}


///////////////////////////////////////////////////////
//
//  call the LAPACK to diagonalize the DVR of H
//
//
void DVR::fulldiag(double *hop)
{


  int info;
  int naux = 5*ngp;
  char job[] = "N";
  
  double *aux = new double[naux];  // aux mem for DSYEV
  double *d   = new double[ngp];    // eigenvalue vector

  dsyev(job, "U", &ngp, hop, &ngp, d, aux, &naux, &info);
  if (info != 0) {
     std::ostringstream error;
     error << "error calling dsyev in fulldiag: info = " << info;
     throw std::runtime_error( error.str() );
  }

  // print the lowest n eigenvalues
  int nev = 20;
  printf("\nEigenvalues in au and meV\n");
  for (int k = 0; k < nev; ++k)
    printf("%3d  %23.15e   %12.4f\n", k+1, d[k], d[k]*AU2MEV);

  delete[] d;
  delete[] aux;
}

