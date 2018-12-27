//
//  set up an harmonic oscillator DVR
//
//  calculates kinetic energy matrix t and grid points x of the DVR basis
//
//   input:
//
//     ndim : no. of HO basis functions
//     w    : mass scaled frequency of the harmonic oscillator (mass is one)
//     ldt  : flag for the output format 
//
//   output:
//
//     x[ndim]      : gridpoints
//     t, ldt       : ldt > 0   t is a full matrix
//                    ldt = 0   t is a lower triangle
//     z, itz       : eigenvector matrix needed to transform and back-transform
//     itz = 0        temporatry allocate z (itz is NOT used as a leading dimension)
//
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lapackblas.h"

#define SQR(a) ((a)*(a))
#define AUWN  219474.0




void ho_dvr (int ndim, double w, double *x, double *t, int ldt, double *zz, int itz, int verbose)
{
  int i, j, k, ic, info;
  int one = 1;
  double zwi, *z, *aux;

  aux = new double[ndim*ndim]; // work space for dsteqr
  if (itz == 0)
    z = new double[ndim*ndim]; // eigenvectors of X
  else 
    z = zz;

  /*  find the grid
   *  by diagonalising the X operator  (X is tridiagonal in the energy eigenfunctions)
   * 
   *                                          FBR
   *  t is used to store the offdiagonal of  X       
   */
  for (i = 0; i < ndim; ++i) {
    x[i] = 0.0;
    t[i] = sqrt((double)(i+1));
  }
  dsteqr("I", &ndim, x, t, z, &ndim, aux, &info);
  if (info != 0) {printf("error calling dsteqr in ho_dvr: info = %i\n", info); exit(1);}

  zwi = 1.0 / sqrt(2.0 * w);
  for (i = 0; i < ndim; ++i)
    x[i] *= zwi;
  if (verbose > 0) 
    printf(" w0  %f (%f cm-1)\n xmax  %f\n", w, w*AUWN, x[ndim-1]);

  // set the phase of the eigenvectors uniquely
  for (i = 0; i < ndim; ++i)
    if (z[i*ndim] < 0)
      for (k = 0; k < ndim; ++k)
	z[i*ndim+k] = -z[i*ndim+k];


  /*  set up T = 0.5 p^2 in t
   *   step #1  T is band diagonal in energy representation
   *   step #2  aux = T * Z
   *   step #3  T   = Z+ * aux
   */              
  for (i = 0; i < ndim; i++) {
    t[i] = 0.5 * w * ((double)i + 0.5);
    t[ndim+i] = -0.25 * w * sqrt((double)((i+1)*(i+2)));
  }
  for (i = 0; i < ndim; i++) 
    for (j = 0; j < ndim; j++) {
      k = i + ndim*j;
      aux[k] = t[i] * z[k];
      if (i < ndim-2) aux[k] += t[ndim+i]   * z[k+2]; 
      if (i > 1)      aux[k] += t[ndim+i-2] * z[k-2]; 
    }
  for (i = 0; i < ndim; i++)
    for (j = 0; j <= i; j++) {
      if (ldt > 0)
	ic = i + ldt*j;
      else
	ic = i*(i+1)/2 + j;
      t[ic] = ddot(&ndim, &z[ndim*i], &one, &aux[ndim*j], &one);  
    }
  if (ldt > 0)
    for (i = 0; i < ndim; i++)
      for (j = 0; j < i; j++)
	t[j + i*ldt] = t[i + j*ldt];
  
  delete[] aux;
  if (itz == 0)
    delete[] z;
}
