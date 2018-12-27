//
//  sine DVR 
//
//  sine DVR is based on the particle-in-a-box wavefunctions
//  this leads to an equidistant grid
//  the no of grid points is N (x1, ... xN), and the wavefunction 
//  vanishes by construction at x0 and x(N+1)
//  
//  the kinetic energy in DVR is known analytically
//  so no matrix diagonalization is needed to set up the DVR
//
//  the formulae are taken from Physics Reports 324, 1 (2000) appendix B
//
//  
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SQR(a) ((a)*(a))
#ifndef PI
#define PI 3.14159265358979
#endif


//////////////////////////////
//
//   sine_dvr
//
//  input:
//    ngp     : no of grid points
//    mu      : mass of the particle
//    x0      : left wall  :  leftmost gridpoint is x1
//    xN+1    : right wall :  rightmost gridpoint is xN
//    ldt     : leading dimension of ekin (if ldt = 0, ekin is stored as a lower triangle)
//    verbose : output level 
//
//  output:
//    x       : grid points (x[0] ... x[ngp-1] are x1 ... xN) 
//    ekin    : DVR of the kinetic energy operator   
//
//
void sine_dvr(int ngp, double mu, double x0, double xNp1, double *x, double *ekin, int ldt, int verbose)
{
  double DeltaX = (xNp1 - x0) / (double)(ngp+1);   // L = x(N+1) - x0

  if (verbose > 0) {
    printf("Sine DVR with %d grid points\n", ngp);
    printf("  mass: %f\n", mu);
    printf("  box walls at %10.6f  -  %10.6f,  L = %10.6f\n", x0, xNp1, xNp1-x0);
    printf("  grid range   %10.6f  -  %10.6f\n", x0+DeltaX, xNp1-DeltaX);
  }

  //
  //  compute the grid points; eq. (B.62)
  //  note that x0 and xN+1 are not grid points
  //
  int j;
  for (j = 0; j < ngp; ++j)
    x[j] = x0 + DeltaX * (double)(j+1);
  if (verbose > 1) {
    printf("grid points:\n");
    for (int j = 0; j < ngp; ++j)
      printf("%3i  %12.6f\n", j+1, x[j]);
  }

  //
  //  compute kinetic energy (eq. (B.65))
  //
  double f1 = 0.5 * PI*PI / (mu * SQR(DeltaX));
  double f2 = PI / (double)(ngp+1);
  double f3 = SQR(ngp+1);
  int a, b, ic;

  //
  //  diagonal terms
  //
  for (a = 0; a < ngp; ++a) {
    if (ldt > 0)
      ic = a + ldt*a;
    else
      ic = a*(a+1)/2 + a;
    ekin[ic] = f1 * ( 1.0/3.0 + 1.0/(6.0*f3) - 0.5/(f3 * SQR(sin((a+1)*f2))) );
  }
  if (verbose > 2) {
    printf("\nTii diagonal terms\n");
    for (a = 0; a < ngp; ++a) {
      if (ldt > 0)
	ic = a + ldt*a;
      else
	ic = a*(a+1)/2 + a;
      printf("%3i  %15.7f\n", a+1, ekin[ic]);
    }
  }

  //
  //   off diagonal terms
  //
  for (a = 0; a < ngp; ++a) {
    for (b = 0; b < a; ++b) {
      if (ldt > 0)
	ic = a + ldt*b;
      else
	ic = a*(a+1)/2 + b;
      ekin[ic] = f1 * 2.0/f3 * sin(f2*(a+1)) * sin(f2*(b+1)) / SQR(cos(f2*(a+1)) - cos(f2*(b+1))) ;
      if ((a-b)%2 == 1)
	ekin[ic] *= -1.0;
    }
  }
  //  fill full matrix if ldt > 0
  if (ldt > 0)
    for (a = 0; a < ngp; ++a)
      for (b = 0; b < a; ++b)
	ekin[b + a*ldt] = ekin[a + b*ldt];

  if (verbose > 3) {
    printf("\nTij off-diagonal terms\n");    
    for (a = 0; a < ngp; ++a) {
      for (b = 0; b < a; ++b) {
	if (ldt > 0)
	  ic = a + ldt*b;
	else
	  ic = a*(a+1)/2 + b;
	printf("%3i %3i %15.7e\n", a+1, b+1, ekin[ic]);
      }
    }
  }
    
}
