#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define SQR(a) ((a)*(a))
#ifndef PI
#define PI 3.14159265358979
#endif

void cm_dvr (int m, double len, double mass, double t[], int ldt);

/*
    compute the kinetic energy matrix and the grid points
    of a Colbert Miller DVR (this should be sine dvr)

    nb : no. of grid points
    mu : mass of the particle
    rl, ru : radial range of the grid
    
*/


void cm_ham (int nb, double mu, double rl, double ru, double *ekin, double *r, int ldt, int verbose)
{
  int i;
  double length = ru - rl;
  double dr = length / (double) nb;

  for (i = 0; i < nb; ++i)
    r[i] = rl + dr * (double)i;

  printf("#  CM-DVR with %d grid points\n", nb);
  printf("#  mass: %f\n", mu);
  printf("#  grid range: %f  -  %f\n", rl, ru);
  cm_dvr(nb, length, mu, ekin, ldt);
}

/*   
     calculate T in Colbert Miller DVR basis
     
     m     : no. of basis fn / grid points
     len   : length of the grid
     mass  : mass in au

     t, ldt : ldt > 0   t is a full matrix
            : ldt = 0   t is a lower triangle

     copied from the fortran by Andreas J.

*/

void cm_dvr (int m, double len, double mass, double t[], int ldt)
{

  int n = m + 1;
  double f1 = SQR(PI) / (4.0 * mass * SQR(len));
  double f2 = PI / (double)(2*n);
  double f3 = (2*n*n + 1.0)/3.0;

  int i, j, ic;

  for (i = 0; i < m; i++) {
    if (ldt > 0)
      ic = i + ldt*i;
    else
      ic = i*(i+1)/2 + i;
    t[ic] = f1 * (f3 - 1.0/SQR(sin(PI*(double)(i+1)/(double)n)));
    for (j = 0; j < i; j++) {
      if (ldt > 0)
	ic = i + ldt*j;
      else
	ic = i*(i+1)/2 + j;
      t[ic] = f1 * (1.0/SQR(sin(f2*(j-i))) - 1.0/SQR(sin(f2*(j+i+2))));
      if ((i-j)%2 == 1)
	t[ic] *= -1.0;
    }
  }

  if (ldt > 0)
    for (i = 0; i < m; i++)
      for (j = 0; j < i; j++)
	t[j + i*ldt] = t[i + j*ldt];

}
