#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>


#include "constants.h"

using namespace std;

/////////////////////////////////////////
//
//  write a cube file
//  only cartesian cubes
//
void WriteCubeFile(const char *fname, int nAtoms, const int *Z, const double *position, const int *npts, const double *x0, const double *dx, const double *cube, const int *verbose)
{
  if (verbose > 0)
    cout << "Writing cube-file " << fname << "\n";

  int nx = npts[0];
  int ny = npts[1];
  int nz = npts[2];

  FILE *cubefile;
  cubefile = fopen(fname,"w+");


  // Gaussian cube file
  int ValuesPerLine = 6;
  double nought = 0.0;
  // two comment lines suitable for cubeint
  fprintf(cubefile, " 5 0  (these 2 lines are for cubeint)\n");
  fprintf(cubefile, " 0.01 0.001 0.0001 0.00001 0.000001\n");
  // header: no of atoms and definition of the grid
  fprintf(cubefile, "%5i  %11.6f  %11.6f  %11.6f\n", nAtoms, x0[0], x0[1], x0[2]);
  fprintf(cubefile, "%5i  %11.6f  %11.6f  %11.6f\n", nx, dx[0], nought, nought);
  fprintf(cubefile, "%5i  %11.6f  %11.6f  %11.6f\n", ny, nought, dx[1], nought);
  fprintf(cubefile, "%5i  %11.6f  %11.6f  %11.6f\n", nz, nought, nought, dx[2]);
  // atoms list: the 2nd number is ignored by most programs and usually 0.0
  // for cubeint it is set to the van der Waals radius
  for (int k = 0; k < nAtoms; ++k) {
    const double *r = position + 3*k;
    double RvdW = 0.0;
    switch (Z[k])
      {
      case 1: RvdW = 1.20; break;
      case 6: RvdW = 1.70; break;
      case 7: RvdW = 1.55; break;
      case 8: RvdW = 1.52; break;
      default: break; // do nothing; 
      }
    fprintf(cubefile, "   %i %11.6f  %11.6f  %11.6f  %11.6f\n", 
	    Z[k], RvdW*Angs2Bohr, r[0]*Angs2Bohr, r[1]*Angs2Bohr, r[2]*Angs2Bohr);
  }


  // here comes the cube
  int icube = 0;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      int iPerLine = 1;
      for (int iz = 0; iz < nz; ++iz) {
	fprintf(cubefile, "%13.6e ", cube[icube]);
	icube ++;
	if (iPerLine % ValuesPerLine == 0)
	  fprintf(cubefile, "\n");
	iPerLine++;
      }
      if (nz % ValuesPerLine != 0)
	fprintf(cubefile, "\n");
    }
  }

  fclose(cubefile);

}
