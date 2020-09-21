//
//  write files with cuts through the potential and through the wavefunctions
//
//  all this works only for 3D,  there are many more cuts in nD
//
//
//
//  for each 3D function we write three files: .XY, .XZ, and .YZ  = cuts in these planes
//
//  the filenames are POTENTIAL, WAVEFN1, ....
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include "constants.h"

using namespace std;

void writewfcuts(int nwf, int ng, int no_dim, int *n_1dbas, int max1db, 
		 double *x_dvr, int *incv, double *v, double *wf, int cubeflag, int nWaters, double *WaterPos)
{
  if (no_dim != 3) {
    cout << "Error in writewfcuts; this is a function for 3D grids only\n";
    exit(1);
  }

  int nx = n_1dbas[0];
  int ny = n_1dbas[1];
  int nz = n_1dbas[2];
  int ix, iy, iz;
  FILE *flout;

  double *xgrid = x_dvr;
  double *ygrid = x_dvr+max1db;
  double *zgrid = x_dvr+2*max1db;


  iz = nz/2; cout << "Cut XY at the " << iz << "-th z-grid point\n"; 
  flout = fopen("POTENTIAL.XY","w+");
  for (iy = 0; iy < ny; ++iy) {
    for (ix = 0; ix < nx; ++ix) {
      fprintf(flout, "%10.7f %10.7f %15.7e\n", xgrid[ix], ygrid[iy], v[ix*incv[0]+iy*incv[1]+iz*incv[2]]);
    }
    fprintf(flout,"\n");
  }
  fclose(flout);

  iy = ny/2; cout << "Cut XZ at the " << iy << "-th y-grid point\n"; 
  flout = fopen("POTENTIAL.XZ","w+");
  for (iz = 0; iz < nz; ++iz) {
    for (ix = 0; ix < nx; ++ix) {
      fprintf(flout, "%10.7f %10.7f %15.7e\n", xgrid[ix], zgrid[iz], v[ix*incv[0]+iy*incv[1]+iz*incv[2]]);
    }
    fprintf(flout,"\n");
  }
  fclose(flout);

  ix = nx/2; cout << "Cut YZ at the " << ix << "-th x-grid point\n"; 
  flout = fopen("POTENTIAL.YZ","w+");
  for (iy = 0; iy < ny; ++iy) {
    for (iz = 0; iz < nz; ++iz) {
      fprintf(flout, "%10.7f %10.7f %15.7e\n", ygrid[iy], zgrid[iz], v[ix*incv[0]+iy*incv[1]+iz*incv[2]]);
    }
    fprintf(flout,"\n");
  }
  fclose(flout);


  // now the wavefunctions
  //
  //  this is a DVR wavefunction, so the volume element is already in the
  //  value of the wavefunction at that grid point
  //  in the grid-integrater we have dV=dxdydz  dx = (xmax-xmin)/(nx-1)
  //
  //  good check: the sum over a DVR grid should always be 1.00000000
  //
  //  multiplication with 1/sqrt(dV) should give the Bohr^(-3/2) unit of the wavefunction
  //

  double dV = (zgrid[nz-1]-zgrid[0])*(ygrid[nz-1]-ygrid[0])*(xgrid[nz-1]-xgrid[0]) / double((nx-1)*(ny-1)*(nz-1));
  double oosqrdv = 1.0/sqrt(dV);

  cout << "Cube normalization factor is " << oosqrdv << "\n";

  for (int iwf = 0; iwf < nwf; ++iwf) { 

    char fname[15];
    double *wfp = wf+iwf*ng;

    // cuts
    if (0) {
    iz = nz/2; 
    sprintf(fname,"WaveFn%i.XY", iwf+1); flout = fopen(fname,"w+");
    for (iy = 0; iy < ny; ++iy) {
      for (ix = 0; ix < nx; ++ix) {
	fprintf(flout, "%10.7f %10.7f %15.7e\n", xgrid[ix], ygrid[iy],  wfp[ix*incv[0]+iy*incv[1]+iz*incv[2]]*oosqrdv);
      }
      fprintf(flout,"\n");
    }
    fclose(flout);
    iy = ny/2;
    sprintf(fname,"WaveFn%i.XZ", iwf+1); flout = fopen(fname,"w+");
    for (iz = 0; iz < nz; ++iz) {
      for (ix = 0; ix < nx; ++ix) {
	fprintf(flout, "%10.7f %10.7f %15.7e\n", xgrid[ix], zgrid[iz],  wfp[ix*incv[0]+iy*incv[1]+iz*incv[2]]*oosqrdv);
      }
      fprintf(flout,"\n");
    }
    fclose(flout);
    ix = nx/2;
    sprintf(fname,"WaveFn%i.YZ", iwf+1); flout = fopen(fname,"w+");
    for (iy = 0; iy < ny; ++iy) {
      for (iz = 0; iz < nz; ++iz) {
	fprintf(flout, "%10.7f %10.7f %15.7e\n", ygrid[iy], zgrid[iz],  wfp[ix*incv[0]+iy*incv[1]+iz*incv[2]]*oosqrdv);
      }
      fprintf(flout,"\n");
    }
    fclose(flout);
    }

    //
    // 1=gOpenMol 2=Gaussian cube file format
    // all this does not work for HO DVR (this should be tested in dvr3d!)
    //
    double rho = 0, intr = 0;
    sprintf(fname,"WaveFn%i.cub", iwf+1); flout = fopen(fname,"w+");
    cout << "Writing the cube file " << fname << " for wave function " << iwf+1 << endl;
    if (cubeflag == 1) 
      {
	// gOpenMol
	fprintf(flout, "3 3\n%i %i %i\n", nz, ny, nx);
	const double B2A = 0.529177249;
	fprintf(flout, "%13.6e %13.6e    %13.6e %13.6e    %13.6e %13.6e\n", 
		zgrid[0]*B2A, zgrid[nz-1]*B2A, ygrid[0]*B2A, ygrid[ny-1]*B2A, xgrid[0]*B2A, xgrid[nx-1]*B2A);
	for (iz = 0; iz < nz; ++iz) {
	  for (iy = 0; iy < ny; ++iy) {
	    for (ix = 0; ix < nx; ++ix) {
	      rho = wfp[ix*incv[0]+iy*incv[1]+iz*incv[2]] * oosqrdv;
	      fprintf(flout, "%13.6e\n", rho);
	      intr += rho*rho;
	    }
	  }
	}
      }
    else 
      {
	// header of a Gaussian cube file
	int ValuesPerLine = 6;
	double nought = 0.0;
	double dx = xgrid[1] - xgrid[0];
	double dy = ygrid[1] - ygrid[0];
	double dz = zgrid[1] - zgrid[0];
	fprintf(flout, " 5 0\n"); // this is for cubeint
	fprintf(flout, " 0.01 0.001 0.0001 0.00001 0.000001\n");
	fprintf(flout, "%5i  %11.6f  %11.6f  %11.6f\n", 3*nWaters, xgrid[0], ygrid[0], zgrid[0]);
	fprintf(flout, "%5i  %11.6f  %11.6f  %11.6f\n", nx, dx, nought, nought);
	fprintf(flout, "%5i  %11.6f  %11.6f  %11.6f\n", ny, nought, dy, nought);
	fprintf(flout, "%5i  %11.6f  %11.6f  %11.6f\n", nz, nought, nought, dz);
	// 1.52 and 1.2 are the van der Waals radii of O and H (relevant only for cubeint, in Gaussian ignored, I think)
	for (int k = 0; k < nWaters; ++k) {
	  int add = 9*k;
	  fprintf(flout, "   8 %11.6f  %11.6f  %11.6f  %11.6f\n", 1.52*Angs2Bohr,
		  WaterPos[add+0]*Angs2Bohr, WaterPos[add+1]*Angs2Bohr, WaterPos[add+2]*Angs2Bohr);
	  fprintf(flout, "   1 %11.6f  %11.6f  %11.6f  %11.6f\n", 1.20*Angs2Bohr,
		  WaterPos[add+3]*Angs2Bohr, WaterPos[add+4]*Angs2Bohr, WaterPos[add+5]*Angs2Bohr);
	  fprintf(flout, "   1 %11.6f  %11.6f  %11.6f  %11.6f\n", 1.20*Angs2Bohr,
		  WaterPos[add+6]*Angs2Bohr, WaterPos[add+7]*Angs2Bohr, WaterPos[add+8]*Angs2Bohr);
	}
	// here is the cube
	for (ix = 0; ix < nx; ++ix) {
	  for (iy = 0; iy < ny; ++iy) {
	    int iPerLine = 1;
	    for (iz = 0; iz < nz; ++iz) {
	      rho = wfp[ix*incv[0]+iy*incv[1]+iz*incv[2]] * oosqrdv;
	      fprintf(flout, "%13.6e ", rho);
	      if (iPerLine % ValuesPerLine == 0)
		fprintf(flout, "\n");
	      iPerLine++;
	      intr += rho*rho;
	    }
	    if (nz % ValuesPerLine != 0)
	      fprintf(flout, "\n");
	  }
	}
      }
    printf("  Int d3r rho(r) = %11.9f\n", intr*dV);
  }
  

}
