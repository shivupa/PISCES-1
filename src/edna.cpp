/*
 *   dressed down version of the cube integrator
 *
 *   it reads an orbital evaluated on a grid and computes numerical integrals of
 *   <r^2>
 *   <r>
 *   <r_nearest>   where r_nearest is the distance to the nearest atom
 *
 *   It compiles with g++ and icpp, but should really compile with any C++ compiler
 *
 *   It is called from the command prompt: edna  input.cube 
 *   where input.cube is a Gaussian-cube file 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>

using namespace std;

const int LineLength = 1024;
const double Bohr2Angs = 0.52917720859;

double r_nearest(double &x, double &y, double &z, int &nnuc, double *xnuc, double *ynuc, double *znuc);


int main (int argc, char *argv[])
{

  int i;
  int verbose = 1;

  printf("\nCube integrator for EDNA\n");
  if (argc < 2 || argc > 3) {
    printf("call: edna cube-file [for_binning.dat]\n");
    exit(0);
  }

  FILE *output = 0;
  if (argc == 3)
    output = fopen(argv[2], "w+");
  else
    output = fopen("edna_bins.dat", "w+");
  fprintf(output,"# r, rho, rho_nearest,  rho, roh_nearest\n# cols 2,3 for plotting, cols 4,5 sum to 1\n");
  
  ifstream input;
  input.open(argv[1]);
  // skip two comment lines
  char buffer[LineLength+1];
  input.getline(buffer, LineLength);
  input.getline(buffer, LineLength);

  // read number of atoms and grid parameters
  int nnuc;
  int nx, ny, nz;
  double xmin, ymin, zmin, xx, yy, zz, dx, dy, dz;
  input >> nnuc >> xmin >> ymin >> zmin;
  input >> nx >> dx >> yy >> zz;
  if (yy != 0 || zz != 0) {cout << "Only cartesian grids please.\n"; exit(1);}
  input >> ny >> xx >> dy >> zz;
  if (xx != 0 || zz != 0) {cout << "Only cartesian grids please.\n"; exit(1);}
  input >> nz >> xx >> yy >> dz;
  if (xx != 0 || yy != 0) {cout << "Only cartesian grids please.\n"; exit(1);}
  long npts = nx * ny * nz;
  double vol = dx * dy * dz;
  if (verbose > 0) {
    cout << "\nGrid parameters:\n";
    cout << nx << " Points along x from " << xmin << " to " << xmin+(nx-1)*dx << " steplength = " << dx << endl;
    cout << ny << " Points along y from " << ymin << " to " << ymin+(ny-1)*dy << " steplength = " << dy << endl;
    cout << nz << " Points along z from " << zmin << " to " << zmin+(nz-1)*dz << " steplength = " << dz << endl;
  }

  // try to guess reasonable values for the width and number of bins
  double xmax = fabs(xmin + (nx-1)*dx);
  double ymax = fabs(ymin + (ny-1)*dx);
  double zmax = fabs(zmin + (nz-1)*dx);
  if (fabs(xmin) > xmax) xmax = fabs(xmin);
  if (fabs(ymin) > ymax) xmax = fabs(ymin);
  if (fabs(zmin) > zmax) xmax = fabs(zmin);
  double rmax = sqrt(xmax*xmax + ymax*ymax + zmax*zmax);
  double bin_width = sqrt(dx*dx+dy*dy+dz*dz);
  int nbins = int(rmax/bin_width); 
  double *rbins = new double[nbins+1];
  double *rnbins = new double[nbins+1];

  // read the nuclear coordinates
  double *rnuc = new double[nnuc];
  double *xnuc = new double[nnuc];
  double *ynuc = new double[nnuc];
  double *znuc = new double[nnuc];
  int charge;
  if (verbose > 0)
    cout << "\nNuclei\n";
  for (i = 0; i < nnuc; ++i) {
    input >> charge >> rnuc[i] >> xnuc[i] >> ynuc[i] >> znuc[i];
    if (verbose > 0)
      printf("%3i  %2i  %11.5f     %11.5f %11.5f %11.5f\n", i+1, charge, rnuc[i], xnuc[i], ynuc[i], znuc[i]); 
  }


  // 
  //  read the grid data and do your sums  
  //  sum either density enclosed in contours (JobFlag > 0)
  //  or density enclosed in spheres with different radii (JobFlag < 0)
  //

  double x, y, z, orbital, density; 
      
  double total_sum = 0;       // sum over the grid 
  double rsqexp = 0;          // sum for <r^2>
  double rabsexp = 0;         // sum for <|r|>
  double x_exp = 0, y_exp = 0, z_exp = 0; // for <r_vec>
  double rnexp = 0;           // sum for <r_nearest>
  double maximum = -1e10;     // maximum density
  double minimum = 1e10;      // minimum density

  long nvaluesread = 0;  
  for (int ix = 0; ix < nx; ++ix) {
    x = xmin + double(ix)*dx;
    for (int iy = 0; iy < ny; ++iy) {
      y = ymin + double(iy)*dy;
      for (int iz = 0; iz < nz; ++iz) {
	z = zmin + double(iz)*dz;
	// get next density
	input >> orbital;
	density = orbital*orbital;
	nvaluesread++;
	//  find maximal and minimal density and maximal edge value
	if (density > maximum)
	  maximum = density;
	if (density < minimum)
	  minimum = density;
	// integrate the density 
	total_sum += density;
	// <r^2> and <|r|>, and <r_nearest> integrals
	double rsq = x*x + y*y + z*z;
	double rpt = sqrt(rsq);
	double rnr = r_nearest(x, y, z, nnuc, xnuc, ynuc, znuc);
	rsqexp  += density * rsq;
	rabsexp += density * rpt;
	rnexp   += density * rnr;
	x_exp   += density * x;
	y_exp   += density * y;
	z_exp   += density * z;
	//
	int i_rbin  = int(rpt / bin_width);
	int i_rnbin = int(rnr / bin_width);
	if (i_rbin > nbins) i_rbin = nbins;
	if (i_rnbin > nbins) i_rnbin = nbins;
	rbins[i_rbin] += density;
	rnbins[i_rnbin] += density;
	// done; next grid point
      }
    }
  }

  // done with reading cube data
  input.close();

  printf("Last value read is %e (no %li)\n", orbital, nvaluesread);
  if (nvaluesread != npts)
    cout << "But there should be " << npts << endl;

  // accounting
  rsqexp *= vol;
  rabsexp *= vol;
  rnexp *= vol;
  total_sum *= vol;
  x_exp *= vol;
  y_exp *= vol;  
  z_exp *= vol;
  
  double r_vec_ex_sq = x_exp * x_exp + y_exp * y_exp + z_exp * z_exp;
  double r_gy  = sqrt(rsqexp - r_vec_ex_sq);
  double r_gy2 = sqrt(rsqexp - rabsexp*rabsexp);

  //  print summary  
  printf("Maximum density = %e\n", maximum);
  printf("Minimum density = %e\n", minimum);
  printf("\nIntegrals:\n");
  printf("  over the grid = %f\n", total_sum);
  printf("  x             = %f Bohr = %f Angs\n", x_exp , x_exp * Bohr2Angs);
  printf("  y             = %f Bohr = %f Angs\n", y_exp , y_exp * Bohr2Angs);
  printf("  z             = %f Bohr = %f Angs\n", z_exp , z_exp * Bohr2Angs);
  printf("  <|r|>         = %f Bohr = %f Angs\n", rabsexp , rabsexp * Bohr2Angs);
  printf("  sqrt(<r^2>)   = %f Bohr = %f Angs\n", sqrt(rsqexp), sqrt(rsqexp)*Bohr2Angs);
  printf("  r_gyration    = %f Bohr = %f Angs\n", r_gy, r_gy * Bohr2Angs);
  printf("  <r_nearest>   = %f Bohr = %f Angs\n", rnexp , rnexp  * Bohr2Angs);

  printf("\nSummary for grep:\n  %f  %f  %f  %f\n", 
	 total_sum, sqrt(rsqexp), rabsexp, rnexp);

  // write bins file
  for (int i = 0; i <= nbins; ++i)
    fprintf(output,"%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n", i*bin_width, 
	    rbins[i]*vol/bin_width, rnbins[i]*vol/bin_width,
	    rbins[i]*vol, rnbins[i]*vol);
  fclose(output);


  exit(0);
}



double r_nearest(double &x, double &y, double &z, int &nnuc, double *xnuc, double *ynuc, double *znuc)
{
  // finds the distance to the nucleus nearest to (x, y, z) 
  double dsqmin = (x-xnuc[0])*(x-xnuc[0]) + (y-ynuc[0])*(y-ynuc[0]) + (z-znuc[0])*(z-znuc[0]);

  for (int i = 1; i < nnuc; ++i){
    double dsq = (x-xnuc[i])*(x-xnuc[i]) + (y-ynuc[i])*(y-ynuc[i]) + (z-znuc[i])*(z-znuc[i]);
    if (dsq < dsqmin)
      dsqmin = dsq;
  }

  return sqrt(dsqmin);
}


