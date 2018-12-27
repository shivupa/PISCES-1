#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>

#include <mpi.h>

#include "vecdefs.h"
#include "Parameters.h"
#include "tsin.h"
#include "GetInput.h"


using namespace std;


////////////////////////////////////////////////////////
//
//  Read, and print parameters read from the following input groups:
//       Job, WaterModel, LocalPolarization, and SoftBox input groups
//
void GetInputParameters(Parameters &P, TSIN Input)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Job
  P.runtype    = Input.GetInt("Job", "runtype", 0);          
  P.nElectron  = Input.GetInt("Job", "nElectron", 1); // can be 0 or 1
  P.CubeFile   = Input.GetInt("Job", "CubeFile", 0); 
  P.WfCuts     = Input.GetInt("Job", "WfCuts", 0); 
  
  // WaterModel group
  P.WMverbose  = Input.GetInt("WaterModel", "Verbose", 1);
  P.nWater     = Input.GetInt("WaterModel", "NoOfWaters", 0);         // no of water monomers
  P.KTFlag     = Input.GetInt("WaterModel", "KTCharges", 0);
  P.CenterFlag = Input.GetInt("WaterModel", "CoMOrigin", 0);

  // ElectronPotential
  if (P.nElectron > 0) {
    P.PotFlag[0] = Input.GetInt("ElectronPotential", "Potential", 3); 
    P.PotFlag[1] = Input.GetInt("ElectronPotential", "DampType", 2);
    P.PotFlag[2] = Input.GetInt("ElectronPotential", "SigmaOFlag", 0);
    P.PotFlag[3] = Input.GetInt("ElectronPotential", "Polarization", 1);
    P.PotFlag[4] = Input.GetInt("ElectronPotential", "AdiabaticPolPot", 0);
    P.PotFlag[5] = Input.GetInt("ElectronPotential", "InternalParam", 0);
    P.PotVerbose = Input.GetInt("ElectronPotential", "Verbose", 0);
    switch (P.PotFlag[0])
      {
      case 1:  // DPP-6S_GTO-P
      case 2:  // DPP-12S_GTO-P
      case 3:  // DPP-4S_STO-P
      case 4:  // DPP-P without explicit repulsive potential
	P.PotPara[0] = Input.GetDouble("ElectronPotential", "Alpha", 9.745);
	P.PotPara[1] = Input.GetDouble("ElectronPotential", "PolDamping", 2.536);
	P.PotPara[2] = Input.GetDouble("ElectronPotential", "RepCoreScale", 4.4837);
	P.PotPara[3] = Input.GetDouble("ElectronPotential", "CoulombDamping", 1.7724);
	P.PotPara[4] = Input.GetDouble("ElectronPotential", "DipoleDamping", 1.7724);
	P.PotPara[5] = Input.GetDouble("ElectronPotential", "ORepCoreScale", 4.4837);
	P.PotPara[6] = Input.GetDouble("ElectronPotential", "EpsDrude", 0.5 / sqrt(P.PotPara[0]));
	P.PotPara[10] = Input.GetDouble("ElectronPotential", "Rtol", 20.0 );
	break;
      case 11:
	P.PotPara[0] = Input.GetDouble("ElectronPotential", "Alpha", 9.745);
	P.PotPara[1] = Input.GetDouble("ElectronPotential", "PolDamping", 0.288);
	P.PotPara[2] = Input.GetDouble("ElectronPotential", "CoulombDamping", 0.6);
	P.PotPara[3] = Input.GetDouble("ElectronPotential", "VRepA1H", 0.41);  // this is so far not used
	P.PotPara[4] = Input.GetDouble("ElectronPotential", "VRepB1O", 0.620);
	P.PotPara[5] = Input.GetDouble("ElectronPotential", "VRepB1H", 0.150);
	P.PotPara[6] = Input.GetDouble("ElectronPotential", "VRepB2O", 1.0);
	P.PotPara[7] = Input.GetDouble("ElectronPotential", "VRepB2H", 0.5);
	P.PotPara[8] = Input.GetDouble("ElectronPotential", "VRepB3O", 0.4);
	P.PotPara[9] = Input.GetDouble("ElectronPotential", "VRepB3H", 0.35);
	break;
      case 1021:
        break ; 
      case 1023:
        break ; 
      default:
	if(rank == 0) cout << "GetInputParameters: unknown potential for the excess electron: " << P.PotFlag[0] << "\n";
	exit(1);
      }
    
    if(P.PotFlag[5] == 1){   //#1 means default parameters are turned on
       P.PotFlag[0] = 3 ;
       P.PotFlag[1] = 2 ;
       if(rank == 0)cout << "InternalParameters uses Damptype= " << P.PotFlag[1] << ", and \n" ;
       if(rank == 0)cout << "                  Repulsive Pot = " << P.PotFlag[1] << "\n" ;
        
       switch (P.PotFlag[3]) {   //polarization model type
         case 1:
           P.PotPara[1] = 5.011 ;
           P.PotPara[4] = P.PotPara[1] ; // thomas-vamsee decided to keep these equal
           P.PotPara[2] = 6.6396 ; 
           break ;
         case 2:
           P.PotPara[1] = 5.1543 ; 
           P.PotPara[4] = P.PotPara[1] ; // thomas-vamsee decided to keep these equal
           P.PotPara[2] = 6.6717 ;
           break ; 
         case 3:
           P.PotPara[1] = 5.005 ;
           P.PotPara[4] = P.PotPara[1] ; // thomas-vamsee decided to keep these equal
           P.PotPara[2] = 6.6382 ;
           break;
         case 4:
           P.PotPara[1] =  4.84;
           P.PotPara[4] = P.PotPara[1] ; // thomas-vamsee decided to keep these equal
           P.PotPara[2] = 6.602 ;
           break;
         default:
         if(rank == 0)cout << "GetInputParameters: this combination has unknown internal-parameters \n" ; 
         }
    }
  }

  // GridDef group
  P.gridverbose = Input.GetInt("GridDef", "Verbose", 0);
  P.DVRType = Input.GetInt("GridDef", "DVRType", 0);
  P.Sampling = Input.GetInt("GridDef", "Sampling", 1);
  Input.GetIntArray("GridDef", "NoOfGridPoints", P.ngrid, 3);
  switch (P.DVRType)
    {
    case 1:
      Input.GetDoubleArray("GridDef", "Omega0", P.gpara, 3);      
      break;
    default:
      Input.GetDoubleArray("GridDef", "Length", P.gpara, 3);
      break;
    }

  // Diag group
  P.diagverbose = Input.GetInt("Diag", "Verbose", 1);
  P.DiagMethod = Input.GetInt("Diag", "Method", 2);// 1=Lanczos  2=Davidson  3= 0th-order Jacobi-Davidson
  P.nStates = Input.GetInt("Diag", "nStates", 1);
  P.maxSub = Input.GetInt("Diag", "maxSubspace", 20);   // max no of micro-iterations
  P.maxIter = Input.GetInt("Diag", "maxIter", 100);      // max no of macro-iterations
  P.ptol = Input.GetInt("Diag", "pTol", 5);          // tolerance = 10^-pTol
  P.istartvec = Input.GetInt("Diag", "StartVector", 1); // this is different for Lanczos + Davidson and needs work

  // Optimize group
  if (P.runtype == 2) {
    P.gtol = Input.GetDouble("Optimize", "Gtol", 1e-4);
    P.optverbose = Input.GetDouble("Optimize", "Verbose", 0);
  }
  // Molecular Dynamics group
  if (P.runtype == 3) {
    P.nsteps = Input.GetInt("MolecularDynamics", "nsteps", 100);
    P.timestep = Input.GetDouble("MolecularDynamics", "timestep", 0.001);
    P.RigidBody = Input.GetInt("MolecularDynamics", "RigidBody", 1);

  }   




  // potfit
  P.nParaOpt =  Input.GetInt("PotFit", "nPara", 0);
  if (P.nParaOpt > 0) {
    P.LowCutOff  = Input.GetDouble("PotFit", "LowCutOff", 2.0); 
    P.HighCutOff = Input.GetDouble("PotFit", "HighCutOff", 15.0); 
    P.wexp = Input.GetDouble("PotFit", "WeightExp", -1.0);
    P.minimizer = Input.GetInt("PotFit", "Method", 0); 
    P.mapping.resize(P.nParaOpt);
    Input.GetIntArray("PotFit", "Mapping", &P.mapping[0], P.nParaOpt);
  } 


   //vkv
   //PolPlot group
  if (P.runtype == 60) {
   Input.GetDoubleArray("PolPlot","Origin", P.origin, 3) ;   
   Input.GetDoubleArray("PolPlot","Direction", P.direction, 3) ;   
   P.polplot_pts = Input.GetInt("PolPlot","NumberOfGridPoints",10);
   P.polplot_stepsize = Input.GetDouble("PolPlot","StepSize",1.889);
  }
  //vkv

  // nice output
  if(rank == 0)cout << "-----------------------------------------------------------------\n";
  if(rank == 0) cout << "Input Parameters used:\n\n";
  P.Print();
  if(rank == 0) cout << "-----------------------------------------------------------------\n";
}






////////////////////////////////////////////////////////////////////////////////////
//
//   takes 3*nw lines and reads the coordinates of nw water monomers
//
//   the order (still) is O, H, H 
//   however, both 8/1 and O/H are accepted
//
//   the coordinates are put into a vector with 9*nw doubles
//
//
void GetWaterCoordinates(int nw, char **coorlines, double *WaterPos)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  int i, add, line = 0;
  for (i = 0; i < nw; ++i) {
    char atomsymbol[2];
    double ox, oy, oz, h1x, h1y, h1z, h2x, h2y, h2z;

    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &ox, &oy, &oz) != 4) {
     if(rank==0) cout << "Error reading input group Waters line" << line << ":"
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '8' && atomsymbol[0] != 'O')
      {printf("Line %i: O or 8 is expected!\n", line);exit(1);}
    line ++;

    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &h1x, &h1y, &h1z) != 4) {
      if(rank==0)cout << "Error reading input group Waters line" << line << ":" 
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '1' && atomsymbol[0] != 'H')
      {printf("Line %i: H or 1 is expected!\n", line);exit(1);}
    line ++;
      
    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &h2x, &h2y, &h2z) != 4) {
      if(rank==0)cout << "Error reading input group Waters line" << line << ":" 
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '1' && atomsymbol[0] != 'H')
      {printf("Line %i: H or 1 is expected!\n", line);exit(1);}
    line ++;
    
    add = i * 9;
    WaterPos[add+0] = ox;  WaterPos[add+1] = oy;  WaterPos[add+2] = oz;
    WaterPos[add+3] = h1x; WaterPos[add+4] = h1y; WaterPos[add+5] = h1z;
    WaterPos[add+6] = h2x; WaterPos[add+7] = h2y; WaterPos[add+8] = h2z;
  }
}



////////////////////////////////////////////////////////////////////////////////////
//
//  process the input group "Waters"
//  put coordinates into WaterPos, return the number of Water monomers found
//  read O, H, H (or 8, 1, 1)
//
int ReadWatersGroup(TSIN Input, dVec &WaterPos)
{

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  char **coorlines;
  int nlines = Input.GetGroup("Waters", &coorlines);
  if (nlines % 3 != 0) 
    {if(rank==0)cout << "error in ReadWaterCoordinates: No of input lines is not a multiple of 3\n"; exit(1);}
  int nwater = nlines / 3;
  
  WaterPos.resize(3*nlines, 0); // x, y, z for each atom

  int i, add, line = 0;
  for (i = 0; i < nwater; ++i) {
    char atomsymbol[2];
    double ox, oy, oz, h1x, h1y, h1z, h2x, h2y, h2z;

    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &ox, &oy, &oz) != 4) {
      if(rank==0)cout << "Error reading input group Waters line" << line << ":"
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '8' && atomsymbol[0] != 'O')
      {printf("Line %i: O or 8 is expected!\n", line);exit(1);}
    line ++;

    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &h1x, &h1y, &h1z) != 4) {
      if(rank==0)cout << "Error reading input group Waters line" << line << ":" 
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '1' && atomsymbol[0] != 'H')
      {printf("Line %i: H or 1 is expected!\n", line);exit(1);}
    line ++;
      
    if (sscanf(coorlines[line]," %s %lf %lf %lf", atomsymbol, &h2x, &h2y, &h2z) != 4) {
      if(rank==0)cout << "Error reading input group Waters line" << line << ":" 
	   << coorlines[line];
      exit(1);
    }
    if (atomsymbol[0] != '1' && atomsymbol[0] != 'H')
      {printf("Line %i: H or 1 is expected!\n", line);exit(1);}
    line ++;
    
    add = i * 9;
    WaterPos[add+0] = ox;  WaterPos[add+1] = oy;  WaterPos[add+2] = oz;
    WaterPos[add+3] = h1x; WaterPos[add+4] = h1y; WaterPos[add+5] = h1z;
    WaterPos[add+6] = h2x; WaterPos[add+7] = h2y; WaterPos[add+8] = h2z;
  }

  return nwater;

}





////////////////////////////////////////////////////////////////////
//  
//  Get a water cluster structure from an xyz file
//
//  read a structure as O, H, H    (no 8s or 1s are accepted)
//
//  natoms is the no of atoms to be read (without the B that may be present)
//
//
//  there may be a B atom (basis: coors of floating center in MCGeometries files)
//  setting nboron==1 ignores this lines
//
//  only for MCGeometry files:
//  the comment-line is analysed and if present the energy, EBE, dipole, and Seymour's Jz 
//  are returned in results
//
//  up to four values may be returned in results[]
//
//  this function is also used to read the PTMC input
//  here the comment-line is MCData  T maxstep
//  these 2 values are returned
//
//
//
int GetStructure(FILE *flin, int natoms, int nboron, double *coordinates, double *results)
{


  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  const int LineLength = 132;
  char buffer[LineLength+1];
  // overread leading empty lines
  do {
    if(fgets(buffer, LineLength, flin) == 0)
      return 1;
  } while (strlen(buffer) < 2);
  // no of atoms:
  if (atol(buffer) != natoms + nboron)
    {if(rank==0)cout << "Error in GetStructure:atoms = " << natoms << " != >>>" << buffer << "<<<\n"; exit(1);}
  // comment line
  // if data are present return them in results
  if(fgets(buffer, LineLength, flin) == 0)
    {if(rank==0)cout << "Error in GetStructure reading the comment line\n"; exit(1);}
  {
    double ene, ebe, dpl, a, b, c; 
    if (sscanf(buffer,"Move %*i %lf %lf %lf", &ene, &ebe, &dpl) == 3) {
      results[0] = ene;
      results[1] = ebe;
      results[2] = dpl;
       if (sscanf(buffer,"Move %*i %*f %*f %*f %lf %lf %lf", &c, &b, &a) == 3)
	 results[3] = sqrt(0.5*(a+b-c));
    }
    else {
      if (sscanf(buffer,"MCData %lf %lf", &a ,&b) == 2) {
	results[0] = a; 
	results[1] = b; 
      }
    }
  }

  // debug: cout << "# " << buffer;

  // first line may contain the position of the basis set as B
  if (nboron == 1) {
    if(fgets(buffer, LineLength, flin)==0)
      {if(rank==0)cout << "Error in GetStructure reading the boron atoms\n"; exit(1);}
    if (buffer[0] != 'B')
      {if(rank==0)cout << "Error in GetStructure: no B found in " << buffer << "\n"; exit(1);}
  }

  // read water monomers
  double x, y, z;
  int nw = natoms / 3;
  for (int i = 0; i < nw; ++i) {
    int add = 9*i;
    if(fgets(buffer, LineLength, flin)==0)
      {if(rank==0)cout << "Error in GetStructure reading O line of water " << i << "\n"; exit(1);} 
    if (sscanf(buffer, "O %lf %lf %lf", &x, &y, &z) != 3)
      {if(rank==0)cout << "Error in GetStructure reading O of water " << i << "\n"; exit(1);} 
    coordinates[add + 0] = x;
    coordinates[add + 1] = y;
    coordinates[add + 2] = z;
    if(fgets(buffer, LineLength, flin)==0)
      {if(rank==0)cout << "Error in GetStructure reading H1 line of water " << i << "\n"; exit(1);} 
    if (sscanf(buffer, "H %lf %lf %lf", &x, &y, &z) != 3)
      {if(rank==0)cout << "Error in GetStructure reading H1 of water " << i << "\n"; exit(1);} 
    coordinates[add + 3] = x;
    coordinates[add + 4] = y;
    coordinates[add + 5] = z;
    if(fgets(buffer, LineLength, flin)==0)
      {if(rank==0)cout << "Error in GetStructure reading H2 line of water " << i << "\n"; exit(1);} 
    if (sscanf(buffer, "H %lf %lf %lf", &x, &y, &z) != 3)
      {if(rank==0)cout << "Error in GetStructure reading H2 of water " << i << "\n"; exit(1);} 
    coordinates[add + 6] = x;
    coordinates[add + 7] = y;
    coordinates[add + 8] = z;
  }
  return(0);
} 



////////////////////////////////////////////////////////////////////
//  
//  Write a water cluster structure to an .xyz file
//
//  atoms as O, H, H
//
//  the comment line is:  #tag index comments[0]  ... comments[ncomments]   
//
//  or, if index < 0,   #tag comments[0]  ... comments[ncomments] 
//
//
//
void PutStructure(FILE *flout, int natoms, double *coordinates, int index, char *tag, int ncomments, double *comments)
{
  int i;
  // 1st line: number of atoms:
  fprintf(flout, "%i\n", natoms);

  // comment line
  fprintf(flout, "%s", tag);
  if (index >= 0)
    fprintf(flout, " %-8i", index);
  for (i = 0; i < ncomments; ++i)
    fprintf(flout, " %f", comments[i]);
  fprintf(flout, "\n");

  // write water monomers
  int nw = natoms / 3;
  double x, y, z;
  for (i = 0; i < nw; ++i) {
    int add = 9*i;
    x = coordinates[add + 0];
    y = coordinates[add + 1];
    z = coordinates[add + 2];
    fprintf(flout, "O  %10.6f  %10.6f  %10.6f\n", x, y, z);
    x = coordinates[add + 3];
    y = coordinates[add + 4];
    z = coordinates[add + 5];
    fprintf(flout, "H  %10.6f  %10.6f  %10.6f\n", x, y, z);
    x = coordinates[add + 6];
    y = coordinates[add + 7];
    z = coordinates[add + 8];
    fprintf(flout, "H  %10.6f  %10.6f  %10.6f\n", x, y, z);
  }
  fflush(flout);
} 

