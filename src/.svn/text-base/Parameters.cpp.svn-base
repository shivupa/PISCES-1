#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "vecdefs.h"
#include "Parameters.h"


using namespace std;

void Parameters::Print()
{

  // Job
  cout << "  runtype   = " << runtype;
  switch (runtype)
    {
    case 1   : cout << " (single point)\n"; break;
    case 2   : cout << " (optimization)\n"; break;
    case 3   : cout << " (MDsimulation)\n"; break;
    case 42  : cout << " (potfit for electron model potentials)\n"; break;
    case 60  : cout << " (plot polarization potential)\n";break;  //vkv
    case 101 : cout << " (NaCl cluster)\n"; break;
    case 19  : break;
    case 102 : cout << " (C60 cluster)\n"; break; //vkv
    case 103 : cout << " (Model spherical potential)\n"; break; //vkv
    default  : cout << " Error in GetInputParameters: Unknown runtype\n"; exit(1);
    }
  cout << "  nElectron = " << nElectron;
  if (nElectron > 0)
    cout << " (Cluster with excess electron)\n";
  else
    cout << " (Neutral Cluster)\n";
  if (runtype == 1) {
    if (CubeFile > 0)
    cout << "  Cube files will be created for all states\n";
    if (WfCuts > 0)
      cout << "  Cuts through the main axis will be created for all states\n";
  }

  // Water model
  if (KTFlag == 1)
    cout << "\n  Water model uses DPP-like parameters with KT charges\n";
  else
    cout << "\n  Standard DPP water model\n";
  if (CenterFlag)
    cout << "  Origin will be shifted to the center of mass.\n";


  // ElectronPotential
  if (nElectron > 0) {
    cout << "\n  ElectronPotential = " << PotFlag[0];
    switch (PotFlag[0])
      {
      case 1:  // DPP-6S-P
      case 2:  // DPP-12S-P
      case 3:  // DPP-12S-P but with 4 STOs instead of 12 GTOs
      case 4:  // DPP-P (no explicit repulsive potential)
	switch (PotFlag[0])
	  {
	  case 1: cout << " (DPP charges, 6S-GTO repulsive core)\n"; break;
	  case 2: cout << " (DPP charges, 12S-GTO repulsive core)\n"; break;
	  case 3: cout << " (DPP charges, 4S-STO repulsive core)\n"; break;
	  case 4: cout << " (DPP charges, no explicit repulsive potential)\n"; break;
	  default: cout << " - unknown Potential for the electron\n"; exit(1);
	  }
	cout << "   PotPara 0: " << PotPara[0] << " (isotropic polarizability alpha)\n";
	cout << "   PotPara 1: " << PotPara[1] << " (damping length for polarizable sites [Bohr])\n";
	cout << "   PotPara 2: " << PotPara[2] << " (sigma for Vrep scaling)\n";
	cout << "   PotPara 3: " << PotPara[3] << " (Coulomb damping length [Bohr])\n";
	cout << "   PotPara 4: " << PotPara[4] << " (Dipole damping length [Bohr])\n";
	cout << "   Polarization model:\n";
	switch (PotFlag[3])
	  {
	  case 1:
	    cout << "     Non-interacting polarizable sites; ";
	    if (PotFlag[4] == 1)
	      cout << "adiabatic polarization potential with EpsDrude = " << PotPara[6] << " Hartree\n";
	    else
	      cout << "standard Sum alpha/R^4 potential\n";
	    break;
	  case 2:
	    cout << "     Non-interacting Molecules with interacting polarizable sites.\n";
	    break;
	  default:
	    cout << "     No polarization included.\n";
	  }
	if (PotFlag[2] == 1)
	  cout << "   PotPara 5: " << PotPara[5] << " (separate Vrep scaling for O centered functions)\n";
	cout << "   DampType = " <<  PotFlag[1];
	if (PotFlag[1] == 1)
	  cout << "   (Gaussian Damping as in the Drude code)\n";
	else if (PotFlag[1] == 2)
	  cout << "   (effective-R damping)\n";	
	else 
	  {cout << " Unknown DampType\n"; exit(1);}	
	break;
      case 11:
	cout << " (DPP charges, single polarizable site, Turi-Borgis repulsive core)\n";
	cout << "   PotPara 0: " << PotPara[0] << " (alpha)\n";
	cout << "   PotPara 1: " << PotPara[1] << " (b for polarization damping)\n";
	cout << "   PotPara 2: " << PotPara[2] << " (A1 parameter for Coulomb damping for O)\n";
	cout << "   PotPara 3: " << PotPara[3] << " (A1 parameter for Coulomb damping for H (not used))\n";
	cout << "   PotPara 4: " << PotPara[4] << " (B1O parameter of Vrep)\n";
	cout << "   PotPara 5: " << PotPara[5] << " (B1H parameter of Vrep)\n";
	cout << "   PotPara 6: " << PotPara[6] << " (B2O parameter of Vrep)\n";
	cout << "   PotPara 7: " << PotPara[7] << " (B2H parameter of Vrep)\n";
	cout << "   PotPara 8: " << PotPara[8] << " (B3O parameter of Vrep)\n";
	cout << "   PotPara 9: " << PotPara[9] << " (B3H parameter of Vrep)\n";
	break;
      case 1021:
        cout << "fuller elec potential, use with caution\n" ;
        break ;
      case 1023:
        cout << " spherical\n" ;
        break ;
      default:
	cout << "GetInputParameters: unknown potential for the excess electron: " << PotFlag[0] << "\n";
	exit(1);
      }
    switch (PotFlag[5])
      {case 0:
         cout << "Warning :User should know that the default parameters or the parameters\n \
                  specified in the input are being used for electron-water potential\n";
         break;
       case 1:
         cout << "Pisces internal parameters of electron potential being used \n" ;
         break;
       default:
         cout << "GetInputParameters: unknown option for internal-parameters in electron-water potential \n" ; 
      }
  }

  // GridDef group
  switch (DVRType)
    {
    case 1:
      cout << "\n  Harmonic oscillator DVR\n    Freqs [au] = ";
      break;
    default:
      cout << "\n  Particle-in-a-box (sine) DVR\n    Length [Bohr] = ";
      break;
    }
  cout << gpara[0] << ", "<<  gpara[1] << ", " <<  gpara[2] << "\n";
  cout << "    no of points = " << ngrid[0] << " x " << ngrid[1] << " x " << ngrid[2] 
       << " = " << ngrid[0] * ngrid[1] * ngrid[2]<< "\n";
  switch (Sampling)
    {
    case 1: cout << "    Simple sampling at DVR points.\n"; break;
    case 2: cout << "    Sampling with double density (8x more calls).\n"; break;
    case 3: cout << "    Sampling with triple density (27x more calls).\n"; break;
    case 4: cout << "    Sampling with 0.2 Bohr step (6x more calls).\n"; break;
    default:
      cout << "    Sampling = " << Sampling << " used as q in Sq smoothing OP\n";
    }


  // Diag group
  cout << "  \n  Diagonalization Parameters\n";
  cout << "    DiagMethod = " <<  DiagMethod;
  switch (DiagMethod)
    {
    case 0: cout << " (Full Diagonalization; only for debugging)\n"; break;
    case 1: cout << " (Lanczos-Arnoldi; ARPACK library)\n"; break;
    case 2: cout << " (Davidson)\n"; break;
    case 3: cout << " (Davidson with 0th-order Jacobi correction)\n"; break;
    case 4: cout << " (Davidson without correction = effective Lanczos-Arnoldi)\n"; break;
    default: cout << " Error in GetInputParameters: unknown diag method\n"; exit(1);
    }
  cout << "    nStates = " << nStates << " (no of eigenpairs to be found)\n";
  cout << "    maxSub = " << maxSub << " (maximal subspace size)\n";
  cout << "    maxIter = " << maxIter << " (maximal no of macro iterations)\n";
  cout << "    ptol = " << ptol << " (iteration tolerance = 10^-pTol)\n";
  cout << "    istartvec = " << istartvec;
  switch (istartvec)
    {
    case 0: cout << " (read from StartVector.cub)\n"; break;
    case 1: cout << " (particle-in-a-box-like startvector)\n"; break;
    case 2: cout << " (use last converged vector if available)\n"; break;
    default:
      cout << "  GetInput: illegal value for startvector flag \n";
      exit(1);
    }

  // Optimize group
  if (runtype == 2) {
    cout << "  \n  Optimization Parameters\n";    
    cout << "  Convergence tolerance = " << gtol << "\n"; 
    cout << "  Optimizer Verbose = " << optverbose << "\n"; 
  }

  // PotFit group
  if (nParaOpt > 0) {
    cout << "  \n  Fitting the electron's potential to reproduce a EOM-NO\n";
    cout << "    The EOM-NO is read from StartVector.cub\n";
    cout << "    Points are weighted according to the distance to the nearest O atom Rnext\n";
    cout << "    LowCutOff  = " << LowCutOff << " (of Rnext)\n";
    cout << "    HighCutOff = " << HighCutOff << " (of Rnext)\n";
    cout << "    Exponent for Rnext = " << wexp <<  "\n";
    cout << "    No. of Parameters used for the fit" << nParaOpt << "\n";
    for (int ip = 0; ip < nParaOpt; ++ip)
      cout << "      no. " << mapping[ip] << " = " << PotPara[mapping[ip]] << "\n";
    cout << "    Minimizer is ";
    switch (minimizer)
      {
      case 0: cout << "Downhill Monte Carlo\n"; break;
      case 1: cout << "Simplex\n"; break;
      case 2: cout << "Powell\n"; break;
      default:cout << "Nothing is minimized, too bad that.\n"; 
      }


  }
  
}
