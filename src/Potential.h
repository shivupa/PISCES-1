//
//  this a a collection of different potentials for an excess electron
//
//  for each potential there is a Setup function that sets all the parameters
//  and site lists, say, as obtained from the neutral cluster calculation
//
//
//  there is only one function that evaluates the potential, which was setup, at r = (x, y, z)
//
#include "vecdefs.h"

#include "DistributedPolarizabilities.h"

class Potential
{

public:

  Potential()
    : verbose(0)
    , nSites(0)
    , nCharges(0)
    , nDipoles(0)
    , nPPS(0)
    , nGauss(0)
    , nContributions(0)
    , ChargeDamping(0)
    , DipoleDamping(0)
    , Alpha(0)
    , PolDamping(0)
    , VRepScale(0)
    , CationDamping(0)
    , AnionDamping(0)
    , nMolPol(0) {}


   void SetVerbose(int v);
   void Setup(const iVec& potflag, int nr, const double *r,
	      int nq, const double *q, const int *iq,
	      int nd, const double *d, const int *id,
	      int np, const double *a, const int *ip, const double *epc,
	      const double *DmuByDR, const double *potpara);

   double Evaluate(const double *r);
   void ReportEnergies(int n, double *e);
   void CalcDDTensorDer(double *TDerivXX ,
                     double *TDerivYY ,
                     double *TDerivZZ ,
                     double *TDerivXY ,
                     double *TDerivYZ ,
                     double *TDerivXZ ,
                     double *Rij, double alpha_i, double alpha_j, double aThole) ;

   void EvaluateGradient(const double *r, dVec& Grad, dVec& mu,dVec& mu_cross_mu, double wavefn, class WaterCluster &WaterN);
   void DerivElecField ( int nSites, double DampParameter, double *Rminus3 , double *R , double *Grad , double *mu) ;
   void FinalGradient( int nAtoms, double *Gradient, double *mu, double *mu_cross_mu,
                               double *dT_x, double *dT_y, double *dT_z, double *dEfield);

   void MuCrossMu (int nAtoms, double *mu ,dVec& mu1,dVec& mu_cross_mu  , double wavefn) ;
   void SubtractWWGradient (int nSites, double *PolGrad, double *Gradient) ;


   void SetupDPPTB(int nr, const double *r,
		   int nq, const double *q, const int *iq,
		   int nd, const double *d, const int *id,
		   const double *DmuByDR, const double *potpara);

  void SetupBloomfield(int nr, const double *r, int nq, const double *q, const int *iq, const double *PotPara);
  void SetupFullerElec(int nr, const double *r, const double *PotPara, const int *nAtomArray, const double *DipoleArray, const double *CenterArray);

  void UpdateChargeDipole(int nr, double *r, int nq, double *q, int *iq, int nd, double *d, int *id);

  void UpdateParameters(const double *PotPara);

  void PrintMinMax();
  
  // for SphericalPot
  void SetupSphericalPotential(const double *PotPara);


private: // METHODS

   void SetSites(int n, const double* sites);
   void SetCharges(int n, const double *q, const int *iq);
   void SetDipoles(int n, const double *d, const int *id);
   void SetGauss(int n, const double *expcoeff, const int *ig);
   void SetPPS(int n, const double *a, const int *ia);


   void SetupDPPGTOP(int nr, const double* r,
		     int nq, const double *q, const int *iq,
		     int nd, const double *d, const int *id,
		     int np, const double *a, const int *ip, const double *epc,
		     const double *DmuByDR, const double *potpara);

   double EvaluateDPPGTOP(const double *x);

   void Set6GTORepCore(int n);
   void Set12GTORepCore(int n);
   void Set4STORepCore(int n);

   double EvaluateChargePotential(const double *x, int DampFlag, double DampParameter);
   double EvaluateDipolePotential(const double *x, int DampFlag, double DampParameter);
   double EvaluateRepulsivePotential(const double *x, int RepCoreFlag);
   double EvaluatePolPot(const double *x, int DampFlag, double DampParameter);
   double EvaluateMolecularPolarizableSites(const double *x, int DampFlag, double DampParameter);
   double SelfConsistentPolarizability(const double *x, int DampFlag, double DampParameter);
   double SelfConsistentPolarizability_33(const double *x, int DampFlag, double DampParameter);

   void EvaluateDPP6SPGradient(const double *x, dVec& Grad);

//   void EvaluateDPPGTOPGradient(const double *x, dVec& Grad, double *dEfield, double *PolGrad,  class WaterCluster &WaterN);
   void EvaluateDPPGTOPGradient(const double *x, dVec& Grad,dVec& mu,dVec& mu_cross_mu, double wavefn, class WaterCluster &WaterN);
   double EvaluateDPPTB(const double *x);
   double EvaluateBloomfield(const double *x);

   void SetupMinMax();
   void CheckMinMax(const double *v);

   //for c60
   void SetPtDipoles( const int *NoAtomsArray, const double *DipoleArray , const double *CenterArray);
   void SetMoleculeFields(const int *NoAtomsArray);
   double EvaluateCsixty(const double *x);
   void SetSitesCsixty(int n, const double* sites);
   double EvaluateCsixtyPolPot(const double *x);
   double EvaluateCsixtyElecPot(const double *x);
   double EvaluateCsixtyRepPot(const double *x, double RepScale, double GaussExp);
   double InterFullerPolPot(const double *EfieldVec);
   // for sphericalPot
   double EvaluateSphericalPotential(const double *x);

private: // DATA MEMBERS

   // Flags  see Parameters.h/Parameters.cpp
   iVec PotFlags;
   int verbose;
   int DampType;           // 1=Thole-like GTO  2=effective distance
   int SigmaOFlag;         // separate scaling for Vrep on O nad H (kind of deprecated)
   int PolType;            // 1=isotropic molecular, 2=interacting atomic, non-interacting molecular, 3=fully interacting
   int AdiabaticPolPot;    // use adiabatic square-root formula (only for PolType 1 and 2)
   int RepCoreType;        // 1:GTOs  2:STOs

   // parameters from the input
   double ChargeDamping;
   double DipoleDamping;
   double PolDamping;
   double VRepScale;
   double VRepScaleO;
   double Alpha;
   double EpsDrude;
   // for NaCl
   double CationDamping;
   double AnionDamping;

   // for C60
   double CEx ; 
   double EsDampCsixty ; 
   double RepScaleCsixty ; 
   double GaussExpCsixty ; 
   double PolDampCsixty ; 
   double InterFullerPol ; 
   int PolFlagCsixty ; 
   double ExtraNucQ ; 
   int RepFlagCsixty ; 
   dVec PtDip ;
   dVec EfieldOnAtoms ; 
   int nFullerenes ;
   iVec nAtomsArray ;
   //iVec nAtomsArray ; 

   // for SphericalPot
   double WellDepth ; double WellWidth ;
   double ForceConstant ;  
   // data
   int nSites;      // number of sites for distance tables
   int nCharges;    // number of point charges
   int nDipoles;    // number of point dipoles
   int nPPS;        // number of point polarizable sites (non-interacting)
   int nGauss;      // number of sites with repulsive potentials
   int nMolPol;     // number of molecules with interacting atomic point polarizable sites

   double VpolWaterWater; // energy offset for self-consistent polarization

   dVec Site;        // 3 * nSites
   dVec Charge;      // nCharges
   iVec ChargeSite;  // mapping Charge -> Site
   dVec Dipole;      // 3 * nDipoles
   iVec DipoleSite;  // mapping Dipole -> Site
   dVec PPS;         // 2 * nPPs  (alpha / damping parameter)
   iVec PPSSite;     // mapping PP -> Site
   iVec PPSType;     // not used
   dVec Gauss;       // 2*nGauss (exp + coef)
   iVec GaussSite;   // mapping again

   struct DistributedPolarizabilities *MolPol;

   dVec Rx;
   dVec Ry;
   dVec Rz;
   dVec R;
   dVec R2;
   dVec Rminus3;

   dVec RVec;
   dVec X;
   dVec Y;
   dVec Z;

   // for MinMax debug output and detailed reporting
   int nContributions;
   int nReturnEnergies;
   dVec ReturnEnergies;
   dVec MaxV;
   dVec MinV;

};
