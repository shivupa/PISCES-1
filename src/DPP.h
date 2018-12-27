#ifndef PISCES_DPP_H_
#define PISCES_DPP_H_

//
//  this class comes from the Drude code, where it is used for the neutral
//  cluster and for the excess electron Hamiltonian
//  this led to a lot of coding problems, and the idea in Pisces is to keep the 
//  water model and the model Hamiltonian as separate as practical
//  so in Pisces WaterCluster is just used for the water model,
//  but there is still lots of fossils from the Drude code

class WaterCluster
{
public:

   WaterCluster(void)
   {
      nwaters = 0;
      waters = 0;      
      nSitesPerMonomer = 4;        // DPP parameter also set in Water::nSites
      nChargesPerMonomer = 3;      // DPP parameter
      DipoleSCFFlag = 0;
      KTFlag = 0;
      //nOscs = 0;
      //Oscs = 0;
      //nRepCores = 0;
      //RepCores = 0;
      //dampingfactor = 0;
      //oscIndex = 0;
      //repCoreIndex = 0;


      //LMPotential = 0;    // Quantum Drude Model is the default 
      //SBE_Flag = 0;
      //ConfinementFlag = 0;
      //needrecalculation=true;
      //nevalue=0;
      //electronbindingenergy=0;
      //oldLJenergy=0;
      //oldLJNRenergy=0;	// LJ + Nuclear Replusion.
      //MomentofInertia[0]=0;
      //MomentofInertia[1]=0;
      //MomentofInertia[2]=0;
      //autodetachment_counter = 0;
      //dissociation_counter = 0;
      //CenterOfMass[0] = CenterOfMass[1] = CenterOfMass[2] = 0;
      //      SavePattern=NULL;
   }

   int ReportNoOfWaters(void);   
   int ReportNoOfSitesPerMonomer(void);
   int ReportNoOfPolSites(void);


   double ReportPolarization(void);
   void ReportNoOfSites(int &nSites, int &nCharges, int &nDipoles);


   void SetStructure(int nW, const double *WCoor, int unit, int CenterFlag, int KTCharges, int verbose = 0);
   int GetStructure(int nW, double *WCoor, double scale, int verbose = 0);

   void SetConfiguration(int nW, const double *Conf, int scale, int verbose = 0);
   int GetConfiguration(int nW, double *WConf, double scale, int verbose = 0);

   double CalcEnergy(int verbose = 0);
   double CalcGradient(int n, double *Grad, double *PolGrad, int verbose = 0);

   void GetLists(int nSites, double *Sites, int nCharges, double *q, int *iqs, 
		 int nDipoles, double *d, int *ids, double *DmuByDR, int verbose);
   void ReportInducedDipoles(int nDipoles, double *d, int *ids, int verbose);
   void ReportPolSitesAndField(int nSites, double *alpha, int *ipp, double *Epc, int verbose);


   void ConvertF (double *Gradient, double *gradEular);
   void GetDInducedDipoles(double *DmuByDR, double *dEfield, double *dTensor );
   void Calc_dInducedDipoles(int verbose, double *DmuByDR, double *dEfield, double *dTensor );
   void Calc_dTenosrDotDipoles(int verbose, double *mu, double *dTensor );



private:

   int nwaters;
   int nSitesPerMonomer;
   int nChargesPerMonomer;
   int DipoleSCFFlag;
   int KTFlag;

   Water **waters;


   //double dampingfactor;

   //int oscIndex, repCoreIndex;

   //   double CenterOfMass[3];
   //   double MomentofInertia[3];

   //int *SavePattern;


   double NuclearRepulsionEnergy; // these three are public member variables of Cluster with a somewhat misleading names
   double LJEnergy;               // the 1st is electrostatics+polarization, the 2nd is the van-der-Waals interaction
   double ElectrostaticEnergy;    // the 3rd is only the electrostatic contribution
   double PolarizationEnergy;    //  water-water polarization energy




   //virtual double CalculateRMSD();
   double ComputeInducedDipolesAndIntermolecularEnergy(int verbose);
   double ComputeInducedDipolesAndIntermolecularEnergyAndGradient(double *Gradient, int verbose);
   void LocateMassCenter();
   //virtual void SetConfiguration(double *configuration);
   double CalculateEnergy(int basissetflag, double *configuration, double &ebe);
   double CalculateEnergy(double *configuration, int NoOfBasisFns[6], double Exponents[9], 
      double FloatingPos[3], double Energy[3], int verbose);
   double CalculateGradient(double *configuration, int NoOfBasisFns[6], double Exponents[9],
      double FloatingPos[3], double Energy[3], double *Gradient, int verbose);
   void Setup(int nOscs, int nRepCores);
   void Setup(void);
   //void AddMolecule(RepulsiveCore *RepCore, DrudeOscillator *Osc);
   void AddMolecule(Water *water);
   void ReadBasisFunctions(std::istream &astream);
   void SetBasisFunctions(int HBasisSelector);
   void SetBasisFunctions(int NoOfBasisFns[6], double Exponents[9], double FloatingPos[3], int NewEl, int verbose);
   void AddDerivativeOfGTOs();
   void FormDerivativeOfSMatrix();
   void FormDerivativeOfKMatrix();
   void FormDerivativeOfVenMatrix();
   void FormDerivativeOfRepulsiveCores();
   void FormDerivativeOfDipoleofOscillator();
   void FormDerivativeOfGTOPotential();
   void FormDerivativeOfDipole(int verbose);
   double StartOverlapIntegral_S(GTO& gto1, GTO& gto2);
   double KineticEnergyIntegral(GTO& gto1, GTO& gto2);
   double StartRepulsivePotential(GTO& gto1, GTO& gto2);
   double StartCoulombicPotential(GTO& gto1, GTO& gto2);
   double StartdCoulombicPotential(GTO& gto1, GTO& gto2, int i, int k, int type);
   double Start_dRepulsivePotential(GTO& gto1, GTO& gto3, int i, int Atomtype, int XYZtype);
   double CalculateDipoleofOscillator(GTO& gto1, GTO& gto2);
   void SaveOscillator_Vec(GTO& gto1, GTO& gto2, double *SaveVec );
   double Calculate_dDipoleofOscillator(GTO& gto1, GTO& gto2, int iosc, int type);
   double CalculateDipoleofOscillator(GTO& gto1, GTO& gto2, double *dipolex, double *dipoley, double *dipolez);
   double CalculateDipoleofOscillator(double *dipolex, double *dipoley, double *dipolez, double *SaveVec);
   double CalcIntermolecularPotential(int verbose = 0);
   double CalcIntermolecularGradient(double *Gradient, double *PolGrad, int verbose = 0);
   double CalculatGTOPotential(GTO& gto1, GTO& gto2 );
   double Calculate_dGTOPotential(GTO& gto1, GTO& gto2, int iosc, int type );
   void ReadPositions(std::istream &astream);
   void ReadWaterClusterProperty(std::istream &astream);

   //void Water_Read(std::istream &astream);
   //void Water_Read_MC(std::istream &astream, int basissetflag);
   //void ProcessStructures(std::istream &astream);
   //void Getdistribution(std::istream &astream);
   //virtual void RestoreMassCenter();
   //virtual void PrintGeometry_InA(ostream& astream=cout);
   //double CalculateMomentsofInertia();
   //void CalculateMomentsofInertia(double moi[3]);
   //bool ReadxyzFormat(std::istream &astream);

   //void Playground(TSIN Input);
   void CalcCurrentEnergy(double &HFEBE, double &CIEBE, double &TotalE, int verbose);
   void LocalPolarization(int grid_flag, double *energy, int verbose);

   void GetRigidConf(double *configuration);
   void SetRigidConf(double *configuration);
   void GetXYZ(double *coordinates);
   void GetXYZM(double *coordinates);
   void ComputeCOM(int nwater, double *pos, double com[3]);

   void SetAutodetachmentCounter(int what);
   int GetAutodetachmentCounter();
   void SetDissociationCounter(int what);
   int GetDissociationCounter();

   //void SetGeometry(TSIN Input);
   void SetGeometry(double *WaterCoordinates);
   //void SetGlobalParameters(TSIN Input);

   void SetGenBasisSet(int nShells, int *l, double *Alpha, double *Center);


   void CalcInducedDipoles(int verbose=0);
   void ReportInducedDipoles();
   void AddUpTTMDipoles(int verbose);
   void SetLocalPotFlag(int value, int verbose=0);
   void SetSoftBoxParameters(double *para);
   void SetConfinementPotential(int flag, double rcut, double eta);


   void PutIntoCenterOfMassSystem(int nw, dVec& pos);
   //private:

   //double ConfinementPotential(double *configuration);

   //bool needrecalculation;
   //int autodetachment_counter;
   //int dissociation_counter;
   //int LMPotential;
   //int SBE_Flag; // flag for soft-box potential in the electron Hamiltonian

   // Confinement of the water monomers
   // 1=confinement on shortest next neighbour distance
   // 2=confinement in a cubic soft-box centered at the origin (=center of mass) 
   //int ConfinementFlag;  
   //double ConfineRcut, ConfineEta;

};

#endif // PISCES_DPP_H_
