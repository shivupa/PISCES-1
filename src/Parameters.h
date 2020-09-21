//TS
//    collecting many of the parameters for WaterCluster in one struct
//    rudimentary documentation is in Parameters::Print() 
//
struct Parameters
{

  // Job group
  int runtype;
  int nElectron;
  int CubeFile;
  int WfCuts;

  // WaterModel group
  int WMverbose;
  int nWater;
  int KTFlag;
  int CenterFlag;

  // ElectronPotential group
  int PotVerbose;
  iVec PotFlag;
  double PotPara[32];

  // GridDef group
  int DVRType;
  int Sampling;
  int gridverbose;
  int ngrid[3];
  double gpara[3];

  // Diag group
  int diagverbose;
  int DiagMethod;
  int nStates;
  int maxSub;
  int maxIter;
  int ptol;
  int istartvec;

  // Optimize group
  int optverbose;
  double gtol;

  // Molecular Dynamics group
  int nsteps;
  int RigidBody;
  double timestep;

  // potfit 
  //   weights for NO from cube file
  double LowCutOff;   // weights R from next O must be bigger than this
  double HighCutOff;  // and smaller than this
  double wexp;        // and will get weighted with R^wexp * sqrt(NO amplitude)
  //   parameters to be optimized
  int nParaOpt;
  iVec mapping;
  int minimizer;
  
//vkv

  // PolPlot parameters
  double origin[3] ; 	 
  double direction[3] ; 	 
  int polplot_pts ;
  double polplot_stepsize;

//vkv

  //
  // Defaults are in GetInputParameters
  //
  Parameters() 
  {
    runtype = -1;
    PotFlag.resize(32); // so you were saying?
  }

  ~Parameters()
  {
  }

  void Print(void);

};
