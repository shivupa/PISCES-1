//TS
//    a Water object does not have any data
//    all (position, orientation, point charges, LJ terms, as well as the RepulsiveCore parameters)
//    are stored in its RepulsiveCore object
//

class Water
{
 public:

  int nSites;

  int nAtomCenters;
  AtomCenter *AtomCenters;  // the atom centers have charges, LJ parameters and basis functions
  double x_origin, y_origin, z_origin;
  double phi, theta, psi;


  Water(double x, double y, double z, int chargeselector, int RepcoreSelector, double divider);
  virtual ~Water();
  
  void Setup(int nAtoms, double x, double y, double z);
  void AddBasisfunctionCenter(int AtomIndex, double x, double y, double z, int nGTOs,
			      double charge, double sigma, double epsilon, double alpha);


  void AddRepulsiveCore(int RepcoreSelector);
  void SetupCharge(int chargeselector);
  void AddRepulsive_6S();
  void AddRepulsiveCCPVTZ();
  void AddRepulsiveCCPVTZ_all();
  void PrintPositions_InA(std::ostream& astream=std::cout);
  void ResetPositions();
  void UpdatePositions();
  void Reposition(double o_x, double o_y, double o_z, 
		  double h1_x, double h1_y, double h1_z, 
		  double h2_x, double h2_y, double h2_z);
  void Reposition_A(double o_x, double o_y, double o_z, 
		    double h1_x, double h1_y, double h1_z, 
		    double h2_x, double h2_y, double h2_z);
  
  int Dispersion; //used as a function selector
  double* DispersionOO;
  //double* DispersionOH;
  //double* DispersionHH;
  
  int Repulsion; //used as a function selector
  double* RepulsionOO;
  double* RepulsionOH;
  double* RepulsionHH;
  
  double* Damping;
  
  double PolSites;
  
 private:

  void ApplyEuler(double angle1, double angle2, double angle3);

  double RefPosOx, RefPosOy, RefPosOz;
  double RefPosHx, RefPosHy, RefPosHz;
  double RefPosMx, RefPosMy, RefPosMz;
  

  
  

};


