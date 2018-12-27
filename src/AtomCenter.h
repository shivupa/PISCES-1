//
//  this should be a "site" class
//  a site can have many things, and some of the original LJ parameters
//  are now being used differently in DPP
//  the GTOs are used for two things:
//  (1) GTOs of AtomCenters in Repulsive core are used for the repulsive potential
//  (2) GTOs of AtomCenters of AtomCenters in ExcessElectron are used as a GTO basis set  
//  So for DPP stand alone there should not be any GTO related objects
//
struct AtomCenter
{
  int identifier; // should equal to atomsite index;
  int nGTOs;
  double x, y, z;
  double dx, dy, dz;
  double charge;
  double sigma;      // for L-J potential
  double epsilon;    // for L-J potential
  double alpha;      // polarizability
  double mx, my, mz; // induced dipole moment

  // derivative of the dipoles of the neutral cluster
  double *dxmx;
  double *dxmy;
  double *dxmz;
  double *dymx;
  double *dymy;
  double *dymz;
  double *dzmx;
  double *dzmy;
  double *dzmz;


  // derivatives of electronic energy in GTO basis
  double *dxSMatrix;
  double *dySMatrix;
  double *dzSMatrix;
  double *dxFMatrix;
  double *dyFMatrix;
  double *dzFMatrix;
  GTO* GTOs; 
  GTO* dxGTOs; 
  GTO* dxGTOs_1; 
  GTO* dyGTOs; 
  GTO* dyGTOs_1; 
  GTO* dzGTOs; 
  GTO* dzGTOs_1; 
  GTO* dxRGTOs; 
  GTO* dxRGTOs_1; 
  GTO* dyRGTOs; 
  GTO* dyRGTOs_1; 
  GTO* dzRGTOs; 
  GTO* dzRGTOs_1; 


  AtomCenter() {
     identifier = 0;
     nGTOs = 0;
     x = y = z = 0;
     dx = dy = dz = 0;
     charge = 0;
     sigma = 0;  
     epsilon = 0;
     alpha = 0; 
     mx = my = mz = 0; 

    GTOs = NULL;
    dxSMatrix =NULL;
    dySMatrix =NULL;
    dzSMatrix =NULL;
    dxFMatrix =NULL;
    dyFMatrix =NULL;
    dzFMatrix =NULL;
    dxmx =NULL;
    dxmy =NULL;
    dxmz =NULL;
    dymx =NULL;
    dymy =NULL;
    dymz =NULL;
    dzmx =NULL;
    dzmy =NULL;
    dzmz =NULL;

    dxGTOs =NULL;
    dxGTOs_1 =NULL;
    dyGTOs =NULL;
    dyGTOs_1 =NULL;
    dzGTOs =NULL;
    dzGTOs_1 =NULL;
    dxRGTOs =NULL;
    dxRGTOs_1 =NULL;
    dyRGTOs =NULL;
    dyRGTOs_1 =NULL;
    dzRGTOs =NULL;
    dzRGTOs_1 =NULL;
  }

  ~AtomCenter()
  {
    if (GTOs) delete[] GTOs; GTOs = NULL;
    if (dxGTOs) delete[] dxGTOs; dxGTOs = NULL;
    if (dxGTOs_1) delete[] dxGTOs_1; dxGTOs_1 = NULL;
    if (dyGTOs) delete[] dyGTOs; dyGTOs = NULL;
    if (dyGTOs_1) delete[] dyGTOs_1; dyGTOs_1 = NULL;
    if (dzGTOs) delete[] dzGTOs; dzGTOs = NULL;
    if (dzGTOs_1) delete[] dzGTOs_1; dzGTOs_1 = NULL;
    if (dxRGTOs) delete[] dxRGTOs; dxRGTOs = NULL;
    if (dxRGTOs_1) delete[] dxRGTOs_1; dxRGTOs_1 = NULL;
    if (dyRGTOs) delete[] dyRGTOs; dyRGTOs = NULL;
    if (dyRGTOs_1) delete[] dyRGTOs_1; dyRGTOs_1 = NULL;
    if (dzRGTOs) delete[] dzRGTOs; dzRGTOs = NULL;
    if (dzRGTOs_1) delete[] dzRGTOs_1; dzRGTOs_1 = NULL;
    if (dxSMatrix) delete[] dxSMatrix; dxSMatrix = NULL;
    if (dySMatrix) delete[] dySMatrix; dySMatrix = NULL;
    if (dzSMatrix) delete[] dzSMatrix; dzSMatrix = NULL;
    if (dxFMatrix) delete[] dxFMatrix; dxFMatrix = NULL;
    if (dyFMatrix) delete[] dyFMatrix; dyFMatrix = NULL;
    if (dzFMatrix) delete[] dzFMatrix; dzFMatrix = NULL;
    if (dxmx) delete[] dxmx; dxmx = NULL;
    if (dxmy) delete[] dxmy; dxmy = NULL;
    if (dxmz) delete[] dxmz; dxmz = NULL;
    if (dymx) delete[] dymx; dymx = NULL;
    if (dymy) delete[] dymy; dymy = NULL;
    if (dymz) delete[] dymz; dymz = NULL;
    if (dzmx) delete[] dzmx; dzmx = NULL;
    if (dzmy) delete[] dzmy; dzmy = NULL;
    if (dzmz) delete[] dzmz; dzmz = NULL;
  }

private :
//  AtomCenter& operator=(AtomCenter const &other){exit(1);};
//  AtomCenter(AtomCenter const &other){exit(1);};
};
