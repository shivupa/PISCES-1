//TS
//  everything about the MO class
//
//  so far I do not know more than the names suggest by themselves
//
struct MO
{
  int state;
  int nGTOs;
  GTO** GTOs;
  double *coefficients;
  double energy;
  bool needdeletecoefficients;
  
  MO() {
    coefficients=NULL;
    GTOs=NULL;
    nGTOs=0;
    needdeletecoefficients=true;
  }

  ~MO() {
    if (needdeletecoefficients && coefficients != 0) 
      {delete[] coefficients; coefficients = 0;}
    // we sometimes do not need to delete it.
    // we have to be very careful.
    GTOs=NULL;	// we do not need to delete it.
                // this is not a memory leak, GTOs should be deleted in AtomCenters.
  }
  
private :
//  MO& operator=(MO const &other){exit(1);};
//  MO(MO const &other){exit(1);};
};
