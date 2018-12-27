// TS
//    
//  this is a data strcuture to be populated by the input parser, and will be given 
//  first to the force field, which may populate dipoles for the electron,
//  and then to the electron model
//
//  so far a molecule can have 
//      sites with point charges, 
//      sites with dipole moments, 
//      sites with point-polarizabilities 
//
//  other sites needed either by the force field or by the electron are not (yet) 
//  taken care of, but I hope can be dealt with by "Type" and "Tag"
//
//
//  Type:  1  rigid 4-site water Tags: O H H M
//
//
//
//

#include "vecdefs.h"

class Molecule
{
 public:

  Molecule()
    : Type(0)
    , nSites(0)
    , nCharges(0)
    , nDipoles(0)
    , nPntPols(0) {}


  int Type;
  int nSites;

  char *Tag;  // one character for each Site

  int nCharges;
  int nDipoles;
  int nPntPols;
  
  dVec Coor;         // Cartesian Corrdinates for each site
  dVec ChargeList;   // List of charges (nCharges)
  dVec DipoleList;   // list of dipoles (3*nDipoles)
  dVec PntPolList;   // list of polarizabilities (nPntPol)

  iVec ChargeSite;   // mapping Charge -> Site
  iVec DipoleSite;   // mapping Dipole -> Site
  iVec PntPolSite;   // mapping PntPol -> Site

  void SetupDPPWater(const double *WaterCoor);

 
 private:

};


