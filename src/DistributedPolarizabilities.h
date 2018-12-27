//
//
// the electron interacts with one or more groups of interacting polarizabilities
//
// each group needs to establish its matrix InvA = (alpha^-1 - T)^-1
//
// in addition to the electron's field, there may be an "external" field, say, 
// due to point charges: Epc
//
#include "vecdefs.h"


struct DistributedPolarizabilities
{
 public:
  
  DistributedPolarizabilities()
    : nAtoms(0) {}

  int nAtoms;
  iVec SiteList;
  dVec Alpha;
  dVec InvA;
  dVec Tensor;
  dVec Epc;
//  dVec dTij ; 
  dVec dT_x ; 
  dVec dT_y ; 
  dVec dT_z ; 
};
