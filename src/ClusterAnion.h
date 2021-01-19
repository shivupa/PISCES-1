#ifndef PISCES_CLUSTERANION_H_
#define PISCES_CLUSTERANION_H_

class ClusterAnion
{
 public:

  ClusterAnion() 
    {
    }

  void SetUpClusterAnion(const dVec WaterPos, const Parameters InP);
  double SinglePoint(const dVec WaterPos);
  double EnergyFromConfiguration(double *WaterConf);
  double PrintConf(int nWater, double *WaterConf) ; 
  void GetAnalGrad( const double *WaterConf, double *analgrad);
  void GetNumGrad( const double *WaterConf, double *analgrad);

  ~ClusterAnion() {};

 private:

  Parameters Para;
  WaterCluster WaterN;
  Molecule *Molecules;
  Potential Vel;
  dVec ebes ;
  int nWater ;
  int nconverged ;
  DVR Hel ;  

};


#endif // PISCES_CLUSTERANION_H_

