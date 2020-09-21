void Optimizer(dVec WaterPos, const Parameters InP);
double GetClusterEnergy(double *smplx);
void GetAnalGrad(const double *configuration, double *analgrad);
void GetNumGrad(double *configuration, double *analgrad);
void DumpGeometry(int nWater, double energy, double *coor);
void polplot(const dVec WaterPos, const Parameters InP);
