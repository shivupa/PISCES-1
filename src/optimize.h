void Optimizer(dVec WaterPos, const Parameters InP);
void MolecularDynamics(dVec WaterPos, const Parameters InP);
double GetClusterEnergy(double *smplx);
void GetGeometry(double *smplx, double *energy);
double GetClusterEnergyMD(const dVec WaterPos);
//double GetClusterEnergy(double *smplx, const Parameters InP);
void GetAnalGrad(const double *configuration, double *analgrad);
void GetRigidbodyForce(const dVec WaterPos, double *force);
//void GetNumGrad(const double *configuration, double *numgrad, const Parameters InP);
void GetNumGrad(const double *configuration, double *numgrad);
void DumpGeometry(int nWater, double energy, double *coor);
void single_point(const dVec WaterPos, const Parameters InP);
void SetUpOptimizer(const dVec WaterPos, const Parameters InP);
void  UpdateCoordinate( dVec WaterPos, double *vel, double *acc, double dt);

// =BEGIN xis19
void WriteSPEnergyGradient(dVec WaterPos, const Parameters InP);
// =END xis19
