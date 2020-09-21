void InvChargeDipolePol(int nSites, const double *R, double *InvA, int nFullerenes, const int *NoAtomsArray);
void BuildPolarizationMatrix(int nSites, const double *R,  double *A, int nFullerenes, const int *NoAtomsArray);
void InvertCDMatrix(int n, double *A);

