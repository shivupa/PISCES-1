void InvChargeDipolePol(int nSites, const double *R, double *InvA, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType);
void BuildPolarizationMatrix(int nSites, const double *R,  double *A, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType);
void BuildPolarizationMatrix2(int nSites, const double *R,  double *A, int nFullerenes, const int *NoAtomsArray);
void BuildPolarizationMatrix3(int nSites, const double *R,  double *A, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType);
void InvertCDMatrix(int n, double *A);
void InverseMatrix(double* Input, double* Output,  int n );
double det(int n, double mat[10][10]);


