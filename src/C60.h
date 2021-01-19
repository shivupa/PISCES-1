#ifndef PISCES_C60_H_
#define PISCES_C60_H_
#include "tsin.h"
void C60_sp(TSIN Input);
int  GetAtomCoordinates(int nLines, char **coorlines, double *IonPos, int *IonType);
void PrintAtomCoordinates(int nIons, double *IonPos, int *IonType, double scale);
#endif // PISCES_C60_H_
