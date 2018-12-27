#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime> 

#include <iomanip>
#include <fstream>
#include <iostream>

using namespace std;
#include "ChargeDipPol.h"
#include "vecdefs.h"
#include "lapackblas.h"
#include "constants.h"
#include "timer.hpp"

//
//  this is for setting up a molecule with interacting atomic polarizability
//  InvA is the matrix needed to compute the induced dipoles from the "external" E field
//

//////////////////////////////////////////////////////////////
//
//   InvA = (alpha^-1 - T)^-1 
//
//   this is the matrix that can be multiplied with the electric field at
//   the polarizable sites to yield the induced dipoles:
//
//   mu = InvA * E
//
//   dim(A) = 3*nSiyes x 3*nSites
//
void InvChargeDipolePol(int nSites, const double *R, double *InvA, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType)
{

 int verbose = 0;
progress_timer t("ComputeInvPolarization", verbose);

  BuildPolarizationMatrix(nSites, R, InvA, nFullerenes, NoAtomsArray, CarbonType, IonType);  
  InvertCDMatrix(4*nSites+nFullerenes, InvA);

//  BuildPolarizationMatrix2(nSites, R, InvA, nFullerenes, NoAtomsArray );  
//  InvertCDMatrix(3*nSites, InvA);

//  BuildPolarizationMatrix3(nSites, R, InvA, nFullerenes, NoAtomsArray, CarbonType, IonType);  
//  InvertCDMatrix(nSites+nFullerenes, InvA);
 cout<<"start InvA ~~~~"<<endl;
//  InvertCDMatrix(nSites*(nSites-1)/2, InvA);

 int n = nSites;
int n4plus1 = 4*nSites+nFullerenes; 
double alpha_xx = 0.0;
double alpha_yy = 0.0;
double alpha_zz = 0.0;
double alpha_xy = 0.0;
double alpha_yx = 0.0;
double alpha_xz = 0.0;
double alpha_zx = 0.0;
double alpha_yz = 0.0;
double alpha_zy = 0.0;
double temp=0.0;
/*

double Q_xx = 0.0;
double Q_yy = 0.0;
double Q_zz = 0.0;
double Q_xy = 0.0;
double Q_xz = 0.0;
double Q_yz = 0.0;

  for (int i = 0; i < n; ++i){
    double Ri2 = R[3*i+0]*R[3*i+0] + R[3*i+1]*R[3*i+1] + R[3*i+2]*R[3*i+2];
    Q_xx += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+0] - Ri2;
    Q_yy += 3.0*InvA[n4plus1*i]*R[3*i+1]*R[3*i+1] - Ri2;
    Q_zz += 3.0*InvA[n4plus1*i]*R[3*i+2]*R[3*i+2] - Ri2;
    Q_xy += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+1];
    Q_xz += 3.0*InvA[n4plus1*i]*R[3*i+0]*R[3*i+2];
    Q_yz += 3.0*InvA[n4plus1*i]*R[3*i+1]*R[3*i+2];
  }
*/

  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){

//For the charge-charge term




       alpha_xx += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+0];
       alpha_xy += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+1];
       alpha_xz += InvA[n4plus1*i+j]*R[3*i+0]*R[3*j+2];
       alpha_yx += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+0];
       alpha_yy += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+1];
       alpha_yz += InvA[n4plus1*i+j]*R[3*i+1]*R[3*j+2];
       alpha_zx += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+0];
       alpha_zy += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+1];
       alpha_zz += InvA[n4plus1*i+j]*R[3*i+2]*R[3*j+2];



// For the dipole-chrage term


 //    cout<<"n4plus1 ="<<n4plus1<<endl;
 //   cout<<" i and j zz " <<i<<" "<<j<<" "<<n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i<<" "<<endl;
 //   cout<<" i and j xx " <<i<<" "<<j<<" "<<n*n4plus1+j*3*n4plus1+ i<<" "<<endl;
  



       alpha_xx -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+0]; 
       alpha_xy -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+1]; 
       alpha_xz -= 2*InvA[n*n4plus1+j*3*n4plus1                  +i]*R[3*i+2]; 
       alpha_yx -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+0]; 
       alpha_yy -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+1]; 
       alpha_yz -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1         +i]*R[3*i+2]; 
       alpha_zx -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+0]; 
       alpha_zy -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+1]; 
       alpha_zz -= 2*InvA[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i]*R[3*i+2]; 



 //       cout<<"InvA[n4plus1*i+j] = " << InvA[n4plus1*i+j]<<endl;

//      cout<<"R["<<3*i+0<<"]="<<R[3*i+0]<<endl;
//      cout<<"R["<<3*i+1<<"]="<<R[3*i+1]<<endl;
//      cout<<"R["<<3*i+2<<"]="<<R[3*i+2]<<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+0] =" <<InvA[i*n4plus1+n+3*j+0] <<endl;
//      cout<<"muti =" <<InvA[i*n4plus1+n+3*j+0] * R[3*i+0] <<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+1] =" <<InvA[i*n4plus1+n+3*j+1] <<endl;
//      cout<<"InvA[i*n4plus1+n+3*j+2] =" <<InvA[i*n4plus1+n+3*j+2] <<endl;

//        cout<<"alpha_xx = "<<alpha_xx<<endl;
//        cout<<"R[3*i+0] R[3*i+1]  = "<<R[3*i+0]<<" "<<R[3*i+1]<<endl;


// For the dipole-dipole term



       alpha_xx += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+0];
       alpha_xy += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+1]; 
       alpha_xz += InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+2]; 
       alpha_yx += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0]; 
       alpha_yy += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1];
       alpha_yz += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2]; 
       alpha_zx += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0]; 
       alpha_zy += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1]; 
       alpha_zz += InvA[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2];


/*
       alpha_xx += InvA[3*i*n*3           +j*3+0];
       alpha_xy += InvA[3*i*n*3           +j*3+1]; 
       alpha_xz += InvA[3*i*n*3           +j*3+2]; 
       alpha_yx += InvA[3*i*n*3 +n*3      +j*3+0]; 
       alpha_yy += InvA[3*i*n*3 +n*3      +j*3+1];
       alpha_yz += InvA[3*i*n*3 +n*3      +j*3+2]; 
       alpha_zx += InvA[3*i*n*3 +n*3+ n*3 +j*3+0]; 
       alpha_zy += InvA[3*i*n*3 +n*3+ n*3 +j*3+1]; 
       alpha_zz += InvA[3*i*n*3 +n*3+ n*3 +j*3+2];
*/

//       cout<<" alpha_xx ["<<3*i*n*3+j*3+0<<"]="<<InvA[3*i*n*3 +j*3+0]<<endl;
//       cout<<" alpha_yy ["<<3*i*n*3+ n*3+j*3+1<<"]="<<InvA[3*i*n*3+n*3 +j*3+1]<<endl;
//       cout<<" alpha_zz ["<<3*i*n*3+ n*3+n*3+j*3+2<<"]="<<InvA[3*i*n*3+n*3+n*3 +j*3+2]<<endl;

//       cout<<"2 alpha_xx ["<< n*n4plus1+3*i*n4plus1+n+j*3+0 <<"]="<<InvA[n*n4plus1+3*i*n4plus1 +                  n+j*3+0]<<endl;
 

       temp += InvA[i*n4plus1+n+3*j+0];
//      cout<<"index =  "<< i*n4plus1+n+3*j+0<<endl;
//      cout<<"temp = "<<temp<<" "<< InvA[i*n4plus1+n+3*j+0]<<endl;

/*
       alpha_xx += InvA[i*n4plus1+n+3*j+0];
       alpha_yy += InvA[i*n4plus1+n+3*j+1];
       alpha_zz += InvA[i*n4plus1+n+3*j+2];
*/


    }
  } 
 cout<<"---Polarizability  -------"<<endl;
 cout<<"alpha_xx = "<<alpha_xx<<endl;
 cout<<"alpha_yy = "<<alpha_yy<<endl;
 cout<<"alpha_zz = "<<alpha_zz<<endl;
 cout<<"alpha_xy = "<<alpha_xy<<endl;
 cout<<"alpha_yx = "<<alpha_yx<<endl;
 cout<<"alpha_xz = "<<alpha_xz<<endl;
 cout<<"alpha_zx = "<<alpha_zx<<endl;
 cout<<"alpha_yz = "<<alpha_yz<<endl;
 cout<<"alpha_zy = "<<alpha_zy<<endl;

/*
 cout<<"---Quadrupole moment  -------"<<endl;
 cout<<"Q_xx = "<<Q_xx<<endl;
 cout<<"Q_yy = "<<Q_yy<<endl;
 cout<<"Q_zz = "<<Q_zz<<endl;
 cout<<"Q_xy = "<<Q_xy<<endl;
 cout<<"Q_xz = "<<Q_xz<<endl;
 cout<<"Q_yz = "<<Q_yz<<endl;
*/

/*

  for (int i = 0; i < n4plus1; ++i){
    for (int j = 0; j < n4plus1; ++j){
     cout<<" "<<InvA[i*n4plus1+j];
    }
    cout<<endl;
  }
*/


}

//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix2(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray)
{

  int n = nSites;
  int n3 = n*3 ;

  for (int k = 0; k < n3*n3 ; ++k)
    A[k] = 0.0;
  
//  double R_const = 0.62*1.889725989 ; 
  double R_const = 0.68620399*1.889725989 ; 
  double chi = sqrt(2.0/PI)/R_const  ;  
  double Rqq = R_const*sqrt(2.0)  ;
  //cout << "R_const , R_qq" << R_const << " "<<Rqq << endl ; 
  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  
  double alpha = 3.0*sqrt(PI/2)*R_const*R_const*R_const ;
  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    for (int i = 0; i < n; ++i){
      if(i==j){
        A[3*i*n3 +        j*3+0] = 1.0/alpha ;
        A[3*i*n3 +n3+     j*3+1] = 1.0/alpha ;
        A[3*i*n3 +n3+ n3+ j*3+2] = 1.0/alpha ;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq ; 
        double Txy, Txz, Tyz ;
        fac = erf(Rij/Rqq)-2.0*Rij*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq) ;
        fac_pq = fac/(Rij*Rij*Rij) ; 
        fac1 = fac_pq/(Rij*Rij) ; 
        fac2 = 4.0*exp(-1.0*(Rij/Rqq)*(Rij/Rqq))/(sqrt(PI)*Rqq*Rqq*Rqq*Rij*Rij) ;
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;
        A[3*i*n3 +        j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[3*i*n3 +n3+     j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[3*i*n3 +n3+ n3+ j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 
        A[3*i*n3 +        j*3+1] = -1.0*Txy ; 
        A[3*i*n3 +n3+     j*3+0] = -1.0*Txy ; 
        A[3*i*n3 +        j*3+2] = -1.0*Txz ; 
        A[3*i*n3 +n3+ n3+ j*3+0] = -1.0*Txz ; 
        A[3*i*n3 +n3+     j*3+2] = -1.0*Tyz ; 
        A[3*i*n3 +n3+ n3+ j*3+1] = -1.0*Tyz ;
         
        //Tpq

      }
     
    }
  }


}


//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType)
{

  int n = nSites;
  int n4 = n*4 ;
  int n4plus1 = 4*n + nFullerenes ; // dimension of super-matrix A and leading dimension of any block B
  int nplus1 = n + nFullerenes ; 

  cout<<"nSites ="<<nSites<<endl;
  int verbose =0;
  progress_timer t("ComputeBuildPolarization", verbose);

//   cout<<"n4plus1 = "<<n4plus1<<endl;
//   cout<<"CarbonType = "<<CarbonType<<endl;
// for (int k= 0; k<n ; ++k )
//   cout<<"IonType["<<k<<"] = "<<IonType[k]<<endl;

  for (int k = 0; k < n4plus1*n4plus1 ; ++k)
    A[k] = 0.0;
  
  double* Matrix= new double[nplus1*nplus1];
  double ChargePen=-1.5;
// double R_const = 0.62*1.889725989 ; 
//  double R_const = 0.68620399*1.889725989 ; 
  double Rqi, Rqj, chi; 
 // double Rqq = R_const*sqrt(2.0)  ;
  //cout << "R_const , R_qq" << R_const << " "<<Rqq << endl ; 
  //Tqq
  
  double R_const;
   if (CarbonType == 2) R_const=0.01945;
   else if (CarbonType == 3) R_const=0.02279;
   else if (CarbonType == 4) R_const=0.02652;
   else if (CarbonType == 5) R_const=0.03029;
   else if (CarbonType == 6) R_const=0.215;
   else if (CarbonType == 7) R_const=0.01838;
   else if (CarbonType == 8) R_const=0.01758;

   double H_R=0.0;

  for (int j = 0; j < n; ++j) {
    //if (IonType[j] == 1) Rqj= 0.67200149*1.889725989 ;
   // if (IonType[j] == 1) Rqj= 0.66329437*1.889725989 ;
    if (IonType[j] == 1) Rqj= 0.68620399*1.889725989 ;
    else if (IonType[j] == 2) Rqj= H_R*1.889725989 ;
    else if (IonType[j] == 3) Rqj= R_const*1.889725989 ;
   // cout << "Rqj " << Rqj << endl ; 

    for (int i = 0; i < n; ++i) {
     // if (IonType[i] == 1) Rqi= 0.67200149*1.889725989 ;
     // if (IonType[i] == 1) Rqi= 0.66329437*1.889725989 ;
      if (IonType[i] == 1) Rqi= 0.68620399*1.889725989 ;
      else if (IonType[i] == 2) Rqi= H_R*1.889725989 ;
      else if (IonType[i] == 3) Rqi= R_const*1.889725989 ;
      chi = sqrt(2.0/PI)/Rqi  ;
     // cout << "Rqi " << Rqi << endl ; 
      if (i == j) {
        A[n4plus1*i+j] = chi; 
        Matrix[nplus1*i+j] = chi;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double CT=4.0*pow(10,4)*exp(-1.5*Rij);
        double Rqq = sqrt(Rqi*Rqi+Rqj*Rqj)  ;

/*
        double C = 1.00;
        double B = 0.43; 
        double w = Rij*B;
        if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
         A[n4plus1*i+j] = erf(Rij/Rqq)/Rij *  C*exp(-w)*(1+w+1.0/3.0*w*w) ; 
        }
        else  A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ; 
*/
         A[n4plus1*i+j] = erf(Rij/Rqq)/Rij ; 
         Matrix[nplus1*i+j] = erf(Rij/Rqq)/Rij ;
         // A[n4plus1*i+j] = erf(Rij/Rqq)/(Rij*CT) ; 
         // cout<<" erf(Rij/Rqq) =  "<< erf(Rij/Rqq) <<endl; 
         //cout<<" i j : "<<i<<" "<<j<<" "<< erf(Rij/Rqq)/Rij<<endl; 
      }
    }
  }

  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  //double Aconst=-341.; // For 8.0 Angstrom of dimer
  // double Aconst=-11.85; // For 9.0 Angstrom of dimer
 //  double Aconst=-1.0; // 
 //  double Aconst=-4.02; // For 9.25 Angstrom of dimer
 //  double Aconst=-2.08; // For 9.50 Angstrom of dimer
 //  double Aconst=-1.9811; // For fit 9.50 Angstrom of dimer
  // double Aconst=-1.03; // For 9.75 Angstrom of dimer
 //  double Aconst=-1.1209; // For fit 9.75 Angstrom of dimer
 //  double Aconst=-0.525; // For 10.05 Angstgrom of dimer
 //  double Aconst=-0.566; // For fit 10.05 Angstgrom of dimer
  // double Aconst=-0.345; // For 10.25 Angstgrom of dimer
  // double Aconst=-0.3589; // For fit 10.25 Angstgrom of dimer
  //double Aconst=-0.21; // For 10.50 Angstgrom of dimer
 //double Aconst=-0.2031; // For fit 10.50 Angstgrom of dimer
 // double Aconst=-0.128; // For 10.75 Angstgrom of dimer
 // double Aconst=-0.079; // For 11.0 Angstrom of dimer 
 // double Aconst=-0.065; // For fit 11.0 Angstrom of dimer 
 // double Aconst=-0.0091; // For 12.0 Angstrom of dimer 
 double Aconst=-0.0067; // For fit 12.0 Angstrom of dimer 
  //double Aconst=-0.0003; // For 13.0 Angstrom of dimer 
//  double Aconst=-0.00000001; // For 14.0 Angstrom of dimer 
  // double Aconst= 0.0000000;



 //double Bconst= -11.85; 
 //  double Bconst=-11.85; // For 9.0 Angstrom of dimer
 //  double Bconst=-4.02; // For 9.25 Angstrom of dimer
 //  double Bconst=-2.08; // For 9.50 Angstrom of dimer
 //  double Bconst=-1.03; // For 9.75 Angstrom of dimer
   double Bconst=-0.525; // For 10.05 Angstgrom of dimer
 //  double Bconst=-0.345; // For 10.25 Angstgrom of dimer
 // double Bconst=-0.079; // For 11.0 Angstrom of dimer 
 // double Bconst=-0.0091; // For 12.0 Angstrom of dimer 
 // double Bconst= 0.000000; 
  for ( int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom      = current_atom ;
    current_atom += NoAtomsArray[iFuller] ;
    for (int i = old_atom; i < current_atom ; ++i) {
      A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows


/*
     if (iFuller ==0) {
      A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[i*n4plus1+n4 + iFuller+1]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      A[n4plus1*n4 + n4plus1*iFuller + i + n4plus1]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[i*nplus1+n + iFuller+1]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[nplus1*n + nplus1*iFuller + i + nplus1]     = 1.0 ; // corresponds to rows of lagrange rows
    }
      if (iFuller ==1) { 
      A[i*n4plus1+n4 + iFuller]               = -1.0 ; // corresponds to columns of lagrange vector
      A[i*n4plus1+n4 + iFuller-1]               = 1.0 ; // corresponds to columns of lagrange vector
      A[n4plus1*n4 + n4plus1*iFuller + i]     = -1.0 ; // corresponds to rows of lagrange rows
      A[n4plus1*n4 + n4plus1*iFuller + i - n4plus1]     = 1.0 ; // corresponds to rows of lagrange rows
      Matrix[i*nplus1+n + iFuller]               = -1.0 ; // corresponds to columns of lagrange vector
      Matrix[i*nplus1+n + iFuller-1]               = 1.0 ; // corresponds to columns of lagrange vector
      Matrix[nplus1*n + nplus1*iFuller + i]     = -1.0 ; // corresponds to rows of lagrange rows
      Matrix[nplus1*n + nplus1*iFuller + i - nplus1]     = 1.0 ; // corresponds to rows of lagrange rows
     }
*/
     // cout<< "A[ "<<i*n4plus1+n4 + iFuller<<"]= "<<A[i*n4plus1+n4 + iFuller]<<endl; 
    //  cout<< "A[ "<<i*n4plus1+n4 + iFuller-1<<"]= "<<A[i*n4plus1+n4 + iFuller-1]<<endl; 
     // cout<< "A[ "<<n4plus1*n4 + n4plus1*iFuller + i<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + i]<<endl; 
     // A[i*n4plus1+n4 + iFuller]               = 0.0 ; // corresponds to columns of lagrange vector
     //  A[n4plus1*n4 + n4plus1*iFuller + i]     = 0.0 ; // corresponds to rows of lagrange rows
    }
    if (iFuller == 0 ) {
      cout<<"nFullerenes="<<nFullerenes<<endl;
      if (nFullerenes == 2 ) {
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Aconst; // corresponds to rows of lagrange rows
      }
      else if (nFullerenes == 3 ) {
       // A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     = Aconst; // corresponds to rows of lagrange rows
       // A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
      //  A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = 0.0 ; // corresponds to rows of lagrange rows
// triangle
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     = Aconst+Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Aconst; // corresponds to rows of lagrange rows
      }
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]<<endl; 
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
    cout<< "0 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }
    else if (iFuller == 1 ) {
       if (nFullerenes == 2 ) {
         A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
         A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = Aconst ; // corresponds to rows of lagrange rows
       }
      else if (nFullerenes == 3 ) {
     //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  -Aconst; // corresponds to rows of lagrange rows
     //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = Aconst+Bconst; // corresponds to rows of lagrange rows
     //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Bconst ; // corresponds to rows of lagrange rows
// triangle
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  -Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = Aconst+Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = -Aconst ; // corresponds to rows of lagrange rows

      }
       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]<<endl; 
       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
       cout<< " 1 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }
    else if (iFuller == 2 ) {
    //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  0.0; // corresponds to rows of lagrange rows
    //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Bconst; // corresponds to rows of lagrange rows
    //   A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = Bconst ; // corresponds to rows of lagrange rows

// triangle
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]     =  -Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]     = -Aconst; // corresponds to rows of lagrange rows
       A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]     = Aconst+Aconst ; // corresponds to rows of lagrange rows


       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-3<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-3]<<endl; 
       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-2<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-2]<<endl; 
       cout<< " 2 A[ "<<n4plus1*n4 + n4plus1*iFuller + n4plus1-1<<"]= "<<A[n4plus1*n4 + n4plus1*iFuller + n4plus1-1]<<endl; 
    }

  }
  

  double alpha, alphaX, alphaY, alphaZ;
  double Rpi, Rpj;
  double alpha_parC, alpha_perpC, alpha_parH, alpha_perpH;

   if (CarbonType == 1) { 
     // R_const=0.67200149;
    //  alpha_parC = 1.2517105;
    //  alpha_perpC =1.2517105; 
    
      alpha_parC = 1.214899;
      alpha_perpC =1.214899; 
      R_const= 0.68620399;
    //  R_const= 0.66329437;
    //  alpha_parC = 1.41 ; 
    //  alpha_perpC = 1.14 ;

    //  alpha_parC = 1.4499997 ; 
    //  alpha_perpC = 0.54496647 ;

    //  alpha_perpC = 1.4499997 ; 
    //  alpha_parC = 0.54496647 ;
   }
   else if (CarbonType == 2) { 
      R_const=0.01945;
      alpha_parC = 1.2885;
      alpha_perpC = 1.2885;
      alpha_parH = 0.3384 ;
      alpha_perpH = 0.3384 ;
   }
   else if (CarbonType == 3) {
      R_const=0.02279;
      alpha_parC = 1.4755;
      alpha_perpC = 2.2611;
      alpha_parH = 0.3057 ;
      alpha_perpH = 0.3057 ;
   }
   else if (CarbonType == 4) { 
      R_const=0.02652;
      alpha_parC = 1.2080;
      alpha_perpC = 1.2080;
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 5) {
       R_const=0.03029;
      alpha_parC = 1.3632;
      alpha_perpC = 2.1671; 
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 6) {
      R_const=0.215;
      alpha_parC = 1.317;
      alpha_perpC = 1.84 ; 
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 7) {
      R_const=0.01838;
      alpha_parC = 1.5324;
      alpha_perpC = 2.2360 ;
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 8) {
      R_const=0.01758;
      alpha_parC = 1.2597;
      alpha_perpC = 1.2597 ;
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }



  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    if (IonType[j] == 1) { 
      Rqj= R_const*1.889725989 ; 
      Rpj=Rqj;
      alphaX= alpha_parC/(0.529177*0.529177*0.529177);
      alphaY= alpha_parC/(0.529177*0.529177*0.529177);
      alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);

    //  alphaX = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
    //  alphaY = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
    //  alphaZ = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
    }
    else if (IonType[j] == 2) {
     // Rqj= 0.00 ;
      Rqj= H_R*1.889725989 ;
      alphaX= alpha_parH/(0.529177*0.529177*0.529177);
      alphaY= alpha_parH/(0.529177*0.529177*0.529177);
      alphaZ= alpha_perpH/(0.529177*0.529177*0.529177);
      Rpj= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
     // Rpj= pow(alphaX/3.0/sqrt(PI/2),1.0/3.0);
    }
    else if (IonType[j] == 3){
      Rqj= R_const*1.889725989 ; 
      alphaX= alpha_parC/(0.529177*0.529177*0.529177);
      alphaY= alpha_parC/(0.529177*0.529177*0.529177);
      alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);
      Rpj= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
     // Rpj= pow(alphaX/3.0/sqrt(PI/2),1.0/3.0);
    }

    for (int i = 0; i < n; ++i){
      if (IonType[i] == 1) {
        Rqi= R_const*1.889725989 ; 
        Rpi=Rqi;
        alphaX= alpha_parC/(0.529177*0.529177*0.529177);
        alphaY= alpha_parC/(0.529177*0.529177*0.529177);
        alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);

      //  alphaX = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
      //  alphaY = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
      //  alphaZ = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
       }
      else if (IonType[i] == 2) {
       // Rqi= 0.00 ;
        Rqi= H_R*1.889725989 ;
        alphaX= alpha_parH/(0.529177*0.529177*0.529177);
        alphaY= alpha_parH/(0.529177*0.529177*0.529177);
        alphaZ= alpha_perpH/(0.529177*0.529177*0.529177);
        Rpi= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
       // Rpi= pow(alphaX/3.0/sqrt(PI/2),1.0/3.0);
      }
      else if (IonType[i] == 3) {
        Rqi= R_const*1.889725989 ;
        alphaX= alpha_parC/(0.529177*0.529177*0.529177);
        alphaY= alpha_parC/(0.529177*0.529177*0.529177);
        alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);
        Rpi= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
       // Rpi= pow(alphaX/3.0/sqrt(PI/2),1.0/3.0);
      }
      if(i==j){
      //   cout <<"alphaX : "<<alphaX<<endl;
      //   cout <<"alphaY : "<<alphaX<<endl;
      //   cout <<"alphaZ : "<<alphaX<<endl;
      


          A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = 1.0/alphaX ;
          A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = 1.0/alphaY ;
          A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = 1.0/alphaZ ;


       // cout<<" i=j 1 : "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+0 <<endl;
      //  cout<<" i=j 2 : "<<n*n4plus1+3*i*n4plus1 + n4plus1+         n+j*3+1 <<endl;
      //  cout<<" i=j 3 : "<<n*n4plus1+3*i*n4plus1 + n4plus1+ n4plus1+n+j*3+2 <<endl;
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq, fac_qp; 
        double Txy, Txz, Tyz, Rpq, Rqp, Rpp ;
        Rpq = sqrt(Rpi*Rpi+Rqj*Rqj)  ;
        Rqp = sqrt(Rqi*Rqi+Rpj*Rpj)  ;
        Rpp = sqrt(Rpi*Rpi+Rpj*Rpj)  ;

      //  cout<<"Rpi Rpj = "<<Rpi<<" "<<Rpj<<endl;
      //  cout<<"Rpp Rqp = "<<Rpp<<" "<<Rqp<<endl;

        fac = erf(Rij/Rpp)-2.0*Rij*exp(-1.0*(Rij/Rpp)*(Rij/Rpp))/(sqrt(PI)*Rpp) ;
        fac_pq = (erf(Rij/Rpq)-2.0*Rij*exp(-1.0*(Rij/Rpq)*(Rij/Rpq))/(sqrt(PI)*Rpq))/(Rij*Rij*Rij) ;
        fac_qp = (erf(Rij/Rqp)-2.0*Rij*exp(-1.0*(Rij/Rqp)*(Rij/Rqp))/(sqrt(PI)*Rqp))/(Rij*Rij*Rij) ;
       // fac_pq = fac/(Rij*Rij*Rij) ; 
       // cout<<"Rij = "<<Rij<<endl;
        //cout<<"Rqq = "<<Rqq<<endl;

        fac1 = fac/(Rij*Rij*Rij*Rij*Rij) ; 

        fac2 = 4.0*exp(-1.0*(Rij/Rpp)*(Rij/Rpp))/(sqrt(PI)*Rpp*Rpp*Rpp*Rij*Rij) ;
        //cout<<"fac_qp = "<<fac_qp<<endl;
 
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;







        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 


        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+1] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+2] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2] = -1.0*Tyz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1] = -1.0*Tyz ;

         

//        cout<<" Txy: "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+1<<endl; 
//        cout<<" Txy: "<<n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0<<endl; 
//        cout<<" Txz: "<<n*n4plus1+3*i*n4plus1 +                  n+j*3+2<<endl; 
//        cout<<" Txz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0<<endl; 
//        cout<<" Tyz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2<<endl; 
//        cout<<" Tyz: "<<n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1<<endl;




        //Tpq
       

        
//        A[i*n4plus1+n+3*j+0] = 0;
//        A[i*n4plus1+n+3*j+1] = 0;
//        A[i*n4plus1+n+3*j+2] = 0;


//        A[n*n4plus1+i*3*n4plus1                  +j] = 0;
//        A[n*n4plus1+i*3*n4plus1+ n4plus1         +j] = 0;
//        A[n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j] = 0;


       


//        if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
//        A[i*n4plus1+n+3*j+0] = -1.0*rij[0]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[i*n4plus1+n+3*j+1] = -1.0*rij[1]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[i*n4plus1+n+3*j+2] = -1.0*rij[2]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        }
//        else {



        A[i*n4plus1+n+3*j+0] = -1.0*rij[0]*fac_qp;
        A[i*n4plus1+n+3*j+1] = -1.0*rij[1]*fac_qp;
        A[i*n4plus1+n+3*j+2] = -1.0*rij[2]*fac_qp;


//        }




//        A[n*n4plus1+i*3*n4plus1                  +j] = rij[0]*fac_pq*-1.0;
//        A[n*n4plus1+i*3*n4plus1+ n4plus1         +j] = rij[1]*fac_pq*-1.0;
//        A[n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j] = rij[2]*fac_pq*-1.0;





//       if ((i < 60 && j >= 60) || (i >= 60 && j < 60) )  {
//        A[n*n4plus1+j*3*n4plus1                  +i] = -1.0*rij[0]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[n*n4plus1+j*3*n4plus1+ n4plus1         +i] = -1.0*rij[1]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        A[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i] = -1.0*rij[2]*fac_pq*(1.0 - exp(ChargePen*Rij));
//        }
//        else {



        A[n*n4plus1+j*3*n4plus1                  +i] = -1.0*rij[0]*fac_qp;
        A[n*n4plus1+j*3*n4plus1+ n4plus1         +i] = -1.0*rij[1]*fac_qp;
        A[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i] = -1.0*rij[2]*fac_qp;



//        }




//        cout<<"first A["<<i*n4plus1+n+3*j+0<<"] = "<<rij[0]*fac_qp*-1.0<<endl;
//        cout<<"first A["<<i*n4plus1+n+3*j+1<<"] = "<<rij[1]*fac_qp*-1.0<<endl;
//        cout<<"first A["<<i*n4plus1+n+3*j+2<<"] = "<<rij[2]*fac_qp*-1.0<<endl;
//        cout<<"second A["<<n*n4plus1+j*3*n4plus1                  +i<<"]="<<rij[0]*fac_qp*-1.0<<endl;
//        cout<<"second A["<<n*n4plus1+j*3*n4plus1+ n4plus1         +i<<"]="<<rij[1]*fac_qp*-1.0<<endl;
//        cout<<"second A["<<n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i<<"]="<<rij[2]*fac_qp*-1.0<<endl;

//        cout<<"second A["<<n*n4plus1+i*3*n4plus1                  +j<<"]="<<rij[0]*fac_qp*-1.0<<endl;
//        cout<<"second A["<<n*n4plus1+i*3*n4plus1+ n4plus1         +j<<"]="<<rij[1]*fac_qp*-1.0<<endl;
//        cout<<"second A["<<n*n4plus1+i*3*n4plus1+ n4plus1+ n4plus1+j<<"]="<<rij[2]*fac_qp*-1.0<<endl;


//        cout<<"Tqx("<<i<<","<<j<<": "<<i*n4plus1+n+3*j+0 <<" :  "<<  rij[0]*fac_pq*-1.0 <<endl;
//        cout<<"Tqy : "<<i*n4plus1+n+3*j+1 <<" :  "<<  rij[1]*fac_pq*-1.0 <<endl;
//        cout<<"Tqz : "<<i*n4plus1+n+3*j+2 <<" :  "<<  rij[2]*fac_pq*-1.0 <<endl;
//        cout<<"Txq : "<<n*n4plus1+j*3*n4plus1                  +i<<" :  "<<  rij[0]*fac_pq*-1.0 <<endl;
//        cout<<"Tyq : "<<n*n4plus1+j*3*n4plus1+ n4plus1         +i<<" :  "<<  rij[1]*fac_pq*-1.0 <<endl;
//        cout<<"Tzq : "<<n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i<<" :  "<<  rij[2]*fac_pq*-1.0 <<endl;
 
      }
     
    }
  }








/*
  for (int i = 0; i < n4plus1*n4plus1; ++i){
     cout<<" A["<<i<<"]= "<< A[i]<<endl;;
  }
*/

// Callculate the Hessian matrix diagonalization
  char jobz = 'V';
  char uplo = 'L';
  double* eigenvectors = new double[nplus1*nplus1];
  double* eigenvalues = new double[nplus1];

/*
  for (int i = 0; i < nplus1; ++i){
    cout<<i<<" : "; 
    for (int j = 0; j < nplus1; ++j){
     cout<< Matrix[i*nplus1+j]<<" ";
    }
    cout<<endl;
  }
*/
/*
  for (int i = 0; i < nplus1; ++i)
    memcpy(eigenvectors+i*(n+1), Matrix+(i+1)*(2*n-i)/2-(n-i),(n-i)*sizeof(double));
  int lwork = 12 * n;
  double* work = new double[lwork];
  int info=0;
  dsyev(&jobz, &uplo, &n, eigenvectors, &n, eigenvalues, work, &lwork, &info);

  for (int i = 0; i < nplus1; ++i){
     cout<<" eigenvalues["<<i<<"]= "<< eigenvalues[i]<<endl;;
  }

  if (info != 0) {
    cout << "DMDCBase::CallDSYEV: DSYEV returned with info = " << info << endl;
    exit(1);
  }
  delete[] work;
*/
  delete[] eigenvalues;
  delete[] eigenvectors;
  delete[] Matrix;






/*
  for (int i = 0; i < n4plus1*n4plus1; ++i){
     A[i]=0.0;
   }
*/

//  int icount = 0;
//   for (int i = 0; i < n4plus1; ++i)
//    for (int j = 0; j < n4plus1; ++j) {
//     cout<<"symmetry check ["<<icount<<"], ("<<i<<","<<j<<")=" << " : "<< n4plus1*i+j<<","<<i+n4plus1*j<<" : "<<A[n4plus1*i+j]- A[i+n4plus1*j]<<endl;
//      icount++;
//    }
}

/////////////////////////////////////////////////
//
//  symmetric matrix A is inverted in place
//
void InvertCDMatrix(int n, double *A)
{

  int N = n;
  char uplo = 'U';
  static iVec ipiv; ipiv.resize(n);    // this is called repeatedly to invert 9x9 matrices, but only during setup, so
  static dVec work; work.resize(n*n);  // static vectors should be fine here as long as nobody calls this in mutiple threads
  int lwork = n*n;
  int info;
 // static dVec B; B.resize(n*n);
  double* B= new double[n*n];


  for (int i = 0; i < n; ++i)
   ipiv[i]=0;

 cout<<"Matrix =============" <<endl;
/*
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
     cout<<" "<< A[n*i+j];
    }
   cout<<endl;
 }
*/


  cout<<"dimension n = "<<n<<endl;
  cout<<"dimension N = "<<N<<endl;
  dgetrf(&N, &N, A, &n, &ipiv[0], &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrf = " << info << "\n";
    exit(1);
  }

  dgetri(&N, A, &n, &ipiv[0], &work[0],  &lwork, &info);
  if (info != 0) {
    cout << "InverseMatrix, Error in dgetrfi = " << info << "\n";
    exit(1);
  }



//
//  for (int i = 0; i < n; ++i)
//   for (int j = 0; j < n; ++j)
//     cout<<i<<" "<<j<<" :InvMatrix  = " << A[n*i+j]<<endl;
    // cout<<"2 symmetry check =" << A[n*i+j]- A[i+n*j]<<endl;
     //if (abs(A[n*i+j]- A[i+n*j]) >= 0.00001 ) cout<<"symmetry check =" << A[n*i+j]- A[i+n*j]<<endl;
  delete [] B;

}



//////////////////////////////////////////////////////////////
//
//  build matrix A = (alpha^-1 - T)
//
//  this is a super matrix consisting of:
//  
//  T(qq) T(qd) T(1)
//  T(qd) T(dd) T(1)
//  T(1)  T(1)  0 
//
void BuildPolarizationMatrix3(int nSites, const double *R, double *A, int nFullerenes, const int *NoAtomsArray, int CarbonType, int *IonType)
{

  int n = nSites;
  int n4 = n*4 ;
  int n4plus1 = 4*n + nFullerenes ; // dimension of super-matrix A and leading dimension of any block B
  int nplus1 = n + nFullerenes ; 
  int ncomb = n*(n-1)/2-1 ; 
  int ncplus1 = ncomb+1 ; 

   cout<<"n4plus1 = "<<n4plus1<<endl;
   cout<<"nFullerenes = "<<nFullerenes<<endl;
   cout<<"CarbonType = "<<CarbonType<<endl;
// for (int k= 0; k<n ; ++k )
//   cout<<"IonType["<<k<<"] = "<<IonType[k]<<endl;

//  for (int k = 0; k < n4plus1*n4plus1 ; ++k)
//    A[k] = 0.0;
  for (int k = 0; k < nplus1*nplus1 ; ++k)
    A[k] = 0.0;
  
 // double* Matrix= new double[ncomb*ncomb];
 // dVec Matrix; Matrix.resize(ncomb*ncomb);
  double Matrix[nplus1*nplus1];
  double ChargePen=-1.5;
  double Rqi, Rqj, chi; 
  
  for (int k = 0; k < nplus1*nplus1 ; ++k)
    Matrix[k] = 0.0;


  //Tqq
  
  double R_const;
   if (CarbonType == 2) R_const=0.01945;
   else if (CarbonType == 3) R_const=0.02279;
   else if (CarbonType == 4) R_const=0.02652;
   else if (CarbonType == 5) R_const=0.03029;
   else if (CarbonType == 6) R_const=0.215;
   else if (CarbonType == 7) R_const=0.01838;
   else if (CarbonType == 8) R_const=0.01758;

  for (int j = 0; j < n; ++j) {
    if (IonType[j] == 1) Rqj= 0.68620399*1.889725989 ;
    else if (IonType[j] == 2) Rqj= 0.01758*1.889725989 ;
    else if (IonType[j] == 3) Rqj= R_const*1.889725989 ;

    for (int i = 0; i < n; ++i) {
      if (IonType[i] == 1) Rqi= 0.68620399*1.889725989 ;
      else if (IonType[j] == 2) Rqi= 0.01758*1.889725989 ;
      else if (IonType[i] == 3) Rqi= R_const*1.889725989 ;
      chi = sqrt(2.0/PI)/Rqi  ;
      if (i == j) {
        A[nplus1*i+j] = chi; 
        Matrix[nplus1*i+j] = chi; 
      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double CT=4.0*pow(10,4)*exp(-1.5*Rij);
        double Rqq = sqrt(Rqi*Rqi+Rqj*Rqj)  ;

         A[nplus1*i+j] = erf(Rij/Rqq)/Rij ; 
         Matrix[nplus1*i+j] = erf(Rij/Rqq)/Rij ; 
      }
    }
  }


 int halfn = n/2.0;
 int icount=0;
 int comb_count=0;
// double *B = A;

/*
  for (int i = 0; i < n; ++i) {
    for (int j=i+1 ; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          for (int l=k+1 ; l < n; ++l) {
          // if (k <60 && l >=60 ) comb_count++ ;
            cout<<"icount ncplus icount%ncplus1:  " <<icount<<" "<<ncplus1<<" "<<icount%ncplus1<<endl;
            if (icount%ncplus1 == ncplus1-1 ){
               if (i <halfn && j >=halfn ) B[icount++]=1.0 ;
               else B[icount++]=0.0 ;
            }
            else B[icount++]= Matrix[nplus1*i+k]-Matrix[nplus1*i+l]-Matrix[nplus1*j+k]+Matrix[nplus1*j+l]; 

            cout<<"break: "<<icount<<" "<<ncplus1*(ncplus1-1) <<endl;
            if (icount >= ncplus1*(ncplus1-1)) {
              cout<<"brake: "<<ncplus1*(ncplus1-1) <<endl; 
              goto done; 
            }
          // A[icount ] = Matrix[nplus1*i+k]-Matrix[nplus1*i+l]-Matrix[nplus1*j+k]+Matrix[nplus1*j+l]; 
          }
        }
    }
  }

done:

 cout<<" icount = "<<icount<<endl;
 cout<<" comb_count = "<<comb_count<<endl;


  for (int k = 0; k < n; ++k) 
    for (int l=k+1 ; l < n; ++l) { 
        cout<<" icount k l "<<icount<<" "<<k<<" "<<l<<endl;
       if (k <halfn && l >=halfn ) B[icount++]=1.0 ;
       else B[icount++]=0.0;
    }
*/


  //T-lagrange
  int current_atom = 0 ; 
  int old_atom     = 0 ;    
  double Aconst= 0.0000000;
  double Bconst= 0.000000; 
  for ( int iFuller = 0 ; iFuller < nFullerenes ; ++iFuller){
    old_atom      = current_atom ;
    current_atom += NoAtomsArray[iFuller] ;
    for (int i = old_atom; i < current_atom ; ++i) {
      A[i*nplus1 +n + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
      A[nplus1*n + nplus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
  // cout<<"i*nplus1 +n + iFuller " <<i*nplus1 +n + iFuller <<endl;
  // cout<<"nplus1*n + nplus1*iFuller + i " <<nplus1*n + nplus1*iFuller + i <<endl;
    //  A[i*n4plus1+n4 + iFuller]               = 1.0 ; // corresponds to columns of lagrange vector
     // A[n4plus1*n4 + n4plus1*iFuller + i]     = 1.0 ; // corresponds to rows of lagrange rows
    }
  }
 



/*
  double alpha, alphaX, alphaY, alphaZ;
  double Rpi, Rpj;
  double alpha_parC, alpha_perpC, alpha_parH, alpha_perpH;


   if (CarbonType == 2) { 
      R_const=0.01945;
      alpha_parC = 1.2885;
      alpha_perpC = 1.2885;
      alpha_parH = 0.3384 ;
      alpha_perpH = 0.3384 ;
   }
   else if (CarbonType == 3) {
      R_const=0.02279;
      alpha_parC = 1.4755;
      alpha_perpC = 2.2611;
      alpha_parH = 0.3057 ;
      alpha_perpH = 0.3057 ;
   }
   else if (CarbonType == 4) { 
      R_const=0.02652;
      alpha_parC = 1.2080;
      alpha_perpC = 1.2080;
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 5) {
       R_const=0.03029;
      alpha_parC = 1.3632;
      alpha_perpC = 2.1671; 
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }
   else if (CarbonType == 6) {
      R_const=0.215;
      alpha_parC = 1.317;
      alpha_perpC = 1.84 ; 
      alpha_parH = 0.4471 ;
      alpha_perpH = 0.4471 ;
   }


  //Tpq and Tpp
  for (int j = 0; j < n; ++j){
    if (IonType[j] == 1) { 
      Rqj= 0.68620399*1.889725989 ; 
      Rpj=Rqj;
      alphaX = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
      alphaY = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
      alphaZ = 3.0*sqrt(PI/2)*Rpj*Rpj*Rpj ;
    }
    else if (IonType[j] == 2) {
      Rqj= 0.00 ;
      alphaX= alpha_parH/(0.529177*0.529177*0.529177);
      alphaY= alpha_parH/(0.529177*0.529177*0.529177);
      alphaZ= alpha_perpH/(0.529177*0.529177*0.529177);
      Rpj= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
    }
    else if (IonType[j] == 3){
      Rqj= R_const*1.889725989 ; 
      alphaX= alpha_parC/(0.529177*0.529177*0.529177);
      alphaY= alpha_parC/(0.529177*0.529177*0.529177);
      alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);
      Rpj= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
    }

    for (int i = 0; i < n; ++i){
      if (IonType[i] == 1) {
        Rqi= 0.68620399*1.889725989 ; 
        Rpi=Rqi;
        alphaX = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
        alphaY = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
        alphaZ = 3.0*sqrt(PI/2)*Rpi*Rpi*Rpi ;
       }
      else if (IonType[i] == 2) {
        Rqi= 0.00 ;
        alphaX= alpha_parH/(0.529177*0.529177*0.529177);
        alphaY= alpha_parH/(0.529177*0.529177*0.529177);
        alphaZ= alpha_perpH/(0.529177*0.529177*0.529177);
        Rpi= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
      }
      else if (IonType[i] == 3) {
        Rqi= R_const*1.889725989 ;
        alphaX= alpha_parC/(0.529177*0.529177*0.529177);
        alphaY= alpha_parC/(0.529177*0.529177*0.529177);
        alphaZ= alpha_perpC/(0.529177*0.529177*0.529177);
        Rpi= pow(sqrt(2.0/PI)/(2.0/alphaX +1.0/alphaZ),1.0/3.0);
      }
      if(i==j){
        // cout <<"alphaX : "<<alphaX<<endl;
       //  cout <<"alphaY : "<<alphaX<<endl;
       //  cout <<"alphaZ : "<<alphaX<<endl;

          A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = 1.0/alphaX ;
          A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = 1.0/alphaY ;
          A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = 1.0/alphaZ ;

      }
      else {
        double rij[3] ;
        rij[0] = R[3*i+0] - R[3*j+0];
        rij[1] = R[3*i+1] - R[3*j+1];
        rij[2] = R[3*i+2] - R[3*j+2];
        double Rij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
        double fac, fac1, fac2, fac_pq, fac_qp; 
        double Txy, Txz, Tyz, Rpq, Rqp, Rpp ;
        Rpq = sqrt(Rpi*Rpi+Rqj*Rqj)  ;
        Rqp = sqrt(Rqi*Rqi+Rpj*Rpj)  ;
        Rpp = sqrt(Rpi*Rpi+Rpj*Rpj)  ;


        fac = erf(Rij/Rpp)-2.0*Rij*exp(-1.0*(Rij/Rpp)*(Rij/Rpp))/(sqrt(PI)*Rpp) ;
        fac_pq = (erf(Rij/Rpq)-2.0*Rij*exp(-1.0*(Rij/Rpq)*(Rij/Rpq))/(sqrt(PI)*Rpq))/(Rij*Rij*Rij) ;
        fac_qp = (erf(Rij/Rqp)-2.0*Rij*exp(-1.0*(Rij/Rqp)*(Rij/Rqp))/(sqrt(PI)*Rqp))/(Rij*Rij*Rij) ;

        fac1 = fac/(Rij*Rij*Rij*Rij*Rij) ; 

        fac2 = 4.0*exp(-1.0*(Rij/Rpp)*(Rij/Rpp))/(sqrt(PI)*Rpp*Rpp*Rpp*Rij*Rij) ;
        //cout<<"fac_qp = "<<fac_qp<<endl;
 
        Txy = 3.0*rij[0]*rij[1]*fac1 - rij[0]*rij[1]*fac2 ;
        Txz = 3.0*rij[0]*rij[2]*fac1 - rij[0]*rij[2]*fac2 ;
        Tyz = 3.0*rij[1]*rij[2]*fac1 - rij[1]*rij[2]*fac2 ;



        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+0] = -1.0*((3.0*rij[0]*rij[0]-Rij*Rij)*fac1 - rij[0]*rij[0]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+1] = -1.0*((3.0*rij[1]*rij[1]-Rij*Rij)*fac1 - rij[1]*rij[1]*fac2) ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+2] = -1.0*((3.0*rij[2]*rij[2]-Rij*Rij)*fac1 - rij[2]*rij[2]*fac2) ; 


        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+1] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+0] = -1.0*Txy ; 
        A[n*n4plus1+3*i*n4plus1 +                  n+j*3+2] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+0] = -1.0*Txz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+          n+j*3+2] = -1.0*Tyz ; 
        A[n*n4plus1+3*i*n4plus1 +n4plus1+ n4plus1+ n+j*3+1] = -1.0*Tyz ;

         
        //Tpq
       

        A[i*n4plus1+n+3*j+0] = -1.0*rij[0]*fac_qp;
        A[i*n4plus1+n+3*j+1] = -1.0*rij[1]*fac_qp;
        A[i*n4plus1+n+3*j+2] = -1.0*rij[2]*fac_qp;


        A[n*n4plus1+j*3*n4plus1                  +i] = -1.0*rij[0]*fac_qp;
        A[n*n4plus1+j*3*n4plus1+ n4plus1         +i] = -1.0*rij[1]*fac_qp;
        A[n*n4plus1+j*3*n4plus1+ n4plus1+ n4plus1+i] = -1.0*rij[2]*fac_qp;

      }
     
    }
  }
*/

/*
  for (int i = 0; i < ncomb; ++i)
    for (int j = i+1; j < ncomb; ++j)
     cout<<"symmetry check =" << A[ncomb*i+j]- A[i+ncomb*j]<<endl;
 */   
     
/*
    for (int i = 0; i < ncomb*ncomb; ++i){
         cout<<" A["<<i<<"]= "<< A[i]<<endl;;
   }
*/


// Callculate the Hessian matrix diagonalization
  char jobz = 'V';
  double* eigenvectors = new double[nplus1*nplus1];
  double* eigenvalues = new double[nplus1];

  delete[] eigenvalues;
  delete[] eigenvectors;
 // delete[] Matrix;


}

void InverseMatrix(double* Input, double* Output,  int n )
{

  cout<<" n ===="<<n<<endl;
  int m =n;
  double **a = new double*[n];
  double **b = new double*[n];
  for (int ii=0; ii<n; ++ii){
    a[ii]= new double[n];
    b[ii]= new double[m];
  }

  for(int i=0; i<n*(n+1)/2; ++i){ 
   printf("Input %10.10f \n", Input[i]);
 }
  int index=0;
  for (int ii=0; ii< n; ++ii){
    for (int ij=ii; ij< n; ++ij){
      a[ii][ij] = Input[index++];
      if (ii != ij) a[ij][ii]=a[ii][ij];
    }
  }

/*

   I'm going to find the Inverse Matrix using Gauss-Jordan Elimination method
   this source code came from Numerical recipes: the art of scientific computing
   by William H. Press
*/

  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv;
  int *indxc, *indxr, *ipiv;
  indxc = new int[n];
  indxr = new int[n];
  ipiv = new int[n];
  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++){
          big=0.0;
          for (j=0;j<n;j++)
                  if (ipiv[j] !=1)
                          for (k=0;k<n;k++){
                                  if (ipiv[k]==0){
                                          if (abs(a[j][k]) >= big){
                                                  big=abs(a[j][k]);
                                                  irow=j;
                                                  icol=k;
                                          }
                                  }
                          }
          ++(ipiv[icol]);
          if( irow != icol){
                  for(l=0;l<n;l++) swap(a[irow][l],a[icol][l]);
                  for(l=0;l<m;l++) swap(b[irow][l],b[icol][l]);
          }
          indxr[i]=irow;
          indxc[i]=icol;
          if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
          pivinv=1.0/a[icol][icol];
          a[icol][icol]=1.0;
          for (l=0;l<n;l++) a[icol][l] *= pivinv;
          for (l=0;l<m;l++) b[icol][l] *= pivinv;
          for (ll=0;ll<n;ll++)
                  if(ll != icol){
                          dum=a[ll][icol];
                          a[ll][icol]=0.0;
                          for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                          for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
                  }

  }
  for(l=n-1; l>=0;l--){
          if (indxr[l] != indxc[l])
                  for (k=0;k<n;k++)
                          swap(a[k][indxr[l]],a[k][indxc[l]]);
  }

  index=0;
  for (int ii=0; ii< n; ++ii){
    for (int ij=ii; ij< n; ++ij){
      Output[index++]=a[ii][ij];
      cout<<ii<<" "<<ij<<" Output =["<<index-1<<"]="<<Output[index-1]<<endl;
    }
  }

 delete [] a;
 delete [] b;

}

double d=0.0;
double det(int n, double mat[10][10])
{
    int c, subi, i, j, subj;
    double submat[10][10];  

 //   for (int i=0; i <n ; i++)
 //    for (int j=0; j <n ; j++)
 //      cout<<"mat["<<i<<"]["<<j<<"] "<<mat[i][j]<<endl;

    if (n == 2) 
    {
        return( (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
    }
    else
    {  
        for(c = 0; c < n; c++)
        {  
            subi = 0;  
            for(i = 1; i < n; i++)
            {  
                subj = 0;
                for(j = 0; j < n; j++)
                {    
                    if (j == c)
                    {
                        continue;
                    }
                    submat[subi][subj] = mat[i][j];
                    subj++;
                }
                subi++;
            }
        d = d + (pow(-1 ,c) * mat[0][c] * det(n - 1 ,submat));
        }
    }
    return d;
}

