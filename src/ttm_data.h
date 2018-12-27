
//
//   here are the parameters of the Dang-Chang water model
//   see JCP 106, 8149 (1997)
// 
//

const double TTM_ROH = 0.9572;      // OH bond length
const double TTM_HOH = 104.52;      // bond angle


//
// point charges on H and on the M-site
//
const double TTM_QH = 0.5742185854046506;
const double TTM_QM = -2 * TTM_QH;

const double KT_QH = 0.56964;     // This value comes from Seymour, and
const double KT_QM = -2 * KT_QH;  // the negative charge goes on the DC M-site


//DC m-site
const double DC_QH = 0.5190;
const double DC_QM = -2*DC_QH;

//
//  This is the choice of reference system
//  this has to be the same as in Water::Reposition()
//
// Center of mass system  [Angs]
const double TTM_Oz =  0.0650980307;    
const double TTM_Hz = -0.5207842453;
const double TTM_Hy =  0.756950327;
//
//  Gaussian convention (center of charge)
//const double TTM_Oz =  0.117176;    
//const double TTM_Hz = -0.468706;
//const double TTM_Hy =  0.756950;
//
//  Oxygen as origin
//const double TTM_Oz =  0.0;    
//const double TTM_Hz = -0.585882276;
//const double TTM_Hy =  0.756950327;
//
// position of the M-site
// 
const double DPP_DeltaM = 0.25;
const double TTM_Mz = TTM_Oz - DPP_DeltaM;
const double DC_Mz  = TTM_Oz - 0.215;

//
//  polarizable sites on O and H
//
const int TTM_NoOfPolSites = 3;
//const double TTM_POL_O = 0.837*pow(Angs2Bohr,3);
//const double TTM_POL_H = 0.496*pow(Angs2Bohr,3);
const double TTM_POL_O = 1.22*pow(Angs2Bohr,3);
const double TTM_POL_H = 0.28*pow(Angs2Bohr,3);
//const double TTM_POL_H = 0.0;
//const double TTM_POL_O = 0.0;

//these constants are for the damping polarizabitlity 
const double TTM_alpha_o = 0.837*pow(Angs2Bohr,3);
const double TTM_alpha_h = 0.496*pow(Angs2Bohr,3);
//const double TTM_POL_H = 0.496;

//DC
const int DC_NoOfPolSites = 1;
const double DC_POL_M = 1.444*pow(Angs2Bohr,3);



// AAD 4-11-06
//	all-atom repulsion and damped dispersion (ttm2-u6)
//
const double u6_delta = 2.23/Angs2Bohr;
const double u6_c6 = -1300.0000000/AU2KCAL*pow(Angs2Bohr,6);

// OO
const double u6_A1 = 150400.0/AU2KCAL;
const double u6_b1 = 4.231/Angs2Bohr;

// HH
const double u6_A2 = 24252.61091436135/AU2KCAL;
const double u6_b2 = 4.7393090112470775/Angs2Bohr;

// OH
const double u6_A3 = 2897.3950768710165/AU2KCAL;
const double u6_b3 = 3.784608449684255/Angs2Bohr;

/*
// AAD 4-11-06
//	all-atom repulsion and damped dispersion
// OO
const double DPP_A1 = 150400.0/627.50956;
const double DPP_b1 = 4.231/1.88972613392;

// HH
const double DPP_A2 = 24252.61091436135/627.50956;
const double DPP_b2 = 4.7393090112470775/1.88972613392;

// OH
const double DPP_A3 = 2897.3950768710165/627.50956;
const double DPP_b3 = 3.784608449684255/1.88972613392;
*/

// AAD 6-6-07
//   REAL DPP parameters
// OO
const double DPP_A1 = 115271.77066806905/AU2KCAL;
const double DPP_b1 = 4.167868773196562/Angs2Bohr;

// HH
const double DPP_A2 = 119085.35806729361/AU2KCAL;
const double DPP_b2 = 5.581698710292979/Angs2Bohr;

// OH
const double DPP_A3 = 1803.662375255084/AU2KCAL;
const double DPP_b3 = 3.5354602029164868/Angs2Bohr;


//DC parameters
const double DC_c6  = -835.14705414204481800768/AU2KCAL*pow(Angs2Bohr,6);
const double DC_c12 = 955439.18087963772412374922/AU2KCAL*pow(Angs2Bohr,12);

const double DC_sigma = 3.2340*Angs2Bohr;
const double DC_eps = 0.1825/AU2KCAL;


//
//  isotropic polarizability of water for the Drude oscillators
//
const double W_ALPHA_EXP = 9.745;
//const double W_ALPHA_EXP = 9.25775;
const double W_ALPHA_HF  = 8.534;
