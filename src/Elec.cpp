

double Potential::EvaluateCoroleneElecPot(const double *x)
{ 

// progress_timer tmr("EvaluateColroleneElecPot:", verbose);

  double Velec = 0.0 ;


  for (int i =0 ; i < nSites; ++i){
    Q00[i] = 
    PtDip[3*i+0]=
    PtDip[3*i+1]=
    PtDip[3*i+2]=
  
  for (int i =0 ; i < nSites; ++i){
    // effective-r damping
    double DampParameter = EsDampCsixty ;
    double Reff = R[i];
    double gij = Rminus3[i] ; 
    if (Reff < DampParameter) {
      double ror0 = Reff / DampParameter;
      Reff = DampParameter * (0.5 + ror0 * ror0 * ror0 * (1.0 - 0.5 * ror0));
      // old not quite so smooth fn: Reff = 0.5 * (Reff*Reff/DampParameter + DampParameter);
      gij = 1.0 / (Reff * Reff * Reff);
    }

    Velec += -1.0*ExtraNucQ/Reff + (Rx[i]*PtDip[3*i+0] + Ry[i]*PtDip[3*i+1] + Rz[i]*PtDip[3*i+2])*gij ;
    
  }
  return Velec ;
   
}

double Potential::EvaluateCsixtyRepPot(const double *x, double GaussExp, double RepScale)
{
// progress_timer tmr("EvaluateCsixtyRepPot", verbose);
  double Vrep = 0.0 ;
  switch (RepFlagCsixty)
    {
      case 0:
        Vrep = 0.0 ; //no repulsion
      break;
      case 1:
        for (int i = 0; i < nSites; ++i) {
        Vrep +=  RepScale*exp(-GaussExp * R2[i]); //gaussians
          }
      break ;
      default:
      printf("EvaluateCsixtyRep: This should never happen\n");
      exit(42);

    }
  return Vrep ; 
}


