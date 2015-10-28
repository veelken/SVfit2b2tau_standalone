#include "TauAnalysis/SVfit2b2tau_standalone/interface/BJetEtResolution.h"

#include <TMath.h>

#include <iostream>
#include "math.h"

using namespace svFitStandalone2b2tau;

BJetEtResolution::BJetEtResolution()
{}

BJetEtResolution::~BJetEtResolution()
{}

namespace
{
  // b-jet energy resolution taken from HHKinFit package:
  //   https://github.com/bvormwald/HHKinFit

  double getBjetEtResolution(double et, double eta)
  {
    //std::cout << "<getBjetEtResolution>: Et = " << et << ", eta = " << eta << std::endl;

    double det = 0.;
    //double de = 10.;

    if ( 0.000 <= fabs(eta) && fabs(eta) < 0.087 ) {
      det = et * (sqrt(0.0686*0.0686 + (1.03/sqrt(et))*(1.03/sqrt(et)) + (1.68/et)*(1.68/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.000+0.087)/2))) * det;
    }
    
    if ( 0.087 <= fabs(eta) && fabs(eta) < 0.174 ) {
      det = et * (sqrt(0.0737*0.0737 + (1.01/sqrt(et))*(1.01/sqrt(et)) + (1.74/et)*(1.74/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.087+0.174)/2))) * det;
    }
    
    if ( 0.174 <= fabs(eta) && fabs(eta) < 0.261 ) {
      det = et * (sqrt(0.0657*0.0657 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (5.16e-06/et)*(5.16e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.174+0.261)/2))) * det;
    }
    
    if ( 0.261 <= fabs(eta) && fabs(eta) < 0.348 ) {
      det = et * (sqrt(0.062*0.062 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (0.000134/et)*(0.000134/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.261+0.348)/2))) * det;
    }
    
    if ( 0.348 <= fabs(eta) && fabs(eta) < 0.435 ) {
      det = et * (sqrt(0.0605*0.0605 + (1.07/sqrt(et))*(1.07/sqrt(et)) + (1.84e-07/et)*(1.84e-07/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.348+0.435)/2))) * det;
    }
    
    if ( 0.435 <= fabs(eta) && fabs(eta) < 0.522 ) {
      det = et * (sqrt(0.059*0.059 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (9.06e-09/et)*(9.06e-09/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.435+0.522)/2))) * det;
    }
    
    if ( 0.522 <= fabs(eta) && fabs(eta) < 0.609 ) {
      det = et * (sqrt(0.0577*0.0577 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (5.46e-06/et)*(5.46e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.522+0.609)/2))) * det;
    }
    
    if ( 0.609 <= fabs(eta) && fabs(eta) < 0.696 ) {
      det = et * (sqrt(0.0525*0.0525 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (4.05e-05/et)*(4.05e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.609+0.696)/2))) * det;
    }
    
    if ( 0.696 <= fabs(eta) && fabs(eta) < 0.783 ) {
      det = et * (sqrt(0.0582*0.0582 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.17e-05/et)*(1.17e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.696+0.783)/2))) * det;
    }
    
    if ( 0.783 <= fabs(eta) && fabs(eta) < 0.870 ) {
      det = et * (sqrt(0.0649*0.0649 + (1.08/sqrt(et))*(1.08/sqrt(et)) + (7.85e-06/et)*(7.85e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.783+0.870)/2))) * det;
    }
    
    if ( 0.870 <= fabs(eta) && fabs(eta) < 0.957 ) {
      det = et * (sqrt(0.0654*0.0654 + (1.1/sqrt(et))*(1.1/sqrt(et)) + (1.09e-07/et)*(1.09e-07/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.870+0.957)/2))) * det;
    }
    
    if ( 0.957 <= fabs(eta) && fabs(eta) < 1.044 ) {
      det = et * (sqrt(0.0669*0.0669 + (1.11/sqrt(et))*(1.11/sqrt(et)) + (1.87e-06/et)*(1.87e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(0.957+1.044)/2))) * det;
    }
    
    if ( 1.044 <= fabs(eta) && fabs(eta) < 1.131 ) {
      det = et * (sqrt(0.0643*0.0643 + (1.15/sqrt(et))*(1.15/sqrt(et)) + (2.76e-05/et)*(2.76e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.044+1.131)/2))) * det;
    }
    
    if ( 1.131 <= fabs(eta) && fabs(eta) < 1.218 ) {
      det = et * (sqrt(0.0645*0.0645 + (1.16/sqrt(et))*(1.16/sqrt(et)) + (1.04e-06/et)*(1.04e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.131+1.218)/2))) * det;
    }
    
    if ( 1.218 <= fabs(eta) && fabs(eta) < 1.305 ) {
      det = et * (sqrt(0.0637*0.0637 + (1.19/sqrt(et))*(1.19/sqrt(et)) + (1.08e-07/et)*(1.08e-07/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.218+1.305)/2))) * det;
    }
    
    if ( 1.305 <= fabs(eta) && fabs(eta) < 1.392 ) {
      det = et * (sqrt(0.0695*0.0695 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5.75e-06/et)*(5.75e-06/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.305+1.392)/2))) * det;
    }
    
    if ( 1.392 <= fabs(eta) && fabs(eta) < 1.479 ) {
      det = et * (sqrt(0.0748*0.0748 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (5.15e-08/et)*(5.15e-08/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.392+1.479)/2))) * det;
    }
    
    if ( 1.479 <= fabs(eta) && fabs(eta) < 1.566 ) {
      det = et * (sqrt(0.0624*0.0624 + (1.23/sqrt(et))*(1.23/sqrt(et)) + (2.28e-05/et)*(2.28e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.479+1.566)/2))) * det;
    }
    
    if ( 1.566 <= fabs(eta) && fabs(eta) < 1.653 ) {
      det = et * (sqrt(0.0283*0.0283 + (1.25/sqrt(et))*(1.25/sqrt(et)) + (4.79e-07/et)*(4.79e-07/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.566+1.653)/2))) * det;
    }
    
    if ( 1.653 <= fabs(eta) && fabs(eta) < 1.740 ) {
      det = et * (sqrt(0.0316*0.0316 + (1.21/sqrt(et))*(1.21/sqrt(et)) + (5e-05/et)*(5e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.653+1.740)/2))) * det;
    }
    
    if ( 1.740 <= fabs(eta) && fabs(eta) < 1.830 ) {
      det = et * (sqrt(2.29e-07*2.29e-07 + (1.2/sqrt(et))*(1.2/sqrt(et)) + (1.71e-05/et)*(1.71e-05/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.740+1.830)/2))) * det;
    }
    
    if ( 1.830 <= fabs(eta) && fabs(eta) < 1.930 ) {
      det = et * (sqrt(5.18e-09*5.18e-09 + (1.14/sqrt(et))*(1.14/sqrt(et)) + (1.7/et)*(1.7/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.830+1.930)/2))) * det;
    }
    
    if ( 1.930 <= fabs(eta) && fabs(eta) < 2.043 ) {
      det = et * (sqrt(2.17e-07*2.17e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (2.08/et)*(2.08/et)));
      //de = 1.0/sin(2 * atan(exp(-(1.930+2.043)/2))) * det;
    }
    
    if ( 2.043 <= fabs(eta) && fabs(eta) < 2.172 ) {
      det = et * (sqrt(3.65e-07*3.65e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.63/et)*(1.63/et)));
      //de = 1.0/sin(2 * atan(exp(-(2.043+2.172)/2))) * det;
    }
    
    if ( 2.172 <= fabs(eta) && fabs(eta) < 2.322 ) {
      det = et * (sqrt(2.02e-07*2.02e-07 + (1.09/sqrt(et))*(1.09/sqrt(et)) + (1.68/et)*(1.68/et)));
      //de = 1.0/sin(2 * atan(exp(-(2.172+2.322)/2))) * det;
    }
    
    if ( 2.322 <= fabs(eta) && fabs(eta) < 2.500 ) {
      det = et * (sqrt(5.27e-07*5.27e-07 + (1.12/sqrt(et))*(1.12/sqrt(et)) + (1.78/et)*(1.78/et)));
      //de = 1.0/sin(2 * atan(exp(-(2.322+2.500)/2))) * det;
    }
    
    //std::cout<< "--> returning det = "<< det << std::endl;
    return det;
  }
}

double BJetEtResolution::operator()(double genEt, double genEta, double recEt) const
{
  //std::cout << "<BJetEtResolution::operator>: genEt = " << genEt << ", genEta = " << genEta << ", recEt = " << recEt << std::endl;
  double sigmaEt = getBjetEtResolution(genEt, genEta);
  double pull = (recEt - genEt)/sigmaEt;
  const double sqrt2Pi = TMath::Sqrt(2.*TMath::Pi());
  double retVal = (1./(sigmaEt*sqrt2Pi))*TMath::Exp(-0.5*pull*pull);
  //std::cout<< "--> returning retVal = "<< retVal << std::endl;
  return retVal;
}
