#ifndef TauAnalysis_SVfitStandalone2b2tau_BJetEtResolution_h
#define TauAnalysis_SVfitStandalone2b2tau_BJetEtResolution_h

namespace svFitStandalone2b2tau
{
  class BJetEtResolution
  {
   public:
    BJetEtResolution();
    ~BJetEtResolution();
    
    double operator()(double, double, double) const;
  };
}

#endif
