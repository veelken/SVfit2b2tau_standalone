#include "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneAlgorithm2b2tau.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"

#include <TGraphErrors.h>

namespace svFitStandalone2b2tau
{
  TH1* makeHistogram(const std::string& histogramName, double xMin, double xMax, double logBinWidth)
  {
    if ( xMin <= 0 ) xMin = 0.1;
    int numBins = 1 + TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
    TArrayF binning(numBins + 1);
    binning[0] = 0.;
    double x = xMin;  
    for ( int iBin = 1; iBin <= numBins; ++iBin ) {
      binning[iBin] = x;
      x *= logBinWidth;
    }  
    TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), numBins, binning.GetArray());
    return histogram;
  }
  void compHistogramDensity(const TH1* histogram, TH1* histogram_density) 
  {
    for ( int iBin = 1; iBin <= histogram->GetNbinsX(); ++iBin ) {
      double binContent = histogram->GetBinContent(iBin);
      double binError = histogram->GetBinError(iBin);
      double binWidth = histogram->GetBinWidth(iBin);
      //if ( histogram == histogramMass_ ) {
      //  TAxis* xAxis = histogram->GetXaxis();
      //  std::cout << "bin #" << iBin << " (x = " << xAxis->GetBinLowEdge(iBin) << ".." << xAxis->GetBinUpEdge(iBin) << ", width = " << binWidth << "):"
      //	      << " " << binContent << " +/- " << binError << std::endl;
      //}
      assert(binWidth > 0.);
      histogram_density->SetBinContent(iBin, binContent/binWidth);
      histogram_density->SetBinError(iBin, binError/binWidth);
    }
  }
  double extractValue(const TH1* histogram, TH1* histogram_density) 
  {
    double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084, mean3sigmaWithinMax, mean5sigmaWithinMax;
    compHistogramDensity(histogram, histogram_density);
    svFitStandalone::extractHistogramProperties(histogram, histogram_density, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084, mean3sigmaWithinMax, mean5sigmaWithinMax);
    //double value = maximum_interpol;
    double value = maximum;
    return value;
  }
  double extractUncertainty(const TH1* histogram, TH1* histogram_density)
  {
    double maximum, maximum_interpol, mean, quantile016, quantile050, quantile084, mean3sigmaWithinMax, mean5sigmaWithinMax;
    compHistogramDensity(histogram, histogram_density);
    svFitStandalone::extractHistogramProperties(histogram, histogram_density, maximum, maximum_interpol, mean, quantile016, quantile050, quantile084, mean3sigmaWithinMax, mean5sigmaWithinMax);
    //double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum_interpol, 2.) + TMath::Power(maximum_interpol - quantile016, 2.)));
    double uncertainty = TMath::Sqrt(0.5*(TMath::Power(quantile084 - maximum, 2.) + TMath::Power(maximum - quantile016, 2.)));
    return uncertainty;  
  }
  double extractLmax(const TH1* histogram, TH1* histogram_density)
  {
    compHistogramDensity(histogram, histogram_density);
    double Lmax = histogram_density->GetMaximum();
    return Lmax;
  }

  void map_xMarkovChain2b2tau(const double* x, 
			      bool l1isLep, bool l2isLep, bool marginalizeVisMass, bool shiftVisMass, bool shiftVisPt, 
			      double diTauMass_unfitted, double mHiggsTauTau, 
			      double diBJetMass_unfitted, double bJet1Et_unfitted, double bJet2Et_unfitted, double mHiggsBB, 
			      double* x_mapped)
  {
    using namespace svFitStandalone;
    int offset1 = 0;
    double shiftVisPt1 = 1.;
    if ( l1isLep ) {
      x_mapped[kXFrac]                 = x[0];
      x_mapped[kMNuNu]                 = x[1];
      x_mapped[kPhi]                   = x[2];
      x_mapped[kVisMassShifted]        = 0.;
      x_mapped[kRecTauPtDivGenTauPt]   = 0.;
      offset1 = 3;
    } else {
      x_mapped[kXFrac]                 = x[0];
      x_mapped[kMNuNu]                 = 0.;
      x_mapped[kPhi]                   = x[1];
      offset1 = 2;
      if ( marginalizeVisMass || shiftVisMass ) {
	x_mapped[kVisMassShifted]      = x[offset1];
	++offset1;
      } else {
	x_mapped[kVisMassShifted]      = 0.;
      }
      if ( shiftVisPt ) {
	x_mapped[kRecTauPtDivGenTauPt] = x[offset1];
	shiftVisPt1 = x_mapped[kRecTauPtDivGenTauPt];
	++offset1;
      } else {
	x_mapped[kRecTauPtDivGenTauPt] = 0.;
      }
    }
    int offset2 = 0;
    double shiftVisPt2 = 1.;
    if ( l2isLep ) {
      if ( mHiggsTauTau <= 0. ) {
	x_mapped[kMaxFitParams + kXFrac]               = x[offset1 + offset2];
	++offset2;
      }
      x_mapped[kMaxFitParams + kMNuNu]                 = x[offset1 + offset2]; 
      ++offset2;
      x_mapped[kMaxFitParams + kPhi]                   = x[offset1 + offset2];
      ++offset2;
      x_mapped[kMaxFitParams + kVisMassShifted]        = 0.;
      x_mapped[kMaxFitParams + kRecTauPtDivGenTauPt]   = 0.;
    } else {
      if ( mHiggsTauTau <= 0. ) {
	x_mapped[kMaxFitParams + kXFrac]               = x[offset1 + offset2];
	++offset2;
      }
      x_mapped[kMaxFitParams + kMNuNu]                 = 0.;
      x_mapped[kMaxFitParams + kPhi]                   = x[offset1 + offset2];
      ++offset2;
      if ( marginalizeVisMass || shiftVisMass ) {
	x_mapped[kMaxFitParams + kVisMassShifted]      = x[offset1 + offset2];
	++offset2;
      } else {
	x_mapped[kMaxFitParams + kVisMassShifted]      = 0.;
      }
      if ( shiftVisPt ) {	
	x_mapped[kMaxFitParams + kRecTauPtDivGenTauPt] = x[offset1 + offset2];
	shiftVisPt2 = x_mapped[kMaxFitParams + kRecTauPtDivGenTauPt];
	++offset2;
      } else {
	x_mapped[kMaxFitParams + kRecTauPtDivGenTauPt] = 0.;
      }
    }
    if ( mHiggsTauTau > 0. ) {
      double x1 = x_mapped[kXFrac];
      if ( x1 > 0. ) {
	double diTauMass_shifted = TMath::Sqrt(shiftVisPt1*shiftVisPt2)*diTauMass_unfitted;
	double x2 = (1./x1)*square(diTauMass_shifted/mHiggsTauTau);
	x_mapped[kMaxFitParams + kXFrac] = x2;
      }
    }
    double bJet1Et_fitted = x[offset1 + offset2];
    x_mapped[2*kMaxFitParams] = bJet1Et_fitted;
    double bJet2Et_fitted;
    if ( mHiggsBB > 0. ) {
      bJet2Et_fitted = bJet2Et_unfitted*(bJet1Et_unfitted/TMath::Max(1., bJet1Et_fitted))*square(mHiggsBB/diBJetMass_unfitted);
    } else {
      bJet2Et_fitted = x[offset1 + offset2 + 1];
    }
    x_mapped[2*kMaxFitParams + 1] = bJet2Et_fitted;
    //std::cout << "<map_xMarkovChain2b2tau>:" << std::endl;
    //for ( int i = 0; i < 6; ++i ) {
    //  std::cout << " x_mapped[" << i << "] = " << x_mapped[i] << std::endl;
    //}
    //std::cout << " x_mapped[" << (2*kMaxFitParams) << "] = " << x_mapped[2*kMaxFitParams] << std::endl;
    //std::cout << " x_mapped[" << (2*kMaxFitParams + 1) << "] = " << x_mapped[2*kMaxFitParams + 1] << std::endl;    
  }
}

SVfitStandaloneAlgorithm2b2tau::SVfitStandaloneAlgorithm2b2tau(const std::vector<svFitStandalone::MeasuredTauLepton>& measuredTauLeptons, double measuredMETx, double measuredMETy, const TMatrixD& covMET, 
							       const std::vector<svFitStandalone2b2tau::MeasuredBJet>& measuredBJets, int verbose) 
  : fitStatus_(-1), 
    verbose_(verbose), 
    maxObjFunctionCalls_(10000),
    mHiggsTauTau_(-1.), // Higgs mass constraint disabled for Higgs that decays into two taus
    mHiggsBB_(125.),
    mcObjectiveFunctionAdapter_(0),
    mcPtEtaPhiMassAdapter_higgs2tau_(0),
    mcPtEtaPhiMassAdapter_higgs2b_(0),
    mcPtEtaPhiMassAdapter_diHiggs_(0),
    integrator2_(0),
    integrator2_nDim_(0),
    isInitialized2_(false),
    maxObjFunctionCalls2_(100000),
    marginalizeVisMass_(false),
    lutVisMassAllDMs_(0),
    shiftVisMass_(false),
    lutVisMassResDM0_(0),
    lutVisMassResDM1_(0),
    lutVisMassResDM10_(0),
    shiftVisPt_(false),
    lutVisPtResDM0_(0),
    lutVisPtResDM1_(0),
    lutVisPtResDM10_(0)
{ 
  //std::cout << "<SVfitStandaloneAlgorithm2b2tau::SVfitStandaloneAlgorithm2b2tau>:" << std::endl;
  //std::cout << " verbose = " << verbose << std::endl;

  // instantiate minuit, the arguments might turn into configurables once
  minimizer_ = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  // instantiate the combined likelihood
  svFitStandalone::Vector measuredMET_rounded(svFitStandalone::roundToNdigits(measuredMETx), svFitStandalone::roundToNdigits(measuredMETy), 0.);
  TMatrixD covMET_rounded(2, 2);
  covMET_rounded[0][0] = svFitStandalone::roundToNdigits(covMET[0][0]);
  covMET_rounded[1][0] = svFitStandalone::roundToNdigits(covMET[1][0]);
  covMET_rounded[0][1] = svFitStandalone::roundToNdigits(covMET[0][1]);
  covMET_rounded[1][1] = svFitStandalone::roundToNdigits(covMET[1][1]);
  nll_ = new svFitStandalone2b2tau::SVfitStandaloneLikelihood2b2tau(measuredTauLeptons, measuredMET_rounded, covMET_rounded, measuredBJets, (verbose_ > 2));
  nll_->setMassHiggsTauTau(mHiggsTauTau_);
  nll_->setMassHiggsBB(mHiggsBB_);
  nllStatus_ = nll_->error();
  
  this->reset();

  clock_ = new TBenchmark();
}

void SVfitStandaloneAlgorithm2b2tau::reset()
{
  diHiggsMass_ = 0.;
  diHiggsMassUncert_ = 0.;
  diHiggsMassLmax_ = 0.; 
  diHiggsPt_ = 0.; 
  diHiggsPtUncert_ = 0.; 
  diHiggsPtLmax_ = 0.; 
  diHiggsEta_ = 0.; 
  diHiggsEtaUncert_ = 0.; 
  diHiggsEtaLmax_ = 0.; 
  diHiggsPhi_ = 0.; 
  diHiggsPhiUncert_ = 0.; 
  diHiggsPhiLmax_ = 0.; 
  diHiggsTransverseMass_ = 0.; 
  diHiggsTransverseMassUncert_ = 0.; 
  diHiggsTransverseMassLmax_ = 0.; 

  higgs2tauMass_ = 0.; 
  higgs2tauMassUncert_ = 0.; 
  higgs2tauMassLmax_ = 0.; 
  higgs2tauPt_ = 0.; 
  higgs2tauPtUncert_ = 0.; 
  higgs2tauPtLmax_ = 0.; 
  higgs2tauEta_ = 0.; 
  higgs2tauEtaUncert_ = 0.; 
  higgs2tauEtaLmax_ = 0.; 
  higgs2tauPhi_ = 0.; 
  higgs2tauPhiUncert_ = 0.; 
  higgs2tauPhiLmax_ = 0.; 
  higgs2tauTransverseMass_ = 0.; 
  higgs2tauTransverseMassUncert_ = 0.; 
  higgs2tauTransverseMassLmax_ = 0.; 

  higgs2bMass_ = 0.; 
  higgs2bMassUncert_ = 0.; 
  higgs2bMassLmax_ = 0.; 
  higgs2bPt_ = 0.; 
  higgs2bPtUncert_ = 0.; 
  higgs2bPtLmax_ = 0.; 
  higgs2bEta_ = 0.; 
  higgs2bEtaUncert_ = 0.; 
  higgs2bEtaLmax_ = 0.; 
  higgs2bPhi_ = 0.; 
  higgs2bPhiUncert_ = 0.; 
  higgs2bPhiLmax_ = 0.; 
  higgs2bTransverseMass_ = 0.; 
  higgs2bTransverseMassUncert_ = 0.; 
  higgs2bTransverseMassLmax_ = 0.; 
}

SVfitStandaloneAlgorithm2b2tau::~SVfitStandaloneAlgorithm2b2tau() 
{
  delete nll_;
  delete minimizer_;
  delete mcObjectiveFunctionAdapter_;
  delete mcPtEtaPhiMassAdapter_higgs2tau_;
  delete mcPtEtaPhiMassAdapter_higgs2b_;
  delete mcPtEtaPhiMassAdapter_diHiggs_;
  delete integrator2_;
  //delete lutVisMassAllDMs_;
  //delete lutVisMassResDM0_;
  //delete lutVisMassResDM1_;
  //delete lutVisMassResDM10_;
  //delete lutVisPtResDM0_;
  //delete lutVisPtResDM1_;
  //delete lutVisPtResDM10_;

  delete clock_;
}

namespace
{
  TH1* readHistogram(TFile* inputFile, const std::string& histogramName)
  {
    TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
    if ( !histogram ) {
      std::cerr << "<readHistogram>: Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
      assert(0);
    }
    return histogram;
  }
}

void 
SVfitStandaloneAlgorithm2b2tau::marginalizeVisMass(bool value, TFile* inputFile)
{
  marginalizeVisMass_ = value;
  if ( marginalizeVisMass_ ) {
    delete lutVisMassAllDMs_;
    TH1* lutVisMassAllDMs_tmp = readHistogram(inputFile, "DQMData/genTauMassAnalyzer/genTauJetMass");
    if ( lutVisMassAllDMs_tmp->GetNbinsX() >= 1000 ) lutVisMassAllDMs_tmp->Rebin(100);
    lutVisMassAllDMs_ = lutVisMassAllDMs_tmp;
  }
}

void 
SVfitStandaloneAlgorithm2b2tau::marginalizeVisMass(bool value, const TH1* lut)
{
  marginalizeVisMass_ = value;
  if ( marginalizeVisMass_ ) {
    lutVisMassAllDMs_ = lut;
  }
}

void 
SVfitStandaloneAlgorithm2b2tau::shiftVisMass(bool value, TFile* inputFile)
{
  shiftVisMass_ = value;
  if ( shiftVisMass_ ) {
    delete lutVisMassResDM0_;
    lutVisMassResDM0_ = readHistogram(inputFile, "recMinusGenTauMass_recDecayModeEq0");
    delete lutVisMassResDM1_;
    lutVisMassResDM1_ = readHistogram(inputFile, "recMinusGenTauMass_recDecayModeEq1");
    delete lutVisMassResDM10_;
    lutVisMassResDM10_ = readHistogram(inputFile, "recMinusGenTauMass_recDecayModeEq10");
  }
}

void 
SVfitStandaloneAlgorithm2b2tau::shiftVisPt(bool value, TFile* inputFile)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    delete lutVisPtResDM0_;
    lutVisPtResDM0_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq0");
    delete lutVisPtResDM1_;
    lutVisPtResDM1_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq1");
    delete lutVisPtResDM10_;
    lutVisPtResDM10_ = readHistogram(inputFile, "recTauPtDivGenTauPt_recDecayModeEq10");
  }
}

void
SVfitStandaloneAlgorithm2b2tau::setup()
{
  using namespace svFitStandalone;

  //if ( verbose_ >= 1 ) {
  //  std::cout << "<SVfitStandaloneAlgorithm2b2tau::setup()>:" << std::endl;
  //}
  for ( size_t idx = 0; idx < nll_->measuredTauLeptons().size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = nll_->measuredTauLeptons()[idx];
    //if ( verbose_ >= 1 ) {
    //  std::cout << " --> upper limit of leg1::mNuNu will be set to "; 
    //  if ( measuredTauLepton.type() == kTauToHadDecay ) { 
    //	  std::cout << "0";
    //  } else {
    //	  std::cout << (svFitStandalone::tauLeptonMass - TMath::Min(measuredTauLepton.mass(), 1.5));
    //  } 
    //  std::cout << std::endl;
    //}
    // start values for xFrac
    if ( idx == 0 || mHiggsTauTau_ <= 0. ) {
      minimizer_->SetLimitedVariable(
        idx*kMaxFitParams + kXFrac, 
        std::string(TString::Format("leg%i::xFrac", (int)idx + 1)).c_str(), 
        0.5, 0.1, 0., 1.);
    } else {
      minimizer_->SetFixedVariable(
        idx*kMaxFitParams + kXFrac, 
        std::string(TString::Format("leg%i::xFrac", (int)idx + 1)).c_str(), 
        0.);				     
    }
    // start values for nunuMass (leptonic tau decays only)
    if ( measuredTauLepton.type() == kTauToHadDecay ) { 
      minimizer_->SetFixedVariable(
        idx*kMaxFitParams + kMNuNu, 
	std::string(TString::Format("leg%i::mNuNu", (int)idx + 1)).c_str(), 
	0.); 
    } else { 
      minimizer_->SetLimitedVariable(
        idx*kMaxFitParams + kMNuNu, 
	std::string(TString::Format("leg%i::mNuNu", (int)idx + 1)).c_str(), 
	0.8, 0.10, 0., svFitStandalone::tauLeptonMass - TMath::Min(measuredTauLepton.mass(), 1.5)); 
    }
    // start values for phi
    minimizer_->SetVariable(
      idx*kMaxFitParams + kPhi, 
      std::string(TString::Format("leg%i::phi", (int)idx + 1)).c_str(), 
      0.0, 0.25);
    // start values for Pt and mass of visible tau decay products (hadronic tau decays only)
    if ( measuredTauLepton.type() == kTauToHadDecay && (marginalizeVisMass_ || shiftVisMass_) ) {
      minimizer_->SetLimitedVariable(
        idx*kMaxFitParams + kVisMassShifted, 
        std::string(TString::Format("leg%i::mVisShift", (int)idx + 1)).c_str(), 
        0.8, 0.10, svFitStandalone::chargedPionMass, svFitStandalone::tauLeptonMass);
    } else {
      minimizer_->SetFixedVariable(
        idx*kMaxFitParams + kVisMassShifted, 
        std::string(TString::Format("leg%i::mVisShift", (int)idx + 1)).c_str(), 
        measuredTauLepton.mass());
    }
    if ( measuredTauLepton.type() == kTauToHadDecay && shiftVisPt_ ) {
      minimizer_->SetLimitedVariable(
        idx*kMaxFitParams + kRecTauPtDivGenTauPt, 
        std::string(TString::Format("leg%i::tauPtDivGenVisPt", (int)idx + 1)).c_str(), 
        0., 0.10, -1., +1.5);
    } else {     
      minimizer_->SetFixedVariable(
        idx*kMaxFitParams + kRecTauPtDivGenTauPt, 
        std::string(TString::Format("leg%i::tauPtDivGenVisPt", (int)idx + 1)).c_str(), 
        0.);
    }
  }
  double maxBJet1Et = 10.*nll_->measuredBJets()[0].Et();
  if ( maxBJet1Et > 5.e+3 ) maxBJet1Et = 5.e+3;
  minimizer_->SetLimitedVariable(
    2*kMaxFitParams, 
    "bJet1Et", 
    nll_->measuredBJets()[0].Et(), 1., 1., maxBJet1Et);
  if ( mHiggsBB_ <= 0. ) {
    double maxBJet2Et = 10.*nll_->measuredBJets()[1].Et();
    if ( maxBJet2Et > 5.e+3 ) maxBJet2Et = 5.e+3;
    minimizer_->SetLimitedVariable(
      2*kMaxFitParams + 1, 
      "bJet2Et", 
      nll_->measuredBJets()[1].Et(), 1., 1., maxBJet2Et);
  } else {
    minimizer_->SetFixedVariable(
      2*kMaxFitParams + 1, 
      "bJet2Et", 
      0.);
  }
}

void
SVfitStandaloneAlgorithm2b2tau::fit()
{
  //if ( verbose_ >= 1 ) {
  //  std::cout << "<SVfitStandaloneAlgorithm2b2tau::fit()>" << std::endl
  //	        << " dimension of fit    : " << nll_->measuredTauLeptons().size()*svFitStandalone::kMaxFitParams << std::endl
  //	        << " maxObjFunctionCalls : " << maxObjFunctionCalls_ << std::endl; 
  //}
  // clear minimizer
  minimizer_->Clear();
  // set verbosity level of minimizer
  minimizer_->SetPrintLevel(-1);
  // setup the function to be called and the dimension of the fit
  ROOT::Math::Functor toMinimize(standaloneObjectiveFunctionAdapterMINUIT_, nll_->measuredTauLeptons().size()*svFitStandalone::kMaxFitParams + 2);
  minimizer_->SetFunction(toMinimize); 
  setup();
  minimizer_->SetMaxFunctionCalls(maxObjFunctionCalls_);
  // set Minuit strategy = 2, in order to get reliable error estimates:
  // http://www-cdf.fnal.gov/physics/statistics/recommendations/minuit.html
  minimizer_->SetStrategy(2);
  // compute uncertainties for increase of objective function by 0.5 wrt. 
  // minimum (objective function is log-likelihood function)
  minimizer_->SetErrorDef(0.5);
  //if ( verbose_ >= 1 ) {
  //  std::cout << "starting ROOT::Math::Minimizer::Minimize..." << std::endl;
  //  std::cout << " #freeParameters = " << minimizer_->NFree() << ","
  //	        << " #constrainedParameters = " << (minimizer_->NDim() - minimizer_->NFree()) << std::endl;
  //}
  // do the minimization
  nll_->addDelta(false);
  nll_->addSinTheta(true);
  nll_->requirePhysicalSolution(false);
  minimizer_->Minimize();
  //if ( verbose_ >= 2 ) { 
  //  minimizer_->PrintResults(); 
  //};
  /* get Minimizer status code, check if solution is valid:
    
     0: Valid solution
     1: Covariance matrix was made positive definite
     2: Hesse matrix is invalid
     3: Estimated distance to minimum (EDM) is above maximum
     4: Reached maximum number of function calls before reaching convergence
     5: Any other failure
  */
  fitStatus_ = minimizer_->Status();
  //if ( verbose_ >=1 ) { 
  //  std::cout << "--> fitStatus = " << fitStatus_ << std::endl; 
  //}
  
  // and write out the result
  this->reset();
  using svFitStandalone::kXFrac;
  using svFitStandalone::kMNuNu;
  using svFitStandalone::kPhi;
  using svFitStandalone::kMaxFitParams;
  // update Higgs -> tautau system with final fit results
  nll_->results2tau(fittedTauLeptons_, minimizer_->X());
  // determine uncertainty of the fitted di-tau mass
  double x1RelErr = minimizer_->Errors()[kXFrac]/minimizer_->X()[kXFrac];
  double x2RelErr = minimizer_->Errors()[kMaxFitParams + kXFrac]/minimizer_->X()[kMaxFitParams + kXFrac];
  // this gives a unified treatment for retrieving the result for integration mode and fit mode
  fittedHiggs2tauSystem_ = fittedTauLeptons_[0] + fittedTauLeptons_[1];
  higgs2tauMass_ = fittedHiggs2tauSystem_.mass();
  higgs2tauMassUncert_ = TMath::Sqrt(0.25*x1RelErr*x1RelErr + 0.25*x2RelErr*x2RelErr)*fittedHiggs2tauSystem().mass();
  // update Higgs -> bb system with final fit results
  nll_->results2tau(fittedBJets_, minimizer_->X());
  fittedHiggs2bSystem_ = fittedBJets_[0] + fittedBJets_[1];
  higgs2bMass_ = fittedHiggs2bSystem_.mass();
  double bJet1EtErr = minimizer_->Errors()[2*kMaxFitParams];
  double bJet2EtErr = minimizer_->Errors()[2*kMaxFitParams + 1];
  using svFitStandalone::square;
  higgs2bMassUncert_ = TMath::Sqrt(0.25*square(bJet1EtErr/nll_->measuredBJets()[0].Et()) + 0.25*square(bJet2EtErr/nll_->measuredBJets()[1].Et()))*fittedHiggs2bSystem_.mass();
  // update di-Higgs system with final fit results
  fittedDiHiggsSystem_ = fittedHiggs2tauSystem_ + fittedHiggs2bSystem_;
  diHiggsMass_ = fittedDiHiggsSystem_.mass();
  diHiggsMassUncert_ = 0.; // CV: not implemented yet
}

void
SVfitStandaloneAlgorithm2b2tau::integrateMarkovChain(const std::string& likelihoodFileName)
{
  if ( verbose_ >= 1 ) {
    std::cout << "<SVfitStandaloneAlgorithm2b2tau::integrateMarkovChain>:" << std::endl;
    clock_->Start("<SVfitStandaloneAlgorithm2b2tau::integrateMarkovChain>");
  }
  if ( isInitialized2_ ) {
    mcPtEtaPhiMassAdapter_higgs2tau_->Reset();
    mcPtEtaPhiMassAdapter_higgs2b_->Reset();
    mcPtEtaPhiMassAdapter_diHiggs_->Reset();
  } else {
    // initialize    
    std::string initMode = "none";
    unsigned numIterBurnin = TMath::Nint(0.10*maxObjFunctionCalls2_);
    unsigned numIterSampling = maxObjFunctionCalls2_;
    unsigned numIterSimAnnealingPhase1 = TMath::Nint(0.02*maxObjFunctionCalls2_);
    unsigned numIterSimAnnealingPhase2 = TMath::Nint(0.06*maxObjFunctionCalls2_);
    double T0 = 15.;
    double alpha = 1.0 - 1.e+2/maxObjFunctionCalls2_;
    unsigned numChains = 7;
    unsigned numBatches = 1;
    unsigned L = 1;
    double epsilon0 = 1.e-2;
    double nu = 0.71;
    int verbose = ( verbose_ >= 3 ) ? verbose_ : -1;
    integrator2_ = new SVfitStandaloneMarkovChainIntegrator(
                         initMode, numIterBurnin, numIterSampling, numIterSimAnnealingPhase1, numIterSimAnnealingPhase2,
			 T0, alpha, numChains, numBatches, L, epsilon0, nu,
			 verbose);
    using namespace svFitStandalone2b2tau;
    mcObjectiveFunctionAdapter_ = new svFitStandalone2b2tau::MCObjectiveFunctionAdapter2b2tau();
    integrator2_->setIntegrand(*mcObjectiveFunctionAdapter_);
    integrator2_nDim_ = 0;
    mcPtEtaPhiMassAdapter_diHiggs_ = new svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau(MCPtEtaPhiMassAdapter2b2tau::kDiHiggs);
    integrator2_->registerCallBackFunction(*mcPtEtaPhiMassAdapter_diHiggs_);
    mcPtEtaPhiMassAdapter_higgs2tau_ = new svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau(MCPtEtaPhiMassAdapter2b2tau::kHiggs2tau);
    integrator2_->registerCallBackFunction(*mcPtEtaPhiMassAdapter_higgs2tau_);
    mcPtEtaPhiMassAdapter_higgs2b_ = new svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau(MCPtEtaPhiMassAdapter2b2tau::kHiggs2b);
    integrator2_->registerCallBackFunction(*mcPtEtaPhiMassAdapter_higgs2b_);
    isInitialized2_ = true;    
  }

  // number of parameters for fit
  int nDim = 0;
  l1isLep_ = false;
  l2isLep_ = false;
  const TH1* l1lutVisMass = 0;
  const TH1* l1lutVisMassRes = 0;
  const TH1* l1lutVisPtRes = 0;
  const TH1* l2lutVisMass = 0;
  const TH1* l2lutVisMassRes = 0;
  const TH1* l2lutVisPtRes = 0;
  for ( size_t idx = 0; idx < nll_->measuredTauLeptons().size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = nll_->measuredTauLeptons()[idx];
    if ( idx == 0 ) {
      idxFitParLeg1_ = 0;
      if ( measuredTauLepton.type() == svFitStandalone::kTauToHadDecay ) { 
	if ( marginalizeVisMass_ ) {
	  l1lutVisMass = lutVisMassAllDMs_;
	}
	if ( shiftVisMass_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    l1lutVisMassRes = lutVisMassResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    l1lutVisMassRes = lutVisMassResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    l1lutVisMassRes = lutVisMassResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisMass is enabled, but leg1 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
		      << " --> disabling shiftVisMass for this event !!" << std::endl;
	  }
        }
	if ( shiftVisPt_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    l1lutVisPtRes = lutVisPtResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    l1lutVisPtRes = lutVisPtResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    l1lutVisPtRes = lutVisPtResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisPt is enabled, but leg1 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
		      << " --> disabling shiftVisPt for this event !!" << std::endl;
	  }
        }
	nDim += 2;
	if ( l1lutVisMass || l1lutVisMassRes ) {
	  ++nDim;
	}
	if ( l1lutVisPtRes ) {
	  ++nDim;
	}
      } else {
	l1isLep_ = true;
	nDim += 3;
      }
    }
    if ( idx == 1 ) {
      idxFitParLeg2_ = nDim;
      if ( measuredTauLepton.type() == svFitStandalone::kTauToHadDecay ) { 
	if ( marginalizeVisMass_ ) {
	  l2lutVisMass = lutVisMassAllDMs_;
	}
	if ( shiftVisMass_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    l2lutVisMassRes = lutVisMassResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    l2lutVisMassRes = lutVisMassResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    l2lutVisMassRes = lutVisMassResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisMass is enabled, but leg2 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
		      << " --> disabling shiftVisMass for this event !!" << std::endl;
	  }
        }
	if ( shiftVisPt_ ) {
	  if ( measuredTauLepton.decayMode() == 0 ) {
	    l2lutVisPtRes = lutVisPtResDM0_;
	  } else if ( measuredTauLepton.decayMode() == 1 || measuredTauLepton.decayMode() == 2 ) {
	    l2lutVisPtRes = lutVisPtResDM1_;
	  } else if ( measuredTauLepton.decayMode() == 10 ) {
	    l2lutVisPtRes = lutVisPtResDM10_;
	  } else {
	    std::cerr << "Warning: shiftVisPt is enabled, but leg2 decay mode = " << measuredTauLepton.decayMode() << " is not supported" 
		      << " --> disabling shiftVisPt for this event !!" << std::endl;
	  }
        }
	if ( mHiggsTauTau_ <= 0. ) nDim += 2;
	else nDim += 1;
	if ( l2lutVisMass || l2lutVisMassRes ) {
	  ++nDim;
	}
	if ( l2lutVisPtRes ) {
	  ++nDim;
	}
      } else {
	l2isLep_ = true;
	if ( mHiggsTauTau_ <= 0. ) nDim += 3;
	else nDim += 2;
      }
    }
  }
  idxFitParBJet1Et_ = nDim;
  nDim += 1;
  if ( mHiggsBB_ <= 0. ) {
    idxFitParBJet2Et_ = nDim;
    nDim += 1;
  } else {
    idxFitParBJet2Et_ = -1;
  }

  if ( nDim != integrator2_nDim_ ) {
    mcObjectiveFunctionAdapter_->SetNDim(nDim);    
    integrator2_->setIntegrand(*mcObjectiveFunctionAdapter_);
    mcPtEtaPhiMassAdapter_diHiggs_->SetNDim(nDim);
    mcPtEtaPhiMassAdapter_higgs2tau_->SetNDim(nDim);
    mcPtEtaPhiMassAdapter_higgs2b_->SetNDim(nDim);
    integrator2_nDim_ = nDim;
  }
  mcObjectiveFunctionAdapter_->SetL1isLep(l1isLep_);
  mcObjectiveFunctionAdapter_->SetL2isLep(l2isLep_);
  if ( marginalizeVisMass_ && shiftVisMass_ ) {
    std::cerr << "Error: marginalizeVisMass and shiftVisMass flags must not both be enabled !!" << std::endl;
    assert(0);
  }  
  mcObjectiveFunctionAdapter_->SetMarginalizeVisMass(marginalizeVisMass_ && (l1lutVisMass || l2lutVisMass));
  mcObjectiveFunctionAdapter_->SetShiftVisMass(shiftVisMass_ && (l1lutVisMassRes || l2lutVisMassRes));
  mcObjectiveFunctionAdapter_->SetShiftVisPt(shiftVisPt_ && (l1lutVisPtRes || l2lutVisPtRes));
  mcObjectiveFunctionAdapter_->SetDiTauMass_unfitted(nll_->diTauMass_unfitted());
  mcObjectiveFunctionAdapter_->SetMassHiggsTauTau(mHiggsTauTau_);
  mcObjectiveFunctionAdapter_->SetDiBJetMass_unfitted(nll_->diBJetMass_unfitted());
  mcObjectiveFunctionAdapter_->SetBJet1Et_unfitted(nll_->bJet1Et_unfitted());
  mcObjectiveFunctionAdapter_->SetBJet2Et_unfitted(nll_->bJet2Et_unfitted());
  mcObjectiveFunctionAdapter_->SetMassHiggsBB(mHiggsBB_);
  // CV: use binning relative to the visible mass to avoid that SVfit mass is "quantizied" at the same discrete values for all events.
  TH1* histogramMass_diHiggs = svFitStandalone2b2tau::makeHistogram("SVfitStandaloneAlgorithmMarkovChain_histogramMass_diHiggs", measuredHiggs2tauSystem().mass()/1.0125, 1.e+4, 1.025);
  TH1* histogramMass_diHiggs_density = (TH1*)histogramMass_diHiggs->Clone(Form("%s_density", histogramMass_diHiggs->GetName()));
  mcPtEtaPhiMassAdapter_diHiggs_->SetHistogramMass(histogramMass_diHiggs, histogramMass_diHiggs_density);
  mcPtEtaPhiMassAdapter_diHiggs_->SetL1isLep(l1isLep_);
  mcPtEtaPhiMassAdapter_diHiggs_->SetL2isLep(l2isLep_);
  mcPtEtaPhiMassAdapter_diHiggs_->SetMarginalizeVisMass(marginalizeVisMass_ && (l1lutVisMass || l2lutVisMass));
  mcPtEtaPhiMassAdapter_diHiggs_->SetShiftVisMass(shiftVisMass_ && (l1lutVisMassRes || l2lutVisMassRes));
  mcPtEtaPhiMassAdapter_diHiggs_->SetShiftVisPt(shiftVisPt_ && (l1lutVisPtRes || l2lutVisPtRes));
  mcPtEtaPhiMassAdapter_diHiggs_->SetDiTauMass_unfitted(nll_->diTauMass_unfitted());
  mcPtEtaPhiMassAdapter_diHiggs_->SetMassHiggsTauTau(mHiggsTauTau_);
  mcPtEtaPhiMassAdapter_diHiggs_->SetDiBJetMass_unfitted(nll_->diBJetMass_unfitted());
  mcPtEtaPhiMassAdapter_diHiggs_->SetBJet1Et_unfitted(nll_->bJet1Et_unfitted());
  mcPtEtaPhiMassAdapter_diHiggs_->SetBJet2Et_unfitted(nll_->bJet2Et_unfitted());
  mcPtEtaPhiMassAdapter_diHiggs_->SetMassHiggsBB(mHiggsBB_);
  TH1* histogramMass_higgs2tau = svFitStandalone2b2tau::makeHistogram("SVfitStandaloneAlgorithmMarkovChain_histogramMass_higgs2tau", measuredHiggs2tauSystem().mass()/1.0125, 1.e+4, 1.025);
  TH1* histogramMass_higgs2tau_density = (TH1*)histogramMass_higgs2tau->Clone(Form("%s_density", histogramMass_higgs2tau->GetName()));
  mcPtEtaPhiMassAdapter_higgs2tau_->SetHistogramMass(histogramMass_higgs2tau, histogramMass_higgs2tau_density);
  mcPtEtaPhiMassAdapter_higgs2tau_->CopySettingsFrom(*mcPtEtaPhiMassAdapter_diHiggs_);
  TH1* histogramMass_higgs2b = svFitStandalone2b2tau::makeHistogram("SVfitStandaloneAlgorithmMarkovChain_histogramMass_higgs2b", measuredHiggs2bSystem().mass()*0.0125, 1.e+4, 1.025);
  TH1* histogramMass_higgs2b_density = (TH1*)histogramMass_higgs2b->Clone(Form("%s_density", histogramMass_higgs2b->GetName()));
  mcPtEtaPhiMassAdapter_higgs2b_->SetHistogramMass(histogramMass_higgs2b, histogramMass_higgs2b_density);
  mcPtEtaPhiMassAdapter_higgs2b_->CopySettingsFrom(*mcPtEtaPhiMassAdapter_diHiggs_);
  
  /* --------------------------------------------------------------------------------------
     lower and upper bounds for integration. Boundaries are defined for each decay channel
     separately. The order is: 
     
     - fully hadronic {xhad1, phihad1, (masshad1, pthad1), xhad2, phihad2, (masshad2, pthad2)}
     - semi  leptonic {xlep, nunuMass, philep, xhad, phihad, (masshad, pthad)}
     - fully leptonic {xlep1, nunuMass1, philep1, xlep2, nunuMass2, philep2}
     
     x0* defines the start value for the integration, xl* defines the lower integation bound, 
     xh* defines the upper integration bound in the following definitions. 
     ATTENTION: order matters here! In the semi-leptonic decay the lepton must go first in 
     the parametrization, as it is first in the definition of integral boundaries. This is
     the reason why the measuredLeptons are eventually re-ordered in the constructor of 
     this class before passing them on to SVfitStandaloneLikelihood.
  */
  std::vector<double> x0(nDim);
  std::vector<double> xl(nDim);
  std::vector<double> xh(nDim);
  if ( l1isLep_ ) {
    x0[idxFitParLeg1_ + 0] = 0.5; 
    xl[idxFitParLeg1_ + 0] = 0.0; 
    xh[idxFitParLeg1_ + 0] = 1.0; 
    x0[idxFitParLeg1_ + 1] = 0.8; 
    xl[idxFitParLeg1_ + 1] = 0.0; 
    xh[idxFitParLeg1_ + 1] = svFitStandalone::tauLeptonMass; 
    x0[idxFitParLeg1_ + 2] = 0.0; 
    xl[idxFitParLeg1_ + 2] = -TMath::Pi(); 
    xh[idxFitParLeg1_ + 2] = +TMath::Pi();
  } else {
    x0[idxFitParLeg1_ + 0] = 0.5; 
    xl[idxFitParLeg1_ + 0] = 0.0; 
    xh[idxFitParLeg1_ + 0] = 1.0; 
    x0[idxFitParLeg1_ + 1] = 0.0; 
    xl[idxFitParLeg1_ + 1] = -TMath::Pi(); 
    xh[idxFitParLeg1_ + 1] = +TMath::Pi();
    int offset1 = 2;
    if ( marginalizeVisMass_ || shiftVisMass_ ) {
      x0[idxFitParLeg1_ + offset1] = 0.8; 
      xl[idxFitParLeg1_ + offset1] = svFitStandalone::chargedPionMass; 
      xh[idxFitParLeg1_ + offset1] = svFitStandalone::tauLeptonMass; 
      ++offset1;
    }
    if ( shiftVisPt_ ) {
      x0[idxFitParLeg1_ + offset1] = 0.0; 
      xl[idxFitParLeg1_ + offset1] = -1.0; 
      xh[idxFitParLeg1_ + offset1] = +1.5;
      ++offset1;
    }
  }
  if ( l2isLep_ ) {
    int offset2 = 0;
    if ( mHiggsTauTau_ <= 0. ) {
      x0[idxFitParLeg2_ + offset2] = 0.5; 
      xl[idxFitParLeg2_ + offset2] = 0.0; 
      xh[idxFitParLeg2_ + offset2] = 1.0; 
      ++offset2;
    }
    x0[idxFitParLeg2_ + offset2] = 0.8; 
    xl[idxFitParLeg2_ + offset2] = 0.0; 
    xh[idxFitParLeg2_ + offset2] = svFitStandalone::tauLeptonMass; 
    ++offset2;
    x0[idxFitParLeg2_ + offset2] = 0.0; 
    xl[idxFitParLeg2_ + offset2] = -TMath::Pi(); 
    xh[idxFitParLeg2_ + offset2] = +TMath::Pi();
    ++offset2;
  } else {
    int offset2 = 0;
    if ( mHiggsTauTau_ <= 0. ) {
      x0[idxFitParLeg2_ + offset2] = 0.5; 
      xl[idxFitParLeg2_ + offset2] = 0.0; 
      xh[idxFitParLeg2_ + offset2] = 1.0; 
      ++offset2;
    }
    x0[idxFitParLeg2_ + offset2] = 0.0; 
    xl[idxFitParLeg2_ + offset2] = -TMath::Pi(); 
    xh[idxFitParLeg2_ + offset2] = +TMath::Pi();
    ++offset2;
    if ( marginalizeVisMass_ || shiftVisMass_ ) {
      x0[idxFitParLeg2_ + offset2] = 0.8; 
      xl[idxFitParLeg2_ + offset2] = svFitStandalone::chargedPionMass; 
      xh[idxFitParLeg2_ + offset2] = svFitStandalone::tauLeptonMass; 
      ++offset2;
    }
    if ( shiftVisPt_ ) {     
      x0[idxFitParLeg2_ + offset2] = 0.0; 
      xl[idxFitParLeg2_ + offset2] = -1.0; 
      xh[idxFitParLeg2_ + offset2] = +1.5;
      ++offset2;
    }
  }
  x0[idxFitParBJet1Et_] = nll_->bJet1Et_unfitted();
  xl[idxFitParBJet1Et_] = 1.;
  double maxBJet1Et = 10.*nll_->bJet1Et_unfitted();
  if ( maxBJet1Et > 5.e+3 ) maxBJet1Et = 5.e+3;
  if ( maxBJet1Et > 2.*nll_->bJet1Et_unfitted() ) maxBJet1Et = 2.*nll_->bJet1Et_unfitted(); // ONLY FOR TESTING !!!
  xh[idxFitParBJet1Et_] = maxBJet1Et;
  if ( mHiggsBB_ <= 0. ) {
    x0[idxFitParBJet2Et_] = nll_->bJet2Et_unfitted();
    xl[idxFitParBJet2Et_] = 1.;
    double maxBJet2Et = 10.*nll_->bJet2Et_unfitted();
    if ( maxBJet2Et > 5.e+3 ) maxBJet2Et = 5.e+3;
    if ( maxBJet2Et > 2.*nll_->bJet2Et_unfitted() ) maxBJet2Et = 2.*nll_->bJet2Et_unfitted(); // ONLY FOR TESTING !!!
    xh[idxFitParBJet2Et_] = maxBJet2Et;
  }

  //std::cout << "nDim = " << nDim << std::endl;
  //std::cout << " idxFitParLeg1 = " << idxFitParLeg1_ << std::endl;
  //std::cout << " idxFitParLeg2 = " << idxFitParLeg2_ << std::endl;
  //std::cout << " idxFitParBJet1Et = " << idxFitParBJet1Et_ << std::endl;
  //std::cout << " idxFitParBJet2Et = " << idxFitParBJet2Et_ << std::endl;
  for ( int i = 0; i < nDim; ++i ) {
    // transform startPosition into interval ]0..1[
    // expected by MarkovChainIntegrator class
    //std::cout << "before transformation: x0[" << i << "] = " << x0[i] << " (xl[" << i << "] = " << xl[i] << ", xh[" << i << "] = " << xh[i] << ")" << std::endl;
    x0[i] = (x0[i] - xl[i])/(xh[i] - xl[i]);
    //std::cout << "after transformation: x0[" << i << "] = " << x0[i] << std::endl;
  }
  integrator2_->initializeStartPosition_and_Momentum(x0);
  nll_->addDelta(true);
  nll_->addSinTheta(false);
  nll_->addPhiPenalty(false);
  nll_->marginalizeVisMass(marginalizeVisMass_ && (l1lutVisMass || l2lutVisMass), l1lutVisMass, l2lutVisMass);
  nll_->shiftVisMass(shiftVisMass_ && (l1lutVisMassRes || l2lutVisMassRes), l1lutVisMassRes, l2lutVisMassRes);
  nll_->shiftVisPt(shiftVisPt_ && (l1lutVisPtRes && l2lutVisPtRes), l1lutVisPtRes, l2lutVisPtRes);
  nll_->requirePhysicalSolution(true);

  double integral = 0.;
  double integralErr = 0.;
  int errorFlag = 0;
  integrator2_->integrate(xl, xh, integral, integralErr, errorFlag);
  fitStatus_ = errorFlag;
  diHiggsPt_ = mcPtEtaPhiMassAdapter_diHiggs_->getPt();
  diHiggsPtUncert_ = mcPtEtaPhiMassAdapter_diHiggs_->getPtUncert();
  diHiggsPtLmax_ = mcPtEtaPhiMassAdapter_diHiggs_->getPtLmax();
  diHiggsEta_ = mcPtEtaPhiMassAdapter_diHiggs_->getEta();
  diHiggsEtaUncert_ = mcPtEtaPhiMassAdapter_diHiggs_->getEtaUncert();
  diHiggsEtaLmax_ = mcPtEtaPhiMassAdapter_diHiggs_->getEtaLmax();
  diHiggsPhi_ = mcPtEtaPhiMassAdapter_diHiggs_->getPhi();
  diHiggsPhiUncert_ = mcPtEtaPhiMassAdapter_diHiggs_->getPhiUncert();
  diHiggsPhiLmax_ = mcPtEtaPhiMassAdapter_diHiggs_->getPhiLmax();
  diHiggsMass_ = mcPtEtaPhiMassAdapter_diHiggs_->getMass();
  diHiggsMassUncert_ = mcPtEtaPhiMassAdapter_diHiggs_->getMassUncert();
  diHiggsMassLmax_ = mcPtEtaPhiMassAdapter_diHiggs_->getMassLmax();
  diHiggsTransverseMass_ = mcPtEtaPhiMassAdapter_diHiggs_->getTransverseMass();
  diHiggsTransverseMassUncert_ = mcPtEtaPhiMassAdapter_diHiggs_->getTransverseMassUncert();
  diHiggsTransverseMassLmax_ = mcPtEtaPhiMassAdapter_diHiggs_->getTransverseMassLmax();
  fittedDiHiggsSystem_ = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(diHiggsPt_, diHiggsEta_, diHiggsPhi_, diHiggsMass_);
  higgs2tauPt_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPt();
  higgs2tauPtUncert_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPtUncert();
  higgs2tauPtLmax_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPtLmax();
  higgs2tauEta_ = mcPtEtaPhiMassAdapter_higgs2tau_->getEta();
  higgs2tauEtaUncert_ = mcPtEtaPhiMassAdapter_higgs2tau_->getEtaUncert();
  higgs2tauEtaLmax_ = mcPtEtaPhiMassAdapter_higgs2tau_->getEtaLmax();
  higgs2tauPhi_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPhi();
  higgs2tauPhiUncert_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPhiUncert();
  higgs2tauPhiLmax_ = mcPtEtaPhiMassAdapter_higgs2tau_->getPhiLmax();
  higgs2tauMass_ = mcPtEtaPhiMassAdapter_higgs2tau_->getMass();
  higgs2tauMassUncert_ = mcPtEtaPhiMassAdapter_higgs2tau_->getMassUncert();
  higgs2tauMassLmax_ = mcPtEtaPhiMassAdapter_higgs2tau_->getMassLmax();
  higgs2tauTransverseMass_ = mcPtEtaPhiMassAdapter_higgs2tau_->getTransverseMass();
  higgs2tauTransverseMassUncert_ = mcPtEtaPhiMassAdapter_higgs2tau_->getTransverseMassUncert();
  higgs2tauTransverseMassLmax_ = mcPtEtaPhiMassAdapter_higgs2tau_->getTransverseMassLmax();
  fittedHiggs2tauSystem_ = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(higgs2tauPt_, higgs2tauEta_, higgs2tauPhi_, higgs2tauMass_);
  higgs2bPt_ = mcPtEtaPhiMassAdapter_higgs2b_->getPt();
  higgs2bPtUncert_ = mcPtEtaPhiMassAdapter_higgs2b_->getPtUncert();
  higgs2bPtLmax_ = mcPtEtaPhiMassAdapter_higgs2b_->getPtLmax();
  higgs2bEta_ = mcPtEtaPhiMassAdapter_higgs2b_->getEta();
  higgs2bEtaUncert_ = mcPtEtaPhiMassAdapter_higgs2b_->getEtaUncert();
  higgs2bEtaLmax_ = mcPtEtaPhiMassAdapter_higgs2b_->getEtaLmax();
  higgs2bPhi_ = mcPtEtaPhiMassAdapter_higgs2b_->getPhi();
  higgs2bPhiUncert_ = mcPtEtaPhiMassAdapter_higgs2b_->getPhiUncert();
  higgs2bPhiLmax_ = mcPtEtaPhiMassAdapter_higgs2b_->getPhiLmax();
  higgs2bMass_ = mcPtEtaPhiMassAdapter_higgs2b_->getMass();
  higgs2bMassUncert_ = mcPtEtaPhiMassAdapter_higgs2b_->getMassUncert();
  higgs2bMassLmax_ = mcPtEtaPhiMassAdapter_higgs2b_->getMassLmax();
  higgs2bTransverseMass_ = mcPtEtaPhiMassAdapter_higgs2b_->getTransverseMass();
  higgs2bTransverseMassUncert_ = mcPtEtaPhiMassAdapter_higgs2b_->getTransverseMassUncert();
  higgs2bTransverseMassLmax_ = mcPtEtaPhiMassAdapter_higgs2b_->getTransverseMassLmax();
  fittedHiggs2bSystem_ = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >(higgs2bPt_, higgs2bEta_, higgs2bPhi_, higgs2bMass_);

  if ( likelihoodFileName != "" ) {
    TFile* likelihoodFile = new TFile(likelihoodFileName.data(), "RECREATE");
    histogramMass_diHiggs->Write();
    histogramMass_diHiggs_density->Write();
    histogramMass_higgs2tau->Write();
    histogramMass_higgs2tau_density->Write();
    histogramMass_higgs2b->Write();
    histogramMass_higgs2b_density->Write();
    delete likelihoodFile;
  }

  delete histogramMass_diHiggs;
  delete histogramMass_diHiggs_density;
  mcPtEtaPhiMassAdapter_diHiggs_->SetHistogramMass(0, 0);
  delete histogramMass_higgs2tau;
  delete histogramMass_higgs2tau_density;
  mcPtEtaPhiMassAdapter_higgs2tau_->SetHistogramMass(0, 0);
  delete histogramMass_higgs2b;
  delete histogramMass_higgs2b_density;
  mcPtEtaPhiMassAdapter_higgs2b_->SetHistogramMass(0, 0);

  if ( verbose_ >= 1 ) {
    std::cout << "--> di-Higgs: Pt = " << diHiggsPt_ << ", eta = " << diHiggsEta_ << ", phi = " << diHiggsPhi_ << ", mass  = " << diHiggsMass_ 
	      << " (mtt = " << higgs2tauMass_ << ", mbb = " << higgs2bMass_ << ")" << std::endl;
    clock_->Show("<SVfitStandaloneAlgorithm2b2tau::integrateMarkovChain>");
  }
}
