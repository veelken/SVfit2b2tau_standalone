#ifndef TauAnalysis_SVfitStandalone_SVfitStandaloneAlgorithm2b2tau_h
#define TauAnalysis_SVfitStandalone_SVfitStandaloneAlgorithm2b2tau_h

#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneMarkovChainIntegrator.h"
#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
#include "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneLikelihood2b2tau.h"

#include <TMath.h>
#include <TArrayF.h>
#include <TString.h>
#include <TH1.h>
#include <TBenchmark.h>

using svFitStandalone::Vector;
using svFitStandalone::LorentzVector;
using svFitStandalone::MeasuredTauLepton;
using svFitStandalone2b2tau::MeasuredBJet;

/**
   \class   ObjectFunctionAdapter SVfitStandaloneAlgorithm2b2tau.h "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneAlgorithm2b2tau.h"
   
   \brief   Function interface to minuit.
   
   This class is an interface, which is used as global function pointer of the combined likelihood as defined in src/SVfitStandaloneLikelihood2b2tau.cc
   to VEGAS or minuit. It is a member of the of the SVfitStandaloneAlgorithm class defined below and is used in SVfitStandalone::fit(), or 
   SVfitStandalone::integrate(), where it is passed on to a ROOT::Math::Functor. The parameters x correspond to the array of fit/integration 
   paramters as defined in interface/SVfitStandaloneLikelihood.h of this package. In the fit mode these are made known to minuit in the function
   SVfitStandaloneAlgorithm::setup. In the integration mode the mapping is done internally in the SVfitStandaloneLikelihood::tansformint. This
   has to be in sync. with the definition of the integration boundaries in SVfitStandaloneAlgorithm::integrate. 
*/

namespace svFitStandalone2b2tau
{
  TH1* makeHistogram(const std::string&, double, double, double);
  void compHistogramDensity(const TH1*, TH1*);
  double extractValue(const TH1*, TH1*);
  double extractUncertainty(const TH1*, TH1*);
  double extractLmax(const TH1*, TH1*);

  // for "fit" (MINUIT) mode
  class ObjectiveFunctionAdapter2b2tauMINUIT
  {
  public:
    double operator()(const double* x) const // NOTE: return value = -log(likelihood)
    {
      double prob = SVfitStandaloneLikelihood2b2tau::gSVfitStandaloneLikelihood2b2tau->prob(x);
      double nll;
      if ( prob > 0. ) nll = -TMath::Log(prob);
      else nll = std::numeric_limits<float>::max();
      return nll;
    }
  };
  // for Markov Chain integration
  void map_xMarkovChain2b2tau(const double*, bool, bool, bool, bool, bool, double, double, double, double, double, double, double*);
  class MCObjectiveFunctionAdapter2b2tau : public ROOT::Math::Functor
  {
   public:
    void SetL1isLep(bool l1isLep) { l1isLep_ = l1isLep; }
    void SetL2isLep(bool l2isLep) { l2isLep_ = l2isLep; }
    void SetMarginalizeVisMass(bool marginalizeVisMass) { marginalizeVisMass_ = marginalizeVisMass; }
    void SetShiftVisMass(bool shiftVisMass) { shiftVisMass_ = shiftVisMass; }
    void SetShiftVisPt(bool shiftVisPt) { shiftVisPt_ = shiftVisPt; }
    void SetNDim(int nDim) { nDim_ = nDim; }
    unsigned int NDim() const { return nDim_; }
    void SetDiTauMass_unfitted(double diTauMass) { diTauMass_unfitted_ = diTauMass; }
    void SetMassHiggsTauTau(double mH) { mHiggsTauTau_ = mH; }
    void SetDiBJetMass_unfitted(double diBJetMass) { diBJetMass_unfitted_ = diBJetMass; }
    void SetBJet1Et_unfitted(double bJet1Et) { bJet1Et_unfitted_ = bJet1Et; }
    void SetBJet2Et_unfitted(double bJet2Et) { bJet2Et_unfitted_ = bJet2Et; }
    void SetMassHiggsBB(double mH) { mHiggsBB_ = mH; }
   private:
    virtual double DoEval(const double* x) const
    {
      //std::cout << "<MCObjectiveFunctionAdapter2b2tau::DoEval(const double*)>:" << std::endl;
      map_xMarkovChain2b2tau(
        x, 
	l1isLep_, l2isLep_, marginalizeVisMass_, shiftVisMass_, shiftVisPt_, 
	diTauMass_unfitted_, mHiggsTauTau_, 
	diBJetMass_unfitted_, bJet1Et_unfitted_, bJet2Et_unfitted_, mHiggsBB_, 
	x_mapped_);
      double prob = SVfitStandaloneLikelihood2b2tau::gSVfitStandaloneLikelihood2b2tau->prob(x_mapped_);
      if ( TMath::IsNaN(prob) ) prob = 0.;
      return prob;
    } 
    mutable double x_mapped_[12];
    int nDim_;
    bool l1isLep_;
    bool l2isLep_;
    bool marginalizeVisMass_;
    bool shiftVisMass_;
    bool shiftVisPt_;
    double diTauMass_unfitted_;
    double mHiggsTauTau_;
    double diBJetMass_unfitted_;
    double bJet1Et_unfitted_;
    double bJet2Et_unfitted_;
    double mHiggsBB_;
  };
  class MCPtEtaPhiMassAdapter2b2tau : public ROOT::Math::Functor
  {
   public:
    enum { kDiHiggs, kHiggs2tau, kHiggs2b };
    MCPtEtaPhiMassAdapter2b2tau(int mode) 
      : mode_(mode)
    {
      // mode parameter:
      //   0 = di-Higgs system
      //   1 = Higgs -> tautau system
      //   2 = Higgs -> bb system
      std::string label;
      if      ( mode == kDiHiggs   ) label = "diHiggs";
      else if ( mode == kHiggs2tau ) label = "higgs2tau";
      else if ( mode == kHiggs2b   ) label = "higgs2b";
      else {
	std::cerr << "Error: mode parameter = " << mode << " invalid !!" << std::endl;	
	assert(0);
      }
      std::string histogramNamePt = Form("SVfitStandaloneAlgorithm2b2tau_histogramPt_%s", label.data());
      histogramPt_ = makeHistogram(histogramNamePt.data(), 1., 1.e+3, 1.025);
      histogramPt_density_ = (TH1*)histogramPt_->Clone(Form("%s_density", histogramPt_->GetName()));
      std::string histogramNameEta = Form("SVfitStandaloneAlgorithm2b2tau_histogramEta_%s", label.data());
      histogramEta_ = new TH1D(histogramNameEta.data(), histogramNameEta.data(), 198, -9.9, +9.9);
      histogramEta_density_ = (TH1*)histogramEta_->Clone(Form("%s_density", histogramEta_->GetName()));
      std::string histogramNamePhi = Form("SVfitStandaloneAlgorithm2b2tau_histogramPhi_%s", label.data());
      histogramPhi_ = new TH1D(histogramNamePhi.data(), histogramNamePhi.data(), 180, -TMath::Pi(), +TMath::Pi());
      histogramPhi_density_ = (TH1*)histogramPhi_->Clone(Form("%s_density", histogramPhi_->GetName()));
      std::string histogramNameMass = Form("SVfitStandaloneAlgorithm2b2tau_histogramMass_%s", label.data());
      histogramMass_ = makeHistogram(histogramNameMass.data(), 1.e+1, 1.e+4, 1.025);
      histogramMass_density_ = (TH1*)histogramMass_->Clone(Form("%s_density", histogramMass_->GetName()));
      std::string histogramNameTransverseMass = Form("SVfitStandaloneAlgorithm2b2tau_histogramTransverseMass_%s", label.data());
      histogramTransverseMass_ = makeHistogram(histogramNameTransverseMass.data(), 1., 1.e+4, 1.025);
      histogramTransverseMass_density_ = (TH1*)histogramTransverseMass_->Clone(Form("%s_density", histogramTransverseMass_->GetName()));
    }      
    ~MCPtEtaPhiMassAdapter2b2tau()
    {
      delete histogramPt_;
      delete histogramPt_density_;
      delete histogramEta_;
      delete histogramEta_density_;
      delete histogramPhi_;
      delete histogramPhi_density_;
      delete histogramMass_;
      delete histogramMass_density_;
      delete histogramTransverseMass_;
      delete histogramTransverseMass_density_;
    }
    void SetHistogramMass(TH1* histogramMass, TH1* histogramMass_density)
    {
      // CV: passing null pointers to the SetHistogramMass function
      //     indicates that the histograms have been deleted by the calling code
      if ( histogramMass != 0 ) delete histogramMass_;
      histogramMass_ = histogramMass;
      if ( histogramMass_density != 0 ) delete histogramMass_density_;
      histogramMass_density_ = histogramMass_density;
    }
    void SetL1isLep(bool l1isLep) { l1isLep_ = l1isLep; }
    void SetL2isLep(bool l2isLep) { l2isLep_ = l2isLep; }
    void SetMarginalizeVisMass(bool marginalizeVisMass) { marginalizeVisMass_ = marginalizeVisMass; }
    void SetShiftVisMass(bool shiftVisMass) { shiftVisMass_ = shiftVisMass; }
    void SetShiftVisPt(bool shiftVisPt) { shiftVisPt_ = shiftVisPt; }
    void SetNDim(int nDim) { nDim_ = nDim; }
    unsigned int NDim() const { return nDim_; }
    void SetDiTauMass_unfitted(double diTauMass) { diTauMass_unfitted_ = diTauMass; }
    void SetMassHiggsTauTau(double mH) { mHiggsTauTau_ = mH; }
    void SetDiBJetMass_unfitted(double diBJetMass) { diBJetMass_unfitted_ = diBJetMass; }
    void SetBJet1Et_unfitted(double bJet1Et) { bJet1Et_unfitted_ = bJet1Et; }
    void SetBJet2Et_unfitted(double bJet2Et) { bJet2Et_unfitted_ = bJet2Et; }
    void SetMassHiggsBB(double mH) { mHiggsBB_ = mH; }
    void CopySettingsFrom(const MCPtEtaPhiMassAdapter2b2tau& mcPtEtaPhiMassAdapter)
    {
      l1isLep_ = mcPtEtaPhiMassAdapter.l1isLep_;
      l2isLep_ = mcPtEtaPhiMassAdapter.l2isLep_;
      marginalizeVisMass_ = mcPtEtaPhiMassAdapter.marginalizeVisMass_;
      shiftVisMass_ = mcPtEtaPhiMassAdapter.shiftVisMass_;
      shiftVisPt_ = mcPtEtaPhiMassAdapter.shiftVisPt_;
      nDim_ = mcPtEtaPhiMassAdapter.nDim_;
      diTauMass_unfitted_ = mcPtEtaPhiMassAdapter.diTauMass_unfitted_;
      mHiggsTauTau_ = mcPtEtaPhiMassAdapter.mHiggsTauTau_;
      diBJetMass_unfitted_ = mcPtEtaPhiMassAdapter.diBJetMass_unfitted_;
      bJet1Et_unfitted_ = mcPtEtaPhiMassAdapter.bJet1Et_unfitted_;
      bJet2Et_unfitted_ = mcPtEtaPhiMassAdapter.bJet2Et_unfitted_;
      mHiggsBB_ = mcPtEtaPhiMassAdapter.mHiggsBB_;
    }
    void Reset()
    {
      histogramPt_->Reset();
      histogramEta_->Reset();
      histogramPhi_->Reset();
      histogramMass_->Reset();
      histogramTransverseMass_->Reset();
    }
    double getPt() const { return extractValue(histogramPt_, histogramPt_density_); }
    double getPtUncert() const { return extractUncertainty(histogramPt_, histogramPt_density_); }
    double getPtLmax() const { return extractLmax(histogramPt_, histogramPt_density_); }
    double getEta() const { return extractValue(histogramEta_, histogramEta_density_); }
    double getEtaUncert() const { return extractUncertainty(histogramEta_, histogramEta_density_); }
    double getEtaLmax() const { return extractLmax(histogramEta_, histogramEta_density_); }
    double getPhi() const { return extractValue(histogramPhi_, histogramPhi_density_); }
    double getPhiUncert() const { return extractUncertainty(histogramPhi_, histogramPhi_density_); }
    double getPhiLmax() const { return extractLmax(histogramPhi_, histogramPhi_density_); }
    double getMass() const { return extractValue(histogramMass_, histogramMass_density_); }
    double getMassUncert() const { return extractUncertainty(histogramMass_, histogramMass_density_); }
    double getMassLmax() const { return extractLmax(histogramMass_, histogramMass_density_); }
    double getTransverseMass() const { return extractValue(histogramTransverseMass_, histogramTransverseMass_density_); }
    double getTransverseMassUncert() const { return extractUncertainty(histogramTransverseMass_, histogramTransverseMass_density_); }
    double getTransverseMassLmax() const { return extractLmax(histogramTransverseMass_, histogramTransverseMass_density_); }
   private:    
    virtual double DoEval(const double* x) const
    {
      map_xMarkovChain2b2tau(
        x, 
	l1isLep_, l2isLep_, marginalizeVisMass_, shiftVisMass_, shiftVisPt_, 
	diTauMass_unfitted_, mHiggsTauTau_, 
	diBJetMass_unfitted_, bJet1Et_unfitted_, bJet2Et_unfitted_, mHiggsBB_, 
	x_mapped_);
      SVfitStandaloneLikelihood2b2tau::gSVfitStandaloneLikelihood2b2tau->results2tau(fittedTauLeptons_, x_mapped_);
      fittedHiggs2tauSystem_ = fittedTauLeptons_[0] + fittedTauLeptons_[1];
      SVfitStandaloneLikelihood2b2tau::gSVfitStandaloneLikelihood2b2tau->results2b(fittedBJets_, x_mapped_);
      fittedHiggs2bSystem_ = fittedBJets_[0] + fittedBJets_[1];
      fittedDiHiggsSystem_ = fittedHiggs2tauSystem_ + fittedHiggs2bSystem_;
      double pt = 0.;
      double eta = 0.;
      double phi = 0.;
      double mass = 0.;
      double mT = 0.;
      if ( mode_ == kDiHiggs ) {
	pt = fittedDiHiggsSystem_.pt();
	eta = fittedDiHiggsSystem_.eta();
	phi = fittedDiHiggsSystem_.phi();
	mass = fittedDiHiggsSystem_.mass();
	mT = TMath::Sqrt(2.*fittedHiggs2tauSystem_.pt()*fittedHiggs2bSystem_.pt()*(1. - TMath::Cos(fittedHiggs2tauSystem_.phi() - fittedHiggs2bSystem_.phi())));
      } else if ( mode_ == kHiggs2tau ) {
	pt = fittedHiggs2tauSystem_.pt();
	eta = fittedHiggs2tauSystem_.eta();
	phi = fittedHiggs2tauSystem_.phi();
	mass = fittedHiggs2tauSystem_.mass();
	mT = TMath::Sqrt(2.*fittedTauLeptons_[0].pt()*fittedTauLeptons_[1].pt()*(1. - TMath::Cos(fittedTauLeptons_[0].phi() - fittedTauLeptons_[1].phi())));
      } else if ( mode_ == kHiggs2b ) {
	pt = fittedHiggs2bSystem_.pt();
	eta = fittedHiggs2bSystem_.eta();
	phi = fittedHiggs2bSystem_.phi();
	mass = fittedHiggs2bSystem_.mass();
	//std::cout << "mbb = " << mass << std::endl;
	//if ( mHiggsBB_ > 100. && mass < 100. ) {
	//  std::cerr << "Warning in <MCPtEtaPhiMassAdapter2b2tau::DoEval(const double*)>:" << std::endl;
	//  std::cerr << " fittedBJet #0: Pt = " << fittedBJets_[0].pt() << ", eta = " << fittedBJets_[0].eta() << ", phi = " << fittedBJets_[0].phi() << ", mass = " << fittedBJets_[0].mass() << std::endl;
	//  std::cerr << " fittedBJet #1: Pt = " << fittedBJets_[1].pt() << ", eta = " << fittedBJets_[1].eta() << ", phi = " << fittedBJets_[1].phi() << ", mass = " << fittedBJets_[1].mass() << std::endl;
	//  std::cerr << "--> mbb = " << mass << ", expected 125 GeV !!" << std::endl;
	//  //assert(0);
	//}
	mT = TMath::Sqrt(2.*fittedBJets_[0].Et()*fittedBJets_[1].Et()*(1. - TMath::Cos(fittedBJets_[0].phi() - fittedBJets_[1].phi())));
      } else assert(0);
      histogramPt_->Fill(pt);
      histogramEta_->Fill(eta);
      histogramPhi_->Fill(phi);
      histogramMass_->Fill(mass);
      histogramTransverseMass_->Fill(mT);
      return 0.;
    } 
   protected:
 //public:
    mutable std::vector<svFitStandalone::LorentzVector> fittedTauLeptons_;
    mutable LorentzVector fittedHiggs2tauSystem_;
    mutable std::vector<svFitStandalone::LorentzVector> fittedBJets_;
    mutable LorentzVector fittedHiggs2bSystem_;
    mutable LorentzVector fittedDiHiggsSystem_;
    mutable TH1* histogramPt_;
    mutable TH1* histogramPt_density_;
    mutable TH1* histogramEta_;
    mutable TH1* histogramEta_density_;
    mutable TH1* histogramPhi_;
    mutable TH1* histogramPhi_density_;
    mutable TH1* histogramMass_;
    mutable TH1* histogramMass_density_;
    mutable TH1* histogramTransverseMass_;
    mutable TH1* histogramTransverseMass_density_;
    mutable double x_mapped_[12];
    int nDim_;
    bool l1isLep_;
    bool l2isLep_;
    bool marginalizeVisMass_;
    bool shiftVisMass_;
    bool shiftVisPt_;
    int mode_;
    double diTauMass_unfitted_;
    double mHiggsTauTau_;
    double diBJetMass_unfitted_;
    double bJet1Et_unfitted_;
    double bJet2Et_unfitted_;
    double mHiggsBB_;
  };
}

/**
   \class   SVfitStandaloneAlgorithm2b2tau SVfitStandaloneAlgorithm2b2tau.h "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneAlgorithm2b2tau.h"
   
   \brief   Standalone version of the SVfitAlgorithm2b2tau.

   This class is a standalone version of the SVfitAlgorithm to perform the full reconstruction of a di-Higgs boson to 2b + 2tau decay. The 
   implementation is supposed to deal with any combination of leptonic or hadronic tau decays. It exploits likelihood functions 
   as defined in interface/LikelihoodFunctions.h of this package, which are combined into a single likelihood function as defined 
   interface/SVfitStandaloneLikelihood2b2tau.h in this package. The combined likelihood function depends on the following variables: 

   \var nunuMass   : the invariant mass of the neutrino system for each decay branch (two parameters)
   \var decayAngle : the decay angle in the restframe of each decay branch (two parameters)
   \var visMass    : the mass of the visible component of the di-tau system (two parameters)

   The actual fit parameters are:

   \var nunuMass   : the invariant mass of the neutrino system for each decay branch (two parameters)
   \var xFrac      : the fraction of the visible energy on the energy of the tau lepton in the labframe (two parameters)
   \var phi        : the azimuthal angle of the tau lepton (two parameters)

   \var bJet1Et    : the transverse energy of the first b-jet (the Et of the second b-jet is determined by the Higgs mass contraint)

   In the fit mode. The azimuthal angle of each tau lepton is not constraint by measurement. It is limited to the physical values 
   from -Math::Phi to Math::Phi in the likelihood function of the combined likelihood class. The parameter nunuMass is constraint 
   to the tau lepton mass minus the mass of the visible part of the decay (which is itself constraint to values below the tau 
   lepton mass) in the setup function of this class. The parameter xFrac is constraint to values between 0. and 1. in the setup 
   function of this class. The invariant mass of the neutrino system is fixed to be zero for hadronic tau lepton decays as only 
   one (tau-) neutrino is involved in the decay. The original number of free parameters of 6 is therefore reduced by one for each 
   hadronic tau decay within the resonance. All information about the negative log likelihood is stored in the SVfitStandaloneLikelihood 
   class as defined in the same package. This class interfaces the combined likelihood to the ROOT::Math::Minuit minimization program. 
   It does setup/initialize the fit parameters as defined in interface/SVfitStandaloneLikelihood.h in this package, initializes the 
   minimization procedure, executes the fit algorithm and returns the fit result. The fit result consists of the fully reconstructed 
   di-tau system, from which also the invariant mass can be derived.

   In the integration mode xFrac for the second leptons is determiend from xFrac of the first lepton for given di-tau mass, thus reducing 
   the number of parameters to be integrated out wrt. to the fit version by one. The di-tau mass is scanned for the highest likelihood 
   starting from the visible mass of the two leptons. The return value is just the di-tau mass. 

   Common usage is: 
   
   // construct the class object from the minimal necessary information
   SVfitStandaloneAlgorithm2b2tau algo(measuredTauLeptons, measuredMET, covMET, measuredBJets);
   // apply customized configurations if wanted (examples are given below)
   //algo.maxObjFunctionCalls(10000); // only applies for fit mode
   //algo.addLogM(false);             // applies for fit and integration mode
   //algo.metPower(1.);               // applies for fit and integration mode
   // run the fit in fit mode
   algo.fit();
   // retrieve the results upon success
   if ( algo.isValidSolution() ) {
     std::cout << algo.diHiggsMass();
   }
   // run the integration in integration mode
   algo.integrate();
   std::cout << algo.diHiggsMass();

   The following optional parameters can be applied after initialization but before running the fit in fit mode: 

   \var metPower : indicating an additional power to enhance the MET likelihood (default is 1.)
   \var addLogM : specifying whether to use the LogM penalty term or not (default is true)     
   \var maxObjFunctionCalls : the maximum of function calls before the minimization procedure is terminated (default is 5000)
*/

class SVfitStandaloneAlgorithm2b2tau
{
 public:
  /// constructor from a minimal set of configurables
  SVfitStandaloneAlgorithm2b2tau(const std::vector<svFitStandalone::MeasuredTauLepton>&, double, double, const TMatrixD&, const std::vector<svFitStandalone2b2tau::MeasuredBJet>&, int verbose = 0);
  /// destructor
  ~SVfitStandaloneAlgorithm2b2tau();

  /// add an additional logM(tau,tau) term to the nll to suppress tails on M(tau,tau) (default is false)
  void addLogM(bool value, double power = 1.) { nll_->addLogM(value, power); }
  /// modify the MET term in the nll by an additional power (default is 1.)
  void metPower(double value) { nll_->metPower(value); }
  /// marginalize unknown mass of hadronic tau decay products (ATLAS case)
  void marginalizeVisMass(bool value, TFile* inputFile);
  void marginalizeVisMass(bool value, const TH1*);    
  /// take resolution on energy and mass of hadronic tau decays into account
  void shiftVisMass(bool value, TFile* inputFile);
  void shiftVisPt(bool value, TFile* inputFile);
  /// set mass of Higgs -> bb system (used as constraint for visible energy fraction X of second tau);
  /// set to -1 in order to disable this constraint
  void setMassHiggsTauTau(double mH) { mHiggsTauTau_ = mH; nll_->setMassHiggsTauTau(mHiggsTauTau_); }
  /// set mass of Higgs -> bb system (used as constraint for Et of second b-jet);
  /// set to -1 in order to disable this constraint
  void setMassHiggsBB(double mH) { mHiggsBB_ = mH; nll_->setMassHiggsBB(mHiggsBB_); }
  /// maximum function calls after which to stop the minimization procedure (default is 5000)
  void maxObjFunctionCalls(double value) { maxObjFunctionCalls_ = value; }

  /// fit to be called from outside
  void fit();
  /// integration by Markov Chain ('integrate' function kept for backwards-compatibility)
  void integrate() { return integrateMarkovChain(""); }
  /// integration by Markov Chain 
  void integrateMarkovChain(const std::string& likelihoodFileName = "");

  /// return status of minuit fit
  /*    
      0: Valid solution
      1: Covariance matrix was made positive definite
      2: Hesse matrix is invalid
      3: Estimated distance to minimum (EDM) is above maximum
      4: Reached maximum number of function calls before reaching convergence
      5: Any other failure
  */
  int fitStatus() { return fitStatus_; }
  /// return whether this is a valid solution or not
  bool isValidSolution() { return (nllStatus_ == 0 && fitStatus_ <= 0); }
  /// return whether this is a valid solution or not
  bool isValidFit() { return fitStatus_ == 0; }
  /// return whether this is a valid solution or not
  bool isValidNLL() { return nllStatus_ == 0; }

  /// return reconstructed mass, pT, eta and phi of the di-Higgs system, their corresponding uncertainties and maxima of the likelihood scan
  double diHiggsMass() const { return diHiggsMass_; }
  double diHiggsMassUncert() const { return diHiggsMassUncert_; }
  double diHiggsMassLmax() const { return diHiggsMassLmax_; }
  double diHiggsPt() const { return diHiggsPt_; }
  double diHiggsPtUncert() const { return diHiggsPtUncert_; }
  double diHiggsPtLmax() const { return diHiggsPtLmax_; }
  double diHiggsEta() const { return diHiggsEta_; }
  double diHiggsEtaUncert() const { return diHiggsEtaUncert_; }
  double diHiggsEtaLmax() const { return diHiggsEtaLmax_; }
  double diHiggsPhi() const { return diHiggsPhi_; }
  double diHiggsPhiUncert() const { return diHiggsPhiUncert_; }
  double diHiggsPhiLmax() const { return diHiggsPhiLmax_; }
  /// return reconstructed transverse mass of the di-Higgs system, the corresponding uncertainty and the maximum of the likelihood scan
  double diHiggsTransverseMass() const { return diHiggsTransverseMass_; }
  double diHiggsTransverseMassUncert() const { return diHiggsTransverseMassUncert_; }
  double diHiggsTransverseMassLmax() const { return diHiggsTransverseMassLmax_; }

  /// return reconstructed mass, pT, eta and phi of the Higgs boson that decays into taus, their corresponding uncertainties and maxima of the likelihood scan
  double higgs2tauMass() const { return higgs2tauMass_; }
  double higgs2tauMassUncert() const { return higgs2tauMassUncert_; }
  double higgs2tauMassLmax() const { return higgs2tauMassLmax_; }
  double higgs2tauPt() const { return higgs2tauPt_; }
  double higgs2tauPtUncert() const { return higgs2tauPtUncert_; }
  double higgs2tauPtLmax() const { return higgs2tauPtLmax_; }
  double higgs2tauEta() const { return higgs2tauEta_; }
  double higgs2tauEtaUncert() const { return higgs2tauEtaUncert_; }
  double higgs2tauEtaLmax() const { return higgs2tauEtaLmax_; }
  double higgs2tauPhi() const { return higgs2tauPhi_; }
  double higgs2tauPhiUncert() const { return higgs2tauPhiUncert_; }
  double higgs2tauPhiLmax() const { return higgs2tauPhiLmax_; }
  /// return reconstructed transverse mass of the Higgs boson that decays into taus, the corresponding uncertainty and the maximum of the likelihood scan
  double higgs2tauTransverseMass() const { return higgs2tauTransverseMass_; }
  double higgs2tauTransverseMassUncert() const { return higgs2tauTransverseMassUncert_; }
  double higgs2tauTransverseMassLmax() const { return higgs2tauTransverseMassLmax_; }

  /// return reconstructed mass, pT, eta and phi of the Higgs boson that decays into b-quarks, their corresponding uncertainties and maxima of the likelihood scan
  double higgs2bMass() const { return higgs2bMass_; }
  double higgs2bMassUncert() const { return higgs2bMassUncert_; }
  double higgs2bMassLmax() const { return higgs2bMassLmax_; }
  double higgs2bPt() const { return higgs2bPt_; }
  double higgs2bPtUncert() const { return higgs2bPtUncert_; }
  double higgs2bPtLmax() const { return higgs2bPtLmax_; }
  double higgs2bEta() const { return higgs2bEta_; }
  double higgs2bEtaUncert() const { return higgs2bEtaUncert_; }
  double higgs2bEtaLmax() const { return higgs2bEtaLmax_; }
  double higgs2bPhi() const { return higgs2bPhi_; }
  double higgs2bPhiUncert() const { return higgs2bPhiUncert_; }
  double higgs2bPhiLmax() const { return higgs2bPhiLmax_; }
  /// return reconstructed transverse mass of the Higgs boson that decays into b-quarks, the corresponding uncertainty and the maximum of the likelihood scan
  double higgs2bTransverseMass() const { return higgs2bTransverseMass_; }
  double higgs2bTransverseMassUncert() const { return higgs2bTransverseMassUncert_; }
  double higgs2bTransverseMassLmax() const { return higgs2bTransverseMassLmax_; }

  /// return 4-vector of the fitted di-Higgs system
  LorentzVector fittedDiHiggsSystem() const { return fittedDiHiggsSystem_; }
  /// return 4-vector of the measured di-Higgs system
  LorentzVector measuredDiHiggsSystem() const { return measuredTauLeptons()[0] + measuredTauLeptons()[1] + measuredBJets()[0] + measuredBJets()[1]; }
  /// return 4-vectors of the fitted tau leptons
  std::vector<LorentzVector> fittedTauLeptons() const { return fittedTauLeptons_; }
  /// return 4-vectors of measured visible decay products of the two taus
  std::vector<LorentzVector> measuredTauLeptons() const; 
  /// return 4-vector of the fitted Higgs -> tautau system
  LorentzVector fittedHiggs2tauSystem() const { return fittedHiggs2tauSystem_; }
  /// return 4-vector of the measured tau1 + tau2 system
  LorentzVector measuredHiggs2tauSystem() const { return measuredTauLeptons()[0] + measuredTauLeptons()[1]; }
  /// return 4-vectors of the fitted tau leptons
  std::vector<LorentzVector> fittedBJets() const { return fittedBJets_; }
  /// return 4-vectors of measured b-jets
  std::vector<LorentzVector> measuredBJets() const; 
  /// return 4-vector of the fitted Higgs -> bb system
  LorentzVector fittedHiggs2bSystem() const { return fittedHiggs2bSystem_; }
  /// return 4-vector of the measured tau1 + tau2 system
  LorentzVector measuredHiggs2bSystem() const { return measuredBJets()[0] + measuredBJets()[1]; }
  /// return spatial vector of the fitted MET
  Vector fittedMET() const { return (fittedHiggs2tauSystem().Vect() - measuredHiggs2tauSystem().Vect()); }
  // return spatial vector of the measured MET
  Vector measuredMET() const { return nll_->measuredMET(); }

 protected:
  /// setup the starting values for the minimization (default values for the fit parameters are taken from src/SVFitParameters.cc in the same package)
  void setup();

  /// reset all pT, eta, phi and mass values
  void reset();

 protected:
  /// return whether this is a valid solution or not
  int fitStatus_;
  /// return whether this is a valid solution or not
  unsigned int nllStatus_;
  /// verbosity level
  unsigned int verbose_;
  /// stop minimization after a maximal number of function calls
  unsigned int maxObjFunctionCalls_;

  /// mass of the Higgs boson that decays into two taus (used as constraint for visible energy fraction X of second tau, unless mHiggsTauTau_ == -1)
  double mHiggsTauTau_;
  /// mass of the Higgs boson that decays into two b-jets (used as constraint for Et of second b-jet, unless mHiggsBB_ == -1)
  double mHiggsBB_;

  /// minuit instance 
  ROOT::Math::Minimizer* minimizer_;
  /// standalone combined likelihood
  svFitStandalone2b2tau::SVfitStandaloneLikelihood2b2tau* nll_;
  /// needed to make the fit function callable from within minuit
  svFitStandalone2b2tau::ObjectiveFunctionAdapter2b2tauMINUIT standaloneObjectiveFunctionAdapterMINUIT_;

  /// reconstructed mass, pT, eta and phi of the di-Higgs system, their corresponding uncertainties and maxima of the likelihood scan
  double diHiggsMass_;
  double diHiggsMassUncert_;
  double diHiggsMassLmax_; 
  double diHiggsPt_; 
  double diHiggsPtUncert_; 
  double diHiggsPtLmax_; 
  double diHiggsEta_; 
  double diHiggsEtaUncert_; 
  double diHiggsEtaLmax_; 
  double diHiggsPhi_; 
  double diHiggsPhiUncert_; 
  double diHiggsPhiLmax_; 
  /// reconstructed transverse mass of the di-Higgs system, the corresponding uncertainty and the maximum of the likelihood scan
  double diHiggsTransverseMass_; 
  double diHiggsTransverseMassUncert_; 
  double diHiggsTransverseMassLmax_; 

  /// reconstructed mass, pT, eta and phi of the Higgs boson that decays into taus, their corresponding uncertainties and maxima of the likelihood scan
  double higgs2tauMass_; 
  double higgs2tauMassUncert_; 
  double higgs2tauMassLmax_; 
  double higgs2tauPt_; 
  double higgs2tauPtUncert_; 
  double higgs2tauPtLmax_; 
  double higgs2tauEta_; 
  double higgs2tauEtaUncert_; 
  double higgs2tauEtaLmax_; 
  double higgs2tauPhi_; 
  double higgs2tauPhiUncert_; 
  double higgs2tauPhiLmax_; 
  /// reconstructed transverse mass of the Higgs boson that decays into taus, the corresponding uncertainty and the maximum of the likelihood scan
  double higgs2tauTransverseMass_; 
  double higgs2tauTransverseMassUncert_; 
  double higgs2tauTransverseMassLmax_; 

  /// reconstructed mass, pT, eta and phi of the Higgs boson that decays into b-quarks, their corresponding uncertainties and maxima of the likelihood scan
  double higgs2bMass_; 
  double higgs2bMassUncert_; 
  double higgs2bMassLmax_; 
  double higgs2bPt_; 
  double higgs2bPtUncert_; 
  double higgs2bPtLmax_; 
  double higgs2bEta_; 
  double higgs2bEtaUncert_; 
  double higgs2bEtaLmax_; 
  double higgs2bPhi_; 
  double higgs2bPhiUncert_; 
  double higgs2bPhiLmax_; 
  /// reconstructed transverse mass of the Higgs boson that decays into b-quarks, the corresponding uncertainty and the maximum of the likelihood scan
  double higgs2bTransverseMass_; 
  double higgs2bTransverseMassUncert_; 
  double higgs2bTransverseMassLmax_; 

  /// fit result for each tau lepton
  std::vector<svFitStandalone::LorentzVector> fittedTauLeptons_;
  /// fit result for each b-jet
  std::vector<svFitStandalone::LorentzVector> fittedBJets_;
  /// fitted di-Higgs system
  svFitStandalone::LorentzVector fittedDiHiggsSystem_;
  /// fitted Higgs --> tautau system
  svFitStandalone::LorentzVector fittedHiggs2tauSystem_;
  /// fitted Higgs --> bb system
  svFitStandalone::LorentzVector fittedHiggs2bSystem_;

  /// needed for markov chain integration
  svFitStandalone2b2tau::MCObjectiveFunctionAdapter2b2tau* mcObjectiveFunctionAdapter_;
  svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau* mcPtEtaPhiMassAdapter_higgs2tau_;
  svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau* mcPtEtaPhiMassAdapter_higgs2b_;
  svFitStandalone2b2tau::MCPtEtaPhiMassAdapter2b2tau* mcPtEtaPhiMassAdapter_diHiggs_;
  SVfitStandaloneMarkovChainIntegrator* integrator2_;
  int integrator2_nDim_;
  bool isInitialized2_;
  unsigned maxObjFunctionCalls2_;

  TBenchmark* clock_;

  /// resolution on Pt and mass of hadronic taus
  bool marginalizeVisMass_;
  const TH1* lutVisMassAllDMs_;
  bool shiftVisMass_;
  const TH1* lutVisMassResDM0_;
  const TH1* lutVisMassResDM1_;
  const TH1* lutVisMassResDM10_;
  bool shiftVisPt_;  
  const TH1* lutVisPtResDM0_;
  const TH1* lutVisPtResDM1_;
  const TH1* lutVisPtResDM10_;

  bool l1isLep_;
  int idxFitParLeg1_;
  bool l2isLep_;
  int idxFitParLeg2_;
  int idxFitParBJet1Et_;
  int idxFitParBJet2Et_;
};

inline
std::vector<svFitStandalone::LorentzVector> 
SVfitStandaloneAlgorithm2b2tau::measuredTauLeptons() const 
{ 
  std::vector<svFitStandalone::LorentzVector> measuredTauLeptons;
  measuredTauLeptons.push_back(nll_->measuredTauLeptons()[0].p4());
  measuredTauLeptons.push_back(nll_->measuredTauLeptons()[1].p4());
  return measuredTauLeptons; 
}

inline
std::vector<svFitStandalone::LorentzVector> 
SVfitStandaloneAlgorithm2b2tau::measuredBJets() const 
{ 
  std::vector<svFitStandalone::LorentzVector> measuredBJets;
  measuredBJets.push_back(nll_->measuredBJets()[0].p4());
  measuredBJets.push_back(nll_->measuredBJets()[1].p4());
  return measuredBJets; 
}

#endif
