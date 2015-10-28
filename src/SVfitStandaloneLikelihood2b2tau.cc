#include "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneLikelihood2b2tau.h"

#include "TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
#include "TauAnalysis/SVfitStandalone/interface/LikelihoodFunctions.h"

using namespace svFitStandalone2b2tau;
using namespace svFitStandalone;

/// global function pointer for minuit or VEGAS
const SVfitStandaloneLikelihood2b2tau* SVfitStandaloneLikelihood2b2tau::gSVfitStandaloneLikelihood2b2tau = 0;
/// indicate first iteration for integration or fit cycle for debugging
static bool FIRST = true;

SVfitStandaloneLikelihood2b2tau::SVfitStandaloneLikelihood2b2tau(const std::vector<MeasuredTauLepton>& measuredTauLeptons, const Vector& measuredMET, const TMatrixD& covMET, 
								 const std::vector<MeasuredBJet>& measuredBJets, bool verbose) 
  : metPower_(1.0), 
    addLogM_(false), 
    powerLogM_(1.),
    addDelta_(true),
    addSinTheta_(false),
    addPhiPenalty_(true),
    verbose_(verbose), 
    idxObjFunctionCall_(0), 
    invCovMET_(2,2),
    mHiggsTauTau_(-1.), // Higgs mass constraint disabled for Higgs that decays into two taus
    mHiggsBB_(125.),
    errorCode_(0),
    requirePhysicalSolution_(false),
    marginalizeVisMass_(false),
    l1lutVisMass_(0),
    l2lutVisMass_(0),
    shiftVisMass_(false),
    l1lutVisMassRes_(0),
    l2lutVisMassRes_(0),
    shiftVisPt_(false),
    l1lutVisPtRes_(0),
    l2lutVisPtRes_(0)
{
  //if ( verbose_ ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::SVfitStandaloneLikelihood2b2tau>:" << std::endl;
  //}
  measuredMET_ = measuredMET;
  // for integration mode the order of lepton or tau matters due to the choice of order in which 
  // way the integration boundaries are defined. In this case the lepton should always go before
  // the tau. For tau-tau or lep-lep the order is irrelevant.
  if ( measuredTauLeptons[0].type() == svFitStandalone::kTauToHadDecay ) {
    measuredTauLeptons_.push_back(measuredTauLeptons[1]);
    measuredTauLeptons_.push_back(measuredTauLeptons[0]);
  } else {
    measuredTauLeptons_= measuredTauLeptons;
  }
  if ( measuredTauLeptons_.size() == 2 ) {
    diTauMass_unfitted_ = (measuredTauLeptons[0].p4() + measuredTauLeptons[1].p4()).mass();
  } else {
    std::cout << " >> ERROR : the number of measured leptons must be 2 but is found to be: " << measuredTauLeptons_.size() << std::endl;
    errorCode_ |= LeptonNumber;
  }
  // determine transfer matrix for MET
  invCovMET_= covMET;
  covDet_ = invCovMET_.Determinant();
  if ( covDet_ != 0 ) { 
    invCovMET_.Invert(); 
  } else{
    std::cout << " >> ERROR: cannot invert MET covariance Matrix (det=0)." << std::endl;
    errorCode_ |= MatrixInversion;
  }
  measuredBJets_ = measuredBJets;
  if ( measuredBJets_.size() == 2 ) {
    diBJetMass_unfitted_ = (measuredBJets[0].p4() + measuredBJets[1].p4()).mass();
    bJet1Et_unfitted_ = measuredBJets[0].Et();
    bJet1Eta_ = measuredBJets[0].eta();
    bJet2Et_unfitted_ = measuredBJets[1].Et();
    bJet2Eta_ = measuredBJets[1].eta();
  } else {
    std::cout << " >> ERROR : the number of measured b-jets must be 2 but is found to be: " << measuredBJets_.size() << std::endl;
    errorCode_ |= BJetNumber;
  }
  // set global function pointer to this
  gSVfitStandaloneLikelihood2b2tau = this;
}

void 
SVfitStandaloneLikelihood2b2tau::marginalizeVisMass(bool value, const TH1* l1lutVisMass, const TH1* l2lutVisMass)
{
  marginalizeVisMass_ = value;
  if ( marginalizeVisMass_ ) {
    l1lutVisMass_ = l1lutVisMass;
    l2lutVisMass_ = l2lutVisMass;
  }
}

void 
SVfitStandaloneLikelihood2b2tau::shiftVisMass(bool value, const TH1* l1lutVisMassRes, const TH1* l2lutVisMassRes)
{
  shiftVisMass_ = value;
  if ( shiftVisMass_ ) {
    l1lutVisMassRes_ = l1lutVisMassRes;
    l2lutVisMassRes_ = l2lutVisMassRes;
  }
}

void 
SVfitStandaloneLikelihood2b2tau::shiftVisPt(bool value, const TH1* l1lutVisPtRes, const TH1* l2lutVisPtRes)
{
  shiftVisPt_ = value;
  if ( shiftVisPt_ ) {
    l1lutVisPtRes_ = l1lutVisPtRes;
    l2lutVisPtRes_ = l2lutVisPtRes;
  }
}

const double*
SVfitStandaloneLikelihood2b2tau::transform(double* xPrime, const double* x) const
{
  //if ( verbose_ ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::transform(double*, const double*)>:" << std::endl;
  //}
  LorentzVector fittedHiggs2tauSystem;
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];

    // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
    // for each tau decay
    double nunuMass, labframeXFrac, labframePhi;
    double visMass_unshifted = measuredTauLepton.mass();
    double visMass = visMass_unshifted; // visible momentum in lab-frame
    double labframeVisMom_unshifted = measuredTauLepton.momentum(); 
    double labframeVisMom = labframeVisMom_unshifted; // visible momentum in lab-frame
    double labframeVisEn  = measuredTauLepton.energy(); // visible energy in lab-frame    
    if ( measuredTauLepton.type() == kTauToElecDecay || measuredTauLepton.type() == kTauToMuDecay ) {
      labframeXFrac = x[idx*kMaxFitParams + kXFrac];
      nunuMass = x[idx*kMaxFitParams + kMNuNu];
      labframePhi = x[idx*kMaxFitParams + kPhi];
    } else {
      labframeXFrac = x[idx*kMaxFitParams + kXFrac];
      nunuMass = 0.;
      labframePhi = x[idx*kMaxFitParams + kPhi];
      if ( marginalizeVisMass_ || shiftVisMass_ ) {
	visMass = x[idx*kMaxFitParams + kVisMassShifted];
      } 
      if ( shiftVisPt_ ) {
	double shiftInv = 1. + x[idx*kMaxFitParams + kRecTauPtDivGenTauPt];
	double shift = ( shiftInv > 1.e-1 ) ?
	  (1./shiftInv) : 1.e+1;
	labframeVisMom *= shift;
	//visMass *= shift; // CV: take mass and momentum to be correlated
	//labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
	labframeVisEn *= shift;
      }
    }
    bool isValidSolution = true;
    // add protection against unphysical mass of visible tau decay products
    if ( visMass < electronMass || visMass > tauLeptonMass ) { 
      std::cout << "--> setting isValidSolution = false, because visMass = " << visMass << std::endl;
      isValidSolution = false;
    }  
    // add protection against unphysical visible energy fractions
    if ( !(labframeXFrac >= 0. && labframeXFrac <= 1.) ) {
      std::cout << "--> setting isValidSolution = false, because labframeXFrac = " << labframeXFrac << std::endl;
      isValidSolution = false;
    }
    // CV: do not spend time on unphysical solutions: returning 0 pointer will lead to 0 evaluation of prob
    if ( !isValidSolution ) {
      return 0;
    }
    double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
    double enTau_lab = labframeVisEn/labframeXFrac;
    if ( (enTau_lab*enTau_lab) < tauLeptonMass2 ) {
      std::cout << "--> setting isValidSolution = false, because enTau_lab = " << enTau_lab << std::endl;
      enTau_lab = tauLeptonMass;
      isValidSolution = false;
    }  
    double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
    double gamma = enTau_lab/tauLeptonMass;
    double beta = TMath::Sqrt(1. - 1./(gamma*gamma));
    double pVis_parl_rf = -beta*gamma*labframeVisEn + gamma*TMath::Cos(gjAngle_lab)*labframeVisMom;
    double pVis_perp = labframeVisMom*TMath::Sin(gjAngle_lab);
    double gjAngle_rf = TMath::ATan2(pVis_perp, pVis_parl_rf);
    Vector p3Tau_unit = motherDirection(measuredTauLepton.direction(), gjAngle_lab, labframePhi);
    LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
    fittedHiggs2tauSystem += p4Tau_lab;
    // fill branch-wise nll parameters
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kNuNuMass1            : svFitStandalone2b2tau::kNuNuMass2            ] = nunuMass;
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kVisMass1             : svFitStandalone2b2tau::kVisMass2             ] = visMass;
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kDecayAngle1          : svFitStandalone2b2tau::kDecayAngle2          ] = gjAngle_rf;
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kDeltaVisMass1        : svFitStandalone2b2tau::kDeltaVisMass2        ] = ( marginalizeVisMass_ ) ? visMass : (visMass_unshifted - visMass);
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kRecTauPtDivGenTauPt1 : svFitStandalone2b2tau::kRecTauPtDivGenTauPt2 ] = ( labframeVisMom > 0. ) ? (labframeVisMom_unshifted/labframeVisMom) : 1.e+3;
    xPrime[ idx == 0 ? svFitStandalone2b2tau::kMaxNLLParams         : (svFitStandalone2b2tau::kMaxNLLParams + 1)   ] = labframeXFrac;
    xPrime[ idx == 0 ? (svFitStandalone2b2tau::kMaxNLLParams + 2)   : (svFitStandalone2b2tau::kMaxNLLParams + 3)   ] = isValidSolution;
  }
 
  Vector fittedMET = fittedHiggs2tauSystem.Vect() - (measuredTauLeptons_[0].p() + measuredTauLeptons_[1].p()); 
  // fill event-wise nll parameters
  xPrime[ svFitStandalone2b2tau::kDMETx ] = measuredMET_.x() - fittedMET.x(); 
  xPrime[ svFitStandalone2b2tau::kDMETy ] = measuredMET_.y() - fittedMET.y();
  if ( mHiggsTauTau_ > 0. ) xPrime[ svFitStandalone2b2tau::kMTauTau ] = mHiggsTauTau_;
  else xPrime[ svFitStandalone2b2tau::kMTauTau ] = fittedHiggs2tauSystem.mass(); 

  double bJet1Et = x[2*kMaxFitParams];
  xPrime[ svFitStandalone2b2tau::kBJet1Et ] = bJet1Et;
  double bJet2Et = x[2*kMaxFitParams + 1];
  xPrime[ svFitStandalone2b2tau::kBJet2Et ] = bJet2Et;

  return xPrime;
}

double
SVfitStandaloneLikelihood2b2tau::prob(const double* x) const 
{
  // in case of initialization errors don't start to do anything
  if ( error() ) { 
    return 0.;
  }
  //if ( verbose_ ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::prob(const double*)>:" << std::endl;
  //}
  ++idxObjFunctionCall_;
  //if ( verbose_ && FIRST ) {
  //  std::cout << " >> ixdObjFunctionCall : " << idxObjFunctionCall_ << std::endl;  
  //}
  // prevent kPhi in the fit parameters (kFitParams) from trespassing the 
  // +/-pi boundaries
  double phiPenalty = 0.;
  if ( addPhiPenalty_ ) {
    for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
      if ( TMath::Abs(idx*kMaxFitParams + x[kPhi]) > TMath::Pi() ) {
	phiPenalty += (TMath::Abs(x[kPhi]) - TMath::Pi())*(TMath::Abs(x[kPhi]) - TMath::Pi());
      }
    }
  }
  // xPrime are the transformed variables from which to construct the nll
  // transform performs the transformation from the fit parameters x to the 
  // nll parameters xPrime. prob is the actual combined likelihood. The
  // phiPenalty prevents the fit to converge to unphysical values beyond
  // +/-pi 
  double xPrime[svFitStandalone2b2tau::kMaxNLLParams + 4];
  const double* xPrime_ptr = transform(xPrime, x);
  if ( xPrime_ptr ) {
    return prob(xPrime_ptr, phiPenalty);
  } else {
    return 0.;
  }
}

double 
SVfitStandaloneLikelihood2b2tau::prob(const double* xPrime, double phiPenalty) const
{
  //if ( verbose_ && FIRST ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::prob(const double*, double)>:" << std::endl;
  //  std::cout << " isValidSolution1 = " << xPrime[svFitStandalone2b2tau::kMaxNLLParams + 2] << std::endl;
  //  std::cout << " isValidSolution2 = " << xPrime[svFitStandalone2b2tau::kMaxNLLParams + 3] << std::endl;
  //}
  if ( requirePhysicalSolution_ && (xPrime[ svFitStandalone2b2tau::kMaxNLLParams + 2 ] < 0.5 || xPrime[ svFitStandalone2b2tau::kMaxNLLParams + 3 ] < 0.5) ) return 0.;
  // start the combined likelihood construction from MET
  double prob = probMET(xPrime[svFitStandalone2b2tau::kDMETx], xPrime[svFitStandalone2b2tau::kDMETy], covDet_, invCovMET_, metPower_, (verbose_&& FIRST));
  //if ( verbose_ && FIRST ) {
  //  std::cout << "probMET = " << prob << std::endl;
  //}
  // add likelihoods for the decay branches
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    switch ( measuredTauLeptons_[idx].type() ) {
    case kTauToHadDecay :
      prob *= probTauToHadPhaseSpace(
                xPrime[idx == 0 ? svFitStandalone2b2tau::kDecayAngle1 : svFitStandalone2b2tau::kDecayAngle2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kNuNuMass1 : svFitStandalone2b2tau::kNuNuMass2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kVisMass1 : svFitStandalone2b2tau::kVisMass2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kMaxNLLParams : (svFitStandalone2b2tau::kMaxNLLParams + 1)], 
		addSinTheta_, 
		(verbose_&& FIRST));
      assert(!(marginalizeVisMass_ && shiftVisMass_));
      if ( marginalizeVisMass_ ) {
	prob *= probVisMass(
                  xPrime[idx == 0 ? svFitStandalone2b2tau::kDeltaVisMass1 : svFitStandalone2b2tau::kDeltaVisMass2], 
		  idx == 0 ? l1lutVisMass_ : l2lutVisMass_,
		  (verbose_&& FIRST));
      }
      if ( shiftVisMass_ ) {
	prob *= probVisMassShift(
                  xPrime[idx == 0 ? svFitStandalone2b2tau::kDeltaVisMass1 : svFitStandalone2b2tau::kDeltaVisMass2], 
		  idx == 0 ? l1lutVisMassRes_ : l2lutVisMassRes_,
		  (verbose_&& FIRST));
      }
      if ( shiftVisPt_ ) {
	prob *= probVisPtShift(
		  xPrime[idx == 0 ? svFitStandalone2b2tau::kRecTauPtDivGenTauPt1 : svFitStandalone2b2tau::kRecTauPtDivGenTauPt2], 
		  idx == 0 ? l1lutVisPtRes_ : l2lutVisPtRes_, 
		  (verbose_&& FIRST));
      }
      //if ( verbose_ && FIRST ) {
      //  std::cout << " *probTauToHad (" << idx << ") = " << prob << std::endl;
      //}
      break;
    case kTauToElecDecay :
    case kTauToMuDecay :
      prob *= probTauToLepMatrixElement(
		xPrime[idx == 0 ? svFitStandalone2b2tau::kDecayAngle1 : svFitStandalone2b2tau::kDecayAngle2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kNuNuMass1 : svFitStandalone2b2tau::kNuNuMass2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kVisMass1 : svFitStandalone2b2tau::kVisMass2], 
		xPrime[idx == 0 ? svFitStandalone2b2tau::kMaxNLLParams : (svFitStandalone2b2tau::kMaxNLLParams + 1)], 
		addSinTheta_, 
		(verbose_&& FIRST));
      //if ( verbose_ && FIRST ) {
      //  std::cout << " *probTauToLep (" << idx << ") = " << prob << std::endl;
      //}
      break;
    default :
      break;
    }
  }
  // add additional logM term if configured such 
  if ( addLogM_ && powerLogM_ > 0. ) {
    if ( xPrime[svFitStandalone2b2tau::kMTauTau] > 0. ) {
      prob *= TMath::Power(1.0/xPrime[svFitStandalone2b2tau::kMTauTau], powerLogM_);
    }
    //if ( verbose_ && FIRST ) {
    //  std::cout << " *logM = " << prob << std::endl;
    //}
  }
  if ( addDelta_ && mHiggsTauTau_ > 0. ) { // delta function derrivative of Higgs mass constraint for H -> tautau system
    prob *= (2.*xPrime[svFitStandalone2b2tau::kMaxNLLParams + 1]/xPrime[svFitStandalone2b2tau::kMTauTau]);
    //if ( verbose_ && FIRST ) {
    //  std::cout << " *delta (tau) = " << prob << std::endl;
    //}
  }
  // add additional phiPenalty in case kPhi in the fit parameters 
  // (kFitParams) trespassed the physical boundaries from +/-pi 
  if ( phiPenalty > 0. ) {
    prob *= TMath::Exp(-phiPenalty);
    //if ( verbose_ && FIRST ) {
    //  std::cout << " *phiPenalty = " << prob << std::endl;
    //}
  }
  // add likelihood for b-jet1 and b-jet2
  double bJet1Et_fitted = xPrime[svFitStandalone2b2tau::kBJet1Et];
  double bJet2Et_fitted = xPrime[svFitStandalone2b2tau::kBJet2Et];
  prob *= bJetEtResolution_(bJet1Et_fitted, bJet1Eta_, bJet1Et_unfitted_);
  //if ( verbose_ && FIRST ) {
  //  std::cout << " *bJetEtResolution (1) = " << prob << std::endl;
  //}
  prob *= bJetEtResolution_(bJet2Et_fitted, bJet2Eta_, bJet2Et_unfitted_);
  //if ( verbose_ && FIRST ) {
  //  std::cout << " *bJetEtResolution (2) = " << prob << std::endl;
  //}
  if ( addDelta_ && mHiggsBB_ > 0. ) { // delta function derrivative of Higgs mass constraint for H -> bb system
    prob *= (1./(2.*TMath::Sqrt(bJet2Et_fitted)))*(TMath::Sqrt(bJet1Et_fitted/(bJet1Et_unfitted_*bJet2Et_unfitted_))*diBJetMass_unfitted_);
    //if ( verbose_ && FIRST ) {
    //  std::cout << " *delta (b-jet) = " << prob << std::endl;
    //}
  }
  // set FIRST to false after the first complete evaluation of the likelihood 
  //assert(FIRST);
  FIRST = false;
  return prob;
}

void
SVfitStandaloneLikelihood2b2tau::results2tau(std::vector<LorentzVector>& fittedTauLeptons, const double* x) const
{
  //if ( verbose_ ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::results2tau(std::vector<LorentzVector>&, const double*)>:" << std::endl;
  //}
  for ( size_t idx = 0; idx < measuredTauLeptons_.size(); ++idx ) {
    const MeasuredTauLepton& measuredTauLepton = measuredTauLeptons_[idx];

    // map to local variables to be more clear on the meaning of the individual parameters. The fit parameters are ayered 
    // for each tau decay
    double nunuMass                 = x[ idx*kMaxFitParams + kMNuNu ];       // nunu inv mass (can be const 0 for had tau decays) 
    double labframeXFrac            = x[ idx*kMaxFitParams + kXFrac ];       // visible energy fraction x in labframe
    double labframePhi              = x[ idx*kMaxFitParams + kPhi   ];       // phi in labframe 
    double visMass                  = measuredTauLepton.mass(); 
    double labframeVisMom_unshifted = measuredTauLepton.momentum(); 
    double labframeVisMom           = labframeVisMom_unshifted; // visible momentum in lab-frame
    double labframeVisEn            = measuredTauLepton.energy(); // visible energy in lab-frame    
    if ( measuredTauLepton.type() == kTauToHadDecay ) {
      if ( marginalizeVisMass_ || shiftVisMass_ ) {
	visMass = x[idx*kMaxFitParams + kVisMassShifted];
      }
      if ( shiftVisPt_ ) {
	double shiftInv = 1. + x[idx*kMaxFitParams + kRecTauPtDivGenTauPt];
        double shift = ( shiftInv > 1.e-1 ) ?
	  (1./shiftInv) : 1.e+1;
        labframeVisMom *= shift;
        //visMass *= shift; // CV: take mass and momentum to be correlated
        //labframeVisEn = TMath::Sqrt(labframeVisMom*labframeVisMom + visMass*visMass);
        labframeVisEn *= shift;
      }
    }
    if ( visMass < 5.1e-4 ) { 
      visMass = 5.1e-4; 
    } 
    bool isValidSolution = true;
    double gjAngle_lab = gjAngleLabFrameFromX(labframeXFrac, visMass, nunuMass, labframeVisMom, labframeVisEn, tauLeptonMass, isValidSolution);
    double enTau_lab = labframeVisEn/labframeXFrac;
    double pTau_lab = TMath::Sqrt(square(enTau_lab) - tauLeptonMass2);
    Vector p3Tau_unit = motherDirection(measuredTauLepton.direction(), gjAngle_lab, labframePhi);
    LorentzVector p4Tau_lab = motherP4(p3Tau_unit, pTau_lab, enTau_lab);
    // tau lepton four vector in labframe
    if ( idx < fittedTauLeptons.size() ) fittedTauLeptons[idx] = p4Tau_lab;
    else fittedTauLeptons.push_back(p4Tau_lab);
  }
}

void
SVfitStandaloneLikelihood2b2tau::results2b(std::vector<LorentzVector>& fittedBJets, const double* x) const
{
  //if ( verbose_ ) {
  //  std::cout << "<SVfitStandaloneLikelihood2b2tau::results2b(std::vector<LorentzVector>&, const double*)>:" << std::endl;
  //}
  for ( size_t idx = 0; idx < measuredBJets_.size(); ++idx ) {
    const MeasuredBJet& measuredBJet = measuredBJets_[idx];

    double bJetEt_fitted = x[ 2*kMaxFitParams + idx ];
    double bJetEt_unfitted = measuredBJet.Et();
    double sf = bJetEt_fitted/bJetEt_unfitted;
    LorentzVector p4BJet_fitted(sf*measuredBJet.px(), sf*measuredBJet.py(), sf*measuredBJet.pz(), sf*measuredBJet.energy());

    if ( idx < fittedBJets.size() ) fittedBJets[idx] = p4BJet_fitted;
    else fittedBJets.push_back(p4BJet_fitted);
  }
}
