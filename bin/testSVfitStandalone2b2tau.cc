
/**
   \class testSVfitStandalone testSVfitStandalone.cc "TauAnalysis/SVfitStandalone/bin/testSVfitStandalone.cc"
   \brief Basic example of the use of the standalone version of SVfit

   This is an example executable to show the use of the standalone version of SVfit 
   from a flat n-tuple or single event.
*/

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/SVfit2b2tau_standalone/interface/SVfitStandaloneAlgorithm2b2tau.h"

#include "TFile.h"
#include "TH1.h"

void singleEvent()
{
  /* 
     This is a single event for testing in the integration mode.
  */

  // define MET
  double measuredMETx = -53.3447;
  double measuredMETy =  25.3619; 
  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] = 670.5;
  covMET[1][0] =   8.9;
  covMET[0][1] =   8.9;
  covMET[1][1] = 558.1;
  // define lepton four vectors
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, 61.2491,  0.2258, -2.6170, 1.1381));
  measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, 51.5053, -0.9034,  2.7945, 1.2575));
  std::vector<svFitStandalone2b2tau::MeasuredBJet> measuredBJets;
  measuredBJets.push_back(svFitStandalone2b2tau::MeasuredBJet(121.528, 0.3238, 1.3196, 13.1555));
  measuredBJets.push_back(svFitStandalone2b2tau::MeasuredBJet(157.978, 0.5676, 0.4814, 22.1124));
  // define algorithm (set the debug level to 3 for testing)
  unsigned verbosity = 2;
  SVfitStandaloneAlgorithm2b2tau algo(measuredTauLeptons, measuredMETx, measuredMETy, covMET, measuredBJets, verbosity);
  algo.setMassHiggsTauTau(-1.);
  algo.setMassHiggsBB(125.);
  //algo.setMassHiggsBB(-1.);
  //edm::FileInPath inputFileName_visPtResolution("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root");
  //TH1::AddDirectory(false);  
  //TFile* inputFile_visPtResolution = new TFile(inputFileName_visPtResolution.fullPath().data());
  //algo.shiftVisPt(true, inputFile_visPtResolution);
  algo.integrateMarkovChain();

  double mass = algo.diHiggsMass(); // return value is in units of GeV
  if ( algo.isValidSolution() ) {
    std::cout << "found mass = " << mass << " (expected value = 485.299)" << std::endl;
  } else {
    std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  }
  
  //delete inputFile_visPtResolution;
}

int main(int argc, char* argv[]) 
{
  singleEvent();
  return 0;
}
