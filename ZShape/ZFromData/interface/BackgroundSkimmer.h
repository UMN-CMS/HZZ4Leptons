// -*- C++ -*-
//
// Package:    McLevelEff
// Class:      BackgroundSkimmer
// 
/**\class BackgroundSkimmer BackgroundSkimmer.cc ZShape/BackgroundSkimmer/src/BackgroundSkimmer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/



// system include files
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/Handle.h"

#include <CLHEP/Vector/LorentzVector.h>


///#include "AnalysisDataFormats/ElectronEfficiency/interface/EmObjectFwd.h"
///#include "AnalysisDataFormats/ElectronEfficiency/interface/TagProbeAssociation.h"
///#include "MyAnalyzer/TagProbeCollectionCuts/interface/HLTisPassed.h"
///#include "MyAnalyzer/TagProbeCollectionCuts/interface/TrackisIsolated.h"



#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

class TProfile;

//
// class decleration
//

class BackgroundSkimmer : public edm::EDFilter {
public:
  explicit BackgroundSkimmer(const edm::ParameterSet&);
  ~BackgroundSkimmer();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  std::vector<edm::InputTag> allProbeCandTags_;
  
};


/** TODO: 

    Efficiencies need to specify dependencies (eta, pT, etc) [current: hack]

    Randomization of efficiencies

 */
