// -*- C++ -*-
//
// Package:    McLevelEff
// Class:      MCAnalyTrial
// 
/**\class MCAnalyTrial MCAnalyTrial.cc ZShape/MCAnalyTrial/src/MCAnalyTrial.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/



// system include files
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/Handle.h"

#include <CLHEP/Vector/LorentzVector.h>

#include "ZShape/Base/interface/ZShapeEvent.h"
#include "ZShape/Base/interface/ZShapeZDef.h"
#include "ZShape/Base/interface/ZShapeStandardCuts.h"
#include "ZShape/EffAcc/interface/EffHistos.h"
#include "ZShape/EffAcc/interface/EfficiencyStore.h"
#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "ZShape/ZFromData/interface/EffExtraHistos.h"

#include "PhysicsTools/TagAndProbe/interface/EffTableLoader.h"

///#include "AnalysisDataFormats/ElectronEfficiency/interface/EmObjectFwd.h"
///#include "AnalysisDataFormats/ElectronEfficiency/interface/TagProbeAssociation.h"
///#include "MyAnalyzer/TagProbeCollectionCuts/interface/HLTisPassed.h"
///#include "MyAnalyzer/TagProbeCollectionCuts/interface/TrackisIsolated.h"

#include "TRandom3.h"
#include "TFile.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

class TProfile;

//
// class decleration
//

class MCAnalyTrial : public edm::EDAnalyzer {
public:
  explicit MCAnalyTrial(const edm::ParameterSet&);
  ~MCAnalyTrial();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  int ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, edm::Handle<reco::CandidateView>& passprobes, bool checkExactOverlap_ = true );

  bool MatchObjects( const reco::Candidate *hltObj,  const reco::CandidateBaseRef& tagObj, bool exact = true );
  

  // ----------member data ---------------------------
  ZShapeEvent evt_;
  ZShapeEvent evtMC_;

  edm::InputTag MCTag_, SCTag_;

  TH1F *Hgood_,*Hall_,*Heff_; 
  TH1F *SCtoMCERatio_[40];
  TH2F *HMCptvsSCpt_;

  //std::map<std::string, EffInfo*> effInfo_;
  edm::InputTag tnpProducer_;

  //Iso and hlt information
  ///HLTisPassed _HLTpassed;
  ///TrackisIsolated _Isopassed;

  TRandom3 randomNum;

  edm::InputTag m_srcTag;
  bool quiet_;
  std::string outFileName_;
  bool        writeHistoConservatively_;
  TFile*      histoFile_;

  // the efficiency objects: strings identify the cuts efficiencies and corresponding cuts
  std::map<std::string, EfficiencyStore*> efficiencies_;
  std::map<std::string, EfficiencyCut*> theCuts_;
  
  // object which applies the eta and pt cuts
  ZShapeStandardCuts stdCuts_;

  // Z Definitions
  std::map<std::string, ZShapeZDef*> zdefs_;
  std::map<std::string, bool> zptorder_;

  // when doing efficiency statistics runs
  struct StatsBox {
    std::string targetEffStat, targetZDefStat;
    int trials;
    std::vector<EfficiencyCut*> alternateCuts;
    std::vector<EffHistos> hists;
    std::vector<TH1*> rawHists;
    TProfile* cutsProfile;
  } statsBox_;
  
  // targets for systematic variations
  std::string targetEffSys_, targetZDefSys_;

  // Z Plots:

  // before any z-def selection
  EffHistos allCase_;
  EffHistos allCaseFirst_;
  EffExtraHistos allCaseExtra_;
  EffExtraHistos allCaseExtraFirst_;

  // at each step of selection for any z-definition
  struct ZPlots {
    ZPlots(int nc) : postCut_(nc), postCutExtra_(nc) { }
    ///ZPlots(int nc) : postCut_(nc) { }
    EffHistos acceptance_;
    EffExtraHistos acceptanceExtra_; //added for the extra histos
    std::vector<EffHistos> postCut_;
    std::vector<EffExtraHistos> postCutExtra_; //added for the extra histos
  };
  std::map<std::string,ZPlots*> zplots_;
  TRandom3 *random_;

  bool extraHistos_;
  
  
  // Matching parameters
  double delRMatchingCut_;
  double delPtRelMatchingCut_;
  
  // Tag probe map tags
  std::vector<edm::InputTag> tagProbeMapTags_;

  // Candidate collection tags
  edm::InputTag genParticlesTag_;
  std::vector<edm::InputTag> tagCandTags_;
  std::vector<edm::InputTag> allProbeCandTags_;
  std::vector<edm::InputTag> passProbeCandTags_;
  
  //Vector of strings to hold the efficiency names
  std::vector<std::string> cutName_;
  
  //Vector of for determining "exact" matching
  std::vector<int> exactMatch_;

  double scEtCut_;
  std::string textBinsFile_;       // This is the name of the file that holds the 2D bin information
  EffTableLoader* effBinsFromTxt_; // This holds the efficiency bins information

  std::string cbBinsFile_; //This file will hold the outputed crystal ball fit information
  bool fromfile_;
  std::string fromfilefile_;
  TFile *myFile_;
  
};


/** TODO: 

    Efficiencies need to specify dependencies (eta, pT, etc) [current: hack]

    Randomization of efficiencies

 */
