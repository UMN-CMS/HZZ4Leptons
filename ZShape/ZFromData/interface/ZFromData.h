// -*- C++ -*-
//
// Package:    McLevelEff
// Class:      ZFromData
// 
/**\class ZFromData ZFromData.cc ZShape/ZFromData/src/ZFromData.cc

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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/Handle.h"

#include <CLHEP/Vector/LorentzVector.h>

#include "ZShape/Base/interface/ZShapeEvent.h"
#include "ZShape/Base/interface/ZShapeZDef.h"
#include "ZShape/Base/interface/ZShapeStandardCuts.h"
#include "ZShape/EffAcc/interface/EffHistos.h"
#include "ZShape/EffAcc/interface/EfficiencyStore.h"
#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "ZShape/ZFromData/interface/EffExtraHistos.h"

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

class ZFromData : public edm::EDAnalyzer {
public:
  explicit ZFromData(const edm::ParameterSet&);
  ~ZFromData();


private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

   struct EffInfo {
    EffInfo(){}
    edm::InputTag trialCollection;
    edm::InputTag passCollection;
    bool needIso;
    bool needHlt;
  };

  void fillMCEvent(const HepMC::GenEvent* evt);
  void applyEfficiencies();
  void createAlternateEfficiencies();
  void createAlternateZDefs(const std::string& targetZDefSys, const std::string& targetEffSys);
  ///bool oneCut(EffInfo*, const edm::Event &,const EgEff::EmObjectRef &, math::PtEtaPhiMLorentzVector*);
  ///bool isolated(const edm::Event &, const EgEff::EmObjectRef &);
  ///bool hlt(const edm::Event &, const EgEff::EmObjectRef &);
  int ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, edm::Handle<reco::CandidateView>& passprobes, bool checkExactOverlap_ = true );

  bool MatchObjects( const reco::Candidate *hltObj,  const reco::CandidateBaseRef& tagObj, bool exact = true );
  

  // ----------member data ---------------------------
  ZShapeEvent evt_;
  ZShapeEvent evtMC_;

 

  std::map<std::string, EffInfo*> effInfo_;
  edm::InputTag tnpProducer_;
  edm::InputTag gsfProducer_; //added to count the number of gsf's above 20

  //Iso and hlt information
  ///HLTisPassed _HLTpassed;
  ///TrackisIsolated _Isopassed;

  TRandom3 randomNum;

  edm::InputTag m_srcTag;
  bool quiet_;
  std::string outFileName_;
  bool        writeHistoConservatively_;
  TFile*      histoFile_;

  
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


  struct RegionalPlots {
    TH1F* ecal_ecal, *ecal_ntrk, *ecal_hf, *ecal_noacc;
    TH1F* ntrk_ntrk, *ntrk_hf, *ntrk_noacc;
    TH1F* hf_hf, *hf_noacc, *noacc_noacc;
  } accHistos_;

  // at each step of selection for any z-definition
  struct ZPlots {
    ZPlots(int nc, int nz) : postCut_(nc), postCutExtra_(nc), zCut_(nz), zCutExtra_(nz) { }
    ///ZPlots(int nc) : postCut_(nc) { }
    EffHistos acceptance_;
    EffExtraHistos acceptanceExtra_; //added for the extra histos
    std::vector<EffHistos> postCut_;
    std::vector<EffExtraHistos> postCutExtra_; //added for the extra histos
    std::vector<EffHistos> zCut_;
    std::vector<EffExtraHistos> zCutExtra_; //added for the extra histos
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
  
};


/** TODO: 

    Efficiencies need to specify dependencies (eta, pT, etc) [current: hack]

    Randomization of efficiencies

 */
