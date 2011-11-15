// -*- C++ -*-
//
// Package:    McLevelEff
// Class:      ZEfficiencyCalculator
// 
/**\class ZEfficiencyCalculator ZEfficiencyCalculator.cc ZShape/EffAcc/src/ZEfficiencyCalculator.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Giovanni FRANZONI
//         Created:  Thu Oct  4 11:30:13 CEST 2007
// $Id: ZEfficiencyCalculator.h,v 1.20 2011/01/27 17:28:25 mansj Exp $
//
//


// system include files
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/Common/interface/Handle.h"

#include "ZShape/Base/interface/ZShapeEvent.h"
#include "ZShape/Base/interface/ZShapeZDef.h"
#include "ZShape/Base/interface/ZShapeStandardCuts.h"
#include "ZShape/EffAcc/interface/EffHistos.h"
#include "ZShape/EffAcc/interface/EfficiencyStore.h"
#include "ZShape/EffAcc/interface/EfficiencyCut.h"
#include "ZShape/Base/interface/Systematics.h"

#include "TRandom3.h"
#include "TFile.h"

class TProfile;

//
// class decleration
//

class ZEfficiencyCalculator : public edm::EDAnalyzer {
public:
  explicit ZEfficiencyCalculator(const edm::ParameterSet&);
  ~ZEfficiencyCalculator();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillEvent(const reco::GenParticleCollection* ZeeParticles, const reco::GenParticleCollection* ZeeTreeLevelParticles); 
  void loadEfficiency(const std::string& name, const std::string& fname);
  void applyEfficiencies();
  void createAlternateEfficiencies(int cycle, TFileDirectory& fd);
  void createAlternateZDefs(const std::string& targetZDefSys, const std::string& targetEffSys);

  double getPDFWeight(const edm::Event&);
  
  // ----------member data ---------------------------
  ZShapeEvent evt_;

  TRandom3 randomNum;

  edm::InputTag m_srcTag;
  edm::InputTag zElectronsTag;
  edm::InputTag zTreeLevelElectronsTag;
  edm::InputTag partonIdTag,partonMomFractionTag,scaleQTag;
  bool doPDFreweight_;
  std::string pdfReweightBaseName, pdfReweightTargetName;
  int pdfReweightBaseId, pdfReweightTargetId;
  bool pdfReweightAddZmass_;
  bool quiet_;
  float zElectronsCone_;

  // the efficiency objects: strings identify the cuts efficiencies and corresponding cuts
  std::map<std::string, EfficiencyStore*> efficiencies_;
  std::map<std::string, EfficiencyCut*> theCuts_;
  
  // object which applies the eta and pt cuts
  ZShapeStandardCuts stdCuts_;

  // Z Definitions
  std::map<std::string, ZShapeZDef*> zdefs_;

  // when doing efficiency statistics runs
  struct StatsBox {
    std::string targetEffStat;
    int trials;
    std::vector<EfficiencyCut*> alternateCuts;
    std::map<std::string, std::vector<EffHistos> > hists;
    std::vector<TH1*> rawHists;
    TProfile* cutsProfile;
  } statsBox_;
  
  // targets for systematic variations
  std::string targetEffSys_, targetZDefSys_;

  // other systematics
  std::string systematicVariation_;
  struct SystematicsBlock {
    zshape::EnergyScale* energyScale;
  } m_systematics;

  // Z Plots:

  // before any z-def selection
  EffHistos allCase_;

  struct RegionalPlots {
    TH1F* ecal_ecal, *ecal_ntrk, *ecal_hf, *ecal_noacc;
    TH1F* ntrk_ntrk, *ntrk_hf, *ntrk_noacc;
    TH1F* hf_hf, *hf_noacc, *noacc_noacc;
    TH1F* weights;
  } accHistos_;
    

  // at each step of selection for any z-definition
  struct ZPlots {
    ZPlots(int nc, int nz) : postCut_(nc),zCut_(nz) { }
    EffHistos acceptance_;
    std::vector<EffHistos> postCut_;
    std::vector<EffHistos> zCut_;
  };
  std::map<std::string,ZPlots*> zplots_;

  // invariant mass cuts
  double massLow_ ;
  double massHigh_;
};


