// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNu
// 
/**\class HeavyNu HeavyNu.cc HeavyNu/AnalyzerModules/src/HeavyNu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNu.cc,v 1.18 2011/01/13 12:40:23 dudero Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// According to
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
// this must be included before Frameworkfwd.h
//
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"


//////////////////////////////////////////////////////////////////
// generic maximum/minimum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;
}
template <class T> const T& min ( const T& a, const T& b ) {
  return (b<a)?b:a;
}

template <class T>
inline std::string int2str(T i) {
  std::ostringstream ss;
  ss << i;
  return ss.str();
}

//============================================================

class compare {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

//============================================================

inline void outputCandidate(const reco::CandidateBaseRef& rc) {
  std::cout<<"pt="<<rc->pt()<<" GeV, eta="<<rc->eta()<<", phi="<<rc->phi();
}

//============================================================

inline void dumpJetCorInfo(const pat::Jet& j) {
  std::cout << "Available JEC sets and levels:\n";
  const std::vector<std::string> jecsets = j.availableJECSets();
  for (size_t i=0; i<jecsets.size(); i++) {
    std::cout << jecsets[i]<<":";
    const std::vector<std::string> jeclevs = j.availableJECLevels(i);
    for (size_t j=0; j<jeclevs.size(); j++)
      std::cout << " " << jeclevs[i];
    std::cout << std::endl;
  }
  std::cout << "current JEC set   : " << j.currentJECSet()    << std::endl;
  std::cout << "current JEC level : " << j.currentJECLevel()  << std::endl;
  std::cout << "current JEC flavor: " << j.currentJECFlavor() << std::endl;
}

//============================================================

// Returns 0=neither, 1=loose or 2=tight, -1 if tight but not loose (!)
int jetID(const pat::Jet& j)
{
  JetIDSelectionFunctor jetIDloose(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::LOOSE);
  JetIDSelectionFunctor jetIDtight(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::TIGHT);

  pat::strbitset ret = jetIDloose.getBitTemplate();
  ret.set(false);  bool loose = jetIDloose(j, ret);
  ret.set(false);  bool tight = jetIDtight(j, ret);
  return (tight ? (loose ? 2 : -1) : (loose ? 1 : 0));
}

inline void labelJetIDaxis(TAxis *ax)
{
  ax->SetBinLabel(1,"Neither");
  ax->SetBinLabel(2,"PURE09 Loose");
  ax->SetBinLabel(3,"PURE09 Tight");
}

//============================================================


static std::string btagName;

class HeavyNu : public edm::EDFilter {
public:
  explicit HeavyNu(const edm::ParameterSet&);
  ~HeavyNu();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual bool isVBTFloose(const pat::Muon& m);
  virtual bool isVBTFtight(const pat::Muon& m);
  virtual TH1 *bookRunHisto(uint32_t runNumber);
  
  inline bool inZmassWindow( double mMuMu ) {
    return (mMuMu <= ZwinMaxGeV_) && (mMuMu >= ZwinMinGeV_);
  }

  edm::InputTag muonTag_;
  edm::InputTag jetTag_;
  edm::InputTag elecTag_;

  double ZwinMinGeV_, ZwinMaxGeV_;

  std::string currentFile_;
  bool dolog_;
  HeavyNu_NNIF *nnif_;
  HeavyNuTrigger *trig_;

  std::map<uint32_t,TH1 *> m_runHistos_;

  // ----------member data ---------------------------
  static const std::string muonQuality[] ; 
  static const int muonQualityFlags ; 

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&) ;
    // fill all histos of the set with the two electron candidates
    void fill(pat::MuonCollection muons, pat::JetCollection jets,bool isMC) ;
    // fill all histos of the set with the two electron candidates
    void fill(const HeavyNuEvent& hne, const std::vector<hNuMassHypothesis>&) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;
    TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet ;
    TH1 *dEtaPhiMu, *dEtaPhiJet ; 
    TH1 *dRminMu1jet, *dRminMu2jet ; 
    TH1 *hptrelMu1, *hptrelMu2 ; 
    TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet ;
    TH2 *jetID2d;

    TH1 *dptMu1gen, *dptMu2gen ; 
    TH1 *dRMu1gen, *dRMu2gen ; 
    TH1 *qualMu1, *qualMu2 ; 

    TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
    TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

    TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
    TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

    TH1 *mWR, *mNuR1, *mNuR2, *mMuMu, *mMuMuZoom, *mJJ ; 
    TH2 *mNuR2D, *jetPtvsNum ; 

    TH1* btagJet1, *btagJet2;

    TH1* czeta_mumu;
    TH1* czeta_mumu_zoom;

    // Jeremy's crazy angles...
    TH1* ctheta_mumu, *cthetaz_mumu;
    TH1* ctheta_jj, *cthetaz_jj;
    TH1* ctheta_mu1_jj, *cthetaz_mu1_jj;
    TH1* ctheta_mu2_jj, *cthetaz_mu2_jj;

    TFileDirectory *mydir;
    TFileDirectory *nndir;

    HeavyNuTrigger::trigHistos_t trigHistos;
  };

  bool init_;

  // gf set of histo for all Z definitions in a stack
  struct HistStruct {
    TH1 *nelec, *nmu, *njet ;
    TH1 *muPt, *muEta, *muPhi ; 
    TH1 *jetPt, *jetEta, *jetPhi, *jetID ; 
    TH2 *jetPtvsNum;
    TFileDirectory *rundir;
    HistPerDef noCuts; 
    HistPerDef LLptCuts;
    HistPerDef MuTightCuts;
    HistPerDef TrigMatches;
    HistPerDef Mu1TrigMatchesInZwin;
    HistPerDef Mu2TrigMatchesInZwin;
    HistPerDef Mu1Mu2TrigMatchesInZwin;
    HistPerDef MuTightInZwin;
    HistPerDef LLJJptCuts;
    HistPerDef diLmassCut;
    HistPerDef mWRmassCut;
    HistPerDef oneBtag;
    HistPerDef twoBtag;
  } hists;

  struct CutsStruct {
    double minimum_mu1_pt;
    double minimum_mu2_pt;
    double minimum_jet_pt;
    double maximum_mu_abseta;
    double maximum_jet_abseta;
    double minimum_mumu_mass;
    double minimum_mWR_mass;
    double minimum_muon_jet_dR;
    double muon_trackiso_limit;
  } cuts;
  
};

//======================================================================

const std::string HeavyNu::muonQuality[] = {"AllGlobalMuons","AllStandAloneMuons","AllTrackerMuons"} ; 
const int HeavyNu::muonQualityFlags = 3 ;
inline std::string nnhistoname(int mwr,int mnu) {
  return ("WR"+int2str<int>(mwr)+"nuRmu"+int2str<int>(mnu));
}

//======================================================================

void
HeavyNu::HistPerDef::book(TFileDirectory *td,
			  const std::string& post,
			  const std::vector<hNuMassHypothesis>& v_masspts )
{
  std::string t; // histogram title string;
  
  TH1::SetDefaultSumw2();

  mydir = td;

  // ----------  Muon histograms  ----------

  t="p_{T}(#mu_{1}) "+post;   ptMu1=td->make<TH1D>("ptMu1",t.c_str(),100,0.,1000.);
  t="p_{T}(#mu_{2}) "+post;   ptMu2=td->make<TH1D>("ptMu2",t.c_str(),100,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1=td->make<TH1D>("etaMu1",t.c_str(),40,-2.5,2.5);
  t="#eta(#mu_{2}) " +post;  etaMu2=td->make<TH1D>("etaMu2",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1=td->make<TH1D>("phiMu1",t.c_str(),30,-3.14159,3.14159);
  t="#phi(#mu_{2}) " +post;  phiMu2=td->make<TH1D>("phiMu2",t.c_str(),30,-3.14159,3.14159);

  // delta angles

  t="#Delta#eta(#mu_{1},#mu_{2}) "  +post;    dEtaMu=td->make<TH1D>("dEtaMu",t.c_str(),40,0,5);
  t="#Delta#phi(#mu_{1},#mu_{2}) "  +post;    dPhiMu=td->make<TH1D>("dPhiMu",t.c_str(),30,0,3.14159);
  t="#Delta p_{T}(#mu_{1},gen) "    +post; dptMu1gen=td->make<TH1D>("dptMu1gen",t.c_str(),50,-0.50,0.50);
  t="#Delta p_{T}(#mu_{2},gen) "    +post; dptMu2gen=td->make<TH1D>("dptMu2gen",t.c_str(),50,-0.50,0.50);
  t="#Delta R(#mu_{1},gen) "        +post;  dRMu1gen=td->make<TH1D>("dRMu1gen", t.c_str(),50,0,0.01);
  t="#Delta R(#mu_{2},gen) "        +post;  dRMu2gen=td->make<TH1D>("dRMu2gen", t.c_str(),50,0,0.01);
  t="#mu #Delta#eta vs. #Delta#phi "+post;
  t+=  ";#Delta#eta; #Delta#phi";          dEtaPhiMu=td->make<TH2D>("dEtaPhiMu",t.c_str(),
								   50,0,5,30,0,3.14159);

  t="Quality (#mu_{1}) "+post; qualMu1=td->make<TH1D>("qualMu1",t.c_str(),muonQualityFlags,0,muonQualityFlags);
  t="Quality (#mu_{2}) "+post; qualMu2=td->make<TH1D>("qualMu2",t.c_str(),muonQualityFlags,0,muonQualityFlags);
  for (int i=0; i<muonQualityFlags; i++) {
    qualMu1->GetXaxis()->SetBinLabel(i+1,muonQuality[i].c_str()) ;
    qualMu2->GetXaxis()->SetBinLabel(i+1,muonQuality[i].c_str()) ;
  }

  // isolation

  t=   "trackIso(#mu_{1}) "+post;    mu1trackIso=td->make<TH1D>("mu1trackIso",   t.c_str(),40,0.,200.);
  t=    "hcalIso(#mu_{1}) "+post;     mu1hcalIso=td->make<TH1D>("mu1hcalIso",    t.c_str(),40,0.,200.);
  t=    "ecalIso(#mu_{1}) "+post;     mu1ecalIso=td->make<TH1D>("mu1ecalIso",    t.c_str(),40,0.,200.);
  t=    "caloIso(#mu_{1}) "+post;     mu1caloIso=td->make<TH1D>("mu1caloIso",    t.c_str(),40,0.,200.);
  t=        "Dxy(#mu_{1}) "+post;          mu1dB=td->make<TH1D>("mu1dB",         t.c_str(),50,-5.,5.);

  t=   "trackIso(#mu_{2}) "+post;    mu2trackIso=td->make<TH1D>("mu2trackIso",   t.c_str(),40,0.,200.);
  t=    "hcalIso(#mu_{2}) "+post;     mu2hcalIso=td->make<TH1D>("mu2hcalIso",    t.c_str(),40,0.,200.);
  t=    "ecalIso(#mu_{2}) "+post;     mu2ecalIso=td->make<TH1D>("mu2ecalIso",    t.c_str(),40,0.,200.);
  t=    "caloIso(#mu_{2}) "+post;     mu2caloIso=td->make<TH1D>("mu2caloIso",    t.c_str(),40,0.,200.);
  t=        "Dxy(#mu_{2}) "+post;          mu2dB=td->make<TH1D>("mu2dB",         t.c_str(),50,-5.,5.);

  t="trackRelIso(#mu_{1}) "+post; mu1trackRelIso=td->make<TH1D>("mu1trackRelIso",t.c_str(),50,0.,5.);
  t= "hcalRelIso(#mu_{1}) "+post;  mu1hcalRelIso=td->make<TH1D>("mu1hcalRelIso", t.c_str(),50,0.,5.);
  t= "ecalRelIso(#mu_{1}) "+post;  mu1ecalRelIso=td->make<TH1D>("mu1ecalRelIso", t.c_str(),50,0.,5.);
  t= "caloRelIso(#mu_{1}) "+post;  mu1caloRelIso=td->make<TH1D>("mu1caloRelIso", t.c_str(),50,0.,5.);

  t="trackRelIso(#mu_{2}) "+post; mu2trackRelIso=td->make<TH1D>("mu2trackRelIso",t.c_str(),50,0.,5.);
  t= "hcalRelIso(#mu_{2}) "+post;  mu2hcalRelIso=td->make<TH1D>("mu2hcalRelIso", t.c_str(),50,0.,5.);
  t= "ecalRelIso(#mu_{2}) "+post;  mu2ecalRelIso=td->make<TH1D>("mu2ecalRelIso", t.c_str(),50,0.,5.);
  t= "caloRelIso(#mu_{2}) "+post;  mu2caloRelIso=td->make<TH1D>("mu2caloRelIso", t.c_str(),50,0.,5.);

  // ----------  Jet histograms ----------

  t="p_{T}(j_{1}) "            +post;     ptJet1=td->make<TH1D>("ptJet1",  t.c_str(),50,0.,500.);
  t="p_{T}(j_{2}) "            +post;     ptJet2=td->make<TH1D>("ptJet2",  t.c_str(),50,0.,500.);
  t= "#eta(j_{1}) "            +post;    etaJet1=td->make<TH1D>("etaJet1", t.c_str(),40,-5,5);
  t= "#eta(j_{2}) "            +post;    etaJet2=td->make<TH1D>("etaJet2", t.c_str(),40,-5,5);
  t= "#phi(j_{1}) "            +post;    phiJet1=td->make<TH1D>("phiJet1", t.c_str(),30,-3.14159,3.14159);
  t= "#phi(j_{2}) "            +post;    phiJet2=td->make<TH1D>("phiJet2", t.c_str(),30,-3.14159,3.14159);

  t="#Delta#eta(j_{1},j_{2}) " +post;    dEtaJet=td->make<TH1D>("dEtaJet", t.c_str(),40,0,5);
  t="#Delta#phi(j_{1},j_{2}) " +post;    dPhiJet=td->make<TH1D>("dPhiJet", t.c_str(),30,0,3.14159);

  t= "btag(j_{1}) "            +post;    btagJet1=td->make<TH1D>("btagJet1",t.c_str(),40,0,5);
  t= "btag(j_{2}) "            +post;    btagJet2=td->make<TH1D>("btagJet2",t.c_str(),40,0,5);

  t ="Jet #Delta#eta vs. #Delta#phi ";
  t+= post+";#Delta#eta; #Delta#phi";  dEtaPhiJet=td->make<TH2D>("dEtaPhiJet",t.c_str(),
								50,0,5,30,0,3.14159);
  t ="Jet ID(j_{2}) vs. ID(j_{1}) ";
  t+=post+"; ID(j_{1}); ID(j_{2})";         jetID2d=td->make<TH2I>("jetID2d",t.c_str(),3,0,3,3,0,3);
  labelJetIDaxis(jetID2d->GetXaxis());
  labelJetIDaxis(jetID2d->GetYaxis());

  // ----------  Mu/Jet histograms  ----------

  t="Minimum #Delta R(#mu_{1},jet) "+post;  dRminMu1jet=td->make<TH1D>("dRminMu1jet",t.c_str(),50,0,5.);
  t="Minimum #Delta R(#mu_{2},jet) "+post;  dRminMu2jet=td->make<TH1D>("dRminMu2jet",t.c_str(),50,0,5.);

  t="p_{T,rel}(#mu_{1},jet)"+post;            hptrelMu1=td->make<TH1D>("ptrelMu1",t.c_str(),50,0,1000.);
  t="p_{T,rel}(#mu_{2},jet)"+post;            hptrelMu2=td->make<TH1D>("ptrelMu2",t.c_str(),50,0,1000.);

  t="p_{T,rel}(#mu_{1},jet) vs #Delta R(#mu_{1},jet)"+post;
  t+="; #Delta R(#mu_{1},jet); p_{T,rel}(#mu_{1},jet)"; ptrelVsdRminMu1jet=td->make<TH2D>("ptrelVsdRminMu1jet",
											 t.c_str(),
											 50,0,5.,50,0,1000);
  t="p_{T,rel}(#mu_{2},jet) vs #Delta R(#mu_{2},jet)"+post;
  t+="; #Delta R(#mu_{2},jet); p_{T,rel}(#mu_{2},jet)"; ptrelVsdRminMu2jet=td->make<TH2D>("ptrelVsdRminMu2jet",
											 t.c_str(),
											 50,0,5.,50,0,1000);
  
  // ----------  Composite histograms  ----------

  t="M(W_{R}) "                    +post;       mWR=td->make<TH1D>("mWR",   t.c_str(),50,0,2500);
  t="M(N_{R}) with #mu_{1} "       +post;     mNuR1=td->make<TH1D>("mNuR1", t.c_str(),50,0,2500);
  t="M(N_{R}) with #mu_{2} "       +post;     mNuR2=td->make<TH1D>("mNuR2", t.c_str(),50,0,1000);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post;    mNuR2D=td->make<TH2D>("mNuR2D",t.c_str(),50,0,2500,50,0,1000);

  t="M(#mu #mu)"                   +post;     mMuMu=td->make<TH1D>("mMuMu",    t.c_str(),50,0,2500);
  t="M(#mu #mu)"                   +post; mMuMuZoom=td->make<TH1D>("mMuMuZoom",t.c_str(),50,0,200);
  t="M(jj)"                        +post;       mJJ=td->make<TH1D>("mJJ",      t.c_str(),50,0,2500);

  t="cZeta(mumu)"                 +post; czeta_mumu=td->make<TH1D>("czMM",    t.c_str(),100,-1,1);
  t="cZeta(mumu) Zoom"       +post; czeta_mumu_zoom=td->make<TH1D>("czMMzoom",t.c_str(),100,-1,-0.9);

  // crazy angles
  t="cT(mumu)"    +post;     ctheta_mumu=td->make<TH1D>("ctMM",   t.c_str(),50,0,1);
  t="cT(jj)"      +post;       ctheta_jj=td->make<TH1D>("ctJJ",   t.c_str(),50,0,1);
  t="cT(mu1-jj)"  +post;   ctheta_mu1_jj=td->make<TH1D>("ctM1JJ", t.c_str(),50,0,1);
  t="cT(mu2-jj)"  +post;   ctheta_mu2_jj=td->make<TH1D>("ctM2JJ", t.c_str(),50,0,1);

  t="cTz(mumu)"   +post;    cthetaz_mumu=td->make<TH1D>("ctzMM",  t.c_str(),50,0,1);
  t="cTz(jj)"     +post;      cthetaz_jj=td->make<TH1D>("ctzJJ",  t.c_str(),50,0,1);
  t="cTz(mu1-jj)" +post;  cthetaz_mu1_jj=td->make<TH1D>("ctzM1JJ",t.c_str(),50,0,1);
  t="cTz(mu2-jj)" +post;  cthetaz_mu2_jj=td->make<TH1D>("ctzM2JJ",t.c_str(),50,0,1);

  // ----------  Neural Net histograms  ----------

  if (v_masspts.size()) {
    nndir= new TFileDirectory(td->mkdir("_NNdata"));
    for (size_t i=0; i<v_masspts.size(); i++) {
      int mwr = v_masspts[i].first;
      int mnu = v_masspts[i].second;
      std::string name = nnhistoname(mwr,mnu);
      nndir->make<TH1D>(name.c_str(),(name+post).c_str(),51,-0.01,1.01);
    }
  }
}// end of book()

//======================================================================

void HeavyNu::HistPerDef::fill(pat::MuonCollection muons,
			       pat::JetCollection  jets,
			       bool isMC)
{  
  std::sort(muons.begin(),muons.end(),compare()) ; 
  std::sort(jets.begin(),jets.end(),compare()) ; 

  reco::Particle::LorentzVector WR ; 

  // Muons 
  ptMu1->Fill(muons.at(0).pt()) ; 
  ptMu2->Fill(muons.at(1).pt()) ; 

  etaMu1->Fill(muons.at(0).eta()) ; 
  etaMu2->Fill(muons.at(1).eta()) ; 

  phiMu1->Fill(muons.at(0).phi()) ; 
  phiMu2->Fill(muons.at(1).phi()) ; 

  dPhiMu->Fill( fabs(deltaPhi(muons.at(0).phi(),muons.at(1).phi())) ) ; 
  dEtaMu->Fill( fabs(muons.at(0).eta() - muons.at(1).eta()) ) ; 
  dEtaPhiMu->Fill(fabs(muons.at(0).eta()-muons.at(1).eta()),
		  fabs(deltaPhi(muons.at(0).phi(),muons.at(1).phi()))) ; 
  
  mu1trackIso->Fill(muons.at(0).trackIso());
  mu1hcalIso ->Fill(muons.at(0).hcalIso());
  mu1ecalIso ->Fill(muons.at(0).ecalIso());
  mu1caloIso ->Fill(muons.at(0).caloIso());
  mu1dB      ->Fill(muons.at(0).dB());
  mu2trackIso->Fill(muons.at(1).trackIso());
  mu2hcalIso ->Fill(muons.at(1).hcalIso());
  mu2ecalIso ->Fill(muons.at(1).ecalIso());
  mu2caloIso ->Fill(muons.at(1).caloIso());
  mu2dB      ->Fill(muons.at(1).dB());
  
  mu1trackRelIso->Fill(muons.at(0).trackIso()/muons.at(0).pt());
  mu1hcalRelIso ->Fill(muons.at(0).hcalIso() /muons.at(0).pt());
  mu1ecalRelIso ->Fill(muons.at(0).ecalIso() /muons.at(0).pt());
  mu1caloRelIso ->Fill(muons.at(0).caloIso() /muons.at(0).pt());
  mu2trackRelIso->Fill(muons.at(1).trackIso()/muons.at(1).pt());
  mu2hcalRelIso ->Fill(muons.at(1).hcalIso() /muons.at(1).pt());
  mu2ecalRelIso ->Fill(muons.at(1).ecalIso() /muons.at(1).pt());
  mu2caloRelIso ->Fill(muons.at(1).caloIso() /muons.at(1).pt());

  if (isMC) {
    for (unsigned int i=0; i<2; i++) { 
      if ( muons.at(i).genLepton() != 0 ) {
	float dpt = muons.at(i).pt()-muons.at(i).genLepton()->pt() ; 
	float dR = deltaR(muons.at(i).eta(),muons.at(i).phi(),
			  muons.at(i).genLepton()->eta(),muons.at(i).genLepton()->phi()) ; 
	if ( i == 0 ) { 
	  dptMu1gen->Fill(dpt/muons.at(i).genLepton()->pt()) ; 
	  dRMu1gen->Fill(dR) ; 
	} else { 
	  dptMu2gen->Fill(dpt/muons.at(i).genLepton()->pt()) ; 
	  dRMu2gen->Fill(dR) ; 
	}
      }
    }
  }
  for (int i=0; i<muonQualityFlags; i++) { 
    if (muons.at(0).muonID(muonQuality[i])) qualMu1->Fill( i ) ; 
    if (muons.at(1).muonID(muonQuality[i])) qualMu2->Fill( i ) ; 
  }

  // Jets 
  ptJet1->Fill(jets.at(0).pt()) ; 
  ptJet2->Fill(jets.at(1).pt()) ; 

  etaJet1->Fill(jets.at(0).eta()) ; 
  etaJet2->Fill(jets.at(1).eta()) ; 

  phiJet1->Fill(jets.at(0).phi()) ; 
  phiJet2->Fill(jets.at(1).phi()) ; 

  btagJet1->Fill(jets.at(0).bDiscriminator(btagName));
  btagJet1->Fill(jets.at(1).bDiscriminator(btagName));

  dPhiJet->Fill( fabs(deltaPhi(jets.at(0).phi(),jets.at(1).phi())) ) ; 
  dEtaJet->Fill( fabs(jets.at(0).eta() - jets.at(1).eta()) ) ; 
  dEtaPhiJet->Fill(fabs(jets.at(0).eta()-jets.at(1).eta()),
		   fabs(deltaPhi(jets.at(0).phi(),jets.at(1).phi()))) ;

  jetID2d->Fill(jetID(jets.at(0)),jetID(jets.at(1)));

  // Muon-Jet plots
  float dRmu1jet1 = deltaR(muons.at(0).eta(),muons.at(0).phi(),jets.at(0).eta(),jets.at(0).phi()) ; 
  float dRmu1jet2 = deltaR(muons.at(0).eta(),muons.at(0).phi(),jets.at(1).eta(),jets.at(1).phi()) ; 
  float dRmu2jet1 = deltaR(muons.at(1).eta(),muons.at(1).phi(),jets.at(0).eta(),jets.at(0).phi()) ; 
  float dRmu2jet2 = deltaR(muons.at(1).eta(),muons.at(1).phi(),jets.at(1).eta(),jets.at(1).phi()) ; 

  const pat::Jet&  j4mu1 = (dRmu1jet1 < dRmu1jet2) ? jets.at(0) : jets.at(1);
  const pat::Jet&  j4mu2 = (dRmu2jet1 < dRmu2jet2) ? jets.at(0) : jets.at(1);

  TVector3 mu1vec(muons.at(0).momentum().X(), muons.at(0).momentum().Y(), muons.at(0).momentum().Z());
  TVector3 mu2vec(muons.at(1).momentum().X(), muons.at(1).momentum().Y(), muons.at(1).momentum().Z());

  TVector3 jt1vec(j4mu1.p4().Vect().X(), j4mu1.p4().Vect().Y(), j4mu1.p4().Vect().Z() );
  TVector3 jt2vec(j4mu2.p4().Vect().X(), j4mu2.p4().Vect().Y(), j4mu2.p4().Vect().Z() );

  double ptrelMu1 = mu1vec.Perp(jt1vec);
  double ptrelMu2 = mu2vec.Perp(jt2vec);

  dRminMu1jet->Fill( min(dRmu1jet1,dRmu1jet2) ) ; 
  dRminMu2jet->Fill( min(dRmu2jet1,dRmu2jet2) ) ; 

  hptrelMu1->Fill( ptrelMu1 );
  hptrelMu2->Fill( ptrelMu2 );

  ptrelVsdRminMu1jet->Fill(min(dRmu1jet1,dRmu1jet2),ptrelMu1);
  ptrelVsdRminMu2jet->Fill(min(dRmu2jet1,dRmu1jet2),ptrelMu2);

  // Composite objects
  reco::Particle::LorentzVector vWR = jets.at(0).p4() + jets.at(1).p4() ; 
  WR = vWR + muons.at(0).p4() + muons.at(1).p4() ; 

  mWR->Fill(WR.M()) ; 
  mNuR1->Fill( (vWR + muons.at(0).p4()).M() ) ; 
  mNuR2->Fill( (vWR + muons.at(1).p4()).M() ) ; 
  mNuR2D->Fill( (vWR + muons.at(0).p4()).M(),(vWR + muons.at(1).p4()).M() ) ; 

  reco::Particle::LorentzVector mumu=muons.at(0).p4()+muons.at(1).p4();
  reco::Particle::LorentzVector jj=jets.at(0).p4()+jets.at(1).p4();

  mMuMu->Fill(mumu.M() );
  mMuMuZoom->Fill(mumu.M() );
  mJJ->Fill(jj.M() );

}// end of fill()

//======================================================================

void
HeavyNu::HistPerDef::fill(const HeavyNuEvent& hne,
			  const std::vector<hNuMassHypothesis>& v_masspts)
{
  // Muons 
  ptMu1->Fill(hne.mu1->pt()) ; 
  ptMu2->Fill(hne.mu2->pt()) ; 

  etaMu1->Fill(hne.mu1->eta()) ; 
  etaMu2->Fill(hne.mu2->eta()) ; 

  phiMu1->Fill(hne.mu1->phi()) ; 
  phiMu2->Fill(hne.mu2->phi()) ; 

  dPhiMu->Fill( fabs(deltaPhi(hne.mu1->phi(),hne.mu2->phi())) ) ; 
  dEtaMu->Fill( fabs(hne.mu1->eta() - hne.mu2->eta()) ) ; 
  dEtaPhiMu->Fill(fabs(hne.mu1->eta()-hne.mu2->eta()),
		  fabs(deltaPhi(hne.mu1->phi(),hne.mu2->phi()))) ; 

  mu1trackIso->Fill(hne.mu1->trackIso());
  mu1hcalIso ->Fill(hne.mu1->hcalIso());
  mu1ecalIso ->Fill(hne.mu1->ecalIso());
  mu1caloIso ->Fill(hne.mu1->caloIso());
  mu1dB      ->Fill(hne.mu1->dB());
  mu2trackIso->Fill(hne.mu2->trackIso());
  mu2hcalIso ->Fill(hne.mu2->hcalIso());
  mu2ecalIso ->Fill(hne.mu2->ecalIso());
  mu2caloIso ->Fill(hne.mu2->caloIso());
  mu2dB      ->Fill(hne.mu2->dB());
  
  mu1trackRelIso->Fill(hne.mu1->trackIso()/hne.mu1->pt());
  mu1hcalRelIso ->Fill(hne.mu1->hcalIso() /hne.mu1->pt());
  mu1ecalRelIso ->Fill(hne.mu1->ecalIso() /hne.mu1->pt());
  mu1caloRelIso ->Fill(hne.mu1->caloIso() /hne.mu1->pt());
  mu2trackRelIso->Fill(hne.mu2->trackIso()/hne.mu2->pt());
  mu2hcalRelIso ->Fill(hne.mu2->hcalIso() /hne.mu2->pt());
  mu2ecalRelIso ->Fill(hne.mu2->ecalIso() /hne.mu2->pt());
  mu2caloRelIso ->Fill(hne.mu2->caloIso() /hne.mu2->pt());
  
  if (hne.isMC) {
    for (unsigned int i=0; i<2; i++) { 
      if ( hne.mu[i]->genLepton() != 0 ) {
	float dpt = hne.mu[i]->pt()-hne.mu[i]->genLepton()->pt() ; 
	float dR  = deltaR(hne.mu[i]->eta(),hne.mu[i]->phi(),
			   hne.mu[i]->genLepton()->eta(),hne.mu[i]->genLepton()->phi()) ; 
	if ( i == 0 ) { 
	  dptMu1gen->Fill(dpt/hne.mu[i]->genLepton()->pt()) ; 
	  dRMu1gen->Fill(dR) ; 
	} else { 
	  dptMu2gen->Fill(dpt/hne.mu[i]->genLepton()->pt()) ; 
	  dRMu2gen->Fill(dR) ; 
	}
      }
    }
  }
  for (int i=0; i<muonQualityFlags; i++) { 
    if (hne.mu1->muonID(muonQuality[i])) qualMu1->Fill( i ) ; 
    if (hne.mu2->muonID(muonQuality[i])) qualMu2->Fill( i ) ; 
  }

  int jet1id = 0;
  int jet2id = 0;

  // Jets 
  if (hne.j1.isAvailable()) {
    jet1id = jetID(*(hne.j1));

    ptJet1->Fill(hne.j1->pt()) ; 
    etaJet1->Fill(hne.j1->eta()) ; 
    phiJet1->Fill(hne.j1->phi()) ; 
    btagJet1->Fill(hne.j1->bDiscriminator(btagName));

    if (hne.j2.isAvailable()) {
      ptJet2->Fill(hne.j2->pt()) ; 
      etaJet2->Fill(hne.j2->eta()) ; 
      phiJet2->Fill(hne.j2->phi()) ; 
      btagJet2->Fill(hne.j2->bDiscriminator(btagName));

      jet2id = jetID(*(hne.j2));

      dPhiJet->Fill( fabs(deltaPhi(hne.j1->phi(),hne.j2->phi())) ) ; 
      dEtaJet->Fill( fabs(hne.j1->eta() - hne.j2->eta()) ) ; 
      dEtaPhiJet->Fill(fabs(hne.j1->eta()-hne.j2->eta()),
		       fabs(deltaPhi(hne.j1->phi(),hne.j2->phi()))) ;

      mWR->Fill   ( hne.mWR   ) ; 
      mNuR1->Fill ( hne.mNuR1 ) ; 
      mNuR2->Fill ( hne.mNuR2 ) ; 
      mNuR2D->Fill( hne.mNuR1, hne.mNuR2 );
      mJJ->Fill   ( hne.mJJ   );

      ctheta_jj->Fill(hne.ctheta_jj);
      ctheta_mu1_jj->Fill(hne.ctheta_mu1_jj);
      ctheta_mu2_jj->Fill(hne.ctheta_mu2_jj);
      cthetaz_jj->Fill(hne.cthetaz_jj);
      cthetaz_mu1_jj->Fill(hne.cthetaz_mu1_jj);
      cthetaz_mu2_jj->Fill(hne.cthetaz_mu2_jj);
    }

    dRminMu1jet->Fill(hne.dRminMu1jet);
    dRminMu2jet->Fill(hne.dRminMu2jet);

    hptrelMu1->Fill(hne.ptrelMu1);
    hptrelMu2->Fill(hne.ptrelMu2);

    ptrelVsdRminMu1jet->Fill(hne.dRminMu1jet,hne.ptrelMu1);
    ptrelVsdRminMu2jet->Fill(hne.dRminMu2jet,hne.ptrelMu2);
  }

  jetID2d->Fill(jet1id,jet2id);

  mMuMu->Fill( hne.mMuMu );
  mMuMuZoom->Fill( hne.mMuMu );

  czeta_mumu->Fill(hne.czeta_mumu);
  czeta_mumu_zoom->Fill(hne.czeta_mumu);
  ctheta_mumu->Fill(hne.ctheta_mumu);
  cthetaz_mumu->Fill(hne.cthetaz_mumu);

  // Neural net histos
  if (v_masspts.size()) {
    TDirectory *nnrootdir = nndir->cd();
    for (size_t i=0; i<v_masspts.size(); i++) {
      int mwr = v_masspts[i].first;
      int mnu = v_masspts[i].second;
      std::string name = nnhistoname(mwr,mnu);
      TH1D *nnh = (TH1D *)nnrootdir->Get(name.c_str());
      assert(nnh);
      nnh->Fill(hne.nnoutputs[i]);
    }
  }

}// end of fill()

//======================================================================

//
// constants, enums and typedefs
//
const std::vector<hNuMassHypothesis> v_null;

//
// static data member definitions
//

//
// constructors and destructor
//
HeavyNu::HeavyNu(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  dolog_=iConfig.getParameter<bool>("DoLog");

  muonTag_ = iConfig.getParameter< edm::InputTag >( "muonTag" );
  jetTag_  = iConfig.getParameter< edm::InputTag >( "jetTag"  );
  elecTag_ = iConfig.getParameter< edm::InputTag >( "electronTag" );

  btagName=iConfig.getParameter<std::string>("BtagName");

  nnif_ = new HeavyNu_NNIF(iConfig);

  trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet>("trigMatchPset"));

  edm::Service<TFileService> fs;
  hists.nelec    = fs->make<TH1D>("nelec", "N(e^{#pm})",10,-0.5,9.5);
  hists.nmu      = fs->make<TH1D>("nmu",   "N(#mu^{#pm})",10,-0.5,9.5);
  hists.njet     = fs->make<TH1D>("njet",  "N(Jet)",50,-0.5,49.5);
  hists.muPt     = fs->make<TH1D>("muPt",  "#mu p_{T} distribution",100,0,2000);
  hists.muEta    = fs->make<TH1D>("muEta", "#mu #eta distribution",50,-2.5,2.5);
  hists.muPhi    = fs->make<TH1D>("muPhi", "#mu #phi distribution",60,-3.14159,3.14159);
  hists.jetPt    = fs->make<TH1D>("jetPt", "jet p_{T} distribution",100,0,2000);
  hists.jetEta   = fs->make<TH1D>("jetEta","jet #eta distribution",50,-5,5);
  hists.jetPhi   = fs->make<TH1D>("jetPhi","jet #phi distribution",60,-3.14159,3.14159);
  hists.jetID    = fs->make<TH1I>("jetID", "Jet ID",3,0,3);
  hists.jetPtvsNum=fs->make<TH2D>("jetPtvsNum","Jet P_{T} vs. Jet # ",11,-0.5,10.5,200,0.,2000.);
  labelJetIDaxis(hists.jetID->GetXaxis());

  hists.noCuts.book      ( new TFileDirectory(fs->mkdir("cut0_none")),    "(no cuts)",                v_null );
  hists.LLptCuts.book    ( new TFileDirectory(fs->mkdir("cut1_LLpt")),    "(dileptons with ptcuts:1)",v_null );
  hists.MuTightCuts.book ( new TFileDirectory(fs->mkdir("cut2_MuTight")), "(Mu tight cuts:2)",        v_null );
  hists.LLJJptCuts.book  ( new TFileDirectory(fs->mkdir("cut4_LLJJpt")),  "(4objects with ptcuts:4)", nnif_->masspts() );
  hists.diLmassCut.book  ( new TFileDirectory(fs->mkdir("cut5_diLmass")), "(mumu mass cut:5)",        nnif_->masspts() );
  hists.mWRmassCut.book  ( new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(mumujj mass cut:6)",      nnif_->masspts() );

  if (trig_->matchingEnabled()) {
    hists.TrigMatches.book            ( new TFileDirectory(fs->mkdir("cut3_TrigMatches")),       "(Trigger match:3)",               v_null );

    hists.MuTightInZwin.book          ( new TFileDirectory(fs->mkdir("MuTightInZwin")),          "(Mu1 tight in Z mass Window)",    v_null );
    hists.Mu1TrigMatchesInZwin.book   ( new TFileDirectory(fs->mkdir("Mu1TrigMatchesInZwin")),   "(#mu1 trigger match in Z mass Window)",v_null );
    hists.Mu2TrigMatchesInZwin.book   ( new TFileDirectory(fs->mkdir("Mu2TrigMatchesInZwin")),   "(#mu2 Trigger match in Z mass Window)",v_null );
    hists.Mu1Mu2TrigMatchesInZwin.book( new TFileDirectory(fs->mkdir("Mu1Mu2TrigMatchesInZwin")),"(#mu1,#mu2 Trigger match in Z mass Window)",v_null );
    trig_->book(*(hists.Mu1TrigMatchesInZwin.mydir), &(hists.Mu1TrigMatchesInZwin.trigHistos));
    trig_->book(*(hists.Mu2TrigMatchesInZwin.mydir), &(hists.Mu2TrigMatchesInZwin.trigHistos));
  }

  hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

  init_=false;

  cuts.minimum_mu1_pt      = iConfig.getParameter<double>("minMu1pt");
  cuts.minimum_mu2_pt      = iConfig.getParameter<double>("minMu2pt");
  cuts.minimum_jet_pt      = iConfig.getParameter<double>("minJetPt");
  cuts.maximum_mu_abseta   = iConfig.getParameter<double>("maxMuAbsEta");
  cuts.maximum_jet_abseta  = iConfig.getParameter<double>("maxJetAbsEta");
  cuts.minimum_mumu_mass   = iConfig.getParameter<double>("minMuMuMass");
  cuts.minimum_mWR_mass    = iConfig.getParameter<double>("min4objMass");
  cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");
  cuts.muon_trackiso_limit = iConfig.getParameter<double>("muonTrackIsoLimitGeV");

  ZwinMinGeV_ = iConfig.getParameter<double>("ZmassWinMinGeV");
  ZwinMaxGeV_ = iConfig.getParameter<double>("ZmassWinMaxGeV");

  // For the record...
  std::cout << "Configurable cut values applied:" << std::endl;
  std::cout << "muonTag        = " << muonTag_                 << std::endl;
  std::cout << "jetTag         = " << jetTag_                  << std::endl;
  std::cout << "electronTag    = " << elecTag_                 << std::endl;
  std::cout << "ZmassWinMinGeV = " << ZwinMinGeV_              << " GeV" << std::endl;
  std::cout << "ZmassWinMaxGeV = " << ZwinMaxGeV_              << " GeV" << std::endl;
  std::cout << "minMu1pt       = " << cuts.minimum_mu1_pt      << " GeV" << std::endl;
  std::cout << "minMu2pt       = " << cuts.minimum_mu2_pt      << " GeV" << std::endl;
  std::cout << "minJetPt       = " << cuts.minimum_jet_pt      << " GeV" << std::endl;
  std::cout << "maxMuAbsEta    = " << cuts.maximum_mu_abseta   << std::endl;
  std::cout << "maxJetAbsEta   = " << cuts.maximum_jet_abseta  << std::endl;
  std::cout << "minMuonJetdR   = " << cuts.minimum_muon_jet_dR << std::endl;
  std::cout << "muonTrackIso   = " << cuts.muon_trackiso_limit << " GeV" << std::endl;
  std::cout << "minMuMuMass    = " << cuts.minimum_mumu_mass   << " GeV" << std::endl;
  std::cout << "min4objMass    = " << cuts.minimum_mWR_mass    << " GeV" << std::endl;
}
  
HeavyNu::~HeavyNu()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool
HeavyNu::isVBTFloose(const pat::Muon& m)
{
  return (m.muonID("AllGlobalMuons") &&
	  (m.numberOfValidHits() > 10));
}

bool
HeavyNu::isVBTFtight(const pat::Muon& m)
{
  if (!isVBTFloose(m)) return false;

  reco::TrackRef gt = m.globalTrack();
  if (gt.isNull()) {
    std::cerr << "Mu global track reference is NULL" << std::endl;
    return false;
  }
  return (m.muonID("AllTrackerMuons") &&
	  (m.dB() < 0.2) &&
	  (m.normChi2() < 10) &&
	  (m.numberOfMatches() > 1) &&
	  (gt->hitPattern().numberOfValidMuonHits()>0) &&
	  (gt->hitPattern().numberOfValidPixelHits()>0) );
}

TH1 *
HeavyNu::bookRunHisto(uint32_t runNumber)
{
  std::string runstr = int2str<uint32_t>(runNumber);
  return hists.rundir->make <TH1I> (runstr.c_str(), runstr.c_str(),1,1,2);
}

// ------------ method called to for each event  ------------
bool
HeavyNu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  HeavyNuEvent hnuEvent;


  hnuEvent.isMC = !iEvent.isRealData();

  if (iEvent.isRealData())
  {
    uint32_t runn = iEvent.id().run();
    std::map<uint32_t,TH1 *>::const_iterator it = m_runHistos_.find(runn);
    TH1 *runh;
    if (it == m_runHistos_.end()) {
      runh = bookRunHisto(runn);
      m_runHistos_[runn] = runh;
    } else
      runh = it->second;
    runh->Fill(1);
  }

  edm::Handle<pat::MuonCollection> pMuons ; 
  iEvent.getByLabel(muonTag_,pMuons) ; 

  edm::Handle<pat::ElectronCollection> pElecs ;
  iEvent.getByLabel(elecTag_, pElecs) ;

  edm::Handle<pat::JetCollection> pJets ;
  iEvent.getByLabel(jetTag_, pJets) ;

  if ( !pElecs.isValid() || 
       !pMuons.isValid() || 
       !pJets.isValid() ) {
    std::cout << "Exiting as valid PAT objects not found" << std::endl ;
    return false; 
  }

  pat::JetRef  iJ;
  pat::MuonRef iM;

  hists.nelec->Fill(pElecs->size()) ;
  hists.nmu->Fill(pMuons->size()) ;
  hists.njet->Fill(pJets->size()) ;

  int njet=1;
  for (size_t iJet=0; iJet<pJets->size(); iJet++) { 
    iJ=pat::JetRef(pJets,iJet);
    hists.jetPt->Fill( (*iJ).pt() ) ; 
    hists.jetEta->Fill( (*iJ).eta() ) ; 
    hists.jetPhi->Fill( (*iJ).phi() ) ; 
    hists.jetPtvsNum->Fill(njet++, (*iJ).pt() ) ; 
    hists.jetID->Fill(jetID(*iJ));
  }
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) { 
    iM=pat::MuonRef(pMuons,iMuon);
    hists.muPt->Fill( (*iM).pt() ) ; 
    hists.muEta->Fill( (*iM).eta() ) ; 
    hists.muPhi->Fill( (*iM).phi() ) ; 
  }

  // Basic selection requirements: Require at least two muons, two jets
  if ( pMuons->size() >= 2 && pJets->size() >= 2 ) {
    hists.noCuts.fill( *pMuons,*pJets, hnuEvent.isMC ) ; 
  } else return false;

  // next, we look for valid muons and jets and put them into the Event
  for (size_t iJet=0; iJet<pJets->size(); iJet++) {
    iJ=pat::JetRef(pJets,iJet);
    float jpt = (*iJ).pt();
    if ((jpt > cuts.minimum_jet_pt)   && // more later!
	(fabs((*iJ).eta())<=cuts.maximum_jet_abseta) ) {
      if        ((hnuEvent.j1.isNull()) || (hnuEvent.j1->pt() < jpt)) {
	hnuEvent.j2=hnuEvent.j1;
	hnuEvent.j1=iJ;
      } else if ((hnuEvent.j2.isNull()) || (hnuEvent.j2->pt() < jpt)) {
	hnuEvent.j2=iJ;
      }
    }
  }

  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    iM=pat::MuonRef(pMuons,iMuon);
    double dr1=(hnuEvent.j1.isNull())?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hnuEvent.j1->eta(),hnuEvent.j1->phi()));
    double dr2=(hnuEvent.j2.isNull())?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hnuEvent.j2->eta(),hnuEvent.j2->phi()));

    if (((*iM).pt()>cuts.minimum_mu2_pt)
	&& isVBTFloose(*iM)
	&& (fabs((*iM).eta()) < cuts.maximum_mu_abseta)
	&& (std::min(dr1,dr2) > cuts.minimum_muon_jet_dR)
	&& ((*iM).trackIso()  < cuts.muon_trackiso_limit)
	) {
      if ( (hnuEvent.mu1.isNull()) ||
	   (hnuEvent.mu1->pt()<(*iM).pt()) ) {
	hnuEvent.mu2=hnuEvent.mu1;
	hnuEvent.mu1=iM;
      } else 	if (hnuEvent.mu2.isNull() ||
		    hnuEvent.mu2->pt()<(*iM).pt()) {
	hnuEvent.mu2=iM;
      }
    }
  }

  // require two "loose" muons first
  if (hnuEvent.mu2.isNull())
    return false;

  hnuEvent.regularize(); // assign internal standards
  hnuEvent.calculateMuMu();
  hists.LLptCuts.fill(hnuEvent,v_null);
  
  // Require mu1.OR.mu2 meets tight requirements
  //
  bool mu1isTight = isVBTFtight(*(hnuEvent.mu1));
  bool mu2isTight = isVBTFtight(*(hnuEvent.mu2));

  if ( !mu1isTight && !mu2isTight )
    return false;

  hists.MuTightCuts.fill( hnuEvent,v_null );

  // split out trigger matching requirement to study trigger eff.
  if ( trig_->matchingEnabled() ) {
    if ( inZmassWindow( hnuEvent.mMuMu ) )
      hists.MuTightInZwin.fill( hnuEvent,v_null );

    // require that one muon be BOTH tight and trigger-matched
    //
    bool mu1matches = mu1isTight &&
      trig_->isTriggerMatched( hnuEvent.mu1, iEvent,
			       &(hists.Mu1TrigMatchesInZwin.trigHistos));

    bool mu2matches = mu2isTight &&
      trig_->isTriggerMatched( hnuEvent.mu2, iEvent,
			       &(hists.Mu2TrigMatchesInZwin.trigHistos));

    if ( !mu1matches && !mu2matches )
      return false;

    hists.TrigMatches.fill( hnuEvent,v_null );

    // histos for trigger efficiency study:
    //
    if ( inZmassWindow( hnuEvent.mMuMu ) ) {
      if ( mu1matches ) {
	hists.Mu1TrigMatchesInZwin.fill     ( hnuEvent,v_null );
	if ( mu2matches ) {
	  hists.Mu2TrigMatchesInZwin.fill   ( hnuEvent,v_null );
	  hists.Mu1Mu2TrigMatchesInZwin.fill( hnuEvent,v_null );
	}
      }
      else // mu2matches
	hists.Mu2TrigMatchesInZwin.fill     ( hnuEvent,v_null );
    }
  }
  
  // require four objects
  // require also the selected jets to pass loose ID,
  // per JetMET recommendation
  //
  if ( (        hnuEvent.mu2.isNull()) ||
       (        hnuEvent.j2.isNull())  ||
       (jetID(*(hnuEvent.j1)) < 1)     || 
       (jetID(*(hnuEvent.j2)) < 1) ) return false;

  hnuEvent.calculate(); // calculate various details

  //dumpJetCorInfo( *(hnuEvent.j1) );

  nnif_->fillvector( hnuEvent );
  nnif_->output( hnuEvent.nnoutputs );

  hists.LLJJptCuts.fill( hnuEvent,nnif_->masspts() );

  if ( hnuEvent.mMuMu<cuts.minimum_mumu_mass ) return false;  // dimuon mass cut
  hists.diLmassCut.fill( hnuEvent,nnif_->masspts() );

  if ( iEvent.isRealData() ) {
    std::cout<<"\t"<<iEvent.id() << std::endl;
    std::cout<<"\tM(W_R)  = "<<hnuEvent.mWR  <<" GeV";
    std::cout<<", M(NuR1) = "<<hnuEvent.mNuR1<<" GeV";
    std::cout<<", M(NuR2) = "<<hnuEvent.mNuR2<<" GeV"<<std::endl;
    std::cout<<"\tJets:   j1 ";outputCandidate( reco::CandidateBaseRef( hnuEvent.j1  ) );
    std::cout<<        ", j2 ";outputCandidate( reco::CandidateBaseRef( hnuEvent.j2  ) ); std::cout<<std::endl;
    std::cout<<"\tMuons: mu1 ";outputCandidate( reco::CandidateBaseRef( hnuEvent.mu1 ) );
    std::cout<<       ", mu2 ";outputCandidate( reco::CandidateBaseRef( hnuEvent.mu2 ) ); std::cout<<std::endl;
  }

  if ( hnuEvent.mWR<cuts.minimum_mWR_mass ) return false;  // 4-object mass cut
  hists.mWRmassCut.fill( hnuEvent,nnif_->masspts() );

  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNu::beginJob() {
  nnif_->beginJob();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNu::endJob() {
  nnif_->endJob();
  trig_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);