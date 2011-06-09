// -*- C++ -*-
//
// Package:    HFZeeVBTF
// Class:      HFZeeVBTF
// 
/**\class HFZeeVBTF HFZeeVBTF.cc Work/HFZeeVBTF/src/HFZeeVBTF.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HFZeeVBTF.cc,v 1.23 2011/04/15 01:45:06 mansj Exp $
//
//

// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ZShape/HFZeeVBTF/interface/SelectElectron.h"


//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}




class HFZeeVBTF : public edm::EDFilter {
public:
  explicit HFZeeVBTF(const edm::ParameterSet&);
  ~HFZeeVBTF();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }

  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string ecalID_, currentFile_;
  bool dolog_;
  std::string myName_;
  double minEtECAL_;
  double minEtHF_;
  std::vector<int> acceptedElectronIDs_;
  std::vector<double> hfIdParams_;
  edm::InputTag hFElectrons_;
  selectElectron w;
  std::string filterId_;
  std::vector<double> massWindow_;

  // ----------member data ---------------------------

  // gf set of histograms per each == Z definition ==
  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td,const std::string&,const std::vector<double>& win);
    // fill all histos of the set with the two electron candidates
    void fill(pat::ElectronCollection::const_iterator ecalE,  const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, 
	      const bool hasEleIDPassed, const edm::ParameterSet& myPs, const int isTight, std::vector<double> hfIdParams);
    TH1* mee, *meeHFP, *meeHFM;
    TH2* mee_vsEta, *mee_vsAbsEta, *meeHFP_vsEta, *meeHFM_vsEta;
    TH1* Yee, *Ptee;
    TH1* ec_eta,       *ec_phi,     *ec_pt;
    TH1* hf_eta,       *hf_phi,     *hf_pt;
    TH1* hfp_eta,      *hfp_phi,    *hfp_pt;
    TH1* hfm_eta,      *hfm_phi,    *hfm_pt;
    TH1* hf_e9e25,     *hf_eCOREe9, *hf_eSeL, *hf_var2d;
    TH1* hf_seed_eSeL, *hf_e1e9;
    TH1* hf_e,         *hf_e1x1,    *hf_core, *hf_e3x3, *hf_e5x5;
    TH1* el_1_eta,     *el_2_eta;
    TH2* el_1_eta_VS_el_2_eta, *el_1_phi_VS_el_2_phi;
    // gf: add hisos for variables N-1 selected
    TH1* combIsoEB_nmo,  *relTkIsoEB_nmo, *relEcIsoEB_nmo, *relHcIsoEB_nmo, *sigiEtaiEtaEB_nmo, *dphiEB_nmo, *detaEB_nmo, *hOeEB_nmo;
    TH1* combIsoEE_nmo,  *relTkIsoEE_nmo, *relEcIsoEE_nmo, *relHcIsoEE_nmo, *sigiEtaiEtaEE_nmo, *dphiEE_nmo, *detaEE_nmo, *hOeEE_nmo;
    TH1* e9e25_nmo, *var2d_nmo, *eCOREe9_nmo, *eSeL_nmo;
    TH2* eSeL_vs_logEl;
    std::vector<double> massWindow_;

  };
  std::vector<std::string> HLT_Names;
  bool init_,skipTrigger_;

  // gf set of histo for all Z definitios in a stack
  struct HistStruct {
    TH1 *nelec,*nhf;
    HistPerDef base, mee90none, zMass, ;
    HistPerDef mee90loose, mee80loose, mee70loose, mee60loose;
    HistPerDef mee90tight, mee80tight, mee70tight, mee60tight;
  } hists;


  edm::ParameterSet robust95relIsoEleIDCutsV04_ps_;
  edm::ParameterSet robust90relIsoEleIDCutsV04_ps_;
  edm::ParameterSet robust85relIsoEleIDCutsV04_ps_;
  edm::ParameterSet robust80relIsoEleIDCutsV04_ps_;
  edm::ParameterSet robust70relIsoEleIDCutsV04_ps_;
  edm::ParameterSet robust60relIsoEleIDCutsV04_ps_;

  edm::ParameterSet robust95cIsoEleIDCutsV04_ps_  ;
  edm::ParameterSet robust90cIsoEleIDCutsV04_ps_  ;
  edm::ParameterSet robust85cIsoEleIDCutsV04_ps_  ;
  edm::ParameterSet robust80cIsoEleIDCutsV04_ps_  ;
  edm::ParameterSet robust70cIsoEleIDCutsV04_ps_  ;
  edm::ParameterSet robust60cIsoEleIDCutsV04_ps_  ;

};


void HFZeeVBTF::HistPerDef::book(TFileDirectory td, const std::string& post, const std::vector<double>& win) {
  massWindow_=win;
  std::string title;
  double minZrange = 40;
  double maxZrange = 120;
  
  title=std::string("M_{ee} ")+post;
  mee=td.make<TH1D>("mee",title.c_str(),120,minZrange,maxZrange);  
  title=std::string("M_{ee,HF+} ")+post;
  meeHFP=td.make<TH1D>("mee-HFP",title.c_str(),120,minZrange,maxZrange);  
  title=std::string("M_{ee,HF-} ")+post;
  meeHFM=td.make<TH1D>("mee-HFM",title.c_str(),120,minZrange,maxZrange);  

  title=std::string("M_{ee} vs eta")+post;
  mee_vsEta=td.make<TH2D>("mee_vsEta",title.c_str(),10,3,5,30,minZrange,maxZrange);
  title=std::string("M_{ee,HF+} vs abs(eta) ")+post;
  mee_vsAbsEta=td.make<TH2D>("mee_vsAbsEta",title.c_str(),10,3,5,30,minZrange,maxZrange);  
  title=std::string("M_{ee,HF+} vs eta ")+post;
  meeHFP_vsEta=td.make<TH2D>("mee-HFP_vsEta",title.c_str(),10,3,5,30,minZrange,maxZrange);  
  title=std::string("M_{ee,HF-} vs eta ")+post;
  meeHFM_vsEta=td.make<TH2D>("mee-HFM_vsEta",title.c_str(),10,-5,-3,30,minZrange,maxZrange);  

  title=std::string("Y_{ee} ")+post;
  Yee=td.make<TH1D>("yee",title.c_str(),50,-4,4);  
  title=std::string("pT_{ee} ")+post;
  Ptee=td.make<TH1D>("ptee",title.c_str(),100,0,100);  

  title=std::string("eta_{ecal} ")+post;
  ec_eta=td.make<TH1D>("etaecal",title.c_str(),30,-3,3);  
  title=std::string("phi_{ecal} ")+post;
  ec_phi=td.make<TH1D>("phiecal",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("pt_{ecal} ")+post;
  ec_pt=td.make<TH1D>("ptecal",title.c_str(),120,0,120);

  title=std::string("eta_{hf} ")+post;
  hf_eta=td.make<TH1D>("etahf",title.c_str(),50,-5,5);  
  title=std::string("phi_{hf} ")+post;
  hf_phi=td.make<TH1D>("phihf",title.c_str(),30,-3.14159,3.14159);
  title=std::string("pt_{hf} ")+post;
  hf_pt=td.make<TH1D>("pthf",title.c_str(),120,0,120);  

  title=std::string("eta_{hf+} ")+post;
  hfp_eta=td.make<TH1D>("etahfp",title.c_str(),50,-5,5);  
  title=std::string("phi_{hf+} ")+post;
  hfp_phi=td.make<TH1D>("phihfp",title.c_str(),30,-3.14159,3.14159);
  title=std::string("pt_{hf+} ")+post;
  hfp_pt=td.make<TH1D>("pthfp",title.c_str(),120,0,120);  

  title=std::string("eta_{hf-} ")+post;
  hfm_eta=td.make<TH1D>("etahfm",title.c_str(),50,-5,5);  
  title=std::string("phi_{hf-} ")+post;
  hfm_phi=td.make<TH1D>("phihfm",title.c_str(),30,-3.14159,3.14159);
  title=std::string("pt_{hf-} ")+post;
  hfm_pt=td.make<TH1D>("pthfm",title.c_str(),120,0,120);  

  title=std::string("iso e9e25 ")+post;
  hf_e9e25=td.make<TH1D>("e9e25",title.c_str(),60,0.5,1.2);
  title=std::string("eldId var2d")+post;
  hf_var2d=td.make<TH1D>("var2d",title.c_str(),75,0,1.5);  
  title=std::string("eCOREe9 ")+post;
  hf_eCOREe9=td.make<TH1D>("eCOREe9",title.c_str(),60,0,1.2);  
  title=std::string("eSeL ")+post;
  hf_eSeL=td.make<TH1D>("eSeL",title.c_str(),75,0,1.5);  

  title=std::string("seed_eSeL ")+post;
  hf_seed_eSeL=td.make<TH1D>("seed_eSeL",title.c_str(),75,0,1.5);  
  title=std::string("e1e9 ")+post;
  hf_e1e9=td.make<TH1D>("e1e9",title.c_str(),60,0,1.2);
  title=std::string("e ")+post;
  hf_e=td.make<TH1D>("e",title.c_str(),150,0,1500);
  title=std::string("e1x1 ")+post;
  hf_e1x1=td.make<TH1D>("e1x1",title.c_str(),150,0,1500);
  title=std::string("e3x3 ")+post;
  hf_e3x3=td.make<TH1D>("e3x3",title.c_str(),150,0,1500);
  title=std::string("e5x5 ")+post;
  hf_e5x5=td.make<TH1D>("e5x5",title.c_str(),150,0,1500);
  title=std::string("core ")+post;
  hf_core=td.make<TH1D>("core",title.c_str(),150,0,1500);

  title=std::string("eta_{el_1} ")+post;
  el_1_eta=td.make<TH1D>("etael1",title.c_str(),50,-5,5);  
  title=std::string("eta_{el_2} ")+post;
  el_2_eta=td.make<TH1D>("etael2",title.c_str(),50,-5,5);  

  title=std::string("eta_{el_1}_VS_eta_{el_2} ")+post;
  el_1_eta_VS_el_2_eta=td.make<TH2D>("etael1-vs-etael2",title.c_str(),50,-5,5,50,-5,5);
  title=std::string("phi_{el_1}_VS_phi_{el_2} ")+post;
  el_1_phi_VS_el_2_phi=td.make<TH2D>("phiel1-vs-phiel2",title.c_str(),30,-3.14159,3.14159,30,-3.14159,3.14159);

  title=std::string("N-1_combIsoEB ")+post;
  combIsoEB_nmo=td.make<TH1D>("N-1_combIsoEB",title.c_str(),100,0,2);
  title=std::string("N-1_relTkIsoEB ")+post;
  relTkIsoEB_nmo=td.make<TH1D>("N-1_relTkIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1_relEcIsoEB ")+post;
  relEcIsoEB_nmo=td.make<TH1D>("N-1_relEcIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1_relHcIsoEB ")+post;
  relHcIsoEB_nmo=td.make<TH1D>("N-1_relHcIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1_sigiEtaiEtaEB ")+post;
  sigiEtaiEtaEB_nmo=td.make<TH1D>("N-1_sigiEtaiEtaEB",title.c_str(),60,0,0.06);  
  title=std::string("N-1_dphiEB ")+post;
  dphiEB_nmo=td.make<TH1D>("N-1_dphiEB",title.c_str(),75,-0.75,0.75);
  title=std::string("N-1_detaEB ")+post;
  detaEB_nmo=td.make<TH1D>("N-1_detaEB",title.c_str(),40,-0.02,0.02);  
  title=std::string("N-1_hOeEB ")+post;
  hOeEB_nmo=td.make<TH1D>("N-1_hOeEB",title.c_str(),40,0,1);
  title=std::string("N-1_combIsoEE ")+post;
  combIsoEE_nmo=td.make<TH1D>("N-1_combIsoEE",title.c_str(),100,0,2); 
  title=std::string("N-1_relTkIsoEE ")+post;
  relTkIsoEE_nmo=td.make<TH1D>("N-1_relTkIsoEE",title.c_str(),40,0,0.5);  
  title=std::string("N-1_relEcIsoEE ")+post;
  relEcIsoEE_nmo=td.make<TH1D>("N-1_relEcIsoEE",title.c_str(),40,0,0.5);
  title=std::string("N-1_relHcIsoEE ")+post;
  relHcIsoEE_nmo=td.make<TH1D>("N-1_relHcIsoEE",title.c_str(),40,0,0.5);  
  title=std::string("N-1_sigiEtaiEtaEE ")+post;
  sigiEtaiEtaEE_nmo=td.make<TH1D>("N-1_sigiEtaiEtaEE",title.c_str(),60,0,0.06);  
  title=std::string("N-1_dphiEE ")+post;
  dphiEE_nmo=td.make<TH1D>("N-1_dphiEE",title.c_str(),75,-0.75,0.75);
  title=std::string("N-1_detaEE ")+post;
  detaEE_nmo=td.make<TH1D>("N-1_detaEE",title.c_str(),40,-0.02,0.02);  
  title=std::string("N-1_hOeEE ")+post;
  hOeEE_nmo=td.make<TH1D>("N-1_hOeEE",title.c_str(),40,0,0.2);

  title=std::string("N-1_HF iso e9e25 ")+post;
  e9e25_nmo=td.make<TH1D>("N-1_e9e25",title.c_str(),40,0,1);
  title=std::string("N-1_HF eldId var2d")+post;
  var2d_nmo=td.make<TH1D>("N-1_var2d",title.c_str(),40,0,1.5);  
  title=std::string("N-1_eCOREe9 ")+post;
  eCOREe9_nmo=td.make<TH1D>("N-1_eCOREe9",title.c_str(),40,0,1);  
  title=std::string("N-1_eSeL ")+post;
  eSeL_nmo=td.make<TH1D>("N-1_eSeL",title.c_str(),75,0,1.5);
  title=std::string("eSeL_vs_logEl ")+post;
  eSeL_vs_logEl=td.make<TH2D>("eSeL_vs_logEl",title.c_str(),25,3,8,75,0,1.5);
}

void HFZeeVBTF::HistPerDef::fill(pat::ElectronCollection::const_iterator ecalE,  
				 const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, 
				 const bool hasEleIDPassed,
				 const edm::ParameterSet& myPs, 
				 const int isTight,
				 std::vector<double> hfIdParams) {

  float e9e25_cut;
  float hf_2d_cut;
  float eCOREe9_cut;
  float eSeL_cut;

  if (isTight==0){     // this is the loose HFeleId
    e9e25_cut   = hfIdParams.at(0);  
    hf_2d_cut   = hfIdParams.at(2);
    eCOREe9_cut = hfIdParams.at(4);
    eSeL_cut    = hfIdParams.at(6);}
  else if (isTight==1){// this is the tight HFeleId
    e9e25_cut   = hfIdParams.at(1);  
    hf_2d_cut   = hfIdParams.at(3);
    eCOREe9_cut = hfIdParams.at(5);
    eSeL_cut    = hfIdParams.at(7);}
  else if (isTight==-1) {
    e9e25_cut=-1;
    hf_2d_cut=-1;
    eCOREe9_cut=-1;
    eSeL_cut=-1;
  } else 
    {assert(0);} // no values otherthan 1 and 0 are supported for the moment
  
  // for debugging
  //  std::cout << "hf_2d_cut: " << hf_2d_cut << " hfIdParams[2]: (loose) " << hfIdParams.at(2) << " hfIdParams[3]: (tight) " << hfIdParams.at(3) << std::endl;
  //  std::cout << "e9e25_cut: " << e9e25_cut << " hfIdParams[0]: (loose) " << hfIdParams.at(0) << " hfIdParams[1]: (tight) " << hfIdParams.at(1) << std::endl;
  //  std::cout << "eCOREe9: " << eCOREe9_cut << "\t" << hfIdParams.at(4) << "\t" << hfIdParams.at(5) << std::endl; 
  //  std::cout << "eSeL: " << eSeL_cut << "\t" << hfIdParams.at(6) << "\t" <<  hfIdParams.at(7) << std::endl; 
  
  

  bool isEb(false);
  int isEe=1;// 0: EB;  1: EE
  if(fabs(ecalE->eta())<1.4442){
    isEb=true;
    isEe=0;
  }
  
  float e9e25      = hfshape.eLong3x3()/hfshape.eLong5x5();
  float var2d      = hfshape.eCOREe9()-(hfshape.eSeL()*9./8.);
  float eCOREe9    = hfshape.eCOREe9();
  float eSeL       = hfshape.eSeL();
  
  reco::Particle::LorentzVector Z(ecalE->p4());
  Z+=hfE.p4();

  // if all selections are passed, fill standard plots
  if (isTight==-1 || (hasEleIDPassed && 
     var2d   > hf_2d_cut   &&
     e9e25   > e9e25_cut   &&
     eCOREe9 > eCOREe9_cut &&
     eSeL    < eSeL_cut    &&
     massWindow_[0] < Z.M() && Z.M() < massWindow_[1])){ // make inv. mass requirement explicit
    
    mee          ->Fill(Z.M());
    mee_vsEta    -> Fill(hfE.p4().eta(),Z.M());
    mee_vsAbsEta -> Fill(fabs(hfE.p4().eta()),Z.M());
    eSeL_vs_logEl-> Fill(log(hfshape.eLong3x3()),hfshape.eSeL());

    if(hfE.p4().eta()>0){
      meeHFP   ->Fill(Z.M());
      hfp_eta  ->Fill(hfE.p4().eta());
      hfp_phi  ->Fill(hfE.p4().phi());
      hfp_pt   ->Fill(hfE.pt());
      meeHFP_vsEta -> Fill(hfE.p4().eta(),Z.M());
    }
    else{
      meeHFM   ->Fill(Z.M());
      hfm_eta  ->Fill(hfE.p4().eta());
      hfm_phi  ->Fill(hfE.p4().phi());
      hfm_pt   ->Fill(hfE.pt());
      meeHFM_vsEta -> Fill(hfE.p4().eta(),Z.M());
    }
    
    Yee   ->Fill(Z.Rapidity());
    Ptee  ->Fill(Z.pt());
    
    ec_eta->Fill(ecalE->eta());
    
    ec_phi->Fill(ecalE->phi());
    ec_pt ->Fill(ecalE->pt());
    hf_eta->Fill(hfE.eta());
    hf_phi->Fill(hfE.phi());
    hf_pt ->Fill(hfE.pt());
    // HF Id variables after the cuts, as a cross check
    hf_e9e25  ->Fill(e9e25);
    hf_var2d  ->Fill(var2d);
    hf_eCOREe9->Fill(hfshape.eCOREe9());
    hf_eSeL   ->Fill(hfshape.eSeL());
    

    hf_seed_eSeL->Fill(hfshape.eShort1x1()/hfshape.eLong1x1());
    hf_e1e9     ->Fill(hfshape.eLong1x1()/hfshape.eLong3x3());
    hf_e        ->Fill(hfshape.eLong3x3()); // GF: check this with Kevin
    hf_e1x1     ->Fill(hfshape.e1x1());
    hf_e3x3     ->Fill(hfshape.e3x3());
    hf_e5x5     ->Fill(hfshape.e5x5());
    hf_core     ->Fill(hfshape.eCore());


    if(ecalE->pt() > hfE.pt()){
      el_1_eta->Fill(ecalE->eta());
      el_2_eta->Fill(hfE.eta());
      el_1_eta_VS_el_2_eta->Fill(ecalE->eta(),hfE.eta());
      el_1_phi_VS_el_2_phi->Fill(ecalE->phi(),hfE.phi());
    }
    else{
      el_1_eta->Fill(hfE.eta());
      el_2_eta->Fill(ecalE->eta());
      el_1_eta_VS_el_2_eta->Fill(hfE.eta(),ecalE->eta());
      el_1_phi_VS_el_2_phi->Fill(hfE.phi(),ecalE->phi());
    }
 
  }// END
  // if all selections are passed, fill standard plots
  
  //////////////////////////////////////////////////////////////////////////////////////
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID
  // here follows filling of histograms for selection variables when the other N-1_have passed
  // which cuts values are used? Depends on which parameter set you pass...
  std::vector<double> ebParam = myPs.getParameter< std::vector<double> >("barrel"); 
  std::vector<double> eeParam = myPs.getParameter< std::vector<double> >("endcap"); 


  float combIsolation_cut[2];
  float trackRel03_cut[2];
  float ecalRel03_cut[2];
  float hcalRel03_cut[2];
  float sigiEtaiEta_cut[2];
  float dPhi_cut[2];
  float dEta_cut[2];
  float HoE_cut[2];


  bool verboseDebug(false);
  
  // EB isolation variables
  combIsolation_cut[0] = ebParam.at(17); if(verboseDebug) std::cout  << "\n\ncombIsolationEB: " << combIsolation_cut[0] << std::endl;
  
  trackRel03_cut[0]    = ebParam.at(11); if(verboseDebug) std::cout  << "trackRel03EB: " << trackRel03_cut[0] << std::endl;
  ecalRel03_cut[0]     = ebParam.at(12); if(verboseDebug) std::cout  << "ecalRel03EB: " << ecalRel03_cut[0] << std::endl;
  hcalRel03_cut[0]     = ebParam.at(13); if(verboseDebug) std::cout  << "hcalRel03EB: " << hcalRel03_cut[0] << std::endl;
  
  // EB eleId variables
  sigiEtaiEta_cut[0]  = ebParam.at(1);  if(verboseDebug) std::cout  << "sigiEetaiEtaEB: " << sigiEtaiEta_cut[0] << std::endl;
  dPhi_cut[0]         = ebParam.at(2);  if(verboseDebug) std::cout  << "DphiEB: " << dPhi_cut[0] << std::endl;
  dEta_cut[0]         = ebParam.at(3);  if(verboseDebug) std::cout  << "DetaEB: " << dEta_cut[0] << std::endl;
  HoE_cut[0]          = ebParam.at(0);  if(verboseDebug) std::cout  << "HoEEB: " << HoE_cut[0] << std::endl;

  
  // EE isolation variables
  combIsolation_cut[1] = eeParam.at(17); if(verboseDebug) std::cout  << "\n\ncombIsolationEE: " << combIsolation_cut[1] << std::endl;
  
  trackRel03_cut[1]    = eeParam.at(11); if(verboseDebug) std::cout  << "trackRel03EE: " << trackRel03_cut[1] << std::endl;
  ecalRel03_cut[1]     = eeParam.at(12); if(verboseDebug) std::cout  << "ecalRel03EE: " << ecalRel03_cut[1] << std::endl;
  hcalRel03_cut[1]     = eeParam.at(13); if(verboseDebug) std::cout  << "hcalRel03EE: " << hcalRel03_cut[1] << std::endl;
  
  // EE eleId variables
  sigiEtaiEta_cut[1]  = eeParam.at(1);  if(verboseDebug) std::cout  << "sigiEetaiEtaEE: " << sigiEtaiEta_cut[1] << std::endl;
  dPhi_cut[1]          = eeParam.at(2);  if(verboseDebug) std::cout  << "DphiEE: " << dPhi_cut[1] << std::endl;
  dEta_cut[1]          = eeParam.at(3);  if(verboseDebug) std::cout  << "DetaEE: " << dEta_cut[1] << std::endl;
  HoE_cut[1]           = eeParam.at(0);  if(verboseDebug) std::cout  << "HoEEE: " << HoE_cut[1] << std::endl;
  

  //  //////////////////////////////////////////////////////////////////////////////////////
  //  // now get hold of the actual variables
  //  //hfE

  float combinedEcalIso;
  if (isEb) combinedEcalIso  = ( ecalE->dr03TkSumPt() + max(0., ecalE->dr03EcalRecHitSumEt() - 1.) + ecalE->dr03HcalTowerSumEt() ) / ecalE->p4().Pt();
  else      combinedEcalIso  = ( ecalE->dr03TkSumPt() + ecalE->dr03EcalRecHitSumEt() + ecalE->dr03HcalTowerSumEt() ) / ecalE->p4().Pt();

  // relative iso variable
  float trackRel03 = ecalE->dr03TkSumPt()/ecalE->p4().Pt();
  float ecalRel03  = ecalE->dr03EcalRecHitSumEt()/ecalE->p4().Pt(); 
  float hcalRel03  = ecalE->dr03HcalTowerSumEt()/ecalE->p4().Pt();

  // eleID variables  
  float sigiEtaiEta = ecalE->scSigmaIEtaIEta();
  float dEta        = ecalE->deltaEtaSuperClusterTrackAtVtx();
  float dPhi        = ecalE->deltaPhiSuperClusterTrackAtVtx();
  float HoE         = ecalE->hadronicOverEm();

  //////////////////////////////////////////////////////////////////////////////////////
  // now establish whether cuts are passed or not; 0-9 ECAL, 10-19 HF
  bool  cut[20];  for (int v =0; v<20; v++) cut[v]=false;
  short passedCuts=0;
  // ECAL electron cuts  
  if (combinedEcalIso < combIsolation_cut[isEe]) {cut[0] = true; passedCuts++; /*std::cout<<"\t c0 ";*/}
  // if (trackRel03      < trackRel03_cut[isEe])   {cut[1] = true; passedCuts++;/*std::cout<<"c1 ";*/}
  // if (ecalRel03       < ecalRel03_cut[isEe])     {cut[2] = true; passedCuts++;/*std::cout<<"c2 ";*/}
  // if (hcalRel03       < hcalRel03_cut[isEe])     {cut[3] = true; passedCuts++;/*std::cout<<"c3 ";*/}
  if (sigiEtaiEta     < sigiEtaiEta_cut[isEe])   {cut[4] = true; passedCuts++;/*std::cout<<"c4 ";*/}
  if (dEta            < dEta_cut[isEe])          {cut[5] = true; passedCuts++;/*std::cout<<"c5 ";*/}
  if (dPhi            < dPhi_cut[isEe])          {cut[6] = true; passedCuts++;/*std::cout<<"c6 ";*/}
  if (HoE             < HoE_cut[isEe])           {cut[7] = true; passedCuts++;/*std::cout<<"c7 ";*/}

  if (e9e25           > e9e25_cut)               {cut[10] = true; passedCuts++;/*std::cout<<"c10 ";*/}
  if (var2d           > hf_2d_cut)               {cut[11] = true; passedCuts++;/*std::cout<<"c11 ";*/}
  if (eSeL            < eSeL_cut )               {cut[12] = true; passedCuts++;/*std::cout<<"c12 ";*/}
  if (eCOREe9         > eCOREe9_cut)             {cut[13] = true; passedCuts++;/*std::cout<<"c13 ";*/}
  /*std::cout<<std::endl;*/
  short numCuts = 9; // was 7 when eSeL and eCOREe9 were not considered sigularly


  //////////////////////////////////////////////////////////////////////////////////////
  // all cuts passed     OR   all passed except present variable
  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[0])) ) {
    if(isEb) combIsoEB_nmo -> Fill(combinedEcalIso);
    else     combIsoEE_nmo -> Fill(combinedEcalIso);
  }

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[4])) ) {
    if(isEb) sigiEtaiEtaEB_nmo -> Fill(sigiEtaiEta);
    else     sigiEtaiEtaEE_nmo -> Fill(sigiEtaiEta);
  }

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[5])) ) {
    if(isEb) detaEB_nmo -> Fill(dEta);
    else     detaEE_nmo -> Fill(dEta);
  }

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[6])) ) {
    if(isEb) dphiEB_nmo -> Fill(dPhi);
    else     dphiEE_nmo -> Fill(dPhi);
  }

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[7])) ) {
    if(isEb) hOeEB_nmo -> Fill(HoE);
    else     hOeEE_nmo -> Fill(HoE);
  }

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[10])) ) {
    e9e25_nmo -> Fill(e9e25);
  }

  // the following three variables are correlated
  // assume that only vard2 or eSeL&eCOREe9 will be used
  // and when not used, give huge values
  
  // 2d var
  if(      ( passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[11]))  )
	   && fabs(hf_2d_cut)<5
	   ) {
    var2d_nmo -> Fill(var2d); eSeL_nmo-> Fill(eSeL); eCOREe9_nmo -> Fill(eCOREe9);
  }
  // eSeL
  if(      ( passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[12]))  )
	   && fabs(eSeL_cut)<5
	   ) {
    var2d_nmo -> Fill(var2d); eSeL_nmo-> Fill(eSeL); eCOREe9_nmo -> Fill(eCOREe9);
  }
  // eCOREe9
  if(      ( passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[13]))  )
	   && fabs(eCOREe9_cut)<5
	   ) {
    var2d_nmo -> Fill(var2d); eSeL_nmo-> Fill(eSeL); eCOREe9_nmo -> Fill(eCOREe9);
  }

  
}// end of fill()

// this value for tight is lower than CMS AN-2009/106 to accont for S/L miscalib in data
// see figure 15 therein: turn-on is sharp
static const double hf_2d_loose   = 0.32;
static const double hf_2d_tight   = 0.45;
static const double e9e25_loose   = 0.90;
static const double e9e25_tight   = 0.94;
//////////////////////////////////////////////////////////////////////////////////////
// FIXME
// gf: bring HF variables to PSet as well, as opposed to having them hard-coded



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HFZeeVBTF::HFZeeVBTF(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  ecalID_              = iConfig.getParameter<std::string>("ECALid");
  dolog_               = iConfig.getParameter<bool>("DoLog");
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  filterId_=iConfig.getUntrackedParameter<std::string>("FILTERid","");
  hfIdParams_          = iConfig.getParameter< std::vector<double> >("hFselParams");
  hFElectrons_         = iConfig.getParameter< edm::InputTag> ("hFElectrons");
  minEtECAL_           = iConfig.getParameter< double >("minEtECAL");
  minEtHF_             = iConfig.getParameter< double >("minEtHF");
  myName_              = iConfig.getParameter<std::string>("myName");
  skipTrigger_         = iConfig.getParameter<bool>("skipTrigger");
  massWindow_          = iConfig.getParameter< std::vector<double> >("Zmass");

  myName_+=std::string("    ");
  
  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    w.add(*it);//      std::cout << "adding: " << (*it) << std::endl;
    }

  edm::Service<TFileService> fs;
  hists.nelec=fs->make<TH1D>("nelec","N_Elec",10,-0.5,9.5);
  hists.nhf=fs  ->make<TH1D>("nhf","N_HF",10,-0.5,9.5);
  hists.base.book(fs->mkdir("base"),"(base)",massWindow_);
  hists.mee90none.book(fs->mkdir("mee90none"),"(mee90none)",massWindow_);
  hists.zMass.book(fs->mkdir("zMass"),"(zMass)",massWindow_);
  hists.mee90loose.book(fs->mkdir("mee90loose"),"(90,loose)",massWindow_);
  hists.mee80loose.book(fs->mkdir("mee80loose"),"(80,loose)",massWindow_);
  hists.mee70loose.book(fs->mkdir("mee70loose"),"(70,loose)",massWindow_);  
  hists.mee60loose.book(fs->mkdir("mee60loose"),"(60,loose)",massWindow_);  
  hists.mee90tight.book(fs->mkdir("mee90tight"),"(90,tight)",massWindow_);  
  hists.mee80tight.book(fs->mkdir("mee80tight"),"(80,tight)",massWindow_);  
  hists.mee70tight.book(fs->mkdir("mee70tight"),"(70,tight)",massWindow_);  
  hists.mee60tight.book(fs->mkdir("mee60tight"),"(60,tight)",massWindow_);  
  init_=false;

  // import parameter set which carry threshold vaues
  robust95relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust95relIsoEleIDCutsV04");
  robust90relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust90relIsoEleIDCutsV04");
  robust85relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust85relIsoEleIDCutsV04");
  robust80relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust80relIsoEleIDCutsV04");
  robust70relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust70relIsoEleIDCutsV04");
  robust60relIsoEleIDCutsV04_ps_ = iConfig.getParameter<edm::ParameterSet>("robust60relIsoEleIDCutsV04");
  
  robust95cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust95cIsoEleIDCutsV04");
  robust90cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust90cIsoEleIDCutsV04");
  robust85cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust85cIsoEleIDCutsV04");
  robust80cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust80cIsoEleIDCutsV04");
  robust70cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust70cIsoEleIDCutsV04");
  robust60cIsoEleIDCutsV04_ps_   = iConfig.getParameter<edm::ParameterSet>("robust60cIsoEleIDCutsV04");

}
  
HFZeeVBTF::~HFZeeVBTF()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool
HFZeeVBTF::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  bool filterOk=false;

  edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
  edm::Handle<reco::SuperClusterCollection> SuperClusters;
  edm::Handle<reco::HFEMClusterShapeAssociationCollection> ClusterAssociation;

  iEvent.getByLabel("hfEMClusters",SuperClusters);
  
  //std::cout <<  "\n\n debug inputTag hFElectrons is: " << hFElectrons_ << std::endl;

  //  ++recoSuperClustersToOnerecoHFEMClusterShapesAssociation "hfEMClusters" "" "RECO" (productId = 2:697)
  iEvent.getByLabel("hfEMClusters",ClusterAssociation);

  //  InputTag:                                     label = , instance = HFElectrons
  //  ++recoRecoEcalCandidates "hfRecoEcalCandidate" "" "RECO" (productId = 2:680)
  //  ++recoRecoEcalCandidates "hfRecoEcalCandidate" "" "hfzeevbtfAnalysFromSkim" (productId = 4:2)
  //  ++recoRecoEcalCandidates "hfRecoEcalCandidateLoose" "" "hfzeevbtfAnalysFromSkim" (productId = 4:3)
  //iEvent.getByLabel("hfRecoEcalCandidate",HFElectrons);
  iEvent.getByLabel(hFElectrons_,HFElectrons);

  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    std::cout << "No electrons found in this event with tag " 
	      << "patElectrons" << std::endl;
    return false; // RETURN if no elecs in the event
  }

  const pat::ElectronCollection& pElecs = *(patElectron.product());
  pat::ElectronCollection::const_iterator i, ecalE=pElecs.end();
  
  for (i=pElecs.begin(); i!=pElecs.end(); i++) {
    reco::SuperClusterRef scr=i->superCluster();
    if (scr.isNull()) {
      std::cout << "I'm NULL!!!\n";
      continue;
    }

    // gf: ECAL eta acceptance cut on supercluster. What about phi acceptance and central crack?
    // explanation of electronID value: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
    double eta_det=scr.get()->eta();
    if (ecalE==pElecs.end() && i->electronID(ecalID_)>=(.0-0.1)   // very loose to let events pass; cuts done later
	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)
	&& HFElectrons->size()>0) {
      const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0];
      reco::SuperClusterRef hfclusRef=hfE.superCluster();
      const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
      // const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;
    }
    
    // ECAL acceptance cut
    if (ecalE==pElecs.end() && w.doesElePass( i->electronID(ecalID_) )  // 3: passes electron iso/Id/convRej depending on configuration; ecalID_ defines working point
    	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)) {             // ECAL acceptance cut
      ecalE=i;
    }// if
    
  }// end loop over electrons
  
  hists.nelec->Fill(pElecs.size());
  hists.nhf  ->Fill(HFElectrons->size());

  if (pElecs.size()>0 && (HFElectrons->size()>0)) {
    unsigned int hfEleWithMaxPt = 0;
    for(unsigned int u=0; u<(HFElectrons->size()); u++){
      if ( (*HFElectrons).at(u).pt() > (*HFElectrons).at(hfEleWithMaxPt).pt() )  hfEleWithMaxPt = u;
    }
    const reco::RecoEcalCandidate& hfE=(*HFElectrons).at(hfEleWithMaxPt);
    
    reco::SuperClusterRef hfclusRef=hfE.superCluster();
    const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
    const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;

    if (ecalE!=pElecs.end()) 
      hists.base.fill(ecalE,hfE,hfshape,true,robust90relIsoEleIDCutsV04_ps_,-1,hfIdParams_);
    else
      hists.base.fill(pElecs.begin(),hfE,hfshape,true,robust90relIsoEleIDCutsV04_ps_,-1,hfIdParams_);
  }

  if ( (pElecs.size()>0) &&  (HFElectrons->size()>0) && (ecalE!=pElecs.end())){
    if(dolog_ && 0) std::cout << myName_  << "size ele coll:" << pElecs.size() << "  size HF coll: " << HFElectrons->size() << std::endl;
    // add here histograms w/o request of Z mass or eleId
    // mee90none is at the very beginning: one HF supercluste and one GSF w/ ID and iso
    // const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0]; // choose first HF candidate - old crude way!
    // find highest pt HF candidate
    unsigned int hfEleWithMaxPt = 0;
    for(unsigned int u=0; u<(HFElectrons->size()); u++){
      if ( (*HFElectrons).at(u).pt() > (*HFElectrons).at(hfEleWithMaxPt).pt() )  hfEleWithMaxPt = u;
    }
    const reco::RecoEcalCandidate& hfE=(*HFElectrons).at(hfEleWithMaxPt);
    
    reco::SuperClusterRef hfclusRef=hfE.superCluster();
    const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
    const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;
    hists.mee90none.fill(ecalE,hfE,hfshape,true,robust90relIsoEleIDCutsV04_ps_,0,hfIdParams_);//gf: for each Z definition, N-1 control plots could be filled here too
  }

  // gf: up to this point requirements are:
  // 1) at least one ECAL electron; ecalE is ==first== ECAL electron with SC within acceptance
  //    && electron ID and Isolation only
  // 2) at least one HF electrons
  // than try match combinatorics
  
  if (ecalE!=pElecs.end() && HFElectrons->size()>0) {
    //const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0]; // choose first HF candidate - old crude way!
    // find highest pt HF candidate
    unsigned int hfEleWithMaxPt = 0;
    for(unsigned int u=0; u<(HFElectrons->size()); u++){
      if ( (*HFElectrons).at(u).pt() > (*HFElectrons).at(hfEleWithMaxPt).pt() )  hfEleWithMaxPt = u;
    }
    const reco::RecoEcalCandidate& hfE=(*HFElectrons).at(hfEleWithMaxPt);
    
     reco::SuperClusterRef           hfclusRef=hfE.superCluster();
     const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
     const reco::HFEMClusterShape&   hfshape=*hfclusShapeRef;
    
    // make a Z
    reco::Particle::LorentzVector Z(ecalE->p4());
    Z+=hfE.p4();

    // gf: 3) loose Z mass cut and min pt
    // ==> these are the CANDIATES
    if (Z.M()>massWindow_[0] && Z.M()<massWindow_[1] && ecalE->pt()>minEtECAL_ && hfE.pt()>minEtHF_) {
      double var2d=hfshape.eCOREe9()-(hfshape.eSeL()*1.125);
      double eta_det=ecalE->superCluster().get()->eta();
      if (dolog_) {
	std::cout << myName_  << "------------------------------------------------------" << std::endl;
	
	std::cout << myName_  << "Candidate!: event " << iEvent.id().event() << " run " << iEvent.id().run();
	
	std::cout << myName_  << "  m_ee=" << Z.M() << " Y_ee=" << Z.Rapidity() << " pT_ee=" << Z.pt() <<  std::endl;
	std::cout << myName_  << " file: " << currentFile_ << std::endl;
	
	// https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID#Electron_ID_Implementation_in_Re
	std::cout << myName_  << "ECAL (pt,eta,phi, eta_det) : " << ecalE->pt() << ", " << ecalE->eta() << ", " << ecalE->phi() << " , (" << eta_det << ")" << std::endl;
	//DELETE THIS std::cout << myName_  << "ECAL (pt,tkIso,phi, eta_det) : " << ecalE->pt() << ", " << ecalE->trackIso() << ", " << ecalE->phi() << " , (" << eta_det << ")" << std::endl;
	// gf: pull up here all the ECALId variables as descrived in twiki here above 
	std::cout << myName_  << "\t ECAL ele // REL ISOLATION -  trackRel03: " << ecalE->dr03TkSumPt()/ecalE->p4().Pt() 
		  << "\t ecalRel03: "  << ecalE->dr03EcalRecHitSumEt()/ecalE->p4().Pt() 
		  << "\t hcalRel03: "  <<  ecalE->dr03HcalTowerSumEt()/ecalE->p4().Pt() 
		  << "\n\t Electron ID - sigIetaIeta: " << ecalE->scSigmaIEtaIEta()
		  << "\t Deta: " << ecalE->deltaEtaSuperClusterTrackAtVtx()
		  << "\t Dphi: " << ecalE->deltaPhiSuperClusterTrackAtVtx()
		  << "\t H/E: "  << ecalE->hadronicOverEm()
		  << std::endl;

	std::cout << myName_  << "  ID (60,70,80,90) : " 
		  << ecalE->electronID("simpleEleId60relIso") << " " 
		  << ecalE->electronID("simpleEleId70relIso") << " " 
		  << ecalE->electronID("simpleEleId80relIso") << " " 
		  << ecalE->electronID("simpleEleId90relIso") << " " 
		  << std::endl;
	
	
	std::cout << myName_  << "  HF (pt, eta, phi): " << hfE.pt() << ", " << hfE.eta() << ", " << hfE.phi() << std::endl;
	std::cout << myName_  << "  Shape (S/L, ec/e9, e1/e9, e9/e25, 2d) : " 
		  << hfshape.eSeL() << " "
		  << hfshape.eCOREe9() << " "
		  << hfshape.e1x1()/hfshape.e3x3() << " "
		  << hfshape.e9e25() << " "
		  << var2d << " "
		  << std::endl;

	if (!skipTrigger_) {

	  std::cout << myName_  << "\n Triggers fired: ";
	
	  Handle<TriggerResults> hltResults ;
	  iEvent.getByLabel( edm::InputTag("TriggerResults","","HLT"), hltResults ) ;
      
	  if (!init_) {
	    init_=true;
	    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);
	    HLT_Names=triggerNames.triggerNames();
	  }
	  
	  int nf=0;
	  // checking against both size of names and of results: found events when they're diferent (someonelse's bug?)
	  for (unsigned int i=0; i<HLT_Names.size() && i<hltResults->size(); i++) {
	    if ( hltResults->accept(i) ) {
	      if ((nf % 5)==0) std::cout << myName_  << "\n    ";
	      std::cout << myName_  << HLT_Names.at(i) << " ";
	      nf++;
	  }
	  }
	  std::cout << myName_  << std::endl;
	}
	std::cout << myName_  << "======================================================" << std::endl;
      }// 4) here require electronId to have passed, for different working points
      //  5) and require HF eleID loose or tight
      //  preselections on HF isolation (0.9) and 2dvar (0.3) are already applied at clustering level.
      //  loosen as much as reasonable here

      hists.zMass.fill(ecalE,hfE,hfshape,true,robust90relIsoEleIDCutsV04_ps_,0,hfIdParams_);
      hists.mee90loose.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId90relIso")),robust90relIsoEleIDCutsV04_ps_,0,hfIdParams_);
      hists.mee80loose.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId80relIso")),robust80relIsoEleIDCutsV04_ps_,0,hfIdParams_);
      hists.mee70loose.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId70relIso")),robust70relIsoEleIDCutsV04_ps_,0,hfIdParams_);
      hists.mee60loose.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId60relIso")),robust60relIsoEleIDCutsV04_ps_,0,hfIdParams_);

      hists.mee90tight.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId90relIso")),robust90relIsoEleIDCutsV04_ps_,1,hfIdParams_);
      hists.mee80tight.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId80relIso")),robust80relIsoEleIDCutsV04_ps_,1,hfIdParams_);
      hists.mee70tight.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId70relIso")),robust70relIsoEleIDCutsV04_ps_,1,hfIdParams_);
      hists.mee60tight.fill(ecalE,hfE,hfshape,w.doesElePass( ecalE->electronID("simpleEleId60relIso")),robust60relIsoEleIDCutsV04_ps_,1,hfIdParams_);

      if (!filterId_.empty()) {
	filterOk=w.doesElePass( ecalE->electronID(filterId_) ) && var2d>hf_2d_loose;
      } 

    }// if mass is in Z window
    
  }// if the is at least on ECAL and one HF candidate

  if (filterId_.empty()) filterOk=true;

  return filterOk;
  
}


// reminder re electronID(""):
//0: fails
//1: passes electron ID only
//2: passes electron Isolation only
//3: passes electron ID and Isolation only   <<< 
//4: passes conversion rejection
//5: passes conversion rejection and ID
//6: passes conversion rejection and Isolation
//7: passes the whole selection
//https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID 

// ------------ method called once each job just before starting event loop  ------------
void 
HFZeeVBTF::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HFZeeVBTF::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFZeeVBTF);
