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
// $Id: HFZeeVBTF.cc,v 1.2 2010/07/20 17:32:11 franzoni Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}

//////////////////////////////////////////////////////////////////
// little helper class to manage selection of electron candidates
class selectElectron {
  int x, y;
  std::vector<int> allowed;
public:
  void add(int);
  bool doesElePass(int);
} rect;

void selectElectron::add(int v)
{ allowed.push_back(v);}

bool selectElectron::doesElePass(int v)
{ 
  if ( find(allowed.begin(),allowed.end(),v) != allowed.end() ) return true;
  else return false;
}


//// original implementation; and how Jeremy submitted MC jobs for bck
//bool doesElePass(int k)
//{
//  if (k<0 || k>7) {
//    std::cout << "value entering doesElePass makes no sense" << std::endl;
//    return false;
//  }
//  else if (k < 3)         return false;
//  else if (k==4)          return false;
//  else                    return true;
//}

// tighter version proposed by Jason: require at least "electron ID and Isolation" (3), possibly with conversion rejection as well (7: passes the whole selection)
//bool doesElePass(int k)
//{
//  if (k<0 || k>7) {
//    std::cout << "value entering doesElePass makes no sense" << std::endl;
//    return false;
//  }
//  //if (k==3 || k==7)         return true;
//  if (k==7)         return true;
//  else                      return false;
//}



class HFZeeVBTF : public edm::EDAnalyzer {
public:
  explicit HFZeeVBTF(const edm::ParameterSet&);
  ~HFZeeVBTF();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string ecalID_, currentFile_;
  bool dolog_;
  std::vector<int> acceptedElectronIDs_;
  selectElectron w;

  // ----------member data ---------------------------

  // gf set of histograms per each == Z definition ==
  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td,const std::string&);
    // fill all histos of the set with the two electron candidates
    void fill(pat::ElectronCollection::const_iterator ecalE,  const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, const edm::ParameterSet& myPs, const float hf_2d_cut);
    TH1* mee, *meeHFP, *meeHFM, *Yee, *Ptee;
    TH1* ec_eta, *ec_phi, *ec_pt;
    TH1* hf_eta, *hf_phi, *hf_pt;
    TH1* el_1_eta, *el_2_eta;
    TH2* el_1_eta_VS_el_2_eta, *el_1_phi_VS_el_2_phi;
    // gf: add hisos for variables N-1 selected
    TH1* combIsoEB_nmo,  *relTkIsoEB_nmo, *relEcIsoEB_nmo, *relHcIsoEB_nmo, *sigiEtaiEtaEB_nmo, *dphiEB_nmo, *detaEB_nmo, *hOeEB_nmo;
    TH1* combIsoEE_nmo,  *relTkIsoEE_nmo, *relEcIsoEE_nmo, *relHcIsoEE_nmo, *sigiEtaiEtaEE_nmo, *dphiEE_nmo, *detaEE_nmo, *hOeEE_nmo;
    TH1* e9e25_nmo, *var2d_nmo, *eCOREe9_nmo, *eSeL_nmo;
  };
std::vector<std::string> HLT_Names;
bool init_;

  // gf set of histo for all Z definitios in a stack
  struct HistStruct {
    TH1 *nelec,*nhf;
    HistPerDef /*base,*/ mee90none, zMass;
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


void HFZeeVBTF::HistPerDef::book(TFileDirectory td, const std::string& post) {
  std::string title;

  title=std::string("M_{ee} ")+post;
  mee=td.make<TH1D>("mee",title.c_str(),90,40,130);  
  title=std::string("M_{ee,HF+} ")+post;
  meeHFP=td.make<TH1D>("mee-HFP",title.c_str(),90,40,130);  
  title=std::string("M_{ee,HF-} ")+post;
  meeHFM=td.make<TH1D>("mee-HFM",title.c_str(),90,40,130);  
  title=std::string("Y_{ee} ")+post;
  Yee=td.make<TH1D>("yee",title.c_str(),50,-4,4);  
  title=std::string("pT_{ee} ")+post;
  Ptee=td.make<TH1D>("ptee",title.c_str(),90,0,90);  

  title=std::string("eta_{ecal} ")+post;
  ec_eta=td.make<TH1D>("etaecal",title.c_str(),30,-3,3);  
  title=std::string("phi_{ecal} ")+post;
  ec_phi=td.make<TH1D>("phiecal",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("pt_{ecal} ")+post;
  ec_pt=td.make<TH1D>("ptecal",title.c_str(),90,0,90);  

  title=std::string("eta_{hf} ")+post;
  hf_eta=td.make<TH1D>("etahf",title.c_str(),50,-5,5);  
  title=std::string("phi_{hf} ")+post;
  hf_phi=td.make<TH1D>("phihf",title.c_str(),30,-3.14159,3.14159);
  title=std::string("pt_{hf} ")+post;
  hf_pt=td.make<TH1D>("pthf",title.c_str(),90,0,90);  


  title=std::string("eta_{el_1} ")+post;
  el_1_eta=td.make<TH1D>("etael1",title.c_str(),50,-5,5);  
  title=std::string("eta_{el_2} ")+post;
  el_2_eta=td.make<TH1D>("etael2",title.c_str(),50,-5,5);  

  title=std::string("eta_{el_1}_VS_eta_{el_2} ")+post;
  el_1_eta_VS_el_2_eta=td.make<TH2D>("etael1-vs-etael2",title.c_str(),50,-5,5,50,-5,5);
  title=std::string("phi_{el_1}_VS_phi_{el_2} ")+post;
  el_1_phi_VS_el_2_phi=td.make<TH2D>("phiel1-vs-phiel2",title.c_str(),30,-3.14159,3.14159,30,-3.14159,3.14159);

  title=std::string("N-1 combIsoEB ")+post;
  combIsoEB_nmo=td.make<TH1D>("N-1 combIsoEB",title.c_str(),100,0,2);
  title=std::string("N-1 relTkIsoEB ")+post;
  relTkIsoEB_nmo=td.make<TH1D>("N-1 relTkIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1 relEcIsoEB ")+post;
  relEcIsoEB_nmo=td.make<TH1D>("N-1 relEcIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1 relHcIsoEB ")+post;
  relHcIsoEB_nmo=td.make<TH1D>("N-1 relHcIsoEB",title.c_str(),40,0,0.5);  
  title=std::string("N-1 sigiEtaiEtaEB ")+post;
  sigiEtaiEtaEB_nmo=td.make<TH1D>("N-1 sigiEtaiEtaEB",title.c_str(),60,0,0.06);  
  title=std::string("N-1 dphiEB ")+post;
  dphiEB_nmo=td.make<TH1D>("N-1 dphiEB",title.c_str(),75,-0.75,0.75);
  title=std::string("N-1 detaEB ")+post;
  detaEB_nmo=td.make<TH1D>("N-1 detaEB",title.c_str(),40,-0.02,0.02);  
  title=std::string("N-1 hOeEB ")+post;
  hOeEB_nmo=td.make<TH1D>("N-1 hOeEB",title.c_str(),40,0,1);
  title=std::string("N-1 combIsoEE ")+post;
  combIsoEE_nmo=td.make<TH1D>("N-1 combIsoEE",title.c_str(),100,0,2); 
  title=std::string("N-1 relTkIsoEE ")+post;
  relTkIsoEE_nmo=td.make<TH1D>("N-1 relTkIsoEE",title.c_str(),40,0,0.5);  
  title=std::string("N-1 relEcIsoEE ")+post;
  relEcIsoEE_nmo=td.make<TH1D>("N-1 relEcIsoEE",title.c_str(),40,0,0.5);
  title=std::string("N-1 relHcIsoEE ")+post;
  relHcIsoEE_nmo=td.make<TH1D>("N-1 relHcIsoEE",title.c_str(),40,0,0.5);  
  title=std::string("N-1 sigiEtaiEtaEE ")+post;
  sigiEtaiEtaEE_nmo=td.make<TH1D>("N-1 sigiEtaiEtaEE",title.c_str(),60,0,0.06);  
  title=std::string("N-1 dphiEE ")+post;
  dphiEE_nmo=td.make<TH1D>("N-1 dphiEE",title.c_str(),75,-0.75,0.75);
  title=std::string("N-1 detaEE ")+post;
  detaEE_nmo=td.make<TH1D>("N-1 detaEE",title.c_str(),40,-0.02,0.02);  
  title=std::string("N-1 hOeEE ")+post;
  hOeEE_nmo=td.make<TH1D>("N-1 hOeEE",title.c_str(),40,0,0.2);

  title=std::string("N-1 HF iso e9e25 ")+post;
  e9e25_nmo=td.make<TH1D>("N-1 e9e25",title.c_str(),40,0,1);
  title=std::string("N-1 HF eldId var2d")+post;
  var2d_nmo=td.make<TH1D>("N-1 var2d",title.c_str(),40,0,1.5);  
  title=std::string("N-1 eCOREe9 ")+post;
  eCOREe9_nmo=td.make<TH1D>("N-1 eCOREe9",title.c_str(),40,0,1);  
  title=std::string("N-1 eSeL ")+post;
  eSeL_nmo=td.make<TH1D>("N-1 eSeL",title.c_str(),40,0,1.5);  

}

void HFZeeVBTF::HistPerDef::fill(pat::ElectronCollection::const_iterator ecalE,  
				 const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, 
				 const edm::ParameterSet& myPs, const float hf_2d_cut) {

  reco::Particle::LorentzVector Z(ecalE->p4());
  Z+=hfE.p4();

  mee   ->Fill(Z.M());

  if(hfE.p4().eta()>0)   meeHFP   ->Fill(Z.M());
  else                   meeHFM   ->Fill(Z.M());

  Yee   ->Fill(Z.Rapidity());
  Ptee  ->Fill(Z.pt());

  ec_eta->Fill(ecalE->eta());
  bool isEb(false);
  int isEe=1;// 0: EB;  1: EE
  if(fabs(ecalE->eta())<1.4442){
    isEb=true;
    isEe=0;
  }
  ec_phi->Fill(ecalE->phi());
  ec_pt ->Fill(ecalE->pt());
  hf_eta->Fill(hfE.eta());
  hf_phi->Fill(hfE.phi());
  hf_pt ->Fill(hfE.pt());
  
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

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID
  // here follows filling of histograms for selection variables when the other N-1 have passed
  // which cuts values are used? Depends on which parameter set you pass...
  std::vector<double> ebParam = myPs.getParameter< std::vector<double> >("barrel"); 
  std::vector<double> eeParam = myPs.getParameter< std::vector<double> >("endcap"); 



  //////////////////////////////////////////////////////////////////////////////////////
  // establish cuts 
  // gf: bring HF variables to PSet as well
  float e9e25_cut    = 0.90;
  // float eSeL_cut     = 0.50;
  // float eCOREe9_cut  = 0.70;
  // float hf_2d_cut    = 0.32; // loose 
  
  // stay super loos for now
  //float e9e25_cut    = 0.40;
  //float eSeL_cut     = 0.75;
  //float eCOREe9_cut  = 0.35;
  //float hf_2d_cut    = 0.20; // loose 
  //float hf_2d_cut    = 0.45; // mid tight
  //float hf_2d_cut    = 0.55; // tight

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
  combIsolation_cut[0] = ebParam[17]; if(verboseDebug) std::cout  << "\n\ncombIsolationEB: " << combIsolation_cut[0] << std::endl;
  
  trackRel03_cut[0]    = ebParam[11]; if(verboseDebug) std::cout  << "trackRel03EB: " << trackRel03_cut[0] << std::endl;
  ecalRel03_cut[0]     = ebParam[12]; if(verboseDebug) std::cout  << "ecalRel03EB: " << ecalRel03_cut[0] << std::endl;
  hcalRel03_cut[0]     = ebParam[13]; if(verboseDebug) std::cout  << "hcalRel03EB: " << hcalRel03_cut[0] << std::endl;
  
  // EB eleId variables
  sigiEtaiEta_cut[0]  = ebParam[1];  if(verboseDebug) std::cout  << "sigiEetaiEtaEB: " << sigiEtaiEta_cut[0] << std::endl;
  dPhi_cut[0]         = ebParam[2];  if(verboseDebug) std::cout  << "DphiEB: " << dPhi_cut[0] << std::endl;
  dEta_cut[0]         = ebParam[3];  if(verboseDebug) std::cout  << "DetaEB: " << dEta_cut[0] << std::endl;
  HoE_cut[0]          = ebParam[0];  if(verboseDebug) std::cout  << "HoEEB: " << HoE_cut[0] << std::endl;

  
  // EE isolation variables
  combIsolation_cut[1] = eeParam[17]; if(verboseDebug) std::cout  << "\n\ncombIsolationEE: " << combIsolation_cut[1] << std::endl;
  
  trackRel03_cut[1]    = eeParam[11]; if(verboseDebug) std::cout  << "trackRel03EE: " << trackRel03_cut[1] << std::endl;
  ecalRel03_cut[1]     = eeParam[12]; if(verboseDebug) std::cout  << "ecalRel03EE: " << ecalRel03_cut[1] << std::endl;
  hcalRel03_cut[1]     = eeParam[13]; if(verboseDebug) std::cout  << "hcalRel03EE: " << hcalRel03_cut[1] << std::endl;
  
  // EE eleId variables
  sigiEtaiEta_cut[1]  = eeParam[1];  if(verboseDebug) std::cout  << "sigiEetaiEtaEE: " << sigiEtaiEta_cut[1] << std::endl;
  dPhi_cut[1]          = eeParam[2];  if(verboseDebug) std::cout  << "DphiEE: " << dPhi_cut[1] << std::endl;
  dEta_cut[1]          = eeParam[3];  if(verboseDebug) std::cout  << "DetaEE: " << dEta_cut[1] << std::endl;
  HoE_cut[1]           = eeParam[0];  if(verboseDebug) std::cout  << "HoEEE: " << HoE_cut[1] << std::endl;
  
  
  //////////////////////////////////////////////////////////////////////////////////////
  // now get hold of the variables
  //hfE
  float e9e25   = hfshape.e9e25();
  float var2d   = hfshape.eCOREe9()-(hfshape.eSeL()*9./8.);
  float eSeL    = hfshape.eSeL();
  float eCOREe9 = hfshape.eCOREe9();

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
  // if (eSeL            < eSeL_cut )               {cut[12] = true; passedCuts++;/*std::cout<<"c12 ";*/}
  // if (eCOREe9         < eCOREe9_cut)             {cut[13] = true; passedCuts++;/*std::cout<<"c13 ";*/}
  /*std::cout<<std::endl;*/
  short numCuts = 7;


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

  if(passedCuts==numCuts || ( (passedCuts==(numCuts-1)) && (!cut[11])) ) {
     var2d_nmo -> Fill(var2d);  eSeL_nmo-> Fill(eSeL); eCOREe9_nmo -> Fill(eCOREe9);
  }

}// end of fill()

static const double hf_2d_loose = 0.32;
// this value for tight is lower than CMS AN-2009/106 to accont for S/L miscalib in data
// see figure 15 therein: turn-on is sharp
static const double hf_2d_tight = 0.45;


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
  ecalID_=iConfig.getParameter<std::string>("ECALid");
  dolog_=iConfig.getParameter<bool>("DoLog");
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  
  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    w.add(*it);//      std::cout << "adding: " << (*it) << std::endl;
    }



  edm::Service<TFileService> fs;
  hists.nelec=fs->make<TH1D>("nelec","N_Elec",10,-0.5,9.5);
  hists.nhf=fs  ->make<TH1D>("nhf","N_HF",10,-0.5,9.5);
  //hists.base.book(fs->mkdir("base"),"(base)");
  hists.mee90none.book(fs->mkdir("mee90none"),"(mee90none)");
  hists.zMass.book(fs->mkdir("zMass"),"(zMass)");
  hists.mee90loose.book(fs->mkdir("mee90loose"),"(90,loose)");
  hists.mee80loose.book(fs->mkdir("mee80loose"),"(80,loose)");
  hists.mee70loose.book(fs->mkdir("mee70loose"),"(70,loose)");  
  hists.mee60loose.book(fs->mkdir("mee60loose"),"(60,loose)");  
  hists.mee90tight.book(fs->mkdir("mee90tight"),"(90,tight)");  
  hists.mee80tight.book(fs->mkdir("mee80tight"),"(80,tight)");  
  hists.mee70tight.book(fs->mkdir("mee70tight"),"(70,tight)");  
  hists.mee60tight.book(fs->mkdir("mee60tight"),"(60,tight)");  
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
void
HFZeeVBTF::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  //std::cout << "ciao! Entering analyze" << std::endl;

  edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
  edm::Handle<reco::SuperClusterCollection> SuperClusters;
  edm::Handle<reco::HFEMClusterShapeAssociationCollection> ClusterAssociation;
  
  iEvent.getByLabel("hfEMClusters",SuperClusters);
  iEvent.getByLabel("hfEMClusters",ClusterAssociation);
  iEvent.getByLabel("hfRecoEcalCandidate",HFElectrons);
  
  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    std::cout << "No electrons found in this event with tag " 
	      << "patElectrons" << std::endl;
    return ; // RETURN if no elecs in the event
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
      // hists.base.fill(i,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_);//gf: for each Z definition, N-1 control plots could be filled here too
      ;
    }
    
    // ECAL acceptance cut
    if (ecalE==pElecs.end() && w.doesElePass( i->electronID(ecalID_) )  // 3: passes electron iso/Id/convRej depending on configuration; ecalID_ defines working point
    	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)) {             // ECAL acceptance cut
      ecalE=i;
    }// if
    
  }// end loop over electrons
  
  hists.nelec->Fill(pElecs.size());
  hists.nhf  ->Fill(HFElectrons->size());

  if ( (pElecs.size()>0) &&  (HFElectrons->size()>0) && (ecalE!=pElecs.end())){
    if(dolog_) std::cout << "size ele coll:" << pElecs.size() << "  size HF coll: " << HFElectrons->size() << std::endl;
    // add here histograms w/o request of Z mass or eleId
    // mee90none is at the very beginning: one HF supercluste and one GSF w/ ID and iso
    const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0];
    reco::SuperClusterRef hfclusRef=hfE.superCluster();
    const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
    const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;
    hists.mee90none.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,0);//gf: for each Z definition, N-1 control plots could be filled here too
  }

  // gf: up to this point requirements are:
  // 1) at least one ECAL electron; ecalE is ==first== ECAL electron with SC within acceptance
  //    && electron ID and Isolation only
  // 2) at least one HF electrons
  // than try match combinatorics
  
  if (ecalE!=pElecs.end() && HFElectrons->size()>0) {
     const reco::RecoEcalCandidate& hfE=(*HFElectrons)[0]; // choose first HF candidate
     reco::SuperClusterRef hfclusRef=hfE.superCluster();
     const reco::HFEMClusterShapeRef hfclusShapeRef=(*ClusterAssociation).find(hfclusRef)->val;
     const reco::HFEMClusterShape& hfshape=*hfclusShapeRef;
    
    // make a Z
    reco::Particle::LorentzVector Z(ecalE->p4());
    Z+=hfE.p4();

    // gf: 3) loose mass cut and Z min pt
    // ==> these are the CANDIATES
    if (Z.M()>40 && Z.M()<130 && ecalE->pt()>15 && hfE.pt()>15) {
      double var2d=hfshape.eCOREe9()-(hfshape.eSeL()*1.125);
      double eta_det=ecalE->superCluster().get()->eta();
      if (dolog_) {
	std::cout << "------------------------------------------------------" << std::endl;
	
	std::cout << "Candidate!: event " << iEvent.id().event() << " run " << iEvent.id().run();
	
	std::cout << "  m_ee=" << Z.M() << " Y_ee=" << Z.Rapidity() << " pT_ee=" << Z.pt() <<  std::endl;
	std::cout << " file: " << currentFile_ << std::endl;
	
	// https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID#Electron_ID_Implementation_in_Re
	std::cout << "ECAL (pt,eta,phi, eta_det) : " << ecalE->pt() << ", " << ecalE->eta() << ", " << ecalE->phi() << " , (" << eta_det << ")" << std::endl;
	//DELETE THIS std::cout << "ECAL (pt,tkIso,phi, eta_det) : " << ecalE->pt() << ", " << ecalE->trackIso() << ", " << ecalE->phi() << " , (" << eta_det << ")" << std::endl;
	// gf: pull up here all the ECALId variables as descrived in twiki here above 
	std::cout << "\t ECAL ele // REL ISOLATION -  trackRel03: " << ecalE->dr03TkSumPt()/ecalE->p4().Pt() 
		  << "\t ecalRel03: "  << ecalE->dr03EcalRecHitSumEt()/ecalE->p4().Pt() 
		  << "\t hcalRel03: "  <<  ecalE->dr03HcalTowerSumEt()/ecalE->p4().Pt() 
		  << "\n\t Electron ID - sigIetaIeta: " << ecalE->scSigmaIEtaIEta()
		  << "\t Deta: " << ecalE->deltaEtaSuperClusterTrackAtVtx()
		  << "\t Dphi: " << ecalE->deltaPhiSuperClusterTrackAtVtx()
		  << "\t H/E: "  << ecalE->hadronicOverEm()
		  << std::endl;

	std::cout << "  ID (60,70,80,90) : " 
		  << ecalE->electronID("simpleEleId60relIso") << " " 
		  << ecalE->electronID("simpleEleId70relIso") << " " 
		  << ecalE->electronID("simpleEleId80relIso") << " " 
		  << ecalE->electronID("simpleEleId90relIso") << " " 
		  << std::endl;
	
	
	std::cout << "  HF (pt, eta, phi): " << hfE.pt() << ", " << hfE.eta() << ", " << hfE.phi() << std::endl;
	std::cout << "  Shape (S/L, ec/e9, e1/e9, e9/e25, 2d) : " 
		  << hfshape.eSeL() << " "
		  << hfshape.eCOREe9() << " "
		  << hfshape.e1x1()/hfshape.e3x3() << " "
		  << hfshape.e9e25() << " "
		  << var2d << " "
		  << std::endl;

	std::cout << "\n Triggers fired: ";
	
	Handle<TriggerResults> hltResults ;
	iEvent.getByLabel( edm::InputTag("TriggerResults","","HLT"), hltResults ) ;
      
	if (!init_) {
	  init_=true;
	  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltResults);
	  HLT_Names=triggerNames.triggerNames();
	}
	
	int nf=0;
	for (unsigned int i=0; i<HLT_Names.size(); i++) {
	  if ( hltResults->accept(i) ) {
	    if ((nf % 5)==0) std::cout << "\n    ";
	    std::cout << HLT_Names.at(i) << " ";
	    nf++;
	  }
	}
	std::cout << std::endl;
	
	std::cout << "======================================================" << std::endl;
      }// 4) here require electronId to have passed, for different working points
      //  5) and require HF eleID loose or tight
      hists.zMass.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,0);//gf: for each Z definition, N-1 control plots could be filled here too
      if (var2d>=hf_2d_loose) {
	if (w.doesElePass( ecalE->electronID("simpleEleId90relIso") )) hists.mee90loose.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,hf_2d_loose);
	if (w.doesElePass( ecalE->electronID("simpleEleId80relIso") )) hists.mee80loose.fill(ecalE,hfE,hfshape,robust80relIsoEleIDCutsV04_ps_,hf_2d_loose);
	if (w.doesElePass( ecalE->electronID("simpleEleId70relIso") )) hists.mee70loose.fill(ecalE,hfE,hfshape,robust70relIsoEleIDCutsV04_ps_,hf_2d_loose);
	if (w.doesElePass( ecalE->electronID("simpleEleId60relIso") )) hists.mee60loose.fill(ecalE,hfE,hfshape,robust60relIsoEleIDCutsV04_ps_,hf_2d_loose);
      }
      if (var2d>=hf_2d_tight) {
	if (w.doesElePass( ecalE->electronID("simpleEleId90relIso") ) ) hists.mee90tight.fill(ecalE,hfE,hfshape,robust90relIsoEleIDCutsV04_ps_,hf_2d_tight);
	if (w.doesElePass( ecalE->electronID("simpleEleId80relIso") ) )
	  { hists.mee80tight.fill(ecalE,hfE,hfshape,robust80relIsoEleIDCutsV04_ps_,hf_2d_tight); std::cout << "passed simpleEleId80relIso - HFtight" << std::endl;}
	if (w.doesElePass( ecalE->electronID("simpleEleId70relIso") ) ) 
	  {hists.mee70tight.fill(ecalE,hfE,hfshape,robust70relIsoEleIDCutsV04_ps_,hf_2d_tight);std::cout << "passed simpleEleI70relIso - HFtight" << std::endl;}	
	if (w.doesElePass( ecalE->electronID("simpleEleId60relIso") ) ) hists.mee60tight.fill(ecalE,hfE,hfshape,robust60relIsoEleIDCutsV04_ps_,hf_2d_tight);
      }//gf: why are electronID all treated the same?
      
    }// if mass is in Z window

  }// if the is at least on ECAL and one HF candidate

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
