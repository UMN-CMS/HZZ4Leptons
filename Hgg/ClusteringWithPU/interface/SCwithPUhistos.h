#ifndef PUHS_INC
#define PUHS_INC

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

// holds histograms for pile up studies

class SCwithPUhistos {

public:
  void Book(TFileDirectory& tfd);
  
  void FillSc(const reco::SuperClusterCollection::const_iterator sc1, HepMC::GenEvent::particle_const_iterator particle1,
	    const EcalRecHitCollection* ebRecHits, const EcalRecHitCollection* eeRecHits);
  
  void FillGamma(const reco::PhotonCollection::const_iterator pho1, HepMC::GenEvent::particle_const_iterator particle1);

  void FillH(const std::vector<reco::Photon> higgsPhotons, const std::vector<reco::Photon> higgsPhotons_trueVtx);
  
  //  void FillDouble(const reco::SuperClusterCollection::const_iterator sc1, const reco::SuperClusterCollection::const_iterator sc2,
  //		  HepMC::GenEvent::particle_const_iterator particle1, HepMC::GenEvent::particle_const_iterator particle2,
  //		  const EcalRecHitCollection* ebRecHits, const EcalRecHitCollection* eeRecHits);
  //  

  float etaTransformation(  float EtaParticle , float Zvertex);
  
  SCwithPUhistos() { }
  
 private:
  
  TH1F *h_scet_barl;// single-photon mult
  TH1F *h_scet_endc;// single-photon mult
  TH1F *h_EoverEtrue_barl;// single-photon mult
  TH1F *h_EoverEtrue_endc;// single-photon mult
  TH1F *h_E5x5overEtrue_barl;// single-photon mult
  TH1F *h_E5x5overEtrue_endc;// single-photon mult
  TH1F *h_E5x5R9overEtrue_barl;// single-photon mult
  TH1F *h_PhoER9overEtrue_barl;// single-photon mult
  TH1F *h_E5x5notR9overEtrue_barl;// single-photon mult
  TH1F *h_PhoEnotR9overEtrue_barl;// single-photon mult
  TH1F *h_PhoEoverEtrue_barl;// single-photon mult
  TH1F *h_PhoEoverEtrue_endc;// single-photon mult
  TH1F *h_E5x5R9overEtrue_endc;// single-photon mult
  TH1F *h_PhoER9overEtrue_endc;// single-photon mult
  TH1F *h_E5x5notR9overEtrue_endc;// single-photon mult
  TH1F *h_PhoEnotR9overEtrue_endc;// single-photon mult

  TH1F *h_mHiggs_EBEB; // diphoton
  TH1F *h_mHiggs_EBEE; // diphoton
  TH1F *h_mHiggs_EEEE; // diphoton
  TH1F *h_mHiggs_EBEB_trueVtx; // diphoton
  TH1F *h_mHiggs_EBEE_trueVtx; // diphoton
  TH1F *h_mHiggs_EEEE_trueVtx; // diphoton
  TH2F *h_E5x5overEtrueVsEphoEtrue_barl;// single-photon mult
  TH2F *h_E5x5overEtrueVsEphoEtrue_endc;// single-photon mult
  TH1F *h_phiWidth;// single-photon mult
  TH2F *h_phiWidthVsE;// single-photon mult
  TH1F *h_phiWidth_endc;// single-photon mult
  TH2F *h_phiWidthVsE_endc;// single-photon mult
  TH1F *h_phiSize;// single-photon mult
  TH2F *h_phiSizeVsE;// single-photon mult
  TH2F *h_phiSizeVsEt;// single-photon mult
  TH1F *h_phiSize_endc;// single-photon mult
  TH2F *h_phiSizeVsE_endc;// single-photon mult
  TH2F *h_phiSizeVsEt_endc;// single-photon mult
  TH1F *h_phiShape_barl; // single-photon mult
  TH1F *h_absPhiShape_barl; // single-photon mult
  TH1F *h_phiShape_endc; // single-photon mult
  TH1F *h_absPhiShape_endc; // single-photon mult
  TH2F *h_phiShapeVsE_barl; // single-photon mult
  TH2F *h_absPhiShapeVsE_barl; // single-photon mult
  TH2F *h_phiShapeVsE_endc; // single-photon mult
  TH2F *h_absPhiShapeVsE_endc; // single-photon mult
  TH1F *h_etaShape_barl;// single-photon mult
  TH1F *h_etaShape_barlPLus;// single-photon mult
  TH1F *h_etaShape_barlMinus;// single-photon mult
  TH1F *h_etaShape_barlSymm;// single-photon mult
  TH2F *h_etaPhiShape_barl;// single-photon mult
  TH2F *h_etaPhiShape_barlPLus;// single-photon mult
  TH2F *h_etaPhiShape_barlMinus;// single-photon mult
  TH2F *h_etaPhiShape_barlSymm;// single-photon mult
  TH1F *h_phiBCminusPhiSeed_barl;// single-photon mult
  TH1F *h_etaBCminusEtaSeed_barl;// single-photon mult
  TH1F *h_absEtaBCminusAbsEtaSeed_barl;// single-photon mult
  TH1F *h_absEaBCminusAbsEtaSeed_barl;// single-photon mult
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;// single-photon mult
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1;// single-photon mult
  TH1F *h_phiBCminusPhiSeed_barlm1;// single-photon mult
  TH1F *h_etaBCminusEtaSeed_barlm1;// single-photon mult
  TH1F *h_absEtaBCminusAbsEtaSeed_barlm1;// single-photon mult
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2;// single-photon mult
  TH1F *h_phiBCminusPhiSeed_barlm2;// single-photon mult
  TH1F *h_etaBCminusEtaSeed_barlm2;// single-photon mult
  TH1F *h_absEtaBCminusAbsEtaSeed_barlm2;// single-photon mult
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3;// single-photon mult
  TH1F *h_phiBCminusPhiSeed_barlm3;// single-photon mult
  TH1F *h_etaBCminusEtaSeed_barlm3;// single-photon mult
  TH1F *h_absEtaBCminusAbsEtaSeed_barlm3;// single-photon mult
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4;// single-photon mult
  TH1F *h_phiBCminusPhiSeed_barlm4;// single-photon mult
  TH1F *h_etaBCminusEtaSeed_barlm4;// single-photon mult
  TH1F *h_absEtaBCminusAbsEtaSeed_barlm4;// single-photon mult
  TH1F *h_maxCryInDomino_barl;// single-photon mult
  TH2F *h_maxCryInDominoVsPhi_barl;// single-photon mult
  TH1F *h_maxCryInLocMax_barlPLus;// single-photon mult
  TH2F *h_maxCryInLocMaxVsPhi_barlPLus;// single-photon mult
  TH1F *h_maxCryInLocMax_barlMinus;// single-photon mult
  TH2F *h_maxCryInLocMaxVsPhi_barlMinus;// single-photon mult
  TH1F *h_maxCryInLocMax_barlSymm;// single-photon mult
  TH2F *h_maxCryInLocMaxVsPhi_barlSymm;// single-photon mult

//  double massWindowLow_,massWindowHigh_;
//  
//  TH1 *mZ_,*YZ_, *ptZ_,*ptZmon_,*polarizationZ_,*cosPolarizationZ_;
//  TH1 *YZmasscut_,*ptZmasscut_;
//  TH1 *YZmasscutTL_,*ptZmasscutTL_;
//  TH1 *YZTLmasscut_,*ptZTLmasscut_;
//  TH1 *YZTL_,*ptZTLmon_,*mZTL_,*ptZTL_;
//  bool booked_;
//  TH1 *e1eta_,*e2eta_,*e1phi_,*e2phi_,*e1pt_,*e2pt_,*eeta_,*ephi_,*hfeta_;
//  TH2 *mZ_Y_,*mZ_pt_,*pt_Y_,*e1eta_YZ_,*e2eta_YZ_,*e1eta_ptZ_,*e2eta_ptZ_,*e1eta_e2eta_,*YZTL_YZ_,*YZTL_YZ_matrix_;
//  TH1 *evt_PVz_, *evt_BSz_, *evt_MET_, *evt_PFMET_, *evt_TCMET_;
//  TH2* ptZTL_ptZ_, *ptZTL_ptZ_zoom_;
//  TH3* mZ_e2pt_e2eta_;
//
//  bool massFinals_;
//  struct {
//    TH2 *e1eta_,*e2eta_,*e1phi_,*e2phi_,*e1pt_,*e2pt_,*hfeta_;    
//  } mZ_binned_finals;

}; 

#endif
