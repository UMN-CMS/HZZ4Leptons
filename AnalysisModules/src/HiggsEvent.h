/** -- C++ -- **/
#ifndef HiggsEvent_h_included
#define HiggsEvent_h_included

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

/** The purpose of this class is contain the key items for
    a HiggsEvent and provide a simple way to pass this information
    between sections of the code.
*/

class HiggsEvent {
public:

  HiggsEvent() ; 

  // void initialize(int mode);

  void regularize();
  void scaleMuE(double mufactor=1.0, double efactor=1.0) ; 
  int getZ1(double minElePt1, double minElePt2, double minMuPt1, double minMuPt2, double minMass, double minHFpt, edm::Handle< std::vector<reco::Vertex> >& vertices) ; 
  int getZ2(double minElePt, double minMuPt, double minMass, double minM4, edm::Handle< std::vector<reco::Vertex> >& vertices) ; 
  void calculate();
  void decayID(const reco::GenParticleCollection& gpc);

  bool isMC ;
  // mc_class=0 (something else), 1=ee, 2=mm, 3=tau tau
  int mc_class;

  std::vector<reco::Muon>              muCands ; 
  std::vector<reco::GsfElectron>       gsfCands ; 
  std::vector<reco::Photon>            ntCands ; 
  std::vector<reco::RecoEcalCandidate> hfCands ; 
  void SetMuonList(std::vector<reco::Muon> MuCandList){muCands = MuCandList;};
  void SetGsfElectronList(std::vector<reco::GsfElectron> GsfCandList){gsfCands=GsfCandList;};
  void SetPhotonList(std::vector<reco::Photon> PhotonCandList){ntCands = PhotonCandList;};
  void SetHFList(std::vector<reco::RecoEcalCandidate> HFCandList){hfCands = HFCandList;};
  
  reco::Muon              mu1, mu2, mu3, mu4, mu[4] ; 
  reco::GsfElectron       e1, e2, e3, e4, e[4] ; 
  reco::Photon            nt1, nt[1] ; 
  reco::RecoEcalCandidate hf1, hf[1] ; 

  int nMuons, nElectrons ; 

  int n_primary_vertex, n_pue;

  // pat::MET met1;

  // separately stored for JEC Uncertainty studies
  // (saves space and time not copying whole jet objects,
  //  particularly during the jet selection)
  //

  float MuScale, ElecScale;

  int cutlevel ; 
  double eventWgt ; 

  int Z1flavor, Z2flavor ; 
  std::pair<unsigned int, unsigned int> Z1idx, Z2idx ; 

  reco::Particle::LorentzVector vZ1, vZ2, vl1, vl2, vl3, vl4;
  reco::Particle::LorentzVector lv_evt;
  reco::Particle::LorentzVector vH;

  double HY, mH, mZ1, mZ2, YZ1, YZ2;
  double l1pt, l1eta, l2pt, l2eta, l3pt, l3eta, l4pt, l4eta ; 
  double ecalSumEt_1, ecalSumEt_2, ecalSumEt_3, ecalSumEt_4;
  double ecalIsoByGSF_1, ecalIsoByGSF_2, ecalIsoByGSF_3, ecalIsoByGSF_4;
  double tkSumPt_1,tkSumPt_2,tkSumPt_3,tkSumPt_4;
  double trackIso_1, trackIso_2, trackIso_3, trackIso_4;
  double hcalIso_1, hcalIso_2, hcalIso_3, hcalIso_4;
  double netIso_1, netIso_2, netIso_3, netIso_4;
  double netIso_12, netIso_13, netIso_14, netIso_23, netIso_24, netIso_34;
  double scTheta_1, scTheta_2, scTheta_3, scTheta_4;
  double e25Max_1, e25Max_2, e25Max_3, e25Max_4;
  double e15_1, e15_2, e15_3, e15_4;
  double e55_1, e55_2, e55_3, e55_4;
  double HoEM_1, HoEM_2, HoEM_3, HoEM_4;
  double sIeIe_1, sIeIe_2, sIeIe_3, sIeIe_4;
  double var2d, e9e25;  
  double dz_1, dz_2, dz_3, dz_4;
  //double HoEMBarrel_1, HoEMBarrel_2, HoEMBarrel_3, HoEMBarrel_4;
  //double HoEMEndcap_1, HoEMEndcap_2, HoEMEndcap_3, HoEMEndcap_4;
  double PfIso04_1, PfIso04_2, PfIso04_3, PfIso04_4;
  double PfIsoNHadron, PfIsoCHadron, PfIsoPhoton;
  
  double cosTheta0, cosTheta1, cosTheta2, Phi, Phi1;
  
  
  bool passCombinedIso;
  
 // double muIso_1, muIso_2, muIso_3, muIso_4;

};

#endif
