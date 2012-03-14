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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
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
  bool getZ1(double minElePt1, double minElePt2, double minMuPt1, double minMuPt2, double minMass) ; 
  bool getZ2(double minElePt, double minMuPt, double minMass, double minM4) ; 
  void calculate();
  void decayID(const reco::GenParticleCollection& gpc);

  bool isMC ;
  // mc_class=0 (something else), 1=ee, 2=mm, 3=tau tau
  int mc_class;

  std::vector<reco::Muon>              muCands ; 
  std::vector<reco::GsfElectron>       gsfCands ; 
  std::vector<reco::Photon>            ntCands ; 
  std::vector<reco::RecoEcalCandidate> hfCands ; 

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

  reco::Particle::LorentzVector vZ1, vZ2;
  reco::Particle::LorentzVector lv_evt;
  reco::Particle::LorentzVector vH;

    double mH, mZ1, mZ2 ; 

};

#endif
