// -*- C++ -*-
//
// Package:    HiggsGenAnalysis
// Class:      HiggsGenAnalysis
// 
/**\class HiggsGenAnalysis HiggsGenAnalysis.cc HeavyNu/AnalysisModules/src/HiggsGenAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bryan Dahmes
//         Created:  Wed Sep 22 04:49:56 CDT 2010
// $Id: HiggsGenAnalysis.cc,v 1.2 2012/03/12 19:02:14 afinkel Exp $
//
// Edited by:   ALexey Finkel
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Candidate/interface/Particle.h"

//
// class declaration
//

class HiggsGenAnalysis : public edm::EDFilter {
public:
  explicit HiggsGenAnalysis(const edm::ParameterSet&);
  ~HiggsGenAnalysis();
  
  struct zboson {
      int mode; // 11, 13
      double mass;
      double l1pt, l2pt ; 
      double l1eta, l2eta ; 
      reco::Particle::LorentzVector PZ, Pl1, Pl2;
  } ;

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  bool validZmumu(zboson theZ, bool offshell=false) ; 
  bool validZee(zboson theZ, bool offshell=false) ;
  bool recZ(zboson theZ);
  void Angular(zboson z1, zboson z2, reco::Particle::LorentzVector H);
  TLorentzVector TLV(reco::Particle::LorentzVector);

  bool farElectronFilter_ ; 

  // ----------member data ---------------------------
  struct HistStruct {
      TH1D* h_higgsY, *h_higgsY_mm, *h_higgsY_ee, *h_higgsY_eFe, *mu1pt, *mu2pt, *mupt, *el1pt, *el2pt, *elpt, *Phi1, *Phi, *CosTheta1, *CosTheta2, *CosTheta0, *HiggsRestP2, *HiggsP2 ; 
  } hists ; 
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HiggsGenAnalysis::HiggsGenAnalysis(const edm::ParameterSet& iConfig) { 
    edm::Service<TFileService> fs;

    farElectronFilter_ = iConfig.getParameter<bool> ("filterFarElectronsOnly") ;

    hists.h_higgsY     = fs->make<TH1D > ("h_higgsY","Higgs rapidity", 100,-5.,5.) ;
    hists.h_higgsY_mm  = fs->make<TH1D > ("h_higgsY_mm","Higgs rapidity", 100,-5.,5.) ; 
    hists.h_higgsY_ee  = fs->make<TH1D > ("h_higgsY_ee","Higgs rapidity", 100,-5.,5.) ; 
    hists.h_higgsY_eFe = fs->make<TH1D > ("h_higgsY_eFe","Higgs rapidity", 100,-5.,5.) ; 
     
    hists.mu1pt = fs->make<TH1D >("mu1pt","Mu1 pt",50,0.,100.);    
    hists.mu2pt = fs->make<TH1D >("mu2pt","Mu2 pt",50,0.,100.);
    hists.mupt = fs->make<TH1D >("mupt","Muon pt",50,0.,100.);
    hists.el1pt = fs->make<TH1D >("el1pt","El1 pt",50,0.,100.);
    hists.el2pt = fs->make<TH1D >("el2pt","El2 pt",50,0.,100.);
    hists.elpt = fs->make<TH1D >("elpt","Electron pt",50,0.,100.);
    
    hists.Phi = fs->make<TH1D>("Phi","#Phi",20, -3.2,3.2);
    hists.Phi1 = fs->make<TH1D>("Phi1","#Phi_{1}",20, -3.2,3.2);
    hists.CosTheta0 = fs->make<TH1D>("CostTheta0","Cos(#theta_{0})",50,-1.,1.);
    hists.CosTheta1 = fs->make<TH1D>("CostTheta1","Cos(#theta_{1})",20,-1.,1.);
    hists.CosTheta2 = fs->make<TH1D>("CostTheta2","Cos(#theta_{2})",20,-1.,1.);
    hists.HiggsRestP2 = fs->make<TH1D>("HiggsRestP2","H rest frame p^2",100,0.,150);
    hists.HiggsP2 = fs->make<TH1D>("HiggsP2","H lab frame p^2",100,0.,1000);
}


HiggsGenAnalysis::~HiggsGenAnalysis()
{

}


//
// member functions
//
bool
HiggsGenAnalysis::validZmumu(zboson theZ, bool offshell) {

    double minLep1pt = 20. ; 
    double minLep2pt = 10. ; 
    if ( offshell ) {
        minLep1pt = 5. ; minLep2pt = 5. ; 
    }
        
    double maxPt = std::max(theZ.l1pt,theZ.l2pt);
    double minPt = std::min(theZ.l1pt,theZ.l2pt);

    if ( maxPt > minLep1pt && minPt > minLep2pt ) {
        if ( (fabs(theZ.l1eta) < 2.1 || fabs(theZ.l2eta) < 2.1) &&
             (fabs(theZ.l1eta) < 2.4 && fabs(theZ.l2eta) < 2.4) )
            return true ;
    }
    return false ;
}

bool
HiggsGenAnalysis::validZee(zboson theZ, bool offshell) {
    
    double minLep1pt = 20. ; 
    double minLep2pt = 10. ; 
    if ( offshell ) {
        minLep1pt = 7. ; minLep2pt = 7. ; 
    }
        
    double maxPt = std::max(theZ.l1pt,theZ.l2pt);
    double minPt = std::min(theZ.l1pt,theZ.l2pt);

    if ( maxPt > minLep1pt && minPt > minLep2pt ) {
      if ( fabs(theZ.l1eta) < 2.5 && fabs(theZ.l2eta) < 2.5 ) return true ;
    }
    return false ;
}

bool HiggsGenAnalysis::recZ(zboson theZ)//this checks if Z daughters could be reconstructed
{
	double minMupt = 5.;
	double minElpt = 7.;
	double maxMuEta = 2.4;
	double maxElEta = 5.;
	
	if (theZ.mode==13)
		return ((fabs(theZ.l1eta)<maxMuEta) && (fabs(theZ.l2eta)<maxMuEta) && (theZ.l1pt>minMupt) && (theZ.l2pt>minMupt));
		
	if (theZ.mode==11)
		return ((fabs(theZ.l1eta)<maxElEta) && (fabs(theZ.l2eta)<maxElEta) && (theZ.l1pt>minElpt) && (theZ.l2pt>minElpt));
		
	return false;	
}

void HiggsGenAnalysis::Angular(zboson z1, zboson z2, reco::Particle::LorentzVector H0)
{
	TLorentzVector Z1, Z2, l11, l12, l21, l22, H, Prot(0,0,1,1);
	Z1 = TLV(z1.PZ);
	l11 = TLV(z1.Pl1);
	l12 = TLV(z1.Pl2);
	Z2 = TLV(z2.PZ);
	l21 = TLV(z2.Pl1);
	l22 = TLV(z2.Pl2);
	H = TLV(H0);
	
	//test to see that boost is working correctly
	TLorentzVector Hclone = H;
	hists.HiggsP2->Fill(Hclone.E());
	Hclone.Boost(-H.BoostVector());
	hists.HiggsRestP2->Fill(Hclone.E());
	
	//find theta-star and phi from Z1 in Higgs frame
	TLorentzVector Z1clone = Z1;
	Hclone = H;
	Prot.Boost(-H.BoostVector());
	Z1clone.Boost(-H.BoostVector());
	hists.Phi1->Fill(Z1clone.Phi()+Hclone.Phi());
	hists.CosTheta0->Fill(Z1clone.Vect()*Prot.Vect()/Z1clone.P()/Prot.P());
	Z1clone.Boost(H.BoostVector());
	
	//go to Z1 frame to get theta1
	TLorentzVector Z2clone = Z2;
	TLorentzVector l11clone = l11;
	Hclone = H;
	Z2clone.Boost(-Z1.BoostVector());
	Hclone.Boost(-Z1.BoostVector());
	l11clone.Boost(-Z1.BoostVector());
	
	double cosTh1 = l11clone.Vect()*Hclone.Vect()/l11clone.P()/Hclone.P();
	hists.CosTheta1->Fill(cosTh1);
	
	//go to Z2 frame and get theta2
	Z1clone = Z1;
	Hclone = H;
	TLorentzVector l21clone = l21;
	Z1clone.Boost(-Z2.BoostVector());
	l21clone.Boost(-Z2.BoostVector());
	
	double cosTh2 = l21clone.Vect()*Z1clone.Vect()/l21clone.P()/Z1clone.P();
	hists.CosTheta2->Fill(cosTh2);
	
	//now try to get the phi angle (between Z decay planes)
	l11clone=l11;
	l21clone=l21;
	Z1clone=Z1;
	Z2clone=Z2;
	l11clone.Boost(-H.BoostVector());
	l21clone.Boost(-H.BoostVector());
	Z1clone.Boost(-H.BoostVector());
	Z2clone.Boost(-H.BoostVector());
	
	TLorentzVector l1perp = l11clone-(l11clone.Vect()*Z1clone.Vect())/(Z1clone.Vect()*Z1clone.Vect())*Z1clone;
	TLorentzVector l2perp = l21clone-(l21clone.Vect()*Z2clone.Vect())/(Z2clone.Vect()*Z2clone.Vect())*Z2clone;
	
	double delPhi = acos( l1perp.Vect()*l2perp.Vect()/(l1perp.P()*l2perp.P()) );
	delPhi *= (l11clone.Phi()>l21clone.Phi())? 1: -1;
	
	hists.Phi->Fill( delPhi );
	
}

TLorentzVector HiggsGenAnalysis::TLV(reco::Particle::LorentzVector V)
{
//converts a reco::Particle::LorentzVector into a TLorentzVector
	return TLorentzVector(V.px(),V.py(),V.pz(),V.E());
}


// ------------ method called on each new Event  ------------
bool
HiggsGenAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> genInfo;
  int nZ = 0 ; bool passEvent = false ; 

  std::cout << "-------------------------------------------------------" << std::endl ; 
  if (iEvent.getByLabel("genParticles",genInfo)) {
      reco::GenParticleCollection::const_iterator iparticle ;
      zboson z1, z2 ; 
      for (iparticle=genInfo->begin(); iparticle!=genInfo->end(); iparticle++) {

          if ( iparticle->pdgId() == 25 && iparticle->numberOfDaughters() > 0 ) { // Higgs
              std::cout << "Higgs with mass " << iparticle->mass()
                        << " and rapidity " << iparticle->y()
                        << " decays to " << std::endl ;
              zboson newZ ;
              for (unsigned int idg=0; idg<iparticle->numberOfDaughters(); idg++) {
                  std::cout << "    " << iparticle->daughter(idg)->pdgId()
                            << " with mass " << iparticle->daughter(idg)->mass()
                            << " rapidity " << iparticle->daughter(idg)->y() 
                            << " and " << iparticle->daughter(idg)->numberOfDaughters()
                            << " daughters" << std::endl ;
                  if ( iparticle->daughter(idg)->pdgId() == 23 ) { // Z boson
                      int nele = 0 ; int nmu  = 0 ;
                      newZ.PZ = iparticle->daughter(idg)->p4();
                      for (unsigned int idg2=0; idg2<iparticle->daughter(idg)->numberOfDaughters(); idg2++) {
                          std::cout << iparticle->daughter(idg)->daughter(idg2)->pdgId() << " " ; 
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 11 ) {
                              nele++ ;
                              if ( nele == 1 ) {
                                  newZ.l1pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l1eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                                  newZ.Pl1	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                              } else {
                                  newZ.l2pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l2eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                                  newZ.Pl2	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                                  hists.el1pt->Fill( std::max(newZ.l1pt,newZ.l2pt) );
								  hists.el2pt->Fill( std::min(newZ.l1pt,newZ.l2pt) );
                              }
                              hists.elpt->Fill( iparticle->daughter(idg)->daughter(idg2)->pt() );
                          }
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 13 ) {
                              nmu++ ; 
                              if ( nmu == 1 ) {
                                  newZ.l1pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l1eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                              } else {
                                  newZ.l2pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l2eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                                  hists.mu1pt->Fill( std::max(newZ.l1pt,newZ.l2pt) );
                              	  hists.mu2pt->Fill( std::min(newZ.l1pt,newZ.l2pt) );
                              }
                              hists.mupt->Fill( iparticle->daughter(idg)->daughter(idg2)->pt() );
                          }
                      }
                      std::cout << std::endl ; 
                      if ( nele == 2 || nmu == 2 ) {
                          nZ++ ; 
                          newZ.mass = iparticle->daughter(idg)->mass() ;
                          newZ.mode = ( (nele > 0) ? 11 : 13 ) ;
                      }

                      if ( nZ == 1 ) { // 1st Z found
                          z1 = newZ ;
                      } else if ( nZ > 1 ) { // last Z found
                          double dm1 = fabs(z1.mass - 91. ) ; 
                          double dm  = fabs(newZ.mass - 91. ) ;
                          if ( dm1 < dm ) z2 = newZ ;
                          else {
                              z2 = z1 ;
                              z1 = newZ ;
                          }
                      }
                      Angular(z1,z2, iparticle->p4());
                  }
              }

              std::cout << "Event information: " << std::endl ;
              if ( nZ != 2 ) std::cout << "Not a valid decay: " << nZ << std::endl ;
              else {
                  std::cout << "Found a Higgs -> ZZ -> 4 lepton decay" << std::endl ;
                  std::cout << "Z1: " << z1.mass << " GeV decays to 2 " << ( z1.mode == 11 ? "electrons" : "muons" ) << std::endl ; 
                  std::cout << "primary:   pt = " << z1.l1pt << " GeV and eta = " << z1.l1eta << std::endl ; 
                  std::cout << "secondary: pt = " << z1.l2pt << " GeV and eta = " << z1.l2eta << std::endl ; 
                  std::cout << "Z2: " << z2.mass << " GeV decays to 2 " << ( z2.mode == 11 ? "electrons" : "muons" ) << std::endl ; 
                  std::cout << "primary:   pt = " << z2.l1pt << " GeV and eta = " << z2.l1eta << std::endl ; 
                  std::cout << "secondary: pt = " << z2.l2pt << " GeV and eta = " << z2.l2eta << std::endl ;
              }
              
              // Determine a category
              int zCat = -1 ;
              if ( nZ == 2 )
              {
                  zCat = 0 ; 
                  if ( z1.mode == 13 )
                  {
                      if ( z1.mass > 50. && z1.mass < 120. && validZmumu(z1) && z2.mass > 12. )
                      {
                          if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 1 ;
                          if ( z2.mode == 11 && validZee(z2,true) )   zCat = 2 ;
                      }
                  }
                  else if ( z1.mode == 11 )
                  {
                      if ( z1.mass > 50. && z1.mass < 120. && z2.mass > 12. ) 
                      {
                          if ( validZee(z1) ) // Requires both electrons in the tracker region
                          {
                              if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 3 ;
                              if ( z2.mode == 11 && validZee(z2,true) )   zCat = 4 ;
                          }
                          else
                          {
                              if ( (fabs(z1.l1eta) < 2.5 || fabs(z1.l2eta) < 2.5) &&
                                   (fabs(z1.l1eta) < 5.0 && fabs(z1.l2eta) < 5.0) )
                              {
                                  if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 5 ;
                                  if ( z2.mode == 11 && validZee(z2,true) )   zCat = 6 ;
                              }
                          }
                      }
                  }
              }
              
              if ( zCat >= 0 ) {
		passEvent = true ; 

		if ( recZ(z1) && recZ(z2) )	hists.h_higgsY->Fill( iparticle->y() ) ;
		if ( zCat == 1 || zCat == 2 ) hists.h_higgsY_mm->Fill( iparticle->y() ) ;
		if ( zCat == 3 || zCat == 4 ) hists.h_higgsY_ee->Fill( iparticle->y() ) ;
		if ( zCat == 5 || zCat == 6 ) hists.h_higgsY_eFe->Fill( iparticle->y() ) ;
		if ( farElectronFilter_ ) { 
		  if ( z1.mode != 11 || z1.mass < 50. || z1.mass > 120. ) return false ; 
		  if ( fabs(z1.l1eta) < 2.5 && fabs(z1.l2eta) < 2.5 ) return false ; 
		  std::cout << "*** This event has an HF/NT electron candidate ***" << std::endl ; 
		}
              }
          }
      }
      
  }
  std::cout << "-------------------------------------------------------" << std::endl ; 

  return passEvent ; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiggsGenAnalysis::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsGenAnalysis::endJob() {

}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HiggsGenAnalysis);
