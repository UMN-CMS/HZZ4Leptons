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
// $Id: HiggsGenAnalysis.cc,v 1.1 2011/06/10 17:41:44 mansj Exp $
//
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
  } ;

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  bool validZmumu(zboson theZ, bool offshell=false) ; 
  bool validZee(zboson theZ, bool offshell=false) ;

  // ----------member data ---------------------------
  struct HistStruct {
      TH1D* h_higgsY, *h_higgsY_mm, *h_higgsY_ee, *h_higgsY_eFe ; 
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

    hists.h_higgsY     = fs->make<TH1D > ("h_higgsY","Higgs rapidity", 100,-5.,5.) ;

    hists.h_higgsY_mm  = fs->make<TH1D > ("h_higgsY_mm","Higgs rapidity", 100,-5.,5.) ; 
    hists.h_higgsY_ee  = fs->make<TH1D > ("h_higgsY_ee","Higgs rapidity", 100,-5.,5.) ; 
    hists.h_higgsY_eFe = fs->make<TH1D > ("h_higgsY_eFe","Higgs rapidity", 100,-5.,5.) ; 
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
        
    double maxPt = (( theZ.l1pt > theZ.l1pt ) ? theZ.l1pt : theZ.l2pt ) ; 
    double minPt = (( theZ.l1pt > theZ.l1pt ) ? theZ.l2pt : theZ.l1pt ) ;
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
        
    double maxPt = (( theZ.l1pt > theZ.l1pt ) ? theZ.l1pt : theZ.l2pt ) ; 
    double minPt = (( theZ.l1pt > theZ.l1pt ) ? theZ.l2pt : theZ.l1pt ) ;
    if ( maxPt > minLep1pt && minPt > minLep2pt ) {
        if ( (fabs(theZ.l1eta) < 2.1 || fabs(theZ.l2eta) < 2.1) &&
             (fabs(theZ.l1eta) < 2.4 && fabs(theZ.l2eta) < 2.4) )
            return true ;
    }
    return false ;
}


// ------------ method called on each new Event  ------------
bool
HiggsGenAnalysis::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::GenParticleCollection> genInfo;
  if (iEvent.getByLabel("genParticles",genInfo)) {
      reco::GenParticleCollection::const_iterator iparticle ;
      zboson z1, z2 ; 
      int nZ = 0 ; 
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
                      for (unsigned int idg2=0; idg2<iparticle->daughter(idg)->numberOfDaughters(); idg2++) {
                          std::cout << iparticle->daughter(idg)->daughter(idg2)->pdgId() << " " ; 
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 11 ) {
                              nele++ ;
                              if ( nele == 1 ) {
                                  newZ.l1pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l1eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                              } else {
                                  newZ.l2pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l2eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                              }
                          }
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 13 ) {
                              nmu++ ; 
                              if ( nmu == 1 ) {
                                  newZ.l1pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l1eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                              } else {
                                  newZ.l2pt  = iparticle->daughter(idg)->daughter(idg2)->pt() ;
                                  newZ.l2eta = iparticle->daughter(idg)->daughter(idg2)->eta() ;
                              }
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
              if ( nZ == 2 ) {
                  zCat = 0 ; 
                  if ( z1.mode == 13 ) {
                      if ( z1.mass > 50. && z1.mass < 120. && validZmumu(z1) && z2.mass > 12. ) {
                          if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 1 ;
                          if ( z2.mode == 11 && validZee(z2,true) )   zCat = 2 ;
                      }
                  } else if ( z1.mode == 11 ) {
                      if ( z1.mass > 50. && z1.mass < 120. && z2.mass > 12. ) {
                          if ( validZee(z1) ) {
                              if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 3 ;
                              if ( z2.mode == 11 && validZee(z2,true) )   zCat = 4 ;
                          }
                          else {
                              if ( (fabs(z1.l1eta) < 2.5 || fabs(z1.l2eta) < 2.5) &&
                                   (fabs(z1.l1eta) < 5.0 && fabs(z1.l2eta) < 5.0) ) {
                                  if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 5 ;
                                  if ( z2.mode == 11 && validZee(z2,true) )   zCat = 6 ;
                              }
                          }
                      }
                  }
              }
              
              if ( zCat >= 0 ) {
                  hists.h_higgsY->Fill( iparticle->y() ) ;
                  if ( zCat == 1 || zCat == 2 ) hists.h_higgsY_mm->Fill( iparticle->y() ) ;
                  if ( zCat == 3 || zCat == 4 ) hists.h_higgsY_ee->Fill( iparticle->y() ) ;
                  if ( zCat == 5 || zCat == 6 ) hists.h_higgsY_eFe->Fill( iparticle->y() ) ;
              }
          }
      }
      
  }

  return false ; 
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
