// -*- C++ -*-
//
// Package:    ZEventProducer
// Class:      ZEventProducer
// 
/**\class ZEventProducer ZEventProducer.cc ZShape/ZEventProducer/src/ZEventProducer.cc

Description: <one line class summary>
Produces pair of electrons (GenParticleCollection) from Z decay as generated with Phytia
Implementation:
<Notes on implementation>
*/
//
// Original Author:  Giovanni FRANZONI
//         Created:  Mon Feb 18 21:18:39 CET 2008
// $Id: ZEventProducer.cc,v 1.5 2008/06/15 16:03:19 franzoni Exp $
//
//


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/GeometryVector/interface/Phi.h"
//
// class decleration
//


class ZEventProducer : public edm::EDProducer {
public:
  explicit ZEventProducer(const edm::ParameterSet&);
  ~ZEventProducer();

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  // ----------member data ---------------------------

  edm::InputTag m_srcTag;
  double deltaR_, minET_;
};


ZEventProducer::ZEventProducer(const edm::ParameterSet& iConfig) : 
    m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source"))),
    deltaR_(iConfig.getParameter<double>("deltaR")),
    minET_(iConfig.getParameter<double>("minET"))
{
  using namespace reco;

  //register your products
  produces<GenParticleCollection>("ZEventEle1");
  produces<GenParticleCollection>("ZEventEle3");
  produces<GenParticleCollection>("ZEventParticles");

}


ZEventProducer::~ZEventProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ZEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm; 
  using namespace std; 
  using namespace reco;
  
  std::auto_ptr<GenParticleCollection> pZeeEle1(  new GenParticleCollection() );
  std::auto_ptr<GenParticleCollection> pZeeEle3(  new GenParticleCollection() );
  std::auto_ptr<GenParticleCollection> pZeeParticles(  new GenParticleCollection() );

  pZeeEle3->clear();
  pZeeEle1->clear();
  pZeeParticles->clear();

  typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;
  // electric charge type
  int q=0;
  // Lorentz vector
  typedef int Charge;
  typedef math::XYZTLorentzVector LorentzVector;
  // Lorentz vector
  typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;
  // point in the space
  typedef math::XYZPoint Point;
  // point in the space
  typedef math::XYZVector Vector;
  

  Handle< HepMCProduct > EvtHandle ; 
  
  // find initial (unsmeared, unfiltered,...) HepMCProduct 
  iEvent.getByLabel(m_srcTag, EvtHandle ) ; 
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ; 
     
  int status;
  int pid;
  std::vector<int> ele1Barcodes;
  
  ////////////////////////
  // looking for electrons
  int myPid = 11;

  // first loop aiming to pre-FSR electrons
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {
    
    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();
    
    // https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookGenParticleCandidate
    // http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
    // status = 1 : an existing entry, which has not decayed or fragmented. This is the main class of entries, which represents the `final state' given by the generator.
    // status = 2 : decayed or fragmented entry (i.e. decayed particle or parton produced in shower.)
    // status = 3 : a documentation line, defined separately from the event history. This could include the two incoming reacting particles, etc
    if (status == 3 && abs(pid) == myPid ){
        
      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();
      
      ele1Barcodes.push_back( (*mcpart)->barcode() );
      
      LorentzVector momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
      
      pZeeEle3->push_back( GenParticle(q,
                                       momentum,
                                       Point(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() ),
                                       pid,
                                       status,
                                       true
                                       )
                           );
}
}// end of loop on particles aiming to type-3 electrons


  
    //GF
    // int num(0);
    

  /////////////////////////////////////////
  // second loop: aiming to type-1 electrons
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {
    
    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();

    // https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookGenParticleCandidate
    // http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
    // status = 1 : an existing entry, which has not decayed or fragmented. This is the main class of entries, which represents the `final state' given by the generator.
    // status = 2 : decayed or fragmented entry (i.e. decayed particle or parton produced in shower.)
    // status = 3 : a documentation line, defined separately from the event history. This could include the two incoming reacting particles, etc

    
    // aiming at electrons  after FSR (status 1) which match mother electrons with status 3
    if (status == 1 && abs(pid) == myPid ){
        
        const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();
        
        if ( vertex_ != 0 ) {
            size_t numberOfMothers = vertex_->particles_in_size();
            if ( numberOfMothers > 0 ) {
                HepMC::GenVertex::particles_in_const_iterator motherIt = vertex_->particles_in_const_begin();
                
                for( ; motherIt != vertex_->particles_in_const_end(); motherIt++) {
                    
                    if ( std::find(ele1Barcodes.begin() ,  ele1Barcodes.end(), (*motherIt)->barcode() ) != ele1Barcodes.end() )
                    {//GF
                        // std::cout << "they match: " << std::endl;
                        //std::cout << "barcode: " << (*motherIt)->barcode() << " py: " << (*motherIt)->momentum().y() << std::endl;
                        //std::cout << "barcode: " << (*mcpart)->barcode() << " py: "  <<  (*mcpart)->momentum().y() << std::endl;
                        //num++;//GF
                        
                        LorentzVector momentum;
                        momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
                        momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
                        
                        pZeeEle1->push_back( GenParticle(q,
                                                         momentum,
                                                         Point(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() ),
                                                         pid,
                                                         status,
                                                         true
                                                         )
                                             );
                        
                    }// end if there's matching
                }// end loop on moters
            }// end if there's mothers
        }// end if production vertex found
    }// end if electron with status 1
  }// end loop on particles, aiming at type-1 electrons 

  //GF
  //if (num!=2) std::cout << "no two matches" << std::endl;

////////////////////////////////////////
      // looping on all the (status 1) particles of the simulated event to collect those falling in a cone around status 1 electrons 
      for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
          mcpart != Evt->particles_end();
          ++ mcpart ) {
          
          // http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
          // status = 1 : an existing entry, which has not decayed or fragmented. This is the main class of entries, which represents the `final state' given by the generator.
          // status = 3 : a documentation line, defined separately from the event history. This could include the two incoming reacting particles, etc
          status = (*mcpart)->status();
          if       (status !=1) continue;
          
          pid    = (*mcpart)->pdg_id();
          
          if     ((*mcpart)->pdg_id() >0) q=1;
          else  q=-1;
          
      LorentzVector momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
      
      // kill very low PT stuff
      //if (momentum.Pt()<minET_) continue;
      
      float DR=9999;
      GenParticleCollection::const_iterator oneElectron;
      for (oneElectron  = pZeeEle1->begin();
           oneElectron != pZeeEle1->end();
           oneElectron ++)
      {
          float dphi = (  Geom::Phi<float>( momentum.phi())  -  Geom::Phi<float>( (*oneElectron).phi())  ).value() ;
          float deta = momentum.Eta()-(*oneElectron).eta();
          float dr   = pow(dphi,2)+pow(deta,2);
          //std::cout <<pZeeEle1->size() << " " <<momentum.e() << " "  << momentum.phi() << "  " <<(*oneElectron).phi()
          //<< "  " <<momentum.Eta() << "  " <<(*oneElectron).eta()  << "  DR: " << DR << std::endl;
          if (dr<DR) DR=dr;
      }
      
      DR=sqrt(DR);
      
      //GF
      //if (DR<0.000001) std::cout << "found myself" << std::endl;
      //GF
      
      if (DR>deltaR_) continue;
      
      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();
      pZeeParticles->push_back( GenParticle(q,
                                            momentum,
                                                Point(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() ),
                                            pid,
                                            status,
                                            true
                                            )
                                );
      
      
  }// end loop to store particles
  
  
  // reorder electrons in pT
  //float pt1 = evt_.elec(0).p4_.Pt();
  //float pt2 = evt_.elec(1).p4_.Pt();
  
  //  std::cout << "collection size: " << pZeeParticles->size() << std::endl;
  iEvent.put(pZeeEle1,"ZEventEle1");
  iEvent.put(pZeeEle3,"ZEventEle3");
  iEvent.put(pZeeParticles,"ZEventParticles");
  
  
}


// ------------ Method called once each job just before starting event loop  ------------
void ZEventProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZEventProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEventProducer);
