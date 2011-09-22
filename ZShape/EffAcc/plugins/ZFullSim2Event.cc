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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
//
// class decleration
//


class ZFullSim2Event : public edm::EDProducer {
public:
  explicit ZFullSim2Event(const edm::ParameterSet&);
  virtual ~ZFullSim2Event() { }

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
      
  // ----------member data ---------------------------
  edm::InputTag m_mcTag;
  bool m_fromPAT;
  double dRmatch_;

};


ZFullSim2Event::ZFullSim2Event(const edm::ParameterSet& iConfig) : 
  m_mcTag(iConfig.getUntrackedParameter<edm::InputTag>("mc",edm::InputTag("generator"))),
  m_fromPAT(iConfig.getUntrackedParameter<bool>("fromPAT",false)),
  dRmatch_(iConfig.getUntrackedParameter<double>("dRmatch",0.15))
{
  using namespace reco;

  //register your products
  produces<GenParticleCollection>("ZEventEle1");
  produces<GenParticleCollection>("ZEventEle3");
  produces<GenParticleCollection>("ZEventParticles");
}

void
ZFullSim2Event::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  iEvent.getByLabel(m_mcTag, EvtHandle ) ; 
  if (!EvtHandle.isValid()) return;
  
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

  /////////////////////////////////////////
  // second loop: aiming to type-1 electrons
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {
    
    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();
    
    // aiming at electrons  after FSR (status 1) which match mother electrons with status 3
    if (status == 1 && abs(pid) == myPid ){
      if (true) {
	
	HepMC::GenEvent::particle_const_iterator ipart=mcpart;
	bool gotMatch=false;
	
	const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();
		
	while (!gotMatch && vertex_!=0) {
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
		  gotMatch=true;
		}// end if there's matching
	    }// end loop on moters
	    
	    if (!gotMatch) {
	      bool foundOfInterest=false;
	      for( motherIt = vertex_->particles_in_const_begin(); motherIt != vertex_->particles_in_const_end(); motherIt++) {
		if (abs((*motherIt)->pdg_id())==myPid && (*motherIt)->status()!=3) {
		  vertex_=(*motherIt)->production_vertex();
		  foundOfInterest=true;
		  break;
		}	      
	      }
	      if (!foundOfInterest) vertex_=0;
	    }// end if there's mothers	    
	  }// end if production vertex found
	}
      } else {
	LorentzVector momentum;
	momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
	momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
        
	pZeeEle1->push_back( GenParticle(q,
					 momentum,
					 Point(0,0,0),
					 pid,
					 status,
					 true
					 )
			     );
      }
    }// end if electron with status 1
  }// end loop on particles, aiming at type-1 electrons 

  if (pZeeEle1->empty()) {
    std::cout << "No matches\n";
    Evt->print(cout);
  }

  

  // now, we prepare for matching to reconstructed objects...
  if (pZeeEle3->size()==2) { // two-and-only-two for this step

    GenParticle bestMatches[2];
    double qualityOfMatch[2];
    qualityOfMatch[0]=1e6;     qualityOfMatch[1]=1e6; 

    Particle::LorentzVector dilep= (*pZeeEle3)[0].p4()+(*pZeeEle3)[1].p4();

    bool dump=false;

    if (dilep.M()>60 && dump) {
      std::cout << dilep.M() << " GeV / ";
      std::cout << (*pZeeEle3)[0].pt() << " " << (*pZeeEle3)[0].eta() << " " <<  (*pZeeEle3)[0].phi() << " / ";
      std::cout << (*pZeeEle3)[1].pt() << " " << (*pZeeEle3)[1].eta() << " " <<  (*pZeeEle3)[1].phi() << std::endl;
    }

    if (m_fromPAT) {
   
      edm::Handle<pat::ElectronCollection> patElectron;
      iEvent.getByLabel("patElectrons", patElectron);
      if ( ! patElectron.isValid()) {
	std::cout << "No electrons found in this event with tag " 
		  << "patElectrons" << std::endl;
	return; // RETURN if no elecs in the event
      }
      
      const pat::ElectronCollection& pElecs = *(patElectron.product());
      pat::ElectronCollection::const_iterator i, ecalE=pElecs.end();
      
      for (i=pElecs.begin(); i!=pElecs.end(); i++) {
	double qf1=reco::deltaR(*i,(*pZeeEle3)[0]);
	double qf2=reco::deltaR(*i,(*pZeeEle3)[1]);
    
	//	std::cout << i->pt() << " " << i->eta() << " " <<  i->phi() << " : " << qf1 << " " << qf2 << std::endl;	
     
      } 
    } else {
      edm::Handle<reco::GsfElectronCollection> gsfElectrons;

      iEvent.getByLabel("gsfElectrons", gsfElectrons);
      if ( ! gsfElectrons.isValid()) {
	std::cout << "No electrons found in this event with tag " 
		  << "gsfElectrons" << std::endl;
	return; // RETURN if no elecs in the event
      }
      
      const reco::GsfElectronCollection& pElecs = *gsfElectrons.product();
      reco::GsfElectronCollection::const_iterator i, ecalE=pElecs.end();
      
      for (i=pElecs.begin(); i!=pElecs.end(); i++) {
	for (int k=0; k<2; k++) {
	  double qf=reco::deltaR(*i,(*pZeeEle3)[k]);

	  if (qf<qualityOfMatch[k]) {
	    bestMatches[k]=GenParticle(i->charge(),
				       i->p4(),
				       i->vertex(),
				       11,
				       100,
				       true);
	    qualityOfMatch[k]=qf;
	  }
	}
      }
    }

    // next try HF!
    edm::Handle<reco::RecoEcalCandidateCollection> HFElectronsH;
    iEvent.getByLabel("hfRecoEcalCandidate",HFElectronsH);

    if (!HFElectronsH.isValid()) {
      std::cerr << "Invalid HF Electrons Handle" << std::endl;
      return;
    }
    
    const reco::RecoEcalCandidateCollection& HFelecs=*HFElectronsH;
    reco::RecoEcalCandidateCollection::const_iterator j;
    for (j=HFelecs.begin(); j!=HFelecs.end(); j++) {
      for (int k=0; k<2; k++) {
	double qf=reco::deltaR(*j,(*pZeeEle3)[k]);
	
	if (qf<qualityOfMatch[k]) {
	  bestMatches[k]=GenParticle(j->charge(),
				     j->p4(),
				     j->vertex(),
				     11,
				     100,
				     true);
	  qualityOfMatch[k]=qf;
	}
      }
    }

    if (dilep.M()>60 && dump) {
      std::cout << bestMatches[0].pt() << " " << bestMatches[0].eta() << " " <<  bestMatches[0].phi() << " : " << qualityOfMatch[0] << " / ";
      std::cout << bestMatches[1].pt() << " " << bestMatches[1].eta() << " " <<  bestMatches[1].phi() << " : " << qualityOfMatch[1] << std::endl;
      std::cout << std::endl;

    }

    if (qualityOfMatch[0]<dRmatch_ && qualityOfMatch[1]<dRmatch_) {
      pZeeParticles->push_back(bestMatches[0]);
      pZeeParticles->push_back(bestMatches[1]);
    }

         
  } else {
    std::cout << "We have " << pZeeEle3->size() << std::endl;
  }

  if (pZeeParticles->empty()) {
    LorentzVector p4(1.0,0.0,1.1e4,sqrt(1.0+1.1e4*1.1e4));
    
    pZeeParticles->push_back(GenParticle(1,
					p4,
					Point(0,0,0),
					11,
					110,
					true));
    pZeeParticles->push_back(GenParticle(-1,
					p4,
					Point(0,0,0),
					11,
					110,
					true));
  }


  iEvent.put(pZeeEle1,"ZEventEle1");
  iEvent.put(pZeeEle3,"ZEventEle3");
  iEvent.put(pZeeParticles,"ZEventParticles");

}
	    

//define this as a plug-in
DEFINE_FWK_MODULE(ZFullSim2Event);
