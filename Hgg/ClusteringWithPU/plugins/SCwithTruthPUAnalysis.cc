// -*- C++ -*-
//
// Package:    SCwithTruthPUAnalysis
// Class:      SCwithTruthPUAnalysis
// 
/**\class SCwithTruthPUAnalysis SCwithTruthPUAnalysis.cc UserCode/SCwithTruthPUAnalysis/src/SCwithTruthPUAnalysis.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  G. Franzoni (previous instance was inspired to D. Futyans's)
//         Created:  Mon Jul 18 14:50:50 CEST 2011
// $Id: SCwithTruthPUAnalysis.cc,v 1.1 2011/07/18 14:50:29 franzoni Exp $
//
//


// system include files
#include <memory>

// user include files
#include <assert.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <TMath.h>
#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "Hgg/ClusteringWithPU/interface/SCwithPUhistos.h"
#include "Hgg/ClusteringWithPU/interface/Utils.h"


//#define PI    3.14159
//#define TWOPI 6.28318539
//
//#define R_ECAL           136.5  // radius of maximum containement
//#define Z_Endcap         328.0
//#define etaBarrelEndcap  1.479
//#define eeCrySize        2.5


//
// class declaration
//

class SCwithTruthPUAnalysis : public edm::EDAnalyzer {
public:
  explicit SCwithTruthPUAnalysis(const edm::ParameterSet&);
  ~SCwithTruthPUAnalysis();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  float  etaTransformation( float a, float b);

  // ----------member data ---------------------------

  bool useRawEnergy_;

  TH1F *h_nVtx;// global
  TH1F *h_dzVtx;// global

  
  struct HistSingleSet {
    
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(edm::Service<TFileService>& td,const std::string&,float xi);
    
    // fill all histos of the set with the two electron candidates
    void fill(const HepMC::GenEvent::particle_const_iterator truthParticle, 
	      const reco::SuperClusterCollection::const_iterator scIt,
	      const reco::VertexCollection * theRecVtxs  );
    
    TH1 *nVertices;
    
    TH1* particleEta,   *particlePhi,   *particleEt, *particleNumBC, *particleFrac;
    TH1* fractions[20];
    TH1* eOverTruth[20];
    TH1* eOverTruthOrig[20];
    
    float theXi;
    
  } theSingleHists;
  



};

//
// constants, enums and typedefs
//


//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}


//////////////////////////////////////////////////////////////////
// convert integer to string
std::string convertInt(int number)
{
  std::stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}



//
// static data member definitions
//

//
// constructors and destructor
//
SCwithTruthPUAnalysis::SCwithTruthPUAnalysis(const edm::ParameterSet& iConfig)
{
  useRawEnergy_ = iConfig.getParameter< bool >("useRawEnergy");
  if(useRawEnergy_) std::cout << "\n[SCwithTruthPUAnalysis] raw energy of SC will be used\n" << std::endl;
  else              std::cout << "\n[SCwithTruthPUAnalysis] corrected energy of SC will be used\n" << std::endl;
}

SCwithTruthPUAnalysis::~SCwithTruthPUAnalysis()
{
}

void SCwithTruthPUAnalysis::HistSingleSet::book(edm::Service<TFileService> &td, const std::string& post, float xi) {
  std::string title;

  title=std::string("num vertices ")+post+std::string("; num vertices");
  nVertices=td->make<TH1D>("num vertices","num vertices; num vertices",41,-0.5,40.5);

  title=std::string("eta_{SC} ")+post+std::string(";#eta_{SC}");
  particleEta=td->make<TH1D>("eta_{elLead}",title.c_str(),60,-3,3);
  title=std::string("phi_{elLead} ")+post+std::string(";#phi_{elLead}");
  particlePhi=td->make<TH1D>("phi_{elLead}",title.c_str(),60,-1*TMath::Pi(),TMath::Pi());
  title=std::string("et_{elLead} ")+post+std::string(";Et_{elLead} [GeV]");
  particleEt=td->make<TH1D>("et_{elLead}",title.c_str(),120,0,120);
  title=std::string("frac_{elLead} ")+post+std::string(";frac_{elLead}");
  particleFrac=td->make<TH1D>("frac_{elLead}",title.c_str(),11000,0,1.1);


  for(uint vertex=1; vertex<=20; vertex++){
    title=std::string("E_{SC} / E_{true} num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))+post+std::string("; E_{SC} / E_{true}");
    eOverTruth[vertex-1]=td->make<TH1D>( (std::string("E_{SC} / E_{true} num. vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),90,0.75,1.2);
    title=std::string("E_{SC} / E_{true} (orig.) num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))+post+std::string("; orig: E_{SC} / E_{true}");
    eOverTruthOrig[vertex-1]=td->make<TH1D>( (std::string("E_{SC} / E_{true} num. vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),90,0.75,1.2);
    fractions[vertex-1]=td->make<TH1D>( (std::string("lead E_{dyn}/E_{SC} num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),11000,0,1.1);
  } // loop over PU bins
  
  // the proportionality parameter for dynamic superclustering 
  theXi = xi;
}


void SCwithTruthPUAnalysis::HistSingleSet::fill(const HepMC::GenEvent::particle_const_iterator truthParticle, 
						const reco::SuperClusterCollection::const_iterator scIt,
						const reco::VertexCollection * theRecVtxs  )
//void SCwithTruthPUAnalysis::HistSet::fill(const pat::ElectronCollection::const_iterator leadEle, const pat::ElectronCollection::const_iterator subEle,
//				 const reco::VertexCollection * theRecVtxs)
{
  // float fraction = removePU(leadEle,theXi);
  // why is this commented out? Ah, needs  be imported... what about  Utils.h  ?

  particleEta -> Fill(scIt->eta());
  particlePhi -> Fill(scIt->phi());
  particleEt  -> Fill(scIt->energy()  / cosh( scIt->eta()) );
  particleNumBC->Fill( scIt->clustersSize() ); 
  
  float particleFrac = 1.;
  //particleFrac-> Fill(fraction);

  if( 0 ) std::cout  <<  " number of vertices: " << theRecVtxs->size() << std::endl;

  float theEOverTruth = (scIt->energy()  / cosh( scIt->eta() ) )  / (*truthParticle)->momentum().e() ;
 
  bool doLogFractions_=false;
  if( theRecVtxs->size()>0 && theRecVtxs->size()<=40) {
    fractions[( (theRecVtxs->size()-1) /2 )]      -> Fill( particleFrac );
    if(   particleFrac <1  ){
      eOverTruth[( (theRecVtxs->size()-1) /2 )]       -> Fill( theEOverTruth * particleFrac);
      eOverTruthOrig[ ( (theRecVtxs->size()-1) /2 ) ] -> Fill( theEOverTruth );
    }
  }
  else {
    fractions[ 19 ] -> Fill( particleFrac );
    //massPlots[ 19 ] -> Fill(Z.M());
    if( particleFrac <1  ){
      eOverTruth[ 19 ]      -> Fill( theEOverTruth * particleFrac);
      eOverTruthOrig[ 19 ]  -> Fill( theEOverTruth );
    }
  }
  //  massVsVertex->Fill( theRecVtxs->size(), Z.M() );

}// end of fill()





//
// member functions
//

// ------------ method called to for each event  ------------
void
SCwithTruthPUAnalysis::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  Handle<vector<SuperCluster> > barrelSCHandle;
  ev.getByLabel("correctedHybridSuperClusters","",barrelSCHandle);
  const SuperClusterCollection * barrelSCCollection = barrelSCHandle.product();

  Handle<vector<SuperCluster> > endcapSCHandle;
  ev.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",endcapSCHandle);
  const SuperClusterCollection * endcapSCCollection = endcapSCHandle.product();

  // do i still need this? remove!!
  Handle<PhotonCollection> photonHandle;
  ev.getByLabel("photons","", photonHandle);
  const PhotonCollection * photonCollection = photonHandle.product();

  Handle<VertexCollection> vertexHandle;
  ev.getByLabel("offlinePrimaryVerticesWithBS", vertexHandle);
  const VertexCollection * vertexCollection = vertexHandle.product();
  math::XYZPoint recoVtx(0.,0.,0.);
  if (vertexCollection->size()>0) recoVtx = vertexCollection->begin()->position();
  
  //  std::cout << "++ SCwithTruthPUAnalysis " << std::endl;  // GF: check heart beat

  // this histogram can be used as a counter of events 
  h_nVtx->Fill(vertexCollection->size());

  Handle< HepMCProduct > hepProd ;
  ev.getByLabel("generator",hepProd) ;
  const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();

  //  vector<Photon> higgsPhotons,higgsPhotons_trueVtx;
  int isEB=0; 
  int hasOneClus(0);   int hasTwoClus(0);   int hasThreeClus(0);   int hasMoreClus(0); 
  bool trueVtxFound=false;
  math::XYZPoint trueVtx(0.,0.,0.);

  // loop over MC-truth particles
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
    // match only to truth of electrons or photons
    if ( !( (fabs((*p)->pdg_id())==11 || (*p)->pdg_id()==22) && (*p)->status()==1 )  )  continue;
    
    HepMC::GenParticle* mother = 0;
    HepMC::GenParticle* mother2 = 0;
    if ( (*p)->production_vertex() )  {
      if ( (*p)->production_vertex()->particles_begin(HepMC::parents) !=
           (*p)->production_vertex()->particles_end(HepMC::parents))
	mother = *((*p)->production_vertex()->particles_begin(HepMC::parents));
    }
    if (mother!=0) mother2 = *(mother->production_vertex()->particles_begin(HepMC::parents));
    if (mother == 0 || (mother2!=0 && mother2->pdg_id()==25)) {

      if (!trueVtxFound && mother2!=0 && mother2->pdg_id()==25) {
	HepMC::ThreeVector vtx = mother2->production_vertex()->point3d();
	trueVtx = math::XYZPoint(vtx.x()/10.,vtx.y()/10.,vtx.z()/10.);
	trueVtxFound=true;
	h_dzVtx->Fill(recoVtx.z()-trueVtx.z());
      }

      float phi_true=(*p)->momentum().phi();
      float eta_true=(*p)->momentum().eta();
      float etaEcal_true = etaTransformation(eta_true, (*p)->production_vertex()->position().z()/10. );
      float et_true = (*p)->momentum().e()/cosh((*p)->momentum().eta());

      // Barrel SuperClusters (inside loop over MC particles)
      for(SuperClusterCollection::const_iterator scIt = barrelSCCollection->begin(); scIt != barrelSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())<1.4442 && scIt->energy()/cosh(scIt->eta())>20.) {

	  // perform SC-MCparticle matching
	  float deltaPhi = scIt->phi()-phi_true;
	  float deltaEta = scIt->eta()-etaEcal_true;
	  float phiWidth = scIt->phiWidth();          // covariance of cluster in phi
	                                              // we want to look at the absolute extension too... phi_MAX-phi_min 
	  // float energySC = scIt->energy();
	  float energySC       = scIt->rawEnergy();
	  float energySCtoCheck=0;

	  float phiSize = getPhiSize(scIt,true);
	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1


	    // do things here
	    
	  }// if SC matches photons
	  //COPY1//
	}// if EB and Et>20
      }// loop over superclusters

      // Endcap SuperClusters
      for(SuperClusterCollection::const_iterator scIt = endcapSCCollection->begin(); scIt != endcapSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())>1.566 && fabs(scIt->eta())<2.5 && scIt->energy()/cosh(scIt->eta())>20.) {
	  float deltaPhi = scIt->phi()-phi_true;
	  float deltaEta = scIt->eta()-etaEcal_true;
	  float phiWidth = scIt->phiWidth();
	  float energySC = scIt->rawEnergy();
	  float energySCtoCheck=0;

	  float phiSize = getPhiSize(scIt,false);

	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1


	    // do things here
	    

	  }// if SC matches one MC particle 
	}// if in EE and Et>20 GeV 
      }// loop over EE Sc's
      
      


      /*
      ///////////////////////////////////////////////// keep for now /////////////////////////////////////////////////
      //////////////////////////////////////////////
      // follows the photon part of the histograms
      // still within loop over MC particles
      // special treatment if current MCparticle is a photon
      // match also to photon object
      if ((*p)->pdg_id()==22) {
	for(PhotonCollection::const_iterator phoIt = photonCollection->begin(); phoIt != photonCollection->end(); phoIt++) {
	  float deltaPhi = phoIt->caloPosition().phi()-phi_true;
	  float deltaEta = phoIt->caloPosition().eta()-etaEcal_true;
	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1 && phoIt->energy()/cosh(phoIt->caloPosition().eta())>20.) {
	    
	    if (mother2!=0 && mother2->pdg_id()==25) { // specifically look out for HIGGses
	      if (phoIt->isEB() || phoIt->isEE()) {
		// use this variable to classify SC's associated to photons
		int numOfBC = phoIt->superCluster()->clustersSize();
		if      (numOfBC==1) hasOneClus++;
		else if (numOfBC==2) hasTwoClus++;
		else if (numOfBC==3) hasThreeClus++;
		else                 hasMoreClus++;
		//		higgsPhotons.push_back(*phoIt);
		//		Photon localPho = Photon(*phoIt);
		//		if (!trueVtxFound) cout << "Error: true vertex not found!" << endl;
		//		localPho.setVertex(trueVtx);
		//		higgsPhotons_trueVtx.push_back(localPho);
		//cout << phoIt->et() << " " << localPho.et() << " " << trueVtx << " " << recoVtx << endl;
	      }
	    } // if mother
	  }// if matching between SC and photon and ET_{pho} > 20
	}// loop over photons
      }// if gamma
      ///////////////////////////////////////////////// keep for now /////////////////////////////////////////////////
      */


    }// if it's MC-truth photon
  }// loop over all MC-truth particles
}

// ------------ method called once each job just before starting event loop  ------------
void 
SCwithTruthPUAnalysis::beginJob()
{
  edm::Service<TFileService> fs;
  std::cout << "[SCwithTruthPUAnalysis] making directories" << std::endl;
  
  TFileDirectory subDir=fs->mkdir("allCase");  
  h_nVtx = fs->make<TH1F>("h_nVtx","no. of primary vertices; num vertices reco",30,0.,30.);  // usig h_nVtx  histogram can be used as a counter of events 
  h_dzVtx = fs->make<TH1F>("h_dzVtx","delta_z for reconstructed PV w.r.t. true PV",200,-5.,5.);

  // add histos to the class, instead!
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCwithTruthPUAnalysis::endJob() {
}

// compute eta at detector, given physics eta and location of primary vertex
float SCwithTruthPUAnalysis::etaTransformation(  float EtaParticle , float Zvertex)  {

  //---Definitions for ECAL
  // const float R_ECAL           = 136.5;  // radius of maximum containement
  // const float Z_Endcap         = 328.0;
  // const float etaBarrelEndcap  = 1.479;

  //---ETA correction

  float Theta = 0.0  ;
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0)    Theta = Theta+PI ;
  float ETA =      - log(tan(0.5*Theta));

  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle);
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+PI ;
      ETA = - log(tan(0.5*Theta));
    }
  //---Return the result
  return ETA;
  //---end
}


#include "FWCore/Framework/interface/MakerMacros.h"  
DEFINE_FWK_MODULE( SCwithTruthPUAnalysis );
