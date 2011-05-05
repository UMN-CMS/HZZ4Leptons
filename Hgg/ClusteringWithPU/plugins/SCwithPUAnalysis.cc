// -*- C++ -*-
//
// Package:    SCwithPUAnalysis
// Class:      SCwithPUAnalysis
// 
/**\class SCwithPUAnalysis SCwithPUAnalysis.cc UserCode/SCwithPUAnalysis/src/SCwithPUAnalysis.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  David Futyan,40 4-B32,+41227671591,
//         Created:  Thu Dec  2 20:20:57 CET 2010
// $Id: SCwithPUAnalysis.cc,v 1.1 2011/05/05 10:02:43 franzoni Exp $
//
//

// analyzer from David Futyan - Imperial College
// imported for usage of UMN G. Franzoni, Y. Kubota, V. Rekovic

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

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

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

class SCwithPUAnalysis : public edm::EDAnalyzer {
public:
  explicit SCwithPUAnalysis(const edm::ParameterSet&);
  ~SCwithPUAnalysis();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  float  etaTransformation( float a, float b);

  // ----------member data ---------------------------

  bool useRawEnergy_;

  SCwithPUhistos allCase_;
  SCwithPUhistos oneCase_, twoCase_, threeCase_, moreCase_;

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
  TH1F *h_nVtx;// global
  TH1F *h_dzVtx;// global
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
SCwithPUAnalysis::SCwithPUAnalysis(const edm::ParameterSet& iConfig)
{
  useRawEnergy_ = iConfig.getParameter< bool >("useRawEnergy");
  if(useRawEnergy_) std::cout << "\n[SCwithPUAnalysis] raw energy of SC will be used\n" << std::endl;
  else              std::cout << "\n[SCwithPUAnalysis] corrected energy of SC will be used\n" << std::endl;
}


SCwithPUAnalysis::~SCwithPUAnalysis()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SCwithPUAnalysis::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
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

  Handle<PhotonCollection> photonHandle;
  ev.getByLabel("photons","", photonHandle);
  const PhotonCollection * photonCollection = photonHandle.product();

  Handle<VertexCollection> vertexHandle;
  ev.getByLabel("offlinePrimaryVerticesWithBS", vertexHandle);
  const VertexCollection * vertexCollection = vertexHandle.product();
  math::XYZPoint recoVtx(0.,0.,0.);
  if (vertexCollection->size()>0) recoVtx = vertexCollection->begin()->position();

  Handle<EcalRecHitCollection> ebRecHitsHandle;
  ev.getByLabel("ecalRecHit", "EcalRecHitsEB", ebRecHitsHandle);
  const EcalRecHitCollection * ebRecHits = ebRecHitsHandle.product();

  Handle<EcalRecHitCollection> eeRecHitsHandle;
  ev.getByLabel("ecalRecHit", "EcalRecHitsEE", eeRecHitsHandle);
  const EcalRecHitCollection * eeRecHits = eeRecHitsHandle.product();

  // this histogram can be used as a counter of events 
  h_nVtx->Fill(vertexCollection->size());

  Handle< HepMCProduct > hepProd ;
  ev.getByLabel("generator",hepProd) ;
  const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();

  vector<Photon> higgsPhotons,higgsPhotons_trueVtx;
  int isEB=0; 
  int hasOneClus(0);   int hasTwoClus(0);   int hasThreeClus(0);   int hasMoreClus(0); 
  bool trueVtxFound=false;
  math::XYZPoint trueVtx(0.,0.,0.);

  // loop over MC-truth particles
  for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
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
	  float energySC = scIt->rawEnergy();
	  float energySCtoCheck=0;

	  float phiSize = getPhiSize(scIt,true);
	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1

	    allCase_.FillSc(scIt,p,ebRecHits,eeRecHits);
	    // use this variable to classify EB SC's
	    int numOfBC = scIt->clustersSize();
	    std::cout << "EB numOfBC is: " << numOfBC << std::endl;
	    if      (numOfBC==0){ std::cout << "ZERO basic cluster found in an EB SC; mess? - Bailing out"; assert (-1);}
	    else if (numOfBC==1){ oneCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else if (numOfBC==2){ twoCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else if (numOfBC==3){ threeCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else                { moreCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}

	    // COPY1 //
	    h_scet_barl->Fill((scIt->energy()/cosh(scIt->eta()))/et_true);   // using: h_scet_barl as counter of EB matched photons
	    h_EoverEtrue_barl->Fill(scIt->energy()/(*p)->momentum().e());    //        supercluster that matches a MC-truth particle 
	    h_phiWidth->Fill(phiWidth);
	    h_phiWidthVsE->Fill(phiWidth,(*p)->momentum().e());	  
	    h_phiSize->Fill(phiSize);
	    h_phiSizeVsE->Fill(phiSize,(*p)->momentum().e());
	    h_phiSizeVsEt->Fill(phiSize, et_true);

	    int   whereIsMaxInDomino[35];
	    float whatIsMaxInDomino[35];
	    for(int u=0; u<35; u++) {
	      whereIsMaxInDomino[u]=-999;
	      whatIsMaxInDomino[u] =-999;
	    } 

	    std::vector< std::pair<DetId, float> >    theHitsAndFractions =  scIt->hitsAndFractions();  
	    // loop over all the components of this supercluster	    
	    for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
		idsIt != theHitsAndFractions.end(); ++idsIt) 
	      {
		float thePhiDistance = getPhiDistance( (*scIt->seed()).seed() , 
						       (*idsIt).first );
		float currEbEnergy    = (ebRecHits->find( (*idsIt).first ))->energy();
		energySCtoCheck       += currEbEnergy;
		h_phiShape_barl       -> Fill(thePhiDistance, currEbEnergy/energySC);
		h_absPhiShape_barl    -> Fill(fabs(thePhiDistance), currEbEnergy/energySC);
		h_phiShapeVsE_barl    -> Fill(thePhiDistance, energySC , currEbEnergy/energySC);   
		h_absPhiShapeVsE_barl -> Fill(fabs(thePhiDistance), energySC , currEbEnergy/energySC);   

		float theEtaDistance = getEtaDistance( (*scIt->seed()).seed() , 
						       (*idsIt).first );
		float theSeedEta = EBDetId ( (*scIt->seed()).seed().rawId() ).ieta();
		// treat EB+ and EB- separately, give the sign of ieta; and make a symmetrized plot
		h_etaShape_barl        -> Fill(theEtaDistance, currEbEnergy/energySC);
		h_etaShape_barlSymm    -> Fill( signum(theSeedEta) *theEtaDistance, currEbEnergy/energySC);
		h_etaPhiShape_barl     -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energySC);
		h_etaPhiShape_barlSymm -> Fill( signum(theSeedEta) *theEtaDistance, thePhiDistance, currEbEnergy/energySC);
		if ( theSeedEta>0 )  {
		  h_etaShape_barlPLus    -> Fill(theEtaDistance, currEbEnergy/energySC);
		  h_etaPhiShape_barlPLus -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energySC);
		}
		else{
		  h_etaShape_barlMinus    -> Fill(theEtaDistance, currEbEnergy/energySC);
		  h_etaPhiShape_barlMinus -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energySC);
		}

		// at each given phi [0... 34], heck where in eta max energy crystal occurs
		// phi=17 is the seed by construction, in order to allow +-17
		int theIntegerPhiDistance = ((int)thePhiDistance) +17;
		int theIntegerEtaDistance = ((int)theEtaDistance);

		if (currEbEnergy > whatIsMaxInDomino[theIntegerPhiDistance]){
		  whatIsMaxInDomino[theIntegerPhiDistance]  = currEbEnergy;
		  whereIsMaxInDomino[theIntegerPhiDistance] = theIntegerEtaDistance;

		}
	      } // loop over crystals of the SC

	    for(int u=0; u<35; u++) {
	      if (whereIsMaxInDomino[u]<-100) continue;
	      h_maxCryInDomino_barl        -> Fill(whereIsMaxInDomino[u]);
	      h_maxCryInDominoVsPhi_barl   -> Fill(u-17, whereIsMaxInDomino[u]);
	    }

	    // looping on basic clusters within the Supercluster
	    for(reco::CaloCluster_iterator bcIt = scIt->clustersBegin(); bcIt!=scIt->clustersEnd(); bcIt++)
	      {
		// fill the plots that follow only if the basic cluster is NOT the seed basic cluster
		if (   EBDetId(   (*bcIt)->seed().rawId())  
		       == EBDetId ( (*scIt->seed()).seed().rawId() )  
		       ) continue;

		EBDetId theBCSeed           = EBDetId(   (*bcIt)->seed().rawId());  // this is always SCseed? FIXme
		EBDetId theSCseed           = EBDetId(   ( *(scIt->seed()) ).seed().rawId()  ) ;

		int     etaOfMaxInsidBC     =-999; 
		int     phiOfMaxInsidBC     =-999; 
		float   energyOfMaxInsideBC =-888;

		std::vector< std::pair<DetId, float> >  theBCHitsAndFractions =  (*bcIt)->hitsAndFractions();
		for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theBCHitsAndFractions.begin(); 
		    idsIt != theBCHitsAndFractions.end(); ++idsIt) {
		  
		  float currCryEnergy    = (ebRecHits->find( (*idsIt).first ))->energy();
		  if (currCryEnergy > energyOfMaxInsideBC){
		    etaOfMaxInsidBC         = (   (int)getEtaDistance ( (*scIt->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
		    phiOfMaxInsidBC         = (   (int)getPhiDistance ( (*scIt->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
		    energyOfMaxInsideBC     = currCryEnergy;
		  } // if energy surpassed
		} // loop over theBCHitsAndFractions of BC
	      
		// fill the histograms with the location of the BC maximum, w.r.t. to SC seed
		if(fabs(theSCseed.ieta()) >2)
		  {
		    h_maxCryInLocMax_barlSymm           ->Fill( etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
		    h_maxCryInLocMaxVsPhi_barlSymm      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
		  }
		if(theSCseed.ieta() >2)
		  {
		    h_maxCryInLocMax_barlPLus           ->Fill( etaOfMaxInsidBC );
		    h_maxCryInLocMaxVsPhi_barlPLus      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
		  } else if (theSCseed.ieta() <-2)
		  {
		    h_maxCryInLocMax_barlMinus          ->Fill( etaOfMaxInsidBC );
		    h_maxCryInLocMaxVsPhi_barlMinus     ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
		  }


		h_phiBCminusPhiSeed_barl  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() );
		h_etaBCminusEtaSeed_barl  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		h_absEtaBCminusAbsEtaSeed_barl  -> Fill( fabs((*bcIt)->eta())   -  fabs((*scIt->seed()).eta())  );
		h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta()) );
		if ( fabs(theSCseed.ieta()) <= 25 ){
		  h_phiBCminusPhiSeed_barlm1  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm1  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta()) );
		  h_absEtaBCminusAbsEtaSeed_barlm1  -> Fill( fabs((*bcIt)->eta())  -  fabs((*scIt->seed()).eta())   );	} 
		else if	    ( 25<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 45 ){
		  h_phiBCminusPhiSeed_barlm2  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm2  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta()) );		
		  h_absEtaBCminusAbsEtaSeed_barlm2  -> Fill( fabs((*bcIt)->eta())  - fabs((*scIt->seed()).eta()) ); } 
		else if	    ( 45<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 65 ){
		  h_phiBCminusPhiSeed_barlm3  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm3  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta()) );		
		  h_absEtaBCminusAbsEtaSeed_barlm3  -> Fill( fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta())  );} 
		else if	    ( 65<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 85 ){
		  h_phiBCminusPhiSeed_barlm4  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm4  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*scIt->seed()).eta()) );
		  h_absEtaBCminusAbsEtaSeed_barlm4  -> Fill( fabs((*bcIt)->eta())  -  fabs((*scIt->seed()).eta()) );} 
		
	      }// loop over the BC's 
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

	    allCase_.FillSc(scIt,p,ebRecHits,eeRecHits);
	    // use this variable to classify EB SC's
	    int numOfBC = scIt->clustersSize();
	    std::cout << "EE numOfBC is: " << numOfBC << std::endl;
	    if      (numOfBC==0){ std::cout << "ZERO basic cluster found in an EE SC; mess? - Bailing out"; assert (-1);}
	    else if (numOfBC==1){ oneCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else if (numOfBC==2){ twoCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else if (numOfBC==3){ threeCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}
	    else                { moreCase_.FillSc(scIt,p,ebRecHits,eeRecHits);}

	    // COPY2 //
	    h_scet_endc->Fill((scIt->energy()/cosh(scIt->eta()))/et_true); // using: h_scet_barl as counter of EE matched photons
	    h_EoverEtrue_endc->Fill(scIt->energy()/(*p)->momentum().e());
	    h_phiWidth_endc->Fill(phiWidth);
	    h_phiWidthVsE_endc->Fill(phiWidth,(*p)->momentum().e());	  
	    h_phiSize_endc->Fill(phiSize);
	    h_phiSizeVsE_endc->Fill(phiSize,(*p)->momentum().e());	  
	    h_phiSizeVsEt_endc->Fill(phiSize,et_true);	  

	    std::vector< std::pair<DetId, float> >    theHitsAndFractions =  scIt->hitsAndFractions();  
	    for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
		idsIt != theHitsAndFractions.end(); ++idsIt) 
	      {
		float thePhiDistance = getPhiDistance( (*scIt->seed()).seed() , 
						       (*idsIt).first );
		float currEeEnergy   = (eeRecHits->find( (*idsIt).first ))->energy();
		energySCtoCheck      += currEeEnergy;
		h_phiShape_endc       -> Fill(thePhiDistance, currEeEnergy/energySC);
		h_absPhiShape_endc    -> Fill(fabs(thePhiDistance), currEeEnergy/energySC);
		h_phiShapeVsE_endc    -> Fill(thePhiDistance, energySC , currEeEnergy/energySC);   
		h_absPhiShapeVsE_endc -> Fill(fabs(thePhiDistance), energySC, currEeEnergy/energySC );   

	      } // loop over crystals of the SC
	    // COPY2 //

	  }// if SC matches one MC particle 
	}// if in EE and Et>20 GeV 
      }// loop over EE Sc's
    

      
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
	    
	    allCase_.FillGamma(phoIt,p);
	    // use this variable to classify the SC's associated to PHOTONS
	    int numOfBC = phoIt->superCluster()->clustersSize();
	    std::cout << "numOfBC for PHOTON is: " << numOfBC << std::endl;
	    if      (numOfBC==0){ std::cout << "ZERO basic cluster found in a PHOTOn; mess? - Bailing out"; assert (-1);}
	    else if (numOfBC==1){ oneCase_.FillGamma(phoIt,p);}
	    else if (numOfBC==2){ twoCase_.FillGamma(phoIt,p);}
	    else if (numOfBC==3){ threeCase_.FillGamma(phoIt,p);}
	    else                { moreCase_.FillGamma(phoIt,p);}

	    // COPY 3 //
	    if (phoIt->isEB()) {
	      h_PhoEoverEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e());
	      h_E5x5overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
	      h_E5x5overEtrueVsEphoEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e(),phoIt->energy()/(*p)->momentum().e());
	      if (phoIt->r9()>0.94) { // r9 cut is 0.94 for EB and 0.95 for EE  
		h_E5x5R9overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoER9overEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e()); // counter r9 EB
	      } else {
		h_E5x5notR9overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoEnotR9overEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e()); // counter !r9 EB
	      }
	      isEB++;
	    }
	    else if (phoIt->isEE()) {
	      h_PhoEoverEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      h_E5x5overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
	      h_E5x5overEtrueVsEphoEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e(),phoIt->energy()/(*p)->momentum().e());
	      if (phoIt->r9()>0.94) { // r9 cut is 0.94 for EB and 0.95 for EE  
		h_E5x5R9overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoER9overEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      } else {
		h_E5x5notR9overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoEnotR9overEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      }
	    }//clustersSize()
	    // COPY 3 //

	    if (mother2!=0 && mother2->pdg_id()==25) { // specifically look out for HIGGses
	      if (phoIt->isEB() || phoIt->isEE()) {
		// use this variable to classify SC's associated to photons
		int numOfBC = phoIt->superCluster()->clustersSize();
		if      (numOfBC==1) hasOneClus++;
		else if (numOfBC==2) hasTwoClus++;
		else if (numOfBC==3) hasThreeClus++;
		else                 hasMoreClus++;
		
		higgsPhotons.push_back(*phoIt);
		Photon localPho = Photon(*phoIt);
		if (!trueVtxFound) cout << "Error: true vertex not found!" << endl;
		localPho.setVertex(trueVtx);
		higgsPhotons_trueVtx.push_back(localPho);
		//cout << phoIt->et() << " " << localPho.et() << " " << trueVtx << " " << recoVtx << endl;
	      }
	    } // if mother
	  }// if matching between SC and photon
	}// loop over photons
	
      }// if it's MC-truth photon

    }
  }// loop over all MC-truth particles

  allCase_.FillH(higgsPhotons, higgsPhotons_trueVtx);
  if(hasOneClus==2)    oneCase_.FillH(higgsPhotons, higgsPhotons_trueVtx);
  if(hasTwoClus==2)    twoCase_.FillH(higgsPhotons, higgsPhotons_trueVtx);
  if(hasThreeClus==2)  threeCase_.FillH(higgsPhotons, higgsPhotons_trueVtx);
  if(hasMoreClus==2)   moreCase_.FillH(higgsPhotons, higgsPhotons_trueVtx);

  //COPY4
  if (higgsPhotons.size()>1) {
    if ((higgsPhotons[0].et()>40. && higgsPhotons[1].et()>30.) ||
	(higgsPhotons[0].et()>30. && higgsPhotons[1].et()>40.)) { // get two highest energy candidates
      double mHiggs = (higgsPhotons[0].p4()+higgsPhotons[1].p4()).M();
      if (isEB==2) {
	h_mHiggs_EBEB->Fill(mHiggs);
      } else if (isEB==1) {
	h_mHiggs_EBEE->Fill(mHiggs);
      } else {
	h_mHiggs_EEEE->Fill(mHiggs);
      }
    }
  }

  if (higgsPhotons_trueVtx.size()>1) {
    //cout << higgsPhotons_trueVtx[0].et() << " " << higgsPhotons_trueVtx[1].et() << endl;
    if ((higgsPhotons_trueVtx[0].et()>40. && higgsPhotons_trueVtx[1].et()>30.) ||
	(higgsPhotons_trueVtx[0].et()>30. && higgsPhotons_trueVtx[1].et()>40.)) {
      double mHiggs_trueVtx = (higgsPhotons_trueVtx[0].p4()+higgsPhotons_trueVtx[1].p4()).M();
      if (isEB==2) {
	h_mHiggs_EBEB_trueVtx->Fill(mHiggs_trueVtx);
      } else if (isEB==1) {
	h_mHiggs_EBEE_trueVtx->Fill(mHiggs_trueVtx);
      } else {
	h_mHiggs_EEEE_trueVtx->Fill(mHiggs_trueVtx);
      }
    }// if ET selections
  }// if there's a true vertex
  //COPY4

}// END of analyze


// ------------ method called once each job just before starting event loop  ------------
void 
SCwithPUAnalysis::beginJob()
{

  allCase_.setRawEnergy(useRawEnergy_);
  oneCase_.setRawEnergy(useRawEnergy_);   twoCase_.setRawEnergy(useRawEnergy_);
  threeCase_.setRawEnergy(useRawEnergy_); moreCase_.setRawEnergy(useRawEnergy_);
  
  edm::Service<TFileService> fs;
  std::cout << "[SCwithPUAnalysis] making directories" << std::endl;
  
  TFileDirectory subDir=fs->mkdir("allCase");  
  allCase_.Book(subDir);
  subDir=fs->mkdir("oneBC");     oneCase_.Book(subDir);
  subDir=fs->mkdir("twoBC");     twoCase_.Book(subDir);
  subDir=fs->mkdir("threeBC");   threeCase_.Book(subDir);
  subDir=fs->mkdir("moreBC");    moreCase_.Book(subDir);
  
  // keeping these ones for COMPARISON and validation
  h_scet_barl = fs->make<TH1F>("h_scet_barl","SC ET over true, barrel; EB: E_{T,superC}/E_{T,true}",90,0.75,1.2);  // using: h_scet_barl as counter of EB matched photons 
  h_scet_endc = fs->make<TH1F>("h_scet_endc","SC ET over true, endcap; EE: E_{T,superC}/E_{T,true}",90,0.75,1.2);  // using: h_scet_endc as counter of EE matched photons 
  h_EoverEtrue_barl = fs->make<TH1F>("h_EoverEtrue_barl","E_SC/Etrue, barrel; EB: E_{superC}/E_{true}",90,0.75,1.2);
  h_EoverEtrue_endc = fs->make<TH1F>("h_EoverEtrue_endc","E_SC/Etrue, endcap; EE: E_{superC}/E_{true}",90,0.75,1.2);
  h_E5x5overEtrue_barl = fs->make<TH1F>("h_E5x5overEtrue_barl","E5x5/Etrue, barrel; EB: E_{5x5}/E_{true}",90,0.75,1.2);
  h_PhoEoverEtrue_barl = fs->make<TH1F>("h_PhoEoverEtrue_barl","E_Photon/Etrue, barrel; EB: E_{photon}/E_{true}",90,0.75,1.2);
  h_E5x5R9overEtrue_barl = fs->make<TH1F>("h_E5x5R9overEtrue_barl","E5x5/Etrue for R9>0.94, barrel; EB: E_{5x5, r9>0.94}/E_{true}",90,0.75,1.2);
  h_PhoER9overEtrue_barl = fs->make<TH1F>("h_PhoER9overEtrue_barl","E_Photon/Etrue for R9>0.94, barrel; EB: E_{photon, r9>0.94}/E_{true}",90,0.75,1.2); // counter r9 EB
  h_E5x5notR9overEtrue_barl = fs->make<TH1F>("h_E5x5notR9overEtrue_barl","E5x5/Etrue for R9<0.94, barrel; EB: E_{5x5, r9<0.94}/E_{true}",90,0.75,1.2);
  h_PhoEnotR9overEtrue_barl = fs->make<TH1F>("h_PhoEnotR9overEtrue_barl","E_Photon/Etrue for R9<0.94, barrel; EB: E_{photon, r9<0.94}/E_{true}",90,0.75,1.2); // counter !r9 EB
  h_E5x5overEtrue_endc = fs->make<TH1F>("h_E5x5overEtrue_endc","E5x5/Etrue, endcap; EE: E_{5x5}/E_{true}",90,0.75,1.2);
  h_PhoEoverEtrue_endc = fs->make<TH1F>("h_PhoEoverEtrue_endc","E_Photon/Etrue, endcap; EE: E_{photon}/E_{true}",90,0.75,1.2);
  h_E5x5R9overEtrue_endc = fs->make<TH1F>("h_E5x5R9overEtrue_endc","E5x5/Etrue for R9>0.95, endcap; EE: E_{5x5, r9>0.95}/E_{true}",90,0.75,1.2);
  h_PhoER9overEtrue_endc = fs->make<TH1F>("h_PhoER9overEtrue_endc","E_Photon/Etrue for R9>0.95, endcap; EE: E_{photon, r9>0.95}/E_{true}",90,0.75,1.2); // counter r9 EE
  h_E5x5notR9overEtrue_endc = fs->make<TH1F>("h_E5x5notR9overEtrue_endc","E5x5/Etrue for R9<0.95, endcap; EE: E_{5x5, r9<0.95}/E_{true}",90,0.75,1.2);
  h_PhoEnotR9overEtrue_endc = fs->make<TH1F>("h_PhoEnotR9overEtrue_endc","E_Photon/Etrue for R9<0.95, endcap; EE: E_{photon, r9<0.95}/E_{true}",90,0.75,1.2); // counter !r9 EE
  h_nVtx = fs->make<TH1F>("h_nVtx","no. of primary vertices; num vertices reco",30,0.,30.);  // usig h_nVtx  histogram can be used as a counter of events 
  h_dzVtx = fs->make<TH1F>("h_dzVtx","delta_z for reconstructed PV w.r.t. true PV",200,-5.,5.);
  h_mHiggs_EBEB = fs->make<TH1F>("h_mHiggs_EBEB","2 photon invariant mass, EBEB; EB-EB: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEE = fs->make<TH1F>("h_mHiggs_EBEE","2 photon invariant mass, EBEE; EB-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EEEE = fs->make<TH1F>("h_mHiggs_EEEE","2 photon invariant mass, EEEE; EE-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEB_trueVtx = fs->make<TH1F>("h_mHiggs_EBEB_trueVtx","2 photon invariant mass, EBEB, true PV; EB-EB: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EBEE_trueVtx = fs->make<TH1F>("h_mHiggs_EBEE_trueVtx","2 photon invariant mass, EBEE, true PV; EB-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EEEE_trueVtx = fs->make<TH1F>("h_mHiggs_EEEE_trueVtx","2 photon invariant mass, EEEE, true PV; EE-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_E5x5overEtrueVsEphoEtrue_barl = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_barl","E5x5/Etrue vs Epho/Etrue, barrel;EB: E_{5x5}/E_{true}; EB: E_{photon}/E_{true}",90,0.75,1.2,90,0.75,1.2);
  h_E5x5overEtrueVsEphoEtrue_endc = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_endc","E5x5/Etrue vs Epho/Etrue, endcap;EE: E_{5x5}/E_{true}; EE: E_{photon}/E_{true}",90,0.75,1.2,90,0.75,1.2);
  h_phiWidth = fs->make<TH1F>("h_phiWidth_barl","phi Width (barrel); EB: i#phi width", 100,0,0.2);
  h_phiWidthVsE = fs->make<TH2F>("h_phiWidthVsE_barl","phi Width Vs. E (Barrel); EB: i#phi width; EB: E [GeV]", 100,0,0.2,100,0,200);
  h_phiWidth_endc = fs->make<TH1F>("h_phiWidth_endc","phi Width (Endcap); EE: i#phi width", 100,0,0.2);
  h_phiWidthVsE_endc = fs->make<TH2F>("h_phiWidthVsE_endc","phi Width Vs. E (Endcap); EB: i#phi width; EB: E [GeV]", 100,0,0.2,100,0,200);
  h_phiSize         = fs->make<TH1F>("h_phiSize_barl","phi Size (barrel); EB: i#phi size", 50,0,50);
  h_phiSizeVsE      = fs->make<TH2F>("h_phiSizeVsE_barl","phi Size Vs. E (Barrel); EB: i#phi size; EB: E [GeV]", 50,0,50.,100,0,200);
  h_phiSizeVsEt     = fs->make<TH2F>("h_phiSizeVsEt_barl","phi Size Vs. E_{T} (Barrel); EB: i#phi size; EB: E_{T} [GeV]", 50,0,50.,100,0,200);
  h_phiSize_endc    = fs->make<TH1F>("h_phiSize_endc","phi Size (Endcap); EE: i#phi size", 50,0,1);
  h_phiSizeVsE_endc = fs->make<TH2F>("h_phiSizeVsE_endc","phi Size Vs. E (Endcap); EB: i#phi size; EE: E [GeV]", 50,0,1.,100,0,200);
  h_phiSizeVsEt_endc= fs->make<TH2F>("h_phiSizeVsEt_endc","phi Size Vs. E_{T} (Endcap); EB: i#phi size; EE: E_{T} [GeV]", 50,0,1.,100,0,200);
  h_phiShape_barl       = fs->make<TH1F>("h_phiShape_barl","phi Shape (barrel); EB i#phi - i#phi_{SCseed}", 35,-17,18); 
  h_absPhiShape_barl    = fs->make<TH1F>("h_absPhiShape_barl","phi AbsShape (barrel); EB  abs(i#phi - i#phi_{SCseed})", 18,0,18); 
  h_phiShape_endc       = fs->make<TH1F>("h_phiShape_endc","phi Shape (endcap) EE i#phi - i#phi_{SCseed}", 35,-17,18); 
  h_absPhiShape_endc    = fs->make<TH1F>("h_absPhiShape_endc","phi AbsShape (endcap); EE  abs(i#phi - i#phi_{SCseed})", 18,0,18); 
  h_phiShapeVsE_barl    = fs->make<TH2F>("h_phiShapeVsE_barl","phi Shape Vs E (barrel); EB i#phi - i#phi_{SCseed}; EB E [GeV]", 35,-17,18,100,0,200); 
  h_absPhiShapeVsE_barl = fs->make<TH2F>("h_absPhiShapeVsE_barl","phi AbsShape Vs E (barrel); EB  abs(i#phi - i#phi_{SCseed}); EB E [GeV]", 18,0,18,100,0,200); 
  h_phiShapeVsE_endc    = fs->make<TH2F>("h_phiShapeVsE_endc","phi Shape Vs E (endcap); EE i#phi - i#phi_{SCseed}; EE E [GeV]", 35,-17,18,100,0,200); 
  h_absPhiShapeVsE_endc = fs->make<TH2F>("h_absPhiShapeVsE_endc","phi AbsShape Vs E (endcap); EE  abs(i#phi - i#phi_{SCseed}); EE E [GeV] ", 18,0,18,100,0,200); 

  h_etaShape_barl       = fs->make<TH1F>("h_etaShape_barl","eta Shape (barrel); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlPLus   = fs->make<TH1F>("h_etaShape_barlPLus","eta Shape (barrel plus); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlMinus  = fs->make<TH1F>("h_etaShape_barlMinus","eta Shape (barrel minus); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlSymm   = fs->make<TH1F>("h_etaShape_barlSymm","eta Shape (barrel symm); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaPhiShape_barl    = fs->make<TH2F>("h_etaPhiShape_barl","eta Shape (barrel); EB i#phi_{BC} - phi_{SCseed}; EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlPLus  = fs->make<TH2F>("h_etaPhiShape_barlPlus","eta Shape (barrel plus); EB+ i#phi_{BC} - phi_{SCseed}; EB+ i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlMinus = fs->make<TH2F>("h_etaPhiShape_barlMinus","eta Shape (barrel minus); EB- i#phi_{BC} - phi_{SCseed}; EB- i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlSymm  = fs->make<TH2F>("h_etaPhiShape_barlSymm","eta Shape (barrel symm); EBsymm i#phi_{BC} - phi_{SCseed}; EBsymm i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 

  // each phi bin is half crystal wide; each eta bin is 1/10 crystal wide
  h_phiBCminusPhiSeed_barl= fs->make<TH1F>("h_phiBCminusPhiSeed_barl","h_phiBCminusPhiSeed_barl; EB #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barl= fs->make<TH1F>("h_etaBCminusEtaSeed_barl","h_etaBCminusEtaSeed_barl; EB #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barl= fs->make<TH1F>("h_absEtaBCminusAbsEtaSeed_barl","h_absEtaBCminusAbsEtaSeed_barl; EB |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EB #phi_{BC} - #phi_{SCseed}; EB |#eta_{BC}| - |#eta_{SCseed}| ", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm1= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm1","h_phiBCminusPhiSeed_barl; EBm1  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm1= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm1","h_etaBCminusEtaSeed_barl; EBm1  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm1= fs->make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm1","h_absEtaBCminusAbsEtaSeed_barl; EBm1  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EBm1 #phi_{BC} - #phi_{SCseed}; EBm1 |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm2= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm2","h_phiBCminusPhiSeed_barlm2; EBm2 #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm2= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm2","h_etaBCminusEtaSeed_barlm2; EBm2 #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm2= fs->make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm2","h_absEtaBCminusAbsEtaSeed_barlm2; EBm2 |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2;  EBm2  #phi_{BC} - #phi_{SCseed}; EBm2  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm3= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm3","h_phiBCminusPhiSeed_barlm3; EBm3  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm3= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm3","h_etaBCminusEtaSeed_barlm3; EBm3  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm3= fs->make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm3","h_absEtaBCminusAbsEtaSeed_barlm3; EBm3  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3;  EBm3  #phi_{BC} - #phi_{SCseed}; EBm3  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm4= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm4","h_phiBCminusPhiSeed_barlm4; EBm4  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm4= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm4","h_etaBCminusEtaSeed_barlm4; EBm4  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm4= fs->make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm4","h_absEtaBCminusAbsEtaSeed_barlm4; EBm4  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4;  EBm4  #phi_{BC} - #phi_{SCseed}; EBm4  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));

  h_maxCryInDomino_barl           = fs->make<TH1F>("h_maxCryInDomino_barl","max Cry In Domino (barrel); i#eta", 5,-2,3); 
  h_maxCryInDominoVsPhi_barl      = fs->make<TH2F>("h_maxCryInDominoVsPhi_barl","max Cry In Domino Vs i#phi (barrel); EB i#eta_{BC} - i#eta_{SCseed}; EB i#phi_{BC} - i#phi_{SCseed}", 35,-17,18,5,-2,3); 

  h_maxCryInLocMax_barlSymm       = fs->make<TH1F>("h_maxCryInLocalMax_barlSymm","max Cry In Local Max (EB); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlSymm  = fs->make<TH2F>("h_maxCryInLocalMaxVsPhi_barlSymm","max Cry In Local Max Vs i#phi (EB); EB i#phi_{BC} - i#phi_{SCseed}; EB  i#eta_{BC} - i#eta_{SCseed};", 35,-17,18,5,-2,3); 
  h_maxCryInLocMax_barlPLus       = fs->make<TH1F>("h_maxCryInLocalMax_barlPLus","max Cry In Local Max (EB+); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlPLus  = fs->make<TH2F>("h_maxCryInLocalMaxVsPhi_barlPlus","max Cry In Local Max Vs i#phi (EB+); EB i#phi_{BC} - i#phi_{SCseed}; EB+  i#eta_{BC} - i#eta_{SCseed}", 35,-17,18,5,-2,3); 
  h_maxCryInLocMax_barlMinus      = fs->make<TH1F>("h_maxCryInLocalMax_barlMinus","max Cry In Local Max (EB-); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlMinus = fs->make<TH2F>("h_maxCryInLocalMaxVsPhi_barlMinus","max Cry In Local Max Vs i#phi (EB-); EB i#phi_{BC} - i#phi_{SCseed}; EB-  i#eta - i#eta_{SCseed}", 35,-17,18,5,-2,3); 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCwithPUAnalysis::endJob() {
}

// compute eta at detector, given physics eta and location of primary vertex
float SCwithPUAnalysis::etaTransformation(  float EtaParticle , float Zvertex)  {

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


//define this as a plug-in
//DEFINE_FWK_MODULE(SCwithPUAnalysis);

#include "FWCore/Framework/interface/MakerMacros.h"  
DEFINE_FWK_MODULE( SCwithPUAnalysis );
