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
// $Id: SCwithPUAnalysis.cc,v 1.13 2011/04/08 14:21:10 franzoni Exp $
//
//

// analyzer from David Futyan - Imperial College
// imported for usage of UMN G. Franzoni, Y. Kubota, V. Rekovic

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
//

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

#define PI    3.14159
#define TWOPI 6.28318539

#define R_ECAL           136.5  // radius of maximum containement
#define Z_Endcap         328.0
#define etaBarrelEndcap  1.479
#define eeCrySize        2.5


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

  TH1F *h_scet_barl;
  TH1F *h_scet_endc;
  TH1F *h_EoverEtrue_barl;
  TH1F *h_EoverEtrue_endc;
  TH1F *h_E5x5overEtrue_barl;
  TH1F *h_E5x5overEtrue_endc;
  TH1F *h_E5x5R9overEtrue_barl;
  TH1F *h_PhoER9overEtrue_barl;
  TH1F *h_E5x5notR9overEtrue_barl;
  TH1F *h_PhoEnotR9overEtrue_barl;
  TH1F *h_PhoEoverEtrue_barl;
  TH1F *h_PhoEoverEtrue_endc;
  TH1F *h_E5x5R9overEtrue_endc;
  TH1F *h_PhoER9overEtrue_endc;
  TH1F *h_E5x5notR9overEtrue_endc;
  TH1F *h_PhoEnotR9overEtrue_endc;
  TH1F *h_nVtx;
  TH1F *h_dzVtx;
  TH1F *h_mHiggs_EBEB;
  TH1F *h_mHiggs_EBEE;
  TH1F *h_mHiggs_EEEE;
  TH1F *h_mHiggs_EBEB_trueVtx;
  TH1F *h_mHiggs_EBEE_trueVtx;
  TH1F *h_mHiggs_EEEE_trueVtx;
  TH2F *h_E5x5overEtrueVsEphoEtrue_barl;
  TH2F *h_E5x5overEtrueVsEphoEtrue_endc;
  TH1F *h_phiWidth;
  TH2F *h_phiWidthVsE;
  TH1F *h_phiWidth_endc;
  TH2F *h_phiWidthVsE_endc;
  TH1F *h_phiSize;
  TH2F *h_phiSizeVsE;
  TH2F *h_phiSizeVsEt;
  TH1F *h_phiSize_endc;
  TH2F *h_phiSizeVsE_endc;
  TH2F *h_phiSizeVsEt_endc;
  TH1F *h_phiShape_barl; 
  TH1F *h_absPhiShape_barl; 
  TH1F *h_phiShape_endc; 
  TH1F *h_absPhiShape_endc; 
  TH2F *h_phiShapeVsE_barl; 
  TH2F *h_absPhiShapeVsE_barl; 
  TH2F *h_phiShapeVsE_endc; 
  TH2F *h_absPhiShapeVsE_endc; 
  TH1F *h_etaShape_barl;
  TH1F *h_etaShape_barlPLus;
  TH1F *h_etaShape_barlMinus;
  TH1F *h_etaShape_barlSymm;
  TH2F *h_etaPhiShape_barl;
  TH2F *h_etaPhiShape_barlPLus;
  TH2F *h_etaPhiShape_barlMinus;
  TH2F *h_etaPhiShape_barlSymm;
  TH1F *h_phiBCminusPhiSeed_barl;
  TH1F *h_etaBCminusEtaSeed_barl;
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1;
  TH1F *h_phiBCminusPhiSeed_barlm1;
  TH1F *h_etaBCminusEtaSeed_barlm1;
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2;
  TH1F *h_phiBCminusPhiSeed_barlm2;
  TH1F *h_etaBCminusEtaSeed_barlm2;
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3;
  TH1F *h_phiBCminusPhiSeed_barlm3;
  TH1F *h_etaBCminusEtaSeed_barlm3;
  TH2F *h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4;
  TH1F *h_phiBCminusPhiSeed_barlm4;
  TH1F *h_etaBCminusEtaSeed_barlm4;
  TH1F *h_maxCryInDomino_barl;
  TH2F *h_maxCryInDominoVsPhi_barl;
  TH1F *h_maxCryInLocMax_barlPLus;
  TH2F *h_maxCryInLocMaxVsPhi_barlPLus;
  TH1F *h_maxCryInLocMax_barlMinus;
  TH2F *h_maxCryInLocMaxVsPhi_barlMinus;
  TH1F *h_maxCryInLocMax_barlSymm;
  TH2F *h_maxCryInLocMaxVsPhi_barlSymm;
};

//
// constants, enums and typedefs
//

float signum(float x) {    return (x>0)?1:((x<0)?-1:0); }

//
// static data member definitions
//

//
// constructors and destructor
//
SCwithPUAnalysis::SCwithPUAnalysis(const edm::ParameterSet& iConfig)

{

}


SCwithPUAnalysis::~SCwithPUAnalysis()
{
}



// returns distance between two crystals, in units of crystal size (both  for EB and EE)
float getEtaDistance(DetId theSeed, DetId theCry ){

  if ( theSeed.subdetId()==EcalBarrel && theCry.subdetId()==EcalBarrel)     // the barrel case
    {
      float theSeedEta = EBDetId ( theSeed.rawId() ).ieta();
      float theCryEta  = EBDetId ( theCry.rawId() ).ieta();
      
      if ( theSeedEta * theCryEta > 0 )     // crystal and seed are on the same side of EB
	{ 
	  return (theCryEta - theSeedEta);
	}
      else	                           // crystal and seed are on opposite sides of EB
	{
	  if (theSeedEta<0) return (theCryEta - theSeedEta -1);
	  else              return (theCryEta - theSeedEta +1);
	}
    }
  else if ( theSeed.subdetId()==EcalEndcap && theCry.subdetId()==EcalEndcap) // the endcap case
    {
      // to be implemented w/ geometry...
      return -9999;
    }

  return -9999;

}


// returns distance between two crystals, in units of crystal size (both  for EB and EE)
float getPhiDistance(DetId theSeed, DetId theCry ){

  if ( theSeed.subdetId()==EcalBarrel && theCry.subdetId()==EcalBarrel)     // the barrel case
    {
      float theSeedPhi = EBDetId ( theSeed.rawId() ).iphi();
      float theCryPhi  = EBDetId ( theCry.rawId() ).iphi();
      if (fabs( theCryPhi - theSeedPhi ) > 180 )                            // this means cluster is crossing the iphi==1 <--> iphi==360 boundary
	{
	  theCryPhi -= 360 * signum( theCryPhi - theSeedPhi );
	}
      return (theCryPhi-theSeedPhi);    // this is in units of EB crystals
    }
  else if ( theSeed.subdetId()==EcalEndcap && theCry.subdetId()==EcalEndcap) // the endcap case
    {
      float theSeedX = EEDetId ( theSeed.rawId() ).ix() -50;
      float theSeedY = EEDetId ( theSeed.rawId() ).iy() -50;
      float theSeedPhi = atan2( theSeedX , theSeedY); 
      float theCryX  = EEDetId ( theCry.rawId() ).ix()  -50;
      float theCryY  = EEDetId ( theCry.rawId() ).iy()  -50;
      float theCryPhi  = atan2( theCryX , theCryY); 
      
      if ( fabs( theCryPhi - theSeedPhi ) > PI ) // this means cluster is crossing the phi==-pi <--> phi==pi boundary
	{
	  theCryPhi -= TWOPI * signum ( theCryPhi - theSeedPhi );
	}
      //        angular span              radius in units of EE crystal size 
      return ( (theCryPhi - theSeedPhi) * pow( theSeedX*theSeedX + theSeedY*theSeedY, 0.5) ) ;
    }
  else // unsupported mixed case
    {
      return -9999;
    }
  
}








// PhiSize is in unit of: crystals for EB 
//                        radiants for EE (change this to crystals too)             
float getPhiSize(reco::SuperClusterCollection::const_iterator scIt , float isBarrel ){
  // scIt is an iterator pointing to a supercluster
  // isBarrel has to be passed from outside - resolve there the clusters across EB and EE
  
  std::vector< std::pair<DetId, float> >  theHitsAndFractions =  (*scIt).hitsAndFractions();  
  
  float phiLow =0;
  float phiHigh=0;
  float phiFirstEB =-9999;  // a sort of seeding is done with the first crystal
  float phiFirstEE =-9999;  // in order to handle the boundaries where phi is discontinuos
  
  for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
      idsIt != theHitsAndFractions.end(); ++idsIt) 
    {
      
      float thePhi=-99999;
      
      if      ( isBarrel && (*idsIt).first.subdetId()==EcalBarrel ) { 

	thePhi = EBDetId(  (*idsIt).first.rawId()  ).iphi() ; 
	if (phiFirstEB<-9998) phiFirstEB = thePhi;
	
	if (fabs( thePhi - phiFirstEB ) > 180 ) // this means cluster is crossing the iphi==1 <--> iphi==360 boundary
	  {
	    //	    std::cout << "====== EB thePhi was: " << thePhi << " and now has become: " << ( thePhi + 360 * signum( thePhi - phiFirstEB ) * -1 )<< std::endl;
	    thePhi -= 360 * signum( thePhi - phiFirstEB );
	  }
      }

      else if ( (!isBarrel) &&  (*idsIt).first.subdetId()==EcalEndcap ) {
	double ix = EEDetId(  (*idsIt).first.rawId()  ).ix() -50;
	double iy = EEDetId(  (*idsIt).first.rawId()  ).iy() -50;
	thePhi = atan2( ix,  iy );

	if (phiFirstEE<-9998) phiFirstEE = thePhi;

	if (fabs( thePhi - phiFirstEE ) > PI ) // this means cluster is crossing the phi==-pi <--> phi==pi boundary
	  {
	    //	    std::cout << "====== EE thePhi was: " << thePhi << " and now has become: " << ( thePhi + TWOPI * signum( thePhi - phiFirstEE ) * -1 )<< std::endl;
	    thePhi -= TWOPI * signum( thePhi - phiFirstEE );
	  }
	
      }
      
      if (phiLow==0 && phiHigh==0)
	{ phiLow  = thePhi;
	  phiHigh = thePhi;	 }
      if       (thePhi < phiLow)  phiLow  = thePhi;
      else if  (thePhi > phiHigh) phiHigh = thePhi;

    }

  float thePhiSize = phiHigh-phiLow;

  // PhiSize is in unit of crystals for EB 
  //                       radiants for EE (change this to crystals too)             
  return thePhiSize;
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


  h_nVtx->Fill(vertexCollection->size());

  Handle< HepMCProduct > hepProd ;
  ev.getByLabel("generator",hepProd) ;
  const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();

  vector<Photon> higgsPhotons,higgsPhotons_trueVtx;
  int isEB=0;
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

      // Barrel SuperClusters
      for(SuperClusterCollection::const_iterator scIt = barrelSCCollection->begin(); scIt != barrelSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())<1.4442 && scIt->energy()/cosh(scIt->eta())>20.) {

	  // std::cout << "\n\n new SC number: " << (scIt-barrelSCCollection->begin()) << std::endl; //GF

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
	    h_scet_barl->Fill((scIt->energy()/cosh(scIt->eta()))/et_true);
	    h_EoverEtrue_barl->Fill(scIt->energy()/(*p)->momentum().e());
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

	    // std::cout << "size of BC clollection for this SC:  " << (scIt->clustersEnd() - scIt->clustersBegin() ) << std::endl; 
	    // looping on basic clusters within the Supercluster
	    for(CaloCluster_iterator bcIt = scIt->clustersBegin(); bcIt!=scIt->clustersEnd(); bcIt++)
	      {
		// std::cout << "BC number: "<< (bcIt-scIt->clustersBegin()) << std::endl;
		// std::cout << "\n\n new BC number: " << ( bcIt-(scIt->clustersBegin()) ) << std::endl;
		EBDetId theBCSeed           = EBDetId(   (*bcIt)->seed().rawId());  // this is always SCseed? FIXme
		//EBDetId theBCSeed           = EBDetId(   bcIt->seed().rawId());  // this is always SCseed? FIXme
		EBDetId theSCseed           = EBDetId(   ( *(scIt->seed()) ).seed().rawId()  ) ;
		// std::cout << " BC seed is: " << theBCSeed.rawId() << " while the SC seed is: " << theSCseed.rawId() << std::endl; 

		int     etaOfMaxInsidBC     =-999; 
		int     phiOfMaxInsidBC     =-999; 
		float   energyOfMaxInsideBC =-888;
		// std::cout << "++ BEFORE etaOfMaxInsidBC: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;

		std::vector< std::pair<DetId, float> >  theBCHitsAndFractions =  (*bcIt)->hitsAndFractions();
		for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theBCHitsAndFractions.begin(); 
		    idsIt != theBCHitsAndFractions.end(); ++idsIt) {
		  
		  float currCryEnergy    = (ebRecHits->find( (*idsIt).first ))->energy();
		  if (currCryEnergy > energyOfMaxInsideBC){
		    etaOfMaxInsidBC         = (   (int)getEtaDistance ( (*scIt->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
		    phiOfMaxInsidBC         = (   (int)getPhiDistance ( (*scIt->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
		    energyOfMaxInsideBC     = currCryEnergy;
		    //std::cout << "eta distance from SC seed is: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;
		  } // if energy surpassed
		} // loop over theBCHitsAndFractions of BC
		// std::cout << "++ AFTER etaOfMaxInsidBC: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;

		// fill the histograms with the location of the BC maximum, w.r.t. to SC seed
		if(fabs(theSCseed.ieta()) >2)
		  {
		    h_maxCryInLocMax_barlSymm           ->Fill( etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
		    h_maxCryInLocMaxVsPhi_barlSymm      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
		    //std::cout << "++ SYMM eta distance from SC seed is: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;
		  }
		if(theSCseed.ieta() >2)
		  {
		    h_maxCryInLocMax_barlPLus           ->Fill( etaOfMaxInsidBC );
		    h_maxCryInLocMaxVsPhi_barlPLus      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
		    //std::cout << "++ EB+ eta distance from SC seed is: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;
		  } else if (theSCseed.ieta() <-2)
		  {
		    h_maxCryInLocMax_barlMinus          ->Fill( etaOfMaxInsidBC );
		    h_maxCryInLocMaxVsPhi_barlMinus     ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
		    //std::cout << "++ EB- eta distance from SC seed is: " << etaOfMaxInsidBC << " phiOfMaxInsidBC  " << phiOfMaxInsidBC << " energy: " << energyOfMaxInsideBC << std::endl;
		  }

		// basic cluster position
		//std::cout << "position of basic cluster number: " << (bcIt-scIt->clustersBegin()) << " eta: " << (*bcIt)->eta() << " phi: " << (*bcIt)->phi() 
		//	  << "\t seed eta: " << (*scIt->seed()).eta() << " seed phi: " << (*scIt->seed()).phi() << std::endl;

		h_phiBCminusPhiSeed_barl  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() );
		h_etaBCminusEtaSeed_barl  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , (*bcIt)->eta() - (*scIt->seed()).eta() );
		if ( fabs(theSCseed.ieta()) <= 25 ){
		  h_phiBCminusPhiSeed_barlm1  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm1  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , (*bcIt)->eta() - (*scIt->seed()).eta() );
		} 
		else if	    ( 25<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 45 ){
		  h_phiBCminusPhiSeed_barlm2  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm2  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , (*bcIt)->eta() - (*scIt->seed()).eta() );		} 
		else if	    ( 45<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 65 ){
		  h_phiBCminusPhiSeed_barlm3  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm3  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , (*bcIt)->eta() - (*scIt->seed()).eta() );		} 
		else if	    ( 65<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 85 ){
		  h_phiBCminusPhiSeed_barlm4  -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() ); h_etaBCminusEtaSeed_barlm4  -> Fill( (*bcIt)->eta() - (*scIt->seed()).eta() );
		  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4   -> Fill( (*bcIt)->phi() - (*scIt->seed()).phi() , (*bcIt)->eta() - (*scIt->seed()).eta() );		} 
		
	      }// loop over the BC's 
	  }// if SC matches photons
	}// if EB and Et>20
      }// loop over superclusters
      
      // Endcap SuperClusters
      for(SuperClusterCollection::const_iterator scIt = endcapSCCollection->begin(); scIt != endcapSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())>1.566 && fabs(scIt->eta())<2.5 && scIt->energy()/cosh(scIt->eta())>20.) {
	  float deltaPhi = scIt->phi()-phi_true;
	  float deltaEta = scIt->eta()-etaEcal_true;
	  float phiWidth = scIt->phiWidth();
	  float energySC = scIt->rawEnergy();//gf
	  float energySCtoCheck=0;

	  float phiSize = getPhiSize(scIt,false);

	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1
	    h_scet_endc->Fill((scIt->energy()/cosh(scIt->eta()))/et_true);
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
	    // std::cout << "EE energySC: " << energySC << " energySCtoCheck: " << energySCtoCheck << std::endl;
	    // unresolved: energySC > energySCtoCheck???
	  }// if matching
	}
      }
      
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
	    if (phoIt->isEB()) {
	      h_PhoEoverEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e());
	      h_E5x5overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
	      h_E5x5overEtrueVsEphoEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e(),phoIt->energy()/(*p)->momentum().e());
	      if (phoIt->r9()>0.94) { // r9 cut is 0.94 for EB and 0.95 for EE  
		h_E5x5R9overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoER9overEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e());
	      } else {
		h_E5x5notR9overEtrue_barl->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoEnotR9overEtrue_barl->Fill(phoIt->energy()/(*p)->momentum().e());
	      }
	      isEB++;
	    }
	    else if (phoIt->isEE()) {
	      h_PhoEoverEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      h_E5x5overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
	      h_E5x5overEtrueVsEphoEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e(),phoIt->energy()/(*p)->momentum().e());
	      if (phoIt->r9()>0.95) { // r9 cut is 0.94 for EB and 0.95 for EE  
		h_E5x5R9overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoER9overEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      } else {
		h_E5x5notR9overEtrue_endc->Fill(phoIt->e5x5()/(*p)->momentum().e());
		h_PhoEnotR9overEtrue_endc->Fill(phoIt->energy()/(*p)->momentum().e());
	      }
	    }
	    if (mother2!=0 && mother2->pdg_id()==25) { // specifically look out for HIGGses
	      if (phoIt->isEB() || phoIt->isEE()) {
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
  
  if (higgsPhotons.size()>1) {
    if ((higgsPhotons[0].et()>40. && higgsPhotons[1].et()>30.) ||
	(higgsPhotons[0].et()>30. && higgsPhotons[1].et()>40.)) {
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

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
SCwithPUAnalysis::beginJob()
{

  edm::Service<TFileService> fs;

  h_scet_barl = fs->make<TH1F>("h_scet_barl","SC ET over true, barrel; EB: E_{T,superC}/E_{T,true}",120,0.6,1.3);
  h_scet_endc = fs->make<TH1F>("h_scet_endc","SC ET over true, endcap; EE: E_{T,superC}/E_{T,true}",120,0.6,1.3);
  h_EoverEtrue_barl = fs->make<TH1F>("h_EoverEtrue_barl","E_SC/Etrue, barrel; EB: E_{superC}/E_{true}",120,0.6,1.3);
  h_EoverEtrue_endc = fs->make<TH1F>("h_EoverEtrue_endc","E_SC/Etrue, endcap; EB: E_{superC}/E_{true}",120,0.6,1.3);
  h_E5x5overEtrue_barl = fs->make<TH1F>("h_E5x5overEtrue_barl","E5x5/Etrue, barrel; EB: E_{5x5}/E_{true}",120,0.6,1.3);
  h_PhoEoverEtrue_barl = fs->make<TH1F>("h_PhoEoverEtrue_barl","E_Photon/Etrue, barrel; EB: E_{photon}/E_{true}",120,0.6,1.3);
  h_E5x5R9overEtrue_barl = fs->make<TH1F>("h_E5x5R9overEtrue_barl","E5x5/Etrue for R9>0.94, barrel; EB: E_{5x5, r9>0.94}/E_{true}",120,0.6,1.3);
  h_PhoER9overEtrue_barl = fs->make<TH1F>("h_PhoER9overEtrue_barl","E_Photon/Etrue for R9>0.94, barrel; EB: E_{photon, r9>0.94}/E_{true}",120,0.6,1.3);
  h_E5x5notR9overEtrue_barl = fs->make<TH1F>("h_E5x5notR9overEtrue_barl","E5x5/Etrue for R9<0.94, barrel; EB: E_{5x5, r9<0.94}/E_{true}",120,0.6,1.3);
  h_PhoEnotR9overEtrue_barl = fs->make<TH1F>("h_PhoEnotR9overEtrue_barl","E_Photon/Etrue for R9<0.94, barrel; EB: E_{photon, r9<0.94}/E_{true}",120,0.6,1.3);
  h_E5x5overEtrue_endc = fs->make<TH1F>("h_E5x5overEtrue_endc","E5x5/Etrue, endcap; EE: E_{5x5}/E_{true}",120,0.6,1.3);
  h_PhoEoverEtrue_endc = fs->make<TH1F>("h_PhoEoverEtrue_endc","E_Photon/Etrue, endcap; EE: E_{photon}/E_{true}",120,0.6,1.3);
  h_E5x5R9overEtrue_endc = fs->make<TH1F>("h_E5x5R9overEtrue_endc","E5x5/Etrue for R9>0.95, endcap; EE: E_{5x5, r9>0.95}/E_{true}",120,0.6,1.3);
  h_PhoER9overEtrue_endc = fs->make<TH1F>("h_PhoER9overEtrue_endc","E_Photon/Etrue for R9>0.95, endcap; EE: E_{photon, r9>0.95}/E_{true}",120,0.6,1.3);
  h_E5x5notR9overEtrue_endc = fs->make<TH1F>("h_E5x5notR9overEtrue_endc","E5x5/Etrue for R9<0.95, endcap; EE: E_{5x5, r9<0.95}/E_{true}",120,0.6,1.3);
  h_PhoEnotR9overEtrue_endc = fs->make<TH1F>("h_PhoEnotR9overEtrue_endc","E_Photon/Etrue for R9<0.95, endcap; EE: E_{photon, r9<0.95}/E_{true}",120,0.6,1.3);
  h_nVtx = fs->make<TH1F>("h_nVtx","no. of primary vertices; num vertices reco",30,0.,30.);
  h_dzVtx = fs->make<TH1F>("h_dzVtx","delta_z for reconstructed PV w.r.t. true PV",200,-5.,5.);
  h_mHiggs_EBEB = fs->make<TH1F>("h_mHiggs_EBEB","2 photon invariant mass, EBEB; EB-EB: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEE = fs->make<TH1F>("h_mHiggs_EBEE","2 photon invariant mass, EBEE; EB-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EEEE = fs->make<TH1F>("h_mHiggs_EEEE","2 photon invariant mass, EEEE; EE-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEB_trueVtx = fs->make<TH1F>("h_mHiggs_EBEB_trueVtx","2 photon invariant mass, EBEB, true PV; EB-EB: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EBEE_trueVtx = fs->make<TH1F>("h_mHiggs_EBEE_trueVtx","2 photon invariant mass, EBEE, true PV; EB-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EEEE_trueVtx = fs->make<TH1F>("h_mHiggs_EEEE_trueVtx","2 photon invariant mass, EEEE, true PV; EE-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_E5x5overEtrueVsEphoEtrue_barl = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_barl","E5x5/Etrue vs Epho/Etrue, barrel;EB: E_{5x5}/E_{true}; EB: E_{photon}/E_{true}",120,0.6,1.3,120,0.6,1.3);
  h_E5x5overEtrueVsEphoEtrue_endc = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_endc","E5x5/Etrue vs Epho/Etrue, endcap;EE: E_{5x5}/E_{true}; EE: E_{photon}/E_{true}",120,0.6,1.3,120,0.6,1.3);
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
  
  h_phiBCminusPhiSeed_barl= fs->make<TH1F>("h_phiBCminusPhiSeed_barl","h_phiBCminusPhiSeed_barl; EB #phi_{BC} - #phi_{SCseed}", 50,-0.4,0.4); 
  h_etaBCminusEtaSeed_barl= fs->make<TH1F>("h_etaBCminusEtaSeed_barl","h_etaBCminusEtaSeed_barl; EB #eta_{BC} - #eta_{SCseed}", 50,-0.05,0.05); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EB #phi_{BC} - #phi_{SCseed}; EB #eta_{BC} - #eta_{SCseed}", 20,-0.4,0.4, 20,-0.05,0.05);
  h_phiBCminusPhiSeed_barlm1= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm1","h_phiBCminusPhiSeed_barl; EBm1  #phi_{BC} - #phi_{SCseed}", 50,-0.4,0.4); 
  h_etaBCminusEtaSeed_barlm1= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm1","h_etaBCminusEtaSeed_barl; EBm1  #eta_{BC} - #eta_{SCseed}", 50,-0.05,0.05); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EBm1 #phi_{BC} - #phi_{SCseed}; EBm1 #eta_{BC} - #eta_{SCseed}", 20,-0.4,0.4, 20,-0.05,0.05);
  h_phiBCminusPhiSeed_barlm2= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm2","h_phiBCminusPhiSeed_barlm2; EBm2 #phi_{BC} - #phi_{SCseed}", 50,-0.4,0.4); 
  h_etaBCminusEtaSeed_barlm2= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm2","h_etaBCminusEtaSeed_barlm2; EBm2 #eta_{BC} - #eta_{SCseed}", 50,-0.05,0.05); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2;  EBm2  #phi_{BC} - #phi_{SCseed}; EBm2  #eta_{BC} - #eta_{SCseed}", 20,-0.4,0.4, 20,-0.05,0.05);
  h_phiBCminusPhiSeed_barlm3= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm3","h_phiBCminusPhiSeed_barlm3; EBm3  #phi_{BC} - #phi_{SCseed}", 50,-0.4,0.4); 
  h_etaBCminusEtaSeed_barlm3= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm3","h_etaBCminusEtaSeed_barlm3; EBm3  #eta_{BC} - #eta_{SCseed}", 50,-0.05,0.05); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3;  EBm3  #phi_{BC} - #phi_{SCseed}; EBm3  #eta_{BC} - #eta_{SCseed}", 20,-0.4,0.4, 20,-0.05,0.05);
  h_phiBCminusPhiSeed_barlm4= fs->make<TH1F>("h_phiBCminusPhiSeed_barlm4","h_phiBCminusPhiSeed_barlm4; EBm4  #phi_{BC} - #phi_{SCseed}", 50,-0.4,0.4); 
  h_etaBCminusEtaSeed_barlm4= fs->make<TH1F>("h_etaBCminusEtaSeed_barlm4","h_etaBCminusEtaSeed_barlm4; EBm4  #eta_{BC} - #eta_{SCseed}", 50,-0.05,0.05); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4 = fs->make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4;  EBm4  #phi_{BC} - #phi_{SCseed}; EBm4  #eta_{BC} - #eta_{SCseed}", 20,-0.4,0.4, 20,-0.05,0.05);

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
