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
// $Id: SCwithPUAnalysis.cc,v 1.1 2011/03/14 22:27:06 franzoni Exp $
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
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

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
	//cout << " " << trueVtx << " " << recoVtx << endl;
      }

      
      float phi_true=(*p)->momentum().phi();
      float eta_true=(*p)->momentum().eta();
      float etaEcal_true = etaTransformation(eta_true, (*p)->production_vertex()->position().z()/10. );
      float et_true = (*p)->momentum().e()/cosh((*p)->momentum().eta());

      // Barrel SuperClusters
      for(SuperClusterCollection::const_iterator scIt = barrelSCCollection->begin(); scIt != barrelSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())<1.4442 && scIt->energy()/cosh(scIt->eta())>20.) {
	  // perform SC-MCparticle matching
	  float deltaPhi = scIt->phi()-phi_true;
	  float deltaEta = scIt->eta()-etaEcal_true;
	  float phiWidth = scIt->phiWidth();          // covariance of cluster in phi
	                                              // we want to look at the absolute extension too... phi_MAX-phi_min 

	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1
	    h_scet_barl->Fill((scIt->energy()/cosh(scIt->eta()))/et_true);
	    h_EoverEtrue_barl->Fill(scIt->energy()/(*p)->momentum().e());
	    h_phiWidth->Fill(phiWidth);
	    h_phiWidthVsE->Fill(phiWidth,(*p)->momentum().e());	  
	  }
	}
      }
      
      // Endcap SuperClusters
      for(SuperClusterCollection::const_iterator scIt = endcapSCCollection->begin(); scIt != endcapSCCollection->end(); scIt++) {
	if (fabs(scIt->eta())>1.566 && fabs(scIt->eta())<2.5 && scIt->energy()/cosh(scIt->eta())>20.) {
	  float deltaPhi = scIt->phi()-phi_true;
	  float deltaEta = scIt->eta()-etaEcal_true;
	  float phiWidth = scIt->phiWidth();

	  if ( deltaPhi > pi ) deltaPhi -= twopi;
	  if ( deltaPhi < -pi) deltaPhi += twopi;
	  float delta = sqrt( deltaPhi*deltaPhi+deltaEta*deltaEta);
	  if ( delta<0.1) { // match if dr<0.1
	    h_scet_endc->Fill((scIt->energy()/cosh(scIt->eta()))/et_true);
	    h_EoverEtrue_endc->Fill(scIt->energy()/(*p)->momentum().e());
	    h_phiWidth_endc->Fill(phiWidth);
	    h_phiWidthVsE_endc->Fill(phiWidth,(*p)->momentum().e());	  
	  }
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

  h_scet_barl = fs->make<TH1F>("h_scet_barl","SC ET, barrel",120,0.6,1.2);
  h_scet_endc = fs->make<TH1F>("h_scet_endc","SC ET, endcap",120,0.6,1.2);
  h_EoverEtrue_barl = fs->make<TH1F>("h_EoverEtrue_barl","E_SC/Etrue, barrel",120,0.6,1.2);
  h_EoverEtrue_endc = fs->make<TH1F>("h_EoverEtrue_endc","E_SC/Etrue, endcap",120,0.6,1.2);
  h_E5x5overEtrue_barl = fs->make<TH1F>("h_E5x5overEtrue_barl","E5x5/Etrue, barrel",120,0.6,1.2);
  h_PhoEoverEtrue_barl = fs->make<TH1F>("h_PhoEoverEtrue_barl","E_Photon/Etrue, barrel",120,0.6,1.2);
  h_E5x5R9overEtrue_barl = fs->make<TH1F>("h_E5x5R9overEtrue_barl","E5x5/Etrue for R9>0.93, barrel",120,0.6,1.2);
  h_PhoER9overEtrue_barl = fs->make<TH1F>("h_PhoER9overEtrue_barl","E_Photon/Etrue for R9>0.93, barrel",120,0.6,1.2);
  h_E5x5notR9overEtrue_barl = fs->make<TH1F>("h_E5x5notR9overEtrue_barl","E5x5/Etrue for R9<0.93, barrel",120,0.6,1.2);
  h_PhoEnotR9overEtrue_barl = fs->make<TH1F>("h_PhoEnotR9overEtrue_barl","E_Photon/Etrue for R9<0.93, barrel",120,0.6,1.2);
  h_E5x5overEtrue_endc = fs->make<TH1F>("h_E5x5overEtrue_endc","E5x5/Etrue, endcap",120,0.6,1.2);
  h_PhoEoverEtrue_endc = fs->make<TH1F>("h_PhoEoverEtrue_endc","E_Photon/Etrue, endcap",120,0.6,1.2);
  h_E5x5R9overEtrue_endc = fs->make<TH1F>("h_E5x5R9overEtrue_endc","E5x5/Etrue for R9>0.93, endcap",120,0.6,1.2);
  h_PhoER9overEtrue_endc = fs->make<TH1F>("h_PhoER9overEtrue_endc","E_Photon/Etrue for R9>0.93, endcap",120,0.6,1.2);
  h_E5x5notR9overEtrue_endc = fs->make<TH1F>("h_E5x5notR9overEtrue_endc","E5x5/Etrue for R9<0.93, endcap",120,0.6,1.2);
  h_PhoEnotR9overEtrue_endc = fs->make<TH1F>("h_PhoEnotR9overEtrue_endc","E_Photon/Etrue for R9<0.93, endcap",120,0.6,1.2);
  h_nVtx = fs->make<TH1F>("h_nVtx","no. of primary vertices",30,0.,30.);
  h_dzVtx = fs->make<TH1F>("h_dzVtx","delta_z for reconstructed PV w.r.t. true PV",200,-5.,5.);
  h_mHiggs_EBEB = fs->make<TH1F>("h_mHiggs_EBEB","2 photon invariant mass, EBEB",120,100.,140.);
  h_mHiggs_EBEE = fs->make<TH1F>("h_mHiggs_EBEE","2 photon invariant mass, EBEE",120,100.,140.);
  h_mHiggs_EEEE = fs->make<TH1F>("h_mHiggs_EEEE","2 photon invariant mass, EEEE",120,100.,140.);
  h_mHiggs_EBEB_trueVtx = fs->make<TH1F>("h_mHiggs_EBEB_trueVtx","2 photon invariant mass, EBEB, true PV",120,100.,140.);
  h_mHiggs_EBEE_trueVtx = fs->make<TH1F>("h_mHiggs_EBEE_trueVtx","2 photon invariant mass, EBEE, true PV",120,100.,140.);
  h_mHiggs_EEEE_trueVtx = fs->make<TH1F>("h_mHiggs_EEEE_trueVtx","2 photon invariant mass, EEEE, true PV",120,100.,140.);
  h_E5x5overEtrueVsEphoEtrue_barl = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_barl","E5x5/Etrue vs Epho/Etrue, barrel",120,0.6,1.2,120,0.6,1.2);
  h_E5x5overEtrueVsEphoEtrue_endc = fs->make<TH2F>("h_E5x5overEtrueVsEphoEtrue_endc","E5x5/Etrue vs Epho/Etrue, endcap",120,0.6,1.2,120,0.6,1.2);
  h_phiWidth = fs->make<TH1F>("h_phiWidth_barrel","phi Width (barrel)", 100,0,0.2);
  h_phiWidthVsE = fs->make<TH2F>("h_phiWidthVsE_barrel","phi Width Vs. E (Barrel)", 100,0,0.2,200,0,200);
  h_phiWidth_endc = fs->make<TH1F>("h_phiWidth_endc","phi Width (Endcap)", 100,0,0.2);
  h_phiWidthVsE_endc = fs->make<TH2F>("h_phiWidthVsE_endc","phi Width Vs. E (Endcap)", 100,0,0.2,200,0,200);




}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCwithPUAnalysis::endJob() {
}

// compute eta at detector, given physics eta and location of primary vertex
float SCwithPUAnalysis::etaTransformation(  float EtaParticle , float Zvertex)  {

  //---Definitions
  const float PI    = 3.1415927;
  //UNUSED const float TWOPI = 2.0*PI;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;  // radius of maximum containement
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479;

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
DEFINE_FWK_MODULE(SCwithPUAnalysis);
