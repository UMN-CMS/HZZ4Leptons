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
// $Id: HiggsGenAnalysis.cc,v 1.5 2012/05/11 02:14:43 afinkel Exp $
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
#include "TH2.h"
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
      double E1, E2;
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
      TH1D* h_higgsY, *PotRecoHiggs, //*h_higgsY_mm, *h_higgsY_ee, *h_higgsY_eFe, *mu1pt, *mu2pt, *mupt, *el1pt, *el2pt, *elpt,
      		*Phi1, *Phi, *CosTheta1, *CosTheta2, *CosTheta0,
      		//*HiggsRestP2, *HiggsP2,
      		*FarEEelEnergy, *HFelEnergy,
      		*FarEEelPt, *HFelPt,
      		*FarEEeta, *HFeta,
      		*Z2e1Pt, *Z2e2Pt, *Z2elPt,
      		*H4Mu, *H4GSFe, *H2mu2GSF, *HGSF_FEE_2mu, *HGSF_HF_2mu, *HGSF_HF_2GSF, *HGSF_FEE_2GSF, *H2GSFe2mu,
      		*Ztypes, *AllForwardE, *AllForwardPt, *AllZ1e1Pt, *AllZ1e2Pt,*AllZ1ElecPt, *AllZ1e1Eta, *AllZ1e2Eta, *AllZ1ElecEta;
      		
      TH2D *e1PtVsEta, *e2PtVsEta, *AllPtVsEta;
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
    hists.PotRecoHiggs = fs->make<TH1D > ("PotRecoHiggs","Potentially reconstructable Higgs rapidity", 100,-5.,5.) ;
    //hists.h_higgsY_mm  = fs->make<TH1D > ("h_higgsY_mm","Higgs rapidity", 100,-5.,5.) ; 
    //hists.h_higgsY_ee  = fs->make<TH1D > ("h_higgsY_ee","Higgs rapidity", 100,-5.,5.) ; 
    //hists.h_higgsY_eFe = fs->make<TH1D > ("h_higgsY_eFe","Higgs rapidity", 100,-5.,5.) ; 
     
    //hists.mu1pt = fs->make<TH1D >("mu1pt","Mu1 pt",50,0.,100.);    
    //hists.mu2pt = fs->make<TH1D >("mu2pt","Mu2 pt",50,0.,100.);
    //hists.mupt = fs->make<TH1D >("mupt","Muon pt",50,0.,100.);
    //hists.el1pt = fs->make<TH1D >("el1pt","El1 pt",50,0.,100.);
    //hists.el2pt = fs->make<TH1D >("el2pt","El2 pt",50,0.,100.);
    //hists.elpt = fs->make<TH1D >("elpt","Electron pt",50,0.,100.);
    //hists.HiggsRestP2 = fs->make<TH1D>("HiggsRestP2","H rest frame p^2",100,0.,150);
    //hists.HiggsP2 = fs->make<TH1D>("HiggsP2","H lab frame p^2",100,0.,1000);
    
    //hists.Z2e1Pt = fs->make<TH1D>("Z2e1Pt","Off-Shell Z GSFel1 Pt", 50, 0, 50);
    //hists.Z2e2Pt = fs->make<TH1D>("Z2e2Pt","Off-Shell Z GSFel2 Pt", 50, 0, 50);
    //hists.Z2elPt = fs->make<TH1D>("Z2elPt","Off-Shell Z GSFel-all Pt", 50, 0, 50);

    
    hists.Phi = fs->make<TH1D>("Phi","#Phi",20, -3.2,3.2);
    hists.Phi1 = fs->make<TH1D>("Phi1","#Phi_{1}",20, -3.2,3.2);
    hists.CosTheta0 = fs->make<TH1D>("CostTheta0","Cos(#theta_{0})",50,-1.,1.);
    hists.CosTheta1 = fs->make<TH1D>("CostTheta1","Cos(#theta_{1})",20,-1.,1.);
    hists.CosTheta2 = fs->make<TH1D>("CostTheta2","Cos(#theta_{2})",20,-1.,1.);
    
    hists.FarEEelEnergy = fs->make<TH1D>("FarEEelEnergy","(On-shell Z) FarEE electron energy", 100, 0., 1000);
    hists.HFelEnergy = fs->make<TH1D>("HFelEnergy","(On-shell Z) HF electron energy", 100, 0., 1000);    
    hists.FarEEelPt = fs->make<TH1D>("FarEEelPt","(On-shell Z) FarEE electron Pt", 50, 0., 150);
    hists.HFelPt = fs->make<TH1D>("HFelPt","(On-shell Z) HF electron Pt", 50, 0., 150);
    hists.FarEEeta = fs->make<TH1D>("FarEEeta", "Far EE electron Eta",100, -5, 5 );
    hists.HFeta = fs->make<TH1D>("HFeta", "HF electron Eta",100, -5, 5 );
    
    hists.H4Mu = fs->make<TH1D>("4Mu","H to 4Mu Y",20, -5., 5.);
    hists.H4GSFe = fs->make<TH1D>("4GSF","H to 4GSFel Y",20, -5., 5.);
    hists.H2mu2GSF = fs->make<TH1D>("2Mu2GSF","H to 2Mu & 2GSFel Y",20, -5., 5.);
    hists.H2GSFe2mu = fs->make<TH1D>("2GSFe2Mu","H to 2GSFel & 2Mu Y",20, -5., 5.);
    hists.HGSF_FEE_2mu = fs->make<TH1D>("GSF_FEE_2mu","H to GSF+FarEE & 2Mu Y",20, -5., 5.);
    hists.HGSF_HF_2mu = fs->make<TH1D>("GSF_HF_2mu","H to GSF+HF & 2Mu Y",20, -5., 5.);
    hists.HGSF_FEE_2GSF = fs->make<TH1D>("GSF_FEE_2GSF","H to GSF+FarEE & 2GSF Y",20, -5., 5.);
    hists.HGSF_HF_2GSF = fs->make<TH1D>("GSF_HF_2GSF","H to GSF+HF & 2GSF Y",20, -5., 5.);
    
    hists.Ztypes = fs->make<TH1D>("Ztypes", "Decay Channels of Z, Z*", 9, -0.5, 8.5 );
    hists.Ztypes->GetXaxis()->SetBinLabel(1, "Other");
    hists.Ztypes->GetXaxis()->SetBinLabel(2, "MuMu");
    hists.Ztypes->GetXaxis()->SetBinLabel(3, "MuE");
    hists.Ztypes->GetXaxis()->SetBinLabel(4, "EMu");
    hists.Ztypes->GetXaxis()->SetBinLabel(5, "EE");
    hists.Ztypes->GetXaxis()->SetBinLabel(6, "FarE-Mu");
    hists.Ztypes->GetXaxis()->SetBinLabel(7, "FarE-E");
    hists.Ztypes->GetXaxis()->SetBinLabel(8, "HF-Mu");
    hists.Ztypes->GetXaxis()->SetBinLabel(9, "HF-E");
    
    hists.AllForwardE = fs->make<TH1D>("AllForwardE","All Forward Electron Energies", 100, 0, 1000);
    hists.AllForwardPt = fs->make<TH1D>("AllForwardPt","All Forward Electron Pt", 50, 0, 150);
    
    hists.AllZ1e1Pt = fs->make<TH1D>("AllZ1e1Pt","All e1s from Z1 with Central Z2 Pt",50,0,150);
    hists.AllZ1e2Pt = fs->make<TH1D>("AllZ1e2Pt","All e2s from Z1 with Central Z2 Pt",50,0,150);
    hists.AllZ1ElecPt = fs->make<TH1D>("AllZ1ElecPt","All electrons from Z1 with Central Z2 Pt",50,0,150);
    hists.AllZ1e1Eta = fs->make<TH1D>("AllZ1e1Eta","All e1s from Z1 with Central Z2 Eta",48,-6,6);
    hists.AllZ1e2Eta = fs->make<TH1D>("AllZ1e2Eta","All e2s from Z1 with Central Z2 Eta",48,-6,6);
    hists.AllZ1ElecEta = fs->make<TH1D>("AllZ1ElecEta","All electrons from Z1 with Central Z2 Eta",48,-6,6);
    
    hists.e1PtVsEta = fs->make<TH2D>("e1PtVsEta","e1 Pt vs. Eta", 25, -5, 5, 50, 0, 150 );
    hists.e2PtVsEta = fs->make<TH2D>("e2PtVsEta","e2 Pt vs. Eta", 25, -5, 5, 50, 0, 150 );
    hists.AllPtVsEta = fs->make<TH2D>("AllPtVsEta","All electron Pt vs. Eta", 25, -5, 5, 50, 0, 150 );
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

    if ( (maxPt > minLep1pt && minPt > minLep2pt ) && 
	   ( fabs(theZ.l1eta) < 2.1 || fabs(theZ.l2eta) < 2.1) &&
       (fabs(theZ.l1eta) < 2.4 && fabs(theZ.l2eta) < 2.4) && theZ.mode==13) return true ;
       
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

    if ( (maxPt > minLep1pt && minPt > minLep2pt ) &&
       ( fabs(theZ.l1eta) < 2.5 && fabs(theZ.l2eta) < 2.5 ) && theZ.mode==11) return true ;
    
    return false ;
}

bool HiggsGenAnalysis::recZ(zboson theZ)//this checks if Z daughters could be reconstructed
{
	double minMupt = 5.;
	double minElpt = 7.;
	double maxMuEta = 2.4;
	double maxElEta = 3.;
	double maxHFeta = 5.;
	double maxGSFeta = 2.5;
	double minHFpt = 15.;
	double minZmass = 12;
	
	if (theZ.mode==13)
		return ((fabs(theZ.l1eta)<maxMuEta) && (fabs(theZ.l2eta)<maxMuEta) && (theZ.l1pt>minMupt) && (theZ.l2pt>minMupt));
		
	if ( (theZ.mode==11) && (theZ.mass>minZmass) && (fabs(theZ.l1eta)<maxGSFeta && theZ.l1pt>minElpt) )//e1 central
	{
		if( fabs(theZ.l2eta)<maxElEta && theZ.l2pt>minElpt ) return true;//e2 is GSF or FEE
		else if( fabs(theZ.l2eta)<maxHFeta && theZ.l2pt>minHFpt ) return true; //or e2 in HF
	}
	else if ( (theZ.mode==11) && (theZ.mass>minZmass) && (fabs(theZ.l2eta)<maxGSFeta && theZ.l2pt>minElpt) )//e2 central
		return ( fabs(theZ.l1eta)<maxHFeta && theZ.l1pt>minElpt ); //since e1 is already NOT central
		
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
	//hists.HiggsP2->Fill(Hclone.E());
	Hclone.Boost(-H.BoostVector());
	//hists.HiggsRestP2->Fill(Hclone.E());
	
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

  //std::cout << "-------------------------------------------------------" << std::endl ; 
  if (iEvent.getByLabel("genParticles",genInfo)) {
      reco::GenParticleCollection::const_iterator iparticle ;
      zboson z1, z2 ; 
      for (iparticle=genInfo->begin(); iparticle!=genInfo->end(); iparticle++) {

          if ( iparticle->pdgId() == 25 && iparticle->numberOfDaughters() > 0 ) { // Higgs
              //std::cout << "Higgs with mass " << iparticle->mass()
                        //<< " and rapidity " << iparticle->y()
                        //<< " decays to " << std::endl ;
              zboson newZ ;
              for (unsigned int idg=0; idg<iparticle->numberOfDaughters(); idg++) {//loop over Higgs products
                  if ( iparticle->daughter(idg)->pdgId() == 23 ) { // Z boson
                      int nele = 0 ; int nmu  = 0 ;
                      newZ.PZ = iparticle->daughter(idg)->p4();
                      for (unsigned int idg2=0; idg2<iparticle->daughter(idg)->numberOfDaughters(); idg2++) {//loop over each Z's products
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 11 ) {
                              nele++ ;
                              if ( nele == 1 ) {
                                  newZ.Pl1	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                              } 
                              else {
                                  newZ.Pl2	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                              }
                          }
                          if ( abs(iparticle->daughter(idg)->daughter(idg2)->pdgId()) == 13 ) {
                              nmu++ ; 
                              if ( nmu == 1 ) {
                                  newZ.Pl1	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                              }
                              else {
                                  newZ.Pl2	 = iparticle->daughter(idg)->daughter(idg2)->p4();
                              }
                          }
                      }
                      if ( nele == 2 || nmu == 2 ) {
                          nZ++ ; 
                          newZ.mass = iparticle->daughter(idg)->mass() ;
                          newZ.mode = ( (nele > 0) ? 11 : 13 ) ;
	                      if(newZ.Pl1.pt()>newZ.Pl2.pt())
	                      {   
	                          newZ.E1 = newZ.Pl1.E();
	                          newZ.E2 = newZ.Pl2.E();
	                          newZ.l1pt = newZ.Pl1.pt();
	                          newZ.l2pt = newZ.Pl2.pt();
	                          newZ.l1eta = newZ.Pl1.eta();
	                          newZ.l2eta = newZ.Pl2.eta();
	                      }
	                      else
	                      {
	                  		  newZ.E1 = newZ.Pl2.E();
	                          newZ.E2 = newZ.Pl1.E();
	                          newZ.l1pt = newZ.Pl2.pt();
	                          newZ.l2pt = newZ.Pl1.pt();
	                          newZ.l1eta = newZ.Pl2.eta();
	                          newZ.l2eta = newZ.Pl1.eta();
	                      }
	                      
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
                      Angular(z1,z2, iparticle->p4()); //angular distributions computed
                  }
              }
              
              // Determine a category
              int zCat = -1 ;
              if ( nZ == 2 )
              {
                  zCat = 0 ; 
                  if ( z1.mode == 13 ) //mu
                  {
                      if ( z1.mass > 50. && z1.mass < 120. && validZmumu(z1) && z2.mass > 12. ) //within cuts
                      {
                          if ( z2.mode == 13 && validZmumu(z2,true) ) zCat = 1 ;
                          if ( z2.mode == 11 && validZee(z2,true) )   zCat = 2 ;
                      }
                  }
                  else if ( z1.mode == 11 ) //el
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
                              if ( //(fabs(z1.l1eta) < 2.5 || fabs(z1.l2eta) < 2.5) 
                              		//&&(fabs(z1.l1eta) < 5.0 && fabs(z1.l2eta) < 5.0) 
                                   	 recZ(z1) )
                              {
                                  if ( z2.mode == 13 && validZmumu(z2,true) && (fabs(z1.l1eta) < 3.0 || fabs(z1.l2eta) < 3.0) ) zCat = 5 ;
                                  if ( z2.mode == 11 && validZee(z2,true) && (fabs(z1.l1eta) < 3.0 || fabs(z1.l2eta) < 3.0 ) )   zCat = 6 ;
                                  if ( z2.mode == 13 && validZmumu(z2,true) && (fabs(z1.l1eta) > 3.0 || fabs(z1.l2eta) > 3.0 ) ) zCat = 7 ;
                                  if ( z2.mode == 11 && validZee(z2,true) && (fabs(z1.l1eta) > 3.0 || fabs(z1.l2eta) > 3.0 ) )   zCat = 8 ;
                              }
                          }
                      }
                  }
              }
             
              if ( zCat >= 0 ) {
		passEvent = true ; 
		
		hists.h_higgsY->Fill( iparticle->y() ) ;
		hists.Ztypes->Fill(zCat);

		if ( recZ(z1) && recZ(z2) )
		{
			hists.PotRecoHiggs->Fill(iparticle->y());
			if( z1.mode==11 )
			{
				hists.e1PtVsEta->Fill(z1.l1eta,z1.l1pt);
				hists.e2PtVsEta->Fill(z1.l2eta,z1.l2pt);
				hists.AllPtVsEta->Fill(z1.l1eta,z1.l1pt);
				hists.AllPtVsEta->Fill(z1.l2eta,z1.l2pt);
			}
		}
		
		switch(zCat)
		{
			case 1:
				hists.H4Mu->Fill( iparticle->y() );
				break;
			case 2:
				hists.H2mu2GSF->Fill( iparticle->y() );
				break;
			case 3:
				hists.H2GSFe2mu->Fill( iparticle->y() ) ;
				break;
			case 4:
				hists.H4GSFe->Fill( iparticle->y() ) ;
				break;
			case 5:
				if( (fabs(z1.l1eta) > 2.5) && (fabs(z1.l1eta) < 3)  && (fabs(z1.l2eta) < 2.5) )
				{
					hists.HGSF_FEE_2mu->Fill( iparticle->y() ) ;
					hists.FarEEelEnergy->Fill(z1.E1);
					hists.FarEEelPt->Fill(z1.l1pt);
					hists.FarEEeta->Fill(z1.l1eta );
				}
				else if( (fabs(z1.l2eta) > 2.5) && (fabs(z1.l2eta) < 3)  && (fabs(z1.l1eta) < 2.5) )
				{
					hists.HGSF_FEE_2mu->Fill( iparticle->y() ) ;
					hists.FarEEelEnergy->Fill(z1.E2);
					hists.FarEEelPt->Fill(z1.l2pt);
					hists.FarEEeta->Fill(z1.l2eta );
				}
				break;
			case 6:
				if( (fabs(z1.l1eta) > 2.5) && (fabs(z1.l1eta) < 3) && (fabs(z1.l2eta) < 2.5) )
				{
					hists.HGSF_FEE_2GSF->Fill( iparticle->y() ) ;
					hists.FarEEelEnergy->Fill(z1.E1);
					hists.FarEEelPt->Fill(z1.l1pt);
					hists.FarEEeta->Fill(z1.l1eta);
				}
				else if( (fabs(z1.l2eta) > 2.5) && (fabs(z1.l2eta) < 3)  && (fabs(z1.l1eta) < 2.5) )
				{
					hists.HGSF_FEE_2GSF->Fill( iparticle->y() ) ;
					hists.FarEEelEnergy->Fill(z1.E2);
					hists.FarEEelPt->Fill(z1.l2pt);
					hists.FarEEeta->Fill(z1.l2eta);
				}
				break;
			case 7:
				if( (fabs(z1.l1eta) > 3) && (fabs(z1.l2eta) < 2.5) && (z1.l1pt > 15) )
				{
					hists.HGSF_HF_2mu->Fill( iparticle->y() ) ;
					hists.HFelEnergy->Fill(z1.E1);
					hists.HFelPt->Fill(z1.l1pt);
					hists.HFeta->Fill(z1.l1eta);
				}
				else if( (fabs(z1.l2eta) > 3) && (fabs(z1.l1eta) < 2.5) && (z1.l2pt > 15) )
				{
					hists.HGSF_HF_2mu->Fill( iparticle->y() ) ;
					hists.HFelEnergy->Fill(z1.E2);
					hists.HFelPt->Fill(z1.l2pt);
					hists.HFeta->Fill(z1.l2eta);
				}
				break;
			case 8:
				if( (fabs(z1.l1eta) > 3) && (fabs(z1.l2eta) < 2.5) && (z1.l1pt > 15) )
				{
					hists.HGSF_HF_2GSF->Fill( iparticle->y() ) ;
					hists.HFelEnergy->Fill(z1.E1);
					hists.HFelPt->Fill(z1.l1pt);
					hists.HFeta->Fill(z1.l1eta );
				}
				else if( (fabs(z1.l2eta) > 3) && (fabs(z1.l1eta) < 2.5) && (z1.l2pt > 15) )
				{
					hists.HGSF_HF_2GSF->Fill( iparticle->y() ) ;
					hists.HFelEnergy->Fill(z1.E2);
					hists.HFelPt->Fill(z1.l2pt);
					hists.HFeta->Fill(z1.l2eta);
				}
				break;
		}
		
		if( (z1.mode == 11) && recZ(z2) )
		{
			hists.AllZ1e1Pt->Fill(z1.l1pt);
			hists.AllZ1e2Pt->Fill(z1.l2pt);
			hists.AllZ1e1Eta->Fill(z1.l1eta);
			hists.AllZ1e2Eta->Fill(z1.l2eta);
		}
		
		if ( farElectronFilter_ ) { 
		  if ( z1.mode != 11 || z1.mass < 50. || z1.mass > 120. ) return false ; 
		  if ( fabs(z1.l1eta) < 2.5 && fabs(z1.l2eta) < 2.5 ) return false ; 
		  //std::cout << "*** This event has an HF/NT electron candidate ***" << std::endl ; 
		}
              }
          }
      }
      
  }
  //std::cout << "-------------------------------------------------------" << std::endl ; 

  return passEvent ; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
HiggsGenAnalysis::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiggsGenAnalysis::endJob() 
{
		hists.AllForwardE->Add(hists.FarEEelEnergy);
		hists.AllForwardE->Add(hists.HFelEnergy);
		hists.AllForwardPt->Add(hists.FarEEelPt);
		hists.AllForwardPt->Add(hists.HFelPt);
		
		hists.AllZ1ElecPt->Add(hists.AllZ1e1Pt);
		hists.AllZ1ElecPt->Add(hists.AllZ1e2Pt);
		hists.AllZ1ElecEta->Add(hists.AllZ1e1Eta);
		hists.AllZ1ElecEta->Add(hists.AllZ1e2Eta);

}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HiggsGenAnalysis);
