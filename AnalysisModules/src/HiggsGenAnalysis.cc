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
// $Id: HiggsGenAnalysis.cc,v 1.6 2012/05/30 01:03:50 afinkel Exp $
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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "HZZ4Leptons/AnalysisModules/src/HiggsCommon.h"

//
// class declaration
//

class HiggsGenAnalysis : public edm::EDFilter {
public:
  explicit HiggsGenAnalysis(const edm::ParameterSet&);
  ~HiggsGenAnalysis();
  
  struct zboson {
      int mode; // 11, 13
      double mass, Y;
      double l1pt, l2pt ; 
      double l1eta, l2eta ; 
      double E1, E2;
      reco::Particle::LorentzVector PZ, Pl1, Pl2;
  } ;
  
  struct HiggsGenEvent
  {
  	double  Hm, HY, Z1m, Z2m, Z1Y, Z2Y,
  			l1pt, l2pt, l3pt, l4pt,
  			l1eta, l2eta, l3eta, l4eta;
	int n_primary_vertex, n_pue;
  } he;
  

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
      TH1D *Phi1, *Phi, *CosTheta1, *CosTheta2, *CosTheta0, *Ztypes, *HY;

  } hists ;
    
  struct HistPerDef
  {
    	public:
        //book histogram set w/ common suffix inside the provided TFileDirectory
        void Book(TFileDirectory *, const std::string&, int type=0) ; 
      	void Fill(const HiggsGenEvent& he) ; 
        
      TH1  *HMass, *HY, *Z1mass, *Z2mass, *Z1Y, *Z2Y, *nVertex,
		   *l1Pt, *l1Eta, *l2Pt, *l2Eta, *l3Pt, *l3Eta, *l4Pt, *l4Eta,
      	   *EcalIso_1, *EcalIso_2, *EcalIsoByGSF_1, *EcalIsoByGSF_2,
      	   *E25Max_1, *E25Max_2, //*E25Max_3, *E25Max_4,
      	   *E15_1, *E15_2, //*E15_3, *E15_4,
      	   *E55_1, *E55_2, //*E55_3, *E55_4;
      	   *FarEEPt, *FarEEeta, *HadrOverEM, *sIeIe,
      	   *HFPt, *HFeta, *e9e25, *var2d,
      	   *FwdElPt, *FwdElEta;
      
      TH2 *AllPtVsEta, *l1PtVsEta, *l2PtVsEta;
  } ;
    
    TFileDirectory *rundir;    
        
    HistPerDef  H4Mu, H4GSFe, H2mu2GSF, HGSF_FEE_2mu, HGSF_HF_2mu, HGSF_HF_2GSF, HGSF_FEE_2GSF, H2GSFe2mu,
    			AllHF, AllFEE, AllFwd, AllPassing;
    			
    TH1 *Ztypes;
};
    
void HiggsGenAnalysis::HistPerDef::Book(TFileDirectory *mydir, const std::string& post, int type) 
  {
    std::string t, T; // histogram title string;
    TH1::SetDefaultSumw2();
    t = post + "_HY";
    T = post + " Gen H Rapidity";
    HY = mydir->make<TH1D> (t.c_str(), T.c_str(), 20, -5, 5 );
    t = post + "_HMass";
    T = post + " Gen H mass";
    std::cout<<"Created titles: "<<t<<", "<<T<<std::endl;        
    HMass = mydir->make<TH1D> (t.c_str(), T.c_str(), 50, 100, 200 );  
    t = post + "_Z1Mass";
    T = post + " Gen Z1 mass";
    Z1mass = mydir->make<TH1D>(t.c_str(), T.c_str(), 70, 50, 120 );  
    t = post + "_Z2Mass";
    T = post + " Gen Z2 mass";
    Z2mass = mydir->make<TH1D>(t.c_str(), T.c_str(), 90, 0, 90 );    
    t = post + "_Z1Y";
    T = post + " Gen Z1 Rapidity";
    Z1Y = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, -5, 5 );  
    t = post + "_Z2Y";
    T = post + " Gen Z2 Rapidity";
    Z2Y = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, -5, 5 );
    t = post + "_l1Pt";
    T = post + " Gen L1 Pt";
    l1Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 150 );    
    t = post + "_l2Pt";
    T = post + " Gen L2 Pt";
    l2Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 100 );
    t = post + "_l3Pt";
    T = post + " Gen L3Pt";
    l3Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 25, 0, 100 );    
    t = post + "_l4Pt";
    T = post + " Gen L4 Pt";
    l4Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 25, 0, 100 );
    t = post + "_l1Eta";
    T = post + " Gen L1 Eta";
    l1Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );    
    t = post + "_l2Eta";
    T = post + " Gen L2 Eta";
    l2Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );
    t = post + "_l3Eta";
    T = post + " Gen L3 Eta";
    l3Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );    
    t = post + "_l4Eta";
    T = post + " Gen L4 Eta";
    l4Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );
    t = post + "_nVertex";
    T = post + " Number of Vertices";
    nVertex = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 50 );
    
    t = post + "_AllPtVsEta";
    T = post + " l1, l2 Pt vs. Eta";
    AllPtVsEta = mydir->make<TH2D>(t.c_str(), T.c_str(), 25, -5, 5, 50, 0, 150 );
    t = post + "_l1PtVsEta";
    T = post + " Highest Pt vs. Eta";
    l1PtVsEta = mydir->make<TH2D>(t.c_str(), T.c_str(), 25, -5, 5, 50, 0, 150 );
    t = post + "_l2PtVsEta";
    T = post + " Second Pt vs. Eta";
    l2PtVsEta = mydir->make<TH2D>(t.c_str(), T.c_str(), 25, -5, 5, 50, 0, 150 );
  }
  
  void HiggsGenAnalysis::HistPerDef::Fill(const HiggsGenEvent& he) 
  {
    HY->Fill(he.HY);
    HMass->Fill(he.Hm);
    Z1mass->Fill(he.Z1m);
	Z2mass->Fill(he.Z2m);
	Z1Y->Fill(he.Z1Y);
	Z2Y->Fill(he.Z2Y);
	l1Pt->Fill(he.l1pt);
	l1Eta->Fill(he.l1eta);
	l2Pt->Fill(he.l2pt);
	l2Eta->Fill(he.l2eta);
	l3Pt->Fill(he.l3pt);
	l3Eta->Fill(he.l3eta);
	l4Pt->Fill(he.l4pt);
	l4Eta->Fill(he.l4eta);
	nVertex->Fill(he.n_primary_vertex);
	
	AllPtVsEta->Fill(he.l1eta,he.l1pt);
	AllPtVsEta->Fill(he.l2eta,he.l2pt);	
	if(he.l1pt>he.l2pt)
	{
		l1PtVsEta->Fill(he.l1eta,he.l1pt);
		l2PtVsEta->Fill(he.l2eta,he.l2pt);
	}
	else
	{	
		l1PtVsEta->Fill(he.l2eta,he.l2pt);
		l2PtVsEta->Fill(he.l1eta,he.l1pt);
	}
  }
  
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
    
    AllHF.Book(new TFileDirectory(fs->mkdir("AllHF")),"AllHF");
	AllFEE.Book(new TFileDirectory(fs->mkdir("AllFEE")),"AllFEE");
	AllFwd.Book(new TFileDirectory(fs->mkdir("AllFwd")),"AllFwd");
	AllPassing.Book(new TFileDirectory(fs->mkdir("AllPassing")),"AllPassing");
    H4Mu.Book(new TFileDirectory(fs->mkdir("4mu")),"4Mu");
    H4GSFe.Book(new TFileDirectory(fs->mkdir("4GSF")),"4GSF");
    H2mu2GSF.Book(new TFileDirectory(fs->mkdir("2Mu_2GSF")),"2Mu_2GSF");
    HGSF_FEE_2mu.Book(new TFileDirectory(fs->mkdir("GSF_FarEE_2Mu")),"GSF_FarEE_2Mu");
    HGSF_FEE_2GSF.Book(new TFileDirectory(fs->mkdir("GSF_FarEE_2GSF")),"GSF_FarEE_2GSF");
    HGSF_HF_2mu.Book(new TFileDirectory(fs->mkdir("GSF_HF_2Mu")),"GSF_HF_2Mu");
    HGSF_HF_2GSF.Book(new TFileDirectory(fs->mkdir("GSF_HF_2GSF")),"GSF_HF_2GSF");
    H2GSFe2mu.Book(new TFileDirectory(fs->mkdir("2GSF_2Mu_Cut0")),"2GSF_2Mu");
    
    hists.Phi = fs->make<TH1D>("Phi","#Phi",20, -3.2,3.2);
    hists.Phi1 = fs->make<TH1D>("Phi1","#Phi_{1}",20, -3.2,3.2);
    hists.CosTheta0 = fs->make<TH1D>("CostTheta0","Cos(#theta_{0})",50,-1.,1.);
    hists.CosTheta1 = fs->make<TH1D>("CostTheta1","Cos(#theta_{1})",20,-1.,1.);
    hists.CosTheta2 = fs->make<TH1D>("CostTheta2","Cos(#theta_{2})",20,-1.,1.);
       
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
    
    hists.HY = fs->make<TH1D > ("HY","Higgs rapidity", 100,-5.,5.) ;
    
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
	int nZ = 0 ; 
	bool passEvent = false ;  
	
	edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    he.n_primary_vertex = higgs::numberOfPrimaryVertices(pvHandle); 

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
                      
                      he.Hm = iparticle->mass();
                      he.HY = iparticle->y();
                      he.Z1m = z1.mass;
                      he.Z1Y = z1.PZ.Rapidity();
                      he.l1pt = z1.l1pt;
                      he.l1eta = z1.l1eta;
                      he.l2pt = z1.l2pt;
                      he.l2eta = z1.l2eta;
                      he.Z2m = z2.mass;
                      he.Z2Y = z2.PZ.Rapidity();
                      he.l3pt = z1.l1pt;
                      he.l3eta = z1.l1eta;
                      he.l4pt = z1.l2pt;
                      he.l4eta = z1.l2eta;
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
		
		hists.HY->Fill( iparticle->y() ) ;
		hists.Ztypes->Fill(zCat);
		AllPassing.Fill( he );

		/* if ( recZ(z1) && recZ(z2) )
		{
			hists.PotRecoHiggs->Fill(iparticle->y());
			if( z1.mode==11 )
			{
				hists.e1PtVsEta->Fill(z1.l1eta,z1.l1pt);
				hists.e2PtVsEta->Fill(z1.l2eta,z1.l2pt);
				hists.AllPtVsEta->Fill(z1.l1eta,z1.l1pt);
				hists.AllPtVsEta->Fill(z1.l2eta,z1.l2pt);
			}
		} */
		
		switch(zCat)
		{
			case 1:
				H4Mu.Fill( he );
				break;
			case 2:
				H2mu2GSF.Fill( he );
				break;
			case 3:
				H2GSFe2mu.Fill( he ) ;
				break;
			case 4:
				H4GSFe.Fill( he ) ;
				break;
			case 5:
				HGSF_FEE_2mu.Fill( he ) ;
				AllFEE.Fill( he );
				AllFwd.Fill( he );
				break;
			case 6:
				HGSF_FEE_2GSF.Fill( he ) ;
				AllFEE.Fill( he );
				AllFwd.Fill( he );
				break;
			case 7:
				HGSF_HF_2mu.Fill( he ) ;
				AllHF.Fill( he );
				AllFwd.Fill( he );
				break;
			case 8:
				HGSF_HF_2GSF.Fill( he ) ;
				AllHF.Fill( he );
				AllFwd.Fill( he );
				break;
		}
		
		/*if( (z1.mode == 11) && recZ(z2) )
		{
			hists.AllZ1e1Pt->Fill(z1.l1pt);
			hists.AllZ1e2Pt->Fill(z1.l2pt);
			hists.AllZ1e1Eta->Fill(z1.l1eta);
			hists.AllZ1e2Eta->Fill(z1.l2eta);
		} */
		
		if ( farElectronFilter_ ) { 
		  if ( z1.mode != 11 || z1.mass < 50. || z1.mass > 120. ) return false ; 
		  if ( fabs(z1.l1eta) < 2.5 && fabs(z1.l2eta) < 2.5 ) return false ; 
		  //std::cout << "*** This event has an HF/NT electron candidate ***" << std::endl ; 
		}
              }
          }
      }
      
  }

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
		/* hists.AllForwardE->Add(hists.FarEEelEnergy);
		hists.AllForwardE->Add(hists.HFelEnergy);
		hists.AllForwardPt->Add(hists.FarEEelPt);
		hists.AllForwardPt->Add(hists.HFelPt);
		
		hists.AllZ1ElecPt->Add(hists.AllZ1e1Pt);
		hists.AllZ1ElecPt->Add(hists.AllZ1e2Pt);
		hists.AllZ1ElecEta->Add(hists.AllZ1e1Eta);
		hists.AllZ1ElecEta->Add(hists.AllZ1e2Eta); */

}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HiggsGenAnalysis);
