// -*- C++ -*-
//
// Package:    SCwithPUData
// Class:      SCwithPUData
// 
/**\class SCwithPUData SCwithPUData.cc ECALTime/SCwithPUData/src/SCwithPUData.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Franzoni,27 2-013,+41227678347,
//         Created:  Mon Jun 20 15:07:58 CEST 2011
// $Id: SCwithPUData.cc,v 1.2 2011/07/08 21:39:11 franzoni Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>

#include <iostream>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include <TMath.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "Hgg/ClusteringWithPU/interface/SelectElectron.h"

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


// gives ratio of PU-cleaned energy over initial energy 
float removePU(const pat::ElectronCollection::const_iterator oneEle, const float xi){
  reco::SuperClusterRef scr=oneEle->superCluster();

  if (scr.isNull()) {        std::cout  << "removePU: reference to superlcuster is NULL: we have a problem" << std::endl; assert(0); }
  
  float seedBCEnergy      = (scr->seed())->energy(); // this should be replaced by the 5x5 around the seed, to reproduce earlier study
  float eSeed             = 0.35;                    // standard eSeed in EB 
  float cumulateRawEnergy = 0;
  
  // looping on basic clusters within the Supercluster
  for(reco::CaloCluster_iterator bcIt = scr->clustersBegin(); bcIt!=scr->clustersEnd(); bcIt++)
    {
      if( (*bcIt)->energy() > sqrt( eSeed*eSeed + xi*xi*seedBCEnergy*seedBCEnergy/cosh((*bcIt)->eta())/cosh((*bcIt)->eta())  ) ) {
	cumulateRawEnergy+=(*bcIt)->energy();
      }
    }

  if(0) std::cout  << "cumulateRawEnergy / scr->rawEnergy() " << cumulateRawEnergy / scr->rawEnergy() << " xi being: " << xi << " num_BC:  " << scr->clustersSize() << " seed/rawEnergy: " << (scr->seed())->energy() << "/" << scr->rawEnergy() << std::endl;
  return (cumulateRawEnergy / scr->rawEnergy() );

}




//
// class declaration
//

class SCwithPUData : public edm::EDAnalyzer {
public:
  explicit SCwithPUData(const edm::ParameterSet&);
  ~SCwithPUData();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  
  edm::InputTag vertexCollection_ ;
  edm::InputTag patElectrons_ ;


  std::string ecalID_, currentFile_;
  bool dolog_;
  std::string myName_;
  float ETCut_, etaMin_, etaMax_;

  std::string workingPoint_;
  std::vector<int> acceptedElectronIDs_;
  SelectElectron elSelector; 

  float xi_;

  // ----------member data ---------------------------
  
  struct HistSet {

    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(edm::Service<TFileService>& td,const std::string&,float xi);

    // fill all histos of the set with the two electron candidates
    void fill(const pat::ElectronCollection::const_iterator leadEle, 
	      const pat::ElectronCollection::const_iterator subEle,
	      const reco::VertexCollection * theRecVtxs  );

    TH1 *nElec, *nVertices;

    TH1* eleLEta,   *eleLPhi,   *eleLEt, *eleLNumBC, *eleLFrac;
    TH1* eleSubEta, *eleSubPhi, *eleSubEt, *eleSubNumBC, *eleSubFrac;
    TH1* massPlots[20];
    TH1* massPlotsPUOrig[20];
    TH1* massPlotsPUMod[20];
    TH2* massVsVertex;

    float theXi;

  } theHists;
  
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
SCwithPUData::SCwithPUData(const edm::ParameterSet& iConfig)
  
{
  //now do what ever initialization is needed
  dolog_              =iConfig.getParameter<bool>("DoLog");
  myName_             =iConfig.getParameter<std::string>("myName");
  workingPoint_       =iConfig.getParameter<std::string>("eleWorkingPoint");

  ETCut_              = iConfig.getParameter<double>("ETCut");
  etaMin_             = iConfig.getParameter<double>("etaMin");
  etaMax_             = iConfig.getParameter<double>("etaMax");
  vertexCollection_   = iConfig.getParameter<edm::InputTag> ("vertexCollection");
  patElectrons_       = iConfig.getParameter<edm::InputTag> ("patElectrons");
  
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    elSelector.add(*it);    }
  
  xi_ = iConfig.getParameter<double>("xi");
  

  edm::Service<TFileService> fs;
  
  // book all the histograms of the structure
  theHists.book( fs, "", xi_ );
  
  std::cout << ">>> this is SCwithPUData instance: " << myName_ << " eleWP: " << workingPoint_ << " xi: " << xi_ << std::endl;

}


SCwithPUData::~SCwithPUData()
{ }


void SCwithPUData::HistSet::book(edm::Service<TFileService> &td, const std::string& post, float xi) {
  std::string title;

  title=std::string("num electrons ")+post+std::string("; num electrons");
  nElec=td->make<TH1D>("num electrons","num electrons; num electrons",10,-0.5,9.5);
  title=std::string("num vertices ")+post+std::string("; num vertices");
  nVertices=td->make<TH1D>("num vertices","num vertices; num vertices",41,-0.5,40.5);

  title=std::string("eta_{elLead} ")+post+std::string(";#eta_{elLead}");
  eleLEta=td->make<TH1D>("eta_{elLead}",title.c_str(),60,-3,3);
  title=std::string("phi_{elLead} ")+post+std::string(";#phi_{elLead}");
  eleLPhi=td->make<TH1D>("phi_{elLead}",title.c_str(),60,-1*TMath::Pi(),TMath::Pi());
  title=std::string("et_{elLead} ")+post+std::string(";Et_{elLead} [GeV]");
  eleLEt=td->make<TH1D>("et_{elLead}",title.c_str(),120,0,120);
  title=std::string("frac_{elLead} ")+post+std::string(";frac_{elLead}");
  eleLFrac=td->make<TH1D>("frac_{elLead}",title.c_str(),1100,0,1.1);

  title=std::string("eta_{elSub} ")+post+std::string(";#eta_{elSub}");
  eleSubEta=td->make<TH1D>("eta_{elSub}",title.c_str(),60,-3,3);
  title=std::string("phi_{elSub} ")+post+std::string(";#phi_{elSub}");
  eleSubPhi=td->make<TH1D>("phi_{elSub}",title.c_str(),60,-1*TMath::Pi(),TMath::Pi());
  title=std::string("et_{elSub} ")+post+std::string(";Et_{elSub} [GeV]");
  eleSubEt=td->make<TH1D>("et_{elSub}",title.c_str(),120,0,120);
  title=std::string("frac_{elSub} ")+post+std::string(";frac_{elSub}");
  eleSubFrac=td->make<TH1D>("frac_{elSub}",title.c_str(),1100,0,1.1);
    
  title=std::string("numBC_{elLead} ")+post+std::string("; numBC_{elLead}");
  eleLNumBC=td->make<TH1D>("numBC_{elLead}",title.c_str(),20,0,20);
  title=std::string("numBC_{elSub} ")+post+std::string("; numBC_{elSub}");
  eleSubNumBC=td->make<TH1D>("numBC_{elSub}",title.c_str(),20,0,20);

  for(uint vertex=1; vertex<=20; vertex++){
    title=std::string("m_{ee} num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))+post+std::string(";m_{ee} [GeV/c^{2}]");
    massPlots[vertex-1]=td->make<TH1D>( (std::string("m_{ee} num. vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),140,50,120);
    title=std::string("m_{ee} PU orig  num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))+post+std::string(";m_{ee} [GeV/c^{2}]");
    massPlotsPUOrig[vertex-1]=td->make<TH1D>( (std::string("m_{ee} PU orig num. vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),140,50,120);
    title=std::string("m_{ee} PU mod num vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))+post+std::string(";m_{ee} [GeV/c^{2}]");
    massPlotsPUMod[vertex-1]=td->make<TH1D>( (std::string("m_{ee} PU mod num. vertices: ")+convertInt(2*(vertex-1))+std::string("-")+convertInt(2*(vertex))).c_str() ,title.c_str(),140,50,120);

  }
  massVsVertex=td->make<TH2D>("m_{ee} Vs NumVertex","m_{ee} Vs NumVertex; m_{ee} [GeV/c^{2}]; NumVertex  ",40,0,40,140,50,120);

  theXi = xi;

}


void SCwithPUData::HistSet::fill(const pat::ElectronCollection::const_iterator leadEle, const pat::ElectronCollection::const_iterator subEle,
				 const reco::VertexCollection * theRecVtxs)
{
  
  float fractionLead = removePU(leadEle,theXi);
  float fractionSub  = removePU(subEle,theXi);
  
  eleLEta -> Fill(leadEle->eta());
  eleLPhi -> Fill(leadEle->phi());
  eleLEt  -> Fill(leadEle->pt());
  eleLNumBC->Fill( (leadEle->superCluster())->clustersSize() ); 
  eleLFrac-> Fill(fractionLead);
  
  eleSubEta -> Fill(subEle->eta());
  eleSubPhi -> Fill(subEle->phi());
  eleSubEt  -> Fill(subEle->pt());
  eleSubNumBC->Fill( (subEle->superCluster())->clustersSize() ); 
  eleSubFrac-> Fill(fractionSub);

  if( 0 ) std::cout  <<  " number of vertices: " << theRecVtxs->size() << std::endl;
  
  reco::Particle::LorentzVector Z( fractionLead * leadEle->p4());
  Z+= ( fractionSub * subEle->p4() );

  reco::Particle::LorentzVector Zraw( leadEle->p4() );
  Zraw+= ( subEle->p4() );



  if( theRecVtxs->size()>0 && theRecVtxs->size()<=40) {
    massPlots[( (theRecVtxs->size()-1) /2 )] -> Fill(Z.M());
    if(   fractionLead * fractionSub <1  ){
      massPlotsPUOrig[( (theRecVtxs->size()-1) /2 )] -> Fill(Zraw.M());
      massPlotsPUMod[( (theRecVtxs->size()-1) /2 )]  -> Fill(Z.M());
    }
  }
  else {
    massPlots[ 19 ] -> Fill(Z.M());
    if(  ( fractionLead * fractionSub ) <1  ){
      massPlotsPUMod[ 19 ]  -> Fill(Zraw.M());
      massPlotsPUOrig[ 19 ] -> Fill(Z.M());
    }
  }
  massVsVertex->Fill( theRecVtxs->size(), Z.M() );

}// end of fill()


//
// member functions
//

// ------------ method called for each event  ------------
void
SCwithPUData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   edm::Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(vertexCollection_, recVtxs);
   const reco::VertexCollection * theRecVtxs = recVtxs.product();
   
   // get collection of pat electrons
   edm::Handle<pat::ElectronCollection> patElectrons;
   iEvent.getByLabel(patElectrons_, patElectrons);
   if ( ! patElectrons.isValid()) {      
     std::cout << myName_ << " No electrons found in this event with tag "  << patElectrons << std::endl;     return ;    }
   else {      if( dolog_) std::cout  << myName_ << " got electrons collection with tag: " << patElectrons_ << "\t"<< patElectrons << std::endl;    }
   const pat::ElectronCollection& pElecs = *(patElectrons.product());
   
   // if there's less than TWO electrons in the collection, bail out
   unsigned int numMinElectrons=2;
   if(pElecs.size()<numMinElectrons)
     { if( dolog_) std::cout  << myName_ << " Too few electrons (less than " << numMinElectrons << "  - skip to next event" << std::endl;        return;}
   else     { if( dolog_) std::cout  << myName_ << " Enough electrons found,  they're " << pElecs.size() << std::endl;}

   
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // looking for leading electron
   pat::ElectronCollection::const_iterator i;
   pat::ElectronCollection::const_iterator theLeadingEle=pElecs.begin();
   bool leadHasBeenFound(false); float leadElePt=0;
   for (i=pElecs.begin(); i!=pElecs.end(); i++) {
     
     reco::SuperClusterRef scr=i->superCluster();
     if (scr.isNull()) {        std::cout  << myName_ << " reference to superlcuster is NULL: we have a problem" << std::endl;        continue;     }
     double eta_det=scr.get()->eta();
     if( dolog_ ) std::cout  << myName_ << " got ONE electron with SC-pt: " << i->pt() << " and with ID : " <<  i->electronID( workingPoint_.c_str()  ) << std::endl; 
     
     // gf: ECAL eta acceptance cut on supercluster. What about phi acceptance and central crack?
     // explanation of electronID value: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
     if (	 !(  ( fabs(eta_det)<1.4442 || fabs(eta_det)>1.560 )    // stay away from cracks
		     && fabs(eta_det)<2.5
		     && i->pt() > ETCut_ 
		     && elSelector.doesElePass( i->electronID( workingPoint_.c_str()  )	)
		     && etaMin_<fabs(eta_det)    // here you limit range from .py file; 
		     && fabs(eta_det)<etaMax_    // e.g. to limit to the EE...
		     )
		 ) continue;
     
     if( i->pt() >= leadElePt ) {theLeadingEle=i; leadElePt=i->pt(); leadHasBeenFound=true;}
     if( dolog_ && leadHasBeenFound) std::cout  << myName_ << " got ONE leading electron passing cuts with SC-pt: " << i->pt() << std::endl; 
     
   }// end loop over electrons to look for leading electron 
   


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // looking for sub-leading electron
   pat::ElectronCollection::const_iterator theSubLEle;//=pElecs.begin();
   bool hasBeenFound(false); float subElePt=0;
   for (i=pElecs.begin(); i!=pElecs.end(); i++) {
     reco::SuperClusterRef scr=i->superCluster();
     if (scr.isNull()) {        
       std::cout  << myName_ << " reference to superlcuster is NULL: we have a problem" << std::endl;        continue;     }
     
     double eta_det=scr.get()->eta();
     if (	 !(  ( fabs(eta_det)<1.4442 || fabs(eta_det)>1.560 )   && fabs(eta_det)<2.5
		     && i->pt() > ETCut_ 
		     && elSelector.doesElePass( i->electronID( workingPoint_.c_str()  ) )
		     && etaMin_<fabs(eta_det)
		     && fabs(eta_det)<etaMax_ 
		     )
		 ) continue;
     
     if( ( i->pt() >= subElePt ) && (i!=theLeadingEle) ) {theSubLEle=i; subElePt=theSubLEle->pt(); hasBeenFound=true;}
     if( dolog_ && hasBeenFound) std::cout  << myName_ << " got ONE sub-ele eelectron passing cuts with SC-pt: " << i->pt() << std::endl; 
   }// end loop over electrons to look for sub-leading electron 
   
   
   // no suitable sub-leading electron found, bail out
   if( (!hasBeenFound)  || (!leadHasBeenFound) ) { 
     if( dolog_) std::cout  << myName_ << " no subleading electron found - skip to next event" << std::endl; 
     return; 
   }
   
   
   if (0) std::cout  << myName_ << " xi: " << xi_ << " leading electron pt: " << theLeadingEle->pt() 
		     << " fraction leading is: " << removePU(theLeadingEle,xi_) 
		     << " subleading pt: " << theSubLEle->pt() << " fraction is: " << removePU(theSubLEle,xi_) << std::endl;

   theHists.nElec    ->Fill(pElecs.size());
   theHists.nVertices->Fill(theRecVtxs->size());
   theHists.fill(theLeadingEle,theSubLEle,theRecVtxs);
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
SCwithPUData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCwithPUData::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
SCwithPUData::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SCwithPUData::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SCwithPUData::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SCwithPUData::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SCwithPUData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCwithPUData);
