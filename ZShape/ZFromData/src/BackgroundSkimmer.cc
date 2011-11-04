#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/ZFromData/interface/BackgroundSkimmer.h"
///#include "ZShape/Base/interface/ZShapeUtils.h"

#include <cmath>
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "PhysicsTools/TagAndProbe/interface/CandidateAssociation.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include <TF1.h>
#include <TMath.h>
#include <math.h>
#include <iostream>
#include <fstream>


Double_t CrystalBall(Double_t *x,Double_t *par);
//
// constructors and destructor
//
BackgroundSkimmer::BackgroundSkimmer(const edm::ParameterSet& iConfig) 
{
  // 
  //now do what ever initialization is needed

  //This little guy reads in the differnt InputTag collections.
  const edm::InputTag dBlankTag("Blank");
  
  std::vector< edm::InputTag > defaultAllProbeCandTags;
  allProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
      "allProbeCandTags",defaultAllProbeCandTags);
   for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
      std::cout << allProbeCandTags_[i] << std::endl;
  
}


BackgroundSkimmer::~BackgroundSkimmer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool BackgroundSkimmer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  using namespace reco;
 
  //if (fromfile_) return;
  
  edm::Handle<reco::CandViewCandViewAssociation> tagprobes;
      if ( !iEvent.getByLabel(allProbeCandTags_[0],tagprobes) ) {
         //std::cout << " DIDN't get the darn tag " << allProbeCandTags_[i] << std::endl;
      }

  bool goodevent = false;

  if (tagprobes.isValid()) if (tagprobes->size() > 0) goodevent = true;

  return goodevent;


}




// ------------ method called once each job just before starting event loop  ------------
void 
BackgroundSkimmer::beginJob(const edm::EventSetup&)
{



  
}




// ------------ method called once each job just after ending the event loop  ------------
void 
BackgroundSkimmer::endJob() {

 


}


// ************************************************************** //

//define this as a plug-in
//DEFINE_ANOTHER_FWK_MODULE(BackgroundSkimmer);
