#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/ZFromData/interface/InclusiveCheck.h"
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
InclusiveCheck::InclusiveCheck(const edm::ParameterSet& iConfig) :
  m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source"))),
  quiet_(iConfig.getUntrackedParameter<bool>("quiet",false))
{
  // 
  //now do what ever initialization is needed


  outFileName_              = iConfig.getUntrackedParameter<std::string>("outHistogramsFile", "");
  writeHistoConservatively_ = iConfig.getUntrackedParameter<bool>("writeHistoBeforeEndJob", false);
  delRMatchingCut_          = iConfig.getUntrackedParameter<double>("dRMatchCut",0.2);
  delPtRelMatchingCut_      = iConfig.getUntrackedParameter<double>("dPtMatchCut",15.0);

  textBinsFile_      = iConfig.getUntrackedParameter< std::string >("EffBinsFile","EffBinsFile.txt");
  effBinsFromTxt_    = new EffTableLoader(textBinsFile_);

  numMatch_          =  iConfig.getUntrackedParameter<std::vector<int > >("NumMatch");

  //This little guy reads in the differnt InputTag collections.
  const edm::InputTag dBlankTag("Blank");
  
  std::vector< edm::InputTag > defaultAllProbeCandTags;
  allProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
      "allProbeCandTags",defaultAllProbeCandTags);
   for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
      std::cout << allProbeCandTags_[i] << std::endl;
  
}


InclusiveCheck::~InclusiveCheck()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void InclusiveCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  using namespace reco;
 
  //if (fromfile_) return;

  ///FIRST FIND THE NUMBER FOR HLT... or ignore it for the time being
  bool goodevent = true;
  for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
  {
      edm::Handle<reco::CandidateView> thecand;
      if ( !iEvent.getByLabel(allProbeCandTags_[i],thecand) ) {
         //std::cout << " DIDN't get the darn tag " << allProbeCandTags_[i] << std::endl;
      }
      //put stuff here to get the TP collection from the event....
      edm::Handle<reco::CandViewCandViewAssociation> tagprobes;
      if ( !iEvent.getByLabel(allProbeCandTags_[i],tagprobes) ) {
         //std::cout << " DIDN't get the darn tag " << allProbeCandTags_[i] << std::endl;
      }
      int thesize = 0;
      if (thecand.isValid())   thesize = thecand->size();
      if (tagprobes.isValid()) thesize = tagprobes->size();

      if ( thesize > numMatch_[i] & goodevent ){ seqHist_->Fill(i); }//Maybe only do this by checking if previous ones existed
      else goodevent = false;
  }



}




// ------------ method called once each job just before starting event loop  ------------
void 
InclusiveCheck::beginJob(const edm::EventSetup&)
{
  // smearing the target efficiency according to statistics 
  ///if (statsBox_.trials>0) createAlternateEfficiencies();
  
  // file to write out the histograms produced by the InclusiveCheck
  ///histoFile_   = new TFile(outFileName_.c_str(),"RECREATE");

  edm::Service<TFileService> fs;
  seqHist_  = fs->make<TH1F>("SeqHist","SeqHist", 10, 0, 10);


  
}




// ------------ method called once each job just after ending the event loop  ------------
void 
InclusiveCheck::endJob() {

 


}



int InclusiveCheck::ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, 
					      edm::Handle<reco::CandidateView>& passprobes, bool checkExactOverlap_ )
{
   int ppass = 0;
   if( passprobes.isValid() )
   {
      for( int ipp=0; ipp<(int)passprobes->size(); ++ipp )
      {
        if (((*passprobes)[ipp]).et() < scEtCut_ ) continue;
	bool isOverlap = MatchObjects(&((*passprobes)[ipp]),  probe,checkExactOverlap_);

	reco::SuperClusterRef probeSC;
	reco::SuperClusterRef passprobeSC; 

	if( not probe.isNull() ) probeSC  =  probe->get<reco::SuperClusterRef>();

	reco::CandidateBaseRef ref = passprobes->refAt(ipp);

	if( not ref.isNull() ) passprobeSC = ref->get<reco::SuperClusterRef>();
	
	//isOverlap = isOverlap && ( probeSC == passprobeSC );
        if( isOverlap ) 
        {
           int myindex = effBinsFromTxt_->GetBandIndex(probe->pt(), probe->eta());
            if ((myindex == 1 || myindex == 3) )  HMCptvsSCpt_->Fill(probe->pt(),(*passprobes)[ipp].et());
           SCtoMCERatio_[myindex]->Fill(((*passprobes)[ipp].et())/(probe->pt()));
        }

	if( isOverlap ) ppass = 1;
      } 
   }

   return ppass;
}
// ************************************************************ //


// ***************** Trigger object matching ******************** //
bool InclusiveCheck::MatchObjects( const reco::Candidate *hltObj, 
				      const reco::CandidateBaseRef& tagObj,
				      bool exact )
{
   double tEta = tagObj->eta();
   double tPhi = tagObj->phi();
   double tPt  = tagObj->pt();
   double hEta = hltObj->eta();
   double hPhi = hltObj->phi();
   double hPt  = hltObj->pt();

   double dRval = deltaR(tEta, tPhi, hEta, hPhi);
   double dPtRel = 999.0;
   if( tPt > 0.0 ) dPtRel = fabs( hPt - tPt )/tPt;

   // If we are comparing two objects for which the candidates should
   // be exactly the same, cut hard. Otherwise take cuts from user.
   if( exact ) return ( dRval < 1e-3 && dPtRel < 1e-3 );
   else        return ( dRval < delRMatchingCut_ && dPtRel < delPtRelMatchingCut_ );
}
// ************************************************************** //

//define this as a plug-in
//DEFINE_ANOTHER_FWK_MODULE(InclusiveCheck);
