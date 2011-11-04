#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/ZFromData/interface/ZFromData.h"
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
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"



//
// constructors and destructor
//
ZFromData::ZFromData(const edm::ParameterSet& iConfig) :
  m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source"))),
  quiet_(iConfig.getUntrackedParameter<bool>("quiet",false))
{
  // 
  //now do what ever initialization is needed

  outFileName_              = iConfig.getUntrackedParameter<std::string>("outHistogramsFile", "");
  writeHistoConservatively_ = iConfig.getUntrackedParameter<bool>("writeHistoBeforeEndJob", false);
  delRMatchingCut_          = iConfig.getUntrackedParameter<double>("dRMatchCut",0.2);
  delPtRelMatchingCut_      = iConfig.getUntrackedParameter<double>("dPtMatchCut",15.0);

  //
  //Extra histos generated from the "From Data trials"
  extraHistos_ = iConfig.getUntrackedParameter<bool>("ExtraFromDataHistos", false);

  //
  // multiple PSets==Zdefinitions: each electron required to pass a set of cuts

  std::vector<edm::ParameterSet> zdefSetups =iConfig.getUntrackedParameter< std::vector<edm::ParameterSet> >("ZDefs");

  for (std::vector<edm::ParameterSet>::iterator i=zdefSetups.begin(); i!=zdefSetups.end(); ++i) {
    std::string name=i->getUntrackedParameter<std::string>("name");
    std::vector<std::string> req1=i->getUntrackedParameter<std::vector<std::string> >("e1");
    std::vector<std::string> req2=i->getUntrackedParameter<std::vector<std::string> >("e2");
    std::vector<std::string> reqZ=i->getUntrackedParameter<std::vector<std::string> >("Z");
    ZShapeZDef* zdef=new ZShapeZDef();
    for (std::vector<std::string>::iterator k=req1.begin(); k!=req1.end(); k++) 
      zdef->addCriterion(ZShapeZDef::crit_E1,*k);
    for (std::vector<std::string>::iterator k=req2.begin(); k!=req2.end(); k++) 
      zdef->addCriterion(ZShapeZDef::crit_E2,*k);
    for (std::vector<std::string>::iterator k=reqZ.begin(); k!=reqZ.end(); k++) 
      zdef->addCriterion(ZShapeZDef::crit_Z,*k);
    zdefs_[name]=zdef;
    zptorder_[name] = i->getUntrackedParameter<bool>("ptorder", false);
  }

  //
  // list of all available efficiencies
  //MAYBE REMOVE THIS... WHAT DO I NEED IT FOR????
  
  std::vector<edm::ParameterSet> effs =iConfig.getUntrackedParameter< std::vector<edm::ParameterSet> >("Effs");

  //for (std::vector<edm::ParameterSet>::iterator i=effs.begin(); i!=effs.end(); ++i) {
  //  std::string name=i->getUntrackedParameter<std::string>("name");
  //  std::string file=i->getUntrackedParameter<std::string>("effFile");
    //std::string var=i->getUntrackedParameter<std::string>("variable");
  //   loadEfficiency(name,file);
  // }


  //
  // Set up the Z from the data initializer
  //This is the Zee association map that we will start the efficiencies on
  tnpProducer_ =  iConfig.getUntrackedParameter<edm::InputTag >("TagProbeProducer");
   
  //Collection of GsfElectrons greater than 20 GeV pt.
  gsfProducer_ =  iConfig.getUntrackedParameter<edm::InputTag >("GsfProducer");
  //The Isolation & HLT  Parameters
  ///_HLTpassed.Init(iConfig);
  ///_Isopassed.Init(iConfig);

  random_ = new TRandom3();
  
  // ********** Tag-Probes ********** //
  //std::vector< edm::InputTag > defaultTagProbeMapTags;
  //tagProbeMapTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
  //    "tagProbeMapTags",defaultTagProbeMapTags);
  //for( int i=0; i<(int)tagProbeMapTags_.size(); ++i )
  //    std::cout << tagProbeMapTags_[i] << std::endl;
  // ******************************** //
  //Just for now allow more than 1 TagandProbe to be read in...
  
  // ********** Candiates ********** //
  //int map_size = (int)tagProbeMapTags_.size();
  const edm::InputTag dBlankTag("Blank");
  
  std::vector< edm::InputTag > defaultAllProbeCandTags;
   allProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
      "allProbeCandTags",defaultAllProbeCandTags);
   for( int i=0; i<(int)allProbeCandTags_.size(); ++i ) 
      std::cout << allProbeCandTags_[i] << std::endl;
   int map_size = (int)allProbeCandTags_.size();
   // Make sure the arrays won't cause a seg fault!
   //if( (int)allProbeCandTags_.size() < map_size )
   //{
   //   std::cout << "Warning: Number of TagProbe maps bigger than number of tag arrays!" << std::endl;
   //   for( int i=0; i<(map_size-(int)allProbeCandTags_.size()); ++i ) 
//	 allProbeCandTags_.push_back(dBlankTag);
   //} 

   std::vector< edm::InputTag > defaultPassProbeCandTags;
   passProbeCandTags_ = iConfig.getUntrackedParameter< std::vector<edm::InputTag> >(
      "passProbeCandTags",defaultPassProbeCandTags);
   for( int i=0; i<(int)passProbeCandTags_.size(); ++i ) 
      std::cout << passProbeCandTags_[i] << std::endl;
   // Make sure the arrays won't cause a seg fault!
   if( (int)passProbeCandTags_.size() < map_size )
   {
      std::cout << "Warning: Number of TagProbe maps bigger than number of probe arrays!" << std::endl;
      for( int i=0; i<(map_size-(int)passProbeCandTags_.size()); ++i ) 
	 passProbeCandTags_.push_back(dBlankTag);
   } 
   // ******************************************* //
   
   
   // **********Cut Names for Candiates************ //
   cutName_ = iConfig.getUntrackedParameter< std::vector<std::string> >("CutNames");
   // ******************************************* //
   
   
   // **********Exact Matching for Candiates************ //
   exactMatch_ = iConfig.getUntrackedParameter< std::vector<int> >("ExactMatch");
   // ******************************************* //
   

   std::cout << " Done with the constructor " << std::endl;
}


ZFromData::~ZFromData()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void ZFromData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  
  //std::cout << " Start 0 " << std::endl;
  Handle< HepMCProduct > EvtHandle ;
   
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  iEvent.getByLabel(m_srcTag, EvtHandle ) ;
  //std::cout << " Start 1 " << std::endl;
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
   //std::cout << " Start 2 " << std::endl;
  fillMCEvent(Evt);
  //evt_.findDataZ(iEvent);
  //int numDataZ = 0;
  //std::cout << " Clearing an Event " << std::endl;
  evt_.clear();
  //std::cout << " Checking an Event " << std::endl;
   
  //First set the number of Gsf Electrons above 20 GeV
  edm::Handle< reco::GsfElectronCollection> electronHandle;
  try
    {
      iEvent.getByLabel(gsfProducer_,electronHandle);
    }
  catch(cms::Exception &ex)
    {
      edm::LogError("ZFromData::GsfElectron ") << "Error! Can't get collection " << 
	gsfProducer_;
      throw ex;
    }
   evt_.n_gsf20 = electronHandle->size();

  
  // Fill some information about the tag & probe collections
  //int nrtp = 0;
  ///for( int itype=0; itype<(int)tagProbeMapTags_.size(); ++itype )
  ///{
  // Tag-Probes
  //int ptype=0;
  edm::Handle<reco::CandViewCandViewAssociation> tagprobes;
  if ( !iEvent.getByLabel(tnpProducer_,tagprobes) ) {
	 LogWarning("ZFromData") << "Could not extract tag-probe map with input tag "  << tnpProducer_;
         std::cout << "Didn't get the input tag " << std::endl;
     }

     if( tagprobes.isValid() )
      {
	  std::cout << "In tag probes " << std::endl;
	  //--Here I pick a random T&Probe combination based on the event number
	  int tpsize = tagprobes->size();
	  int usetp = (tpsize > 0 ) ? (iEvent.id().event())%tpsize : 0;
	  int tpnum = 0;
	  //--End picking a random Tag and Probe combination
	  std::cout << " tpnum " << tpnum << " usetp " << usetp << " tpsize " << tpsize << std::endl;
	  reco::CandViewCandViewAssociation::const_iterator tpItr = tagprobes->begin();
	  for( ; tpItr != tagprobes->end(); ++tpItr,++tpnum)
	  { 
	    std::cout << " tpnum is " << tpnum <<  " usetp is " << usetp << std::endl;
	    if (tpnum != usetp ) continue;
            std::cout << "In tag probes Iterator " << std::endl;
	    const reco::CandidateBaseRef &tag = tpItr->key;
	    vector< pair<reco::CandidateBaseRef,double> > vprobes = (*tagprobes)[tag];

	    math::XYZTLorentzVector tpP4 = tag->p4() + (vprobes[0].first)->p4();//Probably not needed for anything.
	    
	    evt_.elec(0).p4_=tag->p4();
	    evt_.elec(1).p4_=(vprobes[0].first)->p4();
		
		reco::SuperClusterRef SCt = tag->get<reco::SuperClusterRef>();
		reco::SuperClusterRef SCp = (vprobes[0].first)->get<reco::SuperClusterRef>();
		evt_.elec(0).detEta_ = SCt->eta();
		evt_.elec(1).detEta_ = SCp->eta();
		
	//std::cout << "Got the 4-vectors " << std::endl;
		 // If there are two probes with the tag continue
	    if( vprobes.size() > 1 ) { std::cout << " More Than 2 Probes; Z Rapidity "<< tpP4.Rapidity() << " Mass " << sqrt ( tpP4.Dot(tpP4) ) << " PT " << tpP4.Pt()<< std::endl; continue;}
	//std::cout << "Only 2 probes " << std::endl;	
		for( int itype=0; itype<(int)passProbeCandTags_.size(); ++itype ){
	   //std::cout << "Looping over the types" << std::endl;
           // Passing Probe Candidates
           edm::Handle<reco::CandidateView> passprobes;
           if ( !iEvent.getByLabel(passProbeCandTags_[itype],passprobes) ) {
	           LogWarning("ZFromData") << "Could not extract tag cands with input tag "  << passProbeCandTags_[itype];
                   std::cout << "DIDn't get the darn tag..... " << std::endl;
           } 
		
	       ///int ppass = ProbePassProbeOverlap(vprobes[0].first,passprobes);
		   //std::cout << "doing the cuts themselves " << std::endl;
		   evt_.elec(0).cutResult(cutName_[itype], ProbePassProbeOverlap(tag,passprobes,exactMatch_[itype]));
		   evt_.elec(1).cutResult(cutName_[itype], ProbePassProbeOverlap(vprobes[0].first,passprobes,exactMatch_[itype]));
		   evt_.n_elec=2;
	       // Did this tag cause a L1 and/or HLT trigger?
	       //bool l1Trigger = false;
	       //bool hltTrigger = false;
           }
         }
	  }
   ///}
  //std::cout << " Ending this one loop woohoo " << std::endl;
  /* 

          if (!oneevt){
	    if (extraHistos_ ) allCaseExtraFirst_.Fill(evt_.elec(0).p4_, evt_.elec(1).p4_, evtMC_.elec(0).p4_, evtMC_.elec(1).p4_);
             allCaseFirst_.Fill(evt_.elec(0).p4_, evt_.elec(1).p4_);
	  }
          
	  for (std::map<std::string,EffInfo*>::iterator eff=effInfo_.begin(); eff!=effInfo_.end(); eff++) {
	    evt_.elec(0).cutResult(eff->first, true); //For now I pick only "nice" electrons
	    evt_.elec(1).cutResult(eff->first,oneCut(eff->second, iEvent, probes[j], &evt_.elec(1).p4_));// Can be rearranged later
	  
    std::cout << "Returning " << numDataZ << " fellas " << std::endl;
  */
 
  allCaseFirst_.Fill(evtMC_.elec(0), evtMC_.elec(1),evtMC_.elec(0).p4_, evtMC_.elec(1).p4_); //This actually just is the RAW MC information
  if (evt_.n_elec!=2) return; // need 2 and only 2
  evt_.afterLoad(); //added to be consistent 
  std::cout << " There was 1 good pair " << std::endl;
  //
  // fill histograms for before any selection
  allCase_.Fill(evt_.elec(0), evt_.elec(1),evtMC_.elec(0).p4_, evtMC_.elec(1).p4_);
  if (extraHistos_) allCaseExtra_.Fill(evt_.elec(0), evt_.elec(1), evtMC_.elec(0), evtMC_.elec(1));

  // these stages merely _fill_ bits.  They do not apply cuts!
  stdCuts_.acceptanceCuts(evt_.elec(0));
  stdCuts_.acceptanceCuts(evt_.elec(1));
  stdCuts_.ptCuts(evt_.elec(0));
  stdCuts_.ptCuts(evt_.elec(1));
  stdCuts_.dummyCuts(evt_.elec(0));
  stdCuts_.dummyCuts(evt_.elec(1));
  //std::cout << "MYevt Start 3" << std::endl;
  //EfficiencyCut* keep=0;
  //std::cout << "Got the 4-vectors " << std::endl;
 
    

    ////// apply all cuts and store results into ZEle objects
    ////for (std::map<std::string,EfficiencyCut*>::iterator cut=theCuts_.begin(); cut!=theCuts_.end(); cut++) {
    ////  for (int ne=0; ne<2; ne++) {
	////evt_.elec(ne).cutResult(cut->first,cut->second->passesCut(evt_.elec(ne).p4_));
     //// }
    /////}
    //------------------------------------------------------------------//
    // all what is before this point goes to a separate producer
    //                           ++ decouple ++
    // all what is after this point goes to a separate analyzer
    //------------------------------------------------------------------//
    
    bool pairing;
    
      if (!quiet_) evt_.dump();

      // selections done at first pass, no matter whether there are trials or not
      for (std::map<std::string,ZShapeZDef*>::const_iterator q=zdefs_.begin(); q!=zdefs_.end(); ++q) {
	ZPlots* plots=zplots_[q->first];
        bool ptorder = zptorder_[q->first];
        int e1 = 0;
        int e2 = 1;
        if (ptorder) if ( evt_.elec(0).p4_.Pt() < evt_.elec(1).p4_.Pt() ) {e1 = 1; e2 =0;}
        
	// acceptance is always the first cut
	if (!q->second->pass(evt_,1,1,0,&pairing)) continue;
	
        if (pairing) {std::swap(e1,e2); }
        int e1n = e1;
        int e2n = e2;
	// fill standard histograms after acceptance
	plots->acceptance_.Fill(evt_.elec(e1), evt_.elec(e2),evtMC_.elec(e1).p4_, evtMC_.elec(e2).p4_);
	if (extraHistos_) plots->acceptanceExtra_.Fill(evt_.elec(e1), evt_.elec(e2), evtMC_.elec(e1), evtMC_.elec(e2));
	
	// next n-cuts
	bool ok=true;
	for (int j=1; ok && j<q->second->criteriaCount(ZShapeZDef::crit_E1); j++) {
	  ok=q->second->pass(evt_,j+1,j+1,0,&pairing);
          e1n = e1;
          e2n = e2;
          if (pairing) {std::swap(e1n,e2n); }
	  if (ok)
          {
	    plots->postCut_[j-1].Fill(evt_.elec(e1n), evt_.elec(e2n),evtMC_.elec(e1n).p4_, evtMC_.elec(e2n).p4_);
            if (extraHistos_) plots->postCutExtra_[j-1].Fill(evt_.elec(e1n), evt_.elec(e2n), evtMC_.elec(e1n), evtMC_.elec(e2n));
	  }
        }// now the Z cuts
	for (int j=0; ok && j<q->second->criteriaCount(ZShapeZDef::crit_Z); j++) {
	  ok=q->second->pass(evt_,1000,1000,j+1,&pairing);
          e1n = e1;
          e2n = e2;	  
          if (pairing) {std::swap(e1n,e2n); }
          if (ok) 
          {
            plots->zCut_[j].Fill(evt_.elec(e1n), evt_.elec(e2n),evtMC_.elec(e1n).p4_, evtMC_.elec(e2n).p4_);
            if (extraHistos_) plots->zCutExtra_[j].Fill(evt_.elec(e1n), evt_.elec(e2n),evtMC_.elec(e1n), evtMC_.elec(e2n));
          }	
        }// criteria 
        
    }


}


////THIS IS NO LONGER NEEDED, can be deleted
void ZFromData::fillMCEvent(const HepMC::GenEvent* Evt) {
  int status;
  int pid;

  ////////////////////////
  // looking for electrons
  int myPid = 11;
  evtMC_.clear();
  int ne=0;

  //std::cout << " Do I even do the MC stuff? " << std::endl;
  // looping on all the particles of the simulated event
  for(HepMC::GenEvent::particle_const_iterator mcpart = Evt->particles_begin();
      mcpart != Evt->particles_end();
      ++ mcpart ) {

    status = (*mcpart)->status();
    pid    = (*mcpart)->pdg_id();
     
    if (status == 3 && abs(pid) == myPid ){

      const HepMC::GenVertex * vertex_=(*mcpart)->production_vertex();

      if (ne==0) {
	evtMC_.vtx_=::math::XYZPoint(vertex_->position().x(),vertex_->position().y(),vertex_->position().z() );
      }
      
      math::XYZTLorentzVector  momentum;
      momentum.SetPx(  (*mcpart)->momentum().x() );        momentum.SetPy(  (*mcpart)->momentum().y() );
      momentum.SetPz(  (*mcpart)->momentum().z() );        momentum.SetE(   (*mcpart)->momentum().t() );
      
      if (ne<2) {
	evtMC_.elec(ne).p4_=momentum;
        evtMC_.elec(ne).detEta_= evtMC_.elec(ne).detectorEta(evtMC_.vtx_);
      }     
      ne++;

    }

  }
  // end loop on particles
  evtMC_.n_elec=ne;


  //
  // check than we have 2 and only 2 electrons
 
  if (ne!= 2)
    {
      edm::LogWarning("ZFromData") << " we need two electrons (pid: " << myPid << " ) while we have: " << ne;
    //  std::cout << " WHAT is going on here? " << std::endl;
      return;
    }


  //
  // reorder electrons in pT
  float pt1 = evtMC_.elec(0).p4_.Pt();
  float pt2 = evtMC_.elec(1).p4_.Pt();
  
  //std::cout << " The Pt of these are.. " << pt1 << " and " << pt2 << std::endl;

  // flipping order only if necessary
  if (pt2 > pt1){
    std::swap( evtMC_.elec(0), evtMC_.elec(1));
  }

  
}  
  





// ------------ method called once each job just before starting event loop  ------------
void 
ZFromData::beginJob(const edm::EventSetup&)
{
  // smearing the target efficiency according to statistics 
  ///if (statsBox_.trials>0) createAlternateEfficiencies();
  
  // file to write out the histograms produced by the ZFromData
  ///histoFile_   = new TFile(outFileName_.c_str(),"RECREATE");
  
  edm::Service<TFileService> fs;

  // one directory for histos of event before any selection
  TFileDirectory subDir=fs->mkdir("All");

  allCase_.Book(subDir);
  if (extraHistos_){ allCaseExtra_.Book(subDir);}

  TFileDirectory subDirF=fs->mkdir("AllFirst");
  allCaseFirst_.Book(subDirF);
  if (extraHistos_) allCaseExtraFirst_.Book(subDirF);
  
  //
  // one directory for each Z definition
  for (std::map<std::string,ZShapeZDef*>::const_iterator q=zdefs_.begin(); q!=zdefs_.end(); ++q) {
    ZPlots* zplots=new ZPlots(q->second->criteriaCount(ZShapeZDef::crit_E1)-1, q->second->criteriaCount(ZShapeZDef::crit_Z)); // -1 for acceptance

    TFileDirectory sd = fs->mkdir(q->first.c_str());

    // organize in sub-directories in order to have histograms at each step

    TFileDirectory td = sd.mkdir("Acceptance");
    zplots->acceptance_.Book(td);    
    if (extraHistos_) zplots->acceptanceExtra_.Book(td);
    
    // looping over cuts foreseen in given Z definition (for E1 and E2 simultaneously)
    for (int i=1; i<q->second->criteriaCount(ZShapeZDef::crit_E1); ++i) {
      char dirname[1024];
      std::string c1=q->second->criteria(ZShapeZDef::crit_E1,i);
      std::string c2=q->second->criteria(ZShapeZDef::crit_E2,i);
      
      // name is combination of the cut on E1 and E2, when both available
      if (c1==c2 || i>=q->second->criteriaCount(ZShapeZDef::crit_E2)) 
	sprintf(dirname,"C%02d-%s",i,c1.c_str());
      else 
	sprintf(dirname,"C%02d-%s-%s",i,c1.c_str(),c2.c_str());
      //std::cout << " BGDBG 62 " << std::endl;
      // one sub-dir for each step of selection
      td = sd.mkdir(dirname);
      zplots->postCut_[i-1].Book(td);
      
       //if (extraHistos_){td->cd(); zplots->postCutExtra_[i-1].Book();}
      if (extraHistos_){ zplots->postCutExtra_[i-1].Book(td);}
    }

    int iOffset=q->second->criteriaCount(ZShapeZDef::crit_E1);
    for (int i=0; i<q->second->criteriaCount(ZShapeZDef::crit_Z); ++i) {
      char dirname[1024];
      std::string c1=q->second->criteria(ZShapeZDef::crit_Z,i);

      sprintf(dirname,"C%02d-%s",i+iOffset,c1.c_str());

      // one sub-dir for each step of selection
      td = sd.mkdir(dirname);
      zplots->zCut_[i].Book(td);
      if (extraHistos_){ zplots->zCutExtra_[i].Book(td);}
    }

    zplots_[q->first]=zplots;
  }
 
 

}



// ------------ method called once each job just after ending the event loop  ------------
void 
ZFromData::endJob() {


}

void ZFromData::createAlternateEfficiencies() {
  return;
}


void ZFromData::createAlternateZDefs(const std::string& targetZDefSys, const std::string& targetEffSys) {
  
  return;
  
}
/*
bool ZFromData::oneCut(EffInfo *effinfo, const edm::Event &event,const EgEff::EmObjectRef &myprobe, math::PtEtaPhiMLorentzVector *e2p4){
  //Look at the collecion of the trials
  //see if that member of the collection is in the passed.
  //I should add possibly the probe to do some matching....
  //STEPS
  //Get EmObject branch from "passcollection"
  //Match probe to one EmObject from "passcollection"
  bool passEvt = false;
  edm::Handle<EgEff::EmObjectCollection> probeCollection;
  edm::Handle<EgEff::EmObjectCollection> passCollection;
   
  //Get the collections
  //std::cout << " ZFD Do I get here -3 " << std::endl;
  //std::cout << " ZFD Name of Collection " << effinfo->trialCollection << std::endl;
  try
    {
      event.getByLabel((effinfo->trialCollection), probeCollection);
    }
  catch(cms::Exception &ex)
    {
      //edm::LogError("ZShapeEventFromData") << "Error! Can't get collection " << effinfo->trialCollection;
      //throw ex;
      //std::cout << " ZFD SHIT, ain't no probes " << std::endl;
      return false;
    }
  //std::cout << " ZFD Do I get here -4 " << std::endl;
  try
    {
      event.getByLabel((effinfo->passCollection), passCollection);
    }
  catch(cms::Exception &ex)
    {
      //edm::LogError("ZShapeEventFromData") << "Error! Can't get collection " << effinfo->passCollection;
      //throw ex;
      //std::cout << " ZFD SHIT, ain't no passes " << std::endl;
      return false;
    }
   
  //Now, if there is a match in the probe collection we look in the pass collection
  // if there is a match in the pass collection, then we are good.
  //std::cout << " ZFD ProbeCollection Size "<< probeCollection->size() << std::endl;
  for(unsigned int i = 0; i < probeCollection->size(); ++i)
    {
      //std::cout << " ZFD Do I get here -1 " << std::endl;
      const EgEff::EmObjectRef probe = edm::Ref<EgEff::EmObjectCollection>(probeCollection, i);
      //std::cout << " ZFD Do I get here -2 " << std::endl;
      if ((ZShapeUtils::matchEmObjects(myprobe, probe))) for(unsigned int j = 0; j < passCollection->size(); ++j)
	{ 
	  *e2p4 = probe->p4();//Changed from myprobe->p4();
	  std::cout << " MYNEW Now " << *e2p4 << std::endl;
          std::cout << " MYNEW Pro " << myprobe->p4() << std::endl;
          std::cout << " MYNEW Probe " << probe->p4() << std::endl;
	  //std::cout << " ZFD Do I get here 0 " << std::endl;
	  const EgEff::EmObjectRef pass = edm::Ref<EgEff::EmObjectCollection>(passCollection, j);
	  //std::cout << " ZFD Do I get here 1 " << std::endl;
	  if((ZShapeUtils::matchEmObjects(myprobe, pass)))  {passEvt = true;} //To same some time I can put the return in here.   
	  //std::cout << " ZFD Do I get here 2 Passed " << passEvt << std::endl;
	  std::cout << " TYU Passed " << passEvt;
	  if (effinfo->needIso) {if (!isolated(event,pass)) {passEvt = false;}}
	  if (effinfo->needHlt) {if (!hlt(event,pass)) {passEvt = false;}}
	  std::cout << " : now " << passEvt << std::endl;
	}
    }
   

  return passEvt; //returning a passed event
}
*/
/*
bool ZFromData::isolated(const edm::Event &_iEvent, const EgEff::EmObjectRef &probe){
  return _Isopassed.CheckIso(_iEvent,probe);
}

bool ZFromData::hlt(const edm::Event &_iEvent, const EgEff::EmObjectRef &probe){
  return _HLTpassed.CheckTrigger(_iEvent,probe);
}
*/
// ***************** Does this probe pass? ******************** //

int ZFromData::ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, 
					      edm::Handle<reco::CandidateView>& passprobes, bool checkExactOverlap_ )
{
   int ppass = 0;
   if( passprobes.isValid() )
   {
      for( int ipp=0; ipp<(int)passprobes->size(); ++ipp )
      {

	bool isOverlap = MatchObjects(&((*passprobes)[ipp]),  probe,checkExactOverlap_);

	reco::SuperClusterRef probeSC;
	reco::SuperClusterRef passprobeSC; 

	if( not probe.isNull() ) probeSC  =  probe->get<reco::SuperClusterRef>();

	reco::CandidateBaseRef ref = passprobes->refAt(ipp);

	if( not ref.isNull() ) passprobeSC = ref->get<reco::SuperClusterRef>();
	
	isOverlap = isOverlap && ( probeSC == passprobeSC );


	if( isOverlap ) ppass = 1;
      } 
   }

   return ppass;
}
// ************************************************************ //


// ***************** Trigger object matching ******************** //
bool ZFromData::MatchObjects( const reco::Candidate *hltObj, 
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
//DEFINE_ANOTHER_FWK_MODULE(ZFromData);
