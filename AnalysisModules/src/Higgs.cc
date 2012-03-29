// -*- C++ -*-
//
// Package:    Higgs
// Class:      Higgs
// 
/**\class Higgs Higgs.cc HZZ4Leptons/AnalysisModules/src/Higgs.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: Higgs.cc,v 1.3 2012/03/21 16:09:25 afinkel Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID 
// Order valid for 38X only.  Can be moved after Frameworkfwd.h in 39X
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
// Needed for 39X
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
// #include "DataFormats/PatCandidates/interface/GenericParticle.h"
// #include "DataFormats/PatCandidates/interface/Electron.h"
// #include "DataFormats/PatCandidates/interface/Muon.h"
// #include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TRandom.h"

#include "HZZ4Leptons/AnalysisModules/src/HiggsEvent.h"
#include "HZZ4Leptons/AnalysisModules/src/HiggsCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


template <class T>
inline std::string int2str(T i)
{
    std::ostringstream ss;
    ss << i;
    return ss.str();
}

template <class T> void outputCandidate(const T& p)
{
    std::cout << "pt=" << p.pt() << " GeV, eta=" << p.eta() << ", phi=" << p.phi();
}

inline void labelJetIDaxis(TAxis *ax)
{
    ax->SetBinLabel(1, "Neither");
    ax->SetBinLabel(2, "Loose");
    ax->SetBinLabel(3, "Tight");
}

class Higgs : public edm::EDFilter
{
public:
    explicit Higgs(const edm::ParameterSet&);
    ~Higgs();


private:

    virtual void respondToOpenInputFile(edm::FileBlock const& fb)
    {
        currentFile_ = fb.fileName();
    }

    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    virtual TH1 *bookRunHisto(uint32_t runNumber);

    edm::InputTag muonTag_;
  // edm::InputTag trackTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;
    edm::InputTag elecTag_;
  edm::InputTag elecMap_; 
    edm::InputTag photTag_;
    edm::InputTag hfTag_;
  edm::InputTag rhoTag_ ; 


    int evtCounter;

    double ZwinMinGeV_, ZwinMaxGeV_; // for trigger efficiency studies

	int elecCut_ ; 

    int pileupEra_;
    double puShift_ ;
    reweight::PoissonMeanShifter poissonNvtxShifter_ ; 
    
    std::string currentFile_;
    bool dolog_;
    bool firstEvent_;

    JetCorrectionUncertainty *jecuObj_;

    edm::LumiReWeighting MCweightByVertex_;

    std::map<uint32_t, TH1 *> m_runHistos_;

  bool doPDFreweight_;
  std::string pdfReweightBaseName, pdfReweightTargetName;
  int pdfReweightBaseId, pdfReweightTargetId;


    // ----------member data ---------------------------

    struct HistPerDef
    {
        //book histogram set w/ common suffix inside the provided TFileDirectory
        void book(TFileDirectory *, const std::string&) ; 
      // void bookTagProbe(TFileDirectory *, const std::string&);
        // fill all histos of the set with the two electron candidates
      void fill(// pat::MuonCollection muons,
                  // pat::JetCollection jets,
                  // pat::METCollection metc,
                  bool isMC,
                  double wgt,
                  bool pfJets,
		  int nPU,
		  int nPV);
        // fill all histos of the set with the two electron candidates
      void fill(const HiggsEvent& he) ; 
        // Special fill for muon efficiency studies
        // void fill(const pat::Muon& theTag, const pat::Muon& theProbe, const double probeTrkIso, const double wgt) ; 
        // void fill(const pat::Muon& theTag, const pat::GenericParticle& theProbe, const double trkIso, const double wgt) ; 
        

        TFileDirectory *mydir;

    };

    bool init_;

    // gf set of histo for all Z definitions in a stack

    struct HistStruct
    {

        TFileDirectory *rundir;
        TH1 *cutlevel, *HcandMass, *H4muMass, *H4GSFeMass, *H2mu2GSFMass, *HGSF_FEE_2muMass, *HGSF_HF_2muMass, *H2GSFe2muMass, *HGSF_HF_2GSFMass, *HGSF_FEE_2GSFMass, *PassingZ1mass, *PassingZ2mass ;

        HistPerDef noCuts;

    } hists;

    struct CutsStruct
    {
      double minimum_e1_z1_pt ;
      double minimum_e2_z1_pt ;
      double minimum_e1_z2_pt ;
      double minimum_e2_z2_pt ;
      double minimum_eNT_z1_pt ; 
      double minimum_eHF_z1_pt ; 
      double minimum_mu1_z1_pt ;
      double minimum_mu2_z1_pt ;
      double minimum_mu1_z2_pt ;
      double minimum_mu2_z2_pt ;
      double maximum_e_abseta ;
      double minimum_eNT_abseta ;
      double maximum_eNT_abseta ;
      double minimum_eHF_abseta ;
      double maximum_eHF_abseta ;
      double maximum_mu_abseta ;
      double minimum_z1_mass ;
      double maximum_z1_mass ;
      double minimum_z2_mass ;
      double maximum_z2_mass ;
      double minimum_zz_mass;
      double electron_reliso_limit ;
      double muon_reliso_limit ;
    } cuts;

};

void Higgs::HistPerDef::book(TFileDirectory *td, const std::string& post) 
{
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    mydir = td;

    std::cout << "book doing nothing at the moment" << std::endl ; 

}// end of book()

void Higgs::HistPerDef::fill(// pat::MuonCollection muons,
			     // pat::JetCollection jets,
			     // pat::METCollection metc,
                               bool isMC,
                               double wgt,
                               bool pfJets,
			       int nPU, int nPV)
{
  std::cout << "Basic fill not doing anything yet" << std::endl ; 
}

void Higgs::HistPerDef::fill(const HiggsEvent& he) 
{
  std::cout << "fill not doing anything yet" << std::endl ;


}// end of fill()


//
// constants, enums and typedefs
//

//
// static data member definitions
//



//
// constructors and destructor
//

Higgs::Higgs(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //
    dolog_ = iConfig.getParameter<bool>("DoLog");

    muonTag_ = iConfig.getParameter< edm::InputTag > ("muonTag");
    elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");
    photTag_ = iConfig.getParameter< edm::InputTag > ("photonTag");
    hfTag_   = iConfig.getParameter< edm::InputTag > ("hfTag");
    // trackTag_ = iConfig.getParameter< edm::InputTag > ("trackTag");

    elecMap_ = iConfig.getParameter< edm::InputTag > ("electronMap");
    elecCut_ = iConfig.getParameter<int> ("electronCutValue") ;  

    rhoTag_  = iConfig.getParameter< edm::InputTag > ("fjPhotonRho");

    cuts.minimum_mu1_z1_pt  = iConfig.getParameter<double>("minZ1mu1pt");
    cuts.minimum_mu2_z1_pt  = iConfig.getParameter<double>("minZ1mu2pt");
    cuts.minimum_mu2_z2_pt  = iConfig.getParameter<double>("minZ2mu2pt");
    cuts.maximum_mu_abseta  = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.minimum_e1_z1_pt   = iConfig.getParameter<double>("minZ1e1pt");
    cuts.minimum_e2_z1_pt   = iConfig.getParameter<double>("minZ1e2pt");
    cuts.minimum_e2_z2_pt   = iConfig.getParameter<double>("minZ2e2pt");
    cuts.minimum_eNT_z1_pt  = iConfig.getParameter<double>("minNoTrkpt");
    cuts.minimum_eHF_z1_pt  = iConfig.getParameter<double>("minHFpt");
    cuts.maximum_e_abseta   = iConfig.getParameter<double>("maxElecAbsEta");
    cuts.minimum_eNT_abseta = iConfig.getParameter<double>("minNoTrkAbsEta");
    cuts.maximum_eNT_abseta = iConfig.getParameter<double>("maxNoTrkAbsEta");
    cuts.minimum_eHF_abseta = iConfig.getParameter<double>("minHFElecAbsEta");
    cuts.maximum_eHF_abseta = iConfig.getParameter<double>("maxHFElecAbsEta");
    cuts.minimum_z1_mass    = iConfig.getParameter<double>("minZ1Mass");
    cuts.maximum_z1_mass    = iConfig.getParameter<double>("maxZ1Mass");
    cuts.minimum_z2_mass    = iConfig.getParameter<double>("minZ2Mass");
    cuts.maximum_z2_mass    = iConfig.getParameter<double>("maxZ2Mass");
    cuts.minimum_zz_mass    = iConfig.getParameter<double>("min4objMass");
    cuts.electron_reliso_limit = iConfig.getParameter<double>("electronRelIsoLimit");
    cuts.muon_reliso_limit     = iConfig.getParameter<double>("muonRelIsoLimit");

    pileupEra_ = iConfig.getParameter<int>("pileupEra");
    puShift_   = iConfig.getParameter<double>("systPileupShift") ;
    if ( fabs(puShift_) > 0.001 ) poissonNvtxShifter_ = reweight::PoissonMeanShifter( puShift_ ) ; 

    // ==================== Init other members ====================
    //

    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;

    hists.cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 11, -1.5, 9.5);
    hists.cutlevel->GetXaxis()->SetBinLabel(1, "Raw");
    hists.cutlevel->GetXaxis()->SetBinLabel(2, "No cuts");
    hists.cutlevel->GetXaxis()->SetBinLabel(3, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(4, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(5, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(6, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(7, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(8, "Dummy");
    hists.cutlevel->GetXaxis()->SetBinLabel(9, "Dummy");
    hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));
    
    hists.H4muMass = fs->make<TH1D > ("H4muMass", "H to 4 mu ", 50, 100, 150);
    hists.H4GSFeMass = fs->make<TH1D > ("H4GSFeMass", "H to 4 GSF el ", 50, 100, 150);
    hists.H2mu2GSFMass = fs->make<TH1D > ("H2mu2GSFMass", "H to 2 mu with 2 GSF el", 50, 100, 150);
    hists.HGSF_FEE_2muMass = fs->make<TH1D > ("HGSF_FEE_2muMass", "H to GSF+FEE and 2mu", 50, 100, 150);
    hists.HGSF_HF_2muMass = fs->make<TH1D > ("HGSF_HF_2muMass", "H to GSF+HF and 2 mu", 50, 100, 150);
    hists.HcandMass = fs->make<TH1D > ("HcandMass", "Potentially Reconstructable H Masses", 50, 100, 150);
    hists.H2GSFe2muMass = fs->make<TH1D > ("H2GSFe2muMass", "H to 2GSF el and Z*->mumu ", 50, 100, 150);
    hists.HGSF_HF_2GSFMass = fs->make<TH1D > ("HGSF_HF_2GSFMass", "H to GSF+HF and 2GSF", 50, 100, 150);
    hists.HGSF_FEE_2GSFMass = fs->make<TH1D > ("HGSF_FEE_2GSFMass", "H to 2GSF el and HF el ", 50, 100, 150);
    hists.PassingZ1mass = fs->make<TH1D> ("PassingZ1mass", "Passing Z masses", 60, 0., 120);
    hists.PassingZ2mass = fs->make<TH1D> ("PassingZ2mass", "Passing Z* masses", 60, 0., 120);

    init_ = false;

    MCweightByVertex_ = edm::LumiReWeighting(higgs::generate_flat10_mc(pileupEra_),
            higgs::get_standard_pileup_data(pileupEra_));
    
    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "minMu1pt (z1)     = " << cuts.minimum_mu1_z1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt (z1)     = " << cuts.minimum_mu2_z1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt (z2)     = " << cuts.minimum_mu2_z2_pt << " GeV" << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "electronTag       = " << elecTag_ << std::endl;
    std::cout << "minE1pt (z1)      = " << cuts.minimum_e1_z1_pt << " GeV" << std::endl;
    std::cout << "minE2pt (z1)      = " << cuts.minimum_e2_z1_pt << " GeV" << std::endl;
    std::cout << "minE2pt (z2)      = " << cuts.minimum_e2_z2_pt << " GeV" << std::endl;
    std::cout << "maxElecAbsEta     = " << cuts.maximum_e_abseta << std::endl;
    std::cout << "electronTag (NT)  = " << photTag_ << std::endl;
    std::cout << "minNTpt (z1)      = " << cuts.minimum_eNT_z1_pt << " GeV" << std::endl;
    std::cout << "NT Elec AbsEta    = " << cuts.minimum_eNT_abseta << " to " << cuts.maximum_eNT_abseta << std::endl;
    std::cout << "electronTag (HF)  = " << hfTag_ << std::endl;
    std::cout << "HF Elec AbsEta    = " << cuts.minimum_eHF_abseta << " to " << cuts.maximum_eHF_abseta << std::endl;
    std::cout << "minHFpt (z1)      = " << cuts.minimum_eHF_z1_pt << " GeV" << std::endl;
    std::cout << "Z1 window         : " << cuts.minimum_z1_mass << " to " << cuts.maximum_z1_mass << " GeV" << std::endl ; 
    std::cout << "Z2 window         : " << cuts.minimum_z2_mass << " to " << cuts.maximum_z2_mass << " GeV" << std::endl ; 
    std::cout << "electronRelIso    = " << cuts.electron_reliso_limit << std::endl;
    std::cout << "muonRelIso        = " << cuts.muon_reliso_limit << std::endl;
    std::cout << "min4objMass       = " << cuts.minimum_zz_mass << " GeV" << std::endl;

    std::cout << "PU era (shift)    = " << pileupEra_ << " (" << puShift_ << ")" << std::endl;


}

Higgs::~Higgs()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}



//
// member functions
//


TH1 * Higgs::bookRunHisto(uint32_t runNumber)
{
    std::string runstr = int2str<uint32_t > (runNumber);
    return hists.rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
}


// ------------ method called to for each event  ------------

bool Higgs::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HiggsEvent higgsEvent;

    evtCounter++;

    higgsEvent.isMC = !iEvent.isRealData();
    higgsEvent.scaleMuE();

    if(iEvent.isRealData())
    {
        uint32_t runn = iEvent.id().run();
        std::map<uint32_t, TH1 *>::const_iterator it = m_runHistos_.find(runn);
        TH1 *runh;
        if(it == m_runHistos_.end())
        {
            runh = bookRunHisto(runn);
            m_runHistos_[runn] = runh;
        }
        else
            runh = it->second;
        runh->Fill(1);
    }

    edm::Handle<reco::MuonCollection> recoMuons;
    iEvent.getByLabel(muonTag_, recoMuons);

    edm::Handle<reco::GsfElectronCollection> recoElectrons;
    iEvent.getByLabel(elecTag_, recoElectrons);

    edm::Handle< edm::ValueMap<float> > eIDValueMap;
    iEvent.getByLabel(elecMap_, eIDValueMap);

    edm::Handle<reco::PhotonCollection> recoGammas;
    iEvent.getByLabel(photTag_, recoGammas);

    edm::Handle<double> fastjetPhotonRho; 
    iEvent.getByLabel(rhoTag_, fastjetPhotonRho) ; 
    double photonRho = (*(fastjetPhotonRho.product())) ; 

    edm::Handle<reco::RecoEcalCandidateCollection> recoHFElectrons;
    iEvent.getByLabel(hfTag_, recoHFElectrons);

    edm::Handle<reco::HFEMClusterShapeAssociationCollection> clusterAssociation;
    iEvent.getByLabel("hfEMClusters",clusterAssociation);

    if(higgsEvent.isMC)
    {
        edm::Handle<std::vector<PileupSummaryInfo> > pPU;
        iEvent.getByLabel("addPileupInfo", pPU);
        std::pair<float, double> pileup = higgs::pileupReweighting(pPU, MCweightByVertex_);
        higgsEvent.n_pue = pileup.first ; // Will only be used for studies, thus no syst. correction necessary
        higgsEvent.eventWgt *= pileup.second;
        if ( fabs(puShift_) > 0.001 ) higgsEvent.eventWgt *= poissonNvtxShifter_.ShiftWeight( pileup.first ) ; 
	// PDF reweighting
	// if (doPDFreweight_) {
	//   edm::Handle<GenEventInfoProduct> geip;
	//   iEvent.getByLabel("generator",geip);
      
	//   float Q=geip->pdf()->scalePDF;
	//   int id1=geip->pdf()->id.first;
	//   int id2=geip->pdf()->id.second;
	//   float x1=geip->pdf()->x.first;
	//   float x2=geip->pdf()->x.second;

	//   higgsEvent.eventWgt *= higgs::getPDFWeight(Q,id1,x1,id2,x2,
        //                                          doPDFreweight_,
        //                                          pdfReweightBaseId,pdfReweightTargetId);	  
	// }
        
        // generator information
        edm::Handle<reco::GenParticleCollection> genInfo;
        if(iEvent.getByLabel("genParticles", genInfo))
        {
            higgsEvent.decayID(*genInfo);
        }
    }
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    higgsEvent.n_primary_vertex = higgs::numberOfPrimaryVertices(pvHandle);

    if (!recoElectrons.isValid() || !recoMuons.isValid() || !recoGammas.isValid() || !recoHFElectrons.isValid()) 
    {
        /*std::cout << "Exiting as valid PAT objects not found" << std::endl;
        std::cout << "Electrons:    " << recoElectrons.isValid() << std::endl;
        std::cout << "HF Electrons: " << recoElectrons.isValid() << std::endl;
        std::cout << "Photons:      " << recoGammas.isValid() << std::endl;
        std::cout << "Muons:        " << recoMuons.isValid() << std::endl;*/
        return false;
    }

    if(firstEvent_)
    {
        //std::cout << "===============================================" << std::endl;
        firstEvent_ = false;
    }

    hists.cutlevel->Fill(-1.0, higgsEvent.eventWgt);


    
    // Look for valid muons
    std::vector<reco::Muon> muCands =
      higgs::getMuonList(recoMuons, cuts.minimum_mu2_z2_pt, cuts.maximum_mu_abseta, false) ; 
      
    std::vector<reco::GsfElectron> eleCands =
      higgs::getElectronList(recoElectrons, eIDValueMap, cuts.minimum_e2_z2_pt, cuts.maximum_e_abseta, elecCut_) ; 

    std::vector<reco::RecoEcalCandidate> hfEleCands =
      higgs::getElectronList(recoHFElectrons, clusterAssociation, cuts.minimum_e2_z1_pt, cuts.maximum_eHF_abseta) ;

    std::vector<reco::Photon> ntEleCands =
      higgs::getElectronList(recoGammas, photonRho, cuts.minimum_e2_z1_pt, 
			     cuts.minimum_eNT_abseta, cuts.maximum_eNT_abseta) ;

    /*std::cout << "Event has " << muCands.size() << " muons, " 
	      << eleCands.size() << "(gsf) + " 
	      << ntEleCands.size() << "(nt) + " 
	      << hfEleCands.size() << "(hf) = " 
	      << eleCands.size() + ntEleCands.size() + hfEleCands.size() 
	      << " electrons" << std::endl ; */

    higgsEvent.muCands  = muCands ; 
    higgsEvent.gsfCands = eleCands ; 
    higgsEvent.ntCands  = ntEleCands ; 
    higgsEvent.hfCands  = hfEleCands ; 

    if( !higgsEvent.getZ1(cuts.minimum_e1_z1_pt,cuts.minimum_e2_z1_pt,
			   cuts.minimum_mu1_z1_pt,cuts.minimum_mu2_z1_pt,cuts.minimum_z1_mass) ) return false ; 
    if ( !higgsEvent.getZ2(cuts.minimum_e2_z2_pt,cuts.minimum_mu2_z2_pt,
			   cuts.minimum_z2_mass,cuts.minimum_zz_mass) ) return false ;
			   
	higgsEvent.calculate();
    // Basic selection requirements: Require at least two muons, two jets

	int totalEMcands = recoElectrons->size() + recoGammas->size() + recoHFElectrons->size() ; 
    if ( (recoMuons->size() + totalEMcands) > 3 ) 
    {
        hists.cutlevel->Fill(0.0, higgsEvent.eventWgt);
        hists.HcandMass->Fill(higgsEvent.mH);
        hists.PassingZ1mass->Fill(higgsEvent.mZ1);
		hists.PassingZ2mass->Fill(higgsEvent.mZ2);
    }

    else return false;    
    
	//Determine what type of event we got
	
	if( higgsEvent.getZ2(cuts.minimum_e2_z2_pt,cuts.minimum_mu2_z2_pt,
			   cuts.minimum_z2_mass,cuts.minimum_zz_mass)==1 )
	{
		switch( higgsEvent.getZ1(cuts.minimum_e1_z1_pt,cuts.minimum_e2_z1_pt,
			   cuts.minimum_mu1_z1_pt,cuts.minimum_mu2_z1_pt,cuts.minimum_z1_mass) )
		{
			case 1:
				std::cout<<"Case 1 evoked: 4mu"<<std::endl;
				hists.H4muMass->Fill(higgsEvent.mH);
				break;
			case 2:
				std::cout<<"Case 2 evoked: 2GSFel and 2mu"<<std::endl;
				hists.H2GSFe2muMass->Fill(higgsEvent.mH);
				break;
			case 3:
				std::cout<<"Case 3 evoked: GSF + FarEE and 2mu"<<std::endl;
				hists.HGSF_FEE_2muMass->Fill(higgsEvent.mH);
				break;
			case 4:
				std::cout<<"Case 4 evoked: GSF + HF and 2mu"<<std::endl;
				hists.HGSF_HF_2muMass->Fill(higgsEvent.mH);
				break;
		}
	}	
	
	if( higgsEvent.getZ2(cuts.minimum_e2_z2_pt,cuts.minimum_mu2_z2_pt,
			   cuts.minimum_z2_mass,cuts.minimum_zz_mass)==2 )
	{
		switch( higgsEvent.getZ1(cuts.minimum_e1_z1_pt,cuts.minimum_e2_z1_pt,
			   cuts.minimum_mu1_z1_pt,cuts.minimum_mu2_z1_pt,cuts.minimum_z1_mass) )
		{
			case 1:
				std::cout<<"Case 5 evoked: 2mu + 2GSF"<<std::endl;
				hists.H2mu2GSFMass->Fill(higgsEvent.mH);
				break;
			case 2:
				std::cout<<"Case 6 evoked: 4GSF el"<<std::endl;
				hists.H4GSFeMass->Fill(higgsEvent.mH);
				break;
			case 3:
				std::cout<<"Case 7 evoked: GSF + FarEE and 2GSF"<<std::endl;
				hists.HGSF_FEE_2GSFMass->Fill(higgsEvent.mH);
				break;
			case 4:
				std::cout<<"Case 8 evoked: GSF + HF and 2GSF"<<std::endl;
				hists.HGSF_HF_2GSFMass->Fill(higgsEvent.mH);
				break;
		}
	}		

    return true;
}


// ------------ method called once each job just before starting event loop  ------------

void Higgs::beginJob()
{
    firstEvent_ = true;
    evtCounter = 0;

    // if (doPDFreweight_) {
    //   higgs::initPDFSet(1,pdfReweightBaseName);
    //   higgs::initPDFSet(2,pdfReweightTargetName);
    // }

}

// ------------ method called once each job just after ending the event loop  ------------

void Higgs::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(Higgs);



