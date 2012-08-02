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
// $Id: Higgs.cc,v 1.12 2012/06/26 23:43:23 afinkel Exp $
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
#include "TLorentzVector.h"

#include "HZZ4Leptons/AnalysisModules/src/HiggsEvent.h"
#include "HZZ4Leptons/AnalysisModules/src/HiggsCommon.h"
#include "HZZ4Leptons/AnalysisModules/src/HiggsHists.h"

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
    
    TLorentzVector TLV(reco::Particle::LorentzVector V);
    void Angular(HiggsEvent& he);

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
    
    int Z1type, Z2type;
    int e1Cut, e2Cut ;
    bool passIso, passIsoHoEM, passIsoHoEMsIe, passElecsCut_1, passElecsCut_2, passElecsCut_3 ;

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
    
	struct CutLevels //this contains the folders for the sub-channels; one for every cut level
    {   
    	public:    
    	HistPerDef  H4mu, H4GSFe, H2mu2GSF, HGSF_FEE_2mu, HGSF_HF_2mu, HGSF_HF_2GSF, HGSF_FEE_2GSF,	H2GSFe2mu,
    				AllHF, AllFEE, AllFwd, AllEvents,
    				AngularAll, AngularExMu, AngularGSFxMu;
    				
    	TH1 *dummy; //*Phi, *Phi1, *cosTheta0, *cosTheta1, *cosTheta2;

    	void book(TFileDirectory *mydir, const std::string);
    	bool fill(const HiggsEvent& );

    } ;

    bool init_;
    
    TFileDirectory *rundir;
    
    CutLevels NoCuts, ElecsCut_1, ElecsCut_2, ElecsCut_3, IsoCut, IsoHoEMCut, IsoHoEMsIeCut;
    
    TH1 *genPU, *recoPU, *cutlevel;
    
    // gf set of histo for all Z definitions in a stack

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

void Higgs::CutLevels::book(TFileDirectory *mydir, const std::string name)
{
	/*edm::Service<TFileService> fs;
	TFileDirectory *mydir = new TFileDirectory(fs->mkdir(name.c_str()));
	mydir->cd();*/
	
	AllHF.Book(new TFileDirectory(mydir->mkdir("AllHF")),"AllHF",1);
	AllFEE.Book(new TFileDirectory(mydir->mkdir("AllFEE")),"AllFEE",4);
	AllFwd.Book(new TFileDirectory(mydir->mkdir("AllFwd")),"AllFwd",3);
    AllEvents.Book(new TFileDirectory(mydir->mkdir("AllEvents")),"AllEvents");
    H4mu.Book(new TFileDirectory(mydir->mkdir("4mu")),"4Mu");
    H4GSFe.Book(new TFileDirectory(mydir->mkdir("4GSF")),"4GSF",4);
    H2mu2GSF.Book(new TFileDirectory(mydir->mkdir("2Mu_2GSF")),"2Mu_2GSF",4);
    HGSF_FEE_2mu.Book(new TFileDirectory(mydir->mkdir("GSF_FarEE_2Mu")),"GSF_FarEE_2Mu",4);
    HGSF_FEE_2GSF.Book(new TFileDirectory(mydir->mkdir("GSF_FarEE_2GSF")),"GSF_FarEE_2GSF",4);
    HGSF_HF_2mu.Book(new TFileDirectory(mydir->mkdir("GSF_HF_2Mu")),"GSF_HF_2Mu",4);
    HGSF_HF_2GSF.Book(new TFileDirectory(mydir->mkdir("GSF_HF_2GSF")),"GSF_HF_2GSF",4);
    H2GSFe2mu.Book(new TFileDirectory(mydir->mkdir("2GSF_2Mu")),"2GSF_2Mu",4);
    AngularAll.Book(new TFileDirectory(mydir->mkdir("AngularAll")),"AngularAll",10); //Warning: 10 is for the Angular kind only!  
    AngularExMu.Book(new TFileDirectory(mydir->mkdir("AngularExMu")),"AngularExMu",10);
    AngularGSFxMu.Book(new TFileDirectory(mydir->mkdir("AngularGSFxMu")),"AngularGSFxMu",10);
    
    dummy = mydir->make<TH1D>("dummy","Dummy Hist",10,1,0);
}

bool Higgs::CutLevels::fill(const HiggsEvent& he)
{
	//if (!he.Z1flavor || !he.Z2flavor) return false;
	AllEvents.Fill(he);
	AngularAll.Fill(he,10);
    if( he.Z2flavor==1 )
	{
		switch( he.Z1flavor )
		{
			case 1:
				H4mu.Fill(he);
				return true;
			case 2:
				H2GSFe2mu.Fill(he,4);
				AngularGSFxMu.Fill(he,10);
				AngularExMu.Fill(he,10);
				return true;
			case 3:
				HGSF_FEE_2mu.Fill(he,4);
				AllFEE.Fill(he,4);
				AllFwd.Fill(he,3);
				AngularExMu.Fill(he,10);
				return true;
			case 4:
				HGSF_HF_2mu.Fill(he,4);
				AngularExMu.Fill(he,10);
				AllHF.Fill(he,1);
				AllFwd.Fill(he,3);
				return true;
		} 
	}	
	if( he.Z2flavor==2 )
	{
		switch( he.Z1flavor )
		{
			case 1:
				H2mu2GSF.Fill(he,4);
				AngularGSFxMu.Fill(he,10);
				AngularExMu.Fill(he,10);
				return true;
			case 2:
				H4GSFe.Fill(he,4);
				return true;
			case 3:
				HGSF_FEE_2GSF.Fill(he,4);
				AllFEE.Fill(he,4);
				AllFwd.Fill(he,3);
				return true;
			case 4:
				HGSF_HF_2GSF.Fill(he,4);
				AllHF.Fill(he,1);
				AllFwd.Fill(he,3);
				return true;
		}
	}
	return false;
}

void Higgs::Angular(HiggsEvent& he)
{
	TLorentzVector Z1, Z2, l11, l12, l21, l22, H, Prot(0,0,1,1);
	Z1 = TLV(he.vZ1);
	l11 = TLV(he.vl1);
	l12 = TLV(he.vl2);
	Z2 = TLV(he.vZ2);
	l21 = TLV(he.vl3);
	l22 = TLV(he.vl4);
	H = TLV(he.vH);
	
	//find theta-star and phi1 from Z1 in Higgs frame
	TLorentzVector Z1clone = Z1;
	TLorentzVector l11clone = l11;
	TLorentzVector Hclone = H;
	Prot.Boost(-H.BoostVector());
	l11clone.Boost(-H.BoostVector());
	Z1clone.Boost(-H.BoostVector());
	he.cosTheta0 = Z1clone.Vect()*Prot.Vect()/Z1clone.P()/Prot.P();
	
	TVector3 NZ1 = Z1clone.Vect().Cross(l11clone.Vect());
	TVector3 NPZ = Prot.Vect().Cross(Z1clone.Vect());
	
	he.Phi1 = acos( NZ1*NPZ/( NZ1.Mag()*NPZ.Mag() ) );
	he.Phi1 *= (NZ1.Phi()>NPZ.Phi())? 1: -1;
	
	//go to Z1 frame to get theta1
	TLorentzVector Z2clone = Z2;
	l11clone = l11;
	Hclone = H;
	Z2clone.Boost(-Z1.BoostVector());
	Hclone.Boost(-Z1.BoostVector());
	l11clone.Boost(-Z1.BoostVector());
	
	he.cosTheta1 = l11clone.Vect()*(-Z2clone.Vect())/l11clone.P()/Z2clone.P();
	
	//go to Z2 frame and get theta2
	Z1clone = Z1;
	Hclone = H;
	TLorentzVector l21clone = l21;
	Z1clone.Boost(-Z2.BoostVector());
	l21clone.Boost(-Z2.BoostVector());
	
	he.cosTheta2 = l21clone.Vect()*(-Z1clone.Vect())/l21clone.P()/Z1clone.P();
	
	//now try to get the phi angle (between Z decay planes)
	l11clone=l11;
	l21clone=l21;
	Z1clone=Z1;
	Z2clone=Z2;
	l11clone.Boost(-H.BoostVector());
	l21clone.Boost(-H.BoostVector());
	Z1clone.Boost(-H.BoostVector());
	Z2clone.Boost(-H.BoostVector());
	
	TVector3 N1 = Z1clone.Vect().Cross(l11clone.Vect());
	
	TVector3 N2 = Z2clone.Vect().Cross(l21clone.Vect());
	
	he.Phi = acos( N1*N2/( N1.Mag()*N2.Mag() ) );
	he.Phi *= (N1.Phi()>N2.Phi())? 1: -1;
}

TLorentzVector Higgs::TLV(reco::Particle::LorentzVector V)
{
//converts a reco::Particle::LorentzVector into a TLorentzVector
	return TLorentzVector(V.px(),V.py(),V.pz(),V.E());
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
    Z1type = 0;
    Z2type = 0;
    e1Cut = 0;
    e2Cut = 0;
    passIso = passIsoHoEM = passIsoHoEMsIe = passElecsCut_1 = passElecsCut_2 = passElecsCut_3 = false ;
    

    // ==================== Book the histos ====================
    //
    
	edm::Service<TFileService> fs;
    
    NoCuts.book(new TFileDirectory(fs->mkdir("NoCuts")),"NoCuts");
    ElecsCut_1.book(new TFileDirectory(fs->mkdir("ElecsCut_1")),"ElecsCut_1");
    ElecsCut_2.book(new TFileDirectory(fs->mkdir("ElecsCut_2")),"ElecsCut_2");
    ElecsCut_3.book(new TFileDirectory(fs->mkdir("ElecsCut_3")),"ElecsCut_3");
    IsoCut.book(new TFileDirectory(fs->mkdir("IsoCut")),"IsoCut");
    IsoHoEMCut.book(new TFileDirectory(fs->mkdir("IsoHoEMCut")),"IsoHoEMCut");
    IsoHoEMsIeCut.book(new TFileDirectory(fs->mkdir("IsoHoEMsIeCut")),"IsoHoEMsIeCut");
    
    cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 11, -1.5, 9.5);
    cutlevel->GetXaxis()->SetBinLabel(1, "Raw");
    cutlevel->GetXaxis()->SetBinLabel(2, "No cuts");
    cutlevel->GetXaxis()->SetBinLabel(3, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(4, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(5, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(6, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(7, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(8, "Dummy");
    cutlevel->GetXaxis()->SetBinLabel(9, "Dummy");
    rundir = new TFileDirectory(fs->mkdir("RunDir"));
    
    genPU = fs->make<TH1D>("genPU", "Gen-Level Pile-Up", 50, 0, 50);
    recoPU = fs->make<TH1D>("recoPU", "Reco-Level Pile-Up", 50, 0, 50);
    
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
    return rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
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

    if (!recoElectrons.isValid() || !recoMuons.isValid() || !recoGammas.isValid() || !recoHFElectrons.isValid()) return false;

    if(firstEvent_)
    {
        //std::cout << "===============================================" << std::endl;
        firstEvent_ = false;
    }
    
    std::vector< std::pair<double,unsigned int> > eleCutsByPt;

    // Look for valid muons
    std::vector<reco::Muon> muCands =
      higgs::getMuonList(recoMuons, cuts.minimum_mu2_z2_pt, cuts.maximum_mu_abseta, pvHandle, false) ; 
      
    std::vector<reco::GsfElectron> eleCands =
      higgs::getElectronList(recoElectrons, eIDValueMap, cuts.minimum_e2_z2_pt, cuts.maximum_e_abseta, elecCut_, eleCutsByPt, pvHandle) ; 

    std::vector<reco::RecoEcalCandidate> hfEleCands =
      higgs::getElectronList(recoHFElectrons, clusterAssociation, cuts.minimum_e2_z1_pt, cuts.maximum_eHF_abseta, higgsEvent) ;

    std::vector<reco::Photon> ntEleCands =
      higgs::getElectronList(recoGammas, photonRho, cuts.minimum_e2_z1_pt, 
		     cuts.minimum_eNT_abseta, cuts.maximum_eNT_abseta) ;

    higgsEvent.SetMuonList(muCands) ; 
    higgsEvent.SetGsfElectronList(eleCands) ; 
    higgsEvent.SetPhotonList(ntEleCands) ; 
    higgsEvent.SetHFList(hfEleCands) ; 
    
    genPU->Fill(higgsEvent.n_pue);
    
    Z1type = higgsEvent.getZ1(cuts.minimum_e1_z1_pt,cuts.minimum_e2_z1_pt, cuts.minimum_mu1_z1_pt,
    							cuts.minimum_mu2_z1_pt,cuts.minimum_z1_mass,cuts.minimum_eHF_z1_pt, pvHandle);
	Z2type = higgsEvent.getZ2(cuts.minimum_e2_z2_pt,cuts.minimum_mu2_z2_pt,
			   cuts.minimum_z2_mass,cuts.minimum_zz_mass, pvHandle); 

    if( !Z1type ) return false ; 
    if ( !Z2type ) return false ;
    if (Z1type>4 || Z2type>4) return false;
    
    recoPU->Fill(higgsEvent.n_pue);
			   
	higgsEvent.calculate();
	Angular(higgsEvent);
	
	if ( !higgsEvent.passCombinedIso ) return false;
	
	for(std::vector< std::pair<double,unsigned int> >::const_iterator i = eleCutsByPt.begin(); i != eleCutsByPt.end(); i++)
    {
    	if( (fabs(higgsEvent.l1pt - i->first) < 1e-9))
    	{
    		e1Cut = i->second;
    	}
    	if( (fabs(higgsEvent.l2pt - i->first) < 1e-9))
    	{
    		e2Cut = i->second;
    	}    	
    }
	
	passIso = ( (higgsEvent.netIso_1<0.15)&&(higgsEvent.netIso_2<0.15)&&(higgsEvent.netIso_3<0.15)&&(higgsEvent.netIso_4<0.15) );
	passIsoHoEM = ( (passIso)&&(higgsEvent.HoEM_1<0.15)&&(higgsEvent.HoEM_2<0.15)&&(higgsEvent.HoEM_3<0.15)&&(higgsEvent.HoEM_4<0.15) );
	passIsoHoEMsIe = ( (passIsoHoEM)&&(higgsEvent.sIeIe_1<0.03)&&(higgsEvent.sIeIe_2<0.03)&&(higgsEvent.sIeIe_3<0.03)&&(higgsEvent.sIeIe_4<0.03) );
	passElecsCut_1 = ( ( (e1Cut & 1) == 1 )&&( (e2Cut & 1) == 1 ) );
	passElecsCut_2 = ( ( (e1Cut & 2) == 2 )&&( (e2Cut & 2) == 2 ) );
	passElecsCut_3 = ( ( (e1Cut & 3) == 3 )&&( (e2Cut & 3) == 3 ) );
	
    // Basic selection requirements: Require at four leptopns
	int totalEMcands = recoElectrons->size() + recoGammas->size() + recoHFElectrons->size() ; 
    if ( (recoMuons->size() + totalEMcands) < 4 ) return false; 
    //In following original paper and recent constraints, limit the H mass range:
    if( higgsEvent.mH<100 || higgsEvent.mH>140 ) return false;

    cutlevel->Fill(0.0, higgsEvent.eventWgt);
	NoCuts.fill(higgsEvent);
	if(passIso) IsoCut.fill(higgsEvent); 	
	if(passIsoHoEM) IsoHoEMCut.fill(higgsEvent);
	if(passIsoHoEMsIe) IsoHoEMsIeCut.fill(higgsEvent);
	if(passElecsCut_1)  ElecsCut_1.fill(higgsEvent);
	if(passElecsCut_2)  ElecsCut_2.fill(higgsEvent);
	if(passElecsCut_3)  ElecsCut_3.fill(higgsEvent);

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



