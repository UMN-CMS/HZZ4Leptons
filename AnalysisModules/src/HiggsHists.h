#include "TH1.h"
#include "TH2.h"
#include "HiggsEvent.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


class HistPerDef
    {
    	public:
        //book histogram set w/ common suffix inside the provided TFileDirectory
        void Book(TFileDirectory *, const std::string&, int type=0) ; 
    	// void bookTagProbe(TFileDirectory *, const std::string&);
        // fill all histos of the set with the two electron candidates
    	void Fill(// pat::MuonCollection muons,
                  // pat::JetCollection jets,
                  // pat::METCollection metc,
                  bool isMC,
                  double wgt,
                  bool pfJets,
		  		  int nPU,
		  		  int nPV);
        // fill all histos of the set with the two electron candidates
      	void Fill(const HiggsEvent& he,int type=0) ; 
        // Special fill for muon efficiency studies
        // void fill(const pat::Muon& theTag, const pat::Muon& theProbe, const double probeTrkIso, const double wgt) ; 
        // void fill(const pat::Muon& theTag, const pat::GenericParticle& theProbe, const double trkIso, const double wgt) ; 
        
      //private:
      //TFileDirectory *mydir;
        
      TH1  *HMass, *HY, *Z1mass, *Z2mass, *Z1Y, *Z2Y, *nVertex,
		   *l1Pt, *l1Eta, *l2Pt, *l2Eta, *l3Pt, *l3Eta, *l4Pt, *l4Eta,
      	   *ecalSumEt_1, *ecalSumEt_2, *ecalSumEt_3, *ecalSumEt_4,
      	   *hcalIso_1, *hcalIso_2, *hcalIso_3, *hcalIso_4,
      	   *tkSumPt_1, *tkSumPt_2, *tkSumPt_3, *tkSumPt_4,
      	   *EcalIsoByGSF_1, *EcalIsoByGSF_2, *EcalIsoByGSF_3, *EcalIsoByGSF_4,
      	   *trackIso_1, *trackIso_2, *trackIso_3, *trackIso_4,
      	   *netIso_1, *netIso_2, *netIso_3, *netIso_4,
      	   *netIso_12, *netIso_13, *netIso_14, *netIso_23, *netIso_24, *netIso_34,
      	   *E25Max_1, *E25Max_2, //*E25Max_3, *E25Max_4,
      	   *E15_1, *E15_2, //*E15_3, *E15_4,
      	   *E55_1, *E55_2, //*E55_3, *E55_4;
      	   *HadrOverEM_1, *HadrOverEM_2, *HadrOverEM_3, *HadrOverEM_4,
      	   *sIeIe_1, *sIeIe_2, *sIeIe_3, *sIeIe_4,
      	   *FarEEPt, *FarEEeta, 
      	   *HFPt, *HFeta, *e9e25, *var2d,
      	   *FwdElPt, *FwdElEta,
      	   *HoEMBarrel, *HoEMEndcap,
      	   *PfIso04LC, *PfIso04LF, *PfIso04HC, *PfIso04HF,
      	   *PfIsoPhoton, *PfIsoNHadron, *PfIsoCHadron,
      	   *dz_1, *dz_2, *dz_3, *dz_4,
      	   *Phi, *Phi1, *cosTheta0, *cosTheta1, *cosTheta2;
			//*H4muMass, *H4GSFeMass, *H2mu2GSFMass, *HGSF_FEE_2muMass, *HGSF_HF_2muMass, *H2GSFe2muMass, *HGSF_HF_2GSFMass, *HGSF_FEE_2GSFMass,
      
      TH2 *AllPtVsEta, *l1PtVsEta, *l2PtVsEta,
      		*l1PtVsNV, *l1etaVsNV,
      		*e9e25VsNV, *var2dVsNV;
    } ;
