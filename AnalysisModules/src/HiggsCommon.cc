#include "HiggsCommon.h"
#include "LHAPDF/LHAPDF.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeFwd.h"

const float pileup2010A = 1.2;
const float pileup2010B = 2.2;
const float pileup2011A = 5.0;

// Data corrections provided by N. Kypreos.
// MC corrections are inverse of data, and given below
const double muScaleLUTarray100[5] = { 0.379012,0.162107,0.144514,0.0125131,-0.392431 } ; 

namespace higgs {
  bool isVBTFloose(const reco::Muon& m)
  {
    return ( muon::isGoodMuon(m,muon::AllGlobalMuons) && m.innerTrack()->numberOfValidHits() > 10 ) ;  
  }
  
  // Note: No impact parameter requirement is applied!!!
  bool isVBTFtight(const reco::Muon& m)
  {
    if( !isVBTFloose(m) ) return false; // this should already have been checked.
    return ( muon::isGoodMuon(m,muon::GlobalMuonPromptTight) &&
	     m.isTrackerMuon() && 
	     m.numberOfMatches() > 1 && 
	     m.globalTrack()->hitPattern().numberOfValidPixelHits()>0 ) ; 
  }


  double muIsolation(const pat::Muon& m, const double pTscale) {
    double mupt = pTscale * m.pt() ; 
    return (m.trackIso() / mupt) ; 
  }

  void initPDFSet(int i, std::string name) {
      LHAPDF::initPDFSet(i,name) ; 
  }
    
  double getPDFWeight(float Q, int id1, float x1, int id2, float x2,
                      bool doPDFreweight, int pdfReweightBaseId, int pdfReweightTargetId) {
  
      if (!doPDFreweight) return 1.0;
  
      LHAPDF::usePDFMember(1,pdfReweightBaseId);
      double pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
      double pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  
      LHAPDF::usePDFMember(2,pdfReweightTargetId);
      double newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
      double newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
      
      double w=(newpdf1/pdf1*newpdf2/pdf2);
  
      //  printf("My weight is %f\n",w);
  
      return w;
  }
    
  int numberOfPrimaryVertices(edm::Handle<reco::VertexCollection> pvHandle) { 
    int nvertex = 0 ; 
    
    const reco::VertexCollection& vertices = *pvHandle.product();
    static const int minNDOF = 4;
    static const double maxAbsZ = 15.0;
    static const double maxd0 = 2.0;

    for (reco::VertexCollection::const_iterator vit=vertices.begin(); 
	 vit!=vertices.end(); vit++) {
      if (vit->ndof() > minNDOF && 
	  (fabs(vit->z()) <= maxAbsZ) && 
	  (fabs(vit->position().rho()) <= maxd0)) nvertex++;
    }
    return nvertex ; 
  }

 
  std::vector<float> generate_flat10_mc(int era){
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    // const double npu_probs[25] = {0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 0-4
    // 				  0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 5-9
    // 				  0.0698146584,0.0630151648,0.0526654164,0.0402754482,0.0292988928,  // 10-14
    // 				  0.0194384503,0.0122016783,0.0072070420,0.0040036370,0.0020278322,  // 15-19
    // 				  0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 // 20-24 
    // };

    // see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities for PU_S4 samples
    // Using the "spike at zero + smearing distribution" as shown on the twiki and recommended for in-time PU corrections
    const double npu_probs_summer11[35] = { 1.45346E-01,6.42802E-02,6.95255E-02,6.96747E-02,6.92955E-02, // 0-4
					    6.84997E-02,6.69528E-02,6.45515E-02,6.09865E-02,5.63323E-02, // 5-9
					    5.07322E-02,4.44681E-02,3.79205E-02,3.15131E-02,2.54220E-02, // 10-14
					    2.00184E-02,1.53776E-02,1.15387E-02,8.47608E-03,6.08715E-03, // 15-19
					    4.28255E-03,2.97185E-03,2.01918E-03,1.34490E-03,8.81587E-04, // 20-24
					    5.69954E-04,3.61493E-04,2.28692E-04,1.40791E-04,8.44606E-05, // 25-29
					    5.10204E-05,3.07802E-05,1.81401E-05,1.00201E-05,5.80004E-06  // 30-34
    }; 

    const double npu_probs_fall11[50] = { 0.003388501, 0.010357558, 0.024724258, 0.042348605, 0.058279812, // 0-4
					  0.068851751, 0.072914824, 0.071579609, 0.066811668, 0.060672356, // 5-9
					  0.054528356,  0.04919354, 0.044886042, 0.041341896,   0.0384679, // 10-14
					  0.035871463,  0.03341952, 0.030915649, 0.028395374, 0.025798107, // 15-19
					  0.023237445, 0.020602754,   0.0180688, 0.015559693, 0.013211063, // 20-24
					  0.010964293, 0.008920993, 0.007080504, 0.005499239, 0.004187022, // 25-29
					  0.003096474, 0.002237361, 0.001566428, 0.001074149, 0.000721755, // 30-34
					  0.000470838,  0.00030268, 0.000184665, 0.000112883, 6.74043E-05, // 35-39
					  3.82178E-05, 2.22847E-05, 1.20933E-05, 6.96173E-06,  3.4689E-06, // 40-44
					  1.96172E-06, 8.49283E-07, 5.02393E-07, 2.15311E-07, 9.56938E-08  // 45-49
    };

    bool isFall11 = ( era == 20113 || era == 20114 ) ; 
    const double* npu_probs = npu_probs_summer11 ;
    int npt = sizeof(npu_probs_summer11)/sizeof(double);
    if ( isFall11 ) {
        npt = sizeof(npu_probs_fall11)/sizeof(double);
        npu_probs = npu_probs_fall11 ;
    }

    std::vector<float> retval;
    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(npu_probs[i]);
    return retval;
    // double s = 0.0;
    // for(int npu=0; npu<25; ++npu){
    //   double npu_estimated = dataDist[npu];
    //   result[npu] = npu_estimated / npu_probs[npu];
    //   s += npu_estimated;
    // }
    // // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    // for(int npu=0; npu<25; ++npu){
    //   result[npu] /= s;
    // }
    // return result;
  }


  std::vector<float> get_standard_pileup_data(int era) {
    const double default_pd[] = { 100, 100, 100, 0, 0, 0, 0, 0, 0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
    const double may10_json[]= {
      3920760.81, 6081805.28, 13810357.99, 22505758.94, 28864043.84, 30917427.86, 28721324.56, 23746403.90, 17803439.77, 12274902.61, 
      7868110.47, 4729915.40, 2686011.14, 1449831.56, 747892.03, 370496.38, 177039.19, 81929.35, 36852.78, 16164.45, 
      6932.97, 2914.39, 1202.92, 488.15, 194.93, 76.64, 29.66, 11.30, 4.24, 1.56, 
      0.57, 0.20, 0.07, 0.02, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00 };
    const double dec22_json[]= {
      5525914.08, 9064207.39, 8760233.30, 6230731.01, 3615158.71, 1806273.52, 802707.46, 323868.48, 120245.31, 41458.73, 
      13360.27, 4043.62, 1153.89, 311.48, 79.76, 19.43, 4.51, 1.00, 0.21, 0.04, 
      0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00 };

    // Pileup histograms assembled from inputs in this directory: 
    // /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp

    // Note: there is a 36th entry with 0.008, but other ingredients have only 35 entries so omitting
    const double json_2011a[] = {
      1.29654E07, 5.58514E07, 1.29329E08, 2.12134E08, 2.76138E08,
      3.03604E08, 2.93258E08, 2.55633E08,  2.0497E08, 1.53264E08,
      1.07936E08, 7.21006E07,  4.5913E07,   2.797E07, 1.63426E07,
      9.17598E06, 4.95861E06, 2.58239E06,  1.2977E06,     629975,
          295784,     134470,    59260.1,    25343.9,    10530.1,
         4255.05,    1673.95,    641.776,    240.022,    87.6504,
          31.281,    10.9195,    3.73146,    1.24923,   0.602368
    } ; 

    const double json_2011b[] = {
          481142, 3.21393E06, 1.15733E07, 2.91676E07, 5.76072E07,
      9.51074E07, 1.36849E08,  1.7665E08,  2.0885E08, 2.29582E08,
      2.37228E08, 2.32243E08, 2.16642E08, 1.93361E08,  1.6564E08,
      1.36514E08, 1.08455E08, 8.31965E07, 6.17147E07, 4.43296E07,
      3.08733E07, 2.08734E07, 1.37166E07, 8.77106E06, 5.46389E06,
      3.31952E06, 1.96896E06,  1.1414E06,     647299,     359460,
          195642,     104449,    54741.4,    28184.3,    28004.9
    } ; 

    // Pileup histograms for Fall11 assembled from inputs in this directory: 
    // /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp
    // Per instructions, the so-called "true" distributions from data are used
    const double json_2011a_44x[] = { 
             0.0,   252573.0, 5.73861E06, 5.05641E07, 2.68197E08,
      5.12977E08, 5.21466E08, 3.75661E08, 2.41467E08,  1.4315E08,
      5.38318E07, 1.18125E07, 1.29073E06,     111569,    6537.93,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0
    } ; 

    const double json_2011b_44x[] = { 
             0.0,    27267.4,      35590,    74493.3,     574590,
      2.90648E06, 3.36311E07, 9.36661E07, 1.38283E08, 1.87624E08,
      2.15647E08,  2.1173E08, 1.87002E08, 1.46693E08, 9.44372E07,
      4.60317E07, 1.69231E07, 5.18161E06, 1.42805E06,     437008,
          102694,     6516.2,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0
    } ; 


    const double* pileupDist=default_pd;
    int npt = 35;
    if (era == 20111) 
    {
        npt = sizeof(json_2011a)/sizeof(double);
        pileupDist=json_2011a;
    }
    if (era == 20112) 
    {
        npt = sizeof(json_2011b)/sizeof(double);
        pileupDist=json_2011b;
    }
    if (era == 20113) 
    {
        npt = sizeof(json_2011a_44x)/sizeof(double);
        pileupDist=json_2011a_44x;
    }
    if (era == 20114) 
    {
        npt = sizeof(json_2011b_44x)/sizeof(double);
        pileupDist=json_2011b_44x;
    }
    if (era == 20110) 
    {
        npt = sizeof(may10_json)/sizeof(double);
        pileupDist=may10_json;
    }
    if (era>=20100 && era<=20109) 
    {
        npt = sizeof(dec22_json)/sizeof(double);
        pileupDist=dec22_json;
    }

    std::vector<float> retval;

    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(pileupDist[i]);
    return retval;
  }

  std::pair<float,double> pileupReweighting(const edm::Handle<std::vector<PileupSummaryInfo> >& pPU,
					    edm::LumiReWeighting& puWeight) { 

    int   nPileup = -1 ; 
    double weight = 1.0 ; 
    // float  avg_nvtx = -1. ;

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for (PVI = pPU->begin(); PVI != pPU->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();

      if (BX == 0) { 
        nPileup = PVI->getPU_NumInteractions();
        continue;
      }
    }
    weight = puWeight.weight( nPileup );
    

    // std::pair<float,double> pileupInfo = std::make_pair(avg_nvtx,weight) ; 
    std::pair<float,double> pileupInfo = std::make_pair(nPileup,weight) ; 
    return pileupInfo ; 
  }

  float jecTotalUncertainty(float jpt, float jeta,
			    JetCorrectionUncertainty *jecUnc,
			    int correctEra,
			    bool isBjet,
			    bool directionIsUp) {
    float offunc; // the "official" eta/pt dependent uncertainty

    jecUnc->setJetPt((float)jpt);
    jecUnc->setJetEta((float)jeta);
    offunc = jecUnc->getUncertainty(directionIsUp);
        // these are no longer necessary, commented out to ensure they never get applied without express intent
        //    float pileupCorrection = (2 * 0.75 * 0.8) / jpt;
        //    if      (correctEra == 1) pileupCorrection *= pileup2010A;
        //    else if (correctEra == 2) pileupCorrection *= pileup2010B;
        //    else if (correctEra == 3) pileupCorrection *= pileup2011A;
        //    else { // Merge 2010A and 2010B corrections
        //      float pileup2010 = ((3.18 * pileup2010A) + (32.96 * pileup2010B)) / 36.14;
        //      pileupCorrection *= pileup2010;
        //    }
    // Calculations taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC, version 23
    float totalUnc = offunc;//sqrt((offunc*offunc) + (pileupCorrection*pileupCorrection) + (0.025*0.025));
    return totalUnc;
  }


  // double muScaleLUT(pat::Muon& iM) { 

  //   const double etastep   = 2.0 * 2.4 / 5 ; 
  //   const double etavec[6] = {-2.4,(-2.4+1.0*etastep),(-2.4+2.0*etastep),
  // 			      (-2.4+3.0*etastep),(-2.4+4.0*etastep),2.4} ; 
  //   unsigned int ieta = 0 ; 
  //   while (ieta < 5 && iM.eta() > etavec[ieta+1]) ieta++ ;  

  //   // Alignment corrections are given in chg/TeV.
  //   // To get GeV scale, need to divide by 10^3
  //   return muScaleLUTarray100[ieta] / 1000. ; 
  // } 

  std::vector<reco::Muon> getMuonList(edm::Handle<reco::MuonCollection>& recoMuons,
				     double minPt, double maxAbsEta, 
				     bool trackerPt) {

    std::vector<reco::Muon> muonList ; 
    for (unsigned int iMuon = 0; iMuon < recoMuons->size(); iMuon++) {
        reco::Muon iM = recoMuons->at(iMuon) ; 
        if ( fabs(iM.eta()) > maxAbsEta ) continue ; 
        if ( iM.pt() < minPt ) continue ; 
        if ( !isVBTFloose(iM) ) continue ; 
        
        muonList.push_back(iM) ; 
    }
    
    std::sort(muonList.begin(),muonList.end(),pTcompare()) ; 
    /*if ( muonList.size() ) { 
      std::cout << "Sorted muon list: " << std::endl ; 
      for (unsigned int i=0; i<muonList.size(); i++) { 
	std::cout << "Muon " << i+1 << " of " << muonList.size() 
		  << " with pT " << muonList.at(i).pt() << " GeV" << std::endl ; 
      }
    }*/
    return muonList ; 
  }

  std::vector< reco::GsfElectron > getElectronList(edm::Handle<reco::GsfElectronCollection>& recoElecs,
						   							edm::Handle< edm::ValueMap<float> >& valueMap, 
                                                    double minEt, double maxAbsEta, 
						   							int cutlevel, std::vector< std::pair<double,unsigned int> >& eleCutsByPt) { 

    const edm::ValueMap<float> &eIDmap = *valueMap ;
    
    std::pair<double,unsigned int> thePair;
    
    std::vector< reco::GsfElectron > electronList ; 
    for (unsigned int iElectron=0; iElectron < recoElecs->size(); iElectron++) {
        edm::Ref<reco::GsfElectronCollection> eRef(recoElecs,iElectron);
        if ( eRef->pt() < minEt ) continue ;
        if ( fabs(eRef->eta()) > maxAbsEta ) continue ;
        thePair.first = eRef->pt();
        thePair.second = int(eIDmap[eRef]);
        eleCutsByPt.push_back(thePair);
	//std::cout << "Found an electron with pT: " << eRef->pt() << " and eta " << eRef->eta() << std::endl ; 
	double scTheta = (2*atan(exp(-eRef->superCluster()->eta()))) ;
	double e25Max = eRef->e2x5Max() ; 
	double e15 = eRef->e1x5() ; 
	double e55 = eRef->e5x5() ; 
	double ecalIso = eRef->dr03EcalRecHitSumEt() ; 

	//if ( (int(eIDmap[eRef]) & cutlevel) != cutlevel ) continue ; // electron fails ID/iso/IP/conv
        electronList.push_back( *eRef ) ;        
    }
    std::sort(electronList.begin(),electronList.end(),pTcompare()) ;
    std::sort(eleCutsByPt.begin(),eleCutsByPt.end()); 
    /*if ( electronList.size() ) { 
      std::cout << "Sorted GSF electron list: " << std::endl ; 
      for (unsigned int i=0; i<electronList.size(); i++) { 
	std::cout << "Electron " << i+1 << " of " << electronList.size() 
		  << " with pT " << electronList.at(i).pt() << " GeV" << std::endl ; 
      }
    }*/
    return electronList ; 
  }

  // A very basic HF electron ID
  bool passesHFElectronID(const reco::RecoEcalCandidate& electron, 
			  const edm::Handle<reco::HFEMClusterShapeAssociationCollection>& clusterAssociation, HiggsEvent& HE) { 

    reco::SuperClusterRef hfclusRef = electron.superCluster() ;
    const reco::HFEMClusterShapeRef hfclusShapeRef = (*clusterAssociation).find(hfclusRef)->val ;
    const reco::HFEMClusterShape&   hfshape = *hfclusShapeRef ;

    double e9e25      = hfshape.eLong3x3()/hfshape.eLong5x5();
    double var2d      = hfshape.eCOREe9()-(hfshape.eSeL()*9./8.);
    double eCOREe9    = hfshape.eCOREe9();
    double eSeL       = hfshape.eSeL();
    
    HE.var2d = var2d;
    HE.e9e25 = e9e25;

    // Parameters: e9e25_loose, e9e25_tight,  var2d_loose, var2d_tight,  eCOREe9_loose, eCOREe9_tight,  eSeL_loose, eSeL_tight;
    // hFselParams =  cms.vdouble(0.90, 0.94,      0.2, 0.40,    -9999, -9999,     9999, 9999),

    //if ( e9e25 <= 0.94) return false ; 
    //if ( var2d <= 0.40 ) return false ; 
    
    //std::cout << "HF candidate passes selection" << std::endl ; 

    return true ; 
  }

  std::vector< reco::RecoEcalCandidate > getElectronList(edm::Handle<reco::RecoEcalCandidateCollection>& recoElecs,
                                                         edm::Handle<reco::HFEMClusterShapeAssociationCollection>& clusterAssociation, 
							 double minEt,double maxAbsEta, HiggsEvent& HE) {
      
    std::vector< reco::RecoEcalCandidate > electronList ; 
    for (unsigned int iElectron=0; iElectron < recoElecs->size(); iElectron++) {
        edm::Ref<reco::RecoEcalCandidateCollection> eRef(recoElecs,iElectron);
        if ( eRef->pt() < minEt ) continue ;
        if ( fabs(eRef->eta()) > maxAbsEta ) continue ;
	if ( passesHFElectronID(*eRef,clusterAssociation, HE) ) electronList.push_back( *eRef ) ;        
    }
    std::sort(electronList.begin(),electronList.end(),pTcompare()) ; 
    return electronList ; 
  }

  // At the moment, a very rough guess for "No Track" electron ID
  // Starting point: https://twiki.cern.ch/twiki/bin/view/CMS/Vgamma2011PhotonID
  bool passesNoTrackID( const reco::Photon& electron, const double rho ) {

    //if ( electron.hadronicOverEm() >= 0.05 ) return false ; 
    //if ( electron.sigmaIetaIeta() >= 0.03 ) return false ; 

    // Values to compute for isolation 
    double scEt = electron.pt() ; 
    double ecalIsoThreshold = 4.2 +  (0.006 * scEt) + (0.090 * rho) ; 
    double hcalIsoThreshold = 2.2 + (0.0025 * scEt) + (0.180 * rho) ; 

    if ( electron.ecalRecHitSumEtConeDR04() >= ecalIsoThreshold ) return false ; 
    if ( electron.hcalTowerSumEtConeDR04() >= hcalIsoThreshold ) return false ; 
    //std::cout<<"H/EM is "<<electron.hadronicOverEm()<<"; sigmaIetaIeta is "<<electron.sigmaIetaIeta()<<std::endl;

    return true ; 
  }

  std::vector< reco::Photon > getElectronList(edm::Handle<reco::PhotonCollection>& recoElecs,
					      const double rho, double minEt, 
					      double minAbsEta, double maxAbsEta) {
      
    //if ( recoElecs->size() ) std::cout << "*** Pileup rho is: " << rho << std::endl ; 
    std::vector< reco::Photon > electronList ; 
    for (unsigned int iElectron=0; iElectron < recoElecs->size(); iElectron++) {
        edm::Ref<reco::PhotonCollection> eRef(recoElecs,iElectron);
        if ( eRef->pt() < minEt ) continue ;
        if ( fabs(eRef->eta()) < minAbsEta ) continue ;
        if ( fabs(eRef->eta()) > maxAbsEta ) continue ;
        if ( passesNoTrackID( *eRef,rho ) ) electronList.push_back( *eRef ) ;        
    }
    std::sort(electronList.begin(),electronList.end(),pTcompare()) ; 
    return electronList ; 
  }
    
}
