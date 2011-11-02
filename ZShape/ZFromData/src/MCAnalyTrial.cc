#include "TProfile.h"
#include "ZShape/Base/interface/EfficiencyStatistics.h"
#include "ZShape/ZFromData/interface/MCAnalyTrial.h"
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
MCAnalyTrial::MCAnalyTrial(const edm::ParameterSet& iConfig) :
  m_srcTag(iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("source"))),
  quiet_(iConfig.getUntrackedParameter<bool>("quiet",false))
{
  // 
  //now do what ever initialization is needed

  MCTag_ = iConfig.getUntrackedParameter< edm::InputTag>("MCTag");
  SCTag_ = iConfig.getUntrackedParameter< edm::InputTag>("SCTag");

  outFileName_              = iConfig.getUntrackedParameter<std::string>("outHistogramsFile", "");
  writeHistoConservatively_ = iConfig.getUntrackedParameter<bool>("writeHistoBeforeEndJob", false);
  delRMatchingCut_          = iConfig.getUntrackedParameter<double>("dRMatchCut",0.2);
  delPtRelMatchingCut_      = iConfig.getUntrackedParameter<double>("dPtMatchCut",15.0);


  scEtCut_      = iConfig.getUntrackedParameter<double>("scEtCut",5.0);

  textBinsFile_      = iConfig.getUntrackedParameter< std::string >("EffBinsFile","EffBinsFile.txt");
  effBinsFromTxt_    = new EffTableLoader(textBinsFile_);

  cbBinsFile_   = iConfig.getUntrackedParameter< std::string >("CBBinsFile","CBBinsFile.txt");
  fromfile_     = iConfig.getUntrackedParameter<bool>("FromFile", false);
  fromfilefile_ = iConfig.getUntrackedParameter< std::string >("InputFile","histoMC.root");
  
}


MCAnalyTrial::~MCAnalyTrial()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void MCAnalyTrial::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace cms;
  using namespace reco;
 
  if (fromfile_) return;

  Handle<CandidateView> mctags;
  iEvent.getByLabel(MCTag_,mctags);

  Handle<CandidateView> sctags;
  iEvent.getByLabel(SCTag_,sctags);
  std::cout << " the side of the MC tags is " << mctags->size() << std::endl;
  CandidateView::const_iterator f = mctags->begin();
  for( ; f != mctags->end(); ++f )
    {
      unsigned int index = f - mctags->begin();
      CandidateBaseRef mctagRef = mctags->refAt(index);
      if( mctagRef.isNull() ) continue;
      bool good = false;
      good = (bool) ProbePassProbeOverlap(mctagRef,sctags,false);
      int myindex = effBinsFromTxt_->GetBandIndex(mctagRef->pt(), mctagRef->eta());
      if ( good ) Hgood_->Fill(myindex);
      Hall_->Fill(myindex);
    }


}




// ------------ method called once each job just before starting event loop  ------------
void 
MCAnalyTrial::beginJob(const edm::EventSetup&)
{
  // smearing the target efficiency according to statistics 
  ///if (statsBox_.trials>0) createAlternateEfficiencies();
  
  // file to write out the histograms produced by the MCAnalyTrial
  ///histoFile_   = new TFile(outFileName_.c_str(),"RECREATE");
  if (!fromfile_)
  {
     edm::Service<TFileService> fs;

     Hgood_  = fs->make<TH1F>("Good","good", 100, 0, 100);
     Hall_   = fs->make<TH1F>("All","All", 100, 0, 100);
     Heff_   = fs->make<TH1F>("eff","eff", 100, 0, 100);
     HMCptvsSCpt_ = fs->make<TH2F>("MCptVsSCpt","MCptVsSCpt", 50, 0, 20, 50, 0, 20);
     for (int i = 0;  i < effBinsFromTxt_->size(); ++i)
     {
         SCtoMCERatio_[i] = fs->make<TH1F>(Form("zMCtoSC_%d",i),Form("zMCtoSC_%d",i), 100, 0.5, 1.3);
     }
     for ( int i =  effBinsFromTxt_->size(); i < 40; ++i)
     {
         SCtoMCERatio_[i] = fs->make<TH1F>(Form("zzMCtoSC_%d",i),Form("zzMCtoSC_%d",i), 100, 0.5, 1.3);
     }
  }

  
}




// ------------ method called once each job just after ending the event loop  ------------
void 
MCAnalyTrial::endJob() {

  if (!fromfile_) Heff_->Divide(Hgood_,Hall_,1.,1.);
  else
  {
     myFile_ = new TFile(fromfilefile_.c_str());
     //myFile_->cd("ZFromData");
     Hgood_ = (TH1F*) myFile_->Get("ZFromData/Good");
     Hall_ = (TH1F*) myFile_->Get("ZFromData/All");
     Heff_ = (TH1F*) myFile_->Get("ZFromData/eff");
     HMCptvsSCpt_ = (TH2F*) myFile_->Get("ZFromData/MCptVsSCpt");
     for (int i = 0;  i < effBinsFromTxt_->size(); ++i)
     {
         SCtoMCERatio_[i] =(TH1F*) myFile_->Get(Form("ZFromData/zMCtoSC_%d",i));
     }
  }

  std::cout << " I made it to the endjobie " << std::endl;
  //Now, I need to loop over the single histos to output the parameters to a file.
  //Maybe I should add the option to read the histograms from a file, that way
  //  i can play with the fitting settings in a real time type of enviroment.
  ofstream the1DEffFile;
  the1DEffFile.open (cbBinsFile_.c_str(),std::ios::out);
  std::cout << " Starting to write the file " << std::endl;
  the1DEffFile << "#put your comments after a '#'.\n";
  the1DEffFile << "#dimension: 1 \n\n\n";
  the1DEffFile << " " << std::setw(10) << "#PtMin"
               << " " << std::setw(10) << "PtMax"
               << " " << std::setw(10)  << "EtaMin"
               << " " << std::setw(10)  <<  "EtaMax"
               << " " << std::setw(10)  << "NumParms"
               << " " << std::setw(10)  << "alpha"
               << " " << std::setw(10) << "n"
               << " " << std::setw(10) << "Mean"
               << " " << std::setw(10) << "sigma"
               << " " << std::setw(10) << "N"
               << "\n";
  if (the1DEffFile.is_open())
    {
      TF1 *crystal = new TF1("crystal",CrystalBall,0.5,1.3,5); 
      crystal->SetParNames("#alpha","n","Mean","#sigma","N"); 
      std::cout << "Initting the fitting function " << std::endl;
      for (int r=0; r<effBinsFromTxt_->size(); ++r)
        {
          std::vector<std::pair<float, float> > EffInfo = effBinsFromTxt_->GetCellInfo(r);
          std::cout << "Setting the parameters " << std::endl;
          std::cout << " The Max is " << SCtoMCERatio_[r]->GetMaximum() << std::endl;
          //crystal->SetParameters(1.01,3.0,1.0,0.03,SCtoMCERatio_[r]->GetYaxis()->GetXmax());
          crystal->SetParameters(1.01,2.0,1.0,0.03,SCtoMCERatio_[r]->GetMaximum());
          if (r > 0 ) crystal->SetParLimits(1,1.,10.);
          std::cout << "Fitting " << std::endl;
          SCtoMCERatio_[r]->Fit(crystal,"r");
          double npars[5];
          crystal->GetParameters( &npars[0] );
          std::cout << " Done Fitting and starting to write  " << r << std::endl;

          the1DEffFile << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[0].first
                       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[0].second
                       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[1].first
                       << " " << std::setw(10) << std::setprecision(4) <<  EffInfo[1].second
                       << " " << std::setw(10) << std::setprecision(4) <<  5
                       << " " << std::setw(10) << std::setprecision(4) <<  npars[0]
                       << " " << std::setw(10) << std::setprecision(4) <<  npars[1]
                       << " " << std::setw(10) << std::setprecision(4) <<  npars[2]
                       << " " << std::setw(10) << std::setprecision(4) <<  npars[3]
                       << " " << std::setw(10) << std::setprecision(4) <<  npars[4]
                       << "\n";
        }
    }

  the1DEffFile.close();


}

Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;

    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}


int MCAnalyTrial::ProbePassProbeOverlap( const reco::CandidateBaseRef& probe, 
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
bool MCAnalyTrial::MatchObjects( const reco::Candidate *hltObj, 
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
//DEFINE_ANOTHER_FWK_MODULE(MCAnalyTrial);
