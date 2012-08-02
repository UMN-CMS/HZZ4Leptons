#include "TH1.h"
#include "TH2.h"
#include "HiggsEvent.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "HiggsHists.h"
    //HZZ4Leptons/AnalysisModules/src/
//HistPerDef::HistPerDef(){};
    
void HistPerDef::Book(TFileDirectory *mydir, const std::string& post, int type) 
{
    std::string t, T; // histogram title string;
    if(type != 10)
    {
	    TH1::SetDefaultSumw2();
	    t = post + "_HY";
	    T = post + " Reco H Rapidity";
	    HY = mydir->make<TH1D> (t.c_str(), T.c_str(), 20, -5, 5 );
	    t = post + "_HMass";
	    T = post + " Reco H mass";        
	    HMass = mydir->make<TH1D> (t.c_str(), T.c_str(), 60, 90, 150 );  
	    t = post + "_Z1Mass";
	    T = post + " Reco Z1 mass";
	    Z1mass = mydir->make<TH1D>(t.c_str(), T.c_str(), 70, 50, 120 );  
	    t = post + "_Z2Mass";
	    T = post + " Reco Z2 mass";
	    Z2mass = mydir->make<TH1D>(t.c_str(), T.c_str(), 90, 0, 90 );    
	    t = post + "_Z1Y";
	    T = post + " Reco Z1 Rapidity";
	    Z1Y = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, -5, 5 );  
	    t = post + "_Z2Y";
	    T = post + " Reco Z2 Rapidity";
	    Z2Y = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, -5, 5 );
	    t = post + "_l1Pt";
	    T = post + " Reco L1 Pt";
	    l1Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 150 );    
	    t = post + "_l2Pt";
	    T = post + " Reco L2 Pt";
	    l2Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 100 );
	    t = post + "_l3Pt";
	    T = post + " Reco L3Pt";
	    l3Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 25, 0, 100 );    
	    t = post + "_l4Pt";
	    T = post + " Reco L4 Pt";
	    l4Pt = mydir->make<TH1D>(t.c_str(), T.c_str(), 25, 0, 100 );
	    t = post + "_l1Eta";
	    T = post + " Reco L1 Eta";
	    l1Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );    
	    t = post + "_l2Eta";
	    T = post + " Reco L2 Eta";
	    l2Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );
	    t = post + "_l3Eta";
	    T = post + " Reco L3 Eta";
	    l3Eta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );    
	    t = post + "_l4Eta";
	    T = post + " Reco L4 Eta";
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
	    t = post + "_l1PtVsNV";
	    T = post + " Top Pt vs. nVertex";
	    l1PtVsNV = mydir->make<TH2D>(t.c_str(), T.c_str(), 50, 0, 50, 50, 0, 150 );
	    t = post + "_l1etaVsNV";
	    T = post + " First Eta vs. nVertex";
	    l1etaVsNV = mydir->make<TH2D>(t.c_str(), T.c_str(), 50, 0, 50, 20, -5, 5 );
    }
    
    if(type == 1) //All things HF
    {
    	t = post + "_Pt";
    	T = post + " Pt";
    	HFPt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 150 );
    	t = post + "_Eta";
    	T = post + " Eta";
    	HFeta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );
    	t = post + "_e9e25";
    	T = post + " E(3x3)/E(5x5)";
    	e9e25 = mydir->make<TH1D>(t.c_str(), T.c_str(), 30, 0.9, 1.05 );
    	t = post + "_var2d";
    	T = post + " Var2d";
    	var2d = mydir->make<TH1D>(t.c_str(), T.c_str(), 20, 0.2, 1.2 );    	
    	t = post + "_e9e25VsNT";
	    T = post + " E(3x30/E(5x5) vs. nVertex";
	    e9e25VsNV = mydir->make<TH2D>(t.c_str(), T.c_str(), 50, 0, 50, 30, 0.9, 1.05 );
	    t = post + "_var2dVsNT";
	    T = post + " var2d vs. nVertex";
	    var2dVsNV = mydir->make<TH2D>(t.c_str(), T.c_str(), 50, 0, 50, 20, 0.2, 1.2 );
    }
    
    if(type == 3) //All Forward
    {
    	t = post + "_Pt";
    	T = post + " Pt";
    	FwdElPt = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 150 );
    	t = post + "_Eta";
    	T = post + " Eta";
    	FwdElEta = mydir->make<TH1D>(t.c_str(), T.c_str(), 100, -5, 5 );
    }
    
    if(type==4) //All things GSF
    {
    	t = post + "_ecalSumEt_1";
    	T = post + " E1 Ecal Sum Et";
    	ecalSumEt_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 10 );
    	t = post + "_ecalSumEt_2";
    	T = post + " E2 Ecal Sum Et";
    	ecalSumEt_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 10 );
    	t = post + "_ecalSumEt_3";
    	T = post + " E3 Ecal Sum Et";
    	ecalSumEt_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 10 );
    	t = post + "_ecalSumEt_4";
    	T = post + " E4 Ecal Sum Et";
    	ecalSumEt_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 10 );
    	t = post + "_hcalIso_1";
    	T = post + " E1 Hcal Isolation";
    	hcalIso_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );
    	t = post + "_hcalIso_2";
    	T = post + " E2 Hcal Isolation";
    	hcalIso_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );
    	t = post + "_hcalIso_3";
    	T = post + " E3 Hcal Isolation";
    	hcalIso_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );
    	t = post + "_hcalIso_4";
    	T = post + " E4 Hcal Isolation";
    	hcalIso_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );    	
    	t = post + "_tkSumPt_1";
    	T = post + " E1 Tracker Sum Pt";
    	tkSumPt_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 5 );
    	t = post + "_tkSumPt_2";
    	T = post + " E2 Tracker Sum Pt";
    	tkSumPt_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 5 );
    	t = post + "_tkSumPt_3";
    	T = post + " E3 Tracker Sum Pt";
    	tkSumPt_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 5 );
    	t = post + "_tkSumPt_4";
    	T = post + " E4 Tracker Sum Pt";
    	tkSumPt_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 5 );    	
	    t = post + "_ecalIsoByGSF_1";
		T = post + " E1 Ecal Isolation";
		EcalIsoByGSF_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_ecalIsoByGSF_2";
		T = post + " E2 Ecal Isolation";
		EcalIsoByGSF_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_ecalIsoByGSF_3";
		T = post + " E3 Ecal Isolation";
		EcalIsoByGSF_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_ecalIsoByGSF_4";
		T = post + " E4 Ecal Isolation";
		EcalIsoByGSF_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_trackIso_1";
    	T = post + " E1 Tracker Isolation";
    	trackIso_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.25 );
    	t = post + "_trackIso_2";
    	T = post + " E2 Tracker Isolation";
    	trackIso_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.25 );
    	t = post + "_trackIso_3";
    	T = post + " E3 Tracker Isolation";
    	trackIso_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.25 );
    	t = post + "_trackIso_4";
    	T = post + " E4 Tracker Isolation";
    	trackIso_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.25 );		
		t = post + "_netIso_1";
		T = post + " E1 Net Isolation";
		netIso_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_netIso_2";
		T = post + " E2 Net Isolation";
		netIso_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_netIso_3";
		T = post + " E3 Net Isolation";
		netIso_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		t = post + "_netIso_4";
		T = post + " E4 Net Isolation";
		netIso_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 0.6 );
		
		t = post + "_netIso_12";
		T = post + " Net Isolation e1+e2";
		netIso_12 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		t = post + "_netIso_13";
		T = post + " Net Isolation e1+e3";
		netIso_13 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		t = post + "_netIso_14";
		T = post + " Net Isolation e1+e4";
		netIso_14 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		t = post + "_netIso_23";
		T = post + " Net Isolation e2+e3";
		netIso_23 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		t = post + "_netIso_24";
		T = post + " Net Isolation e2+e4";
		netIso_24 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		t = post + "_netIso_34";
		T = post + " Net Isolation e3+e4";
		netIso_34 = mydir->make<TH1D>(t.c_str(), T.c_str(), 60, 0, 1. );
		
		t = post + "_sIeIe_1";
		T = post + " E1 SigmaIetaIeta";
		sIeIe_1 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.1);
		t = post + "_sIeIe_2";
		T = post + " E2 SigmaIetaIeta";
		sIeIe_2 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.1);
		t = post + "_sIeIe_3";
		T = post + " E3 SigmaIetaIeta";
		sIeIe_3 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.1);
		t = post + "_sIeIe_4";
		T = post + " E4 SigmaIetaIeta";
		sIeIe_4 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.1);
		t = post + "_HadrOverEM_1";
		T = post + " E1 HardonicOverEM";
		HadrOverEM_1 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.25);
		t = post + "_HadrOverEM_2";
		T = post + " E2 HardonicOverEM";
		HadrOverEM_2 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.25);
		t = post + "_HadrOverEM_3";
		T = post + " E3 HardonicOverEM";
		HadrOverEM_3 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.25);
		t = post + "_HadrOverEM_4";
		T = post + " E4 HardonicOverEM";
		HadrOverEM_4 = mydir->make<TH1D>(t.c_str(),T.c_str(), 50, 0, 0.25);    	
    	t = post + "_E25Max_1";
    	T = post + " el-1 E 2x5 Max";
    	E25Max_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 20, 0, 1 );    	
    	t = post + "_E25Max_2";
    	T = post + " el-2 E 2x5 Max";
    	E25Max_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 20, 0, 1 );    	
    	t = post + "_E15_1";
    	T = post + " el-1 E 1x5";
    	E15_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 20, 0, 1 );    	
    	t = post + "_E15_2";
    	T = post + " el-2 E 1x5";
    	E15_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1 );   	
    	t = post + "_E55_1";
    	T = post + " el-1 E 5x5";
    	E55_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1 );    	
    	t = post + "_E55_2";
    	T = post + " el-2 E 5x5";
    	E55_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1 );
    	
    	t = post + "_HoEMBarrel";
    	T = post + " H/EM Barrel";
    	HoEMBarrel = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );
    	t = post + "_HoEMEndcap";
    	T = post + " H/EM Endcap";
    	HoEMEndcap = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 0.5 );
    	t = post + "_PfIso04LC";
    	T = post + " PfIsolation Low Pt, Central";
    	PfIso04LC = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_PfIso04LF";
    	T = post + " PfIsolation Low Pt, Forward";
    	PfIso04LF = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_PfIso04HC";
    	T = post + " PfIsolation High Pt, Central";
    	PfIso04HC = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_PfIso04HF";
    	T = post + " PfIsolation High Pt, Forward";
    	PfIso04HF = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	
    	t = post + "_PfIsoPhoton";
    	T = post + " PfIsolation Photon, e1";
    	PfIsoPhoton = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_PfIsoNHadron";
    	T = post + " PfIsolation Neutral Hadron, e1";
    	PfIsoNHadron = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_PfIsoCHadron";
    	T = post + " PfIsolation Charged Hadron, e1";
    	PfIsoCHadron = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	
    	t = post + "_dz_1";
    	T = post + " e1 |dz|";
    	dz_1 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_dz_2";
    	T = post + " e1 |dz|";
    	dz_2 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_dz_3";
    	T = post + " e1 |dz|";
    	dz_3 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    	t = post + "_dz_4";
    	T = post + " e1 |dz|";
    	dz_4 = mydir->make<TH1D>(t.c_str(), T.c_str(), 50, 0, 1. );
    }
    
    if(type == 10)
    {
	    Phi = mydir->make<TH1D>("Phi", "Phi (Angle between Z decay planes)", 50, -3.1, 3.1 );
	    Phi1 = mydir->make<TH1D>("Phi1", "Phi1 (Angle between Z1 plane and beam)", 50, -3.1, 3.1 );
	    cosTheta0 = mydir->make<TH1D>("cosTheta0","Cos(#theta_{0}) (Angle between PZ1 and beam)",50,-1,1);
	    cosTheta1 = mydir->make<TH1D>("cosTheta1","Cos(#theta_{1}) (Angle between Pl1 and PZ1)",50,-1,1);
	    cosTheta2 = mydir->make<TH1D>("cosTheta2","Cos(#theta_{2}) (Angle between Pl3 and PZ2)",50,-1,1);
    }
}// end of HistPerDef::Book()

void HistPerDef::Fill(const HiggsEvent& he,int type) 
{
	if(type != 10)
	{
	    HY->Fill(he.HY);
	    HMass->Fill(he.mH);
	    Z1mass->Fill(he.mZ1);
		Z2mass->Fill(he.mZ2);
		Z1Y->Fill(he.YZ1);
		Z2Y->Fill(he.YZ2);
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
			l1PtVsNV->Fill(he.n_primary_vertex,he.l1pt);
			l1etaVsNV->Fill(he.n_primary_vertex,he.l1eta);
		}
		else
		{	
			l1PtVsEta->Fill(he.l2eta,he.l2pt);
			l2PtVsEta->Fill(he.l1eta,he.l1pt);
			l1PtVsNV->Fill(he.n_primary_vertex,he.l2pt);
			l1etaVsNV->Fill(he.n_primary_vertex,he.l2eta);
		}
	}
   
	//Determine what type of event we got
	if(type==1)//HF
	{
		HFPt->Fill(he.l2pt);
		HFeta->Fill(he.l2eta);
		e9e25->Fill(he.e9e25);
		var2d->Fill(he.var2d);
		e9e25VsNV->Fill(he.n_primary_vertex,he.e9e25);
		var2dVsNV->Fill(he.n_primary_vertex,he.var2d);
	}
	
	if(type==3) //All Forward
	{
		FwdElPt->Fill(he.l2pt);
		FwdElEta->Fill(he.l2eta);
	}
	
	if(type==4)//GSF
	{
		ecalSumEt_1->Fill(he.ecalSumEt_1);
		ecalSumEt_2->Fill(he.ecalSumEt_2);
		ecalSumEt_3->Fill(he.ecalSumEt_3);
		ecalSumEt_4->Fill(he.ecalSumEt_4);
		hcalIso_1->Fill(he.hcalIso_1);
		hcalIso_2->Fill(he.hcalIso_2);
		hcalIso_3->Fill(he.hcalIso_3);
		hcalIso_4->Fill(he.hcalIso_4);
		tkSumPt_1->Fill(he.tkSumPt_1);
		tkSumPt_2->Fill(he.tkSumPt_2);
		tkSumPt_3->Fill(he.tkSumPt_3);
		tkSumPt_4->Fill(he.tkSumPt_4);	
		EcalIsoByGSF_1->Fill(he.ecalIsoByGSF_1);
		EcalIsoByGSF_2->Fill(he.ecalIsoByGSF_2);
		EcalIsoByGSF_3->Fill(he.ecalIsoByGSF_3);
		EcalIsoByGSF_4->Fill(he.ecalIsoByGSF_4);
		trackIso_1->Fill(he.trackIso_1);
		trackIso_2->Fill(he.trackIso_2);
		trackIso_3->Fill(he.trackIso_3);
		trackIso_4->Fill(he.trackIso_4);
		netIso_1->Fill(he.netIso_1);
		netIso_2->Fill(he.netIso_2);
		netIso_3->Fill(he.netIso_3);
		netIso_4->Fill(he.netIso_4);
		
		netIso_12->Fill(he.netIso_12);
		netIso_13->Fill(he.netIso_13);
		netIso_14->Fill(he.netIso_14);
		netIso_23->Fill(he.netIso_23);
		netIso_24->Fill(he.netIso_24);
		netIso_34->Fill(he.netIso_34);
		
		HadrOverEM_1->Fill(he.HoEM_1);
		HadrOverEM_2->Fill(he.HoEM_2);
		HadrOverEM_3->Fill(he.HoEM_3);
		HadrOverEM_4->Fill(he.HoEM_4);
		sIeIe_1->Fill(he.sIeIe_1);
		sIeIe_2->Fill(he.sIeIe_2);
		sIeIe_3->Fill(he.sIeIe_3);
		sIeIe_4->Fill(he.sIeIe_4);		
		E25Max_1->Fill(he.e25Max_1);
		E25Max_2->Fill(he.e25Max_2);		
		E15_1->Fill(he.e15_1);
		E15_2->Fill(he.e15_2);		
		E55_1->Fill(he.e55_1);
		E55_2->Fill(he.e55_2);
		
		PfIsoPhoton->Fill(he.PfIsoPhoton);
		PfIsoNHadron->Fill(he.PfIsoNHadron);
		PfIsoCHadron->Fill(he.PfIsoCHadron);
		
		if(he.l1pt<20. && he.l1eta<1.48)
    	{
    		PfIso04LC->Fill(he.PfIso04_1);
    		HoEMBarrel->Fill(he.HoEM_1);
    	}
    	else if(he.l1pt<20. && he.l1eta>1.48)
    	{
    		PfIso04LF->Fill(he.PfIso04_1);
    		HoEMEndcap->Fill(he.HoEM_1);
    	}
    	else if(he.l1pt>20. && he.l1eta<1.48)
    	{
    		PfIso04HC->Fill(he.PfIso04_1);
    		HoEMBarrel->Fill(he.HoEM_1);
    	}
    	else if(he.l1pt>20. && he.l1eta>1.48)
    	{
    		PfIso04HF->Fill(he.PfIso04_1);
    		HoEMEndcap->Fill(he.HoEM_1);
    	}
    	
    	if(he.l2pt<20. && he.l2eta<1.48)
    	{
    		PfIso04LC->Fill(he.PfIso04_2);
    		HoEMBarrel->Fill(he.HoEM_2);
    	}
    	else if(he.l2pt<20. && he.l2eta>1.48)
    	{
    		PfIso04LF->Fill(he.PfIso04_2);
    		HoEMEndcap->Fill(he.HoEM_2);
    	}
    	else if(he.l2pt>20. && he.l2eta<1.48)
    	{
    		PfIso04HC->Fill(he.PfIso04_2);
    		HoEMBarrel->Fill(he.HoEM_2);
    	}
    	else if(he.l2pt>20. && he.l2eta>1.48)
    	{
    		PfIso04HF->Fill(he.PfIso04_2);
    		HoEMEndcap->Fill(he.HoEM_2);
    	}
    	
    	if(he.l3pt<20. && he.l3eta<1.48)
    	{
    		PfIso04LC->Fill(he.PfIso04_3);
    		HoEMBarrel->Fill(he.HoEM_3);
    	}
    	else if(he.l3pt<20. && he.l3eta>1.48)
    	{
    		PfIso04LF->Fill(he.PfIso04_3);
    		HoEMEndcap->Fill(he.HoEM_3);
    	}
    	else if(he.l3pt>20. && he.l3eta<1.48)
    	{
    		PfIso04HC->Fill(he.PfIso04_3);
    		HoEMBarrel->Fill(he.HoEM_3);
    	}
    	else if(he.l3pt>20. && he.l3eta>1.48)
    	{
    		PfIso04HF->Fill(he.PfIso04_3);
    		HoEMEndcap->Fill(he.HoEM_3);
    	}
    	
    	if(he.l4pt<20. && he.l4eta<1.48)
    	{
    		PfIso04LC->Fill(he.PfIso04_4);
    		HoEMBarrel->Fill(he.HoEM_4);
    	}
    	else if(he.l4pt<20. && he.l4eta>1.48)
    	{
    		PfIso04LF->Fill(he.PfIso04_4);
    		HoEMEndcap->Fill(he.HoEM_4);
    	}
    	else if(he.l4pt>20. && he.l4eta<1.48)
    	{
    		PfIso04HC->Fill(he.PfIso04_4);
    		HoEMBarrel->Fill(he.HoEM_4);
    	}
    	else if(he.l4pt>20. && he.l4eta>1.48)
    	{
    		PfIso04HF->Fill(he.PfIso04_4);
    		HoEMEndcap->Fill(he.HoEM_4);
    	}
	}
	
	if(type == 10)
	{
		Phi->Fill(he.Phi);
		Phi1->Fill(he.Phi1);
		cosTheta0->Fill(he.cosTheta0);	
		cosTheta1->Fill(he.cosTheta1);
		cosTheta2->Fill(he.cosTheta2);
	}

}// end of HistPerDef:Fill()
















