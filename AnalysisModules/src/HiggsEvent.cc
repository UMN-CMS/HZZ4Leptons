#include "HZZ4Leptons/AnalysisModules/src/HiggsEvent.h"
#include "TVector3.h"

HiggsEvent::HiggsEvent() { 

  n_pue = -1;
  eventWgt = 1.0 ; 
  cutlevel = -1 ; 
  nMuons = 0 ; 
  nElectrons = 0 ;

  Z1flavor = 0 ; 
  Z2flavor = 0 ;
  
  ecalIsoByGSF_1 = -1;
  ecalIsoByGSF_2 = -1;
  ecalIsoByGSF_3 = -1;
  ecalIsoByGSF_4 = -1;
  ecalSumEt_1 = -1;
  ecalSumEt_2 = -1;
  ecalSumEt_3 = -1;
  ecalSumEt_4 = -1;
  HoEM_1 = -1;
  HoEM_2 = -1;
  HoEM_3 = -1;
  HoEM_4 = -1;
  sIeIe_1 = -1;
  sIeIe_2 = -1;
  sIeIe_3 = -1;
  sIeIe_4 = -1;
  
  e25Max_1 = -1;
  e15_1 = -1;
  e55_1 = -1;
  e25Max_2 = -1;
  e15_2 = -1;
  e55_2 = -1;
  
  netIso_1 = -1;
  netIso_2 = -1;
  netIso_3 = -1;
  netIso_4 = -1;

  // Protection: scale factors set to 1 by default
  MuScale = 1.0 ;
  ElecScale = 1.0 ; 
}

void HiggsEvent::regularize() {
    mu[0] = mu1 ; mu[1] = mu2 ; mu[2] = mu3 ; mu[3] = mu4 ; 
    e[0] = e1 ; e[1] = e2 ; e[2] = e3 ; e[3] = e4 ;
}

void HiggsEvent::scaleMuE(double mufactor, double efactor) {
  MuScale   = mufactor ; 
  ElecScale = efactor ; 
}

// This routine determines the best Z candidate, among other things...
int HiggsEvent::getZ1(double minElePt1, double minElePt2, double minMuPt1, double minMuPt2, double minMass, double minHFpt) { 

  const double mZpdg = 91.1876 ; 
  double mZdiff = mZpdg ; 
  double deltaM = 0;
  reco::Particle::LorentzVector zCand, l1Cand, l2Cand;
  //std::cout << "Looking for the primary Z" << std::endl ; 

  if ( muCands.size() ) { //Muon channel
    for (unsigned int i=0; i<muCands.size()-1; i++) { 
      if (muCands.at(i).pt() < minMuPt1) break ; // List is pT ordered
      for (unsigned int j=i+1; j<muCands.size(); j++) { 
		if (muCands.at(j).pt() < minMuPt2) break ; // List is pT ordered
		if ((muCands.at(i).charge()*muCands.at(j).charge()) > 0) continue ;
		 
		zCand = muCands.at(i).p4() + muCands.at(j).p4() ;
		if (zCand.M() < minMass) continue ; 
		deltaM = fabs(zCand.M() - mZpdg) ; 
		if (deltaM < mZdiff) { // New "best" Z candidate
		  mZdiff = deltaM ; 
		  vZ1 = zCand ; 
		  Z1flavor = 1 ; 
		  Z1idx.first = i ; Z1idx.second = j ;
		  l1Cand = muCands.at(i).p4();
		  l2Cand = muCands.at(j).p4(); 
		}
      }
    }
  }

  for (unsigned int i=0; i<gsfCands.size(); i++) { //tracked Electron channel
    if (gsfCands.at(i).pt() < minElePt1) break ; // List is pT ordered
    // Second electron is a GSF electron
    for (unsigned int j=i+1; j<gsfCands.size(); j++) { 
      if (gsfCands.at(j).pt() < minElePt2) break ; // List is pT ordered
      if ((gsfCands.at(i).charge()*gsfCands.at(j).charge()) > 0) continue ;
       
      zCand = gsfCands.at(i).p4() + gsfCands.at(j).p4() ; 
      if (zCand.M() < minMass) continue ; 
      deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
		mZdiff = deltaM ; 
		vZ1 = zCand ; 
		Z1flavor = 2 ; 
		Z1idx.first = i ; Z1idx.second = j ; 
		vl1 =  gsfCands.at(i).p4();
		ecalSumEt_1 = gsfCands.at(i).dr03EcalRecHitSumEt();
		tkSumPt_1 = gsfCands.at(i).dr03TkSumPt();
		hcalSumEt_1 = gsfCands.at(i).dr03HcalTowerSumEt();
		ecalIsoByGSF_1 = ecalSumEt_1/vl1.Et();
		trackIso_1 = tkSumPt_1/vl1.Pt();
		e25Max_1 = gsfCands.at(i).e2x5Max()/gsfCands.at(i).superCluster()->energy();
		e15_1 = gsfCands.at(i).e1x5()/gsfCands.at(i).superCluster()->energy();
		e55_1 = gsfCands.at(i).e5x5()/gsfCands.at(i).superCluster()->energy();
		HoEM_1 = gsfCands.at(i).hcalOverEcal();
		sIeIe_1 = gsfCands.at(i).sigmaIetaIeta();
		netIso_1 = ( tkSumPt_1 + ecalSumEt_1 + hcalSumEt_1 )/vl1.Pt();
		
		l2Cand = gsfCands.at(j).p4();		
		ecalSumEt_2 = gsfCands.at(j).dr03EcalRecHitSumEt();
		tkSumPt_2 = gsfCands.at(j).dr03TkSumPt();
		hcalSumEt_2 = gsfCands.at(j).dr03HcalTowerSumEt();
		ecalIsoByGSF_2 = ecalSumEt_2/l2Cand.Et();
		trackIso_2 = tkSumPt_2/l2Cand.Pt();
		e25Max_2 = gsfCands.at(j).e2x5Max()/gsfCands.at(j).superCluster()->energy();
		e15_2 = gsfCands.at(j).e1x5()/gsfCands.at(j).superCluster()->energy();
		e55_2 = gsfCands.at(j).e5x5()/gsfCands.at(j).superCluster()->energy();
		HoEM_2 = gsfCands.at(j).hcalOverEcal();
		sIeIe_2 = gsfCands.at(j).sigmaIetaIeta();
		netIso_2 = ( tkSumPt_2 + ecalSumEt_2 + hcalSumEt_2 )/l1Cand.Pt();
      }
    }

    // Second electron is from the far ECAL region (no tracker)
    for (unsigned int j=0; j<ntCands.size(); j++) { 
      if (ntCands.at(j).pt() < minElePt2) break ; // List is pT ordered
       
      zCand = gsfCands.at(i).p4() + ntCands.at(j).p4() ; 
      if (zCand.M() < minMass) continue ; 
      double deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
		mZdiff = deltaM ; 
		vZ1 = zCand ; 
		Z1flavor = 3 ; 
		Z1idx.first = i ; Z1idx.second = j ;
		vl1 =  gsfCands.at(i).p4();
		ecalSumEt_1 = gsfCands.at(i).dr03EcalRecHitSumEt();
		tkSumPt_1 = gsfCands.at(i).dr03TkSumPt();
		hcalSumEt_1 = gsfCands.at(i).dr03HcalTowerSumEt();
		ecalIsoByGSF_1 = ecalSumEt_1/vl1.Et();
		trackIso_1 = tkSumPt_1/vl1.Pt();
		e25Max_1 = gsfCands.at(i).e2x5Max()/gsfCands.at(i).superCluster()->energy();
		e15_1 = gsfCands.at(i).e1x5()/gsfCands.at(i).superCluster()->energy();
		e55_1 = gsfCands.at(i).e5x5()/gsfCands.at(i).superCluster()->energy();
		HoEM_1 = gsfCands.at(i).hcalOverEcal();
		sIeIe_1 = gsfCands.at(i).sigmaIetaIeta();
		netIso_1 = ( tkSumPt_1 + ecalSumEt_1 + hcalSumEt_1 )/vl1.Pt();
		
		l2Cand = ntCands.at(j).p4();
		HoEM_2 = ntCands.at(j).hadronicOverEm();
		sIeIe_2 = ntCands.at(j).sigmaIetaIeta();
		//ecalSumEt_2 = ntCands.at(j).dr03EcalRecHitSumEt();
      }
    }

    // Second electron is from HF 
    for (unsigned int j=0; j<hfCands.size(); j++) { 
      if (hfCands.at(j).pt() < minHFpt) break ; // List is pT ordered
      
      zCand = gsfCands.at(i).p4() + hfCands.at(j).p4() ; 
      if (zCand.M() < minMass) continue ; 
      double deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
	mZdiff = deltaM ; 
	vZ1 = zCand ; 
	Z1flavor = 4 ; 
	Z1idx.first = i ; Z1idx.second = j ; 
	vl1 =  gsfCands.at(i).p4();
	ecalSumEt_1 = gsfCands.at(i).dr03EcalRecHitSumEt();
	tkSumPt_1 = gsfCands.at(i).dr03TkSumPt();
	hcalSumEt_1 = gsfCands.at(i).dr03HcalTowerSumEt();
	ecalIsoByGSF_1 = ecalSumEt_1/vl1.Et();
	trackIso_1 = tkSumPt_1/vl1.Pt();
	e25Max_1 = gsfCands.at(i).e2x5Max()/gsfCands.at(i).superCluster()->energy();
	e15_1 = gsfCands.at(i).e1x5()/gsfCands.at(i).superCluster()->energy();
	e55_1 = gsfCands.at(i).e5x5()/gsfCands.at(i).superCluster()->energy();
	HoEM_1 = gsfCands.at(i).hcalOverEcal();
	sIeIe_1 = gsfCands.at(i).sigmaIetaIeta();
	netIso_1 = ( tkSumPt_1 + ecalSumEt_1 + hcalSumEt_1 )/vl1.Pt();
		
	l2Cand = hfCands.at(j).p4();
	ecalSumEt_2 = -1;
	e25Max_2 = -1;
	e15_2 = -1;
	e55_2 = -1;
	
	//std::cout << "New best Z candidate: " << i << "," << j << std::endl ; 
      }
    }
  }

  vl1 = l1Cand;
  vl2 = l2Cand;
  return Z1flavor; 
}

int HiggsEvent::getZ2(double minElePt, double minMuPt, double minMass, double minM4) { 

  double minSumPt = 0. ;
  double sumPt;
  reco::Particle::LorentzVector zCand, l3Cand, l4Cand;
  
  if ( muCands.size() ) { //muon candidates
    for (unsigned int i=0; i<muCands.size()-1; i++) { 
      if ( Z1flavor == 1 && 
	   ( (i == Z1idx.first) || (i == Z1idx.second) ) ) continue ; 
      if (muCands.at(i).pt() < minMuPt) break ; // List is pT ordered
      for (unsigned int j=i+1; j<muCands.size(); j++) { 
	if ( Z1flavor == 1 && 
	     ( (j == Z1idx.first) || (j == Z1idx.second) ) ) continue ; 
	if (muCands.at(j).pt() < minMuPt) break ; // List is pT ordered
	if ((muCands.at(i).charge()*muCands.at(j).charge()) > 0) continue ; 
	
	zCand = muCands.at(i).p4() + muCands.at(j).p4() ;
	if (zCand.M() < minMass) continue ; 
	l4Cand = zCand + vZ1 ; 
	
	if (l4Cand.M() < minM4) continue ; 
	sumPt = muCands.at(i).pt() + muCands.at(j).pt() ; 
	if (sumPt > minSumPt) { // New "best" Z candidate
	  minSumPt = sumPt ; 
	  vZ2 = zCand ; 
	  vH = l4Cand ; 
	  Z2flavor = 1 ; 
	  Z2idx.first = i ; Z2idx.second = j ; 
	  vl3 = muCands.at(i).p4();
	  vl4 = muCands.at(j).p4();
		}
      }
    }
  }

  for (unsigned int i=0; i<gsfCands.size(); i++) { //GSF candidates
    if ( Z1flavor == 2 && ( (i == Z1idx.first) || (i == Z1idx.second) ) ) continue ; 
    if ( ((Z1flavor == 3) || (Z1flavor == 4)) && (i == Z1idx.first) ) continue ; 
    if (gsfCands.at(i).pt() < minElePt) break ; // List is pT ordered
    for (unsigned int j=i+1; j<gsfCands.size(); j++) { 
      if ( Z1flavor == 2 && 
	   ( (j == Z1idx.first) || (j == Z1idx.second) ) ) continue ; 
      if (gsfCands.at(j).pt() < minElePt) break ; // List is pT ordered
      if ((gsfCands.at(i).charge()*gsfCands.at(j).charge()) > 0) continue ;
       
      zCand = gsfCands.at(i).p4() + gsfCands.at(j).p4() ;
      if (zCand.M() < minMass) continue ; 
      if (zCand.M() > vZ1.M()) continue ; // Z1 should be mostly "on shell"      
      l4Cand = zCand + vZ1 ;       
      if (l4Cand.M() < minM4) continue ; 
      sumPt = gsfCands.at(i).pt() + gsfCands.at(j).pt() ; 
      if (sumPt > minSumPt) { // New "best" Z candidate
	minSumPt = sumPt ; 
	vZ2 = zCand ; 
	vH = l4Cand ; 
	Z2flavor = 2 ; 
	Z2idx.first = i ; Z2idx.second = j ; 
	vl3 =  gsfCands.at(i).p4();
	ecalSumEt_3 = gsfCands.at(i).dr03EcalRecHitSumEt();
	tkSumPt_3 = gsfCands.at(i).dr03TkSumPt();
	hcalSumEt_3 = gsfCands.at(i).dr03HcalTowerSumEt();
	ecalIsoByGSF_3 = ecalSumEt_3/vl3.Et();
	trackIso_3 = tkSumPt_3/vl3.Pt();
	e25Max_3 = gsfCands.at(i).e2x5Max()/gsfCands.at(i).superCluster()->energy();
	e15_3 = gsfCands.at(i).e1x5()/gsfCands.at(i).superCluster()->energy();
	e55_3 = gsfCands.at(i).e5x5()/gsfCands.at(i).superCluster()->energy();
	HoEM_3 = gsfCands.at(i).hcalOverEcal();
	sIeIe_3 = gsfCands.at(i).sigmaIetaIeta();
	netIso_3 = ( tkSumPt_3 + ecalSumEt_3 + hcalSumEt_3 )/vl3.Pt();
		
	vl4 =  gsfCands.at(j).p4();
	ecalSumEt_4 = gsfCands.at(j).dr03EcalRecHitSumEt();
	tkSumPt_4 = gsfCands.at(j).dr03TkSumPt();
	hcalSumEt_4 = gsfCands.at(j).dr03HcalTowerSumEt();
	ecalIsoByGSF_4 = ecalSumEt_4/vl4.Et();
	trackIso_4 = tkSumPt_4/vl4.Pt();
	e25Max_4 = gsfCands.at(j).e2x5Max()/gsfCands.at(j).superCluster()->energy();
	e15_4 = gsfCands.at(j).e1x5()/gsfCands.at(j).superCluster()->energy();
	e55_4 = gsfCands.at(j).e5x5()/gsfCands.at(j).superCluster()->energy();
	HoEM_4 = gsfCands.at(j).hcalOverEcal();
	sIeIe_4 = gsfCands.at(j).sigmaIetaIeta();
	netIso_4 = ( tkSumPt_4 + ecalSumEt_4 + hcalSumEt_4 )/vl4.Pt();
      }
    }
  }
  return Z2flavor; 
}

void HiggsEvent::calculate() {
	
	mZ1 = vZ1.M();
	mZ2 = vZ2.M();
	YZ1 = vZ1.Rapidity();
	YZ2 = vZ2.Rapidity();
	vH = vZ1 + vZ2 ; 
	mH = vH.M() ;
	HY = vH.Rapidity();
	l1pt = vl1.Pt();
	l1eta = vl1.Eta();
	l2pt = vl2.Pt();
	l2eta = vl2.Eta();
	l3pt = vl3.Pt();
	l3eta = vl3.Eta();
	l4pt = vl4.Pt();
	l4eta = vl4.Eta();	

  /*
  reco::Particle::LorentzVector j1p4 = j1.p4();
  reco::Particle::LorentzVector j2p4 = j2.p4();

  reco::Particle::LorentzVector lep1p4 = mu1.p4();
  reco::Particle::LorentzVector lep2p4 = (mode == CLO) ? mu1.p4() : ((mode == TOP) ? e1.p4() : mu2.p4()) ;  

  reco::Particle::Vector lep1mom = lep1p4.Vect();
  reco::Particle::Vector lep2mom = lep2p4.Vect();

  vMuMu  = lep1p4+lep2p4;
  vJJ    = j1p4+j2p4;
  lv_evt = vMuMu+vJJ;

  // LorentzVector of just the Z deboost.
  reco::Particle::LorentzVector deboostz(0,0,-lv_evt.pz(),lv_evt.pz());
  
  reco::Particle::LorentzVector lep1z=lep1p4+deboostz;
  reco::Particle::LorentzVector lep2z=lep2p4+deboostz;
  reco::Particle::LorentzVector j1z=j1p4+deboostz;
  reco::Particle::LorentzVector j2z=j2p4+deboostz;
  
  cthetaz_mumu   = planeCosAngle(lep1z.Vect(),lep2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_jj     = planeCosAngle(j1z.Vect(),j2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_mu1_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),lep1z.Vect());
  cthetaz_mu2_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),lep2z.Vect());

  float dRlep1jet1 = deltaR( mu1.eta(), mu1.phi(), j1.eta(), j1.phi() ) ; 
  float dRlep1jet2 = deltaR( mu1.eta(), mu1.phi(), j2.eta(), j2.phi() ) ; 
  float dRlep2jet1 = (( mode == TOP ) ? 
		      (deltaR( e1.eta(), e1.phi(), j1.eta(), j1.phi() )) :  
		      (( mode == CLO ) ? 
		       dRlep1jet1 : 
		       (deltaR( mu2.eta(), mu2.phi(), j1.eta(), j1.phi()))) ) ;  
  float dRlep2jet2 = (( mode == TOP ) ? 
		      (deltaR( e1.eta(), e1.phi(), j2.eta(), j2.phi() )) : 
		      (( mode == CLO ) ? 
		       dRlep1jet2 : 
		       (deltaR( mu2.eta(), mu2.phi(), j2.eta(), j2.phi()))) ) ; 
  
  // find the closest jets
  dRminMu1jet = std::min( dRlep1jet1,dRlep1jet2 );
  dRminMu2jet = std::min( dRlep2jet1,dRlep2jet2 );

  // what are the muon transverse momenta relative to the closest jets?
  reco::Particle::Vector jmom4lep1 = (dRminMu1jet == dRlep1jet1) ? j1mom : j2mom;
  reco::Particle::Vector jmom4lep2 = (dRminMu2jet == dRlep2jet1) ? j1mom : j2mom;

  TVector3 lep1vec( lep1mom.X(), lep1mom.Y(), lep1mom.Z() );
  TVector3 lep2vec( lep2mom.X(), lep2mom.Y(), lep2mom.Z() );

  TVector3 jt4lep1vec( jmom4lep1.X(), jmom4lep1.Y(), jmom4lep1.Z() );
  TVector3 jt4lep2vec( jmom4lep2.X(), jmom4lep2.Y(), jmom4lep2.Z() );

  ptrelMu1 = lep1vec.Perp( jt4lep1vec );
  ptrelMu2 = lep2vec.Perp( jt4lep2vec );

  // Composite objects
  mJJ   = vJJ.M();
  mMuMu = vMuMu.M();

  mWR   = lv_evt.M();

  mNuR1 = (vJJ + lep1p4).M();
  mNuR2 = (vJJ + lep2p4).M();
  */

}
/*
void HiggsEvent::decayID(const HepMC::GenEvent& genE) {

  HepMC::GenEvent::vertex_const_iterator vtex;
  HepMC::GenVertex::particles_out_const_iterator Pout;

  mc_class=0;

  for (vtex=genE.vertices_begin();vtex!=genE.vertices_end();vtex++){
    if(((*vtex)->particles_in_size())==1){ // We want a decay, not collision
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==9900024){ // Is a WR
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 1;
            break;
	    
          case 13:	// mu
            mc_class = 2;
            break;
	    
          case 15:  // tau
            mc_class = 3;
            break;
	    
          default:	// else
	    break;
	  }
	  if (mc_class!=0) break;
	
	}
	break; //end of id loop
      }
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==23){ // Is a Z
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 11;
            break;
	    
          case 13:	// mu
            mc_class = 12;
            break;
	    
          case 15:  // tau
            mc_class = 13;
            break;
	    
          default:	// else
            break;
	  }
	  if (mc_class!=0) break;	
	  
	}
	break; //end of id loop
      }
    }

  }
}
*/
void HiggsEvent::decayID(const reco::GenParticleCollection& gpp) {
  mc_class=0;

  reco::GenParticleCollection::const_iterator i;
  for (i=gpp.begin(); i!=gpp.end(); i++) {
    if (abs(i->pdgId())==9900024 && i->numberOfDaughters()>=2) {
      int apid=abs(i->daughter(0)->pdgId());
      if (apid==11) mc_class=1;
      if (apid==13) mc_class=2;
      if (apid==15) mc_class=3;
      break;
    }
    if (abs(i->pdgId())==23 && i->numberOfDaughters()>=2) {
      int apid=abs(i->daughter(0)->pdgId());
      if (apid==11) mc_class=11;
      if (apid==13) mc_class=12;
      if (apid==15) mc_class=13;
      break;
      /*
      std::cout << i->pdgId() << " " << i->numberOfDaughters() << "  ";
      for (unsigned int j=0; j<i->numberOfDaughters(); j++)
	std::cout  << i->daughter(j)->pdgId() << " " ;
      std::cout << std::endl;
      */
    }
  }
}
