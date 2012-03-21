#include "HZZ4Leptons/AnalysisModules/src/HiggsEvent.h"
#include "TVector3.h"

HiggsEvent::HiggsEvent() { 
  eventWgt = 1.0 ; 
  cutlevel = -1 ; 
  nMuons = 0 ; 
  nElectrons = 0 ;

  Z1flavor = 0 ; 
  Z2flavor = 0 ; 

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
int HiggsEvent::getZ1(double minElePt1, double minElePt2, double minMuPt1, double minMuPt2, double minMass) { 

  const double mZpdg = 91.1876 ; 
  double mZdiff = mZpdg ; 
  double deltaM = 0;
  reco::Particle::LorentzVector zCand;
  
  //std::cout << "Looking for the primary Z" << std::endl ; 

  if ( muCands.size() ) { //Muon channel
    for (unsigned int i=0; i<muCands.size()-1; i++) { 
      if (muCands.at(i).pt() < minMuPt1) break ; // List is pT ordered
      //std::cout << "Examining primary muon candidate with pT = " << muCands.at(i).pt() 
		//<< " GeV and charge " << muCands.at(i).charge() << std::endl ; 
      for (unsigned int j=i+1; j<muCands.size(); j++) { 
		if (muCands.at(j).pt() < minMuPt2) break ; // List is pT ordered
		//std::cout << "Examining secondary muon candidate with pT = " << muCands.at(j).pt() 
			  //<< " GeV and charge " << muCands.at(j).charge() << std::endl ; 
		if ((muCands.at(i).charge()*muCands.at(j).charge()) > 0) continue ;
		 
		zCand = muCands.at(i).p4() + muCands.at(j).p4() ; 
		
		//std::cout << "Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
		if (zCand.M() < minMass) continue ; 
		deltaM = fabs(zCand.M() - mZpdg) ; 
		if (deltaM < mZdiff) { // New "best" Z candidate
		  mZdiff = deltaM ; 
		  vZ1 = zCand ; 
		  Z1flavor = 1 ; 
		  Z1idx.first = i ; Z1idx.second = j ; 
		  //std::cout << "New best Z candidate: " << i << "," << j << std::endl ; 
		}
      }
    }
  }

  for (unsigned int i=0; i<gsfCands.size(); i++) { //tracked Electron channel
    if (gsfCands.at(i).pt() < minElePt1) break ; // List is pT ordered
    //std::cout << "Examining primary GSF candidate with pT = " << gsfCands.at(i).pt() 
	      //<< " GeV and charge " << gsfCands.at(i).charge() << std::endl ; 

    // Second electron is a GSF electron
    for (unsigned int j=i+1; j<gsfCands.size(); j++) { 
      if (gsfCands.at(j).pt() < minElePt2) break ; // List is pT ordered
      //std::cout << "Examining secondary GSF candidate with pT = " << gsfCands.at(j).pt() 
		//<< " GeV and charge " << gsfCands.at(j).charge() << std::endl ; 
      if ((gsfCands.at(i).charge()*gsfCands.at(j).charge()) > 0) continue ;
       
      zCand = gsfCands.at(i).p4() + gsfCands.at(j).p4() ; 
      
      //std::cout << "Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
      if (zCand.M() < minMass) continue ; 
      deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
		mZdiff = deltaM ; 
		vZ1 = zCand ; 
		Z1flavor = 2 ; 
		Z1idx.first = i ; Z1idx.second = j ; 
		//std::cout << "New best Z candidate: " << i << "," << j << std::endl ; 
      }
    }

    // Second electron is from the far ECAL region (no tracker)
    for (unsigned int j=0; j<ntCands.size(); j++) { 
      if (ntCands.at(j).pt() < minElePt2) break ; // List is pT ordered
      //std::cout << "Examining secondary NT candidate with pT = " << ntCands.at(j).pt() << std::endl ;
       
      zCand = gsfCands.at(i).p4() + ntCands.at(j).p4() ; 
      
      //std::cout << "Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
      if (zCand.M() < minMass) continue ; 
      double deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
	mZdiff = deltaM ; 
	vZ1 = zCand ; 
	Z1flavor = 3 ; 
	Z1idx.first = i ; Z1idx.second = j ; 
	//std::cout << "New best Z candidate: " << i << "," << j << std::endl ; 
      }
    }

    // Second electron is from HF 
    for (unsigned int j=0; j<hfCands.size(); j++) { 
      if (hfCands.at(j).pt() < minElePt2) break ; // List is pT ordered
      //std::cout << "Examining secondary HF candidate with pT = " << hfCands.at(j).pt() << std::endl ; 
      
      zCand = gsfCands.at(i).p4() + hfCands.at(j).p4() ; 
      
      //std::cout << "Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
      if (zCand.M() < minMass) continue ; 
      double deltaM = fabs(zCand.M() - mZpdg) ; 
      if (deltaM < mZdiff) { // New "best" Z candidate
	mZdiff = deltaM ; 
	vZ1 = zCand ; 
	Z1flavor = 4 ; 
	Z1idx.first = i ; Z1idx.second = j ; 
	//std::cout << "New best Z candidate: " << i << "," << j << std::endl ; 
      }
    }
  }

  return Z1flavor; 
}

int HiggsEvent::getZ2(double minElePt, double minMuPt, double minMass, double minM4) { 

  //std::cout << "Looking for the Z*" << std::endl ; 

  double minSumPt = 0. ;
  double sumPt;
  reco::Particle::LorentzVector zCand, l4Cand;
  
  if ( muCands.size() ) { 
    for (unsigned int i=0; i<muCands.size()-1; i++) { 
      if ( Z1flavor == 1 && 
	   ( (i == Z1idx.first) || (i == Z1idx.second) ) ) continue ; 
      if (muCands.at(i).pt() < minMuPt) break ; // List is pT ordered
      //std::cout << "(2) Examining primary muon candidate with pT = " << muCands.at(i).pt() 
		//<< " GeV and charge " << muCands.at(i).charge() << std::endl ; 
      for (unsigned int j=i+1; j<muCands.size(); j++) { 
	if ( Z1flavor == 1 && 
	     ( (j == Z1idx.first) || (j == Z1idx.second) ) ) continue ; 
	if (muCands.at(j).pt() < minMuPt) break ; // List is pT ordered
	//std::cout << "(2) Examining secondary muon candidate with pT = " << muCands.at(j).pt() 
		  //<< " GeV and charge " << muCands.at(j).charge() << std::endl ; 
	if ((muCands.at(i).charge()*muCands.at(j).charge()) > 0) continue ; 
	
	zCand = muCands.at(i).p4() + muCands.at(j).p4() ; 
	
	//std::cout << "(2) Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
	if (zCand.M() < minMass) continue ; 
	l4Cand = zCand + vZ1 ; 
	
	//std::cout << "4 lepton candidate mass is " << l4Cand.M() << " GeV" << std::endl ; 
	if (l4Cand.M() < minM4) continue ; 
	sumPt = muCands.at(i).pt() + muCands.at(j).pt() ; 
	if (sumPt > minSumPt) { // New "best" Z candidate
	  minSumPt = sumPt ; 
	  vZ2 = zCand ; 
	  vH = l4Cand ; 
	  Z2flavor = 1 ; 
	  Z2idx.first = i ; Z2idx.second = j ; 
	  //std::cout << "(2) New best Z candidate: " << i << "," << j << std::endl ; 
		}
      }
    }
  }

  for (unsigned int i=0; i<gsfCands.size(); i++) { 
    if ( Z1flavor == 2 && ( (i == Z1idx.first) || (i == Z1idx.second) ) ) continue ; 
    if ( ((Z1flavor == 3) || (Z1flavor == 4)) && (i == Z1idx.first) ) continue ; 
    if (gsfCands.at(i).pt() < minElePt) break ; // List is pT ordered
    //std::cout << "(2) Examining primary GSF candidate with pT = " << gsfCands.at(i).pt() 
	      //<< " GeV and charge " << gsfCands.at(i).charge() << std::endl ; 
    for (unsigned int j=i+1; j<gsfCands.size(); j++) { 
      if ( Z1flavor == 2 && 
	   ( (j == Z1idx.first) || (j == Z1idx.second) ) ) continue ; 
      if (gsfCands.at(j).pt() < minElePt) break ; // List is pT ordered
      //std::cout << "Examining secondary GSF candidate with pT = " << gsfCands.at(j).pt() 
		//<< " GeV and charge " << gsfCands.at(j).charge() << std::endl ; 
      if ((gsfCands.at(i).charge()*gsfCands.at(j).charge()) > 0) continue ;
       
      zCand = gsfCands.at(i).p4() + gsfCands.at(j).p4() ;
       
      //std::cout << "(2) Z candidate mass is " << zCand.M() << " GeV" << std::endl ; 
      if (zCand.M() < minMass) continue ; 
      if (zCand.M() > vZ1.M()) continue ; // Z1 should be mostly "on shell"
      
      l4Cand = zCand + vZ1 ; 
      
      //std::cout << "4 lepton candidate mass is " << l4Cand.M() << " GeV" << std::endl ; 
      if (l4Cand.M() < minM4) continue ; 
      sumPt = gsfCands.at(i).pt() + gsfCands.at(j).pt() ; 
      if (sumPt > minSumPt) { // New "best" Z candidate
	minSumPt = sumPt ; 
	vZ2 = zCand ; 
	vH = l4Cand ; 
	Z2flavor = 2 ; 
	Z2idx.first = i ; Z2idx.second = j ; 
	//std::cout << "(2) New best Z candidate: " << i << "," << j << std::endl ; 
      }
    }
  }

  return Z2flavor; 
}

void HiggsEvent::calculate() {

  vH = vZ1 + vZ2 ; 
  mH = vH.M() ;
  
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
