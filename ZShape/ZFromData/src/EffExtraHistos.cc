#include "ZShape/ZFromData/interface/EffExtraHistos.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffExtraHistos::Book(TFileDirectory& tdf) {

  // add here all extremes, tidily
  float pi       = 3.141593;
  float maxPt    = 400;
  float maxY     = 5.5;
  float minZmass = 50;
  float maxZmass = 140;


  float detetaBins[104] = {-5., -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4., -3.9, -3.8, -3.7, -3.6,  //14
                       -3.5, -3.4, -3.3,  -3.2, -3.1, -3., //6
                       -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.89, -1.78, -1.67,  //13
                       -1.56, -1.4442, -1.35, -1.257, -1.163, -1.127,  -1.02, -0.913, -0.806, -0.770,  //10
                       -0.667, -.564, -0.461, -0.423, -0.32, -0.22, -0.12, -0.018, 0.018, 0.12, 0.22, 0.32,  0.423, 0.461, 0.564,0.667, //16
                       0.770, 0.806, 0.913, 1.02, 1.127, 1.163,  1.257, 1.35, 1.4442, //9
                       1.56, 1.67, 1.78, 1.89, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, //14
                       3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, //10
                       4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5. }; //11
 
 
  DmZ_  = tdf.make<TH1F>("Z0_mass_Delta","Z0_mass_Delta;Delta m_{Z0} (GeV/c^{2})", 50, -60, 60);  
  DYZ_  = tdf.make<TH1F>("Z0_Y_Delta","Z0_Y_Delta;Y_{Z0}", int((maxY*2)*4), -maxY, maxY);  
  DptZ_ = tdf.make<TH1F>("Z0_Pt_Delta","Z0_Pt_Delta;p_{T,Z0}", 100, -50., 50.); 
  De1eta_ = tdf.make<TH1F>("e1_eta_Delta","e1_eta_Delta;detector #eta_{e1}", 103, detetaBins);  
  De2eta_ = tdf.make<TH1F>("e2_eta_Delta","e2_eta_Delta;detector #eta_{e2}", 103, detetaBins);
  De1phi_ = tdf.make<TH1F>("e1_phi_Delta","e1_phi_Delta;#phi_{e1}", 100, -pi, pi);  
  De2phi_ = tdf.make<TH1F>("e2_phi_Delta","e2_phi_Delta;#phi_{e2}", 100,  -pi, pi);  
  De1pt_  = tdf.make<TH1F>("e1_P_t_Delta","e1_P_t_Delta;p_{T,e1}", 200, -50., 50.);  
  De2pt_  = tdf.make<TH1F>("e2_P_t_Delta","e2_P_t_Delta;p_{T,e2}", 200, -50., 50.);

  MCmZ_  = tdf.make<TH1F>("Z0_mass_MC","Z0_mass_MC;m_{Z0} (GeV/c^{2})", 60, 35., 150.);  
  MCYZ_  = tdf.make<TH1F>("Z0_Y_MC","Z0_Y_MC;Y_{Z0}", int((maxY*2)*4), -maxY, maxY);  
  MCptZ_ = tdf.make<TH1F>("Z0_Pt_MC","Z0_Pt_MC;p_{T,Z0}", 200, 0., maxPt); 
  MCe1eta_ = tdf.make<TH1F>("e1_eta_MC","e1_eta_MC;detector #eta_{e1}", 103, detetaBins);  
  MCe2eta_ = tdf.make<TH1F>("e2_eta_MC","e2_eta_MC;detector #eta_{e2}", 103, detetaBins);
  MCe1phi_ = tdf.make<TH1F>("e1_phi_MC","e1_phi_MC;#phi_{e1}", 100, -pi, pi);  
  MCe2phi_ = tdf.make<TH1F>("e2_phi_MC","e2_phi_MC;#phi_{e2}", 100,  -pi, pi);  
  MCe1pt_  = tdf.make<TH1F>("e1_P_t_MC","e1_P_t_MC;p_{T,e1}", 200, 0., maxPt);  
  MCe2pt_  = tdf.make<TH1F>("e2_P_t_MC","e2_P_t_MC;p_{T,e2}", 200, 0., maxPt);
}

void EffExtraHistos::Fill(const ZShapeElectron& e1, const ZShapeElectron& e2, const ZShapeElectron& em1, const ZShapeElectron& em2) {

  XYZTLorentzVector p1(e1.p4_);
  XYZTLorentzVector p2(e2.p4_);

  XYZTLorentzVector pZ = p1 + p2 ;

  XYZTLorentzVector pm1(em1.p4_);
  XYZTLorentzVector pm2(em2.p4_);

  XYZTLorentzVector pmZ = pm1 + pm2 ;
  //  debug
//  std::cout << "this is e1:"  << p1 << std::endl;
//  std::cout << "this is e2:"  << p2 << std::endl;
//  std::cout << "this is Z:"  << pZ << std::endl;
//  float rapidity = 0.5 * log( ( pZ.E()+pZ.Pz() ) / ( pZ.E()-pZ.Pz() ) );  
//  std::cout << "rapidity calculated: " << rapidity  << " while from method: " << pZ.Rapidity() << std::endl;
  
  float zMass = sqrt ( pZ.Dot(pZ) );
  float zY    = pZ.Rapidity();
  float zPt   = pZ.Pt();

  float e1eta = e1.detEta_;
  float e1phi = p1.Phi();
  float e1Pt  = p1.Pt();
  
  float e2eta = e2.detEta_;
  float e2phi = p2.Phi();
  float e2Pt  = p2.Pt();

  float zmMass = sqrt ( pmZ.Dot(pmZ) );
  float zmY    = pmZ.Rapidity();
  float zmPt   = pmZ.Pt();

  float em1eta = em1.detEta_;
  float em1phi = pm1.Phi();
  float em1Pt  = pm1.Pt();
  
  float em2eta = em2.detEta_;
  float em2phi = pm2.Phi();
  float em2Pt  = pm2.Pt();

  DmZ_  -> Fill(zMass-zmMass);
  DYZ_  -> Fill(zY-zmY);
  DptZ_ -> Fill(zPt-zmPt);

  De1eta_ -> Fill(e1eta-em1eta);
  De1pt_  -> Fill(e1Pt-em1Pt);
  De1phi_ -> Fill(e1phi-em1phi);

  De2eta_ -> Fill(e2eta-em2eta);
  De2pt_  -> Fill(e2Pt-em2Pt);
  De2phi_ -> Fill(e2phi-em2phi);

  MCmZ_  -> Fill(zmMass);
  MCYZ_  -> Fill(zmY);
  MCptZ_ -> Fill(zmPt);

  MCe1eta_ -> Fill(em1eta);
  MCe1pt_  -> Fill(em1Pt);
  MCe1phi_ -> Fill(em1phi);

  MCe2eta_ -> Fill(em2eta);
  MCe2pt_  -> Fill(em2Pt);
  MCe2phi_ -> Fill(em2phi);


}
