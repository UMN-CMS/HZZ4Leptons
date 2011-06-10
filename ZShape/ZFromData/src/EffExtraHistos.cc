#include "ZShape/ZFromData/interface/EffExtraHistos.h"
#include "ZShape/Base/interface/ZShapeBinning.h"

typedef math::XYZTLorentzVector XYZTLorentzVector;

// for reference
// $ROOTSYS/include/Math/GenVector/LorentzVector.h

void EffExtraHistos::Book(TFileDirectory& tdf) {
 
  // add here all extremes, tidily
  float pi       = 3.141593;
  float maxPt    = 400;
  float minZmass = 50;
  float maxZmass = 140;
  const int   binsZmass  = 50;


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
  DYZ_  = tdf.make<TH1F>("Z0_Y_Delta","Z0_Y_Delta;Y_{Z0}", int((zshape::y_max*2)*4), -zshape::y_max, zshape::y_max);  
  DptZ_ = tdf.make<TH1F>("Z0_Pt_Delta","Z0_Pt_Delta;p_{T,Z0}", 100, -50., 50.); 



  De1eta_ = tdf.make<TH1F>("e1_eta_Delta","e1_eta_Delta;detector #eta_{e1}", 50, -5, 5);  
  De2eta_ = tdf.make<TH1F>("e2_eta_Delta","e2_eta_Delta;detector #eta_{e2}", 50, -5, 5);
  De1phi_ = tdf.make<TH1F>("e1_phi_Delta","e1_phi_Delta;#phi_{e1}", 50, -pi, pi);  
  De2phi_ = tdf.make<TH1F>("e2_phi_Delta","e2_phi_Delta;#phi_{e2}", 50,  -pi, pi);  
  De1pt_  = tdf.make<TH1F>("e1_P_t_Delta","e1_P_t_Delta;p_{T,e1}", 200, -50., 50.);  
  De2pt_  = tdf.make<TH1F>("e2_P_t_Delta","e2_P_t_Delta;p_{T,e2}", 200, -50., 50.);

  atZ_ = tdf.make<TH1F>("wZ0_at","Z0_at;a_{T,Z0}", 200, 0, 150.);
  alZ_ = tdf.make<TH1F>("wZ0_al","Z0_al;a_{T,Z0}", 200, 0, 150.);
  ptZ_ = tdf.make<TH1F>("wZ0_pt","Z0_Pt;p_{T,Z0}", 200, 0, maxPt);

  
  MCatZ_ = tdf.make<TH1F>("Z0_at_MC","Z0_at_MC;a_{T,Z0}", 200, 0, 150.);
  MCalZ_ = tdf.make<TH1F>("Z0_al_MC","Z0_al_MC;a_{T,Z0}", 200, 0, 150.);


  MCmZ_  = tdf.make<TH1F>("Z0_mass_MC","Z0_mass_MC;m_{Z0} (GeV/c^{2})", 60, 35., 150.);  
  MCYZ_  = tdf.make<TH1F>("Z0_Y_MC","Z0_Y_MC;Y_{Z0}", int((zshape::y_max*2)*4), -zshape::y_max, zshape::y_max);  
  MCptZ_ = tdf.make<TH1F>("Z0_Pt_MC","Z0_Pt_MC;p_{T,Z0}", 200, 0., maxPt); 
  MCe1eta_ = tdf.make<TH1F>("e1_eta_MC","e1_eta_MC;detector #eta_{e1}", 50, -5, 5);  
  MCe2eta_ = tdf.make<TH1F>("e2_eta_MC","e2_eta_MC;detector #eta_{e2}", 50, -5, 5);
  MCeeta_ = tdf.make<TH1F>("e_eta_MC","e_eta_MC;detector #eta_{e}", 50, -5, 5);  
  MCe1phi_ = tdf.make<TH1F>("e1_phi_MC","e1_phi_MC;#phi_{e1}", 50, -pi, pi);  
  MCe2phi_ = tdf.make<TH1F>("e2_phi_MC","e2_phi_MC;#phi_{e2}", 50,  -pi, pi);  
  MCephi_ = tdf.make<TH1F>("e_phi_MC","e_phi_MC;#phi_{e}", 50,  -pi, pi);  
  MCe1pt_  = tdf.make<TH1F>("e1_P_t_MC","e1_P_t_MC;p_{T,e1}", 200, 0., maxPt);  
  MCe2pt_  = tdf.make<TH1F>("e2_P_t_MC","e2_P_t_MC;p_{T,e2}", 200, 0., maxPt);

  sc_mZ_  = tdf.make<TH1F>("Z0_mass_sc","Z0_mass_sc;m_{Z} (GeV/c^{2})", binsZmass, minZmass, maxZmass);
  sc_Y_  = tdf.make<TH1F>("Z0_Y_sc","Z0_Y_sc;Y_{Z0}", zshape::y_bins, -zshape::y_max, zshape::y_max);
  sc_pt_  = tdf.make<TH1F>("Z0_Pt_sc","Z0_Pt_sc;p_{T,Z0}", zshape::pt_bins, zshape::pt_binning);
  sc_ptMon_ = tdf.make<TH1F>("Z0_PtMon_sc","Z0_PtMon_sc;p_{T,Z0}", 200, 0, maxPt);


}

void EffExtraHistos::Fill(const ZShapeElectron& e1, const ZShapeElectron& e2, const ZShapeElectron& em1, const ZShapeElectron& em2,double wgt, bool doMC) {

  XYZTLorentzVector p1(e1.p4_);
  XYZTLorentzVector p2(e2.p4_);

  XYZTLorentzVector pZ = p1 + p2 ;

  math::PtEtaPhiMLorentzVector sce1(e1.scEnergy_/cosh(e1.detEtaVtx_),e1.detEtaVtx_,e1.detPhi_,0);
  math::PtEtaPhiMLorentzVector sce2(e2.scEnergy_/cosh(e2.detEtaVtx_),e2.detEtaVtx_,e2.detPhi_,0);

  math::PtEtaPhiMLorentzVector scZ = sce1 + sce2;

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
  
  //double mzY    = pmZ.Rapidity();
  //double mzPt   = pmZ.Pt();

 
  //calculate at, al, t

  
  math::XYZVector pt1(p1.Vect().X(),p1.Vect().Y(),0);
  math::XYZVector pt2(p2.Vect().X(),p2.Vect().Y(),0);
  math::XYZVector Pt=pt1+pt2;
  math::XYZVector t=(pt1-pt2);
  t*=1.0/(t.r());
  double at=(Pt.Cross(t)).r();
  double al=Pt.Dot(t);

  //calculate at, al, t for MC, (i.e. atm, alm, tm)
 
  math::XYZVector ptm1(pm1.Vect().X(),pm1.Vect().Y(),0);
  math::XYZVector ptm2(pm2.Vect().X(),pm2.Vect().Y(),0);
  math::XYZVector Ptm=ptm1+ptm2;
  math::XYZVector tm=(ptm1-ptm2);
  tm*=1.0/(tm.r());

  double atm=(Ptm.Cross(tm)).r();
  double alm=Ptm.Dot(tm);

 

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


  sc_mZ_ -> Fill(scZ.M(),wgt);
  sc_Y_ -> Fill(scZ.Rapidity(),wgt);
  sc_pt_ -> Fill(scZ.Pt(),wgt);
  sc_ptMon_ -> Fill(scZ.Pt(),wgt);

  if (doMC) 
    {
      DmZ_  -> Fill(zMass-zmMass,wgt);
      DYZ_  -> Fill(zY-zmY,wgt);
      DptZ_ -> Fill(zPt-zmPt,wgt);
  

      De1eta_ -> Fill(e1eta-em1eta,wgt);
      De1pt_  -> Fill(e1Pt-em1Pt,wgt);
      De1phi_ -> Fill(e1phi-em1phi,wgt);

      De2eta_ -> Fill(e2eta-em2eta,wgt);
      De2pt_  -> Fill(e2Pt-em2Pt,wgt);
      De2phi_ -> Fill(e2phi-em2phi,wgt);
 
      MCmZ_  -> Fill(zmMass,wgt);
      MCYZ_  -> Fill(zmY,wgt);
      MCptZ_ -> Fill(zmPt,wgt);

      MCe1eta_ -> Fill(em1eta,wgt);
      MCe1pt_  -> Fill(em1Pt,wgt);
      MCe1phi_ -> Fill(em1phi,wgt);

      MCe2eta_ -> Fill(em2eta,wgt);
      MCe2pt_  -> Fill(em2Pt,wgt);
      MCe2phi_ -> Fill(em2phi,wgt);

      MCeeta_ -> Fill(em1eta,wgt);
      MCeeta_ -> Fill(em2eta,wgt);
      MCephi_ -> Fill(em1phi,wgt);
      MCephi_ -> Fill(em2phi,wgt);
 
      MCatZ_->Fill(atm,wgt); 
      MCalZ_->Fill(alm,wgt);
    }
  atZ_->Fill(at,wgt); 
  alZ_->Fill(al,wgt);
  ptZ_->Fill(zPt,wgt);
 

}
