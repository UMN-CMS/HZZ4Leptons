#include "Hgg/ClusteringWithPU/interface/SCwithPUhistos.h"
#include "Hgg/ClusteringWithPU/interface/Utils.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"


SCwithPUhistos::SCwithPUhistos(bool useRawEnergy) { 
  useRawEnergy_ = useRawEnergy;
}

SCwithPUhistos::SCwithPUhistos() { 
  useRawEnergy_ = false;
}

void SCwithPUhistos::setRawEnergy(bool useRawEnergy){
  useRawEnergy_ = useRawEnergy;
}

void SCwithPUhistos::Book(TFileDirectory& fs) {
  
  if(useRawEnergy_) std::cout << "\n[SCwithPUhistos] raw energy of SC will be used\n" << std::endl;
  else              std::cout << "\n[SCwithPUhistos] corrected energy of SC will be used\n" << std::endl;

  h_scet_barl = fs.make<TH1F>("h_scet_barl","SC ET over true, barrel; EB: E_{T,superC}/E_{T,true}",90,0.75,1.2);  // using: h_scet_barl as counter of EB matched photons 
  h_scet_endc = fs.make<TH1F>("h_scet_endc","SC ET over true, endcap; EE: E_{T,superC}/E_{T,true}",90,0.75,1.2);  // using: h_scet_endc as counter of EE matched photons 
  h_EoverEtrue_barl = fs.make<TH1F>("h_EoverEtrue_barl","E_SC/Etrue, barrel; EB: E_{superC}/E_{true}",90,0.75,1.2);
  h_EoverEtrue_endc = fs.make<TH1F>("h_EoverEtrue_endc","E_SC/Etrue, endcap; EE: E_{superC}/E_{true}",90,0.75,1.2);
  h_E5x5overEtrue_barl = fs.make<TH1F>("h_E5x5overEtrue_barl","E5x5/Etrue, barrel; EB: E_{5x5}/E_{true}",90,0.75,1.2);
  h_PhoEoverEtrue_barl = fs.make<TH1F>("h_PhoEoverEtrue_barl","E_Photon/Etrue, barrel; EB: E_{photon}/E_{true}",90,0.75,1.2);
  h_E5x5R9overEtrue_barl = fs.make<TH1F>("h_E5x5R9overEtrue_barl","E5x5/Etrue for R9>0.94, barrel; EB: E_{5x5, r9>0.94}/E_{true}",90,0.75,1.2);
  h_PhoER9overEtrue_barl = fs.make<TH1F>("h_PhoER9overEtrue_barl","E_Photon/Etrue for R9>0.94, barrel; EB: E_{photon, r9>0.94}/E_{true}",90,0.75,1.2); // counter r9 EB
  h_E5x5notR9overEtrue_barl = fs.make<TH1F>("h_E5x5notR9overEtrue_barl","E5x5/Etrue for R9<0.94, barrel; EB: E_{5x5, r9<0.94}/E_{true}",90,0.75,1.2);
  h_PhoEnotR9overEtrue_barl = fs.make<TH1F>("h_PhoEnotR9overEtrue_barl","E_Photon/Etrue for R9<0.94, barrel; EB: E_{photon, r9<0.94}/E_{true}",90,0.75,1.2); // counter !r9 EB
  h_E5x5overEtrue_endc = fs.make<TH1F>("h_E5x5overEtrue_endc","E5x5/Etrue, endcap; EE: E_{5x5}/E_{true}",90,0.75,1.2);
  h_PhoEoverEtrue_endc = fs.make<TH1F>("h_PhoEoverEtrue_endc","E_Photon/Etrue, endcap; EE: E_{photon}/E_{true}",90,0.75,1.2);
  h_E5x5R9overEtrue_endc = fs.make<TH1F>("h_E5x5R9overEtrue_endc","E5x5/Etrue for R9>0.95, endcap; EE: E_{5x5, r9>0.95}/E_{true}",90,0.75,1.2);
  h_PhoER9overEtrue_endc = fs.make<TH1F>("h_PhoER9overEtrue_endc","E_Photon/Etrue for R9>0.95, endcap; EE: E_{photon, r9>0.95}/E_{true}",90,0.75,1.2); // counter r9 EE
  h_E5x5notR9overEtrue_endc = fs.make<TH1F>("h_E5x5notR9overEtrue_endc","E5x5/Etrue for R9<0.95, endcap; EE: E_{5x5, r9<0.95}/E_{true}",90,0.75,1.2);
  h_PhoEnotR9overEtrue_endc = fs.make<TH1F>("h_PhoEnotR9overEtrue_endc","E_Photon/Etrue for R9<0.95, endcap; EE: E_{photon, r9<0.95}/E_{true}",90,0.75,1.2); // counter !r9 EE
  h_mHiggs_EBEB = fs.make<TH1F>("h_mHiggs_EBEB","2 photon invariant mass, EBEB; EB-EB: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEE = fs.make<TH1F>("h_mHiggs_EBEE","2 photon invariant mass, EBEE; EB-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EEEE = fs.make<TH1F>("h_mHiggs_EEEE","2 photon invariant mass, EEEE; EE-EE: m_{#gamma#gamma} GeV/c^{2}",120,100.,140.);
  h_mHiggs_EBEB_trueVtx = fs.make<TH1F>("h_mHiggs_EBEB_trueVtx","2 photon invariant mass, EBEB, true PV; EB-EB: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EBEE_trueVtx = fs.make<TH1F>("h_mHiggs_EBEE_trueVtx","2 photon invariant mass, EBEE, true PV; EB-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_mHiggs_EEEE_trueVtx = fs.make<TH1F>("h_mHiggs_EEEE_trueVtx","2 photon invariant mass, EEEE, true PV; EE-EE: m_{#gamma#gamma} GeV/c^{2} (true vtx)",120,100.,140.);
  h_E5x5overEtrueVsEphoEtrue_barl = fs.make<TH2F>("h_E5x5overEtrueVsEphoEtrue_barl","E5x5/Etrue vs Epho/Etrue, barrel;EB: E_{5x5}/E_{true}; EB: E_{photon}/E_{true}",90,0.75,1.2,90,0.75,1.2);
  h_E5x5overEtrueVsEphoEtrue_endc = fs.make<TH2F>("h_E5x5overEtrueVsEphoEtrue_endc","E5x5/Etrue vs Epho/Etrue, endcap;EE: E_{5x5}/E_{true}; EE: E_{photon}/E_{true}",90,0.75,1.2,90,0.75,1.2);
  h_phiWidth = fs.make<TH1F>("h_phiWidth_barl","phi Width (barrel); EB: i#phi width", 100,0,0.2);
  h_phiWidthVsE = fs.make<TH2F>("h_phiWidthVsE_barl","phi Width Vs. E (Barrel); EB: i#phi width; EB: E [GeV]", 100,0,0.2,100,0,200);
  h_phiWidth_endc = fs.make<TH1F>("h_phiWidth_endc","phi Width (Endcap); EE: i#phi width", 100,0,0.2);
  h_phiWidthVsE_endc = fs.make<TH2F>("h_phiWidthVsE_endc","phi Width Vs. E (Endcap); EB: i#phi width; EB: E [GeV]", 100,0,0.2,100,0,200);
  h_phiSize         = fs.make<TH1F>("h_phiSize_barl","phi Size (barrel); EB: i#phi size", 50,0,50);
  h_phiSizeVsE      = fs.make<TH2F>("h_phiSizeVsE_barl","phi Size Vs. E (Barrel); EB: i#phi size; EB: E [GeV]", 50,0,50.,100,0,200);
  h_phiSizeVsEt     = fs.make<TH2F>("h_phiSizeVsEt_barl","phi Size Vs. E_{T} (Barrel); EB: i#phi size; EB: E_{T} [GeV]", 50,0,50.,100,0,200);
  h_phiSize_endc    = fs.make<TH1F>("h_phiSize_endc","phi Size (Endcap); EE: i#phi size", 50,0,1);
  h_phiSizeVsE_endc = fs.make<TH2F>("h_phiSizeVsE_endc","phi Size Vs. E (Endcap); EB: i#phi size; EE: E [GeV]", 50,0,1.,100,0,200);
  h_phiSizeVsEt_endc= fs.make<TH2F>("h_phiSizeVsEt_endc","phi Size Vs. E_{T} (Endcap); EB: i#phi size; EE: E_{T} [GeV]", 50,0,1.,100,0,200);
  h_phiShape_barl       = fs.make<TH1F>("h_phiShape_barl","phi Shape (barrel); EB i#phi - i#phi_{SCseed}", 35,-17,18); 
  h_absPhiShape_barl    = fs.make<TH1F>("h_absPhiShape_barl","phi AbsShape (barrel); EB  abs(i#phi - i#phi_{SCseed})", 18,0,18); 
  h_phiShape_endc       = fs.make<TH1F>("h_phiShape_endc","phi Shape (endcap) EE i#phi - i#phi_{SCseed}", 35,-17,18); 
  h_absPhiShape_endc    = fs.make<TH1F>("h_absPhiShape_endc","phi AbsShape (endcap); EE  abs(i#phi - i#phi_{SCseed})", 18,0,18); 
  h_phiShapeVsE_barl    = fs.make<TH2F>("h_phiShapeVsE_barl","phi Shape Vs E (barrel); EB i#phi - i#phi_{SCseed}; EB E [GeV]", 35,-17,18,100,0,200); 
  h_absPhiShapeVsE_barl = fs.make<TH2F>("h_absPhiShapeVsE_barl","phi AbsShape Vs E (barrel); EB  abs(i#phi - i#phi_{SCseed}); EB E [GeV]", 18,0,18,100,0,200); 
  h_phiShapeVsE_endc    = fs.make<TH2F>("h_phiShapeVsE_endc","phi Shape Vs E (endcap); EE i#phi - i#phi_{SCseed}; EE E [GeV]", 35,-17,18,100,0,200); 
  h_absPhiShapeVsE_endc = fs.make<TH2F>("h_absPhiShapeVsE_endc","phi AbsShape Vs E (endcap); EE  abs(i#phi - i#phi_{SCseed}); EE E [GeV] ", 18,0,18,100,0,200); 

  h_etaShape_barl       = fs.make<TH1F>("h_etaShape_barl","eta Shape (barrel); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlPLus   = fs.make<TH1F>("h_etaShape_barlPLus","eta Shape (barrel plus); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlMinus  = fs.make<TH1F>("h_etaShape_barlMinus","eta Shape (barrel minus); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaShape_barlSymm   = fs.make<TH1F>("h_etaShape_barlSymm","eta Shape (barrel symm); EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3); 
  h_etaPhiShape_barl    = fs.make<TH2F>("h_etaPhiShape_barl","eta Shape (barrel); EB i#phi_{BC} - phi_{SCseed}; EB i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlPLus  = fs.make<TH2F>("h_etaPhiShape_barlPlus","eta Shape (barrel plus); EB+ i#phi_{BC} - phi_{SCseed}; EB+ i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlMinus = fs.make<TH2F>("h_etaPhiShape_barlMinus","eta Shape (barrel minus); EB- i#phi_{BC} - phi_{SCseed}; EB- i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 
  h_etaPhiShape_barlSymm  = fs.make<TH2F>("h_etaPhiShape_barlSymm","eta Shape (barrel symm); EBsymm i#phi_{BC} - phi_{SCseed}; EBsymm i#eta_{BC} - i#eta_{SCseed}", 7,-3,3,35,-17,18); 

  // each phi bin is half crystal wide; each eta bin is 1/10 crystal wide
  h_phiBCminusPhiSeed_barl= fs.make<TH1F>("h_phiBCminusPhiSeed_barl","h_phiBCminusPhiSeed_barl; EB #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barl= fs.make<TH1F>("h_etaBCminusEtaSeed_barl","h_etaBCminusEtaSeed_barl; EB #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barl= fs.make<TH1F>("h_absEtaBCminusAbsEtaSeed_barl","h_absEtaBCminusAbsEtaSeed_barl; EB |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl = fs.make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EB #phi_{BC} - #phi_{SCseed}; EB |#eta_{BC}| - |#eta_{SCseed}| ", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm1= fs.make<TH1F>("h_phiBCminusPhiSeed_barlm1","h_phiBCminusPhiSeed_barl; EBm1  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm1= fs.make<TH1F>("h_etaBCminusEtaSeed_barlm1","h_etaBCminusEtaSeed_barl; EBm1  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm1= fs.make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm1","h_absEtaBCminusAbsEtaSeed_barl; EBm1  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1 = fs.make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl;  EBm1 #phi_{BC} - #phi_{SCseed}; EBm1 |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm2= fs.make<TH1F>("h_phiBCminusPhiSeed_barlm2","h_phiBCminusPhiSeed_barlm2; EBm2 #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm2= fs.make<TH1F>("h_etaBCminusEtaSeed_barlm2","h_etaBCminusEtaSeed_barlm2; EBm2 #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm2= fs.make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm2","h_absEtaBCminusAbsEtaSeed_barlm2; EBm2 |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2 = fs.make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2;  EBm2  #phi_{BC} - #phi_{SCseed}; EBm2  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm3= fs.make<TH1F>("h_phiBCminusPhiSeed_barlm3","h_phiBCminusPhiSeed_barlm3; EBm3  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm3= fs.make<TH1F>("h_etaBCminusEtaSeed_barlm3","h_etaBCminusEtaSeed_barlm3; EBm3  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm3= fs.make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm3","h_absEtaBCminusAbsEtaSeed_barlm3; EBm3  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3 = fs.make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3;  EBm3  #phi_{BC} - #phi_{SCseed}; EBm3  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));
  h_phiBCminusPhiSeed_barlm4= fs.make<TH1F>("h_phiBCminusPhiSeed_barlm4","h_phiBCminusPhiSeed_barlm4; EBm4  #phi_{BC} - #phi_{SCseed}", 80,-(6.28/360*20),(6.28/360*20)); 
  h_etaBCminusEtaSeed_barlm4= fs.make<TH1F>("h_etaBCminusEtaSeed_barlm4","h_etaBCminusEtaSeed_barlm4; EBm4  #eta_{BC} - #eta_{SCseed}", 60,-(6.28/360*3),(6.28/360*3)); 
  h_absEtaBCminusAbsEtaSeed_barlm4= fs.make<TH1F>("h_absEtaBCminusAbsEtaSeed_barlm4","h_absEtaBCminusAbsEtaSeed_barlm4; EBm4  |#eta_{BC}| - |#eta_{SCseed}|", 60,-(6.28/360*3),(6.28/360*3)); 
  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4 = fs.make<TH2F>("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4","h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4;  EBm4  #phi_{BC} - #phi_{SCseed}; EBm4  |#eta_{BC}| - |#eta_{SCseed}|", 20,-(6.28/360*20),(6.28/360*20), 20,-(6.28/360*3),(6.28/360*3));

  h_maxCryInDomino_barl           = fs.make<TH1F>("h_maxCryInDomino_barl","max Cry In Domino (barrel); i#eta", 5,-2,3); 
  h_maxCryInDominoVsPhi_barl      = fs.make<TH2F>("h_maxCryInDominoVsPhi_barl","max Cry In Domino Vs i#phi (barrel); EB i#eta_{BC} - i#eta_{SCseed}; EB i#phi_{BC} - i#phi_{SCseed}", 35,-17,18,5,-2,3); 

  h_maxCryInLocMax_barlSymm       = fs.make<TH1F>("h_maxCryInLocalMax_barlSymm","max Cry In Local Max (EB); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlSymm  = fs.make<TH2F>("h_maxCryInLocalMaxVsPhi_barlSymm","max Cry In Local Max Vs i#phi (EB); EB i#phi_{BC} - i#phi_{SCseed}; EB  i#eta_{BC} - i#eta_{SCseed};", 35,-17,18,5,-2,3); 
  h_maxCryInLocMax_barlPLus       = fs.make<TH1F>("h_maxCryInLocalMax_barlPLus","max Cry In Local Max (EB+); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlPLus  = fs.make<TH2F>("h_maxCryInLocalMaxVsPhi_barlPlus","max Cry In Local Max Vs i#phi (EB+); EB i#phi_{BC} - i#phi_{SCseed}; EB+  i#eta_{BC} - i#eta_{SCseed}", 35,-17,18,5,-2,3); 
  h_maxCryInLocMax_barlMinus      = fs.make<TH1F>("h_maxCryInLocalMax_barlMinus","max Cry In Local Max (EB-); i#eta_{BC} - i#eta_{SCseed}", 5,-2,3); 
  h_maxCryInLocMaxVsPhi_barlMinus = fs.make<TH2F>("h_maxCryInLocalMaxVsPhi_barlMinus","max Cry In Local Max Vs i#phi (EB-); EB i#phi_{BC} - i#phi_{SCseed}; EB-  i#eta - i#eta_{SCseed}", 35,-17,18,5,-2,3); 




}

void SCwithPUhistos::FillSc(const reco::SuperClusterCollection::const_iterator sc1, HepMC::GenEvent::particle_const_iterator particle1,
			    const EcalRecHitCollection* ebRecHits, const EcalRecHitCollection* eeRecHits)
{
  std::cout << "Filling histograms - ciao " << std::endl;
  

  // initial photon true variables 
  float phi_true=(*particle1)->momentum().phi();
  float eta_true=(*particle1)->momentum().eta();
  float etaEcal_true = etaTransformation(eta_true, (*particle1)->production_vertex()->position().z()/10. );
  float et_true = (*particle1)->momentum().e()/cosh((*particle1)->momentum().eta());
  
  // choose which energy to use in SC
  float energyScRAW  = sc1->rawEnergy();
  float energyScCorr = sc1->energy();
  
  // this is a SC in the EB case
  if (fabs(sc1->eta())<1.4442 && energyScCorr/cosh(sc1->eta())>20.) {
    
    // perform SC-MCparticle matching
    float deltaPhi = sc1->phi()-phi_true;
    float phiWidth = sc1->phiWidth();          // covariance of cluster in phi
    float energySCtoCheck=0;
    
    float phiSize = getPhiSize(sc1,true);
    if ( deltaPhi > PI ) deltaPhi -= TWOPI;
    if ( deltaPhi < -PI) deltaPhi += TWOPI;
  
    h_scet_barl->Fill((energyScCorr/cosh(sc1->eta()))/et_true);   // using: h_scet_barl as counter of EB matched photons
    h_EoverEtrue_barl->Fill(energyScCorr/(*particle1)->momentum().e());    //        supercluster that matches a MC-truth particle 
    h_phiWidth->Fill(phiWidth);
    h_phiWidthVsE->Fill(phiWidth,(*particle1)->momentum().e());	  
    h_phiSize->Fill(phiSize);
    h_phiSizeVsE->Fill(phiSize,(*particle1)->momentum().e());
    h_phiSizeVsEt->Fill(phiSize, et_true);
  
    int   whereIsMaxInDomino[35];
    float whatIsMaxInDomino[35];
    for(int u=0; u<35; u++) {
      whereIsMaxInDomino[u]=-999;
      whatIsMaxInDomino[u] =-999;
    } 
  
    std::vector< std::pair<DetId, float> >    theHitsAndFractions =  sc1->hitsAndFractions();  
    // loop over all the components of this supercluster	    
    for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
	idsIt != theHitsAndFractions.end(); ++idsIt) 
      {
	float thePhiDistance = getPhiDistance( (*sc1->seed()).seed() , 
					       (*idsIt).first );
	float currEbEnergy    = (ebRecHits->find( (*idsIt).first ))->energy();
	energySCtoCheck       += currEbEnergy;
	h_phiShape_barl       -> Fill(thePhiDistance, currEbEnergy/energyScRAW);
	h_absPhiShape_barl    -> Fill(fabs(thePhiDistance), currEbEnergy/energyScRAW);
	h_phiShapeVsE_barl    -> Fill(thePhiDistance, energyScRAW , currEbEnergy/energyScRAW);   
	h_absPhiShapeVsE_barl -> Fill(fabs(thePhiDistance), energyScRAW , currEbEnergy/energyScRAW);   
      
	float theEtaDistance = getEtaDistance( (*sc1->seed()).seed() , 
					       (*idsIt).first );
	float theSeedEta = EBDetId ( (*sc1->seed()).seed().rawId() ).ieta();
	// treat EB+ and EB- separately, give the sign of ieta; and make a symmetrized plot
	h_etaShape_barl        -> Fill(theEtaDistance, currEbEnergy/energyScRAW);
	h_etaShape_barlSymm    -> Fill( signum(theSeedEta) *theEtaDistance, currEbEnergy/energyScRAW);
	h_etaPhiShape_barl     -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energyScRAW);
	h_etaPhiShape_barlSymm -> Fill( signum(theSeedEta) *theEtaDistance, thePhiDistance, currEbEnergy/energyScRAW);
	if ( theSeedEta>0 )  {
	  h_etaShape_barlPLus    -> Fill(theEtaDistance, currEbEnergy/energyScRAW);
	  h_etaPhiShape_barlPLus -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energyScRAW);
	}
	else{
	  h_etaShape_barlMinus    -> Fill(theEtaDistance, currEbEnergy/energyScRAW);
	  h_etaPhiShape_barlMinus -> Fill(theEtaDistance, thePhiDistance, currEbEnergy/energyScRAW);
	}
      
	// at each given phi [0... 34], heck where in eta max energy crystal occurs
	// phi=17 is the seed by construction, in order to allow +-17
	int theIntegerPhiDistance = ((int)thePhiDistance) +17;
	int theIntegerEtaDistance = ((int)theEtaDistance);
      
	if (currEbEnergy > whatIsMaxInDomino[theIntegerPhiDistance]){
	  whatIsMaxInDomino[theIntegerPhiDistance]  = currEbEnergy;
	  whereIsMaxInDomino[theIntegerPhiDistance] = theIntegerEtaDistance;
	
	}
      } // loop over crystals of the SC
    
    for(int u=0; u<35; u++) {
      if (whereIsMaxInDomino[u]<-100) continue;
      h_maxCryInDomino_barl        -> Fill(whereIsMaxInDomino[u]);
      h_maxCryInDominoVsPhi_barl   -> Fill(u-17, whereIsMaxInDomino[u]);
    }
  
    // looping on basic clusters within the Supercluster    for(reco::CaloCluster_iterator bcIt = sc1->clustersBegin(); bcIt!=sc1->clustersEnd(); bcIt++)
    for(reco::CaloCluster_iterator bcIt = sc1->clustersBegin(); bcIt!=sc1->clustersEnd(); bcIt++)
      {
	// fill the plots that follow only if the basic cluster is NOT the seed basic cluster
	if (   EBDetId(   (*bcIt)->seed().rawId())  
	       == EBDetId ( (*sc1->seed()).seed().rawId() )  
	       ) continue;
	
	EBDetId theBCSeed           = EBDetId(   (*bcIt)->seed().rawId());  // this is always SCseed? FIXme
	EBDetId theSCseed           = EBDetId(   ( *(sc1->seed()) ).seed().rawId()  ) ;
		
	int     etaOfMaxInsidBC     =-999; 
	int     phiOfMaxInsidBC     =-999; 
	float   energyOfMaxInsideBC =-888;
		
	std::vector< std::pair<DetId, float> >  theBCHitsAndFractions =  (*bcIt)->hitsAndFractions();
	for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theBCHitsAndFractions.begin(); 
	    idsIt != theBCHitsAndFractions.end(); ++idsIt) {
		  
	  float currCryEnergy    = (ebRecHits->find( (*idsIt).first ))->energy();
	  if (currCryEnergy > energyOfMaxInsideBC){
	    etaOfMaxInsidBC         = (   (int)getEtaDistance ( (*sc1->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
	    phiOfMaxInsidBC         = (   (int)getPhiDistance ( (*sc1->seed()).seed() , EBDetId( (*idsIt).first.rawId() ) )  );
	    energyOfMaxInsideBC     = currCryEnergy;
	  } // if energy surpassed
	} // loop over theBCHitsAndFractions of BC
	      
	// fill the histograms with the location of the BC maximum, w.r.t. to SC seed
	if(fabs(theSCseed.ieta()) >2)
	  {
	    h_maxCryInLocMax_barlSymm           ->Fill( etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
	    h_maxCryInLocMaxVsPhi_barlSymm      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC * signum( theSCseed.ieta() ) );
	  }
	if(theSCseed.ieta() >2)
	  {
	    h_maxCryInLocMax_barlPLus           ->Fill( etaOfMaxInsidBC );
	    h_maxCryInLocMaxVsPhi_barlPLus      ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
	  } else if (theSCseed.ieta() <-2)
	  {
	    h_maxCryInLocMax_barlMinus          ->Fill( etaOfMaxInsidBC );
	    h_maxCryInLocMaxVsPhi_barlMinus     ->Fill( phiOfMaxInsidBC , etaOfMaxInsidBC );
	  }


	h_phiBCminusPhiSeed_barl  -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() );
	h_etaBCminusEtaSeed_barl  -> Fill( (*bcIt)->eta() - (*sc1->seed()).eta() );
	h_absEtaBCminusAbsEtaSeed_barl  -> Fill( fabs((*bcIt)->eta())   -  fabs((*sc1->seed()).eta())  );
	h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl   -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta()) );
	if ( fabs(theSCseed.ieta()) <= 25 ){
	  h_phiBCminusPhiSeed_barlm1  -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() ); h_etaBCminusEtaSeed_barlm1  -> Fill( (*bcIt)->eta() - (*sc1->seed()).eta() );
	  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1   -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta()) );
	  h_absEtaBCminusAbsEtaSeed_barlm1  -> Fill( fabs((*bcIt)->eta())  -  fabs((*sc1->seed()).eta())   );	} 
	else if	    ( 25<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 45 ){
	  h_phiBCminusPhiSeed_barlm2  -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() ); h_etaBCminusEtaSeed_barlm2  -> Fill( (*bcIt)->eta() - (*sc1->seed()).eta() );
	  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2   -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta()) );		
	  h_absEtaBCminusAbsEtaSeed_barlm2  -> Fill( fabs((*bcIt)->eta())  - fabs((*sc1->seed()).eta()) ); } 
	else if	    ( 45<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 65 ){
	  h_phiBCminusPhiSeed_barlm3  -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() ); h_etaBCminusEtaSeed_barlm3  -> Fill( (*bcIt)->eta() - (*sc1->seed()).eta() );
	  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3   -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta()) );		
	  h_absEtaBCminusAbsEtaSeed_barlm3  -> Fill( fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta())  );} 
	else if	    ( 65<fabs(theSCseed.ieta()) &&  fabs(theSCseed.ieta()) <= 85 ){
	  h_phiBCminusPhiSeed_barlm4  -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() ); h_etaBCminusEtaSeed_barlm4  -> Fill( (*bcIt)->eta() - (*sc1->seed()).eta() );
	  h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4   -> Fill( (*bcIt)->phi() - (*sc1->seed()).phi() , fabs((*bcIt)->eta()) - fabs((*sc1->seed()).eta()) );
	  h_absEtaBCminusAbsEtaSeed_barlm4  -> Fill( fabs((*bcIt)->eta())  -  fabs((*sc1->seed()).eta()) );} 
		
      }// loop over the BC's 
  }// end EB case


  

  // this is a SC in the EE case
  if (fabs(sc1->eta())>1.566 && fabs(sc1->eta())<2.5 && energyScCorr/cosh(sc1->eta())>20.) {
    

    float deltaPhi = sc1->phi()-phi_true;
    //float deltaEta = sc1->eta()-etaEcal_true;
    float phiWidth = sc1->phiWidth();
    float energyScRAW = sc1->rawEnergy();
    float energySCtoCheck=0;
    
    float phiSize = getPhiSize(sc1,false);
    
    if ( deltaPhi > pi ) deltaPhi -= twopi;
    if ( deltaPhi < -pi) deltaPhi += twopi;
	  
    h_scet_endc->Fill((energyScCorr/cosh(sc1->eta()))/et_true); // using: h_scet_barl as counter of EE matched photons
    h_EoverEtrue_endc->Fill(energyScCorr/(*particle1)->momentum().e());
    h_phiWidth_endc->Fill(phiWidth);
    h_phiWidthVsE_endc->Fill(phiWidth,(*particle1)->momentum().e());	  
    h_phiSize_endc->Fill(phiSize);
    h_phiSizeVsE_endc->Fill(phiSize,(*particle1)->momentum().e());	  
    h_phiSizeVsEt_endc->Fill(phiSize,et_true);	  
	  
    std::vector< std::pair<DetId, float> >    theHitsAndFractions =  sc1->hitsAndFractions();  
    for(std::vector< std::pair<DetId, float> >::const_iterator idsIt = theHitsAndFractions.begin(); 
	idsIt != theHitsAndFractions.end(); ++idsIt) 
      {
	float thePhiDistance = getPhiDistance( (*sc1->seed()).seed() , 
					       (*idsIt).first );
	float currEeEnergy   = (eeRecHits->find( (*idsIt).first ))->energy();
	energySCtoCheck      += currEeEnergy;
	h_phiShape_endc       -> Fill(thePhiDistance, currEeEnergy/energyScRAW);
	h_absPhiShape_endc    -> Fill(fabs(thePhiDistance), currEeEnergy/energyScRAW);
	h_phiShapeVsE_endc    -> Fill(thePhiDistance, energyScRAW , currEeEnergy/energyScRAW);   
	h_absPhiShapeVsE_endc -> Fill(fabs(thePhiDistance), energyScRAW, currEeEnergy/energyScRAW );     }// end EE case

  }// end ENDCAP case
}



void SCwithPUhistos::FillGamma(const reco::PhotonCollection::const_iterator pho1, HepMC::GenEvent::particle_const_iterator particle1){
  
  float energyPho    = pho1->energy();
  float energyPhoRAW = pho1->superCluster()->rawEnergy();

  if (pho1->isEB()) {
    h_PhoEoverEtrue_barl->Fill(energyPho/(*particle1)->momentum().e());
    h_E5x5overEtrue_barl->Fill(pho1->e5x5()/(*particle1)->momentum().e());
    h_E5x5overEtrueVsEphoEtrue_barl->Fill(pho1->e5x5()/(*particle1)->momentum().e(),energyPho/(*particle1)->momentum().e());
    if (pho1->r9()>0.94) { // r9 cut is 0.94 for EB and 0.95 for EE  
      h_E5x5R9overEtrue_barl->Fill(pho1->e5x5()/(*particle1)->momentum().e());
      h_PhoER9overEtrue_barl->Fill(energyPho/(*particle1)->momentum().e()); // counter r9 EB
    } else {
      h_E5x5notR9overEtrue_barl->Fill(pho1->e5x5()/(*particle1)->momentum().e());
      h_PhoEnotR9overEtrue_barl->Fill(energyPho/(*particle1)->momentum().e()); // counter !r9 EB
    }
    //isEB++;
  }// if photon in EB
  else if (pho1->isEE()) {
    h_PhoEoverEtrue_endc->Fill(energyPho/(*particle1)->momentum().e());
    h_E5x5overEtrue_endc->Fill(pho1->e5x5()/(*particle1)->momentum().e());
    h_E5x5overEtrueVsEphoEtrue_endc->Fill(pho1->e5x5()/(*particle1)->momentum().e(),energyPho/(*particle1)->momentum().e());
    if (pho1->r9()>0.94) { // r9 cut is 0.94 for EB and 0.95 for EE  
      h_E5x5R9overEtrue_endc->Fill(pho1->e5x5()/(*particle1)->momentum().e());
      h_PhoER9overEtrue_endc->Fill(energyPho/(*particle1)->momentum().e());
    } else {
      h_E5x5notR9overEtrue_endc->Fill(pho1->e5x5()/(*particle1)->momentum().e());
      h_PhoEnotR9overEtrue_endc->Fill(energyPho/(*particle1)->momentum().e());
    }
  }// if photon in EE
}





void SCwithPUhistos::FillH(const std::vector<reco::Photon> higgsPhotons, const std::vector<reco::Photon> higgsPhotons_trueVtx)
{
  std::cout << "H Filling histograms - ciao " << std::endl;

  int isEB(0);
  for(std::vector< reco::Photon>::const_iterator phoIt = higgsPhotons.begin();
      phoIt !=  higgsPhotons.end(); phoIt++ ) {
    if (phoIt->isEB()) isEB++;
  }
  

  if (higgsPhotons.size()>1) {
    if ((higgsPhotons[0].et()>40. && higgsPhotons[1].et()>30.) ||
	(higgsPhotons[0].et()>30. && higgsPhotons[1].et()>40.)) { // get two highest energy candidates
      double mHiggs = (higgsPhotons[0].p4()+higgsPhotons[1].p4()).M();
      if (isEB==2) {
	h_mHiggs_EBEB->Fill(mHiggs);
      } else if (isEB==1) {
	h_mHiggs_EBEE->Fill(mHiggs);
      } else {
	h_mHiggs_EEEE->Fill(mHiggs);
      }
    }
  }

  if (higgsPhotons_trueVtx.size()>1) {
    if ((higgsPhotons_trueVtx[0].et()>40. && higgsPhotons_trueVtx[1].et()>30.) ||
	(higgsPhotons_trueVtx[0].et()>30. && higgsPhotons_trueVtx[1].et()>40.)) {
      double mHiggs_trueVtx = (higgsPhotons_trueVtx[0].p4()+higgsPhotons_trueVtx[1].p4()).M();
      if (isEB==2) {
	h_mHiggs_EBEB_trueVtx->Fill(mHiggs_trueVtx);
      } else if (isEB==1) {
	h_mHiggs_EBEE_trueVtx->Fill(mHiggs_trueVtx);
      } else {
	h_mHiggs_EEEE_trueVtx->Fill(mHiggs_trueVtx);
      }
    }// if ET selections
  }// if there's a true vertex

}




// compute eta at detector, given physics eta and location of primary vertex
float SCwithPUhistos::etaTransformation(  float EtaParticle , float Zvertex)  {

  //---Definitions for ECAL
  // const float R_ECAL           = 136.5;  // radius of maximum containement
  // const float Z_Endcap         = 328.0;
  // const float etaBarrelEndcap  = 1.479;

  //---ETA correction

  float Theta = 0.0  ;
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0)    Theta = Theta+PI ;
  float ETA =      - log(tan(0.5*Theta));

  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle);
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+PI ;
      ETA = - log(tan(0.5*Theta));
    }
  //---Return the result
  return ETA;
  //---end
}
