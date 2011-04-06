#include "readStandardFile.C"

double anXSec(TH1* data, TH1* bkgd, TH1* numer, TH1* denom) {
  double ndata=data->Integral();
  double nbkgd=bkgd->Integral();
  double nnumer=numer->Integral();
  double ndenom=denom->Integral();

  double effacc=(nnumer/ndenom);

  printf("%.0f %.1f %.3f\n",ndata,nbkgd,effacc);

  return (ndata-nbkgd)/effacc;


}

double anError(TH1* data, TH1* bkgd) {
  double ndata=data->Integral();
  double nbkgd=bkgd->Integral();
  return 1.0/sqrt(ndata-nbkgd);
}

void estimateXSecs(TFile* effacc, double lumi=36.0) {

  TH1* y_denom=effacc->Get("mcEff/All/Z0_YTL");

  TH1* ee_y=readStandardFile("ee_y","ECAL80-ECAL95_RAW.csv");
  TH1* ee_y_bkg=readStandardFile("ee_y_bkg","background_all_y_ee.csv");
  TH1* ee_y_numer=effacc->Get("mcEff/ECAL80-ECAL95/C07-HLT-GSF/Z0_Y_masscut");

  double ee_y_xsec=1.0/lumi*anXSec(ee_y,ee_y_bkg,ee_y_numer,y_denom);

  TH1* ef_y=readStandardFile("ef_y","ECAL80-HF_RAW.csv");
  TH1* ef_y_bkg=readStandardFile("ef_y_bkg","background_all_y_ef.csv");
  TH1* ef_numer=effacc->Get("mcEff/ECAL80-HF/C07-HLT-GSF/Z0_Y_masscut");
  
  double ef_xsec=1.0/lumi*anXSec(ef_y,ef_y_bkg,ef_numer,y_denom);


  //  TH1* pt_denom=effacc->Get("mcEff/MUON-MUON/Acceptance/Z0_PtTL");
  TH1* pt_denom=effacc->Get("mcEff/All/Z0_PtTL");

  TH1* ee_pt=readStandardFilePt("ee_pt","ECAL80-ECAL95_Pt_RAW.csv");
  TH1* ee_pt_bkg=readStandardFilePt("ee_pt_bkg","background_all_qt.csv");
  TH1* ee_pt_numer=effacc->Get("mcEff/ECAL80-ECAL95-MUO/C07-HLT-GSF/Z0_Pt_masscut");

  double ee_pt_xsec=1.0/lumi*anXSec(ee_pt,ee_pt_bkg,ee_pt_numer,pt_denom);

  double ee_y_error=anError(ee_y,ee_y_bkg);
  double ef_y_error=sqrt(pow(anError(ef_y,ef_y_bkg),2)+pow(0.04,2));
  double ee_pt_error=anError(ee_pt,ee_pt_bkg);

  printf("ECAL80-ECAL95 Y analysis : %.3f nb +/- %.3f\n",ee_y_xsec/1000,ee_y_error*ee_y_xsec/1000);
  printf("ECAL80-HF Y analysis : %.3f nb +/- %.3f\n",ef_xsec/1000,ef_y_error*ef_xsec/1000);
  printf("ECAL80-ECAL95-MUON Pt analysis : %.3f nb +/- %.3f\n",ee_pt_xsec/1000,ee_pt_error*ee_pt_xsec/1000);
}
