void extractData(TFile* f0) {

  TH1* hqt=f0->Get("ZFromData/ECAL80-ECAL95-MUO/C08-m(60,120)/Z0_Pt_masscut")->Clone("hqt");

  FILE* record;
  record=fopen("data_qt.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=18; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hqt->GetXaxis()->GetBinLowEdge(i), hqt->GetXaxis()->GetBinUpEdge(i));
    double terr=sqrt(hqt->GetBinContent(i));
    fprintf(record,"%9.2f %9.2f\n",hqt->GetBinContent(i),terr);
  }

  fclose(record);
  hqt->Draw();

  TH1* hye=f0->Get("ZFromData/ECAL80-ECAL95/C08-m(60,120)/Z0_Y_masscut")->Clone("hye");

  FILE* record;
  record=fopen("data_y_ee.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=100; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hye->GetXaxis()->GetBinLowEdge(i), hye->GetXaxis()->GetBinUpEdge(i));
    double terr=sqrt(hye->GetBinContent(i));
    fprintf(record,"%9.2f %9.2f\n",hye->GetBinContent(i),terr);
  }

  fclose(record);
  hye->Draw();

  TH1* hyf=f0->Get("ZFromData/ECAL80-HF/C08-m(60,120)/Z0_Y_masscut")->Clone("hyf");

  FILE* record;
  record=fopen("data_y_ef.csv","wt");
  fprintf(record,"#    ptBin     ptMin     ptMax   Value     Error\n");

  for (int i=1; i<=100; i++) {
    fprintf(record,"%8d %9.1f %9.1f",i,hyf->GetXaxis()->GetBinLowEdge(i), hyf->GetXaxis()->GetBinUpEdge(i));
    // total error includes 11% luminosity uncertainty
    double terr=sqrt(hyf->GetBinContent(i));
    fprintf(record,"%9.2f %9.2f\n",hyf->GetBinContent(i),terr);
  }

  fclose(record);
  hyf->Draw();

  



}
