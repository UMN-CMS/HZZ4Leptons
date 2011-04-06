#include "readStandardFile.C"
#include "zrapidityStandard.C"
#include "tdrstyle.C"

void finalBackground(int mode) {
  setTDRStyle();
  int nbins;
  FILE *outf;
  TH1* qcd,*ttbar,*tau, *vv;

  switch (mode) {
  case (1) : 
    nbins=100;
    qcd=readStandardFile("qcd","background_qcd_y_ee.csv");
    ttbar=readStandardFile("ttbar","background_ttbar_y_ee.csv");
    tau=readStandardFile("tau","background_tautau_y_ee.csv");
    vv=readStandardFile("vv","background_vv_y_ee.csv");
    outf=fopen("background_all_y_ee.csv","wt");
    break;
  case (2) : 
    nbins=100;
    qcd=readStandardFile("qcd","background_qcd_y_ef.csv");
    ttbar=readStandardFile("ttbar","background_ttbar_y_ef.csv");
    tau=readStandardFile("tau","background_tautau_y_ef.csv");
    vv=readStandardFile("vv","background_vv_y_ef.csv");
    outf=fopen("background_all_y_ef.csv","wt");
    break;
  case (3) : 
    nbins=18;
    qcd=readStandardFilePt("qcd","background_qcd_qt.csv");
    ttbar=readStandardFilePt("ttbar","background_ttbar_qt.csv");
    tau=readStandardFilePt("tau","background_tautau_qt.csv");
    vv=readStandardFilePt("vv","background_vv_qt.csv");

    outf=fopen("background_all_qt.csv","wt");
    break;
  }

  for (int i=1; i<=nbins; i++) {
    fprintf(outf,"%3d %8.1f %8.1f ",i,qcd->GetXaxis()->GetBinLowEdge(i),
	    qcd->GetXaxis()->GetBinUpEdge(i));

    double tback=qcd->GetBinContent(i)+ttbar->GetBinContent(i)+tau->GetBinContent(i)+vv->GetBinContent(i);
    double terr=sqrt(pow(qcd->GetBinError(i),2)+pow(ttbar->GetBinError(i),2)+pow(tau->GetBinError(i),2)+pow(vv->GetBinError(i),2));

    fprintf(outf,"%8.2f %8.2f\n",tback,terr);
  }

  fclose(outf);

  if (mode==1 || mode==2) {
    qcd->GetXaxis()->SetRangeUser(-3.5,3.5);
    qcd->SetMinimum(0.01);
    if (mode==2)     qcd->SetMaximum(50);
    if (mode==1)     qcd->SetMaximum(100);
    qcd->GetXaxis()->SetTitle(y_xaxis_label);
    qcd->GetXaxis()->CenterTitle();
    qcd->GetYaxis()->SetTitle("Events/0.1 Units of Rapidity");
  }
  if (mode==3) {
    qcd->GetXaxis()->SetRangeUser(0.5,500);
    qcd->SetMinimum(0.005);
    qcd->SetMaximum(200);

    qcd->GetXaxis()->SetTitle(qt_xaxis_label);
    qcd->GetXaxis()->CenterTitle();
    qcd->GetYaxis()->SetTitle("Events/Bin");
  }

  qcd->SetMarkerStyle(20);
  ttbar->SetMarkerStyle(24);
  tau->SetMarkerStyle(26);
vv->SetMarkerStyle(30);
  TCanvas* c1=new TCanvas("c1","c1",700,700);
  c1->SetLogy();

  qcd->Draw("E1");
  ttbar->Draw("E1 SAME");
  tau->Draw("E1 SAME");
vv->Draw("E1 SAME");

  if (mode==1) {
    TLegend* tl=new TLegend(0.20,0.80,0.51,0.92,"ECAL-ECAL Channel");
    tl->AddEntry(qcd,"QCD Background ","P");
    tl->AddEntry(ttbar,"t#bar{t} Background","P");
    tl->AddEntry(tau,"Z/#gamma*#rightarrow#tau#tau Background","P");
    tl->AddEntry(vv,"VV Background","P");

    tl->Draw();
    zrap_Prelim(0.75,0.9,0.8,0.8);
    zrap_Lumi(0.86,0.86,36);
    c1->Print("background_breakdown_ee.eps");
  }
  if (mode==2) {
    TLegend* tl=new TLegend(0.35,0.75,0.75,0.85,"ECAL-HF Channel");
    tl->AddEntry(qcd,"QCD Background ","P");
    tl->AddEntry(ttbar,"t#bar{t} Background","P");
    tl->AddEntry(tau,"Z/#gamma*#rightarrow#tau#tau Background","P");
    tl->AddEntry(vv,"VV Background","P");

    tl->Draw();
    zrap_Prelim(0.75,0.9,0.8,0.8);
    zrap_Lumi(0.86,0.86,36);
    c1->Print("background_breakdown_ef.eps");
  }

  if (mode==3) {
    c1->SetLogx();
    c1->SetBottomMargin(0.16);
    qcd->GetXaxis()->SetTitleOffset(1.1);
    TLegend* tl=new TLegend(0.20,0.81,0.51,0.90);
    tl->AddEntry(qcd,"QCD Background ","P");
    tl->AddEntry(ttbar,"t#bar{t} Background","P");
    tl->AddEntry(tau,"Z/#gamma*#rightarrow#tau#tau Background","P");
    tl->AddEntry(vv,"VV Background","P");

    tl->Draw();
    zrap_Prelim(0.75,0.9,0.8,0.8);
    zrap_Lumi(0.86,0.86,36);
    c1->Print("background_breakdown_qt.eps");
  }

}
