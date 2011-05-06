#ifndef ZRAPIDITY_STANDARD_INCLUDED
#define ZRAPIDITY_STANDARD_INCLUDED

#include "TText.h"
#include "TLatex.h"
#include "TDatime.h"

static bool dotime=false;
static const char* y_xaxis_label="y_{ee}";
static const char* qt_xaxis_label="p_{T,ee} [GeV/c]";
static const char* qt_xaxis_ll_label="p_{T} [GeV/c]";

TH1* zpt_rebinForPlot(TH1* original) {
  char name[125];
  sprintf(name,"%s_4plot",original->GetName());
  int    pt_bins=19;
  double pt_binning[]={ 0.7, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
  TH1* rebinned=new TH1F(name,name,pt_bins,pt_binning);

  for(int i=1;i<=pt_bins;i++){
    rebinned->SetBinContent(i,original->GetBinContent(i));
    rebinned->SetBinError(i,original->GetBinError(i));
  }
  return rebinned;
}

void zrap_colors() {
  const Int_t NRGBs = 2;
  const Int_t NCont = 30;
  
  Double_t stops[NRGBs] = { 0.00, 1.0 };
  Double_t red[NRGBs]   = { 1.00, 0.0 };
  Double_t green[NRGBs] = { 1.00, 0.0 };
  Double_t blue[NRGBs]  = { 1.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}


void zrap_Prelim(double px, double py, double tx=-1, double ty=-1) {
  const char* time_;
  TDatime mytime;
  time_ = mytime.AsString();

  TText *plabel = new TText();
  plabel-> SetNDC();
  plabel -> SetTextFont(42);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.040);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);

  
  //Then for each plot, pick a nice spot and draw
  plabel -> DrawText(px, py, "CMS preliminary"); // 'ala VBTF'
  if (tx>=0 && ty>=0 && dotime) {
    TText *tlabel = new TText();
    tlabel-> SetNDC();
    tlabel -> SetTextFont(42);
    tlabel -> SetTextColor(1);
    tlabel -> SetTextSize(0.02);
    tlabel -> SetTextAlign(22);
    tlabel -> SetTextAngle(0);
    
    tlabel -> DrawText(tx, ty, Form("%s",time_));
  }
}

void zrap_Lumi(double px, double py, int val) {
  char text[40];
  sprintf(text,"%d pb^{-1} at #sqrt{s}=7 TeV",val);

  TLatex *plabel = new TLatex(px,py,text);
  plabel-> SetNDC();
  plabel -> SetTextFont(42);
  plabel -> SetTextColor(1);
  plabel -> SetTextSize(0.035);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);


  //Then for each plot, pick a nice spot and draw
  plabel -> Draw();

}

void zqt_Cut(double px, double py, const char* lept="l",double sizef=1.0) {
  char text[40];
  sprintf(text,"|#eta_{%s}|<2.1, p_{T,%s}>20 GeV/c",lept,lept);

  TLatex *plabel = new TLatex(px,py,text);
  plabel-> SetNDC();
  plabel -> SetTextFont(42);
  plabel -> SetTextColor(kBlue);
  plabel -> SetTextSize(0.045*sizef);
  plabel -> SetTextAlign(22);
  plabel -> SetTextAngle(0);


  //Then for each plot, pick a nice spot and draw
  plabel -> Draw();

}
#endif // ZRAPIDITY_STANDARD_INCLUDED
