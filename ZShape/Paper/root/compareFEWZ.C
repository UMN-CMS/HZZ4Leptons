#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

int colorForId(int i) {
  const int colors[] = {kCyan,kRed, kBlue, kGreen+2,
			kBlack, kGray
  };
  return colors[i];
}


void compareFEWZ(const char* models,const char* outgoing=0) {
  setTDRStyle();
  bool limiNormal=false;
  

  const char* knownModels[] = { 
    "D6T","Pythia Tune D6T","mcFile_D6T_TL.txt",
    "P0","Pythia P0","mcFile_P0_TL.txt",
    "ProPT0","Pythia ProPT0","mcFile_ProPT0_TL.txt",
    "ProQ20","Pythia ProQ20","mcFile_ProQ20_TL.txt",
    "Z2","Pythia Tune Z2 ","mcFile_Z2_TL.txt",
    "FEWZ-1","FEWZ NLO","fewz_v1.txt",
    "Powheg","Powheg","powheg_z2_scaleUnc.txt",
    "Resbos-P","Resbos-P","ptShape_Resbos_a0_80.txt",
    "Resbos-CP","ResBos","resbos-cp.txt",
    "R0.50","Resbos #alpha=0.50","ptShape_Resbos_a0_50.txt",
    "R0.70","Resbos #alpha=0.70","ptShape_Resbos_a0_70.txt",
    "R0.75","Resbos #alpha=0.75","ptShape_Resbos_a0_75.txt",
    "R0.80","Resbos #alpha=0.80","ptShape_Resbos_a0_80.txt",
    "R0.85","Resbos #alpha=0.85","ptShape_Resbos_a0_85.txt",
    "R0.90","Resbos #alpha=0.90","ptShape_Resbos_a0_90.txt",
    "R0.95","Resbos #alpha=0.95","ptShape_Resbos_a0_95.txt",
    "R1.00","Resbos #alpha=1.00","ptShape_Resbos_a1_0.txt",
    "R1.25","Resbos #alpha=1.25","ptShape_Resbos_a1_25.txt",
    "R1.50","Resbos #alpha=1.50","ptShape_Resbos_a1_5.txt",
    "R1.75","Resbos #alpha=1.75","ptShape_Resbos_a1_75.txt",
    "R2.00","Resbos #alpha=2.00","ptShape_Resbos_a2_0.txt",
    "R2.25","Resbos #alpha=2.25","ptShape_Resbos_a2_25.txt",
    "R2.50","Resbos #alpha=2.50","ptShape_Resbos_a2_5.txt",
    "R2.75","Resbos #alpha=2.75","ptShape_Resbos_a2_75.txt",
    "R3.00","Resbos #alpha=3.00","ptShape_Resbos_a3_0.txt",
    "R3.25","Resbos #alpha=3.25","ptShape_Resbos_a3_25.txt",
    "R3.50","Resbos #alpha=3.50","ptShape_Resbos_a3_5.txt",
    "R3.75","Resbos #alpha=3.75","ptShape_Resbos_a3_75.txt",
    "R4.00","Resbos #alpha=4.00","ptShape_Resbos_a4_0.txt",
    "PWG_D6T","Powheg + Pythia Tune D6T","PtTL_pwg_PythiaUED6T.txt",
    "PWG_P0","Powheg + Pythia P0","PtTL_pwg_PythiaUEP0.txt",
    "PWG_ProPT0","POWHEG + Pythia ProPT0","PtTL_pwg_PythiaUEProPT0.txt",
    "PWG_ProQ20","Powheg + Pythia ProQ20","PtTL_pwg_PythiaUEProQ20.txt",
    "PWG_Z2","POWHEG + Pythia Tune Z2","PtTL_pwg_PythiaUEZ2.txt",
    0,0,0
    };
  int nhists=0;
  TH1* hists[10],*histsDelta[10];
  int imodels[10];
  char names[10][50];
  double chi2[10];
  
  // find the models
  int nmatch=sscanf(models,"%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,],%[^,]",
		    names[0],
		    names[1],
		    names[2],
		    names[3],
		    names[4],
		    names[5],
  		    names[6],
		    names[7],
		    names[8],
		    names[9]);
  for (int i=0; i<nmatch; i++) {
    int jmodel=0;
    for (jmodel=0; knownModels[jmodel*3]!=0; jmodel++)
      if (!strcmp(knownModels[jmodel*3],names[i])) break;
    if (knownModels[jmodel*3]!=0) {
      imodels[nhists]=jmodel;
      nhists++;
    } else {
      printf("Didn't know: '%s'\n",names[i]);
    }
  }


  
  TCanvas* c1=new TCanvas("c1","c1",750,1000);
  c1->SetLogx()
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.20);
  TPad* p1=new TPad("p1","p1",0,0.38,1,0.98,0,0);
  p1->SetNumber(1);
  p1->SetFillStyle(0);
  p1->SetLogx();
  p1->SetLogy();
  p1->Draw();
  p1->SetBottomMargin(0.02);
  p1->SetRightMargin(0.02);
  TPad* p2=new TPad("p2","p2",0,00,1,0.38,0,0);
  p2->SetBottomMargin(0.25);
  p2->SetLogx();
  p2->SetNumber(2);
  p2->SetRightMargin(0.02);
  p2->Draw();
  p1->cd();

  TH1* dummya=new TH1F("dummy0","",20,20.0,600.0);

  TH1* data=readStandardFilePt("data","../data/dsdqt_combined.csv");
  dummya->SetMaximum(0.08);
  dummya->SetMinimum(1e-6);
  dummya->GetXaxis()->SetRangeUser(0.7,599.0);
  dummya->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dummya->GetXaxis()->SetLabelSize(0.0);
  dummya->GetXaxis()->SetTitleOffset(1.0);
  dummya->GetYaxis()->SetTitleOffset(1.0);
  dummya->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T} [(GeV/c)^{-1} ]");
  dummya->GetXaxis()->CenterTitle();
  dummya->Draw("HIST");

  char* Leg="";
  //  if (limiNormal) Leg="MC Normalized over 0 < q_{T} < 30";

  TLegend* tl=new TLegend(0.16,0.18,0.63,0.35,Leg);
  tl->SetFillColor(kWhite);
  tl->SetBorderSize(0);
  TLegend* tl2=new TLegend(0.23,0.77,0.95,0.92,Leg);
  tl->AddEntry(data,"data (e + #mu combined)","p");
  tl2->SetFillColor(kWhite);
  tl2->SetBorderSize(0);

  //  tl2->AddEntry(data,"Combined e + #mu","p");

  int ndof=-1;
  
  
  for (int i=0; i<nhists; i++) {
    char fullfname[1024];
    sprintf(fullfname,"../theory/%s",knownModels[imodels[i]*3+2]);
    printf("%s\n",fullfname);
    hists[i]=readStandardFilePt(knownModels[imodels[i]*3],fullfname);

    if (strstr(fullfname,"Resbos")) {
      double sum=hists[i]->Integral();
      for (int j=1; j<=18; j++) 
	hists[i]->SetBinContent(j,hists[i]->GetBinContent(j)/hists[i]->GetBinWidth(j)/sum);
    }

    if (strstr(fullfname,"fewz") || limiNormal) {
      double sum=data->Integral(9,18)/hists[i]->Integral(9,18);
      for (int j=9; j<=18; j++)  {
	hists[i]->SetBinContent(j,hists[i]->GetBinContent(j)*sum);
	hists[i]->SetBinError(j,hists[i]->GetBinError(j)*sum);
      }
      /*
      hists[i]->SetBinContent(8,hists[i]->GetBinContent(9));
      hists[i]->SetBinError(8,hists[i]->GetBinError(9));
      */
    }



    sprintf(fullfname,"%s-delta",knownModels[imodels[i]*3+2]);
    histsDelta[i]=(TH1*)(hists[i]->Clone(fullfname));

    chi2[i]=0;
    
    for (int j=1; j<=18; j++) {
      if (hists[i]->GetBinContent(j)==0) continue;
      /*
      histsDelta[i]->SetBinContent(j,data->GetBinContent(j)/hists[i]->GetBinContent(j)-1);
      histsDelta[i]->SetBinError(j,data->GetBinError(j)/data->GetBinContent(j));
      */
      /*
      histsDelta[i]->SetBinContent(j,(data->GetBinContent(j)-hists[i]->GetBinContent(j))/data->GetBinError(j));
      histsDelta[i]->SetBinError(j,hists[i]->GetBinError(j)/data->GetBinError(j));
      */
      histsDelta[i]->SetBinContent(j,(data->GetBinContent(j)-hists[i]->GetBinContent(j))/hists[i]->GetBinError(j));
      printf("%d %f %f %f %f\n",j,data->GetBinContent(j),hists[i]->GetBinContent(j),hists[i]->GetBinError(j),histsDelta[i]->GetBinContent(j));

      histsDelta[i]->SetBinError(j,data->GetBinError(j)/hists[i]->GetBinError(j));

      if (j>=9 && j<=18) {
	double binerror2=pow(data->GetBinError(j),2) + pow(hists[i]->GetBinError(j),2);
	chi2[i]+=pow(data->GetBinContent(j)-hists[i]->GetBinContent(j),2)/binerror2;

	if (i==0) ndof++;
      }
    }

    hists[i]=zpt_rebinForPlot(hists[i]);
    histsDelta[i]=zpt_rebinForPlot(histsDelta[i]);
    
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colorForId(i));
    hists[i]->SetFillColor(colorForId(i));

    histsDelta[i]->SetLineWidth(2);
    histsDelta[i]->SetLineColor(colorForId(i));

    tl->AddEntry(hists[i],knownModels[imodels[i]*3+1],"LF");
    tl2->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");

    hists[i]->GetXaxis()->SetRangeUser(20,600);
    hists[i]->SetMarkerStyle(1);
    hists[i]->SetMarkerSize(0.0001);
    hists[i]->SetFillStyle(1001);
    hists[i]->Draw("E3 SAME");
    
  }
  tl->Draw();
  c1->SetLogx();
  c1->SetLogy();
  data->Draw("E1 SAME");

  c1->cd();
  c1->Draw();

  zrap_Prelim(0.85,0.975,0.4,0.17);
  zrap_Lumi(0.81,0.91,36);
  zqt_Cut(0.33,0.43,"l",0.75);

  p2->cd();

  //  TCanvas* c2=new TCanvas("c2","c2",800,400);
  //  TH2* dummy=new TH2F("dummy","dummy",20,0.0,30.0,30,-0.35,0.35);
  TH2* dummy=new TH2F("dummy","",20,20,600.0,30,-2.5,2.5);
  dummy->Draw();

  dummy->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dummy->GetYaxis()->SetTitleOffset(0.5);
  dummy->GetYaxis()->SetTitle("(data - theory)/#sigma_{theory}");
  dummy->GetYaxis()->SetTitleSize(0.09);
  dummy->GetXaxis()->SetTitleSize(0.12);
  dummy->GetXaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetLabelSize(0.08);
  dummy->GetYaxis()->SetNdivisions(503);
  dummy->GetXaxis()->SetNdivisions(310);


  dummy->GetXaxis()->SetNoExponent(true);
  dummy->GetXaxis()->SetMoreLogLabels(true);
  
  

  /*
  static const int magicn=17;
  double xv[(magicn+2)*2],yv[(magicn+2)*2];

  for (int i=1; i<=magicn; i++) {
    xv[i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[i]=histsDelta[0]->GetBinError(i);
  }
  xv[0]=0; yv[0]=yv[1];
  xv[magicn+1]=600; yv[magicn+1]=yv[magicn];
  for (int i=1; i<=magicn; i++) {
    xv[(magicn+2)*2-1-i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[(magicn+2)*2-1-i]=-histsDelta[0]->GetBinError(i);
  }
  xv[magicn+2]=xv[magicn+1]; yv[magicn+2]=yv[magicn+3];
  xv[(magicn+2)*2-1]=0; yv[(magicn+2)*2-1]=yv[(magicn+2)*2-2];

  TGraph* deltaBand=new TGraph((magicn+2)*2,xv,yv);
  */
  const double x[4]={0,600,600,0};
  const double y1[4]={1,1,-1,-1};
  const double y2[4]={2,2,-2,-2};
  TGraph* band1s=new TGraph(4,x,y1);
  TGraph* band2s=new TGraph(4,x,y2);
  band1s->SetFillColor(kGreen-7);
  band2s->SetFillColor(kYellow);
  // band2s->Draw("F SAME");
    band1s->Draw("F SAME");


  for (int i=0; i<nhists; i++) {
    histsDelta[i]->SetLineColor(kBlack);
    histsDelta[i]->Draw("SAME ][ E1");
  }
  /*
  tl2->SetNColumns(2);

  tl2->Draw();
  */
  //  zrap_Prelim(0.18,0.97,0.4,0.17);
  //  zrap_Lumi(0.88,0.97,36);


  char outx[1024];
  if (outgoing!=0) {
    sprintf(outx,"%s.png",outgoing);
    //    c1->Print(outx);
    sprintf(outx,"%s.eps",outgoing);
    c1->Print(outx);
    sprintf(outx,"%s.C",outgoing);
    c1->Print(outx);
    sprintf(outx,"%s_delta.png",outgoing);
    //c2->Print(outx);
    sprintf(outx,"%s_delta.eps",outgoing);
    c2->Print(outx);
    sprintf(outx,"%s_delta.C",outgoing);
    c2->Print(outx);
  }
  for (int i=0; i<nhists; i++) {
    printf(" Comparison %s : %f/%d %f \n",knownModels[imodels[i]*3+1],chi2[i],ndof,TMath::Prob(chi2[i],ndof));

  }
}
