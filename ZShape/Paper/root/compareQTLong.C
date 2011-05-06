#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"
#include "../theory/database.C"

int colorForId(int i) {
  const int colors[] = {kRed, kBlue, kGreen+2,
			kBlack, kGray
  };
  return colors[i];
}


void compareQTLong(const char* models,const char* outgoing=0) {
  setTDRStyle();
  bool limiNormal=true;

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
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.22);
  TPad* p1=new TPad("p1","p1",0,0.38,1,0.98,0,0);
  p1->SetNumber(1);
  p1->SetFillStyle(0);
  p1->Draw();
  p1->SetLeftMargin(0.15);
  p1->SetBottomMargin(0.02);
  TPad* p2=new TPad("p2","p2",0,00,1,0.38,0,0);
  p2->SetBottomMargin(0.25);  
  p2->SetLeftMargin(0.15);
  p2->SetNumber(2);
  p2->Draw();
  p1->cd();
  
  TH1* dummya=new TH1F("dummy0","dummy0",20,0.7,600.0);

  TH1* data=readStandardFilePt("data","../data/dsdqt_combined.csv");
  dummya->SetMaximum(0.08);
  dummya->SetMinimum(1e-6);
  dummya->GetXaxis()->SetRangeUser(0.7,599.0);
  dummya->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dummya->GetXaxis()->SetTitleOffset(1.2);
  dummya->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T} [(GeV/c)^{-1}]");
  dummya->GetXaxis()->CenterTitle();
  dummya->Draw("HIST");

  char* Leg="";
  //  if (limiNormal) Leg="MC Normalized over 0 < q_{T} < 30";

  TLegend* tl=new TLegend(0.18,0.22,0.74,0.40,Leg);
  tl->SetFillStyle(0);
  TLegend* tl2=new TLegend(0.23,0.77,0.95,0.92,Leg);
  tl2->SetFillStyle(0);
  tl->AddEntry(data,"Combined e + #mu","p");
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

      double binerror=sqrt(pow(data->GetBinError(j),2)+pow(hists[i]->GetBinError(j),2));
      histsDelta[i]->SetBinContent(j,(data->GetBinContent(j)-hists[i]->GetBinContent(j))/binerror);
      histsDelta[i]->SetBinError(j,hists[i]->GetBinError(j)/data->GetBinError(j));

      histsDelta[i]->SetBinError(j,hists[i]->GetBinError(j)/data->GetBinError(j));
      if (j>=9 && j<=18) {
	chi2[i]+=pow(data->GetBinContent(j)-hists[i]->GetBinContent(j),2)/
	  pow(data->GetBinError(j),2);
	if (i==0) ndof++;
      }
    }

    hists[i]=zpt_rebinForPlot(hists[i]);
    histsDelta[i]=zpt_rebinForPlot(histsDelta[i]);
    
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colorForId(i));

    histsDelta[i]->SetLineWidth(2);
    histsDelta[i]->SetLineColor(colorForId(i));

    tl->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");
    tl2->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");


    hists[i]->Draw("HIST ][ SAME");
  }
  tl->SetNColumns(2);
  tl->Draw();
  p1->SetLogx();
  p1->SetLogy();
  data->Draw("E SAME");

  //  TH2* dummy=new TH2F("dummy","dummy",20,0.0,30.0,30,-0.35,0.35);
  p2->cd();
  TH2* dummy=new TH2F("dummy","dummy",20,0.7,600.0,30,-4.0,6.0);
  dummy->Draw();
  p2->SetLogx();

  dummy->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dummy->GetYaxis()->SetTitleOffset(0.5);
  dummy->GetYaxis()->SetTitle("(data-theory)/Error");
  dummy->GetYaxis()->SetTitleSize(0.09);


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
  band2s->Draw("F SAME");
  band1s->Draw("F SAME");

  for (int i=0; i<nhists; i++) {
    histsDelta[i]->Draw("SAME ][ LHISt");
  }

  tl2->SetNColumns(2);

  //  tl2->Draw();

  c1->cd();
  c1->Draw();
  zrap_Prelim(0.80,0.98,0.4,0.17);
  zrap_Lumi(0.75,0.90,36);
  zqt_Cut(0.34,0.43,"l",0.75);

  char outx[1024];
  if (outgoing!=0) {
    sprintf(outx,"%s.png",outgoing);
    //    c1->Print(outx);
    sprintf(outx,"%s.eps",outgoing);
    c1->Print(outx);
    //sprintf(outx,"%s_delta.png",outgoing);
    //p2->Print(outx);
    //sprintf(outx,"%s_delta.eps",outgoing);
    //c2->Print(outx);
  }
  for (int i=0; i<nhists; i++) {
    printf(" Comparison %s : %f/%d %f \n",knownModels[imodels[i]*3+1],chi2[i],ndof,TMath::Prob(chi2[i],ndof));

  }
}
