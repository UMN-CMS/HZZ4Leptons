#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"
#include "../theory/database.C"
#include <algorithm>

TMatrixD readCov(const char* covFile);

int colorForId(int i) {
  const int colors[] = {kCyan,kRed, kBlue, kGreen+2,
			kBlack, kGray
  };
  return colors[i];
}


void compareFEWZ(const char* models,const char* outgoing=0,bool limiNormal=true) {
  setTDRStyle();
  
  int nhists=0;
  TGraphAsymmErrors* hists[10],*histsDelta[10];
  int imodels[10];
  char names[10][50];
  double chi2[10];
  TMatrixD deltas[10];
  
  double terror[20];

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
  c1->SetLogx();
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

  TH1* dummya=new TH1F("dummy0","",20,19.5,600.0);

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

  for (int i=0; i<10; i++) {
    deltas[i].ResizeTo(1,10,1,1);
  }  
  
  for (int i=0; i<nhists; i++) {
    char fullfname[1024];
    sprintf(fullfname,"../theory/%s",knownModels[imodels[i]*3+2]);
    printf("%s\n",fullfname);
    const int iskip=8;
    hists[i]=readStandardFilePtTGE(fullfname,iskip);

    if (limiNormal) {
      double sum=0;
      for (int j=9-iskip-1; j<=18-iskip-1; j++)  
	sum+=hists[i]->GetY()[j];
      sum=data->Integral(9,18)/sum;

      for (int j=9-iskip-1; j<=18-iskip-1; j++)  {
	hists[i]->SetPoint(j,hists[i]->GetX()[j],sum*hists[i]->GetY()[j]);
	hists[i]->SetPointEYlow(j,hists[i]->GetEYlow()[j]*sum);
	hists[i]->SetPointEYhigh(j,hists[i]->GetEYhigh()[j]*sum);
      }
      
      //      hists[i]->SetBinContent(8,hists[i]->GetBinContent(9));
      //  hists[i]->SetBinError(8,hists[i]->GetBinError(9));
      
    }
   
    sprintf(fullfname,"%s-delta",knownModels[imodels[i]*3+2]);
    //    histsDelta[i]=(TH1*)(hists[i]->Clone(fullfname));

    chi2[i]=0;
    histsDelta[i]=new TGraphAsymmErrors(hists[i]->GetN());
        
    double px,py;
    for (int j=1; j<=18; j++) {
      int which=j-1-iskip;
      hists[i]->GetPoint(which,px,py);
      if (py==0) continue;
      double eyl=hists[i]->GetErrorYlow(which);
      double eyh=hists[i]->GetErrorYhigh(which);
      double exl=hists[i]->GetErrorXlow(which);
      double exh=hists[i]->GetErrorXhigh(which);
      
      /*
      histsDelta[i]->SetBinContent(j,data->GetBinContent(j)/hists[i]->GetBinContent(j)-1);
      histsDelta[i]->SetBinError(j,data->GetBinError(j)/data->GetBinContent(j));
      */
      /*
      histsDelta[i]->SetBinContent(j,(data->GetBinContent(j)-hists[i]->GetBinContent(j))/data->GetBinError(j));
      histsDelta[i]->SetBinError(j,hists[i]->GetBinError(j)/data->GetBinError(j));
      */

      if (data->GetBinContent(j)>py) {
	histsDelta[i]->SetPoint(which,px,(data->GetBinContent(j)-py)/(eyh));
	histsDelta[i]->SetPointError(which,exl,exh,data->GetBinError(j)/eyh,data->GetBinError(j)/eyh);
      } else {
	histsDelta[i]->SetPoint(which,px,(data->GetBinContent(j)-py)/(eyl));
	histsDelta[i]->SetPointError(which,exl,exh,data->GetBinError(j)/eyl,data->GetBinError(j)/eyl);
      }

      printf("%d %e %d %e %f\n",j,data->GetBinContent(j), j-iskip-1, py, histsDelta[i]->GetY()[which]);

      //      histsDelta[i]->SetBinError(j,data->GetBinError(j)/hists[i]->GetBinError(j));
      
      if (j>=9 && j<=18) {
	double binerror2=pow(data->GetBinError(j),2) + pow(std::max(fabs(eyh),fabs(eyl)),2);
	deltas[i](1,j-8)=(data->GetBinContent(j)-py);
	chi2[i]+=pow(data->GetBinContent(j)-py,2)/binerror2;
	if (i==0) terror[j]=std::max(fabs(eyh),fabs(eyl));

	if (i==0) ndof++;
      }
      
    }

    //    hists[i]=zpt_rebinForPlot(hists[i]);
    //    histsDelta[i]=zpt_rebinForPlot(histsDelta[i]);
    
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colorForId(i));
    hists[i]->SetFillColor(colorForId(i));
    
    histsDelta[i]->SetLineWidth(2);
    histsDelta[i]->SetLineColor(colorForId(i));
    
    tl->AddEntry(hists[i],knownModels[imodels[i]*3+1],"LF");
    tl2->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");

    histsDelta[i]->SetMarkerStyle(20+i);

    hists[i]->GetXaxis()->SetRangeUser(20,600);
    hists[i]->SetMarkerStyle(1);
    hists[i]->SetMarkerSize(0.0001);
    hists[i]->SetFillStyle(1001);
    //    hists[i]->Draw("L SAME");
    hists[i]->Draw("E3 SAME");
    
  }
  tl->Draw();
  c1->SetLogx();
  c1->SetLogy();
  data->SetMarkerStyle(20);
  data->Draw("E1 SAME");

  c1->cd();
  c1->Draw();

  zrap_Prelim(0.85,0.975,0.4,0.17);
  zrap_Lumi(0.81,0.91,36);
  zqt_Cut(0.33,0.43,"l",0.75);

  p2->cd();

  //  TCanvas* c2=new TCanvas("c2","c2",800,400);
  //  TH2* dummy=new TH2F("dummy","dummy",20,0.0,30.0,30,-0.35,0.35);
  TH2* dummy=new TH2F("dummy","",20,19.5,600.0,30,-5.0,5.0);
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
  
  TMatrixD comboCovFull=readCov("dsdqt_combined_covariance.dat");
  TMatrixD comboCov(1,10,1,10);
  for (int i=9; i<=18; i++) 
    for (int j=9; j<=18; j++) {
      if (i==j) 
	comboCov(i-8,j-8)=comboCovFull(i,j)+pow(terror[i],2);
      else 
	comboCov(i-8,j-8)=comboCovFull(i,j);
    }

  //comboCov.Print();

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
    histsDelta[i]->Draw("SAME P E1");
  }
  /*
  tl2->SetNColumns(2);

  tl2->Draw();
  */
  //  zrap_Prelim(0.18,0.97,0.4,0.17);
  //  zrap_Lumi(0.88,0.97,36);


  char outx[1024];
  if (outgoing!=0) {
    sprintf(outx,"%s",outgoing);
    zrapPrint(c1,outx);
  }

  comboCov.Invert();

  for (int i=0; i<nhists; i++) {

    double chi2mat=0;
    
    TMatrix deltaT(deltas[i]);
    deltaT.Transpose(deltaT);
    //    deltaT.Print();

    //    comboCov.Print();

    TMatrix a1=comboCov*deltas[i];

    //    a1.Print();
    
    TMatrix a2=deltaT*a1;

    //    a2.Print();

    chi2mat=a2(1,1);

    printf(" Comparison %s : %f/%d %f | %f/%d %f \n",knownModels[imodels[i]*3+1],
	   chi2[i],ndof,TMath::Prob(chi2[i],ndof),
	   chi2mat,ndof,TMath::Prob(chi2mat,ndof)
	   );

  }
}

TMatrixD readCov(const char* covFile) {
  TMatrixD m(1,18,1,18);

  FILE* f=fopen(covFile,"rt");
  int i,j;
  float val;
  char buffer[256];
  while (f!=0 && !feof(f)) {
    buffer[0]=0;
    fgets(buffer,200,f);
    int igot=sscanf(buffer," %d %d %f ",&i,&j,&val);
    if (igot==3) {
      m(i,j)=val;
    }
  }
  fclose(f);
  return m;
}
