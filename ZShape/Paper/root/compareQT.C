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

TMatrixD readCov(const char* covFile);

int colorForId(int i) {
  const int colors[] = {kRed, kBlue, kGreen+2,
			kViolet, kGray
  };
  return colors[i];
}

void draw_axis_labels( float pos, float size, int iV );

void compareQT(const char* models,const char* outgoing=0) {
  setTDRStyle();
  bool limiNormal=false;

  int nhists=0;
  TH1* hists[10],*histsDelta[10];
  int imodels[10];
  char names[10][50];
  double chi2[10];
  TMatrixD deltas[10];
  
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
  //  c1->Divide(1,2);
  // c1->cd(1);

  TMatrixD comboCovFull=readCov("dsdqt_combined_covariance.dat");
  TMatrixD comboCov(1,9,1,9);
  for (int i=1; i<=9; i++) 
    for (int j=1; j<=9; j++)
      comboCov(i,j)=comboCovFull(i,j);

  //  comboCov.Print();

  for (int i=0; i<10; i++) {
    deltas[i].ResizeTo(1,9,1,1);
  }

  TH1* data=readStandardFilePt("data","../data/dsdqt_combined.csv");
  data->SetMaximum(0.08);
  data->SetMinimum(0.00);
  data->SetTitle(0);
  data->GetXaxis()->SetRangeUser(0.0,29.0);
  //  data->GetXaxis()->SetTitle(qt_xaxis_ll_label);
  data->GetXaxis()->SetTitle("");
  data->GetXaxis()->SetLabelSize(0.0);
  data->GetXaxis()->SetNdivisions(310);
  data->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T} [(GeV/c)^{-1} ]");
  data->GetYaxis()->SetTitleOffset(1.2);
  data->GetXaxis()->CenterTitle();
  
  data->Draw("E1");

  char* Leg="";
  if (limiNormal) Leg="MC Normalized over 0 < p_{T} < 30";
  const int ndof=8;

  TLegend* tl=new TLegend(0.45,0.45,0.9,0.8,Leg);
  tl->SetTextFont(42);
  tl->AddEntry(data,"data (e + #mu combined)","P");
  tl->SetFillColor(kWhite);
  tl->SetBorderSize(0);
  TLegend* tl2=new TLegend(0.4,0.75,0.9,0.9,Leg);
  tl2->SetFillColor(kWhite);
  tl2->SetBorderSize(0);
  //  tl->AddEntry(data,"Combined e + #mu","p");

  
  for (int i=0; i<nhists; i++) {
    char fullfname[1024];
    sprintf(fullfname,"../theory/%s",knownModels[imodels[i]*3+2]);
    printf("%s\n",fullfname);
    hists[i]=readStandardFilePt(knownModels[imodels[i]*3],fullfname);

    if (strstr(fullfname,"Resbos")) {
      double sum=hists[i]->Integral();
      for (int j=1; j<18; j++) 
	hists[i]->SetBinContent(j,hists[i]->GetBinContent(j)/hists[i]->GetBinWidth(j)/sum);
    }
    sprintf(fullfname,"%s-delta",knownModels[imodels[i]*3+2]);
    histsDelta[i]=(TH1*)(hists[i]->Clone(fullfname));

    chi2[i]=0;
    
    for (int j=1; j<18; j++) {
      histsDelta[i]->SetBinContent(j,(data->GetBinContent(j)-hists[i]->GetBinContent(j))/data->GetBinError(j));
      histsDelta[i]->SetBinError(j,data->GetBinError(j)/hists[i]->GetBinContent(j));
      if (j<10) {
	chi2[i]+=pow(data->GetBinContent(j)-hists[i]->GetBinContent(j),2)/
	  pow(data->GetBinError(j),2);
	deltas[i](1,j)=data->GetBinContent(j)-hists[i]->GetBinContent(j);
      }
    }
    
    hists[i]->SetLineWidth(2);
    hists[i]->SetLineColor(colorForId(i));

    histsDelta[i]->SetLineWidth(2);
    histsDelta[i]->SetLineColor(colorForId(i));

    tl->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");
    tl2->AddEntry(hists[i],knownModels[imodels[i]*3+1],"L");

    if (limiNormal) 
      hists[i]->Scale(data->Integral(1,9,"width")/hists[i]->Integral(1,9,"width"));

    hists[i]->SetFillStyle(0);
    histsDelta[i]->SetFillStyle(0);
    hists[i]->Draw("HIST SAME ][");
  }
  tl->Draw();

  data->Draw("E1 SAME");

  c1->cd();
  c1->Draw();
  zrap_Prelim(0.80,0.98,0.4,0.17);
  zrap_Lumi(0.75,0.90,36);
  zqt_Cut(0.34,0.43,"l",0.75);

  p2->cd();

  //  TCanvas* c2=new TCanvas("c2","c2",800,400);
  //  c2->SetRightMargin(0.05);
  //  TH2* dummy=new TH2F("dummy","dummy",20,0.0,30.0,30,-0.35,0.35);
  TH2* dummy=new TH2F("dummy","",20,0.0,30.0,30,-5.0,5.0);
  dummy->Draw();
  dummy->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  dummy->GetYaxis()->SetTitleOffset(0.6);
  dummy->GetYaxis()->SetTitleSize(0.10);
  dummy->GetXaxis()->SetLabelSize(0.09);
  dummy->GetYaxis()->SetLabelSize(0.09);
  dummy->GetYaxis()->SetNdivisions(503);
  dummy->GetXaxis()->SetNdivisions(310);
  dummy->GetXaxis()->SetTitleSize(0.10);
  dummy->GetYaxis()->SetTitle("(data - theory)/#sigma_{data}");

  const double x[4]={0,600,600,0};
  const double y1[4]={1,1,-1,-1};
  const double y2[4]={2,2,-2,-2};
  TGraph* band1s=new TGraph(4,x,y1);
  TGraph* band2s=new TGraph(4,x,y2);
  band1s->SetFillColor(kGreen-7);
  band2s->SetFillColor(kYellow);
  band2s->Draw("F SAME");
  band1s->Draw("F SAME");

  const double xx[2]={0,600};
  const double yy[2]={0,0};
  TGraph* zeroLine=new TGraph(2,xx,yy);
  zeroLine->SetLineStyle(3);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineWidth(1);
  zeroLine->Draw("L SAME");

  gPad->GetFrame()->Draw();

  /*
  double xv[11*2],yv[11*2];

  for (int i=1; i<=9; i++) {
    xv[i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[i]=histsDelta[0]->GetBinError(i);
  }
  xv[0]=0; yv[0]=yv[1];
  xv[10]=30; yv[10]=yv[9];
  for (int i=1; i<=9; i++) {
    xv[21-i]=histsDelta[0]->GetXaxis()->GetBinCenter(i);
    yv[21-i]=-histsDelta[0]->GetBinError(i);
  }
  xv[11]=30; yv[11]=yv[12];
  xv[21]=0; yv[21]=yv[20];

  TGraph* deltaBand=new TGraph(11*2,xv,yv);
  deltaBand->SetFillColor(kYellow);
  deltaBand->Draw("F SAME");
  */

  for (int i=0; i<nhists; i++) {
    histsDelta[i]->Draw("SAME HIST L");
  }

  //  tl2->SetNColumns(2);

  //  tl2->Draw();

  //  zrap_Prelim(0.18,0.97,0.4,0.17);
  //  zrap_Lumi(0.88,0.90,36);

  char outx[1024];
  if (outgoing!=0) {
    sprintf(outx,"%s",outgoing);
    zrapPrint(c1,outx);
    sprintf(outx,"%s_delta",outgoing);
    //    zrapPrint(c2,outx);
  }

  comboCov.Invert();

  for (int i=0; i<nhists; i++) {
    double chi2mat=0;
    
    TMatrix deltaT(deltas[i]);
    deltaT.Transpose(deltaT);

    TMatrix a2=deltaT*comboCov*deltas[i];

    chi2mat=a2(1,1);
    

    printf(" Comparison %s : %f/%d %f | %f/%d %f\n",knownModels[imodels[i]*3+1],
	   chi2[i],ndof,TMath::Prob(chi2[i],ndof),
	   chi2mat,ndof,TMath::Prob(chi2mat,ndof)
	   );

  }

}

void draw_axis_labels( float pos, float size, int iV )
{
  TLatex* txt_ = new TLatex();
  // draw axis labels
  //  float xlab_ = 65;
  txt_->SetTextFont(42);
  txt_->SetTextSize(size);
  txt_->SetTextAngle(0);
  if( iV==1 )
    {
      txt_->SetTextAlign(23);
    }
  else
    {
      txt_->SetTextAlign(32);
    }
  float mass[7] = {20,30,60,100,200,300,600};
  for( size_t ii=0; ii<7; ii++ )
    {
      float x_, y_;
      if( iV==1 )
	{
	  x_ = mass[ii];
	  y_ = pos;
	}
      else
	{
	  y_ = mass[ii];
	  x_ = pos;
	}
      txt_->SetText(x_,y_,Form("%-2.0f",mass[ii])); txt_->DrawClone();
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
