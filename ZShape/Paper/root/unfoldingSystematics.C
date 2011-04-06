#include <iostream>
#include <iomanip>
#include <stdio.h>


#include <TDatime.h>

#include "tdrstyle.C"
#include "zrapidityStandard.C"

#include "readUnfoldingMatrices.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

int unfoldingSystematics(std::string variable){

  // schema:
  // 0. get distribution of 'experimental' smeared rapidities, i.e. after the eff Acc trials (folders "ECAL-ECAL" and "ECAL-HF")
  //    divide it by eff acc   ==> the 'data' distribution obtained, I call it eaD
  //      \--> obtain it from the ratio of: [smeared rapidity in (ECAL-ECAL + ECAL-HF)] / [smeared in (All)]
  // 1 . unfold eaD with the matrix obtained from all the 'base' parameters
  //     \--> here one has to verify closure against true initial distribution (in a |Y| range, I guess)  
  // 2 . unfold edD with parameter-varied matrices
  // 3 . compute the ratios of point 2 over point 1
  // 4 . sum in quadrature the positive deviations and the negative deviations
  //     (this is adopted from the pdf studies - unsatisfactory but can do for now)

  bool isRapidity(true);
  if (variable==std::string("Y")){
	std::cout << "\n\tyou have chosen matrix systematics for Y\n\n" << std::endl;
  }  
  else if (variable==std::string("pt")){
	std::cout << "\n\tyou have chosen matrix systematics for pt\n\n" << std::endl;
	isRapidity=false;
  }
  else{
	std::cout << "\n\nnot done a valid choice of variable. You can choose Y or pt. Bailing out\n\n" << std::endl;
	return 0;
  }


  setTDRStyle();
  zrap_colors();

  //std::string filename("matrices/matrix-BASE.py.hadded.root ");  GF
  std::string filename("");
  if(isRapidity){  filename = std::string("matrices-Y/matrix-BASE.py.hadded.root");}
  else          {  filename = std::string("matrices-pt/matrix-BASE.py.hadded.root");}
  TFile *theBASEmatrix = new TFile(filename.c_str(),"read");

  // I need to divide this distribution by eff acc 
  std::string fileEffAcc("effAcc-V00-05-XX/BASE.py.hadded.root"); 
  TFile *theEffAcc = new TFile(fileEffAcc.c_str(),"read");

  std::string theInitTLString("");
  std::string theInitSmearedSring("");
  if(isRapidity){
    theInitTLString = std::string("mcEff/All/Z0_YTL");
    theInitSmearedSring = std::string("mcEff/All/Z0_Y");
    theInit2dSring = std::string("mcEff/All/YZTL_vs_YZ");  //     you can also get from the 2d plot (already including gen-level mass cut!)
  } else    {
    theInitTLString = std::string("mcEff/All/Z0_PtTL");
    theInitSmearedSring = std::string("mcEff/All/Z0_Pt");
    theInit2dSring = std::string("mcEff/All/PtTL_vs_PtZ"); //     you can also get from the 2d plot (already including gen-level mass cut!)
  }
  
  TH2 *  the2dPlot = theEffAcc->Get(theInit2dSring.c_str());
  TH1D * theInitialTreeLevelHisto = the2dPlot->ProjectionY();
  TH1D * theInitialHisto          = the2dPlot->ProjectionX();

  //TH1F * theInitialTreeLevelHisto = theEffAcc->Get(theInitTLString.c_str());         // this is the tree level Y before Eff Acc trials
  //TH1F * theInitialHisto          = theEffAcc->Get(theInitSmearedSring.c_str());     // this is the SMEARED Y before Eff Acc trials

  std::cout << "\tThere are: " << theInitialTreeLevelHisto->GetNbinsX() << " bins in the game\n\n" << std::endl;
  std::cout << "\tFirst bin left: " << theInitialTreeLevelHisto->GetBinLowEdge(1) << "\t" << theInitialTreeLevelHisto->GetBinContent(1)  << std::endl;
  std::cout << "\tFirst bin left: " << theInitialHisto->GetBinLowEdge(1)  << "\t" << theInitialHisto->GetBinContent(1) << std::endl;
  std::string theECALstring("");
  std::string theHFstring("");
  if(isRapidity){
    theECALstring = std::string("mcEff/ECAL80-ECAL95/C08-m(60,120)/Z0_Y");
    theHFstring = std::string("mcEff/ECAL80-HF/C08-m(60,120)/Z0_Y");
  } else    {
    theECALstring = std::string("mcEff/ECAL80-ECAL95-MUO/C08-m(60,120)/Z0_Pt");
    theHFstring = std::string("mcEff/ECAL80-HF/C08-m(60,120)/Z0_Pt");
  }

  
  // here looking at events which pass the eff acceptance and are smeared
  TH1F * theECAL  = theEffAcc->Get(theECALstring.c_str());
  TH1F * theHF    = theEffAcc->Get(theHFstring.c_str());
  if(!isRapidity) theHF->Reset();
  // theTotal is the smeared Y distribution after the eff Acc trials
  TH1F * theTotal = theHF->Clone("theTotal");
  theTotal->Reset();
  theTotal->Add(theECAL,theHF);

  // theTotal is events which pass the eff acceptance and are smeared
  TH1F * theEffAccHistoBase = theTotal->Clone();
  theTotal->Sumw2();
  // theEffAccHistoBase is the SMEARED Y before Eff Acc trials   
  theEffAccHistoBase->Divide(theInitialHisto);    // here is the efficiency X acc histogram
  theEffAccHistoBase->SetTitle("theEffAccHistoBase");
  TCanvas * theEffAccCanvas = new TCanvas("theEffAccCanvas","theEffAccCanvas",100,100,750,600);  
  theEffAccHistoBase->Draw();


  // now prepare the base rapidity histogram to be passed in to the parameter-varied matrices 
  TH1F *   theEffAccCorrBaseRapidity = theTotal->Clone();
  theEffAccCorrBaseRapidity->Sumw2();
  // now eff-ACC correct the base 
  theEffAccCorrBaseRapidity->Divide(theEffAccHistoBase);
  TH1F    * baseUnfoldedRapidity = unfold(theEffAccCorrBaseRapidity,filename.c_str());
  TCanvas * theCheck = new TCanvas("check the closure","check the closure",150,150, 750,600);   theCheck->cd();
  baseUnfoldedRapidity->Draw("");
  theInitialTreeLevelHisto->Draw("same");
  TLegend* leg = new TLegend(0.1,0.7,0.30,0.9);
  leg->SetHeader("closure, isn't it");
  if(isRapidity){
  leg->AddEntry(baseUnfoldedRapidity,"Y_{Z,reco} after effAcc and unsmearing","lp");
  leg->AddEntry(theInitialTreeLevelHisto,"Y_{Z,tree level}","pl");
  } else {
    leg->AddEntry(baseUnfoldedRapidity,"p_{T,ee reco} after effAcc and unsmearing","lp");
    leg->AddEntry(theInitialTreeLevelHisto,"q_{T,ee tree level}","pl");
  }
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); 
  leg->Draw();  



  TCanvas * CanvasTLMinusRecoOverTL = new TCanvas("TL-Reco/TL","TL-Reco/TL",150,150,800,650);
  if(!isRapidity) CanvasTLMinusRecoOverTL->SetLogx();
  TH1F * TLMinusRecoOverTL = (TH1F*) theInitialTreeLevelHisto->Clone();
  TLMinusRecoOverTL->Sumw2();
  theInitialHisto->Sumw2();
  TLMinusRecoOverTL->Add(theInitialHisto,-1);
  TLMinusRecoOverTL->Divide(theInitialTreeLevelHisto);
  TLMinusRecoOverTL->Divide(theInitialHisto);
  TH1F* TLMinusRecoOverTLBis = (TH1F*)TLMinusRecoOverTL->Clone();
  TLMinusRecoOverTLBis->GetXaxis()->SetTitle("Y_{Z0}");
  TLMinusRecoOverTLBis->GetXaxis()->CenterTitle();
  TLMinusRecoOverTLBis->GetYaxis()->SetTitle("Evt(true) - Evt(smear) / Evt(true)");
  TLMinusRecoOverTLBis->GetYaxis()->CenterTitle();
  TLMinusRecoOverTLBis->SetTitleOffset(1.2);
  TLMinusRecoOverTLBis->SetLineColor(kBlack);
  TLMinusRecoOverTLBis->Draw("hpx e HIST");
  TLMinusRecoOverTL->GetXaxis()->SetTitle("Y_{Z0}");
  TLMinusRecoOverTL->SetLineColor(kBlue);
  CanvasTLMinusRecoOverTL->cd();
  TLMinusRecoOverTL->Draw("hpx HIST same");

  TCanvas * CanvasRecoOverTL = new TCanvas("Reco/TL","Reco/TL",0,0, 1050,750); CanvasRecoOverTL->cd();
  TH1F * RecoOverTL = (TH1F*) theInitialTreeLevelHisto->Clone();
  theInitialHisto->Sumw2();
  RecoOverTL->Sumw2();
  RecoOverTL->Divide(theInitialHisto);
  RecoOverTL->SetLineColor(kBlack);
  TH1F*  RecoOverTLBis = RecoOverTL->Clone(); 
  RecoOverTL->GetYaxis()->SetTitle("#rho=Evt(true)/Evt(reco)");
  if(isRapidity)    {     
    RecoOverTL->GetXaxis()->SetTitle(y_axis_label); }
  else              {
    RecoOverTL->GetXaxis()->SetTitle(qt_xaxis_label);
    RecoOverTLBis->GetXaxis()->SetTitle(qt_xaxis_label); 
    RecoOverTL->SetTitleOffset(1.2);
    RecoOverTLBis->SetTitleOffset(1.2);
    CanvasRecoOverTL->SetLogx();
  }
  RecoOverTL->GetXaxis()->SetTitleOffset(1.3);
  RecoOverTL->GetXaxis()->CenterTitle();
  RecoOverTLBis->SetLineColor(kBlack);
  RecoOverTLBis->Draw("hpx e HIST");
  RecoOverTL->GetXaxis()->SetTitle("Y_{Z0}");
  RecoOverTLBis->SetLineColor(kBlue);

  TH1D* theInitialTreeLevelHistoKevin = (TH1D*)zpt_rebinForPlot(theInitialTreeLevelHisto);
  TH1D* theInitialHistoKevin = (TH1D*)zpt_rebinForPlot(theInitialHisto);
  //theInitialTreeLevelHistoKevin->Draw("hpx e HIST");
  //theInitialHistoKevin->Draw("hpx HIST same");

  TH1D* RecoOverTLKevin;
  TH1D* RecoOverTLBisKevin;
  
  if(isRapidity)    {     
    RecoOverTLKevin = (TH1D*)RecoOverTL;
    RecoOverTLBisKevin =(TH1D*)RecoOverTLBis;
    RecoOverTLKevin->GetXaxis()->SetTitle(y_xaxis_label); }
  else              {
    RecoOverTLKevin  = (TH1D*)zpt_rebinForPlot(RecoOverTL);
    RecoOverTLBisKevin=  (TH1D*)zpt_rebinForPlot(RecoOverTLBis);
    RecoOverTLKevin->GetXaxis()->SetTitle(qt_xaxis_label);
    RecoOverTLBisKevin->GetXaxis()->SetTitle(qt_xaxis_label); 
    RecoOverTLKevin->GetXaxis()->SetTitleOffset(0.98);
    RecoOverTLBisKevin->GetXaxis()->SetTitleOffset(0.98);
    RecoOverTLBisKevin->GetXaxis()->CenterTitle();
    RecoOverTLBisKevin->GetYaxis()->CenterTitle();
    CanvasRecoOverTL->SetLogx();
  }
  RecoOverTLKevin->GetYaxis()->SetTitle("#rho=Evt(true)/Evt(reco)");
  RecoOverTLKevin->GetXaxis()->CenterTitle();
  RecoOverTLBisKevin->SetLineColor(kBlue);
  RecoOverTLKevin->Draw("hx e HIST");
  RecoOverTLBisKevin->Draw("hx HIST same");


  TCanvas * theBaseCanvas = new TCanvas("(unsmeared-tree)/tree","(unsmeared-tree)/tree",200,200,1000,800);  
  TH1D * h_effectUnsmearingBase = new TH1D("Y unsmearing Base","Y unsmearing Base",theInitialTreeLevelHisto->GetNbinsX(),theInitialTreeLevelHisto->GetXaxis()->GetXmin(),theInitialTreeLevelHisto->GetXaxis()->GetXmax()); 
  h_effectUnsmearingBase->Add(baseUnfoldedRapidity, theInitialTreeLevelHisto, 1, -1);
  h_effectUnsmearingBase->Divide(theInitialTreeLevelHisto);
  theBaseCanvas->cd();
  h_effectUnsmearingBase->Draw();

  // what I want is: unfold the smeared BASE after it's been corrected by effAcc 
  // wtih unfold matrices obtained varying the parameters
  // compare the outcome to the treeLevel of the base



  TH1D* theVariations[100];
  int index=0;

  TCanvas * theCanvases[10];
  TH2D * theEarlyDummy;
  if(isRapidity)    { 
    theEarlyDummy = new TH2D("EarlyDummy","EarlyDummy",100,-3.6,3.6,2,-0.01,0.01); 
    theEarlyDummy->GetXaxis()->SetTitle(y_xaxis_label);
    theEarlyDummy->GetXaxis()->CenterTitle();
  } else {
    theEarlyDummy = new TH2D("EarlyDummy","EarlyDummy",100,0,3500,2,-0.025,0.025);
    theEarlyDummy->GetXaxis()->SetTitle(qt_xaxis_label);
    theEarlyDummy->GetXaxis()->CenterTitle();
  }
  theEarlyDummy->GetYaxis()->SetTitle("(events_{varied unfolding}/events_{base unfolding}) -1");
  theEarlyDummy->SetTitle(0);

  TLegend* leg;
  char buffer[10];

  if(isRapidity)  filename = std::string("matrices-Y/matrix-BASE.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-BASE.py.hadded.root");

  theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,baseUnfoldedRapidity,filename.c_str());

  
  
  
  if(isRapidity) filename = std::string("matrices-Y/matrix-a-EB-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-a-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  if(isRapidity) filename = std::string("matrices-Y/matrix-a-EB-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-a-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
  theCanvases[index/2-1]->cd(); if(!isRapidity) theCanvases[index/2-1]->SetLogx();
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  
  
  if(isRapidity) filename = std::string("matrices-Y/matrix-a-EE-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-a-EE-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  if(isRapidity) filename = std::string("matrices-Y/matrix-a-EE-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-a-EE-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
  theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  
  
  if(isRapidity) filename = std::string("matrices-Y/matrix-alpha-EB-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-alpha-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  if(isRapidity) filename = std::string("matrices-Y/matrix-alpha-EB-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-alpha-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
  theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  
  

  if(isRapidity) filename = std::string("matrices-Y/matrix-alpha-EE-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-alpha-EE-down.py.hadded.root");
 index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
 if(isRapidity)  filename = std::string("matrices-Y/matrix-alpha-EE-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-alpha-EE-up.py.hadded.root");
 index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
  theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());


  
  if(isRapidity) filename = std::string("matrices-Y/matrix-c-EB-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  if(isRapidity) filename = std::string("matrices-Y/matrix-c-EB-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
  theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());


  if(isRapidity) filename = std::string("matrices-Y/matrix-c-HFm-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-HFm-down.py.hadded.root");
  if(isRapidity) {index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());}
  if(isRapidity)  filename = std::string("matrices-Y/matrix-c-HFm-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-HFm-up.py.hadded.root");
  // the HF is not included in the pt analysis
  if(isRapidity) {
    index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
    theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
    theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
    theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
    leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
    leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
    leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  }
  
  if(isRapidity) filename = std::string("matrices-Y/matrix-c-HFp-down.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-HFp-down.py.hadded.root");
  if(isRapidity) {index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());}
  if(isRapidity) filename = std::string("matrices-Y/matrix-c-HFp-up.py.hadded.root");
  else            filename = std::string("matrices-pt/matrix-c-HFp-up.py.hadded.root");
  // the HF is not included in the pt analysis
  if(isRapidity){
    index++; theVariations[index] = computeRelativeDifference(theEffAccCorrBaseRapidity,theInitialTreeLevelHisto,filename.c_str());
    theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100,100,700,500);
    theCanvases[index/2-1]->cd();  if(!isRapidity) theCanvases[index/2-1]->SetLogx();
    theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0); theVariations[index-1]->Draw(" hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
    leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
    leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
    leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  }


  // loop over the parameters variation 

  TH1D * theEffectPositive = theInitialTreeLevelHisto->Clone("positiveEff");
  TH1D * theEffectNegative = theInitialTreeLevelHisto->Clone("negativeEff");
  theEffectPositive->Reset();
  theEffectNegative->Reset();
  
  for(int bin=1; bin<=theInitialTreeLevelHisto->GetNbinsX(); bin++){

    float cumulativePositive=0;
    float cumulativeNegative=0;
    
    for(int v=1; v<=index; v++){
      
      float theValue = theVariations[v]->GetBinContent(bin);
      if(theValue>0) {
	cumulativePositive    += theValue*theValue;
      }
      else{
	cumulativeNegative    += theValue*theValue;
      }
      
    }

    theEffectPositive->SetBinContent(bin,sqrt(cumulativePositive));
    theEffectNegative->SetBinContent(bin,-1*sqrt(cumulativeNegative));

  }



  TCanvas * theControlCanvasBis = new TCanvas(filename.c_str(),filename.c_str(),0,0, 1050,750); theControlCanvasBis->cd();
  theControlCanvasBis-> UseCurrentStyle();
  theControlCanvasBis->cd();

  TH2D * theDummy;
  if(isRapidity)    { 
    theDummy = new TH2D("EarlyDummy","EarlyDummy",100,-3.6,3.6,2,-0.01,0.01); 
    theDummy->GetXaxis()->SetTitle("Y_{Z0}");
    theDummy->GetXaxis()->CenterTitle();
  } else {
    TH2D* theDummy = new TH2D("EarlyDummy","EarlyDummy",100,0.7,600,2,-0.01,0.01);
    //    theDummy=(TH2)zpt_rebinForPlot(theDummyPrelim);
    theDummy->GetXaxis()->SetTitle(qt_xaxis_label);
    //theDummy->GetXaxis()->CenterTitle();
    theDummy->GetXaxis()->SetTitleOffset(0.96);
    theControlCanvasBis->SetLogx(); 
  }

  theDummy->SetTitle(0);
  theDummy->GetYaxis()->SetTitle("(events_{varied unfolding}/events_{base}) -1");
  theDummy->GetXaxis()->CenterTitle();
  theControlCanvasBis->SetGridy();
  theDummy->Draw("");
  theEffectPositive->SetLineColor(kBlue);
  theEffectNegative->SetLineColor(kBlue);
  theEffectNegative->SetMarkerStyle(0);
  if(isRapidity){
  theEffectPositive->Draw("same hx HIST");
  theEffectNegative->  Draw("same hx HIST");
  } else {
  // use Kevin's tool to conform binning (and go around log probls with x==0)
  TH1D* theEffectPositiveKevin, *theEffectNegativeKevin;
  theEffectPositiveKevin=(TH1D*)zpt_rebinForPlot(theEffectPositive);
  theEffectPositiveKevin->SetLineColor(kBlue);
  theEffectPositiveKevin->SetMarkerStyle(0);
  theEffectNegativeKevin=(TH1D*)zpt_rebinForPlot(theEffectNegative);
  theEffectNegativeKevin->SetLineColor(kBlue);
  theEffectNegativeKevin->SetMarkerStyle(0);
  // removing 'p' from drawing options remove markers
  theEffectNegativeKevin->Draw("same hx HIST");
  theEffectPositiveKevin->Draw("same  hx HIST");
  }


 TH1F * histo;
 TH1D * histoKevin;
 for(int v=1; v<=index; v++){
   
   histo = (TH1F *) theVariations[v]->Clone();
   if(isRapidity){
     histoKevin = (TH1D*) histo;}
   else{
     histoKevin=(TH1D*)zpt_rebinForPlot(histo);
   }
   histoKevin->SetLineStyle(2);
   histoKevin->SetMarkerStyle(1);
   // removing 'p' from drawing options remove markers
   histoKevin->Draw("same hx HIST");
 }

  // http://root.cern.ch/root/html/THistPainter.html#HP01a
  //  TCanvas * theControlCanvasBisWIP = new TCanvas("pappo",filename.c_str(),0,0, 1050,750);  
  //  theVariations[1]->Draw("hx HIST");

  theDummy->Fill(0,0);
  theDummy->SetLineStyle(2);
  theDummy->SetMarkerStyle(0);
  theDummy->SetLineColor(kRed);
  theDummy->Draw("same");

  TLegend* legFinal = new TLegend(0.55,0.7,0.95,0.9);
  if(isRapidity)    { 
    legFinal= new TLegend(0.55,0.7,0.95,0.9);}
  else {
    legFinal= new TLegend(0.55,0.72,0.95,0.91);}

   legFinal->SetHeader("Unfolding systematics");
   legFinal->SetFillColor(0);
   legFinal->SetShadowColor(0);
   legFinal->SetTextSize(0.03);
   legFinal->AddEntry(theEffectPositive,"cumulative effect","l");
   legFinal->AddEntry(histo,"unsmearing variations","l");
   legFinal->Draw();


   //new way of getting the time 
   const char* time_;
   TDatime *mytime = new TDatime();
   time_ = mytime->AsString();
   TText *plabel = new TText();
   plabel-> SetNDC();
   plabel -> SetTextFont(1);
   plabel -> SetTextColor(1);
   plabel -> SetTextSize(0.04);
   plabel -> SetTextAlign(22);
   plabel -> SetTextAngle(0);
   TText *tlabel = new TText();
   tlabel-> SetNDC();
   tlabel -> SetTextFont(1);
   tlabel -> SetTextColor(1);
   tlabel -> SetTextSize(0.02);
   tlabel -> SetTextAlign(22);
   tlabel -> SetTextAngle(0);
   plabel -> DrawText(0.33, 0.91, "CMS preliminary 2010");
   tlabel -> DrawText(0.33, 0.88, Form("%s",time_));
   std::cout << "\n\n The local time is: " << time_ << "\n\n" << std::endl;

   theControlCanvasBis->Print("relative-errors-unfolding.png");
   theControlCanvasBis->Print("relative-errors-unfolding.eps");
   //   theControlCanvasBis->Print("relative-errors-unfolding.pdf");
   
   CanvasRecoOverTL->cd();
   plabel -> DrawText(0.55, 0.33, "CMS preliminary 2010");
   tlabel -> DrawText(0.55, 0.30, Form("%s",time_));
  

   std::string   systematicsFileName("");
   if(isRapidity){  systematicsFileName = std::string("relative-errors-unfolding-matrix-Y.txt");}
   else          {  systematicsFileName = std::string("relative-errors-unfolding-matrix-pt.txt");}
   
   ofstream myfile;
   myfile.open (systematicsFileName.c_str());
   myfile<< "#put your comments after a '#'.\n";
   myfile<<"#Variable name: unfolding_error\n";
   if(isRapidity)  myfile<<"RapBin   RapMin   RapMax    Value    Error\n";
   else            myfile<<"PtBin     PtMin      PtMax    Value    Error\n";
   for(int bin=1; bin<=theInitialTreeLevelHisto->GetNbinsX(); bin++){
     myfile.setf(ios::left);

     //     myfile
     myfile << bin
       	    << "\t" << theEffectPositive->GetBinLowEdge(bin)
       	    << "\t" << (theEffectPositive->GetBinLowEdge(bin) + theEffectPositive->GetBinWidth(bin));
     myfile.setf(ios::scientific);
     myfile << "\t"  << std::setprecision(3) << std::fabs((theEffectPositive->GetBinContent(bin) + fabs(theEffectNegative->GetBinContent(bin))  ) /2. );
     myfile.unsetf(ios::scientific);
     myfile << "\t0\n";
   }
   myfile.close();


   return 0;
}


TH1D * computeRelativeDifference(const TH1* effAccCorrectedBase, const TH1* theInitialTreeLevelHisto, const char* filename){
  
  // effAccCorrectedBase is simply effxacc- corrected Y distribution
  
  // unsmear the effAccCorrectedBase with the current matrix (so-called "UNsmeared_BAD") 
  TH1D * theUnfoldedDistribution = (TH1D *) unfold(effAccCorrectedBase,filename);

  // { UNsmeared_BAD[effxacc corrected base] - UNsmeared_RIGHT(effxacc corrected base)   } / UNsmeared_RIGHT(effxacc corrected base) 
  //TH1D * theRelativeDifference = new TH1D("Y unsmearing BH","Y unsmearing BH",100,-5,5); 
  TH1D * theRelativeDifference = theInitialTreeLevelHisto ->Clone();
  theRelativeDifference->Reset();
  
  //now I need the based unfoled...!~
  
  theRelativeDifference->Add(theUnfoldedDistribution,theInitialTreeLevelHisto,1,-1);
  theRelativeDifference->Divide(theInitialTreeLevelHisto);

  std::string returnedHisto = std::string("return_")+std::string(filename);
  
  TH1D * toBeReturned = theRelativeDifference->Clone(returnedHisto.c_str());
  toBeReturned->SetLineStyle(2);
  toBeReturned->SetMarkerStyle(0);
  return (TH1D *) toBeReturned;

}



int unfoldingSystematicsRatios(std::string variable){

  // schema:
  // 0. + get distribution of 'experimental' smeared rapidities, i.e. after the eff Acc trials (folders "ECAL-ECAL" and "ECAL-HF")
  //      divide it by eff acc  ==> the 'data' distribution obtained, I call it eaD
  //       \--> obtain it from the ratio of: [smeared rapidity in (ECAL-ECAL + ECAL-HF)] / [smeared in (All)]
  //            NumCandidates_TL_k is what comes from the truth
  //            NumCandidates_smeared_k is eff-Acc corrected
  //    + define rho_k = ( NumCandidates_TL_k / NumCandidates_smeared_k)
  // 1 . compute {rho_k} for the BASE
  // 2 . compute {rho_k} for each parameter variation
  // 3 . compute the ratios of point 2 over point 1
  // 4 . sum in quadrature the positive deviations and the negative deviations
  //     (this is adopted from the pdf studies - unsatisfactory but can do for now)

  bool isRapidity(true);
  if (variable==std::string("Y")){
	std::cout << "\n\tyou have chosen ratios systematics for Y\n\n" << std::endl;
  }  
  else if (variable==std::string("pt")){
	std::cout << "\n\tyou have chosen ratios systematics for pt\n\n" << std::endl;
	isRapidity=false;
  }
  else{
	std::cout << "\n\nnot done a valid choice of variable. You can choose Y or pt. Bailing out\n\n" << std::endl;
	return 0;
  }


  setTDRStyle();
  zrap_colors();


  std::string fileEffAcc("effAcc-V00-05-XX/BASE.py.hadded.root"); 
  TFile* theEffAcc = new TFile(fileEffAcc.c_str(),"read");
  
  
  std::string theInitTLString("");
  std::string theInitSmearedSring("");
  std::string theInit2dSring("");
  if(isRapidity){
    theInitTLString = std::string("mcEff/All/Z0_YTL");
    theInitSmearedSring = std::string("mcEff/All/Z0_Y");
    theInit2dSring = std::string("mcEff/All/YZTL_vs_YZ");
  } else    {
    theInitTLString = std::string("mcEff/All/Z0_PtTL");  
    theInitSmearedSring = std::string("mcEff/All/Z0_Pt");
    theInit2dSring = std::string("mcEff/All/PtTL_vs_PtZ");
  }

  TH2  * the2dPlot = theEffAcc->Get(theInit2dSring.c_str());
//  TH1D * theInitialTreeLevelHisto = the2dPlot->ProjectionY();
//  TH1D * theInitialHisto          = the2dPlot->ProjectionX();
  
  TH1F * theInitialTreeLevelHisto = theEffAcc->Get(theInitTLString.c_str());         // this is the tree level Y before Eff Acc trials
  TH1F * theInitialHisto          = theEffAcc->Get(theInitSmearedSring.c_str());     // this is the SMEARED Y before Eff Acc trials
  
  //  TH1F *  theInitialTreeLevelHisto = theEffAcc->Get("mcEff/All/Z0_YTL");   // this is the tree level Y before Eff Acc trials
  //  TH1F *  theInitialHisto          = theEffAcc->Get("mcEff/All/Z0_Y");     // this is the SMEARED Y before Eff Acc trials
  //  TH1F * theECAL  = theEffAcc->Get("mcEff/ECAL80-ECAL95/C08-m(70,110)/Z0_Y");
  //  TH1F * theHF    = theEffAcc->Get("mcEff/ECAL80-HF/C08-m(70,110)/Z0_Y");
  
  TH1F*    theVariations[100];
  TCanvas* theCanvases[100];
  int index=0;

  TH2D * theEarlyDummy;
  if(isRapidity)    { 
    theEarlyDummy = new TH2D("EarlyDummy","EarlyDummy",100,-3.6,3.6,2,-0.01,0.01); 
    theEarlyDummy->GetXaxis()->SetTitle(y_axis_label);
  } else {
    theEarlyDummy = new TH2D("EarlyDummy","EarlyDummy",100,0,3500,2,-0.025,0.025);
    theEarlyDummy->GetXaxis()->SetTitle(qt_xaxis_label);
  }
  theEarlyDummy->GetYaxis()->SetTitle("(#rho_{varied smearing}/#rho_{base}) -1");
  theEarlyDummy->GetXaxis()->CenterTitle();
  theEarlyDummy->SetTitle(0);
  theEarlyDummy->GetXaxis()->CenterTitle();



  TLegend* leg;
  char buffer[10];


  filename = std::string("effAcc-V00-05-XX/BASE.py.hadded.root");
  TH1F* rho_BASE = theInitialTreeLevelHisto->Clone();
  rho_BASE->Sumw2(); theInitialHisto->Sumw2();
  rho_BASE->Divide(theInitialHisto);  // this is the first rho I compute
  rho_BASE->GetXaxis()->CenterTitle();
  theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity); // since I am computing ratios or rho, this will be constatn 1 
  TCanvas * rhoForBase = new TCanvas("rhoForBase","rhoForBase",100,100,700,500);  // keep it anyway to set the structure
  rhoForBase->SetGridy();
  theVariations[index]->Draw();
  

  filename = std::string("effAcc-V00-05-XX/a-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/a-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());

  filename = std::string("effAcc-V00-05-XX/alpha-EE-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/alpha-EE-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());

  if(isRapidity) {
  filename = std::string("effAcc-V00-05-XX/c-HFp-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/c-HFp-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  }

  filename = std::string("effAcc-V00-05-XX/alpha-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/alpha-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());

  if(isRapidity) {
  filename = std::string("effAcc-V00-05-XX/c-HFm-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/c-HFm-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());
  }

  filename = std::string("effAcc-V00-05-XX/a-EE-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/a-EE-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());

  filename = std::string("effAcc-V00-05-XX/c-EB-up.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  filename = std::string("effAcc-V00-05-XX/c-EB-down.py.hadded.root");
  index++; theVariations[index] = computeRelativeDifferenceRho(rho_BASE,filename.c_str(),isRapidity);
  theCanvases[index/2-1] = new TCanvas(filename.c_str(),filename.c_str(),100+10*index,100+10*index,700+10*index,500+10*index);
  theCanvases[index/2-1]->cd(); if(!isRapidity) {theCanvases[index/2-1]->SetLogx();}
  theEarlyDummy->Draw(); theVariations[index]->SetLineStyle(0);
  theVariations[index-1]->Draw("hpx HIST same"); theVariations[index]->Draw(" hpx HIST same"); 
  leg = new TLegend(0.35,0.8,0.9,0.93);   leg->SetHeader(filename.c_str());
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.03); leg->AddEntry(theVariations[index],"+1#sigma","l"); leg->AddEntry(theVariations[index-1],"-1#sigma","l");
  leg->Draw(); sprintf (buffer, "%d", (index/2)); ; filename = std::string( buffer ) + std::string(".png"); theCanvases[index/2-1]->Print(filename.c_str());


  TH1D * theEffectPositive = theInitialTreeLevelHisto->Clone("positive");
  TH1D * theEffectNegative = theInitialTreeLevelHisto->Clone("negative");
  theEffectPositive->Reset();
  theEffectNegative->Reset();

  //  for(int bin=1; bin<=100; bin++){
  for(int bin=1; bin<=theInitialTreeLevelHisto->GetNbinsX(); bin++){

    float cumulativePositive=0;
    float cumulativeNegative=0;
    
    for(int v=1; v<=index; v++){
      
      float theValue = theVariations[v]->GetBinContent(bin);
      if(theValue>0) {
	cumulativePositive  += theValue*theValue;
      }
      else{
	cumulativeNegative  += theValue*theValue;
      }
      
    }

    theEffectPositive->SetBinContent(bin,sqrt(cumulativePositive));
    theEffectNegative->SetBinContent(bin,-1*sqrt(cumulativeNegative));

  }
 


  TCanvas * theControlCanvasBis = new TCanvas(filename.c_str(),filename.c_str(),0,0, 1050,750); theControlCanvasBis->cd();
  theControlCanvasBis-> UseCurrentStyle();
  TH2D * theDummy;
  if(isRapidity)    { 
    theDummy = new TH2D("EarlyDummy","EarlyDummy",100,-3.6,3.6,2,-0.01,0.01); 
    theDummy->GetXaxis()->SetTitle(y_xaxis_label);
    theDummy->GetXaxis()->CenterTitle();
  } else {
    theDummy = new TH2D("EarlyDummy","EarlyDummy",100,0.7,600,2,-0.025,0.025);
    theDummy->GetXaxis()->SetTitle(qt_xaxis_label);
    theDummy->GetXaxis()->CenterTitle();
    theControlCanvasBis->SetLogx();
  }
  theDummy->GetYaxis()->SetTitle("(#rho_{vaied smearing}/#rho_{base}) -1");
  theDummy->GetYaxis()->CenterTitle();
  theDummy->GetXaxis()->CenterTitle();
  theDummy->SetTitle(0);

  //theDummy->GetXaxis()->SetTitle("Y_{Z0}");
  //theDummy->GetYaxis()->SetTitle("(d#sigma/dY_{Z0,base unfolding}/d#sigma/dY_{Z0,varied unfolding}) -1");
  theControlCanvasBis->SetGridy();
  theControlCanvasBis->cd();
  theDummy->Draw("");
  theEffectPositive->SetLineColor(kBlue);
  theEffectNegative->SetLineColor(kBlue);
  theEffectNegative->SetMarkerStyle(0);
  // theEffectPositive->Draw("same");
  //  theEffectNegative->Draw("same");
  TH1D* theEffectPositiveKevin, *theEffectNegativeKevin;
  if(isRapidity){
    theEffectPositive->Draw("same hx HIST");
    theEffectNegative->Draw("same hx HIST");
  } else {
    // use Kevin's tool to conform binning (and go around log probls with x==0)
    theEffectPositiveKevin=(TH1D*)zpt_rebinForPlot(theEffectPositive);
    theEffectPositiveKevin->SetLineColor(kBlue);
    theEffectPositiveKevin->SetMarkerStyle(0);
    theEffectNegativeKevin=(TH1D*)zpt_rebinForPlot(theEffectNegative);
    theEffectNegativeKevin->SetLineColor(kBlue);
    theEffectNegativeKevin->SetMarkerStyle(0);
    // removing 'p' from drawing options remove markers
    theEffectNegativeKevin->Draw("same hx HIST");
    theEffectPositiveKevin->Draw("same  hx HIST");
  }



  TH1F * histo;
  TH1D * histoKevin;
  for(int v=1; v<=index; v++){
    
    histo = (TH1F *) theVariations[v]->Clone();
   if(isRapidity){
     histoKevin = (TH1D*) histo;}
   else{
     histoKevin=(TH1D*)zpt_rebinForPlot(histo);
   }
   histoKevin->SetLineStyle(2);
   histoKevin->SetMarkerStyle(1);
   // removing 'p' from drawing options remove markers
   histoKevin->Draw("same hx HIST");

  }




  // http://root.cern.ch/root/html/THistPainter.html#HP01a
  //  TCanvas * theControlCanvasBisWIP = new TCanvas("pappo",filename.c_str(),0,0, 1050,750);  
  //  theVariations[1]->Draw("hx HIST");

  theDummy->Fill(0,0);
  theDummy->SetLineStyle(2);
  theDummy->SetMarkerStyle(0);
  theDummy->SetLineColor(kRed);
  theDummy->Draw("same");


  TLegend* legFinal = new TLegend(0.55,0.7,0.95,0.9);
  if(isRapidity)    { 
    legFinal= new TLegend(0.55,0.7,0.95,0.9);}
  else {
    legFinal= new TLegend(0.55,0.15,0.95,0.35);}

   legFinal->SetHeader("Unfolding systematics");
   legFinal->SetFillColor(0);
   legFinal->SetShadowColor(0);
   legFinal->SetTextSize(0.03);
   legFinal->AddEntry(theEffectPositive,"cumulative effect","l");
   legFinal->AddEntry(histo,"unsmearing variations","l");
   legFinal->Draw();

   //new way of getting the time 
   const char* time_;
   TDatime *mytime = new TDatime();
   time_ = mytime->AsString();
   TText *plabel = new TText();
   plabel-> SetNDC();
   plabel -> SetTextFont(1);
   plabel -> SetTextColor(1);
   plabel -> SetTextSize(0.04);
   plabel -> SetTextAlign(22);
   plabel -> SetTextAngle(0);
   TText *tlabel = new TText();
   tlabel-> SetNDC();
   tlabel -> SetTextFont(1);
   tlabel -> SetTextColor(1);
   tlabel -> SetTextSize(0.02);
   tlabel -> SetTextAlign(22);
   tlabel -> SetTextAngle(0);
   plabel -> DrawText(0.33, 0.91, "CMS preliminary 2010");
   tlabel -> DrawText(0.33, 0.88, Form("%s",time_));
   std::cout << "\n\n The local time is: " << time_ << "\n\n" << std::endl;


   std::string   systematicsFileName("");
   if(isRapidity){  systematicsFileName = std::string("relative-errors-unfolding-average-Y.txt");}
   else          {  systematicsFileName = std::string("relative-errors-unfolding-average-pt.txt");}
   
   ofstream myfile;
   myfile.open (systematicsFileName.c_str());

   myfile<< "#put your comments after a '#'.\n";
   myfile<<"#Variable name: unfolding_error\n";
   if(isRapidity)  myfile<<"RapBin   RapMin   RapMax    Value    Error\n";
   else            myfile<<"PtBin     PtMin      PtMax    Value    Error\n";
   //   for(int bin=1; bin<=100; bin++){
   for(int bin=1; bin<=theInitialTreeLevelHisto->GetNbinsX(); bin++){
     myfile.setf(ios::left);
     
     //     myfile
     myfile << bin
       	    << "\t" << theEffectPositive->GetBinLowEdge(bin)
       	    << "\t" << (theEffectPositive->GetBinLowEdge(bin) + theEffectPositive->GetBinWidth(bin));
     myfile.setf(ios::scientific);
     myfile << "\t"  << std::setprecision(3) << std::fabs((theEffectPositive->GetBinContent(bin) + fabs(theEffectNegative->GetBinContent(bin) ) ) /2. );
     myfile.unsetf(ios::scientific);
     myfile << "\t0\n";
   }
   myfile.close();

   theControlCanvasBis->Print("unfoldSystematicsRatios.png");
   theControlCanvasBis->Print("unfoldSystematicsRatios.eps");
   //   theControlCanvasBis->Print("unfoldSystematicsRatios.pdf");

   return 0;
}


TH1F * computeRelativeDifferenceRho(const TH1* rho_BASE,  const char* filename, bool isRapidity){

  TFile* theEffAcc = new TFile(filename,"read");

  //  TH1F * theInitialTreeLevelHisto;
  //  TH1F * theInitialSmearedHisto;


  if(isRapidity){
    theInitialTreeLevelHisto = (TH1F *) theEffAcc->Get("mcEff/All/Z0_YTL");   // this is the tree level Y (actually always the same)
    theInitialSmearedHisto   = (TH1F *) theEffAcc->Get("mcEff/All/Z0_Y");     // this is the SMEARED Y before Eff Acc trials (for the present parameter variation)
  } else     {
    theInitialTreeLevelHisto = (TH1F *) theEffAcc->Get("mcEff/All/Z0_PtTL");   // this is the tree level pt (actually always the same)
    theInitialSmearedHisto   = (TH1F *) theEffAcc->Get("mcEff/All/Z0_Pt");     // this is the SMEARED pt before Eff Acc trials (for the present parameter variation)
  }
  
  /*
  std::string theInit2dString("");
  if(isRapidity){
    theInit2dString = std::string("mcEff/All/YZTL_vs_YZ");  // this has full info and gen-level mass cut
  } else    {
    theInit2dString = std::string("mcEff/All/PtTL_vs_PtZ"); // this has full info and gen-level mass cut
  }
  
  TH2  * the2dPlot                 = theEffAcc->Get(theInit2dString.c_str());
  TH1F * theInitialTreeLevelHisto  = (TH1F*) the2dPlot->ProjectionY();
  TH1F * theInitialSmearedHisto    = (TH1F*) the2dPlot->ProjectionX();
  */

  TH1F * theRhos   = theInitialTreeLevelHisto->Clone();
  theRhos->Sumw2();  theInitialSmearedHisto->Sumw2();
  theRhos->Divide(theInitialSmearedHisto);

  TH1F * theRatioOfRhos = theRhos->Clone(); 
  theRatioOfRhos->Sumw2(); rho_BASE->Sumw2();
  theRatioOfRhos->Add(rho_BASE,-1);
  theRatioOfRhos->Divide(rho_BASE);
  theRatioOfRhos->SetLineStyle(2);
  theRatioOfRhos->SetMarkerStyle(0);

  return theRatioOfRhos;
  
}
