#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPaveStats.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>


#include "Validation/EcalClusters/test/macro/effSigma.C"

int  scwithpu(char *);
void makeAplot( std::vector<TFile*> , std::vector<std::string> theLabels, std::vector<int> theColors, std::vector<int> theFillStyle, std::string myPlot, 
		ofstream &myfile, 
		//bool doLog=false,
		std::string path=std::string("./"));
void make2dAplot(TFile*, TFile* , std::vector<std::string> theLabels, std::vector<int> theColors, std::vector<int> theFillStyle,
		 std::string myPlot, 
		 ofstream &myfile, 
		 //bool doLog=false,
		 std::string path=std::string("./"));




int scwithpu(char myword[100]){
  
  std::string theDirectory = std::string(myword) + std::string("/");
  std::string theDir       = std::string(theDirectory);
  std::string theCommand   = std::string("mkdir ") + theDir; 
  std::cout << "locating output in: " << theDir << std::endl; 
  system(theCommand.c_str());

  std::string htmlFileName = theDir + std::string("index.html"); 
  ofstream myfile;
  myfile.open (htmlFileName.c_str()); 
  myfile << "<HEAD><TITLE>Superclusters with pile up studies</TITLE></HEAD>"<<std::endl;; 
  myfile << "<h1><A name=\"EB\"><FONT color=\"Blue\">ECAL Superclusters with pile up studies</FONT></A><BR></h1>"<<std::endl;; 

  gStyle->SetOptTitle(0);
  
  std::vector<TFile*>      theFiles;
  std::vector<std::string> theLabels;
  std::vector<int>         theColors;
  std::vector<int>         theFillStyle;

  TFile *file1 = new TFile("PU-0.root","READ");
  theFiles.push_back(file1);   theLabels.push_back(std::string("PU-0 #xi=0")); theColors.push_back(kBlue); theFillStyle.push_back(3395);
  TFile *file2 = new TFile("PU-20.root","READ");
  theFiles.push_back(file2);   theLabels.push_back(std::string("PU-20 #xi=0")); theColors.push_back(kRed); theFillStyle.push_back(3365);

  TFile *file11 = new TFile("PU-20-0.02.root","READ");
  theFiles.push_back(file11);   theLabels.push_back(std::string("PU-0 #xi=0.02")); theColors.push_back(kAzure+10); theFillStyle.push_back(3365);
  TFile *file12 = new TFile("PU-20-0.04_Su.root","READ");
  theFiles.push_back(file12);   theLabels.push_back(std::string("PU-20 #xi=0.04")); theColors.push_back(kGreen); theFillStyle.push_back(3365);
  TFile *file13 = new TFile("PU-20-0.06_Su.root","READ");
  theFiles.push_back(file13);   theLabels.push_back(std::string("PU-20 #xi=0.06")); theColors.push_back(kOrange); theFillStyle.push_back(3365);


  //TFile *file11 = new TFile("PU-0-deltaEta.root","READ");
  //theFiles.push_back(file11);   theLabels.push_back(std::string("PU-0 clus-w-#Delta#eta")); theColors.push_back(kAzure+10); theFillStyle.push_back(3365);
  //TFile *file12 = new TFile("PU-20-deltaEta.root","READ");
  //theFiles.push_back(file12);   theLabels.push_back(std::string("PU-20 clus-w-#Delta#eta")); theColors.push_back(kGreen); theFillStyle.push_back(3365);

  

  TFile *file0 = new TFile("PU-0.root","READ");
  theFiles.push_back(file0);   theLabels.push_back(std::string("PU-0 #xi=0")); theColors.push_back(kBlue); theFillStyle.push_back(3395);

  std::string myPlot;
  std::vector<std::string> allMyPlots;
  myPlot = std::string("h_phiShape_barl"); allMyPlots.push_back(myPlot);
  
  myPlot = std::string("h_phiSize_barl"); allMyPlots.push_back(myPlot);

  myPlot = std::string("h_absEtaBCminusAbsEtaSeed_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_absEtaBCminusAbsEtaSeed_barlm1"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_absEtaBCminusAbsEtaSeed_barlm2"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_absEtaBCminusAbsEtaSeed_barlm3"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_absEtaBCminusAbsEtaSeed_barlm4"); allMyPlots.push_back(myPlot);

  myPlot = std::string("h_phiBCminusPhiSeed_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_phiBCminusPhiSeed_barlm1"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_phiBCminusPhiSeed_barlm2"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_phiBCminusPhiSeed_barlm3"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_phiBCminusPhiSeed_barlm4"); allMyPlots.push_back(myPlot);

  myPlot = std::string("h_EoverEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoEoverEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoER9overEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoEnotR9overEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_mHiggs_EBEB_trueVtx"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_mHiggs_EEEE_trueVtx"); allMyPlots.push_back(myPlot);
  

  // heading to the web page
  std::string theLineForHtml = std::string("<h2 id=\"ECAL superclusters with pile up - ECAL95-ECAL95\"><A name=\"EB\"><FONT color=\"Black\">G. Franzoni and Y. Kubota") +   std::string("</FONT></A><BR></h2>");   myfile << theLineForHtml << std::endl;
  //  theLineForHtml = std::string("<h1><A name=\"EB\"><FONT color=\"Black\">Total number of events:") + std::string("") + std::string("</FONT></A><BR></h1>");   myfile << theLineForHtml << std::endl;

  std::cout << "\n\n" << std::endl;
  std::vector< TFile* >::const_iterator constIterator;
  for ( constIterator = theFiles.begin();
	constIterator != (theFiles.end()-1); constIterator++ ) {
    int numEvents = ( (TH1F*)  (*constIterator)->Get( (std::string("scwithpuanalyzer/h_nVtx")).c_str() ) ) ->GetEntries();
    int scMatchedEB = ( (TH1F*)  (*constIterator)->Get( (std::string("scwithpuanalyzer/h_scet_barl")).c_str() ) ) ->GetEntries(); 
    int scMatchedEE = ( (TH1F*)  (*constIterator)->Get( (std::string("scwithpuanalyzer/h_scet_endc")).c_str() ) ) ->GetEntries();
    int R9EB           = ( (TH1F*)  (*constIterator)->Get( (std::string("scwithpuanalyzer/h_PhoER9overEtrue_barl")).c_str() ) ) ->GetEntries();
    int notR9EB        = ( (TH1F*)  (*constIterator)->Get( (std::string("scwithpuanalyzer/h_PhoEnotR9overEtrue_barl")).c_str() ) ) ->GetEntries();


    std::cout << "Filename: " << (*constIterator)->GetName() << "  \tevts: " << numEvents
      //<< " MC-matched SC: " << (scMatchedEB+scMatchedEE)    // supercluster that matches a MC-truth particle 
	      << " MC-matched photons in EB: " << (R9EB+notR9EB)
      //<< " in EE: " << scMatchedEE
	      << "\t EB r9: " << R9EB << " ( " << (1.*R9EB)/(R9EB+notR9EB) << " ) " 
	      << "\t EB !r9: " << notR9EB << " ( " << (1.*notR9EB)/(R9EB+notR9EB) << " ) " 
	      << std::endl;
  }
  std::cout << "\n\n" << std::endl;
  
  // loop where the plots get actually made!
  for(std::vector<std::string>::iterator iter=allMyPlots.begin(); iter!=allMyPlots.end(); iter++){
    makeAplot( theFiles, theLabels, theColors, theFillStyle, (*iter) , 
	       myfile, 
	       /*true,*/ theDirectory);
  }


  theLabels.clear();
  theLabels.push_back(std::string("PU-0 #xi=0"));   theLabels.push_back(std::string("PU-20 #xi=0"));
  make2dAplot(file1 , file2, theLabels, theColors, theFillStyle, std::string("h_phiSizeVsE_barl"),myfile,/*true,*/ theDirectory);
  make2dAplot(file1 , file2, theLabels, theColors, theFillStyle, std::string("h_phiSizeVsEt_barl"),myfile,/*true,*/ theDirectory);
  

  myfile << "<hr>" << std::endl; 

  return 0;

}




void makeAplot( std::vector<TFile*> theFiles, std::vector<std::string> theLabels, std::vector<int> theColors, std::vector<int> theFillStyle, 
		std::string myPlot, 
		ofstream &myfile,
		//bool doLog , 
		std::string thePath)
{
  
  std::string theBasePath("scwithpuanalyzer/");
  std::string drawOption("he");
  int         countPlots(0);


  TCanvas *theCanvas = new TCanvas(myPlot.c_str(),myPlot.c_str(),150,10,990,760);
  theCanvas->cd();

  TLegend * leg = new TLegend(0.70,0.70,0.99,0.99);
  leg->SetHeader(myPlot.c_str());

  // loop over the files 
  std::vector< TFile* >::const_iterator constIterator;
  for ( constIterator = theFiles.begin();
	constIterator != (theFiles.end()-1); constIterator++ ) {

    std::string thePlot = theBasePath + myPlot;
    TH1F * aHisto = (TH1F*)  (*constIterator)->Get(thePlot.c_str());
    aHisto->SetNormFactor(1);
    aHisto->SetLineWidth(2);
    aHisto->SetLineColor(theColors.at(countPlots));
    //aHisto->SetFillColor(theColors.at(countPlots));
    //aHisto->SetFillStyle(theFillStyle.at(countPlots));

    aHisto->SetStats(0);
    aHisto->Draw(drawOption.c_str());
    //}
    
    std::string sigEffLabel("");
    std::vector<std::string> theHistosNeedingEffSigna;
    theHistosNeedingEffSigna.push_back("h_EoverEtrue_barl");
    theHistosNeedingEffSigna.push_back("h_EoverEtrue_endc");
    theHistosNeedingEffSigna.push_back("h_PhoEnotR9overEtrue_barl");
    theHistosNeedingEffSigna.push_back("h_mHiggs_EBEB_trueVtx");
    theHistosNeedingEffSigna.push_back("h_mHiggs_EEEE_trueVtx");
    theHistosNeedingEffSigna.push_back("h_PhoEoverEtrue_barl");
    theHistosNeedingEffSigna.push_back("h_PhoER9overEtrue_barl");
    theHistosNeedingEffSigna.push_back("h_PhoEnotR9overEtrue_barl");
    for(std::vector<std::string>::iterator iter=theHistosNeedingEffSigna.begin(); iter!=theHistosNeedingEffSigna.end(); iter++){
      size_t found = myPlot.find(   (*iter)   );
      if( found!=string::npos ){
	std::cout << "eff sigma for plot: " << myPlot << " is: " << effSigma(aHisto) << std::endl;
	char buffer [50]; 
	sprintf (buffer, "%.2e ", effSigma(aHisto) );
	sigEffLabel = std::string("    #sigma_{eff} = ") + std::string( buffer )  + std::string("");
      }
    }// loop over hist needing sigmaeff
    std::string theFullLabel = theLabels.at(countPlots) + sigEffLabel;
    leg->AddEntry(aHisto,theFullLabel.c_str());


    drawOption = std::string("hesames");
    countPlots++;

  }
  // loop over the files 

  leg->Draw();
  theCanvas->Print( (thePath+myPlot+std::string(".png")) .c_str() );
  std::string theLineForHtml = std::string("<A HREF=") +myPlot+std::string(".png") + std::string("> <img height=\"300\" src=\"")+myPlot+std::string(".png") + std::string("\"> </A>"); 
  myfile << theLineForHtml << std::endl; 
  
  //if (doLog){
  //    theCanvas->SetLogy();
  //    theCanvas->Print( (thePath+myPlot+std::string("_log")+std::string(".png")) .c_str() );
  //  }

  theLineForHtml = std::string("<A HREF=") +myPlot+std::string("_log")+std::string(".png")  + std::string("> <img height=\"300\" src=\"") +myPlot+std::string("_log")+std::string(".png")  + std::string("\"> </A>"); 
  myfile << theLineForHtml << std::endl;   

}// end makeAplot




void make2dAplot( TFile *file1 , TFile *file2, std::vector<std::string> theLabels, std::vector<int> theColors, std::vector<int> theFillStyle, 
		std::string myPlot, 
		ofstream &myfile,
		  //bool doLog , 
		std::string thePath)
{
  
  std::string theBasePath("scwithpuanalyzer/");
  std::string drawOption("he");
  //int         countPlots(0);
  std::string thePlot = theBasePath + myPlot;

  TCanvas *theCanvas = new TCanvas(myPlot.c_str(),myPlot.c_str(),150,10,990,760);
  theCanvas->Divide(1,3);
  theCanvas->cd(1);

  TLegend * leg = new TLegend(0.05,0.85,0.20,0.99);
  //leg->SetHeader(myPlot.c_str());

  TH2F * aHisto       = (TH2F*)  file1->Get(thePlot.c_str());
  TH2F * anotherHisto = (TH2F*)  file2->Get(thePlot.c_str());

  aHisto        ->SetLineColor(theColors.at(0));
  anotherHisto  ->SetLineColor(theColors.at(1));

  theCanvas->cd(1);  aHisto->Draw("box");
  bool DoFitting(true);
  if(true){  TProfile *a = ((TProfile*)aHisto->ProfileX()); a->SetLineColor(kBlack); a->Draw("same");}
  else{
    TH2F* theClone = (TH2F*) aHisto->Clone("pippo");
    //   aHisto->FitSlicesY();
    theClone->FitSlicesY();
    std::string theNameOfFitted = std::string("pippo") + std::string("_1");
    //    std::string theNameOfFitted = myPlot + std::string("_1");
    TH1D* fittedMean = (TH1D*)gDirectory->Get(theNameOfFitted.c_str());
    //    fittedMean = (TH1D*)fittedMean->Clone("pippo");
    fittedMean->Draw("same");
    theNameOfFitted = std::string("pippo") + std::string("_2");
    //    theNameOfFitted = std::string("pippo") + std::string("_2");
    TH1D* fittedSigma = (TH1D*)gDirectory->Get(theNameOfFitted.c_str());
    fittedSigma->GetYaxis()->SetTitle("#sigma( |#eta_{BC}| - |#eta_{SCseed}|)");
    fittedSigma->SetLineColor(theColors.at(0));
    theCanvas->cd(3);
    fittedSigma->SetStats(0);
    fittedSigma->Draw();
  }


  theCanvas->cd(2);  anotherHisto ->Draw("box");
  if(true){    TProfile *b = ((TProfile*)anotherHisto->ProfileX()); b->SetLineColor(kBlack); b->Draw("same");}
  else{
    TH2F* theClone = (TH2F*) anotherHisto->Clone("pluto");
    //anotherHisto->FitSlicesY();
    theClone->FitSlicesY();
    //std::string theNameOfFitted = myPlot + std::string("_1");
    std::string theNameOfFitted = std::string("pluto") + std::string("_1");
    TH1D* fittedMean = (TH1D*)gDirectory->Get(theNameOfFitted.c_str());
    fittedMean->Draw("same");
    theNameOfFitted = std::string("pluto") + std::string("_2");
    TH1D* fittedSigma = (TH1D*)gDirectory->Get(theNameOfFitted.c_str());
    fittedSigma->GetYaxis()->SetTitle("#sigma( |#eta_{BC}| - |#eta_{SCseed}|)");
    //fittedSigma->SetLineColor(theColors.at(1));
    fittedSigma->SetLineColor(theColors.at(1));
    theCanvas->cd(3);
    fittedSigma->SetStats(0);
    fittedSigma->Draw("same");
  }


  theCanvas->cd(1);
  leg->AddEntry(aHisto,theLabels.at(0).c_str());
  leg->AddEntry(anotherHisto,theLabels.at(1).c_str());
  leg->Draw("");

  //theCanvas->SetLogz(1);
  //theCanvas->SetLogz(2);
  theCanvas->Print( (thePath+myPlot+std::string(".png")) .c_str() );
  std::string theLineForHtml = std::string("<A HREF=") +myPlot+std::string(".png") + std::string("> <img height=\"300\" src=\"")+myPlot+std::string(".png") + std::string("\"> </A>"); 
  myfile << theLineForHtml << std::endl; 


}// end make2dAplot




//  theLabels.clear();
//  theLabels.push_back(std::string("PU-0 #xi=0"));   theLabels.push_back(std::string("PU-0 clus-w-#Delta#eta"));
//  make2dAplot(file1 , file11, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl"),myfile,true, theDirectory);
//  make2dAplot(file1 , file11, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1"),myfile,true, theDirectory);
//  make2dAplot(file1 , file11, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2"),myfile,true, theDirectory);
//  make2dAplot(file1 , file11, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3"),myfile,true, theDirectory);
//  make2dAplot(file1 , file11, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4"),myfile,true, theDirectory);
//
//  theLabels.clear();
//  theLabels.push_back(std::string("PU-20 #xi=0"));   theLabels.push_back(std::string("PU-20 clus-w-#Delta#eta"));
//  make2dAplot(file2 , file12, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barl"),myfile,true, theDirectory);
//  make2dAplot(file2 , file12, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm1"),myfile,true, theDirectory);
//  make2dAplot(file2 , file12, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm2"),myfile,true, theDirectory);
//  make2dAplot(file2 , file12, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm3"),myfile,true, theDirectory);
//  make2dAplot(file2 , file12, theLabels, theColors, theFillStyle, std::string("h_phiBCminusPhiSeed_VS_etaBCminusEtaSeed_barlm4"),myfile,true, theDirectory);

