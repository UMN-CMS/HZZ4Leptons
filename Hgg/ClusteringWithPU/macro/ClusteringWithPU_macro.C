#include "TFile.h"
#include "TH1F.h"
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
		bool doLog=false, std::string path=std::string("./"));


//int scwithpu(std::string theDirectory){
int scwithpu(char myword[100]){
  
  std::string theDirectory = std::string(myword) + std::string("/");
  std::string theDir       = std::string(theDirectory) + std::string("/"); 
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

  TFile *file2 = new TFile("PU-20.root","READ");
  theFiles.push_back(file2);   theLabels.push_back(std::string("PU-20 #xi=0")); theColors.push_back(kRed); theFillStyle.push_back(3365);
  TFile *file4 = new TFile("PU-20-0.01-mar29.root","READ");
  theFiles.push_back(file4);   theLabels.push_back(std::string("PU-20 #xi=0.01")); theColors.push_back(kGreen); theFillStyle.push_back(3354);
  TFile *file3 = new TFile("PU-20-0.02.root","READ");
  theFiles.push_back(file3);   theLabels.push_back(std::string("PU-20 #xi=0.02")); theColors.push_back(kViolet); theFillStyle.push_back(3356);
  //TFile *file3 = new TFile("PU-20-0.03.root","READ");
  //theFiles.push_back(file3);   theLabels.push_back(std::string("PU-20 #xi=0.03")); theColors.push_back(kViolet); theFillStyle.push_back(3356);
  //TFile *file4 = new TFile("PU-20-0.04.root","READ");
  //theFiles.push_back(file4);   theLabels.push_back(std::string("PU-20 #xi=0.04")); theColors.push_back(kGreen); theFillStyle.push_back(3354);
  //TFile *file5 = new TFile("PU-20-0.06.root","READ");
  //theFiles.push_back(file5);   theLabels.push_back(std::string("PU-20 #xi=0.06")); theColors.push_back(kAzure+10); theFillStyle.push_back(3354);
  TFile *file1 = new TFile("PU-0.root","READ");
  theFiles.push_back(file1);   theLabels.push_back(std::string("PU-0 #xi=0")); theColors.push_back(kBlue); theFillStyle.push_back(3395);
  //TFile *file7 = new TFile("PU-0-0.02.root","READ");
  //theFiles.push_back(file7);   theLabels.push_back(std::string("PU-0 #xi=0.02")); theColors.push_back(kAzure+10); theFillStyle.push_back(3356);
  TFile *file7 = new TFile("PU-0-0.01-mar29.root","READ");
  theFiles.push_back(file7);   theLabels.push_back(std::string("PU-0 #xi=0.01")); theColors.push_back(kAzure+10); theFillStyle.push_back(3356);

  std::string myPlot;
  std::vector<std::string> allMyPlots;
  myPlot = std::string("h_phiShape_barl"); allMyPlots.push_back(myPlot);
  //makeAplot( theFiles, theLabels, theColors, theFillStyle, myPlot , myfile, true, theDirectory);
  
  //myPlot = std::string("h_etaShape_barl"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_etaShape_barlPLus"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_etaShape_barlMinus"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_etaShape_barlSymm"); allMyPlots.push_back(myPlot);

  //myPlot = std::string("h_maxCryInDomino_barl"); allMyPlots.push_back(myPlot);

  //myPlot = std::string("h_maxCryInLocalMax_barlSymm"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_maxCryInLocalMax_barlPLus"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_maxCryInLocalMax_barlMinus"); allMyPlots.push_back(myPlot);


  //myPlot = std::string("h_phiShape_endc"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_phiSize_barl"); allMyPlots.push_back(myPlot);
  //myPlot = std::string("h_phiSize_endc"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_EoverEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_EoverEtrue_endc"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoEoverEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoER9overEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_PhoEnotR9overEtrue_barl"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_mHiggs_EBEB_trueVtx"); allMyPlots.push_back(myPlot);
  myPlot = std::string("h_mHiggs_EEEE_trueVtx"); allMyPlots.push_back(myPlot);
  

  
  std::string theLineForHtml = std::string("<h2 id=\"ECAL superclusters with pile up - ECAL95-ECAL95\"><A name=\"EB\"><FONT color=\"Black\">G. Franzoni and Y. Kubota") +   std::string("</FONT></A><BR></h2>"); 
  myfile << theLineForHtml << std::endl;
  for(std::vector<std::string>::iterator iter=allMyPlots.begin(); iter!=allMyPlots.end(); iter++){
    makeAplot( theFiles, theLabels, theColors, theFillStyle, (*iter) , 
	       myfile, 
	       true, theDirectory);
  }
  
  myfile << "<hr>" << std::endl; 

  return 0;
}




void makeAplot( std::vector<TFile*> theFiles, std::vector<std::string> theLabels, std::vector<int> theColors, std::vector<int> theFillStyle, 
		std::string myPlot, 
		ofstream &myfile,
		bool doLog , 
		std::string thePath)
{
  
  std::string theBasePath("scwithpuanalyzer/");
  std::string drawOption("he");
  int         countPlots(0);

  TCanvas *theCanvas = new TCanvas(myPlot.c_str(),myPlot.c_str(),150,10,990,760);
  theCanvas->cd();

  TLegend * leg = new TLegend(0.10,0.73,0.37,0.97);
  leg->SetHeader(myPlot.c_str());

  // loop over the files 
  std::vector< TFile* >::const_iterator constIterator;
  for ( constIterator = theFiles.begin();
	constIterator != theFiles.end(); constIterator++ ) {

    std::string thePlot = theBasePath + myPlot;
    TH1F * aHisto = (TH1F*)  (*constIterator)->Get(thePlot.c_str());
    aHisto->SetNormFactor(1);
    aHisto->SetLineWidth(2);
    aHisto->SetLineColor(theColors.at(countPlots));
    //aHisto->SetFillColor(theColors.at(countPlots));
    //aHisto->SetFillStyle(theFillStyle.at(countPlots));

    aHisto->Draw(drawOption.c_str());

    //std::cout << "option is : " << drawOption.c_str() << std::endl;
    //if(countPlots>0){
    //      std::cout << "any list? " << aHisto->GetListOfFunctions() << std::endl;
    //      TPaveStats * st = (TPaveStats*)aHisto->GetListOfFunctions()->FindObject("TPaveStats::stats");
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
	sigEffLabel = std::string("       #sigma_{eff} = ") + std::string( buffer )  + std::string(";");
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

  if (doLog){
    theCanvas->SetLogy();
    theCanvas->Print( (thePath+myPlot+std::string("_log")+std::string(".png")) .c_str() );
  }

  theLineForHtml = std::string("<A HREF=") +myPlot+std::string("_log")+std::string(".png")  + std::string("> <img height=\"300\" src=\"") +myPlot+std::string("_log")+std::string(".png")  + std::string("\"> </A>"); 
  myfile << theLineForHtml << std::endl;   



}// end makeAplot
