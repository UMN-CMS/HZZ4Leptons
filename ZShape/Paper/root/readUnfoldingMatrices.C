#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include "TH1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "../root/tdrstyle.C"
#include "../root/zrapidityStandard.C"
#include <math.h>

// In this example, you pick up an unfolding matrix (made with makeUnfoldingMatrices.C) from a file (unfoldingMatrix_theOutPut.root)
// and apply to a rapidity/pt distribution you get from an effAcc file (effAccSource.root)

// the function:
//  readUnfoldingMatrices(std::string unfoldingMatrixFileInputFile, std::string effAccFileForTest, std::string variable) 
// works as well for Rapidity and transverse momentum, you have to choose when you launch it
// the functions:
//    TH1* unfold(const TH1* source, TMatrixD * theUnfoldingMatrix)
//    TH1* unfold(const TH1* source, const char* unfoldingMatrixIF)
// also as well for Rapidity and transverse momentum, but need not to know whether you're unfolding Y or py

// note if you want to _read_ the code: for how this code developed (written initially for rapidity)
// somewhere the 'variable' which is to be unsmeared is referred to as 'rapidity',
// while the macro consistently works for both Pt and Y


// to execute this example
//      wget http://homepages.spa.umn.edu/~franzoni/reps/11/jan12/effAccSource.root
//      (getting this file is only needed to run the example)
//      pb-d-128-141-34-30:2011-01-12-migMatrixReleaseFRESH franzoni$ root -l 
//      root [0] .L readUnfoldingMatrices.C                                                  
//      root [1] readUnfoldingMatrices("matrices-Y/matrix-BASE.py.hadded.root","effAcc-V00-05-05/BASE.py.hadded.root","Y")  
//      root [2] readUnfoldingMatrices("matrices-pt/matrix-BASE.py.hadded.root","effAcc-V00-05-05/BASE.py.hadded.root","pt")

TH1* unfold(const TH1* source, TMatrixD * theUnfoldingMatrix);

TH1* unfold(const TH1* source, const char* unfoldingMatrixIF, bool isFoldedRap=false)  {
  // load the matrix...

  TFile theunfoldingMatrixInputFile(unfoldingMatrixIF,"read");

  TMatrixD * theUnfoldingMatrix;
  theUnfoldingMatrix = (TMatrixD*)theunfoldingMatrixInputFile.Get("unsmearMatrices/unfoldingMatrixTotal");
  if (isFoldedRap) {
    theUnfoldingMatrix = (TMatrixD*)theunfoldingMatrixInputFile.Get("unsmearMatrices/unfoldingMatrixTotalFolded");
    std::cout << "you've chosen to unfold |Y|" << std::endl;
  }
  return unfold(source, theUnfoldingMatrix);
}



TH1* unfold(const TH1* source, TMatrixD * theUnfoldingMatrix)  {
  Double_t arrayRapiditySmeared[100]; // this needs be hardcoded...

  for(int bin=1; bin<=source->GetNbinsX() ; bin++)
    {
      arrayRapiditySmeared[bin-1]=source->GetBinContent(bin);
      if (0) std::cout << source->GetBinContent(bin) << std::endl;
    }
  
  
  //setTDRStyle();
  //zrap_colors();
  
  
  TVectorD vectorRapiditySmeared; vectorRapiditySmeared.Use(source->GetNbinsX(),arrayRapiditySmeared);
  
  // this is the multiplication for clusure test
  TVectorD vectorRapidityUNSmeared = (*theUnfoldingMatrix) * vectorRapiditySmeared;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  char name[1024];
  sprintf(name,"%s_unfolded",source->GetName());
  TH1D * h_RapidityUNSmeared = (TH1D*)(source->Clone(name));
  h_RapidityUNSmeared->SetDirectory(0);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // turning the unsmeared vector into a histogram
  for(int bin=1; bin<=source->GetNbinsX() ; bin++)
    {
      h_RapidityUNSmeared->SetBinContent(bin,vectorRapidityUNSmeared[bin-1]);
      if(fabs(vectorRapiditySmeared[bin-1]) > 1.e-9) h_RapidityUNSmeared->SetBinError(bin,vectorRapidityUNSmeared[bin-1]/vectorRapiditySmeared[bin-1]*source->GetBinError(bin));
      else                                           h_RapidityUNSmeared->SetBinError(bin,source->GetBinError(bin));
    }
  return h_RapidityUNSmeared;
}

void makeCovarianceMatrix(const char* file,  TMatrixD * theUnfoldingMatrix)  {
  
  setTDRStyle();
  zrap_colors();

  // get unfolding matrix and make its transpose

  TMatrixD unfoldingMatrixRaw           = *((TMatrixD*) theUnfoldingMatrix->Clone("theUnfoldingMatrix")) ;
  // necessary only for pt since results have one less bin 
  //unfoldingMatrixRaw.ResizeTo(18,18);
  TMatrixD unfoldingMatrix = (TMatrixD)unfoldingMatrixRaw;
  // necessary only for pt since results have one less bin 
  //unfoldingMatrix.ResizeTo(18,18);
  TMatrixD transpUnfoldingMatrix     =  (        theUnfoldingMatrix->Transpose((*theUnfoldingMatrix))        );
  // necessary only for pt since results have one less bin 
  //transpUnfoldingMatrix.ResizeTo(18,18);

  ofstream myfileData;
  myfileData.open ("originalDataTextFile.txt");
  myfileData<<"#bin\tqt_min\tqt_max\tvalue\tstat error" << std::endl;
  int progressive=1;

  // import total errors
  Double_t errorsArray[1000];
  Double_t systErrosArray[1000];
  int   errorsCounter(0);
  FILE* f=fopen(file,"r");
  if (f==0) {
    fprintf(stderr,"Cannot open '%s'\n",file);
    return;
  }
  char line[4096];
  while (!feof(f)) {
    line[0]=0;
    fgets(line,4095,f);
    if (strlen(line)<5 || line[0]=='#') continue;
    int i;
      float a,b,c,d,e;
      //      int ate, atetotal=0;
      int found=sscanf(line," %d %f %f %f %f %f",&i,&a,&b,&c,&d,&e);

      errorsArray[errorsCounter] = d*d; // use statistical to compute covariance matrix
      systErrosArray[errorsCounter] = d*d+e*e;
      std::cout << "stat errors2 directly from unfolded distribution - bin : " << i  <<"\t" << " err2: " << sqrt(errorsArray[errorsCounter])  << std::endl;    
      myfileData << scientific << progressive << "\t"
		 << a << "\t" //bin_min
		 << b << "\t"  //bin_max
		 << c << "\t" //value
	//<< sqrt(d*d) //staterr
	//<<   sqrt(e*e) //systerr
		 <<   sqrt(e*e+d*d) //TOTterr
		 << std::endl;
      progressive++;
      errorsCounter++;
      if (found!=3) continue;
  }
  fclose(f);
  // this is the file contains the erros of the smeared results that I start from 
  myfileData.close();

  TVectorD totalErrorsVector; totalErrorsVector.Use(errorsCounter,errorsArray);
  TMatrixD errorsOnDiagonalMatrix(errorsCounter,errorsCounter);  
  TMatrixDDiag diag1(errorsOnDiagonalMatrix); diag1 = totalErrorsVector;
  for(int i=0; i<errorsCounter; i++){
    for(int j=0; j<errorsCounter; j++){
      if( errorsOnDiagonalMatrix(i,j) !=0){       std::cout << "errorsOnDiagonal: " << i << "\t" << j << "\t" << errorsOnDiagonalMatrix(i,j) << std::endl;       }
    }}
  

  gStyle->SetPalette(1);
  // manipulate erros and matrices to get correlation matrix
  TCanvas * theErrorsCanvas = new TCanvas("theErrors^2 from unfolded results","theErrors^2 from unfolded results",250,250,1050,850);  
  theErrorsCanvas->cd();
  errorsOnDiagonalMatrix.Draw("colz");
  errorsOnDiagonalMatrix.Draw("text same");

  TCanvas * theUnfoldingCanvas = new TCanvas("the Unfolding matrix Canvas","the Unfolding matrix Canvas",210,210,1000,800);  
  theUnfoldingCanvas->cd();
  unfoldingMatrix.Draw("colz");

  TCanvas * theDebuggingCanvas = new TCanvas("the debugging Canvas","the debugging Canvas",280,280,1080,880);  
  theDebuggingCanvas->cd();

  std::cout << "sizes:   errorsOnDiagonalMatrix is " << errorsOnDiagonalMatrix.GetNcols() << "\t transpUnfoldingMatrix is: " << transpUnfoldingMatrix.GetNcols() << std::endl;
  TMatrixD tmp(errorsOnDiagonalMatrix,TMatrixD::kMult,transpUnfoldingMatrix);
  tmp.Draw("colz");
  TMatrixD covarianceMatrix(unfoldingMatrix,TMatrixD::kMult,tmp);
  TCanvas * theCovarianceCanvas = new TCanvas("theCovarianceCanvas","theCovarianceCanvas",200,200,1000,800);
  theCovarianceCanvas->cd();
  covarianceMatrix.Draw("colz");
  //covarianceMatrix.Draw("text same");

  
  // write out text output
  // by default systematic errors are added in quadrature to the diagonal
  // of the covariance matrix, which has been obtained using only stat errors
  ofstream myfile;
  myfile.open ("covarianceTextFile.txt");

  std::cout << "\n\n\n covariance matrix (num errors in file) " << errorsCounter << " x " << errorsCounter << std::endl;
  std::cout << " unfolding  matrix (from .root) " << unfoldingMatrix.GetNcols() << " x " << unfoldingMatrix.GetNrows() << std::endl;
  
  bool addSystematicsOnDiagonal(false);
  if(addSystematicsOnDiagonal) std::cout << " adding systematic errors in quadrature to the diagonal of the cov matrix\n" << std::endl;
  else                         std::cout << " NOT adding systematic errors in quadrature to the diagonal of the cov matrix\n" << std::endl;
  for(int i=0; i<errorsCounter; i++){
    for(int j=0; j<errorsCounter; j++){
      //std::cout << (i+1) << "\t" << (j+1) << "\t" << covarianceMatrix(i,j) << std::endl;     
      if(addSystematicsOnDiagonal) {
	float theMatrixElement = covarianceMatrix(i,j);
	if(i==j) theMatrixElement+= systErrosArray[i]*systErrosArray[i];
	  myfile << scientific << (i+1) << "\t" << (j+1) << "\t" << theMatrixElement << std::endl;
	  if(i==j) std::cout << "cov matrix diagonal (sig^2): " << theMatrixElement << "\t sig: " << sqrt(theMatrixElement) << std::endl;
      }
      else{
	myfile << scientific << (i+1) << "\t" << (j+1) << "\t" << covarianceMatrix(i,j) << std::endl;
	if(i==j) std::cout << "cov matrix diagonal (sig^2): " << covarianceMatrix(i,j) << "\t sig: " << sqrt(covarianceMatrix(i,j)) << std::endl;
      }
    }// loop on matrixIndex
  }// loop on matrixIndex
 
  myfile.close();
  
  // export formulas from here:
  // http://en.wikipedia.org/wiki/Propagation_of_uncertainty#Linear_combinations

}




int readUnfoldingMatrices(std::string unfoldingMatrixFileInputFile, std::string effAccFileForTest, std::string variable) 

{

  bool isRapidity(true);
  if (variable==std::string("Y")){
	std::cout << "\n\tyou have chosen Y" << std::endl;
  }  
  else if (variable==std::string("pt")){
	std::cout << "\n\tyou have chosen pt" << std::endl;
	isRapidity=false;
  }
  else{
	std::cout << "not done a valid choice of variable. Bailing out" << std::endl;
	return 0;
  }


  gStyle->SetOptStat(110011);
  gStyle->SetPaintTextFormat(".2f");

  std::cout << "\n\n\tReading test histogram  from: " << effAccFileForTest          << std::endl;
  std::cout << "\tReading unfolding matrices in:  "   << unfoldingMatrixFileInputFile <<"\n\n"<< std::endl;

  TFile *theEffAccInpuntFile = new TFile(effAccFileForTest.c_str(),"read");

  TFile *theunfoldingMatrixInputFile = new TFile(unfoldingMatrixFileInputFile.c_str(),"read");
  theunfoldingMatrixInputFile->cd();
  std::cout << "The file: " << theunfoldingMatrixInputFile->GetName() << " has been opened"<< std::endl;
    
  // this is the matrix you want to use for to unfold
  TMatrixD * theUnfoldingMatrix = (TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal");


  //////////////////////////////////////////////////////
  // to be consistent with the definitions 
  // of AN-10-367 one needs to use the YZTL_vs_YZ before Eff and Acc
  // are applied ==> move to the "All" folder

  std::string prefix("mcEff/");
  std::string definition=std::string("All");
  std::string cut("");
  std::string plot("ciao");
  if(isRapidity){
    plot=std::string("Z0_Y"); }// smeared Y before Eff&Acc
  else{
    plot=std::string("Z0_Pt"); 
  }
  std::string plotName = prefix + definition;
  plotName = plotName + std::string("/");    plotName = plotName + cut;
  plotName = plotName + std::string("/");    plotName = plotName + plot;
  if(isRapidity){
    std::cout << "\tUsing directly the total smeared Y plot:     "<< plotName << std::endl;}
  else {
    std::cout << "\tUsing directly the total smeared Pt plot:     "<< plotName << std::endl;  }
  TH1F * h_fullTest_smeared = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());


  plotName = prefix + definition;
  plotName = plotName + std::string("/");    plotName = plotName + cut;
  if(isRapidity){
    plot=std::string("Z0_YTL"); // tree-level Y before Eff&Acc
    std::cout << "\tComparing against the total tree-Level Y plot:     "<< plotName << std::endl;}
  else{
    plot=std::string("Z0_PtTL"); 
    std::cout << "\tComparing against the total tree-Level Pt plot:     "<< plotName << std::endl;  }
  plotName = plotName + std::string("/");    plotName = plotName + plot;


  TH1F * h_fullTest_treeLevel = (TH1F*) theEffAccInpuntFile->Get(plotName.c_str());

  TH1D * h_RapidityUNSmeared;
  h_RapidityUNSmeared = (TH1D*) unfold(h_fullTest_smeared, theUnfoldingMatrix );

  gStyle->SetOptStat(0);
  TCanvas * theBaseCanvas = new TCanvas("check of closure","check of closure",200,200,1000,800);  
  theBaseCanvas->cd();
  if(!isRapidity){
   theBaseCanvas->SetLogx();   theBaseCanvas->SetLogy();
  }
  
  h_RapidityUNSmeared->Draw();

  h_fullTest_smeared->SetLineColor(kGreen);
  h_fullTest_smeared->Draw("same");

  h_fullTest_treeLevel->SetLineColor(kRed);
  h_fullTest_treeLevel->SetLineStyle(8);
  h_fullTest_treeLevel->Draw("same");
  
  TLegend* leg = new TLegend(0.7,0.7,0.95,0.9);
  leg->SetHeader("closure: isn't it?");
  leg->AddEntry(h_fullTest_treeLevel,"tree level","lp");
  leg->AddEntry(h_RapidityUNSmeared,"unsmeared: does it overlap?","lp");
  leg->AddEntry(h_fullTest_smeared,"smeared","lp");
  leg->SetFillColor(0); leg->SetShadowColor(0); leg->SetTextSize(0.02);
  leg->Draw();  
  

  std::cout << "\n\n\n\t If you see tree plots and two are overlayed (red and black), you've applied unfolding matrix succesfully" << std::endl;

  return 0;
}
