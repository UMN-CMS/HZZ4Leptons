#include <fstream>
#include <iomanip>

#include <string.h>
#include <vector>
#include <iostream>
#include <algorithm>

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TKey.h>
#include <TText.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TMath.h>

#include "ZShape/EffAcc/interface/SignalBackgroundExtraction.h"
#include "ZShape/Base/interface/EffTableLoader.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include <boost/algorithm/string/replace.hpp>
#include "PhysicsTools/TagAndProbe/interface/ZGeneratorLineShape.h"

SignalBackgroundExtraction::SignalBackgroundExtraction()
{
  theFile_ =0;
  sigefftable_ =0;
  bkgefftable_ =0;

  doInit();

  sigtextFileName_="";
  bkgtextFileName_="";

  
}


SignalBackgroundExtraction::SignalBackgroundExtraction(TFile * file, std::string EffName, std::string EffType, bool doSig, bool doBkg, bool doSimpleOut , double masslo, double masshi)
{ 
  doInit();
  effName_ = EffName;
  effType_ = EffType;
  doSig_ = doSig;
  doBkg_ = doBkg;
  doSimpleOut_ = doSimpleOut;
  sigtextFileName_="";
  bkgtextFileName_="";
  masslo_ = masslo;
  masshi_ = masshi;
  std::cout << "class SignalBackgroundExtraction created with root file: " << file->GetName() << std::endl;
  isPt_ = strstr(EffType.c_str(),"pt");
  isBoth_ = (bool (strstr(EffType.c_str(),"pt")) && bool (strstr(EffType.c_str(),"rap")));
  isRap_ = strstr(EffType.c_str(),"rap");
  setRootFile(file);

}

SignalBackgroundExtraction::SignalBackgroundExtraction(const std::string & sigtextFileName  , const std::string & bkgtextFileName )
{ 
  theFile_ =0;
  doInit();
  if ( !(sigtextFileName == "0")) {doSig_ = true; }
  if ( !(bkgtextFileName == "0")) {doBkg_ = true; }

  if ( doSig_)
    {
      sigtextFileName_ = sigtextFileName;
      sigefftable_ = new EffTableLoader(sigtextFileName_);
      std::cout << "class SignalBackgroundExtraction created with Signal text file: " << sigtextFileName_ << std::endl;
    }
  if ( doBkg_)
    {
      bkgtextFileName_ = bkgtextFileName;
      bkgefftable_ = new EffTableLoader(bkgtextFileName_);
      std::cout << "class SignalBackgroundExtraction created with Background text file: " << bkgtextFileName_ << std::endl;
    }
  
  produceHistograms(0,0);
  produceHistograms(0,1);

}


void SignalBackgroundExtraction::doInit(void)
{
  values1DHisto_=0;   systPlus1DHisto_=0;
  systMinus1DHisto_=0;  denominator1DHisto_=0;
  values2DHisto_=0;  systPlus2DHisto_=0;
  systMinus2DHisto_=0; denominator2DHisto_=0;
  doSig_ = false;
  doBkg_ = false;
}
  

//-----------------------------------------------------------------------------------//
  
void SignalBackgroundExtraction::setRootFile(TFile * file){
  
  theFile_ = file;
  if ( (!theFile_) &&
       theFile_->IsOpen() == kFALSE
       )
    {
      std::cout << "invalid root file. Returning." << std::endl;
    }

  std::cout<<"object SignalBackgroundExtraction opened file: " << theFile_->GetName() << std::endl;
    
    
  // getting list of keys of all objects in the file
  TIter next(theFile_->GetListOfKeys());
  TKey * key;
    
  values1DHisto_=0;     systPlus1DHisto_=0;
  systMinus1DHisto_=0;  denominator1DHisto_=0;
    
  values2DHisto_=0;     systPlus2DHisto_=0;
  systMinus2DHisto_=0;  denominator2DHisto_=0;
    
  
  // looping on objects in the file aiming at the 4 needed histograms 
    
  while((key=((TKey*) next()) )){
    if (strcmp(key->GetClassName(),"TDirectoryFile")==0)
      {
	TDirectory *mdir = (TDirectory*) theFile_->Get(key->GetName());
	TIter next1(mdir->GetListOfKeys());
	TKey * key1;
	while((key1=((TKey*) next1())) ){
	  //at this step I need to check if I want eta or pt or whatever...
	  //std::cout << " key1 " << key1->GetName() << " type " << key1->GetClassName() << std::endl;
	  if (strcmp(key1->GetClassName(),"TDirectoryFile")==0 && strcmp(key1->GetName(),effType_.c_str())==0 )
	    {	    
	      TDirectory *mdir1 = (TDirectory*) mdir->Get(key1->GetName());
	      TIter next2(mdir1->GetListOfKeys());
	      TKey * key2;
	      
	      while((key2=((TKey*) next2())) ){
		//std::cout << " key2 " << key2->GetName() << " type " << key2->GetClassName() << std::endl;
		if (strcmp(key2->GetClassName(),"TDirectoryFile")==0 && strstr(key2->GetName(),"bin") )
		  {
		    TDirectory *mdir2 = (TDirectory*) mdir1->Get(key2->GetName());
		    TIter next3(mdir2->GetListOfKeys());
		    TKey * key3;
		    
		    while((key3=((TKey*) next3())) ){
		      if ( strcmp(key3->GetName(), "fit_canvas")==0 )
			{
			  std::string dirname = effName_ + key2->GetName();
			  TCanvas *tcan = (TCanvas* ) mdir2->Get(key3->GetName());
			  tcan->cd();
			  TCanvas *myc2 = new TCanvas("myc2","myc2", 1200, 800);
			  myc2->cd();

			  TPad *mypad = (TPad*)tcan->GetPrimitive("fit_canvas_3");
			  TIter mnext (mypad->GetListOfPrimitives());
			  TObject *mkey;
			  while((mkey=((TObject*) mnext()) )){
			    if (strcmp(mkey->ClassName(),"TH1D")==0)
			      {	
				((TH1D*) mypad->GetPrimitive(mkey->GetName()))->SetTitle("");
				((TH1D*) mypad->GetPrimitive(mkey->GetName()))->GetXaxis()->SetTitle("m_{Z} (GeV/c^{2})");
				((TH1D*) mypad->GetPrimitive(mkey->GetName()))->GetXaxis()->CenterTitle(1);
				((TH1D*) mypad->GetPrimitive(mkey->GetName()))->GetYaxis()->CenterTitle(1);
				((TH1D*)mypad->GetPrimitive(mkey->GetName()))->GetYaxis()->SetRangeUser(0.01,75);
				
			      }
			  }
			  TPad *mypl =   (TPad*)tcan->GetPrimitive("fit_canvas_4");
 			  TPaveText *mypav = (TPaveText*) mypl->GetPrimitive("simPdf_paramBox");

			  TPaveText *nmypav = new TPaveText(0.6,0.78,0.90,0.90,"brNDC");
			  
			  TText *l1 = nmypav->AddText(mypav->GetLine(2)->GetTitle());
			  TText *l2 = nmypav->AddText(mypav->GetLine(1)->GetTitle());
			  nmypav->SetLineColor(kBlack);
			  nmypav->SetTextColor(kBlack);
			  nmypav->SetFillColor(kWhite);

			  myc2->cd();
			  mypad->Draw();
 

			  mypad->SetPad(0,0,1,1);
			  mypad->cd();
			  mypad->SetLogy(1);
			  mypad->SetTitle("");
			  nmypav->Draw();
			  myc2->Modified();
			  myc2->Update();
			  mypad->Modified();


			  myc2->Print(Form("%s.png",dirname.c_str()));
			  myc2->Print(Form("%s.eps",dirname.c_str()));
			  //This then prints out all the cute little plots... or it should.
			}
		      //std::cout << " key3 " << key3->GetName() << " type " << key3->GetClassName() << std::endl;
		      //here is where I get the w's!!!
		      if ( strcmp(key3->GetName(), "w")==0 )
			{
			  std::cout << " Found " << key3->GetName() << " in " << key2->GetName() << " in "<< key1->GetName() << " in "<< key->GetName() << std::endl;
			  std::string dirname = key2->GetName();
			  std::replace(dirname.begin(),dirname.end(),'_',' ');
			  char myvartype[1000];
			  char restofinfo[1000];
			  int rbinval=-1,pbinval=-1;

			  sscanf(dirname.c_str(), "pair %s bin%d %s",myvartype,  &rbinval, restofinfo);
			  //std::cout << " bin " << binval << " vartype " << myvartype << " rest " << restofinfo << std::endl;
			  if (isPt_ && !isBoth_) {std::swap(pbinval,rbinval); }
			  char nrestofinfo[1000];
			  std::string newstr = Form("pair_rapidity_%s_bin%d_%%s",myvartype,rbinval);

			  sscanf(key2->GetName(), newstr.c_str(), nrestofinfo);

			  std::string darest = nrestofinfo;
			  std::replace(darest.begin(),darest.end(),'_',' ');
			  //std::cout << " nrest of info " << nrestofinfo << " darest " << darest << std::endl;
			  int obinval = -1 ;

			  char omyvartype[1000];
			  char orestofinfo[1000];
			  int good = sscanf(darest.c_str(), " pair %s bin%d %s",omyvartype,  &obinval, orestofinfo);
			  //if (good) std::cout << " newstuff obinval " << obinval << " " << omyvartype << std::endl;
			  if (!good) obinval = -1;
			  else pbinval = obinval;
			  //if ( strstr(myvartype,"eta") ) std::swap(binval,obinval);
			  
			  
 			  RooWorkspace* w = (RooWorkspace*) mdir2->Get(key3->GetName());
			  //w->Print();
			  RooRealVar *var1  = (isRap_) ? w->var("pair_rapidity") : w->var("pair_pt") ; //dear anyone reading this, it is EITHER et or eta!!
			  RooRealVar *var2 = isBoth_ ? w->var("pair_pt"): 0 ;
			  //Here I will one day need to do 'two-D'... that is why I have the bin1 and bin2
			  RooRealVar *eff = w->var("efficiency");
			  RooRealVar *nS  = w->var("numSignalAll");
			  if (! var1 ) std::cout << " no var1? " << std::endl;
			  if (! eff ) std::cout << " no eff? " << std::endl;
			  if (! nS ) std::cout << " no nS? " << std::endl;
			  if (! nS || ! var1 ) continue; 
			  int hbinnum =var1->getBinning().binNumber(var1->getVal());
			  Double_t hval = var1->getBinning().binHigh(hbinnum);
			  Double_t lval = var1->getBinning().binLow(hbinnum);

			  if (isBoth_)
			  {
			    int hbinnumeta =var2->getBinning().binNumber(var2->getVal());
			    Double_t hvaleta = var2->getBinning().binHigh(hbinnumeta);
			    Double_t lvaleta = var2->getBinning().binLow(hbinnumeta); 
			    bin1Min_.push_back( lval );     bin1Max_.push_back( hval);
			    bin2Min_.push_back( lvaleta );  bin2Max_.push_back( hvaleta );
			  }
			  else 
			  {
			    if (isPt_) {bin2Min_.push_back( lval );  bin2Max_.push_back( hval );}
			    else { bin1Min_.push_back( lval );  bin1Max_.push_back( hval);}			  
			  }

			  ptbin_.push_back(pbinval);
			  rapbin_.push_back(rbinval);

			  int mappos =  rapbin_.size()-1;
			  //systPlus1D_.push_back(errorHi);
			  //systMinus1D_.push_back(eff->getErrorLo());
			  std::cout << " mass low " << masslo_ << " mass hi " << masshi_ << std::endl;
			  w->var("mass")->setRange("tight",masslo_,masshi_);
			  
			  //Background
			  RooRealVar *numBackgroundPass = w->var("numBackgroundPass");
			  Double_t NBPScale =  w->pdf("backgroundPass")->createIntegral(*w->var("mass"),RooFit::NormSet(*w->var("mass")),RooFit::Range("tight"))->getVal();
			  Double_t NBP   = numBackgroundPass->getVal()*NBPScale;
			  Double_t NBPLo = numBackgroundPass->getErrorLo()*NBPScale;
			  Double_t NBPHi = numBackgroundPass->getErrorHi()*NBPScale;

			  Double_t NSPScale =  w->pdf("signal")->createIntegral(*w->var("mass"),RooFit::NormSet(*w->var("mass")),RooFit::Range("tight"))->getVal();
			  Double_t NSP   = nS->getVal()* NSPScale;
			  Double_t NSPLo = nS->getErrorLo()*NSPScale;
			  Double_t NSPHi = nS->getErrorHi()*NSPScale;
			  //nS->getVal()
			  std::cout << "Original values " << nS->getVal() << " bkg " <<  numBackgroundPass->getVal() << std::endl;
			  std::cout << " scales " << NBPScale << " sig " << NSPScale << std::endl;
			  std::cout << "Values are " << NSP << " back " << NBP << " hi " << NSPHi << " bhi " << NBPHi << std::endl;
			  values1D_.push_back( NSP );
			  denominator1D_.push_back( NBP );

			  systPlus1D_.push_back(NSPHi);
			  systMinus1D_.push_back(NBPHi);
			  
			
			}
		    }
		  }
	      }
	    }
	}
      }
  }
    
  vecToHist();//So one can use the histograms if one wishes...

}


//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::vecToHist(void){
  //Simple internal program converts the vector quantities into the histogram ones.
  //Maybe I should keep some 'bools' around to test if the vec's or hist's are around.....

  int binCounter = values1D_.size();

  efficiencyName_ = "empty";

  std::string title = std::string("values_")+efficiencyName_ + "EtaDet_Pt"; /// JUST A HACK for the MOMENT 
  values1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicMinus_")+efficiencyName_;
  systMinus1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicPlus_")+efficiencyName_;
  systPlus1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("denominator_")+efficiencyName_;
  denominator1DHisto_ = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);

  for (int i=0; i<binCounter; ++i)
    {
      values1DHisto_->SetBinContent(i, values1D_[i] );
      systPlus1DHisto_->SetBinContent(i, systPlus1D_[i]);
      systMinus1DHisto_->SetBinContent(i, systMinus1D_[i]);
      denominator1DHisto_->SetBinContent(i, denominator1D_[i]);
    }
  
}

//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::histToVec(void){
  
}


//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::setVarMinMax(double vmin, double vmax){
  vmax_=vmax;
  vmin_=vmin;
  if (isBoth_) return;
  for (uint j = 0; j < values1D_.size(); ++j)
    {
      if (isPt_) 
	{
	  bin1Min_.push_back(vmin);
	  bin1Max_.push_back(vmax);
	}
      else
	{
	  bin2Min_.push_back(vmin);
	  bin2Max_.push_back(vmax);
	}
    }
}
  
  

//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::produceTxtFile(std::string &textFileName, std::string sbtype){
  /*    
  // check if there are 4 histograms of the same type and with the same number of bins;
  bool histosAre1D =false;
  if(  ( values1DHisto_ && systPlus1DHisto_ &&  systMinus1DHisto_  &&  denominator1DHisto_ ) ) histosAre1D = true;
  if(  ( values1DHisto_ && denominator1DHisto_ ) ) histosAre1D = true; //To allow for now systematic histos
  bool histosAre2D =false;
  if(  ( values2DHisto_ && systPlus2DHisto_ &&  systMinus2DHisto_  &&  denominator2DHisto_ ) ) histosAre2D = true;
  if(  ( values2DHisto_ &&  denominator2DHisto_ ) ) histosAre2D = true; //To allow for now systematic histos

  if ( (!histosAre1D)   && (!histosAre2D))
    {
      std::cout << "there need be 4 1-d or 2-d histograms which were not found in file: " 
		<< theFile_->GetName() << ". Text file cannot be produced" << std::endl;
      return;
    }



  textFileName_ = textFileName;

  if(histosAre1D ) produceTxtFile1D();
  if(histosAre2D ) produceTxtFile2D();
  */
  textFileName_ = textFileName;
  produceTxtFile1D(sbtype);
}


  
//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::produceHistograms(std::string rootFileName){
  TFile * rootFile = new TFile(rootFileName.c_str(),"RECREATE");
  produceHistograms(rootFile,0);
  produceHistograms(rootFile,1);

}    


  
//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::produceHistograms(TFile * rootFile, bool issig){
    
  if ( issig && !doSig_ )  {std::cout << " no signal set"     << std::endl; return;}
  if ( !issig && !doBkg_ ) {std::cout << " no background set" << std::endl; return;}

  std::string textFileName = "0";
  std::vector<std::string> nameVec;

  if ( issig ) {textFileName = sigtextFileName_;nameVec = signameVec_; } 
  else {textFileName = bkgtextFileName_;nameVec = bkgnameVec_;}

  // check if the text file is available as needed for input
  if (textFileName == ""){
    std::cout << "no inpunt txt file give, cannot produce histograms." <<std::endl;
    return;
  }
 
  dimension_=0;
  char str [1000];    char effName[1000];
  ifstream is(textFileName.c_str());

  while ( 	is.getline(str,1000)){
    // std::cout << str <<"\n";
    if ( strstr(str,"#dimension:") )           sscanf(str, "#dimension:%d",  &dimension_);
    if ( strstr(str,"#efficiency name:") )  sscanf(str, "#efficiency name: %s",  effName);
    if ( strstr(str,"#Vars:") )
      {
	char * tstr;	
	tstr = strtok (str," ,.-");
	while (tstr != NULL)
	  {
	    if ( strstr(tstr,"Vars")) {	 tstr = strtok (NULL, " ,.-");  continue;}
	    printf ("%s\n",tstr);
	    nameVec.push_back(tstr);	 
	    tstr = strtok (NULL, " ,.-");
	    
	  }
	  
      }
  } 
  is.close();
  efficiencyName_ = std::string(effName);
  if ( issig ) { signameVec_ = nameVec; } 
  else {bkgnameVec_=nameVec;}

  // call dedicated private method depending on the dimension of the histogram
  //  std::cout<<"histogram dimension is: " << dimension_ << std::endl;
  if(dimension_==1) produce1DHistograms(rootFile, issig);
  //  else if(dimension_==2) produce2DHistograms(rootFile);
  else  {
    std::cout<<"*** histogram dimension is: " << dimension_ << " but only 1 and 2 are valid.\n" <<std::endl;
  }
    
  return;

}

  


//-----------------------------------------------------------------------------------//
void SignalBackgroundExtraction::produceTxtFile1D(std::string sbtype){


  varmap MyVars_; //I probably should be using pointers here... but oh well :)
  std::vector<std::string> nameVec;
  if (sbtype == "signal" && doSig_ ) {MyVars_ = SigVars_;nameVec = signameVec_; }
  else if (sbtype == "background" && doBkg_ ) {MyVars_ = BkgVars_;nameVec = bkgnameVec_;}
  else {std::cout << " There are no " << sbtype << " parameters set " << std::endl; return;}

  ofstream the1DEffFile;
  the1DEffFile.open (Form("%s%s",sbtype.c_str(),textFileName_.c_str()),std::ios::out);


  the1DEffFile << "#put your comments after a '#'.\n";

  if ( doSimpleOut_)
    {
      the1DEffFile << "#Variable name: " << effName_  <<"\n";
      the1DEffFile << "#" << std::setw(10) << "RapBin"
		   << " " << std::setw(10) << "RapMin"
		   << " " << std::setw(10) << "RapMax"
		   << " " << std::setw(10) << "Value"
		   << " " << std::setw(10) << "Error\n";
      if (the1DEffFile.is_open())
	{
	  
	  for (uint r=0; r< bin1Min_.size(); r++) 
	    {
	      bool doem = false ;
	      if ( sbtype == "signal" && doSig_) doem = true;
	      if ( !doem) if ((!doBkg_) || (!(sbtype == "background"))) continue;
	      double myval = ((doem) ? values1D_[r] :  denominator1D_[r]);
	      //std::vector<std::pair<float, float> > EffInfo = efftable_->GetCellInfo(r);
	      the1DEffFile << " " << std::setw(10) << std::setprecision(4) <<  ((isRap_) ? rapbin_[r]+1 : ptbin_[r]+1)
			   << " " << std::setw(10) << std::setprecision(4) <<  ((isRap_) ? bin1Min_[r] :bin2Min_[r] )  
			   << " " << std::setw(10) << std::setprecision(4) <<  ((isRap_) ? bin1Max_[r] :bin2Max_[r] )
			   << " " << std::setw(10) << std::setprecision(4) <<  myval
			   << " " << std::setw(10) << std::setprecision(4) <<  ((!doem) ? systMinus1D_[r] :systPlus1D_[r])  << "\n";
	    }
	}
    }
  else
    {
      if ( MyVars_.size() < 1 ) {std::cout << " no variables " << std::endl; return;}
      int totalvars = MyVars_[0].size()*4 +3;
      the1DEffFile << "#efficiency name: " << effName_  <<"\n";
      the1DEffFile << "#Vars:"; 
      for ( uint u = 0; u < nameVec.size(); ++u)
	{
	  the1DEffFile << " " << nameVec[u]; 
	}
      the1DEffFile << "\n";

      the1DEffFile << "#dimension: 1 \n\n\n";
      the1DEffFile << " " << std::setw(10) << "#PtMin"
		   << " " << std::setw(10) << "PtMax"
		   << " " << std::setw(10)  << "EtaMin"
		   << " " << std::setw(10)  <<  "EtaMax"
		   << " " << std::setw(10)  << "NumParms"
		   << " " << std::setw(10) << "ptbin"
		   << " " << std::setw(10) << "etabin"
		   << " " << std::setw(10) << "rapbin"
	//<< " " << std::setw(10)  << "parname"
		   << " " << std::setw(10) << "parval"
		   << " " << std::setw(10) << "parplus"
		   << " " << std::setw(10) << "parminus"
		   << " " << std::setw(10) << "parconst"
	
		   << "\n";
     
      if (the1DEffFile.is_open())
	{
	  
	  for (uint r=0; r< bin1Min_.size(); r++) 
	    {
	      //std::vector<std::pair<float, float> > EffInfo = efftable_->GetCellInfo(r);
	      the1DEffFile << " " << std::setw(10) << std::setprecision(4) <<  bin1Min_[r]  
			   << " " << std::setw(10) << std::setprecision(4) <<  bin1Max_[r]
			   << " " << std::setw(10) << std::setprecision(4) <<  bin2Min_[r]
			   << " " << std::setw(10) << std::setprecision(4) <<  bin2Max_[r]
			   << " " << std::setw(10) << std::setprecision(4) <<  totalvars
			   << " " << std::setw(10) << std::setprecision(4) <<  ptbin_[r]
			   << " " << std::setw(10) << std::setprecision(4) <<  etabin_[r]
			   << " " << std::setw(10) << std::setprecision(4) <<  rapbin_[r];

	      for ( varit myit = MyVars_[r].begin(); myit != MyVars_[r].end(); ++myit)
		{
		  the1DEffFile //<< " " << std::setw(10) <<  (*myit).first
		    << " " << std::setw(10) << std::setprecision(4) <<  (*myit).second.val
		    << " " << std::setw(10) << std::setprecision(4) <<  (*myit).second.max
		    << " " << std::setw(10) << std::setprecision(4) <<  (*myit).second.min
		    << " " << std::setw(10) << std::setprecision(4) <<  (*myit).second.con;
		  
		}
	      the1DEffFile << "\n";
	    }
	}
    }
  the1DEffFile.close();
    
}

  
  

void SignalBackgroundExtraction::produce1DHistograms(TFile * rootFile, bool issig){

  if ( issig && !doSig_ )  {std::cout << " no signal set"     << std::endl; return;}
  if ( !issig && !doBkg_ ) {std::cout << " no background set" << std::endl; return;}

  EffTableLoader *efftable = 0;
  varmap MyVars_;
  std::vector<std::string> nameVec;

  if ( issig ) {efftable = sigefftable_;nameVec = signameVec_; }
  else  {efftable =bkgefftable_;nameVec = bkgnameVec_;}


  int binCounter = efftable->size();

  if (rootFile!=0) rootFile->cd();
  std::string title = std::string("values_")+efficiencyName_ + "EtaDet_Pt"; /// JUST A HACK for the MOMENT 
  TH1F * valuesHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicMinus_")+efficiencyName_;
  TH1F * systematicMHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("systematicPlus_")+efficiencyName_;
  TH1F * systematicPHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);
  title = std::string("denominator_")+efficiencyName_;
  TH1F * denominatorHisto = new TH1F(title.c_str(),title.c_str(), binCounter, 0, binCounter);

  if (rootFile==0) {
    valuesHisto->SetDirectory(0);
    systematicMHisto->SetDirectory(0);
    systematicPHisto->SetDirectory(0);
    denominatorHisto->SetDirectory(0);
    values1DHisto_=valuesHisto;
    systPlus1DHisto_=systematicMHisto;
    systMinus1DHisto_=systematicPHisto;
    denominator1DHisto_=denominatorHisto;
  }


  for (int u=0; u<binCounter; u++)
    {  
      //HERE IS THE MEAT OF THE EFF FIT FROM TXT

      std::vector<std::pair<float, float> > BinInfo = efftable->GetCellInfo(u);
      std::vector<float> BinParms = efftable->correctionEff(u);
      int binsize =  BinParms.size();
      bin1Min_.push_back(BinInfo[0].first);
      bin1Max_.push_back(BinInfo[0].second);
      bin2Min_.push_back(BinInfo[1].first);
      bin2Max_.push_back(BinInfo[1].second);
      ptbin_.push_back(BinParms[1]);
      etabin_.push_back(BinParms[2]);
      rapbin_.push_back(BinParms[3]);

      //std::cout << " bin1 " << BinInfo[0].first << " " << BinInfo[0].second << " " << BinInfo[1].first << " " <<  BinInfo[1].second << " " << BinParms[0] << " " <<BinParms[1] << std::endl; 
	    
      for ( int b = 0; b < (binsize-3)/4; ++b)
	{

	  Vars_t davar;
	  davar.val    = BinParms[b*4+4];
	  davar.max    = BinParms[b*4+5];
	  davar.min    = BinParms[b*4+6];
	  davar.con    = BinParms[b*4+7];
	  MyVars_[u][nameVec[b]]=davar;
	  //std::cout << " namevec " << nameVec[b] <<" val " <<davar.val <<std::endl;
	}

    }

  if ( issig ) {SigVars_=MyVars_;}
  else  { BkgVars_=MyVars_;}

  if (rootFile!=0) {
    valuesHisto->Write();systematicMHisto ->Write();
    systematicPHisto->Write(); denominatorHisto->Write();
    rootFile->Close();
  }
}
  

//void produce2DHistograms(TFile * rootFile){
//  std::cout<<"produce2DHistograms not yet implemented" << std::endl;
//}
  
void SignalBackgroundExtraction::produceTxtFile2D(){
  std::cout<<"produceTxtFile2D not yet implemented" << std::endl;
}

void SignalBackgroundExtraction::produceConfig(std::string myconfig, std::string mytemplate,std::string sbtype){
  std::cout<<"Producing Configuration File" << myconfig << std::endl;
  
  //I need to write an output condifugration file based on the template
  // the output needs to match the destination (Python) format 
  // it also needs to be in the SignalBackgroundExtractionting
  
  // DETAILS

  //LOOP over the different variables
  // for each set the param and or the errors based on the templates.
  
  //Different cases... one with errors... one without? 
  // Start with the no error case for now
  //
  varmap MyVars_; //I probably should be using pointers here... but oh well :)
  bool doErrors = false;
  std::vector<std::string> nameVec;
  if (sbtype == "signal" && doSig_ ) {MyVars_ = SigVars_;nameVec = signameVec_; }
  else if (sbtype == "background" && doBkg_ ) {MyVars_ = BkgVars_;nameVec = bkgnameVec_;}
  else {std::cout << " Missing something signal or background sbtype " << sbtype << " doSig_ " << doSig_ << " doBkg_ " << doBkg_ << std::endl; return;}
  
  //OPEN THE OUTPUT FILE
  ofstream theConfFile;
  theConfFile.open (myconfig.c_str(),std::ios::out);

  std::string mapconf = "map"+myconfig;
  ofstream theMapConfFile;
  theMapConfFile.open (mapconf.c_str(),std::ios::out);
  std::string bintopdfmap="";


  for (uint r=0; r< bin1Min_.size(); r++) 
    {
      std::string etastr="", ptstr="", rapstr="",bincomb="";
      if (ptbin_[r] > -1 ) ptstr=Form("pt%d",ptbin_[r]);  
      if (etabin_[r] > -1 ) etastr=Form("eta%d",etabin_[r]);
      if (rapbin_[r] > -1 ) rapstr=Form("rap%d",rapbin_[r]);
      bincomb=ptstr+etastr+rapstr;

      std::string metastr="", mptstr="",mbincomb="";
      if (ptbin_[r] > -1 ) mptstr=Form("et_bin%d",ptbin_[r]);   
      if (etabin_[r] > -1 ) metastr=Form("eta_bin%d",etabin_[r]);
      if (rapbin_[r] > -1 ) metastr=Form("rapidity_bin%d",rapbin_[r]);
      mbincomb="*"+mptstr+metastr+"}*"; //This is the magic needed
      if ( ptbin_[r] > -1 && etabin_[r] > -1)  mbincomb="*"+mptstr+";*"+metastr+"}*"; //This is the magic needed...

      std::string mbinname = "";

      //Open the input file
      ifstream is(mytemplate.c_str());
      char str [1000];
      while ( 	is.getline(str,1000)){
	// std::cout << str <<"\n";
	std::string mystr = str; //me likes strings, easy to use predefined algos to manipulate!

	//find the xxbinxx thing first, and take the first part
	std::string::size_type pos = mystr.find("xxbinxx");
	if (pos != string::npos) 
	  {
	    mbinname = mystr.substr(0,pos)+ bincomb;
	    if ( r == 0 ) bintopdfmap += "\""+mystr.substr(0,pos)+"\"";
	  }

	//first, always replace the bin value
	boost::algorithm::replace_all(mystr,"xxbinxx",bincomb);

	for ( varit myit = MyVars_[r].begin(); myit != MyVars_[r].end(); ++myit)
	{
	  std::string thevar = Form("xx%sxx",((*myit).first).c_str());
	  float val = (*myit).second.val;
	  float min = (*myit).second.min+val;
	  float max = (*myit).second.max+val;
	  std::string theval = (doErrors && !((*myit).second.con)  )?  Form("%f,%f,%f",val,min,max ): Form("%f",val);
	  boost::algorithm::replace_all(mystr,thevar,theval);
	  
	}
	theConfFile << mystr << "\n";
	      

      }
      bintopdfmap += ",\""+mbincomb + "\",\""+mbinname+ "\"" ;
      //std::cout << " bin to pdfmap " << bintopdfmap << std::endl;
    }
  theConfFile.close();
  theMapConfFile << bintopdfmap << std::endl;
  theMapConfFile.close();
    
  /*
   voigtPlusBKGxxbinxx = cms.vstring(
             "Voigtian::signal(mass, mean[91.2,85,96], width[2.495], sigma[2.3,1,9])",
             "EXPR::backgroundPass('(RooMath::erfc((alpha - mass) * beta))*(( (mass-peak)*gamma < -70  )* (1e20) + ( (mass-peak)*gamma > 70 )* (0.0) + ( (mass-peak)*gamma >= -70 && (mass-peak)*gamma <= 70) * exp(-(mass-peak)*gamma))',mass,peak[91.1876],gamma[xxgammaxx],beta[xxbetaxx],alpha[xxalphaxx])",
             "EXPR::backgroundFail('(RooMath::erfc((alpha - mass) * beta))*(( (mass-peak)*gamma < -70  )* (1e20) + ( (mass-peak)*gamma > 70 )* (0.0) + ( (mass-peak)*gamma >= -70 && (mass-peak)*gamma <= 70) * exp(-(mass-peak)*gamma))',mass,peak,gamma,beta,alpha)",
             "efficiency[0.9,0,1]",
             "signalFractionInPassing[0.9]"
        ),
  */
}
