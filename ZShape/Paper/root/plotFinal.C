#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH2.h"
#include "readUnfoldingMatrices.C"

static const int BINCOUNT=100;
static const double xlow=-5.0;
static const double xhigh=5.0;

struct DataSeries {

  DataSeries() {
    clear();
  }
  void clear() {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=xlow+(xhigh-xlow)/BINCOUNT*(i+0.5);
      xwidth[i]=(xhigh-xlow)/BINCOUNT/2;
      y[i]=0;
      ey[i]=0;
    }
  }
  void clearContent() {
    for (int i=0; i<BINCOUNT; i++) {
      y[i]=0;
      ey[i]=0;
    }
  }
  DataSeries(const DataSeries& ds) {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=xlow+(xhigh-xlow)/BINCOUNT*(i+0.5);
      xwidth[i]=(xhigh-xlow)/BINCOUNT/2;
      y[i]=ds.y[i];
      ey[i]=ds.ey[i];
    }
  }
  DataSeries(const char* fname, int iy=0);
  DataSeries(int lumi, const char* fnameformat, int iy=0);

  TGraphErrors* makeTGE(int bins, int startFrom) {
    return new TGraphErrors(bins,xave+startFrom,y+startFrom,xwidth+startFrom,ey+startFrom);
  }
  TGraph* makeTG(int bins, int startFrom) {
    return new TGraph(bins,xave+startFrom,y+startFrom);
  }

  TH1* makeTH1(const char* name, int bins, int startFrom) {
    printf("%s %d %f %f \n",name,bins,xave[startFrom]-xwidth[startFrom],xave[startFrom+bins-1]+xwidth[startFrom+bins-1]);
    TH1F* hist=new TH1F(name,name,bins,xave[startFrom]-xwidth[startFrom],xave[startFrom+bins-1]+xwidth[startFrom+bins-1]);
    for (int i=1; i<=bins; i++) {
      hist->SetBinContent(i,y[startFrom+i-1]);
      hist->SetBinError(i,ey[startFrom+i-1]);
      //      printf("%s %d %f %f \n",name,i+startFrom,y[startFrom+i-1],ey[startFrom+i-1]);
    }
    hist->SetDirectory(0);
    return hist;
  }

  void create(const char* file, int iy);

  double xave[BINCOUNT];
  double xwidth[BINCOUNT];
  double y[BINCOUNT];
  double ey[BINCOUNT];
};  


#include "tdrstyle.C"
#include "zrapidityStandard.C"

DataSeries::DataSeries(int lumi, const char* file, int iy) {
  char fullname[1024];
  sprintf(fullname,file,lumi);
  create(fullname,iy);
}
DataSeries::DataSeries(const char* file, int iy) {
  create(file,iy);
}
void DataSeries::create(const char* file, int iy) {
  clear();

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
    float a,b,q[20];
    int ate, atetotal=0;
    int found=sscanf(line," %d %f %f %n",&i,&a,&b,&ate);
    if (found!=3) continue;

    atetotal+=ate;
    for (int iq=0; iq<=iy && found!=0; iq++) {
      found=sscanf(line+atetotal," %f %n ",&q[iq],&ate);
      //      printf("%d %d %d %f\n",i,iq,ate,q[iq]);
      atetotal+=ate;
    }
    if (found!=1) continue;
    if (i<1 || i>BINCOUNT) continue;
    y[i-1]=q[iy];
  }
  fclose(f);
}

void plotFinal(TFile* mctruth, int mode=1) {
  int lumi=36;
  setTDRStyle();

  TH1* truth=(TH1*)(mctruth->Get("mcEff/All/Z0_Y_masscut")->Clone("truth"));
  truth->SetDirectory(0);
  
  DataSeries data_ee("ECAL80-ECAL95_RAW.csv");
  DataSeries data_ef("ECAL80-HF_RAW.csv"), data_all;
  DataSeries data_bkgd(data_all);
  
  DataSeries effAcc;

  const char* postfix="";
  const char* label="";

  if (mode==1 || mode==3) {
    effAcc=DataSeries("effAcc_combined.csv");
    postfix="_smeared";
    label="No Bin Migration Correction";
  }
  if (mode==2) {
    effAcc=DataSeries("effAcc_smearOverTrue.csv");
    postfix="_avemig";
    label="Average Bin Migration Correction";
  }
  if (mode==3) {
    effAcc=DataSeries("effAcc_combined.csv");
    postfix="_matrix";
    label="Matrix-based Unsmearing";
  }
  if (mode==4) {
    effAcc=DataSeries("effAcc_combined.csv");
    postfix="_foldMatrix";
    //    label="Folded Matrix-based Unsmearing";
  }
  DataSeries ea_statErr("effStatError.csv");

  const int firsti=16;
  const int lasti=85;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;

  /* 
     // study HF efficiency shifts
  for (int i=firsti-1; i<=lasti-1; i++)
    if (i<26 || i>74) effAcc.y[i]*=0.95;
  */

  for (int i=0; i<BINCOUNT; i++) {
    data_bkgd.y[i]=data_ef.y[i]+data_ee.y[i];
  }


  DataSeries dataStatError(data_all);

  DataSeries backgroundEE("background_all_y_ee.csv",0);
  DataSeries backgroundHF("background_all_y_ef.csv",0);

  DataSeries backgroundEEUnc("background_all_y_ee.csv",1);
  DataSeries backgroundHFUnc("background_all_y_ef.csv",1);

  DataSeries backgroundAll(backgroundEE), backgroundAllUnc(backgroundEEUnc),
    backgroundUncFrac(backgroundEEUnc);


  for (int i=0; i<BINCOUNT; i++) {  
    backgroundAll.y[i]=backgroundEE.y[i]+backgroundHF.y[i];
    backgroundAllUnc.y[i]=sqrt(pow(backgroundEEUnc.y[i],2)+pow(backgroundHFUnc.y[i],2));

    data_all.y[i]=data_bkgd.y[i]-backgroundAll.y[i];    
    data_all.ey[i]=data_bkgd.ey[i];
  }    

  DataSeries corrData(data_all);
  DataSeries totalSyst(data_all);


  for (int i=0; i<BINCOUNT; i++) {  
    if   (effAcc.y[i]==0) corrData.y[i]=0;
    else                  {
      corrData.y[i]/=effAcc.y[i];
      corrData.ey[i]=sqrt(data_bkgd.y[i])/effAcc.y[i];
    }
  }

  // unsmear
  DataSeries corrDataClone(corrData);


  bool       doMatrixUnsmearing(mode==3); // use this to activate matrix (type-3) unsmearing + associated errors or leave it off  
  if (doMatrixUnsmearing) {
    TH1* data_corr_smeared=corrData.makeTH1("data_corr_smeared",100,0);
    TH1* data_corr_unsmeared=unfold(data_corr_smeared,"../root/unfoldingMatrix_theOutPut.root");
    TFile*    theunfoldingMatrixInputFile = new TFile("../root/unfoldingMatrix_theOutPut.root","read");
    TMatrixD* theUnfoldingMatrix          = (TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal"); // indices [0.. 99] X [0.. 99]
    double errorCumul=0;
    for (int i=0; i<BINCOUNT; i++) {  
      corrData.y[i]=data_corr_unsmeared->GetBinContent(i+1);
      for (int s=0; s<BINCOUNT; s++) {
	errorCumul+= pow( corrDataClone.ey[s] * (*theUnfoldingMatrix)(s,i) , 2);
      }
      if (doMatrixUnsmearing) corrData.ey[i]=sqrt(errorCumul);
      errorCumul=0;
    }
  }
  DataSeries corrDataSyst(corrData);
  DataSeries corrDataBkgd(corrData);

  DataSeries pdfPos("pdfErrsPos_Y_TL.txt");
  DataSeries pdfNeg("pdfErrsNeg_Y_TL.txt");
  DataSeries pdfTotal(data_all);
  DataSeries pdfFrac(data_all);

  DataSeries energyScale("energyScaleError.csv");

  DataSeries unfoldSystMatrix("unfolding-relative-errors.txt");
  DataSeries unfoldSystAverage("averageUnfolding-relative-errors.txt");
  DataSeries unfoldSyst(unfoldSystAverage);
  // there's a switch doMatrixUnsmearing  to activate matrix (type 3) unsmearing + associated errors or leave it  off (e.g. when type-2 is done)
  if(doMatrixUnsmearing){
    for (int i=0; i<BINCOUNT; i++) {
      unfoldSyst.y[i]=unfoldSystMatrix.y[i];
    }
  }

  DataSeries effSystGSF("effSystGSFError.csv");
  DataSeries effSystWP80("effSystWP80Error.csv");
  DataSeries effSystWP95("effSystWP95Error.csv");
  DataSeries effSystHFEID("effSystHFEIDError.csv");
  DataSeries effSystHLT("effSystHLTError.csv");

  DataSeries effSystErr;

  for (int i=firsti-1; i<lasti; i++) {
    effSystErr.y[i]=sqrt(pow(effSystGSF.y[i],2)+
			 pow(effSystWP80.y[i],2)+
			 pow(effSystWP95.y[i],2)+
			 pow(effSystHFEID.y[i],2)+
			 pow(effSystHLT.y[i],2));
  }  

  double sumtotal=0;
  
  for (int i=firsti-1; i<lasti; i++) {
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      continue;
    }
    // what is this for??
    //corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars


    tt+=truth->GetBinContent(i+1);
    td+=corrData.y[i];

    backgroundUncFrac.y[i]=backgroundAllUnc.y[i]/std::max(1.0,data_all.y[i]);
    // GF
    //corrDataBkgd.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i];
    
    dataStatError.y[i]= corrDataBkgd.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    printf("%d %f\n",i+1,backgroundUncFrac.y[i]);

    //    corrDataSyst.y[i]/=effAcc.y[i];

    
    pdfTotal.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2*corrDataBkgd.y[i];
    pdfFrac.y[i]=(pdfPos.y[i]+(-pdfNeg.y[i]))/2;
    
    printf("%d %f %f %f %f\n",i+1,corrData.ey[i],corrData.y[i]*ea_statErr.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrDataBkgd.ey[i],2)+
			    pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			    pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			    pow(energyScale.y[i]*corrData.y[i],2)+
			    pow(unfoldSyst.y[i]*corrData.y[i],2)+
			    pow(effSystErr.y[i]*corrData.y[i],2)+
			    pow(pdfTotal.y[i],2)
			    );

    totalSyst.ey[i]=sqrt( pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			  pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			  pow(energyScale.y[i]*corrData.y[i],2)+
			  pow(unfoldSyst.y[i]*corrData.y[i],2)+
			  pow(effSystErr.y[i]*corrData.y[i],2)+
			  pow(pdfTotal.y[i],2)
			  );

    totalSyst.y[i]=totalSyst.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    sumtotal+=corrDataSyst.y[i];
  }
  char fnamework[1024];

  
  DataSeries corrDataFold(corrData),corrDataSystFold(corrData);

  sprintf(fnamework,"ZRapidity_final_fold_table%s.tex",postfix);

  FILE* ftable=fopen(fnamework,"w");

  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n");
  fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Measurement & Statistical Error & Systematic Error \\\\ \\hline\n");
  
  const int izero=51;
  const double binwidth=0.1;
  //  double binwidth=1.0;
  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;
    corrDataFold.y[j]=(corrDataBkgd.y[j]+corrDataBkgd.y[jadd])/std::max(1.0,sumtotal);


    //    printf("%d %d (%f) %d (%f) \n",i,j,corrDataBkgd.y[j],jadd,corrDataBkgd.y[jadd]);

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    double pdferr=(pdfFrac.y[j]+pdfFrac.y[jadd])/2*corrDataFold.y[j];
    double bkgderr=sqrt(pow(backgroundAllUnc.y[j],2)+pow(backgroundAllUnc.y[jadd],2))/ea_Ave/sumtotal;
    double energyerr=(energyScale.y[j]+energyScale.y[jadd])/2*corrDataFold.y[j];
    double unfolder=(unfoldSyst.y[j]+unfoldSyst.y[jadd])/2*corrDataFold.y[j];

    //corrDataFold.ey[j]=sqrt(data_all.y[j]+data_all.y[jadd]+backgroundAll.y[j]+backgroundAll.y[jadd])/ea_Ave/sumtotal;//FG
    corrDataFold.ey[j]=sqrt(data_bkgd.y[j]+data_bkgd.y[jadd])/ea_Ave/sumtotal;
    //corrData.ey[j]*effAcc.y[j],2) + pow(corrData.ey[jadd]*effAcc.y[j],2) )/ea_Ave/sumtotal;

    corrDataSystFold.y[j]=corrDataFold.y[j];
    corrDataSystFold.ey[j]=sqrt(pow(corrDataFold.ey[j],2)+
				pow(corrDataFold.y[j]*(ea_statErr.y[j]+ea_statErr.y[jadd])/2,2)+
				pow(corrDataFold.y[j]*(effSystErr.y[j]+effSystErr.y[jadd])/2,2)+
				pow(pdferr,2)+
				pow(energyerr,2)+
				pow(unfolder,2)+
				pow(bkgderr,2)
				);

    corrDataFold.y[j]/=binwidth;
    corrDataFold.ey[j]/=binwidth;
    corrDataSystFold.y[j]/=binwidth;
    corrDataSystFold.ey[j]/=binwidth;


    fprintf(ftable,"%7.2f & %7.2f & %7.4f & %7.4f & %7.4f \\\\ \n",
	    corrDataFold.xave[j]-corrDataFold.xwidth[j],
	    corrDataFold.xave[j]+corrDataFold.xwidth[j],
	    corrDataFold.y[j],
	    corrDataFold.ey[j],
	    sqrt(pow(corrDataSystFold.ey[j],2)-pow(corrDataFold.ey[j],2))
	    );
  }

  fprintf(ftable,"\\hline\n\\end{tabular}");
  fclose(ftable);
 
  truth->Scale(td/tt*1.0);
  
  TGraph* rawd=data_all.makeTG(npoint,firsti-1);
  TGraph* rawdb=data_bkgd.makeTG(npoint,firsti-1);
  TGraph* corrd=corrData.makeTGE(npoint,firsti-1);
  TGraph* corrdsys=corrDataSyst.makeTGE(npoint,firsti-1);
  TGraph* corrdbkgd=corrDataBkgd.makeTGE(npoint,firsti-1);
  // this is special
  TGraph* corrdsysb=new TGraphErrors(npoint,corrDataBkgd.xave+firsti-1,corrDataBkgd.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  TH2* dummy=new TH2F("dummy","dummy",20,-3.5,3.5,30,0,700);

  dummy->GetYaxis()->SetTitle("Events/0.1 Units of Rapidity");
  dummy->GetXaxis()->SetTitle(y_xaxis_label);
  dummy->GetXaxis()->CenterTitle();
  dummy->GetYaxis()->CenterTitle();
  
  TCanvas* c1=new TCanvas("finalZRap","finalZRap",800,600);
  dummy->Draw();

  dummy->SetDirectory(0);

  rawd->SetMarkerStyle(24);
  rawd->Draw("PSAME");
  rawdb->SetMarkerStyle(25);
  rawdb->Draw("PSAME");

  //  corrd->SetMarkerStyle(20);
  //  corrd->Draw("PSAME");

  // corrdsys->SetMarkerStyle(0);
  // corrdsys->Draw("PSAME");

  corrdbkgd->SetMarkerStyle(20);
  corrdbkgd->Draw("PSAME");

  corrdsysb->SetMarkerStyle(0);
  corrdsysb->Draw("PSAME");

  truth->Draw("SAMEHIST");

  TLegend* tl=new TLegend(0.4,0.16,0.72,0.34,label);
  tl->SetFillStyle(0);
  tl->AddEntry(corrd,"Corrected Data","P");
  tl->AddEntry(rawdb,"Raw Data","P");
  tl->AddEntry(rawd,"Background-Subtracted Data","P");
  tl->AddEntry(truth,"POWHEG+CT10 Prediction","l");
  tl->Draw();

  zrap_Prelim(0.85,0.975,0.82,0.80);
  zrap_Lumi(0.81,0.89,lumi);

  sprintf(fnamework,"ZRapidity_final%s.eps",postfix);
  c1->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final%s.png",postfix);
  c1->Print(fnamework);


  TCanvas* c2=new TCanvas("finalZRapFold","finalZRapFold",800,600);
  TH2* dummy2=new TH2F("dummy2","dummy2",10,0.0,3.5,30,0,0.4);
  dummy2->SetDirectory(0);

  TH1* truth2=(TH1*)(truth->Clone("truth2"));
  TH1* truthalt=(TH1*)(truth->Clone("truthalt"));
  truth2->SetDirectory(0);
  truthalt->SetDirectory(0);


  // this is where you unfold |Y|, if mode==4
  DataSeries corrDataFoldClone(corrDataFold);
  DataSeries corrDataSystFoldClone(corrDataSystFold);
    // use this to activate matrix (type-3) unsmearing for the folded rapidity distribution |Y|
  if (mode==4) {
    TH1*      data_corr_fold_smeared      =corrDataFold.makeTH1("data_corr_fold_smeared",50,50);
    TH1*      data_corr_unsmeared         =unfold(data_corr_fold_smeared,"../root/unfoldingMatrix_theOutPut.root",true);              // unfolding for |Y|
    TFile*    theunfoldingMatrixInputFile =new TFile("../root/unfoldingMatrix_theOutPut.root","read");
    TMatrixD* theUnfoldingMatrix          =(TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotalFolded"); // indices [0.. 49] X [0.. 49]
    double errorCumul=0;
    double errorCumulSyst=0;
    for (int i=0; i<=(lasti-firsti)/2; i++) { 
      int j=izero+i-1;
      corrDataFold.y[j]=data_corr_fold_smeared->GetBinContent(i+1);
      // now do error propagation
      for (int s=0; s<=(lasti-firsti)/2; s++) { 
	int js=izero+s-1;
	errorCumul+=     pow( corrDataFoldClone.ey[js]     * (*theUnfoldingMatrix)(s,i) , 2);
	errorCumulSyst+= pow( corrDataSystFoldClone.ey[js] * (*theUnfoldingMatrix)(s,i) , 2);
	//std::cout << "js: " << js << "\t" << corrDataFold.y[js] << "\t"<< corrDataFoldClone.ey[js] 
	//	  << "\t" << corrDataSystFoldClone.ey[js] << "\t"<< (corrDataFoldClone.ey[js]/corrDataSystFoldClone.ey[js]) << std::endl;
      }
      //std::cout << "errorCumul: "<< errorCumul << "\t" << errorCumulSyst << std::endl;
      corrDataFold.ey[j]    =sqrt(errorCumul*1);
      corrDataSystFold.ey[j]=sqrt(errorCumulSyst*1);
      errorCumul=0;
      errorCumulSyst=0;
    }
  }
  

  double chi2=0;
  int ndof=0;
  for (int i=0; i<=(lasti-firsti)/2+2; i++) {

    int j=izero+i-1;
    int jadd=izero-i-2;

    truth2->SetBinContent(j+1,(truth2->GetBinContent(j+1)+truth2->GetBinContent(jadd+1))/sumtotal/binwidth);

    double toadd=pow(truth2->GetBinContent(j+1)-corrDataFold.y[j],2)/pow(corrDataSystFold.ey[j],2);
    printf("%d %f\n",i,toadd);
    if (corrDataFold.y[j]>0) {
      ndof++;
      chi2+=toadd;
    }
  }
  ndof--; // because of normalization


  c2->cd();
  dummy2->GetYaxis()->SetTitle("1/#sigma d#sigma/d|y|");
  dummy2->GetXaxis()->SetTitle("|y_{ee}|");
  dummy2->GetXaxis()->CenterTitle();
  dummy2->GetYaxis()->CenterTitle();

  dummy2->Draw();

  TGraph* corrdfold=new TGraphErrors(npoint/2,corrDataFold.xave+izero-1,corrDataFold.y+izero-1,corrDataFold.xwidth+izero-1,corrDataFold.ey+izero-1);
  corrdfold->Draw("PSAME");
  TGraph* corrdfoldsys=new TGraphErrors(npoint/2,corrDataFold.xave+izero-1,corrDataFold.y+izero-1,corrDataSystFold.xwidth+izero-1,corrDataSystFold.ey+izero-1);
  corrdfold->Draw("PSAME");
  corrdfoldsys->SetMarkerStyle(0);
  corrdfoldsys->Draw("PSAME"); // is this stat-only error?

  truth2->Draw("SAMEHIST");
  truthalt->SetLineColor(kRed);
  truthalt->Draw("SAMEHIST");

  zrap_Prelim(0.85,0.975,0.32,0.15);
  zrap_Lumi(0.82,0.89,lumi);

  if (strlen(label)>2) {
      TLatex *plabel = new TLatex(0.7,0.17,label);
      plabel-> SetNDC();
      plabel -> SetTextFont(42);
      plabel -> SetTextColor(1);
      plabel -> SetTextSize(0.040);
      plabel -> SetTextAlign(22);
      plabel -> SetTextAngle(0);
      plabel -> Draw();
  }

  sprintf(fnamework,"ZRapidity_final_fold%s.eps",postfix);
  c2->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final_fold%s.png",postfix);
  c2->Print(fnamework);

  TCanvas* c3=new TCanvas("finalZRapErrors","finalZRapErrors",800,600);
  TH1* dummy3=new TH1F("dummy3","dummy3",10,-3.5,3.5);
  dummy3->SetDirectory(0);
  c3->SetLogy();
  dummy3->SetMaximum(1.0);
  dummy3->SetMinimum(0.001);
  dummy3->GetYaxis()->SetTitle("Fractional Error");
  dummy3->GetXaxis()->SetTitle(y_xaxis_label);
  dummy3->GetXaxis()->CenterTitle();

  
  //  TH1* binmiggr=binMigPos.makeTH1("bin migration",npoint,firsti-1);
  TH1* effacc_stat_gr=ea_statErr.makeTH1("eff stat",npoint,firsti-1);
  TH1* effacc_syst_gr=effSystErr.makeTH1("eff syst",npoint,firsti-1);
  TH1* gr_dataStatError=dataStatError.makeTH1("data stats err",npoint,firsti-1);
  TH1* gr_bkgd_unc=backgroundUncFrac.makeTH1("bkgd unc",npoint,firsti-1);
  TH1* gr_pdf=pdfFrac.makeTH1("pdf_err",npoint,firsti-1);
  TH1* gr_energyscale=energyScale.makeTH1("es_unc",npoint,firsti-1);
  TH1* gr_unfold=unfoldSyst.makeTH1("unfold unc",npoint,firsti-1);
  TH1* gr_systotal=totalSyst.makeTH1("total_syst",npoint,firsti-1);
  dummy3->Draw();

  //  binmiggr->SetLineWidth(2);
  //  binmiggr->SetLineColor(kBlue);

  effacc_stat_gr->SetLineWidth(2);
  effacc_stat_gr->SetLineColor(kGreen+2);

  effacc_syst_gr->SetLineWidth(2);
  effacc_syst_gr->SetLineColor(kOrange+8);

  gr_dataStatError->SetLineWidth(2);
  gr_dataStatError->SetLineColor(kBlack);
  gr_dataStatError->SetLineStyle(2);

  gr_bkgd_unc->SetLineWidth(2);
  gr_bkgd_unc->SetLineColor(kRed-4);

  gr_energyscale->SetLineWidth(2);
  gr_energyscale->SetLineColor(kBlue);


  gr_unfold->SetLineWidth(2);
  gr_unfold->SetLineStyle(3);
  gr_unfold->SetLineColor(kViolet);
  
  gr_pdf->SetLineWidth(2);
  gr_pdf->SetLineColor(kViolet-5);

  gr_systotal->SetLineWidth(2);
  

  effacc_stat_gr->Draw("HSAME");
  effacc_syst_gr->Draw("HSAME");
  gr_dataStatError->Draw("HSAME");
  gr_bkgd_unc->Draw("HSAME");
  gr_energyscale->Draw("HSAME");
  if (mode==2) 
    gr_unfold->Draw("HSAME");

  gr_pdf->Draw("HSAME");
  gr_systotal->Draw("HISTSAME");


  tl=new TLegend(0.18,0.73,0.48,0.92);
  tl->SetFillStyle(0);
  TLegend* tl2=new TLegend(0.47,0.73,0.73,0.92,label);
  tl2->SetFillStyle(0);
   //  tl->AddEntry(binmiggr,"Bin Migration","L");
  tl->AddEntry(gr_dataStatError,"Data Statistics","L");
  tl->AddEntry(gr_systotal,"Total Systematic Error","L");
  tl->AddEntry(gr_bkgd_unc,"Background Uncertainty","L");
  tl->AddEntry(effacc_stat_gr,"Efficiency Statistics","L");
  tl2->AddEntry(gr_energyscale,"Energy Scale Uncertainty","L");
  tl2->AddEntry(gr_pdf,"PDF (#epsilon #times A) Uncertainty","L");
  if (mode==2)
    tl2->AddEntry(gr_unfold,"Unfolding Uncertainty","L");
  tl2->AddEntry(effacc_syst_gr,"Efficiency Systematics","L");
  tl->Draw();
  tl2->Draw();
  
  zrap_Prelim(0.85,0.975,0.3,0.16);
  zrap_Lumi(0.84,0.90,lumi);

  sprintf(fnamework,"ZRapidity_final_errors%s.eps",postfix);
  c3->Print(fnamework);
  sprintf(fnamework,"ZRapidity_final_errors%s.png",postfix);
  c3->Print(fnamework);

  sprintf(fnamework,"ZRapidity_final_fold_errtable%s.tex",postfix);
  ftable=fopen(fnamework,"w");

  sprintf(fnamework,"ZRapidity_fold_zee%s.txt",postfix);
  FILE* ftxt=fopen(fnamework,"w");
  
  for (int i=0; i<=(lasti-firsti)/2; i++) {
    int j=izero+i-1;
    int jadd=izero-i-2;

    double ea_Ave=(effAcc.y[j]+effAcc.y[jadd])/2;

    //    double errmig=(fabs(binMig.y[j])*corrData.y[j]+fabs(binMig.y[jadd])*corrData.y[jadd])/sumtotal;
    double pdferr=(pdfFrac.y[j]+pdfFrac.y[jadd])/2;//*corrDataFold.y[j]*binwidth/sumtotal;
    double bkgderr=sqrt(pow(backgroundAllUnc.y[j],2)+pow(backgroundAllUnc.y[jadd],2))/ea_Ave/sumtotal;
    double staterr=(ea_statErr.y[j]+ea_statErr.y[jadd])/2;
    double systerr=(effSystErr.y[j]+effSystErr.y[jadd])/2;
    double energyerr=(energyScale.y[j]+energyScale.y[jadd])/2;
    double unfolderr=(unfoldSyst.y[j]+unfoldSyst.y[jadd])/2;

    //    errmig/=binwidth*corrDataFold.y[j];
    //    pdferr/=binwidth*corrDataFold.y[j];
    bkgderr/=binwidth*corrDataFold.y[j];
    //    staterr/=binwidth*corrDataFold.y[j];
    //    systerr/=binwidth*corrDataFold.y[j];


    if (i==0) {
      fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|c|c|}\\hline\n");
      fprintf(ftable,"            &             &  Background & Efficiency & Energy Scale & PDF \\effacc & Unfolding  \\\\ \n");
      fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Estimation & Determination & Systematics & Error & Uncertainty \\\\ \\hline \n");
      fprintf(ftxt,"# Zee rapidity results\n");
      fprintf(ftxt,"# bin  y_min  y_max value staterr systerr allerr \n");
    }
      fprintf(ftable," %7.2f & %7.2f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n",
             corrDataFold.xave[j]-corrDataFold.xwidth[j],
             corrDataFold.xave[j]+corrDataFold.xwidth[j],
	      bkgderr,sqrt(staterr*staterr+systerr*systerr),
	      energyerr,pdferr,unfolderr);
      fprintf(ftxt," %3d  %7.1f %7.1f %8.4f %8.4f %8.4f %8.4f\n",j+1,
	      corrDataFold.xave[j]-corrDataFold.xwidth[j],
	      corrDataFold.xave[j]+corrDataFold.xwidth[j],
	      corrDataFold.y[j],corrDataFold.ey[j],
	      sqrt(pow(corrDataSystFold.ey[j],2)-pow(corrDataFold.ey[j],2)),
	      corrDataSystFold.ey[j]);
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);
  fclose(ftxt);

  printf(" Consistency: %f/%d %f \n",chi2,ndof,TMath::Prob(chi2,ndof));

}
