#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH2.h"
#include "readUnfoldingMatrices.C"
#include "TMath.h"

static const int    pt_bins=19;
static const double pt_binning[]={ 0.0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
static const double pt_binning_vis[]={ 0.7, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };

static const int BINCOUNT=pt_bins;


struct DataSeries {

  DataSeries() {
    clear();
  }
  void clear() {
    for (int i=0; i<BINCOUNT; i++) {
      xave[i]=(pt_binning[i]+pt_binning[i+1])/2;
      xwidth[i]=(pt_binning[i+1]-pt_binning[i])/2;
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
      xave[i]=(pt_binning[i]+pt_binning[i+1])/2;
      xwidth[i]=(pt_binning[i+1]-pt_binning[i])/2;
      y[i]=ds.y[i];
      ey[i]=ds.ey[i];
    }
  }
  DataSeries(const char* fname, int iy=0);
  DataSeries(int lumi, const char* fnameformat, int iy=0);

  TGraphErrors* makeTGE() {
    return new TGraphErrors(BINCOUNT,xave,y,xwidth,ey);
  }
  TGraph* makeTG() {
    return new TGraph(BINCOUNT,xave,y);
  }

  TH1* makeTH1(const char* name) {
    TH1F* hist=new TH1F(name,name,BINCOUNT,pt_binning_vis);
    for (int i=1; i<=BINCOUNT; i++) {
      hist->SetBinContent(i,y[i-1]);
      hist->SetBinError(i,ey[i-1]);
      //      printf("%s %d %f %f \n",name,i+startFrom,y[startFrom+i-1],ey[startFrom+i-1]);
    }
    hist->SetDirectory(0);
    return hist;
  }

  void create(const char* file, int iy);

  void normalizeToWidth(double td=1.0) {
    for (int i=0; i<BINCOUNT; i++) {
      y[i]/=(xwidth[i]*2)*td;
      ey[i]/=(xwidth[i]*2)*td;
    }
  }

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

void plotFinalQt(TFile* mctruth, int mode) {
  int lumi=36;
  setTDRStyle();

  //  TH2* truth2d=(TH1*)(mctruth->Get("mcEff/ECAL80-ECAL95-MUO/Acceptance/PtTL_vs_PtZ"));
  //  TH1* proj=truth2d->ProjectionY();
  TH1* truth=0;
  
  DataSeries data_ee("ECAL80-ECAL95_Pt_RAW.csv");
  DataSeries data_all;
  DataSeries data_bkgd(data_all);
  
  DataSeries effAcc;
  DataSeries ea_statErr("effStatPtError.csv");

  const char* postfix="";

  if (mode==1) {
    effAcc=DataSeries("effAcc_pt_massCutdenom.txt");
    postfix="_smear";
    //    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON/Acceptance/Z0_Pt_masscut")->Clone("truth"));
    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON-G20/C01-GPT20/Z0_Pt_masscut")->Clone("truth"));
  } else if (mode==2) {
    effAcc=DataSeries("effAcc_pt_TLdenom.txt");
    postfix="_avemig";
    //    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON/Acceptance/Z0_PtTL_masscut")->Clone("truth"));
    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON-G20/C01-GPT20/Z0_PtTL_masscut")->Clone("truth"));
  } else if (mode==3) {
    effAcc=DataSeries("effAcc_pt_massCutdenom.txt");
    postfix="_matrix";
    //    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON/Acceptance/Z0_PtTL_masscut")->Clone("truth"));
    truth=(TH1*)(mctruth->Get("mcEff/MUON-MUON-G20/C01-GPT20/Z0_PtTL_masscut")->Clone("truth"));
  } else {
    std::cout << "the value you've input: " << mode << " is not valid" << std::endl;
    return;
  }

  truth->SetDirectory(0);

  const int firsti=1;
  const int lasti=18;
  const int npoint=lasti-firsti+1;

  double tt=0, td=0;
 
  for (int i=0; i<BINCOUNT; i++) {
    data_bkgd.y[i]=data_ee.y[i];
  }


  DataSeries dataStatError(data_all);

  DataSeries backgroundAll("background_all_qt.csv",0);
  DataSeries backgroundAllUnc("background_all_qt.csv",1);

  DataSeries backgroundUncFrac(backgroundAllUnc);

  for (int i=0; i<BINCOUNT; i++) {  
    data_all.y[i]=data_bkgd.y[i]-backgroundAll.y[i];    
    data_all.ey[i]=data_bkgd.ey[i];
  }    

  DataSeries corrData(data_all);
  DataSeries totalSyst(data_all);

  for (int i=0; i<BINCOUNT; i++) {  
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      continue;
    }
    //    printf("%2d %f %f \n",i,corrData.y[i],effAcc.y[i]);
    corrData.y[i]/=effAcc.y[i];
  }
  // unsmear


  for (int i=0; i<BINCOUNT; i++) {  
    if   (effAcc.y[i]==0) corrData.y[i]=0;
    else                  corrData.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i]; // this is statistical err, post effxacc corrections
  }
  DataSeries corrDataClone(corrData);

  bool doUnsmearing(mode==3); // use this to activate unsmearing + associated errors or leave it completely off
  if (doUnsmearing) {
    TH1* data_corr_smeared=corrData.makeTH1("data_corr_smeared");
    TH1* data_corr_unsmeared=unfold(data_corr_smeared,"../root/unfoldingMatrix-pt_theOutPut.root");
    TFile*    theunfoldingMatrixInputFile = new TFile("../root/unfoldingMatrix-pt_theOutPut.root","read");
    TMatrixD* theUnfoldingMatrix          = (TMatrixD*)theunfoldingMatrixInputFile->Get("unsmearMatrices/unfoldingMatrixTotal");
    double errorCumul=0;
    for (int i=0; i<BINCOUNT; i++) {  
      corrData.y[i]=data_corr_unsmeared->GetBinContent(i+1);
      // unfolded bin = linear combination ( folded bins, theUnfoldingMatrix-column)
      // propagate errors accordingly, assuming uncorrelated variables; errors on corrData are statistical only 
      for (int s=0; s<BINCOUNT; s++) {
	errorCumul+= pow( corrDataClone.ey[s] * (*theUnfoldingMatrix)(s,i) , 2);
      }
      if (doUnsmearing) corrData.ey[i]=sqrt(errorCumul);
      errorCumul=0;
    }
  }
  DataSeries corrDataSyst(corrData);
  DataSeries corrDataBkgd(corrData);

  DataSeries pdfPos("pdfErrsPos_QT.txt");
  DataSeries pdfNeg("pdfErrsNeg_QT.txt");
  DataSeries pdfTotal(data_all);
  DataSeries pdfFrac(data_all);

  DataSeries energyScale("energyScaleErrorPt.csv");

  DataSeries effSystGSF("effSystGSFPtError.csv");
  DataSeries effSystWP80("effSystWP80PtError.csv");
  DataSeries effSystWP95("effSystWP95PtError.csv");
  DataSeries effSystHLT("effSystHLTPtError.csv");
  
  // pick unfolding systematic consistent with unfolding performed
  DataSeries unfoldSyst;
  if(mode==2) {
    DataSeries unfoldSystTmp("averageUnfolding-relative-errors-pt.txt");
    for (int i=firsti-1; i<lasti; i++) {    unfoldSyst.y[i]=unfoldSystTmp.y[i];   }
  }
  if(mode==3){ 
    DataSeries unfoldSystTmp("unfolding-relative-errors-pt.txt");
    for (int i=firsti-1; i<lasti; i++) {  unfoldSyst.y[i]=unfoldSystTmp.y[i];  } 
  }
  
  DataSeries effSystErr;

  for (int i=firsti-1; i<lasti; i++) {
    effSystErr.y[i]=sqrt(pow(effSystGSF.y[i],2)+
			 pow(effSystWP80.y[i],2)+
			 pow(effSystWP95.y[i],2)+
			 pow(effSystHLT.y[i],2)); //GF add unfolding here
  }  

  double sumtotal=0;
  
  for (int i=firsti-1; i<lasti; i++) {
    if (effAcc.y[i]==0) {
      corrData.y[i]=0;
      printf("hih\n");
      continue;
    }
    // what is this for??
    //corrData.ey[i]=sqrt(data_all.y[i])/effAcc.y[i];
    //    corrDataSyst.xave[i]+=0.02; // uncomment to offset error bars

    backgroundUncFrac.y[i]=backgroundAllUnc.y[i]/std::max(1.0,data_bkgd.y[i]);
    // GF
    //corrDataBkgd.ey[i]=sqrt(data_all.y[i]+backgroundAll.y[i])/effAcc.y[i];
    
    dataStatError.y[i]= corrDataBkgd.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    //    corrDataSyst.y[i]/=effAcc.y[i];

    
    pdfTotal.y[i]=(pdfPos.y[i]+fabs(pdfNeg.y[i]))/2*corrDataBkgd.y[i];
    pdfFrac.y[i]=(pdfPos.y[i]+fabs(pdfNeg.y[i]))/2;
    
    //  printf("%d %f %f %f %f\n",i,corrData.ey[i],corrData.y[i]*ea_statErr.y[i]);
    corrDataSyst.ey[i]=sqrt(pow(corrDataBkgd.ey[i],2)+
			    pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			    pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			    pow(energyScale.y[i]*corrData.y[i],2)+
			    pow(unfoldSyst.y[i]*corrData.y[i],2)+
			    pow(effSystErr.y[i]*corrData.y[i],2)+
			    0
			    //			    pow(pdfTotal.y[i],2)
			    );
    
    totalSyst.ey[i]=sqrt( pow(corrDataBkgd.y[i]*ea_statErr.y[i],2)+
			  pow(backgroundAllUnc.y[i]/effAcc.y[i],2)+
			  pow(energyScale.y[i]*corrData.y[i],2)+
			  pow(unfoldSyst.y[i]*corrData.y[i],2)+
			  pow(effSystErr.y[i]*corrData.y[i],2)+
			  0
			  //			  pow(pdfTotal.y[i],2)
			  );

    totalSyst.y[i]=totalSyst.ey[i]/std::max(1.0,corrDataBkgd.y[i]);

    sumtotal+=corrDataSyst.y[i];
  }
  char fnamework[1024];

  for (int i=0; i<BINCOUNT; i++) {
    td+=corrData.y[i];
  }  
 
  data_all.normalizeToWidth(td);
  data_bkgd.normalizeToWidth(td);
  corrData.normalizeToWidth(td);
  corrDataSyst.normalizeToWidth(td);
  corrDataBkgd.normalizeToWidth(td);


  TH1* truth_vis=new TH1F("tvis","tvis",BINCOUNT,pt_binning_vis);
  double td2=0;

  for (int i=0; i<BINCOUNT; i++) {
    truth_vis->SetBinContent(i+1,truth->GetBinContent(i+1)/(2*corrData.xwidth[i]));

    tt+=truth_vis->GetBinContent(i+1);
    td2+=corrData.y[i];
  }

  truth_vis->Scale(td2/tt*1.0);

  sprintf(fnamework,"ZQt_final_table%s.tex",postfix);

  FILE* ftable=fopen(fnamework,"w");

  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|}\n\\hline\n");
  fprintf(ftable,"$|q_{T,min}|$ & $|q_{T,max}|$ & Measurement & Statistical Error & Systematic Error \\\\ \\hline\n");
  for (int j=firsti-1; j<lasti; j++) {
    fprintf(ftable,"%7.2f & %7.2f & ",
	    corrData.xave[j]-corrData.xwidth[j],
	    corrData.xave[j]+corrData.xwidth[j]);
    if (corrData.y[j]>0.01)
      fprintf(ftable,"%7.4f & %7.4f & %7.4f \\\\ \n",
	      corrData.y[j],
	      corrData.ey[j],                                       // stat_error (post-unfolding) since corrData includes only statistical
	      sqrt(pow(corrDataSyst.ey[j],2)-pow(corrData.ey[j],2)) // syst_error := (tot_error quadMinus stat_err) 
	      );
    else     if (corrData.y[j]>0.001)

      fprintf(ftable,"%7.5f & %7.5f & %7.5f \\\\ \n",
	      corrData.y[j],
	      corrData.ey[j],
	      sqrt(pow(corrDataSyst.ey[j],2)-pow(corrData.ey[j],2))
	      );
    else
      fprintf(ftable,"%7.6f & %7.6f & %7.6f \\\\ \n",
	      corrData.y[j],
	      corrData.ey[j],
	      sqrt(pow(corrDataSyst.ey[j],2)-pow(corrData.ey[j],2))
	      );
  }
  
  fprintf(ftable,"\\hline\n\\end{tabular}");
  fclose(ftable);
 



  
  TGraph* rawd=data_all.makeTG();
  TGraph* rawdb=data_bkgd.makeTG();
  TGraph* corrd=corrData.makeTGE();
  //  TGraph* corrdsys=corrDataSyst.makeTGE();
  TGraph* corrdbkgd=corrDataBkgd.makeTGE();
  // this is special
  TGraph* corrdsysb=new TGraphErrors(npoint,corrDataBkgd.xave+firsti-1,corrDataBkgd.y+firsti-1,corrDataSyst.xwidth+firsti-1,corrDataSyst.ey+firsti-1);

  double ybinage[21];
  for (int i=0; i<21; i++) {
    ybinage[i]=1e-6+0.07/20.0*i;
  }

  TH2* dummy=new TH2F("dummy","",18,pt_binning_vis,20,ybinage);

  dummy->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T} [GeV^{-1}]");
  dummy->GetXaxis()->SetTitle(qt_xaxis_label);
  dummy->GetXaxis()->SetTitleOffset(1.2);
  dummy->GetXaxis()->CenterTitle();
  
  TCanvas* c1=new TCanvas("finalZQt","finalZQt",800,600);
  c1->SetLogx();
  c1->SetLogy(false);
  c1->SetBottomMargin(0.16);
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

  truth_vis->Draw("SAMEHIST ][");

  TLegend* tl=new TLegend(0.55,0.45,0.94,0.70);
  tl->SetFillStyle(0);
  tl->AddEntry(corrd,"Corrected Data","P");
  tl->AddEntry(rawdb,"Raw Data","P");
  tl->AddEntry(rawd,"Background-Subtracted Data","P");
  tl->AddEntry(truth_vis,"POWHEG+CT10 Prediction","l");
  tl->Draw();


  zrap_Prelim(0.22,0.97,0.82,0.17);
  zrap_Lumi(0.85,0.97,lumi);

  sprintf(fnamework,"ZQt_final%s.eps",postfix);
  c1->Print(fnamework);
  sprintf(fnamework,"ZQt_final%s.C",postfix);
  c1->Print(fnamework);

  TCanvas* c1l=new TCanvas("finalZQtLog","finalZQtLog",800,600);
  c1l->cd();
  c1l->SetLogx();
  c1l->SetLogy();
  c1l->SetBottomMargin(0.16);

  TH1* dummyl=(TH2*)(dummy->Clone("logdummy"));
  dummyl->Draw();
  dummyl->SetDirectory(0);

  TH1* rawdl=(TH1*)(rawd->Clone("raw log"));
  rawdl->Draw("PSAME");
  TH1* rawdbl=(TH1*)(rawdb->Clone("rawsub log"));
  rawdbl->Draw("PSAME");

  //  corrd->SetMarkerStyle(20);
  //  corrd->Draw("PSAME");

  // corrdsys->SetMarkerStyle(0);
  // corrdsys->Draw("PSAME");
  
  TH1* corrdbkgdl=(TH1*)(corrdbkgd->Clone("final log"));
  corrdbkgdl->Draw("PSAME");
  
  TH1* corrdsysbl=(TH1*)(corrdsysb->Clone("finalsys log"));
  corrdsysbl->Draw("PSAME");
  
  TH1* truth_visl=(TH1*)(truth_vis->Clone("truth log"));
  truth_visl->Draw("SAMEHIST ][");
  
  
  TLegend* tl2=new TLegend(0.2,0.20,0.62,0.38);
  tl2->SetFillStyle(0);

  tl2->AddEntry(corrdbkgdl,"Corrected Data","P");
  tl2->AddEntry(rawdbl,"Raw Data","P");
  tl2->AddEntry(rawdl,"Background-Subtracted Data","P");
    tl2->AddEntry(truth_visl,"POWHEG+CT10 Prediction","l");
  tl2->Draw();
  
  zrap_Prelim(0.22,0.97,0.82,0.17);
  zrap_Lumi(0.85,0.97,lumi);

  sprintf(fnamework,"ZQt_final-log%s.eps",postfix);
  c1l->Print(fnamework);
  sprintf(fnamework,"ZQt_final-log%s.C",postfix);
  c1l->Print(fnamework);


  TCanvas* c3=new TCanvas("finalZQtErrors","finalZQtErrors",800,600);
  TH1* dummy3=new TH1F("dummy3","",10,0.7,600);
  c3->SetBottomMargin(0.15);
  dummy3->SetDirectory(0);
  c3->SetLogy();
  c3->SetLogx();
  dummy3->SetMaximum(1.0);
  dummy3->SetMinimum(0.001);
  dummy3->GetYaxis()->SetTitle("Fractional Error");
  dummy3->GetXaxis()->SetTitle(qt_xaxis_label);
  dummy3->GetXaxis()->CenterTitle();
  dummy3->GetXaxis()->SetTitleOffset(1.05);

  
  TH1* effacc_stat_gr=ea_statErr.makeTH1("eff stat");
  TH1* effacc_syst_gr=effSystErr.makeTH1("eff syst");
  TH1* gr_dataStatError=dataStatError.makeTH1("data stats err");
  TH1* gr_bkgd_unc=backgroundUncFrac.makeTH1("bkgd unc");
  TH1* gr_pdf=pdfFrac.makeTH1("pdf_err");
  TH1* gr_energyscale=energyScale.makeTH1("es_unc");
  TH1* gr_unfolding=unfoldSyst.makeTH1("unfold unc");
  TH1* gr_systotal=totalSyst.makeTH1("total_syst");
  dummy3->Draw();

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

  gr_unfolding->SetLineWidth(2);
  gr_unfolding->SetLineColor(kViolet);

  gr_pdf->SetLineWidth(3);
  gr_pdf->SetLineStyle(3);
  gr_pdf->SetLineColor(kViolet-5);

  gr_systotal->SetLineWidth(2);
  

  effacc_stat_gr->Draw("HSAME ][");
  effacc_syst_gr->Draw("HSAME ][");
  gr_dataStatError->Draw("HSAME ][");
  gr_bkgd_unc->Draw("HSAME ][");
  gr_energyscale->Draw("HSAME ][");
  gr_pdf->Draw("HSAME ][");
  gr_unfolding->Draw("HSAME ][");
  gr_systotal->Draw("HISTSAME ][");

  tl=new TLegend(0.15,0.75,0.75,0.94,"");
  tl->SetFillStyle(0);
  tl->SetNColumns(2);
  tl->AddEntry(gr_dataStatError,"Data Statistics","L");
  tl->AddEntry(gr_systotal,"Total Systematic Error","L");
  tl->AddEntry(gr_bkgd_unc,"Background Uncertainty","L");
  tl->AddEntry(effacc_stat_gr,"Efficiency Statistics","L");
  tl->AddEntry(gr_pdf,"PDF (#epsilon #times A) Uncertainty","L");
  tl->AddEntry(gr_unfolding,"Unfolding Uncertainty","L");
  tl->AddEntry(gr_energyscale,"Energy Scale Uncertainty","L");
  tl->AddEntry(effacc_syst_gr,"Efficiency Systematics","L");
  tl->Draw();
  
  zrap_Prelim(0.85,0.975,0.82,0.17);
  zrap_Lumi(0.25,0.18,lumi);

  sprintf(fnamework,"ZQt_final_errors%s.eps",postfix);
  c3->Print(fnamework);
  //  sprintf(fnamework,"ZQt_final_errors%s.png",postfix);
  // c3->Print(fnamework);
  sprintf(fnamework,"ZQt_final_errors%s.C",postfix);
  c3->Print(fnamework);

  sprintf(fnamework,"ZQt_final_errtable%s.tex",postfix);
  ftable=fopen(fnamework,"w");

  for (int i=0; i<lasti-1; i++)
  {

    if (i==0) {
      fprintf(ftable,"\\begin{tabular}{|c|c||c|c|c|c|c|}\\hline\n");
      fprintf(ftable,"            &             & Efficiency &   Bin     & Energy Scale  & Background & PDF\\\\ \n");
      fprintf(ftable,"$|q_{T,min}|$ & $|q_{T,max}|$ & Uncertainties & Migration & Uncertainty & Estimation & Uncertainty \\\\ \\hline \n");
    }

    double staterr=ea_statErr.y[i];
    double systerr=effSystErr.y[i];
    double efferr=corrData.y[i]*sqrt(staterr*staterr+systerr*systerr);
    double energyerr=energyScale.y[i]*corrData.y[i];
    double bkgderr=backgroundUncFrac.y[i]*corrData.y[i];
    double migerr=unfoldSyst.y[i]*corrData.y[i];
    double pdferr=pdfFrac.y[i]*corrData.y[i];
    
    //    fprintf(ftable," %7.2f & %7.2f & %.5f & %.5f & %.5f & %.5f \\\\ \n",
    fprintf(ftable," %7.2f & %7.2f & $%.3f \\times 10^{-4}$ &  $%.3f \\times 10^{-4}$ & $%.3f \\times 10^{-4}$ & $%.3f \\times 10^{-4}$ & $%.3f \\times 10^{-4}$ \\\\ \n",
	    corrData.xave[i]-corrData.xwidth[i],
	    corrData.xave[i]+corrData.xwidth[i],
	    efferr*10000,migerr*10000,energyerr*10000,bkgderr*10000,
	    pdferr*10000);
    
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);


  double chi2=0;
  int ndof=0;

  sprintf(fnamework,"dsdqt_zee_final%s.csv",postfix);
  FILE* outFile = fopen(fnamework, "wt");

  fprintf(outFile, "#%9s %9s %9s %9s %12s %12s\n", "bin", "lower", "upper", "data", "stat error", "syst error");
  for(int i = 0; i < BINCOUNT-1; i++)
    {
      fprintf(outFile, "%9i %9.1f %9.1f %9f %12.10f %12.10f\n", i+1, pt_binning[i], pt_binning[i+1], corrData.y[i], corrData.ey[i], sqrt(pow(corrDataSyst.ey[i],2)-pow(corrData.ey[i],2)));

      double toadd=pow(corrData.y[i]-truth_vis->GetBinContent(i+1),2)/pow(corrDataSyst.ey[i],2);
      //      printf("%d %f %f %f %f\n",i,corrData.y[i],truth_vis->GetBinContent(i+1),corrDataSyst.ey[i],toadd);
      if (i>5) {
	chi2+=toadd;
	ndof++;
      }
    }

  ndof--;
  fclose(outFile);
  

  printf(" Comparison: %f/%d %f \n",chi2,ndof,TMath::Prob(chi2,ndof));

}
