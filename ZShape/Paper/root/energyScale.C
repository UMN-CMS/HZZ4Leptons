#include "readStandardFile.C"
#include "readUnfoldingMatrices.C"

TH1* extractDistribution(const char* name, TFile* f, int mode) {

  
  TH1* effAcc;
  if (mode==2) effAcc=readStandardFile("effAcc","../data/effAcc_smearOverTrue.csv");
  else effAcc=readStandardFile("effAcc","../data/effAcc_combined.csv");
  
  const char* h1name="mcEff/ECAL80-ECAL95/C08-m(60,120)/Z0_Y";
  const char* h2name="mcEff/ECAL80-HF/C08-m(60,120)/Z0_Y";

  TH1* sig1=f->Get(h1name)->Clone(name);
  TH1* sig2=f->Get(h2name);
  sig1->Sumw2();

  sig1->Add(sig2);
  sig1->Divide(effAcc);

  delete effAcc;

  return sig1;
}

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void energyScale(TFile* base, TFile* vpe, TFile* vme,  
		 TFile* vpf, TFile* vmf, 
		 TFile* vpv, TFile* vmv, int mode,
		 const char* errfilename="energyScaleError.csv" ) {

  const char* smode=0;
  switch (mode) {
  case (1) : smode="No unsmearing, EffAcc is smear/smear"; break;
  case (2) : smode="No unsmearing, EffAcc is smear/true"; break;
  case (3) : smode="Unsmearing"; break;
  default : return;
  }

  FILE* errfile=0;
  if (errfilename!=0) {
    errfile=fopen(errfilename,"wt");
    TDatime now;
    fprintf(errfile,"# Energy Scale Error # %s\n",now.AsString());
    fprintf(errfile,"#Bin Y_Min Y_Max Value Uncertainty\n");
  }

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TH1* hb=extractDistribution("base",base,mode);
  TH1* valep=extractDistribution("vpe",vpe,mode);
  TH1* valem=extractDistribution("vme",vme,mode);
  TH1* valfp=extractDistribution("vpf",vpf,mode);
  TH1* valfm=extractDistribution("vmf",vmf,mode);
  TH1* valvp=extractDistribution("vpv",vpv,mode);
  TH1* valvm=extractDistribution("vmv",vmv,mode);
  
  valep->Scale(hb->Integral(16,85)/valep->Integral(16,85));
  valem->Scale(hb->Integral(16,85)/valem->Integral(16,85));

  valfp->Scale(hb->Integral(16,85)/valfp->Integral(16,85));
  valfm->Scale(hb->Integral(16,85)/valfm->Integral(16,85));

  valvp->Scale(hb->Integral(16,85)/valvp->Integral(16,85));
  valvm->Scale(hb->Integral(16,85)/valvm->Integral(16,85));

  if (mode==3) {
    hb=unfold(hb,"unfoldingMatrix_theOutPut.root");
    valep=unfold(valep,"unfoldingMatrix_theOutPut.root");
    valem=unfold(valem,"unfoldingMatrix_theOutPut.root");
    valfp=unfold(valfp,"unfoldingMatrix_theOutPut.root");
    valfm=unfold(valfm,"unfoldingMatrix_theOutPut.root");
    valvp=unfold(valvp,"unfoldingMatrix_theOutPut.root");
    valvm=unfold(valvm,"unfoldingMatrix_theOutPut.root");
  }

  valep->Add(hb,-1);
  valem->Add(hb,-1);

  valep->Divide(hb);
  valem->Divide(hb);
  
  valfp->Add(hb,-1);
  valfm->Add(hb,-1);
  
  valfp->Divide(hb);
  valfm->Divide(hb);

  valvp->Add(hb,-1);
  valvm->Add(hb,-1);
  
  valvp->Divide(hb);
  valvm->Divide(hb);

  TH1* valm=valem->Clone("valm");
  TH1* valp=valep->Clone("valp");
  TH1* valave=valep->Clone("valave");

  for (int i=1; i<=valep->GetXaxis()->GetNbins(); i++) {
    valm->SetBinContent(i,sqrt(pow(valem->GetBinContent(i),2)+
			       pow(valvm->GetBinContent(i),2)+
			       pow(valfm->GetBinContent(i),2)));
    valp->SetBinContent(i,sqrt(pow(valep->GetBinContent(i),2)+
			       pow(valvp->GetBinContent(i),2)+
			       pow(valfp->GetBinContent(i),2)));
    valave->SetBinContent(i,0.5*(
				 valm->GetBinContent(i)+
				 valp->GetBinContent(i)
				 ));
    valave->SetBinError(i,0.5*sqrt(
			       pow(valm->GetBinError(i),2)+
			       pow(valp->GetBinError(i),2)
				 ));

    if (i>=16 && i<=85) {
      if (errfile!=0)
	fprintf(errfile,"%4d %5.2f %5.2f %6.4f %6.4f\n",i,
		valave->GetXaxis()->GetBinLowEdge(i),
		valave->GetXaxis()->GetBinUpEdge(i),
		valave->GetBinContent(i),
		valave->GetBinError(i)
		);
    }

  }

  TCanvas* c1=0;
  if (mode==1 || mode==2) {
    c1=new TCanvas("c1","c1",800,800);
    c1->SetLeftMargin(0.15);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetBottomMargin(0.10);
  }

  //  valm->SetLineColor(kBlue);
  //  valp->SetLineColor(kRed);

  //  valm->SetLineWidth(2);
  //  valp->SetLineWidth(2);


  valave->SetLineWidth(2);
  valave->SetMaximum(0.05);
  valave->SetMinimum(0.0);
  
  valave->SetTitle(0);
  valave->SetStats(0);
  valave->GetXaxis()->SetRangeUser(-3.5,3.5);
  valave->GetXaxis()->SetTitle("y_{Z}");
  valave->GetXaxis()->CenterTitle(true);
  valave->GetYaxis()->SetTitle("Fractional Error from Energy Scale");
  valave->GetYaxis()->SetTitleOffset(1.8);
  valave->GetYaxis()->CenterTitle(true);
  if (mode==1 || mode==2)   valave->Draw("HIST");
  else   valave->Draw("HISTSAME");
  //  valm->Draw("HIST SAME");

  zrap_Prelim(0.6,0.9,0.6,0.87);
  if (c1!=0) c1->Print("energyScale.eps");

  if (errfile!=0) fclose(errfile);

}


