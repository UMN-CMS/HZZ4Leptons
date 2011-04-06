EColor colorFor(int i) {
  switch (i) {
  case (0): return kBlack;
  case (1): return kRed;
  case (2): return kGreen+1;
  case (3): return kBlue;
  case (4): return kOrange;
  case (5): return kViolet-5;
  case (6): return kRed-6;
  default: return i+1;
  };

}

const char* labelFor(int i) {
  switch (i) {
  case(0) : return "Geometrical Acceptance";
  case(1) : return "Supercluster";
  case(2) : return "Supercluster with E_{T}>20 GeV";
  case(3) : return "Electron Reconstruction (Track)";
  case(4) : return "SKIPME";
  case(5) : return "Electron ID";
  case(6) : return "HLT";
  default : return "Huh?";
  }
}

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void pt_plotEffAcc(TFile* f, int mode=1) {

  const char* var="Z0_Pt_masscut";
  TDirectory* basedir=f->Get("mcEff");

  TH1* base=0, *base0=0;

  if (mode==1) {
    base=(TH1*)(((TDirectory*)basedir->Get("MUON-MUON-G20/C01-GPT20"))->Get(var));
    base0=(TH1*)(((TDirectory*)basedir->Get("MUON-MUON-G20/Acceptance"))->Get(var));
  } else
    base=(TH1*)(((TDirectory*)basedir->Get("MUON-MUON/Acceptance"))->Get("Z0_PtTL_masscut"));

  TDirectory* zdir=basedir->Get("ECAL80-ECAL95-MUO");
 
  TList* cutz=zdir->GetListOfKeys();
 

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TCanvas* c1=new TCanvas("c1","c1",700,600);
  
  c1->cd();

  TLegend* tl;
  tl=new TLegend(0.15,0.2,0.75,0.45);
  tl->SetFillStyle(0);

  TLegend* tl2=new TLegend(0.50,0.2,0.92,0.32);
  tl2->SetFillStyle(0);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.15);
  //  c1->SetRightMargin(0.05);

  TH1* final=0,*finalt=0;

  int n=0;
  for (int i=0; i<cutz->GetSize(); i++) {
    TH1* stage=((TDirectory*)zdir->Get(cutz->At(i)->GetName()))->Get(var);

    char cleanname[256];
    char* pt=strchr(cutz->At(i)->GetName(),'-');
    if (pt!=0) strcpy(cleanname,pt+1);
    else strcpy(cleanname,cutz->At(i)->GetName());

    pt=strchr(cleanname,'-'); if (pt!=0) *pt=0;

    if (!strcmp(cleanname,"PT10")) continue;
    if (strstr(cleanname,"m(")!=0) continue;
    // if (!strcmp(cleanname,"PT20")) continue;
    printf("%d %s %s %p %p\n",n,cutz->At(i)->GetName(),cleanname,stage,pt);
    TH1* work=zpt_rebinForPlot(stage);
    TH1* rebase=zpt_rebinForPlot(base);
    final=(TH1*)(stage->Clone());
    if (i==0)  {
      final->Divide(final,base0);
      rebase=zpt_rebinForPlot(base0);
    } else {
      final->Divide(final,base);
    }



    TH1* work2=base->Clone();

    double w2f=0;
    int w2fn=0;
    for (int iii=-3; iii<=4; iii++) {
      w2f+=work2->GetBinContent(work2->GetNbinsX()/2+iii);
      w2fn++;
    }
    w2f/=w2fn;

    work2->Scale(1.0/w2f);
    work->Divide(work,rebase);
    work->SetLineColor(colorFor(n));
    work->SetLineWidth(2);
    if (i==0) {
      work->GetYaxis()->SetTitle("Efficiency #times Acceptance");
      work->GetXaxis()->SetTitle(qt_xaxis_label);
      work->GetYaxis()->SetTitleOffset(1.05);
      work->SetTitle(0);
      work->SetStats(false);
      work->SetMaximum(1.00);
      work->SetMinimum(0.5);
      if (strstr(var,"Y")!=0) {
	work->GetXaxis()->SetRangeUser(-4.5,4.5);
      } else {
	work->GetXaxis()->SetRangeUser(0.0,500.0);
      }
      work->GetXaxis()->CenterTitle(true);
      work->GetYaxis()->CenterTitle(true);
      work->GetXaxis()->SetNoExponent(true);
      work->GetXaxis()->SetMoreLogLabels(false);
      work->Draw("HIST ][");
    } else if (n!=4 && n!=1) 
      work->Draw("SAMEHIST ][");

    if (n==4 || n==1) { 
    } else if (n<3) 
      tl->AddEntry(work,labelFor(n),"L");
    else
      tl->AddEntry(work,labelFor(n),"L");
    
    n++;
    finalt=(TH1*)(work->Clone());
  }
  tl->SetFillColor(0);
  tl->Draw("SAME");
  //  tl2->SetFillColor(0);
  //tl2->Draw("SAME");

  c1->SetLogx();

  zrap_Prelim(0.84,0.975,0.8,0.74);
  zqt_Cut(0.73,0.23,"e");

  char stuff[1024];
  //  sprintf(stuff,"eff_pt_bycut.png");
  //c1->Print(stuff);
  if (mode==1) sprintf(stuff,"eff_pt_bycut_smear.eps");
  if (mode==2) sprintf(stuff,"eff_pt_bycut_denomtl.eps");
  c1->Print(stuff);

  /*
  c2->cd();

  final->SetLineColor(kBlue);
  final2->SetLineColor(kRed);
  finalt->SetLineColor(kBlack);
  final->SetLineWidth(2);
  final2->SetLineWidth(2);
  finalt->SetLineWidth(2);

  final->GetYaxis()->SetTitle("#epsilon");
  final->GetXaxis()->SetTitle("Y_{Z}");

  final->GetYaxis()->SetTitleOffset(1.2);
  final->SetTitle(0);
  final->SetStats(false);
  final->SetMaximum(0.65);
  final->SetMinimum(0.0);
  if (strstr(var,"Y")!=0) {
    final->GetXaxis()->SetRangeUser(-4.5,4.5);
  } else {
    final->GetXaxis()->SetRangeUser(0.0,500.0);
  }
  
  final->GetXaxis()->CenterTitle(true);
  final->GetYaxis()->CenterTitle(true);
  
  final->Draw("HIST");
  final2->Draw("SAMEHIST");
  finalt->Draw("SAMEHIST");

  tl3=new TLegend(0.65,0.75,0.9,0.9);
  tl3->SetFillStyle(0);
  tl3->AddEntry(final,"ECAL-ECAL");
  tl3->AddEntry(final2,"ECAL-HF");
  tl3->AddEntry(finalt,"Total");
  tl3->Draw();

  zrap_Prelim(0.3,0.9,0.3,0.87);

  if (post==0) sprintf(stuff,"effacc_bytype.png");
  else sprintf(stuff,"effacc_bytype_%s.png",post);
  c2->Print(stuff);
  if (post==0) sprintf(stuff,"effacc_bytype.eps");
  else sprintf(stuff,"effacc_bytype_%s.eps",post);
  c2->Print(stuff);

  */


  if (mode==1) {
    sprintf(stuff,"effAcc_pt_massCutdenom.txt");
  } else if (mode==2) {
    sprintf(stuff,"effAcc_pt_TLdenom.txt");
  } else return;
  FILE* stats=fopen(stuff,"wt");

  TDatime now;

  fprintf(stats,"#Bin pt_min pt_max Eff %s\n",now.AsString());

  for (int i=1; i<=finalt->GetNbinsX(); i++) {
    fprintf(stats," %2d %5.2f %5.2f %.5f\n",i,finalt->GetBinLowEdge(i),finalt->GetBinWidth(i)+finalt->GetBinLowEdge(i),finalt->GetBinContent(i)); 
  }
  fclose(stats);

}
