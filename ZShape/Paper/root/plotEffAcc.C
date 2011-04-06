EColor colorFor(int i) {
  switch (i) {
  case (0): return kBlack;
  case (1): return kRed;
  case (2): return kBlue;
  case (3): return kRed;
  case (4): return kOrange;
  case (5): return kViolet-5;
  case (6): return kOrange;
  default: return i+1;
  };

}

const char* labelFor(int i) {
  switch (i) {
  case(0) : return "Geometrical Acceptance";
  case(1) : return "Supercluster";
  case(2) : return "Supercluster P_{T}>20 GeV";
  case(3) : return "Electron Reconstruction (Track)";
  case(5) : return "Electron Id";
  case(6) : return "HLT";
    //  case(6) : return "HLT";
  default : return "Huh?";
  }
}

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void plotEffAcc(TFile* f, int mode=0, const char* var="Z0_Y_masscut") {
  TDirectory* basedir=f->Get("mcEff");

  TH1* base;
  const char* post=0;
  if (mode==0) {
    if (strstr(var,"Y")!=0) base=(TH1*)((TDirectory*)basedir->Get("All"))->Get("Z0_Y");
  } else {
    post="smeartrue";
    if (strstr(var,"Y")!=0) base=(TH1*)((TDirectory*)basedir->Get("All"))->Get("Z0_YTL");
  }
  if (strstr(var,"Pt")!=0) {
    base=(TH1*)((TDirectory*)basedir->Get("All"))->Get("Z0_Pt");
  }

  TH1* baseclone=base->Clone("BaseClone");
  baseclone->Scale(0.64/(0.5*(baseclone->GetBinContent(50)+baseclone->GetBinContent(51))));

  TDirectory* zdir=basedir->Get("ECAL80-ECAL95");
  TDirectory* zdir2=basedir->Get("ECAL80-HF");

  TList* cutz=zdir->GetListOfKeys();
  TList* cutz2=zdir2->GetListOfKeys();

  gROOT->SetStyle("Plain");
  setTDRStyle();

  TCanvas* c1=new TCanvas("c1","c1",700,600);
  TCanvas* c2=new TCanvas("c2","c2",700,600);
  
  c1->cd();

  TLegend* tl;
  tl=new TLegend(0.13,0.8,0.47,0.92);
  tl->SetFillStyle(0);

  TLegend* tl2=new TLegend(0.50,0.8,0.95,0.92);
  tl2->SetFillStyle(0);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.02);
  c1->SetLeftMargin(0.10);
  c1->SetBottomMargin(0.10);
  c2->SetTopMargin(0.05);
  c2->SetRightMargin(0.02);
  c2->SetLeftMargin(0.10);
  c2->SetBottomMargin(0.10);

  TH1* final=0,*final2=0,*finalt=0;

  int n=0;
  for (int i=0; i<cutz->GetSize(); i++) {
    TH1* stage=((TDirectory*)zdir->Get(cutz->At(i)->GetName()))->Get(var);
    TH1* stage2=((TDirectory*)zdir2->Get(cutz2->At(i)->GetName()))->Get(var);
    char cleanname[256];
    char* fullname=cutz->At(i)->GetName();
    char* pt=strchr(cutz->At(i)->GetName(),'-');
    if (pt!=0) strcpy(cleanname,pt+1);
    else strcpy(cleanname,cutz->At(i)->GetName());

    pt=strchr(cleanname,'-'); if (pt!=0) *pt=0;

    if (!strcmp(cleanname,"PT10")) continue;
    if (!strcmp(fullname,"C05-PT20")) continue;
    //    if (strstr(fullname,"-m(")!=0) continue;
    // if (!strcmp(cleanname,"PT20")) continue;
    printf("%d %s %s %p %p\n",n,fullname,cleanname,stage,pt);
    TH1* work=stage->Clone();
    work->Sumw2();
    final=(TH1*)(stage->Clone());
    final->Divide(final,base);

    final2=(TH1*)(stage2->Clone());
    final2->Divide(final2,base);

    work->Add(stage2);
    TH1* work2=base->Clone();

    double w2f=0;
    int w2fn=0;
    for (int iii=-3; iii<=4; iii++) {
      w2f+=work2->GetBinContent(work2->GetNbinsX()/2+iii);
      w2fn++;
    }
    w2f/=w2fn;

    work2->Scale(1.0/w2f);
    work->Divide(work,base);
    work->SetLineColor(colorFor(n));
    work->SetLineWidth(2);
    if (i==0) {
      work->GetYaxis()->SetTitle("Efficiency #times Acceptance");
      work->GetXaxis()->SetTitle("y_{ee}");
      work->GetYaxis()->SetTitleOffset(1.2);
      work->SetTitle(0);
      work->SetStats(false);
      work->SetMaximum(1.2);
      work->SetMinimum(0.0);
      if (strstr(var,"Y")!=0) {
	work->GetXaxis()->SetRangeUser(-4.5,4.5);
      } else {
	work->GetXaxis()->SetRangeUser(0.0,150.0);
      }
      work->GetXaxis()->CenterTitle(true);
      work->GetYaxis()->CenterTitle(true);
      work->Draw("HIST");
    } else if (n!=4 && n!=7 && n!=1) work->Draw("SAMEHIST");

    if (n!=4 && n!=7 && n!=1) {
      if (n<3) 
	tl->AddEntry(work,labelFor(n));
      else
	tl2->AddEntry(work,labelFor(n));
    }
    n++;
    finalt=(TH1*)(work->Clone());
  }
  tl->SetFillColor(0);
  tl->Draw("SAME");
  tl2->SetFillColor(0);
  tl2->Draw("SAME");

  zrap_Prelim(0.87,0.975,0.3,0.87);

  char stuff[1024];
  if (post==0) sprintf(stuff,"effacc_bycut.png");
  else sprintf(stuff,"effacc_bycut_%s.png",post);
  c1->Print(stuff);
  if (post==0) sprintf(stuff,"effacc_bycut.eps");
  else sprintf(stuff,"effacc_bycut_%s.eps",post);
  c1->Print(stuff);

  c2->cd();

  final->SetLineColor(kBlue);
  final2->SetLineColor(kRed);
  finalt->SetLineColor(kBlack);
  final->SetLineWidth(2);
  final2->SetLineWidth(2);
  finalt->SetLineWidth(2);

  final->GetYaxis()->SetTitle("Efficiency #times Acceptance");
  final->GetXaxis()->SetTitle("y_{ee}");

  final->GetYaxis()->SetTitleOffset(1.2);
  final->SetTitle(0);
  final->SetStats(false);
  final->SetMaximum(0.90);
  final->SetMinimum(0.0);
  if (strstr(var,"Y")!=0) {
    final->GetXaxis()->SetRangeUser(-4.5,4.5);
  } else {
    final->GetXaxis()->SetRangeUser(0.0,120.0);
  }
  final->GetXaxis()->CenterTitle(true);
  final->GetYaxis()->CenterTitle(true);
  
  final->Draw("HIST");
  final2->Draw("SAMEHIST");
  finalt->Draw("SAMEHIST");

  tl3=new TLegend(0.12,0.74,0.95,0.95);
  tl3->SetNColumns(2);
  tl3->SetFillStyle(0);
  tl3->AddEntry(final,"ECAL-ECAL");
  tl3->AddEntry(final2,"ECAL-HF");
  tl3->AddEntry(finalt,"Total");
  tl3->Draw();

  zrap_Prelim(0.87,0.975,0.3,0.87);

  if (post==0) sprintf(stuff,"effacc_bytype.png");
  else sprintf(stuff,"effacc_bytype_%s.png",post);
  c2->Print(stuff);
  if (post==0) sprintf(stuff,"effacc_bytype.eps");
  else sprintf(stuff,"effacc_bytype_%s.eps",post);
  c2->Print(stuff);

  baseclone->SetLineWidth(2);
  baseclone->SetLineStyle(2);
  baseclone->SetLineColor(kGreen+1);
  tl3->AddEntry(baseclone,"POWHEG+CT10 Prediction");

  baseclone->Draw("SAMEHIST");

  if (post==0) sprintf(stuff,"effacc_bytypetheory.eps");
  else sprintf(stuff,"effacc_bytypetheory_%s.eps",post);
  c2->Print(stuff);
  if (post==0) sprintf(stuff,"effacc_bytypetheory.C");
  else sprintf(stuff,"effacc_bytypetheory_%s.C",post);
  c2->Print(stuff);

  


  if (post==0) sprintf(stuff,"effacc.txt");
  else sprintf(stuff,"effacc_%s.txt",post);
  FILE* stats=fopen(stuff,"wt");

  TDatime now;

  if (mode==0)  fprintf(stats,"# Smear/Smear\n");
  else  fprintf(stats,"# Smear/True\n");
  fprintf(stats,"#Bin Y_min Y_max Eff*Acc MCErr#  %s \n",now.AsString());

  for (int i=1; i<=finalt->GetNbinsX(); i++) {
    fprintf(stats," %2d %5.2f %5.2f %7.5f %7.5f\n",i,finalt->GetBinLowEdge(i),finalt->GetBinWidth(i)+finalt->GetBinLowEdge(i),finalt->GetBinContent(i),finalt->GetBinError(i)); 
  }
  fclose(stats);

}
