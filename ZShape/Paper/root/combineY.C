#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"
#include "TMath.h"

void combineY(TFile* theory,int mode=0) {
  TH1* truth=(TH1*)(theory->Get("mcEff/All/Z0_YTL")->Clone("truth"));
  truth->SetDirectory(0);  

  setTDRStyle();
  TH1* pdf_sens_low=readStandardFile("pdf_frac_low","pdfSensitivitiesNeg_Y.txt");
  TH1* pdf_sens_high=readStandardFile("pdf_frac_high","pdfSensitivitiesPos_Y.txt");

  TH1* elec_stat;
  TH1* elec_syst;
  if (mode==0) {

    elec_stat=readStandardFile("elec_stat","ZRapidity_fold_zee_foldMatrix.txt",0);
    elec_syst=readStandardFile("elec_all","ZRapidity_fold_zee_foldMatrix.txt",1);
  } else {
    elec_stat=readStandardFile("elec_stat","ZRapidity_fold_zee_avemig.txt",0);
    elec_syst=readStandardFile("elec_all","ZRapidity_fold_zee_avemig.txt",1);
  }

  TH1* mu_stat=readStandardFile("mu_stat","dsdy_muons.txt",0);
  TH1* mu_syst=readStandardFile("mu_syst","dsdy_muons.txt",1);

  TH1* elec_all=elec_syst->Clone("mu_all");
  TH1* mu_all=mu_syst->Clone("mu_all");
  TH1* comb_all=elec_syst->Clone("comb_all");

  TH1* pull_elec=new TH1F("pull_elec","pull_elec",50,-3,3);
  TH1* pull_mu=new TH1F("pull_mu","pull_mu",50,-3,3);

  TH1* pull_elec_y=new TH1F("pull_elec_Y","pull_elec_y",35,-0.01,3.49);
  TH1* pull_mu_y=new TH1F("pull_mu_Y","pull_mu_y",35,0.01,3.51);

  double td=0,tt=0;
  for (int i=51; i<=85; i++) {
    double total_err_elec=sqrt(pow(elec_stat->GetBinError(i),2)+pow(elec_syst->GetBinError(i),2));
    elec_all->SetBinError(i,total_err_elec);
    double bestEstimate,sumerr;
    double pelec=0,pmu=0;
    if (mu_stat->GetBinContent(i)>0.001) {
      double total_err_mu=sqrt(pow(mu_stat->GetBinError(i),2)+pow(mu_syst->GetBinError(i),2));
      if (mode==1) total_err_mu/=(1.05); // temporary hack for average-corrections
      mu_all->SetBinError(i,total_err_mu);      
      sumerr=1.0/(total_err_elec*total_err_elec)+1.0/(total_err_mu*total_err_mu);
    
      bestEstimate=(elec_stat->GetBinContent(i)/pow(total_err_elec,2)+mu_stat->GetBinContent(i)/pow(total_err_mu,2))/sumerr;

      
      pelec=(elec_stat->GetBinContent(i)-bestEstimate)/total_err_elec;
      pull_elec->Fill(pelec);
		   
      pmu=((mu_stat->GetBinContent(i)-bestEstimate)/total_err_mu);   
      pull_mu->Fill(pmu);

    } else {
      sumerr=1.0/(total_err_elec*total_err_elec);
    
      bestEstimate=elec_stat->GetBinContent(i);
    }
    sumerr=sqrt(1.0/sumerr);

    printf("%2d %.3f %.3f %.3f %.3f %.1f %.1f\n",i,elec_stat->GetBinContent(i),mu_stat->GetBinContent(i),bestEstimate,sumerr,pelec,pmu);
    comb_all->SetBinContent(i,bestEstimate);
    comb_all->SetBinError(i,sumerr);

    pull_elec_y->SetBinContent(i-50, elec_stat->GetBinContent(i)-bestEstimate);
    pull_elec_y->SetBinError(i-50,total_err_elec);

    if (mu_stat->GetBinContent(i)>0.001) {
      pull_mu_y->SetBinContent(i-50, mu_stat->GetBinContent(i)-bestEstimate);
      pull_mu_y->SetBinError(i-50,total_err_mu);      
    }

    td+=bestEstimate;
    tt+=truth->GetBinContent(i);
  }

  truth->Scale(1.0/(tt*0.1));

  // chi2
  double chi2=0;
  int ndof=-1; // for normalization
  for (int i=51;i<=85; i++) {
    chi2+=pow(comb_all->GetBinContent(i)-truth->GetBinContent(i),2)/pow(comb_all->GetBinError(i),2);
    ndof++;
  }

  static const int npoints=(85-51+1+2)*2;
  double terrorsx[npoints];
  double terrorsy[npoints];

  for (int i=51; i<=85; i++) {
    terrorsx[i-51+1]=truth->GetBinCenter(i);
    terrorsy[i-51+1]=truth->GetBinContent(i)*(1+pdf_sens_high->GetBinContent(i));

    terrorsx[123-i]=truth->GetBinCenter(i);
    terrorsy[123-i]=truth->GetBinContent(i)*(1+pdf_sens_low->GetBinContent(i));    
  }

  terrorsx[0]=0; terrorsy[0]=terrorsy[1];
  terrorsx[36]=3.5; terrorsy[36]=terrorsy[35];
  terrorsx[37]=3.5; terrorsy[37]=terrorsy[38];
  terrorsx[73]=0; terrorsy[73]=terrorsy[72];

  TGraph* terrors=new TGraph(npoints,terrorsx,terrorsy);
  terrors->SetFillColor(kCyan);
  terrors->SetFillStyle(1001);

  elec_all->GetXaxis()->SetRangeUser(0.0,3.45);
  comb_all->GetXaxis()->SetRangeUser(0.0,3.45);

  elec_all->SetMarkerStyle(25);
  mu_all->SetMarkerStyle(26);
  comb_all->SetMarkerStyle(21);

  TCanvas* c1=new TCanvas("c1","c1",800,800);


  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.17);
  comb_all->GetXaxis()->SetTitle("|y|");
  comb_all->GetXaxis()->CenterTitle();

  comb_all->GetYaxis()->SetTitle("1/#sigma d#sigma/d|y|");
  comb_all->GetYaxis()->SetTitleOffset(1.3);
  comb_all->GetYaxis()->CenterTitle();

  comb_all->Draw("E");
  terrors->Draw("F");

  comb_all->Draw("SAME E");

  //  truth->Draw("SAME HIST");

  zrap_Prelim(0.82,0.97,0.75,0.82);
  zrap_Lumi(0.77,0.90,36);

  TLegend* tl=new TLegend(0.20,0.25,0.75,0.40);
  tl->SetFillStyle(0);
  tl->AddEntry(comb_all,"data (e and #mu combined)","P");
  tl->AddEntry(terrors,"POWHEG + CT10","F");

  tl->Draw();

  c1->Print("dsday_combined.png");
  c1->Print("dsday_combined.eps");

  /*
  elec_all->Draw("E");
  mu_all->Draw("E SAME");
  comb_all->Draw("E SAME");
  */

  FILE* ftable=fopen("dsday_combined.tex","wt");
  fprintf(ftable,"\\begin{tabular}{|c|c||c|c|}\n\\hline\n");
  fprintf(ftable,"$|Y_{min}|$ & $|Y_{max}|$ & Measurement & Uncertainty  \\\\ \\hline\n");

  for (int i=51; i<=85; i++) {
    fprintf(ftable,"%7.1f & %7.1f & %7.4f & %7.4f \\\\ \n",
	    comb_all->GetXaxis()->GetBinLowEdge(i),
	    comb_all->GetXaxis()->GetBinUpEdge(i),
	    comb_all->GetBinContent(i),
	    comb_all->GetBinError(i));
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);

  TCanvas* c2=new TCanvas("c2","c2",1000,600);

  c2->Divide(2,1);
  c2->cd(1);
  pull_elec->Draw("HIST");
  c2->cd(2);
  pull_mu->Draw("HIST");

  TCanvas* c3=new TCanvas("c3","c3",800,550);
  c3->SetLeftMargin(0.14);

  pull_elec_y->SetMaximum(0.05);
  pull_elec_y->SetMinimum(-0.05);
  pull_elec_y->SetMarkerStyle(20);
  pull_elec_y->GetXaxis()->SetTitle("|y|");
  pull_elec_y->GetXaxis()->CenterTitle();

  pull_elec_y->GetYaxis()->SetTitle("|y_{ll}|-|y_{combined}|");
  pull_elec_y->GetYaxis()->CenterTitle();  
  pull_elec_y->GetYaxis()->SetTitleOffset(1.08);
  
  pull_mu_y->SetMarkerStyle(24);

  pull_elec_y->Draw("E");
  pull_mu_y->Draw("E SAME");

  TLegend* tlp=new TLegend(0.15,0.12,0.6,0.26);
  tlp->SetFillStyle(0);
  tlp->AddEntry(pull_elec_y,"e channel","P");
  tlp->AddEntry(pull_mu_y,"#mu channel","P");
  tlp->SetNColumns(2);
  tlp->Draw();

  zrap_Prelim(0.85,0.97,0.75,0.82);
  zrap_Lumi(0.55,0.90,36);
  
  c3->Print("dsday_combined_deltas.C");
  c3->Print("dsday_combined_deltas.eps");

  printf("Theory check: %f/%d (%f)\n",chi2,ndof,TMath::Prob(chi2,ndof));
}
