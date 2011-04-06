#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

void combineQt(int which=0) {
  TH1* truth=readStandardFilePt("theory","../theory/powheg_pt_etapt_scaleUnc.txt");
  truth->SetDirectory(0);  



  setTDRStyle();
  TH1* elec_stat,*elec_syst;

  if (which==1) {
    elec_stat=readStandardFilePt("elec_stat","dsdqt_zee_final_avemig.csv",0);
    elec_syst=readStandardFilePt("elec_syst","dsdqt_zee_final_avemig.csv",1);
  } else {
    elec_stat=readStandardFilePt("elec_stat","dsdqt_zee_final_matrix.csv",0);
    elec_syst=readStandardFilePt("elec_syst","dsdqt_zee_final_matrix.csv",1);
  }
  TH1* mu_all=readStandardFilePt("mu_syst","dsdqt_muons.txt",0);

  TH1* elec_all=elec_syst->Clone("mu_all");
  //  TH1* mu_all=mu_syst->Clone("mu_all");
  TH1* comb_all=elec_syst->Clone("comb_all");

  TH1* pull_elec=new TH1F("pull_elec","pull_elec",20,-3,3);
  TH1* pull_mu=new TH1F("pull_mu","pull_mu",20,-3,3);


  TGraphErrors* pull_elec_qt=new TGraphErrors(18);
  TGraphErrors* pull_mu_qt=new TGraphErrors(18);


  double td=0,tt=0;
  for (int i=1; i<=18; i++) {
    double total_err_elec=sqrt(pow(elec_stat->GetBinError(i),2)+pow(elec_syst->GetBinError(i),2));
    elec_all->SetBinError(i,total_err_elec);
    double bestEstimate,sumerr;
    double pelec=0,pmu=0;

    double total_err_mu=mu_all->GetBinError(i);
    sumerr=1.0/(total_err_elec*total_err_elec)+1.0/(total_err_mu*total_err_mu);
    
    bestEstimate=(elec_all->GetBinContent(i)/pow(total_err_elec,2)+mu_all->GetBinContent(i)/pow(total_err_mu,2))/sumerr;

    sumerr=sqrt(1.0/sumerr);
    
    pelec=(elec_all->GetBinContent(i)-bestEstimate)/total_err_elec;
    pull_elec->Fill(pelec);

    pull_elec_qt->SetPoint(i-1,
			   (elec_all->GetBinCenter(i)*1.02),
			   (elec_all->GetBinContent(i)/bestEstimate)-1);

    pull_elec_qt->SetPointError(i-1,
				0,
				total_err_elec/bestEstimate);
    
    pmu=((mu_all->GetBinContent(i)-bestEstimate)/total_err_mu);   
    pull_mu->Fill(pmu);

    pull_mu_qt->SetPoint(i-1,
			 (mu_all->GetBinCenter(i)/1.02),
			   (mu_all->GetBinContent(i)/bestEstimate)-1);

    pull_mu_qt->SetPointError(i-1,
			      0,
			      total_err_mu/bestEstimate);

    /*
    pull_mu_qt->SetBinContent(i,(mu_all->GetBinContent(i)/bestEstimate)-1);
    pull_mu_qt->SetBinError(i,(total_err_mu)/bestEstimate);
    */


    printf("%2d %.5e %.5e %.5e %.5e %.1f %.1f\n",i,elec_stat->GetBinContent(i),mu_all->GetBinContent(i),bestEstimate,sumerr,pelec,pmu);
    comb_all->SetBinContent(i,bestEstimate);
    comb_all->SetBinError(i,sumerr);
    td+=bestEstimate;
    tt+=truth->GetBinContent(i);
  }

  //  truth->Scale(1.0/(tt));
  /*
  for (int i=1; i<=18; i++) {
    truth->SetBinContent(i,truth->GetBinContent(i)/truth->GetBinWidth(i));
  }
  */
  TH1* truth2=truth->Clone("truth2");
  TGraphErrors* truth3=new TGraphErrors(18);

  for (int i=1; i<=18; i++) {
    truth2->SetBinContent(i,truth->GetBinContent(i)/comb_all->GetBinContent(i)-1);
    truth2->SetBinError(i,comb_all->GetBinError(i));

    truth3->SetPoint(i-1,truth2->GetBinCenter(i),truth2->GetBinContent(i));
    truth3->SetPointError(i-1,truth2->GetBinWidth(i)/2,truth->GetBinError(i)/truth->GetBinContent(i));

  }

  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.17);  
  c1->SetBottomMargin(0.16);
  comb_all=zpt_rebinForPlot(comb_all);
  comb_all->SetMaximum(0.08);
  elec_all->GetXaxis()->SetRangeUser(0.7,500);
  comb_all->GetXaxis()->SetRangeUser(0.7,500);

  elec_all->SetMarkerStyle(25);
  mu_all->SetMarkerStyle(20);
  comb_all->SetMarkerStyle(21);
  comb_all->SetTitle(0);
  comb_all->GetYaxis()->SetTitleOffset(1.3);
  comb_all->GetXaxis()->SetTitle(qt_xaxis_ll_label);
  comb_all->GetXaxis()->CenterTitle();
  comb_all->GetXaxis()->SetTitleOffset(1.2);

  comb_all->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T,ll}");
  comb_all->GetYaxis()->CenterTitle();

  comb_all->Draw("E");
  truth=zpt_rebinForPlot(truth);

  truth->SetMarkerSize(0);
  truth->SetLineColor(kGreen+2);
  truth->SetLineWidth(1);

  truth->Draw("SAME][ HIST");
  /*
  TH1* t3=truth->Clone("t3");
  t3->SetFillStyle(1001);
  t3->SetFillColor(kGreen);
  t3->Draw("E2");

  TH1* t4=truth->Clone("t4");
  t4->SetFillStyle(3256);
  t4->SetFillColor(kGreen+2);
  t4->Draw("E2");
  */
  zrap_Prelim(0.29,0.97,0.75,0.82);
  zrap_Lumi(0.80,0.97,36);

  TLegend* tl=new TLegend(0.20,0.25,0.75,0.40);
  tl->SetFillStyle(0);
  //  tl->AddEntry(comb_all,"Z#rightarrow e^{+}e^{-} and Z#rightarrow #mu^{+}#mu^{-} combined","P");
  tl->AddEntry(comb_all,"e^{+}e^{-} and #mu^{+}#mu^{-} combined","P");
  tl->AddEntry(truth,"POWHEG + CT10","L");

  tl->Draw();

  c1->SetLogx();

  c1->Print("dsdqt_combined_linear.C");
  c1->Print("dsdqt_combined_linear.eps");

  
  c1->SetLogy();

  c1->Print("dsdqt_combined_log.C");
  c1->Print("dsdqt_combined_log.eps");
  

  /*
  elec_all->Draw("E");
  mu_all->Draw("E SAME");
  comb_all->Draw("E SAME");
  */

  int ndof=-1;
  double chi2=0;
  
  FILE* ftable=fopen("dsdqt_combined.tex","wt");
  FILE* fcsv=fopen("dsdqt_combined.csv","wt");
  fprintf(ftable,"\\begin{tabular}{|c|c||l|l|}\n\\hline\n");
  fprintf(ftable,"$p_{\\mathrm{T,min}}$ & $q_{\\mathrm{T,max}}$ & Measurement & Uncertainty  \\\\ \\hline\n");
  TDatime now;
  fprintf(fcsv,"# combined result electron/muon as of %s\n",now.AsString());
  fprintf(fcsv,"#bin qt_min qt_max value error\n");
  for (int i=1; i<=18; i++) {
    double bl=comb_all->GetXaxis()->GetBinLowEdge(i);
    double bh=comb_all->GetXaxis()->GetBinUpEdge(i);
    if (i==1) bl=0;
    int prec=int(-log10(comb_all->GetBinError(i)))+2;
    fprintf(ftable,"%7.1f & %7.1f & %7.*f & %7.*f \\\\ \n",
	    bl,bh,
	    prec,comb_all->GetBinContent(i),
	    prec,comb_all->GetBinError(i));
    fprintf(fcsv, "%4d %5.1f %5.1f %.9f %.9f\n",
	    i,bl,bh,
	    comb_all->GetBinContent(i),
	    comb_all->GetBinError(i));
    ndof++;

    double toadd=pow(comb_all->GetBinContent(i)-truth->GetBinContent(i),2)/pow(comb_all->GetBinError(i),2);
    chi2+=toadd;
    printf("%d %f\n",i,toadd);
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);
  fclose(fcsv);

  TCanvas* c2=new TCanvas("c2","c2",1000,600);
  truth2=zpt_rebinForPlot(truth2);

  c2->Divide(2,1);
  c2->cd(1);
  pull_elec->Draw("HIST");
  c2->cd(2);
  pull_mu->Draw("HIST");

  TCanvas* c3=new TCanvas("c3","c3",800,800);
  c3->SetLeftMargin(0.18);
  c3->SetBottomMargin(0.16);

  truth2->SetMaximum(0.9);
  truth2->SetMinimum(-0.9);
  truth2->SetMarkerStyle(20);
  truth2->GetXaxis()->SetTitle();
  truth2->SetLineColor(kBlue);
  truth2->SetLineWidth(2);
  truth3->SetLineColor(kBlue);
  truth3->SetLineWidth(2);
  truth2->GetXaxis()->SetRangeUser(0.7,500.0);
  truth2->GetXaxis()->CenterTitle();
  truth2->SetFillColor(kCyan);
  c3->SetLogx();

  truth2->GetYaxis()->SetTitle("(p_{T,ee/#mu#mu}/p_{T,ll}) - 1");
  truth2->GetXaxis()->SetTitle(qt_xaxis_ll_label);
  truth2->GetYaxis()->CenterTitle();  
  truth2->GetXaxis()->SetTitleOffset(1.2);
  truth2->GetYaxis()->SetTitleOffset(1.4);
  
  pull_mu_qt->SetMarkerStyle(20);
  pull_elec_qt->SetMarkerStyle(25);
  
  truth2->GetXaxis()->SetRangeUser(0.7,500);
  truth2->SetTitle(0);
  truth2->Draw("HIST ][");
  
  truth3->SetFillColor(kCyan);
  truth3->Draw("2 SAME");
  
  pull_elec_qt->Draw("P E1 SAME");
  pull_mu_qt->Draw("P E1 SAME");

  truth2->Draw("HIST SAME ][");

  TLegend* tlp=new TLegend(0.20,0.75,0.75,0.85);
  tlp->SetTextFont(42);
  tlp->SetFillStyle(0);
  tlp->AddEntry(pull_elec_qt,"data (e)","P");
  tlp->AddEntry(pull_mu_qt,"data (#mu)","P");
  tlp->SetNColumns(2);
  tlp->Draw();

  TLegend* tlp2=new TLegend(0.25,0.66,0.70,0.75);
  tlp2->SetFillStyle(0);
  tlp2->SetTextFont(42);
  tlp2->AddEntry(truth3,"POWHEG / p_{T,ll} - 1","LF");
  tlp2->Draw();

  zrap_Prelim(0.82,0.975,0.35,0.82);
  zrap_Lumi(0.35,0.90,36);
  zqt_Cut(0.4,0.23);
  c3->Print("dsdqt_combined_deltas.eps");
  c3->Print("dsdqt_combined_deltas.C");

  printf(" Comparison : %f/%d %f \n",chi2,ndof,TMath::Prob(chi2,ndof));

}
