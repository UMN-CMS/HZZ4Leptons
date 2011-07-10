#include "readStandardFile.C"
#include "tdrstyle.C"
#include "zrapidityStandard.C"

void combineQt(int which=0) {
  //  TH1* truth=readStandardFilePt("theory","../theory/powheg_pt_etapt_scaleUnc.txt");
  TH1* truth=readStandardFilePt("theory","../theory/PtTL_pwg_PythiaUEZ2.txt");
  truth->SetDirectory(0);  



  setTDRStyle();

  TH1* mu_all=readStandardFilePt("mu_all","dsdqt_muons.txt",0);

  TH1* elec_all=readStandardFilePt("elec_all","dsdqt_zee_final_matrix.csv",0);

  //  TH1* mu_all=mu_syst->Clone("mu_all");
  TH1* combo_all=readStandardFilePt("combo_all","dsdqt_combined.csv",0);

  TH1* pull_elec=new TH1F("pull_elec","pull_elec",20,-3,3);
  TH1* pull_mu=new TH1F("pull_mu","pull_mu",20,-3,3);


  TGraphErrors* pull_elec_qt=new TGraphErrors(18);
  TGraphErrors* pull_mu_qt=new TGraphErrors(18);


  double td=0,tt=0;
  for (int i=1; i<=18; i++) {
    double bestEstimate,sumerr;
    double pelec=0,pmu=0;

    double total_err_elec=elec_all->GetBinError(i);
    double total_err_mu=mu_all->GetBinError(i);

    sumerr=1.0/(total_err_elec*total_err_elec)+1.0/(total_err_mu*total_err_mu);
    
    bestEstimate=combo_all->GetBinContent(i);

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

    printf("%2d %.5e %.5e %.5e %.5e %.1f %.1f\n",i,elec_all->GetBinContent(i),mu_all->GetBinContent(i),bestEstimate,sumerr,pelec,pmu);
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
    truth2->SetBinContent(i,truth->GetBinContent(i)/combo_all->GetBinContent(i)-1);
    truth2->SetBinError(i,combo_all->GetBinError(i));

    truth3->SetPoint(i-1,truth2->GetBinCenter(i),truth2->GetBinContent(i));
    truth3->SetPointError(i-1,truth2->GetBinWidth(i)/2,truth->GetBinError(i)/truth->GetBinContent(i));

  }

  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetLeftMargin(0.17);  
  c1->SetBottomMargin(0.16);
  combo_all=zpt_rebinForPlot(combo_all);
  combo_all->SetMaximum(0.1);
  elec_all->GetXaxis()->SetRangeUser(0.7,500);
  combo_all->GetXaxis()->SetRangeUser(0.7,500);

  elec_all->SetMarkerStyle(25);
  mu_all->SetMarkerStyle(20);
  combo_all->SetMarkerStyle(21);
  combo_all->SetTitle(0);
  combo_all->GetYaxis()->SetTitleOffset(1.3);
  combo_all->GetXaxis()->SetTitle(qt_xaxis_ll_label);
  combo_all->GetXaxis()->CenterTitle();
  combo_all->GetXaxis()->SetTitleOffset(1.2);

  combo_all->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T,ll}");
  combo_all->GetYaxis()->CenterTitle();

  combo_all->Draw("E");
  truth=zpt_rebinForPlot(truth);

  truth->SetMarkerSize(0);
  truth->SetLineColor(kGreen+2);
  truth->SetLineWidth(1);

  TH1* t3=truth->Clone("t3");
  t3->SetFillStyle(1001);
  t3->SetFillColor(kGreen);
  t3->Draw("E2 SAME");

  truth->Draw("SAME ][ HIST");
  combo_all->Draw("E SAME");
  /*
  TH1* t4=truth->Clone("t4");
  t4->SetFillStyle(3256);
  t4->SetFillColor(kGreen+2);
  t4->Draw("E2");
  */
  zrap_Prelim(0.82,0.975,0.75,0.82);
  zrap_Lumi(0.80,0.90,36);

  TLegend* tl=new TLegend(0.17,0.25,0.86,0.40);
  TLegend* tl2=new TLegend(0.21,0.68,0.90,0.84);
  tl->SetFillStyle(0);
  tl2->SetFillStyle(0);
  //  tl->AddEntry(combo_all,"Z#rightarrow e^{+}e^{-} and Z#rightarrow #mu^{+}#mu^{-} combined","P");
  tl->AddEntry(combo_all,"e^{+}e^{-} and #mu^{+}#mu^{-} combined","P");
  tl->AddEntry(truth,"POWHEG + PYTHIA(Z2) + CT10","L");
  tl2->AddEntry(combo_all,"e^{+}e^{-} and #mu^{+}#mu^{-} combined","P");
  tl2->AddEntry(truth,"POWHEG + PYTHIA(Z2) + CT10","L");

  tl2->Draw();

  c1->SetLogx();

  zrapPrint(c1,"dsdqt_combined_linear");
  
  c1->SetLogy();

  delete tl2;

  tl->Draw();


  zrapPrint(c1,"dsdqt_combined_log");

  int ndof=-1;
  double chi2=0;
  
  FILE* ftable=fopen("dsdqt_combined.tex","wt");
  fprintf(ftable,"\\begin{tabular}{|c|c||l|l|}\n\\hline\n");
  fprintf(ftable,"$p_{\\mathrm{T,min}}$ & $p_{\\mathrm{T,max}}$ & Measurement & Uncertainty  \\\\ \\hline\n");
  for (int i=1; i<=18; i++) {
    double bl=combo_all->GetXaxis()->GetBinLowEdge(i);
    double bh=combo_all->GetXaxis()->GetBinUpEdge(i);
    if (i==1) bl=0;
    int prec=int(-log10(combo_all->GetBinError(i)))+2;
    fprintf(ftable,"%7.1f & %7.1f & %7.*f & %7.*f \\\\ \n",
	    bl,bh,
	    prec,combo_all->GetBinContent(i),
	    prec,combo_all->GetBinError(i));
    ndof++;

    double toadd=pow(combo_all->GetBinContent(i)-truth->GetBinContent(i),2)/pow(combo_all->GetBinError(i),2);
    chi2+=toadd;
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);

  ftable=fopen("dsdqt_combined_all.tex","wt");
  fprintf(ftable,"\\begin{tabular}{|c|c||ll|ll|ll|}\n\\hline\n");
  fprintf(ftable,"$p_{\\mathrm{T,min}}$ & $p_{\\mathrm{T,max}}$ & \\multicolumn{2}{c|}{Electron Channel} & \\multicolumn{2}{c|}{Muon Channel} & \\multicolumn{2}{c|}{Combination} \\\\ \\hline\n");
  for (int i=1; i<=18; i++) {
    double bl=combo_all->GetXaxis()->GetBinLowEdge(i);
    double bh=combo_all->GetXaxis()->GetBinUpEdge(i);
    if (i==1) bl=0;
    int prec=int(-log10(combo_all->GetBinError(i)))+2;
    fprintf(ftable,"%7.1f & %7.1f ",bl,bh);

    fprintf(ftable,"& %7.*f & $\\pm$%7.*f ",
	    prec,elec_all->GetBinContent(i),
	    prec,elec_all->GetBinError(i));

    fprintf(ftable,"& %7.*f & $\\pm$%7.*f ",
	    prec,mu_all->GetBinContent(i),
	    prec,mu_all->GetBinError(i));

    fprintf(ftable,"& %7.*f & $\\pm$%7.*f ",
	    prec,combo_all->GetBinContent(i),
	    prec,combo_all->GetBinError(i));

    fprintf(ftable,"\\\\ \n");
  }
  fprintf(ftable,"\\hline\n\\end{tabular}\n");
  fclose(ftable);

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
  zrapPrint(c3,"dsdqt_combined_deltas");

  printf(" Comparison : %f/%d %f \n",chi2,ndof,TMath::Prob(chi2,ndof));

}
