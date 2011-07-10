#include "readStandardFile.C"
//#include "readUnfoldingMatrices.C"

TH1* extractDistribution(const char* name, TFile* f) {

    TH1* effAcc = readStandardFilePt("effAcc", "../data/eff_pt.txt");

    const char* h1name = "mcEff/ECAL80-ECAL95/C08-m(60,120)/Z0_Pt";

    TH1* sig1 = f->Get(h1name)->Clone(name);
    sig1->Sumw2();

    sig1->Divide(effAcc);

    delete effAcc;
    return sig1;
}

#include "tdrstyle.C"
#include "zrapidityStandard.C"

void energyScalePt(TFile* base, TFile* vpe, TFile* vme, TFile* tpe, TFile* tme, const char* errfilename = 0) {

    FILE* errfile = 0;
    if (errfilename != 0) {
        errfile = fopen(errfilename, "wt");
        TDatime now;
        fprintf(errfile, "# Energy Scale Error # %s\n", now.AsString());
        fprintf(errfile, "#Bin Y_Min Y_Max Value Uncertainty\n");
    }

    gROOT->SetStyle("Plain");
    setTDRStyle();

    TH1* hb = extractDistribution("base", base);
    TH1* valep = extractDistribution("vpe", vpe);
    TH1* valem = extractDistribution("vme", vme);
    TH1* valtep = extractDistribution("tpe", tpe);
    TH1* valtem = extractDistribution("tme", tme);

    valep->Scale(hb->Integral(1, 12) / valep->Integral(1, 12));
    valem->Scale(hb->Integral(1, 12) / valem->Integral(1, 12));

    valtep->Scale(hb->Integral(1, 12) / valtep->Integral(1, 12));
    valtem->Scale(hb->Integral(1, 12) / valtem->Integral(1, 12));

    //hb = unfold(hb, "unfoldingMatrix_theOutPut.root");
    //valep = unfold(valep, "unfoldingMatrix_theOutPut.root");
    //valem = unfold(valem, "unfoldingMatrix_theOutPut.root");

    valep->Add(hb, -1);
    valem->Add(hb, -1);

    valtep->Add(hb, -1);
    valtem->Add(hb, -1);

    valep->Divide(hb);
    valem->Divide(hb);

    valtep->Divide(hb);
    valtem->Divide(hb);

    TH1* valm = valem->Clone("valm");
    TH1* valp = valep->Clone("valp");
    TH1* valtm = valtem->Clone("valtm");
    TH1* valtp = valtep->Clone("valtp");
    TH1* valave = valep->Clone("valave");

    for (int i = 1; i <= valep->GetXaxis()->GetNbins(); i++) {
        valm->SetBinContent(i, sqrt(pow(valem->GetBinContent(i),2)+pow(valtem->GetBinContent(i),2)));
        valp->SetBinContent(i, sqrt(pow(valem->GetBinContent(i),2)+pow(valtem->GetBinContent(i),2)));
        valave->SetBinContent(i, 0.5 * (
                valm->GetBinContent(i) +
                valp->GetBinContent(i)
                ));
        valave->SetBinError(i, 0.5 * sqrt(
                pow(valm->GetBinError(i), 2) +
                pow(valp->GetBinError(i), 2)
                ));

	//        if (i >= 15 && i <= 86) {
            if (errfile != 0)
                fprintf(errfile, "%4d %5.2f %5.2f %6.4f %6.4f\n", i,
                    valave->GetXaxis()->GetBinLowEdge(i),
                    valave->GetXaxis()->GetBinUpEdge(i),
                    valave->GetBinContent(i),
                    valave->GetBinError(i)
                    );
	    //        }

    }

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetLeftMargin(0.17);
    c1->SetTopMargin(0.03);
    c1->SetRightMargin(0.02);
    c1->SetBottomMargin(0.15);

    //  valm->SetLineColor(kBlue);
    //  valp->SetLineColor(kRed);

    //  valm->SetLineWidth(2);
    //  valp->SetLineWidth(2);


    valave=zpt_rebinForPlot(valave);

    c1->SetLogx();

    valave->SetLineWidth(2);
    valave->SetMaximum(0.1);
    valave->SetMinimum(0.0);

    valave->SetTitle(0);
    valave->SetStats(0);
    valave->GetXaxis()->SetTitle("p_{T,ee} [GeV/c]");
    valave->GetXaxis()->SetTitleOffset(1.2);
    valave->GetXaxis()->SetRangeUser(0.7, 500);
    valave->GetXaxis()->CenterTitle(true);
    valave->GetYaxis()->SetTitle("Fractional Error from Energy Scale");
    valave->GetYaxis()->SetTitleOffset(1.4);
    valave->GetYaxis()->CenterTitle(true);
    valave->Draw("HIST");
    //  valm->Draw("HIST SAME");

    zrap_Prelim(0.6, 0.9, 0.6, 0.87);
    c1->Print("energyScalePt.eps");

 const int pxbins=19;
  const double pxbinEdges[]={0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190,250, 600,3500};

  for (int kk=1;kk<=19;kk++){
    printf("%9d %9.1f %9.1f %9.4f %9.4f \n",kk,pxbinEdges[kk-1],pxbinEdges[kk],valave->GetBinContent(kk),sqrt(valave->GetBinContent(kk)));

  }


    if (errfile != 0) fclose(errfile);

}

