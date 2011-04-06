#include "readStandardFile.C"

void csvToTGAE(const char* csvf) {
  TH1* h1=readStandardFilePt("csv",csvf);

  printf("  TGraphAsymmErrors *grae = new TGraphAsymmErrors(18);\n");
  for (int i=1; i<=18; i++) {
    double bc=h1->GetBinCenter(i);
    double bw=h1->GetBinWidth(i)/2;
    if (i==1) {
      bc=1.25;
      bw=1.25;
    }
    printf("  grae->SetPoint(%d,%f,%e);\n",i-1,bc,h1->GetBinContent(i));
    printf("  grae->SetPointError(%d,%f,%f,%e,%e);\n",i-1,bw,bw,h1->GetBinError(i),h1->GetBinError(i));
  }

}
