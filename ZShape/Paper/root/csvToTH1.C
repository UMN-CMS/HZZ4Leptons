#include "readStandardFile.C"

void csvToTH1(const char* csvf) {
  TH1* h1=readStandardFilePt("csv",csvf);

  printf("  TH1 *myh1 = new TH1F(\"SomeSuch\",\"\",18, xAxis1);\n");
  for (int i=1; i<=18; i++) {
    printf("  myh1->SetBinContent(%d,%e);\n",i,h1->GetBinContent(i));
    printf("  myh1->SetBinError(%d,%e);\n",i,h1->GetBinError(i));
  }

}
