TH1* readStandardFile(const char* name, const char* fname, int ec=0, bool isAbs=false) {
  TH1* dist=new TH1F(name,name,100,-5.0,5.0);
  dist->SetDirectory(0);
  dist->Sumw2();

  FILE* af=fopen(fname,"rt");
  char buffer[1024];

  int irow;
  float flow,fhi,val,err[4];
  while (!feof(af)) {
    buffer[0]=0;
    fgets(buffer,1023,af);
    if (strchr(buffer,'#')!=0) *(strchr(buffer,'#'))=0;
    
    int matched=sscanf(buffer,"%d %f %f %f %f %f %f",&irow,&flow,&fhi,&val,&(err[0]),&(err[1]),&(err[2]));

    if (isAbs) irow+=50;

    dist->SetBinContent(irow,val);
    dist->SetBinError(irow,err[ec]);
  }
  fclose(af);

  return dist;
}

TH1* readStandardFilePt(const char* name, const char* fname, int ec=0) {
  static const int    pt_bins=19;
  static const double pt_binning[]={ 0.0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
  TH1* dist=new TH1F(name,name,pt_bins,pt_binning);
  dist->SetDirectory(0);
  dist->Sumw2();

  FILE* af=fopen(fname,"rt");
  if (af==0) {
    fprintf(stderr,"Unable to open '%s'\n",fname);
    return 0;
  }
  char buffer[1024];

  int irow;
  float flow,fhi,val,err[5];
  while (!feof(af)) {
    buffer[0]=0;
    fgets(buffer,1023,af);
    if (strchr(buffer,'#')!=0) *(strchr(buffer,'#'))=0;

    int matched=sscanf(buffer,"%d %f %f %f %f %f %f",&irow,&flow,&fhi,&val,&(err[0]),&(err[1]),&(err[2]));

    dist->SetBinContent(irow,val);
    dist->SetBinError(irow,err[ec]);
  }
  fclose(af);

  return dist;
}


TGraphAsymmErrors* readStandardFilePtTGE(const char* fname, int skip=0) {
  static const int    pt_bins=18;
  static const double pt_binning[]={ 0.0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600, 3500 };
  TGraphAsymmErrors* dist=new TGraphAsymmErrors(pt_bins-skip);
  //  dist->SetDirectory(0);

  FILE* af=fopen(fname,"rt");
  if (af==0) return 0;
  char buffer[1024];

  int irow;
  float flow,fhi,val,err[5];
  while (!feof(af)) {
    buffer[0]=0;
    fgets(buffer,1023,af);
    if (strchr(buffer,'#')!=0) *(strchr(buffer,'#'))=0;

    int matched=sscanf(buffer,"%d %f %f %f %f %f %f",&irow,&flow,&fhi,&val,&(err[0]),&(err[1]),&(err[2]));

    dist->SetPoint(irow-1-skip,(flow+fhi)/2.0,val);
    dist->SetPointError(irow-1-skip,(fhi-flow)/2.0,(fhi-flow)/2.0,fabs(err[0]),err[1]);
  }
  fclose(af);

  return dist;
}
