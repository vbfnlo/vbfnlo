#include "roothist.h"

namespace roothists {
  int max_jets;
  TObjArray* Hlist;
  TString rootfile;
  TFile *hfile;
  vector<TH1F*> lth1f;
  vector<TH2F*> lth2f;
}

void openrootsession_(char* file, int len)
{
  
  rootfile.Resize(0);
  rootfile.Append(file);
  int siz = rootfile.Index(' ');
  rootfile.Resize(siz);
  lth1f.clear();
  lth2f.clear();

  Hlist = new TObjArray(0);
  hfile = new TFile(rootfile,"recreate");
  hfile->cd();
  return;
}

void closerootsession_()
{
  hfile->cd();
  if (Hlist->GetEntries() > 0 )  {
    Hlist->Write();
    hfile->Close();
  }
  return;
}

void createroothisto_(int* TID, int* ID, char* title, int* bins, double* min, double* max, char* type, int len, int len2)
{
  char name[20];
  TString htitle(title,len);
  htitle = htitle.Strip();
  TString htype(type,len2);
  htype = htype.Strip();
  sprintf( name, "H_%d_%s", *TID, htype.Data() );
  hfile->cd();
  TH1F* histo = new TH1F(name, htitle, *bins, *min, *max);
  Hlist->Add(histo);
  lth1f.resize(std::max(static_cast<int>(lth1f.size()),*ID));
  lth1f[*ID-1]=histo;
  return;
}

void createroothisto2d_(int* TID, int* ID, char* title, int* binx, double* minx, double* maxx, int* biny, double* miny, double* maxy, char* type, int len, int len2)
{
  char name[20];
  TString htitle(title,len);
  htitle = htitle.Strip();
  TString htype(type,len2);
  htype = htype.Strip();
  sprintf( name, "H2d_%d_%s", *TID, htype.Data() );
  hfile->cd();
  TH2F* histo = new TH2F(name, htitle, *binx, *minx, *maxx, *biny, *miny, *maxy);
  Hlist->Add(histo);
  lth2f.resize(std::max(static_cast<int>(lth2f.size()),*ID));
  lth2f[*ID-1]=histo;
  return;
}

void fillroothisto_(int* ID, double* value, double* weight, double* error)
{
  int binnumber;
  hfile->cd();
  TH1F* histo = lth1f[*ID-1];
  if (histo != 0) {
    histo->Fill(*value, *weight);
    binnumber = histo->FindBin(*value);
    histo->SetBinError(binnumber, *error);
  }
}

void fillroothisto2d_(int* ID, double* valuex, double* valuey, double* weight, double* error)
{
  int binnumber;
  hfile->cd();
  TH2F* histo = lth2f[*ID-1];
  if (histo != 0) {
    histo->Fill(*valuex, *valuey, *weight);
    binnumber = histo->FindBin(*valuex, *valuey);
    histo->SetBinError(binnumber, *error);
  }
}



