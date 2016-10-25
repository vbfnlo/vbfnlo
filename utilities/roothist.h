#include <iostream>
#include "stdio.h"
#include "TH1F.h" 
#include "TFile.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TH2F.h"
#include "TMatrix.h"
#include "TString.h"
#include "TRegexp.h"

using namespace std;

extern "C" 
{
  void openrootsession_(char*, int);
  void closerootsession_();
  void initroothists_(void);
  void createroothisto_(int*, int*, char*, int*, double*, double*, char*, int, int);
  void createroothisto2d_(int*, int*, char*, int*, double*, double*, int*, double*, double*, char*, int, int);
  void fillroothisto_(int*, double*, double*, double*);
  void fillroothisto2d_(int*, double*, double*, double*, double*);
  void fillroothists_(double*,double*,int*,double*,int*,double*,int*,double*,int*,int*);
  void closeroothists_();
}

namespace roothists {
  extern int max_jets;
  extern TObjArray* Hlist;
  extern TString rootfile;
  extern TFile *hfile;
  extern vector<TH1F*> lth1f;
  extern vector<TH2F*> lth2f;
}

using namespace roothists;

