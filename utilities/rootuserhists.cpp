/*
 *  Book and fill ROOT histograms in VBFNLO
 *  =======================================
 *
 * - VBFNLO has to be compiled with "with-root" flag
 * - Book histos in initroothists_()
 * - fillroothists_() is called for each event. Fill histos here
 * - final processing can be done in closeroothists_()
 * - histos are written to the root file defined in vbfnlo.dat
 *   together with the histos defined in the FORTRAN code
 * - the passed weight is dsig*weight in routine histograms.F
 * - concerning error estimation for NLO distributions:
 *   at NLO there are several histogram entries for each
 *   phase space point which are correlated. Therefore,
 *   the automatically calculated error estimate from ROOT
 *   is not valid, but instead vastly overestimates the error.
 *   For correct error estimates the histograms defined in
 *   the FORTRAN file histograms.F (which also have a ROOT output)
 *   should be used. They take the correlation into account and 
 *   provide valid error estimates.
 *
 * 
 * File contributed by: Christian Zeitnitz, 05 Feb 2011
 * Last modification: Michael Rauch, 18 Feb 2011
 */

#include "roothist.h"

#define DefineHist(t,x) t *x ## _LO, *x ## _NLO;
#define CreateHist(t,x,name,...) { \
  string fnameext, fname = name; \
  fnameext = fname+"_LO"; \
  x ## _LO  = new t(fnameext.c_str(),__VA_ARGS__); \
  fnameext = fname+"_NLO"; \
  x ## _NLO = new t(fnameext.c_str(),__VA_ARGS__); \
  Hlist->Add(x ## _LO ); \
  Hlist->Add(x ## _NLO ); \
}
#define FillHist(x,...) { \
  if (*nlo==0) { \
    x ## _LO ->Fill(__VA_ARGS__); \
  } else if (*nlo!=0) { \
    x ## _NLO ->Fill(__VA_ARGS__); \
  } \
}

/*  
 * Whats stored in jets[7] and leptons[8]:
 * 
 * index    content
 *   0        E
 *   1        Px
 *   2        Py
 *   3        Pz
 *   4        mass
 *   5        Pt 
 *   6        eta
 *   7        phi
 *   8        lepton type (11-13)
 */

namespace roothists {

  vector<TLorentzVector> jet;
  vector<TLorentzVector> lepton;
  vector<TLorentzVector> invisible;
  vector<TLorentzVector> photon;
  
  DefineHist(TH1F,h_ptj);
  DefineHist(TH1F,h_ptmaxl);
  DefineHist(TH2F,h_etajj_mjj);
}
using namespace roothists;


// Initialization of ROOT histograms
// =================================
// - define a pointer to the histo in the above namespace
// - called by InitHistograms() (histograms.F)
// 
void initroothists_() {  
  hfile->cd();
  
  CreateHist(TH1F,h_ptj,"H_ptj","d#sigma/d{p_T}_j (fb/GeV)",100,0.,250.);
  CreateHist(TH1F,h_ptmaxl,"H_ptmaxl" ,"d#sigma/d{p_T,max}_l (fb/GeV)",100,0.,500.);
  CreateHist(TH2F,h_etajj_mjj,"H2_etajj_mjj" ,"d^2#sigma/d#eta_{jj} dm_{jj} (fb/GeV)", 25,0.,6.,25,0.,800.);

}


// Fill ROOT histograms 
// ====================
// - called for each event (several times in case of the real emission part of NLO calculations)
// - called by HistogramEvent() (histograms.F)
//
void fillroothists_(double *weight,double *jets,int *nj,double *leptons,int *nl,double* invisibles,int *ninv,double *photons,int *nph,int *nlo) { 

// store input momenta into LorentzVectors
  jet.resize(*nj);
  lepton.resize(*nl);
  invisible.resize(*ninv);
  photon.resize(*nph);

// copy jets, leptons, invisibles, photons to TLorentzVectors 
  for (Int_t i=0; i < *nj; i++) {
  	jet[i].SetPxPyPzE(jets[8*i+1],jets[8*i+2],jets[8*i+3],jets[8*i+0]);
  }
  for (Int_t i=0; i < *nl; i++) {
  	lepton[i].SetPxPyPzE(leptons[9*i+1],leptons[9*i+2],leptons[9*i+3],leptons[9*i+0]);
  }
  for (Int_t i=0; i < *ninv; i++) {
  	invisible[i].SetPxPyPzE(invisibles[9*i+1],invisibles[9*i+2],invisibles[9*i+3],invisibles[9*i+0]);
  }
  for (Int_t i=0; i < *nph; i++) {
  	photon[i].SetPxPyPzE(photons[8*i+1],photons[8*i+2],photons[8*i+3],photons[8*i+0]);
  }

// Fill histograms

// pt_j
  for (Int_t i=0; i < *nj; i++) {
    FillHist(h_ptj,jet[i].Pt(),*weight/(static_cast<double>(*nj)));
  }

// ptmax_l
  double ptmax_l = 0.;
  for (Int_t i=0; i < *nl; i++) {
    if (lepton[i].Pt() > ptmax_l) ptmax_l = lepton[i].Pt();
  }
  FillHist(h_ptmaxl,ptmax_l,*weight);

// h_etajj_mjj
  if (*nj >= 2) {
    FillHist(h_etajj_mjj,fabs(jet[0].Rapidity()-jet[1].Rapidity()),
                         (jet[0]+jet[1]).M(),*weight);
  }

}

// Close ROOT Histos
// =================
// - endprocessing of root histos
// - file is written in roothist.cpp, so don't do it here !!
// - called by WriteHistograms() (histograms.F)
//

void closeroothists_() {

// correct dsig/dX for bin width of histogram
//	cout << endl << "Correct ROOT histograms for bin width " << endl;
// correct weight 
  TIter Hiter(Hlist);
  TObject* hnext;

  TH1F *h1d;
  TH2F *h2d;
  while( hnext = Hiter.Next() ) { 
    if (h2d = dynamic_cast<TH2F *>(hnext) ) {
      double dx = h2d->GetXaxis()->GetBinWidth(0);
      double dy = h2d->GetYaxis()->GetBinWidth(0);
      if(dx*dy>0.) h2d->Scale(1./(dx*dy));
    } else if (h1d = dynamic_cast<TH1F *>(hnext) ) {
      double dx = h1d->GetBinWidth(0);
      if(dx >0.) h1d->Scale(1./dx);
    } else {
      std::cerr << "Internal error in rootuserhists.cpp: " 
                << "Unconvertible object in Hlist." << std::endl;
    }
  }	

}

