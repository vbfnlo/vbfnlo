#include <iostream>
#include <stdio.h>
#include "TString.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TTree.h"
// #include "TChain.h"

extern "C" 
{
  void initrootoutput_(int*, int*, char*, int);
  void closerootoutput_();
  void writerootevent_(int*, double*, int*, double*, double*);
  void writerootevent2_(int*, double*, int*, double*, double*);
  void writerootevent3_(int*, double*, int*, double*, double*, int*);
}

namespace rootevents {
  TString eventfile;
  TTree *eventtree;
  int num_jets;
  int num_leptons;
  int tau_chi[2];
  TClonesArray *jets; 
  TClonesArray *leptons; 
  double evtweight;
  TFile *eventf;
}

using namespace rootevents;

void initrootoutput_(int* nj, int* nl, char* file, int len)
{
  eventfile.Resize(0);
  eventfile.Append(file);
  int siz = eventfile.Index(' ');
  eventfile.Resize(siz);

  eventf = new TFile(eventfile,"RECREATE");  
  eventf->cd();
  eventtree = new TTree("events","vbf tree output");
  //eventchain = new TChain("events","");
  jets = new TClonesArray("TLorentzVector",*nj);
  leptons = new TClonesArray("TLorentzVector",*nl);

  eventtree->Branch("jets","TClonesArray",&jets,8000,0);
  eventtree->Branch("leptons","TClonesArray",&leptons,8000,0);
  eventtree->Branch("weight",&evtweight,"evtweight/D");
  eventtree->Branch("nj",&num_jets,"num_jets/I");
  eventtree->Branch("nl",&num_leptons,"num_leptons/I");
  eventtree->Branch("tauchi",&tau_chi,"tau_chi[2]/I");
  return;
}

void closerootoutput_()
{
  //  eventf->Write();
  gFile->Write();
  //  eventf->Close();
  gFile->Close();
  return;
}

void writerootevent_(int* nj, double* ajets, int* nl, double* aleptons, double* weight)
{
  if (*weight>0.0) {
    eventf->cd();
    num_jets = *nj;
    num_leptons = *nl; 
    for (Int_t i=0;i<num_jets;i++) { new(jets->AddrAt(i)) TLorentzVector();}
    for (Int_t i=0;i<num_leptons;i++) { new(leptons->AddrAt(i)) TLorentzVector();}

    for (Int_t i=0; i<*nj; i++){
      new(jets->AddrAt(i)) TLorentzVector(ajets[8*i+1],ajets[8*i+2],ajets[8*i+3],ajets[8*i+0]);

    }
    for (Int_t i=0; i<*nl; i++){
      new(leptons->AddrAt(i)) TLorentzVector(aleptons[8*i+1],aleptons[8*i+2],aleptons[8*i+3],aleptons[8*i+0]);

    }
    evtweight=*weight;
    eventtree->Fill();
    jets->Delete();
    leptons->Delete();
  }
  return;
}

void writerootevent2_(int* np, double* p, int* nv, double* v, double* weight)
{
  if (*weight>0.0) {
    eventf->cd();
    num_jets = *np;
    num_leptons = *nv; 
    for (Int_t i=0;i<num_jets;i++) { new(jets->AddrAt(i)) TLorentzVector();}
    for (Int_t i=0;i<num_leptons;i++) { new(leptons->AddrAt(i)) TLorentzVector();}

    for (Int_t i=0; i<*np; i++){
      new(jets->AddrAt(i)) TLorentzVector(p[4*i+1],p[4*i+2],p[4*i+3],p[4*i+0]);

    }
    for (Int_t i=0; i<*nv; i++){
      new(leptons->AddrAt(i)) TLorentzVector(v[4*i+1],v[4*i+2],v[4*i+3],v[4*i+0]);

    }
    evtweight=*weight;
    eventtree->Fill();
    jets->Delete();
    leptons->Delete();
  }
  return;
}

void writerootevent3_(int* nj, double* ajets, int* nl, double* aleptons, double* weight, int tauchi[2])
{
  if (*weight>0.0) {
    eventf->cd();
    num_jets = *nj;
    num_leptons = *nl; 
    for (Int_t i=0;i<num_jets;i++) { new(jets->AddrAt(i)) TLorentzVector();}
    for (Int_t i=0;i<num_leptons;i++) { new(leptons->AddrAt(i)) TLorentzVector();}

    for (Int_t i=0; i<*nj; i++){
      new(jets->AddrAt(i)) TLorentzVector(ajets[8*i+1],ajets[8*i+2],ajets[8*i+3],ajets[8*i+0]);

    }
    for (Int_t i=0; i<*nl; i++){
      new(leptons->AddrAt(i)) TLorentzVector(aleptons[8*i+1],aleptons[8*i+2],aleptons[8*i+3],aleptons[8*i+0]);

    }
    evtweight=*weight;
    tau_chi[0] = tauchi[0];
    tau_chi[1] = tauchi[1];
 
    eventtree->Fill();
    jets->Delete();
    leptons->Delete();
  }
  return;
}
