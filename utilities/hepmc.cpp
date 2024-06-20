/* 
 * File for writing out HepMC event information
 *
 * Author: Michael Rauch <rauch@particle.uni-karlsruhe.de>
 * adapted from the LHEF -> HEPMC converter by Andy Buckley, Les Houches 2011
 * Last modified: Jun 2011
 */

#include <iostream>
#include <sstream>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

extern "C" {
  void hepmcheader_(char*,int);
  void hepmcheader2_(char*,int);
  void hepmcheader3_(int*,int*,
                     double*,double*,
                     int*,int*,int*,int*,int*,int*,
                     double*,double*,double*,
                     int*);
  void hepmcevent_(double*,double*,double*,double*);
  void hepmceventp_(int*,int*,
                     double*,double*,double*,
                     double*,double*);
  void hepmcevente_();
  void hepmcend_();
}

using namespace std;
using namespace HepMC;

IO_GenEvent* writer = NULL;
stringstream ss(stringstream::in|stringstream::out);
GenCrossSection cs;
GenEvent* evt;
GenVertex* v;
int idb1, idb2;
double ebm1, ebm2;

string mytrim(const char* c ,const int l) {
  string s(c,l);
  size_t t = s.find_last_not_of(" ");
  s = s.substr(0,t+1);
  return s;
}

void hepmcheader_(char* file ,int l) {
  writer = new IO_GenEvent(mytrim(file,l));
}

void hepmcheader2_(char* line ,int l) {
  ss << mytrim(line,l) << endl;
}

void hepmcheader3_(int* idbmup1,int* idbmup2,
                   double* ebmup1,double* ebmup2,
                   int* pdfgup1,int* pdfgup2,int* pdfsup1,int* pdfsup2,int* idwtup,int* nprup,
                   double* xsecup,double* xerrup,double* xmaxup,
                   int* lprup) {
  if (writer) {
    writer->write_comment(ss.str());
    cs.set_cross_section(*xsecup,*xerrup);
    idb1 = *idbmup1;
    idb2 = *idbmup2;
    ebm1 = *ebmup1;
    ebm2 = *ebmup2;
  }
}

void hepmcevent_(double* scalup,double* aqcdup,double* aqedup,double *xwgtup) {
  evt = new GenEvent();
  evt->use_units(Units::GEV, Units::MM);
  v = new GenVertex();
  evt->add_vertex(v);
  FourVector beam1(0, 0, ebm1, ebm1);
  GenParticle* gp1 = new GenParticle(beam1,idb1, 4);
  v->add_particle_in(gp1);
  FourVector beam2(0, 0, ebm2, ebm2);
  GenParticle* gp2 = new GenParticle(beam2,idb2, 4);
  v->add_particle_in(gp2);
  evt->set_beam_particles(gp1, gp2);
  evt->set_event_scale(*scalup);
  evt->set_alphaQCD(*aqcdup);
  evt->set_alphaQED(*aqedup);
  evt->weights().push_back(*xwgtup);
}

void hepmceventp_(int* idup,int* istup,
                  double* p1,double* p2,double* p3,
                  double* p4,double* m) {
  FourVector p(*p1, *p2, *p3, *p4);
  GenParticle* gp = new GenParticle(p, *idup, *istup);
  gp->set_generated_mass(*m);
  v->add_particle_out(gp);
}

void hepmcevente_() {
  writer->write_event(evt);
  delete evt;
}

void hepmcend_() {
  delete writer;
}



