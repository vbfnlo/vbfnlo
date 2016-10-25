/*
 * Public Interface for the Binoth Les Houches Accord
 *
 * Author: Michael Rauch <michael.rauch@kit.edu>
 * Initial version: Nov 2013
 * Last modified: Mar 2015
 *
 */

#include <vector>
#include <string>
#include "VBFNLOConfig.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VBFNLO-internal Fortran routines
 */
void FC_FUNC_(olp_setparameter_vbfnlo,OLP_SETPARAMETER_VBFNLO)(char* line, double* re, double* im, int* ierr, int nline);
void FC_FUNC_(olp_getparameter_vbfnlo,OLP_GETPARAMETER_VBFNLO)(char* line, double* re, double* im, int* ierr, int nline);
void FC_FUNC_(olp_evalsubprocess_vbfnlo,OLP_EVALSUBPROCESS_VBFNLO)(int* i, double* pp, double* mu, double* alphas, double* rval);
void FC_FUNC_(olp_evalsubprocess2_vbfnlo,OLP_EVALSUBPROCESS2_VBFNLO)(int* i, double* pp, double* mu, double* rval, double* acc);
void FC_FUNC_(vbfnlo_setupprocess,VBFNLO_SETUPPROCESS)(int* nparticles,int* pdgprocess,int* orderAlphas,int* orderAlpha,int* AmplitudeType,int* procok);
void FC_FUNC_(olp_phasespacepoint,OLP_PHASESPACEPOINT)(int* proc, double* rpsnum, double* r, double* p, double* weight);
void FC_FUNC_(olp_polvec,OLP_POLVEC)(double* p, double* q, double* eps);
void FC_FUNC_(olp_info,OLP_INFO)(char* version, char* message, int nversion=15, int nmessage=255);
void FC_FUNC_(blha_initialize,BLHA_INITIALIZE)();
void FC_FUNC_(blha_start,BLHA_START)();
void FC_FUNC_(blha_setdr,BLHA_SETDR)(int* dr);
void FC_FUNC_(blha_setnf,BLHA_SETNF)(int* nf);
void FC_FUNC_(blha_setnc,BLHA_SETNC)(int* nc);
void FC_FUNC_(blha_setcouplingsoff,BLHA_SETCOUPLINGSOFF)();
void FC_FUNC_(blha_setewrenormscheme,BLHA_SETEWRENORMSCHEME)(int* scheme);

#ifdef __cplusplus
}
#endif

/* C++-internal functions follow here */

#ifdef __cplusplus
namespace VBFNLO {
#endif

class BLHA {

public:
  static BLHA& instance() {
     static BLHA _instance;
     return _instance;
  }
  ~BLHA() {};
  void ParseOrderLine(std::string& line, std::string& result, int* ierr);
  void DetProcess(std::vector<std::string> tokens, std::string& result, int* ierr);
  bool stringtoint(const std::string str, int& i);
  bool stringtodouble(const std::string str, double& d);
  std::string stringlower(const std::string str);
  bool isparton(int pdg);

protected:
// dis-allow constructor and copy constructor
  BLHA() {
    BLHAtype=1;
    orderAlpha=0;
    orderAlphas=0;
    AmplitudeType=-1; // -1=undet
  };           
  BLHA( const BLHA& ); 
  BLHA & operator = (const BLHA &); 
  int BLHAtype;
  int orderAlpha;
  int orderAlphas;
  int AmplitudeType; // -1=undet, 0=tree, 1=loop, 2=cctree, 3=sctree

};

#ifdef __cplusplus
}
#endif

