/*
 * Public Interface for the Binoth Les Houches Accord 
 *
 * Author: Michael Rauch <michael.rauch@kit.edu>
 * Initial version: Nov 2013
 * Last modified: Dec 2013
 *
 * at the moment: BLHA2
 *
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "BLHAinterface.h"
#include "BLHAclass.h"

using namespace std;

#ifdef WITH_CUSTOM_BLHA
void OLP_Order_VBFNLO(char* inname, char* outname, int* ierr)
#else 
void OLP_Order(char* inname, char* outname, int* ierr)
#endif
{

  *ierr = 1;
 
  ifstream orderfile(inname);
  if (! orderfile.is_open()) {
    cerr << "VBFNLO OLP_Order: Cannot open order file " << inname << endl;
    *ierr = -1;
    return;
  }
  orderfile.exceptions(ifstream::badbit);

  ofstream contractfile(outname);
  if (! contractfile.is_open()) {
    cerr << "VBFNLO OLP_Order: Cannot open contract file " << outname << endl;
    *ierr = -1;
    return;
  }
  
  VBFNLO::BLHA& blha = VBFNLO::BLHA::instance();

  FC_FUNC_(blha_initialize,BLHA_INITIALIZE)();

  string line;

  try {
    while (! orderfile.eof() ) {
      getline(orderfile,line);
      contractfile << line;
      // handle line
      int lineerr = 1;
      string result;
      blha.ParseOrderLine(line,result,&lineerr);
      if (result!="") contractfile << " | " << result;
      if (lineerr!=1) *ierr = -2; 
      contractfile << endl;
    }
  } catch (ifstream::failure e) {
    cerr << "VBFNLO OLP_Order: Error reading order file " << inname << endl;
    *ierr = -1;
    return;
  }
  
  return;
}

#ifdef WITH_CUSTOM_BLHA
void OLP_Start_VBFNLO(char* fname, int* ierr)
#else
void OLP_Start(char* fname, int* ierr)
#endif
{

  *ierr = 1;
 
  ifstream contractfile(fname);
  if (! contractfile.is_open()) {
    cerr << "VBFNLO OLP_Start: Cannot open contract file " << fname << endl;
    *ierr = -1;
    return;
  }
  contractfile.exceptions(ifstream::badbit);
  
  VBFNLO::BLHA& blha = VBFNLO::BLHA::instance();

  FC_FUNC_(blha_initialize,BLHA_INITIALIZE)();
  FC_FUNC_(blha_start,BLHA_START)();

  string line;

  try {
    while (! contractfile.eof() ) {
      getline(contractfile,line);
      // handle line
      int lineerr = 1;
      string result;
      blha.ParseOrderLine(line,result,&lineerr);
      if (lineerr!=1) {
        cerr << "VBFNLO OLP_Start: Error in contract file " << fname << " :" << endl;
        cerr << "  " << line << endl;
        cerr << " -> " << result << endl;
        *ierr = -1;
      }
    }
  } catch (ifstream::failure e) {
    cerr << "VBFNLO OLP_Start: Error reading contract file " << fname << endl;
    *ierr = -1;
    return;
  }
  
  return;
}

#ifdef WITH_CUSTOM_BLHA
void OLP_Info_VBFNLO(char olp_name[15], char olp_version[15], char message[255])
#else 
void OLP_Info(char olp_name[15], char olp_version[15], char message[255])
#endif
 {
  strncpy(olp_name,"VBFNLO",15); 
  FC_FUNC_(olp_info,OLP_INFO)(olp_version,message);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_GetParameter_VBFNLO(char* line, double* re, double* im, int* ierr)
#else
void OLP_GetParameter(char* line, double* re, double* im, int* ierr)
#endif 
{
  FC_FUNC_(olp_getparameter_vbfnlo,OLP_GETPARAMETER_VBFNLO)(line, re, im, ierr, strlen(line));
  *ierr = abs(*ierr);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_SetParameter_VBFNLO(char* line, double* re, double* im, int* ierr)
#else
void OLP_SetParameter(char* line, double* re, double* im, int* ierr)
#endif
{
  FC_FUNC_(olp_setparameter_vbfnlo,OLP_SETPARAMETER_VBFNLO)(line, re, im, ierr, strlen(line));
  *ierr = abs(*ierr);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_EvalSubProcess_VBFNLO(int i, double* pp, double mu, double* alphas, double* rval)
#else
void OLP_EvalSubProcess(int i, double* pp, double mu, double* alphas, double* rval)
#endif
{
  FC_FUNC_(olp_evalsubprocess_vbfnlo,OLP_EVALSUBPROCESS_VBFNLO)(&i, pp, &mu, alphas, rval);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_EvalSubProcess2_VBFNLO(int* i, double* pp, double* mu, double* rval, double* acc)
#else
void OLP_EvalSubProcess2(int* i, double* pp, double* mu, double* rval, double* acc) 
#endif
{
  FC_FUNC_(olp_evalsubprocess2_vbfnlo,OLP_EVALSUBPROCESS2_VBFNLO)(i, pp, mu, rval, acc);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_PhaseSpacePoint_VBFNLO(int* proc, double* rpsnum, double* r, double* p, double* weight)
#else
void OLP_PhaseSpacePoint(int* proc, double* rpsnum, double* r, double* p, double* weight)
#endif
{
  FC_FUNC_(olp_phasespacepoint,OLP_PHASESPACEPOINT)(proc, rpsnum, r, p, weight);
}

#ifdef WITH_CUSTOM_BLHA
void OLP_Polvec_VBFNLO(double* p, double* q, double* eps)
#else
void OLP_Polvec(double* p, double* q, double* eps)
#endif
{
  FC_FUNC_(olp_polvec,OLP_POLVEC)(p, q, eps);
}