/*
 * Public Interface for the Binoth Les Houches Accord
 *
 * Author: Michael Rauch <michael.rauch@kit.edu>
 * Initial version: Nov 2013
 * Last modified: Dec 2013
 *
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "BLHAinterface.h"
#include "BLHAclass.h"

using namespace std;

namespace VBFNLO {

void BLHA::ParseOrderLine(string& line, string& result, int* ierr) {

  string whitespace = " \t\r\n";
  size_t pos;

  // remove comments and own output
  pos = line.find_first_of("#|");
  if (pos < string::npos) line.erase(pos);
  // ignore empty lines and remove leading white space
  pos = line.find_first_not_of(whitespace);
  if (pos == string::npos) return;
  string linechomp = line.substr(pos,string::npos);
  // tokenize
  vector<string> tokens;
  while (linechomp != "") {
    pos = linechomp.find_first_of(whitespace);
    tokens.push_back(linechomp.substr(0,pos));
    pos = linechomp.find_first_not_of(whitespace,pos);
    if (pos == string::npos) break;
    linechomp.erase(0,pos); 
  }
  // now parse all statements
  if (tokens.size()<=1) {
    result = "Error: check syntax - value missing";
    *ierr = -1;
    return;
  }
  
  result = "OK";
  int pdgId;

  if (stringlower(tokens[0]) == "interfaceversion") {
    if (stringlower(tokens[1]) != "blha2") {
      result = "Error: unsupported flag\n# BLHA2 is supported";
      *ierr = -1;
    }
    BLHAtype=2;
  }
  else if (stringlower(tokens[0]) == "model") {
    if (stringlower(tokens[1]) != "sm") {
      result = "Error: unsupported flag\n# SM is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "matrixelementsquaretype") {
    if (stringlower(tokens[1]) != "chsummed") {
      result = "Error: unsupported flag\n# CHSummed is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "correctiontype") {
    if (stringlower(tokens[1]) != "qcd") {
      result = "Error: unsupported flag\n# QCD is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "irregularisation") {
    if (stringlower(tokens[1]) == "cdr") {
      int dr=0;
      FC_FUNC_(blha_setdr,BLHA_SETDR)(&dr);
    } else if (stringlower(tokens[1]) == "dred") {
      int dr=1;
      FC_FUNC_(blha_setdr,BLHA_SETDR)(&dr);
    } else {
      result = "Error: unsupported flag\n# CDR or DRED is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "irsubtractionmethod") {
    if (stringlower(tokens[1]) != "none") {
      result = "Error: unsupported flag\n# None is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "extra") {
    if (tokens.size() < 3) {
      result = "Error: check syntax\n# number of arguments too low";
      *ierr = -1;
    }
    if (stringlower(tokens[1]) == "helavginitial") {
      if (stringlower(tokens[2]) != "no") {
        result = "Error: unsupported flag\n# no is supported";
        *ierr = -1;
      }
    }
    else if (stringlower(tokens[1]) == "colavginitial") {
      if (stringlower(tokens[2]) != "no") {
        result = "Error: unsupported flag\n# no is supported";
        *ierr = -1;
      }
    }
    else if (stringlower(tokens[1]) == "mcsymmetrizefinal") {
      if (stringlower(tokens[2]) != "no") {
        result = "Error: unsupported flag\n# no is supported";
        *ierr = -1;
      }
    }
    else if (stringlower(tokens[1]) == "nf") {
      int Nf;
      if (stringtoint(tokens[2], Nf) && (Nf>=0) ) {
        FC_FUNC_(blha_setnf,BLHA_SETNF)(&Nf);
      } else {
        result = "Error: unsupported flag\n# any non-negative integer is supported";
        *ierr = -1;
      }
    }
    else if (stringlower(tokens[1]) == "nc") {
      int Nc;
      if (stringlower(tokens[2]) == "infinity" ) {
        Nc=-1;
        FC_FUNC_(blha_setnc,BLHA_SETNC)(&Nc);
      } else if (stringtoint(tokens[2], Nc) && (Nc>0) ) {
        FC_FUNC_(blha_setnc,BLHA_SETNC)(&Nc);
      } else {
        result = "Error: unsupported flag\n# any non-negative integer or \"Infinity\" is supported";
        *ierr = -1;
      }
    }
    else {
      result = "Error: unsupported flag\n# HelAvgInitial, ColAvgInitial, MCSymmetrizeFinal, Nf, Nc is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "modelfile") {
    result = "Error: flag not yet supported\n# Use OLP_SetParameter instead";
    *ierr = -1;
  }
  else if (stringlower(tokens[0]) == "operationmode") {
    if (stringlower(tokens[1]) != "couplingsstrippedoff") {
      result = "Error: unsupported flag\n# CouplingsStrippedOff is supported";
      *ierr = -1;
    }
    FC_FUNC_(blha_setcouplingsoff,BLHA_SETCOUPLINGSOFF)();
  }
  else if (stringlower(tokens[0]) == "subdividesubprocess") {
    if (stringlower(tokens[1]) != "no") {
      result = "Error: unsupported flag\n# no is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "alphaspower") {
    if (!stringtoint(tokens[1],orderAlphas)) {
      result = "Error: invalid syntax\n# Value not an integer?";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "alphapower") {
    if (!stringtoint(tokens[1],orderAlpha)) {
      result = "Error: invalid syntax\n# Value not an integer?";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "resonancetreatment") {
    if (stringlower(tokens[1]) != "fixedwidthscheme") {
      result = "Error: unsupported flag\n# FixedWidthScheme is supported";
      *ierr = -1;
    }
  }
  else if (stringlower(tokens[0]) == "ewrenormalisationscheme") {
    int scheme = 0;
    if (stringlower(tokens[1]) == "alpha0") {
      scheme = 1;
    } else if (stringlower(tokens[1]) == "alphamz") {
      scheme = 2;
    } else if (stringlower(tokens[1]) == "alphaGF") {
      scheme = 3;
    } else {
      result = "Error: unsupported flag\n# alpha0, alphaMZ, alphaGF is supported";
      *ierr = -1;
    }
    FC_FUNC_(blha_setewrenormscheme,BLHA_SETEWRENORMSCHEME)(&scheme);
  }
  else if (stringlower(tokens[0]) == "amplitudetype") {
    if (stringlower(tokens[1]) == "tree") {
      AmplitudeType = 0;
    } else if (stringlower(tokens[1]) == "loop") {
      AmplitudeType = 1;
    } else if (stringlower(tokens[1]) == "cctree") {
      AmplitudeType = 2;
    } else if (stringlower(tokens[1]) == "sctree") {
      AmplitudeType = 3;
    } else {
      result = "Error: unsupported flag\n# tree, loop, cctree or sctree is supported";
      *ierr = -1;
    }
  }
  else if (stringtoint(tokens[0],pdgId)) {
    DetProcess(tokens,result,ierr);
  }
  else {
    int iparamerr;
    double reval=0, imval=0;
    if (stringtodouble(tokens[1],reval)) {
      OLP_SetParameter(const_cast<char*>(tokens[0].c_str()), &reval, &imval, &iparamerr);
      if (iparamerr==1 || iparamerr==2) {
      // everything ok or ignored
      } else if (iparamerr==0) {
        result = "Error: bad option";
        *ierr = -1;
      } else {
        result = "Error: unknown option";
        *ierr = -1;
      }
    } else {
      result = "Error: unknown option";
      *ierr = -1;
    }
  }

}

void BLHA::DetProcess(vector<string> tokens, string& result, int* ierr) {

  if (AmplitudeType < 0) {
    if (BLHAtype == 1) {
      AmplitudeType = -1; // BLHA1 - loop
    } else {
      result = "Error: Loop-level of process unknown";
      *ierr = -1;
      return;
    }
  }
  if (tokens.size() <= 4) {
    result = "Error: unsupported process";
    *ierr = -1;
    return;
  } 
  if ( tokens[2] != "->" ) {
    result = "Error: unsupported process\n# only 2->n processes available";
    *ierr = -1;
    return;
  }
  tokens.erase(tokens.begin()+2);
  int nparticles = tokens.size();
  int pdgprocess[nparticles];
  for (int i=0;i<nparticles;i++) {
    if (!stringtoint(tokens[i],pdgprocess[i])) {
      result = "Error: unsupported token in process: "+tokens[i];
      *ierr = -1;
      break;
    }
  }
  // Now try to infer from the PDG codes what process we want
  // currently only VBF-H implemented
  int procok;
  FC_FUNC_(vbfnlo_setupprocess,VBFNLO_SETUPPROCESS)(&nparticles,pdgprocess,&orderAlphas,&orderAlpha,&AmplitudeType,&procok);
  if (procok==0) {
    result = "Error: unsupported process";
    *ierr = -1;
    return;
  } else if (procok==-1) {
    result = "Error: invalid process";
    *ierr = -1;
    return;
  } else if (procok==-2) {
    result = "Error: too many subprocesses - increase max_blhaproc in VBFNLO";
    *ierr = -1;
    return;
  } else {
    ostringstream rs;
    rs << "1 " << procok;
    result = rs.str();
  }
}

bool BLHA::stringtoint(const string str, int& i) {
  char* endptr;
  
  i = strtol(str.c_str(),&endptr,0);
  return ((*endptr) == '\0');
}

bool BLHA::stringtodouble(const string str, double& d) {
  char* endptr;
  
  d = strtod(str.c_str(),&endptr);
  return ((*endptr) == '\0');
}

string BLHA::stringlower(const string str) {
  string out=str;
  transform(out.begin(),out.end(),out.begin(),::tolower); 
  return out;
}

}

