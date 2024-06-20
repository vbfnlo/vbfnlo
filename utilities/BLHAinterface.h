/*
 * Public Interface for the Binoth Les Houches Accord
 *
 * Author: Michael Rauch <michael.rauch@kit.edu>
 * Initial version: Nov 2013
 * Last modified: Mar 2015
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(WITH_CUSTOM_BLHA) || defined(VBFNLO3)
/*
 * External function definitions from arXiv:1308.3462
 */
void OLP_Start_VBFNLO(char* fname, int* ierr);
void OLP_Info_VBFNLO(char olp_name[15], char olp_version[15], char message[255]);
void OLP_SetParameter_VBFNLO(char* line, double* re, double* im, int* ierr);
void OLP_Polvec_VBFNLO(double* p, double* q, double* eps); // optional
void OLP_EvalSubProcess_VBFNLO(int i, double* pp, double mu, double* alphas, double* rval);
void OLP_EvalSubProcess2_VBFNLO(int* i, double* pp, double* mu, double* rval, double* acc);

/* 
 * Extra function for order-contract phase
 */
void OLP_Order_VBFNLO(char* inname, char* outname, int* ierr);
/* 
 * Extra functions for phasespace
 */
void OLP_GetParameter_VBFNLO(char* line, double* re, double* im, int* ierr);
void OLP_PhaseSpacePoint_VBFNLO(int* proc, double* rpsnum, double* r, double* p, double* weight);
#else
/*
 * External function definitions from arXiv:1308.3462
 */
void OLP_Start(char* fname, int* ierr);
void OLP_Info(char olp_name[15], char olp_version[15], char message[255]);
void OLP_SetParameter(char* line, double* re, double* im, int* ierr);
void OLP_Polvec(double* p, double* q, double* eps); // optional
void OLP_EvalSubProcess(int i, double* pp, double mu, double* alphas, double* rval);
void OLP_EvalSubProcess2(int* i, double* pp, double* mu, double* rval, double* acc);

/* 
 * Extra function for order-contract phase
 */
void OLP_Order(char* inname, char* outname, int* ierr);
/* 
 * Extra functions for phasespace
 */
void OLP_GetParameter(char* line, double* re, double* im, int* ierr);
void OLP_PhaseSpacePoint(int* proc, double* rpsnum, double* r, double* p, double* weight);
#endif

#ifdef __cplusplus
}
#endif

