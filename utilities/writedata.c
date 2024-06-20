#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

FILE* datafile = NULL;
int oldtype = 0;
int seed = 0;
char** filelist = NULL;
int* ifilelist = NULL;
int ifilelisttot = 0;
int nfilelist = 0;
int nofile=0;

void opendatafile_(int* type);
void closedatafile_();
void writedata_(int*, double*, double*, int*, double*, int*);
void readdata_(int*, double*, double*, int*, double*);
void numdata_(long*, int*, int*, int*);
void getfilelist(int);
void nextfile();

void opendatafile_(int* type) {
  char filename[200];
  int ftype;
  if (datafile == NULL) {
    seed = *type;
    ftype = 1;
  } else {
    closedatafile_();
    ftype = *type;
  }
  sprintf(filename,"binarydata_%03d.%d.out",ftype,seed);
  datafile = fopen(filename,"w");
}

void closedatafile_() {
  fclose(datafile);
  datafile=NULL;
}

void writedata_(int *numentries, double *ran, double *weight, int *nampl, double *ampl, int *type) {
  int i;
  if (oldtype != *type) {
    opendatafile_(type);
    oldtype = *type;
  }
  float *fran = malloc(*numentries*sizeof(float));
  for (i=0;i<*numentries;i++) { fran[i] = ran[i]; }
  float fweight = *weight;
  float *fampl = malloc(*nampl*sizeof(float));
  for (i=0;i<*nampl;i++) { fampl[i] = ampl[i]; }
  fwrite(fran,sizeof(float),*numentries,datafile);
  fwrite(&fweight,sizeof(float),1,datafile);
  fwrite(fampl,sizeof(float),*nampl,datafile);
  free(fran);
  free(fampl);
}

void readdata_(int *numentries, double *ran, double *weight, int *nampl, double *ampl) {
  int i;
  float *fran = malloc(*numentries*sizeof(float));
  float fweight;
  float *fampl = malloc(*nampl*sizeof(float));
  while (fread(fran,sizeof(float),*numentries,datafile)==0) {
    nextfile();
  }
  fread(&fweight,sizeof(float),1,datafile);
  fread(fampl,sizeof(float),*nampl,datafile);
  for (i=0;i<*numentries;i++) { ran[i] = fran[i]; }
  *weight = fweight*ifilelist[nofile]/(double)ifilelisttot;
  if (*nampl>1) {
    ampl[0]=0;
    for (i=0;i<*nampl;i++) { 
      ampl[i+1] = fampl[i]; 
      ampl[0] += ampl[i+1];
   }
  } else { 
    ampl[0] = fampl[0]; 
  }
  free(fran);
  free(fampl);
}

void numdata_(long *numpts, int *numentries, int *nampl, int *type) {
  char filename[200];
  int i;
  getfilelist(*type);
  ifilelisttot = 0;
  for (i=0;i<nfilelist;i++) {
    if (datafile != NULL) closedatafile_();
    datafile = fopen(filelist[i],"r");
    if (datafile==NULL) {
      printf("Error: Cannot open datafile %s\n",filename);
      exit(1);
    }
    fseek(datafile,0,SEEK_END);
    ifilelist[i] = ftell(datafile)/(sizeof(float)*(*numentries+1+*nampl));
    ifilelisttot += ifilelist[i];
  }
  *numpts = ifilelisttot;
  nofile=-1;
  nextfile();
}

void getfilelist(int ftype) {
  DIR *dir;
  struct dirent *ent;
  char filename[200];
  int i;

  if (filelist!=NULL) {
    for (i=0;i<nfilelist;i++) { free(filelist[i]); };
    filelist = realloc(filelist,0);
    ifilelist = realloc(ifilelist,0);
    nfilelist = 0;
  }
  dir = opendir (".");
  if (dir == NULL) { exit(1); };

  sprintf(filename,"binarydata_%03d.",ftype);
  while ((ent = readdir (dir)) != NULL) {
    if(strncmp(ent->d_name,filename,14)==0) {
      nfilelist++;
      filelist = realloc(filelist,nfilelist*sizeof(char*));
      ifilelist = realloc(ifilelist,nfilelist*sizeof(int));
      filelist[nfilelist-1] = malloc(100*sizeof(char));
      strncpy(filelist[nfilelist-1],ent->d_name,100);
    }
  }
  closedir(dir);
}

void nextfile() {
  if (datafile != NULL) closedatafile_();
  nofile++;
  datafile = fopen(filelist[nofile],"r");
}
