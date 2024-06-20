#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "gsl/gsl_integration.h"


extern "C" 
{
    void calckkcoupl_( double*, double*, double* , double* , int* , int* , int* );
}


using namespace std;

FILE *kkdat;


// R(rdwn) and R-tilde(rup)

double R1( double M,double R ){
       return
       y1(M*R)/j1(M*R);
}

double R0( double M,double R ){
       return
       y0(M*R)/j0(M*R);
}

// Z-Boson solution as a function of -b^(B)_k

double aka1( double M,double rup, double rdwn,double g5, double g5t){
       return
       R0(M,rup);
} 

double akzl1( double M,double rup, double rdwn,double g5, double g5t){ 
       double k=g5/g5t;
       return
       (k*R0(M,rdwn)*(-R0(M,rup) + R1(M,rup)))/(2.*(R0(M,rdwn) - R1(M,rup)));
}

double  bkzl1( double M,double rup, double rdwn,double g5,double g5t){
       double k=g5/g5t;
       return
       (k*(R0(M,rup) - R1(M,rup)))/(2.*(R0(M,rdwn) - R1(M,rup)));
}

double  akzr1( double M,double rup, double rdwn,double g5, double g5t){
       double k=g5/g5t;
       return
       -(k*(-2*R0(M,rup)*R1(M,rup) + R0(M,rdwn)*(R0(M,rup) + R1(M,rup))))/(2.*(R0(M,rdwn) - R1(M,rup)));      
}

double  bkzr1( double M,double rup, double rdwn,double g5, double g5t){
       double k=g5/g5t;
       return
       -(k*(-2.*R0(M,rdwn) + R0(M,rup) + R1(M,rup)))/(2.*(R0(M,rdwn) - R1(M,rup)));

}

// Gamma-Boson solution as a function of -b^(B)_k

double aka2( double M,double rup, double rdwn,double g5, double g5t){
       return
       R0(M,rdwn);
}

double  akzl2( double M,double rup, double rdwn,double g5, double g5t){   //Z-Wavefunction coeff sim Photon bb
       double k1=g5t/g5;
       return
       R0(M,rdwn)*k1;
}

double  akzr2( double M,double rup, double rdwn,double g5, double g5t){
       double k1=g5t/g5;
       return
       R0(M,rdwn)*k1;      
}

double  bkzl2( double M,double rup, double rdwn,double g5,double g5t){
       return
       -g5t/g5;
}

double  bkzr2( double M,double rup, double rdwn,double g5, double g5t){
//       double k=g5/g5t;
       return
       -g5t/g5;
}

// W - Boson solution as a function of -b^(R+/-)_k

double akls( double M,double rup, double rdwn){
     return
     R0(M,rup) + ((R1(M,rdwn) - R0(M,rup))*(R0(M,rup) + R1(M,rup)))/(R0(M,rup) - R1(M,rup));
}  

double bkls( double M, double rup, double rdwn){
     return
     (-2.*R1(M,rdwn) + R0(M,rup) + R1(M,rup))/(R0(M,rup) - R1(M,rup));     
}

double akrs( double M, double rup, double rdwn){
     return
     R1(M,rdwn);
}

// Wave functions, mass defining functions and stuff  

double waven( double z, void* params){ // double M, double ar,  double br,  double R1){
     double *alpha = (double *) params;
     return
     alpha[3]/z*(alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z));
}

double wave( double z, void* params){ //  double M, double ar,  double br){
     double *alpha = (double *) params;
     return
     (alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z));
}

double wavew(double z, void* params){ //  double M, double ar,  double br){
     double *alpha = (double *) params;
     return
     alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z)+
     alpha[3]*z*j1(alpha[0]*z) + alpha[4]*z*y1(alpha[0]*z);
}

double wavez(double z, void* params){ //  double M, double ar,  double br){
     double *alpha = (double *) params;
     return
     alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z)+
     alpha[3]*z*j1(alpha[0]*z) + alpha[4]*z*y1(alpha[0]*z)+
     alpha[5]*z*j1(alpha[0]*z) + alpha[6]*z*y1(alpha[0]*z) ;
}

double wavesqrt( double z, void* params ) { //double M, double ar,  double br,  double R1){
     double *alpha = (double *) params;
     return
     alpha[3]/z*
     pow(wave(z, params),2);
}
     
double wavesqrtw( double z, void* params ) { //double M, double ar,  double br,  double R1){
     double *alpha = (double *) params;
     return
     alpha[5]/z*pow(wavew(z, params),2);
}

double wavesqrtz( double z, void* params ) { //double M, double ar,  double br,  double R1){
     double *alpha = (double *) params;
     return
     alpha[7]/z*pow(wavez(z, params),2);
}
     
double Zeros( double M, double rup, double rdwn) {
     double x,y;
     x=M*rdwn;
     y=M*rup;
     return 
     -(j1(x)*(j1(y)*y0(x)*y0(y) + (j0(y)*y0(x) - 2.*j0(x)*y0(y))*y1(y))) - 
     y1(x)*(j0(x)*j1(y)*y0(y) + j0(y)*(-2.*j1(y)*y0(x) + j0(x)*y1(y)));
}

double zZeros( double M, double rup, double rdwn,double g52, double g5t2) {
        double y,x;
       y=M*rup;
       x=M*rdwn;
       return
       -(j1(x)*((g52 + 2.*g5t2)*j1(y)*y0(x)*y0(y) + 
        (g52*j0(y)*y0(x) - 2.*(g52 + g5t2)*j0(x)*y0(y))*y1(y))) 
       - y1(x)*(g52*j0(x)*j1(y)*y0(y) + 
       j0(y)*(-2.*(g52 + g5t2)*j1(y)*y0(x) + (g52 + 2.*g5t2)*j0(x)*y1(y)));
}     
  
double gZeros( double M, double rup, double rdwn) {
       return
       (-(j0(M*rup)*y0(M*rdwn)) + j0(M*rdwn)*y0(M*rup));
}

double coupling(double M, double rup, double rdwn, double g52) {
       return
       -g52*((R0(M,rup)-R0(M,rdwn))*(R1(M,rup)-R1(M,rdwn))+ (R1(M,rup)-R0(M,rdwn))*(R0(M,rup)-R1(M,rdwn)))/
        2./(R1(M,rup)-R0(M,rdwn))/(R0(M,rup)-R1(M,rdwn));
}

double Vertex(double z, void *params) { //params1: alpha (0-2),beta (a 3-5) ,
      double *alpha = (double *) params;  //         gamma (6-8),lambda (a 9-11) 
                                          //         epsilon(12-14), rdwn (15)
      return
      alpha[18]/z*(
        (alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z))*                        //Wave_{kL}
        (alpha[4]*z*j1(alpha[3]*z) + alpha[5]*z*y1(alpha[3]*z))*                        //Wave_{lL}
        (alpha[7]*z*j1(alpha[6]*z) + alpha[8]*z*y1(alpha[6]*z))+                        //Wave_{mL}
        (alpha[10]*z*j1(alpha[9]*z) + alpha[11]*z*y1(alpha[9]*z))*                        //Wave_{kR}
        (alpha[13]*z*j1(alpha[12]*z) + alpha[14]*z*y1(alpha[12]*z))*                      //Wave_{lR}
        (alpha[16]*z*j1(alpha[15]*z) + alpha[17]*z*y1(alpha[15]*z)));                     //Wave_{mR}
}


double Vertex4(double z, void *params) { //params1: alpha (0-2),beta (a 3-5) ,
      double *alpha = (double *) params;  //         gamma (6-8),lambda (a 9-11) 
                                          //         epsilon(12-14), rdwn (15)
      return
       alpha[24]/z*(
       (alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z))*                       //Wave_{kL}
       (alpha[4]*z*j1(alpha[3]*z) + alpha[5]*z*y1(alpha[3]*z))*                       //Wave_{lL}
       (alpha[7]*z*j1(alpha[6]*z) + alpha[8]*z*y1(alpha[6]*z))*                       //Wave_{mL}
       (alpha[10]*z*j1(alpha[9]*z) + alpha[11]*z*y1(alpha[9]*z))+                     //Wave_{nL}
       (alpha[13]*z*j1(alpha[12]*z) + alpha[14]*z*y1(alpha[12]*z))*                       //Wave_{kR}
       (alpha[16]*z*j1(alpha[15]*z) + alpha[17]*z*y1(alpha[15]*z))*                       //Wave_{lR}
       (alpha[19]*z*j1(alpha[18]*z) + alpha[20]*z*y1(alpha[18]*z))*                       //Wave_{mR}
       (alpha[22]*z*j1(alpha[21]*z) + alpha[23]*z*y1(alpha[21]*z)));                      //Wave_{nR}
}

//------------------------------------------------------------------------------------------------------
double testorth(double z, void *params) { //params1: alpha (0-2),beta (a 3-5) ,
      double *alpha = (double *) params;  //         gamma (6-8),lambda (a 9-11) 
                                          //         epsilon(12-14), rdwn (15)
      return
       alpha[6]/z*((alpha[1]*z*j1(alpha[0]*z) + alpha[2]*z*y1(alpha[0]*z))*
                     (alpha[4]*z*j1(alpha[3]*z) + alpha[5]*z*y1(alpha[3]*z)));
}

void calckkcoupl_( double *ggrf, double *wmassrf, 
       double *zmassrf, double *rdwnrf, int *kkmaxwrf, int *kkmaxzrf, int *kkmaxgrf){
   
   double rup,norm1,norm2,norm3,norm4,norm5,norm6,norm7,norm8,normph;
   double a,c,b,d,g52,g5t2,g5,g5t,gwwa2,
          az,cz,bz,dz,gg,wmass,zmass,rdwn,
          ag,cg,bg,dg,fg;
   int kkmaxw,kkmaxz,kkmaxg;
   double const pi(3.14159265358979323);
   double alpha[7],beta[4],gamma[4],lambda[4],epsilon[4],delta[4],omega[4],sigma[4],vert[19];
   double wwz1[20][20][20],wwz2[20][20][20];
   double fvert[25];
   double massw[20],massz[20],massg[20];
   double normw[20],normz1[20],normz2[20];
   double result1,result2;
   double wwfvrt[20][20][20][20],wzfvrt[20][20][20][20];
   double wgfvrt[20][20][20][20],wgzfvrt[20][20][20][20];
   double cutoff;
   double error,sum1,sum1w,sum2w,sum2z,sum1z;
   int i,k,cnt;
   bool choose;

   gsl_function F,verfu;
   gsl_integration_workspace * w   = gsl_integration_workspace_alloc (1000000);
   
   sum1=0;
   sum2w=0;
   c=0;
   cg=0;
   bg=0;
   ag=0;

   wmass=*wmassrf;
   zmass=*zmassrf;
   rdwn=*rdwnrf;
   gg=*ggrf;
   kkmaxw = *kkmaxwrf;
   kkmaxz = *kkmaxzrf;
   kkmaxg = *kkmaxgrf;
   gwwa2=pow(gg,2); 
   choose=false;
   cutoff=3000.;

   kkdat=fopen("kk_coupl_inp.dat", "wt");
//----------------------------------------------------------------------------------------------------------------
// find appropriate rup as a function of wmass and rdwn
   a=0.001;
   b=0.001;
   do {b+=1e-7;} while (Zeros(wmass,b,rdwn)*Zeros(wmass,a,rdwn)>0);
   do {
   d=c;
   c=a-(b-a)/(Zeros(wmass,b,rdwn)-Zeros(wmass,a,rdwn))*Zeros(wmass,a,rdwn);
   if(Zeros(wmass,a,rdwn)*Zeros(wmass,c,rdwn) > 0) {a=c;}
   if(Zeros(wmass,b,rdwn)*Zeros(wmass,c,rdwn) > 0) {b=c;}
   } while (fabs(d-c) > 1e-10);
   rup=d;  


//---------------------------------------------------------------------------------------------------

// calculate g5 by matching g5* \int rdwn/z*\psi_w**2 \psi_z = g4


// physical photon normalization and fixed WWA -coupling

   g52=(-2*gwwa2*rdwn*(R0(zmass,rdwn) - R0(zmass,rup))*(R1(zmass,rdwn) - R1(zmass,rup))*log(rup/rdwn))
         /(-2.*R0(zmass,rup)*R1(zmass,rup) + R1(zmass,rdwn)*(R0(zmass,rup) + R1(zmass,rup)) 
        + R0(zmass,rdwn)*(-2.*R1(zmass,rdwn) + R0(zmass,rup) + R1(zmass,rup)));
   g5t2=(g52*gwwa2*rdwn*log(rup/rdwn))/(g52 - 2.*gwwa2*rdwn*log(rup/rdwn));
   g5=sqrt(g52);
   g5t=sqrt(g5t2);  
 
   normph=1./sqrt( (g52+2.*g5t2)*rdwn*log(rup/rdwn) );


// Left Stm-W normalization

   alpha[0]=wmass;
   alpha[1]=akls(wmass,rup,rdwn);
   alpha[2]=bkls(wmass,rup,rdwn);   
   alpha[3]=rdwn;

   F.function = &wavesqrt;
   F.params = &alpha; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm1, &error);

// Right Stm-W normalization

   beta[0]=wmass;
   beta[1]=akrs(wmass,rup,rdwn);
   beta[2]=-1.;   
   beta[3]=rdwn;

   F.params = &beta; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm2, &error);   

// Left Stm-Z normalization

   gamma[0]=zmass;
   gamma[1]=akzl1(zmass,rup,rdwn,g5,g5t);
   gamma[2]=bkzl1(zmass,rup,rdwn,g5,g5t);   
   gamma[3]=rdwn;

   F.params = &gamma; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm3, &error);

// Right Stm-Z normalization

   lambda[0]=zmass;
   lambda[1]=akzr1(zmass,rup,rdwn,g5,g5t);
   lambda[2]=bkzr1(zmass,rup,rdwn,g5,g5t);   
   lambda[3]=rdwn;

   F.params = &lambda; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm4, &error);  

// B contribution to Z1

   epsilon[0]=zmass;
   epsilon[1]=aka1(zmass,rup,rdwn,g5,g5t);
   epsilon[2]=-1.;   
   epsilon[3]=rdwn;

   F.params = &epsilon;

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm5, &error); 

   vert[0]=wmass;
   vert[1]=akls(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[2]=bkls(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[3]=wmass;
   vert[4]=akls(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[5]=bkls(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[6]=zmass;
   vert[7]=akzl1(zmass,rup,rdwn,g5,g5t)/sqrt(norm3+norm4+norm5);
   vert[8]=bkzl1(zmass,rup,rdwn,g5,g5t)/sqrt(norm3+norm4+norm5);

   vert[9]=wmass;
   vert[10]=akrs(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[11]=(-1.)/sqrt(norm1+norm2); 
   vert[12]=wmass;
   vert[13]=akrs(wmass,rup,rdwn)/sqrt(norm1+norm2);
   vert[14]=(-1.)/sqrt(norm1+norm2);
   vert[15]=zmass;
   vert[16]=akzr1(zmass,rup,rdwn,g5,g5t)/sqrt(norm3+norm4+norm5);
   vert[17]=bkzr1(zmass,rup,rdwn,g5,g5t)/sqrt(norm3+norm4+norm5);
   vert[18]=rdwn; 

   verfu.function = &Vertex;
   verfu.params = &vert;

   gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-10, 1000000,w, &result1, &error); 

//-------------------------------------------------------------------------------------------

// W BOSONS

// Regula falsi for zeros using Bessel function Expansion

  for (i=0;i<10;i++){
  if (i==0) {a=77;}
  else {a=massw[i-1]+50.;}
  b=a;
  do {b=b+10.;} while (Zeros(b,rup,rdwn)*Zeros(a,rup,rdwn)>0);
  do {
  d=c;
  c=a-(b-a)/(Zeros(b,rup,rdwn)-Zeros(a,rup,rdwn))*Zeros(a,rup,rdwn);
  if(Zeros(a,rup,rdwn)*Zeros(c,rup,rdwn) > 0) {a=c;}
  if(Zeros(b,rup,rdwn)*Zeros(c,rup,rdwn) > 0) {b=c;}
  } while (fabs(d-c) > 1e-14);
  massw[i]=c;
   
  if (fabs(massw[0]-wmass)>1e-5) {cout << "W mass insufficiently reconstructed! " << endl; return;}

// Left W normalization
   alpha[0]=massw[i];
   alpha[1]=akls(massw[i],rup,rdwn);
   alpha[2]=bkls(massw[i],rup,rdwn);   
   alpha[3]=rdwn;

   F.function = &wavesqrt;
   F.params = &alpha; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm1, &error);

// Right W normalization

   beta[0]=massw[i];
   beta[1]=akrs(massw[i],rup,rdwn);
   beta[2]=-1.;   
   beta[3]=rdwn;

   F.params = &beta; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-10, 1000000,w, &norm2, &error);   

   normw[i]=norm1+norm2;
}


//-------------------------------------------------------------------------------------

// Z BOSONS

   for (k=0;k<12;k++){
   if(k==0) {az=89;}
   if(k==1) {az=pi/2*(k+0.5)/rup;}
   else {az=massz[k-1]+50;}
   bz=az;
   cz=az;
   do {bz=bz+10;} while (zZeros(bz,rup,rdwn,g52,g5t2)*zZeros(az,rup,rdwn,g52,g5t2)>0);
   az=bz-10.;
   do 
     { dz=cz;
       cz=az-(bz-az)/(zZeros(bz,rup,rdwn,g52,g5t2)-zZeros(az,rup,rdwn,g52,g5t2))*zZeros(az,rup,rdwn,g52,g5t2);
       if(zZeros(az,rup,rdwn,g52,g5t2)*zZeros(cz,rup,rdwn,g52,g5t2) > 0) {az=cz;}
       if(zZeros(bz,rup,rdwn,g52,g5t2)*zZeros(cz,rup,rdwn,g52,g5t2) > 0) {bz=cz;}
     } 
   while (fabs(dz-cz) > 1e-16);
   massz[k]=cz;

   }
   if (fabs(massz[0]-zmass)>1e-10) {cout << "Z mass insufficiently reconstructed! " << endl; return;}
    
   
   for (k=1;k<12;k++){
   if(k==1){ ag=pi/2*(k+0.5)/rup;bg=ag;cg=ag;}
   if(k>1){ ag=massg[k-1]+50.;bg=ag;cg=ag;}
   do {bg=bg+10;} while (gZeros(bg,rup,rdwn)*gZeros(ag,rup,rdwn)>0);
   ag=bg-10;
   do 
   {
   dg=cg;
   cg=ag-(bg-ag)/(gZeros(bg,rup,rdwn)-gZeros(ag,rup,rdwn))*gZeros(ag,rup,rdwn);
   if(gZeros(ag,rup,rdwn)*gZeros(cg,rup,rdwn) > 0) {ag=cg;}
   if(gZeros(bg,rup,rdwn)*gZeros(cg,rup,rdwn) > 0) {bg=cg;}
   } while (fabs(dg-cg) > 1e-16);
   fg=dg;
   massg[k]=dg;
   }
   

   for (k=0;k<12;k++){

// Left Z normalization

   gamma[0]=massz[k];
   gamma[1]=akzl1(massz[k],rup,rdwn,g5,g5t);
   gamma[2]=bkzl1(massz[k],rup,rdwn,g5,g5t);   
   gamma[3]=rdwn;

   F.function = &wavesqrt;
   F.params = &gamma; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7, 1000000,w, &norm3, &error);
 
// Right Z normalization

   lambda[0]=massz[k];
   lambda[1]=akzr1(massz[k],rup,rdwn,g5,g5t);
   lambda[2]=bkzr1(massz[k],rup,rdwn,g5,g5t);   
   lambda[3]=rdwn;

   F.params = &lambda; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7, 1000000,w, &norm4, &error);   


// Photon contribution

   epsilon[0]=massz[k];
   epsilon[1]=aka1(massz[k],rup,rdwn,g5,g5t);
   epsilon[2]=-1.;   
   epsilon[3]=rdwn;

   F.params = &epsilon;

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7,1000000,w, &norm5, &error); 


   normz1[k]=norm3+norm4+norm5;

//-------------------------------------------------------------

// Left Z normalization

   if(k>0){
   sigma[0]=massg[k];
   sigma[1]=akzl2(massg[k],rup,rdwn,g5,g5t);
   sigma[2]=bkzl2(massg[k],rup,rdwn,g5,g5t);   
   sigma[3]=rdwn;

   F.params = &sigma; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7, 1000000,w, &norm6, &error);


// Right Z normalization

   omega[0]=massg[k];
   omega[1]=akzr2(massg[k],rup,rdwn,g5,g5t);
   omega[2]=bkzr2(massg[k],rup,rdwn,g5,g5t);   
   omega[3]=rdwn;

   F.params = &omega; 

   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7, 1000000,w, &norm7, &error);   

   
// Photon contribution

   delta[0]=massg[k];
   delta[1]=aka2(massg[k],rup,rdwn,g5,g5t);
   delta[2]=-1.;   
   delta[3]=rdwn;

   F.params = &delta;
 
   gsl_integration_qags (&F, rdwn,rup, 0, 1e-7, 1000000,w, &norm8, &error); 

   normz2[k]=norm6+norm7+norm8;
   }
}   

//-----------------------------------------------------------------------------------

// find cut-off related kk-index
  
  if(choose){
    cnt=0;
    do {cnt=cnt+1;} while (massw[cnt]<=cutoff);
    kkmaxw=cnt-1;
    cout << cutoff << " "<< " " << cnt+1 << " " << massw[kkmaxw] << endl;
    cnt=0;
    do {cnt=cnt+1;} while (massz[cnt]<=cutoff);
    kkmaxz=cnt-1;
    cout << cutoff << " "<< " " << cnt+1 << " " << massz[kkmaxz]  << endl;
    cnt=0;
    do {cnt=cnt+1;} while (massg[cnt]<=cutoff);
    kkmaxg=cnt-1;
    cout << cutoff << " "<< " " << cnt+1 << " " << massg[kkmaxg] << endl;
   } 

//-------------------------------------------------------------------------------------

//  4-Vertices

// W^4 - vertices

   for (k=0;k<=kkmaxw;k++) {
    for (int l=0;l<=k;l++) {
     for (int m=0;m<=l;m++) {
      for (int n=0;n<=m;n++) { 

   	fvert[0]=massw[k];
   	fvert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]); 
   	fvert[3]=massw[l];
   	fvert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[6]=massw[m];
   	fvert[7]=akls(massw[m],rup,rdwn)/sqrt(normw[m]);
   	fvert[8]=bkls(massw[m],rup,rdwn)/sqrt(normw[m]);
   	fvert[9]=massw[n];
   	fvert[10]=akls(massw[n],rup,rdwn)/sqrt(normw[n]);
   	fvert[11]=bkls(massw[n],rup,rdwn)/sqrt(normw[n]);  
   	fvert[12]=massw[k];
   	fvert[13]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[14]=(-1.)/sqrt(normw[k]);
   	fvert[15]=massw[l];
   	fvert[16]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[17]=(-1.)/sqrt(normw[l]);
   	fvert[18]=massw[m];
   	fvert[19]=akrs(massw[m],rup,rdwn)/sqrt(normw[m]);
   	fvert[20]=(-1.)/sqrt(normw[m]);
   	fvert[21]=massw[n];
   	fvert[22]=akrs(massw[n],rup,rdwn)/sqrt(normw[n]);
   	fvert[23]=(-1.)/sqrt(normw[n]);
   	fvert[24]=rdwn;

   	verfu.function= &Vertex4;
   	verfu.params = &fvert;
   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-10, 1000000,w, &result2, &error);
   	wwfvrt[k][l][m][n]=fabs(result2*pow(g5,2));
        }
       }
      }
     }

//  W^2 Z^2 -vertices

//   W^2 Z^2 vertices ( WWZZ, WW-more-Z more-Z )

   for (k=0;k<=kkmaxw;k++) {
    for (int l=0;l<=k;l++) {
     for (int m=0;m<=kkmaxz;m++) {
      for (int n=0;n<=m;n++) { 

   	fvert[0]=massw[k];
   	fvert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]); 
   	fvert[3]=massw[l];
   	fvert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[6]=massz[m];
   	fvert[7]=akzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	fvert[8]=bkzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
        fvert[9]=massz[n];
   	fvert[10]=akzl1(massz[n],rup,rdwn,g5,g5t)/sqrt(normz1[n]);
   	fvert[11]=bkzl1(massz[n],rup,rdwn,g5,g5t)/sqrt(normz1[n]);
   	fvert[12]=massw[m];
   	fvert[13]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[14]=(-1.)/sqrt(normw[k]);
   	fvert[15]=massw[l];
   	fvert[16]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[17]=(-1.)/sqrt(normw[l]);	
        fvert[18]=massz[m];
   	fvert[19]=akzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	fvert[20]=bkzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
        fvert[21]=massz[n];
   	fvert[22]=akzr1(massz[n],rup,rdwn,g5,g5t)/sqrt(normz1[n]);
   	fvert[23]=bkzr1(massz[n],rup,rdwn,g5,g5t)/sqrt(normz1[n]);
   	fvert[24]=rdwn;

   	verfu.function= &Vertex4;
   	verfu.params = &fvert;
   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-10, 1000000,w, &result2, &error);
   	wzfvrt[k][l][m][n]=fabs(result2*pow(g5,2));

  
        if (m>=1 and n>=1 and m<=kkmaxg and n<=kkmaxg) {
 
        fvert[0]=massw[k];
   	fvert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]); 
   	fvert[3]=massw[l];
   	fvert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[6]=massg[m];
   	fvert[7]=akzl2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	fvert[8]=bkzl2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
        fvert[9]=massg[n];
   	fvert[10]=akzl2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[11]=bkzl2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[12]=massw[m];
   	fvert[13]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[14]=(-1.)/sqrt(normw[k]);
   	fvert[15]=massw[l];
   	fvert[16]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[17]=(-1.)/sqrt(normw[l]);	
        fvert[18]=massg[m];
   	fvert[19]=akzr2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	fvert[20]=bkzr2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
        fvert[21]=massg[n];
   	fvert[22]=akzr2(massz[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[23]=bkzr2(massz[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[24]=rdwn;

   	verfu.function= &Vertex4;
   	verfu.params = &fvert;
   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-10, 1000000,w, &result2, &error);
   	wgfvrt[k][l][m][n]=fabs(result2*pow(g5,2));
        }
      }
     }
    }
   }


//  W^2 Z^2 Vertex ( WW more-B-Z Z )

   for (k=0;k<=kkmaxw;k++) {
    for (int l=0;l<=k;l++) {
     for (int m=0;m<=kkmaxz;m++) {
      for (int n=1;n<=kkmaxg;n++) { 

   	fvert[0]=massw[k];
   	fvert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]); 
   	fvert[3]=massw[l];
   	fvert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[6]=massz[m];
   	fvert[7]=akzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	fvert[8]=bkzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
        fvert[9]=massg[n];
   	fvert[10]=akzl2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[11]=bkzl2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[12]=massw[m];
   	fvert[13]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	fvert[14]=(-1.)/sqrt(normw[k]);
   	fvert[15]=massw[l];
   	fvert[16]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	fvert[17]=(-1.)/sqrt(normw[l]);	
        fvert[18]=massz[m];
   	fvert[19]=akzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	fvert[20]=bkzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
        fvert[21]=massg[n];
   	fvert[22]=akzr2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[23]=bkzr2(massg[n],rup,rdwn,g5,g5t)/sqrt(normz2[n]);
   	fvert[24]=rdwn;

   	verfu.function= &Vertex4;
   	verfu.params = &fvert;
   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-10, 1000000,w, &result2, &error);
   	wgzfvrt[k][l][m][n]=fabs(result2*pow(g5,2));
       }
       }
      }
     }

//---------------------------------------------------------------------------------------
//  3-Vertices    
   
   for (k=0;k<=kkmaxw;k++){
    for (int l=0;l<=k;l++){
     for (int m=0;m<=kkmaxz;m++){
   	vert[0]=massw[k];
   	vert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[3]=massw[l];
   	vert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[6]=massz[m];
   	vert[7]=akzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	vert[8]=bkzl1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	vert[9]=massw[k];
   	vert[10]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[11]=(-1.)/sqrt(normw[k]); 
   	vert[12]=massw[l];
   	vert[13]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[14]=(-1.)/sqrt(normw[l]);
   	vert[15]=massz[m];
   	vert[16]=akzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	vert[17]=bkzr1(massz[m],rup,rdwn,g5,g5t)/sqrt(normz1[m]);
   	vert[18]=rdwn; 

   	verfu.function = &Vertex;
   	verfu.params = &vert;

   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-7, 1000000,w, &result1, &error);
   	wwz1[k][l][m]=fabs(g5*result1);
        }
        }
        }


   for (k=0;k<=kkmaxw;k++){
    for (int l=0;l<=k;l++){
     for (int m=1;m<=kkmaxg;m++){
       

        vert[0]=massw[k];
   	vert[1]=akls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[2]=bkls(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[3]=massw[l];
   	vert[4]=akls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[5]=bkls(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[6]=massg[m];
   	vert[7]=akzl2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	vert[8]=bkzl2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	vert[9]=massw[k];
   	vert[10]=akrs(massw[k],rup,rdwn)/sqrt(normw[k]);
   	vert[11]=(-1.)/sqrt(normw[k]); 
   	vert[12]=massw[l];
   	vert[13]=akrs(massw[l],rup,rdwn)/sqrt(normw[l]);
   	vert[14]=(-1.)/sqrt(normw[l]);
   	vert[15]=massg[m];
   	vert[16]=akzr2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	vert[17]=bkzr2(massg[m],rup,rdwn,g5,g5t)/sqrt(normz2[m]);
   	vert[18]=rdwn; 

   	verfu.function = &Vertex;
   	verfu.params = &vert;

   	gsl_integration_qags (&verfu, rdwn,rup, 0, 1e-7, 1000000,w, &result1, &error);
   	wwz2[k][l][m]=fabs(g5*result1);
     
   	}
      }
      }


//  physical Photon coupling unchanged due to unbroken U(1)

   for (int l=0;l<=kkmaxw;l++) {
     for (int q=0;q<=l;q++) {
       if (l==q)  {
       wgfvrt[l][q][0][0]=pow(gg,2);}
       else
       {wgfvrt[l][q][0][0]=0.;}     // due to normalization, checked numerically
     }
    }

   for (int l=0;l<=kkmaxw;l++) {
     for (int q=0;q<=l;q++) { 
       for (int h=0;h<=kkmaxz;h++) { 
         wgzfvrt[l][q][h][0]=wwz1[l][q][h]*gg;
      }
     }
    }

    for (int l=0;l<=kkmaxw;l++) {
     for (int q=0;q<=l;q++) {  
       if(l==q){
       wwz2[l][q][0]=gg;}
       else {wwz2[l][q][0]=0.;} // due to normalization, checked numerically
     }
   }
    for (int l=0;l<=kkmaxw;l++) {
     for (int q=0;q<=l;q++) { 
       for (int h=1;h<=kkmaxg;h++) { 
         wgfvrt[l][q][h][0]=wwz2[l][q][h]*gg;
      }
     }
    }

//------------------------------------------------------------

//     Generate vbfnlo input
   fprintf(kkdat , "%s\t%d\n", "KKMAXW = ", kkmaxw );
   fprintf(kkdat , "%s\t%d\n", "KKMAXZ = ", kkmaxz );
   fprintf(kkdat , "%s\t%d\n", "KKMAXG = ", kkmaxg ); 

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "---Kaluza-Klein masses---\n" );
   for (int h=1;h<=kkmaxw; h++) { 
          if (h == 1){fprintf(kkdat , "%s" , "KK_MASSW =   ");}
     fprintf(kkdat , "%.9f%s" , massw[h], "   "); 
          if (h == kkmaxw) {fprintf(kkdat , "%s%d\n" , "!W2 to W" , h+1);}
   }
   for (int h=1;h<=kkmaxz; h++) { 
          if (h == 1){fprintf(kkdat , "%s" , "KK_MASSZ =   ");}
     fprintf(kkdat , "%.9f%s" , massz[h], "   "); 
          if (h == kkmaxz) {fprintf(kkdat , "%s%d\n" , "!Z2 to Z" , h+1);} 
   }
   for (int h=1;h<=kkmaxg; h++) { 
          if (h == 1){fprintf(kkdat , "%s" , "KK_MASSG =  ");}
     fprintf(kkdat , "%.9f%s" , massg[h], "   "); 
          if (h == kkmaxg) {fprintf(kkdat , "%s%d\n" , "!Z'1 to Z'" , h);} 
	  //   fprintf(kkdat , "%s%.9f\t%s%d\n" , "KK_MASSZ' = ", massg[h], "!Z'" , h); 
   }

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--RECALL: W0 = SM W, Z0 = SM Z, Z'0 = SM photon --\n" );
   fprintf(kkdat , "--WWZ coupl----\n" ); 

   for (int h=0;h<=kkmaxz; h++) { 
          if (h == 0){fprintf(kkdat , "%s" , "CPL_W0-W0-ZX = ");}
	  fprintf(kkdat , "%.9f%s", wwz1[0][0][h],"    "); 
          if (h == kkmaxz){fprintf(kkdat,"%s%d\n", "!W0-W0- Z0 to Z" , h );}
	  //          if (h == kkmaxz){fprintf(kkdat,"%s%d\n", "!W1-W1- Z1 to Z" , h+1 );}
   //   fprintf(kkdat , "%.9f\t%s%d%s%d%s%d\n", wwz1[0][0][h] , "!W" , 1 , "-W" , 1 , "-Z", h+1 );
   }    
   for (int h=1;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       fprintf(kkdat , "%s%d%s%d%s" , "CPL_W", h, "-W", l, "-ZX = ");
       //       fprintf(kkdat , "%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-ZX = ");
       for (int j=0;j<=kkmaxz; j++) {     
         fprintf(kkdat , "%.9f%s", wwz1[h][l][j] , "    " );
	 //         fprintf(kkdat , "%.9f\t%s%d%s%d%s%d\n", wwz1[h][l][j] , "!W" , h+1 , "-W" , l+1 , "-Z", j+1 );
         if (j == kkmaxz) {fprintf(kkdat , "%s%d%s%d%s%d\n", "!W", h, "-W", l, "- Z0 to Z", j);}
	 //         if (j == kkmaxz) {fprintf(kkdat , "%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "- Z1 to Z", j+1);}
       }
     }
   }

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--WWZ' coupl----\n" ); 
   for (int h=1;h<=kkmaxg; h++) { 
          if (h == 1){fprintf(kkdat , "%s" , "CPL_W0-W0-GX = ");}
   fprintf(kkdat , "%.9f%s", wwz2[0][0][h] , "    ");  
          if (h == kkmaxg){fprintf(kkdat,"%s%d\n", "!W0-W0- Z'0 to Z'" , h );}  
   //   fprintf(kkdat , "%.9f\t%s%d%s%d%s%d\n", wwz2[0][0][h] , "!W" , 1 , "-W" , 1 , "-Z'", h );    
   }
   for (int h=1;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       fprintf(kkdat , "%s%d%s%d%s" , "CPL_W", h, "-W", l, "-GX = ");
       //       fprintf(kkdat , "%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-GX = ");
       for (int j=0;j<=kkmaxg; j++) {     
         fprintf(kkdat , "%.9f%s", wwz2[h][l][j] , "    ");
	 //         fprintf(kkdat , "%.9f\t%s%d%s%d%s%d\n", wwz2[h][l][j] , "!W" , h+1 , "-W" , l+1 , "-Z'", j );
         if (j == kkmaxg) {fprintf(kkdat , "%s%d%s%d%s%d\n", "!W", h, "-W", l, "- Z'0 to Z'", j);}
	 //         if (j == kkmaxg) {fprintf(kkdat , "%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "- Z'0 to Z'", j);}
       }
     }
   }

   sum1w=pow(gg,2);
   for (k=0;k<=kkmaxz;k++) {
   sum1w+= pow(wwz1[0][0][k],2);
   sum2w+= 3*pow(massz[k],2)*pow(wwz1[0][0][k],2);
   }
   for (k=1;k<=kkmaxg;k++){
   sum1w+= pow(wwz2[0][0][k],2);
   sum2w+= 3*pow(massg[k],2)*pow(wwz2[0][0][k],2);
   }

   sum1w=pow(gg,2);
   for (k=0;k<=kkmaxz;k++) {
   sum1w+= pow(wwz1[0][0][k],2);
   sum2w+= 3*pow(massz[k],2)*pow(wwz1[0][0][k],2);
   }
   for (k=1;k<=kkmaxg;k++){
   sum1w+= pow(wwz2[0][0][k],2);
   sum2w+= 3*pow(massg[k],2)*pow(wwz2[0][0][k],2);
   }
   sum2w+=pow(wmass,2)*sqrt(fabs((4.*pow(wmass,2)*sum1w-sum2w)/pow(wmass,2)));


   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--WWWW coupl----\n" ); 
   for (int h=0;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       for (int j=0;j<=l; j++) { 
	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h, "-W", l, "-W", j, "-WX =  ");
	 //	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-W", j+1, "-WX =  ");
         for (int m=0;m<=j; m++) {
           if( fabs(wwfvrt[0][0][0][0]-sum1w)/wwfvrt[0][0][0][0]*100 > 0.1) { wwfvrt[0][0][0][0]=sum1w;}     
          fprintf(kkdat , "%.9f%s", 
           sqrt(fabs(wwfvrt[h][l][j][m])) , "    ");
	  //          fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	  //           sqrt(fabs(wwfvrt[h][l][j][m])) , "!W" , h+1 , "-W" , l+1 , "-W", j+1, "-W",m+1 );
	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h, "-W", l, "-W", j, "- W0 to W", m);}
	  //	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "-W", j+1, "- W1 to W", m+1);}
         }
       }
     }
   }

  sum1z=0.;
  sum2z=0.;
  for (k=0;k<=kkmaxw;k++) {
  sum1z+= pow(wwz1[k][0][0],2);
  sum2z+= 3*pow(massw[k],2)*pow(wwz1[k][0][0],2);
  }
 
   if( fabs(wzfvrt[0][0][0][0]-sum1z)/wzfvrt[0][0][0][0]*100 > 0.1) { wzfvrt[0][0][0][0]=sum1z;}     

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--WWZZ coupl----\n" ); 
   for (int h=0;h<=kkmaxz;h++) {
     fprintf(kkdat , "%s%d%s" , "CPL_W0-W0-Z", h, "-ZX =  ");
     //     fprintf(kkdat , "%s%d%s" , "CPL_W1-W1-Z", h+1, "-ZX =  ");
     for (int l=0;l<=h;l++) {
          fprintf(kkdat , "%.9f%s", 
           (wzfvrt[0][0][h][l]) , "    ");
	  //          fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	  //           (wzfvrt[0][0][h][l]) , "!W" , 1 , "-W" , 1 , "-Z", h+1, "-Z",l+1 );
	  if (l == h) {fprintf(kkdat , "%s%d%s%d\n", "!W0-W0-Z", h, " -Z0 to Z", l);}
	  //	  if (l == h) {fprintf(kkdat , "%s%d%s%d\n", "!W1-W1-Z", h+1, " -Z1 to Z", l+1);}
           }
   }
   for (int h=1;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       for (int j=0;j<=kkmaxz; j++) { 
	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h, "-W", l, "-Z", j, "-ZX =  ");
	 //	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-Z", j+1, "-ZX =  ");
         for (int m=0;m<=j; m++) {
         fprintf(kkdat , "%.9f%s", 
           (wzfvrt[h][l][j][m]) , "    " );
	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h, "-W", l, "-Z", j, "- Z0 to Z", m);}
	  //	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "-Z", j+1, "- Z1 to Z", m+1);}
	 //         fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	 //           (wzfvrt[h][l][j][m]) , "!W" , h+1 , "-W" , l+1 , "-Z", j+1, "-Z",m+1 );
         }
       }
     }
   }

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--WWZ'Z' coupl----\n" ); 
   for (int h=0;h<=kkmaxg;h++) {
     fprintf(kkdat , "%s%d%s" , "CPL_W0-W0-G", h, "-GX =  ");
     //     fprintf(kkdat , "%s%d%s" , "CPL_W1-W1-G", h, "-GX =  ");
     for (int l=0;l<=h;l++) {
          fprintf(kkdat , "%.9f%s", 
           (wgfvrt[0][0][h][l]) , "    ");
	  //          fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	  //           (wgfvrt[0][0][h][l]) , "!W" , 1 , "-W" , 1 , "-Z'", h, "-Z'",l );
	  if (l == h) {fprintf(kkdat , "%s%d%s%d\n", "!W0-W0-Z'", h, " -Z'0 to Z'", l);}
	  //	  if (l == h) {fprintf(kkdat , "%s%d%s%d\n", "!W1-W1-Z'", h, " -Z'0 to Z'", l);}
       }
   }
   for (int h=1;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       for (int j=0;j<=kkmaxg; j++) { 
	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h, "-W", l, "-G", j, "-GX =  ");
	 //	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-G", j, "-GX =  ");
         for (int m=0;m<=j; m++) {      
           fprintf(kkdat , "%.9f%s", 
           (wgfvrt[h][l][j][m]) , "    " );
	   //           fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	   //           (wgfvrt[h][l][j][m]) , "!W" , h+1 , "-W" ,l+1 , "-Z'", j, "-Z'",m );
	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h, "-W", l, "-Z'", j, "- Z'0 to Z'", m);}
	  //	  if (m == j) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "-Z'", j, "- Z'0 to Z'", m);}
         }
       }
     }
   }

   fprintf(kkdat , "--------------------------------------------------\n" );
   fprintf(kkdat , "--WWZ'Z coupl----\n" ); 
   for (int h=0;h<=kkmaxw; h++) {
     for (int l=0;l<=h; l++) {
       for (int j=0;j<=kkmaxg; j++) { 
	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h, "-W", l, "-G", j, "-ZX =  ");
	 //	 fprintf(kkdat , "%s%d%s%d%s%d%s" , "CPL_W", h+1, "-W", l+1, "-G", j, "-ZX =  ");
         for (int m=0;m<=kkmaxz; m++) { 
            fprintf(kkdat , "%.9f%s", 
		    (wgzfvrt[h][l][m][j]) , "    " );
	    //            fprintf(kkdat , "%.9f\t%s%d%s%d%s%d%s%d\n", 
	    //		    (wgzfvrt[h][l][m][j]) , "!W" , h+1 , "-W" ,l+1 , "-Z'", j, "-Z",m+1 );
	  if (m == kkmaxz) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h, "-W", l, "-Z'", j, "- Z1 to Z", m);}
	  //	  if (m == kkmaxz) {fprintf(kkdat , "%s%d%s%d%s%d%s%d\n", "!W", h+1, "-W", l+1, "-Z'", j, "- Z1 to Z", m+1);}
         }
       }
     }
   }

  fclose(kkdat);	
  gsl_integration_workspace_free(w);
  return;
}

