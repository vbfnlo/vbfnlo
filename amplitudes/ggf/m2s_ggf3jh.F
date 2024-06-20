
c	Michael Kubocz
c       Based on 'm2s_ggfh.F' 
c	Initial version:  October 2008 
c
      real*8 function m2s_ggf3jh(
     &                   xi,	!in:  Feynman x parameters of incoming fermions
     &                   p,	!in:  momenta associated with external fermions
     &                   ph,	!in:  momenta of Higgs boson
     &                   rand   !in:  random number array
     &                        )
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "ggf_process.inc"
#include "VBFNLO/utilities/lha.inc"
c... MadGraph check
#include "ggfloMG.inc"
c      logical MGdebug
c      PARAMETER (MGdebug=.true.)
c      common/MGgauge/MGdebug,lgauge2
c      PARAMETER (MGdebug=.false.)
c
c declare input/output variables
c
*******************************************************************
*******************************************************************
c      integer counter(10),gf_iteration
c      logical Singular
c      COMMON/DetCount/Singular,counter,gf_iteration
      integer InterpolCall(4),ICallF(4),SingularF(4)
      COMMON/ICall/InterpolCall,ICallF,SingularF

      integer GCcounter(4)
      logical GaugeFail
      COMMON/GaugeCHeck/GCcounter,GaugeFail

      character*17 wr1,wr2
      CHARACTER*50 wr3,wr4
      Parameter(wr4='(A10,I1,A1,I2,A2,E28.20,E28.20)',
     &          wr1='(E28.20,E28.20)',
     &          wr2='(A10,E28.20)',wr3='(A5,I1,A2,E28.20)')

*******************************************************************
*******************************************************************

      real*8 xi(nx), p(0:3,max_p), ph(0:3), m2s
      real*8 rand(max_PS_dim)

*******************************************************************
*******************************************************************
      logical strange
      COMMON/kinem/strange
*******************************************************************
*******************************************************************

C  variables in the commom blocks    
      REAL*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH,  
     1       CLR,XM2,XMG,B,V,A
      real*8 BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     1       BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM
      real*8 BRFAC,FKF
      integer isubproc, icf

      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM
      COMMON /BCONST/BRFAC,FKF
      COMMON /COUNTER/ isubproc

C  Arguments for call of Matrixelements
      INTEGER FSIGN(4), gsign(5), hsign
      REAL*8 PBAR(0:3,4), KBAR(0:4,5), PHbar(0:4)
      
*----------  Variables for checks (Lorentz check) ----------
      logical lboost
      PARAMETER (lboost=.false.)
c  Boost variables
      real*8 pboost(0:3),pb(0:4,5),phb(0:4)
      integer check,nr
*---------------------------------------------------------------------

c local variables 
      real*8 pdf1(-6:6), pdf2(-6:6), x1p,x2p,w_ps,dsigjj,phi
      integer gpol(5)
      REAL*8 W1,W,RES(11),reslha(9), SHAT, rh(0:2),rh2, srh(0:3), eps,
     1       U1,D1,S1,C1,B1,G1, UB1,DB1,SB1,CB1,BB1,
     2       U2,D2,S2,C2,B2,G2, UB2,DB2,SB2,CB2,BB2,
     3       qqqqs, qqqqt, qqqqid
      DOUBLE PRECISION parton1(-6:6), parton2(-6:6), wsubproc, weight,
     $     rnumb
      complex * 16 amp(3), ci_pent(3),ci_tri(3),ci_box(3),ci_tot(3)
      complex*16 zamp(3,2,2,2,2),iz, samp(3), oamp(3)
      parameter(iz=(0d0,1d0))
      integer i1,i2,i3,i4
      REAL*8 K1(0:3),K2(0:3),PJ(0:3,3),L1(0:3),L2(0:3),QV(0:4)
      REAL*8 QSQDF1,QSQDF2,qfac,POLCOL, WQQ,WQg,Wgg,WTOT
      real*8 xmb, betah, aux
      real*8 rnaux, qqqq0, cth1, sth1, cth2, sth2, rat0, facdec
      real*8 RN8,rwgt 
      double precision rescale

      INTEGER I,J,K,sigma,mu,IRES,I1CALL, icall, fpol(2)

      common/topmass/m
      real * 8 m
      
      real*8 ran2

      REAL*8 alphas5,alphas,RandomNumber
      EXTERNAL alphas5,alphas,RandomNumber

      DATA I1CALL /0/, icall /0/
      save icall, facdec

      INTEGER PERMCALL          !!! for test purposes of the cached C0-function
      DATA PERMCALL/0/
      SAVE PERMCALL

      real*8 maxwgt
      logical  ldebug, lgram, llhadebug, INTERPOL,debugH3j
      debugH3j  = .false.
      ldebug    = .false.
      lgram     = .false.
      llhadebug = .false.
      INTERPOL  = .false.

      Singular=.false.
      
C  Determine the Yukawa coupling Higgs bb from BR(H--->bb) (given in percent)
      if (icall.eq.0) then
         maxwgt = 0
         m = xmt
         eps = 0.005
 
         CALL C_1M_INIT_CACHE() 
         CALL D_1M_INIT_CACHE()
      endif
      icall = icall + 1

      CALL C_1M_CLEAR_CACHE()  
      CALL D_1M_CLEAR_CACHE()

      
c$$$      if (lboost) then
c$$$         pboost(0) =  p(0,1)+p(0,2)
c$$$         pboost(1) =-(p(1,1)+p(1,2))
c$$$         pboost(2) =-(p(2,1)+p(2,2))
c$$$         pboost(3) =-(p(3,1)+p(3,2))
c$$$         do i = 1,5
c$$$            CALL boostn(p(0,i), pboost, pb(0,i))
c$$$c            print*,'p(mu,',i,') boost',(p(mu,i),mu=0,3)
c$$$         enddo
c$$$         pb(4,5) = 0d0
c$$$         CALL boostn(ph,pboost,phb)
c$$$      endif

c$$$      print*,'-------- in m2s_ggf3jh.F --------'
c$$$      write(*,wr2)"p1(0)",P(0,1)
c$$$      write(*,wr2)"p1(1)",P(1,1)
c$$$      write(*,wr2)"p1(2)",P(2,1)
c$$$      write(*,wr2)"p1(3)",P(3,1)
c$$$      print*,''
c$$$      write(*,wr2)"p2(0)",P(0,2)
c$$$      write(*,wr2)"p2(1)",P(1,2)
c$$$      write(*,wr2)"p2(2)",P(2,2)
c$$$      write(*,wr2)"p2(3)",P(3,2)
c$$$      print*,''
c$$$      write(*,wr2)"p3(0)",P(0,3)
c$$$      write(*,wr2)"p3(1)",P(1,3)
c$$$      write(*,wr2)"p3(2)",P(2,3)
c$$$      write(*,wr2)"p3(3)",P(3,3)
c$$$      print*,''
c$$$      write(*,wr2)"p4(0)",P(0,4)
c$$$      write(*,wr2)"p4(1)",P(1,4)
c$$$      write(*,wr2)"p4(2)",P(2,4)
c$$$      write(*,wr2)"p4(3)",P(3,4)
c$$$      print*,''
c$$$      write(*,wr2)"p5(0)",P(0,5)
c$$$      write(*,wr2)"p5(1)",P(1,5)
c$$$      write(*,wr2)"p5(2)",P(2,5)
c$$$      write(*,wr2)"p5(3)",P(3,5)
c$$$      print*,'-----------------'

c$$$      p(0, 1 )=
c$$$      p(1, 1 )=
c$$$      p(2, 1 )=
c$$$      p(3, 1 )=
c$$$      
c$$$      p(0, 2 )=
c$$$      p(1, 2 )=
c$$$      p(2, 2 )=
c$$$      p(3, 2 )=
c$$$      
c$$$      p(0, 3 )=
c$$$      p(1, 3 )=
c$$$      p(2, 3 )=
c$$$      p(3, 3 )=
c$$$      
c$$$      p(0, 4 )= 
c$$$      p(1, 4 )= 
c$$$      p(2, 4 )= 
c$$$      p(3, 4 )= 
c$$$      
c$$$      p(0, 5 )= 
c$$$      p(1, 5 )= 
c$$$      p(2, 5 )= 
c$$$      p(3, 5 )= 

cc start check
c$$$      p(0, 1 )=40.389178910740803d0
c$$$      p(1, 1 )=0.0d0
c$$$      p(2, 1 )=0.0d0
c$$$      p(3, 1 )=40.389178910740803d0
c$$$      
c$$$      p(0, 2 )=36.250754749116794d0
c$$$      p(1, 2 )=28.560504677935327d0
c$$$      p(2, 2 )=15.095858606222116d0
c$$$      p(3, 2 )=16.447791504125952d0
c$$$      
c$$$      p(0, 3 )=2547.9829137573079d0
c$$$      p(1, 3 )=0.0d0
c$$$      p(2, 3 )=0.0d0
c$$$      p(3, 3 )=-2547.9829137573079d0
c$$$      
c$$$      p(0, 4 )=225.07255477451184d0
c$$$      p(1, 4 )= -6.8946602390592293d0
c$$$      p(2, 4 )=28.103403191707393d0 
c$$$      p(3, 4 )=-223.20465340569828d0 
c$$$      
c$$$      p(0, 5 )=160.95599823064785d0
c$$$      p(1, 5 )=-75.920164736055199d0 
c$$$      p(2, 5 )=-29.652751722621897d0 
c$$$      p(3, 5 )=-138.79364635368256d0 
cc end check

c  Define the internal momenta
c  Convert momenta to old ggf-notation
      DO mu=0,3
         K1(mu)=P(mu,1)
         K2(mu)=P(mu,2)
         PJ(mu,1)=P(mu,3)
         PJ(mu,2)=P(mu,4)
         PJ(mu,3)=P(mu,5)
         phbar(mu) = ph(mu)
      ENDDO
      phbar(4) = phbar(0)**2-phbar(1)**2-phbar(2)**2-phbar(3)**2
      hsign = 1

cmk      if(MGdebug) then
cmk         do i = 0,3
cmk            MGpbar(i,1)=k1(i)
cmk            MGpbar(i,2)=k2(i) 
cmk            MGpbar(i,3)=pj(i,1)
cmk            MGpbar(i,4)=pj(i,2)
cmk            MGpbar(i,5)=pj(i,3)
cmk            MGpbar(i,6)=phbar(i)
cmk         enddo
cmkc         do i=1,6
cmkc            print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
cmkc         enddo
cmk      endif


cc start check
c      print*,'in m2s'
c      print*,'PH(mu)',(PH(mu),mu=0,3)
c      print*,'m_Hsq',PH(0)**2-PH(1)**2-PH(2)**2-PH(3)**2
cc end check

c$$$      do i=1,5
c$$$         print*,'P(mu,',i,')',(P(mu,i),mu=0,3)
c$$$      enddo
c$$$      print*,'ph(mu)',(ph(mu),mu=0,3)
c$$$      print*,''

c$$$      write(*,wr2)"phbar(0)",phbar(0)
c$$$      write(*,wr2)"phbar(1)",phbar(1)
c$$$      write(*,wr2)"phbar(2)",phbar(2)
c$$$      write(*,wr2)"phbar(3)",phbar(3)
c$$$      print*,'-----------------'

c calculate factorization scale and running alpha_s
      call calcGGFScales3j(p)

      qfac = sqrt(mufsq(1,1))
      CALL pdfproton( xi(1), qfac, pdf1 )
      qfac = sqrt(mufsq(2,1))
      CALL pdfproton( xi(2), qfac, pdf2 )

      do i=-6,6
         pdf1(i) = pdf1(i)/xi(1)
         pdf2(i) = pdf2(i)/xi(2)
      enddo

      DO i = -5,5
         parton1(i) = pdf1(i*sign1)
         parton2(i) = pdf2(i*sign2)
      ENDDO
      
      g1 = pdf1(0)
      db1 = pdf1(-1*sign1)
      d1 = pdf1(1*sign1)
      ub1 = pdf1(-2*sign1)
      u1 = pdf1(2*sign1)
      s1 = pdf1(3*sign1)
      sb1 = pdf1(-3*sign1)
      c1 = pdf1(4*sign1)
      cb1 = pdf1(-4*sign1)
      b1 = pdf1(5*sign1)
      bb1 = pdf1(-5*sign1)
      
      g2 = pdf2(0)
      db2 = pdf2(-1*sign2)
      d2 = pdf2(1*sign2)
      ub2 = pdf2(-2*sign2)
      u2 = pdf2(2*sign2)
      sb2 = pdf2(-3*sign2)
      s2 = pdf2(3*sign2)
      cb2 = pdf2(-4*sign2)
      c2 = pdf2(4*sign2)
      bb2 = pdf2(-5*sign2)
      b2 = pdf2(5*sign2)
C     
C  Now sum over the subprocesses contributing to Higgs production
C
      RN8=RandomNumber()
              
      if (llhadebug) then
         do i=1,9
            reslha(i)=0d0
         enddo
      endif

      res(1) = 0.0d0
      res(2) = 0.0d0
      res(3) = 0.0d0
      res(4) = 0.0d0
      res(5) = 0.0d0
      res(6) = 0.0d0
      res(7) = 0.0d0
      res(8) = 0.0d0
      res(9) = 0.0d0
      res(10)= 0.0d0
      res(11)= 0.0d0
      WQQ    = 0.0d0
      WQg    = 0.0d0
      Wgg    = 0.0d0

**********************************************************************
********************* QUARK - QUARK  SCATTERING ********************** 
**********************************************************************

*-------------------- q1  q3  ---> q2  q4  g H -----------------------
*-------------------- qb1 qb3 ---> qb2 qb4 g H (by C conjugation) ----
      if (.not.subprQQ) go to 101

      polcol = 1d0 /4d0 /9d0    ! 1/4 for pol., 1/9 for initial colors
      do i = 0,3
         pbar(i,1) = k1(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = k2(i)
         pbar(i,4) = pj(i,2)
         kbar(i,1) = pj(i,3)  
      enddo
      kbar(4,1)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1
      gsign(1) = 1    ! outcoming gluon

      if (lboost) then
         do mu=0,3
            pbar(mu,1) = pb(mu,1)
            pbar(mu,2) = pb(mu,3)
            pbar(mu,3) = pb(mu,2)
            pbar(mu,4) = pb(mu,4)
            kbar(mu,1)   = pb(mu,5) ! outcoming gluon
         enddo    
         print*,''
         print*,'---------------- with BOOST ----------------'
         CALL Higgs3jQQ(PBAR,FSIGN,KBAR(0,1),GSIGN(1),phb,hsign,rh)
      endif

cfc to comment in the call of res(1) and res(2)
cmk    print*, "mguu",mguu
cmk    print*, "mguub",mguub
cmk    print*, "mggu",mggu
cmk    print*,"GGFLOampSq1",GGFLOampSq1

      if(MGdebug) then
         do i = 0,3
            MGpbar(i,1)=pbar(i,1) 
            MGpbar(i,2)=pbar(i,3)
            MGpbar(i,3)=pbar(i,2)
            MGpbar(i,4)=pbar(i,4)
            MGpbar(i,5)=kbar(i,1) !gluon
            MGpbar(i,6)=phbar(i)
         enddo
c         do i=1,6
c            print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
c         enddo
         if(mguu) then
cmk         CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         if(equal_f) then
         CALL Higgs3juuuugH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         else
         CALL Higgs3jucucgH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         endif
         res(1) = ( (u1+d1+s1+c1+b1)*(u2+d2+s2+c2+b2)*rh(1) 
     1           + (u1*u2+d1*d2+s1*s2+c1*c2+b1*b2)*(rh(0)/2-rh(1)) 
     2           ) * polcol
         res(2) = ((ub1+db1+sb1+cb1+bb1)*(ub2+db2+sb2+cb2+bb2)*rh(1)+
     1           (ub1*ub2+db1*db2+sb1*sb2+cb1*cb2+bb1*bb2)*(rh(0)/2-rh(1)) 
     2           ) * polcol
         endif
      else
         CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         res(1) = ( (u1+d1+s1+c1+b1)*(u2+d2+s2+c2+b2)*rh(1) 
     1        + (u1*u2+d1*d2+s1*s2+c1*c2+b1*b2)*(rh(0)/2-rh(1)) 
     2        ) * polcol
         res(2) = ((ub1+db1+sb1+cb1+bb1)*(ub2+db2+sb2+cb2+bb2)*rh(1)+
     1        (ub1*ub2+db1*db2+sb1*sb2+cb1*cb2+bb1*bb2)*(rh(0)/2-rh(1)) 
     2        ) * polcol
cfc            res(1) =( (u1+d1+s1+c1+b1)*(u2+d2+s2+c2+b2)*rh(0) )* polcol
      endif
cfc  Check rh(0),rh(1) vs madgraph at the diamgram by diagram level and 
cfc  cross section. Numericall agreement 10^-10. The differences probably are 
cfc  due to input variables.
cfc  to check it set in ggflo.dat
cfc  mguu =true or false 
cfc  mguu = true      rh(0) checked
cfc  mguu = false      rh(1) checked
cfc  MG_DEBUG =  true 
cfc  MG_PRINT =  true 
cfc  SUBPRQQ      = true      ! quark-quark scattering
cfc  SUBPRQG      = false     ! quark-gluon scattering
cfc  SUBPRGG      = false     ! gluon-gluon scattering
cfc  EFF_TH_CORR = false 
cfc  IF EFF_TH_CORR = true, Then one has to go to mgcheck and define series_t series_p.
cfc one can do this commenting the lines in or out.
cfc  Q_LOOP       = EFF_TH_CORR = false  0 
cfc  For different flavour mguu=false. Only 12 diagrams are not zero.
cfc  Comment out the next line:
cmk      if(MGdebug) then
cmk         if(mguu) CALL mgCheck()
cmk         print*,'mguu, STOP in m2s_ggf3jh.F'
cmk         STOP
cmk      endif
       
cfc  Comparison  with the new MG_5 subroutine uuhuug_MG5 ucucgh also work
cfc  file are : matrixuuuuugH_MG5.F
cfc             matrixucucgH_MG5.F
cfc I had to install a new helas routie jvshxx.F 
cfc This is put in helas directory.
cfc in mgcheck, one can select the old MG2 or MG5 version commenting in or out the lines.
cfc As default now is the MG2 used for same flavour uu> uugh
cfc with MG5, the individual diagrams are not properly re-ordered.
cfc One will need to change this by hand.


*-------------------- q1  qb4 ---> q2  qb3 g H ------------------------
*-------------------- qb1 q4  ---> qb2 q3  g H ------------------------

      do i = 0,3
         pbar(i,1) = k1(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = pj(i,2)
         pbar(i,4) = k2(i)
         kbar(i,1) = pj(i,3)  !gluon
      end do
      kbar(4,1)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1
      gsign(1) = 1  ! outcoming gluon

cmk GF--MG: pbar1--1, pbar2--3, pbar3--2, pbar4--4, kbar1--5
      if(MGdebug) then
         if(mguub) then
            do i = 0,3
               MGpbar(i,1)=pbar(i,1) 
               MGpbar(i,2)=pbar(i,2) !4 (in MG: 2<-->3)
               MGpbar(i,3)=pbar(i,4) !2
               MGpbar(i,4)=pbar(i,3)
               MGpbar(i,5)=kbar(i,1) !gluon
               MGpbar(i,6)=phbar(i)
            enddo
cfc         do i=1,6
cfc            print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
cfc         enddo
cfc         print*
cmk           CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
           if(equal_f) then
           CALL Higgs3juuxuuxgH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
           else
           CALL Higgs3juuxccxgH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
           endif
           res(3) = ( (u1+d1+s1+c1+b1)*(ub2+db2+sb2+cb2+bb2)*rh(1)+(
     1         u1*ub2+d1*db2+s1*sb2+c1*cb2+b1*bb2)*(4*rh(2)+rh(0)-rh(1))
     2          ) * polcol
           res(4) = ( (ub1+db1+sb1+cb1+bb1)*(u2+d2+s2+c2+b2)*rh(1)+(
     1         ub1*u2+db1*d2+sb1*s2+cb1*c2+bb1*b2)*(4*rh(2)+rh(0)-rh(1))
     2          ) * polcol
         endif
      else
         CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         res(3) = ( (u1+d1+s1+c1+b1)*(ub2+db2+sb2+cb2+bb2)*rh(1)+(
     1        u1*ub2+d1*db2+s1*sb2+c1*cb2+b1*bb2)*(4*rh(2)+rh(0)-rh(1))
     2        ) * polcol
         res(4) = ( (ub1+db1+sb1+cb1+bb1)*(u2+d2+s2+c2+b2)*rh(1)+(
     1        ub1*u2+db1*d2+sb1*s2+cb1*c2+bb1*b2)*(4*rh(2)+rh(0)-rh(1))
     2        ) * polcol
      endif
cfc    res(3) = ( (u1)*(ub2)*rh(2)) * polcol
cmk      if(MGdebug) then
cmk         if (mguub)  CALL mgCheck()
cmk         print*,'mguub, STOP in m2s_ggf3jh.F'
cmk         STOP
cmk      endif
cfc


*-------------------- g q1  ---> q2  qb3 q4  H ----------------------
*-------------------- g qb1 ---> qb2 q3  qb4 H (by C conjugation) ---

      polcol = 1d0 /4d0 /24d0 !1/24 for initial colors --> 3 quarks * 8 gluons
                                   
      do i = 0,3
         pbar(i,1) = k2(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = pj(i,2)
         pbar(i,4) = pj(i,3)
         kbar(i,1) = k1(i)   !gluon
      enddo
      kbar(4,1)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = 1
      gsign(1) = -1  ! incoming gluon

      if(MGdebug) then
         if(mggu) then
            do i = 0,3
cmk               MGpbar(i,1)=kbar(i,1) !gluon
cmk               MGpbar(i,2)=pbar(i,3)
cmk               MGpbar(i,3)=pbar(i,4)
cmk               MGpbar(i,4)=pbar(i,2)
cmk               MGpbar(i,5)=pbar(i,1) 
cmk               MGpbar(i,6)=phbar(i)

cmk               MGpbar(i,1)=pbar(i,4)
cmk               MGpbar(i,2)=kbar(i,1) !gluon
cmk               MGpbar(i,3)=pbar(i,3)
cmk               MGpbar(i,4)=pbar(i,2)
cmk               MGpbar(i,5)=pbar(i,1) 
cmk               MGpbar(i,6)=phbar(i)
               
               MGpbar(i,1)=kbar(i,1)
               MGpbar(i,2)=pbar(i,1)
               MGpbar(i,3)=pbar(i,2)
               MGpbar(i,4)=pbar(i,4)
               MGpbar(i,5)=pbar(i,3)
               MGpbar(i,6)=phbar(i)
            enddo
cfc            do i=1,6
cfc               print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
cfc            enddo
cfc            print*
            CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
cfc          if(equal_f) then
cfc          CALL Higgs3jguuuuxH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
cfc          else
cfc          CALL Higgs3jguuccxH(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
cfc          endif
          res(5) = g1*(u2+d2+s2+c2+b2+ub2+db2+sb2+cb2+bb2)
     &           *(4d0*rh(1)+rh(0)/2)*polcol
c     &           *(rh(0))*polcol
         endif
      else
         CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         res(5) = g1*(u2+d2+s2+c2+b2+ub2+db2+sb2+cb2+bb2)
     &        *(4d0*rh(1)+rh(0)/2)*polcol
      endif
      if(MGdebug) then
         if (mggu)    CALL mgCheck()
c         print*,'mggu, STOP in m2s_ggf3jh.F'
c         STOP
      endif

c flipped beams:
c   q1 g --->  q2 qb3 q4  H
c  qb1 g --->  qb2 q3 qb4 H    (via C invariance of amplitude)

      do i = 0,3
         pbar(i,1) = k1(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = pj(i,2)
         pbar(i,4) = pj(i,3)
         kbar(i,1)   = k2(i)   !gluon
      enddo
      kbar(4,1)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = 1
      gsign(1) = -1  ! incoming gluon

      if(MGdebug) then
         if(mggu) then
            do i = 0,3
               MGpbar(i,1)=k1(i) 
               MGpbar(i,2)=k2(i)
               MGpbar(i,5)=pj(i,2)
               MGpbar(i,3)=pj(i,1)
               MGpbar(i,4)=pj(i,3)
               MGpbar(i,6)=phbar(i)
            enddo
cfc            do i=1,6
cfc               print*,'MGpbar',i,'=',(MGpbar(j,i),j=0,3)
cfc            enddo
cfc            print*
            CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
            res(6) = g2*(u1+d1+s1+c1+b1+ub1+db1+sb1+cb1+bb1)
     &           *(4d0*rh(1)+rh(0)/2d0)*polcol
         endif
      else
         CALL Higgs3jQQ(PBAR,FSIGN,KBAR,GSIGN(1),phbar,hsign,rh)
         res(6) = g2*(u1+d1+s1+c1+b1+ub1+db1+sb1+cb1+bb1)
     &        *(4d0*rh(1)+rh(0)/2d0)*polcol
cfc         res(6) = g2*(u1)*rh(0)/2d0*polcol
      endif
      if(MGdebug) then
         if (mggu)    CALL mgCheck()
      endif

cfc To remove the ZERO 

cfc   if(mgdebug) then
         if(mguu) then
            res(3:6)=0d0
         endif
         if(mguub) then
            res(1:2)=0d0
            res(5:6)=0d0
         endif
         if (mggu) then
            res(1:4)=0d0
         endif

         IF (sub_QQ_PART.EQ.1) THEN
            res(5:6)=0d0
         elseif (sub_QQ_PART.EQ.2) THEN
            res(1:4)=0d0
         endif



cfc      else
cfc         res(1)= res(1)
cfc         res(2)= res(2)
cfc         res(3)= res(3)
cfc         res(4)= res(4)
cfc         res(5)= res(5) 
cfc         res(6)= res(6)
cfc      endif
**********************************************************************
********************* QUARK - GLUON SCATTERING  **********************
**********************************************************************

*-------------------- q  g ---> q  g g H -----------------------------
*-------------------- qb g ---> qb g g H -----------------------------

 101  if (.not.subprQG) go to 102

c  first select random fermion pol. for all ggqqH amplitudes

      rwgt = 2d0           !factor 2 is for random fermion helicity summation
      if (rn8.gt.0.5) then
         sigma = 1
      else
         sigma = -1
      endif
c
c  q1 + g ---> q2 + g + g+ H   (q1 = u,d,s,c,b)
c  q1bar + g ---> q2bar + g + g + H  (use C invariance of amplitude)
      do mu = 0,3
         pbar(mu,1) = k2(mu)
         pbar(mu,2) = pj(mu,2)

         kbar(mu,1) = k1(mu) 
         kbar(mu,2) = pj(mu,1)
         kbar(mu,3) = pj(mu,3)
      enddo
      kbar(4,1)   = 0.0d0
      kbar(4,2)   = 0.0d0
      kbar(4,3)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1 
      gsign(1) = -1            ! incoming gluon  
      gsign(2) = 1             ! outgoing gluon
      gsign(3) = 1             ! outgoing gluon

      polcol = rwgt /4d0 /24d0 /2d0 ! 1/2 for 2 identical gluons 
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma       

      CALL Higgs3jQG(PBAR(0,1),FSIGN,FPOL,KBAR(0,1),GSIGN,PHBAR,HSIGN,
     &               RH2)

      res(7)=g1*(u2+d2+s2+c2+b2+ub2+db2+sb2+cb2+bb2)*rh2*polcol

*-------------------- flipped beams ---------------------------------
*-------------------- g q1  ---> g q2  g H (q1 = u,d,s,c,b) ---------
*-------------------- g qb1 ---> g qb2 g H (via C conjugation) ------

      do mu = 0,3
         pbar(mu,1) = k1(mu)
         pbar(mu,2) = pj(mu,1)

         kbar(mu,1) = k2(mu)
         kbar(mu,2) = pj(mu,2)
         kbar(mu,3) = pj(mu,3)
      enddo
      kbar(4,1)   = 0.0d0
      kbar(4,2)   = 0.0d0
      kbar(4,3)   = 0.0d0

      fsign(1) = 1
      fsign(2) = 1
      gsign(1) = -1            ! incoming gluon
      gsign(2) = 1             ! outgoing gluon
      gsign(3) = 1             ! outgoing gluon

      polcol = rwgt /4d0 /24d0 /2d0 ! 1/2 for 2 identical gluons
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL Higgs3jQG(PBAR(0,1),FSIGN,FPOL,KBAR(0,1),GSIGN,PHBAR,HSIGN,RH2)

      res(8) = g2*(u1+d1+s1+c1+b1+ub1+db1+sb1+cb1+bb1)*rh2*polcol
c      print*,'res(8)',rh2*polcol

*-------------------- g g ---> qb1 q2 g H (q1 = u,d,s,c,b) ----------

      do mu = 0,3
         pbar(mu,1) = pj(mu,2)
         pbar(mu,2) = pj(mu,1)

         kbar(mu,1) = k1(mu)
         kbar(mu,2) = k2(mu)
         kbar(mu,3) = pj(mu,3)
      enddo
      fsign(1) = -1
      fsign(2) = 1
      gsign(1) = -1            ! incoming gluons
      gsign(2) = -1            ! incoming gluons
      gsign(3) = 1             ! outgoing gluon

      polcol = rwgt /4d0 /64d0 * 5d0      ! five final state quark flavors
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL Higgs3jQG(PBAR(0,1),FSIGN,FPOL,KBAR(0,1),GSIGN,PHBAR,HSIGN,RH2)

      res(9)=g1*g2*rh2*polcol
      

*-------------------- q1  qb2 ---> g g g H (q1 = u,d,s,c,b) ----------
*-------------------- qb1 q2  ---> g g g H (by C conjugation) -------- 

      do mu = 0,3
         kbar(mu,1) = pj(mu,1)
         kbar(mu,2) = pj(mu,2)
         kbar(mu,3) = pj(mu,3)
         pbar(mu,1) = k1(mu)
         pbar(mu,2) = k2(mu)
      enddo
      fsign(1) = 1
      fsign(2) = -1
      gsign(1) = 1           ! outgoing gluon
      gsign(2) = 1           ! outgoing gluon
      gsign(2) = 1           ! outgoing gluon
c
      polcol = rwgt /4d0 /9d0 /6d0   ! 1/6 for 3 identical gluons 
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL Higgs3jQG(PBAR(0,1),FSIGN,FPOL,KBAR(0,1),GSIGN,PHBAR,HSIGN,RH2)

      res(10)=(u1*ub2+d1*db2+s1*sb2+c1*cb2+b1*bb2+
     &         ub1*u2+db1*d2+sb1*s2+cb1*c2+bb1*b2)*rh2*polcol  

      res(7) =res(7) 
      res(8) =res(8) 
      res(9) =res(9) 
      res(10)=res(10)

      if (sub_QG_PART.EQ.1) THEN
      res(7:8)=0d0
      res(10)=0d0
      elseif (sub_QG_PART.EQ.2) THEN
       res(9:10)=0d0
      elseif (sub_QG_PART.EQ.3) THEN
       res(7:9)=0d0     
      endif 






**********************************************************************
********************* GLUON - GLUON SCATTERING  **********************
**********************************************************************

 102  if (.not.subprGG) go to 103

      do i = 0,3
         kbar(i,1)=k1(i)
         kbar(i,2)=k2(i) 
         kbar(i,3)=pj(i,1)
         kbar(i,4)=pj(i,2)
         kbar(i,5)=pj(i,3)
      enddo

      if(MGdebug) then
         do i = 0,3
            MGpbar(i,1)=k1(i) 
            MGpbar(i,2)=k2(i)
            MGpbar(i,3)=pj(i,1)
            MGpbar(i,4)=pj(i,2)
            MGpbar(i,5)=pj(i,3)
            MGpbar(i,6)=phbar(i)
         enddo
      endif
      
      
      gsign(1) = -1         ! incoming gluons
      gsign(2) = -1         ! incoming gluons 
      gsign(3) = 1          ! outgoing gluons
      gsign(4) = 1          ! outgoing gluons
      gsign(5) = 1          ! outgoing gluons

      rwgt = 32d0  !!!!!!!!!!!!! CHANGED 16 ----> 32 !!!!!!!!!!!!!
      do i = 1,5
         if (rn8.gt.0.5) then
            gpol(i) = 2
            rn8 = 2*rn8 - 1
         else
            gpol(i) = 1
            rn8 = rn8*2
         endif
      enddo
      
      polcol = rwgt /4d0 /64d0 /6d0
      ran2=rand(PS_dimension) ! last entry <=> azimuthal angle


      CALL Higgs3jGG(KBAR,GPOL,GSIGN,RH2)
      res(11)=g1*g2*rh2*polcol

******************************************************************************


* Sum the matrixelements

 103  WQQ = res(1)+res(2)+res(3)+res(4)+res(5)+res(6)
      WQg = res(7)+res(8)+res(9)+res(10)
c      WQg = 0d0
      Wgg = res(11)
c      Wgg = 0d0

cc start check
c$$$      print*,'res(1)',res(1)
c$$$      print*,'res(2)',res(2)
c$$$      print*,'res(3)',res(3)
c$$$      print*,'res(4)',res(4)
c$$$      print*,'res(5)',res(5)
c$$$      print*,'res(6)',res(6)
c$$$      print*,''
c$$$      STOP
cc end check

      if(.not.lgram) then
         ICallF(gf_iteration) = ICallF(gf_iteration) + 1
      endif
      if(.not.Singular) then
         SingularF(gf_iteration) = SingularF(gf_iteration) + 1
      endif

c$$$      if(GaugeFail) then
c$$$         GCcounter(gf_iteration)=GCcounter(gf_iteration)+1
c$$$         WTOT = WQQ+WQg+Wgg
c$$$c         WTOT=0.0d0
c$$$      else
c$$$         WTOT = WQQ+WQg+Wgg
c$$$      endif

      if(singular) then
         counter(gf_iteration)=counter(gf_iteration)+1
c         WTOT = WQQ+WQg+Wgg
c         WTOT=0.0d0
c         print*,"counter(",gf_iteration,")",counter(gf_iteration)
         WTOT = WQQ+WQg+Wgg
      else
         WTOT = WQQ+WQg+Wgg
      endif

* rescale by running alpha_s
      if(resalpha) then
      rescale=(als(1,1)*als(2,1)*als(3,1)*(als(4,1)**2))/(ALFAS**5)
      else
 !!! CHANGE !!!
c      rescale=1d0 !!! for the first time...
      rescale=alphareplace**5/(ALFAS**5)
c      rescale=(als(1,1)*als(2,1)*als(3,1)*(als(4,1)**2))/(ALFAS**5) !!! CHANGE !!!
      endif

      WTOT=WTOT*rescale

      m2s_ggf3jh= WTOT
cc start check
c      print*,'m2s_ggf3jh', m2s_ggf3jh
c      STOP
       if(print_amp) then
         print*, "ms2",WTOT,RH2
      endif
cc end check

      if(MGdebug) CALL mgCheck()



      RETURN
      END

