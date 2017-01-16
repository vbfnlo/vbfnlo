c...  Gunnar Klaemke 
c...  Based on 'Fgf_QQH' 
c...  Initial version:  2005 July 15

c... modified: Michael Kubocz
c... date    :

      real*8 FUNCTION m2s_ggfh(
     &     xi,     !in:  Feynman x parameters of incoming fermions
     &     p,      !in:  momenta associated with external fermions
     &     ph,     !in:  momenta of Higgs boson
     &     rand    !in:  random number array
     &                        )
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "ggf_process.inc"
#include "VBFNLO/utilities/lha.inc"

c... stability flag for Cijs and Dijs
      logical Singular 
      COMMON/DetCount/Singular

c... declare input/output variables
      real*8 xi(nx), p(0:3,max_p), ph(0:3)
      real*8 rand(max_PS_dim)

c... variables in the commom blocks    
      REAL*8 ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo, 
     &       CLR,XM2,XMG,B,V,A
      real*8 BRFAC,FKF
      integer isubproc

      COMMON /BKOPIN/ ALFAS,XMT,ALFA,XMZ,XMW,SIN2W,XMH, gf, alfas_lo
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      COMMON /BCONST/BRFAC,FKF
      COMMON /COUNTERc/ isubproc

c... arguments for call of Matrixelements
      INTEGER FSIGN(4), gsign(4), hsign
      REAL*8 PBAR(0:3,4), QBAR(0:3,4), PHbar(0:4)
      
c... local variables 
      real*8 pdf1(-6:6), pdf2(-6:6), x1p,x2p
      integer gpol(4)
      REAL*8 RES(9),reslha(9), rh(0:3), srh(0:3), eps,
     &       U1,D1,S1,C1,B1,G1, UB1,DB1,SB1,CB1,BB1,
     &       U2,D2,S2,C2,B2,G2, UB2,DB2,SB2,CB2,BB2,
     &       qqqqs, qqqqt, qqqqid
      DOUBLE PRECISION parton1(-6:6), parton2(-6:6), wsubproc, weight,
     &     rnumb
      complex * 16 amp(3), ci_pent(3),ci_tri(3),ci_box(3)
      complex*16 iz, samp(3), oamp(3)
      parameter(iz=(0d0,1d0))

      REAL*8 K1(0:3),K2(0:3),PJ(0:3,2)
      REAL*8 qfac,POLCOL, WQQ,WQg,Wgg,WTOT


      real*8 RN8,rwgt 
      double precision rescale

      INTEGER I,J,K,sigma,mu,I1CALL, icall, fpol(2)

      common/topmass/m
      real*8 m
      
      real*8 ran2

      REAL*8 RandomNumber
      EXTERNAL RandomNumber

      DATA I1CALL /0/, icall /0/
      save icall

      INTEGER PERMCALL   ! for test purposes of the cached C0-function
      DATA PERMCALL/0/
      SAVE PERMCALL

c... variables for boost test
      real*8 pbarb(0:3,4), pb(0:3), qbarb(0:3,4), phbarb(0:4), qqqqb
      real*8 maxwgt
      logical lboost, lboostg, ldebug, lgram, llhadebug, INTERPOL
      parameter (lboost=.false.,lboostg=.false.,INTERPOL=.false. )
      ldebug    = .false.
      lgram     = .false.
      llhadebug = .false.
c... determine the Yukawa coupling Higgs bb from BR(H--->bb) 
c... (given in percent)
      if (icall.eq.0) then
         maxwgt = 0.0d0
         m = xmt
         eps = 0.005d0

c... initialize cache for c- and d-functions
         CALL C_1M_INIT_CACHE() 
         CALL D_1M_INIT_CACHE()
      endif
      icall = icall + 1

c... clear cache for c- and d-functions
      CALL C_1M_CLEAR_CACHE()  
      CALL D_1M_CLEAR_CACHE()

c... convert momenta to old ggf-notation
      DO mu=0,3
         K1(mu)=P(mu,1)
         K2(mu)=P(mu,2)
         PJ(mu,1)=P(mu,3)
         PJ(mu,2)=P(mu,4)
         phbar(mu) = ph(mu)
      ENDDO

      qfac = dsqrt(mufsq(1,1))
      CALL pdfproton( xi(1), qfac, pdf1 )
      qfac = dsqrt(mufsq(2,1))
      CALL pdfproton( xi(2), qfac, pdf2 )

      do i=-6,6
         pdf1(i) = pdf1(i)/xi(1)
         pdf2(i) = pdf2(i)/xi(2)
      enddo

      DO i = -5,5
         parton1(i) = pdf1(i*sign1)
         parton2(i) = pdf2(i*sign2)
      ENDDO


      if (nfl .ge. 5) then
         b1 = pdf1(5*sign1)
         bb1 = pdf1(-5*sign1)
         bb2 = pdf2(-5*sign2)
         b2 = pdf2(5*sign2)
      else
         b1 = 0d0
         bb1 = 0d0
         bb2 = 0d0
         b2 = 0d0
      end if
      if (nfl .ge. 4) then
         c1 = pdf1(4*sign1)
         cb1 = pdf1(-4*sign1)
         cb2 = pdf2(-4*sign2)
         c2 = pdf2(4*sign2)
      else
         c1 = 0d0
         cb1 = 0d0
         cb2 = 0d0
         c2 = 0d0
      end if
      if (nfl .ge. 3) then
         s1 = pdf1(3*sign1)
         sb1 = pdf1(-3*sign1)
         sb2 = pdf2(-3*sign2)
         s2 = pdf2(3*sign2)
      else
         s1 = 0d0
         sb1 = 0d0
         sb2 = 0d0
         s2 = 0d0
      end if
      if (nfl .ge. 2) then
         ub1 = pdf1(-2*sign1)
         u1 = pdf1(2*sign1)
         ub2 = pdf2(-2*sign2)
         u2 = pdf2(2*sign2)
      else
         ub1 = 0d0
         u1 = 0d0
         ub2 = 0d0
         u2 = 0d0
      end if
      g1 = pdf1(0)
      db1 = pdf1(-1*sign1)
      d1 = pdf1(1*sign1)
      g2 = pdf2(0)
      db2 = pdf2(-1*sign2)
      d2 = pdf2(1*sign2)

c... now sum over the subprocesses contributing to Higgs production
c... first set up the H decay products
      phbar(4) = phbar(0)**2-phbar(1)**2-phbar(2)**2-phbar(3)**2
      hsign = 1
      polcol = 1./4. /9.      ! 1/4 for pol., 1/9 for initial colors

      RN8=RandomNumber()
              
c... initialize counter for subprocesses in fcput
      call fcput(0,0,0,0,0,0d0)

      if (llhadebug) then
         do i=1,9
            reslha(i)=0d0
         enddo
      endif

C*****************  QUARK - QUARK SCATTERING  ***********************

C*******************  q1 q3 ---> q2 q4 H   **************************
C*******************  qb1 qb3 ---> qb2 qb4 H   (by charge conj.)*****
C
      res(1)=0.0d0
      res(2)=0.0d0
      res(3)=0.0d0
      res(4)=0.0d0

      if (.not.subprQQ) go to 101
    
      do i = 0,3
         pbar(i,1) = k1(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = k2(i)
         pbar(i,4) = pj(i,2)
      end do

      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = 1
      fsign(4) = 1

      CALL HiggsQQ(pbar,fsign, phbar,hsign, qqqqt, qqqqs, qqqqid)

      if((lha.or.hepmc))then
         do i=1,5
            do j=1,5
               if(i.ne.j)then
                  wsubproc=parton1(i)*parton2(j)*qqqqt*polcol
                  call fcput(i,j,i,j,1,wsubproc/2d0)
                  call fcput(i,j,j,i,1,wsubproc/2d0)
               else
                  wsubproc=parton1(i)*parton2(j)*qqqqid/2*polcol
                  if((qqqqid+(qqqqs+qqqqt)/36d0).ne.0d0)then
                     call fcput(i,j,i,j,1,wsubproc*
     &                    (qqqqt/2d0-17d0*qqqqs/36d0+qqqqid/2d0)/(qqqqid+(qqqqs
     &                    +qqqqt)/36d0))
                     call fcput(i,j,i,j,2,wsubproc*
     &                    (qqqqs/2d0-17d0*qqqqt/36d0+qqqqid/2d0)/(qqqqid+(qqqqt
     &                    +qqqqs)/36d0))
                  endif
               endif
               if(llhadebug) reslha(1)=reslha(1)+wsubproc
            enddo
         enddo
         
         do i=-5,-1
            do j=-5,-1
               if(i.ne.j)then
                  wsubproc=parton1(i)*parton2(j)*qqqqt*polcol
                  call fcput(i,j,i,j,1,wsubproc/2d0)
                  call fcput(i,j,j,i,1,wsubproc/2d0)
               else
                  wsubproc=parton1(i)*parton2(j)*qqqqid/2*polcol
                  call fcput(i,j,i,j,1,wsubproc*
     &                 (qqqqt/2d0-17d0*qqqqs/36d0+qqqqid/2d0)/(qqqqid+(qqqqs
     &                 +qqqqt)/36d0))
                  call fcput(i,j,i,j,2,wsubproc*
     &                 (qqqqs/2d0-17d0*qqqqt/36d0+qqqqid/2d0)/(qqqqid+(qqqqt
     &                 +qqqqs)/36d0))
               endif
               if(llhadebug) reslha(2)=reslha(2)+wsubproc
            enddo
         enddo
      endif
      res(1)=((u1+d1+s1+c1+b1)*(u2+d2+s2+c2+b2)*qqqqt 
     &      +(u1*u2+d1*d2+s1*s2+c1*c2+b1*b2)*(qqqqid/2d0-qqqqt) 
     &       )*polcol
      res(2)=((ub1+db1+sb1+cb1+bb1)*(ub2+db2+sb2+cb2+bb2)*qqqqt
     &      +(ub1*ub2+db1*db2+sb1*sb2+cb1*cb2+bb1*bb2)*(qqqqid/2d0
     &      -qqqqt))*polcol
      
* --------------------------------------------
*      res(1) = u1*u2*qqqqid/2 * polcol
*      res(2) = 0.0d0
      if (ldebug) then
         print*,' '
         print 71,' t, s for qqqq: ',qqqqt, qqqqs, qqqqid
 71      format(a,3g18.10)
         do i = 0,3
            pbar(i,1) = k1(i)
            pbar(i,4) = pj(i,1)
            pbar(i,3) = k2(i)
            pbar(i,2) = pj(i,2)
         end do
         CALL HiggsQQ(pbar,fsign, phbar,hsign, qqqqt, qqqqs, qqqqid)
         print 71,' flipped qqqq:  ',qqqqs, qqqqt, qqqqid
      endif
C
C*******************  q1 qb4 ---> q2 qb3 H   **********************
C*******************  qb1 q4 ---> qb2 q3 H   **********************
C
      do i = 0,3
         pbar(i,1) = k1(i)
         pbar(i,2) = pj(i,1)
         pbar(i,3) = pj(i,2)
         pbar(i,4) = k2(i)
      end do
      fsign(1) = 1
      fsign(2) = 1
      fsign(3) = -1
      fsign(4) = -1
      CALL HiggsQQ(pbar,fsign, phbar,hsign, qqqqt, qqqqs, qqqqid)

      if((lha.or.hepmc))then
         do i=1,5
            do j=-5,-1
               if(i.ne.-j)then
                  wsubproc=parton1(i)*parton2(j)*qqqqt*polcol                
                  call fcput(i,j,i,j,1,wsubproc/2d0)
                  call fcput(i,j,j,i,1,wsubproc/2d0)
               else
                  wsubproc=parton1(i)*parton2(-i)*(4*qqqqs+qqqqid)
     &                    *polcol
                  call fcput(i,i,i,i,1,wsubproc)
               endif
               if(llhadebug) reslha(3)=reslha(3)+wsubproc
            enddo
         enddo
         
         do i=-5,-1
            do j=1,5
               if(i.ne.-j)then
                  wsubproc=parton1(i)*parton2(j)*qqqqt*polcol
                  call fcput(i,j,i,j,1,wsubproc/2d0)
                  call fcput(i,j,j,i,1,wsubproc/2d0)
               else
                  wsubproc=parton1(i)*parton2(-i)*(4*qqqqs+qqqqid)
     &                    *polcol
                  call fcput(i,i,i,i,1,wsubproc)
               endif
               if(llhadebug) reslha(4)=reslha(4)+wsubproc
            enddo
         enddo
      endif   
      
      res(3)=((u1+d1+s1+c1+b1)*(ub2+db2+sb2+cb2+bb2)*qqqqt+(
     &      u1*ub2+d1*db2+s1*sb2+c1*cb2+b1*bb2)*(4d0*qqqqs+qqqqid-qqqqt)
     &      )*polcol
      res(4)=((ub1+db1+sb1+cb1+bb1)*(u2+d2+s2+c2+b2)*qqqqt+(
     &      ub1*u2+db1*d2+sb1*s2+cb1*c2+bb1*b2)*(4d0*qqqqs+qqqqid-qqqqt)
     &      )*polcol
      
* --------------------------------------------
*      res(3) = 0.0d0
*      res(4) = 0.0d0
      if (ldebug) then
         print71,' t, s for qQqQ: ',qqqqt, qqqqs, qqqqid
         do i = 0,3
            pbar(i,1) = k1(i)
            pbar(i,4) = pj(i,1)
            pbar(i,2) = k2(i)
            pbar(i,3) = pj(i,2)
         end do
         fsign(2) = -1
         fsign(4) = 1
         CALL HiggsQQ(pbar,fsign, phbar,hsign, qqqqt, qqqqs, qqqqid)
         print 71,' flipped qQqQ:  ',qqqqs, qqqqt, qqqqid
      endif

C*****************  QUARK - GLUON SCATTERING  ***********************

c
c---------------  q  g ---> q  g H --------------------------------------
c---------------  qb g ---> qb g H --------------------------------------


 101  res(5)=0.0d0
      res(6)=0.0d0
      res(7)=0.0d0
      res(8)=0.0d0

      if (.not.subprQG) go to 102

c  first select random fermion pol. for all ggqqH amplitudes

      rwgt = 2            !factor 2 is for random fermion helicity summation
      if (rn8.gt.0.5) then
         sigma = 1
      else
         sigma = -1
      endif
c
c  q1 + g ---> q2 + g + H   (q1 = u,d,s,c,b)
c  q1bar + g ---> q2bar + g + H  (use C invariance of amplitude)
      do mu = 0,3
         pbar(mu,1) = k2(mu)
         pbar(mu,2) = pj(mu,2)
         qbar(mu,1) = k1(mu)
         qbar(mu,2) = pj(mu,1)
      enddo
      fsign(1) = 1
      fsign(2) = 1
      gsign(1) = -1            ! incoming gluon
      gsign(2) = 1

      polcol = rwgt/4. /24.
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma       

      CALL HiggsQG(PBAR,FSIGN,FPOL,QBAR,GSIGN,PHBAR,HSIGN,RH)

      if((lha.or.hepmc))then
         if(rh(1)+rh(2).ne.0d0)then
            do i=-5,-1
               do j=1,2
                  wsubproc=parton1(0)*parton2(i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(21,i,21,i,j,wsubproc/2d0)
                  call fcput(21,i,i,21,j,wsubproc/2d0)
                  if(llhadebug) reslha(5)=reslha(5)+wsubproc
               enddo
            enddo
            do i=1,5
               do j=1,2
                  wsubproc=parton1(0)*parton2(i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(21,i,21,i,j,wsubproc/2d0)
                  call fcput(21,i,i,21,j,wsubproc/2d0)
                  if(llhadebug) reslha(5)=reslha(5)+wsubproc
               enddo
            enddo
         endif
      endif
      res(5) = g1*(u2 +d2 +s2 +c2 +b2
     1     +ub2+db2+sb2+cb2+bb2)* rh(0) * polcol
      
* --------------------------------------------
c      res(5) = g1*(u2)* rh(0) * polcol

c flipped beams:
c  g + q1 ---> g + q2 + H   (q1 = u,d,s,c,b)
c  g + q1bar ---> g + q2bar + H    (via C invariance of amplitude)
      do mu = 0,3
         pbar(mu,1) = k1(mu)
         pbar(mu,2) = pj(mu,1)
         qbar(mu,1) = k2(mu)
         qbar(mu,2) = pj(mu,2)
      enddo
      fsign(1) = 1
      fsign(2) = 1
      gsign(1) = -1            ! incoming gluon
      gsign(2) = 1
c
      polcol = rwgt/4. /24.
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL HiggsQG(PBAR,FSIGN,FPOL,QBAR,GSIGN,PHBAR,HSIGN,RH)
      
      if((lha.or.hepmc))then
         if(rh(1)+rh(2).gt.0d0)then
            do i=-5,-1
               do j=1,2
                  wsubproc=parton2(0)*parton1(i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(i,21,i,21,j,wsubproc/2d0)
                  call fcput(i,21,21,i,j,wsubproc/2d0)
                  if(llhadebug) reslha(6)=reslha(6)+wsubproc
               enddo
            enddo
            do i=1,5
               do j=1,2
                  wsubproc=parton2(0)*parton1(i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(i,21,i,21,j,wsubproc/2d0)
                  call fcput(i,21,21,i,j,wsubproc/2d0)
                  if(llhadebug) reslha(6)=reslha(6)+wsubproc
               enddo
            enddo
         endif
      endif
      res(6) = g2*(u1 +d1 +s1 +c1 +b1
     1     +       ub1+db1+sb1+cb1+bb1)* rh(0) * polcol
      
* --------------------------------------------
c      res(6) = g2*(u1)* rh(0) * polcol
c      res(6) = 0.0d0

c 
c  g + g ---> q1bar q2 + H   (q1 = u,d,s,c,b)
      do mu = 0,3
         pbar(mu,1) = pj(mu,2)
         pbar(mu,2) = pj(mu,1)
         qbar(mu,1) = k1(mu)
         qbar(mu,2) = k2(mu)
      enddo
      fsign(1) = -1
      fsign(2) = 1
      gsign(1) = -1            ! incoming gluons
      gsign(2) = -1
c
      polcol = rwgt/4. /64. * 5.      ! five final state quark flavors
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL HiggsQG(PBAR,FSIGN,FPOL,QBAR,GSIGN,PHBAR,HSIGN,RH)

      if((lha.or.hepmc))then
         if(rh(1)+rh(2).gt.0d0)then
            do i=-5,-1
               do j=1,2
                  wsubproc=parton1(0)*parton2(0)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(21,21,i,-i,j,wsubproc/10d0)
                  if(llhadebug) reslha(7)=reslha(7)+wsubproc/10
               enddo
            enddo
            do i=1,5
               do j=1,2
                  wsubproc=parton1(0)*parton2(0)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(21,21,i,-i,j,wsubproc/10)
                  if(llhadebug) reslha(7)=reslha(7)+wsubproc/10
               enddo
            enddo
         endif
      endif
      res(7) = g1*g2* rh(0) * polcol
      
* -------------------------------------------- 
c      res(7)=0.0d0

c 
c  q1 + q2bar ---> g + g + H   (q1 = u,d,s,c,b)
c  qbar1 + q2 ---> g + g + H   by C conjugation 
      do mu = 0,3
         qbar(mu,1) = pj(mu,1)
         qbar(mu,2) = pj(mu,2)
         pbar(mu,1) = k1(mu)
         pbar(mu,2) = k2(mu)
      enddo
      fsign(1) = 1
      fsign(2) = -1
      gsign(1) = 1
      gsign(2) = 1
c
      polcol = rwgt/4. /9. /2.      ! 1/2 for 2 identical gluons
      fpol(1) = fsign(1)*sigma
      fpol(2) = fsign(2)*sigma

      CALL HiggsQG(PBAR,FSIGN,FPOL,QBAR,GSIGN,PHBAR,HSIGN,RH)

      if (lboost .and. rh(0).ne.0) then
         qqqqt = rh(0)
         pb(0) = 10.
         pb(1) = 3.172
         pb(2) = 8.76
         pb(3) = 2.976
         do i = 1,2
            call boostn(pbar(0,i), pb, pbarb(0,i))
            call boostn(qbar(0,i), pb, qbarb(0,i))
         enddo
         call boostn(phbar(0), pb, phbarb(0))
         phbarb(4) = phbar(4)

         call HiggsQG(PBARb,FSIGN,FPOL,QBARb,GSIGN,PHBARb,HSIGN,RH)

         qqqqb = rh(0)
         if(abs(qqqqt/qqqqb-1).gt.1d-12) then
            print*,' before boost: ',qqqqt,' gamma = ',phbarb(0)/phbar(0)
            print*,' after boost:  ',qqqqb,' bef./after ',qqqqt/qqqqb
         endif
      endif


      if((lha.or.hepmc))then
         if(rh(1)+rh(2).gt.0d0)then
            do i=-5,-1
               do j=1,2
                  wsubproc=parton1(i)*parton2(-i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(i,-i,21,21,j,wsubproc)
                  if(llhadebug) reslha(8)=reslha(8)+wsubproc
               enddo
            enddo
            do i=1,5
               do j=1,2
                  wsubproc=parton1(i)*parton2(-i)*rh(j)*rh(0)
     $                 /(rh(1)+rh(2))*polcol
                  call fcput(i,-i,21,21,j,wsubproc)
                  if(llhadebug) reslha(8)=reslha(8)+wsubproc
               enddo
            enddo
         endif
      endif
      res(8) = ( u1*ub2+d1*db2+s1*sb2+c1*cb2+b1*bb2 +
     1        ub1*u2+db1*d2+sb1*s2+cb1*c2+bb1*b2 ) * rh(0) * polcol  
      
* -------------------------------------------- 
c      res(8)=0.0d0


C*****************  GLUON - GLUON SCATTERING  ***********************

c
C-------  Pentagon contributions: gg-->ggH ------------------
c 
c  g + g ---> g +g + H

 102  res(9)=0.0d0

      if (.not.subprGG) go to 103

      do mu = 0,3
         qbar(mu,1) = k1(mu)
         qbar(mu,2) = k2(mu)
         qbar(mu,3) = pj(mu,1)
         qbar(mu,4) = pj(mu,2)
      enddo
      gsign(1) = -1            ! incoming gluons
      gsign(2) = -1
      gsign(3) = 1            ! outgoing gluons
      gsign(4) = 1

      rwgt = 16.
      do i = 1,4
         if (rn8.gt.0.5) then
            gpol(i) = 2
            rn8 = 2*rn8 - 1
         else
            gpol(i) = 1
            rn8 = rn8*2
         endif
      enddo
      
      polcol = rwgt/4. /64. /2.
      ran2=rand(PS_dimension) ! last entry <=> azimuthal angle

      IF(INTERPOL) then
      if (.not.(Heff.or.Aeff)) then
         eps= 0.025d0
         if (abs(rand(PS_dimension)-0.25).lt.eps ) then
            ran2 = 0.25d0 - eps 
            lgram = .true.
         elseif (abs(rand(PS_dimension)-0.75).lt.eps ) then
            ran2 = 0.75d0 - eps 
            lgram = .true.
         endif
       
         if (lgram) then
            CALL C_1M_CLEAR_CACHE()  
            CALL D_1M_CLEAR_CACHE()
            if (ldebug) then
               call HiggsGG(QBAR,GPOL,GSIGN,
     1              rh,amp,ci_tri,ci_box,ci_pent)
               x1p = rh(0)
               oamp(1) = amp(1)+amp(2)+amp(3)
               oamp(2) = amp(1)+amp(2)-2*amp(3)
               oamp(3) = amp(1)-amp(2)
               print*,' ' 
               print*,' rand = ',
     &                rand(PS_dimension),' original momenta: '
               print *,' q1  ', (qbar(j,1),j=0,3)
               print *,' q2  ', (qbar(j,2),j=0,3)
               print *,' q3  ', (qbar(j,3),j=0,3)
               print *,' q4  ', (qbar(j,4),j=0,3)
               print *,' pH  ',phbar
               print *,' gpol',gpol
               print *,'amps : ',(oamp(i),i=1,3),'rh = ',rh(0)
            endif
            
            CALL LIPSB(ran2,
     1           qbar(0,1),qbar(0,2),phbar,qbar(0,3), 
     2           qbarb(0,1),phbarb,qbarb(0,3))

            call HiggsGG(QBARb,GPOL,GSIGN,
     1           rh,amp,ci_tri,ci_box,ci_pent)
            
         else
           
            call HiggsGG(QBAR,GPOL,GSIGN,
     1           rh,amp,ci_tri,ci_box,ci_pent)
            
         endif
         if(lgram) then
            oamp(1) = amp(1)+amp(2)+amp(3)
            oamp(2) = amp(1)+amp(2)-2*amp(3)
            oamp(3) = amp(1)-amp(2)
            do i=1,3
               amp(i) = oamp(i)
            enddo
            ran2 = ran2 + 2d0*eps

            CALL LIPSB(ran2,
     1           qbar(0,1),qbar(0,2),phbar,qbar(0,3), 
     2           qbarb(0,1),phbarb,qbarb(0,3))

            call HiggsGG(QBARb,GPOL,GSIGN,
     1           srh,samp,ci_tri,ci_box,ci_pent)

            oamp(1) = samp(1)+samp(2)+samp(3)
            oamp(2) = samp(1)+samp(2)-2d0*samp(3)
            oamp(3) = samp(1)-samp(2)
            if (ldebug) then
               print *,' q1  ', (qbar(j,1),j=0,3)
               print *,' q1b ',(qbarb(j,1),j=0,3)
               print *,' q2  ', (qbar(j,2),j=0,3)
               print *,' q2b ',(qbarb(j,2),j=0,3)
               print *,' q3  ', (qbar(j,3),j=0,3)
               print *,' q3b ',(qbarb(j,3),j=0,3)
               print *,' q4  ', (qbar(j,4),j=0,3)
               print *,' q4b ',(qbarb(j,4),j=0,3)
               print *,' pH  ',phbar
               print *,' pHb ',phbarb
               print *,gpol,(amp(i),i=1,3),rh(0)
               print *,gpol,(oamp(i),i=1,3),srh(0)
            endif
            do i = 1,3
               amp(i) = (oamp(i)+amp(i))/2d0 + 
     1              (rand(PS_dimension)-ran2+eps)/(2d0*eps)*(oamp(i)-amp(i))
            enddo
            x2p = rh(0)
            rh(0) = ( dreal(amp(1))**2+dimag(amp(1))**2 )*1.6666666666666d0
     1           + ( dreal(amp(2))**2+dimag(amp(2))**2 )*1.5d0
     2           + ( dreal(amp(3))**2+dimag(amp(3))**2 )*4.5d0

            if((lha.or.hepmc)) then
c     Reconstruct amplitudes according to different color factors
               oamp(1)=amp(1)/3+amp(2)/6+amp(3)/2
               oamp(2)=amp(1)/3+amp(2)/6-amp(3)/2
               oamp(3)=amp(1)/3-amp(2)/3
c     Calculate squared amplitudes and split matrix element in different color flows
               rh(1) = (dreal(oamp(1))**2 + dimag(oamp(1))**2)
               rh(2) = (dreal(oamp(2))**2 + dimag(oamp(2))**2)
               rh(3) = (dreal(oamp(3))**2 + dimag(oamp(3))**2)
            endif

            if(ldebug) then 
               print *,gpol,(amp(i),i=1,3),rh(0)
               print*,' true/interpol. (amp) ',x1p/rh(0)
               rh(0) = x2p
               rh(0) = (srh(0)+rh(0))/2d0 + 
     1              (rand(PS_dimension)-ran2+eps)/(2d0*eps)*(srh(0)-rh(0))
               print*,' true/interpol. (rh)  ',x1p/rh(0)
            endif
            ran2 = rand(PS_dimension)
         endif

      else
         
         if (lboostg) then
* when checking lorentz invariance, we have to sum over the polarizations
            qqqqt=0.0d0
            qqqqb=0.0d0
            do i=1,16
*     loop over 16 gluon polarizations
*     define gluon polarization
               j=i
               gpol(1)=Mod(j,2)+1
               j=(j-gpol(1)+1)/2
               gpol(2)=Mod(j,2)+1
               j=(j-gpol(2)+1)/2
               gpol(3)=Mod(j,2)+1
               j=(j-gpol(3)+1)/2
               gpol(4)=Mod(j,2)+1
            
               print *,i,'. polarization :',gpol

               call HiggsGG(QBAR,GPOL,GSIGN,
     1              rh,amp,ci_tri,ci_box,ci_pent)

               qqqqt=qqqqt+rh(0)
 
**** boost
               pb(0) = 10.0d0
               pb(1) = 3.172d0
               pb(2) = 8.76d0
               pb(3) = 2.976d0
               do k = 1,4
                  call boostn(qbar(0,k), pb, qbarb(0,k))
               enddo
         

               call HiggsGG(QBARb,GPOL,GSIGN,
     1              rh,amp,ci_tri,ci_box,ci_pent)
               qqqqb = qqqqb+rh(0)
*     the boosted amp^2
            enddo
            print*,' before boost: ',qqqqt,' gamma = ',qbarb(0,1)/qbar(0,1)
            print*,' after boost:  ',qqqqb,' bef./after ',qqqqt/qqqqb
         else
           
            call HiggsGG(QBAR,GPOL,GSIGN,
     1           rh,amp,ci_tri,ci_box,ci_pent)    
         endif
      endif

      res(9)=g1*g2*rh(0)*polcol
      ELSE
*------------- without interpolation -------------
         
         call HiggsGG(QBAR,GPOL,GSIGN,rh,amp,ci_tri,ci_box,ci_pent)
         res(9)=g1*g2*rh(0)*polcol
*-------------------------------------------------
      ENDIF

      if((lha.or.hepmc))then
         if(rh(1)+rh(2)+rh(3).ne.0d0)then
            do i=1,3
               wsubproc=parton1(0)*parton2(0)*rh(i)*rh(0)
     $              /(rh(1)+rh(2)+rh(3))*polcol
               call fcput(21,21,21,21,2*i-1,wsubproc/2d0)
               call fcput(21,21,21,21,2*i,wsubproc/2d0)
               if(llhadebug) reslha(9)=reslha(9)+wsubproc
            enddo         
         endif
      endif 
c---------------------------------------------------------------------------

* Sum the matrixelements

103   WQQ = res(1) + res(2) + res(3) + res(4)
      WQg = res(5) + res(6) + res(7) + res(8)
      Wgg = res(9)

      if(Singular) then
         WTOT=0.0d0
      else
         WTOT = WQQ+WQg+Wgg
      endif

* rescale by running alpha_s

      rescale=(als(1,1)*als(2,1)*(als(3,1)**2))/(ALFAS**4)
      WTOT=WTOT*rescale

      m2s_ggfh= WTOT

      rnumb=RandomNumber()
c     Les Houches interface - the most probable subprocess 
      if ((lha.or.hepmc)) then
         do i=1,isubproc
            lxwgtup(i)=lxwgtup(i)*rescale
         enddo
         i=0
         weight=0.d0
         do while((i.lt.isubproc).and.(weight.le.rnumb*wtot))
            i=i+1
            weight=weight+lxwgtup(i)
            iprocess=i
         enddo
         
         if(llhadebug)then
            do i=1,9
               if(2*(res(i)-reslha(i))/(res(i)+reslha(i)) .gt. 1d-14
     $              .or.2*(reslha(i)-res(i))/(res(i)+reslha(i)) .gt. 1d-15 )then
                  print*,'difference:'
                  print*,'res(',i,')=       ',res(i)
                  print*,'reslha(',i,')=    ',reslha(i)
                  print*,'rel. error(',i,')=',2*(res(i)-reslha(i))
     $                 /(res(i)+reslha(i))
               endif
            enddo
         endif
      endif
           
      RETURN
      END


********************************************************************************
********************************************************************************

      subroutine fcput(flav1, flav2, flav3, flav4, icf, weight)

      implicit none

c... modified: Michael Kubocz
c... date    : 03.03.2011

#include "VBFNLO/utilities/lha.inc"
#include "ggf_process.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      integer flav1, flav2, flav3, flav4, icf, isubproc,i,j,k,
     $     ind(4),gcount
      common /counterc/ isubproc

      double precision weight
      logical iflavid

c     This subroutine fills the LHA HEPEUP common block
c     with flavor and color information about the incoming/outgoing 
c     partons and the outgoing Higgs.


c     The icf tag is used to determine the desired colorflow.
c     In the processes with external gluons:
c     icf=
c     1+2 means: particle_(1,2,3,4)<->color_(1,2,3,4):
c
c     2 -- 3
c     |    |
c     |    |
c     1 -- 4
c     
c     3+4 means: particle_(1,2,3,4)<->color_(2,1,3,4):
c
c     1 -- 3
c     |    |
c     |    |
c     2 -- 4
c     
c     5+6 means: particle_(1,2,3,4)<->color_(2,3,1,4):
c
c     3 -- 1
c     |    |
c     |    |
c     2 -- 4
c     
c     with odd icf meaning clockwise arrows and even icf meaning counter clockwise arrows.

c     Note that in gggg case, all six labels appear whereas in ggqq case,
c     there are only two cases, because there is one color line less than
c     in the gggg case and the orientation of the color flow is given by the
c     fermion line.
c     In the qqqq process with four equal quarks, there exist two different 
c     color structures, resulting by Pauli interference.

c     Initialize process counter
      if(icf.eq.0)then
         isubproc=0
         goto 21101
      endif

c     Count subprocesses:
      isubproc=isubproc+1

c     number of particles in the event
      SELECT CASE(process)
      CASE(GFHjj)
         lnup(isubproc)=5
      CASE(GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar)
         lnup(isubproc) = 7
      CASE(GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
         lnup(isubproc) = 10
      END SELECT

c     put weight of subprocess
      lxwgtup(isubproc)=weight

c     put particle ID of q/g
      lidup(1,isubproc)=flav1
      lidup(2,isubproc)=flav2
      lidup(3,isubproc)=flav3
      lidup(4,isubproc)=flav4

c     put particle ID of Higgs
      SELECT CASE(process)
      CASE(GFHjj,GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar)
         if(model.eq.1) then
            lidup(5,isubproc)=25 ! SM
         endif
         if(model.eq.2) then
            if(HiggsType.eq.1) then
               lidup(5,isubproc)=25 ! h
            endif
            if(HiggsType.eq.2) then
               lidup(5,isubproc)=35 ! H
            endif
            if(HiggsType.eq.3) then
               lidup(5,isubproc)=36 ! A
            endif
         endif
      END SELECT

c     put status code
      listup(1,isubproc)=-1
      listup(2,isubproc)=-1
      listup(3,isubproc)=1
      listup(4,isubproc)=1
      listup(5,isubproc)=1

c     put mothers
      do i=1,2
         do j=3,5
            lmothup(i,j,isubproc)=i
         enddo
      enddo

c     decide whether this is a qqqq-, qqgg- or a gggg-process
      gcount=0
      j=1
      do i=1,4
         if (lidup(i,isubproc).eq.21) then
            gcount=gcount+1
            ind(j)=i
            j=j+1
         endif
      enddo 

c     Higgs wears no color:
      licolup(1,5,isubproc)=0
      licolup(2,5,isubproc)=0

c     For four gluons, we can easily fill the tags:
      if(gcount.eq.4)then
         if(icf.eq.1)then
            licolup(1,1,isubproc)=501
            licolup(2,1,isubproc)=504
            licolup(1,2,isubproc)=502
            licolup(2,2,isubproc)=501
            licolup(1,3,isubproc)=502
            licolup(2,3,isubproc)=503
            licolup(1,4,isubproc)=503
            licolup(2,4,isubproc)=504
         elseif(icf.eq.2)then
            licolup(2,1,isubproc)=501
            licolup(1,1,isubproc)=504
            licolup(2,2,isubproc)=502
            licolup(1,2,isubproc)=501
            licolup(2,3,isubproc)=502
            licolup(1,3,isubproc)=503
            licolup(2,4,isubproc)=503
            licolup(1,4,isubproc)=504
         elseif(icf.eq.3)then
            licolup(1,1,isubproc)=502
            licolup(2,1,isubproc)=501
            licolup(1,2,isubproc)=501
            licolup(2,2,isubproc)=504
            licolup(1,3,isubproc)=502
            licolup(2,3,isubproc)=503
            licolup(1,4,isubproc)=503
            licolup(2,4,isubproc)=504
         elseif(icf.eq.4)then
            licolup(2,1,isubproc)=502
            licolup(1,1,isubproc)=501
            licolup(2,2,isubproc)=501
            licolup(1,2,isubproc)=504
            licolup(2,3,isubproc)=502
            licolup(1,3,isubproc)=503
            licolup(2,4,isubproc)=503
            licolup(1,4,isubproc)=504
         elseif(icf.eq.5)then
            licolup(1,1,isubproc)=502
            licolup(2,1,isubproc)=501
            licolup(1,2,isubproc)=503
            licolup(2,2,isubproc)=502
            licolup(1,3,isubproc)=504
            licolup(2,3,isubproc)=501
            licolup(1,4,isubproc)=503
            licolup(2,4,isubproc)=504
         elseif(icf.eq.6)then
            licolup(2,1,isubproc)=502
            licolup(1,1,isubproc)=501
            licolup(2,2,isubproc)=503
            licolup(1,2,isubproc)=502
            licolup(2,3,isubproc)=504
            licolup(1,3,isubproc)=501
            licolup(2,4,isubproc)=503
            licolup(1,4,isubproc)=504
         endif


c     This is the two-gluon case:
      elseif(gcount.eq.2)then
         do i=1,4
            if(lidup(i,isubproc).ne.21) then
               if(lidup(i,isubproc).lt.0) then   !if the particle is an antiquark,
                  licolup(1,i,isubproc)=0        !then it carries no color
                  licolup(2,i,isubproc)=500+i    !but anticolor
               else
                  licolup(1,i,isubproc)=500+i    !else vice versa
                  licolup(2,i,isubproc)=0
               endif
c     To make sure that no errors occur due to remains of
c     an earlier call of this subroutine:
            else
               licolup(1,i,isubproc)=0
               licolup(2,i,isubproc)=0
            endif
         enddo
         do i=1,4
            if(lidup(i,isubproc).ne.21)then !find q or qbar
               if(icf.eq.1)then
                  if(lidup(i,isubproc)*listup(i,isubproc).lt.0)then   !if incoming q or outgoing qbar
                     if(listup(ind(1),isubproc)*
     $                    listup(i,isubproc).lt.0)then   !if first g has another status, then
                        licolup(1,ind(1),isubproc)=      !connect its color and anticolor lines
     $                       licolup(1,i,isubproc)       !with the color and anticolor lines of the 
                        licolup(2,ind(1),isubproc)=      !fermion. One of the two tags is zero.
     $                       licolup(2,i,isubproc)
                     elseif(listup(ind(1),isubproc)*
     $                       listup(i,isubproc).gt.0)then !if the first gluon has the same status,
                        licolup(1,ind(1),isubproc)=       !then its color and anticolor lines must
     $                       licolup(2,i,isubproc)        !be connected with the anticolor and color
                        licolup(2,ind(1),isubproc)=       !lines of the fermion, i.e. the gluon carries
     $                       licolup(1,i,isubproc)        !the swapped tags of the fermion
                     endif
                  elseif(lidup(i,isubproc)*listup(i,isubproc).gt.0)then    !if outgoing q or incoming qbar.
                     if(listup(ind(2),isubproc)*                           !Note that for each diagram, both
     $                    listup(i,isubproc).lt.0)then                     !conditions are fulfilled once.
                        licolup(1,ind(2),isubproc)=
     $                       licolup(1,i,isubproc)
                        licolup(2,ind(2),isubproc)=    !Same as above, but with the other gluon
     $                       licolup(2,i,isubproc)
                     elseif(listup(ind(2),isubproc)*
     $                       listup(i,isubproc).gt.0)then
                        licolup(1,ind(2),isubproc)=
     $                       licolup(2,i,isubproc)
                        licolup(2,ind(2),isubproc)=
     $                       licolup(1,i,isubproc)
                     endif
                  endif
               elseif(icf.eq.2)then   !Second possible colorflow, which means exchanged gluons
                  if(lidup(i,isubproc)*listup(i,isubproc).lt.0)then
                     if(listup(ind(2),isubproc)*
     $                    listup(i,isubproc).lt.0)then
                        licolup(1,ind(2),isubproc)=
     $                       licolup(1,i,isubproc)
                        licolup(2,ind(2),isubproc)=
     $                       licolup(2,i,isubproc)
                     elseif(listup(ind(2),isubproc)*
     $                       listup(i,isubproc).gt.0)then
                        licolup(1,ind(2),isubproc)=
     $                       licolup(2,i,isubproc)
                        licolup(2,ind(2),isubproc)=
     $                       licolup(1,i,isubproc)
                     endif
                  elseif(lidup(i,isubproc)*listup(i,isubproc).gt.0)then
                     if(listup(ind(1),isubproc)*
     $                    listup(i,isubproc).lt.0)then
                        licolup(1,ind(1),isubproc)=
     $                       licolup(1,i,isubproc)
                        licolup(2,ind(1),isubproc)=
     $                       licolup(2,i,isubproc)
                     elseif(listup(ind(1),isubproc)*
     $                       listup(i,isubproc).gt.0)then
                        licolup(1,ind(1),isubproc)=
     $                       licolup(2,i,isubproc)
                        licolup(2,ind(1),isubproc)=
     $                       licolup(1,i,isubproc)
                     endif
                  endif               
               endif
            endif
         enddo
c     After connecting with fermions, two gluon color tags are still zero.
c     The following loop sets them to 505.
         do i=1,4
            if(lidup(i,isubproc).eq.21)then
               do j=1,2
                  if(licolup(j,i,isubproc).eq.0)then
                     licolup(j,i,isubproc)=505
                  endif
               enddo
            endif
         enddo


c     The qqqq-process:
      else
         iflavid=.true.
         i=1
         k=1
         do while(i.le.4)
            j=i+1
            do while(j.le.4)
               if(lidup(i,isubproc).ne.lidup(j,isubproc))then
                  iflavid=.false.
c     If two quarks with identical flavors are found, their indices are
c     stored in ind(1) and ind(3) or ind(2) and ind(4)
               else
                  if (k.le.2) then
                    ind(k)=i
                    ind(k+2)=j
                  endif
                  k=k+1
               endif
               j=j+1
            enddo
            i=i+1
         enddo
c     Four identical flavors, which is accompanied by two different colorflows:
         if(iflavid)then
            if(lidup(1,isubproc).gt.0)then
               if(icf.eq.1)then
                  licolup(1,1,isubproc)=501
                  licolup(2,1,isubproc)=0
                  licolup(1,2,isubproc)=502
                  licolup(2,2,isubproc)=0
                  licolup(1,3,isubproc)=501
                  licolup(2,3,isubproc)=0
                  licolup(1,4,isubproc)=502
                  licolup(2,4,isubproc)=0
               elseif(icf.eq.2)then
                  licolup(1,1,isubproc)=501
                  licolup(2,1,isubproc)=0
                  licolup(1,2,isubproc)=502
                  licolup(2,2,isubproc)=0
                  licolup(1,3,isubproc)=502
                  licolup(2,3,isubproc)=0
                  licolup(1,4,isubproc)=501
                  licolup(2,4,isubproc)=0
               endif
            elseif(lidup(1,isubproc).lt.0)then
               if(icf.eq.1)then
                  licolup(1,1,isubproc)=0
                  licolup(2,1,isubproc)=501
                  licolup(1,2,isubproc)=0
                  licolup(2,2,isubproc)=502
                  licolup(1,3,isubproc)=0
                  licolup(2,3,isubproc)=501
                  licolup(1,4,isubproc)=0
                  licolup(2,4,isubproc)=502
               elseif(icf.eq.2)then
                  licolup(1,1,isubproc)=0
                  licolup(2,1,isubproc)=501
                  licolup(1,2,isubproc)=0
                  licolup(2,2,isubproc)=502
                  licolup(1,3,isubproc)=0
                  licolup(2,3,isubproc)=502
                  licolup(1,4,isubproc)=0
                  licolup(2,4,isubproc)=501
               endif
            endif
         else
            k=0
c     Set equal color tags to quarks with different flavors, which are
c     identified by ind(i) array.
            if(listup(ind(1),isubproc)*listup(ind(2),isubproc)*
     $           lidup(ind(1),isubproc)*lidup(ind(2),isubproc).lt.0)then
               ! connect 1 and 2, 3 and 4
               if(lidup(ind(1),isubproc).gt.0)then
                   licolup(1,ind(1),isubproc)=501
                   licolup(2,ind(1),isubproc)=0
                   if(lidup(ind(2),isubproc).gt.0)then
                      licolup(1,ind(2),isubproc)=501
                      licolup(2,ind(2),isubproc)=0
                   else
                      licolup(1,ind(2),isubproc)=0
                      licolup(2,ind(2),isubproc)=501
                   endif
                else
                   licolup(1,ind(1),isubproc)=0
                   licolup(2,ind(1),isubproc)=501
                   if(lidup(ind(2),isubproc).gt.0)then
                      licolup(1,ind(2),isubproc)=501
                      licolup(2,ind(2),isubproc)=0
                   else
                      licolup(1,ind(2),isubproc)=0
                      licolup(2,ind(2),isubproc)=501
                   endif
                endif
                if(lidup(ind(3),isubproc).gt.0)then
                   licolup(1,ind(3),isubproc)=502
                   licolup(2,ind(3),isubproc)=0
                   if(lidup(ind(4),isubproc).gt.0)then
                      licolup(1,ind(4),isubproc)=502
                      licolup(2,ind(4),isubproc)=0
                   else
                      licolup(1,ind(4),isubproc)=0
                      licolup(2,ind(4),isubproc)=502
                   endif
                else
                   licolup(1,ind(3),isubproc)=0
                   licolup(2,ind(3),isubproc)=502
                   if(lidup(ind(4),isubproc).gt.0)then
                      licolup(1,ind(4),isubproc)=502
                      licolup(2,ind(4),isubproc)=0
                   else
                      licolup(1,ind(4),isubproc)=0
                      licolup(2,ind(4),isubproc)=502
                   endif
                endif
            elseif(listup(ind(1),isubproc)*listup(ind(4),isubproc)*
     $              lidup(ind(1),isubproc)*lidup(ind(4),isubproc).lt.0)then
                                ! connect 1 and 4, 3 and 2
               if(lidup(ind(1),isubproc)*listup(ind(1),isubproc).lt.0)then
                   licolup(1,ind(1),isubproc)=501
                   licolup(2,ind(1),isubproc)=0
                   if(lidup(ind(4),isubproc).gt.0)then
                      licolup(1,ind(4),isubproc)=501
                      licolup(2,ind(4),isubproc)=0
                   else
                      licolup(1,ind(4),isubproc)=0
                      licolup(2,ind(4),isubproc)=501
                   endif
                else
                   licolup(1,ind(1),isubproc)=0
                   licolup(2,ind(1),isubproc)=501
                   if(lidup(ind(4),isubproc).gt.0)then
                      licolup(1,ind(4),isubproc)=501
                      licolup(2,ind(4),isubproc)=0
                   else
                      licolup(1,ind(4),isubproc)=0
                      licolup(2,ind(4),isubproc)=501
                   endif
                endif
                if(lidup(ind(3),isubproc).gt.0)then
                   licolup(1,ind(3),isubproc)=502
                   licolup(2,ind(3),isubproc)=0
                   if(lidup(ind(2),isubproc).gt.0)then
                      licolup(1,ind(2),isubproc)=502
                      licolup(2,ind(2),isubproc)=0
                   else
                      licolup(1,ind(2),isubproc)=0
                      licolup(2,ind(2),isubproc)=502
                   endif
                else
                   licolup(1,ind(3),isubproc)=0
                   licolup(2,ind(3),isubproc)=502
                   if(lidup(ind(2),isubproc).gt.0)then
                      licolup(1,ind(2),isubproc)=502
                      licolup(2,ind(2),isubproc)=0
                   else
                      licolup(1,ind(2),isubproc)=0
                      licolup(2,ind(2),isubproc)=502
                   endif
                endif
            endif
         endif
      endif

21101 return
      end



