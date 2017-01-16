c  Adapted from ps_qqVVj by Vera Hankele
c
c  Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c  Initial version: November 2012
c  Last modified: November 2012
c
c  phase space generator for pp --> V1 [+ 1 Jet]
c                                  V1 --> l1 l2
c

      subroutine def_ps_qqVj(bos)
          use globalvars, only: ldoblha


      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/utilities/process.inc"
      integer bos, i

c set up variables in common block BLIPSIqqVj for call of LIPSVj

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      COMMON /BLIPSIqqVj/ RM2,RMG,RM2MIN,RM2MAX,
     2                S,YJMIN,YJMAX,PTJMIN

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs


      logical ldebug
      parameter (ldebug=.false.)
      
      if (ldebug) then
         print*, "xm2 from ps_qqvj"
         do i=1,6
            print*, xm2(i)
         enddo
      endif


      s = ecm**2
      
      rm2 = xm2(bos)        ! mass squared
      rmg = xmg(bos)        ! mass * decay width
      rm2min = 1d-3        ! minimal/maximal value of mass squared
!      rm2min = 200**2 ! for comparison with W* -> WH which is "just" a heavy higgs
      rm2max = s/2d0          ! for tan mapping
      if (ldebug) then
         print*, rm2
      endif

c      ptjmin = 0.0d0
      ptjmin = ptmin_def(1)

 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE Vj------------------"
        write(6,*)
           write(6,72) sqrt(rm2),rmg/sqrt(rm2)
           write(6,70) bos,sqrt(rm2min),sqrt(rm2max)
        write(6,*)
        write(6,60) sqrt(s)
        write(6,*)
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end


C**********************  LIPSVj   ****************************************
C
      SUBROUTINE LIPSVj(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)
      IMPLICIT NONE
#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=(1+3)*3, TPI=2D0*PI )
      real*8 RD(*),RN
      real*8 K1(0:3),K2(0:3),V(0:3,max_v,max_kin),X1,X2,W
      real*8 PJET(0:3,N)

C Phase space generator for production of 2 bosons + 1 jet at hadron colliders
C
C  INPUT:   RD        random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C           N         0: no jet
C                     1: one jet      
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4,2)  4-vectors for the vector boson resonances V1 and V2,
C                     the 4th component contains the mass**2
C           V(0:3,3)  4-vectors for the 2 decay leptons, assumed massless
C                     V1(Q(*,1)) --> D(*,1)+D(*,2)
C                     and the Higgs
C           P(0:3)    4-vector for the jet, assumed massless
C           W         the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation


c  Example input for p p -> W -> e ve
c      N = 0, NW = 0 (default)

C  INPUT INFORMATION on the phase space region is transferred via the

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      COMMON /BLIPSIqqVj/ RM2,RMG,RM2MIN,RM2MAX,
     2                S,YJMIN,YJMAX,PTJMIN
C
C   RM2            resonance mass squared for resonance
C   RMG            resonance total width * mass for resonance
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   YJMIN;YJMAX    rapidity range of the individual jets:
C                  YJMIN(I) < ABS(Y(I)) < YJMAX(I)
C   PTJMIN         minimal PT values      REAL*8 Y,PT,PHI,PNJ

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs


C  local variables

      INTEGER I,MU

      real*8 test
      real*8 tau, q(0:4)

      double precision dotrr
      external dotrr

      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)


      logical ldebug
      parameter (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY
      double precision RandomNumber
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RandomNumber

c--------------------------------------------------------------------

      if (ldebug) then
         print*,"--- enter LIPSVj ---"
         print*,"NW = ", NW
      endif


      W = 1d0             !init of weightfactor

      if (NW .eq. 1) then
      ! narrow width: get first masses of vector bosons for lower limit of s-channel
      ! this part for full calculation done afterwards
        q(4) = rm2
        W = W * rmg/2.0d0
      endif

      if (NW .eq. 0) then
         if (ldebug) print*,rm2, rmg, rm2min, rm2max, rd(1), W
         if (ldebug) print*,"do Resonance"
c q^2 of V is generated via tan mapping out of rd(1)
         if (.not.Resonance(rm2, rmg, rm2min, rm2max, 
     1       rd(1), W, q(4))) RETURN 
      endif

c-------------------------------------------------------------------

      if (ldebug) print*,"do twotoone or twotojetplusX"
      if (N.ge.1) then
c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd(4-)
         call TwoToJetsPlusX(N, rd(4), rn, ecm, pTjmin, q(4),
     1        k1(0), k2(0), x1, x2, q(0), pJet, W)
         if (w.eq.0) return
      else if (N.eq.0) then
         call TwoToOne(rd(4),ecm,q(4),k1(0),k2(0),x1,x2,q(0),W)
         if (W.eq.0) return
      endif

      if (ldebug) print*,"do V-decay"
c generate phase space for decay of V1 -> l1 + l2 out of rd(2,3)
      call TwoBodyDecay0(rd(2),rd(3),q(0),q(4),
     1     v(0,1,1),v(0,2,1),W)
 
c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         print*, rd(10), rd(11), rd(12)
c         print*, "y     = ", y
c         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
c         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "------ Momenta: -------"
         

         print*, "ve", v(0,1,1), v(1,1,1), v(2,1,1), v(3,1,1)
         print*, "e", v(0,2,1), v(1,2,1), v(2,2,1), v(3,2,1)
!         print*, "jet ", pjet(0,1), pjet(1,1), pjet(2,1), pjet(3,1)

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-q(mu)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
            test = q(mu)-v(mu,1,1)-v(mu,2,1)
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1,1)-v(mu,2,1)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
         enddo
         print*,"---tau*s = (q)^2 ?---"
         test = tau*s-dotrr(q(0),q(0))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         do i=1,N
           print*,"pjet(",i,")^2 = ",dotrr(pjet(0,i),pjet(0,i))
         enddo
         print*,"q1^2 = ",dotrr(q(0),q(0))
         print*,"v1^2 = m_e^2 = ",dotrr(v(0,1,1),v(0,1,1))
         print*,"v2^2 = m_v^2 = ",dotrr(v(0,2,1),v(0,2,1))
         print*, "-----------------------"

      endif

      return
      end subroutine

