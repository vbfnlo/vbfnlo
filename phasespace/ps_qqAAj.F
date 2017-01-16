c*****************************************************************************
c
c    begin subroutine DEF_PS_QQAAj
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: February 2011
c       Last modified by Michael Rauch: February 2011
c	Based on DEF_PS_QQVVAj by Giuseppe Bozzi
c
c   initialization of phase space generators for 
c         pp --> A + A + 1 Jet
c
c****************************************************************************

      SUBROUTINE DEF_PS_QQAAj(bos1,bos2)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(2), i, bos1, bos2

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      COMMON /BLIPSIAA/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX, RM2MAXX, M2MIN,
     2                S,YJMIN,YJMAX,PTJMIN

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs

      s = ecm**2

      bos(1) = bos1        ! A = 1
      bos(2) = bos2        ! A = 1

      do i = 1,2
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.1) then
            rm2min(i) = 1.0d-3
            rm2max(i) = 5000.0**2
         else
            rm2min(i) = 1.0d-3      ! minimal/maximal value of mass squared
            rm2max(i) = s/2.        ! for tan mapping
         endif
      enddo

      rm2minx = 0d0
      rm2maxx = s                 ! for intermediate particle X -> V1 V2

      if ((floops .eq. 2).and.inFLoops) then 
        ihiggs = 2   ! exclusively Higgs peak
      elseif ((floops .eq. 3).and.inFLoops) then
        ihiggs = 1   ! everything
      else
        ihiggs = 0   ! only continuum
      endif

      xm2h = xm2(6)
      xmgh = xmg(6)

      xxm2 = 200d0**2
      xxmwidth = sqrt(xxm2)*100d0
      if (xm2h.le.600d0**2) then
        call bwborders(xxm2,xxmwidth,xm2h,xmgh,rm2minh,rm2maxh) 
        splitrh = 0.75d0
      else
c everything well mapped by continuum
        ihiggs = 0
      endif

c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10.0d0)

 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
 80   format ( " configuration ", i3 )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE--------------------"
        write(6,*)
        do i = 1,2
           if (bos(i).eq.1) then
              write(6,72) sqrt(rm2(i)),rmg(i)
              write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
           else
              write(6,72) sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
              write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
           endif
        enddo
        write(6,*)
        write(6,60) sqrt(s)
        write(6,*)
        write(6,80) ihiggs
        write(6,*)
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end

c*****************************************************************************
c
c    end subroutine DEF_PS_QQAAj
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSAAJ
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: February 2011
c       Last modified by Michael Rauch: February 2011
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 2 bosons (+ 1 jet)
c at hadron colliders
c
c  INPUT:   RD         random number array (between 0 and 1)
c           RN         additional random number which will not be driven by
c                      VEGAS
c
c  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
c           X1,X2      corresponding values of Feynman x
c           V(0:3,2)   4-vectors for the photons
c                      V(*,1) = photon
c                      V(*,2) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSAAJ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=10, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3,N)

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTC
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      COMMON /BLIPSIAA/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX, RM2MAXX, M2MIN,
     2                S,YJMIN,YJMAX,PTC

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs

C  local variables

      INTEGER I, MU

      REAL*8 Y, Y2


      REAL*8 test, rf
      REAL*8 tau, q(0:4,3)  

      REAL*8 pY(0:4)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,1)=0d0
      q(4,2)=0d0

c generate q^2_Y from rd(1)
c pY^2 of intermediate particle Y -> A1 + A2 out of rd(1)
c for LO/NLO we only need continuum, if FERMIONLOOP=2 or 3 (also) Higgs peak
      if (ihiggs .eq. 0) then
        if (.not. Resonance(xxm2, xxmwidth, rm2minx, rm2maxx,
     1       rd(1), W, pY(4))) return
      else if (ihiggs .eq. 2) then
        if (.not. Resonance(xm2h, xmgh, rm2minx, rm2maxx,
     1       rd(1), W, pY(4))) return
      else
        if (rd(1) .le. splitrh) then
          if (.not. Resonance(xxm2, xxmwidth, 
     &         rm2minx, rm2maxx-(rm2maxh-rm2minh),
     1         rd(1)/splitrh, W, pY(4))) return
          if (pY(4).gt.rm2minh) pY(4) = pY(4) + (rm2maxh-rm2minh)
          W = W/splitrh
        else 
          if (.not. Resonance(xm2h, xmgh, rm2minh, rm2maxh,
     1         (rd(1)-splitrh)/(1d0-splitrh), W, pY(4))) return
          W = W/(1d0-splitrh)
        endif
      endif


      IF(N.GE.2) THEN ! 2 -> JETs+Y
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(4),rd(5),rd(6)
         rf = rd(3)
         if(.not. TwoToJetsPlusX(N,rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0,1), W)) W=0

         if (W.eq.0) return
      ELSEIF(N.EQ.1) THEN ! 2 -> JET+Y
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(4),rd(5),rd(6)
         rf = rd(3)
         call TwoToJetPlusX(rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0,1), W)

         if (W.eq.0) return

      ELSE ! 2 -> Y

         rf = rn
         call TwoToOne(rd(3),ecm,pY(4),k1(0),k2(0),x1,x2,pY(0),W)
         if (W.eq.0) return

      ENDIF

c generate phase space for decay of Y -> A1 + A2 out of rd(3-rf/4)
      call TwoBodyDecay0(rd(2),rf,pY(0),pY(4),
     1     v(0,1),v(0,2),W)

c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4)
         if (N.eq.1) print*, rd(5), rd(6), rd(7)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-sum(Pjet(mu,1:N))-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-sum(Pjet(mu,1:N))-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?--- "
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - dotrr(q(0,3),q(0,3)) - 2*dotrr(q(0,1),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,3)) - 2*dotrr(q(0,2),q(0,3))
         do i=1,N
         test = test
     &     - dotrr(Pjet(0,i),Pjet(0,i)) - 2*dotrr(Pjet(0,i),q(0,1))
     &     - 2*dotrr(Pjet(0,i),q(0,2)) - 2*dotrr(Pjet(0,i),q(0,3))
         enddo
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         do i=1,N
          print*,"Pjet^2 = ",dotrr(Pjet(0,i),Pjet(0,i))
         enddo
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSAAAJ
c
c*****************************************************************************

