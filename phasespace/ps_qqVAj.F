c*****************************************************************************
c
c    begin subroutine DEF_PS_QQVAj
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: June 2011
c       Last modified by Michael Rauch: June 2012
c	Based on DEF_PS_QQVAAj by Michael Rauch
c
c   initialization of phase space generators for 
c         pp --> V1 + A + 1 Jet
c                V1 --> l1 l2
c
c****************************************************************************

      SUBROUTINE DEF_PS_QQVAj(bos1,bos2)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(2), i, bos1, bos2

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs

      s = ecm**2

      bos(1) = bos1        ! Z = 2, W+ = 3, W- = 4
      bos(2) = bos2        ! A = 1

      do i = 1,2
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.1) then       ! this if-statement can be removed
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3      ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0        ! for tan mapping
         endif
      enddo

      rm2miny = 1d-3
      rm2maxy = s/2d0

      rm2minx = 1d-3
      rm2maxx = s/2d0   ! for intermediate particle X -> V1 V2

      rm2max1 = (sqrt(rm2(1))-20d0)**2
      rm2max2 = (sqrt(rm2(1))+20d0)**2

      if ((floops .eq. 2).and.inFLoops) then 
        ihiggs = 2   ! exclusively Higgs peak
      elseif ((floops .eq. 3).and.inFLoops) then
        ihiggs = 1   ! everything
      else
        ihiggs = 0   ! only continuum
      endif

      xm2h = xm2(6)
      xmgh = xmg(6)

      xxm2 = 500d0**2
      xxmwidth = sqrt(xxm2)*500d0
      if (xm2h.le.600d0**2) then
        call bwborders(xxm2,xxmwidth,xm2h,xmgh,rm2minh,rm2maxh) 
        splitrh = 0.75d0
      else
c everything well mapped by continuum
        ihiggs = 0
      endif

c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10d0)

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
           if (rmg(i).eq.0d0) then
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
c    end subroutine DEF_PS_QQVAj
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAJ1
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: June 2011
c       Last modified by Michael Rauch: June 2011
c       Based on LIPSVAAJ by Michael Rauch
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
c           V(0:3,4)   4-vectors for the 2 decay leptons, assumed massless
c                      Q(*,1) = V(*,1) + V(*,2)
c                      V(*,3) = photon
c                      V(*,4) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation
c
c****************************************************************************

      SUBROUTINE LIPSVAJ1(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,3), X1, X2, W
      REAL*8 Pjet(0:3,N)

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

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

      REAL*8 Y


      REAL*8 test
      REAL*8 tau, q(0:4,3)  

      REAL*8 pY(0:4)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVA(1) "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,2)=0d0


c generate q^2_Y from rd(1)
      if (ihiggs .eq. 0) then
        if (NW .eq. 1) then
c need at least mass of the vector boson
          if (.not. NonResonance(rm2(1), rm2max(1),
     1                    rd(1), W, pY(4))) return
        else
c model flat part plus vector boson peak
          if ( rd(1) .le. 0.1d0 ) then
            if (.not. NonResonance(rm2miny, rm2max1,
     1                       rd(1)/0.1d0, W, pY(4))) return
             W = W/0.1d0
          else 
             if (.not. Resonance(rm2(1), rmg(1), rm2max1, rm2max(1),
     1                       (rd(1)-0.1d0)/0.9d0, W, pY(4))) return
             W = W/0.9d0
          endif
        endif
      else if (ihiggs .eq. 2) then
        if (.not. Resonance(xm2h, xmgh, rm2miny, rm2maxy,
     1       rd(1), W, pY(4))) return
      else
        if (rd(1) .le. splitrh) then
          if (.not. Resonance(xxm2, xxmwidth, 
     &         rm2miny, rm2maxy-(rm2maxh-rm2minh),
     1         rd(1)/splitrh, W, pY(4))) return
          if (pY(4).gt.rm2minh) pY(4) = pY(4) + (rm2maxh-rm2minh)
          W = W/splitrh
        else 
          if (.not. Resonance(xm2h, xmgh, rm2minh, rm2maxh,
     1         (rd(1)-splitrh)/(1d0-splitrh), W, pY(4))) return
          W = W/(1d0-splitrh)
        endif
      endif

      IF(N.GE.1) THEN ! 2 -> JETS+Y
c generate 2 -> 1+N process for k1 k2 -> q_Y Pjet from rd(7+)
         if (.not.TwoToJetsPlusX(N, rd(7), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pjet, W)) return

      ELSE ! 2 -> Y

         call TwoToOne(rd(7),ecm,pY(4),k1(0),k2(0),x1,x2,pY(0),W)
         if (W.eq.0) return

      ENDIF

c q^2 of V1 is generated via tan mapping out of rd(2)
      if (NW .eq. 1) then
c q^2_1 in narrow width approximation
        q(4,1) = rm2(1)
        W = W * rmg(1)/2.0d0
      else
        if (.not.NonResonance(rm2min(1), pY(4), rd(2), W, q(4,1))) RETURN 
      endif

c-------------------------------------------------------------------
c-------------------------------------------------------------------

c generate phase space for decay of Y -> V1 + A out of rd(3-4)
      if (.not. TwoBodyDecay(rd(3), rd(4), pY(0), pY(4), q(4,1), q(4,2),
     1     q(0,1),q(0,2),W)) return

c generate phase space for decay of V1 -> l1 + l2 out of rd(5-6)
      call TwoBodyDecay0(rd(5),rd(6),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c set photons
      do mu=0,3
         v(mu,3)=q(mu,2)
      enddo

c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         if (N.eq.1) print*, rd(10), rd(11), rd(12)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-q(mu,1)-q(mu,2)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)-v(mu,3)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         do i=1,N
           print*,"pjet(",i,")^2 = ",dotrr(pjet(0,i),pjet(0,i))
         enddo
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAJ1
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAJ2
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: June 2011 
c       Last modified by Michael Rauch: June 2011
c       Based on LIPSVAAJ by Michael Rauch
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
c           V(0:3,3)   4-vectors for the 2 decay leptons, assumed massless
c                      Q(*,1) = V(*,1) + V(*,2)
c                      V(*,3) = photon
c           Pjets(0:3,N)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation
c
c****************************************************************************

      SUBROUTINE LIPSVAJ2(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,4), X1, X2, W
      REAL*8 Pjet(0:3,N)
      REAL*8 Pjets(0:3,N+1)

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

C  local variables

      INTEGER I, MU

      REAL*8 Y


      REAL*8 test
      REAL*8 tau, q(0:4,2)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVA(2) "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,2)=0d0


c generate q^2_W from rd(1)
      if (NW .eq. 1) then
c need at least mass of the vector boson
        q(4,1) = rm2(1)
        W = W * rmg(1)/2.0d0
      else
        if (.not. Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1),
     1                      rd(1), W, q(4,1))) return
      endif

c generate 2 -> 1+1+N process for k1 k2 -> Pjets q_W from rd(4-9/12)
       if (.not.TwoToJetsPlusX(N+1, rd(4), 0d0, ecm, pTc, q(4,1),
     1         k1(0), k2(0), x1, x2, q(0,1), pjets(0,1), W)) return

       if (N.ge.1) then
         do i=1,N
           do mu=0,3
             pjet(mu,i) = pjets(mu,i)
           enddo
         enddo
       endif

       do mu=0,3
         q(mu,2) = pjets(mu,N+1)
       enddo

c-------------------------------------------------------------------

c generate phase space for decay of V1 -> l1 + l2 out of rd(2-3)
      call TwoBodyDecay0(rd(2),rd(3),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c set photons
      do mu=0,3
         v(mu,3)=q(mu,2)
      enddo

c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         if (N.eq.1) print*, rd(10), rd(11), rd(12)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-q(mu,1)-q(mu,2)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)-v(mu,3)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,2))
     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
     &     - 2*dotrr(Pjet(0,1),q(0,2)) 
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         do i=1,N
           print*,"pjet(",i,")^2 = ",dotrr(pjet(0,i),pjet(0,i))
         enddo
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAJ2
c
c*****************************************************************************
