c
c phase space generator for pp --> V1 + V2 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version:  December 2005
c       Last modified: June 2013, Bastian Feigl
c


      subroutine def_ps_qqVVj(bos1,bos2)
          use globalvars, only: ldoblha


      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/utilities/process.inc"
      integer bos(2), i, bos1, bos2

c set up variables in common block BLIPSIqqVVj for call of LIPSVVJ

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      COMMON /BLIPSIqqVVj/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX, RM2MAXX, M2MIN, 
     2                S,YJMIN,YJMAX,PTJMIN

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      INTEGER ihiggs
      DOUBLE PRECISION xm2h, xmgh, rm2minh, rm2maxh, xxm2, xxmwidth, 
     &                 splitrh, splitvv, rm2mid(2)
      COMMON /BLIPSIFLOOPS/ xm2h, xmgh, rm2minh, rm2maxh, rm2mid,
     &                      xxm2, xxmwidth, splitrh, splitvv, ihiggs
      logical is_Z(2),is_W(2)
      common /phspvv/is_Z,is_W
      double precision M2_i(max_res),MG_i(max_res),F_i(max_res)
      integer N_res

      s = ecm**2

      bos(1) = bos1
      bos(2) = bos2
      is_Z(1)=(bos1.eq.2)
      is_Z(2)=(bos2.eq.2)
      is_W(1)=(bos1.eq.3 .or. bos1.eq.4)
      is_W(2)=(bos2.eq.3 .or. bos2.eq.4)

      do i = 1,2
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         rm2min(i) = 1d-3        ! minimal/maximal value of mass squared
         rm2max(i) = s/2d0          ! for tan mapping
      enddo

      rm2minx = 1d-3
      rm2maxx = s                 ! for intermediate particle X -> V1 V2

      if ((floops .eq. 2).and.inFLoops) then 
        if (narrowwidth.eq.1) then
          ihiggs = 2   ! exclusively Higgs peak
c need at least both masses
          rm2minx = (sqrt(rm2(1))+sqrt(rm2(2)))**2
        else
          ihiggs = 1   ! need continuum for gg->H*->VV
        endif
      elseif ((floops .eq. 3).and.inFLoops) then
        if (narrowwidth.eq.1) then
c need at least both masses
          rm2minx = (sqrt(rm2(1))+sqrt(rm2(2)))**2
ccc if Higgs peak on threshold, this is not well-mapped at the moment
          if (xm2(6) .lt. rm2minx) then ! Higgs peak kinematically inaccessible
            ihiggs = 0   ! only continuum
          else
            ihiggs = 1   ! everything
          endif
        else
          ihiggs = 1   ! everything
        endif
      else
        ihiggs = 0   ! only continuum
      endif


      ! select relevant resonances for VV invariant mass
      N_res = 0
      xm2h = xm2(6)
      xmgh = xmg(6)
      if (xm2h.gt.600d0**2 .and. ihiggs.gt.0) ihiggs = 0    ! higgs well mapped by continuum
      xxm2 = (sqrt(xm2(bos1))+sqrt(xm2(bos2)))**2
      xxmwidth = sqrt(xxm2)*500d0
      ! for ZZ: Z->4l resonance
      if (is_Z(1).and.is_Z(2)) then
        N_res = N_res+1
        M2_i(N_res)=xm2(2)
        MG_i(N_res)=xmg(2)
        F_i(N_res)=0.4
      endif
      ! for WZ: W->3l+nu resonance
      if (is_Z(1).and.is_W(2) .or. is_Z(2).and.is_W(1)) then
        N_res = N_res+1
        M2_i(N_res)=xm2(3)
        MG_i(N_res)=xmg(3)
        F_i(N_res)=0.4
      endif
      if (ihiggs .eq. 2) then                      ! only higgs resonance
        call MultiResonancesInit(1, xm2h, xmgh, 1d0, M2_i, MG_i, F_i, 0)
      elseif (ihiggs .eq. 1) then                  ! continuum + higgs (+ Z->4l / W->3l)
        N_res = N_res+1
        M2_i(N_res)=xm2h
        MG_i(N_res)=xmgh
        F_i(N_res)=0.6
        call MultiResonancesInit(1, xxm2, xxmwidth, 1d0, M2_i, MG_i, F_i, N_res)
      else                                         ! only continuum (+ Z->4l / W->3l)
        call MultiResonancesInit(1, xxm2, xxmwidth, 1d0, M2_i, MG_i, F_i, N_res)
      endif


c cutoff for numerical purposes

      m2min = 0.01d0**2

c      ptjmin = 0.0d0
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
           write(6,72) sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
           write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
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


C**********************  LIPSVVJ   ****************************************
C
      SUBROUTINE LIPSVVJ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)
      IMPLICIT NONE
#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=(1+3)*3, TPI=2D0*PI )
      real*8 RD(*),RN
      real*8 K1(0:3),K2(0:3),V(0:3,4),X1,X2,W
      real*8 PJET(0:3,N)

C Phase space generator for production of 2 bosons + 1 jet at hadron colliders
C
C  INPUT:   RD        random number array (between 0 and 1)
C           RN        additional random number which will not be driven by
C                     VEGAS
C
C  OUTPUT:  K1,K2     4-vectors for 2 incoming partons
C           X1,X2     corresponding values of Feynman x
C           Q(0:4,2)  4-vectors for the vector boson resonances V1 and V2,
C                     the 4th component contains the mass**2
C           V(0:3,4)  4-vectors for the 4 decay leptons, assumed massless
C                     V1(Q(*,1)) --> D(*,1)+D(*,2)
C                     V2(Q(*,2)) --> D(*,3)+D(*,4)
C           P(0:3)    4-vector for the jet, assumed massless
C           W         the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation


C  INPUT INFORMATION on the phase space region is transferred via the
C    common block /BLIPSIqqVVj/

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      COMMON /BLIPSIqqVVj/ RM2(2),RMG(2),RM2MIN(2),RM2MAX(2),
     1                RM2MINX, RM2MAXX, M2MIN,
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
      logical is_Z(2),is_W(2)
      common /phspvv/is_Z,is_W


C  OUTPUT INFORMATION is transferred via the common block /DIBOSLIPSO/


C  local variables

      INTEGER I,MU


      real*8 test, rf
      real*8 y, tau, q(0:4,2), y2
      real*8 pX(0:4)


      double precision dotrr
      external dotrr


      logical ldebug
      parameter (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, multiresonance
      double precision RandomNumber
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RandomNumber, multiresonance

      double precision rx

c--------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVV "

      W = 1d0             !init of weightfactor

c pX^2 of intermediate particle X -> V1 + V2 out of rd(1)
c for LO/NLO we only need continuum, if FERMIONLOOP=2 or 3 (also) Higgs peak

      if (NW .eq. 1) then
c narrow width: get first masses of vector bosons for lower limit of s-channel
c this part for full calculation done afterwards
        q(4,1) = rm2(1)
        W = W * rmg(1)/2.0d0
        q(4,2) = rm2(2)
        W = W * rmg(2)/2.0d0
      endif


      ! q^2 of intermediate particle X
      if (.not. MultiResonance(1, rd(1), rm2minx, rm2maxx, W, pX(4))) return


      if (NW .eq. 0) then

        if (RandomNumber() .lt. 0.5d0) then
          i=1
        else 
          i=2
        endif

c q^2 of Vi is generated via tan mapping out of rd(2)
        rm2max(i) = pX(4)
        if (is_Z(i)) then
          if (rm2max(i).gt.rm2(i)-8d0*rmg(i)) then
            if (rd(2).lt.0.2d0) then
                W=W/0.2d0
                if (.not.NonResonance(rm2min(i), rm2(i)-8d0*rmg(i),
     1                               rd(2)/0.2d0, W, q(4,i))) RETURN 
            else
                W=W/0.8d0
                rx = (rd(2)-0.2d0)/0.8d0
                if (rm2max(i).gt.rm2(i)+25d0*rmg(i)) then
                  if (rx.lt.0.7d0) then
                    W=W/0.7d0
                    if (.not.Resonance(rm2(i), rmg(i), rm2(i)-8d0*rmg(i),
     1                     rm2(i)+25d0*rmg(i), rx/0.7d0, W, q(4,i))) RETURN 
                  else
                    W=W/0.3d0
                    if (.not.NonResonance(rm2(i)+25d0*rmg(i), rm2max(i),
     1                               (rx-0.7d0)/0.3d0, W, q(4,i))) RETURN 
                  endif
                else
                  if (.not.Resonance(rm2(i), rmg(i), rm2(i)-8d0*rmg(i), rm2max(i), 
     1                rx, W, q(4,i))) RETURN 
                endif
            endif
          else
            if (.not.NonResonance(rm2min(i), rm2max(i),
     1                               rd(2), W, q(4,i))) RETURN 
          endif
        else
          if (rm2max(i).gt.rm2(i)+25d0*rmg(i)) then
            if (rd(2).lt.0.7d0) then
                W=W/0.7d0
                if (.not.Resonance(rm2(i), rmg(i), rm2min(i),
     1                     rm2(i)+25d0*rmg(i), rd(2)/0.7d0, W, q(4,i))) RETURN 
            else
                W=W/0.3d0
                if (.not.NonResonance(rm2(i)+25d0*rmg(i), rm2max(i),
     1                               (rd(2)-0.7d0)/0.3d0, W, q(4,i))) RETURN 
            endif
          else
            if (.not.Resonance(rm2(i), rmg(i), rm2min(i), rm2max(i), 
     1         rd(2), W, q(4,i))) RETURN 
          endif
        endif

c q^2 of V(3-i) is generated via tan mapping out of rd(3)
        rm2max(3-i) = (sqrt(pX(4))-sqrt(q(4,i)))**2
        ! generate m_V(i) larger than m_V(3-i)
        if (rm2min(1).eq.rm2min(2)) then
           rm2max(3-i) = min( rm2max(3-i), q(4,i) )
           w = w*2d0
        endif
        if (is_Z(3-i)) then
          if (rm2max(3-i).gt.rm2(i)-8d0*rmg(i)) then
            if (rd(3).lt.0.3d0) then
                W=W/0.3d0
                if (.not.NonResonance(rm2min(3-i), rm2(3-i)-8d0*rmg(3-i),
     1                               rd(3)/0.3d0, W, q(4,3-i))) RETURN 
            else
                W=W/0.7d0
                if (.not.Resonance(rm2(3-i), rmg(3-i), rm2(3-i)-8d0*rmg(3-i),
     1                  rm2max(3-i), (rd(3)-0.3d0)/0.7d0, W, q(4,3-i))) RETURN 
            endif
          else
            if (.not.NonResonance(rm2min(3-i), rm2max(3-i),
     1                               rd(3), W, q(4,3-i))) RETURN 
          endif
        else
          if (.not.Resonance(rm2(3-i), rmg(3-i), rm2min(3-i), rm2max(3-i), 
     1        rd(3), W, q(4,3-i))) RETURN 
        endif
      endif

c-------------------------------------------------------------------

      if (N.ge.1) then
c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd(10-12+)

         rf = rd(9)
         call TwoToJetsPlusX(N, rd(10), rn, ecm, pTjmin, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJet, W)
         if (w.eq.0) return

      else if (N.eq.0) then

         rf = rn
         call TwoToOne(rd(9),ecm,pX(4),k1(0),k2(0),x1,x2,pX(0),W)
         if (W.eq.0) return

      endif

c generate phase space for decay of X -> V1 + V2 out of rd(4), rd(rn/9)
      if (.not. TwoBodyDecay(rd(4), rf, pX(0), pX(4), q(4,1), q(4,2), 
     1     q(0,1), q(0,2), W)) return

c generate phase space for decay of V1 -> l1 + l2 out of rd(5-6)
      call TwoBodyDecay0(rd(5),rd(6),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c generate phase space for decay of V2 -> l3 + l4 out of rd(7-8)
      call TwoBodyDecay0(rd(7),rd(8),q(0,2),q(4,2),
     1     v(0,3),v(0,4),W)
      
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
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
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
            test = q(mu,2)-v(mu,3)-v(mu,4)
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)-v(mu,3)-v(mu,4)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*, test
         enddo
         print*,"---tau*s = (q1 + q2)^2 ?---"
         test = tau*s-dotrr(q(0,1),q(0,1))-dotrr(q(0,2),q(0,2))
     &                -2*dotrr(q(0,1),q(0,2))
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
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"

      endif

      return
      end subroutine

