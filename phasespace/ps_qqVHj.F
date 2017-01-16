c  Adapted from ps_qqVVj by Vera Hankele
c
c  Robin Roth <robin@particle.physik.uni-karlsruhe.de>
c  Initial version: November 2012
c  Last modified: November 2012
c
c  phase space generator for pp --> V1 + H [+ 1 Jet]
c                                   V1 --> l1 l2, H stable
c

      subroutine def_ps_qqVHj(bos1,bos2,bos21,bos22)
          use globalvars, only: ldoblha


      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
#include "VBFNLO/utilities/process.inc"
      integer bos(4), i, bos1, bos2, bos21, bos22, numbos

c set up variables in common block BLIPSIqqVHj for call of LIPSVHj

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      COMMON /BLIPSIqqVHj/ RM2(4),RMG(4),RM2MIN(4),RM2MAX(4),
     1                RM2MINX, RM2MAXX, M2MIN, 
     2                S,YJMIN,YJMAX,PTJMIN,numbos

      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

      DOUBLE PRECISION xxm2, xxmwidth 
      COMMON /BLIPSIqqVHj2/ xxm2, xxmwidth


      logical ldebug
      parameter (ldebug=.false.)
      
      if (ldebug) then
         print*, "xm2 from ps_qqvhj"
         do i=1,6
            print*, xm2(i)
         enddo
      endif


      s = ecm**2
      
      bos(1) = bos1
      bos(2) = bos2
      bos(3) = bos21 
      bos(4) = bos22
      if (bos21 .ge. 2 .and. bos21 .le. 4) then
         numbos = 4
      else
         numbos = 2
      endif

      do i = 1,numbos
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         rm2min(i) = 1d-3        ! minimal/maximal value of mass squared
         rm2max(i) = s/2d0          ! for tan mapping
      enddo
      if (ldebug) then
         if (bos(2) .ne. 6) print*, 'bos2 should be higgs here'
         do i=1,2
            print*, rm2(i)
         enddo
      endif

      !max,min: allow for large higgsmass and large higgswidth to still work
      rm2min(2) = (dmax1(rm2(2)-1000*rmg(2),5d0  )/sqrt(rm2(2)))**2 
      rm2max(2) = (dmin1(rm2(2)+1000*rmg(2),s/2)/sqrt(rm2(2)))**2 

      ! debugging for broken large higgsmasses
!      print*, rm2(2), rmg(2)
!      print*, 'maxmin', rm2min(2), rm2max(2)

!      rm2minx = 1d-3
      !rm2minx gets set on a phasespacedependent basis, lower limit also in lipsvhj
      rm2maxx = s/2d0                 ! for intermediate particle X -> V1 V2

!      xxm2 = 500d0**2!(sqrt(rm2(1))+sqrt(rm2(2)))**2
!      xxmwidth = sqrt(xxm2)*500d0 !xxm2 

!      call bwborders(xxm2,xxmwidth,xm2h,xmg,rm2min(2),rm2max(2)) 

      ! tweaked version:
      xxm2 = (sqrt(rm2(1))+sqrt(rm2(2)))**2 ! approx 100d0**2
      xxmwidth = xxm2 

c      ptjmin = 0.0d0
      ptjmin = ptmin_def(1)

 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE VHj------------------"
        write(6,*)
        do i = 1,numbos
           write(6,72) sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
           write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
        enddo
        write(6,*)
        write(6,60) sqrt(s)
        write(6,*)
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end


C**********************  LIPSVHj   ****************************************
C
      SUBROUTINE LIPSVHj(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW,NWH)
      IMPLICIT NONE
#include "VBFNLO/utilities/global.inc"
c
      INTEGER N,NRANDL,NW,NWH
      REAL*8 TPI
      PARAMETER (nrandl=(1+3)*3, TPI=2D0*PI )
      real*8 RD(*),RN
      real*8 K1(0:3),K2(0:3),V(0:3,max_v,max_kin),X1,X2,W
      real*8 PJET(0:3,N)

C Phase space generator for production of V + H + 1 jet at hadron colliders
C Higgs decays either to 2 fermions or to 4 fermions (via 2 bosons)
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
C           V(0:3,4)  4-vectors for the 2 decay leptons, assumed massless
C                     V1(Q(*,1)) --> D(*,1)+D(*,2)
C                     and the Higgs decay products
C           P(0:3)    4-vector for the jet, assumed massless
C           W         the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation for the W,
c                      NW = 0 means Full calculation
c           NWH        like NW, but for Higgs


c  Example input for p p -> W -> W H -> e ve H
c      N = 0, NW = 0 (default)


C  INPUT INFORMATION on the phase space region is transferred via the

      REAL*8 RM2,RMG,RM2MIN,RM2MAX, S,YJMIN,YJMAX,PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN
      integer numbos
      COMMON /BLIPSIqqVHj/ RM2(4),RMG(4),RM2MIN(4),RM2MAX(4),
     1                RM2MINX, RM2MAXX, M2MIN,
     2                S,YJMIN,YJMAX,PTJMIN,numbos
C
C   RM2            resonance mass squared for resonance
C   RMG            resonance total width * mass for resonance
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   YJMIN;YJMAX    rapidity range of the individual jets:
C                  YJMIN(I) < ABS(Y(I)) < YJMAX(I)
C   PTJMIN         minimal PT values      REAL*8 Y,PT,PHI,PNJ

      DOUBLE PRECISION xxm2, xxmwidth 
      COMMON /BLIPSIqqVHj2/ xxm2, xxmwidth


C  local variables

      INTEGER I,MU,rdoffset,rdoffset2

      real*8 test, rf
      real*8 q(0:4,4) !numbos
      real*8 pX(0:4)
      real*8 rm2max1, rm2maxw2
      double precision dotrr
      external dotrr

      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)


      logical lprintmomenta
      parameter (lprintmomenta=.false.)
      logical ldebug
      parameter (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY
      double precision RandomNumber
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RandomNumber

c--------------------------------------------------------------------

      if (ldebug) then
         print*,"--- enter LIPSVHj ---"
         print*,"NW = ", NW
         print*,"NWH = ", NWH
      endif

      ! some initialisation
      do i=1,4 !numbos
         do mu=0,3
            q(mu,i)=0d0
         enddo
      enddo
      do i=1,n_v
         do mu=0,3
            v(mu,i,1)=0d0
         enddo
      enddo

      W = 1d0             !init of weightfactor

c pX^2 of intermediate particle X -> V1 + V2 out of rd(1)


      rm2minx = 0d0

      if (NW .eq. 1) then
         ! narrow width: get first masses of vector bosons for lower limit of s-channel
         ! this part for full calculation done afterwards
         q(4,1) = rm2(1)
         rm2minx = rm2minx + sqrt(q(4,1))
         W = W * rmg(1)/2.0d0
      endif
      if (NWH .eq. 1) then
         ! Higgs is on-shell q^2 fixed:
         q(4,2) = rm2(2)
         rm2minx = rm2minx + sqrt(q(4,2))
         W = W * rmg(2)/2.0d0 ! this only if higgs decaying 
      endif

      ! NW=0: rm2minx = 0
      ! NW=1 and/or NWH=1: rm2minx sum of masses
      ! need at least both masses
      rm2minx = max(1d-3,rm2minx**2)
      !rm2minx = 1d-3


      if (ldebug) print*,xxm2, xxmwidth, rm2minx, rm2maxx, rd(1), W
      if (ldebug) print*,"do Resonance fullwidth"

      ! wide resonce of W* (going to WH)
      if (.not. Resonance(xxm2, xxmwidth, rm2minx, rm2maxx,
     1       rd(1), w, pX(4))) return
      !this is pretty equal to
!      if (.not. NonResonance(rm2minx, rm2maxx, rd(1), W, pX(4))) return

      if (NWH .eq. 0) then
         if (ldebug) print*,"do Higgs resonance"
         if (ldebug) print*,rm2(2), rmg(2), rm2min(2), rm2max(2), pX(4), q(4,1)
         ! Higgs resonance, offshell
         if (.not.Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2),
     1       rd(3), w, q(4,2))) RETURN 
      endif
      if (NW .eq. 0) then
         if (ldebug) print*,rm2(1), rmg(1), rm2min(1), px(4), rd(2), W
         if (ldebug) print*,"do Resonance"
c q^2 of V is generated via tan mapping out of rd(2)
      ! resonance of W
         rm2max1=min(rm2max(1),(sqrt(pX(4))-sqrt(q(4,2)))**2 )
         !rm2max1 = rm2max(1)
         if (.not.Resonance(rm2(1), rmg(1), rm2min(1), rm2max1, 
     1       rd(2), w, q(4,1))) RETURN 
      endif !NW
      
      rdoffset = 0
      rdoffset2 = 0
      if (numbos .eq. 4) then
         ! randomly select which of the two bosons gets generated first
         ! the second one has a limited mass range and is far off-shell
         if (RandomNumber() .lt. 0.5d0) then
           i=3
         else 
           i=4
         endif
        if (ldebug) print*,"do Boson resonances"
c q^2 of Vi is generated via tan mapping out of rd(4)
        if (.not.Resonance(rm2(i), rmg(i), rm2min(i), q(4,2), 
     1       rd(4), w, q(4,i))) RETURN 
c q^2 of V(3-i) is generated via tan mapping out of rd(5)
      
      ! this cuts away half of phasespace (forcing mv2 < mv1)
      ! compensate this by w=w*2
         rm2maxw2 = min ( (sqrt(q(4,2))-sqrt(q(4,i)))**2 , q(4,i))
        if (.not.Resonance(rm2(7-i), rmg(7-i), rm2min(7-i), rm2maxw2, 
     1       rd(5), w, q(4,7-i))) RETURN 
        W=W*2
        ! these offsets allow both H->2l and H->4l to use same rd
        rdoffset = 2
        rdoffset2 = 6
      endif !numbos 4
       if(ldebug) then
          print*, 'numbos ', numbos
          print*, 'rdoffsets ',rdoffset, rdoffset2
       endif
        
c-------------------------------------------------------------------

      if (ldebug) print*,"do twotoone or twotojetplusX"
      if (N.ge.1) then
c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd
         ! rdoffset2 is either 0 (2l case) or 6 (4l case)
         ! rf gets used later as phi-angle in V-> V H
         ! non-vegas phi-dependence in first jet here
         rf = rd(9+rdoffset2)
         call TwoToJetsPlusX(N, rd(10+rdoffset2), rn, ecm, pTjmin, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJet, W)
         !if (W.eq.0) return
         if (W.eq.0) return

      else if (N.eq.0) then
         ! rf gets used later as phi-angle in V-> V H
         ! gets driven by non-vegas random number
         rf = rn
         call TwoToOne(rd(9+rdoffset2),ecm,pX(4),k1(0),k2(0),x1,x2,pX(0),W)
         if (W.eq.0) return
         !if (W.eq.0) return

      endif


c generate phase space for decay of X -> V1 + V2 out of rd(4), rd(rn/9)
      !if (.not. TwoBodyDecay(rd(3), rf, pX(0), pX(4), q(4,1), q(4,2), 
      if (.not. TwoBodyDecay(rd(4+rdoffset), rf, pX(0), pX(4), q(4,1), q(4,2), 
     1     q(0,1), q(0,2), W)) return
      
      ! force boosted higgs, e.g. p_T > 200 GeV
!      if (q(1,2)**2 + q(2,2)**2 .lt. 200**2) then
!         W=0d0
!         return
!      endif

      if (ldebug) print*,"do V-decay"
c generate phase space for decay of V1 -> l1 + l2 out of rd(5-6)
      call TwoBodyDecay0(rd(5+rdoffset),rd(6+rdoffset),q(0,1),q(4,1),
     1     v(0,1,1),v(0,2,1),W)
      if (ldebug) print*,"V-decay done"

 
      if (numbos .eq. 2) then
c generate phase space for decay of H -> l3 + l4 out of rd(7-8)
         ! as there is no angular correlation those random numbers have a flat distribution
         if (ldebug) print*,"do H-decay"
         call TwoBodyDecay0(rd(7),rd(8),q(0,2),q(4,2),
     1     v(0,3,1),v(0,4,1),w) ! w*= 1/(8*pi)
         if (ldebug) print*,"H-decay done"
      else if (numbos .eq. 4) then
         ! rdoffset2 above gurantees that the random numbers are used only here
         !   offset = 2 hardcoded below
         ! Decay H->VV
         ! fails if m(VV)>m(H) (should be eliminated by rm2max)
         if (ldebug) print*, 'trying decay H->VV'
         if (ldebug) print*, rd(9),rd(10),q(0,2),q(4,2),q(4,3),q(4,4)
      if (.not. TwoBodyDecay(rd(9),rd(10),q(0,2),q(4,2),q(4,3),q(4,4),
     1     q(0,3),q(0,4),w)) return
         if(ldebug) print*, 'decay H->VV passed'
         ! Decay V1 -> l l
      call TwoBodyDecay0(rd(11),rd(12),q(0,3),q(4,3),
     1     v(0,3,1),v(0,4,1),w) 
         ! Decay V2 -> l l
      call TwoBodyDecay0(rd(13),rd(14),q(0,4),q(4,4),
     1     v(0,5,1),v(0,6,1),w) 
      else
         print*, 'numbos in ps_qqVHj.F broken'
         stop
      end if


c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         do i=1,PS_dimension+6,3
            print*, rd(i), rd(i+1), rd(i+2)
         enddo
c         print*, "y     = ", y
c         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
c         print*, "y2    = ", y2

      endif
      if(lprintmomenta) then
         print*, "wgt   = ", w
         print*, "------ Momenta: -------"
         print 100, "ve", v(0,1,1), v(1,1,1), v(2,1,1), v(3,1,1)
         print 100, "e", v(0,2,1), v(1,2,1), v(2,2,1), v(3,2,1)
         print 100, "Hdec1", v(0,3,1), v(1,3,1), v(2,3,1), v(3,3,1)
         print 100, "Hdec2", v(0,4,1), v(1,4,1), v(2,4,1), v(3,4,1)
         if (numbos.eq.4) then
            print 100, "Hdec3", v(0,5,1), v(1,5,1), v(2,5,1), v(3,5,1)
            print 100, "Hdec4", v(0,6,1), v(1,6,1), v(2,6,1), v(3,6,1)
         endif

         print 100, "H", v(0:3,3,1)+v(0:3,4,1)+v(0:3,5,1)+v(0:3,6,1)
         print*, 'numjets', N
         do i=1,N
            print 110, "j",i, pjet(0,i), pjet(1,i), pjet(2,i), pjet(3,i)
         enddo

         print*, 'intermediate particles'
         do i=1,numbos
            print 110, "q",i, q(0,i), q(1,i), q(2,i), q(3,i)
110      FORMAT (2x, A1, I1, 2x, 4(F10.3, 1x))
100      FORMAT (A5, 1x, 4(F10.3, 1x))
         enddo
      endif
      if(ldebug) then
         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-q(mu,1)-q(mu,2)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*,'overall cons with bos: ', test
            test = q(mu,1)-v(mu,1,1)-v(mu,2,1)
            if (test.gt.1E-9) print*,'W dec: ', test
            if (numbos.eq.4) then
               test = q(mu,3)-v(mu,3,1)-v(mu,4,1)
               if (test.gt.1E-9) print*,'V1 dec: ', test
               test = q(mu,4)-v(mu,5,1)-v(mu,6,1)
               if (test.gt.1E-9) print*,'V2 dec: ', test
               test = q(mu,2)-v(mu,3,1)-v(mu,4,1)-v(mu,5,1)-v(mu,6,1)
               if (test.gt.1E-9) print*,'H dec: ', test
            else
               test = q(mu,2)-v(mu,3,1)-v(mu,4,1)
               if (test.gt.1E-9) print*,'H dec: ', test
            endif
            test = k1(mu)+k2(mu)-v(mu,1,1)-v(mu,2,1)-v(mu,3,1) -v(mu,4,1)-v(mu,5,1)-v(mu,6,1)
            do i=1,N
              test = test - pjet(mu,i)
            enddo
            if (test.gt.1E-9) print*,'overall cons with lep: ', test
         enddo
!         print*,"---tau*s = (q1 + q2)^2 ?---"
!         test = tau*s-dotrr(q(0,1),q(0,1))-dotrr(q(0,2),q(0,2))
!     &                -2*dotrr(q(0,1),q(0,2))
!         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         do i=1,N
           print*,"pjet(",i,")^2 = ",dotrr(pjet(0,i),pjet(0,i))
         enddo
         print 200,"sqrt(q1^2) = ",sqrt(dotrr(q(0,1),q(0,1)))
         print 200,"sqrt(q2^2) = ",sqrt(dotrr(q(0,2),q(0,2)))
         print 200,"v1^2 = m_e^2 = ",dotrr(v(0,1,1),v(0,1,1))
         print 200,"v2^2 = m_v^2 = ",dotrr(v(0,2,1),v(0,2,1))
         print 200,"v3^2 = m_v^2 = ",dotrr(v(0,3,1),v(0,3,1))
         print 200,"v4^2 = m_v^2 = ",dotrr(v(0,4,1),v(0,4,1))
         print 200,"v5^2 = m_v^2 = ",dotrr(v(0,5,1),v(0,5,1))
         print 200,"v6^2 = m_v^2 = ",dotrr(v(0,6,1),v(0,6,1))
         print 200,"sqrt((v3+v4)^2) = ",sqrt(dotrr(v(0:3,3,1)+v(0:3,4,1),v(0:3,3,1)+v(0:3,4,1)))
         print 200,"sqrt((v5+v6)^2) = ",sqrt(dotrr(v(0:3,5,1)+v(0:3,6,1),v(0:3,5,1)+v(0:3,6,1)))
         print 201,"mass V1, V2 = ",sqrt(q(4,3)), sqrt(q(4,4))
         print *, "-----------------------"
200      format (A,1x,f10.3)       
201      format (A,1x,2(f10.3))

      endif

      return
      end subroutine

