!
! This file contains the general multi-channel phasespaces for VAA production,
! including the case of anomalous couplings with extra Higgs contributions.
!


c*****************************************************************************
c
c    begin subroutine DEF_PS_QQVAAj
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c	Based on DEF_PS_QQVVAj by Giuseppe Bozzi
c
c   initialization of phase space generators for 
c         pp --> V1 + A + A + 1 Jet
c                V1 --> l1 l2
c
c****************************************************************************

      SUBROUTINE DEF_PS_QQVAAj(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(3), i, bos1, bos2, bos3
      common /psvaabos/bos

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      double precision M2_i(max_res),MG_i(max_res),F_i(max_res)
      integer N_res
      DOUBLE PRECISION xxwidth, xxm

      s = ecm**2

      bos(1) = bos1        ! Z = 2, W+ = 3, W- = 4
      bos(2) = bos2        ! A = 1
      bos(3) = bos3        ! A = 1

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.1) then
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

      ! select relevant resonances for ZA invariant mass in case of anomalous couplings
      if (with_anom .and. (bos(1).eq.2 .or. bos(2).eq.2)) then
        xxm = xm2(2)
        xxwidth = 500d0
        N_res = 1
        M2_i(N_res) = xm2(2)
        MG_i(N_res) = xmg(2)
        F_i(N_res)  = 5d0
        N_res = 2
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        F_i(N_res)  = 3d0
        call MultiResonancesInit(1, xxm**2, xxm*xxwidth, 2d0, M2_i, MG_i, F_i, N_res)
      endif


c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10.0d0)

 60   format ( " machine energy = ", f10.1, "  GeV" )
 70   format ( " integration range for boson ",i3," = ",f10.3,
     1         "  to ",f12.3," GeV" )
 72   format ( " boson has mass ",f8.2," GeV and width ",
     1         f10.4," GeV" )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) "--------------------PHASE SPACE--------------------"
        write(6,*)
        do i = 1,3
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
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end

c*****************************************************************************
c
c    end subroutine DEF_PS_QQVAAj
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAAJ1
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 3 bosons (+ 1 jet)
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
c           Pjets(0:3,N)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation
c
c****************************************************************************

      SUBROUTINE LIPSVAAJ1(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3,N)
      REAL*8 Pjets(0:3,N+2)

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER bos(3)
      common /psvaabos/bos

C  local variables

      INTEGER I, MU

      REAL*8 Y


      REAL*8 test
      REAL*8 tau, q(0:4,3)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos


      double precision phirn, cosrn, sinrn, tmpx, tmpy

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,2)=0d0
      q(4,3)=0d0


c generate q^2_W from rd(1)
      if (NW .eq. 1) then
c need at least mass of the vector boson
        q(4,1) = rm2(1)
        W = W * rmg(1)/2.0d0
      else
        if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, Q(4,1))) RETURN
      endif

!c generate 2 -> 1+2+N process for k1 k2 -> Pjets q_W from rd(4-9/12)
       if (.not.TwoToJetsPlusX(N+2, rd(4), 0d0, ecm, pTc, q(4,1),
     1         k1(0), k2(0), x1, x2, q(0,1), pjets(0,1), W)) return

*       if (W.eq.0) return

       if (N.ge.1) then
         do i=1,N
           do mu=0,3
             pjet(mu,i) = pjets(mu,i)
           enddo
         enddo
       endif

       do mu=0,3
         q(mu,2) = pjets(mu,N+1)
         q(mu,3) = pjets(mu,N+2)
       enddo

c-------------------------------------------------------------------

c generate phase space for decay of V1 -> l1 + l2 out of rd(2-3)
      call TwoBodyDecay0(rd(2),rd(3),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c set photons
      do mu=0,3
         v(mu,3)=q(mu,2)
         v(mu,4)=q(mu,3)
      enddo

c rotate by rn
      phirn = TPI*rn
      cosrn = cos(phirn)
      sinrn = sin(phirn)
      do i=1,4
        tmpx =  cosrn*v(1,i)+sinrn*v(2,i)
        tmpy = -sinrn*v(1,i)+cosrn*v(2,i)
        v(1,i) = tmpx
        v(2,i) = tmpy
      enddo
      do i=1,N
        tmpx =  cosrn*pjet(1,i)+sinrn*pjet(2,i)
        tmpy = -sinrn*pjet(1,i)+cosrn*pjet(2,i)
        pjet(1,i) = tmpx
        pjet(2,i) = tmpy
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
            test = k1(mu)+k2(mu)-Pjet(mu,1)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - dotrr(q(0,3),q(0,3)) - 2*dotrr(q(0,1),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,3)) - 2*dotrr(q(0,2),q(0,3))
     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
     &     - 2*dotrr(Pjet(0,1),q(0,2)) - 2*dotrr(Pjet(0,1),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0,1),Pjet(0,1))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAAJ1
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAAJ2
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified by Vera Hankele: January 2009
c
c Alternative phase-space for pp --> V1 (--> l1 l2 A) V2
c
c****************************************************************************

      SUBROUTINE LIPSVAAJ2(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
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

      INTEGER bos(3)
      common /psvaabos/bos

C  local variables

      INTEGER I, MU

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos, MultiResonance
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos, MultiResonance

      double precision phirn, cosrn, sinrn, tmpx, tmpy

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVAA "

      W = 1d0             !init of weightfactor

c for narrow width just use PS1 here too
      if (NW .eq. 1) then
        call LIPSVAAJ1(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)
        return
      endif

c photons are massless
      q(4,2)=0d0
      q(4,3)=0d0

c p^2 of intermediate particle V1 -> l + l~ + A2 out of rd(1)
       if (bos(1).eq.2 .and. with_anom) then
             if (.not. MultiResonance(1, rd(1), rm2min(1), rm2max(1), W, Q(4,1))) return
       else
             if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, Q(4,1))) RETURN
       endif

      IF(N.GE.1) THEN ! 2 -> JET+Y
!c generate 2 -> 3 process for k1 k2 -> Pjet q_Y A from rd(7-12)
         if (.not.TwoToJetsPlusX(N+1, rd(7), 0d0, ecm, pTc, q(4,1),
     1           k1(0), k2(0), x1, x2, q(0,1), pjets(0,1), W)) return


         do mu=0,3
           do i=1,N
             pjet(mu,i) = pjets(mu,i)
           enddo
           q(mu,3) = pjets(mu,N+1)
         enddo


      ELSE ! 2 -> X + A3
c generate 2 -> 2 process for k1 k2 -> q_X q_A out of rd(7-9)
         call TwoToJetPlusX(rd(7), 0d0, ecm, pTc, q(4,1),
     1           k1(0), k2(0), x1, x2, q(0,1), q(0,3), W)

         if (W.eq.0) return

      ENDIF

c generate phase space for decay of V1 -> l1 + l2 + A2 out of rd(2-6)
      call ThreeBodyDecay0(rd(2),q(0,1),q(4,1),
     1     v(0,1),v(0,2),q(0,2),W)
      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

c set photons
      do mu=0,3
         v(mu,3)=q(mu,2)
         v(mu,4)=q(mu,3)
         v(mu,5)=0d0
         v(mu,6)=0d0
      enddo

c rotate by rn
      phirn = TPI*rn
      cosrn = cos(phirn)
      sinrn = sin(phirn)
      do i=1,4
        tmpx =  cosrn*v(1,i)+sinrn*v(2,i)
        tmpy = -sinrn*v(1,i)+cosrn*v(2,i)
        v(1,i) = tmpx
        v(2,i) = tmpy
      enddo
      do i=1,N
        tmpx =  cosrn*pjet(1,i)+sinrn*pjet(2,i)
        tmpy = -sinrn*pjet(1,i)+cosrn*pjet(2,i)
        pjet(1,i) = tmpx
        pjet(2,i) = tmpy
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
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu,1)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - dotrr(q(0,3),q(0,3)) - 2*dotrr(q(0,1),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,3)) - 2*dotrr(q(0,2),q(0,3))
     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
     &     - 2*dotrr(Pjet(0,1),q(0,2)) - 2*dotrr(Pjet(0,1),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0,1),Pjet(0,1))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAAJ2
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAAJ3
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified by Vera Hankele: January 2009
c
c Alternative phase-space for pp --> V1 (--> l1 l2 A) V2
c
c****************************************************************************

      SUBROUTINE LIPSVAAJ3(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
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

      INTEGER bos(3)
      common /psvaabos/bos

C  local variables

      INTEGER I, MU

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos, MultiResonance
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos, MultiResonance

      double precision phirn, cosrn, sinrn, tmpx, tmpy

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVAA "

      W = 1d0             !init of weightfactor

c for narrow width just use PS1 here too
      if (NW .eq. 1) then
        call LIPSVAAJ1(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)
        return
      endif

c photons are massless
      q(4,2)=0d0
      q(4,3)=0d0

c p^2 of intermediate particle V1 -> l + l~ + A2 out of rd(1)
       if (bos(1).eq.2 .and. with_anom) then
             if (.not. MultiResonance(1, rd(1), rm2min(1), rm2max(1), W, Q(4,1))) return
       else
             if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, Q(4,1))) RETURN
       endif

      IF(N.GE.1) THEN ! 2 -> JET+Y
c generate 2 -> 3 process for k1 k2 -> Pjet q_Y A from rd(7-12)
         if(.not.TwoToJetsPlusX(N+1, rd(7), 0d0, ecm, pTc, q(4,1),
     1           k1(0), k2(0), x1, x2, q(0,1), pjets(0,1), W)) return

         do mu=0,3
           do i=1,N
             pjet(mu,i) = pjets(mu,i)
           enddo
           q(mu,3) = pjets(mu,N+1)
         enddo


      ELSE ! 2 -> X + A3
c generate 2 -> 2 process for k1 k2 -> q_X q_A out of rd(7-9)
         call TwoToJetPlusX(rd(7), 0d0, ecm, pTc, q(4,1),
     1           k1(0), k2(0), x1, x2, q(0,1), q(0,3), W)

         if (W.eq.0) return

      ENDIF

c generate phase space for decay of V1 -> l1 + l2 + A2 out of rd(2-6)
      call ThreeBodyDecay0(rd(2),q(0,1),q(4,1),
     1     v(0,1),v(0,2),q(0,2),W)
      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

c set photons
      do mu=0,3
         v(mu,3)=q(mu,3)
         v(mu,4)=q(mu,2)
         v(mu,5)=0d0
         v(mu,6)=0d0
      enddo

c rotate by rn
      phirn = TPI*rn
      cosrn = cos(phirn)
      sinrn = sin(phirn)
      do i=1,4
        tmpx =  cosrn*v(1,i)+sinrn*v(2,i)
        tmpy = -sinrn*v(1,i)+cosrn*v(2,i)
        v(1,i) = tmpx
        v(2,i) = tmpy
      enddo
      do i=1,N
        tmpx =  cosrn*pjet(1,i)+sinrn*pjet(2,i)
        tmpy = -sinrn*pjet(1,i)+cosrn*pjet(2,i)
        pjet(1,i) = tmpx
        pjet(2,i) = tmpy
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
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu,1)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - dotrr(q(0,3),q(0,3)) - 2*dotrr(q(0,1),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,3)) - 2*dotrr(q(0,2),q(0,3))
     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
     &     - 2*dotrr(Pjet(0,1),q(0,2)) - 2*dotrr(Pjet(0,1),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0,1),Pjet(0,1))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAAJ3
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVAAJ4
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 3 bosons (+ 1 jet)
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

      SUBROUTINE LIPSVAAJ4(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=12, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W, K(0:4)
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

      INTEGER bos(3)
      common /psvaabos/bos

C  local variables

      INTEGER MU

      REAL*8 Y


      REAL*8 test
      REAL*8 y_border, tau, q(0:4,3)  
      REAL*8 pX(0:4)
      REAL*8 pY(0:4)


      DOUBLE PRECISION dotrr
      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX, ResonanceVbos

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,2)=0d0
      q(4,3)=0d0


c generate q^2_Y from rd(1)
      if (NW .eq. 1) then
c need at least mass of the vector boson
        if (.not. NonResonance(rm2(1), rm2max(1),
     1                  rd(1), W, pY(4))) return
      else
c model flat part plus vector boson peak (inside ResonanceVbos)
        if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, pY(4))) RETURN
      endif

      IF(N.GE.2) THEN ! 2 -> JETS+Y
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(10),rd(11),rd(12)
         if (.not.TwoToJetsPlusX(N, rd(10), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pjet(0,1), W)) return

      ELSE IF(N.EQ.1) THEN ! 2 -> JET+Y
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(10),rd(11),rd(12)
         call TwoToJetPlusX(rd(10), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0,1), W)

         if (W.eq.0) return

      ELSE ! 2 -> Y
c Generate y, tau -> x1, x2, shat from RD(1) and RD(2) first **************
         tau = pY(4)/s  !rd(1)
         W = W / s 
         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border !(prop density)^-1 of y

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)

c generate incoming momenta in lab frame out of x1 x2
         k1(0) = 0.5d0 * x1 * ecm
         k1(1) = 0d0
         k1(2) = 0d0
         k1(3) = 0.5d0 * x1 * ecm

         k2(0) = 0.5d0 * x2 * ecm
         k2(1) = 0d0
         k2(2) = 0d0
         k2(3) = - 0.5d0 * x2 * ecm

         do mu = 0, 3
            K(mu) = k1(mu)+k2(mu)
         enddo

         K(4) = K(0)**2 - K(1)**2 - K(2)**2 - K(3)**2

         do mu =0,4
            pY(mu) = K(mu)
         enddo

         W = W * 1.d0/2d0/K(4)        ! Flux
         W = W * TPI                  ! (2*pi)^(4-3) from delta- and ps-int
         W = W * 3.89379304d11        ! conversion into fb

      ENDIF

c pX^2 of intermediate particle X -> V1 + A3 out of rd(3)
      if (NW .eq. 1) then
        if (.not. NonResonance(rm2(1), pY(4),
     1       rd(3), W, pX(4))) return
      else
        if (.not. NonResonance(rm2minx, pY(4),
     1       rd(3), W, pX(4))) return
      endif

c q^2 of V1 is generated via tan mapping out of rd(4)
      if (NW .eq. 1) then
c q^2_1 in narrow width approximation
        q(4,1) = rm2(1)
        W = W * rmg(1)/2.0d0
      else
        if (.not.NonResonance(rm2min(1), pX(4), rd(4), W, q(4,1))) RETURN 
      endif

c-------------------------------------------------------------------
c-------------------------------------------------------------------

c generate phase space for decay of Y -> X + A2 
      IF(N.GE.1) THEN ! 2 -> JET+Y
         if (.not. TwoBodyDecay(rd(5), rd(2), pY(0), pY(4), pX(4), q(4,2),
     1        pX(0),q(0,2),W)) return
      else
         if (.not. TwoBodyDecay(rd(5), 0d0, pY(0), pY(4), pX(4), q(4,2),
     1        pX(0),q(0,2),W)) return
      endif

c generate phase space for decay of X -> V1 + A3 out of rd(6-7)
      if (.not. TwoBodyDecay(rd(6), rd(7), pX(0), pX(4), q(4,1), q(4,3),
     1     q(0,1),q(0,3),W)) return

c generate phase space for decay of V1 -> l1 + l2 out of rd(8-9)
      call TwoBodyDecay0(rd(8),rd(9),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c set photons
      do mu=0,3
         v(mu,3)=q(mu,2)
         v(mu,4)=q(mu,3)
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
            test = k1(mu)+k2(mu)-Pjet(mu,1)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + q3 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - dotrr(q(0,3),q(0,3)) - 2*dotrr(q(0,1),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,3)) - 2*dotrr(q(0,2),q(0,3))
     &     - dotrr(Pjet(0,1),Pjet(0,1)) - 2*dotrr(Pjet(0,1),q(0,1))
     &     - 2*dotrr(Pjet(0,1),q(0,2)) - 2*dotrr(Pjet(0,1),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0,1),Pjet(0,1))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAAJ4
c
c*****************************************************************************


c*****************************************************************************
c
c    begin subroutine LIPSVAAJ5_anomal
c
c*****************************************************************************
c
c       Oliver Schlimpert, <olivers@particle.uni-karlsruhe.de>
c       Initial version: June 2012
c
c phase space generator for the production of 3 bosons (+ 1 jet)
c at hadron colliders
c
c q q -> Z H (j) ; H -> gamma gamma
c
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
c           Pjets(0:3,N)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c           NW         NW = 1 corresponds to Narrow Width Approximation,
c                      NW = 0 means Full calculation
c
c****************************************************************************

      SUBROUTINE LIPSVAAJ5_anomal(N,RD,RN,
     1                     K1,K2,V,Pjet,
     2                     X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL,NW
      REAL*8 TPI
      PARAMETER (nrandl=15, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
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

C  local variables

      INTEGER MU


      REAL*8 q(0:4,3)  
      REAL*8 pX(0:4)


      REAL*8 rm2minh, rm2maxh

      DOUBLE PRECISION V1mg, xxwidth, xxm
      DOUBLE PRECISION dotrr

      PARAMETER ( V1mg = 60d1) !broader V1 resonance
      PARAMETER ( xxwidth = 36d1, xxm = 17d1 )

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TWOTOJETSPLUSX


c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,2)=0d0
      q(4,3)=0d0


c q^2 of V is generated via tan mapping out of rd(1)
      if (.not.Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1), 
     1     rd(1), W, q(4,1))) RETURN 

c q^2 of higgs is generated via tan mapping out of rd(2)
      rm2minh = ((xm2(6)-15d0*xmg(6))/sqrt(xm2(6)) )**2 
      rm2maxh = ((xm2(6)+15d0*xmg(6))/sqrt(xm2(6)) )**2
      if (.not.Resonance(xm2(6), xmg(6), rm2minh, rm2maxh, 
     1     rd(2), W, q(4,2))) RETURN 
      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

      if (N.ge.2) then

c pX^2 of intermediate particle X -> V1 + V2 out of rd(13)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2minx, rm2maxx,
     1        rd(13), W, pX(4))) return

c generate 2 -> 3 process for k1 k2 -> Pjet1 Pjet2 q_X from rd(3-5)
         if (.not.TwoToJetsPlusX(N, rd(3), 0d0, ecm, pTc, pX(4),
     1           k1(0), k2(0), x1, x2, pX(0), pjet(0,1), W)) return

c generate phase space for decay of X -> V1 + V2 out of rd(14), rd(15)
         if (.not. TwoBodyDecay(rd(14), rd(15), pX(0), pX(4), q(4,1), q(4,2), 
     1        q(0,1), q(0,2), W)) return


      elseif (N.eq.1) then

c pX^2 of intermediate particle X -> V1 + V2 out of rd(13)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2minx, rm2maxx,
     1        rd(13), W, pX(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd(3-5)
         call TwoToJetPlusX(rd(3), rn, ecm, pTc, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJet(0,1), W)
         if (w.eq.0) return

c generate phase space for decay of X -> V1 + V2 out of rd(14), rd(15)
         if (.not. TwoBodyDecay(rd(14), rd(15), pX(0), pX(4), q(4,1), q(4,2), 
     1        q(0,1), q(0,2), W)) return

      else if (N.eq.0) then

c generate 2 -> 2 process for k1 k2 -> q_V1 q_V2 out of rd(3-5)
         call TwoToTwo(rd(3), rn, ecm, q(4,1), q(4,2),
     1        k1(0), k2(0), x1, x2, q(0,1), q(0,2), W)
         pjet(0,1) = 0.0d0
         pjet(1,1) = 0.0d0
         pjet(2,1) = 0.0d0
         pjet(3,1) = 0.0d0
      endif

c generate phase space for decay of V1 -> l1 + l2 out of rd(6-7)
      call TwoBodyDecay0(rd(6),rd(7),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c generate phase space for decay of H -> A + A out of rd(8-12)
      call TwoBodyDecay0(rd(8),rd(9),q(0,2),q(4,2),
     1     v(0,3),v(0,4),W)

c set q_A=v(5)
      do mu=0,3
         v(mu,5)=0d0
         v(mu,6)=0d0
      enddo


      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVAAJ5_anomal
c
c*****************************************************************************



