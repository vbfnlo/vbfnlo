c*****************************************************************************
c
c    begin subroutine DEF_PS_QQVVAj
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified by Michael Rauch: April 2011
c
c   initialization of phase space generators for 
c         pp --> V1 + V2 + A + 1 Jet
c                V1 --> l1 l2, V2 --> l3 l4
c
c****************************************************************************

      SUBROUTINE DEF_PS_QQVVAj(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(3), i, bos1, bos2, bos3
      common /psvvabos/bos

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
      bos(2) = bos2        ! Z = 2, W+ = 3, W- = 4
      bos(3) = bos3        ! A = 1

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.1) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3         ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0        ! for tan mapping
         endif
      enddo

      rm2miny = 1d-3
      rm2maxy = s/2d0

      ! for intermediate particle X -> V1 V2
      rm2minx = 1d-3
      rm2maxx = s/2d0
      xxm = 180d0
      xxwidth = 350d0

      ! select relevant resonances for VV invariant mass
      N_res = 0
      ! Higgs resonance in case of WWA/ZZA and anomalous couplings
      if (with_anom .and. xm2(6).lt.600d0**2 .and. (bos(1)+bos(2).eq.7 .or. bos(1)+bos(2).eq.4)) then
        N_res = N_res + 1
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        F_i(N_res)  = 0.5
      endif
      call MultiResonancesInit(1, xxm**2, xxm*xxwidth, 1d0, M2_i, MG_i, F_i, N_res)


      ! select relevant resonances for ZA invariant mass in case of anomalous couplings
      if (with_anom .and. (bos(1).eq.2 .or. bos(2).eq.2)) then
        xxm = sqrt(xm2(2))
        xxwidth = 500d0
        N_res = 1
        M2_i(N_res) = xm2(2)
        MG_i(N_res) = xmg(2)
        F_i(N_res)  = 5d0
        N_res = 2
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        F_i(N_res)  = 3d0
        call MultiResonancesInit(2, xxm**2, xxm*xxwidth, 2d0, M2_i, MG_i, F_i, N_res)
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
           if (i.eq.3.and.bos(3).eq.1) then
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
c    end subroutine DEF_PS_QQVVAj
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVVAJ
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified by Michael Rauch: April 2011
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
c           V(0:3,5)   4-vectors for the 4 decay leptons, assumed massless
c                      Q(*,1) = V(*,1) + V(*,2)
c                      Q(*,2) = V(*,3) + V(*,4)
c                      V(*,5) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSVVAJ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=15, TPI=2D0*PI )
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
      common /psvvabos/bos

C  local variables

      INTEGER I, MU

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  
      REAL*8 pX(0:4)


      REAL*8 RM2MAXV1, RM2MAXV2

!      DOUBLE PRECISION xxwidth, xxm
      DOUBLE PRECISION dotrr

!      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSVVA "

      W = 1d0             !init of weightfactor

      q(4,3)=0d0

c pX^2 of intermediate particle X -> V1 + V2 out of rd(3)
      if (.not. MultiResonance(1, rd(3), rm2minx, rm2maxx, W, pX(4))) return


c q^2 of V1 is generated via tan mapping out of rd(1)
      rm2maxv1 = pX(4)
      if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2maxv1, 16d0, bos(1), rd(1), W, Q(4,1))) RETURN


c q^2 of V2 is generated via tan mapping out of rd(2)
      rm2maxv2 = pX(4)+q(4,1)-2*sqrt(pX(4)*q(4,1))
      if (.not.ResonanceVbos(RM2(2), RMG(2), RM2MIN(2), rm2maxv2, 16d0, bos(2), rd(2), W, Q(4,2))) RETURN

      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y out of rd(10-12-15)
      if (.not.TwoToJetsPlusX(N+1, rd(10), rn, ecm, pTc, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJets(0,1), W)) return

      do mu=0,3
        do i=1,N
          pjet(mu,i) = pjets(mu,i)
        enddo
        q(mu,3) = pjets(mu,N+1)
      enddo

c generate phase space for decay of X -> V1 + V2 out of rd(4-5)
      if (.not. TwoBodyDecay(rd(4), rd(5), pX(0), pX(4), q(4,1), q(4,2),
     1     q(0,1),q(0,2),W)) return

c generate phase space for decay of V1 -> l1 + l2 out of rd(6-7)
      call TwoBodyDecay0(rd(6),rd(7),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c generate phase space for decay of V2 -> l3 + l4 out of rd(8-9)
      call TwoBodyDecay0(rd(8),rd(9),q(0,2),q(4,2),
     1     v(0,3),v(0,4),W)

c set q_A=v(5)
      do mu=0,3
         v(mu,5)=q(mu,3)
         v(mu,6)=0d0
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
         print*, rd(10), rd(11), rd(12)
         if (N.eq.1) print*, rd(13), rd(14), rd(15)
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
            test = q(mu,2)-v(mu,3)-v(mu,4) 
             if (test.gt.1E-9) print*, test
            test = q(mu,3)-v(mu,5)-v(mu,6) 
             if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu,1)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6) 
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
         print*,"v5^2 = ",dotrr(v(0,5),v(0,5))
         print*,"v6^2 = ",dotrr(v(0,6),v(0,6))
         print*, "-----------------------"
         read(*,*)
      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVVAJ
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVVAJ1
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

      SUBROUTINE LIPSVVAJ1(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"

      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=15, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3)
      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)
      INTEGER MU
      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  
      REAL*8 pX(0:4)


      INTEGER bos(3)
      common /psvvabos/bos

      DOUBLE PRECISION xxwidth, xxm
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSWWA1 "

      W = 1d0             !init of weightfactor

      q(4,3)=0d0

c q^2 of V1 is generated via tan mapping out of rd(1), for anom Z case including new Higgs contributions
      if (bos(1).eq.2 .and. with_anom) then
            if (.not. MultiResonance(2, rd(1), rm2min(1), rm2max(1), W, Q(4,1))) return
      else
            if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, Q(4,1))) RETURN
      endif


c q^2 of V2 is generated via tan mapping out of rd(2)
      if (.not.ResonanceVbos(RM2(2), RMG(2), RM2MIN(2), rm2max(2), 16d0, bos(2), rd(2), W, Q(4,2))) RETURN

      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

      if (N.eq.1) then

c pX^2 of intermediate particle X -> V1 + V2 out of rd(13)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2minx, rm2maxx,
     1        rd(13), W, pX(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd(3-5)
         call TwoToJetPlusX(rd(3), rn, ecm, pTc, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJet(0), W)
         if (w.eq.0) return

c generate phase space for decay of X -> V1 + V2 out of rd(14), rd(15)
         if (.not. TwoBodyDecay(rd(14), rd(15), pX(0), pX(4), q(4,1), q(4,2), 
     1        q(0,1), q(0,2), W)) return

      else if (N.eq.0) then

c generate 2 -> 2 process for k1 k2 -> q_V1 q_V2 out of rd(3-5)
         call TwoToTwo(rd(3), rn, ecm, q(4,1), q(4,2),
     1        k1(0), k2(0), x1, x2, q(0,1), q(0,2), W)
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0
      endif

c generate phase space for decay of V1 -> l1 + l2 + A out of rd(6-10)
      call ThreeBodyDecay0(rd(6),q(0,1),q(4,1),
     1     v(0,1),v(0,2),q(0,3),W)
      
c generate phase space for decay of V2 -> l3 + l4 out of rd(11-12)
      call TwoBodyDecay0(rd(12),rd(11),q(0,2),q(4,2),
     1     v(0,3),v(0,4),W)

c set q_A=v(5)
      do mu=0,3
         v(mu,5)=q(mu,3)
         v(mu,6)=0d0
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
         print*, rd(10), rd(11), rd(12)
         if (N.eq.1) print*, rd(13), rd(14), rd(15)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)
            if (test.gt.1E-9) print*,test
            test = q(mu,1)-v(mu,1)-v(mu,2)-v(mu,5) 
            if (test.gt.1E-9) print*,test
            test = q(mu,2)-v(mu,3)-v(mu,4) 
            if (test.gt.1E-9) print*,test
            test = q(mu,3)-v(mu,5)-v(mu,6) 
            if (test.gt.1E-9) print*,test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
     &           -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6) 
            if (test.gt.1E-9) print*,test
         enddo
         tau = x1*x2
         s = ecm**2
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,2))
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*,"v5^2 = ",dotrr(v(0,5),v(0,5))
         print*,"v6^2 = ",dotrr(v(0,6),v(0,6))
         print*, "-----------------------"
         read(*,*)

      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVVAJ1
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSVVAJ2
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified by Vera Hankele: January 2009
c
c Alternative phase-space for pp --> V1 V2 (--> l1 l2 A)
c
c****************************************************************************

      SUBROUTINE LIPSVVAJ2(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=15, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3)
      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)
      INTEGER MU
      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  
      REAL*8 pX(0:4)


      INTEGER bos(3)
      common /psvvabos/bos

      DOUBLE PRECISION xxwidth, xxm
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, RESONANCEVBOS, multiresonance, TWOTOJETSPLUSX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSWWA2 "

      W = 1d0             !init of weightfactor

      q(4,3)=0d0

c q^2 of V1 is generated via tan mapping out of rd(1)
      if (.not.ResonanceVbos(RM2(1), RMG(1), RM2MIN(1), rm2max(1), 16d0, bos(1), rd(1), W, Q(4,1))) RETURN


c q^2 of V2 is generated via tan mapping out of rd(2), for anom Z case including new Higgs contributions
      if (bos(2).eq.2 .and. with_anom) then
            if (.not. MultiResonance(2, rd(2), rm2min(2), rm2max(2), W, Q(4,2))) return
      else
            if (.not.ResonanceVbos(RM2(2), RMG(2), RM2MIN(2), rm2max(2), 16d0, bos(2), rd(2), W, Q(4,2))) RETURN
      endif

      
c-------------------------------------------------------------------
c-------------------------------------------------------------------

      if (N.eq.1) then

c pX^2 of intermediate particle X -> V1 + V2 out of rd(13)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2minx, rm2maxx,
     1        rd(13), W, pX(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_X out of rd(3-5)
         call TwoToJetPlusX(rd(3), rn, ecm, pTc, pX(4),
     1        k1(0), k2(0), x1, x2, pX(0), pJet(0), W)
         if (w.eq.0) return

c generate phase space for decay of X -> V1 + V2 out of rd(14), rd(15)
         if (.not. TwoBodyDecay(rd(14), rd(15), pX(0), pX(4), q(4,1), q(4,2), 
     1        q(0,1), q(0,2), W)) return

      else if (N.eq.0) then

c generate 2 -> 2 process for k1 k2 -> q_V1 q_V2 out of rd(3-5)
         call TwoToTwo(rd(3), rn, ecm, q(4,1), q(4,2),
     1        k1(0), k2(0), x1, x2, q(0,1), q(0,2), W)
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0
      endif

c generate phase space for decay of V1 -> l1 + l2 out of rd(6-7)
      call TwoBodyDecay0(rd(6),rd(7),q(0,1),q(4,1),
     1     v(0,1),v(0,2),W)
      
c generate phase space for decay of V2 -> l3 + l4 + A out of rd(8-12)
      call ThreeBodyDecay0(rd(8),q(0,2),q(4,2),
     1     v(0,3),v(0,4),q(0,3),W)

c set q_A=v(5)
      do mu=0,3
         v(mu,5)=q(mu,3)
         v(mu,6)=0d0
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
         print*, rd(10), rd(11), rd(12)
         if (N.eq.1) print*, rd(13), rd(14), rd(15)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = q(mu,2)-v(mu,3)-v(mu,4)-v(mu,5) 
             if (test.gt.1E-9) print*, test
            test = q(mu,3)-v(mu,5)-v(mu,6) 
             if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6) 
             if (test.gt.1E-9) print*, test
         enddo
         print*,"tau*s = (k1+k2)^2"
         test = tau*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,"---tau*s = (q1 + q2 + pjet)^2 ?---"
         test = tau*s - dotrr(q(0,1),q(0,1)) - dotrr(q(0,2),q(0,2))
     &     - 2*dotrr(q(0,1),q(0,2))
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
         print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         print*,"v2^2 = ",dotrr(v(0,2),v(0,2))
         print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         print*,"v5^2 = ",dotrr(v(0,5),v(0,5))
         print*,"v6^2 = ",dotrr(v(0,6),v(0,6))
         print*, "-----------------------"
         read(*,*)

      endif

      return
      end

c*****************************************************************************
c
c    end subroutine LIPSVVAJ2
c
c*****************************************************************************
