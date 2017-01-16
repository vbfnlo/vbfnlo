c*****************************************************************************
c
c    begin subroutine DEF_PS_QQAAAj
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Oliver Schlimpert: July 2012
c	Based on DEF_PS_QQVVAj by Giuseppe Bozzi
c
c   initialization of phase space generators for 
c         pp --> A + A + A + 1 Jet
c
c****************************************************************************

      SUBROUTINE DEF_PS_QQAAAj(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(3), i, bos1, bos2, bos3

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      REAL*8 RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      s = ecm**2

      bos(1) = bos1        ! A = 1
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
        write(6,*) "---------------------------------------------------"
        write(6,*)
      endif

      end

c*****************************************************************************
c
c    end subroutine DEF_PS_QQAAAj
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSAAAJ
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
c           V(0:3,4)   4-vectors for the photons
c                      V(*,1) = photon
c                      V(*,2) = photon
c                      V(*,3) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSAAAJ(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=10, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W, K(0:4)
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

C  local variables

      INTEGER MU

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 y_border, tau, q(0:4,3)  

      REAL*8 pY(0:4)


      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 100d0, xxm = 200d0 )
      PARAMETER ( xywidth = 200d0, xym = 80d0 ) !broader resonance and new mass

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSAAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,1)=0d0
      q(4,2)=0d0
      q(4,3)=0d0

cc alternative phase space for testing
c      ptmin=20d0
c 
c      if (TwoToJetsPlusX(N+2,RD,RN,ECM,PTmin,q(4,3),K1,K2,X1,X2,
c     &      q(0,3),v(0,1),W)) continue
c 
c      if (N.eq.1) then
c        do mu=0,3
c         pjet(mu) = v(mu,3)
c        enddo
c      else
c        do mu=0,3
c         pjet(mu) = 0.0d0
c        enddo
c      endif
c      do mu=0,3
c        v(mu,3) = q(mu,3)
c      enddo
c      return

c generate q^2_Y from rd(1)
      if (.not. Resonance(xym**2, xym*xywidth, rm2miny, rm2maxy,
     1                 rd(1), W, pY(4))) return

      IF(N.EQ.1) THEN ! 2 -> JET+Y
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(8),rd(9),rd(10)
         call TwoToJetPlusX(rd(8), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

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

c assign jet momenta
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0

         do mu =0,4
            pY(mu) = K(mu)
         enddo

         W = W * 1.d0/2d0/K(4)        ! Flux
         W = W * TPI                  ! (2*pi)^(4-3) from delta- and ps-int
         W = W * 3.89379304d11        ! conversion into fb

      ENDIF

c generate phase space for decay of Y -> A1 + A2 + A3 out of rd(3-7)
      call ThreeBodyDecay0(rd(3),pY(0),pY(4),
     1     v(0,1),v(0,2),v(0,3),W)

c -----------------------------------------------------

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7)
         if (N.eq.1) print*, rd(8), rd(9), rd(10)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
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
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2)) - 2*dotrr(Pjet(0),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
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
c    end subroutine LIPSAAAJ
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine LIPSAAAJ2_anomal
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 3 bosons (+ 1 jet)
c at hadron colliders, q q -> H A3 -> A1 A2 A3
c
c  INPUT:   RD         random number array (between 0 and 1)
c           RN         additional random number which will not be driven by
c                      VEGAS
c
c  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
c           X1,X2      corresponding values of Feynman x
c           V(0:3,4)   4-vectors for the photons
c                      V(*,1) = photon
c                      V(*,2) = photon
c                      V(*,3) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSAAAJ2_anomal(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=10, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3)
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

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  

      REAL*8 pY(0:4)


      REAL*8 rm2minh,rm2maxh

      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 100d0, xxm = 200d0 )
      PARAMETER ( xywidth = 200d0, xym = 80d0 ) 

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSAAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,1)=0d0
      q(4,2)=0d0
      q(4,3)=0d0


c generate q^2_Y from rd(1)
      rm2minh = ((xm2(6)-15d0*xmg(6))/sqrt(xm2(6)) )**2 
      rm2maxh = ((xm2(6)+15d0*xmg(6))/sqrt(xm2(6)) )**2
      if (.not. Resonance(xm2(6), xmg(6), rm2minh, rm2maxh,
     1                 rd(1), W, pY(4))) return

      IF(N.GE.1) THEN ! 2 -> Y + A + JET

         if (.not.TwoToJetsPlusX(N+1, rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pjets(0,1), W)) return


         do mu=0,3
           do i=1,N
             pjet(mu) = pjets(mu,i)
           enddo
           q(mu,1) = pjets(mu,N+1)
         enddo

      ELSE ! 2 -> Y + A

c generate 2 -> 2 process for k1 k2 -> q_A q_Y from rd(4),rd(5),rd(6)
         call TwoToJetPlusX(rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), q(0,1), W)
c         print*, rd(4), pY(0), v(0,1), W
         if (W.eq.0) return

      ENDIF
      

c generate phase space for decay of Y -> A1 + A2 + A3 out of rd(3-7)
      call TwoBodyDecay0(rd(2),rd(3),pY(0),pY(4),
     1     q(0,2),q(0,3),W)

c set photons
      do mu=0,3
         v(mu,1)=q(mu,2)
         v(mu,2)=q(mu,3)
         v(mu,3)=q(mu,1)
         v(mu,5)=0d0
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
         print*, rd(7)
         if (N.eq.1) print*, rd(8), rd(9), rd(10)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
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
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2)) - 2*dotrr(Pjet(0),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
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
c    end subroutine LIPSAAAJ2_anomal
c
c*****************************************************************************


c*****************************************************************************
c
c    begin subroutine LIPSAAAJ3_anomal
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 3 bosons (+ 1 jet)
c at hadron colliders, q q -> H A2 -> A1 A3 A2
c
c  INPUT:   RD         random number array (between 0 and 1)
c           RN         additional random number which will not be driven by
c                      VEGAS
c
c  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
c           X1,X2      corresponding values of Feynman x
c           V(0:3,4)   4-vectors for the photons
c                      V(*,1) = photon
c                      V(*,2) = photon
c                      V(*,3) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSAAAJ3_anomal(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=10, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3)
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

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  

      REAL*8 pY(0:4)


      REAL*8 rm2minh,rm2maxh

      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 100d0, xxm = 200d0 )
      PARAMETER ( xywidth = 200d0, xym = 80d0 ) 

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSAAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,1)=0d0
      q(4,2)=0d0
      q(4,3)=0d0


c generate q^2_Y from rd(1)
      rm2minh = ((xm2(6)-15d0*xmg(6))/sqrt(xm2(6)) )**2 
      rm2maxh = ((xm2(6)+15d0*xmg(6))/sqrt(xm2(6)) )**2
      if (.not. Resonance(xm2(6), xmg(6), rm2minh, rm2maxh,
     1                 rd(1), W, pY(4))) return

      IF(N.GE.1) THEN ! 2 -> Y + A + JET

         if (.not.TwoToJetsPlusX(N+1, rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pjets(0,1), W)) return


         do mu=0,3
           do i=1,N
             pjet(mu) = pjets(mu,i)
           enddo
           q(mu,1) = pjets(mu,N+1)
         enddo


      ELSE ! 2 -> Y + A

c generate 2 -> 2 process for k1 k2 -> q_A q_Y from rd(4),rd(5),rd(6)
         call TwoToJetPlusX(rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), q(0,1), W)
c         print*, rd(4), pY(0), v(0,1), W
         if (W.eq.0) return

      ENDIF
      

c generate phase space for decay of Y -> A1 + A2 + A3 out of rd(3-7)
      call TwoBodyDecay0(rd(2),rd(3),pY(0),pY(4),
     1     q(0,2),q(0,3),W)

c set photons
      do mu=0,3
         v(mu,1)=q(mu,2)
         v(mu,2)=q(mu,1)
         v(mu,3)=q(mu,3)
         v(mu,5)=0d0
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
         print*, rd(7)
         if (N.eq.1) print*, rd(8), rd(9), rd(10)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
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
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2)) - 2*dotrr(Pjet(0),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
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
c    end subroutine LIPSAAAJ3_anomal
c
c*****************************************************************************


c*****************************************************************************
c
c    begin subroutine LIPSAAAJ4_anomal
c
c*****************************************************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version: November 2009
c       Last modified by Michael Rauch: November 2009
c       Based on LIPSVVAJ by Giuseppe Bozzi
c
c phase space generator for the production of 3 bosons (+ 1 jet)
c at hadron colliders, q q -> H A1 -> A2 A3 A1
c
c  INPUT:   RD         random number array (between 0 and 1)
c           RN         additional random number which will not be driven by
c                      VEGAS
c
c  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
c           X1,X2      corresponding values of Feynman x
c           V(0:3,4)   4-vectors for the photons
c                      V(*,1) = photon
c                      V(*,2) = photon
c                      V(*,3) = photon
c           Pjet(0:3)  4-vector for the jet, assumed massless
c           W          the weight of the phase space point
c
c****************************************************************************

      SUBROUTINE LIPSAAAJ4_anomal(N,RD,RN,
     1                    K1,K2,V,Pjet,
     2                    X1,X2,W)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
c
      INTEGER N,NRANDL
      REAL*8 TPI
      PARAMETER (nrandl=10, TPI=2D0*PI )
      REAL*8 RD(nrandl),RN
      REAL*8 K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
      REAL*8 Pjet(0:3)
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

      REAL*8 Y, Y2


      REAL*8 test
      REAL*8 tau, q(0:4,3)  

      REAL*8 pY(0:4)


      REAL*8 rm2minh,rm2maxh

      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 100d0, xxm = 200d0 )
      PARAMETER ( xywidth = 200d0, xym = 80d0 ) 

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX
      EXTERNAL NONRESONANCE, RESONANCE, TWOBODYDECAY, TwoToJetsplusX

c--------------------------------------------------------------------
c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSAAA "

      W = 1d0             !init of weightfactor

c photons are massless
      q(4,1)=0d0
      q(4,2)=0d0
      q(4,3)=0d0


c generate q^2_Y from rd(1)
      rm2minh = ((xm2(6)-15d0*xmg(6))/sqrt(xm2(6)) )**2 
      rm2maxh = ((xm2(6)+15d0*xmg(6))/sqrt(xm2(6)) )**2
      if (.not. Resonance(xm2(6), xmg(6), rm2minh, rm2maxh,
     1                 rd(1), W, pY(4))) return

      IF(N.GE.1) THEN ! 2 -> Y + A + JET

         if (.not.TwoToJetsPlusX(N+1, rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pjets(0,1), W)) return


         do mu=0,3
           do i=1,N
             pjet(mu) = pjets(mu,i)
           enddo
           q(mu,1) = pjets(mu,N+1)
         enddo


      ELSE ! 2 -> Y + A

c generate 2 -> 2 process for k1 k2 -> q_A q_Y from rd(4),rd(5),rd(6)
         call TwoToJetPlusX(rd(4), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), q(0,1), W)
c         print*, rd(4), pY(0), v(0,1), W
         if (W.eq.0) return

      ENDIF
      

c generate phase space for decay of Y -> A1 + A2 + A3 out of rd(3-7)
      call TwoBodyDecay0(rd(2),rd(3),pY(0),pY(4),
     1     q(0,2),q(0,3),W)

c set photons
      do mu=0,3
         v(mu,1)=q(mu,1)
         v(mu,2)=q(mu,2)
         v(mu,3)=q(mu,3)
         v(mu,5)=0d0
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
         print*, rd(7)
         if (N.eq.1) print*, rd(8), rd(9), rd(10)
         print*, "y     = ", y
         print*, "tau   = ", tau
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "y2    = ", y2
         print*, "wgt   = ", w
         print*, "---------------------------"

         print*, "--- momenta conservation ? ---"
         do mu=0,3
            test = k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3)
            if (test.gt.1E-9) print*, test
            test = q(mu,1)-v(mu,1)-v(mu,2) 
            if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
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
     &     - dotrr(Pjet(0),Pjet(0)) - 2*dotrr(Pjet(0),q(0,1))
     &     - 2*dotrr(Pjet(0),q(0,2)) - 2*dotrr(Pjet(0),q(0,3))
         print*, test
         print*, "--- mass equal to zero ? ---"
         print*,"p1^2 = ",dotrr(k1(0),k1(0))
         print*,"p2^2 = ",dotrr(k2(0),k2(0))
         print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
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
c    end subroutine LIPSAAAJ4_anomal
c
c*****************************************************************************

