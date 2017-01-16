c
c phase space generator for pp --> V1 + V2 + V3 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4, V3 --> l5 l6
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: January 2009
c       Last modified: June 2007

      SUBROUTINE DEF_PS_ZZZj(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"
#include "VBFNLO/utilities/cuts.inc"

      INTEGER bos(3), i, bos1, bos2, bos3

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      DOUBLE PRECISION RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      s = ecm**2

      bos(1) = bos1        ! Z = 2
      bos(2) = bos2        ! Z = 2
      bos(3) = bos3        ! Z = 2

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.2) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2D0     
         else
            rm2min(i) = 1d-3         ! minimal/maximal value of mass squared
            rm2max(i) = s/2D0        ! for tan mapping
         endif
      enddo

      rm2miny = 1d-3
      rm2maxy = s/2D0

      rm2minx = 1d-3
      rm2maxx = s/2D0   ! for intermediate particle X -> V1 V2

c      rm2max1 = (sqrt(xm2(6)-30d0*xmg(6)))**2 !Dieters proposal
c      rm2max2 = (sqrt(xm2(6)+30d0*xmg(6)))**2 !Dieters proposal

      rm2max1 = ((xm2(6)-15d0*xmg(6))/sqrt(xm2(6)) )**2
      rm2max2 = ((xm2(6)+15d0*xmg(6))/sqrt(xm2(6)) )**2

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

C**********************  LIPSZZZJ   ****************************************
C
      SUBROUTINE LIPSZZZJ(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W, NW, PS_number)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"

      integer ps_number
c
      INTEGER N,NRANDL, NW
      DOUBLE PRECISION TPI
      PARAMETER (nrandl=19, TPI=2D0*PI )
      DOUBLE PRECISION RD(nrandl),RN
      DOUBLE PRECISION K1(0:3), K2(0:3), V(0:3,6), X1, X2, W, K(0:4)
      DOUBLE PRECISION Pjet(0:3)

C Phase space generator for the production of 3 bosons (+ 1 jet)
C at hadron colliders
C
C  INPUT:   RD         random number array (between 0 and 1)
C           RN         additional random number which will not be driven by
C                      VEGAS
C           NW         NW = 1 corresponds to Narrow Width Approximation,
C                      NW = 0 means Full calculation
C
C  OUTPUT:  K1,K2      4-vectors for 2 incoming partons
C           X1,X2      corresponding values of Feynman x
C           V(0:3,6)   4-vectors for the 4 decay leptons, assumed massless
C                      Q(*,1) = V(*,1) + V(*,2)
C                      Q(*,2) = V(*,3) + V(*,4)
C                      Q(*,3) = V(*,5) + V(*,6)
C           Pjet(0:3)  4-vector for the jet, assumed massless
C           W          the weight of the phase space point
C


      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
C
C   RM2            resonance mass squared
C   RMG            resonance total width * mass
C   RM2_MIN;_MAX   min & max values of the inv. mass**2 of resonance
C   PTJMIN         minimal PT value of the jet

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

C  local variables

      INTEGER MU
      DOUBLE PRECISION m1sq, m2sq, m3sq

      DOUBLE PRECISION test
      DOUBLE PRECISION q(0:4,3)  
      DOUBLE PRECISION pY(0:4)
      DOUBLE PRECISION tauc
      DOUBLE PRECISION q1(0:3), q2(0:3), q3(0:3)

      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION dotrr

      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )
      PARAMETER ( xywidth = 25d1, xym = 45d1 )

      EXTERNAL dotrr
      DOUBLE PRECISION A1, A2, tau, y, y_border

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL RESONANCE, TWOBODYDECAY, ThreeBodyDecay, ThreeBodyDecayPlain
      LOGICAL ThreeBodyDecay3Higgs
      EXTERNAL RESONANCE, TWOBODYDECAY, ThreeBodyDecay, ThreeBodyDecayPlain
      EXTERNAL ThreeBodyDecay3Higgs

c data to be used in histogram routine
      double precision feynmanx1, feynmanx2
      Common /histdata2/ tau, y, feynmanx1, feynmanx2
c--------------------------------------------------------------------

c---------------------------------------------------------------------

      if (ldebug) print*," enter LIPSWWZ "

      W = 1d0             !init of weightfactor

      if (sqrt(xm2(6)) .le. 2.0d0*sqrt(xm2(2))) then ! m_Higgs < 2 * m_Z
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     m_Higgs < 2 * m_Z -> Narrow Higgs resonances
C     Phase space 1: non-Higgs contributions
C     Phase space 2: Higgs contributions with help of msq_min
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c*********************************************************************
c
c Phase space 1: all non-Higgs-contributions
c                and p p -> Z2 H -> Z2 Z1 Z3
c
c*********************************************************************

      if (ps_number.eq.1) then

      IF(N.EQ.1) THEN ! 2 -> JET+Y

c generate q^2_Y which then decays into three massive particles from rd(1)
         if (.not. Resonance(xym**2, xym*xywidth, 1.0d-6, s,
     1                               rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(17),rd(18),rd(19)
         call TwoToJetPlusX(rd(17), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

         if (w.eq.0) return

      ELSE ! 2 -> Y
c Generate y, tau -> x1, x2, shat from RD(1) and RD(2) first **************
         tauc = 0.0001d0
         A1 = 0.2d0/tauc
         A2 = -0.8d0/log(tauc)
         if ( RD(1) .le. 0.2d0 ) then
            tau = RD(1)/A1
            W = W/A1
         else
            tau = (1.0d0/tauc)**((RD(1)-1)/0.8d0)
            W = W/A2*tau
         endif

         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border !(prop density)^-1 of y

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)

         feynmanx1 = x1
         feynmanx2 = x2

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

      ENDIF

c generate phase space for decay of P -> V + V + V out of rd(3):rd(10)
c also in NW
         if (.not. ThreeBodyDecayPlain(rd(3), pY(0), q1(0), q2(0),
     1                     q3(0), m1sq, m2sq, m3sq, W, NW) ) return

         if (W.eq.0) return


c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(11),rd(12),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(13),rd(14),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(15),rd(16),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************
      if(N.EQ.1) then
         W=W
      else
         W = W * 1.d0/2d0/K(4)       ! Flux
         W = W * 3.89379304d11        !conversion into fb
      endif

c*********************************************************************
c
c Phase space 2: p p -> Z1 H -> Z1 Z2 Z3 and
c                p p -> Z3 H -> Z3 Z1 Z2
c
c*********************************************************************

      elseif (ps_number.eq.2) then

         if(N.eq.1) then ! 2 -> 3 + jet
c generate q^2_Y which then decays into three massive particles from rd(1)
            if (.not. Resonance(xym**2, xym*xywidth, 1.0d-6, s,
     1                               rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(17), rd(18),rd(19)
            call TwoToJetPlusX(rd(17), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

            if (w.eq.0) return

         else ! 2 -> 3
c Generate y, tau -> x1, x2, shat from RD(1) and RD(2) first **************
            if (sqrt(xm2(6)).le.100.0d0) then
               tauc = 0.0008d0
            else if (sqrt(xm2(6)).le.180.0d0) then
               tauc = 0.0004d0
            else
               tauc = 0.0008d0
            endif
            A1 = 0.2d0/tauc
            A2 = -0.8d0/log(tauc)
            if ( RD(1) .le. 0.2d0 ) then
               tau = RD(1)/A1
               W = W/A1
            else
               tau = (1.0d0/tauc)**((RD(1)-1)/0.8d0)
               W = W/A2*tau
            endif

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

         endif !N

c generate phase space for decay of P -> V + V + V out of rd(3):rd(10)

         if (NW.eq.0) then
            ! in full calculation: narrow higgs resonance with separate phase space
            if (.not. ThreeBodyDecay3Higgs(rd(3), pY(0), q1(0), q2(0), q3(0), 
     1                    m1sq, m2sq, m3sq,W, NW)) return
         else
            ! in NW: higgs resonance would require off-shell Z, therefore no narrow
            !        resonance in process and default phase space works
            if (.not. ThreeBodyDecayPlain(rd(3), pY(0), q1(0), q2(0),
     1                     q3(0), m1sq, m2sq, m3sq, W, NW) ) return
         endif

         if (W.eq.0) return

c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
         call TwoBodyDecay0(rd(11),rd(12),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
         call TwoBodyDecay0(rd(13),rd(14),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
         call TwoBodyDecay0(rd(15),rd(16),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************
         if(N.EQ.1) then
            W=W
         else
            W = W * 1.d0/2d0/K(4)       ! Flux
            W = W * 3.89379304d11        !conversion into fb
         endif

      endif !ps_number

c************ END PHASE SPACE NUMBER 2 ***********************************


      else ! m_Higgs > 2 * m_Z
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     m_Higgs > 2 * m_Z -> Higgs resonances are broader
C     non-Higgs contributions and Higgs contributions are sampled
C     with one phase space generator
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(N.EQ.1) THEN ! 2 -> JET+Y

c generate q^2_Y which then decays into three massive particles from rd(1)
         if (.not. Resonance(xym**2, xym*xywidth, 1.0d-6, s,
     1                               rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(17), rd(18),rd(19)
         call TwoToJetPlusX(rd(17), rn, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

         if (w.eq.0) return

      ELSE ! 2 -> Y
c Generate y, tau -> x1, x2, shat from RD(1) and RD(2) first **************
         tauc = 0.0001d0
         A1 = 0.2d0/tauc
         A2 = -0.8d0/log(tauc)
         if ( RD(1) .le. 0.2d0 ) then
            tau = RD(1)/A1
            W = W/A1
         else
            tau = (1.0d0/tauc)**((RD(1)-1)/0.8d0)
            W = W/A2*tau
         endif

         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border !(prop density)^-1 of y

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)

         feynmanx1 = x1
         feynmanx2 = x2

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

      ENDIF


c generate phase space for decay of P -> V + V + V out of rd(3):rd(10)
c also in NW
         if (.not. ThreeBodyDecayPlain(rd(3), pY(0), q1(0), q2(0),
     1                     q3(0), m1sq, m2sq, m3sq, W, NW) ) return

         if (W.eq.0) return


c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(11),rd(12),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(13),rd(14),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(15),rd(16),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************
      if(N.EQ.1) then
         W=W
      else
         W = W * 1.d0/2d0/K(4)       ! Flux
         W = W * 3.89379304d11        !conversion into fb
      endif

      ENDIF ! selection of phase space according to m_Higgs

c***************************** DEBUGGING *********************************

      if (ldebug) then
         print*, "--------- new event -----------"
         print*, ""
         print*, "random numbers"
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         print*, rd(10), rd(11), rd(12)
         print*, rd(13), rd(14), rd(15)
         if (N.eq.1) print*, rd(16), rd(17), rd(18)
         print*, "x1    = ", x1
         print*, "x2    = ", x2
         print*, "wgt   = ", w
         print*, "---------------------------"
         print*, "--- momenta conservation ? ---"
         do mu=0,3
            q(mu,1) = q1(mu)
            q(mu,2) = q2(mu)
            q(mu,3) = q3(mu)
            test = abs(k1(mu)+k2(mu)-Pjet(mu)-q(mu,1)-q(mu,2)-q(mu,3))
            if (test.gt.1E-9) print*,"k1+k2-q1-q2-q3 = ", test
            test = abs(q(mu,1)-v(mu,1)-v(mu,2))
            if (test.gt.1E-9) print*,"q1-v1-v2 = ", test
            test = abs(q(mu,2)-v(mu,3)-v(mu,4))
             if (test.gt.1E-9) print*,"q2-v3-v4 = ", test
            test = abs(q(mu,3)-v(mu,5)-v(mu,6))
             if (test.gt.1E-9) print*,"q3-v5-v6 = ", test
            test = abs(k1(mu)+k2(mu)-Pjet(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6))
             if (test.gt.1E-9) print*,"k1+k2-v1-..-v6 = ", test
         enddo
         print*, "--- mass equal to zero ? ---"
         test = dotrr(k1(0),k1(0))
         if (test.gt.1E-7) print*,"p1^2 = ",dotrr(k1(0),k1(0))
         test = dotrr(k2(0),k2(0))
         if (test.gt.1E-7) print*,"p2^2 = ",dotrr(k2(0),k2(0))
         test = dotrr(Pjet(0),Pjet(0))
         if (test.gt.1E-7) print*,"Pjet^2 = ",dotrr(Pjet(0),Pjet(0))
         test = dotrr(v(0,1),v(0,1))
         if (test.gt.1E-7) print*,"v1^2 = ",dotrr(v(0,1),v(0,1))
         test = dotrr(v(0,2),v(0,2))
         if (test.gt.1E-7) print*,"v2^2 = ",dotrr(v(0,2),v(0,2)) 
         test = dotrr(v(0,3),v(0,3))
         if (test.gt.1E-7) print*,"v3^2 = ",dotrr(v(0,3),v(0,3))
         test = dotrr(v(0,4),v(0,4))
         if (test.gt.1E-7) print*,"v4^2 = ",dotrr(v(0,4),v(0,4))
         test = dotrr(v(0,5),v(0,5))
         if (test.gt.1E-7) print*,"v5^2 = ",dotrr(v(0,5),v(0,5))
         test = dotrr(v(0,6),v(0,6))
         if (test.gt.1E-7) print*,"v6^2 = ",dotrr(v(0,6),v(0,6))
         print*, "-----------------------"
      endif

c*********************** END OF DEBUGGING *******************************

      return
      end

