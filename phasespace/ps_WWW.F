c
c phase space generator for pp --> V1 + V2 + V3 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4, V3 --> l5 l6
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified: 

      SUBROUTINE DEF_PS_WWW(bos1,bos2,bos3)
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
      DOUBLE PRECISION RMLIMIT
      COMMON /TRIBOSWWWBLIPSI/ RMLIMIT(2)

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      double precision M2_i(max_res),MG_i(max_res),F_i(max_res)
      integer N_res

      DOUBLE PRECISION xxwidth, xxm
!      PARAMETER ( xxwidth = 200d0, xxm = 170d0 )
      PARAMETER ( xxwidth = 300d0, xxm = 170d0 )

      s = ecm**2

      bos(1) = bos1        ! W+ = 3
      bos(2) = bos2        ! W- = 4
      bos(3) = bos3        ! Z  = 2

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

      rm2minx = 1d-3          ! m_12
      rm2maxx = s/2D0

      rm2miny = 1d-3          ! m_23
      rm2maxy = s/2D0

      ! select relevant resonances for VV invariant mass
      N_res = 0
      if (xm2(6) .lt. 600d0**2) then
        N_res = N_res + 1
        M2_i(N_res) = xm2(6)
        MG_i(N_res) = xmg(6)
        F_i(N_res)  = 0.4
      endif
      call MultiResonancesInit(1, xxm**2, xxm*xxwidth, 1d0, M2_i, MG_i, F_i, N_res)


c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10.0d0)

 60   format ( ' machine energy = ', f10.1, '  GeV' )
 70   format ( ' integration range for boson ',i3,' = ',f10.3,
     1         '  to ',f12.3,' GeV' )
 72   format ( ' boson has mass ',f8.2,' GeV and width ',
     1         f10.4,' GeV' )
      if (.not.ldoblha) then
        write(6,*)
        write(6,*) '--------------------PHASE SPACE--------------------'
        write(6,*)
        do i = 1,3
           write(6,72) sqrt(rm2(i)),rmg(i)/sqrt(rm2(i))
           write(6,70) bos(i),sqrt(rm2min(i)),sqrt(rm2max(i))
        enddo
        write(6,*)
        write(6,60) sqrt(s)
        write(6,*)
        write(6,*) '---------------------------------------------------'
        write(6,*)
      endif

      end

C**********************  LIPSWWW   ****************************************
C
      SUBROUTINE LIPSWWW(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities//global.inc"

c
      INTEGER N,NRANDL,NW
      DOUBLE PRECISION TPI
      PARAMETER (nrandl=18, TPI=2D0*PI )
      DOUBLE PRECISION RD(nrandl),RN
      DOUBLE PRECISION K1(0:3), K2(0:3), V(0:3,6), X1, X2, W
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
C                      Q1(*) = V(*,1) + V(*,2)
C                      Q2(*) = V(*,3) + V(*,4)
C                      Q3(*) = V(*,5) + V(*,6)
C           W          the weight of the phase space point
C

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
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

      INTEGER I, MU

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION q1(0:3), q2(0:3), q3(0:3), K(0:4),pY(0:4), pX(0:4)
      DOUBLE PRECISION tau, y, y_border, tauc

      DOUBLE PRECISION dotrr, test
      DOUBLE PRECISION A1, A2
      DOUBLE PRECISION vdump, phi, cosphi, sinphi
      DOUBLE PRECISION rt

      PARAMETER ( xxwidth = 200d0, xxm = 170d0 )
      PARAMETER ( xywidth = 200d0, xym = 170d0 )


      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL RESONANCE, TwoBodyDecay, ThreeBodyDecayWWW,ThreeBodyDecay3Higgs
      EXTERNAL RESONANCE, TwoBodyDecay, ThreeBodyDecayWWW,ThreeBodyDecay3Higgs

c--------------------------------------------------------------------

      if (ldebug) print*,' enter LIPSWWW '

      W = 1d0             !init of weightfactor

c------------------------------------------------------------------------

      IF (NW.eq.1) THEN
c q^2 of W1 in narrow width approximation
         m1sq = rm2(1)
         W = W * rmg(1)/2.0d0

c q^2 of W2 in narrow width approximation
         m2sq = rm2(2)
         W = W * rmg(2)/2.0d0

c q^2 of W3 in narrow width approximation
         m3sq = rm2(3)
         W = W * rmg(3)/2.0d0

c pX^2 of intermediate particle X -> W1 + W2 out of rd(3)
         rm2max12 = (sqrt(s) - sqrt(m3sq))**2
         rm2min12 = (sqrt(m1sq) + sqrt(m2sq))**2
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min12, rm2max12,
     1                              rd(3), W, pX(4))) return
      endif !NW

      IF(N.EQ.1) THEN ! 2 -> JET+Y

         IF (NW.eq.1) THEN

c q^2_Y generated out of rd(2)
            rm2min23 = (sqrt(pX(4)) + sqrt(m3sq))**2
            if (.not. Resonance(xym**2, xym*xywidth, rm2min23, rm2max23,
     1                              rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(16), rd(17),rd(18)
            call TwoToJetPlusX(rd(16), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

            if (w .eq. 0.d0) return

         ELSE !NW.ne.1

c generate q^2_Y which then decays into three massive particles from rd(1)
            if (.not. Resonance(xym**2, xym*xywidth, 1.0d-6, s,
     1                               rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(16), rd(17),rd(18)
            call TwoToJetPlusX(rd(16), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

            if (w.eq.0) return

      ENDIF !NW

      ELSE ! 2 -> Y
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

! c generate x1 and x2
!       x1 = RD(1)
!       x2 = RD(2)
!       K(4) = tau*s

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

      IF (NW.eq.1) THEN

c generate phase space for decay of Y -> X + Z out of rd(17), rd(18)
         if (.not. TwoBodyDecay(rd(6), rd(7), pY(0), pY(4),
     1                 pX(4), m3sq, pX(0), q3(0), W)) return

c generate phase space for decay of X -> W+ + W- out of rd(4), rd(5)
         if (.not. TwoBodyDecay(rd(4), rd(5), pX(0), pX(4), m1sq, m2sq,
     &                            q1(0),q2(0),W)) return

      ELSE !NW.ne.1

c generate phase space for decay of P -> V + V + V out of rd(3):rd(9), rn
         IF(N.EQ.1) THEN ! 2 -> JET+Y
            rt = rd(2)
         ELSE
            rt = 0d0
         ENDIF
         if (.not. ThreeBodyDecayWWW(rd(3), rt, pY(0), q1(0), q2(0), q3(0), 
     1                                       m1sq, m2sq, m3sq,W, NW)) return
!          if (.not. ThreeBodyDecay3Higgs(rd(2),pY(0), q1(0), q2(0), q3(0), 
!      1                                       m1sq, m2sq, m3sq,W, NW)) return
         ! ThreeBodyDecayWWW can do NW as well, but invariant mass limits
         ! not set up properly...
      ENDIF

c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(10),rd(11),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(12),rd(13),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(14),rd(15),q3(0),m3sq,v(0,5),v(0,6),W)

c phi angle rotation
c  phasespace factors of 2*pi already added in ThreeBodyDecayWWW and 
c  TwoToJetPlusX for WWW and WWW+Jet, respectively

      phi = rn * TPI
      sinphi = sin(phi)
      cosphi = cos(phi)
      do i=1,6
         vdump  = v(1,i)
         v(1,i) = cosphi*v(1,i) + sinphi*v(2,i)
         v(2,i) = cosphi*v(2,i) - sinphi*vdump
      enddo

      if (N.eq.1) then
         vdump   = pjet(1)
         pjet(1) = cosphi*pjet(1) + sinphi*pjet(2)
         pjet(2) = cosphi*pjet(2) - sinphi*vdump
      endif


c Additional weight factors ******************************************
      if(N.EQ.1) then
         W=W
      else
         W = W * 1.d0/2d0/K(4)       ! Flux
         W = W * 3.89379304d11        !conversion into fb
      endif
c-----------------------------------------------------------
c end of phase space generation
c ----------------------------------------------------------
c ----------------------------------------------------------

      if (ldebug) then
         print*, '--------- new event -----------'
         print*, ''
         print*, 'random numbers'
         print*, rn
         print*, rd(1), rd(2), rd(3)
         print*, rd(4), rd(5), rd(6)
         print*, rd(7), rd(8), rd(9)
         print*, rd(10), rd(11), rd(12)
         print*, rd(13), rd(14), rd(15)
         print*, '---------------------------'

         print*,'m1sq + m2sq + m3sq + M^2 = m12sq + m23sq + m13sq ?'
         test = m1sq+m2sq+m3sq+x1*x2*s-m12sq-m23sq-2d0*dotrr(q1(0),q3(0))
     1         -dotrr(q1(0),q1(0))-dotrr(q3(0),q3(0))
         print*, test
         print*,'m1sq = q1^2 ?'
         test = m1sq-dotrr(q1(0),q1(0))
         print*, test
         print*,'m2sq = q2^2 ?'
         test = m2sq-dotrr(q2(0),q2(0))
         print*, test
         print*,'m3sq = q3^2 ?'
         test = m3sq-dotrr(q3(0),q3(0))
         print*, test
         print*, '--- momenta conservation ? (jet from RE not included!) ---'
         do mu=0,3
            test = k1(mu)+k2(mu)-q1(mu)-q2(mu)-q3(mu)
            if (test.gt.1E-9) print*, test
            test = q1(mu)-v(mu,1)-v(mu,2)
            if (test.gt.1E-9) print*, test
            test = q2(mu)-v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*, test
            test = q3(mu)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*, test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*, test
         enddo
         print*,'x1*x2*s = (k1+k2)^2'
         test = x1*x2*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
         print*, test
         print*,'---x1*x2*s = (q1 + q2 + q3)^2 ?---'
         test = x1*x2*s - dotrr(q1(0),q1(0)) - dotrr(q2(0),q2(0))
     &     - dotrr(q3(0),q3(0)) - 2*dotrr(q1(0),q2(0))
     &     - 2*dotrr(q1(0),q3(0)) - 2*dotrr(q2(0),q3(0))
         print*, test
         print*, '--- mass equal to zero ? ---'
         print*,'p1^2 = ',dotrr(k1(0),k1(0))
         print*,'p2^2 = ',dotrr(k2(0),k2(0))
         print*,'v1^2 = ',dotrr(v(0,1),v(0,1))
         print*,'v2^2 = ',dotrr(v(0,2),v(0,2))
         print*,'v3^2 = ',dotrr(v(0,3),v(0,3))
         print*,'v4^2 = ',dotrr(v(0,4),v(0,4))
         print*,'v5^2 = ',dotrr(v(0,5),v(0,5))
         print*,'v6^2 = ',dotrr(v(0,6),v(0,6))
         print*, '-----------------------'

         ! stop
      endif

      return
      end

