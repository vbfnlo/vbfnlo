cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c					               		            c
c	This file contains subroutines for the phasespace of                c
c	Triboson production in the KK scenario				    c
c									    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c
c phase space generator for pp --> V1 + V2 + V3 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4, V3 --> l5 l6
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified: 

      SUBROUTINE DEF_PS_WWZ_kk(bos1,bos2,bos3)
          use globalvars, only: ldoblha

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
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
*	print*,' enter LIPSWWZ_kk ' !MB DEBUG
      s = ecm**2

      bos(1) = bos1        ! W+ = 3
      bos(2) = bos2        ! W- = 4
      bos(3) = bos3        ! Z  = 2

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.2) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3         ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0        ! for tan mapping
         endif
      enddo

      rm2minx = 1d-3
      rm2maxx = s/2d0

      rm2miny = 1d-3
      rm2maxy = s/2d0

c      rm2max1 = (sqrt(xm2(6)-30d0*xmg(6)))**2 !Dieters proposal
c      rm2max2 = (sqrt(xm2(6)+30d0*xmg(6)))**2 !Dieters proposal

c coupling routine handles over 10 times width times mass

      rm2max1 = ((xm2(6)-0.4d0*xmg(6))/sqrt(xm2(6)) )**2
      rm2max2 = ((xm2(6)+0.7d0*xmg(6))/sqrt(xm2(6)) )**2

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

C**********************  LIPSWWZ_kk   ****************************************
C 
      SUBROUTINE LIPSWWZ_kk(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"

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

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq, rt
      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION q1(0:3), q2(0:3), q3(0:3), K(0:4),pY(0:4)
      DOUBLE PRECISION tau, y, y_border

      DOUBLE PRECISION dotrr, test, vdummy(0:3,2)
	Double precision vdump, phi, cosphi, sinphi
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c LO parameter for WWZ inv mass
      PARAMETER ( xxwidth = 300d0, xxm = 490d0 )
      PARAMETER ( xywidth = 6000d0, xym = 480d0 )
c NLO real emission parameter for WWZ invariant mass

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_wwz_kk
      EXTERNAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_wwz_kk
	logical nonresonance, nonresonance_kk
	external nonresonance, nonresonance_kk
c symmetrisation of q1 and q3
      INTEGER symm

      data symm /1/
      save symm    
c DEBUG	
	logical ampout
	double precision mom(0:4)
	common/MB3/ m12sq, m23sq
	common/MB/ ampout

	Double precision expo3v1, m3v1, g3v1
	Double precision expo3v2, m3v2, g3v2		
	PARAMETER (expo3v1=1.5d0, m3v1=200d0, g3v1=850d0)
	PARAMETER (expo3v2=4.26d0, m3v2=1.9d0, g3v2=55500d0)
c--------------------------------------------------------------------

      if (ldebug) print*,' enter LIPSWWZ_kk '
      W = 1d0             !init of weightfactor

c------------------------------------------------------------------------

      IF(N.EQ.1) THEN ! 2 -> JET+Y

c generate q^2_Y which then decays into three massive particles from rd(1)


c	if (rd(1).le. 0.6d0) then
c	      if (.not. Resonance(520d0**2, 520d0*300d0, 
c     1	 1d-6, (1.1*sqrt(xm2(6)))**2, rd(1)/0.6d0, W, pY(4))) return
c	      W=W/0.6d0
c	else
c      	if (.not. Resonance(1550d0**2, 1550d0*1560d0, 
c     1	 (1.1*sqrt(xm2(6)))**2, s, (rd(1)-0.6d0)/0.4d0, W, pY(4))) return
c		W=W/0.4d0
c	endif

	if (rd(1).le. 0.05) then
		if(.not. nonresonance(1d-6,250d0**2,rd(1)/0.05d0,W,pY(4)))return
		W=W/0.05
	else
		if(.not. NonResonance_kk(156d0**2, 156d0*880d0, 2d0, 250d0**2,
     &	   s, (rd(1)-0.05d0)/0.95, W, pY(4))) return
		W=W/0.95
	endif

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(16), rd(17),rd(18)
         call TwoToJetPlusX(rd(16), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

         if (w.eq.0) return

      ELSE ! 2 -> Y

c generate m_wwz^2 out of rd(1) via tan mapping and transform it to tau
	if (rd(1) .le. 0.5) then
         if (.not. Resonance(xxm**2, xxwidth*xxm, 1d-6,
     1 (1.2d0*sqrt(xm2(6)))**2, rd(1)/0.5d0, W, K(4))) return
         W = W/0.5d0*2d0*pi
      else 
         if (.not. Resonance(xym**2, xywidth*xym, 
     1 (1.2d0*sqrt(xm2(6)))**2, s, (rd(1)-0.5d0)/0.5d0, W, K(4))) return
         W = W/0.5d0*2d0*pi
      endif	

c generate tau= K(4)/s
		tau=k(4)/s
		W = W/s
c generate y out of rd(2)

         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border

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

c assign jet momenta
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0

         do mu =0,4
            pY(mu) = K(mu)
         enddo

      ENDIF

		if (n .eq. 0) then
		rt=0d0
		else
		rt=rd(2)
		endif
      
c generate phase space for decay of P -> V + V + V out of rd(3):rd(9), rn
         if (.not. ThreeBodyDecay_wwz_kk(rd(3), rt, pY(0), q1(0), q2(0), q3(0), 
     1                                       m1sq, m2sq, m3sq,W)) return
      

c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(10),rd(11),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(12),rd(13),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(14),rd(15),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************


	phi=twopi*rn
	sinphi=sin(phi)
	cosphi=cos(phi)
	do i=1,6
		vdump=v(1,i)
		v(1,i)=cosphi*v(1,i)+sinphi*v(2,i)
		v(2,i)=cosphi*v(2,i)-sinphi*vdump
	enddo
	if (n .eq. 1) then
		vdump=pjet(1)
		pjet(1)=cosphi*pjet(1)+sinphi*pjet(2)
		pjet(2)=cosphi*pjet(2)-sinphi*vdump
	endif

      if(N.EQ.1) then
         W=W
      else
         W = W * 1.d0/2d0/K(4)       ! Flux
         W = W * 3.89379304d11        !conversion into fb
      endif

c exchange v3,v4 with v5,v6 for Matrixelement in WWZ kk scenario

	do i=0,3
	vdummy(i,1)=v(i,3)
	vdummy(i,2)=v(i,4)
	v(i,3)=v(i,5)
	v(i,4)=v(i,6)
	v(i,5)=vdummy(i,1)
	v(i,6)=vdummy(i,2)
	enddo


c-----------------------------------------------------------
c end of phase space generation
c ----------------------------------------------------------
c ----------------------------------------------------------

      if (ldebug) then
c         print*, '--------- new event -----------'
c         print*, ''
c         print*, 'random numbers'
c         print*, rn
c         print*, rd(1), rd(2), rd(3)
c         print*, rd(4), rd(5), rd(6)
c         print*, rd(7), rd(8), rd(9)
c         print*, rd(10), rd(11), rd(12)
c         print*, rd(13), rd(14), rd(15)
c         print*, '---------------------------'
c
c         print*,'m1sq + m2sq + m3sq + M^2 = m12sq + m23sq + m13sq ?'
c         test = m1sq+m2sq+m3sq+x1*x2*s-m12sq-m23sq-2d0*dotrr(q1(0),q3(0))
c     1         -dotrr(q1(0),q1(0))-dotrr(q3(0),q3(0))
c         print*, test
c	print*,"******************************************"
c	if ((sqrt(m12sq) .le. (sqrt(m1sq)+sqrt(m2sq))) .and. 
c     & (sqrt(m23sq) .le. (sqrt(m2sq)+sqrt(m3sq)))) print*,
c     & "XXXXXXXX  ERROR XXXXXXX" 
c	   
c         print*,'m1sq = q1^2 ?'
c         test = m1sq-dotrr(q1(0),q1(0))
c         print*, test
c         print*,'m2sq = q2^2 ?'
c         test = m2sq-dotrr(q2(0),q2(0))
c         print*, test
c         print*,'m3sq = q3^2 ?'
c         test = m3sq-dotrr(q3(0),q3(0))
c         print*, test
         print*, '--- momenta conservation ? ---'
         do mu=0,3
            test = k1(mu)+k2(mu)-q1(mu)-q2(mu)-q3(mu)
            if (test.gt.1E-9) print*,"1", test
            test = q1(mu)-v(mu,1)-v(mu,2)
           if (test.gt.1E-9) print*,"2", test
            test = q2(mu)-v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*,"3", test
            test = q3(mu)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*,"4", test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*,"5",mu, test
         enddo
c         print*,'x1*x2*s = (k1+k2)^2'
c         test = x1*x2*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
c         print*, test
c         print*,'---x1*x2*s = (q1 + q2 + q3)^2 ?---'
c         test = x1*x2*s - dotrr(q1(0),q1(0)) - dotrr(q2(0),q2(0))
c     &     - dotrr(q3(0),q3(0)) - 2*dotrr(q1(0),q2(0))
c     &     - 2*dotrr(q1(0),q3(0)) - 2*dotrr(q2(0),q3(0))
c         print*, test
c         print*, '--- mass equal to zero ? ---'
c         print*,'p1^2 = ',dotrr(k1(0),k1(0))
c         print*,'p2^2 = ',dotrr(k2(0),k2(0))
c         print*,'v1^2 = ',dotrr(v(0,1),v(0,1))
c         print*,'v2^2 = ',dotrr(v(0,2),v(0,2))
c         print*,'v3^2 = ',dotrr(v(0,3),v(0,3))
c         print*,'v4^2 = ',dotrr(v(0,4),v(0,4))
c         print*,'v5^2 = ',dotrr(v(0,5),v(0,5))
c         print*,'v6^2 = ',dotrr(v(0,6),v(0,6))
c         print*, '-----------------------'
      endif

      return
      end

***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay_wwz_kk(R, RN, P_in, q1_out, q2_out, 
     &                  q3_out,m1sq, m2sq, m3sq, w)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(8), rn       : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecay         : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE
c to get all resonances, include:
#include "VBFNLO/utilities/kkcoupl.inc"
      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)

      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(7), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 40d1, xxm = 28d1 )
      DOUBLE PRECISION xywidth, xym
      PARAMETER ( xywidth = 180d1, xym = 85d1 )

      LOGICAL RESONANCE, nonresonance_kk, nonresonance
      EXTERNAL RESONANCE, nonresonance_kk, nonresonance

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION QL1(0:3), QL3(0:3), msqlocal1, msqlocal3
      data symm /1/
      save symm    

c DEBUG
	common/MB3/ m12sq, m23sq	
	double precision m2v1,g2v1,expo2v1
	parameter (m2v1=118d0, g2v1=235d0, expo2v1=2d0)
      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first

c m_12^2 = (q1 + q2)^2 out of R(1)

      if ( r(1) .le. 0.1d0 ) then
         if (.not. Resonance(150d0**2, 150d0*80d0, rm2min12, 150d0**2d0,
     1                              r(1)/0.1d0, W, m12sq)) return
         W = W/0.1d0

	else if ( r(1) .le.  0.55d0 ) then
	   if (.not. nonResonance_kk( m2v1**2, m2v1*g2v1, expo2v1, 
     1 	150d0**2, rm2max1, (r(1)-0.1d0)/0.45d0, W, m12sq) ) return
	   W = W/0.45d0

      else if ( r(1) .le.  0.8d0 ) then
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*w1width_k(1),
     1	      rm2max1, rm2max2,(r(1)-0.55d0)/0.25d0, W, m12sq)) return
         W = W/0.25d0	

      else
         if (.not. Resonance(xym**2, xym*xywidth, rm2max2, rm2max12,
     1                      (r(1)-0.8d0)/0.2d0, W, m12sq)) return
         W = W/0.2d0
      endif

      if (m12sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c m_23^2 = (q2 + q3)^2 out of r(2)

      if ( r(2) .le. 0.1d0 ) then
         if (.not. Resonance(150d0**2, 80d0*150d0,rm2min12, 150d0**2d0,
     1                              r(2)/0.1d0, W, m23sq)) return
         W = W/0.1d0

	else if ( r(2) .le.  0.55d0 ) then
	   if (.not. nonResonance_kk( m2v1**2, m2v1*g2v1, expo2v1, 
     1 	150d0**2, rm2max1, (r(2)-0.1d0)/0.45d0, W, m23sq) ) return
	   W = W/0.45d0

      else if ( r(2) .le.  0.8d0 ) then
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*w1width_k(1),
     1 	rm2max1, rm2max2,(r(2)-0.55d0)/0.25d0, W, m23sq)) return
         W = W/0.25d0	

      else
         if (.not. Resonance(xym**2, xym*xywidth, rm2max2, rm2max23,
     1                      (r(2)-0.8d0)/0.2d0, W, m23sq)) return
         W = W/0.2d0
      endif

      if (m23sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c now generate m1^2 and m3^2 of the massive vector bosons **********

c m_1^2 is generated via tan mapping out of r(3)
      rm2max(1) = min(p_in(4)-m23sq,p_in(4)+m23sq-2.0d0*sqrt(p_in(4)*m23sq))
      if (rm2max(1).le.0) then
         W = 0.0d0
         return
      endif

         if (.not. Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1),
     1                              r(3), W, m1sq)) return


c m_3^2 is generated via tan mapping out of r(4)
      rm2max(3) = min(p_in(4)-m12sq, p_in(4)+m12sq-2.0d0*sqrt(p_in(4)*m12sq), 
     1                m12sq+m23sq-m1sq)
      if (rm2max(3).le.0) then
         W = 0.0d0
         return
      endif
      if (.not. Resonance(rm2(3), rmg(3), rm2min(3), rm2max(3),
     1                              r(4), W, m3sq)) return

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2/M
      E3 = (p_in(4) + m3sq - m12sq)/2/M
      E2 = M - E1 - E3
      Abs1 = sqrt(E1**2 - m1sq)
      Abs3 = sqrt(E3**2 - m3sq)

c m_2^2 is generated via tan mapping out of r(5)
      rm2min(2) = max(-Abs1**2-Abs3**2-2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1            1.0d-5)
      rm2max(2) = min(-Abs1**2-Abs3**2+2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1          E2**2, P_in(4))
      if (rm2max(2).le.0) then
         W = 0.0d0
         return
      endif
      
         if (.not. Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2),
     1                                 r(5), W, m2sq)) return
      

c now generate Euler angles
      alpha = R(6) * twopi
      W = W*twopi

      cosb = -1.0d0 + R(7)*2.0d0
      W = W*2.0d0

      gamma = RN * twopi
      W = W*twopi

      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      sinb = sqrt(1 - cosb**2)
      sing = sin(gamma)

c calculate now C and D
      D = -(m1sq-m2sq+m3sq+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)/2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
         print*,'ps_tribos_kk: D**2 >1'
         return
      endif

c Momenta of V1
      ql1(0) = E1
      ql1(1) = Abs1 * (cosa * cosg - cosb * sina * sing)
      ql1(2) = Abs1 * (-cosb * cosg * sina - cosa * sing)
      ql1(3) = Abs1 * (sina * sinb)

c Momenta of V3
      ql3(0) = E3
      ql3(1) = Abs3 * C * (cosg * sina + cosa * cosb * sing) 
     1      - Abs3 * D * (cosa * cosg - cosb * sina * sing)
      ql3(2) = Abs3 * C * (cosa * cosb * cosg - sina * sing)
     1      - Abs3 * D * (-cosb * cosg * sina - cosa * sing)
      ql3(3) = Abs3 * C * (- cosa * sinb) - Abs3 * D * (sina * sinb)

      call boostn( ql1(0), P_in(0), ql1(0) )
      call boostn( ql3(0), P_in(0), ql3(0) )

c Momenta of V2
      q2_out(0) = P_in(0) - ql1(0) - ql3(0)
      q2_out(1) = P_in(1) - ql1(1) - ql3(1)
      q2_out(2) = P_in(2) - ql1(2) - ql3(2)
      q2_out(3) = P_in(3) - ql1(3) - ql3(3)


c symmetrisation of q1 and q3
      msqlocal1 = m1sq
      msqlocal3 = m3sq
      symm = mod(symm,2) + 1
      if (symm.eq.1) then !interchange q1 <-> q3
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (symm.eq.2) then
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q3_out(mu) = ql1(mu)
         enddo
         m1sq = msqlocal3
         m3sq = msqlocal1
      endif


c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      ThreeBodyDecay_wwz_kk = .true.

      RETURN
      END


c----------------------------------------------------------------------------
c phase space generator for pp --> V1 + V2 + V3 + 1 Jet
c                                  V1 --> l1 l2, V2 --> l3 l4, V3 --> l5 l6
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c       Initial version: January 2008
c       Last modified: 

      SUBROUTINE DEF_PS_WWW_kk(bos1,bos2,bos3)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"
c to get all resonances, include:
#include "VBFNLO/utilities/kkcoupl.inc"

      INTEGER bos(3), i, bos1, bos2, bos3

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN, x
      DOUBLE PRECISION RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      DOUBLE PRECISION RM2MAX1, RM2MAX2, RM2MAX3, RM2MAX4
	DOUBLE PRECISION XXM1,XXM2,XXW1,XXW2
      COMMON /TRIBOSBLIPSI_WWW/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2, RM2MAX3, RM2MAX4,
     3		    XXM1,XXM2,XXW1,XXW2

      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)
*	print*,' enter LIPSWWW_kk ' !MB DEBUG
      s = ecm**2

      bos(1) = bos1        ! W+ = 3
      bos(2) = bos2        ! W- = 4
      bos(3) = bos3        ! Z  = 2
      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.2) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3         ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0        ! for tan mapping
         endif
      enddo
      rm2minx = 1d-3
      rm2maxx = s/2d0

      rm2miny = 1d-3
      rm2maxy = s/2d0
c set up resonance masses, widths and the range of the resonances for 
c tan mapping	
      if (mhlmswitch .or. (kkha.eq.0)) then ! both resonances in 
                                            ! ps-generator are mapped onto Z'
         rm2max1 = (z1mass_k(1)-6d0*z1width_k(1))**2
         rm2max3 = (z1mass_k(1)-6d0*z1width_k(1))**2
         rm2max2 = (z1mass_k(1)+9d0*z1width_k(1))**2
         rm2max4 = (z1mass_k(1)+9d0*z1width_k(1))**2
         xxm1=z1mass_k(1)
         xxw1=z1width_k(1)
         xxm2=xxm1
         xxw2=xxw1
      else                      ! one mapped onto Z_1 and one onto A_1
         if (z1mass_k(1).lt. a1mass_k(1)) then
            xxm1=z1mass_k(1)
            xxm2=a1mass_k(1)
            xxw1=z1width_k(1)
            xxw2=a1width_k(1)
            x=(a1mass_k(1)-z1mass_k(1))/(a1width_k(1)+z1width_k(1))
            rm2max1 = (z1mass_k(1)-3d0*z1width_k(1))**2
            rm2max2 = ( z1mass_k(1)+x*z1width_k(1) )**2
            rm2max3 =rm2max2
            rm2max4 = (a1mass_k(1)+6d0*a1width_k(1))**2
         else 
            xxm1=a1mass_k(1)
            xxm2=z1mass_k(1)
            xxw1=a1width_k(1)
            xxw2=z1width_k(1)
            x=(z1mass_k(1)-a1mass_k(1))/(a1width_k(1)+z1width_k(1))
            rm2max1 = (a1mass_k(1)-3d0*z1width_k(1))**2
            rm2max2 = (a1mass_k(1)+x*a1width_k(1) )**2
            rm2max3 =rm2max2
            rm2max4 = (z1mass_k(1)+6d0*a1width_k(1))**2
         endif
      endif
      
*	print*, z1mass_k(1), z1width_k(1)	
*	print*, xxm1, xxw1
*	print*, xxm2, xxw2
	

c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10.0d0)

 60   format ( ' machine energy = ', f10.1, '  GeV' )
 70   format ( ' integration range for boson ',i3,' = ',f10.3,
     1         '  to ',f12.3,' GeV' )
 72   format ( ' boson has mass ',f8.2,' GeV and width ',
     1         f10.4,' GeV' )
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

      end

C**********************  LIPSWWW_kk   ****************************************
C 
      SUBROUTINE LIPSWWW_kk(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"

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
      DOUBLE PRECISION RM2MAX1, RM2MAX2,RM2MAX3, RM2MAX4
	DOUBLE PRECISION XXM1,XXM2,XXW1,XXW2
      COMMON /TRIBOSBLIPSI_WWW/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2, RM2MAX3, RM2MAX4,
     3		    XXM1,XXM2,XXW1,XXW2

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

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq, rt
      DOUBLE PRECISION xxwidth, xxm, xywidth, xym
      DOUBLE PRECISION xzwidth, xzm, xwwidth, xwm
      DOUBLE PRECISION q1(0:3), q2(0:3), q3(0:3), K(0:4),pY(0:4)
      DOUBLE PRECISION tau, y, y_border

      DOUBLE PRECISION dotrr, test
	Double precision vdump, phi, cosphi, sinphi
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
c LO + virt parameters
      PARAMETER ( xxwidth = 220d0, xxm = 430d0 )
      PARAMETER ( xywidth = 2200d0, xym = 1100d0 )
c real parameters
      PARAMETER ( xwwidth = 280d0, xwm = 505d0 )
      PARAMETER ( xzwidth = 1373d0, xzm = 1440d0 )	

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_www_kk
      EXTERNAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_www_kk
	logical nonresonance_kk, nonresonance
	external nonresonance_kk, nonresonance

c DEBUG	
	logical ampout
	double precision mom(0:4)
	common/MB3/ m12sq, m23sq
	common/MB/ ampout
c--------------------------------------------------------------------

      if (ldebug) print*,' enter LIPSWWW '
      W = 1d0             !init of weightfactor

c------------------------------------------------------------------------

      IF(N.EQ.1) THEN ! 2 -> JET+Y


c generate q^2_Y which then decays into three massive particles from rd(1)
	if (rd(1) .le. 0.5d0) then
         if (.not. Resonance(xwm**2, xwm*xwwidth, 1.0d-6, (1.2d0*sqrt(xm2(6)))**2,
     1                               rd(1)/0.5d0, W, pY(4))) return
	W=W/0.5d0
	else
         if (.not. Resonance(xzm**2, xzm*xzwidth, (1.2d0*sqrt(xm2(6)))**2, s,
     1                               (rd(1)-0.5)/0.5d0, W, pY(4))) return
	W=W/0.5d0
	endif

c	if (rd(1).le. 0.3d0)then
c	 if (.not. Resonance(330d0**2,330d0*120d0, 1d-6,
c     1   540d0**2, rd(1)/0.3,W,pY(4))) return
c	W=W/0.3
c	else
c         if (.not. Resonance(880d0**2, 880d0*1500d0, 540d0**2, s,
c     1                               (rd(1)-0.3)/0.7, W, pY(4))) return
c	W=W/0.7
c	endif
c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(16), rd(17),rd(18)
         call TwoToJetPlusX(rd(16), 0d0, ecm, pTc, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)


         if (w.eq.0) return

      ELSE ! 2 -> Y

c generate m_www^2 out of rd(1) via tan mapping
      if ( rd(1) .le. 0.6d0 ) then
         if (.not. Resonance(xxm**2, xxm*xxwidth, 0.001d0, (1.2d0*sqrt(xm2(6)))**2,
     1                              rd(1)/0.6d0, W, K(4))) return
         W = W/0.6d0
      else 
         if (.not. Resonance(xym**2, xywidth*xym, (1.2d0*sqrt(xm2(6)))**2, s,
     1                      (rd(1)-0.6d0)/0.4d0, W, K(4))) return
         W = W/0.4d0
      endif	
		
c generate tau= K(4)/s
		tau=k(4)/s
		W = W/s*2*pi

c generate y out of rd(2)
         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border

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

      ENDIF

		if (n .eq. 0) then
		rt=0d0
		else
		rt=rd(2)
		endif
      
c generate phase space for decay of P -> V + V + V out of rd(3):rd(9), rn
         if (.not. ThreeBodyDecay_www_kk(rd(3), rt, pY(0), q1(0), q2(0), q3(0), 
     1                                       m1sq, m2sq, m3sq,W)) return
      

c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(10),rd(11),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(12),rd(13),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(14),rd(15),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************


	phi=twopi*rn
	sinphi=sin(phi)
	cosphi=cos(phi)
	do i=1,6
		vdump=v(1,i)
		v(1,i)=cosphi*v(1,i)+sinphi*v(2,i)
		v(2,i)=cosphi*v(2,i)-sinphi*vdump
	enddo

	if (n .eq. 1) then
		vdump=pjet(1)
		pjet(1)=cosphi*pjet(1)+sinphi*pjet(2)
		pjet(2)=cosphi*pjet(2)-sinphi*vdump
	endif

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
c         print*, '--------- new event -----------'
c         print*, ''
c         print*, 'random numbers'
c         print*, rn
c         print*, rd(1), rd(2), rd(3)
c         print*, rd(4), rd(5), rd(6)
c         print*, rd(7), rd(8), rd(9)
c         print*, rd(10), rd(11), rd(12)
c         print*, rd(13), rd(14), rd(15)
c         print*, '---------------------------'
c
c         print*,'m1sq + m2sq + m3sq + M^2 = m12sq + m23sq + m13sq ?'
c         test = m1sq+m2sq+m3sq+x1*x2*s-m12sq-m23sq-2d0*dotrr(q1(0),q3(0))
c     1         -dotrr(q1(0),q1(0))-dotrr(q3(0),q3(0))
c         print*, test
c	print*,"******************************************"
c	if ((sqrt(m12sq) .le. (sqrt(m1sq)+sqrt(m2sq))) .and. 
c     & (sqrt(m23sq) .le. (sqrt(m2sq)+sqrt(m3sq)))) print*,
c     & "XXXXXXXX  ERROR XXXXXXX" 
c	   
c         print*,'m1sq = q1^2 ?'
c         test = m1sq-dotrr(q1(0),q1(0))
c         print*, test
c         print*,'m2sq = q2^2 ?'
c         test = m2sq-dotrr(q2(0),q2(0))
c         print*, test
c         print*,'m3sq = q3^2 ?'
c         test = m3sq-dotrr(q3(0),q3(0))
c         print*, test
         print*, '--- momenta conservation ? ---'
         do mu=0,3
            test = k1(mu)+k2(mu)-q1(mu)-q2(mu)-q3(mu)
            if (test.gt.1E-9) print*,"1", test
            test = q1(mu)-v(mu,1)-v(mu,2)
           if (test.gt.1E-9) print*,"2", test
            test = q2(mu)-v(mu,3)-v(mu,4)
             if (test.gt.1E-9) print*,"3", test
            test = q3(mu)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*,"4", test
            test = k1(mu)+k2(mu)-v(mu,1)-v(mu,2)
     &             -v(mu,3)-v(mu,4)-v(mu,5)-v(mu,6)
             if (test.gt.1E-9) print*,"5",mu, test
         enddo
c         print*,'x1*x2*s = (k1+k2)^2'
c         test = x1*x2*s - dotrr(k1,k1)-dotrr(k2,k2)-2*dotrr(k1,k2)
c         print*, test
c         print*,'---x1*x2*s = (q1 + q2 + q3)^2 ?---'
c         test = x1*x2*s - dotrr(q1(0),q1(0)) - dotrr(q2(0),q2(0))
c     &     - dotrr(q3(0),q3(0)) - 2*dotrr(q1(0),q2(0))
c     &     - 2*dotrr(q1(0),q3(0)) - 2*dotrr(q2(0),q3(0))
c         print*, test
c         print*, '--- mass equal to zero ? ---'
c         print*,'p1^2 = ',dotrr(k1(0),k1(0))
c         print*,'p2^2 = ',dotrr(k2(0),k2(0))
c         print*,'v1^2 = ',dotrr(v(0,1),v(0,1))
c         print*,'v2^2 = ',dotrr(v(0,2),v(0,2))
c         print*,'v3^2 = ',dotrr(v(0,3),v(0,3))
c         print*,'v4^2 = ',dotrr(v(0,4),v(0,4))
c         print*,'v5^2 = ',dotrr(v(0,5),v(0,5))
c         print*,'v6^2 = ',dotrr(v(0,6),v(0,6))
c         print*, '-----------------------'
      endif

      return
      end
***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay_www_kk(R, RN, P_in, q1_out, 
     &                  q2_out, q3_out, m1sq, m2sq, m3sq, w)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(8), rn       : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_(1.4d0*sqrt(xm2(6)))**2out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecay         : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE
c to get all resonances, include:
#include "VBFNLO/utilities/kkcoupl.inc"

      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      REAL*8 RM2MAX1, RM2MAX2, RM2MAX3, RM2MAX4
	DOUBLE PRECISION XXM1,XXM2,XXW1,XXW2
      COMMON /TRIBOSBLIPSI_WWW/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2, RM2MAX3, RM2MAX4,
     3		    XXM1,XXM2,XXW1,XXW2
      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)

      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(7), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 20d1, xxm = 18d1 )
      DOUBLE PRECISION xywidth, xym
      PARAMETER ( xywidth = 200d1, xym = 718d0 )

      LOGICAL RESONANCE, nonresonance, nonresonance_kk
      EXTERNAL RESONANCE, nonresonance, nonresonance_kk

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION QL1(0:3), QL3(0:3), msqlocal1, msqlocal3
      data symm /1/
      save symm    

c DEBUG
	common/MB3/ m12sq, m23sq	

      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first

c m_12^2 = (q1 + q2)^2 out of R(1)

      if ( r(1) .le. 0.01d0 ) then
         if (.not. nonResonance(rm2min12, 100d0**2d0,
     1                              r(1)/0.01d0, W, m12sq)) return
         W = W/0.01d0
	else if ( r(1) .le.  0.6d0 ) then
	   if (.not. nonResonance_kk( 98d0**2, 98d0*75d0, 12d0, 
     1 	100d0**2, rm2max1, (r(1)-0.01d0)/0.59d0, W, m12sq) ) return
	   W = W/0.59d0

      else if ( r(1) .le.  0.7d0 ) then
         if (.not. Resonance(xxm1**2, xxm1*xxw1,
     1		rm2max1, rm2max2,(r(1)-0.6d0)/0.1d0, W, m12sq)) return
         W = W/0.1d0
	   if (mhlmswitch) W=W/2d0

      else if ( r(1) .le.  0.8d0 ) then
         if (.not. Resonance(xxm2**2, xxm2*xxw2,
     1	      rm2max3, rm2max4,(r(1)-0.7d0)/0.1d0, W, m12sq)) return
         W = W/0.1d0	
	   if (mhlmswitch) W=W/2d0
      else
         if (.not. Resonance(xym**2, xym*xywidth, rm2max4, rm2max12,
     1                      (r(1)-0.8d0)/0.2d0, W, m12sq)) return
         W = W/0.2d0
      endif

      if (m12sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif


c m_23^2 = (q2 + q3)^2 out of r(2)

      if ( r(2) .le. 0.01d0 ) then
         if (.not. nonResonance(rm2min12, 100d0**2d0,
     1                              r(2)/0.01d0, W, m23sq)) return
         W = W/0.01d0

	else if ( r(2) .le.  0.6d0 ) then
	   if (.not. nonResonance_kk( 98d0**2, 98d0*75d0, 12d0, 
     1 	100d0**2, rm2max1, (r(2)-0.01d0)/0.59d0, W, m23sq) ) return
	   W = W/0.59d0

      else if ( r(2) .le.  0.7d0 ) then
         if (.not. Resonance(xxm1**2, xxm1*xxw1,
     1		rm2max1, rm2max2,(r(2)-0.6d0)/0.1d0, W, m23sq)) return
         W = W/0.1d0
	   if (mhlmswitch) W=W/2d0
      else if ( r(2) .le.  0.8d0 ) then
         if (.not. Resonance(xxm2**2, xxm2*xxw2,
     1 rm2max3, rm2max4,(r(2)-0.7d0)/0.1d0, W, m23sq)) return
         W = W/0.1d0	
	   if (mhlmswitch) W=W/2d0
      else
         if (.not. Resonance(xym**2, xym*xywidth, rm2max4, rm2max23,
     1                      (r(2)-0.8d0)/0.2d0, W, m23sq)) return
         W = W/0.2d0
      endif

      if (m23sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c now generate m1^2 and m3^2 of the massive vector bosons **********

c m_1^2 is generated via tan mapping out of r(3)
      rm2max(1) = min(p_in(4)-m23sq,p_in(4)+m23sq-2.0d0*sqrt(p_in(4)*m23sq))
      if (rm2max(1).le.0) then
         W = 0.0d0
         return
      endif

         if (.not. Resonance(rm2(1), rmg(1), rm2min(1), rm2max(1),
     1                              r(3), W, m1sq)) return


c m_3^2 is generated via tan mapping out of r(4)
      rm2max(3) = min(p_in(4)-m12sq, p_in(4)+m12sq-2.0d0*sqrt(p_in(4)*m12sq), 
     1                m12sq+m23sq-m1sq)
      if (rm2max(3).le.0) then
         W = 0.0d0
         return
      endif
      if (.not. Resonance(rm2(3), rmg(3), rm2min(3), rm2max(3),
     1                              r(4), W, m3sq)) return

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2/M
      E3 = (p_in(4) + m3sq - m12sq)/2/M
      E2 = M - E1 - E3
      Abs1 = sqrt(E1**2 - m1sq)
      Abs3 = sqrt(E3**2 - m3sq)

c m_2^2 is generated via tan mapping out of r(5)
      rm2min(2) = max(-Abs1**2-Abs3**2-2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1            1.0d-5)
      rm2max(2) = min(-Abs1**2-Abs3**2+2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1          E2**2, P_in(4))
      if (rm2max(2).le.0) then
         W = 0.0d0
         return
      endif
      
         if (.not. Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2),
     1                                 r(5), W, m2sq)) return
      

c now generate Euler angles
      alpha = R(6) * twopi
      W = W*twopi

      cosb = -1.0d0 + R(7)*2.0d0
      W = W*2.0d0

      gamma = RN * twopi
      W = W*twopi

      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      sinb = sqrt(1 - cosb**2)
      sing = sin(gamma)

c calculate now C and D
      D = -(m1sq-m2sq+m3sq+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)/2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
         print*,'D**2 >1'
         return
      endif

c Momenta of V1
      ql1(0) = E1
      ql1(1) = Abs1 * (cosa * cosg - cosb * sina * sing)
      ql1(2) = Abs1 * (-cosb * cosg * sina - cosa * sing)
      ql1(3) = Abs1 * (sina * sinb)

c Momenta of V3
      ql3(0) = E3
      ql3(1) = Abs3 * C * (cosg * sina + cosa * cosb * sing) 
     1      - Abs3 * D * (cosa * cosg - cosb * sina * sing)
      ql3(2) = Abs3 * C * (cosa * cosb * cosg - sina * sing)
     1      - Abs3 * D * (-cosb * cosg * sina - cosa * sing)
      ql3(3) = Abs3 * C * (- cosa * sinb) - Abs3 * D * (sina * sinb)

      call boostn( ql1(0), P_in(0), ql1(0) )
      call boostn( ql3(0), P_in(0), ql3(0) )

c Momenta of V2
      q2_out(0) = P_in(0) - ql1(0) - ql3(0)
      q2_out(1) = P_in(1) - ql1(1) - ql3(1)
      q2_out(2) = P_in(2) - ql1(2) - ql3(2)
      q2_out(3) = P_in(3) - ql1(3) - ql3(3)


c symmetrisation of q1 and q3
      msqlocal1 = m1sq
      msqlocal3 = m3sq
      symm = mod(symm,2) + 1
      if (symm.eq.1) then !interchange q1 <-> q3
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (symm.eq.2) then
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q3_out(mu) = ql1(mu)
         enddo
         m1sq = msqlocal3
         m3sq = msqlocal1
      endif


c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      ThreeBodyDecay_www_kk = .true.

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE DEF_PS_ZZW_kk(bos1,bos2,bos3)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/kkcoupl.inc"
#include "VBFNLO/utilities/cuts.inc"
      INTEGER bos(3), i, bos1, bos2, bos3

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTJMIN
      DOUBLE PRECISION RM2MINX, RM2MAXX, M2MIN, RM2MINY, RM2MAXY
      DOUBLE PRECISION RM2MAX1, RM2MAX2
	Double precision rd1min,rd3min,rd3max,cwidth
	common /mbrdborders/ rd1min,rd3min,rd3max,cwidth
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTJMIN, RM2MAX1, RM2MAX2


      DOUBLE PRECISION CLR, XM2, XMG, B, V, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), V(4,5), A(4,5)

      s = ecm**2

      bos(1) = bos1        ! W+ = 3
      bos(2) = bos2        ! W- = 4
      bos(3) = bos3        ! Z  = 2

      do i = 1,3
         rm2(i) = xm2(bos(i))        ! mass squared
         rmg(i) = xmg(bos(i))        ! mass * decay width
         if (bos(i).eq.2) then
            rm2min(i) = 1d-3
            rm2max(i) = s/2d0
         else
            rm2min(i) = 1d-3         ! minimal/maximal value of mass squared
            rm2max(i) = s/2d0        ! for tan mapping
         endif
      enddo

      rm2minx = 1d-3
      rm2maxx = s/2d0

      rm2miny = 1d-3
      rm2maxy = s/2d0

c      rm2max1 = (sqrt(xm2(6)-30d0*xmg(6)))**2 !Dieters proposal
c      rm2max2 = (sqrt(xm2(6)+30d0*xmg(6)))**2 !Dieters proposal

	if (sqrt(xm2(6)).lt.400d0)then
	rd1min=0.2d0
	else if (sqrt(xm2(6)).lt.1200d0) then 
	rd1min=0.2d0+0.000625d0*(sqrt(xm2(6))-400d0)
	else
	rd1min=0.7
	endif

	if (sqrt(xm2(6)).lt.400d0)then
	rd3min=0.2d0
	rd3max=0.5d0
	else if (sqrt(xm2(6)).lt.1000d0) then
	rd3min=0.2d0+0.0008333d0*(sqrt(xm2(6))-400d0)
	rd3max=0.5d0+0.000666d0*(sqrt(xm2(6))-400d0)
	else
	rd3min=0.7d0
	rd3max=0.9d0
	endif

	cwidth=1200d0

c coupling routine handles over 10 times width times mass
	if(mhlmswitch) then
      rm2max1 = ((xm2(6)-0.4d0*xmg(6))/sqrt(xm2(6)) )**2
      rm2max2 = ((xm2(6)+0.4d0*xmg(6))/sqrt(xm2(6)) )**2
	else
      rm2max1 = ((xm2(6)-0.3d0*xmg(6))/sqrt(xm2(6)) )**2
      rm2max2 = ((xm2(6)+0.4d0*xmg(6))/sqrt(xm2(6)) )**2
	endif


c cutoff for numerical purposes

      m2min = 0.01d0**2
      ptjmin = max(ptmin_def(1),10.0d0)

 60   format ( ' machine energy = ', f10.1, '  GeV' )
 70   format ( ' integration range for boson ',i3,' = ',f10.3,
     1         '  to ',f12.3,' GeV' )
 72   format ( ' boson has mass ',f8.2,' GeV and width ',
     1         f10.4,' GeV' )
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

      end


C**********************  LIPSZZW_kk   ****************************************
C 
      SUBROUTINE LIPSZZW_kk(N,RD,RN,
     1                   K1,K2,V,Pjet,
     2                   X1,X2,W,NW)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/cuts.inc"

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
	Double precision rd1min,rd3min,rd3max,cwidth
	common /mbrdborders/ rd1min,rd3min,rd3max,cwidth

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
      DOUBLE PRECISION xxwidth, xxm, xywidth, xym, rt
      DOUBLE PRECISION q1(0:3), q2(0:3), q3(0:3), K(0:4),pY(0:4)
      DOUBLE PRECISION tau, y, y_border

      DOUBLE PRECISION dotrr, vdummy(0:3,2)
	Double precision vdump, phi, cosphi, sinphi
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c LO parameter for WWZ inv mass
      PARAMETER ( xxwidth = 280d0, xxm = 470d0 )
      PARAMETER ( xywidth = 2200d0, xym = 1050d0 )
c NLO real emission parameter for WWZ invariant mass

      EXTERNAL dotrr

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

      LOGICAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_zzw_kk
      EXTERNAL RESONANCE, TwoBodyDecay, ThreeBodyDecay_zzw_kk
	logical nonresonance, nonresonance_kk
	external nonresonance, nonresonance_kk
c symmetrisation of q1 and q3
      INTEGER symm

      data symm /1/
      save symm    
c DEBUG	
	logical ampout
	double precision mom(0:4)
	common/MB3/ m12sq, m23sq
	common/MB/ ampout

	Double precision expo3v1, m3v1, g3v1
	Double precision expo3v2, m3v2, g3v2		
	PARAMETER (expo3v1=1.5d0, m3v1=200d0, g3v1=850d0)
	PARAMETER (expo3v2=4.26d0, m3v2=1.9d0, g3v2=55500d0)
c--------------------------------------------------------------------

      if (ldebug) print*,' enter LIPSWWW '
      W = 1d0             !init of weightfactor

c------------------------------------------------------------------------

      IF(N.EQ.1) THEN ! 2 -> JET+Y

c generate q^2_Y which then decays into three massive particles from rd(1)


c	if (rd(1).le. 0.6d0) then
c	      if (.not. Resonance(530d0**2, 530d0*340d0, 
c     1	 1d-6, (1.2*sqrt(xm2(6)))**2, rd(1)/0.6d0, W, pY(4))) return
c	      W=W/0.6d0
c	else
c      	if (.not. Resonance(1550d0**2, 1550d0*1500d0, 
c     1	 (1.2*sqrt(xm2(6)))**2, s, (rd(1)-0.6d0)/0.4d0, W, pY(4))) return
c		W=W/0.4d0
c	endif


c	if (rd(1).le. 0.05) then
c	      if (.not. nonResonance(1d-6, 130d0**2, rd(1)/0.05,
c     &					W, pY(4))) return
c	W=W/0.05
c	else
c	      if (.not. nonResonance_kk(100d0**2, 100d0*1400d0, 2d0, 
c     1	 130d0**2, s, (rd(1)-0.05)/0.95, W, pY(4))) return
c	W=W/0.95
c	endif
      	if (.not. Resonance(470d0**2, 470d0*900d0, 
     1	 1d-6, s, rd(1), W, pY(4))) return

c generate 2 -> 2 process for k1 k2 -> Pjet q_Y from rd(16), rd(17),rd(18)
         call TwoToJetPlusX(rd(16), 0d0, ecm, 20d0, pY(4),
     1           k1(0), k2(0), x1, x2, pY(0), pJet(0), W)

         if (w.eq.0) return

      ELSE ! 2 -> Y

c generate m_wwz^2 out of rd(1) via tan mapping and transform it to tau
	if (rd(1) .le. 0.03) then

         if (.not. nonResonance(1d-6,
     1   	250d0**2, rd(1)/0.03d0, W, K(4))) return
         W = W/0.03d0*2d0*pi

	else if (rd(1) .le. rd1min) then

	  if (.not. Resonance(430d0**2,430d0*500d0,250d0**2,(1.1d0*sqrt(xm2(6)))**2,
     &		(rd(1)-0.03)/(rd1min-0.03d0),W,K(4)))return
c	   if(.not. nonresonance_kk(157d0**2,157d0*840d0, 2d0,250d0**2,
c     & 		(1.15d0*sqrt(xm2(6)))**2, (rd(1)-0.03)/(rd1min-0.03d0),W,K(4)))return
	  	W=W/(rd1min-0.01d0)*2d0*pi

      else 
         if (.not. Resonance((1.5d0*sqrt(xm2(6)))**2, cwidth*(1.5d0*sqrt(xm2(6))), 
     1 (1.1d0*sqrt(xm2(6)))**2, s, (rd(1)-rd1min)/(1d0-rd1min), W, K(4))) return
         W = W/(1d0-rd1min)*2d0*pi
      endif	

c generate tau= K(4)/s
		tau=k(4)/s
		W = W/s
c generate y out of rd(2)

         y_border = log(1/tau)
         y = y_border*( RD(2) - 0.5d0 )
         W = W * y_border

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

c assign jet momenta
         pjet(0) = 0.0d0
         pjet(1) = 0.0d0
         pjet(2) = 0.0d0
         pjet(3) = 0.0d0

         do mu =0,4
            pY(mu) = K(mu)
         enddo

      ENDIF

		if (n .eq. 0) then
		rt=0d0
		else
		rt=rd(2)
		endif
      
c generate phase space for decay of P -> V + V + V out of rd(3):rd(9), rn
         if (.not. ThreeBodyDecay_zzw_kk(rd(3), rt, pY(0), q1(0), q2(0), q3(0), 
     1                                       m1sq, m2sq, m3sq,W,N)) return
      

c now let the vector bosons decay ************************************
c generate phase space for decay of V1 -> l1 + l2 out of rd(10), rd(11)
      call TwoBodyDecay0(rd(10),rd(11),q1(0),m1sq,v(0,1),v(0,2),W)

c generate phase space for decay of V2 -> l3 + l4 out of rd(12), rd(13)
      call TwoBodyDecay0(rd(12),rd(13),q2(0),m2sq,v(0,3),v(0,4),W)

c generate phase space for decay of V3 -> l5 + l6 out of rd(14), rd(15)
      call TwoBodyDecay0(rd(14),rd(15),q3(0),m3sq,v(0,5),v(0,6),W)

c Additional weight factors ******************************************

	phi=twopi*rn
	sinphi=sin(phi)
	cosphi=cos(phi)
	do i=1,6
		vdump=v(1,i)
		v(1,i)=cosphi*v(1,i)+sinphi*v(2,i)
		v(2,i)=cosphi*v(2,i)-sinphi*vdump
	enddo
	if (n .eq. 1) then
		vdump=pjet(1)
		pjet(1)=cosphi*pjet(1)+sinphi*pjet(2)
		pjet(2)=cosphi*pjet(2)-sinphi*vdump
	endif

      if(N.EQ.1) then
         W=W
      else
         W = W * 1.d0/2d0/K(4)       ! Flux
         W = W * 3.89379304d11        !conversion into fb
      endif

c exchange v3,v4 with v5,v6 for Matrixelement in ZZW kk scenario

	do i=0,3
	vdummy(i,1)=v(i,3)
	vdummy(i,2)=v(i,4)
	v(i,3)=v(i,5)
	v(i,4)=v(i,6)
	v(i,5)=vdummy(i,1)
	v(i,6)=vdummy(i,2)
	enddo

      return
      end

***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay_zzw_kk(R, RN, P_in, q1_out, 
     &                  q2_out, q3_out, m1sq, m2sq, m3sq, w, N)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(8), rn       : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecay         : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE
c to get all resonances, include:
#include "VBFNLO/utilities/kkcoupl.inc"
      REAL*8 RM2, RMG, RM2MIN, RM2MAX, S, PTC
      REAL*8 RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      REAL*8 RM2MAX1, RM2MAX2
	
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
	Double precision rd1min,rd3min,rd3max,cwidth
	common /mbrdborders/ rd1min,rd3min,rd3max,cwidth
      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU,N

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(7), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

      DOUBLE PRECISION xywidth, xym,rm2mid2,rm2mid1
      PARAMETER ( xywidth = 400d1, xym = 700d0)

      LOGICAL RESONANCE, nonresonance_kk, nonresonance
      EXTERNAL RESONANCE, nonresonance_kk, nonresonance

c variables for ZW invariant mass:
	double precision expo2v1, m2v1, g2v1,g2v2
	parameter( 	expo2v1	=	2d0		)
	parameter( 	m2v1		=	60d0		)

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION QL1(0:3), QL3(0:3), msqlocal1, msqlocal3
      data symm /1/
      save symm    

	if (n.eq.0) then
	g2v1=690d0
	else
	g2v1=600d0
	if (w1mass_k(1).le.500) then
        rd3min=0.4d0
        rd3max=0.7d0
	else
        rd3min=0.65d0
        rd3max=0.8d0
	endif
	endif

      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first

c m_12^2 = (q1 + q2)^2 out of R(1)

      if ( r(1) .le. 0.05d0 ) then
         if (.not. nonResonance(rm2min12, 66d0**2d0,
     1                              r(1)/0.05d0, W, m12sq)) return
         W = W/0.05d0

	else if ( r(1) .le.  rd3min ) then
	   if (.not. nonResonance_kk( m2v1**2, m2v1*g2v1, expo2v1, 
     1 	66d0**2, rm2max1, (r(1)-0.05d0)/(rd3min-0.05d0), W, m12sq) ) return
	   W = W/(rd3min-0.05d0)

      else if ( r(1) .le.  rd3max ) then
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*w1width_k(1),
     1	      rm2max1, rm2max2,(r(1)-rd3min)/(rd3max-rd3min), W, m12sq)) return
         W = W/(rd3max-rd3min)

      else
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*cwidth, rm2max2, rm2max12,
     1                      (r(1)-rd3max)/(1d0-rd3max), W, m12sq)) return
         W = W/(1d0-rd3max)
      endif


      if (m12sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c m_23^2 = (q2 + q3)^2 out of r(2)

      if ( r(2) .le. 0.05d0 ) then
         if (.not. nonResonance(rm2min12, 66d0**2d0,
     1                              r(2)/0.05d0, W, m23sq)) return
         W = W/0.05d0

	else if ( r(2) .le.  rd3min ) then
	   if (.not. nonResonance_kk( m2v1**2, m2v1*g2v1, expo2v1, 
     1 	66d0**2, rm2max1, (r(2)-0.05d0)/(rd3min-0.05), W, m23sq) ) return
	   W = W/(rd3min-0.05)

      else if ( r(2) .le.  rd3max ) then
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*w1width_k(1),
     1 	rm2max1, rm2max2,(r(2)-rd3min)/(rd3max-rd3min), W, m23sq)) return
         W = W/(rd3max-rd3min)

      else
         if (.not. Resonance(w1mass_k(1)**2, w1mass_k(1)*cwidth, rm2max2, rm2max23,
     1                      (r(2)-rd3max)/(1d0-rd3max), W, m23sq)) return
         W = W/(1d0-rd3max)
      endif

      if (m23sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c now generate m1^2 and m3^2 of the massive vector bosons **********

c m_1^2 is generated via tan mapping out of r(3)
      rm2max(1) = min(p_in(4)-m23sq,p_in(4)+m23sq-2.0d0*sqrt(p_in(4)*m23sq))
	rm2mid2=130d0
	rm2mid1=70d0

      if (rm2max(1).le.0) then
         W = 0.0d0
         return
      elseif(rm2max(1).gt.rm2mid2**2)then
		if (r(3).lt.0.1)then
 	        if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(1), rm2mid1**2,
     1                                 r(3)/0.1d0, W, m1sq)) return
		W=W/0.1d0
		elseif (r(3).le.0.9d0)then
 	        if (.not. Resonance(rm2(1), rmg(1), rm2mid1**2, rm2mid2**2,
     1                                 (r(3)-0.1d0)/0.8d0, W, m1sq)) return
		W=W/0.8d0
		else
 	        if (.not. Resonance(rm2(1), 91d0*1800d0, rm2mid2**2, rm2max(1),
     1                                 (r(3)-0.9d0)/0.1d0, W, m1sq)) return
		W=W/0.1d0
		endif
	elseif(rm2max(1).gt.rm2mid1**2)then
		rm2mid2=sqrt(rm2max(1))
		if (r(3).lt.0.1)then
	         if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(1), rm2mid1**2,
     1                                 r(3)/0.1d0, W, m1sq)) return
		W=W/0.1d0
		else
	         if (.not. Resonance(rm2(1), rmg(1), rm2mid1**2, rm2mid2**2,
     1                                 (r(3)-0.1d0)/0.9d0, W, m1sq)) return
		W=W/0.9d0
		endif
	else
		rm2mid1=sqrt(rm2max(1))
	         if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(1), rm2mid1**2,
     1                                 r(3)/0.1d0, W, m1sq)) return
		W=W/0.1d0
	endif
	
c m_3^2 is generated via tan mapping out of r(4)
      rm2max(3) = min(p_in(4)-m12sq, p_in(4)+m12sq-2.0d0*sqrt(p_in(4)*m12sq), 
     1                m12sq+m23sq-m1sq)
	rm2mid2=130d0
	rm2mid1=70d0

      if (rm2max(3).le.0) then
         W = 0.0d0
         return
	elseif(rm2max(3).gt.rm2mid2**2)then
		if (r(4).lt. 0.1d0) then
	         if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(3), rm2mid1**2,
     1                                 r(4)/0.1d0, W, m3sq)) return
		W=W/0.1d0
		elseif (r(4).le.0.9d0)then
	         if (.not. Resonance(rm2(3), rmg(3), rm2mid1**2, rm2mid2**2,
     1                                 (r(4)-0.1d0)/0.8d0, W, m3sq)) return
		W=W/0.8d0
		else
	         if (.not. Resonance(rm2(3), 91d0*1800d0, rm2mid2**2, rm2max(3),
     1                                 (r(4)-0.9d0)/0.1d0, W, m3sq)) return
		W=W/0.1d0
		endif
	elseif(rm2max(3).gt.rm2mid1**2)then
		rm2mid2=sqrt(rm2max(3))
		if (r(4).lt. 0.1d0) then
	         if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(3), rm2mid1**2,
     1                                 r(4)/0.1d0, W, m3sq)) return
		W=W/0.1d0
		else
	         if (.not. Resonance(rm2(3), rmg(3), rm2mid1**2, rm2mid2**2,
     1                                 (r(4)-0.1d0)/0.9d0, W, m3sq)) return
		W=W/0.9d0
		endif
	else
		rm2mid1=sqrt(rm2max(3))	
	         if (.not. Resonance(91d0**2, 91d0*100d0, rm2min(3), rm2mid1**2,
     1                                 r(4), W, m3sq)) return
	endif

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2/M
      E3 = (p_in(4) + m3sq - m12sq)/2/M
      E2 = M - E1 - E3

      Abs1 = sqrt(E1**2 - m1sq)
      Abs3 = sqrt(E3**2 - m3sq)

c m_2^2 is generated via tan mapping out of r(5)
      rm2min(2) = max(-Abs1**2-Abs3**2-2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1            1.0d-5)
      rm2max(2) = min(-Abs1**2-Abs3**2+2.d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1          E2**2, P_in(4))

      if (rm2max(2).le.0) then
         W = 0.0d0
         return
      endif
      
         if (.not. Resonance(rm2(2), rmg(2), rm2min(2), rm2max(2),
     1                                 r(5), W, m2sq)) return
      
c now generate Euler angles
      alpha = R(6) * twopi
      W = W*twopi

      cosb = -1.0d0 + R(7)*2.0d0
      W = W*2.0d0

      gamma = RN * twopi
      W = W*twopi

      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      sinb = sqrt(1 - cosb**2)
      sing = sin(gamma)
c calculate now C and D
      D = -(m1sq-m2sq+m3sq+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)/2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
c         print*,'D**2 >1'
         return
      endif

c Momenta of V1
      ql1(0) = E1
      ql1(1) = Abs1 * (cosa * cosg - cosb * sina * sing)
      ql1(2) = Abs1 * (-cosb * cosg * sina - cosa * sing)
      ql1(3) = Abs1 * (sina * sinb)

c Momenta of V3
      ql3(0) = E3
      ql3(1) = Abs3 * C * (cosg * sina + cosa * cosb * sing) 
     1      - Abs3 * D * (cosa * cosg - cosb * sina * sing)
      ql3(2) = Abs3 * C * (cosa * cosb * cosg - sina * sing)
     1      - Abs3 * D * (-cosb * cosg * sina - cosa * sing)
      ql3(3) = Abs3 * C * (- cosa * sinb) - Abs3 * D * (sina * sinb)

      call boostn( ql1(0), P_in(0), ql1(0) )
      call boostn( ql3(0), P_in(0), ql3(0) )

c Momenta of V2
      q2_out(0) = P_in(0) - ql1(0) - ql3(0)
      q2_out(1) = P_in(1) - ql1(1) - ql3(1)
      q2_out(2) = P_in(2) - ql1(2) - ql3(2)
      q2_out(3) = P_in(3) - ql1(3) - ql3(3)


c symmetrisation of q1 and q3
      msqlocal1 = m1sq
      msqlocal3 = m3sq
      symm = mod(symm,2) + 1
      if (symm.eq.1) then !interchange q1 <-> q3
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (symm.eq.2) then
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q3_out(mu) = ql1(mu)
         enddo
         m1sq = msqlocal3
         m3sq = msqlocal1
      endif


c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

c	if ((e1*e2*e3).lt.0d0)then
c	print*, "ENERGY LOWER THAN 0d0!!!", E1, E2, E3
c	print*, r, P_in
c	stop
c	endif
      ThreeBodyDecay_zzw_kk = .true.

      RETURN
      END

