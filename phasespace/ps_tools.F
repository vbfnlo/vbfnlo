****************************************************************************
*   This file contains useful functions and subroutines for building a 
*   phasespace generator.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c  FUNCTION Resonance(xm2, xmg, rm2min, rm2max, R, weight, qsq)
c  FUNCTION NonResonance(rm2min, rm2max, R, weight, qsq)
c  FUNCTION TwoBodyDecay(R1,R2,P_in,M2_in,M2_out1,M2_out2,P_out1,P_out2,weight)
c  FUNCTION TwoToJetsPlusX(N,R,RN,ECM,PTmin,XQ2,K1,K2,X1,X2,Q,P,weight)
c  SUBROUTINE TwoBodyDecay0(R1, R2, P_in, M2_in, P_out1, P_out2, weight)
c  SUBROUTINE ThreeBodyDecay0(R, P_in, M2_in, P_1, P_2, P_3, weight)
c  SUBROUTINE TwoToTwo(R,RN,ECM,Qsqr1,Qsqr2,K1,K2,X1,X2,Q1,Q2,weight)
c  SUBROUTINE TwoToOne(R,ECM,Qsqr,K1,K2,X1,X2,Q,weight)
c
c***************************************************************************


****************************************************************************
      LOGICAL FUNCTION Resonance(xm2, xmg, rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       xmg       : mass*width of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, xmg, rm2min, rm2max, R, weight, qsq, tanx
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      DOUBLE PRECISION xmin, xmax, deltax, x 

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN
      ENDIF
      
      xmin = ATAN( (rm2min - xm2)/xmg )
      xmax = ATAN( (rm2max - xm2)/xmg )
      deltax = xmax - xmin
      x = xmin + deltax * R
      tanx = tan(x)
      qsq = xm2 + tanx * xmg

      
      IF (qsq.LE.0 .or. deltax.le.0) THEN
         weight = 0d0 
         Resonance = .false.
         RETURN 
      ENDIF

      weight = weight * deltax / twopi
      weight = weight * xmg * (1 + tanx**2 )

      Resonance = .true.

      RETURN 
      END


****************************************************************************
      LOGICAL FUNCTION Resonance2(xm2, xmg, rm2min, rm2max, nowidth, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance with a high-energy 1/q^2 tail
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       xmg       : mass*width of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       nowidth   : transition from BW to 1/q^2 happens at M + crossing * Gamma
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, xmg, rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      DOUBLE PRECISION nogammamin, nogammamax, rm2split, split
      DOUBLE PRECISION norm, nowidth
      LOGICAL RESONANCE
      EXTERNAL RESONANCE

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         Resonance2 = .false.
         RETURN
      ENDIF
      
      ! nowidth = 20d0

c number of width above resonance
      nogammamin = (sqrt(rm2min)-sqrt(xm2))/(xmg/sqrt(xm2))
      nogammamax = (sqrt(rm2max)-sqrt(xm2))/(xmg/sqrt(xm2))
c adjust split linearly if between resonance and nowidth width
      if (nogammamax .gt. nowidth) then
        split = (nowidth-nogammamin)/nowidth*0.8d0
        if (split .le. 0.1d0 .and. split.gt.0d0) split = 0.1d0
        if (split .le. 0d0) split = 0d0
        if (split .ge. 0.8d0) split = 0.8d0
        rm2split = max((sqrt(xm2)+nowidth*xmg/sqrt(xm2))**2,rm2min)
      else
        split = 1d0
        rm2split = rm2max
      endif

      ! print*, nogammamin, nogammamax, nowidth, split
      if (R .lt. split) then
c normal Breit-Wigner resonance
        Resonance2 = Resonance(xm2, xmg, rm2min, rm2split, R/split, weight, qsq)
        weight = weight/split
      else
c 1/q^2 behaviour
        norm=log((rm2max-xm2)/(rm2split-xm2))
        qsq = xm2 + (rm2split-xm2)*exp((R-split)/(1d0-split)*norm)
        weight = weight/(1d0-split)*(qsq-xm2)*norm/twopi
      endif

      Resonance2 = .true.

      RETURN 
      END


****************************************************************************
      LOGICAL FUNCTION ResonanceQsq(xm2, rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance using a 1/q^2 shape
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      DOUBLE PRECISION norm

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         ResonanceQsq = .false.
         RETURN
      ENDIF

      IF ((xm2.ge.rm2min).and.(xm2.le.rm2max)) THEN
         weight = 0d0 
         ResonanceQsq = .false.
         print *,"Error in ResonanceQsq: trying to integrate over divergence"
         print *," at mass peak: ", sqrt(xm2), " in [", sqrt(rm2min), ",", sqrt(rm2max), "]"
      ENDIF
      
      norm=log((rm2max-xm2)/(rm2min-xm2))
      qsq = xm2 + (rm2min-xm2)*exp(R*norm)
      weight = weight*(qsq-xm2)*norm/twopi

      ResonanceQsq = .true.

      RETURN 
      END


****************************************************************************
      LOGICAL FUNCTION ResonanceVbos(xm2, xmg, rm2min, rm2max, nowidth, bos, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a resonance using a BW for the main part,
*     a 1/q^2 shape at high energies and a 1/q^2 shape (constant) at low energies
*     for a Z (W) boson
*
*    INPUT
*       xm2       : mass**2 of the resonance
*       rm2min    : lowerBound**2 of the resonance
*       rm2max    : upperBound**2 of the resonance
*       nowidth   : transition from BW happens at M +- nowidth * Gamma
*       R         : a random number between 0 and 1
*       weight    : phasespace weight
*
*    OUTPUT
*       weight    : new phasespace weight (multiplied by Jacobi factor)
*       qsq       : the Q**2 
*       Resonance : .true.,  if everything is ok 
*                   .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION xm2, xmg, rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION nowidth
      integer bos
      logical nonresonance, resonanceqsq, resonance2
      external nonresonance, resonanceqsq, resonance2

      ResonanceVbos = .false.

      if (R .lt. 0.3d0) then
         weight = weight / 0.3d0
         if (bos.eq.2) then
           if (.not.ResonanceQsq(0d0,rm2min, min(XM2-nowidth*XMG,rm2max), R/0.3d0, weight, qsq)) RETURN 
         else
           if (.not.NonResonance(RM2MIN, min(XM2-nowidth*XMG,rm2max), R/0.3d0, weight, qsq)) RETURN 
         endif
      else
         weight = weight / 0.7d0
         if (.not.Resonance2(XM2, XMG, XM2-nowidth*XMG, rm2max, nowidth, (R-0.3d0)/0.7d0, weight, qsq)) RETURN
      endif

      ResonanceVbos = .true.

      RETURN 
      END


****************************************************************************
      LOGICAL FUNCTION NonResonance(rm2min, rm2max, R, weight, qsq)
****************************************************************************
*     Calculates the Q-square for a non-resonance (a flat distribution)
*
*    INPUT
*       rm2min       : lowerBound**2 of the resonance
*       rm2max       : upperBound**2 of the resonance
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2 
*       NonResonance : .true.,  if everything is ok 
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION rm2min, rm2max, R, weight, qsq
      DOUBLE PRECISION pi
      PARAMETER (pi = 3.141592653589793d0)
      DOUBLE PRECISION minMass, maxMass, deltaMass, mass

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0 
         NonResonance = .false.
         RETURN 
      ENDIF
      
      minMass = sqrt(rm2min)
      maxMass = sqrt(rm2max)
      deltaMass = maxMass - minMass
      mass = minMass + deltaMass * R
      qsq = mass**2

      weight = weight * deltaMass * mass / pi

      NonResonance = .true.
      RETURN 
      END

****************************************************************************
      LOGICAL FUNCTION NonResonance_kk(xm2, xmg, expo, rm2min, rm2max,
     &                                                R, weight, qsq)
****************************************************************************
*    Calculates q^2 for atan((q^2-xm2)/xmg)**expo / ((q^2-xm2)^2+xmg^2) distribution
*    !!!!!!!!!!!   xm2 needs to be lower than rm2min  !!!!!!!!!!!!!!!!!!!!!
*    INPUT
*       rm2min       : lowerBound**2 of the distribution
*       rm2max       : upperBound**2 of the distribution
*       R            : a random number between 0 and 1
*       weight       : phasespace weight
*       xm2, xmg, expo: parameters of distribution
*    OUTPUT
*       weight       : new phasespace weight (multiplied by Jacobi factor)
*       qsq          : the Q**2
*       NonResonance : .true.,  if everything is ok     
*                      .false., if kinematics not ok -> weight = 0
****************************************************************************

      IMPLICIT NONE

      double precision rm2min, rm2max, R, weight, qsq, xmg, xm2, expo
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      double precision xmax, x, tanx, deltax, xmin, n

      IF (rm2max.LT.rm2min) THEN
         weight = 0d0
         NonResonance_kk = .false.
         RETURN
      ENDIF

      IF (xm2 .ge. rm2min) THEN
         weight = 0d0
         NonResonance_kk = .false.
         RETURN
      ENDIF

      n=1d0/(expo+1d0)

      xmin = ATAN( (rm2min-xm2)/xmg )**(1/n)
      xmax = ATAN( (rm2max-xm2)/xmg )**(1/n)
      deltax= xmax -xmin
      x =xmin + deltax * R

      tanx = tan(x**n)

      qsq = xm2 + tanx * xmg

      weight = weight * deltax / twopi
      weight = weight * xmg * (1 + tanx**2 )
      weight = weight * n * x**(n-1)

      NonResonance_kk = .true.

      RETURN
      END

****************************************************************************
      LOGICAL FUNCTION TwoBodyDecay(R1, R2, P_in, M2_in, M2_out1, M2_out2, 
     &                              P_out1, P_out2, weight)
****************************************************************************
*     Calculates the two-body decay P_in -> P_out1 + P_out2
*
*    INPUT
*       R1, R2         : two random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       M_out1, M_out2 : mass**2 of daughter particles
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_out1, P_out2 : 4-momenta of daughter particles
*       TwoBodyDecay   : .true.,  if everything is ok 
*                        .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION R1, R2, P_in(0:3), M2_in, M2_out1, M2_out2
      DOUBLE PRECISION P_out1(0:3), P_out2(0:3), weight
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)
      INTEGER*4 mu
      double precision Q(0:3), CTH, STH, PHID, M_in, M_out1, M_out2, E1, P1

      M_in = sqrt(M2_in)
      M_out1 = sqrt(M2_out1)
      M_out2 = sqrt(M2_out2)

      IF (M_in .LT. (M_out1 + M_out2) ) THEN
         weight = 0d0 
         TwoBodyDecay = .false.
         RETURN 
      ENDIF

c generate outgoing momenta in CMS and boost back
      CTH = -1d0 + 2d0 * R1
      if ((1D0 - CTH**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      STH = sqrt(1d0 - CTH**2 )
      PHID = twopi * R2 - twopi/2.0d0
      
      E1 = (M2_in - M2_out2 + M2_out1) / (2d0*M_in)
      if ((E1**2 - M2_out1) .le. 1D-10) then
         weight = 0d0
         TwoBodyDecay = .false.
         return
      end if
      P1 = sqrt(E1**2 - M2_out1)
      Q(0) = E1
      Q(1) = P1 * STH*COS(PHID)
      Q(2) = P1 * STH*SIN(PHID)
      Q(3) = P1 * CTH

      call boostn(Q,P_in,P_out1)

      do mu = 0,3
         P_out2(mu) = P_in(mu) - P_out1(mu)
      enddo

      weight = weight * P1/(M_in * 2d0 * twopi) 
      
      TwoBodyDecay = .true.
      RETURN
      END

****************************************************************************
      SUBROUTINE TwoBodyDecay0(R1, R2, P_in, M2_in, P_out1, P_out2, weight)
****************************************************************************
*     Calculates the two-body decay P_in -> P_out1 + P_out2 
*     for MASSLESS decay particles 
*
*    INPUT
*       R1, R2         : two random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_out1, P_out2 : 4-momenta of daughter particles
****************************************************************************
      IMPLICIT NONE
      double precision R1, R2, P_in(0:3), M2_in
      double precision P_out1(0:3), P_out2(0:3), weight
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      integer*4 mu
      double precision Q(0:3), CTH, STH, PHID, M_in, E1

      M_in = sqrt(M2_in)

c generate outgoing momenta in CMS and boost back
      CTH = -1d0 + 2d0 * R1
      if ((1D0 - CTH**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      STH = sqrt(1d0 - CTH**2 )
      PHID = twopi * R2 - twopi/2d0
      
      E1 = 0.5d0 * M_in
      Q(0) = E1
      Q(1) = E1 * STH*COS(PHID)
      Q(2) = E1 * STH*SIN(PHID)
      Q(3) = E1 * CTH

      call boostn(Q,P_in,P_out1)
      do mu = 0,3
         P_out2(mu) = P_in(mu) - P_out1(mu)
      enddo

      weight = weight / (4d0 * twopi) 
      
      RETURN
      END

****************************************************************************
      SUBROUTINE ThreeBodyDecay0(R, P_in, M2_in, P_1, P_2, P_3, weight)
****************************************************************************
*     Calculates the three-body decay P_in -> P_1 + P_2 + p_3
*     for MASSLESS decay particles via Euler angles
*
*    INPUT
*       R              : array of 5 random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       M_in           : mass**2 of mother particle
*       weight         : phasespace weight
*
*    OUTPUT
*       weight         : new phasespace weight (multiplied by Jacobi factor)
*       P_1, P_2, P_3  : 4-momenta of daughter particles
****************************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION R(*)
      double precision P_in(0:3), M2_in
      double precision P_1(0:3), P_2(0:3), P_3(0:3), weight
      double precision twopi
      parameter (twopi = 2d0*3.141592653589793d0)
      integer*4 mu
      double precision alpha, sin_alpha, cos_alpha
      double precision gamma, sin_gamma, cos_gamma
      double precision cos_beta, sin_beta, s1, c1
      double precision E1, E2, M
      
      alpha = twopi*R(1)
      cos_beta = 1-2*R(2)
      gamma = twopi*R(3)
      
      if((1-cos_beta**2) .le. 1D-10) then
	weight = 0D0
	return 
      else	
	sin_beta = sqrt(1-cos_beta**2)
      end if
      sin_alpha = sin(alpha)
      cos_alpha = cos(alpha)
      sin_gamma = sin(gamma)
      cos_gamma = cos(gamma)
      
      M = sqrt(M2_in)
      E1 = 0.5d0*R(4)*M
      E2 = 0.5d0*M*(1-R(4)+R(5)*R(4))
      
      P_1(0) = E1
      P_1(1) = E1*(  cos_alpha*cos_gamma - cos_beta*sin_alpha*sin_gamma )
      P_1(2) = E1*( -cos_alpha*sin_gamma - cos_beta*cos_gamma*sin_alpha )
      P_1(3) = E1*sin_alpha*sin_beta

      c1 = -(2d0*E1*E2 - 2d0*(E1 + E2)*M + M2_in)/(2d0*E1*E2)
      if ((1D0-c1**2) .le. 1D-10) then
         weight = 0D0
         return
      end if
      s1 = sqrt(1-c1**2)

      P_2(0) = E2
      P_2(1) = E2*( s1*(cos_gamma*sin_alpha + cos_alpha*cos_beta*sin_gamma) - 
     &              c1*(cos_gamma*cos_alpha - sin_alpha*cos_beta*sin_gamma) )
      P_2(2) = E2*( s1*( cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma) - 
     &              c1*(-sin_alpha*cos_beta*cos_gamma - cos_alpha*sin_gamma) )
      P_2(3) = E2*sin_beta*(-s1*cos_alpha - c1*sin_alpha)


      CALL BOOSTN(P_1, P_in, P_1)  
      CALL BOOSTN(P_2, P_in, P_2)

      do mu = 0,3
         P_3(mu) = P_in(mu)-P_1(mu)-P_2(mu)
      enddo

      weight = weight*(1d0/twopi**3)*(M2_in*R(4))/(16d0)
      
      RETURN
      END


***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay(R, RN, P_in, m1sq, m2sq, m3sq, 
     &                  ql1, ql2, ql3, w)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(4), RN         : five random numbers between 0 and 1
*       P_in             : 4-momentum of mother particle
*       m1sq, m2sq, m3sq : mass squared of outgoing particles
*       w                : phasespace weight
*
*    OUTPUT
*       w                 : new phasespace weight
*       ql1, ql2, ql3     : 4-momenta of daughter particles
*       ThreeBodyDecay    : .true.,  if everything is ok 
*                           .false., if kinematics not ok -> weight = 0
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)


      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), QL1(0:3), QL2(0:3), QL3(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(4), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c data to be useed in histogram routine
      double precision mZ1sq, mZ2sq, mZ3sq
      Common /histdata1/ m12sq, m23sq, mZ1sq, mZ2sq, mZ3sq, alpha, cosb, gamma

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )
      DOUBLE PRECISION x1width, x1m
      PARAMETER ( x1width = 70.0d0, x1m = 90.0d0 )
      DOUBLE PRECISION x2width, x2m
      PARAMETER ( x2width = 70.0d0, x2m = 90.0d0 )

      LOGICAL RESONANCE, NonResonance_KK, NonResonance
      EXTERNAL RESONANCE, NonResonance_KK, NonResonance

      ThreeBodyDecay = .false.

      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first **********

c generation of m12sq ---------------------------------
      if (.not. NonResonance(1.0d-5, p_in(4), r(1), W, m12sq)) return

c generation of m23sq next ---------------------------------
      if (.not. NonResonance(1.0d-5, p_in(4), r(2), W, m23sq)) return

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2.0d0/M
      E3 = (p_in(4) + m3sq - m12sq)/2.0d0/M
      E2 = M - E1 - E3
      Abs1 = sqrt(E1**2 - m1sq)
      Abs3 = sqrt(E3**2 - m3sq)

c now generate Euler angles (rd(3), rd(4), rn)
      alpha = R(3) * twopi
      W = W*twopi
      cosb = -1.0d0 + R(4)*2.0d0
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
      ql2(0) = P_in(0) - ql1(0) - ql3(0)
      ql2(1) = P_in(1) - ql1(1) - ql3(1)
      ql2(2) = P_in(2) - ql1(2) - ql3(2)
      ql2(3) = P_in(3) - ql1(3) - ql3(3)

c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      ThreeBodyDecay = .true.

      RETURN
      END


***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecayWWW(R, RN, P_in, q1_out, q2_out, q3_out,
     &                  m1sq, m2sq, m3sq, w, NW)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(7), rn       : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*       NW             : Narrow width (NW = 1) or full calculation (NW = 0)
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecayWWW      : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2
      DOUBLE PRECISION RMLIMIT
      COMMON /TRIBOSWWWBLIPSI/ RMLIMIT(2)

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU, NW

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq, mw, m2
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(7), RN
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION QL1(0:3), QL3(0:3), msqlocal1, msqlocal3
      data symm /1/
      save symm    

      LOGICAL RESONANCE, RESONANCE2, multiresonance
      EXTERNAL RESONANCE, RESONANCE2, multiresonance

      M = sqrt(p_in(4))

      ThreeBodyDecayWWW = .false.


c generate two invariant masses (qi + qj)^2  = mij^2 first **********

c m_12^2 = (q1 + q2)^2 out of R(1)
      if (.not. MultiResonance(1, r(1), rm2min12, rm2max12, W, m12sq)) return

      if (m12sq.ge.p_in(4)) then
         W = 0.0d0
         return
      endif

c m_23^2 = (q2 + q3)^2 out of r(2)
      if (.not. MultiResonance(1, r(2), rm2min23, rm2max23, W, m23sq)) return

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

      if (NW.eq.1) then
            if (rm2max(1).lt.rm2(1)) then
               W = 0.0d0
               return
            endif
            m1sq = rm2(1)
            W = W * rmg(1)/2.0d0
      else !NW.ne.1
         if(rm2max(1).ge.rm2(1)-16d0*rmg(1)) then
            if (r(3).lt.0.1d0) then
               W=W/0.1d0
               m2=max(sqrt(xm2(6))-sqrt(rm2(1)),0.0d0)
               mw=rmg(1)*10.0d0*m2
               m2 = m2**2
               if (.not.Resonance(m2, mw, rm2min(1), rm2(1)-16d0*rmg(1),
     1                              r(3)/0.1d0, W, m1sq)) RETURN 
            else
               W=W/0.9d0
               if (.not.Resonance2(rm2(1), rmg(1), rm2(1)-16d0*rmg(1),
     1                 rm2max(1),8d0, (r(3)-0.1d0)/0.9d0, W, m1sq)) RETURN 
            endif
         else
            m2=max(sqrt(xm2(6))-sqrt(rm2(1)),0.0d0)
            mw=rmg(1)*10.0d0*m2
            m2 = m2**2
            if (.not. Resonance(m2, mw, rm2min(1), rm2max(1),
     1                              r(3), W, m1sq)) return
         endif
      endif !NW

c m_3^2 is generated via tan mapping out of r(4)
      rm2max(3) = min(p_in(4)-m12sq, p_in(4)+m12sq-2.0d0*sqrt(p_in(4)*m12sq), 
     1                m12sq+m23sq-m1sq)
      if (rm2max(3).le.0) then
         W = 0.0d0
         return
      endif

      if (NW.eq.1) then
            if (rm2max(3).lt.rm2(3)) then
               W = 0.0d0
               return
            endif
            m3sq = rm2(3)
            W = W * rmg(3)/2.0d0
      else !NW.ne.1
         if(rm2max(3).ge.rm2(3)-16d0*rmg(3)) then
            if (r(4).lt.0.1d0) then
               W=W/0.1d0
               m2=max(sqrt(xm2(6))-sqrt(rm2(3)),0.0d0)
               mw=rmg(3)*10.0d0*m2
               m2 = m2**2
               if (.not.Resonance(m2, mw, rm2min(3), rm2(3)-16d0*rmg(3),
     1                              r(4)/0.1d0, W, m3sq)) RETURN 
            else
               W=W/0.9d0
               if (.not.Resonance2(rm2(3), rmg(3), rm2(3)-16d0*rmg(3),
     1                 rm2max(3), 8d0, (r(4)-0.1d0)/0.9d0, W, m3sq)) RETURN 
            endif
         else
            m2=max(sqrt(xm2(6))-sqrt(rm2(3)),0.0d0)
            mw=rmg(3)*10.0d0*m2
            m2 = m2**2
            if (.not. Resonance(m2, mw, rm2min(3), rm2max(3),
     1                              r(4), W, m3sq)) return
         endif
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

      if (NW.eq.1) then
         if (rm2max(2).lt.rm2(2) .or. rm2min(2).gt.rm2(2)) then
            W = 0.0d0
            return
         endif
         m2sq = rm2(2)
         W = W * rmg(2)/2.0d0
      else !NW.ne.1
         if ((rm2max(2).gt.rm2(2)-16d0*rmg(2)) .and. 
     1       (rm2min(2).lt.rm2(2)-16d0*rmg(2)) ) then
            if (r(5).lt.0.1d0) then
               W=W/0.1d0
               m2=max(sqrt(xm2(6))-sqrt(rm2(2)),0.0d0)
               mw=rmg(2)*10.0d0*m2
               m2 = m2**2
               if (.not.Resonance(m2, mw, rm2min(2), rm2(2)-16d0*rmg(2),
     1                                 r(5)/0.1d0, W, m2sq)) RETURN 
            else
               W=W/0.9d0
               if (.not.Resonance2(rm2(2), rmg(2), rm2(2)-16d0*rmg(2),
     1              rm2max(2),8d0, (r(5)-0.1d0)/0.9d0, W, m2sq)) RETURN 
            endif
         elseif (rm2min(2).lt.rm2(2)-16d0*rmg(2)) then   ! and rm2max is < rm2-16rmg
            m2=max(sqrt(xm2(6))-sqrt(rm2(2)),0.0d0)
            mw=rmg(2)*10.0d0*m2
            m2 = m2**2
            if (.not. Resonance(m2, mw, rm2min(2), rm2max(2),
     1                                 r(5), W, m2sq)) return
         else ! rm2min is > rm2-16rmg
            if (.not. Resonance2(rm2(2), rmg(2), rm2min(2), rm2max(2),8d0,
     1                                 r(5), W, m2sq)) return
         endif
      endif!NW

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

      ThreeBodyDecayWWW = .true.

      RETURN
      END


***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecay3Higgs(R, P_in, q1_out, q2_out,
     &                 q3_out, m1sq, m2sq, m3sq, w, NW)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out with
*     a Higgs resonance in all combinations (12, 13, 23)
*
*    INPUT
*       R(8)           : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*       NW             : Narrow width (NW = 1) or full calculation (NW = 0)
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecay3Higgs   : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************
      IMPLICIT NONE

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MINX, RM2MAXX, RM2MINY, RM2MAXY
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MINX, RM2MAXX, RM2MINY, RM2MAXY, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU, NW

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(8)
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c ordering of mi^2
      DOUBLE PRECISION mlight, mmiddle, mheavy
      INTEGER ident
      DOUBLE PRECISION QL1(0:3), QL2(0:3), QL3(0:3)

c data to be useed in histogram routine
      double precision mZ1sq, mZ2sq, mZ3sq
      Common /histdata1/ m12sq, m23sq, mZ1sq, mZ2sq, mZ3sq, alpha, cosb, gamma

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )

      LOGICAL RESONANCE, NONRESONANCE
      EXTERNAL RESONANCE, NONRESONANCE

      M = sqrt(p_in(4))

c generate m1^2, m2_2 and m3^2 of the massive vector bosons **********

c m_1^2 is generated via tan mapping out of r(1)
      rm2min(1) = 1.0d-5
      rm2max(1) = p_in(4)
      if (NW.eq.1) then
         if (rm2max(1).lt.rm2(1)) then
            W = 0.0d0
            return
         endif
         m1sq = rm2(1)
         W = W * rmg(1)/2.0d0
      else !NW.ne.1
         if (r(1).lt.0.3d0) then
            W=W/0.3d0
            if (.not.NonResonance(rm2min(1), rm2(1)-8d0*rmg(1),
     1                              r(1)/0.3d0, W, m1sq)) RETURN 
         else
            W=W/0.7d0
            if (.not.Resonance(rm2(1), rmg(1), rm2(1)-8d0*rmg(1),
     1                 rm2max(1), (r(1)-0.3d0)/0.7d0, W, m1sq)) RETURN 
         endif
      endif !NW

c m_2^2 is generated via tan mapping out of r(2)
      rm2min(2) = 1.0d-5
      rm2max(2) = (sqrt(p_in(4)) - sqrt(m1sq))**2
      if (NW.eq.1) then
         if (rm2max(2).lt.rm2(2)) then
            W = 0.0d0
            return
         endif
         m2sq = rm2(2)
         W = W * rmg(2)/2.0d0
      else !NW.ne.1
         if (r(2).lt.0.3d0) then
            W=W/0.3d0
            if (.not.NonResonance(rm2min(2), rm2(2)-8d0*rmg(2),
     1                              r(2)/0.3d0, W, m2sq)) RETURN 
         else
            W=W/0.7d0
            if (.not.Resonance(rm2(2), rmg(1), rm2(2)-8d0*rmg(2),
     1                 rm2max(2), (r(2)-0.3d0)/0.7d0, W, m2sq)) RETURN 
         endif
      endif

c m_3^2 is generated via tan mapping out of r(3)
      rm2min(3) = 1.0d-5
      rm2max(3) = (sqrt(p_in(4)) - sqrt(m1sq) - sqrt(m2sq))**2
      if (NW.eq.1) then
         if (rm2max(3).lt.rm2(3)) then
            W = 0.0d0
            return
         endif
         m3sq = rm2(3)
         W = W * rmg(3)/2.0d0
      else !NW.ne.1
         if (r(3).lt.0.3d0) then
            W=W/0.3d0
            if (.not.NonResonance(rm2min(3), rm2(3)-8d0*rmg(3),
     1                              r(3)/0.3d0, W, m3sq)) RETURN 
         else
            W=W/0.7d0
            if (.not.Resonance(rm2(3), rmg(3), rm2(3)-8d0*rmg(3),
     1                 rm2max(3), (r(3)-0.3d0)/0.7d0, W, m3sq)) RETURN 
         endif
      endif

c Select lightest mi^2:
      if (m1sq.lt.m2sq .and. m1sq.lt.m3sq) then
         mlight = m1sq
         if(m2sq.lt.m3sq) then
            mmiddle = m2sq
            mheavy = m3sq
            ident = 123
         else
            mmiddle = m3sq
            mheavy = m2sq
            ident = 132
         endif
      else if (m2sq.lt.m1sq .and. m2sq.lt.m3sq) then
         mlight = m2sq
         if(m1sq.lt.m3sq) then
            mmiddle = m1sq
            mheavy = m3sq
            ident = 213
         else
            mmiddle = m3sq
            mheavy = m1sq
            ident = 231
         endif
      else if (m3sq.lt.m1sq .and. m3sq.lt.m2sq) then
         mlight = m3sq
         if(m1sq.lt.m2sq) then
            mmiddle = m1sq
            mheavy = m2sq
            ident = 312
         else
            mmiddle = m2sq
            mheavy = m1sq
            ident = 321
         endif
      endif

c generate two invariant masses (qi + qj)^2  = mij^2 **********
c m_12^2 = (qlight + qmiddle)^2 out of r(4)
      rm2min12 = max(mmiddle+mlight,1.0d-5)
      rm2max12 = min(p_in(4)-mheavy,p_in(4)+mheavy-2.0d0*sqrt(p_in(4)*mheavy))
      if ( r(4) .le. 0.1d0 ) then
         rm2max12 = min(rm2max12, rm2max1)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min12, rm2max12,
     1                              r(4)/0.1d0, W, m12sq)) return
         W = W/0.1d0
      else if ( r(4) .le.  0.7d0 ) then
         rm2min12 = max(rm2min12, rm2max1 )
         rm2max12 = min(rm2max12, rm2max2)
         if (.not. Resonance(xm2(6), xmg(6), rm2min12, rm2max12,
     1                      (r(4)-0.1d0)/0.6d0, W, m12sq)) return
         W = W/0.6d0
      else
         rm2min12 = max(rm2min12, rm2max2 )
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min12, rm2max12,
     1                      (r(4)-0.7d0)/0.3d0, W, m12sq)) return
         W = W/0.3d0
      endif

c m_23^2 = (qlight + qheavy)^2 out of r(5)
      rm2max23 = min(p_in(4)-mmiddle,
     1               p_in(4)+mmiddle-2.0d0*sqrt(p_in(4)*mmiddle) )
      rm2min23 = max(mheavy+mmiddle-m12sq+2.0d0*sqrt(p_in(4)*mlight),
     1               mheavy+mlight, 1.0d-5)! ->2,3,5)
      if ( r(5) .le. 0.1d0 ) then
         rm2max23 = min(rm2max23, rm2max1)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min23, rm2max23,
     1                              r(5)/0.1d0, W, m23sq)) return
         W = W/0.1d0
      else if ( r(5) .le.  0.4d0 ) then
         rm2max23 = min(rm2max23, rm2max2)
         rm2min23 = max(rm2min23, rm2max1)
         if (.not. Resonance(xm2(6), xmg(6), rm2min23, rm2max23,
     1                      (r(5)-0.1d0)/0.3d0, W, m23sq)) return
         W = W/0.3d0
      else
         rm2min23 = max(rm2min23, rm2max2)
         if (.not. Resonance(xxm**2, xxm*xxwidth, rm2min23, rm2max23,
     1                      (r(5)-0.4d0)/0.6d0, W, m23sq)) return
         W = W/0.6d0
      endif

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + mmiddle - m23sq)/2/M
      E3 = (p_in(4) + mheavy - m12sq)/2/M
      E2 = M - E1 - E3
      if (((E1**2 - mmiddle) .le. 1D-10) .or. 
     -     ((E3**2 - mheavy) .le. 1D-10)) then
         w = 0D0
         return
      else
         Abs1 = sqrt(E1**2 - mmiddle)
         Abs3 = sqrt(E3**2 - mheavy)
      end if

c calculate now C and D
      D = -(mmiddle-mlight+mheavy+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)
     1                                      /2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
         return
      endif

c now generate Euler angles (rd(6), rd(7), rd(8))
      alpha = R(6) * twopi
      W = W*twopi
      cosb = -1.0d0 + R(7)*2.0d0
      W = W*2.0d0
      gamma = r(8) * twopi
      W = W*twopi
      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      if ((1 - cosb**2) .le. 1D-10) then
         w = 0D0
         return
      else
         sinb = sqrt(1 - cosb**2)
      end if
      sing = sin(gamma)

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
      ql2(0) = P_in(0) - ql1(0) - ql3(0)
      ql2(1) = P_in(1) - ql1(1) - ql3(1)
      ql2(2) = P_in(2) - ql1(2) - ql3(2)
      ql2(3) = P_in(3) - ql1(3) - ql3(3)

c order q1, q2, q3, m1, m2, m3 in original order:
      if (ident.eq.123) then
         do mu = 0,3
            q1_out(mu) = ql2(mu)
            q2_out(mu) = ql1(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (ident.eq.132) then
         do mu = 0,3
            q1_out(mu) = ql2(mu)
            q2_out(mu) = ql3(mu)
            q3_out(mu) = ql1(mu)
         enddo
      else if (ident.eq.213) then
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q2_out(mu) = ql2(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (ident.eq.231) then
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q2_out(mu) = ql2(mu)
            q3_out(mu) = ql1(mu)
         enddo
      else if (ident.eq.312) then
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q2_out(mu) = ql3(mu)
            q3_out(mu) = ql2(mu)
         enddo
      else if (ident.eq.321) then
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q2_out(mu) = ql1(mu)
            q3_out(mu) = ql2(mu)
         enddo
      endif

c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      mZ1sq = m1sq
      mZ2sq = m2sq
      mZ3sq = m3sq

      ThreeBodyDecay3Higgs = .true.

      RETURN
      END


***************************************************************************
      LOGICAL FUNCTION ThreeBodyDecayPlain(R, P_in, q1_out, q2_out, 
     &                  q3_out, m1sq, m2sq, m3sq, w, NW)
****************************************************************************
*     Calculates the three-body decay P_in -> q1_out + q2_out + q3_out
*
*    INPUT
*       R(8)           : eight random numbers between 0 and 1
*       P_in           : 4-momentum of mother particle
*       w              : phasespace weight
*       NW             : Narrow width (NW = 1) or full calculation (NW = 0)
*
*    OUTPUT
*       w                      : new phasespace weight
*       q1_out, q2_out, q3_out : 4-momenta of daughter particles
*       m1sq, m2sq, m3sq       : mass squared of outgoing particles
*       ThreeBodyDecayPlain    : .true.,  if everything is ok 
*                                .false., if kinematics not ok -> weight = 0
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION RM2, RMG, RM2MIN, RM2MAX, S, PTC
      DOUBLE PRECISION RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN
      DOUBLE PRECISION RM2MAX1, RM2MAX2
      COMMON /TRIBOSBLIPSI/ RM2(3), RMG(3), RM2MIN(3), RM2MAX(3),
     1                RM2MIN12, RM2MAX12, RM2MIN23, RM2MAX23, M2MIN,
     2                S, PTC, RM2MAX1, RM2MAX2

      DOUBLE PRECISION CLR, XM2, XMG, B, VE, A
      COMMON /BKOPOU/ CLR(4,5,-1:1), XM2(6), XMG(6),
     1                B(6,6,6), VE(4,5), A(4,5)

      INTEGER MU, NW

      DOUBLE PRECISION m1sq, m2sq, m3sq, m12sq, m23sq
      DOUBLE PRECISION P_in(0:4), q1_out(0:3), q2_out(0:3), q3_out(0:3)


      DOUBLE PRECISION M, E1, E2, E3, Abs1, Abs3, C, D
      DOUBLE PRECISION alpha, gamma, cosa, cosb, cosg, sina, sinb, sing

      DOUBLE PRECISION R(8)
      DOUBLE PRECISION w
      DOUBLE PRECISION twopi
      PARAMETER (twopi = 2d0*3.141592653589793d0)

c data to be useed in histogram routine
      double precision mZ1sq, mZ2sq, mZ3sq
      Common /histdata1/ m12sq, m23sq, mZ1sq, mZ2sq, mZ3sq, alpha, cosb, gamma

c symmetrisation of q1 and q3
      INTEGER symm
      DOUBLE PRECISION QL1(0:3), QL2(0:3), QL3(0:3)
      DOUBLE PRECISION msqlocal1, msqlocal2, msqlocal3
      data symm /1/
      save symm    

      DOUBLE PRECISION xxwidth, xxm
      PARAMETER ( xxwidth = 10d1, xxm = 45d1 )
      DOUBLE PRECISION x1width, x1m
      PARAMETER ( x1width = 70.0d0, x1m = 90.0d0 )
      DOUBLE PRECISION x2width, x2m
      PARAMETER ( x2width = 70.0d0, x2m = 90.0d0 )

      LOGICAL RESONANCE, NonResonance_KK, NonResonance
      EXTERNAL RESONANCE, NonResonance_KK, NonResonance

      ThreeBodyDecayPlain = .false.

      M = sqrt(p_in(4))

c generate two invariant masses (qi + qj)^2  = mij^2 first **********

c generation of m12sq ---------------------------------
      rm2max12 = p_in(4)
      rm2min12 = 1.0d-5
      if (r(1).le. 0.1d0) then
        if (.not. Resonance(rm2(1), rmg(1), rm2min12, rm2(1)+10d0*rmg(1),
     1                               r(1)/0.1d0, W, m12sq)) return
         W = W/0.1d0
      else if (r(1) .le. 0.3d0) then
         if (.not. NonResonance_KK(x1m**2, x1m*x1width, 2.0d0, rm2(1)+10d0*rmg(1),
     1                 4.d0*rm2(1)-2.d0*rmg(1), (r(1)- 0.1d0)/0.2d0, W, m12sq)) return
         W = W/0.2d0
      else
         if (.not. NonResonance_KK(x2m**2, x2m*x2width, 2.0d0, 4.d0*rm2(1)-2.d0*rmg(1),
     1                  rm2max12, (r(1)- 0.3d0)/0.7d0, W, m12sq)) return
         W = W/0.7d0
      endif

c generation of m23sq next ---------------------------------
      rm2max23 = p_in(4)
      rm2min23 = 1.0d-5
      if (r(2).le. 0.1d0) then
        if (.not. Resonance(rm2(1), rmg(1), rm2min23, rm2(1)+10d0*rmg(1),
     1                               r(2)/0.1d0, W, m23sq)) return
         W = W/0.1d0
      else if (r(2) .le. 0.3d0) then
         if (.not. NonResonance_KK(x1m**2, x1m*x1width, 2.0d0, rm2(1)+10d0*rmg(1),
     1                 4.d0*rm2(1)-2.d0*rmg(1), (r(2)- 0.1d0)/0.2d0, W, m23sq)) return
         W = W/0.2d0
      else
         if (.not. NonResonance_KK(x2m**2, x2m*x2width, 2.0d0, 4.d0*rm2(1)-2.d0*rmg(1),
     1                  rm2max23, (r(2)- 0.3d0)/0.7d0, W, m23sq)) return
         W = W/0.7d0
      endif

c now generate m1^2 and m3^2 of the massive vector bosons

c m_1^2 is generated via tan mapping out of r(3)
      rm2min(1) = 1.0d-5
      rm2max(1) = min(p_in(4)-m23sq,p_in(4)+m23sq-2.0d0*sqrt(p_in(4)*m23sq),
     1                 m12sq)
      if (rm2max(1).le.0) then
         W = 0.0d0
         return
      endif

      if (NW.eq.1) then
         if (rm2max(1).lt.rm2(1)) then
            W = 0.0d0
            return
         endif
         m1sq = rm2(1)
         W = W * rmg(1)/2.0d0
      else !NW.ne.1

         if(rm2max(1).ge.rm2(1)-16d0*rmg(1)) then
            if (r(3).lt.0.2d0) then
               W=W/0.2d0
               if (.not.NonResonance(rm2min(1), rm2(1)-16d0*rmg(1),
     1                              r(3)/0.2d0, W, m1sq)) RETURN 
            elseif (r(3).lt.0.8d0) then
               W=W/0.6d0
               if (.not.Resonance(rm2(1), rmg(1), rm2(1)-16d0*rmg(1),
     1                 rm2(1)+4d0*rmg(1), (r(3)-0.2d0)/0.6d0, W, m1sq)) RETURN 
            elseif (r(3).lt.0.9d0) then
               W=W/0.1d0
               if (.not.Resonance(rm2(1), rmg(1)*10d00, rm2(1)+4d0*rmg(1),
     1                 rm2(1)+100d0*rmg(1), (r(3)-0.8d0)/0.1d0, W, m1sq)) RETURN 
            else
               W=W/0.1d0
               if (.not.Resonance(rm2(1), rmg(1)*150d0, rm2(1)+100d0*rmg(1),
     1                 rm2max(1), (r(3)-0.9d0)/0.1d0, W, m1sq)) RETURN 
            endif
         else
            if (.not. NonResonance(rm2min(1), rm2max(1),
     1                              r(3), W, m1sq)) return
         endif

      endif !NW

c m_3^2 is generated via tan mapping out of r(4)
      rm2min(3) = 1.0d-5
      rm2max(3) = min(p_in(4)-m12sq, p_in(4)+m12sq-2.0d0*sqrt(p_in(4)*m12sq), 
     1                m12sq+m23sq-m1sq, m23sq, p_in(4)-m1sq)

      if (rm2max(3).le.0) then
         W = 0.0d0
         return
      endif

      if (NW.eq.1) then
         if (rm2max(3).lt.rm2(3)) then
            W = 0.0d0
            return
         endif
         m3sq = rm2(3)
         W = W * rmg(3)/2.0d0
      else !NW.ne.1

         if(rm2max(3).ge.rm2(3)-16d0*rmg(3)) then
            if (r(4).lt.0.2d0) then
               W=W/0.2d0
               if (.not.NonResonance(rm2min(3), rm2(3)-16d0*rmg(3),
     1                              r(4)/0.2d0, W, m3sq)) RETURN 
            elseif (r(4).lt.0.8d0) then
               W=W/0.6d0
               if (.not.Resonance(rm2(3), rmg(3), rm2(3)-16d0*rmg(3),
     1                 rm2(3)+4d0*rmg(3), (r(4)-0.2d0)/0.6d0, W, m3sq)) RETURN 
            elseif (r(4).lt.0.9d0) then
               W=W/0.1d0
               if (.not.Resonance(rm2(3), rmg(3)*10, rm2(3)+4d0*rmg(3),
     1                 rm2(3)+100d0*rmg(3), (r(4)-0.8d0)/0.1d0, W, m3sq)) RETURN 
            else
               W=W/0.1d0
               if (.not.Resonance(rm2(3), rmg(3)*150d0, rm2(3)+100d0*rmg(3),
     1                 rm2max(3), (r(4)-0.9d0)/0.1d0, W, m3sq)) RETURN 
            endif
         else
            if (.not. NonResonance(rm2min(3), rm2max(3),
     1                              r(4), W, m3sq)) return
         endif

      endif

c caculate E1, E2, E3, Abs1, Abs2 from that 
      E1 = (p_in(4) + m1sq - m23sq)/2.0d0/M
      E3 = (p_in(4) + m3sq - m12sq)/2.0d0/M
      E2 = M - E1 - E3
      if (((E1**2 - m1sq) .le. 1D-10) .or. 
     -     ((E3**2 - m3sq) .le. 1D-10)) then
         w = 0D0
         return
      else
         Abs1 = sqrt(E1**2 - m1sq)
         Abs3 = sqrt(E3**2 - m3sq)
      end if

c m_2^2 is generated via tan mapping out of r(5)
      rm2min(2) = max(-Abs1**2-Abs3**2-2.0d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1            1.0d-5)
      rm2max(2) = min(-Abs1**2-Abs3**2+2.0d0*sqrt(Abs1**2*Abs3**2)+E2**2,
     1          E2**2, P_in(4)-m1sq-m3sq, m12sq-m1sq, m23sq-m3sq)

      if (rm2max(2)-rm2min(2).le.0) then
         W = 0.0d0
         return
      endif

      if (NW.eq.1) then
         if (rm2max(2).lt.rm2(2) .or. rm2min(2).gt.rm2(2)) then
            W = 0.0d0
            return
         endif
         m2sq = rm2(2)
         W = W * rmg(2)/2.0d0
      else !NW.ne.1

         if(rm2max(2).ge.rm2(2)-16d0*rmg(2)) then
            if (r(5).lt.0.2d0) then
               W=W/0.2d0
               if (.not.NonResonance(rm2min(2), rm2(2)-16d0*rmg(2),
     1                              r(5)/0.2d0, W, m2sq)) RETURN 
            elseif (r(5).lt.0.8d0) then
               W=W/0.6d0
               if (.not.Resonance(rm2(2), rmg(2), rm2(2)-16d0*rmg(2),
     1                 rm2(2)+4d0*rmg(2), (r(5)-0.2d0)/0.6d0, W, m2sq)) RETURN 
            elseif (r(5).lt.0.9d0) then
               W=W/0.1d0
               if (.not.Resonance(rm2(2), rmg(2)*10d0, rm2(2)+4d0*rmg(2),
     1                 rm2(2)+100d0*rmg(2), (r(5)-0.8d0)/0.1d0, W, m2sq)) RETURN 
            else
               W=W/0.1d0
               if (.not.Resonance(rm2(2), rmg(2)*150d0, rm2(2)+100d0*rmg(2),
     1                 rm2max(2), (r(5)-0.9d0)/0.1d0, W, m2sq)) RETURN 
            endif
         else
            if (.not. NonResonance(rm2min(2), rm2max(2),
     1                              r(5), W, m2sq)) return
         endif

      endif!NW

c now generate Euler angles (rd(6), rd(7), rd(8))
      alpha = R(6) * twopi
      W = W*twopi
      cosb = -1.0d0 + R(7)*2.0d0
      W = W*2.0d0
      gamma = R(8) * twopi
      W = W*twopi

      cosa = cos(alpha)
      cosg = cos(gamma)
      sina = sin(alpha)
      if ((1D0 - cosb**2) .le. 1D-10) then
         w = 0D0
         return
      else
         sinb = sqrt(1 - cosb**2)
      end if
      sing = sin(gamma)

c calculate now C and D
      D = -(m1sq-m2sq+m3sq+p_in(4)-2.0d0*M*(E1+E3)+2.0d0*E1*E3)/2/Abs1/Abs3
      if (D**2.lt.1) then
         C = sqrt(1 - D**2)
      else
         W = 0.0d0
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
      ql2(0) = P_in(0) - ql1(0) - ql3(0)
      ql2(1) = P_in(1) - ql1(1) - ql3(1)
      ql2(2) = P_in(2) - ql1(2) - ql3(2)
      ql2(3) = P_in(3) - ql1(3) - ql3(3)


c symmetrisation of q1, q2, q3
      msqlocal1 = m1sq
      msqlocal2 = m2sq
      msqlocal3 = m3sq
      symm = mod(symm,6)+1
      if (symm.eq.1) then ! 123
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q2_out(mu) = ql2(mu)
            q3_out(mu) = ql3(mu)
         enddo
      else if (symm.eq.2) then !132
         do mu = 0,3
            q1_out(mu) = ql1(mu)
            q2_out(mu) = ql3(mu)
            q3_out(mu) = ql2(mu)
         enddo
         m2sq = msqlocal3
         m3sq = msqlocal2
      else if (symm.eq.3) then !213
         do mu = 0,3
            q1_out(mu) = ql2(mu)
            q2_out(mu) = ql1(mu)
            q3_out(mu) = ql3(mu)
         enddo
         m1sq = msqlocal2
         m2sq = msqlocal1
      else if (symm.eq.4) then !231
         do mu = 0,3
            q1_out(mu) = ql2(mu)
            q2_out(mu) = ql3(mu)
            q3_out(mu) = ql1(mu)
         enddo
         m1sq = msqlocal2
         m2sq = msqlocal3
         m3sq = msqlocal1
      else if (symm.eq.5) then !312
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q2_out(mu) = ql1(mu)
            q3_out(mu) = ql2(mu)
         enddo
         m1sq = msqlocal3
         m2sq = msqlocal1
         m3sq = msqlocal2
      else if (symm.eq.6) then !321
         do mu = 0,3
            q1_out(mu) = ql3(mu)
            q2_out(mu) = ql2(mu)
            q3_out(mu) = ql1(mu)
         enddo
         m1sq = msqlocal3
         m3sq = msqlocal1
      endif

c Additional weight factors ******************************************
      W = W/twopi**3
      W = W/32.0d0/P_in(4)

      mZ1sq = m1sq
      mZ2sq = m2sq
      mZ3sq = m3sq

      ThreeBodyDecayPlain = .true.

      RETURN
      END


****************************************************************************
      LOGICAL FUNCTION TwoToJetsPlusX(N,R,RN,ECM,PTmin,XQ2,
     1                                K1,K2,X1,X2,Q,P,weight)
****************************************************************************
*     Calculates the phasespace for the 2 -> N massless jets + X
*     X is a pseudo-particle with arbitrary mass>0
*
*    INPUT
*       N      : number of massless jets      
*       R      : array with 3*N random numbers between 0 and 1
*       RN     : one additional random number 
*       ECM    : center-of-mass energy square
*       PTmin  : minimum PT of jets
*       XQ2    : Q**2 of X-particle
*       weight : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of X-particle
*       P      : 4-vectors for the N massless jets
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************
      IMPLICIT NONE
      INTEGER*4 N, N1
      PARAMETER (N1 = 5)
      DOUBLE PRECISION TPI,PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )
      DOUBLE PRECISION R(*),RN,ECM,PTmin,XQ2,
     1       K1(0:3),K2(0:3),Q(0:3),P(0:3,N1),
     2       X1,X2,weight
      
      INTEGER*4 I,J
      DOUBLE PRECISION S,FAC,DELY,XMIN,XMAX,X,DELX,YCM,YMIN,YMAX,YBAR
      DOUBLE PRECISION PCM(0:3), PNJ(0:4), Y(N1), PT(N1), PHI(N1)
      DOUBLE PRECISION SINHY, COSHY, M2, M2min, QV, XMJ, PJT2, SHA, ACM
      DOUBLE PRECISION SUMPST, Q0ST, RSHAT
      PARAMETER (m2min = 0.1d0**2)

      double precision dotrr
      external dotrr

      S = ECM**2
      FAC = 1d0/(16d0*PI**3)
      weight = weight*TPI

C  Reset local variables and incoming momenta      
      DO I=0,3
         PCM(I) = 0d0
         PNJ(I) = 0d0
         K1(I) = 0d0
         K2(I) = 0d0
      ENDDO
      PNJ(4) = 0d0

C  Get 4-momenta of the N massless jets
      DO I = 1,N
         weight = weight*FAC
         J = 3*(I-1)
         
         if (R(J+2).gt.0.2) then
            XMIN = 2D0/ECM
            XMAX = 1/max(PTmin,10d0)
            DELX = XMAX-XMIN
            X = XMIN + (R(J+2)-0.2)/0.8d0*DELX
            PT(I) = 1/X
            weight = weight/0.8d0*DELX*PT(I)**3
         else
            xmin = 0
            xmax = max(PTmin,10d0)
            DELX = XMAX-XMIN
            PT(i) = XMIN + R(J+2)/0.2d0*DELX
            weight = weight/0.2d0*DELX*PT(I)
         endif
         YMAX = ECM/(2d0*PT(I))
         YMAX = LOG(YMAX+SQRT(YMAX**2-1d0))

         DELY = 2*YMAX

         IF (I.LE.3) THEN
            Y(I) = -YMAX+R(J+1)*DELY
         ELSE
C  fix random number to 0 or 1 for delta-y = 0 with previous jet: R=0 or 1 
C  corresponds to lowest dijet invariant mass
            IF (Y(I-1).LT.0D0) THEN
               IF (Y(I-1).LT.-YMAX) THEN
                  Y(I) = -YMAX + R(J+1)*DELY
               ELSEIF (Y(I-1).LT.0d0) THEN
                  Y(I) = Y(I-1) + R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ELSE
                  Y(I) =  - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ENDIF
            ELSE
               IF (Y(I-1).GT.YMAX) THEN
                  Y(I) = YMAX - R(J+1)*DELY
               ELSEIF (Y(I-1).GT.0D0) THEN
                  Y(I) = Y(I-1) - R(J+1)*DELY
                  IF (Y(I).LT.-YMAX) Y(I) = 2*YMAX + Y(I)
               ELSE
                  Y(I) = R(J+1)*DELY
                  IF (Y(I).GT.YMAX) Y(I) = -2*YMAX + Y(I)
               ENDIF
            ENDIF
         ENDIF

         SINHY = SINH(Y(I))
         COSHY = SQRT(1+SINHY**2)
         weight = weight*DELY

         P(0,I) = PT(I)*COSHY

         IF (I.EQ.1) THEN
            PHI(I) = PI*(2*RN-1)
         ELSE
            PHI(I) = PHI(I-1) + TPI*R(J+3)
            if (phi(i).gt.pi) phi(i) = phi(i)-tpi
         ENDIF

         weight = weight*TPI

         P(1,I) = PT(I)*COS(PHI(I))
         P(2,I) = PT(I)*SIN(PHI(I))
         P(3,I) = PT(I)*SINHY

         DO J = 0,3
            PNJ(J) = PNJ(J)+P(J,I)
         ENDDO
      ENDDO

C  determine N-jet invariant mass
      IF (N.GT.1) THEN
         PNJ(4) = PNJ(0)**2-PNJ(1)**2-PNJ(2)**2-PNJ(3)**2
      ELSE
         PNJ(4) = 1D-18
      ENDIF
      IF (PNJ(4).GE.S) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF

c  check invariant masses of all parton pairs: eliminate events which 
c  may cause numerical problems 
      do I = 1,N-1
         do J = I+1,N
            m2 = 2*dotrr(p(0,I),p(0,J))
            if (m2.lt.m2min) then
               weight = 0d0
               TwoToJetsPlusX = .false.
               RETURN
            endif
         enddo
      enddo
      
C  Next get the center of mass rapidity YCM
      QV = SQRT(XQ2)
      XMJ = SQRT(PNJ(4))
      if (PNJ(0).le.PNJ(3)) then
        ybar = 100
      elseif (PNJ(0).le.-PNJ(3)) then
        ybar = -100
      else
        YBAR = 0.5D0*LOG( (PNJ(0)+PNJ(3))/(PNJ(0)-PNJ(3)) )
      endif
      PJT2 = PNJ(1)**2+PNJ(2)**2
      SHA = (S-(QV+XMJ)**2-4*PJT2) / (4*(PNJ(4)+PJT2)) 
      IF (SHA.LE.0D0) THEN
        weight=0d0
        TwoToJetsPlusX = .false.
        RETURN
      ENDIF
      SHA = SQRT(SHA)
      ACM = LOG( SHA+SQRT(SHA**2+1) ) 
      YMIN = YBAR - ACM
      YMAX = YBAR + ACM

      DELY = YMAX-YMIN
      YCM = YMIN + R(3)*DELY
      SINHY = SINH(YCM)
      COSHY = SQRT(1+SINHY**2)
      weight = weight*DELY

C  now get Q and SHAT etc.
      SUMPST = PNJ(1)**2+PNJ(2)**2+(PNJ(3)*COSHY-PNJ(0)*SINHY)**2
      Q0ST = SQRT(XQ2+SUMPST)
      RSHAT = Q0ST + SQRT(PNJ(4)+SUMPST)
      PCM(0) = RSHAT*COSHY
      PCM(3) = RSHAT*SINHY

      DO J = 0,3
        Q(J) = PCM(J)-PNJ(J)
      ENDDO

      X1 = (PCM(0)+PCM(3))/ECM
      X2 = (PCM(0)-PCM(3))/ECM
      IF (X1.GE.1.0D0 .OR. X2.GE.1.0D0) THEN
         weight=0d0
         TwoToJetsPlusX = .false.
         RETURN
      ENDIF
      K1(0) =  X1*ECM/2
      K1(3) =  K1(0)
      K2(0) =  X2*ECM/2
      K2(3) = -K2(0)

c  check initial-final momentum transfer for all partons: eliminate events 
c  which may cause numerical problems 
      do I = 1,N         
         m2 = 2*min( p(0,I)*k1(0)-p(3,I)*k1(3),
     1               p(0,I)*k2(0)-p(3,I)*k2(3) )
         if (m2.lt.m2min) then
            weight = 0d0
            TwoToJetsPlusX = .false.
            RETURN
         endif
      enddo

C  insert Jacobian for dtau * delta(PCM(0)-Ql(0,0)-PNJ(0))/(2 Ql(0,0))
      weight = weight * RSHAT/(S*Q0ST)
      
C  finally put in flux factor for the parton collisions and conversion to fb
      weight = weight / (2*RSHAT**2) * 3.89379304D+11

      TwoToJetsPlusX = .true.
      RETURN
      END

****************************************************************************
      SUBROUTINE TwoToOne(R,ECM,Qsqr,K1,K2,X1,X2,Q,weight)
****************************************************************************
*     Calculates the phasespace for 2 -> 1 massive
*
*    INPUT
*       R            : random number
*       ECM          : center-of-mass energy square
*       Qsqr         : Q**2 of the outgoing particle
*       weight       : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of outgoing particle
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION R, ECM, Qsqr
      DOUBLE PRECISION K1(0:3), K2(0:3), Q(0:3)
      DOUBLE PRECISION X1, X2, weight
      DOUBLE PRECISION TPI, PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )

c local variables
      INTEGER*4 mu
      DOUBLE PRECISION s, tau, y, y_border

      S = ecm**2

c Generate y, tau -> x1, x2, shat from R first **************
      tau = Qsqr/s  
      weight = weight / s 
      y_border = log(1/tau)
      y = y_border*( R - 0.5d0 )
      weight = weight * y_border !(prop density)^-1 of y

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

      Q(1) = 0d0
      Q(2) = 0d0
      do mu = 0, 3, 3
         Q(mu) = k1(mu)+k2(mu)
      enddo

      weight = weight * 0.5d0/Qsqr        ! Flux
      weight = weight * TPI                  ! (2*pi)^(4-3) from delta- and ps-int
      weight = weight * 3.89379304d11        ! conversion into fb

      return
      end


****************************************************************************
      SUBROUTINE TwoToTwo(R,RN,ECM,Qsqr1,Qsqr2,K1,K2,X1,X2,Q1,Q2,weight)
****************************************************************************
*     Calculates the phasespace for 2 -> 2 massive
*
*    INPUT
*       R            : array with 3 random numbers between 0 and 1
*       RN           : one additional random number 
*       ECM          : center-of-mass energy square
*       Qsqr1, Qsqr2 : Q**2 of the two outgoing particles
*       weight       : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q1,Q2  : 4-vectors of outgoing particles
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************

      IMPLICIT NONE

      DOUBLE PRECISION R(*), RN, ECM, Qsqr1, Qsqr2
      DOUBLE PRECISION K1(0:3), K2(0:3), Q1(0:3), Q2(0:3)
      DOUBLE PRECISION X1, X2, weight
      DOUBLE PRECISION TPI, PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )

c local variables
      INTEGER*4 mu
      DOUBLE PRECISION k(0:3)
      DOUBLE PRECISION s, ptmax, y2_border, temp, pt, ET, y_border, y, y2
      DOUBLE PRECISION phi, sinhy2, coshy2, abs_p, pL, tau
      DOUBLE PRECISION A, ptc

c test variables
      COMMON /ptdistr/ pt

      S = ecm**2

c generate y, p_T, y*, phi from plain distributed R(1),R(2),R(3), rn

      phi = rn * tpi - tpi/2d0

c ptc is the 'changeover point' of the phasespace generation - 
c below we have a constant value, and above it's more complicated
      ptc = 10.0d0
      pTmax = ecm/2.0d0
c      ptmax = 0.5d0*sqrt((s-(Qsqr1+Qsqr2))**2 - 4d0*Qsqr1*Qsqr2)/ecm
      A = pTmax/ptc*1.0d0/(2*pTmax-ptc)
      if (r(1) .le. A * ptc) then
         pT = r(1)*1.0d0/A
         Weight = Weight/A
      else
         pT = pTmax*ptc*(1-A*ptc)/(pTmax-A*ptc**2-r(1)
     &               *(pTmax-ptc))
         Weight = Weight/A/ptc**2*pT**2
      endif
c pT-distribution not continuous------------------------------
c      if (r(1).le.0.2) then
c         pTmax = ptc
c         pTmin = 0.0d0       
c         pT = r(1)*(pTmax-ptmin)/0.2d0 + ptmin
c         Weight= Weight/0.2d0
c      else
c         pTmin = 2.0d0/ecm
c         pTmax = 1/ptc
c         pT = pTmin + (r(1)-0.2)/0.8d0*(pTmax-pTmin)
c         pT = 1/pT
c         Weight = Weight/0.8d0*pT**2
c      endif
c      Weight = Weight*(ptmax-ptmin)


      ET = sqrt(pT**2+Qsqr1)
      temp = (s+Qsqr1-Qsqr2)/(2*ecm*ET)

c calculate y*
      if ((temp**2 - 1D0) .gt. 0D0) then
         y2_border = log( temp+sqrt(temp**2-1) )
         y2 = y2_border*2d0*( R(2) - 0.5d0 ) ! y*
         weight = weight * y2_border * 2d0 !(prop density)^-1 of y2
         
         sinhy2 = sinh(y2)
         coshy2 = sqrt(sinhy2**2+1)
         abs_p = sqrt(pT**2*coshy2**2+Qsqr1*sinhy2**2)
         pL = sinhy2*ET
         tau=((coshy2*ET+sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))**2)/s
         weight = weight * pT/(sqrt(tau*s))
         weight = weight/(sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))
         
         y_border = log(1/tau)
         y = y_border*( R(3) - 0.5d0 ) !y
         weight = weight * y_border !(prop density)^-1 of y
         
         
c------------------------------------------------------------------------
c generate incoming momenta in lab frame + x1 x2

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)
         
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
         
c Momenta of V1

         q1(0) = sqrt( pT**2 + pL**2 + Qsqr1 )
         q1(1) = pT * cos(phi)
         q1(2) = pT * sin(phi)
         q1(3) = pL
         
         
c Momenta of V2

         q2(0) = sqrt( pT**2 + pL**2 + Qsqr2 )
         q2(1) = - pT * cos(phi)
         q2(2) = - pT * sin(phi)
         q2(3) = - pL
         
c boost in lab frame

         call boostn( q1(0), K, q1(0) )
         call boostn( q2(0), K, q2(0) )

         weight = weight * 1d0/tpi/4d0/s
         weight = weight * 3.89379304d11 !conversion into fb

      else
         weight = 0D0
      end if


      RETURN
      END


****************************************************************************
c$$$      SUBROUTINE TwoToTwo2(R,RN,ECM,Qsqr1,Qsqr2,K1,K2,X1,X2,Q1,Q2,weight)
c$$$****************************************************************************
c$$$*     Calculates the phasespace for 2 -> 2 massive
c$$$*
c$$$*    INPUT
c$$$*       R            : array with 3 random numbers between 0 and 1
c$$$*       RN           : one additional random number 
c$$$*       ECM          : center-of-mass energy square
c$$$*       Qsqr1, Qsqr2 : Q**2 of the two outgoing particles
c$$$*       weight       : phasespace weight
c$$$*
c$$$*    OUTPUT
c$$$*       K1,K2  : 4-vectors for 2 incoming partons
c$$$*       X1,X2  : corresponding values of Feynman x
c$$$*       Q1,Q2  : 4-vectors of outgoing particles
c$$$*       weight : new phasespace weight (multiplied by Jacobi factor)
c$$$****************************************************************************
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION R(*),RN,ECM,Qsqr1,Qsqr2,
c$$$     1       K1(0:3),K2(0:3),Q1(0:3),Q2(0:3),
c$$$     2       X1,X2,weight
c$$$      DOUBLE PRECISION TPI,PI, ptc
c$$$      PARAMETER (PI=3.141592653589793d0)
c$$$      PARAMETER (TPI=2D0*PI )
c$$$
c$$$      integer*4 mu
c$$$      double precision k(0:3)
c$$$      double precision s,ptmax,y2_border,temp,pt,ET,y_border,y,y2
c$$$      double precision phi,sinhy2,coshy2,abs_p,pL,tau
c$$$      double precision A
c$$$
c$$$c test variables
c$$$      common /ptdistr/ pt
c$$$
c$$$
c$$$      S = ecm**2
c$$$
c$$$c generate y, p_T, y*, phi from plain distributed R(1),R(2),R(3), rn
c$$$
c$$$      phi = rn * tpi - tpi/2d0
c$$$
c$$$      pTmax = ecm/2.0d0
c$$$c      ptmax = 0.5d0*sqrt((s-(Qsqr1+Qsqr2))**2 - 4d0*Qsqr1*Qsqr2)/ecm
c$$$      pTC = 10.0d0 ! should be the same as above in twototwo()
c$$$c pT distribution ----constant - 1/pt^2 decrease
c$$$c      A = pTmax/pTc*1.0d0/(2*pTmax-pTc)
c$$$c      if (r(1) .le. A * ptc) then
c$$$c         pT = r(1)*1.0d0/A
c$$$c         Weight = Weight/A
c$$$c      else
c$$$c         pT = pTmax*pTc*(1-A*pTc)/(pTmax-A*pTc**2-r(1)        !1
c$$$c     &               *(pTmax-pTc))
c$$$c         pT = pTmax*pTc/(2.0d0*pTmax-r(1)*(2.0d0*pTmax-pTc))  !2
c$$$c         Weight = Weight/A/pTc**2*pT**2
c$$$c      endif
c$$$c pT distribution ----linear increase - 1/pT^2 decrease
c$$$c      A = ptmax/(3/2.d0*ptmax-ptc)/ptc**2
c$$$c      if ( r(1) .le. A*ptc**2/2.d0 ) then
c$$$c         pT = sqrt(r(1)*2.d0/A)
c$$$c         Weight = Weight/A/pT
c$$$c      else
c$$$c         pT = pTmax*pTc/(3/2.d0*pTmax-r(1)*(3/2.d0*pTmax-pTc))
c$$$c         Weight = Weight/A/ptc**3*pT**2
c$$$c      endif
c$$$c pT distribution ----linear increase - 1/pT decrease
c$$$c      A = 1.d0/ptc**2/(0.5d0+log(pTmax/pTc))
c$$$c      if ( r(1) .le. ptc**2*A/2.d0 ) then
c$$$c         pT = sqrt(r(1)*2.d0/A)
c$$$c         Weight = Weight/A/pT
c$$$c      else
c$$$c         pT = pTc*(ptmax/ptc)**r(1)*exp(0.5d0*(r(1)-1))
c$$$c         Weight = Weight/A/ptc**2*pT
c$$$c      endif
c$$$c pT distribution ----constant - 1/pt decrease
c$$$      A = 1.d0/ptc/(1+log(pTmax/ptc))
c$$$      if ( r(1) .le. A*ptc ) then
c$$$         pT = r(1)/A
c$$$         Weight = Weight/A
c$$$      else
c$$$         pT = pTc*(ptmax/ptc)**r(1)*exp(r(1)-1)
c$$$         Weight = Weight/A/ptc*pT
c$$$      endif
c$$$
c$$$      ET = sqrt(pT**2+Qsqr1)
c$$$      temp = (s+Qsqr1-Qsqr2)/(2*ecm*ET)
c$$$
c$$$c calculate y*
c$$$      if ((temp**2-1) .gt. 0D0) then
c$$$         y2_border = log( temp+sqrt(temp**2-1) )
c$$$         y2 = y2_border*2d0*( R(2) - 0.5d0 )     ! y*
c$$$         weight = weight * y2_border * 2d0       !(prop density)^-1 of y2
c$$$
c$$$         sinhy2 = sinh(y2)
c$$$         coshy2 = sqrt(sinhy2**2+1)
c$$$         abs_p = sqrt(pT**2*coshy2**2+Qsqr1*sinhy2**2)
c$$$         pL = sinhy2*ET
c$$$         tau=((coshy2*ET+sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))**2)/s
c$$$         weight = weight * pT/(sqrt(tau*s))
c$$$         weight = weight/(sqrt(coshy2**2*ET**2-Qsqr1+Qsqr2))
c$$$     
c$$$
c$$$         y_border = log(1/tau)
c$$$         y = y_border*( R(3) - 0.5d0 )            !y
c$$$         weight = weight * y_border               !(prop density)^-1 of y
c$$$         
c$$$c------------------------------------------------------------------------
c$$$c generate incoming momenta in lab frame + x1 x2
c$$$
c$$$         x1 = sqrt(tau) * exp(y)
c$$$         x2 = sqrt(tau) * exp(-y)
c$$$
c$$$         k1(0) = 0.5d0 * x1 * ecm
c$$$         k1(1) = 0d0
c$$$         k1(2) = 0d0
c$$$         k1(3) = 0.5d0 * x1 * ecm
c$$$
c$$$         k2(0) = 0.5d0 * x2 * ecm
c$$$         k2(1) = 0d0
c$$$         k2(2) = 0d0
c$$$         k2(3) = - 0.5d0 * x2 * ecm
c$$$
c$$$         do mu = 0, 3
c$$$            K(mu) = k1(mu)+k2(mu)
c$$$         enddo
c$$$         
c$$$c Momenta of V1
c$$$         q1(0) = sqrt( pT**2 + pL**2 + Qsqr1 )
c$$$         q1(1) = pT * cos(phi)
c$$$         q1(2) = pT * sin(phi)
c$$$         q1(3) = pL
c$$$      
c$$$
c$$$c Momenta of V2
c$$$         q2(0) = sqrt( pT**2 + pL**2 + Qsqr2 )
c$$$         q2(1) = - pT * cos(phi)
c$$$         q2(2) = - pT * sin(phi)
c$$$         q2(3) = - pL
c$$$
c$$$c boost in lab frame
c$$$         call boostn( q1(0), K, q1(0) )
c$$$         call boostn( q2(0), K, q2(0) )
c$$$
c$$$         weight = weight * 1d0/tpi/4d0/s
c$$$         weight = weight * 3.89379304d11 !conversion into fb
c$$$
c$$$      else
c$$$         weight = 0D0
c$$$      end if
c$$$
c$$$      
c$$$      RETURN
c$$$      END


****************************************************************************
      Subroutine TwoToJetPlusX(R,RN,ECM,PTC,XQ2,K1,K2,X1,X2,Q,P,weight)
****************************************************************************
*     Calculates the phasespace for the 2 -> 1 massless Jet + X
*     X is a pseudo-particle with arbitrary mass>0
*
*    INPUT
*       R      : array with 3*N random numbers between 0 and 1
*       RN     : one additional random number 
*       ECM    : center-of-mass energy square
*       PTC    : minimum PT of jets
*       XQ2    : Q**2 of X-particle
*       weight : phasespace weight
*
*    OUTPUT
*       K1,K2  : 4-vectors for 2 incoming partons
*       X1,X2  : corresponding values of Feynman x
*       Q      : 4-vector of X-particle
*       P      : 4-vectors for the N massless jets
*       weight : new phasespace weight (multiplied by Jacobi factor)
****************************************************************************

      IMPLICIT NONE
      DOUBLE PRECISION R(*),RN,ECM,XQ2,
     1       K1(0:3),K2(0:3),Q(0:3),P(0:3),
     2       X1,X2,weight, ptc
      DOUBLE PRECISION TPI,PI
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (TPI=2D0*PI )

c test common block for plots
      double precision ptdistr
      common /ptdistribution/ ptdistr

c local variables
      integer*4 mu
      double precision k(0:3), m2min
      double precision s,ptmax,y2_border,temp,pt,ET,y_border,y,y2
      double precision phi,sinhy2,coshy2,abs_p,pL,tau

      double precision z,lambda
      parameter(lambda = 1.7d0)
      double precision gamma, kappa
      parameter(gamma = 2.5d0)

      kappa = 0.5d0*gamma/atan(gamma)


c cutoff for numerical purposes
      m2min = 0.01d0**2

      S = ecm**2

c generate y, p_T, y*, phi from plain distributed R(1),R(2),R(3), rn

      phi = rn * tpi - tpi/2d0

      pTmax = ecm/2.0d0
c      ptmax = 0.5d0*sqrt((s-(xq2))**2 - 4d0*xq2)/ecm
c      pTC = 40.0d0
c old pt distribution in my PS generator
      z = pTmax/pTc*1.0d0/(2*pTmax-pTc)
      if (r(1) .le. z * ptc) then !constant pT
         pT = r(1)*1.0d0/z
         Weight = Weight/z
      else
!         pT = pTmax*pTc*(1-z*pTc)/(pTmax-z*pTc**2-r(1)
!     &               *(pTmax-pTc))
         pT = pTmax*pTc/(2.0d0*pTmax-r(1)*(2.0d0*pTmax-pTc))
         Weight = Weight/z/pTc**2*pT**2
      endif

c      z = ptmax/(3/2.d0*ptmax-ptc)/ptc**2
c      if ( r(1) .le. z*ptc**2/2.d0 ) then !linear increase
c         pT = sqrt(r(1)*2.d0/z)
c         Weight = Weight/z/pT
c      else                                ! 1/pT**2
c         pT = pTmax*pTc/(3/2.d0*pTmax-r(1)*(3/2.d0*pTmax-pTc))
c         Weight = Weight/z/ptc**3*pT**2
c      endif

      ptdistr = pt

      ET = sqrt(pT**2)
      temp = (s-xq2)/(2*ecm*ET)


c calculate y*
      if ((temp**2-1D0) .gt. 0D0) then
         y2_border = log( temp+sqrt(temp**2-1) )
         y2 = y2_border*2d0*( R(2) - 0.5d0 )        ! y*
         weight = weight * y2_border * 2d0          !(prop density)^-1 of y2

         sinhy2 = sinh(y2)
         if (((sinhy2**2+1D0) .le. 1D-10) .or. 
     -        ((pT**2*coshy2**2+sinhy2**2) .le. 1D-10)) then
            weight = 0D0
            return
         end if
         coshy2 = sqrt(sinhy2**2+1)
         if (((coshy2**2*ET**2+xq2) .le. 1D-10) .or. 
     -        ((pT**2*coshy2**2+sinhy2**2) .le. 1D-10)) then
            weight = 0D0
            return
         end if
         abs_p = sqrt(pT**2*coshy2**2+sinhy2**2)
         pL = sinhy2*ET
         tau=((coshy2*ET+sqrt(coshy2**2*ET**2+xq2))**2)/s
         weight = weight * pT/(sqrt(tau*s))
         weight = weight/(sqrt(coshy2**2*ET**2+xq2))


         y_border = log(1/tau)
         y = y_border*( R(3) - 0.5d0 ) !y
         weight = weight * y_border !(prop density)^-1 of y
         

c------------------------------------------------------------------------
c generate incoming momenta in lab frame + x1 x2

         x1 = sqrt(tau) * exp(y)
         x2 = sqrt(tau) * exp(-y)
         
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
         
c Momenta of V1
         p(0) = sqrt( pT**2 + pL**2 )
         p(1) = pT * cos(phi)
         p(2) = pT * sin(phi)
         p(3) = pL

c Momenta of V2
         q(0) = sqrt( pT**2 + pL**2 + xq2 )
         q(1) = - pT * cos(phi)
         q(2) = - pT * sin(phi)
         q(3) = - pL

c boost in lab frame
         call boostn( p(0), K, p(0) )
         call boostn( q(0), K, q(0) )

         if (2.0d0*min( p(0)*k1(0)-p(3)*k1(3),
     &        p(0)*k2(0)-p(3)*k2(3) ).lt.m2min) then
            Weight = 0
            return
         endif

         weight = weight * 1d0/tpi/4d0/s
         weight = weight * 3.89379304d11 !conversion into fb
      else
         weight = 0D0
      end if

      
      RETURN
      END


      SUBROUTINE BWBORDERS(M12,MG1,M22,MG2,lower,upper)
      IMPLICIT NONE

      DOUBLE PRECISION M12,MG1,M22,MG2,lower,upper
      DOUBLE PRECISION tmpsqrt, ratio, rootp, rootm

      ratio = MG2/MG1
      tmpsqrt = ratio*(MG1**2+MG2**2+(M12-M22)**2)
      tmpsqrt = sqrt(tmpsqrt)

      rootp = (M22 - ratio*M12 + tmpsqrt)/(1-ratio)
      rootm = (M22 - ratio*M12 - tmpsqrt)/(1-ratio)

      if (rootp .le. rootm) then
        lower = rootp
        upper = rootm
      else
        lower = rootm
        upper = rootp
      endif

      return
      end


      subroutine MultiResonancesInit(MultiResNumber, M2_X_in, MG_X_in, F_X_in, M2_i_in, MG_i_in, F_i_in, N_res_in)
      IMPLICIT NONE

      ! Resonance-Function depending on Q^2, M^2, width of function is M*G

      ! Several "MultiResonances" can be used simultaneously with MultiResNumber as index.
      ! X is broad resonance for continuum
      ! i are the sharp resonances on top
      ! F is the fraction of the random number this resonance should get

      ! assumption: Width_i << Width_X
      ! is checked, resonance i removed when Width_i >= Width_X

#include "VBFNLO/utilities/global.inc"
#include "multires.inc"

      integer          N_res_in, i, j, iMin, MultiResNumber
      logical          resonance_erased
      double precision M2_X_in, MG_X_in, F_X_in, M2_i_in(1:max_res), MG_i_in(1:max_res), F_i_in(1:max_res)
      double precision swap, Fsum, c_low, c_high

      logical ldebug
      parameter (ldebug = .false.)

      if (N_res_in.gt.max_res) then
         print*, "Error: Too many resonances to map onto one random number in subroutine MultiResonancesInit."
         print*, "Please increase parameter max_res!"
         STOP
      endif

      if (MultiResNumber.gt.max_multires) then
         print*, "Error: Too many different calls of MultiResonancesInit."
         print*, "Please increase parameter max_multires!"
         STOP
      endif

      if (N_res_in.lt.0 .or. MultiResNumber.le.0) then
         print*, "Error in parameters of MultiResonancesInit:"
         print*, "N_res_in has to be >= 0 and MultiResNumber has to be > 0 !"
         stop
      endif

      ! remove resonances which are wider than continuum resonance
      j = 0
      do i=1,N_res_in
         if (MG_i_in(i) .lt. MG_X_in) then
           j = j + 1
           M2_i(j, MultiResNumber) = M2_i_in(i)
           MG_i(j, MultiResNumber) = MG_i_in(i)
           F_i(i, MultiResNumber)  = F_i_in(i)
         endif
      enddo
      N_res(MultiResNumber) = j

      ! sort list of resonances ascending by M2
      do i=1,N_res(MultiResNumber)-1
         imin = i
         do j=i+1,N_res(MultiResNumber)
            if (M2_i(j, MultiResNumber) .lt. M2_i(iMin, MultiResNumber)) iMin = j
         enddo
         if (iMin .ne. i) then
            swap = M2_i(i, MultiResNumber)
            M2_i(i, MultiResNumber) = M2_i(iMin, MultiResNumber)
            M2_i(iMin, MultiResNumber) = swap
            swap = MG_i(i, MultiResNumber)
            MG_i(i, MultiResNumber) = MG_i(iMin, MultiResNumber)
            MG_i(iMin, MultiResNumber) = swap
            swap = F_i(i, MultiResNumber)
            F_i(i, MultiResNumber) = F_i(iMin, MultiResNumber)
            F_i(iMin, MultiResNumber) = swap
         endif
      enddo

      ! F_i are normalized, so Sum F_i = 1
      Fsum = F_X_in
      do i=1,N_res(MultiResNumber)
         Fsum = Fsum + F_i_in(i)
      enddo
      F_int(0, MultiResNumber)=0.d0
      do i=1,N_res(MultiResNumber)
         F_i(i, MultiResNumber)   = F_i(i,MultiResNumber) / Fsum
         F_int(i, MultiResNumber) = F_int(i-1, MultiResNumber) + F_i(i, MultiResNumber)
      enddo
      ! F_int is for random number splitting
      ! F_int(i) = sum j=1..i F_i(j)
      ! range F_int(max) .. 1 is for continuum

      F_X(MultiResNumber) = F_X_in / Fsum
      M2_X(MultiResNumber) = M2_X_in
      MG_X(MultiResNumber) = MG_X_in

      ! calculate crossing points for breit-wigner with continuum resonance
      ! if two adjacent resonances are so close that they overlap,
      ! calculate border between them
      do
         resonance_erased = .false.
         if (N_res(MultiResNumber).gt.0) call BWBORDERS(M2_X(MultiResNumber),MG_X(MultiResNumber),M2_i(1, MultiResNumber),
     &                  MG_i(1, MultiResNumber),lower(1, MultiResNumber),upper(1, MultiResNumber))
         do i=1,N_res(MultiResNumber)-1
            call BWBORDERS(M2_X(MultiResNumber),MG_X(MultiResNumber),M2_i(i+1, MultiResNumber),
     &                     MG_i(i+1, MultiResNumber),lower(i+1, MultiResNumber),upper(i+1, MultiResNumber))
            if (upper(i, MultiResNumber) .gt. lower(i+1, MultiResNumber)) then
               call BWBORDERS(M2_i(i, MultiResNumber),MG_i(i, MultiResNumber),M2_i(i+1, MultiResNumber),
     &                        MG_i(i+1, MultiResNumber),c_low,c_high)
               if (M2_i(i, MultiResNumber) .le. c_low .and. c_low .lt. M2_i(i+1, MultiResNumber)) then
                  upper(i, MultiResNumber) = c_low
                  lower(i+1, MultiResNumber) = c_low
               elseif (M2_i(i, MultiResNumber) .le. c_high .and. c_high .lt. M2_i(i+1, MultiResNumber)) then
                  upper(i, MultiResNumber) = c_high
                  lower(i+1, MultiResNumber) = c_high
               else
                  ! if no crossing between both maxima:
                  ! use new resonance at position i
                  ! * with M2 of narrow resonance (evaluate first)
                  if (MG_i(i, MultiResNumber) .gt. MG_i(i+1, MultiResNumber)) 
     &                  M2_i(i, MultiResNumber) = M2_i(i+1, MultiResNumber)
                  ! * with width: geometric mean
                  MG_i(i, MultiResNumber) = sqrt(MG_i(i, MultiResNumber)*MG_i(i+1, MultiResNumber))
                  ! adjust fractions
                  F_i(i, MultiResNumber) = F_i(i, MultiResNumber) + F_i(i+1, MultiResNumber)
                  F_int(i, MultiResNumber) = F_int(i+1, MultiResNumber)
                  ! resonance i+1 gets erased. move everything down and restart calculation
                  do j = i+1,N_res(MultiResNumber)-1
                     M2_i(j, MultiResNumber) = M2_i(j+1, MultiResNumber)
                     MG_i(j, MultiResNumber) = MG_i(j+1, MultiResNumber)
                     F_i(j, MultiResNumber)  = F_i(j+1, MultiResNumber)
                     F_int(j, MultiResNumber)= F_int(j+1, MultiResNumber)
                  enddo
                  resonance_erased = .true.
                  exit
               endif
            endif
         enddo
         if (resonance_erased) then
            N_res(MultiResNumber) = N_res(MultiResNumber) - 1
         else
            exit        ! normal termination
         endif
      enddo

      ! print out resonances
      if (ldebug) then
         print*, ""
         print*, ""
         print*, "MultiResonance number", MultiResNumber
         print*, "**************************************"
         print*, ""
         print*, "Variables M^2, M*Gamma and random number fractions:"
         print*, ""
         print*,    "continuum resonance :  ", M2_X(MultiResNumber), MG_X(MultiResNumber), F_X(MultiResNumber)
         do i = 1, N_res(MultiResNumber)
            print*, "resonance", i, ":", M2_i(i, MultiResNumber), MG_i(i, MultiResNumber), F_i(i, MultiResNumber)
            print*, "borders :              ", max(lower(i, MultiResNumber),0.0), upper(i, MultiResNumber)
         enddo
         print*, ""
         print*, ""
         print*, "Variables M, Gamma and random number fractions:"
         print*, ""
         print*,    "continuum resonance :  ", sqrt(M2_X(MultiResNumber)), 
     &            MG_X(MultiResNumber)/sqrt(M2_X(MultiResNumber)), F_X(MultiResNumber)
         do i = 1, N_res(MultiResNumber)
            print*, "resonance", i, ":", sqrt(M2_i(i, MultiResNumber)), 
     &            MG_i(i, MultiResNumber)/sqrt(M2_i(i, MultiResNumber)), F_i(i, MultiResNumber)
            print*, "borders :              ", max(sqrt(lower(i, MultiResNumber)),0.0), sqrt(upper(i, MultiResNumber))
         enddo
         print*, ""
         print*, ""
      endif

      return
      end


      LOGICAL function MultiResonance(MultiResNumber, rd, rm2minx, rm2maxx, W, QSQ)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "multires.inc"

      integer i, resnum, MultiResNumber
      double precision rd, QSQ, W, rm2minx, rm2maxx, upper_temp, lower_temp
      LOGICAL RESONANCE
      EXTERNAL RESONANCE

      ! note: the borders do not match exactly for 2 reasons:
      ! * restriction to rm2minx, rm2maxx
      ! * continuum resonance gets shifted in ranges of other resonances, so "holes" are "filled".
      !   better than keeping "holes" in the random number
      !   possible improvement: shift up for high mass resonance and down for low-mass resonances,
      !                         keeping M^2 of base resonance unchanged.

      MultiResonance = .false.

      if (rd .gt. F_int(N_res(MultiResNumber), MultiResNumber)) then      ! continuum

         upper_temp = rm2maxx
         do i = 1, N_res(MultiResNumber)                  ! subtract "holes" due to sharp resonances from continuum
            if (upper(i, MultiResNumber).le.rm2maxx .and. upper(i, MultiResNumber).gt.rm2minx) then
              upper_temp = upper_temp - (upper(i, MultiResNumber)-max(lower(i, MultiResNumber),rm2minx))
            elseif (lower(i, MultiResNumber).lt.rm2maxx .and. upper(i, MultiResNumber).gt.rm2minx) then
              upper_temp = upper_temp - (rm2maxx-max(lower(i, MultiResNumber),rm2minx))
            endif
         enddo
         lower_temp = rm2minx
         if (.not. Resonance(M2_X(MultiResNumber), MG_X(MultiResNumber), lower_temp, upper_temp, 
     &                       (rd-F_int(N_res(MultiResNumber), MultiResNumber))/F_X(MultiResNumber), W, QSQ)) return
         W = W / F_X(MultiResNumber)

!         ! shift QSQ back due to "holes"
         do i = 1, N_res(MultiResNumber)
            if (QSQ .gt. lower(i, MultiResNumber) .and. upper(i, MultiResNumber).gt.rm2minx) 
     &          QSQ = QSQ + ( min(upper(i, MultiResNumber),rm2maxx) - max(lower(i, MultiResNumber),rm2minx) )
         enddo

      else                                ! one of the sharp resonances

         do resnum = 1, N_res(MultiResNumber)
            if (F_int(resnum, MultiResNumber).gt.rd) exit
         enddo
         if (lower(resnum, MultiResNumber).lt.rm2minx) then
            lower_temp = rm2minx
         else
            lower_temp = lower(resnum, MultiResNumber)
         endif
         if (upper(resnum, MultiResNumber).gt.rm2maxx) then
            upper_temp = rm2maxx
         else
            upper_temp = upper(resnum, MultiResNumber)
         endif
         if (.not. Resonance(M2_i(resnum, MultiResNumber), MG_i(resnum, MultiResNumber), lower_temp, upper_temp, 
     &                       (rd-F_int(resnum-1, MultiResNumber))/F_i(resnum, MultiResNumber), W, QSQ)) return
         W = W / F_i(resnum, MultiResNumber)

      endif

      MultiResonance = .true.

      return
      end


