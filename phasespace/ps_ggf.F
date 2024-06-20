** Phasespace generator used for gluon fusion


C
C**********************  LIPSB   ****************************************
C
      SUBROUTINE LIPSB(R,
     1                 K1,K2,Q,P,
     2                 kn,qn,pn)
      IMPLICIT NONE

      REAL*8 PI,TPI
      PARAMETER (PI=3.14159 26535 89793d0,TPI=2D0*PI)

      REAL*8 R,
     1       K1(0:3),K2(0:3), Q(0:4), P(0:3,2),
     2       kn(0:3,2),      qn(0:4),pn(0:3,2)
C
C Phase space generator for production of W/Z,gamma + N jets at hadron 
c colliders
C
C  INPUT:   R         random number (between 0 and 1)
C
C           K1,K2     4-vectors for 2 incoming partons
C           Q(0:4)    4-vector for the vector boson resonance V, the 4th
C                     component contains the mass**2
C           P(*,N)    4-vectors for the N jets, assumed massless
C
C  OUTPUT:  kn, qn, pn  = new 4-vectors with azimuthal angle difference
c                     of pn(*,1) and pn(*,2) determined by R
      REAL*8 PT(2),PHI(2), pcm(0:3)
      INTEGER MU

C set first jet momentum  = p(*,1)
      do mu = 0,3 
         pn(mu,1) = p(mu,1)
      enddo

      PT(2) = dsqrt(p(1,2)**2+p(2,2)**2)
      phi(1) = atan2(p(2,1),p(1,1))

      PHI(2) = PHI(1) + TPI*(R+0.25)
      if (phi(2).gt.pi) phi(2) = phi(2)-tpi
      if (phi(2).gt.pi) phi(2) = phi(2)-tpi

      PN(1,2) = PT(2)*COS(PHI(2))
      PN(2,2) = PT(2)*SIN(PHI(2))
      PN(3,2) = p(3,2)
      PN(0,2) = dsqrt(pn(1,2)**2+pn(2,2)**2+pn(3,2)**2)
      

      QN(1) = - PN(1,1) - PN(1,2)
      QN(2) = - PN(2,1) - PN(2,2)
      QN(3) = Q(3)
      QN(4) = Q(4)
      QN(0) = dSQRT(QN(4)+qn(1)**2+qn(2)**2+qn(3)**2)

      pcm(0) = PN(0,1)+PN(0,2)+QN(0)
      pcm(3) = PN(3,1)+PN(3,2)+QN(3)

      kn(0,1) = (pcm(0) + pcm(3))/2
      kn(3,1) = kn(0,1)
      kn(0,2) = (pcm(0) - pcm(3))/2
      kn(3,2) = -kn(0,2)
      do mu = 1,2
        kn(mu,1) = 0
        kn(mu,2) = 0
      enddo
      END
