***********************************************************************
      FUNCTION H2G_Loop(q1,q2,J1,J2,mass)
*----------------------------------------------------------------------
c   Calculates the Hgg massive quark-loop tensor
c   and contracts it with currents J1,J2.
c   Gluon momenta are outgoing 
c   Coupling factors are not included.
c     
c         J1,q1
c            \
c             *   
c             | \
c          q  |  *-- H
c             | /
c             *
c            / 
c         J2,q2
c     
c   Gluon momenta q1,q2 must contain the invariant mass 
c   in 4th entry.
***********************************************************************

      IMPLICIT NONE
      
      complex*16 H2G_Loop
      
      real*8 q1(0:4),q2(0:4),mass
      real*8 q1q2,mq1,mq2,mh

      complex*16 J1(0:3),J2(0:3)

      complex*16 B0tM,b0tq1,b0tq2,b0tq12,C0,C0_1M

      real*8 lambda
      complex*16 J, L,Delta,A4,F1
      complex*16 TTjj,TLjj,dotjj

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc,B0tM,C0_1M

c   exploit conservation of two currents, i.e.
c     q1.j1 = 0 ,  q2.j2 = 0  etc.

      mq1=q1(4)
      mq2=q2(4)
      q1q2=dotrr(q1,q2)
      mh=2*q1q2+mq1+mq2

c   calculate tensor coefficients

      b0tq1 = B0tM(mass,mq1) 
      b0tq2 = B0tM(mass,mq2)
      b0tq12= B0tM(mass,mh)
      C0 = C0_1M(mass,mq1,mq2,mh)
      
      lambda = q1q2**2-mq1*mq2
      if ( lambda.le.0d0 ) then 
         print*,' singular lambda in H2G_Loop',lambda
         stop
      endif

      L = 2*q1q2*C0 - 2*b0tq12 + b0tq1 + b0tq2 
      Delta = (mq1 - mq2)*(b0tq1 - b0tq2) 
      J = -(4*mass**2-mh) * C0

      F1 = 2*lambda*(2+J) + L*(2*lambda - 3*q1q2*mh) - 3*q1q2*Delta
      F1 = -f1/(8*lambda**2)

      A4 = -2*q1q2*f1 +(mh*L + Delta)/(2*lambda)

      dotjj=dotcc(J1,J2)
      TTjj=q1q2*dotJJ-dotrc(q1,J2)*dotrc(q2,J1)
      TLjj=2*dotJJ * mq1 * mq2

      H2G_Loop=(A4*TTjj + F1*TLjj) 

      RETURN
      END

***********************************************************************
      FUNCTION H2G_effective(q1,q2,J1,J2)
*----------------------------------------------------------------------
c   Calculates the eff. Hgg Vertex for m_top -> infinity
c   and contracts it with two currents J1,J2.
c   Gluon momenta are either both outgoing or both ingoing.
c   Coupling factors are not included.
c     
c           J1,q1
c              \
c               *-- H
c              /
c           J2,q2
***********************************************************************

      IMPLICIT NONE

      complex*16 H2G_effective

      real*8 q1(0:3),q2(0:3)

      complex*16 J1(0:3),J2(0:3)

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc

      H2G_effective=dotrr(q1,q2)*dotcc(J1,J2)-dotrc(q1,J2)*dotrc(q2,J1)

      END

***********************************************************************
       FUNCTION A2G_Loop(q1,q2,J1,J2,mass) 
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 07.06.2006
*---------------------------------------------------------------------- 
c   Calculates the Agg massive quark loop tensor
c   and contracts it with currents J1,J2.
c   Gluon momenta are outgoing 
c   Coupling factors are not included.
     
c         J1,q1
c            \
c             *   
c             | \
c          q  |  *-- A
c             | /
c             *
c            / 
c         J2,q2
     
c   Gluon momenta q1,q2 must contain the invariant mass 
c   in the 4th entry.
***********************************************************************

      IMPLICIT NONE

      complex*16 A2G_Loop
      
      real*8 q1(0:4),q2(0:4),mass
      real*8 mq1,mq2,mA,t2ps(0:3,0:3)

      complex*16 J1(0:3),J2(0:3)

      complex*16 C0_1M

      real*8 dotrr

      EXTERNAL C0_1M,dotrr 
      
c   get the antisymmetric second rank tensor 
      CALL EPSRR (t2ps,q1,q2)
      
C   exploit conservation of the two currents, i.e.
C     q1.j1 = 0 ,  q2.j2 = 0  etc.

      mq1=q1(4)                
      mq2=q2(4)   
      mA=2*dotrr(q1,q2)+mq1+mq2  
     
      A2G_Loop=(-t2ps(0,1)*J1(0)*J2(1) 
     &          -t2ps(0,2)*J1(0)*J2(2)
     &          -t2ps(0,3)*J1(0)*J2(3)
     &          +t2ps(1,2)*J1(1)*J2(2) 
     &	        +t2ps(1,3)*J1(1)*J2(3) 
     &	        +t2ps(2,3)*J1(2)*J2(3) 
     &          -t2ps(1,0)*J1(1)*J2(0) 
     &	        -t2ps(2,0)*J1(2)*J2(0) 
     &          -t2ps(3,0)*J1(3)*J2(0) 
     &          +t2ps(2,1)*J1(2)*J2(1) 
     &          +t2ps(3,1)*J1(3)*J2(1) 
     &          +t2ps(3,2)*J1(3)*J2(2))*C0_1M(mass,mq1,mq2,mA)

      RETURN
      END

***********************************************************************
      FUNCTION A2G_effective(q1,q2,J1,J2)
*----------------------------------------------------------------------
c   Calculates the eff. Agg Vertex for m_top -> infinity
c   and contracts it with the currents J1,J2.
c   Gluon momenta are either both outgoing or both ingoing.
c   Coupling factors are not included.
c     
c           J1,q1
c              \
c               *-- A0
c              /
c           J2,q2
***********************************************************************

      IMPLICIT NONE

      complex*16 A2G_effective

      real*8 q1(0:3),q2(0:3)

      complex*16 J1(0:3),J2(0:3)

      complex*16 Jaux(0:3)

      complex*16 dotcc

      EXTERNAL dotcc
      
      CALL EPSCRR(Jaux,j2,q2,q1) 

      A2G_effective=dotcc(Jaux,J1)      

      END

***********************************************************************
      SUBROUTINE EpsH_Loop(K,virtual,eps,p_H,Kp_H,eps_H,mass)
*----------------------------------------------------------------------
c... modified: 24.08.2011, Michael Kubocz
*----------------------------------------------------------------------
c... Calculates the effective gluon polarization vector for
c... higgs radiation off an external gluon for the cp-even higgs.

c... Input:
c...      K       - the 4-momentum of the gluon
c...      virtual - if K is virtual too
c...      eps     - the COMPLEX!! polarization vector of the gluon
c...      p_H     - the higgs 4-momentum
c ..      mass    - for quark loop

c... Output:
c ..      eps_H  - the eff. polarization vector (complex)
c...      Kp_H   - the momentum K+p_H

c...           massive quark 
c...   K,eps --<--*---*--<-- (K+p_H),eps_H
c ..               \ /
c...                *
c...                |
c...                V
c ..               p_H
***********************************************************************

      IMPLICIT NONE

      real*8 K(0:4),p_H(0:3),Kp_H(0:3)
      real*8 mKp_H,mh,mk,KKp_H,mass,q1q2
      logical virtual

      complex*16 eps(0:3),eps_H(0:3)
      
      complex*16 B0tq1,B0tq2,B0tp,B0tM,C0_1M,C0

      integer mu
      real*8 lambda
      complex*16 J,L,Delta,A4,F1
      
      real*8 dotrr
      complex*16 dotrc
      
      EXTERNAL dotrr,dotrc,B0tM,C0_1M

      mK=0d0
      do mu=0,3
         Kp_H(mu)=(K(mu)+p_H(mu))  
      enddo

      mKp_H= dotrr(Kp_H,Kp_H)
      mH=dotrr(p_H,p_H)
      KKp_H=dotrr(K,Kp_H)

      b0tq1 = B0tM(mass,mk) 
      b0tq2 = B0tM(mass,mKp_H)
      b0tp  = B0tM(mass,mh)
      C0 = C0_1M(mass,mk,mKp_H,mh)

      q1q2 = 0.5d0*(mH-mk-mKp_H)
      lambda = q1q2**2
      if ( lambda.le.0d0 ) then 
         print*,' singular lambda in  EpsH_Loop',lambda
         stop
      endif

      L = 2*q1q2*C0 - 2*b0tp + b0tq1 + b0tq2 
      Delta = (-mKp_H)*(b0tq1 - b0tq2) 
      J = -(4*mass**2-mh) * C0

      F1 = 2*lambda*(2+J) + L*(2*lambda - 3*q1q2*mh) - 3*q1q2*Delta
      F1 = -f1/(8*lambda**2)

      A4 = -2*q1q2*f1 +(mh*L + Delta)/(2*lambda)

      if(virtual) then
         do mu=0,3
            eps_H(mu)=A4*(dotrc(p_H,eps)*K(mu)-eps(mu)
     &               *dotrr(K,Kp_H))/mKp_h+F1*K(4)*eps(mu)
         enddo
      else
         do mu=0,3
            eps_H(mu)=A4*(dotrc(p_H,eps)*K(mu)-eps(mu)
     &               *dotrr(K(0),p_h(0)))/mKp_h
         enddo
      endif

      END

***********************************************************************
      SUBROUTINE EpsH_effective(K,virtual,eps,p_H,K_H,eps_H)
*----------------------------------------------------------------------
c... modified: 24.08.2011, Michael Kubocz
*----------------------------------------------------------------------
c... Calculates the effective gluon polarization vector for higgs
c... (CP-even) radiation off an external gluon in the eff. theory 
c... limit.

c... Input:
c...    K       - the 4-momentum of the gluon
c...    virtual - if K is virtual
c...    eps     - the COMPLEX!! polarization vector of the gluon
c...    p_H     - the higgs 4-momentum
c...
c... Output:
c...    eps_H   - the eff. polarization vector (complex)
c...    K_H     - the momentum K+p_H

c...    K,eps --<--*--<-- (K+p_H),eps_H
c...               |
c...               V
c...              P_H
***********************************************************************

      IMPLICIT NONE

      real*8 K(0:3),p_H(0:3),K_H(0:3)
      logical virtual
      real*8 qsqr,KK_h

      complex*16 eps(0:3),eps_H(0:3)

      integer mu

      real*8 dotrr
      complex*16 dotrc

      EXTERNAL dotrr,dotrc
      
      do mu=0,3
         K_H(mu)=K(mu)+p_H(mu)
      enddo
      if(virtual) then
         KK_h=dotrr(K,K_h)        
      else
         KK_h=dotrr(K,p_h)
      endif
      qsqr=dotrr(k_h,k_h)

      do mu=0,3
         eps_H(mu)=(dotrc(p_h,eps)*K(mu)-eps(mu)*KK_h)/qsqr
      enddo

      END

***********************************************************************
      SUBROUTINE EpsA_Loop(k,eps,p_A,Kp_A,eps_A,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 07.06.2006 
*----------------------------------------------------------------------
c   Calculates the effective gluon polarization vector for
c   higgs radiation off an external gluon for the cp-even higgs.
c   Input:
c        K      - the 4-momentum of the gluon
c        eps    - the COMPLEX!! polarization vector of the gluon
c        p_A    - the higgs 4-momentum
c        mass   - for quark loop

c   Output:
c        eps_A  - the eff. polarization vector (complex)
c        Kp_A   - the momentum K+p_A

c             massive quark 
c     K,eps --<--*---*--<-- (K+p_),eps_A
c                 \ /
c                  *
c                  |
c                  V
c                 P_A
***********************************************************************

      IMPLICIT NONE

      real*8 k(0:3),p_A(0:3),Kp_A(0:3)
      real*8 mKp_A,mA,mk,mass
      
      complex*16 eps(0:3),eps_A(0:3)
      
      complex*16  C0_1M

      integer mu

      real*8 dotrr

      EXTERNAL dotrr,C0_1M
      
      call EPSCRR(eps_A,eps,k,p_A)
		
      do mu=0,3
         Kp_A(mu)=K(mu)+p_A(mu)
      enddo

      mk=0d0
      mKp_A=dotrr(Kp_A,Kp_A)
      mA=dotrr(p_A,p_A)
      
      do mu=0,3
         eps_A(mu)=-eps_A(mu)*C0_1M(mass,mk,mKp_A,mA)/mKp_A   
      enddo      

      END

***********************************************************************
      SUBROUTINE EpsA_effective(K,eps,p_H,K_H,eps_H)
*----------------------------------------------------------------------
c   Calculates the effective gluon polarization vector for
c   higgs radiation off an external gluon for the cp-odd higgs
c   and loop m_quark -> infinity.
c   Input:
c        K      - the 4-momentum of the gluon
c        eps    - the COMPLEX!! polarization vector of the gluon
c        p_H    - the higgs 4-momentum
c
c   Output:
c        eps_H  - the eff. polarization vector (complex)
c        K_H    - the momentum K+p_H
c
c     K,eps --<--*--<-- (K+p_H),eps_H
c                |
c                V
c               P_H  
***********************************************************************

      IMPLICIT NONE

      real*8 K(0:3),p_H(0:3),K_H(0:3)
      real*8 qsqr

      complex*16 eps(0:3),eps_H(0:3)

      integer mu

      real*8 dotrr

      external dotrr
      
      Do mu=0,3
         K_H(mu)=K(mu)+p_H(mu)
      enddo

      qsqr=dotrr(k_h,k_h)

      CALL EPSCRR(eps_H,eps,p_h,K)

      do mu=0,3
         eps_H(mu)= eps_H(mu)/qsqr
      enddo

      END

***********************************************************************
      FUNCTION H3G_Loop(q1,q2,q3,J1,J2,J3,mass)
*----------------------------------------------------------------------
c   Calculates the Hggg Tensor for a massive quark loop
c   and contracts it with the currents J1,J2,J3.
c   Gluon momenta are outgoing.
c   Coupling factors are not included.
c     
c               J1,q1
c                |
c     J2,q2 --*--*
c             |  |
c             *--*-- H
c             |
c           J3,q3
c
c   the gluon momenta q1,q2,q3 must contain the invariant mass 
c   in the 4th entry.
***********************************************************************

      IMPLICIT NONE

#include "d_coefficients.inc"

      complex*16 H3G_Loop

      real*8 q1(0:4),q2(0:4),q3(0:4)
      real*8 q1sq,q2sq,q3sq,q1q2,q1q3,q2q3
      real*8 q12sq, q23sq, q13sq,q123sq,mass

      complex*16 J1(0:3),J2(0:3),J3(0:3)
      complex*16 q1J3,q2J1,q3J2,q2J3,q1J2,q3J1
      complex*16 J12,J13,J23

      complex*16 C0q1q23,C0q3q12,C0q2q13,C0_1M, D_1M_VAL
      integer key123,key231,key312,key213,key132,key321
      integer D_1M_GET

      complex*16 Ba123, Ba231, Ba312, Ba213, Ba132, Ba321,
     &           Bb123, Bb231, Bb312, Bb213, Bb132, Bb321,
     &           Bc123, Bc213

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrc,dotcc, dotrr, C0_1M, D_1M_VAL, D_1M_GET

*     calculate dot products
      q1sq=q1(4)
      q2sq=q2(4)
      q3sq=q3(4)
      q1q2=dotrr(q1,q2)
      q1q3=dotrr(q1,q3)
      q2q3=dotrr(q2,q3)
      q12sq=2*q1q2+q1sq+q2sq
      q23sq=2*q2q3+q2sq+q3sq
      q13sq=2*q1q3+q1sq+q3sq
      q123sq=q13sq+q2sq+2*q1q2+2*q2q3
      q1J3=dotrc(q1,J3)
      q2J3=dotrc(q2,J3)
      q1J2=dotrc(q1,J2)
      q3J2=dotrc(q3,J2)
      q2J1=dotrc(q2,J1)
      q3J1=dotrc(q3,J1)
      J12=dotcc(J1,J2)
      J13=dotcc(J1,J3)
      J23=dotcc(J2,J3)

*     get some C0-functions
      C0q1q23=C0_1M(mass,q1sq,q23sq,q123sq)
      C0q3q12=C0_1M(mass,q3sq,q12sq,q123sq)
      C0q2q13=C0_1M(mass,q2sq,q13sq,q123sq)

*     get keys for Dij-functions (for 6 permutations of momenta)
      key123 = D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
      key321 = D_1M_GET(mass, q3sq, q2sq, q1sq, q2q3, q1q3, q1q2)
      key231 = D_1M_GET(mass, q2sq, q3sq, q1sq, q2q3, q1q2, q1q3)
      key132 = D_1M_GET(mass, q1sq, q3sq, q2sq, q1q3, q1q2, q2q3)
      key312 = D_1M_GET(mass, q3sq, q1sq, q2sq, q1q3, q2q3, q1q2)
      key213 = D_1M_GET(mass, q2sq, q1sq, q3sq, q1q2, q2q3, q1q3)
      
*     get coefficients Ba, Bb, Bc for all permutations
*     Ba :
       Ba123=0.5*q2q3*(D_1M_VAL(dd0 ,key123) + D_1M_VAL(dd0 ,key231) + D_1M_VAL(dd0 ,key312))
     &          -q1q2*(D_1M_VAL(dd13,key231) + D_1M_VAL(dd12,key312) - D_1M_VAL(dd13,key321))
     &          -C0q1q23
     &             -4*(D_1M_VAL(dd313,key231)+ D_1M_VAL(dd312,key312)- D_1M_VAL(dd313,key321))

       Ba321=0.5*q1q2*(D_1M_VAL(dd0 ,key321) + D_1M_VAL(dd0 ,key213) + D_1M_VAL(dd0 ,key132))
     &          -q2q3*(D_1M_VAL(dd13,key213) + D_1M_VAL(dd12,key132) - D_1M_VAL(dd13,key123))
     &          -C0q3q12
     &             -4*(D_1M_VAL(dd313,key213)+ D_1M_VAL(dd312,key132)- D_1M_VAL(dd313,key123))
       
       Ba231=0.5*q1q3*(D_1M_VAL(dd0 ,key231) + D_1M_VAL(dd0 ,key312) + D_1M_VAL(dd0 ,key123))
     &          -q2q3*(D_1M_VAL(dd13,key312) + D_1M_VAL(dd12,key123) - D_1M_VAL(dd13,key132))
     &          -C0q2q13
     &             -4*(D_1M_VAL(dd313,key312)+ D_1M_VAL(dd312,key123)- D_1M_VAL(dd313,key132))
       
       Ba132=0.5*q2q3*(D_1M_VAL(dd0 ,key132) + D_1M_VAL(dd0 ,key321) + D_1M_VAL(dd0 ,key213))
     &          -q1q3*(D_1M_VAL(dd13,key321) + D_1M_VAL(dd12,key213) - D_1M_VAL(dd13,key231))
     &          -C0q1q23
     &             -4*(D_1M_VAL(dd313,key321)+ D_1M_VAL(dd312,key213)- D_1M_VAL(dd313,key231))

       Ba312=0.5*q1q2*(D_1M_VAL(dd0 ,key312) + D_1M_VAL(dd0 ,key123) + D_1M_VAL(dd0 ,key231))
     &          -q1q3*(D_1M_VAL(dd13,key123) + D_1M_VAL(dd12,key231) - D_1M_VAL(dd13,key213))
     &          -C0q3q12
     &             -4*(D_1M_VAL(dd313,key123)+ D_1M_VAL(dd312,key231)- D_1M_VAL(dd313,key213))

       Ba213=0.5*q1q3*(D_1M_VAL(dd0 ,key213) + D_1M_VAL(dd0 ,key132) + D_1M_VAL(dd0 ,key321))
     &          -q1q2*(D_1M_VAL(dd13,key132) + D_1M_VAL(dd12,key321) - D_1M_VAL(dd13,key312))
     &          -C0q2q13
     &             -4*(D_1M_VAL(dd313,key132)+ D_1M_VAL(dd312,key321)- D_1M_VAL(dd313,key312))
*     Bb :
       Bb123=   D_1M_VAL(dd13,key123) + D_1M_VAL(dd12,key231) - D_1M_VAL(dd13,key213)
     &      +4*(D_1M_VAL(dd37,key123) + D_1M_VAL(dd23,key123) + D_1M_VAL(dd38,key231) 
     &         +D_1M_VAL(dd26,key231) - D_1M_VAL(dd39,key213) - D_1M_VAL(dd23,key213) )

       Bb321=   D_1M_VAL(dd13,key321) + D_1M_VAL(dd12,key213) - D_1M_VAL(dd13,key231)
     &      +4*(D_1M_VAL(dd37,key321) + D_1M_VAL(dd23,key321) + D_1M_VAL(dd38,key213) 
     &         +D_1M_VAL(dd26,key213) - D_1M_VAL(dd39,key231) - D_1M_VAL(dd23,key231) )

       Bb231=   D_1M_VAL(dd13,key231) + D_1M_VAL(dd12,key312) - D_1M_VAL(dd13,key321)
     &      +4*(D_1M_VAL(dd37,key231) + D_1M_VAL(dd23,key231) + D_1M_VAL(dd38,key312) 
     &         +D_1M_VAL(dd26,key312) - D_1M_VAL(dd39,key321) - D_1M_VAL(dd23,key321) )

       Bb132=   D_1M_VAL(dd13,key132) + D_1M_VAL(dd12,key321) - D_1M_VAL(dd13,key312)
     &      +4*(D_1M_VAL(dd37,key132) + D_1M_VAL(dd23,key132) + D_1M_VAL(dd38,key321) 
     &         +D_1M_VAL(dd26,key321) - D_1M_VAL(dd39,key312) - D_1M_VAL(dd23,key312) )

       Bb312=   D_1M_VAL(dd13,key312) + D_1M_VAL(dd12,key123) - D_1M_VAL(dd13,key132)
     &      +4*(D_1M_VAL(dd37,key312) + D_1M_VAL(dd23,key312) + D_1M_VAL(dd38,key123) 
     &         +D_1M_VAL(dd26,key123) - D_1M_VAL(dd39,key132) - D_1M_VAL(dd23,key132) )

       Bb213=   D_1M_VAL(dd13,key213) + D_1M_VAL(dd12,key132) - D_1M_VAL(dd13,key123)
     &      +4*(D_1M_VAL(dd37,key213) + D_1M_VAL(dd23,key213) + D_1M_VAL(dd38,key132) 
     &         +D_1M_VAL(dd26,key132) - D_1M_VAL(dd39,key123) - D_1M_VAL(dd23,key123) )
*     Bc
       Bc123= -0.5*(D_1M_VAL(dd0 ,key123) + D_1M_VAL(dd0 ,key231) + D_1M_VAL(dd0 ,key312))
     &         + 4*(D_1M_VAL(dd26,key123) + D_1M_VAL(dd26,key231) + D_1M_VAL(dd26,key312) 
     &             +D_1M_VAL(dd310,key123)+ D_1M_VAL(dd310,key231)+ D_1M_VAL(dd310,key312))

       Bc213= -0.5*(D_1M_VAL(dd0 ,key213) + D_1M_VAL(dd0 ,key132) + D_1M_VAL(dd0 ,key321))
     &         + 4*(D_1M_VAL(dd26,key213) + D_1M_VAL(dd26,key132) + D_1M_VAL(dd26,key321) 
     &             +D_1M_VAL(dd310,key213)+ D_1M_VAL(dd310,key132)+ D_1M_VAL(dd310,key321))

*     calculate the box
      H3G_Loop = J12*q1J3*Ba123 + J23*q2J1*Ba231 + J13*q3J2*Ba312
     &         - J12*q2J3*Ba213 - J13*q1J2*Ba132 - J23*q3J1*Ba321
     &         + q3J1*q3J2*q1J3*Bb123 + q2J1*q1J2*q1J3*Bb231
     &         + q2J1*q3J2*q2J3*Bb312 - q3J1*q3J2*q2J3*Bb213
     &         - q2J1*q1J2*q2J3*Bb132 - q3J1*q1J2*q1J3*Bb321
     &         + q2J1*q3J2*q1J3*Bc123 - q3J1*q1J2*q2J3*Bc213
      
      RETURN 
      END

***********************************************************************
      FUNCTION H3G_effective(q1,q2,q3,J1,J2,J3)
*----------------------------------------------------------------------
c     Calculates the eff. Hggg Vertex for m_top -> infinity
c     and contracts it with the currents J1,J2,J3.
c     Gluon momenta are outgoing.
c     Coupling factors are not included.
c     
c           J1,q1
c             |
c     J2,q2 --*-- H
c             |
c           J3,q3
***********************************************************************

      IMPLICIT NONE

      complex*16 H3G_effective

      real*8 q1(0:3),q2(0:3),q3(0:3)
      real*8 q13(0:3),q32(0:3),q21(0:3)

      complex*16 J1(0:3),J2(0:3),J3(0:3)

      integer mu

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc
      
      do mu=0,3
         q13(mu)=q1(mu)-q3(mu)
         q32(mu)=q3(mu)-q2(mu)
         q21(mu)=q2(mu)-q1(mu)
      enddo

      H3G_effective = dotrc(q13,J2)*dotcc(J1,J3)
     &               +dotrc(q32,J1)*dotcc(J3,J2)
     &               +dotrc(q21,J3)*dotcc(J2,J1)

      END

***********************************************************************
      FUNCTION A3G_Loop(q1,q2,q3,J1,J2,J3,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 07.06.2006
*----------------------------------------------------------------------
c   Calculates the gggA^0 Tensor of the massive quark loop
c   and contracts it with currents J1,J2,J3.
c   Gluon momenta are outgoing.
c   Coupling factors are not included.
     
c               J1,q1
c                |
c     J2,q2 --*--*
c             |  |
c             *--*-- A^0
c             |
c           J3,q3

c   Gluon momenta q1,q2,q3 must contain the invariant mass
c   in the 4th entry.
***********************************************************************

      IMPLICIT NONE

#include "d_coefficients.inc"

      complex*16 A3G_Loop
      
      real*8 q1(0:4),q2(0:4),q3(0:4),mass
      real*8 q1sq, q2sq, q3sq
      real*8 q1q2,q1q3,q2q3
      real*8 q12sq, q13sq, q23sq, q123sq

      complex*16 J1(0:3),J2(0:3),J3(0:3)
      complex*16 q1J1,q1J2,q1J3,q2J1,q2J2,q2J3,q3J1,q3J2,q3J3
      complex*16 J1J2,J1J3,J2J3

      complex*16 C0q3q12,C0q1q23,C0q2q13,C0_1M,D_1M_VAL
      integer D_1M_GET,DijkeyA1,DijkeyA2,DijkeyA3

      complex*16 EJ1J2J3(0:3),
     &           Eq2J1J2(0:3),Eq2J1J3(0:3),Eq2J2J3(0:3),Eq3J1J2(0:3),
     &           Eq3J1J3(0:3),Eq2q3J1(0:3),Eq2q3J2(0:3),Eq2q3J3(0:3),
     &           Eq3J2J3(0:3)
      complex*16 q1EJ1J2J3,q2EJ1J2J3,q3EJ1J2J3,
     &           q1Eq2J1J2,q1Eq2J1J3,q1Eq2J2J3,q1Eq3J1J2,
     &           q1Eq3J1J3,q1Eq3J2J3,q2Eq3J1J2,q2Eq3J1J3,q2Eq3J2J3,
     &           q1Eq2q3J1,q1Eq2q3J2,q1Eq2q3J3

      complex*16 BoxA1,BoxA2,BoxA3

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc,C0_1M,D_1M_VAL,D_1M_GET

c  dot products between gluon momenta
      q1sq=q1(4)
      q2sq=q2(4)
      q3sq=q3(4)
      q1q2=dotrr(q1,q2)
      q1q3=dotrr(q1,q3)
      q2q3=dotrr(q2,q3)
      q12sq=2*q1q2+q1sq+q2sq
      q13sq=2*q1q3+q1sq+q3sq
      q23sq=2*q2q3+q2sq+q3sq
      q123sq=2*q1q2+2*q1q3+2*q2q3+q1sq+q2sq+q3sq

c  dot products between gluon momenta and currents
      q1J1=dotrc(q1,J1)
      q1J2=dotrc(q1,J2)
      q1J3=dotrc(q1,J3)
      q2J1=dotrc(q2,J1)
      q2J2=dotrc(q2,J2)
      q2J3=dotrc(q2,J3)
      q3J1=dotrc(q3,J1)
      q3J2=dotrc(q3,J2)
      q3J3=dotrc(q3,J3)

c dot products between currents
      J1J2=dotcc(J1,J2)
      J1J3=dotcc(J1,J3)
      J2J3=dotcc(J2,J3)

c  contractions of Ji, q2,q3 with the Levi-Civita symbol
      CALL EPSCCC(EJ1J2J3,J1,J2,J3)
      CALL EPSCCR(Eq2J1J2,J2,J1,q2) ! sign
      CALL EPSCCR(Eq2J1J3,J3,J1,q2) ! sign
      CALL EPSCCR(Eq2J2J3,J3,J2,q2) ! sign
      CALL EPSCCR(Eq3J1J2,J2,J1,q3) ! sign
      CALL EPSCCR(Eq3J1J3,J3,J1,q3) ! sign
      CALL EPSCCR(Eq3J2J3,J3,J2,q3) ! sign

      CALL EPSCRR(Eq2q3J1,J1,q3,q2) ! sign
      CALL EPSCRR(Eq2q3J2,J2,q3,q2) ! sign
      CALL EPSCRR(Eq2q3J3,J3,q3,q2) ! sign

c  remaining dot products between the Levi-Civita symbol and qi 
      q1EJ1J2J3=dotrc(q1,EJ1J2J3) 
      q2EJ1J2J3=dotrc(q2,EJ1J2J3)
      q3EJ1J2J3=dotrc(q3,EJ1J2J3)
	
      q1Eq2J1J2=-dotrc(q1,Eq2J1J2) 
      q1Eq2J1J3=-dotrc(q1,Eq2J1J3) 
      q1Eq2J2J3=-dotrc(q1,Eq2J2J3)
      q1Eq3J1J2=-dotrc(q1,Eq3J1J2) 
      q1Eq3J1J3=-dotrc(q1,Eq3J1J3) 
      q1Eq3J2J3=-dotrc(q1,Eq3J2J3)
      q2Eq3J1J2=-dotrc(q2,Eq3J1J2)
      q2Eq3J1J3=-dotrc(q2,Eq3J1J3)
      q2Eq3J2J3=-dotrc(q2,Eq3J2J3)

      q1Eq2q3J1=-dotrc(q1,Eq2q3J1) 
      q1Eq2q3J2=-dotrc(q1,Eq2q3J2)
      q1Eq2q3J3=-dotrc(q1,Eq2q3J3)

c  scalar 3 point functions
      C0q3q12=C0_1M(mass,q3sq,q12sq,q123sq)
      C0q1q23=C0_1M(mass,q1sq,q23sq,q123sq)
      C0q2q13=C0_1M(mass,q2sq,q13sq,q123sq)


c box loop with gluon permutation (123)
      DijkeyA1=D_1M_GET(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3)

      BoxA1=(-1)*
     &  (q3EJ1J2J3*q1q2-q2EJ1J2J3*q1q3-(2*q1Eq3J1J3+q2Eq3J1J3)*q1J2            
     &  +q2Eq3J1J2*q1J3+q1EJ1J2J3*q2q3+q1Eq3J2J3*q2J1     
     &  -q1Eq3J1J2*q2J3-2*q1Eq2J1J2*(q1J3+q2J3)-q1Eq2J2J3*(q3J1)                
     &  +q1Eq2J1J3*q3J2+q1Eq2q3J3*J1J2-q1Eq2q3J2*J1J3                           
     &  +q1Eq2q3J1*J2J3-q2Eq3J2J3*q1J1-2*q1Eq3J1J3*q2J2                         
     &  -q1Eq2J1J2*q3J3)*D_1M_VAL(dd0,DijkeyA1)          
     &  +q3EJ1J2J3*C0q3q12+q1EJ1J2J3*C0q1q23      
     &  +D_1M_VAL(dd11,DijkeyA1)*2
     &  *(q2Eq3J2J3*q1J1+q1Eq3J1J3*q1J2+q1Eq2J1J2*q1J3)
     &  +D_1M_VAL(dd12,DijkeyA1)*2
     &  *(q2Eq3J2J3*q2J1+q1Eq3J1J3*q2J2+q1Eq2J1J2*q2J3)
     &  +D_1M_VAL(dd13,DijkeyA1)*2
     &  *(q2Eq3J2J3*q3J1+q1Eq3J1J3*q3J2+q1Eq2J1J2*q3J3)

c box loop with permutation (231)
      DijkeyA2=D_1M_GET(mass,q2sq,q3sq,q1sq,q2q3,q1q2,q1q3)

      BoxA2=(-1)*
     &  (q1EJ1J2J3*q2q3-q3EJ1J2J3*q1q2-(2*q1Eq2J1J2+q1Eq3J1J2)*q2J3            
     &  -q1Eq3J2J3*q2J1+q2EJ1J2J3*q1q3+q1Eq2J1J3*q3J2            
     &  +q1Eq2J2J3*q3J1-2*q2Eq3J2J3*(q2J1+q3J1)+q2Eq3J1J3*(q1J2)                
     &  -q2Eq3J1J2*q1J3+q1Eq2q3J1*J2J3-q1Eq2q3J3*J1J2                           
     &  +q1Eq2q3J2*J1J3-q1Eq3J1J3*q2J2-2*q1Eq2J1J2*q3J3                         
     &  -q2Eq3J2J3*q1J1)*D_1M_VAL(dd0,DijkeyA2)                                
     &  +q1EJ1J2J3*C0q1q23+q2EJ1J2J3*C0q2q13      
     &  +D_1M_VAL(dd11,DijkeyA2)*2
     &  *(q1Eq3J1J3*q2J2+q1Eq2J1J2*q2J3+q2Eq3J2J3*q2J1)
     &  +D_1M_VAL(dd12,DijkeyA2)*2
     &  *(q1Eq3J1J3*q3J2+q1Eq2J1J2*q3J3+q2Eq3J2J3*q3J1)
     &  +D_1M_VAL(dd13,DijkeyA2)*2
     &  *(q1Eq3J1J3*q1J2+q1Eq2J1J2*q1J3+q2Eq3J2J3*q1J1)

c box loop with permutation (312)
      DijkeyA3=D_1M_GET(mass,q3sq,q1sq,q2sq,q1q3,q2q3,q1q2)

      BoxA3=(-1)*
     &  (q2EJ1J2J3*q1q3-q1EJ1J2J3*q2q3-(2*q2Eq3J2J3-q1Eq2J2J3)*q3J1            
     &  -q1Eq2J1J3*q3J2+q3EJ1J2J3*q1q2-q2Eq3J1J2*q1J3            
     &  -q2Eq3J1J3*q1J2-2*q1Eq3J1J3*(q3J2+q1J2)+q1Eq3J1J2*(q2J3)                
     &  +q1Eq3J2J3*q2J1+q1Eq2q3J2*J1J3-q1Eq2q3J1*J2J3                           
     &  +q1Eq2q3J3*J1J2-q1Eq2J1J2*q3J3-2*q2Eq3J2J3*q1J1                         
     &  -q1Eq3J1J3*q2J2)*D_1M_VAL(dd0,DijkeyA3)                                
     &  +q2EJ1J2J3*C0q2q13+q3EJ1J2J3*C0q3q12      
     &  +D_1M_VAL(dd11,DijkeyA3)*2
     &  *(q1Eq2J1J2*q3J3+q2Eq3J2J3*q3J1+q1Eq3J1J3*q3J2)
     &  +D_1M_VAL(dd12,DijkeyA3)*2
     &  *(q1Eq2J1J2*q1J3+q2Eq3J2J3*q1J1+q1Eq3J1J3*q1J2)
     &  +D_1M_VAL(dd13,DijkeyA3)*2
     &  *(q1Eq2J1J2*q2J3+q2Eq3J2J3*q2J1+q1Eq3J1J3*q2J2)

c  sum of the 3 box loop permutaions
      A3G_Loop=(BoxA1+BoxA2+BoxA3)*0.5d0
c  factor 0.5d0  --> Furry`s Theorem

      END

***********************************************************************
      FUNCTION A3G_effective(p_H,J1,J2,J3)
*----------------------------------------------------------------------
c     Calculates the eff. A0ggg Vertex for m_top -> infinity
c     and contracts it with the currents J1,J2,J3.
c     Gluon momenta are outgoing.
c     Coupling factors are not included.
c     
c             J1
c             |
c        J2 --*-- p_H, A0
c             |
c            J3
***********************************************************************

      IMPLICIT NONE

      complex*16 A3G_effective

      real*8 p_H(0:3)

      complex*16 J1(0:3),J2(0:3),J3(0:3)

      complex*16 Jaux(0:3)

      complex*16 dotrc

      EXTERNAL dotrc

      CALL EPSCCC(Jaux,J2,J1,J3)

      A3G_effective = dotrc(p_H,Jaux)

      END

***********************************************************************
      SUBROUTINE GGG_current(K1,eps1,K2,eps2,K12,Jgg)
*----------------------------------------------------------------------
c     Calculates the Gluon-Gluon current or rather the 
c     effective 2-Gluon polarization vector.
c
c                  Jgg
c                   v
c                   |
c       q1,eps1 -<--*-->- q2,eps2
c     
***********************************************************************

      IMPLICIT NONE

      real*8 K1(0:3),K2(0:3)
      real*8 k12(0:3),qsqr

      complex*16 eps1(0:3),eps2(0:3),Jgg(0:3)
      
      integer mu

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc

      do mu=0,3
         k12(mu)=k1(mu)+k2(mu)
      enddo

      Qsqr=dotrr(K12,K12)

      do mu=0,3
         Jgg(mu)=(-2*eps1(mu)*dotrc(K1,eps2)
     &            +dotcc(eps1,eps2)*(K1(mu)-K2(mu))
     &            +2*eps2(mu)*dotrc(K2,eps1))/qsqr
      enddo

      END

***********************************************************************
      FUNCTION GGG_vertex(q1,q2,q3,J1,J2,J3)
*----------------------------------------------------------------------
c     Calculates the 3 gluon vertex
c     All gluon momenta are outgoing.
c     Coupling factors are not included.
c     
c             J1,q1
c              /
c     J2,q2 --*
c              \
c             J3,q3
c
***********************************************************************

      IMPLICIT NONE

      complex*16 GGG_vertex

      real*8 q1(0:3),q2(0:3),q3(0:3)
      real*8 q13(0:3),q32(0:3),q21(0:3)

      complex*16 J1(0:3),J2(0:3),J3(0:3)

      integer mu

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc
      
      do mu=0,3
         q13(mu)=q1(mu)-q3(mu)
         q32(mu)=q3(mu)-q2(mu)
         q21(mu)=q2(mu)-q1(mu)
      enddo

      GGG_vertex = dotrc(q13,J2)*dotcc(J1,J3)
     &            +dotrc(q32,J1)*dotcc(J3,J2)
     &            +dotrc(q21,J3)*dotcc(J2,J1)

      END
