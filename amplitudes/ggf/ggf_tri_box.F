cfc      FUNCTION H2G_Loop(q1,q2,J1,J2,mass)
cfc      FUNCTION H3G_Loop(q1,q2,q3,J1,J2,J3,mass)
cfc      SUBROUTINE EpsH_Loop(K,virtual,eps,p_H,Kp_H,eps_H,mass)
*     
cfc     Other functions containing pentagons and hexagons are found in:
*      
c                    ggf_Hpen.F
cfc   and      
cfc               ggf_tri_box_pen_hex.F      
*
*      
cfc    Code for GF Hjj slower with the QUAD system. Need to add IF, for quad precission test
cfc    Maybe use a IF CASE to switch QUAD only for Hjjj and not for Hjj or to use other
cfc    limit to control it.
  
cfc     #ifdef WITH_QUAD
cfc     #endif

c     cfc  difference with last beta version before merge VBFNLO3.0
c     cfc   command to check files regress:
c     cfc   git diff a9d9ebca HEAD amplitudes/ggf/gf_higgs_ME.F
c     cfc  HEAD is equal to 43a9c92658ce791e87c12ca0f051fa023dd03550
c     cfc       
c     cfc   	modified:   amplitudes/hjjj/qqhqqgg_mgfast.F
c     cfc	modified:   loops/ggf_loops/loops_aux.F
c     cfc	modified:   amplitudes/ggf/gf_higgs_ME.F
c     cfc	modified:   amplitudes/ggf/ggf_tri_box.F : QUAD can be desactivated
c     cfc    routine for Hjj is much slower now. To be improve or reduce effect of rescue system.
c     cfc    check the toZERO in the tensor_routines.
c     cfc   How was this implemented before the rescue system in QUAD precision     
cfc Replicas for other H particles or fermion loops.
      
cfc      FUNCTION H2G_effective(q1,q2,J1,J2)
cfc      FUNCTION A2G_Loop(q1,q2,J1,J2,mass)
cfc      FUNCTION A2G_effective(q1,q2,J1,J2)
cfc      FUNCTION H1G_SqLoop(q1,q2,J1,J2,mass)
cfc      FUNCTION H2G_SqLoop(q1,q2,J1,J2,mass)
cfc      FUNCTION H3G_effective(q1,q2,q3,J1,J2,J3)
cfc      FUNCTION A3G_Loop(q1,q2,q3,J1,J2,J3,mass)
cfc      FUNCTION A3G_effective(p_H,J1,J2,J3)
cfc      FUNCTION H3G_SqLoop(q1,q2,q3,J1,J2,J3,mass)
cfc      FUNCTION GGG_vertex(q1,q2,q3,J1,J2,J3)
cfc      FUNCTION H4G_effective2(J1,J2,J3,J4)
cfc      SUBROUTINE EpsH_effective(K,virtual,eps,p_H,K_H,eps_H)
cfc      SUBROUTINE EpsA_Loop(k,eps,p_A,Kp_A,eps_A,mass)
cfc      SUBROUTINE EpsA_effective(K,eps,p_H,K_H,eps_H)
cfc      SUBROUTINE Eps_SqBlobLoop(k,eps,p_H,Kp_H,eps_H,mass)
cfc      SUBROUTINE EpsH_SqTriLoop(k,eps,p_H,Kp_H,eps_H,mass)
cfc      SUBROUTINE GGG_current(K1,eps1,K2,eps2,K12,Jgg)
cfc      SUBROUTINE Vertex4g(J1,J2,J3,Vert4G)
      
      
cmodified: Francisco Campanario
c date : 5.11.2012
c add: Variables for Gram determinants and QUAD precision stuff
c Call to Ward-based Routines for Even Higgs and QUAD precision
c     
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

      real*8 p1(0:3),p2(0:3),ph1(0:3)
      integer jtemp
      complex*16 result
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 GRAMCTO0_old
      real*8 GRAMCTO0,GRAMDTO0
      common/GRAMTO/GRAMCTO0,GRAMDTO0
       logical SingularC(2)
       common/singluonver/SingularC
	logical verdbleSin, verQUADSin
	common/versing/verdbleSin, verQUADSin
        integer counter(10),gf_iteration,counterC(10),counterD(10)
	logical singular
        COMMON /DetCount/singular,counter,gf_iteration,counterC,counterD
        logical QUADVERT,dblevert
cfcc   exploit conservation of two currents, i.e.
cfcc     q1.j1 = 0 ,  q2.j2 = 0  etc.
cfc
cfc      mq1=q1(4)
cfc      mq2=q2(4)
cfc      q1q2=dotrr(q1,q2)
cfc      mh=2*q1q2+mq1+mq2
cfc
cfcc   calculate tensor coefficients
cfc
cfc      b0tq1 = B0tM(mass,mq1) 
cfc      b0tq2 = B0tM(mass,mq2)
cfc      b0tq12= B0tM(mass,mh)
cfc      C0 = C0_1M(mass,mq1,mq2,mh)
cfc      
cfc      lambda = q1q2**2-mq1*mq2
cfc      if ( lambda.le.0d0 ) then 
cfc         print*,' singular lambda in H2G_Loop',lambda
cfc         stop
cfc      endif
cfc
cfc      L = 2*q1q2*C0 - 2*b0tq12 + b0tq1 + b0tq2 
cfc      Delta = (mq1 - mq2)*(b0tq1 - b0tq2) 
cfc      J = -(4*mass**2-mh) * C0
cfc
cfc      F1 = 2*lambda*(2+J) + L*(2*lambda - 3*q1q2*mh) - 3*q1q2*Delta
cfc      F1 = -f1/(8*lambda**2)
cfc
cfc      A4 = -2*q1q2*f1 +(mh*L + Delta)/(2*lambda)
cfc
cfc      dotjj=dotcc(J1,J2)
cfc      TTjj=q1q2*dotJJ-dotrc(q1,J2)*dotrc(q2,J1)
cfc      TLjj=2*dotJJ * mq1 * mq2
cfc
cfc      H2G_Loop=(A4*TTjj + F1*TLjj) 
cfc
c
C     
C     
cfc default value GRAMCTO0=1-d6.
cfc if QUAD implemented, to aim for more accuracy the limit where to put Dble to zero is
cfc lower: Control this with #WITH_QUAD statements for the release.

      SingularC(1)=.false.
      SingularC(2)=.false.
      dbleVERT=.true.              

#ifdef WITH_QUAD
      GRAMCTO0_old=GRAMCTO0
      if(GRAMCTO0.lt.1d-4)GRAMCTO0=1d-4
      if(mass.gt.500.or.mass.lt.100) then
         QUADVERT=.true.
         dbleVERT=.false.
      else
         QUADVERT=.false.
      endif
#endif
      
      do jtemp=0,3
         p1(jtemp)=q1(jtemp)
         p2(jtemp)=q2(jtemp)
         ph1(jtemp)=-(p1(jtemp)+p2(jtemp))
      enddo

       if(dble1.and.dblevert) then
      call gluonverHiggsEven(mass,p1,p2,ph1,j1,j2,1,result)
      result=-result/(4d0*mass)
c       print*, "vert", H2G_Loop/result-1d0
      H2G_Loop=result


cfc TEST for NAN
      if(abs(H2G_Loop).ne.abs(H2G_Loop))then
         H2G_Loop=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif
      if(abs(H2G_Loop)+1d0.eq.abs(H2G_Loop))then
         H2G_Loop=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif
      
c     fc  Set to true in tensor reduction subroutine depending of the limit
c     fc set by GRAMCTO0      
      if(SingularC(1)) then
         verdbleSin =.true.
         H2G_Loop = 0d0
       endif

      endif
#ifdef WITH_QUAD      
       If(QUAD.or. SingularC(1).or.QUADVERT) then
 
c       print*, "mas QUAD", mass
c       stop
        call gluonverHiggsEven_QUAD(mass,p1,p2,ph1,j1,j2,1,result)
       result=-result/(4d0*mass)
c       print*, "vert", H2G_Loop/result-1d0
       H2G_Loop=result

cfc    Test for NAN       
      if(abs(H2G_Loop).ne.abs(H2G_Loop))then
         H2G_Loop=0d0
         verQUADSin =.true.
      endif
      if(abs(H2G_Loop)+1d0.eq.abs(H2G_Loop))then
         H2G_Loop=0d0
         verQUADSin =.true.
      endif
cfc
      

       if(Singular.or.SingularC(2)) then
         verQUADSin =.true.
         H2G_Loop=0d0
       endif

       GRAMCTO0= GRAMCTO0_old

       endif
     
#endif

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
ccccccccccccccccccccccccccccccccccccccccccc
      real*8 p1(0:3),p2(0:3),ph1(0:3)
      integer jtemp
      complex*16 result
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
cccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 GRAMCTO0_old
      real*8 GRAMCTO0,GRAMDTO0
      common/GRAMTO/GRAMCTO0,GRAMDTO0
       logical SingularC(2)
       common/singluonver/SingularC
	logical verdbleSin, verQUADSin
	common/versing/verdbleSin, verQUADSin
        integer counter(10),gf_iteration,counterC(10),counterD(10)
	logical singular
        COMMON /DetCount/singular,counter,gf_iteration,counterC,counterD
        logical QUADVERT,dblevert      
cfcc   get the antisymmetric second rank tensor 
cfc      CALL EPSRR (t2ps,q1,q2)
cfc      
cfcC   exploit conservation of the two currents, i.e.
cfcC     q1.j1 = 0 ,  q2.j2 = 0  etc.
cfc
cfc      mq1=q1(4)                
cfc      mq2=q2(4)   
cfc      mA=2*dotrr(q1,q2)+mq1+mq2  
cfc     
cfc      A2G_Loop=(-t2ps(0,1)*J1(0)*J2(1) 
cfc     &          -t2ps(0,2)*J1(0)*J2(2)
cfc     &          -t2ps(0,3)*J1(0)*J2(3)
cfc     &          +t2ps(1,2)*J1(1)*J2(2) 
cfc     &	        +t2ps(1,3)*J1(1)*J2(3) 
cfc     &	        +t2ps(2,3)*J1(2)*J2(3) 
cfc     &          -t2ps(1,0)*J1(1)*J2(0) 
cfc     &	        -t2ps(2,0)*J1(2)*J2(0) 
cfc     &          -t2ps(3,0)*J1(3)*J2(0) 
cfc     &          +t2ps(2,1)*J1(2)*J2(1) 
cfc     &          +t2ps(3,1)*J1(3)*J2(1) 
cfc     &          +t2ps(3,2)*J1(3)*J2(2))*C0_1M(mass,mq1,mq2,mA)

c      print*, "A2G_loop", A2G_Loop


c     fc default value:
c     fc        GRAMCTO0=1-d6
c     fc  Not changed since rank is lower for the CP-odd case      

      SingularC(1)=.false.
      SingularC(2)=.false.

      dbleVERT=.true.
        
      do jtemp=0,3
         p1(jtemp)=q1(jtemp)
         p2(jtemp)=q2(jtemp)
         ph1(jtemp)=-(p1(jtemp)+p2(jtemp))
      enddo


#ifdef WITH_QUAD            
      if(mass.gt.500.or.mass.lt.100) then
      QUADVERT=.true.
      dbleVERT=.false.
      else
      QUADVERT=.false.
       endif
#endif

       if(dble1.and.dblevert) then
       call gluonverHiggsOdd(mass,p1,p2,ph1,j1,j2,1,result)
       result=-result/(4d0*mass)
       A2G_Loop=result

CFC   TEST FOR NAN's       
      if(abs(A2G_Loop).ne.abs(A2G_Loop))then
         A2G_Loop=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif
      if(abs(A2G_Loop)+1d0.eq.abs(A2G_Loop))then
         A2G_Loop=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif

c     fc  Set to true in tensor reduction subroutine depending of the limit
c     fc set by GRAMCTO0      
      if(SingularC(1)) then
         verdbleSin =.true.
         A2G_Loop = 0d0
       endif

      endif

#ifdef WITH_QUAD
      
       If(QUAD.or. SingularC(1).or.QUADVERT) then
          call gluonverHiggsOdd_QUAD(mass,p1,p2,ph1,j1,j2,1,result)
          result=-result/(4d0*mass)
          A2G_Loop=result
      
c TEST FOR NAN'S
      if(abs(A2G_Loop).ne.abs(A2G_Loop))then
         A2G_Loop=0d0
         verQUADSin =.true.
      endif
      if(abs(A2G_Loop)+1d0.eq.abs(A2G_Loop))then
         A2G_Loop=0d0
         verQUADSin =.true.
      endif

       if(Singular.or.SingularC(2)) then
         verQUADSin =.true.
         A2G_Loop=0d0
       endif

       endif

#endif
       
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
      FUNCTION H1G_SqLoop(q1,q2,J1,J2,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 01.06.2007
*---------------------------------------------------------------------- 
c   Calculates the H^0/h^0 gg massive Squark-loop tensor
c   and contracts it with currents J1,J2.
c   Gluon momenta are outgoing 
c   Coupling factors are not included.
c     
c      J1,q1
c        \
c         \   . .
c          *.     .*--- H^0 / h^0
c         /   . .
c        /       
c      J2,q2

c   Gluon momenta q1,q2 must contain the invariant mass 
c   in 4th entry.
***********************************************************************
      
      IMPLICIT NONE
      
      complex*16 H1G_SqLoop

      real*8 q1(0:4),q2(0:4),mass
      complex*16 J1(0:3),J2(0:3)
      complex*16 B0tM,B0tq12

      complex*16 J1J2
      real*8 q1sq,q2sq,q1q2,q12sq
      
      real*8 dotrr
      complex*16 dotcc

      EXTERNAL dotrr,dotcc,B0tM

      q1sq=q1(4)
      q2sq=q2(4)
      q1q2=dotrr(q1,q2)
      q12sq=(2*q1q2+q1sq+q2sq)
      J1J2=dotcc(J1,J2)
      B0tq12=B0tM(mass,q12sq)

      H1G_SqLoop=J1J2*B0tq12

      RETURN
      END

***********************************************************************
      FUNCTION H2G_SqLoop(q1,q2,J1,J2,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 01.06.2007
*---------------------------------------------------------------------- 
c   Calculates the H^0/h^0 gg massive Squark-loop tensor
c   and contracts it with currents J1,J2.
c   Gluon momenta are outgoing 
c   Coupling factors are not included.
c     
c         J1,q1
c            \
c             *.   
c             . .
c         Sq  .  .*-- H^0 / h^0
c             . .
c             *.
c            / 
c         J2,q2
c     
c   Gluon momenta q1,q2 must contain the invariant mass 
c   in 4th entry.
***********************************************************************

      IMPLICIT NONE

#include "c_coefficients.inc"

      complex*16 H2G_SqLoop

      real*8 q1(0:4),q2(0:4),mass
      complex*16 J1(0:3),J2(0:3)

      complex*16 J1J2
      complex*16 q1J1,q1J2,q2J1,q2J2
      real*8 q1sq,q2sq,q1q2,q12sq

      complex*16 C0q1q2,C0_1M,C_1M_VAL
      integer key12
      integer C_1M_GET

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc,C0_1M,C_1M_VAL,C_1M_GET

      q1sq=q1(4)
      q2sq=q2(4)
      q1q2=dotrr(q1,q2)
      q12sq=2*q1q2+q1sq+q2sq
      q1J1=dotrc(q1,J1)
      q1J2=dotrc(q1,J2)
      q2J1=dotrc(q2,J1)
      q2J2=dotrc(q2,J2)
      J1J2=dotcc(J1,J2)
      
c  get C0-function      
      C0q1q2=C0_1M(mass,q1sq,q2sq,q12sq) !!!!!!!!!!!! sign ?????
c  get Ci- and Cij- functions
      key12=C_1M_GET(mass,q1sq,q2sq,q12sq)

      H2G_SqLoop=
     &   -4*J1J2*C_1M_VAL(cc24,key12)+q1J1*(q2J2*(C0q1q2
     &   +2*(C_1M_VAL(cc11,key12)+C_1M_VAL(cc12,key12)
     &   +2*C_1M_VAL(cc23,key12)))+2*q1J2*(C0q1q2+3*C_1M_VAL(cc11,key12)
     &   +2*C_1M_VAL(cc21,key12)))+2*q2J1*(2*q1J2*(C_1M_VAL(cc12,key12)
     &   +C_1M_VAL(cc23,key12))+q2J2*(C_1M_VAL(cc12,key12)
     &   +2*C_1M_VAL(cc22,key12)))

      RETURN
      END

***********************************************************************
      SUBROUTINE EpsH_Loop(K,virtual,eps,p_H,Kp_H,eps_H,mass)
*----------------------------------------------------------------------
c...  modified: 24.08.2011, Michael Kubocz
c..   modified: F. Campanario to include Quad precision and gluon_subroutines
c..      
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
      real*8 p1(0:3),p2(0:3),ph1(0:3)
      integer jtemp
      complex*16 result,v1(0:3),eps_H1(0:3)
      external delta1
      real*8 delta1
      integer comp
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
       real*8 GRAMCTO0_old
c default values from GramDets (GGFLO)
      real*8 GRAMCTO0,GRAMDTO0
      common/GRAMTO/GRAMCTO0,GRAMDTO0
       logical SingularC(2)
       common/singluonver/SingularC
        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
	logical verdbleSin, verQUADSin
	common/versing/verdbleSin, verQUADSin
        integer counter(10),gf_iteration,counterC(10),counterD(10)
	logical singular
        COMMON /DetCount/singular,counter,gf_iteration,counterC,counterD
        logical QUADVERT,dblevert
cfc      mK=0d0 !org
cfcc      mK=dotrr(K,K)
cfc      do mu=0,3
cfc         Kp_H(mu)=(K(mu)+p_H(mu))  
cfc      enddo
cfc
cfc      mKp_H= dotrr(Kp_H,Kp_H)
cfc      mH=dotrr(p_H,p_H)
cfc      KKp_H=dotrr(K,Kp_H)
cfc
cfc      b0tq1 = B0tM(mass,mk) 
cfc      b0tq2 = B0tM(mass,mKp_H)
cfc      b0tp  = B0tM(mass,mh)
cfc      C0 = C0_1M(mass,mk,mKp_H,mh)
cfc
cfc      q1q2 = 0.5d0*(mH-mk-mKp_H)
cfc      lambda = q1q2**2
cfc      if ( lambda.le.0d0 ) then 
cfc         print*,' singular lambda in  EpsH_Loop',lambda
cfc         stop
cfc      endif
cfc
cfc      L = 2*q1q2*C0 - 2*b0tp + b0tq1 + b0tq2 
cfc      Delta = (-mKp_H)*(b0tq1 - b0tq2) 
cfc      J = -(4*mass**2-mh) * C0
cfc
cfc      F1 = 2*lambda*(2+J) + L*(2*lambda - 3*q1q2*mh) - 3*q1q2*Delta
cfc      F1 = -f1/(8*lambda**2)
cfc
cfc      A4 = -2*q1q2*f1 +(mh*L + Delta)/(2*lambda)
cfc
cfc      if(virtual) then
cfc         do mu=0,3
cfc            eps_H(mu)=A4*(dotrc(p_H,eps)*K(mu)-eps(mu)
cfc     &               *dotrr(K,Kp_H))/mKp_h+F1*K(4)*eps(mu)
cfc         enddo
cfc      else
cfc         do mu=0,3
cfc            eps_H(mu)=A4*(dotrc(p_H,eps)*K(mu)-eps(mu)
cfc     &               *dotrr(K(0),p_h(0)))/mKp_h
cfc         enddo
cfc      endif
cfccfc      do jtemp=0,3
cfccfc         print*, "eps_H(temp)",eps_H(jtemp)
cfccfc       enddo

cfc default value GRAMCTO0=1-d6.
cfc if QUAD implemented, to aim for more accuracy the limit where to put Dble to zero is
cfc lower: Control this with #WITH_QUAD statements for the release.

      SingularC(1)=.false.
      SingularC(2)=.false.
      dbleVERT=.true.

#ifdef WITH_QUAD        
        GRAMCTO0_old=GRAMCTO0
        if(GRAMCTO0.lt.1d-4)GRAMCTO0=1d-4
   
      if(mass.gt.500.or.mass.lt.100) then
      QUADVERT=.true.
      dbleVERT=.false.
      else
      QUADVERT=.false.
       endif
#endif

       do jtemp=0,3
          p1(jtemp)=K(jtemp)
         ph1(jtemp)=p_h(jtemp)
         p2(jtemp)=-(p1(jtemp)+p_h(jtemp))
         Kp_H(jtemp)=-p2(jtemp)
       enddo

       mKp_H= dotrr(p2,p2)
c
       comp=1
       do jtemp=0,3

        v1(0)=delta1(jtemp,0)
        v1(1)=-delta1(jtemp,1)
        v1(2)=-delta1(jtemp,2)
        v1(3)=-delta1(jtemp,3)

       if(dble1.and.dblevert) then
       call gluonverHiggsEven(mass,p1,p2,ph1,eps,v1,comp,result)
       eps_H1(jtemp)=-result/(4d0*mass)/mKp_H
       eps_H(jtemp)=eps_H1(jtemp)

cfc TEST FOR NAN's       
      if(abs(eps_H(jtemp)).ne.abs(eps_H(jtemp)))then
         eps_H(jtemp)=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif
      if(abs(eps_H(jtemp))+1d0.eq.abs(eps_H(jtemp)))then
         eps_H(jtemp)=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif

      
c     fc  Set to true in tensor reduction subroutine depending of the limit
c     fc set by GRAMCTO0      
      if(SingularC(1)) then
         verdbleSin =.true.
        eps_H(jtemp)=0d0
       endif

       endif
#ifdef WITH_QUAD
       if(QUAD.or. SingularC(1).or.QUADVERT) then
       call gluonverHiggsEven_QUAD(mass,p1,p2,ph1,eps,v1,comp,result)
       eps_H1(jtemp)=-result/(4d0*mass)/mKp_H
       eps_H(jtemp)=eps_H1(jtemp)

cfc  TEST for NAN
      if(abs(eps_H(jtemp)).ne.abs(eps_H(jtemp)))then
         eps_H(jtemp)=0d0
         verQUADSin =.true.
      endif
      if(abs(eps_H(jtemp))+1d0.eq.abs(eps_H(jtemp)))then
         eps_H(jtemp)=0d0
         verQUADSin =.true.
      endif
    

       if(Singular.or.SingularC(2)) then
         verQUADSin =.true.
         eps_H(jtemp)=0d0
      endif
      
      endif
#endif
      
       if(comp.eq.1) comp=-1
       enddo

#ifdef WITH_QUAD       
       GRAMCTO0= GRAMCTO0_old
#endif
       END

       function delta1(Int_1,Int_2)
       implicit none
       real*8 delta1
       integer Int_1, Int_2
       If (int_1.eq.Int_2) then
          delta1=1d0
       else
          delta1=0d0
        endif
       end function


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

cccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 p1(0:3),p2(0:3),pA1(0:3)
      integer jtemp
      complex*16 result,v1(0:3),eps_A1(0:3)
      external delta1
      real*8 delta1
      integer comp
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
       real*8 GRAMCTO0_old
c default values from GramDets (GGFLO)
      real*8 GRAMCTO0,GRAMDTO0
      common/GRAMTO/GRAMCTO0,GRAMDTO0
       logical SingularC(2)
       common/singluonver/SingularC
        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
	logical verdbleSin, verQUADSin
	common/versing/verdbleSin, verQUADSin
        integer counter(10),gf_iteration,counterC(10),counterD(10)
	logical singular
        COMMON /DetCount/singular,counter,gf_iteration,counterC,counterD
        logical QUADVERT,dblevert
cfc      call EPSCRR(eps_A,eps,k,p_A)
cfc		
cfc      do mu=0,3
cfc         Kp_A(mu)=K(mu)+p_A(mu)
cfc      enddo
cfc
cfc      mk=0d0
cfc      mKp_A=dotrr(Kp_A,Kp_A)
cfc      mA=dotrr(p_A,p_A)
cfc
cfc      do mu=0,3
cfc         eps_A(mu)=-eps_A(mu)*C0_1M(mass,mk,mKp_A,mA)/mKp_A   
cfc      enddo      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                 NEW
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     fc default value:
c     fc        GRAMCTO0=1-d6
cfc  Not modified as in the CP-even because of the lower rank        

      SingularC(1)=.false.
      SingularC(2)=.false.
    
      dbleVERT=.true.

#ifdef WITH_QUAD      
      if(mass.gt.500.or.mass.lt.100) then
      QUADVERT=.true.
      dbleVERT=.false.
      else
      QUADVERT=.false.
       endif
#endif

       do jtemp=0,3
          p1(jtemp)=K(jtemp)
         pA1(jtemp)=p_A(jtemp)
         p2(jtemp)=-(p1(jtemp)+p_A(jtemp))
         Kp_A(jtemp)=-p2(jtemp)
       enddo

       mKp_A= dotrr(p2,p2)
c
       comp=1
       do jtemp=0,3

        v1(0)=delta1(jtemp,0)
        v1(1)=-delta1(jtemp,1)
        v1(2)=-delta1(jtemp,2)
        v1(3)=-delta1(jtemp,3)

       if(dble1.and.dblevert) then
       call gluonverHiggsOdd(mass,p1,p2,pA1,eps,v1,comp,result)
       eps_A1(jtemp)=-result/(4d0*mass)/mKp_A
       eps_A(jtemp)=eps_A1(jtemp)

cfc Check for NAN       
      if(abs(eps_A(jtemp)).ne.abs(eps_A(jtemp)))then
         eps_A(jtemp)=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif
      if(abs(eps_A(jtemp))+1d0.eq.abs(eps_A(jtemp)))then
         eps_A(jtemp)=0d0
         verdbleSin =.true.
         QUADVERT=.true.
      endif

   
      if(SingularC(1)) then
         verdbleSin =.true.
        eps_A(jtemp)=0d0
       endif

       endif

#ifdef WITH_QUAD       
       if(QUAD.or. SingularC(1).or.QUADVERT) then
       call gluonverHiggsOdd_QUAD(mass,p1,p2,pA1,eps,v1,comp,result)
       eps_A1(jtemp)=-result/(4d0*mass)/mKp_A
       eps_A(jtemp)=eps_A1(jtemp)

cFC  TEST FOR NAN'S
       if(abs(eps_A(jtemp)).ne.abs(eps_A(jtemp)))then
         eps_A(jtemp)=0d0
         verQUADSin =.true.
      endif
       if(abs(eps_A(jtemp))+1D0.eq.abs(eps_A(jtemp)))then
         eps_A(jtemp)=0d0
         verQUADSin =.true.
      endif
      

       if(Singular.or.SingularC(2)) then
         verQUADSin =.true.
         eps_A(jtemp)=0d0
       endif

       endif
#endif
       
       if(comp.eq.1) comp=-1
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
      SUBROUTINE Eps_SqBlobLoop(k,eps,p_H,Kp_H,eps_H,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 01.08.07 
*----------------------------------------------------------------------
c   Calculates the effective gluon polarization vector for higgs radiation 
c   off an external gluon for the MSSM cp-even higgs with squark loop
c  
c   Input:
c        K      - the 4-momentum of the gluon
c        eps    - the COMPLEX!! polarization vector of gluon
c        p_s/bH - the higgs 4-momentum
c        mass   - for squark loop

c   Output:
c        eps_H  - the eff. polarization vector (complex)
c        Kp_H   - the momentum K+p_H

c             massive quark 
c     K,eps --<--*O*--<-- (K+p_),eps_s/bH
c                 *
c                 |
c                 V
c               P_s/bH
***********************************************************************
      
      IMPLICIT NONE

      integer mu
      real*8 k(0:3),p_H(0:3),Kp_H(0:3)
      real*8 mKp_H,mH,mk,mass
            
      complex*16 eps(0:3),eps_H(0:3)
      complex*16 B0tM
      
      real*8 dotrr

      EXTERNAL dotrr,B0tM

      do mu=0,3
         Kp_H(mu)=(K(mu)+p_H(mu))
      enddo

      mKp_H=dotrr(Kp_H,Kp_H)
      mH=dotrr(p_H,p_H)

      do mu=0,3
         eps_H(mu)=eps(mu)*B0tM(mass,mH)/mKp_H 
      enddo      

      END

***********************************************************************
      SUBROUTINE EpsH_SqTriLoop(k,eps,p_H,Kp_H,eps_H,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 01.08.07 
*----------------------------------------------------------------------
c   Calculates the effective gluon polarization vector for higgs radiation 
c   off an external gluon for the MSSM cp-even higgs with squark loop
c  
c   Input:
c        K      - the 4-momentum of the gluon
c        eps    - the COMPLEX!! polarization vector of gluon
c        p_s/bH - the higgs 4-momentum
c        mass   - for squark loop

c   Output:
c        eps_H  - the eff. polarization vector (complex)
c        Kp_H   - the momentum K+p_H

c             massive quark 
c     K,eps --<--*O*--<-- (K+p_),eps_s/bH
c                 *
c                 |
c                 V
c               P_s/bH
***********************************************************************
      
      IMPLICIT NONE

#include "c_coefficients.inc"

      integer mu
      real*8 k(0:3),p_H(0:3),Kp_H(0:3)
      real*8 mKp_H,mH,mk,mass
            
      complex*16 eps(0:3),eps_H(0:3),epsK,epsKp_H
      complex*16 C0KKp_H,C0_1M,C_1M_VAL
      
      integer key12
      integer C_1M_GET

      real*8 dotrr
      complex*16 dotrc

      EXTERNAL dotrr,dotrc,C0_1M,C_1M_VAL,C_1M_GET

      do mu=0,3
         Kp_H(mu)=(K(mu)+p_H(mu))
      enddo

      mk=0d0
      mKp_H=dotrr(Kp_H,Kp_H)
      mH=dotrr(p_H,p_H)

      epsK=-dotrc(K,eps)
      epsKp_H=-dotrc(Kp_H,eps)


c  get Ci- and Cij- functions
      key12=C_1M_GET(mass,mk,mKp_H,mH)
c  get C0-function      
c      C0KKp_H=C0_1M(mass,mk,mKp_H,mH)
      C0KKp_H= C_1M_VAL(cc0,key12)

      do mu=0,3
         eps_H(mu)=(-4*C_1M_VAL(cc24,key12)*eps(mu)+4*epsKp_H
     &    *(C_1M_VAL(cc12,key12)+C_1M_VAL(cc23,key12))*K(mu)+2*epsK
     &    *(C0KKp_H+3*C_1M_VAL(cc11,key12)+2*C_1M_VAL(cc21,key12))*K(mu)
     &    +epsK*(C0KKp_H+2*(C_1M_VAL(cc11,key12)+C_1M_VAL(cc12,key12)+2
     &    *C_1M_VAL(cc23,key12)))*Kp_H(mu)+2*epsKp_H
     &    *(C_1M_VAL(cc12,key12)+2*C_1M_VAL(cc22,key12))*Kp_H(mu))/mKp_H
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
      real*8 q1t(0:3),q2t(0:3),q3t(0:3)
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
      real*8 ph(0:3)
      integer int
      complex*16 sumbox123, box1,box2,box3
ccccc
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
	logical hexeva,peneva,boxeva
	common/evaswitch/hexeva,peneva,boxeva

        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
	logical boxdbleSin, boxQUADSin
	common/boxsing/boxdbleSin, boxQUADSin

c       dble1=.true.
cfc*     calculate dot products
cfc      q1sq=q1(4)
cfc      q2sq=q2(4)
cfc      q3sq=q3(4)
cfc      q1q2=dotrr(q1,q2)
cfc      q1q3=dotrr(q1,q3)
cfc      q2q3=dotrr(q2,q3)
cfc      q12sq=2*q1q2+q1sq+q2sq
cfc      q23sq=2*q2q3+q2sq+q3sq
cfc      q13sq=2*q1q3+q1sq+q3sq
cfc      q123sq=q13sq+q2sq+2*q1q2+2*q2q3
cfc      q1J3=dotrc(q1,J3)
cfc      q2J3=dotrc(q2,J3)
cfc      q1J2=dotrc(q1,J2)
cfc      q3J2=dotrc(q3,J2)
cfc      q2J1=dotrc(q2,J1)
cfc      q3J1=dotrc(q3,J1)
cfc      J12=dotcc(J1,J2)
cfc      J13=dotcc(J1,J3)
cfc      J23=dotcc(J2,J3)

      do int=0,3
         q1t(int)=q1(int)
         q2t(int)=q2(int)
         q3t(int)=q3(int)
         ph(int)=-q1(int)-q2(int)-q3(int)
      enddo
      
      if(boxeva) then
         boxSingular(1)=.false.
         boxSingular(2)=.false.
      if(dble1) then
      CALL gluonboxHiggsEven_T(mass,q1t,q2t,q3t,ph,
     &            J1,J2,J3,Box1,1)

      CALL gluonboxHiggsEven_T(mass,q2t,q3t,q1t,ph,
     &            J2,J3,J1,Box2,1)

      CALL gluonboxHiggsEven_T(mass,q3t,q1t,q2t,ph,
     &     J3,J1,J2,Box3,1)
      
      sumbox123= -((Box1+Box2+Box3)/2d0)/(12d0*mass*mass)
      
      if(boxSingular(1)) then
         boxdbleSin =.true.
         sumbox123=0d0   
c     print*, "boxdbleSin0",boxdbleSin
      endif
      endif
#ifdef WITH_QUAD
      if(QUAD.or.boxSingular(1)) then
       CALL gluonboxHiggsEven_QUADT(mass,q1t,q2t,q3t,ph,
     &            J1,J2,J3,Box1,1)

      CALL gluonboxHiggsEven_QUADT(mass,q2t,q3t,q1t,ph,
     &            J2,J3,J1,Box2,1)

      CALL gluonboxHiggsEven_QUADT(mass,q3t,q1t,q2t,ph,
     &     J3,J1,J2,Box3,1)
      
      sumbox123= -((Box1+Box2+Box3)/2d0)/(12d0*mass*mass)
      
       if(boxSingular(2)) then
          boxquadSin =.true.
          sumbox123=0d0   
       endif
      endif
#endif
      else  ! boxeva
       sumbox123=0d0   
      endif
      
      H3G_Loop=sumbox123

c        Print*, "new",sumbox123
c$$$*     get some C0-functions
c$$$      C0q1q23=C0_1M(mass,q1sq,q23sq,q123sq)
c$$$      C0q3q12=C0_1M(mass,q3sq,q12sq,q123sq)
c$$$      C0q2q13=C0_1M(mass,q2sq,q13sq,q123sq)
c$$$
c$$$*     get keys for Dij-functions (for 6 permutations of momenta)
c$$$      key123 = D_1M_GET(mass, q1sq, q2sq, q3sq, q1q2, q1q3, q2q3)
c$$$      key321 = D_1M_GET(mass, q3sq, q2sq, q1sq, q2q3, q1q3, q1q2)
c$$$      key231 = D_1M_GET(mass, q2sq, q3sq, q1sq, q2q3, q1q2, q1q3)
c$$$      key132 = D_1M_GET(mass, q1sq, q3sq, q2sq, q1q3, q1q2, q2q3)
c$$$      key312 = D_1M_GET(mass, q3sq, q1sq, q2sq, q1q3, q2q3, q1q2)
c$$$      key213 = D_1M_GET(mass, q2sq, q1sq, q3sq, q1q2, q2q3, q1q3)
c$$$      
c$$$*     get coefficients Ba, Bb, Bc for all permutations
c$$$*     Ba :
c$$$       Ba123=0.5d0*q2q3*(D_1M_VAL(dd0 ,key123)+D_1M_VAL(dd0 ,key231) 
c$$$     &     + D_1M_VAL(dd0 ,key312))-q1q2*(D_1M_VAL(dd13,key231) 
c$$$     &     + D_1M_VAL(dd12,key312)-D_1M_VAL(dd13,key321))-C0q1q23
c$$$     &     -4*(D_1M_VAL(dd313,key231)+D_1M_VAL(dd312,key312)
c$$$     &     -D_1M_VAL(dd313,key321))
c$$$
c$$$       Ba321=0.5d0*q1q2*(D_1M_VAL(dd0 ,key321)+D_1M_VAL(dd0 ,key213)
c$$$     &      +D_1M_VAL(dd0 ,key132))-q2q3*(D_1M_VAL(dd13,key213)
c$$$     &      +D_1M_VAL(dd12,key132)-D_1M_VAL(dd13,key123))-C0q3q12
c$$$     &      -4*(D_1M_VAL(dd313,key213)+D_1M_VAL(dd312,key132)
c$$$     &      -D_1M_VAL(dd313,key123))
c$$$       
c$$$       Ba231=0.5d0*q1q3*(D_1M_VAL(dd0 ,key231)+D_1M_VAL(dd0 ,key312)
c$$$     &      +D_1M_VAL(dd0 ,key123))-q2q3*(D_1M_VAL(dd13,key312)
c$$$     &      +D_1M_VAL(dd12,key123)-D_1M_VAL(dd13,key132))-C0q2q13
c$$$     &      -4*(D_1M_VAL(dd313,key312)+D_1M_VAL(dd312,key123)
c$$$     &      -D_1M_VAL(dd313,key132))
c$$$       
c$$$       Ba132=0.5d0*q2q3*(D_1M_VAL(dd0 ,key132)+D_1M_VAL(dd0 ,key321)
c$$$     &      +D_1M_VAL(dd0 ,key213))-q1q3*(D_1M_VAL(dd13,key321)
c$$$     &      +D_1M_VAL(dd12,key213)-D_1M_VAL(dd13,key231))-C0q1q23
c$$$     &      -4*(D_1M_VAL(dd313,key321)+D_1M_VAL(dd312,key213)
c$$$     &      -D_1M_VAL(dd313,key231))
c$$$
c$$$       Ba312=0.5d0*q1q2*(D_1M_VAL(dd0 ,key312)+D_1M_VAL(dd0 ,key123)
c$$$     &      +D_1M_VAL(dd0 ,key231))-q1q3*(D_1M_VAL(dd13,key123)
c$$$     &      +D_1M_VAL(dd12,key231)-D_1M_VAL(dd13,key213))-C0q3q12
c$$$     &      -4*(D_1M_VAL(dd313,key123)+D_1M_VAL(dd312,key231)
c$$$     &      -D_1M_VAL(dd313,key213))
c$$$
c$$$       Ba213=0.5d0*q1q3*(D_1M_VAL(dd0 ,key213)+D_1M_VAL(dd0 ,key132)
c$$$     &      +D_1M_VAL(dd0 ,key321))-q1q2*(D_1M_VAL(dd13,key132)
c$$$     &      +D_1M_VAL(dd12,key321)-D_1M_VAL(dd13,key312))-C0q2q13
c$$$     &      -4*(D_1M_VAL(dd313,key132)+D_1M_VAL(dd312,key321)
c$$$     &      -D_1M_VAL(dd313,key312))
c$$$*     Bb :
c$$$       Bb123=D_1M_VAL(dd13,key123)+D_1M_VAL(dd12,key231) 
c$$$     &      -D_1M_VAL(dd13,key213)+4*(D_1M_VAL(dd37,key123)
c$$$     &      +D_1M_VAL(dd23,key123)+D_1M_VAL(dd38,key231) 
c$$$     &      +D_1M_VAL(dd26,key231)-D_1M_VAL(dd39,key213)
c$$$     &      -D_1M_VAL(dd23,key213))
c$$$
c$$$       Bb321=D_1M_VAL(dd13,key321)+D_1M_VAL(dd12,key213)
c$$$     &      -D_1M_VAL(dd13,key231)+4*(D_1M_VAL(dd37,key321)
c$$$     &      +D_1M_VAL(dd23,key321)+D_1M_VAL(dd38,key213) 
c$$$     &      +D_1M_VAL(dd26,key213)-D_1M_VAL(dd39,key231)
c$$$     &      -D_1M_VAL(dd23,key231))
c$$$
c$$$       Bb231=D_1M_VAL(dd13,key231)+D_1M_VAL(dd12,key312)
c$$$     &      -D_1M_VAL(dd13,key321)+4*(D_1M_VAL(dd37,key231)
c$$$     &      +D_1M_VAL(dd23,key231)+D_1M_VAL(dd38,key312)
c$$$     &      +D_1M_VAL(dd26,key312)-D_1M_VAL(dd39,key321)
c$$$     &      -D_1M_VAL(dd23,key321))
c$$$
c$$$       Bb132=D_1M_VAL(dd13,key132)+D_1M_VAL(dd12,key321)
c$$$     &      -D_1M_VAL(dd13,key312)+4*(D_1M_VAL(dd37,key132)
c$$$     &      +D_1M_VAL(dd23,key132)+D_1M_VAL(dd38,key321) 
c$$$     &      +D_1M_VAL(dd26,key321)-D_1M_VAL(dd39,key312)
c$$$     &      -D_1M_VAL(dd23,key312))
c$$$
c$$$       Bb312=D_1M_VAL(dd13,key312)+D_1M_VAL(dd12,key123)
c$$$     &      -D_1M_VAL(dd13,key132)+4*(D_1M_VAL(dd37,key312)
c$$$     &      +D_1M_VAL(dd23,key312)+D_1M_VAL(dd38,key123)
c$$$     &      +D_1M_VAL(dd26,key123)-D_1M_VAL(dd39,key132)
c$$$     &      -D_1M_VAL(dd23,key132))
c$$$
c$$$       Bb213=D_1M_VAL(dd13,key213)+D_1M_VAL(dd12,key132)
c$$$     &      -D_1M_VAL(dd13,key123)+4*(D_1M_VAL(dd37,key213) 
c$$$     &      +D_1M_VAL(dd23,key213)+D_1M_VAL(dd38,key132) 
c$$$     &      +D_1M_VAL(dd26,key132)-D_1M_VAL(dd39,key123)
c$$$     &      -D_1M_VAL(dd23,key123))
c$$$*     Bc
c$$$       Bc123=-0.5d0*(D_1M_VAL(dd0 ,key123)+D_1M_VAL(dd0 ,key231)
c$$$     &      +D_1M_VAL(dd0 ,key312))+4*(D_1M_VAL(dd26,key123)
c$$$     &      +D_1M_VAL(dd26,key231)+D_1M_VAL(dd26,key312) 
c$$$     &      +D_1M_VAL(dd310,key123)+D_1M_VAL(dd310,key231)
c$$$     &      +D_1M_VAL(dd310,key312))
c$$$
c$$$       Bc213=-0.5d0*(D_1M_VAL(dd0 ,key213)+D_1M_VAL(dd0 ,key132)
c$$$     &      +D_1M_VAL(dd0 ,key321))+4*(D_1M_VAL(dd26,key213)
c$$$     &      +D_1M_VAL(dd26,key132)+D_1M_VAL(dd26,key321) 
c$$$     &      +D_1M_VAL(dd310,key213)+D_1M_VAL(dd310,key132)
c$$$     &      +D_1M_VAL(dd310,key321))
c$$$
c$$$*     calculate the box
c$$$      H3G_Loop = J12*q1J3*Ba123 + J23*q2J1*Ba231 + J13*q3J2*Ba312
c$$$     &         - J12*q2J3*Ba213 - J13*q1J2*Ba132 - J23*q3J1*Ba321
c$$$     &         + q3J1*q3J2*q1J3*Bb123 + q2J1*q1J2*q1J3*Bb231
c$$$     &         + q2J1*q3J2*q2J3*Bb312 - q3J1*q3J2*q2J3*Bb213
c$$$     &         - q2J1*q1J2*q2J3*Bb132 - q3J1*q1J2*q1J3*Bb321
c$$$     &         + q2J1*q3J2*q1J3*Bc123 - q3J1*q1J2*q2J3*Bc213


cc start check       
c      print*
c      Print*, "ratio boxes",sumbox123/H3G_Loop-1d0
c      Print*, "new",sumbox123
c      Print*, "old",H3G_Loop
c      print*
cc en check

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
      real*8 q1t(0:3),q2t(0:3),q3t(0:3)
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

      complex*16 sumBox123,BoxA1,BoxA2,BoxA3

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc,C0_1M,D_1M_VAL,D_1M_GET

ccccc
      real*8 ph(0:3)
      integer int
      common/dbleswitch/dble1,quad
      logical dble1, quad
       real*8 gaugelimit 
       common/gauge/gaugelimit
	logical hexeva,peneva,boxeva
	common/evaswitch/hexeva,peneva,boxeva

        logical HexSingular(2),pensingular(2),boxSingular(2),verSingular(2)
        common/hexSin/HexSingular,pensingular,boxSingular,verSingular
	logical boxdbleSin, boxQUADSin
	common/boxsing/boxdbleSin, boxQUADSin



cfcc        Print*, "new_Q:AB",A3G_Loop
cfc
cfcc  dot products between gluon momenta
cfc      q1sq=q1(4)
cfc      q2sq=q2(4)
cfc      q3sq=q3(4)
cfc      q1q2=dotrr(q1,q2)
cfc      q1q3=dotrr(q1,q3)
cfc      q2q3=dotrr(q2,q3)
cfc      q12sq=2*q1q2+q1sq+q2sq
cfc      q13sq=2*q1q3+q1sq+q3sq
cfc      q23sq=2*q2q3+q2sq+q3sq
cfc      q123sq=2*q1q2+2*q1q3+2*q2q3+q1sq+q2sq+q3sq
cfc
cfcc  dot products between gluon momenta and currents
cfc      q1J1=dotrc(q1,J1)
cfc      q1J2=dotrc(q1,J2)
cfc      q1J3=dotrc(q1,J3)
cfc      q2J1=dotrc(q2,J1)
cfc      q2J2=dotrc(q2,J2)
cfc      q2J3=dotrc(q2,J3)
cfc      q3J1=dotrc(q3,J1)
cfc      q3J2=dotrc(q3,J2)
cfc      q3J3=dotrc(q3,J3)
cfc
cfcc dot products between currents
cfc      J1J2=dotcc(J1,J2)
cfc      J1J3=dotcc(J1,J3)
cfc      J2J3=dotcc(J2,J3)
cfc
cfcc  contractions of Ji, q2,q3 with the Levi-Civita symbol
cfc      CALL EPSCCC(EJ1J2J3,J1,J2,J3)
cfc      CALL EPSCCR(Eq2J1J2,J2,J1,q2) ! sign
cfc      CALL EPSCCR(Eq2J1J3,J3,J1,q2) ! sign
cfc      CALL EPSCCR(Eq2J2J3,J3,J2,q2) ! sign
cfc      CALL EPSCCR(Eq3J1J2,J2,J1,q3) ! sign
cfc      CALL EPSCCR(Eq3J1J3,J3,J1,q3) ! sign
cfc      CALL EPSCCR(Eq3J2J3,J3,J2,q3) ! sign
cfc
cfc      CALL EPSCRR(Eq2q3J1,J1,q3,q2) ! sign
cfc      CALL EPSCRR(Eq2q3J2,J2,q3,q2) ! sign
cfc      CALL EPSCRR(Eq2q3J3,J3,q3,q2) ! sign
cfc
cfcc  remaining dot products between the Levi-Civita symbol and qi 
cfc      q1EJ1J2J3=dotrc(q1,EJ1J2J3) 
cfc      q2EJ1J2J3=dotrc(q2,EJ1J2J3)
cfc      q3EJ1J2J3=dotrc(q3,EJ1J2J3)
cfc	
cfc      q1Eq2J1J2=-dotrc(q1,Eq2J1J2) 
cfc      q1Eq2J1J3=-dotrc(q1,Eq2J1J3) 
cfc      q1Eq2J2J3=-dotrc(q1,Eq2J2J3)
cfc      q1Eq3J1J2=-dotrc(q1,Eq3J1J2) 
cfc      q1Eq3J1J3=-dotrc(q1,Eq3J1J3) 
cfc      q1Eq3J2J3=-dotrc(q1,Eq3J2J3)
cfc      q2Eq3J1J2=-dotrc(q2,Eq3J1J2)
cfc      q2Eq3J1J3=-dotrc(q2,Eq3J1J3)
cfc      q2Eq3J2J3=-dotrc(q2,Eq3J2J3)
cfc
cfc      q1Eq2q3J1=-dotrc(q1,Eq2q3J1) 
cfc      q1Eq2q3J2=-dotrc(q1,Eq2q3J2)
cfc      q1Eq2q3J3=-dotrc(q1,Eq2q3J3)
cfc
cfcc  scalar 3 point functions
cfc      C0q3q12=C0_1M(mass,q3sq,q12sq,q123sq)
cfc      C0q1q23=C0_1M(mass,q1sq,q23sq,q123sq)
cfc      C0q2q13=C0_1M(mass,q2sq,q13sq,q123sq)
cfc
cfc
cfcc box loop with gluon permutation (123)
cfc      DijkeyA1=D_1M_GET(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3)
cfc
cfc      BoxA1=(-1)*
cfc     &  (q3EJ1J2J3*q1q2-q2EJ1J2J3*q1q3-(2*q1Eq3J1J3+q2Eq3J1J3)*q1J2            
cfc     &  +q2Eq3J1J2*q1J3+q1EJ1J2J3*q2q3+q1Eq3J2J3*q2J1                            
cfc     &  -q1Eq3J1J2*q2J3-2*q1Eq2J1J2*(q1J3+q2J3)-q1Eq2J2J3*(q3J1)                
cfc     &  +q1Eq2J1J3*q3J2+q1Eq2q3J3*J1J2-q1Eq2q3J2*J1J3                           
cfc     &  +q1Eq2q3J1*J2J3-q2Eq3J2J3*q1J1-2*q1Eq3J1J3*q2J2                         
cfc     &  -q1Eq2J1J2*q3J3)*D_1M_VAL(dd0,DijkeyA1)                                               
cfc     &  +q3EJ1J2J3*C0q3q12+q1EJ1J2J3*C0q1q23      
cfc     &  +D_1M_VAL(dd11,DijkeyA1)*2
cfc     &  *(q2Eq3J2J3*q1J1+q1Eq3J1J3*q1J2+q1Eq2J1J2*q1J3)
cfc     &  +D_1M_VAL(dd12,DijkeyA1)*2
cfc     &  *(q2Eq3J2J3*q2J1+q1Eq3J1J3*q2J2+q1Eq2J1J2*q2J3)
cfc     &  +D_1M_VAL(dd13,DijkeyA1)*2
cfc     &  *(q2Eq3J2J3*q3J1+q1Eq3J1J3*q3J2+q1Eq2J1J2*q3J3)
cfc
cfcc box loop with permutation (231)
cfc      DijkeyA2=D_1M_GET(mass,q2sq,q3sq,q1sq,q2q3,q1q2,q1q3)
cfc
cfc      BoxA2=(-1)*
cfc     &  (q1EJ1J2J3*q2q3-q3EJ1J2J3*q1q2-(2*q1Eq2J1J2+q1Eq3J1J2)*q2J3            
cfc     &  -q1Eq3J2J3*q2J1+q2EJ1J2J3*q1q3+q1Eq2J1J3*q3J2                             
cfc     &  +q1Eq2J2J3*q3J1-2*q2Eq3J2J3*(q2J1+q3J1)+q2Eq3J1J3*(q1J2)                
cfc     &  -q2Eq3J1J2*q1J3+q1Eq2q3J1*J2J3-q1Eq2q3J3*J1J2                           
cfc     &  +q1Eq2q3J2*J1J3-q1Eq3J1J3*q2J2-2*q1Eq2J1J2*q3J3                         
cfc     &  -q2Eq3J2J3*q1J1)*D_1M_VAL(dd0,DijkeyA2)                                
cfc     &  +q1EJ1J2J3*C0q1q23+q2EJ1J2J3*C0q2q13      
cfc     &  +D_1M_VAL(dd11,DijkeyA2)*2
cfc     &  *(q1Eq3J1J3*q2J2+q1Eq2J1J2*q2J3+q2Eq3J2J3*q2J1)
cfc     &  +D_1M_VAL(dd12,DijkeyA2)*2
cfc     &  *(q1Eq3J1J3*q3J2+q1Eq2J1J2*q3J3+q2Eq3J2J3*q3J1)
cfc     &  +D_1M_VAL(dd13,DijkeyA2)*2
cfc     &  *(q1Eq3J1J3*q1J2+q1Eq2J1J2*q1J3+q2Eq3J2J3*q1J1)
cfc
cfcc box loop with permutation (312)
cfc      DijkeyA3=D_1M_GET(mass,q3sq,q1sq,q2sq,q1q3,q2q3,q1q2)
cfc
cfc      BoxA3=(-1)*
cfc     &  (q2EJ1J2J3*q1q3-q1EJ1J2J3*q2q3-(2*q2Eq3J2J3-q1Eq2J2J3)*q3J1            
cfc     &  -q1Eq2J1J3*q3J2+q3EJ1J2J3*q1q2-q2Eq3J1J2*q1J3                             
cfc     &  -q2Eq3J1J3*q1J2-2*q1Eq3J1J3*(q3J2+q1J2)+q1Eq3J1J2*(q2J3)                
cfc     &  +q1Eq3J2J3*q2J1+q1Eq2q3J2*J1J3-q1Eq2q3J1*J2J3                           
cfc     &  +q1Eq2q3J3*J1J2-q1Eq2J1J2*q3J3-2*q2Eq3J2J3*q1J1                         
cfc     &  -q1Eq3J1J3*q2J2)*D_1M_VAL(dd0,DijkeyA3)                                
cfc     &  +q2EJ1J2J3*C0q2q13+q3EJ1J2J3*C0q3q12      
cfc     &  +D_1M_VAL(dd11,DijkeyA3)*2
cfc     &  *(q1Eq2J1J2*q3J3+q2Eq3J2J3*q3J1+q1Eq3J1J3*q3J2)
cfc     &  +D_1M_VAL(dd12,DijkeyA3)*2
cfc     &  *(q1Eq2J1J2*q1J3+q2Eq3J2J3*q1J1+q1Eq3J1J3*q1J2)
cfc     &  +D_1M_VAL(dd13,DijkeyA3)*2
cfc     &  *(q1Eq2J1J2*q2J3+q2Eq3J2J3*q2J1+q1Eq3J1J3*q2J2)
cfc
cfc      Print*, " BoxA1_M",BoxA1
cfc       Print*, "BoxA2_M",BoxA2
cfc       Print*, "BoxA3_M",BoxA3
cfc  sum of the 3 box loop permutaions
cfc      A3G_Loop=(BoxA1+BoxA2+BoxA3)*0.5d0
cfc  factor 0.5d0  --> Furry`s Theorem
cfc

      do int=0,3
         q1t(int)=q1(int)
         q2t(int)=q2(int)
         q3t(int)=q3(int)
         ph(int)=-q1(int)-q2(int)-q3(int)
      enddo

      if(boxeva) then
         boxSingular(1)=.false.
         boxSingular(2)=.false.

         if(dble1) then
      CALL gluonboxHiggsOdd_T(mass,q1t,q2t,q3t,ph,
     &            J1,J2,J3,BoxA1,1)

      CALL gluonboxHiggsOdd_T(mass,q2t,q3t,q1t,ph,
     &            J2,J3,J1,BoxA2,1)

      CALL gluonboxHiggsOdd_T(mass,q3t,q1t,q2t,ph,
     &            J3,J1,J2,BoxA3,1)
      sumbox123= -((BoxA1+BoxA2+BoxA3)/2d0)/(12d0*mass*mass)
      
      if(boxSingular(1)) then
          boxdbleSin =.true.
          sumbox123=0d0   
        endif

      endif
#ifdef WITH_QUAD
      if(QUAD.or.boxSingular(1)) then
       CALL gluonboxHiggsOdd_QUADT(mass,q1t,q2t,q3t,ph,
     &            J1,J2,J3,BoxA1,1)

      CALL gluonboxHiggsOdd_QUADT(mass,q2t,q3t,q1t,ph,
     &            J2,J3,J1,BoxA2,1)

      CALL gluonboxHiggsOdd_QUADT(mass,q3t,q1t,q2t,ph,
     &     J3,J1,J2,BoxA3,1)
      
      sumbox123= -((BoxA1+BoxA2+BoxA3)/2d0)/(12d0*mass*mass)
      
      if(boxSingular(2)) then
          boxquadSin =.true.
          sumbox123=0d0   
       endif
      endif
#endif      
      else   ! boxeva
       sumbox123=0d0   
       endif
       A3G_Loop= sumbox123
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
      FUNCTION H3G_SqLoop(q1,q2,q3,J1,J2,J3,mass)
*----------------------------------------------------------------------
c   author: Michael Kubocz
c   date  : 01.08.2007
*----------------------------------------------------------------------
c   Calculates the gggH^0 / h^0 Tensor of the massive squark loop
c   and contracts it with currents J1,J2,J3.
c   Gluon momenta are outgoing.
c   Coupling factors are not included.
     
c               J1,q1
c                |
c     J2,q2 --*..*
c             :  :
c             *..*-- H^0 / h^0
c             |
c           J3,q3

c   Gluon momenta q1,q2,q3 must contain the invariant mass
c   in the 4th entry.
***********************************************************************

      IMPLICIT NONE

#include "d_coefficients.inc"

      complex*16 H3G_SqLoop

      real*8 q1(0:4),q2(0:4),q3(0:4),mass
      real*8 q1sq, q2sq, q3sq
      real*8 q1q2,q1q3,q2q3
      real*8 q12sq, q13sq, q23sq, q123sq

      complex*16 J1(0:3),J2(0:3),J3(0:3)
      complex*16 q1J1,q1J2,q1J3,q2J1,q2J2,q2J3,q3J1,q3J2,q3J3
      complex*16 J1J2,J1J3,J2J3

      complex*16 D_1M_VAL
      integer D_1M_GET,DijkeyH1,DijkeyH2,DijkeyH3

      complex*16 BoxH1,BoxH2,BoxH3

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc,D_1M_VAL,D_1M_GET

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
      q123sq=q13sq+q2sq+2*q1q2+2*q2q3

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

c box loop with gluon permutation (123)
      DijkeyH1=D_1M_GET(mass,q1sq,q2sq,q3sq,q1q2,q1q3,q2q3)

      BoxH1=(2*(q1J3+q2J3)+q3J3)*q1J1*(2*q1J2+q2J2)
     & *D_1M_VAL(dd0,DijkeyH1)+4*q1J1*q1J2*(q1J3*D_1M_VAL(dd11,DijkeyH1)
     & +q2J3*D_1M_VAL(dd12,DijkeyH1)+q3J3*D_1M_VAL(dd13,DijkeyH1))
     & +2*q1J1*q2J2*(q1J3*D_1M_VAL(dd11,DijkeyH1)+q2J3
     & *D_1M_VAL(dd12,DijkeyH1)+q3J3*D_1M_VAL(dd13,DijkeyH1))+2*(2
     & *(q1J3+q2J3)+q3J3)*(2*q1J2+q2J2)*(q1J1*D_1M_VAL(dd11,DijkeyH1)
     & +q2J1*D_1M_VAL(dd12,DijkeyH1)+q3J1*D_1M_VAL(dd13,DijkeyH1))+2*(2
     & *(q1J3+q2J3)+q3J3)*q1J1*(q1J2*D_1M_VAL(dd11,DijkeyH1)+q2J2
     & *D_1M_VAL(dd12,DijkeyH1)+q3J2  *D_1M_VAL(dd13,DijkeyH1))+8*q1J2
     & *(-(J1J3)*D_1M_VAL(dd27,DijkeyH1)
     & +q1J3*(q1J1*D_1M_VAL(dd21,DijkeyH1)+q2J1*D_1M_VAL(dd24,DijkeyH1)
     & +q3J1*D_1M_VAL(dd25,DijkeyH1))+q3J3*(q3J1*D_1M_VAL(dd23,DijkeyH1)
     & +q1J1*D_1M_VAL(dd25,DijkeyH1)+q2J1*D_1M_VAL(dd26,DijkeyH1))+q2J3
     & *(q2J1*D_1M_VAL(dd22,DijkeyH1)+q1J1*D_1M_VAL(dd24,DijkeyH1)+q3J1
     & *D_1M_VAL(dd26,DijkeyH1)))+4*q2J2*(-(J1J3)
     & *D_1M_VAL(dd27,DijkeyH1)+q1J3*(q1J1*D_1M_VAL(dd21,DijkeyH1)
     & +q2J1*D_1M_VAL(dd24,DijkeyH1)
     & +q3J1*D_1M_VAL(dd25,DijkeyH1))+q3J3*(q3J1*D_1M_VAL(dd23,DijkeyH1)
     & +q1J1*D_1M_VAL(dd25,DijkeyH1)+q2J1*D_1M_VAL(dd26,DijkeyH1))+q2J3
     & *(q2J1*D_1M_VAL(dd22,DijkeyH1)+q1J1*D_1M_VAL(dd24,DijkeyH1)+q3J1
     & *D_1M_VAL(dd26,DijkeyH1)))+4*q1J1*(-(J2J3)
     & *D_1M_VAL(dd27,DijkeyH1)+q1J3*(q1J2*D_1M_VAL(dd21,DijkeyH1)
     & +q2J2*D_1M_VAL(dd24,DijkeyH1)
     & +q3J2*D_1M_VAL(dd25,DijkeyH1))+q3J3*(q3J2*D_1M_VAL(dd23,DijkeyH1)
     & +q1J2*D_1M_VAL(dd25,DijkeyH1)+q2J2*D_1M_VAL(dd26,DijkeyH1))+q2J3
     & *(q2J2*D_1M_VAL(dd22,DijkeyH1)+q1J2*D_1M_VAL(dd24,DijkeyH1)+q3J2
     & *D_1M_VAL(dd26,DijkeyH1)))+4*(2*(q1J3+q2J3)+q3J3)*(-(J1J2)
     & *D_1M_VAL(dd27,DijkeyH1)+q1J1*(q1J2*D_1M_VAL(dd21,DijkeyH1)+q2J2
     & *D_1M_VAL(dd24,DijkeyH1)+q3J2*D_1M_VAL(dd25,DijkeyH1))+q3J1*(q3J2
     & *D_1M_VAL(dd23,DijkeyH1)+q1J2*D_1M_VAL(dd25,DijkeyH1)+q2J2
     & *D_1M_VAL(dd26,DijkeyH1))+q2J1*(q2J2*D_1M_VAL(dd22,DijkeyH1)+q1J2
     & *D_1M_VAL(dd24,DijkeyH1)+q3J2*D_1M_VAL(dd26,DijkeyH1)))+8*(q2J3
     & *q1J2*q3J1*D_1M_VAL(dd310,DijkeyH1)+q2J3*q1J1*q3J2
     & *D_1M_VAL(dd310,DijkeyH1)-J2J3*q1J1*D_1M_VAL(dd311,DijkeyH1)-J1J3
     & *q1J2*D_1M_VAL(dd311,DijkeyH1)-J2J3*q2J1*D_1M_VAL(dd312,DijkeyH1)
     & -J1J3*q2J2*D_1M_VAL(dd312,DijkeyH1)-q2J3*J1J2
     & *D_1M_VAL(dd312,DijkeyH1)-J2J3*q3J1*D_1M_VAL(dd313,DijkeyH1)
     & -J1J3*q3J2*D_1M_VAL(dd313,DijkeyH1)+q2J3*q2J1*q2J2
     & *D_1M_VAL(dd32,DijkeyH1)+q2J3*q1J1*q1J2*D_1M_VAL(dd34,DijkeyH1)
     & +q2J3*q1J2*q2J1*D_1M_VAL(dd36,DijkeyH1)+q2J3*q1J1*q2J2
     & *D_1M_VAL(dd36,DijkeyH1)+q1J3*((q2J2*q3J1+q2J1*q3J2)
     & *D_1M_VAL(dd310,DijkeyH1)-J1J2*D_1M_VAL(dd311,DijkeyH1)+q1J2
     & *(q2J1*D_1M_VAL(dd34,DijkeyH1)+q3J1*D_1M_VAL(dd35,DijkeyH1))
     & +q1J1*(q1J2*D_1M_VAL(dd31,DijkeyH1)+q2J2*D_1M_VAL(dd34,DijkeyH1)
     & +q3J2*D_1M_VAL(dd35,DijkeyH1))+q2J1*q2J2*D_1M_VAL(dd36,DijkeyH1)
     & +q3J1*q3J2*D_1M_VAL(dd37,DijkeyH1))+q2J3*q2J2*q3J1
     & *D_1M_VAL(dd38,DijkeyH1)+q2J3*q2J1*q3J2*D_1M_VAL(dd38,DijkeyH1)
     & +q2J3*q3J1*q3J2*D_1M_VAL(dd39,DijkeyH1)+q3J3*(-(J1J2)
     & *D_1M_VAL(dd313,DijkeyH1)+q3J1*q3J2*D_1M_VAL(dd33,DijkeyH1)
     & +q1J2*(q2J1*D_1M_VAL(dd310,DijkeyH1)+q1J1*D_1M_VAL(dd35,DijkeyH1)
     & +q3J1*D_1M_VAL(dd37,DijkeyH1))+q1J1
     & *(q2J2*D_1M_VAL(dd310,DijkeyH1)+q3J2*D_1M_VAL(dd37,DijkeyH1))
     & +q2J1*q2J2*D_1M_VAL(dd38,DijkeyH1)+(q2J2*q3J1+q2J1*q3J2)
     & *D_1M_VAL(dd39,DijkeyH1)))

c box loop with permutation (231)
      DijkeyH2=D_1M_GET(mass,q2sq,q3sq,q1sq,q2q3,q1q2,q1q3)

      BoxH2=(2*(q2J1+q3J1)+q1J1)*q2J2*(2*q2J3+q3J3)
     & *D_1M_VAL(dd0,DijkeyH2)+4*q2J2*q2J3*(q2J1*D_1M_VAL(dd11,DijkeyH2)
     & +q3J1*D_1M_VAL(dd12,DijkeyH2)+q1J1*D_1M_VAL(dd13,DijkeyH2))
     & +2*q2J2*q3J3*(q2J1*D_1M_VAL(dd11,DijkeyH2)+q3J1
     & *D_1M_VAL(dd12,DijkeyH2)+q1J1*D_1M_VAL(dd13,DijkeyH2))+2*(2*(q2J1
     & +q3J1)+q1J1)*(2*q2J3+q3J3)*(q2J2*D_1M_VAL(dd11,DijkeyH2)+q3J2
     & *D_1M_VAL(dd12,DijkeyH2)+q1J2*D_1M_VAL(dd13,DijkeyH2))+2*(2*(q2J1
     & +q3J1)+q1J1)*q2J2*(q2J3*D_1M_VAL(dd11,DijkeyH2)+q3J3
     & *D_1M_VAL(dd12,DijkeyH2)+q1J3*D_1M_VAL(dd13,DijkeyH2))+8*q2J3
     & *(-(J1J2)*D_1M_VAL(dd27,DijkeyH2)
     & +q2J1*(q2J2*D_1M_VAL(dd21,DijkeyH2)+q3J2*D_1M_VAL(dd24,DijkeyH2)
     & +q1J2*D_1M_VAL(dd25,DijkeyH2))+q1J1*(q1J2*D_1M_VAL(dd23,DijkeyH2)
     & +q2J2*D_1M_VAL(dd25,DijkeyH2)+q3J2*D_1M_VAL(dd26,DijkeyH2))+q3J1
     & *(q3J2*D_1M_VAL(dd22,DijkeyH2)+q2J2*D_1M_VAL(dd24,DijkeyH2)+q1J2
     & *D_1M_VAL(dd26,DijkeyH2)))+4*q3J3*(-(J1J2)
     & *D_1M_VAL(dd27,DijkeyH2)
     & +q2J1*(q2J2*D_1M_VAL(dd21,DijkeyH2)+q3J2*D_1M_VAL(dd24,DijkeyH2)
     & +q1J2*D_1M_VAL(dd25,DijkeyH2))+q1J1*(q1J2*D_1M_VAL(dd23,DijkeyH2)
     & +q2J2*D_1M_VAL(dd25,DijkeyH2)+q3J2*D_1M_VAL(dd26,DijkeyH2))+q3J1
     & *(q3J2*D_1M_VAL(dd22,DijkeyH2)+q2J2*D_1M_VAL(dd24,DijkeyH2)+q1J2
     & *D_1M_VAL(dd26,DijkeyH2)))+4*q2J2*(-(J1J3)
     & *D_1M_VAL(dd27,DijkeyH2)
     & +q2J1*(q2J3*D_1M_VAL(dd21,DijkeyH2)+q3J3*D_1M_VAL(dd24,DijkeyH2)
     & +q1J3*D_1M_VAL(dd25,DijkeyH2))+q1J1*(q1J3*D_1M_VAL(dd23,DijkeyH2)
     & +q2J3*D_1M_VAL(dd25,DijkeyH2)+q3J3*D_1M_VAL(dd26,DijkeyH2))+q3J1
     & *(q3J3*D_1M_VAL(dd22,DijkeyH2)+q2J3*D_1M_VAL(dd24,DijkeyH2)+q1J3
     & *D_1M_VAL(dd26,DijkeyH2)))+4*(2*(q2J1+q3J1)+q1J1)*(-(J2J3)
     & *D_1M_VAL(dd27,DijkeyH2)+q2J2*(q2J3*D_1M_VAL(dd21,DijkeyH2)+q3J3
     & *D_1M_VAL(dd24,DijkeyH2)+q1J3*D_1M_VAL(dd25,DijkeyH2))+q1J2*(q1J3
     & *D_1M_VAL(dd23,DijkeyH2)+q2J3*D_1M_VAL(dd25,DijkeyH2)+q3J3
     & *D_1M_VAL(dd26,DijkeyH2))+q3J2*(q3J3*D_1M_VAL(dd22,DijkeyH2)+q2J3
     & *D_1M_VAL(dd24,DijkeyH2)+q1J3*D_1M_VAL(dd26,DijkeyH2)))+8*(q3J1
     & *q2J3*q1J2*D_1M_VAL(dd310,DijkeyH2)+q3J1*q2J2*q1J3
     & *D_1M_VAL(dd310,DijkeyH2)-J1J3*q2J2*D_1M_VAL(dd311,DijkeyH2)-J1J2
     & *q2J3*D_1M_VAL(dd311,DijkeyH2)-J1J3*q3J2*D_1M_VAL(dd312,DijkeyH2)
     & -J1J2*q3J3*D_1M_VAL(dd312,DijkeyH2)-q3J1*J2J3
     & *D_1M_VAL(dd312,DijkeyH2)-J1J3*q1J2*D_1M_VAL(dd313,DijkeyH2)
     & -J1J2*q1J3*D_1M_VAL(dd313,DijkeyH2)+q3J1*q3J2*q3J3
     & *D_1M_VAL(dd32,DijkeyH2)+q3J1*q2J2*q2J3*D_1M_VAL(dd34,DijkeyH2)
     & +q3J1*q2J3*q3J2*D_1M_VAL(dd36,DijkeyH2)+q3J1*q2J2*q3J3
     & *D_1M_VAL(dd36,DijkeyH2)+q2J1*((q3J3*q1J2+q3J2*q1J3)
     & *D_1M_VAL(dd310,DijkeyH2)-J2J3*D_1M_VAL(dd311,DijkeyH2)+q2J3
     & *(q3J2*D_1M_VAL(dd34,DijkeyH2)+q1J2*D_1M_VAL(dd35,DijkeyH2))+q2J2
     & *(q2J3*D_1M_VAL(dd31,DijkeyH2)+q3J3*D_1M_VAL(dd34,DijkeyH2)+q1J3
     & *D_1M_VAL(dd35,DijkeyH2))+q3J2*q3J3*D_1M_VAL(dd36,DijkeyH2)+q1J2
     & *q1J3*D_1M_VAL(dd37,DijkeyH2))+q3J1*q3J3*q1J2
     & *D_1M_VAL(dd38,DijkeyH2)+q3J1*q3J2*q1J3*D_1M_VAL(dd38,DijkeyH2)
     & +q3J1*q1J2*q1J3*D_1M_VAL(dd39,DijkeyH2)
     & +q1J1*(-(J2J3)*D_1M_VAL(dd313,DijkeyH2)+q1J2*q1J3 
     & *D_1M_VAL(dd33,DijkeyH2)+q2J3*(q3J2*D_1M_VAL(dd310,DijkeyH2)
     & +q2J2*D_1M_VAL(dd35,DijkeyH2)+q1J2*D_1M_VAL(dd37,DijkeyH2))
     & +q2J2*(q3J3*D_1M_VAL(dd310,DijkeyH2)
     & +q1J3*D_1M_VAL(dd37,DijkeyH2))+q3J2*q3J3*D_1M_VAL(dd38,DijkeyH2)
     & +(q3J3*q1J2+q3J2*q1J3)*D_1M_VAL(dd39,DijkeyH2)))

c box loop with permutation (312)
      DijkeyH3=D_1M_GET(mass,q3sq,q1sq,q2sq,q1q3,q2q3,q1q2)

      BoxH3=(2*(q3J2+q1J2)+q2J2)*q3J3*(2*q3J1+q1J1)
     & *D_1M_VAL(dd0,DijkeyH3)+4*q3J3*q3J1*(q3J2*D_1M_VAL(dd11,DijkeyH3)
     & +q1J2*D_1M_VAL(dd12,DijkeyH3)+q2J2*D_1M_VAL(dd13,DijkeyH3))
     & +2*q3J3*q1J1*(q3J2*D_1M_VAL(dd11,DijkeyH3)
     & +q1J2*D_1M_VAL(dd12,DijkeyH3)+q2J2*D_1M_VAL(dd13,DijkeyH3))+2*(2
     & *(q3J2+q1J2)+q2J2)*(2*q3J1+q1J1)*(q3J3*D_1M_VAL(dd11,DijkeyH3)
     & +q1J3*D_1M_VAL(dd12,DijkeyH3)+q2J3*D_1M_VAL(dd13,DijkeyH3)) 
     & +2*(2*(q3J2+q1J2)+q2J2)*q3J3*(q3J1*D_1M_VAL(dd11,DijkeyH3)+q1J1
     & *D_1M_VAL(dd12,DijkeyH3)+q2J1*D_1M_VAL(dd13,DijkeyH3))+8*q3J1
     & *(-(J2J3)*D_1M_VAL(dd27,DijkeyH3)
     & +q3J2*(q3J3*D_1M_VAL(dd21,DijkeyH3)+q1J3*D_1M_VAL(dd24,DijkeyH3)
     & +q2J3*D_1M_VAL(dd25,DijkeyH3))+q2J2*(q2J3*D_1M_VAL(dd23,DijkeyH3)
     & +q3J3*D_1M_VAL(dd25,DijkeyH3)+q1J3*D_1M_VAL(dd26,DijkeyH3))+q1J2
     & *(q1J3*D_1M_VAL(dd22,DijkeyH3)+q3J3*D_1M_VAL(dd24,DijkeyH3)+q2J3
     & *D_1M_VAL(dd26,DijkeyH3)))+4*q1J1*(-(J2J3)
     & *D_1M_VAL(dd27,DijkeyH3)
     & +q3J2*(q3J3*D_1M_VAL(dd21,DijkeyH3)+q1J3*D_1M_VAL(dd24,DijkeyH3)
     & +q2J3*D_1M_VAL(dd25,DijkeyH3))+q2J2*(q2J3*D_1M_VAL(dd23,DijkeyH3)
     & +q3J3*D_1M_VAL(dd25,DijkeyH3)+q1J3*D_1M_VAL(dd26,DijkeyH3))+q1J2
     & *(q1J3*D_1M_VAL(dd22,DijkeyH3)+q3J3*D_1M_VAL(dd24,DijkeyH3)+q2J3
     & *D_1M_VAL(dd26,DijkeyH3)))+4*q3J3*(-(J1J2)
     & *D_1M_VAL(dd27,DijkeyH3)
     & +q3J2*(q3J1*D_1M_VAL(dd21,DijkeyH3)+q1J1*D_1M_VAL(dd24,DijkeyH3)
     & +q2J1*D_1M_VAL(dd25,DijkeyH3))+q2J2*(q2J1*D_1M_VAL(dd23,DijkeyH3)
     & +q3J1*D_1M_VAL(dd25,DijkeyH3)+q1J1*D_1M_VAL(dd26,DijkeyH3))+q1J2
     & *(q1J1*D_1M_VAL(dd22,DijkeyH3)+q3J1*D_1M_VAL(dd24,DijkeyH3)+q2J1
     & *D_1M_VAL(dd26,DijkeyH3)))+4*(2*(q3J2+q1J2)+q2J2)*(-(J1J3)
     & *D_1M_VAL(dd27,DijkeyH3)+q3J3*(q3J1*D_1M_VAL(dd21,DijkeyH3)+q1J1
     & *D_1M_VAL(dd24,DijkeyH3)+q2J1*D_1M_VAL(dd25,DijkeyH3))+q2J3*(q2J1
     & *D_1M_VAL(dd23,DijkeyH3)+q3J1*D_1M_VAL(dd25,DijkeyH3)+q1J1
     & *D_1M_VAL(dd26,DijkeyH3))+q1J3*(q1J1*D_1M_VAL(dd22,DijkeyH3)+q3J1
     & *D_1M_VAL(dd24,DijkeyH3)+q2J1*D_1M_VAL(dd26,DijkeyH3)))+8*(q1J2
     & *q3J1*q2J3*D_1M_VAL(dd310,DijkeyH3)+q1J2*q3J3*q2J1
     & *D_1M_VAL(dd310,DijkeyH3)-J1J2*q3J3*D_1M_VAL(dd311,DijkeyH3)-J2J3
     & *q3J1*D_1M_VAL(dd311,DijkeyH3)-J1J2*q1J3*D_1M_VAL(dd312,DijkeyH3)
     & -J2J3*q1J1*D_1M_VAL(dd312,DijkeyH3)-q1J2*J1J3
     & *D_1M_VAL(dd312,DijkeyH3)-J1J2*q2J3*D_1M_VAL(dd313,DijkeyH3)
     & -J2J3*q2J1*D_1M_VAL(dd313,DijkeyH3)+q1J2*q1J3*q1J1
     & *D_1M_VAL(dd32,DijkeyH3)+q1J2*q3J3*q3J1*D_1M_VAL(dd34,DijkeyH3)
     & +q1J2*q3J1*q1J3*D_1M_VAL(dd36,DijkeyH3)+q1J2*q3J3*q1J1
     & *D_1M_VAL(dd36,DijkeyH3)+q3J2*((q1J1*q2J3+q1J3*q2J1)
     & *D_1M_VAL(dd310,DijkeyH3)-J1J3*D_1M_VAL(dd311,DijkeyH3)+q3J1
     & *(q1J3*D_1M_VAL(dd34,DijkeyH3)+q2J3*D_1M_VAL(dd35,DijkeyH3))
     & +q3J3*(q3J1*D_1M_VAL(dd31,DijkeyH3)+q1J1*D_1M_VAL(dd34,DijkeyH3)
     & +q2J1*D_1M_VAL(dd35,DijkeyH3))+q1J3*q1J1*D_1M_VAL(dd36,DijkeyH3)
     & +q2J3*q2J1*D_1M_VAL(dd37,DijkeyH3))+q1J2*q1J1*q2J3
     & *D_1M_VAL(dd38,DijkeyH3)+q1J2*q1J3*q2J1*D_1M_VAL(dd38,DijkeyH3)
     & +q1J2*q2J3*q2J1*D_1M_VAL(dd39,DijkeyH3)+q2J2*(-(J1J3)
     & *D_1M_VAL(dd313,DijkeyH3)+q2J3*q2J1*D_1M_VAL(dd33,DijkeyH3)
     & +q3J1*(q1J3*D_1M_VAL(dd310,DijkeyH3)+q3J3*D_1M_VAL(dd35,DijkeyH3)
     & +q2J3*D_1M_VAL(dd37,DijkeyH3))+q3J3*(q1J1
     & *D_1M_VAL(dd310,DijkeyH3)+q2J1*D_1M_VAL(dd37,DijkeyH3))+q1J3*q1J1
     & *D_1M_VAL(dd38,DijkeyH3)+(q1J1*q2J3+q1J3*q2J1)
     & *D_1M_VAL(dd39,DijkeyH3)))

c  sum of the 3 box loop permutations
      H3G_SqLoop=-(BoxH1+BoxH2+BoxH3)*0.5d0

c  factor 0.5d0  --> "Furry`s Theorem"

c  the color factor of two diagrams with sq-sq-g-g vertex and opposite 
c  loop momentum ist zero !!!

      RETURN
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

c      do mu=0,3
c         Jgg(mu)=(-2*eps1(mu)*dotrc(K1,eps2)
c     &            +dotcc(eps1,eps2)*(K1(mu)-K2(mu))
c     &            +2*eps2(mu)*dotrc(K2,eps1))/qsqr
c      enddo

      do mu=0,3
         Jgg(mu)=(-2*eps1(mu)*dotrc(K1,eps2)
     &            +dotcc(eps1,eps2)*(K1(mu)-K2(mu))
     &            +2*eps2(mu)*dotrc(K2,eps1))/qsqr
      enddo

cc start check
c      print*,'qsqr',qsqr
cc end check


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

cc start check
c      print*,'dotrc(q13,J2)',dotrc(q13,J2)
c      print*,'dotcc(J1,J3)',dotcc(J1,J3)
c      print*,'dotrc(q32,J1)',dotrc(q32,J1)
c      print*,'dotcc(J3,J2)',dotcc(J3,J2)
c      print*,'dotrc(q21,J3)',dotrc(q21,J3)
c      print*,'dotcc(J2,J1)',dotcc(J2,J1)
c      print*,'GGG_vertex',GGG_vertex
c      print*,''
cc end check      

      END

***********************************************************************
      FUNCTION H4G_effective2(J1,J2,J3,J4)
*----------------------------------------------------------------------
c  Calculates the eff. Hgggg Vertex for m_top -> infinity,
c  i1,i2,i3,i4 gives the permutation of the 4 gluons.
c  Input gluon momenta and polarization vectors are stored via 
c  the needed dot products in q_dot_ep,ep_dot_ep and q_dot_q.
c  Coupling factors are not included.

***********************************************************************

      IMPLICIT NONE

#include "ggf_dotproducts.inc"

      complex*16 H4G_effective2
      complex*16 J1(0:3),J2(0:3),J3(0:3),J4(0:3)
      complex*16 J1J2,J1J3,J1J4,J2J3,J2J4,J3J4
      
      complex*16 dotcc
      
      EXTERNAL dotcc

      J1J2 = dotcc(J1,J2)
      J1J3 = dotcc(J1,J3)
      J1J4 = dotcc(J1,J4)
      J2J3 = dotcc(J2,J3)
      J2J4 = dotcc(J2,J4)
      J3J4 = dotcc(J3,J4)

      H4G_effective2=-(J1J3*J2J4-J1J4*J2J3)
c      H4G_effective2=(J1J2*J3J4-J1J4*J2J3)
c      H4G_effective2=(J1J2*J3J4-J1J3*J2J4)

      RETURN
      END

***********************************************************************
      SUBROUTINE Vertex4g(J1,J2,J3,Vert4G)
*----------------------------------------------------------------------
c  Calculates the eff. Hgggg Vertex for m_top -> infinity,
c  i1,i2,i3,i4 gives the permutation of the 4 gluons.
c  Input gluon momenta and polarization vectors are stored via 
c  the needed dot products in q_dot_ep,ep_dot_ep and q_dot_q.
c  Coupling factors are not included.

***********************************************************************

      IMPLICIT NONE

#include "ggf_dotproducts.inc"

      complex*16 Vert4G(0:3)
      integer mu
      complex*16 J1(0:3),J2(0:3),J3(0:3),J4(0:3)
      complex*16 J1J2,J1J3,J1J4,J2J3,J2J4,J3J4
      
      complex*16 dotcc
      
      EXTERNAL dotcc

      J1J3 = dotcc(J1,J3)
      J2J3 = dotcc(J2,J3)

      do mu=0,3
         Vert4G(mu)=-(J1J3*J2(mu)-J1(mu)*J2J3)
      enddo

      END
