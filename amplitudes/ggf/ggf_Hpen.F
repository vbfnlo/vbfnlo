***********************************************************************
      SUBROUTINE H4G_LoopDD(i1,i2,i3,i4,m,pent)

      IMPLICIT NONE
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... date     : 01.08.08
c... modified : 30.01.11
c... author   : Francisco Campanario
c...  date     : 20/09/21
c...  Sustitution to pentagons with the Ward Identities to control instabilities
*----------------------------------------------------------------------
c... This pentagon subroutine uses the Denner-Dittmaier tensor 
c... decomposition.      
c... (i1,i2,i3,i4) gives the permutation of the attached four gluons.
c... Input gluon momenta and polarization vectors are stored via the
c... dot products in q_dot_q, q_dot_ep, ep_dot_ep and explicitly 
c... in qi and eps
c... Normalization: this subroutine gives the pentagon contributions
c... without a factor 4*i*mt, the QCD couplings and Yukawa coupling 
c... h_t 
c... The pentagon is divided into three parts: ggf_Hpen_a,b,c.F

c...                  * -- i1    
c...                /   \
c...         i2 -- *     \
c...               |      * -- A^0
c...         i3 -- *     /
c...                \   /
c...                  * -- i4

***********************************************************************
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
#include "VBFNLO/amplitudes/ggf/d_coefficients.inc"
#include "VBFNLO/amplitudes/ggf/ggf_dotproducts.inc"

      integer i1,i2,i3,i4
      real*8 m
      complex*16 pent

c... Tensor oefficients
      complex*16 C0_1M,C034,
     &           D0_1M,D0234,Dij(3,13),  
     &           E0L,Eij(4,46)

c... local dot products 
      real*8 q1sq,q2sq,q3sq,q4sq,
     &       q1q2,q1q3,q1q4,q2q3,q2q4,q3q4

      complex*16 q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &           q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
     &           J1J2,J1J3,J1J4,J2J3,J2J4,J3J4

c... parts of the pentagon
      complex*16 a,b,c	

c... other variables
      integer i,j


      integer init,ie0(4,4,4,4)
      data init/0/, ie0/256*0/
      save ie0,init

      complex*16 H4G_LoopDD2a,H4G_LoopDD2b,H4G_LoopDD2c
      real*8 dotrr

      EXTERNAL dotrr,C0_1M,D0_1M, H4G_LoopDD2a,H4G_LoopDD2b,
     &	       H4G_LoopDD2c


cfc To replace for pengluon routines where Ward Identities are used to control the instabilities
      integer mu
      real*8 Kpen(0:4,4), PH(0:4)
      complex*16 epsgPen(0:3,4), pen1
      Common/Forpent/Kpen,epsgPen

      Common/ForpentCheck/count,countBAD
      integer count,countbad

      complex*16 epvec(0:3,4),J_ij(0:3,4,4) 
      real*8 K(0:4,5),K_ij(0:4,4,4),K_H(0:3,4)
      complex*16 eps_effH(0:3,4),eps_effA(0:3,4),
     &     eps_tH(0:3,4),eps_bH(0:3,4),eps_tA(0:3,4),eps_bA(0:3,4),
     &     eps_st1bH(0:3,4),eps_st1tH(0:3,4),
     &     eps_st2bH(0:3,4),eps_st2tH(0:3,4),
     &     eps_sb1bH(0:3,4),eps_sb1tH(0:3,4),
     &     eps_sb2bH(0:3,4),eps_sb2tH(0:3,4)
      COMMON /AMPINFO/epvec,J_ij,
     &       eps_effH,eps_effA,eps_tH,eps_bH,eps_tA,eps_bA,
     &       eps_st1bH,eps_st1tH,eps_st2bH,eps_st2tH,
     &       eps_sb1bH,eps_sb1tH,eps_sb2bH,eps_sb2tH,K,K_ij,K_H
ccccccccccccccccccccccccccccccccccccccccccccccc
       logical QUAD,Dble1
       integer jtemp
 
       real*8 gaugelimit 
       common/gauge/gaugelimit
       common/dbleswitch/dble1,quad
cfc
cfc      if (init.eq.0) then   
cfc         init = 1
cfc         ie0(4,2,1,3) = 1   ! this is the map from i1,i2,i3,i4 
cfc         ie0(2,1,3,4) = 2   ! to the index of the E0 variables
cfc         ie0(3,2,1,4) = 3
cfc         ie0(2,1,4,3) = 4
cfc         ie0(4,3,2,1) = 5
cfc         ie0(3,4,2,1) = 6
cfc         ie0(2,3,1,4) = 7
cfc         ie0(3,1,4,2) = 8
cfc         ie0(4,2,3,1) = 9
cfc         ie0(2,4,3,1) = 10
cfc         ie0(3,2,4,1) = 11
cfc         ie0(2,3,4,1) = 12
cfc      endif
cfc
cfc      if (m.eq.m_top) then
cfc         E0L=-E0t(ie0(i1,i2,i3,i4))
cfc         do i=1,4
cfc           do j=1,46 
cfc               Eij(i,j)=DDEijtop(i,j,ie0(i1,i2,i3,i4))
cfc            enddo
cfc         enddo  
cfc      endif
cfc      if (m.eq.m_bottom) then
cfc         E0L =-E0b(ie0(i1,i2,i3,i4))
cfc         do i=1,4
cfc            do j=1,46 
cfc               Eij(i,j)=DDEijbottom(i,j,ie0(i1,i2,i3,i4))
cfc            enddo
cfc         enddo   
cfc      endif
cfc
cfc      q1sq=q_dot_q(i1,i1)
cfc      q2sq=q_dot_q(i2,i2)
cfc      q3sq=q_dot_q(i3,i3)
cfc      q4sq=q_dot_q(i4,i4)
cfc
cfc      q1q2 = q_dot_q(i1,i2)
cfc      q1q3 = q_dot_q(i1,i3)
cfc      q1q4 = q_dot_q(i1,i4)
cfc      q2q3 = q_dot_q(i2,i3)
cfc      q2q4 = q_dot_q(i2,i4)
cfc      q3q4 = q_dot_q(i3,i4)
cfc
cfc      q1J1 = q_dot_ep(i1,i1)
cfc      q2J1 = q_dot_ep(i2,i1)
cfc      q3J1 = q_dot_ep(i3,i1)
cfc      q4J1 = q_dot_ep(i4,i1)
cfc
cfc      q1J2 = q_dot_ep(i1,i2)
cfc      q2J2 = q_dot_ep(i2,i2)
cfc      q3J2 = q_dot_ep(i3,i2)
cfc      q4J2 = q_dot_ep(i4,i2)
cfc
cfc      q1J3 = q_dot_ep(i1,i3)
cfc      q2J3 = q_dot_ep(i2,i3)
cfc      q3J3 = q_dot_ep(i3,i3)
cfc      q4J3 = q_dot_ep(i4,i3)
cfc
cfc      q1J4 = q_dot_ep(i1,i4)
cfc      q2J4 = q_dot_ep(i2,i4)
cfc      q3J4 = q_dot_ep(i3,i4)
cfc      q4J4 = q_dot_ep(i4,i4)
cfc
cfc      J1J2 = ep_dot_ep(i1,i2)
cfc      J1J3 = ep_dot_ep(i1,i3)
cfc      J1J4 = ep_dot_ep(i1,i4)
cfc      J2J3 = ep_dot_ep(i2,i3)
cfc      J2J4 = ep_dot_ep(i2,i4)
cfc      J3J4 = ep_dot_ep(i3,i4)
cfc
cfc      C034=-C0_1M(m,q3sq,q4sq,q3sq+q4sq+2d0*q3q4) 
cfc
cfc      D0234=D0_1M(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4)
cfc      CALL D_1M_ALL2(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4,Dij)
cfc      Dij(2,7) =-Dij(2,7)
cfc      Dij(3,11)=-Dij(3,11)
cfc      Dij(3,12)=-Dij(3,12)
cfc      Dij(3,13)=-Dij(3,13)
cfc
cfc      a = H4G_LoopDD2a(q1sq,q2sq,q3sq,q4sq,
cfc     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
cfc     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
cfc     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
cfc     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
cfc     &                 C034,D0234,Dij,E0L,Eij,m)
cfc
cfc      b = H4G_LoopDD2b(q1sq,q2sq,q3sq,q4sq,
cfc     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
cfc     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
cfc     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
cfc     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
cfc     &                 C034,D0234,Dij,E0L,Eij,m)
cfc
cfc      c = H4G_LoopDD2c(q1sq,q2sq,q3sq,q4sq,
cfc     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
cfc     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
cfc     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
cfc     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
cfc     &                 C034,D0234,Dij,E0L,Eij,m)
cfc
cfc      pent=1d0/(4d0*m)*(a+b+c)
cfccfc TO DO exchange to the new pentagons routines
       pensingular(2)=.false.
       pensingular(1)=.false.
 
       do mu=0,4
         Kpen(mu,i1)=K(mu,i1)
         Kpen(mu,i2)=K(mu,i2)
         Kpen(mu,i3)=K(mu,i3)
         Kpen(mu,i4)=K(mu,i4)
         PH(mu)=-(K(mu,i1)+K(mu,i2)+K(mu,i3)+K(mu,i4))
         if(mu.eq.4) exit
         epsgPen(mu,i1)=epvec(mu,i1)
         epsgPen(mu,i2)=epvec(mu,i2)
         epsgPen(mu,i3)=epvec(mu,i3) 
         epsgPen(mu,i4)=epvec(mu,i4) 
      enddo
        
       
        CALL gluonpenHiggsEven_T(m,Kpen(0,i1),Kpen(0,i2),
     &        Kpen(0,i3),Kpen(0,i4),PH,epsgPen(0,i1),epsgPen(0,i2),
     &     epsgPen(0,i3),epsgPen(0,i4),pen1)
        
        pen1=-pen1/(12d0*m*m)
        
        if (pensingular(1)) then
           PendbleSin=.true.
           pen1=0d0
        endif
        
#ifdef WITH_QUAD
        if(QUAD.or.penSingular(1)) then
           CALL gluonpenHiggsEven_QUADT(m,Kpen(0,i1),Kpen(0,i2),
     &          Kpen(0,i3),Kpen(0,i4),PH,epsgPen(0,i1),epsgPen(0,i2),
     &        epsgPen(0,i3),epsgPen(0,i4),pen1)

           pen1=-pen1/(12d0*m*m)
           
         if (pensingular(2)) then
            PenquadSin=.true.
            pen1=0d0
         endif
         
      endif
#endif
      
cfc       if ( abs(pen1/pent-1d0).gt.1d-2) Then
cfc           countbad=countbad+1
c           print*, i1,i2,i3,i4
cfc           print*, "pen1D",abs(pen1/pent-1d0), pen1, pent, abs((countbad*1d0)/(count*1d0))*100, "%", count, countBad
cfc           print*, "gaugelimit", gaugelimit, pensingular(1)
cfc        else
cfc        count=count+1   
cfc        endif
        pent=pen1
       
      END

***********************************************************************
      FUNCTION H4G_SqTriLoop(i1,i2,i3,i4,mass)
*----------------------------------------------------------------------
c... "Triangle" pentagon (has a color structure of a pentagon)
c... author: Michael Kubocz
c... date  : 01.08.07
*----------------------------------------------------------------------

      IMPLICIT NONE
      
#include "ggf_process.inc"
#include "ggf_dotproducts.inc"

      complex*16  H4G_SqTriLoop
      integer i1,i2,i3,i4
      real*8 mass

      real*8 q12sq,q34sq,q1234sq
      complex*16 C0_1M,C0q1234

      EXTERNAL C0_1M

      q12sq=2d0*q_dot_q(i1,i2)
      q34sq=2d0*q_dot_q(i3,i4)
      q1234sq=2d0*(q_dot_q(i1,i2)+q_dot_q(i1,i3)+q_dot_q(i1,i4)
     &          +q_dot_q(i2,i3)+q_dot_q(i2,i4)+q_dot_q(i3,i4))

c... get C0-function      
      C0q1234=C0_1M(mass,q12sq,q34sq,q1234sq) 

      H4G_SqTriLoop=ep_dot_ep(i1,i2)*ep_dot_ep(i3,i4)*C0q1234
      
      RETURN
      END

***********************************************************************
      FUNCTION H4G_SqBoxLoop(i1,i2,i3,i4,mass)
*----------------------------------------------------------------------
c... "Box" Pentagon (has a color structure of a pentagon)
c... author: Michael Kubocz
c... date  : 01.08.07
*----------------------------------------------------------------------

      IMPLICIT NONE
      
#include "d_coefficients.inc"
#include "ggf_dotproducts.inc"

      complex*16 H4G_SqBoxLoop,Box1,Box2,Box3
      integer i1,i2,i3,i4
      real*8 mass

      real*8 q12sq,q3sq,q4sq,q12q3,q12q4,q3q4

      complex*16 D_1M_VAL
      integer D_1M_GET,DijkeyP1,DijkeyP2,DijkeyP3

      EXTERNAL D_1M_VAL,D_1M_GET

      q12sq=2d0*q_dot_q(i1,i2)
      q3sq=q_dot_q(i3,i3)
      q4sq=q_dot_q(i4,i4)

      q12q3=q_dot_q(i1,i3)+q_dot_q(i2,i3)
      q12q4=q_dot_q(i1,i4)+q_dot_q(i2,i4)
      q3q4 =q_dot_q(i3,i4)

c... get D0- and Dij-functions
      DijkeyP1=D_1M_GET(mass,q12sq,q3sq,q4sq,q12q3,q12q4,q3q4)
      DijkeyP2=D_1M_GET(mass,q3sq,q4sq,q12sq,q3q4,q12q3,q12q4)
      DijkeyP3=D_1M_GET(mass,q4sq,q12sq,q3sq,q12q4,q3q4,q12q3)

      Box1=
     & 4d0*ep_dot_ep(i1,i2)*(q_dot_ep(i2,i3)*q_dot_ep(i2,i4)
     & *D_1M_VAL(dd0,DijkeyP1)+q_dot_ep(i2,i3)*q_dot_ep(i3,i4)
     & *D_1M_VAL(dd0,DijkeyP1)+2d0*q_dot_ep(i2,i3)*q_dot_ep(i2,i4)
     & *D_1M_VAL(dd11,DijkeyP1)+q_dot_ep(i2,i3)*q_dot_ep(i3,i4)
     & *D_1M_VAL(dd11,DijkeyP1)+q_dot_ep(i2,i3)*q_dot_ep(i3,i4)
     & *D_1M_VAL(dd12,DijkeyP1)+q_dot_ep(i2,i4)*q_dot_ep(i4,i3)
     & *D_1M_VAL(dd13,DijkeyP1)+q_dot_ep(i3,i4)*q_dot_ep(i4,i3)
     & *D_1M_VAL(dd13,DijkeyP1)-ep_dot_ep(i3,i4)*D_1M_VAL(dd27,DijkeyP1)
     & +q_dot_ep(i2,i3)*q_dot_ep(i2,i4)*D_1M_VAL(dd21,DijkeyP1)
     & +q_dot_ep(i2,i3)*q_dot_ep(i3,i4)*D_1M_VAL(dd24,DijkeyP1)
     & +q_dot_ep(i1,i3)*(q_dot_ep(i1,i4)*(D_1M_VAL(dd0,DijkeyP1)+2d0
     & *D_1M_VAL(dd11,DijkeyP1)+D_1M_VAL(dd21,DijkeyP1))+q_dot_ep(i2,i4)
     & *(D_1M_VAL(dd0,DijkeyP1)+2d0*D_1M_VAL(dd11,DijkeyP1)
     & +D_1M_VAL(dd21,DijkeyP1))+q_dot_ep(i3,i4)*(D_1M_VAL(dd0,DijkeyP1)
     & +D_1M_VAL(dd11,DijkeyP1)+D_1M_VAL(dd12,DijkeyP1)
     & +D_1M_VAL(dd24,DijkeyP1)))+q_dot_ep(i2,i4)*q_dot_ep(i4,i3)
     & *D_1M_VAL(dd25,DijkeyP1)+q_dot_ep(i1,i4)*(q_dot_ep(i2,i3)
     & *(D_1M_VAL(dd0,DijkeyP1)+2d0*D_1M_VAL(dd11,DijkeyP1)
     & +D_1M_VAL(dd21,DijkeyP1))+q_dot_ep(i4,i3)
     & *(D_1M_VAL(dd13,DijkeyP1)+D_1M_VAL(dd25,DijkeyP1)))
     & +q_dot_ep(i3,i4)*q_dot_ep(i4,i3)*D_1M_VAL(dd26,DijkeyP1))

      Box2=
     & 4d0*ep_dot_ep(i1,i2)*(q_dot_ep(i3,i4)*(q_dot_ep(i4,i3)
     & *(D_1M_VAL(dd12,DijkeyP2)+D_1M_VAL(dd24,DijkeyP2))
     & +(q_dot_ep(i1,i3)+q_dot_ep(i2,i3))*(D_1M_VAL(dd13,DijkeyP2)
     & +D_1M_VAL(dd25,DijkeyP2)))+(q_dot_ep(i1,i4)+q_dot_ep(i2,i4))
     & *((q_dot_ep(i1,i3)+q_dot_ep(i2,i3))*D_1M_VAL(dd23,DijkeyP2)
     & +q_dot_ep(i4,i3)*D_1M_VAL(dd26,DijkeyP2))
     & -ep_dot_ep(i3,i4)*D_1M_VAL(dd27,DijkeyP2))

      Box3=
     & 4d0*ep_dot_ep(i1,i2)*((q_dot_ep(i1,i3)+q_dot_ep(i2,i3)
     & +q_dot_ep(i4,i3))*((q_dot_ep(i1,i4)+q_dot_ep(i2,i4))
     & *D_1M_VAL(dd12,DijkeyP3)
     & +q_dot_ep(i3,i4)*D_1M_VAL(dd13,DijkeyP3))+(q_dot_ep(i1,i4)
     & +q_dot_ep(i2,i4))*q_dot_ep(i4,i3)*D_1M_VAL(dd24,DijkeyP3)
     & +q_dot_ep(i3,i4)*q_dot_ep(i4,i3)*D_1M_VAL(dd25,DijkeyP3)
     & +(q_dot_ep(i1,i3)+q_dot_ep(i2,i3))*((q_dot_ep(i1,i4)
     & +q_dot_ep(i2,i4))*D_1M_VAL(dd22,DijkeyP3)+q_dot_ep(i3,i4)
     & *D_1M_VAL(dd26,DijkeyP3))
     & -ep_dot_ep(i3,i4)*D_1M_VAL(dd27,DijkeyP3))

      H4G_SqBoxLoop=Box1+Box2+Box3

      RETURN
      END

***********************************************************************
      FUNCTION H4G_SqLoop(i1,i2,i3,i4,m)
*----------------------------------------------------------------------
c... Pentagon
c... author: Michael Kubocz
c... date  : 01.08.07
*----------------------------------------------------------------------

      IMPLICIT NONE
      
#include "ggf_process.inc"
#include "ggf_dotproducts.inc"
#include "d_coefficients.inc"

      complex*16 H4G_SqLoop
      
      integer i1,i2,i3,i4
      real*8 m
  
      integer init,ie0(4,4,4,4)
      data init/0/, ie0/256*0/
      save ie0,init

      real*8 q1sq,q2sq,q3sq,q4sq
      real*8 q1q2,q1q3,q1q4,q2q3,q2q4,q3q4
      complex*16 q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &           q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4
      complex*16 J1J2,J1J3,J1J4,J2J3,J2J4,J3J4

      complex*16 coeffEmu,coeffEmunu,coeffEmunurho,coeffEmunurhosigma
      complex*16 Eij(4,46)

      integer i,j

      if (init.eq.0) then   
         init = 1
         ie0(4,2,1,3) = 1   ! this is the map from i1,i2,i3,i4 
         ie0(2,1,3,4) = 2   ! to the index of the Eij variables
         ie0(3,2,1,4) = 3
         ie0(2,1,4,3) = 4
         ie0(4,3,2,1) = 5
         ie0(3,4,2,1) = 6
         ie0(2,3,1,4) = 7
         ie0(3,1,4,2) = 8
         ie0(4,2,3,1) = 9
         ie0(2,4,3,1) = 10
         ie0(3,2,4,1) = 11
         ie0(2,3,4,1) = 12
      endif
         
      if (m.eq.m_stop(1)) then
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijstop1(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo  
      endif
      if (m.eq.m_stop(2)) then
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijstop2(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo  
      endif
      if (m.eq.m_sbottom(1)) then
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijsbottom1(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo  
      endif
      if (m.eq.m_sbottom(2)) then
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijsbottom2(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo  
      endif

      q1sq=0d0
      q2sq=0d0
      q3sq=0d0
      q4sq=0d0

      q1q2 = q_dot_q(i1,i2)
      q1q3 = q_dot_q(i1,i3)
      q1q4 = q_dot_q(i1,i4)
      q2q3 = q_dot_q(i2,i3)
      q2q4 = q_dot_q(i2,i4)
      q3q4 = q_dot_q(i3,i4)

      q1J1 = q_dot_ep(i1,i1)
      q2J1 = q_dot_ep(i2,i1)
      q3J1 = q_dot_ep(i3,i1)
      q4J1 = q_dot_ep(i4,i1)

      q1J2 = q_dot_ep(i1,i2)
      q2J2 = q_dot_ep(i2,i2)
      q3J2 = q_dot_ep(i3,i2)
      q4J2 = q_dot_ep(i4,i2)

      q1J3 = q_dot_ep(i1,i3)
      q2J3 = q_dot_ep(i2,i3)
      q3J3 = q_dot_ep(i3,i3)
      q4J3 = q_dot_ep(i4,i3)

      q1J4 = q_dot_ep(i1,i4)
      q2J4 = q_dot_ep(i2,i4)
      q3J4 = q_dot_ep(i3,i4)
      q4J4 = q_dot_ep(i4,i4)

      J1J2 = ep_dot_ep(i1,i2)
      J1J3 = ep_dot_ep(i1,i3)
      J1J4 = ep_dot_ep(i1,i4)
      J2J3 = ep_dot_ep(i2,i3)
      J2J4 = ep_dot_ep(i2,i4)
      J3J4 = ep_dot_ep(i3,i4)

c... coeffE0=0 for qi*Ji=0

      coeffEmu=16d0*q1J2*(q1J3+q2J3)*(q1J4+q2J4+q3J4)*(q2J1*Eij(1,2)
     & +q3J1*Eij(1,3)+q4J1*Eij(1,4))

      coeffEmunu=16d0*(J1J4*q1J2*q1J3*Eij(2,11)+J1J4*q1J2*q2J3*Eij(2,11)
     & +J1J3*q1J2*q1J4*Eij(2,11)+J1J2*q1J3*q1J4*Eij(2,11)+J1J2*q2J3*q1J4
     & *Eij(2,11)+J1J3*q1J2*q2J4*Eij(2,11)+J1J2*q1J3*q2J4*Eij(2,11)+J1J2
     & *q2J3*q2J4*Eij(2,11)+J1J3*q1J2*q3J4*Eij(2,11)+J1J2*q1J3*q3J4
     & *Eij(2,11)+J1J2*q2J3*q3J4*Eij(2,11)+q2J1*q1J2*q2J3*q1J4*Eij(2,2)
     & +q2J1*q1J2*q1J3*q2J4*Eij(2,2)+2d0*q2J1*q1J2*q2J3*q2J4*Eij(2,2)
     & +q2J1*q1J2*q2J3*q3J4*Eij(2,2)+3d0*q2J1*q1J2*q1J3*q1J4*Eij(2,5)
     & +2d0*q2J1*q1J2*q2J3*q1J4*Eij(2,5)+2d0*q2J1*q1J2*q1J3*q2J4
     & *Eij(2,5)+q2J1*q1J2*q2J3*q2J4*Eij(2,5)+2d0*q2J1*q1J2*q1J3*q3J4
     & *Eij(2,5)+q2J1*q1J2*q2J3
     & *q3J4*Eij(2,5)+q2J1*q3J2*q1J3*q1J4*Eij(2,8)+q2J1*q3J2*q2J3*q1J4
     & *Eij(2,8)+q2J1*q3J2*q1J3*q2J4*Eij(2,8)+q2J1*q3J2*q2J3*q2J4
     & *Eij(2,8)+q2J1*q1J2*q1J3*q3J4*Eij(2,8)+q2J1*q3J2*q1J3*q3J4
     & *Eij(2,8)+q2J1*q1J2*q2J3*q3J4*Eij(2,8)+q2J1*q3J2*q2J3*q3J4
     & *Eij(2,8)+q3J1*(q4J2*(q1J3+q2J3)*(q1J4+q2J4+q3J4)*Eij(2,10)+q3J2
     & *(q1J3+q2J3)*(q1J4+q2J4+q3J4)*Eij(2,3)+q1J2*(q4J3*(q1J4+q2J4
     & +q3J4)*Eij(2,10)+q1J3*(3d0*q1J4*Eij(2,6)+q3J4*(Eij(2,3)+2d0
     & *Eij(2,6))+q2J4*(2d0*Eij(2,6)+Eij(2,8)))+q2J3*(q3J4*(Eij(2,3)
     & +Eij(2,6)+Eij(2,8))+q1J4*(2d0*Eij(2,6)
     & +Eij(2,8))+q2J4*(Eij(2,6)+2d0
     & *Eij(2,8)))))+q2J1*(q4J2*(q1J3+q2J3)+q1J2*q4J3)*(q1J4+q2J4+q3J4)
     & *Eij(2,9)+q4J1*(q3J2*(q1J3+q2J3)*(q1J4+q2J4+q3J4)*Eij(2,10)+q4J2
     & *(q1J3+q2J3)*(q1J4+q2J4+q3J4)*Eij(2,4)+q1J2*(q4J3*(q1J4+q2J4
     & +q3J4)*Eij(2,4)+q1J3*(3d0*q1J4*Eij(2,7)+q3J4*(Eij(2,10)+2d0
     & *Eij(2,7))+q2J4*(2d0*Eij(2,7)+Eij(2,9)))+q2J3*(q3J4*(Eij(2,10)
     & +Eij(2,7)+Eij(2,9))+q1J4*(2d0*Eij(2,7)+Eij(2,9))+q2J4*(Eij(2,7)
     & +2d0*Eij(2,9))))))

      coeffEmunurho= 16d0*(2d0*q3J1*q3J2*q1J3*q1J4*Eij(3,11)+q3J1*q3J2
     & *q2J3*q1J4*Eij(3,11)+q3J1*q3J2*q1J3
     & *q2J4*Eij(3,11)+2d0*q3J1*q1J2*q1J3
     & *q3J4*Eij(3,11)+q3J1*q3J2*q1J3*q3J4*Eij(3,11)+q3J1*q1J2*q2J3*q3J4
     & *Eij(3,11)+q3J1*q3J2*q2J3*q1J4*Eij(3,12)+q3J1*q3J2*q1J3*q2J4
     & *Eij(3,12)+2d0*q3J1*q3J2*q2J3*q2J4*Eij(3,12)+q3J1*q1J2*q2J3*q3J4
     & *Eij(3,12)+q3J1*q3J2*q2J3*q3J4*Eij(3,12)+q3J1*q3J2*q4J3*q1J4
     & *Eij(3,13)+q3J1*q3J2*q4J3*q2J4*Eij(3,13)+q3J1*q4J2*q1J3*q3J4
     & *Eij(3,13)+q3J1*q4J2*q2J3*q3J4*Eij(3,13)+q3J1*q1J2*q4J3*q3J4
     & *Eij(3,13)+q3J1*q3J2*q4J3*q3J4*Eij(3,13)+q3J1*q4J2*q4J3*q1J4
     & *Eij(3,16)+q3J1*q4J2*q4J3*q2J4*Eij(3,16)+q3J1*q4J2*q4J3*q3J4
     & *Eij(3,16)+2d0*q3J1*q1J2*q2J3*q1J4
     & *Eij(3,17)+2d0*q3J1*q1J2*q1J3*q2J4
     & *Eij(3,17)+2d0*q3J1*q1J2*q2J3*q2J4*Eij(3,17)+q3J1*q1J2*q2J3*q3J4
     & *Eij(3,17)+2d0*q3J1*q4J2*q1J3*q1J4*Eij(3,19)+q3J1*q4J2*q2J3*q1J4
     & *Eij(3,19)+2d0*q3J1*q1J2*q4J3*q1J4*Eij(3,19)+q3J1*q4J2*q1J3*q2J4
     & *Eij(3,19)+q3J1*q1J2*q4J3*q2J4*Eij(3,19)+q3J1*q4J2*q1J3*q3J4
     & *Eij(3,19)+q3J1*q1J2*q4J3*q3J4*Eij(3,19)+q3J1*q4J2*q2J3*q1J4
     & *Eij(3,20)+q3J1*q4J2*q1J3*q2J4*Eij(3,20)+2d0*q3J1*q4J2*q2J3*q2J4
     & *Eij(3,20)+q3J1*q1J2*q4J3*q2J4*Eij(3,20)+q3J1*q4J2*q2J3*q3J4
     & *Eij(3,20)+2d0*J1J4*q1J2*q1J3*Eij(3,21)
     & +J1J4*q1J2*q2J3*Eij(3,21)+2d0
     & *J1J3*q1J2*q1J4*Eij(3,21)+2d0*J1J2*q1J3*q1J4*Eij(3,21)+J1J2*q2J3
     & *q1J4*Eij(3,21)+J1J3*q1J2*q2J4*Eij(3,21)+J1J2*q1J3*q2J4*Eij(3,21)
     & +J1J3*q1J2*q3J4*Eij(3,21)+J1J2*q1J3*q3J4*Eij(3,21)+J1J4*q1J2*q2J3
     & *Eij(3,22)+J1J2*q2J3*q1J4*Eij(3,22)+J1J3*q1J2*q2J4*Eij(3,22)+J1J2
     & *q1J3*q2J4*Eij(3,22)+2d0*J1J2*q2J3*q2J4*Eij(3,22)+J1J2*q2J3*q3J4
     & *Eij(3,22)+q3J1*q1J2*J3J4*Eij(3,23)+q3J1*J2J4*q1J3*Eij(3,23)+J1J4
     & *q3J2*q1J3*Eij(3,23)+q3J1*J2J4*q2J3*Eij(3,23)+J1J4*q3J2*q2J3
     & *Eij(3,23)+q3J1*J2J3*q1J4*Eij(3,23)+J1J3*q3J2*q1J4*Eij(3,23)+q3J1
     & *J2J3*q2J4*Eij(3,23)+J1J3*q3J2*q2J4*Eij(3,23)+q3J1*J2J3*q3J4
     & *Eij(3,23)+J1J3*q1J2*q3J4*Eij(3,23)+J1J3*q3J2*q3J4*Eij(3,23)+J1J2
     & *q1J3*q3J4*Eij(3,23)+J1J2*q2J3*q3J4*Eij(3,23)+J1J4*q4J2*q1J3
     & *Eij(3,24)+J1J4*q4J2*q2J3*Eij(3,24)+J1J4*q1J2*q4J3*Eij(3,24)+J1J3
     & *q4J2*q1J4*Eij(3,24)+J1J2*q4J3*q1J4*Eij(3,24)+J1J3*q4J2*q2J4
     & *Eij(3,24)+J1J2*q4J3*q2J4*Eij(3,24)+J1J3*q4J2*q3J4*Eij(3,24)+J1J2
     & *q4J3*q3J4*Eij(3,24)+q3J1*q3J2*q1J3*q3J4*Eij(3,3)+q3J1*q3J2*q2J3
     & *q3J4*Eij(3,3)+3d0*q3J1*q1J2*q1J3
     & *q1J4*Eij(3,6)+q3J1*q1J2*q2J3*q1J4
     & *Eij(3,6)+q3J1*q1J2*q1J3*q2J4*Eij(3,6)+q3J1*q1J2*q1J3*q3J4
     & *Eij(3,6)+q4J1*(2d0*q4J2*q1J3*q1J4*Eij(3,14)+q4J2*q2J3*q1J4
     & *Eij(3,14)+q4J2*q1J3*q2J4*Eij(3,14)+q4J2*q1J3*q3J4*Eij(3,14)+q4J2
     & *q2J3*q1J4*Eij(3,15)+q4J2*q1J3*q2J4*Eij(3,15)+2d0*q4J2*q2J3*q2J4
     & *Eij(3,15)+q4J2*q2J3*q3J4*Eij(3,15)+q4J2*q1J3*q3J4*Eij(3,16)+q4J2
     & *q2J3*q3J4*Eij(3,16)+q3J2*(q4J3*(q1J4+q2J4+q3J4)*Eij(3,16)+q2J3
     & *(2d0*q2J4*Eij(3,20)+q3J4*(Eij(3,13)+Eij(3,20))+q1J4*(Eij(3,19)
     & +Eij(3,20)))+q1J3*(2d0*q1J4*Eij(3,19)+q3J4*(Eij(3,13)+Eij(3,19))
     & +q2J4*(Eij(3,19)+Eij(3,20))))+J2J4*q1J3*Eij(3,24)+J2J4*q2J3
     & *Eij(3,24)+J2J3*q1J4*Eij(3,24)+J2J3*q2J4*Eij(3,24)+J2J3*q3J4
     & *Eij(3,24)+q4J2*q4J3*(q1J4+q2J4
     & +q3J4)*Eij(3,4)+q1J2*(q4J3*(2d0*q1J4
     & *Eij(3,14)+q2J4*(Eij(3,14)+Eij(3,15))+q3J4*(Eij(3,14)+Eij(3,16)))
     & +J3J4*Eij(3,24)+q2J3*(q2J4*(Eij(3,10)+2d0*Eij(3,18))+q3J4
     & *(Eij(3,18)+Eij(3,19)+Eij(3,20))+q1J4*(2d0*Eij(3,18)+Eij(3,7)))
     & +q1J3*(3d0*q1J4*Eij(3,7)+q2J4*(2d0*Eij(3,18)+Eij(3,7))+q3J4*(2d0
     & *Eij(3,19)+Eij(3,7)))))+q3J1*q1J2*q2J3*q2J4*Eij(3,9)+q2J1*(q3J2
     & *q1J3*q2J4*Eij(3,17)+q3J2*q1J3*q3J4*Eij(3,17)+q3J2*((q1J3+q2J3)
     & *q3J4*Eij(3,12)+(2d0*q1J3+q2J3)*q1J4*Eij(3,17))+q3J2*q4J3*q1J4
     & *Eij(3,20)+q3J2*q4J3*q2J4*Eij(3,20)+q3J2*q4J3*q3J4*Eij(3,20)+q4J2
     & *(q4J3*(q1J4+q2J4+q3J4)*Eij(3,15)+q2J3*(2d0*q2J4*Eij(3,10)+q1J4
     & *(Eij(3,10)+Eij(3,18))+q3J4*(Eij(3,10)+Eij(3,20)))+q1J3*(2d0*q1J4
     & *Eij(3,18)+q2J4*(Eij(3,10)+Eij(3,18))+q3J4*(Eij(3,18)
     & +Eij(3,20))))+J2J4*q1J3*Eij(3,22)+J2J4*q2J3*Eij(3,22)+J2J3*q1J4
     & *Eij(3,22)+J2J3*q2J4*Eij(3,22)+J2J3*q3J4*Eij(3,22)+q3J2*(q1J3
     & *q2J4+q2J3*(q1J4+2d0*q2J4+q3J4))*Eij(3,9)+q1J2*(q4J3*(2d0*q1J4
     & *Eij(3,18)+q2J4*(Eij(3,10)+Eij(3,18))+q3J4*(Eij(3,18)+Eij(3,20)))
     & +J3J4*Eij(3,22)+q1J3*(3d0*q1J4*Eij(3,5)+q3J4*(2d0*Eij(3,17)
     & +Eij(3,5))+q2J4*(Eij(3,5)+2d0*Eij(3,8)))+q2J3*(q2J4*(Eij(3,2)
     & +2d0*Eij(3,8))+q1J4*(Eij(3,5)+2d0*Eij(3,8))+q3J4*(Eij(3,17)
     & +Eij(3,8)+Eij(3,9))))))

      coeffEmunurhosigma=16d0*(q4J1*q4J2*q4J3*q1J4*Eij(4,14)+q4J1*q4J2
     & *q4J3*q2J4*Eij(4,15)+q4J1*q4J2*q4J3*q3J4*Eij(4,16)+q4J1*q4J2*q1J3
     & *q1J4*Eij(4,19)+q4J1*q1J2*q4J3*q1J4*Eij(4,19)+q4J1*q4J2*q2J3*q2J4
     & *Eij(4,21)+q4J1*q3J2*q4J3*q3J4*Eij(4,22)+q4J1*q1J2*q2J3*q1J4
     & *Eij(4,24)+q4J1*q1J2*q1J3*q2J4*Eij(4,24)+q4J1*q3J2*q1J3*q1J4
     & *Eij(4,25)+q4J1*q1J2*q1J3*q3J4*Eij(4,25)+q4J1*q1J2*q2J3*q2J4
     & *Eij(4,27)+q4J1*q3J2*q2J3*q2J4*Eij(4,28)+q4J1*q3J2*q1J3*q3J4
     & *Eij(4,30)+q4J1*q3J2*q2J3*q3J4*Eij(4,31)+q4J1*q4J2*q2J3*q1J4
     & *Eij(4,32)+q4J1*q4J2*q1J3*q2J4*Eij(4,32)+q4J1*q1J2*q4J3*q2J4
     & *Eij(4,32)+q4J1*q3J2*q4J3*q1J4*Eij(4,33)+q4J1*q4J2*q1J3*q3J4
     & *Eij(4,33)+q4J1*q1J2*q4J3*q3J4*Eij(4,33)+q4J1*q3J2*q4J3*q2J4
     & *Eij(4,34)+q4J1*q4J2*q2J3*q3J4*Eij(4,34)+q4J1*q3J2*q2J3*q1J4
     & *Eij(4,35)+q4J1*q3J2*q1J3*q2J4*Eij(4,35)+q4J1*q1J2*q2J3*q3J4
     & *Eij(4,35)+J1J4*q1J2*q1J3*Eij(4,36)+J1J3*q1J2*q1J4*Eij(4,36)+J1J2
     & *q1J3*q1J4*Eij(4,36)+J1J2*q2J3*q2J4*Eij(4,37)+J1J3*q3J2*q3J4
     & *Eij(4,38)+q4J1*q4J2*J3J4*Eij(4,39)+q4J1*J2J4*q4J3*Eij(4,39)
     & +J1J4*q4J2*q4J3*Eij(4,39)+J1J4*q1J2*q2J3*Eij(4,40)+J1J2*q2J3*q1J4
     & *Eij(4,40)+J1J3*q1J2*q2J4*Eij(4,40)+J1J2*q1J3*q2J4*Eij(4,40)+J1J4
     & *q3J2*q1J3*Eij(4,41)+J1J3*q3J2*q1J4*Eij(4,41)+J1J3*q1J2*q3J4
     & *Eij(4,41)+J1J2*q1J3*q3J4*Eij(4,41)+q4J1*q1J2*J3J4*Eij(4,42)+q4J1
     & *J2J4*q1J3*Eij(4,42)+J1J4*q4J2*q1J3*Eij(4,42)+J1J4*q1J2*q4J3
     & *Eij(4,42)+q4J1*J2J3*q1J4*Eij(4,42)+J1J3*q4J2*q1J4*Eij(4,42)+J1J2
     & *q4J3*q1J4*Eij(4,42)+J1J4*q3J2*q2J3*Eij(4,43)+J1J3*q3J2*q2J4
     & *Eij(4,43)+J1J2*q2J3*q3J4*Eij(4,43)+q4J1*J2J4*q2J3*Eij(4,44)+J1J4
     & *q4J2*q2J3*Eij(4,44)+q4J1*J2J3*q2J4*Eij(4,44)+J1J3*q4J2*q2J4
     & *Eij(4,44)+J1J2*q4J3*q2J4*Eij(4,44)+q4J1*q3J2*J3J4*Eij(4,45)+J1J4
     & *q3J2*q4J3*Eij(4,45)+q4J1*J2J3*q3J4*Eij(4,45)+J1J3*q4J2*q3J4
     & *Eij(4,45)+J1J2*q4J3*q3J4*Eij(4,45)+J1J4*J2J3*Eij(4,46)+J1J3*J2J4
     & *Eij(4,46)+J1J2*J3J4*Eij(4,46)+q3J1*(q4J2*q4J3*q3J4*Eij(4,22)
     & +q4J2*q1J3*q1J4*Eij(4,25)+q4J2*q2J3*q2J4*Eij(4,28)+q4J2*q1J3*q3J4
     & *Eij(4,30)+q4J2*q2J3*q3J4*Eij(4,31)+q4J2*q4J3*q1J4*Eij(4,33)+q4J2
     & *q4J3*q2J4*Eij(4,34)+q4J2*q2J3*q1J4*Eij(4,35)+q4J2*q1J3*q2J4
     & *Eij(4,35)+J2J3*q3J4*Eij(4,38)+q3J2*(q2J3*(q3J4*Eij(4,12)+q2J4
     & *Eij(4,20)+q1J4*Eij(4,29))+q1J3*(q3J4*Eij(4,11)+q1J4*Eij(4,18)
     & +q2J4*Eij(4,29))+q4J3*(q3J4*Eij(4,13)+q1J4*Eij(4,30)+q2J4
     & *Eij(4,31))+J3J4*Eij(4,38))+J2J4*q1J3*Eij(4,41)+J2J3*q1J4
     & *Eij(4,41)+J2J4*q2J3*Eij(4,43)+J2J3*q2J4*Eij(4,43)+(q4J2*J3J4
     & +J2J4*q4J3)*Eij(4,45)+q1J2*(q2J3*(q1J4*Eij(4,23)+q2J4*Eij(4,26)
     & +q3J4*Eij(4,29))+q4J3*(q1J4*Eij(4,25)+q3J4*Eij(4,30)+q2J4
     & *Eij(4,35))+J3J4*Eij(4,41)+q1J3*(q3J4*Eij(4,18)+q2J4*Eij(4,23)
     & +q1J4*Eij(4,6))))+q4J1*q1J2*q1J3*q1J4*Eij(4,7)+q2J1*(q3J2*q2J3
     & *q3J4*Eij(4,20)+q3J2*q1J3*q1J4*Eij(4,23)+q3J2*q2J3*q1J4*Eij(4,26)
     & +q3J2*q1J3*q2J4*Eij(4,26)+q3J2*q4J3*q2J4*Eij(4,28)+q3J2*q1J3*q3J4
     & *Eij(4,29)+q3J2*q4J3*q3J4*Eij(4,31)+q3J2*q4J3*q1J4*Eij(4,35)+J2J4
     & *q2J3*Eij(4,37)+J2J3*q2J4*Eij(4,37)+J2J4*q1J3*Eij(4,40)+J2J3*q1J4
     & *Eij(4,40)+q3J2*J3J4*Eij(4,43)+J2J3*q3J4*Eij(4,43)+J2J4*q4J3
     & *Eij(4,44)+q4J2*(q2J3*(q2J4*Eij(4,10)+q1J4*Eij(4,27)+q3J4
     & *Eij(4,28))+q4J3*(q2J4*Eij(4,21)+q1J4*Eij(4,32)+q3J4*Eij(4,34))
     & +q1J3*(q1J4*Eij(4,24)+q2J4*Eij(4,27)+q3J4*Eij(4,35))+J3J4
     & *Eij(4,44))+q1J2*(q4J3*(q1J4*Eij(4,24)+q2J4*Eij(4,27)+q3J4
     & *Eij(4,35))+J3J4*Eij(4,40)+q1J3*(q2J4*Eij(4,17)+q3J4*Eij(4,23)
     & +q1J4*Eij(4,5))+q2J3*(q1J4*Eij(4,17)+q3J4*Eij(4,26)+q2J4
     & *Eij(4,8)))+q3J2*q2J3*q2J4*Eij(4,9)))

      H4G_SqLoop=-(coeffEmu+coeffEmunu+coeffEmunurho+coeffEmunurhosigma)

      RETURN
      END
