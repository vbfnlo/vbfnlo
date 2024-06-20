***********************************************************************
      FUNCTION A4G_LoopDD(i1,i2,i3,i4,m)
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... date     : 01.08.08
c... modified : 30.01.11 
*----------------------------------------------------------------------
c... This pentagon subroutine uses the Denner-Dittmaier tensor 
c... decomposition.      
c... (i1,i2,i3,i4) gives the permutation of the attached four gluons.
c... Input gluon momenta and polarization vectors are stored via the
c... dot products in q_dot_q, q_dot_ep, ep_dot_ep and explicitly 
c... in qi and eps
c... Normalization: this subroutine gives the pentagon contributions
c... without a factor 4*i*mt, the QCD couplings and Yukawa coupling 
c... h_t (gamma_5). 

c...                  * -- i1    
c...                /   \
c...         i2 -- *     \
c...               |      * -- A^0
c...         i3 -- *     /
c...                \   /
c...                  * -- i4

***********************************************************************

      IMPLICIT NONE

#include "VBFNLO/amplitudes/ggf/ggf_process.inc"
#include "VBFNLO/amplitudes/ggf/d_coefficients.inc"
#include "VBFNLO/amplitudes/ggf/ggf_dotproducts.inc"

      complex*16 A4G_LoopDD

c... four vectors
      real*8 eps,qi,m
      COMMON /bh4g_mom_eps/ qi(0:3,4),eps(0:3,4)

c... explicit momenta for call of eps tensor
      real*8  q1(0:3),q2(0:3),q3(0:3),q4(0:3),
     &        mu1(0:3),mu2(0:3),mu3(0:3),mu4(0:3),
     &     pA1(0:3)

      real*8 Kpen(0:4,4)
      complex*16 epsgPen(0:3,4)
      Common/Forpent/Kpen,epsgPen


c... Passarino-Veltman/Denner-Dittmaier form factors
      complex*16 C0_1M,C034,
     &           D0_1M,D0234,Dij(3,13), 
     &           E0L,Eij(4,46)

c... local dot products 
      real*8 q1sq,q2sq,q3sq,q4sq,q12sq,q34sq,q1234sq,
     &       q1q2,q1q3,q1q4,q2q3,q2q4,q3q4
      
      complex*16 q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &           q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
     &           J1J2,J1J3,J1J4,J2J3,J2J4,J3J4
      
c... other variables
      integer i1,i2,i3,i4
      integer i,j,mu


      integer init,ie0(4,4,4,4)
      data init/0/, ie0/256*0/
      save ie0,init

      complex*16 coeffEmu,coeffEmunu,coeffD0,coeffDmu,coeffE0C0,coeffE0

      real*8 dotrr,epsrrrr

      EXTERNAL dotrr,epsrrrr,C0_1M,D0_1M
      complex*16 pen1

C     BEGIN

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
 
       real*8 gaugelimit 
       common/gauge/gaugelimit
       common/dbleswitch/dble1,quad
cfc
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
cfc            do j=1,46 
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
cfc      do mu = 0,3
cfc         q1(mu) = qi(mu,i1)
cfc         q2(mu) = qi(mu,i2)
cfc         q3(mu) = qi(mu,i3)
cfc         q4(mu) = qi(mu,i4)
cfc         pa1(mu) = -(q1(mu)+q2(mu)+q3(mu)+q4(mu))
cfc         
cfc
cfc         mu1(mu) = eps(mu,i1)
cfc         mu2(mu) = eps(mu,i2)
cfc         mu3(mu) = eps(mu,i3)
cfc         mu4(mu) = eps(mu,i4)
cfc      enddo
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
cfc      q1sq=q_dot_q(i1,i1)
cfc      q2sq=q_dot_q(i2,i2)
cfc      q3sq=q_dot_q(i3,i3)
cfc      q4sq=q_dot_q(i4,i4)
cfc
cfc      q12sq=2d0*q1q2
cfc      q34sq=2d0*q3q4
cfc      q1234sq=2*(q1q2+q1q3+q1q4+q2q3+q2q4+q3q4)
cfc
cfc      C034=C0_1M(m,q3sq,q4sq,q34sq) 
cfc
cfc      D0234 = D0_1M(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4)
cfc      CALL D_1M_ALL2(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4,Dij)
cfc
cfc      coeffD0=-4*m*((q1q2+q1q4-2*q1sq+2*q2q3+q2sq+q3q4+q3sq)*
cfc     &    epsrrrr(mu1,mu2,mu3,mu4)-(q2J4+2*q3J4)*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   q1J4*epsrrrr(mu1,mu2,mu3,q2)-2*(q2J4+q3J4)*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   (q1J4+q3J4)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   (-q2J3+q4J3)*epsrrrr(mu1,mu2,mu4,q1)+q1J3*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   q4J3*epsrrrr(mu1,mu2,mu4,q3)+(q1J3-2*q2J3)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   J3J4*(epsrrrr(mu1,mu2,q1,q2)+epsrrrr(mu1,mu2,q1,q4)+
cfc     &      epsrrrr(mu1,mu2,q3,q4))-q1J2*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   q4J2*(epsrrrr(mu1,mu3,mu4,q1)+epsrrrr(mu1,mu3,mu4,q3))+
cfc     &   (-q1J2+q3J2)*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   J2J4*(epsrrrr(mu1,mu3,q1,q2)+epsrrrr(mu1,mu3,q1,q4)+
cfc     &      epsrrrr(mu1,mu3,q3,q4))+
cfc     &   J2J3*(epsrrrr(mu1,mu4,q1,q2)+epsrrrr(mu1,mu4,q1,q4)+
cfc     &      epsrrrr(mu1,mu4,q3,q4))+(q2J1+q4J1)*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   q4J1*epsrrrr(mu2,mu3,mu4,q3)-q3J1*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   J1J4*(epsrrrr(mu2,mu3,q1,q2)+epsrrrr(mu2,mu3,q1,q4)+
cfc     &      epsrrrr(mu2,mu3,q3,q4))-
cfc     &   J1J3*(epsrrrr(mu2,mu4,q1,q2)+epsrrrr(mu2,mu4,q1,q4)+
cfc     &      epsrrrr(mu2,mu4,q3,q4))+
cfc     &   J1J2*(epsrrrr(mu3,mu4,q1,q2)+epsrrrr(mu3,mu4,q1,q4)+
cfc     &      epsrrrr(mu3,mu4,q3,q4)))*D0234
cfc
cfc      coeffDmu=8*m*(-2*q1q2*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q2q3*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q2sq*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q2q3*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q3sq*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q2q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q3q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)+
cfc     &   q2J4*Dij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q3J4*Dij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*Dij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q3J4*Dij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q2J3*Dij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q4J3*Dij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q3J2*Dij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q4J2*Dij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   (q2J1*Dij(1,1)+q3J1*Dij(1,2)+q4J1*Dij(1,3))*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   (q2J1*Dij(1,1)+q3J1*Dij(1,2)+q4J1*Dij(1,3))*epsrrrr(mu2,mu3,mu4,q4))
cfc
cfc      coeffE0C0=4*m*(-C034+E0L*m**4)*epsrrrr(mu1,mu2,mu3,mu4)
cfc
cfc      coeffE0=-4*E0L*m*((m**4+2*q1q2**2+q1sq**2+q1q4*q2q3-q1q3*q2q4+q1q4*q2sq+
cfc     &      q1sq*(3*q1q2+2*q1q3+q1q4+2*q2q3+q2sq+q3q4+q3sq)+
cfc     &      q1q2*(2*q1q3+2*q1q4+2*q2q3+q2sq+q3q4+q3sq))*
cfc     &    epsrrrr(mu1,mu2,mu3,mu4)-4*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q1q2*q2J4*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q1)-q1sq*q2J4*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q1)-q2J4*q2sq*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   4*q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q2q3*q3J4*epsrrrr(mu1,mu2,mu3,q1)-q2q4*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   2*q2sq*q3J4*epsrrrr(mu1,mu2,mu3,q1)+q2J4*q3q4*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q3sq*epsrrrr(mu1,mu2,mu3,q1)-2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q2)+2*q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   2*q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q2)-q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   2*q1q3*q3J4*epsrrrr(mu1,mu2,mu3,q2)+q1q4*q3J4*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   q1J4*q3q4*epsrrrr(mu1,mu2,mu3,q2)-q1J4*q3sq*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1q2*q2J4*epsrrrr(mu1,mu2,mu3,q3)-q1q4*q2J4*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1sq*q2J4*epsrrrr(mu1,mu2,mu3,q3)+q1J4*q2q4*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q3)-
cfc     &   2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q4)-q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q4)-q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q4)-
cfc     &   q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q4)-q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q4)-
cfc     &   q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q4)-2*q1q2*q2J3*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   2*q1q3*q2J3*epsrrrr(mu1,mu2,mu4,q1)-q1sq*q2J3*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   2*q2J3*q2q3*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   2*q1J3*q2q4*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   2*q2J3*q2q4*epsrrrr(mu1,mu2,mu4,q1)-q2J3*q2sq*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q2J3*q3q4*epsrrrr(mu1,mu2,mu4,q1)-q2J3*q3sq*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   2*q1q2*q4J3*epsrrrr(mu1,mu2,mu4,q1)+q1sq*q4J3*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2q3*q4J3*epsrrrr(mu1,mu2,mu4,q1)+q2sq*q4J3*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   2*q1J3*q1q2*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   2*q1J3*q1q3*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   2*q1J3*q1q4*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q1sq*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   2*q1q4*q2J3*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   2*q1J3*q2q3*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q2sq*epsrrrr(mu1,mu2,mu4,q2)+
cfc     &   q1J3*q3q4*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q3sq*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1q3*q4J3*epsrrrr(mu1,mu2,mu4,q2)+q1q4*q2J3*epsrrrr(mu1,mu2,mu4,q3)-
cfc     &   q1J3*q2q4*epsrrrr(mu1,mu2,mu4,q3)+q1q2*q4J3*epsrrrr(mu1,mu2,mu4,q3)+
cfc     &   q1sq*q4J3*epsrrrr(mu1,mu2,mu4,q3)-q1J3*q1sq*epsrrrr(mu1,mu2,mu4,q4)-
cfc     &   2*q1q2*q2J3*epsrrrr(mu1,mu2,mu4,q4)-q1q3*q2J3*epsrrrr(mu1,mu2,mu4,q4)-
cfc     &   2*q1sq*q2J3*epsrrrr(mu1,mu2,mu4,q4)+q1J3*q2q3*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q2sq*epsrrrr(mu1,mu2,mu4,q4)-4*q1J3*q1J4*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*J3J4*q1q2*epsrrrr(mu1,mu2,q1,q2)+2*J3J4*q1q3*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   J3J4*q1sq*epsrrrr(mu1,mu2,q1,q2)-4*q1J4*q2J3*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   4*q1J3*q2J4*epsrrrr(mu1,mu2,q1,q2)-4*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*J3J4*q2q3*epsrrrr(mu1,mu2,q1,q2)+J3J4*q2sq*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   4*q1J3*q3J4*epsrrrr(mu1,mu2,q1,q2)-4*q2J3*q3J4*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   J3J4*q3q4*epsrrrr(mu1,mu2,q1,q2)+J3J4*q3sq*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   q3J4*q4J3*epsrrrr(mu1,mu2,q1,q2)-2*q1J4*q2J3*epsrrrr(mu1,mu2,q1,q3)-
cfc     &   2*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q3)-J3J4*q2q4*epsrrrr(mu1,mu2,q1,q3)-
cfc     &   2*q2J3*q3J4*epsrrrr(mu1,mu2,q1,q3)+q2J4*q4J3*epsrrrr(mu1,mu2,q1,q3)+
cfc     &   2*J3J4*q1q2*epsrrrr(mu1,mu2,q1,q4)+J3J4*q1sq*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   2*q1J3*q2J4*epsrrrr(mu1,mu2,q1,q4)-2*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q4)+
cfc     &   J3J4*q2q3*epsrrrr(mu1,mu2,q1,q4)+J3J4*q2sq*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   q2J3*q3J4*epsrrrr(mu1,mu2,q1,q4)+2*q1J3*q1J4*epsrrrr(mu1,mu2,q2,q3)+
cfc     &   J3J4*q1q4*epsrrrr(mu1,mu2,q2,q3)+2*q1J3*q2J4*epsrrrr(mu1,mu2,q2,q3)+
cfc     &   2*q1J3*q3J4*epsrrrr(mu1,mu2,q2,q3)-q1J4*q4J3*epsrrrr(mu1,mu2,q2,q3)+
cfc     &   2*q1J3*q1J4*epsrrrr(mu1,mu2,q2,q4)-J3J4*q1q3*epsrrrr(mu1,mu2,q2,q4)+
cfc     &   2*q1J4*q2J3*epsrrrr(mu1,mu2,q2,q4)+q1J3*q3J4*epsrrrr(mu1,mu2,q2,q4)+
cfc     &   J3J4*q1q2*epsrrrr(mu1,mu2,q3,q4)+J3J4*q1sq*epsrrrr(mu1,mu2,q3,q4)+
cfc     &   q1J4*q2J3*epsrrrr(mu1,mu2,q3,q4)-q1J3*q2J4*epsrrrr(mu1,mu2,q3,q4)-
cfc     &   4*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   4*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   4*q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q1)+q2q4*q3J2*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1J2*q3q4*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1J2*q3sq*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1q2*q4J2*epsrrrr(mu1,mu3,mu4,q1)-q1sq*q4J2*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   q2q3*q4J2*epsrrrr(mu1,mu3,mu4,q1)-q2sq*q4J2*epsrrrr(mu1,mu3,mu4,q1)-
cfc     &   2*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   2*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   2*q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   q1q4*q3J2*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q3q4*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   q1J2*q3sq*epsrrrr(mu1,mu3,mu4,q2)+q1q3*q4J2*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   2*q1J2*q1q4*epsrrrr(mu1,mu3,mu4,q3)+q1J2*q2q4*epsrrrr(mu1,mu3,mu4,q3)-
cfc     &   q1q2*q4J2*epsrrrr(mu1,mu3,mu4,q3)-q1sq*q4J2*epsrrrr(mu1,mu3,mu4,q3)-
cfc     &   2*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   2*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q4)-q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q4)-q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1q2*q3J2*epsrrrr(mu1,mu3,mu4,q4)+q1sq*q3J2*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   2*J2J4*q1q2*epsrrrr(mu1,mu3,q1,q2)-2*J2J4*q1q3*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   J2J4*q1sq*epsrrrr(mu1,mu3,q1,q2)-2*J2J4*q2q3*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   J2J4*q2sq*epsrrrr(mu1,mu3,q1,q2)+2*q1J4*q3J2*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   2*q2J4*q3J2*epsrrrr(mu1,mu3,q1,q2)+2*q3J2*q3J4*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   J2J4*q3q4*epsrrrr(mu1,mu3,q1,q2)-J2J4*q3sq*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   q3J4*q4J2*epsrrrr(mu1,mu3,q1,q2)-4*q1J2*q1J4*epsrrrr(mu1,mu3,q1,q3)-
cfc     &   4*q1J2*q2J4*epsrrrr(mu1,mu3,q1,q3)+J2J4*q2q4*epsrrrr(mu1,mu3,q1,q3)-
cfc     &   4*q1J2*q3J4*epsrrrr(mu1,mu3,q1,q3)-q2J4*q4J2*epsrrrr(mu1,mu3,q1,q3)-
cfc     &   2*J2J4*q1q2*epsrrrr(mu1,mu3,q1,q4)-J2J4*q1sq*epsrrrr(mu1,mu3,q1,q4)-
cfc     &   J2J4*q2q3*epsrrrr(mu1,mu3,q1,q4)-J2J4*q2sq*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   q2J4*q3J2*epsrrrr(mu1,mu3,q1,q4)-2*q1J2*q3J4*epsrrrr(mu1,mu3,q1,q4)-
cfc     &   2*q1J2*q1J4*epsrrrr(mu1,mu3,q2,q3)-J2J4*q1q4*epsrrrr(mu1,mu3,q2,q3)-
cfc     &   2*q1J2*q2J4*epsrrrr(mu1,mu3,q2,q3)-2*q1J2*q3J4*epsrrrr(mu1,mu3,q2,q3)+
cfc     &   q1J4*q4J2*epsrrrr(mu1,mu3,q2,q3)+J2J4*q1q3*epsrrrr(mu1,mu3,q2,q4)-
cfc     &   q1J4*q3J2*epsrrrr(mu1,mu3,q2,q4)-q1J2*q3J4*epsrrrr(mu1,mu3,q2,q4)+
cfc     &   2*q1J2*q1J4*epsrrrr(mu1,mu3,q3,q4)-J2J4*q1q2*epsrrrr(mu1,mu3,q3,q4)-
cfc     &   J2J4*q1sq*epsrrrr(mu1,mu3,q3,q4)+q1J2*q2J4*epsrrrr(mu1,mu3,q3,q4)+
cfc     &   2*J2J3*q1q2*epsrrrr(mu1,mu4,q1,q2)+2*J2J3*q1q3*epsrrrr(mu1,mu4,q1,q2)+
cfc     &   J2J3*q1sq*epsrrrr(mu1,mu4,q1,q2)+2*J2J3*q2q3*epsrrrr(mu1,mu4,q1,q2)+
cfc     &   J2J3*q2sq*epsrrrr(mu1,mu4,q1,q2)+J2J3*q3q4*epsrrrr(mu1,mu4,q1,q2)+
cfc     &   J2J3*q3sq*epsrrrr(mu1,mu4,q1,q2)+2*q1J3*q4J2*epsrrrr(mu1,mu4,q1,q2)+
cfc     &   2*q2J3*q4J2*epsrrrr(mu1,mu4,q1,q2)-q3J2*q4J3*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q2q4*epsrrrr(mu1,mu4,q1,q3)+q2J3*q4J2*epsrrrr(mu1,mu4,q1,q3)+
cfc     &   2*q1J2*q4J3*epsrrrr(mu1,mu4,q1,q3)-4*q1J2*q1J3*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*J2J3*q1q2*epsrrrr(mu1,mu4,q1,q4)+J2J3*q1sq*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   4*q1J2*q2J3*epsrrrr(mu1,mu4,q1,q4)+J2J3*q2q3*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   J2J3*q2sq*epsrrrr(mu1,mu4,q1,q4)-q2J3*q3J2*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   J2J3*q1q4*epsrrrr(mu1,mu4,q2,q3)-q1J3*q4J2*epsrrrr(mu1,mu4,q2,q3)+
cfc     &   q1J2*q4J3*epsrrrr(mu1,mu4,q2,q3)-2*q1J2*q1J3*epsrrrr(mu1,mu4,q2,q4)-
cfc     &   J2J3*q1q3*epsrrrr(mu1,mu4,q2,q4)-2*q1J2*q2J3*epsrrrr(mu1,mu4,q2,q4)+
cfc     &   q1J3*q3J2*epsrrrr(mu1,mu4,q2,q4)-2*q1J2*q1J3*epsrrrr(mu1,mu4,q3,q4)+
cfc     &   J2J3*q1q2*epsrrrr(mu1,mu4,q3,q4)+J2J3*q1sq*epsrrrr(mu1,mu4,q3,q4)-
cfc     &   q1J2*q2J3*epsrrrr(mu1,mu4,q3,q4)-2*J2J3*q1J4*epsrrrr(mu1,q1,q2,q3)-
cfc     &   2*J2J3*q2J4*epsrrrr(mu1,q1,q2,q3)-2*J2J3*q3J4*epsrrrr(mu1,q1,q2,q3)-
cfc     &   J3J4*q4J2*epsrrrr(mu1,q1,q2,q3)+J2J4*q4J3*epsrrrr(mu1,q1,q2,q3)-
cfc     &   2*J2J4*q1J3*epsrrrr(mu1,q1,q2,q4)-2*J2J4*q2J3*epsrrrr(mu1,q1,q2,q4)+
cfc     &   J3J4*q3J2*epsrrrr(mu1,q1,q2,q4)-J2J3*q3J4*epsrrrr(mu1,q1,q2,q4)-
cfc     &   2*J3J4*q1J2*epsrrrr(mu1,q1,q3,q4)-J2J4*q2J3*epsrrrr(mu1,q1,q3,q4)+
cfc     &   J2J3*q2J4*epsrrrr(mu1,q1,q3,q4)-J3J4*q1J2*epsrrrr(mu1,q2,q3,q4)+
cfc     &   J2J4*q1J3*epsrrrr(mu1,q2,q3,q4)-J2J3*q1J4*epsrrrr(mu1,q2,q3,q4)+
cfc     &   2*q1q2*q2J1*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   2*q1q3*q2J1*epsrrrr(mu2,mu3,mu4,q1)+q1sq*q2J1*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   2*q2J1*q2q3*epsrrrr(mu2,mu3,mu4,q1)+q2J1*q2sq*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2q4*q3J1*epsrrrr(mu2,mu3,mu4,q1)+q2J1*q3q4*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   q2J1*q3sq*epsrrrr(mu2,mu3,mu4,q1)+2*q1q2*q4J1*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   q1sq*q4J1*epsrrrr(mu2,mu3,mu4,q1)+q2q3*q4J1*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   q2sq*q4J1*epsrrrr(mu2,mu3,mu4,q1)+q1q4*q3J1*epsrrrr(mu2,mu3,mu4,q2)-
cfc     &   q1q3*q4J1*epsrrrr(mu2,mu3,mu4,q2)-q1q4*q2J1*epsrrrr(mu2,mu3,mu4,q3)+
cfc     &   q1q2*q4J1*epsrrrr(mu2,mu3,mu4,q3)+q1sq*q4J1*epsrrrr(mu2,mu3,mu4,q3)+
cfc     &   q1q3*q2J1*epsrrrr(mu2,mu3,mu4,q4)-q1q2*q3J1*epsrrrr(mu2,mu3,mu4,q4)-
cfc     &   q1sq*q3J1*epsrrrr(mu2,mu3,mu4,q4)+2*J1J4*q1q2*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   2*J1J4*q1q3*epsrrrr(mu2,mu3,q1,q2)+J1J4*q1sq*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   2*J1J4*q2q3*epsrrrr(mu2,mu3,q1,q2)+J1J4*q2sq*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   2*q1J4*q3J1*epsrrrr(mu2,mu3,q1,q2)-2*q2J4*q3J1*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   2*q3J1*q3J4*epsrrrr(mu2,mu3,q1,q2)+J1J4*q3q4*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   J1J4*q3sq*epsrrrr(mu2,mu3,q1,q2)-q3J4*q4J1*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   2*q1J4*q2J1*epsrrrr(mu2,mu3,q1,q3)+2*q2J1*q2J4*epsrrrr(mu2,mu3,q1,q3)-
cfc     &   J1J4*q2q4*epsrrrr(mu2,mu3,q1,q3)+2*q2J1*q3J4*epsrrrr(mu2,mu3,q1,q3)+
cfc     &   q2J4*q4J1*epsrrrr(mu2,mu3,q1,q3)+2*J1J4*q1q2*epsrrrr(mu2,mu3,q1,q4)+
cfc     &   J1J4*q1sq*epsrrrr(mu2,mu3,q1,q4)+J1J4*q2q3*epsrrrr(mu2,mu3,q1,q4)+
cfc     &   J1J4*q2sq*epsrrrr(mu2,mu3,q1,q4)-q2J4*q3J1*epsrrrr(mu2,mu3,q1,q4)+
cfc     &   q2J1*q3J4*epsrrrr(mu2,mu3,q1,q4)+J1J4*q1q4*epsrrrr(mu2,mu3,q2,q3)-
cfc     &   q1J4*q4J1*epsrrrr(mu2,mu3,q2,q3)-J1J4*q1q3*epsrrrr(mu2,mu3,q2,q4)+
cfc     &   q1J4*q3J1*epsrrrr(mu2,mu3,q2,q4)+J1J4*q1q2*epsrrrr(mu2,mu3,q3,q4)+
cfc     &   J1J4*q1sq*epsrrrr(mu2,mu3,q3,q4)-q1J4*q2J1*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   2*J1J3*q1q2*epsrrrr(mu2,mu4,q1,q2)-2*J1J3*q1q3*epsrrrr(mu2,mu4,q1,q2)-
cfc     &   J1J3*q1sq*epsrrrr(mu2,mu4,q1,q2)-2*J1J3*q2q3*epsrrrr(mu2,mu4,q1,q2)-
cfc     &   J1J3*q2sq*epsrrrr(mu2,mu4,q1,q2)-J1J3*q3q4*epsrrrr(mu2,mu4,q1,q2)-
cfc     &   J1J3*q3sq*epsrrrr(mu2,mu4,q1,q2)-2*q1J3*q4J1*epsrrrr(mu2,mu4,q1,q2)-
cfc     &   2*q2J3*q4J1*epsrrrr(mu2,mu4,q1,q2)+q3J1*q4J3*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q2q4*epsrrrr(mu2,mu4,q1,q3)-q2J3*q4J1*epsrrrr(mu2,mu4,q1,q3)-
cfc     &   q2J1*q4J3*epsrrrr(mu2,mu4,q1,q3)-2*J1J3*q1q2*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   J1J3*q1sq*epsrrrr(mu2,mu4,q1,q4)+2*q1J3*q2J1*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   2*q2J1*q2J3*epsrrrr(mu2,mu4,q1,q4)-J1J3*q2q3*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   J1J3*q2sq*epsrrrr(mu2,mu4,q1,q4)+q2J3*q3J1*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   J1J3*q1q4*epsrrrr(mu2,mu4,q2,q3)+q1J3*q4J1*epsrrrr(mu2,mu4,q2,q3)+
cfc     &   J1J3*q1q3*epsrrrr(mu2,mu4,q2,q4)-q1J3*q3J1*epsrrrr(mu2,mu4,q2,q4)-
cfc     &   J1J3*q1q2*epsrrrr(mu2,mu4,q3,q4)-J1J3*q1sq*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   q1J3*q2J1*epsrrrr(mu2,mu4,q3,q4)+2*J1J3*q1J4*epsrrrr(mu2,q1,q2,q3)+
cfc     &   2*J1J3*q2J4*epsrrrr(mu2,q1,q2,q3)+2*J1J3*q3J4*epsrrrr(mu2,q1,q2,q3)+
cfc     &   J3J4*q4J1*epsrrrr(mu2,q1,q2,q3)-J1J4*q4J3*epsrrrr(mu2,q1,q2,q3)+
cfc     &   2*J1J4*q1J3*epsrrrr(mu2,q1,q2,q4)+2*J1J4*q2J3*epsrrrr(mu2,q1,q2,q4)-
cfc     &   J3J4*q3J1*epsrrrr(mu2,q1,q2,q4)+J1J3*q3J4*epsrrrr(mu2,q1,q2,q4)+
cfc     &   J3J4*q2J1*epsrrrr(mu2,q1,q3,q4)+J1J4*q2J3*epsrrrr(mu2,q1,q3,q4)-
cfc     &   J1J3*q2J4*epsrrrr(mu2,q1,q3,q4)-J1J4*q1J3*epsrrrr(mu2,q2,q3,q4)+
cfc     &   J1J3*q1J4*epsrrrr(mu2,q2,q3,q4)+2*J1J2*q1q2*epsrrrr(mu3,mu4,q1,q2)+
cfc     &   2*J1J2*q1q3*epsrrrr(mu3,mu4,q1,q2)+J1J2*q1sq*epsrrrr(mu3,mu4,q1,q2)+
cfc     &   2*J1J2*q2q3*epsrrrr(mu3,mu4,q1,q2)+J1J2*q2sq*epsrrrr(mu3,mu4,q1,q2)+
cfc     &   J1J2*q3q4*epsrrrr(mu3,mu4,q1,q2)+J1J2*q3sq*epsrrrr(mu3,mu4,q1,q2)+
cfc     &   q3J2*q4J1*epsrrrr(mu3,mu4,q1,q2)-q3J1*q4J2*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q2q4*epsrrrr(mu3,mu4,q1,q3)-2*q1J2*q4J1*epsrrrr(mu3,mu4,q1,q3)+
cfc     &   q2J1*q4J2*epsrrrr(mu3,mu4,q1,q3)+2*J1J2*q1q2*epsrrrr(mu3,mu4,q1,q4)+
cfc     &   J1J2*q1sq*epsrrrr(mu3,mu4,q1,q4)+J1J2*q2q3*epsrrrr(mu3,mu4,q1,q4)+
cfc     &   J1J2*q2sq*epsrrrr(mu3,mu4,q1,q4)+2*q1J2*q3J1*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   q2J1*q3J2*epsrrrr(mu3,mu4,q1,q4)+J1J2*q1q4*epsrrrr(mu3,mu4,q2,q3)-
cfc     &   q1J2*q4J1*epsrrrr(mu3,mu4,q2,q3)-J1J2*q1q3*epsrrrr(mu3,mu4,q2,q4)+
cfc     &   q1J2*q3J1*epsrrrr(mu3,mu4,q2,q4)+J1J2*q1q2*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   J1J2*q1sq*epsrrrr(mu3,mu4,q3,q4)-q1J2*q2J1*epsrrrr(mu3,mu4,q3,q4)-
cfc     &   2*J1J2*q1J4*epsrrrr(mu3,q1,q2,q3)-2*J1J2*q2J4*epsrrrr(mu3,q1,q2,q3)-
cfc     &   2*J1J2*q3J4*epsrrrr(mu3,q1,q2,q3)-J2J4*q4J1*epsrrrr(mu3,q1,q2,q3)+
cfc     &   J1J4*q4J2*epsrrrr(mu3,q1,q2,q3)+J2J4*q3J1*epsrrrr(mu3,q1,q2,q4)-
cfc     &   J1J4*q3J2*epsrrrr(mu3,q1,q2,q4)-J1J2*q3J4*epsrrrr(mu3,q1,q2,q4)+
cfc     &   2*J1J4*q1J2*epsrrrr(mu3,q1,q3,q4)-J2J4*q2J1*epsrrrr(mu3,q1,q3,q4)+
cfc     &   J1J2*q2J4*epsrrrr(mu3,q1,q3,q4)+J1J4*q1J2*epsrrrr(mu3,q2,q3,q4)-
cfc     &   J1J2*q1J4*epsrrrr(mu3,q2,q3,q4)+J2J3*q4J1*epsrrrr(mu4,q1,q2,q3)-
cfc     &   J1J3*q4J2*epsrrrr(mu4,q1,q2,q3)+J1J2*q4J3*epsrrrr(mu4,q1,q2,q3)-
cfc     &   2*J1J2*q1J3*epsrrrr(mu4,q1,q2,q4)-2*J1J2*q2J3*epsrrrr(mu4,q1,q2,q4)-
cfc     &   J2J3*q3J1*epsrrrr(mu4,q1,q2,q4)+J1J3*q3J2*epsrrrr(mu4,q1,q2,q4)-
cfc     &   2*J1J3*q1J2*epsrrrr(mu4,q1,q3,q4)+J2J3*q2J1*epsrrrr(mu4,q1,q3,q4)-
cfc     &   J1J2*q2J3*epsrrrr(mu4,q1,q3,q4)-J1J3*q1J2*epsrrrr(mu4,q2,q3,q4)+
cfc     &   J1J2*q1J3*epsrrrr(mu4,q2,q3,q4)+
cfc     &   (J1J4*J2J3-J1J3*J2J4+J1J2*J3J4)*epsrrrr(q1,q2,q3,q4))
cfc
cfc      coeffEmu=8*m*(-(q1q2**2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4))-
cfc     &   q1q2*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q1q4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   4*q1q2*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   3*q1q3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1sq**2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1sq*q2q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3q4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   3*q1q2**2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q1q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q1q4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q1sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   3*q1q2*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3q4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   3*q1q2*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3**2*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q1q4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*q1sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q2sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q3q4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   3*q1q2*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4**2*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q4*q1sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q4*q2q3*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q3sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)+
cfc     &   4*q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   3*q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q2*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   q1q3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1sq*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q2q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q2sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q2*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   3*q1q2*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1sq*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q2sq*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q3*q2J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q2*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   3*q2q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   q2J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q4*q2J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   2*q2q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)-
cfc     &   q2J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   q1q3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)-
cfc     &   q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
cfc     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   3*q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1sq*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   3*q1q2*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1sq*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1q3*q2J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q2*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1q4*q2J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   2*q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
cfc     &   q1q2*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1q3*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1sq*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1J3*q2q4*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1q2*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1sq*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1q2*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q2q4*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1q2*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q2sq*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1q3*q2J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q2q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q1q4*q2J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)+
cfc     &   q2J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)-
cfc     &   q1J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1q4*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1q4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
cfc     &   q1sq*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q3)-
cfc     &   q1q2*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q3)-
cfc     &   q1q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q3)-
cfc     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q3)+
cfc     &   2*q1J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   2*q1sq*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   3*q1q2*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(q1sq*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)-
cfc     &   q1J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   2*q1q3*q2J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)-
cfc     &   q1J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   2*q1q4*q2J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)-
cfc     &   q1J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q2*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1sq*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   4*q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J4*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J3*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   4*q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q2J3*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q2J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)-
cfc     &   J3J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q1J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q2J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   2*q3J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J4*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q3)+
cfc     &   q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q3)+
cfc     &   q2J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q3)-
cfc     &   J3J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)+
cfc     &   q1J3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)+
cfc     &   q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   J3J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)+
cfc     &   q2J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)-
cfc     &   q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q2,q3)-
cfc     &   q1J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q2,q3)-
cfc     &   q1J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q2,q3)-
cfc     &   q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q2,q4)-
cfc     &   q1J4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,q2,q4)-
cfc     &   q1J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q2,q4)-
cfc     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q3,q4)-
cfc     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q3,q4)-
cfc     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q3,q4)-
cfc     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q3,q4)+
cfc     &   4*q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   4*q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   3*q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q2q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q2sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q3q4*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q3sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1sq*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2sq*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1q2*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1q3*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1sq*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q2q3*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2sq*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3J2*q3q4*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(q3J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1J2*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1q2*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q1q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   2*q2q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)+
cfc     &   q1J2*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)-
cfc     &   q1J2*q1q4*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q3)+
cfc     &   q1sq*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q3)+
cfc     &   q1q2*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q3)-
cfc     &   q1q4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q3)+
cfc     &   q1q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q3)+
cfc     &   q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   q1sq*q3J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   q1q2*q3J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)-
cfc     &   q1q4*q3J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   q1q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
cfc     &   J2J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   q1J4*q3J2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   q2J4*q3J2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)-
cfc     &   q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   J2J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
cfc     &   4*q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J2*q2J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J2*q3J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J2*q2J4*Eij(1,2)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q2J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q1J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q2J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   2*q3J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
cfc     &   J2J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   q1J2*q3J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(J2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   q3J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4)+
cfc     &   q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q2,q3)+
cfc     &   q1J2*q2J4*Eij(1,2)*epsrrrr(mu1,mu3,q2,q3)+
cfc     &   q1J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q2,q3)-
cfc     &   q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q3,q4)+
cfc     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q3,q4)+
cfc     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q3,q4)+
cfc     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q3,q4)-
cfc     &   q1J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q3,q4)+
cfc     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q3,q4)-
cfc     &   q1J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q3,q4)-
cfc     &   J2J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q1q3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   q1J3*q4J2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   q2J3*q4J2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   J2J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   q4J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
cfc     &   q1J2*q4J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q3)-
cfc     &   q3J2*q4J3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q3)-
cfc     &   q4J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q3)+
cfc     &   4*q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q1J2*q2J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q1J2*q2J3*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q1J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q2J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)-
cfc     &   J2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q1J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q2J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   2*q1J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
cfc     &   q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q2,q4)+
cfc     &   q1J2*q2J3*Eij(1,2)*epsrrrr(mu1,mu4,q2,q4)+
cfc     &   q1J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q2,q4)+
cfc     &   q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q3,q4)-
cfc     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q3,q4)-
cfc     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q3,q4)-
cfc     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q3,q4)+
cfc     &   q1J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q3,q4)-
cfc     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q3,q4)+
cfc     &   q1J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q3,q4)+
cfc     &   J2J3*q1J4*Eij(1,1)*epsrrrr(mu1,q1,q2,q3)+
cfc     &   J2J3*q2J4*Eij(1,2)*epsrrrr(mu1,q1,q2,q3)+
cfc     &   J2J3*q3J4*Eij(1,3)*epsrrrr(mu1,q1,q2,q3)+
cfc     &   J2J4*q1J3*Eij(1,1)*epsrrrr(mu1,q1,q2,q4)+
cfc     &   J2J4*q2J3*Eij(1,2)*epsrrrr(mu1,q1,q2,q4)+
cfc     &   J2J4*q4J3*Eij(1,4)*epsrrrr(mu1,q1,q2,q4)+
cfc     &   J3J4*q1J2*Eij(1,1)*epsrrrr(mu1,q1,q3,q4)+
cfc     &   J3J4*q3J2*Eij(1,3)*epsrrrr(mu1,q1,q3,q4)+
cfc     &   J3J4*q4J2*Eij(1,4)*epsrrrr(mu1,q1,q3,q4)-
cfc     &   q1q2*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q3*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1sq*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q2*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1sq*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q2*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2sq*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q3*q2J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q4*q2J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2J1*q3q4*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q1q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
cfc     &   q2q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)+
cfc     &   2*q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q1q3*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q3q4*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q3sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q1q2*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q1q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q2q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3J1*q3q4*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3J1*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q1q2*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q1q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   2*q2q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)-
cfc     &   q1sq*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q2J1*q2q4*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q1q2*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q2q4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q1q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q1q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q3)-
cfc     &   q2q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q3)+
cfc     &   q1sq*q3J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   2*q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(q1sq*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q3J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   2*q1q2*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q4*q3J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   2*q1q2*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)-
cfc     &   J1J4*q1q2*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q1q3*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   q1J4*q3J1*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   q2J4*q3J1*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)+
cfc     &   q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   J1J4*q3q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
cfc     &   q1J4*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,q1,q3)-
cfc     &   q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q1,q3)-
cfc     &   q2J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q1,q3)-
cfc     &   J1J4*q1q2*Eij(1,1)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q4)-
cfc     &   J1J4*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q4)+
cfc     &   2*q1J4*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q2J1*q3J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q1J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q2J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q1J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q2J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   2*q3J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q2J1*q3J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q4)+
cfc     &   q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q2,q4)+
cfc     &   q3J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q4)-
cfc     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   q2J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q3,q4)-
cfc     &   q2J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q3,q4)+
cfc     &   J1J3*q1q2*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1q3*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   q1J3*q4J1*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q2q3*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q2sq*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   q2J3*q4J1*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(J1J3*q2q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q3sq*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q2q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q3q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   q4J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
cfc     &   J1J3*q1q2*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   q1J3*q2J1*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q2sq*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q2q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)+
cfc     &   J1J3*q2q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   q2J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)-
cfc     &   q2J1*q4J3*Eij(1,2)*epsrrrr(mu2,mu4,q2,q3)-
cfc     &   q3J1*q4J3*Eij(1,3)*epsrrrr(mu2,mu4,q2,q3)-
cfc     &   q4J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q2,q3)+
cfc     &   2*q1J3*q2J1*Eij(1,2)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   2*q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   2*q1J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   2*q2J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   2*q1J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   2*q2J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   q2J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q3,q4)+
cfc     &   q2J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q3,q4)-
cfc     &   J1J3*q1J4*Eij(1,1)*epsrrrr(mu2,q1,q2,q3)-
cfc     &   J1J3*q2J4*Eij(1,2)*epsrrrr(mu2,q1,q2,q3)-
cfc     &   J1J3*q3J4*Eij(1,3)*epsrrrr(mu2,q1,q2,q3)-
cfc     &   J1J4*q1J3*Eij(1,1)*epsrrrr(mu2,q1,q2,q4)-
cfc     &   J1J4*q2J3*Eij(1,2)*epsrrrr(mu2,q1,q2,q4)-
cfc     &   J1J4*q4J3*Eij(1,4)*epsrrrr(mu2,q1,q2,q4)+
cfc     &   J3J4*q2J1*Eij(1,2)*epsrrrr(mu2,q2,q3,q4)+
cfc     &   J3J4*q3J1*Eij(1,3)*epsrrrr(mu2,q2,q3,q4)+
cfc     &   J3J4*q4J1*Eij(1,4)*epsrrrr(mu2,q2,q3,q4)-
cfc     &   J1J2*q1q2*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q1q3*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q2q3*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q2sq*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q2q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q3sq*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q2q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)-
cfc     &   J1J2*q3q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)+
cfc     &   q1J2*q4J1*Eij(1,1)*epsrrrr(mu3,mu4,q1,q3)+
cfc     &   q3J2*q4J1*Eij(1,3)*epsrrrr(mu3,mu4,q1,q3)+
cfc     &   q4J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q1,q3)-
cfc     &   J1J2*q1q2*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   q1J2*q3J1*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q2sq*Eij(1,2)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q2q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   q3J1*q3J2*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   J1J2*q2q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)-
cfc     &   q3J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)+
cfc     &   q2J1*q4J2*Eij(1,2)*epsrrrr(mu3,mu4,q2,q3)+
cfc     &   q3J1*q4J2*Eij(1,3)*epsrrrr(mu3,mu4,q2,q3)+
cfc     &   q4J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q2,q3)-
cfc     &   q2J1*q3J2*Eij(1,2)*epsrrrr(mu3,mu4,q2,q4)-
cfc     &   q3J1*q3J2*Eij(1,3)*epsrrrr(mu3,mu4,q2,q4)-
cfc     &   q3J2*q4J1*Eij(1,4)*epsrrrr(mu3,mu4,q2,q4)-
cfc     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q3,q4)-
cfc     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   2*q1J2*q2J1*Eij(1,2)*epsrrrr(mu3,mu4,q3,q4)-
cfc     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q3,q4))
cfc
cfc      coeffEmu=coeffEmu + 8*m*(2*q1J2*q3J1*Eij(1,3)*epsrrrr(mu3,mu4,q3,q4)-
cfc     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   2*q1J2*q4J1*Eij(1,4)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   J1J2*q1J4*Eij(1,1)*epsrrrr(mu3,q1,q2,q3)+
cfc     &   J1J2*q2J4*Eij(1,2)*epsrrrr(mu3,q1,q2,q3)+
cfc     &   J1J2*q3J4*Eij(1,3)*epsrrrr(mu3,q1,q2,q3)-
cfc     &   J1J4*(q1J2*Eij(1,1)+q3J2*Eij(1,3)+q4J2*Eij(1,4))*
cfc     &    epsrrrr(mu3,q1,q3,q4)-J2J4*q2J1*Eij(1,2)*epsrrrr(mu3,q2,q3,q4)-
cfc     &   J2J4*q3J1*Eij(1,3)*epsrrrr(mu3,q2,q3,q4)-
cfc     &   J2J4*q4J1*Eij(1,4)*epsrrrr(mu3,q2,q3,q4)+
cfc     &   J1J2*q1J3*Eij(1,1)*epsrrrr(mu4,q1,q2,q4)+
cfc     &   J1J2*q2J3*Eij(1,2)*epsrrrr(mu4,q1,q2,q4)+
cfc     &   J1J2*q4J3*Eij(1,4)*epsrrrr(mu4,q1,q2,q4)+
cfc     &   J1J3*(q1J2*Eij(1,1)+q3J2*Eij(1,3)+q4J2*Eij(1,4))*
cfc     &    epsrrrr(mu4,q1,q3,q4)+J2J3*q2J1*Eij(1,2)*epsrrrr(mu4,q2,q3,q4)+
cfc     &   J2J3*q3J1*Eij(1,3)*epsrrrr(mu4,q2,q3,q4)+
cfc     &   J2J3*q4J1*Eij(1,4)*epsrrrr(mu4,q2,q3,q4))
cfc
cfc      coeffEmunu=16*m*(-(q1q2*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4))-
cfc     &   q1q3*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq**2*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2**2*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2q3*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2sq*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3**2*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q2q3*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q3sq*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4**2*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2q4*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q3q4*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2**2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q1q3*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q1sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q3*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3**2*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*q1sq*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3sq*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q4*q1sq*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q2q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*q3q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q1q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q2sq*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3sq*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q2*q1q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2q3*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q2q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2sq*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*q3q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   2*q1q3*q1q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q2q3*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q2q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*q3q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q4*q3sq*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q2*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1q3*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)-
cfc     &   q1sq*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)+
cfc     &   q1J4*q1q2*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q2*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q2sq*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q2*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1sq*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q2sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q2q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q2*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1sq*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q2q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q3*q2J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q2*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2sq*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q4*q2J4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*q2q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q2J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q1)+
cfc     &   q1J4*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q2*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J4*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1sq*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J4*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1sq*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J4*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q3*q2J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q2*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q4*q2J4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q3)+
cfc     &   q1J3*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q2*q2J3*Eij(2,2)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q4*q4J3*Eij(2,4)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1sq*q2J3*Eij(2,5)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1sq*q4J3*Eij(2,7)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q3*q2J3*Eij(2,8)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q4*q2J3*Eij(2,9)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q2*q4J3*Eij(2,9)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1q3*q4J3*Eij(2,10)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*Eij(2,11)*epsrrrr(mu1,mu2,mu4,q4)+
cfc     &   q1J3*q1J4*Eij(2,1)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q2J3*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J4*q2J3*Eij(2,5)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J3*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J3*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J4*q4J3*Eij(2,7)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q2J3*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q2J4*q4J3*Eij(2,9)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q3J4*q4J3*Eij(2,10)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   J3J4*Eij(2,11)*epsrrrr(mu1,mu2,q1,q2)+
cfc     &   q1J2*q1q2*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1q3*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1sq*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q3*q3J2*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q3*q3J2*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3J2*q3sq*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1q2*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q2q3*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q2sq*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1q3*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q2q3*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q3*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1sq*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q3sq*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q1q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q2q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*q3q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q3*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1sq*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q3*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2sq*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q2*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q3*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2sq*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q4*q3J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q4*q3J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3J2*q3q4*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1q3*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q2q3*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3sq*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q1J2*Eij(2,11)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   q3J2*Eij(2,11)*epsrrrr(mu1,mu3,mu4,q1)+
cfc     &   (q3J2*q3J4*Eij(2,3)+q1J4*q3J2*Eij(2,6)+
cfc     &      q1J2*(q1J4*Eij(2,1)+q2J4*Eij(2,5)+q3J4*Eij(2,6))+
cfc     &      q2J4*q3J2*Eij(2,8)+q4J2*
cfc     &       (q1J4*Eij(2,7)+q2J4*Eij(2,9)+q3J4*Eij(2,10))+J2J4*Eij(2,11))*
cfc     &    epsrrrr(mu1,mu3,q1,q3)+
cfc     &   (q4J2*q4J3*Eij(2,4)+q1J3*q4J2*Eij(2,7)+
cfc     &      q1J2*(q1J3*Eij(2,1)+q2J3*Eij(2,5)+q4J3*Eij(2,7))+
cfc     &      q2J3*q4J2*Eij(2,9)+q3J2*
cfc     &       (q1J3*Eij(2,6)+q2J3*Eij(2,8)+q4J3*Eij(2,10))+J2J3*Eij(2,11))*
cfc     &    epsrrrr(mu1,mu4,q1,q4)+q1q2*q2J1*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q2q3*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q2sq*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3J1*q3sq*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1sq*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q2J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q2q3*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q3*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2sq*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q3sq*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q4*q2J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q2q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*q3q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q3*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2sq*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3J1*q3q4*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3sq*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q2J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q3J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q2)+
cfc     &   q1q2*q2J1*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2sq*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1sq*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1sq*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1sq*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q3*q2J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2q3*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2sq*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q4*q2J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q2*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2sq*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q1q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q4)+
cfc     &   q2J1*q2J4*Eij(2,2)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q3J1*q3J4*Eij(2,3)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q1J4*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q1J4*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q1J4*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q2J4*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q2J1*q3J4*Eij(2,8)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q2J4*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q3J4*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   J1J4*Eij(2,11)*epsrrrr(mu2,mu3,q2,q3)+
cfc     &   q2J1*q2J3*Eij(2,2)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q4J1*q4J3*Eij(2,4)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q1J3*q2J1*Eij(2,5)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q1J3*q3J1*Eij(2,6)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q1J3*q4J1*Eij(2,7)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q2J3*q3J1*Eij(2,8)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q2J3*q4J1*Eij(2,9)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q2J1*q4J3*Eij(2,9)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q3J1*q4J3*Eij(2,10)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   J1J3*Eij(2,11)*epsrrrr(mu2,mu4,q2,q4)+
cfc     &   q3J1*q3J2*Eij(2,3)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q4J1*q4J2*Eij(2,4)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q1J2*q2J1*Eij(2,5)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q1J2*q3J1*Eij(2,6)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q1J2*q4J1*Eij(2,7)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q2J1*q3J2*Eij(2,8)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q2J1*q4J2*Eij(2,9)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q3J2*q4J1*Eij(2,10)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   q3J1*q4J2*Eij(2,10)*epsrrrr(mu3,mu4,q3,q4)+
cfc     &   J1J2*Eij(2,11)*epsrrrr(mu3,mu4,q3,q4))
cfc
cfc      A4G_LoopDD=(coeffD0+coeffDmu+coeffE0C0+coeffE0+coeffEmu
cfc     &	         +coeffEmunu)/(4d0*m)
       pensingular(2)=.false.
       pensingular(1)=.false.
     

      do mu=0,3
         Kpen(mu,i1)=K(mu,i1)
         Kpen(mu,i2)=K(mu,i2)
         Kpen(mu,i3)=K(mu,i3)
         Kpen(mu,i4)=K(mu,i4)
         PA1(mu)=-(K(mu,i1)+K(mu,i2)+K(mu,i3)+K(mu,i4))
         if(mu.eq.4) exit
         epsgPen(mu,i1)=epvec(mu,i1)
         epsgPen(mu,i2)=epvec(mu,i2)
         epsgPen(mu,i3)=epvec(mu,i3) 
         epsgPen(mu,i4)=epvec(mu,i4) 
      enddo
         
        CALL gluonpenHiggsOdd_T(m,Kpen(0,i1),Kpen(0,i2),
     &        Kpen(0,i3),Kpen(0,i4),PA1,epsgPen(0,i1),epsgPen(0,i2),
     &     epsgPen(0,i3),epsgPen(0,i4),pen1)
        
        pen1=-pen1/(12d0*m*m)
        
        if (pensingular(1)) then
           PendbleSin=.true.
           pen1=0d0
        endif
#ifdef WITH_QUAD
        if(QUAD.or.penSingular(1)) then
        CALL gluonpenHiggsOdd_QUADT(m,Kpen(0,i1),Kpen(0,i2),
     &        Kpen(0,i3),Kpen(0,i4),PA1,epsgPen(0,i1),epsgPen(0,i2),
     &          epsgPen(0,i3),epsgPen(0,i4),pen1)

      pen1=-pen1/(12d0*m*m)
      
        if (pensingular(2)) then
           PenQuadSin=.true.
           pen1=0d0
        endif
      endif
#endif        

        A4G_LoopDD=pen1
        
cfc        if(ABS(1D0-(pen1)/A4G_LoopDD).gt.1d-3) then
cfc         print*, i1,i2,i3,i4,m
cfcc         print*, "pen1",pen1
cfcc         PRINT*, "A4G_LoopDD",A4G_LoopDD
cfc         print*, "PEn/PEN_OLD",ABS(1D0-(pen1)/A4G_LoopDD)
cfcc     fc         PAUSE
cfc         endif
      RETURN
      END
