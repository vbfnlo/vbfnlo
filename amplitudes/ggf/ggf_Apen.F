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
     &        mu1(0:3),mu2(0:3),mu3(0:3),mu4(0:3)

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

      if (init.eq.0) then   
         init = 1
         ie0(4,2,1,3) = 1   ! this is the map from i1,i2,i3,i4 
         ie0(2,1,3,4) = 2   ! to the index of the E0 variables
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

      if (m.eq.m_top) then
         E0L=-E0t(ie0(i1,i2,i3,i4))
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijtop(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo  
      endif
      if (m.eq.m_bottom) then
         E0L =-E0b(ie0(i1,i2,i3,i4))
         do i=1,4
            do j=1,46 
               Eij(i,j)=DDEijbottom(i,j,ie0(i1,i2,i3,i4))
            enddo
         enddo   
      endif

      do mu = 0,3
         q1(mu) = qi(mu,i1)
         q2(mu) = qi(mu,i2)
         q3(mu) = qi(mu,i3)
         q4(mu) = qi(mu,i4)

         mu1(mu) = eps(mu,i1)
         mu2(mu) = eps(mu,i2)
         mu3(mu) = eps(mu,i3)
         mu4(mu) = eps(mu,i4)
      enddo

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

      q1sq=q_dot_q(i1,i1)
      q2sq=q_dot_q(i2,i2)
      q3sq=q_dot_q(i3,i3)
      q4sq=q_dot_q(i4,i4)

      q12sq=2d0*q1q2
      q34sq=2d0*q3q4
      q1234sq=2*(q1q2+q1q3+q1q4+q2q3+q2q4+q3q4)

      C034=C0_1M(m,q3sq,q4sq,q34sq) 

      D0234 = D0_1M(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4)
      CALL D_1M_ALL2(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4,Dij)

      coeffD0=-4*m*((q1q2+q1q4-2*q1sq+2*q2q3+q2sq+q3q4+q3sq)*
     &    epsrrrr(mu1,mu2,mu3,mu4)-(q2J4+2*q3J4)*epsrrrr(mu1,mu2,mu3,q1)-
     &   q1J4*epsrrrr(mu1,mu2,mu3,q2)-2*(q2J4+q3J4)*epsrrrr(mu1,mu2,mu3,q3)-
     &   (q1J4+q3J4)*epsrrrr(mu1,mu2,mu3,q4)+
     &   (-q2J3+q4J3)*epsrrrr(mu1,mu2,mu4,q1)+q1J3*epsrrrr(mu1,mu2,mu4,q2)+
     &   q4J3*epsrrrr(mu1,mu2,mu4,q3)+(q1J3-2*q2J3)*epsrrrr(mu1,mu2,mu4,q4)+
     &   J3J4*(epsrrrr(mu1,mu2,q1,q2)+epsrrrr(mu1,mu2,q1,q4)+
     &      epsrrrr(mu1,mu2,q3,q4))-q1J2*epsrrrr(mu1,mu3,mu4,q2)-
     &   q4J2*(epsrrrr(mu1,mu3,mu4,q1)+epsrrrr(mu1,mu3,mu4,q3))+
     &   (-q1J2+q3J2)*epsrrrr(mu1,mu3,mu4,q4)-
     &   J2J4*(epsrrrr(mu1,mu3,q1,q2)+epsrrrr(mu1,mu3,q1,q4)+
     &      epsrrrr(mu1,mu3,q3,q4))+
     &   J2J3*(epsrrrr(mu1,mu4,q1,q2)+epsrrrr(mu1,mu4,q1,q4)+
     &      epsrrrr(mu1,mu4,q3,q4))+(q2J1+q4J1)*epsrrrr(mu2,mu3,mu4,q1)+
     &   q4J1*epsrrrr(mu2,mu3,mu4,q3)-q3J1*epsrrrr(mu2,mu3,mu4,q4)+
     &   J1J4*(epsrrrr(mu2,mu3,q1,q2)+epsrrrr(mu2,mu3,q1,q4)+
     &      epsrrrr(mu2,mu3,q3,q4))-
     &   J1J3*(epsrrrr(mu2,mu4,q1,q2)+epsrrrr(mu2,mu4,q1,q4)+
     &      epsrrrr(mu2,mu4,q3,q4))+
     &   J1J2*(epsrrrr(mu3,mu4,q1,q2)+epsrrrr(mu3,mu4,q1,q4)+
     &      epsrrrr(mu3,mu4,q3,q4)))*D0234

      coeffDmu=8*m*(-2*q1q2*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q2q3*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q2sq*Dij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q2q3*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q3sq*Dij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q2q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q3q4*Dij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)+
     &   q2J4*Dij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q3J4*Dij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*Dij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q3J4*Dij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q2J3*Dij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q4J3*Dij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q3J2*Dij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q4J2*Dij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   (q2J1*Dij(1,1)+q3J1*Dij(1,2)+q4J1*Dij(1,3))*epsrrrr(mu2,mu3,mu4,q2)+
     &   (q2J1*Dij(1,1)+q3J1*Dij(1,2)+q4J1*Dij(1,3))*epsrrrr(mu2,mu3,mu4,q4))

      coeffE0C0=4*m*(-C034+E0L*m**4)*epsrrrr(mu1,mu2,mu3,mu4)

      coeffE0=-4*E0L*m*((m**4+2*q1q2**2+q1sq**2+q1q4*q2q3-q1q3*q2q4+q1q4*q2sq+
     &      q1sq*(3*q1q2+2*q1q3+q1q4+2*q2q3+q2sq+q3q4+q3sq)+
     &      q1q2*(2*q1q3+2*q1q4+2*q2q3+q2sq+q3q4+q3sq))*
     &    epsrrrr(mu1,mu2,mu3,mu4)-4*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q1q2*q2J4*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q1)-q1sq*q2J4*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q1)-q2J4*q2sq*epsrrrr(mu1,mu2,mu3,q1)-
     &   4*q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q2q3*q3J4*epsrrrr(mu1,mu2,mu3,q1)-q2q4*q3J4*epsrrrr(mu1,mu2,mu3,q1)-
     &   2*q2sq*q3J4*epsrrrr(mu1,mu2,mu3,q1)+q2J4*q3q4*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q3sq*epsrrrr(mu1,mu2,mu3,q1)-2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q2)-
     &   q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q2)+2*q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q2)-
     &   2*q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q2)-q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q2)+
     &   2*q1q3*q3J4*epsrrrr(mu1,mu2,mu3,q2)+q1q4*q3J4*epsrrrr(mu1,mu2,mu3,q2)-
     &   q1J4*q3q4*epsrrrr(mu1,mu2,mu3,q2)-q1J4*q3sq*epsrrrr(mu1,mu2,mu3,q2)-
     &   2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1q2*q2J4*epsrrrr(mu1,mu2,mu3,q3)-q1q4*q2J4*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1sq*q2J4*epsrrrr(mu1,mu2,mu3,q3)+q1J4*q2q4*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q3)-
     &   2*q1J4*q1q2*epsrrrr(mu1,mu2,mu3,q4)-q1J4*q1sq*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1q3*q2J4*epsrrrr(mu1,mu2,mu3,q4)-q1J4*q2q3*epsrrrr(mu1,mu2,mu3,q4)-
     &   q1J4*q2sq*epsrrrr(mu1,mu2,mu3,q4)-q1q2*q3J4*epsrrrr(mu1,mu2,mu3,q4)-
     &   q1sq*q3J4*epsrrrr(mu1,mu2,mu3,q4)-2*q1q2*q2J3*epsrrrr(mu1,mu2,mu4,q1)-
     &   2*q1q3*q2J3*epsrrrr(mu1,mu2,mu4,q1)-q1sq*q2J3*epsrrrr(mu1,mu2,mu4,q1)-
     &   2*q2J3*q2q3*epsrrrr(mu1,mu2,mu4,q1)-
     &   2*q1J3*q2q4*epsrrrr(mu1,mu2,mu4,q1)-
     &   2*q2J3*q2q4*epsrrrr(mu1,mu2,mu4,q1)-q2J3*q2sq*epsrrrr(mu1,mu2,mu4,q1)-
     &   q2J3*q3q4*epsrrrr(mu1,mu2,mu4,q1)-q2J3*q3sq*epsrrrr(mu1,mu2,mu4,q1)+
     &   2*q1q2*q4J3*epsrrrr(mu1,mu2,mu4,q1)+q1sq*q4J3*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2q3*q4J3*epsrrrr(mu1,mu2,mu4,q1)+q2sq*q4J3*epsrrrr(mu1,mu2,mu4,q1)+
     &   2*q1J3*q1q2*epsrrrr(mu1,mu2,mu4,q2)+
     &   2*q1J3*q1q3*epsrrrr(mu1,mu2,mu4,q2)+
     &   2*q1J3*q1q4*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q1sq*epsrrrr(mu1,mu2,mu4,q2)+
     &   2*q1q4*q2J3*epsrrrr(mu1,mu2,mu4,q2)+
     &   2*q1J3*q2q3*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q2sq*epsrrrr(mu1,mu2,mu4,q2)+
     &   q1J3*q3q4*epsrrrr(mu1,mu2,mu4,q2)+q1J3*q3sq*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1q3*q4J3*epsrrrr(mu1,mu2,mu4,q2)+q1q4*q2J3*epsrrrr(mu1,mu2,mu4,q3)-
     &   q1J3*q2q4*epsrrrr(mu1,mu2,mu4,q3)+q1q2*q4J3*epsrrrr(mu1,mu2,mu4,q3)+
     &   q1sq*q4J3*epsrrrr(mu1,mu2,mu4,q3)-q1J3*q1sq*epsrrrr(mu1,mu2,mu4,q4)-
     &   2*q1q2*q2J3*epsrrrr(mu1,mu2,mu4,q4)-q1q3*q2J3*epsrrrr(mu1,mu2,mu4,q4)-
     &   2*q1sq*q2J3*epsrrrr(mu1,mu2,mu4,q4)+q1J3*q2q3*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q2sq*epsrrrr(mu1,mu2,mu4,q4)-4*q1J3*q1J4*epsrrrr(mu1,mu2,q1,q2)+
     &   2*J3J4*q1q2*epsrrrr(mu1,mu2,q1,q2)+2*J3J4*q1q3*epsrrrr(mu1,mu2,q1,q2)+
     &   J3J4*q1sq*epsrrrr(mu1,mu2,q1,q2)-4*q1J4*q2J3*epsrrrr(mu1,mu2,q1,q2)-
     &   4*q1J3*q2J4*epsrrrr(mu1,mu2,q1,q2)-4*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q2)+
     &   2*J3J4*q2q3*epsrrrr(mu1,mu2,q1,q2)+J3J4*q2sq*epsrrrr(mu1,mu2,q1,q2)-
     &   4*q1J3*q3J4*epsrrrr(mu1,mu2,q1,q2)-4*q2J3*q3J4*epsrrrr(mu1,mu2,q1,q2)+
     &   J3J4*q3q4*epsrrrr(mu1,mu2,q1,q2)+J3J4*q3sq*epsrrrr(mu1,mu2,q1,q2)-
     &   q3J4*q4J3*epsrrrr(mu1,mu2,q1,q2)-2*q1J4*q2J3*epsrrrr(mu1,mu2,q1,q3)-
     &   2*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q3)-J3J4*q2q4*epsrrrr(mu1,mu2,q1,q3)-
     &   2*q2J3*q3J4*epsrrrr(mu1,mu2,q1,q3)+q2J4*q4J3*epsrrrr(mu1,mu2,q1,q3)+
     &   2*J3J4*q1q2*epsrrrr(mu1,mu2,q1,q4)+J3J4*q1sq*epsrrrr(mu1,mu2,q1,q4)-
     &   2*q1J3*q2J4*epsrrrr(mu1,mu2,q1,q4)-2*q2J3*q2J4*epsrrrr(mu1,mu2,q1,q4)+
     &   J3J4*q2q3*epsrrrr(mu1,mu2,q1,q4)+J3J4*q2sq*epsrrrr(mu1,mu2,q1,q4)-
     &   q2J3*q3J4*epsrrrr(mu1,mu2,q1,q4)+2*q1J3*q1J4*epsrrrr(mu1,mu2,q2,q3)+
     &   J3J4*q1q4*epsrrrr(mu1,mu2,q2,q3)+2*q1J3*q2J4*epsrrrr(mu1,mu2,q2,q3)+
     &   2*q1J3*q3J4*epsrrrr(mu1,mu2,q2,q3)-q1J4*q4J3*epsrrrr(mu1,mu2,q2,q3)+
     &   2*q1J3*q1J4*epsrrrr(mu1,mu2,q2,q4)-J3J4*q1q3*epsrrrr(mu1,mu2,q2,q4)+
     &   2*q1J4*q2J3*epsrrrr(mu1,mu2,q2,q4)+q1J3*q3J4*epsrrrr(mu1,mu2,q2,q4)+
     &   J3J4*q1q2*epsrrrr(mu1,mu2,q3,q4)+J3J4*q1sq*epsrrrr(mu1,mu2,q3,q4)+
     &   q1J4*q2J3*epsrrrr(mu1,mu2,q3,q4)-q1J3*q2J4*epsrrrr(mu1,mu2,q3,q4)-
     &   4*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q1)-
     &   4*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q1)-
     &   4*q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q1)+q2q4*q3J2*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1J2*q3q4*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1J2*q3sq*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1q2*q4J2*epsrrrr(mu1,mu3,mu4,q1)-q1sq*q4J2*epsrrrr(mu1,mu3,mu4,q1)-
     &   q2q3*q4J2*epsrrrr(mu1,mu3,mu4,q1)-q2sq*q4J2*epsrrrr(mu1,mu3,mu4,q1)-
     &   2*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q2)-
     &   2*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q2)-
     &   2*q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q2)-
     &   q1q4*q3J2*epsrrrr(mu1,mu3,mu4,q2)-q1J2*q3q4*epsrrrr(mu1,mu3,mu4,q2)-
     &   q1J2*q3sq*epsrrrr(mu1,mu3,mu4,q2)+q1q3*q4J2*epsrrrr(mu1,mu3,mu4,q2)+
     &   2*q1J2*q1q4*epsrrrr(mu1,mu3,mu4,q3)+q1J2*q2q4*epsrrrr(mu1,mu3,mu4,q3)-
     &   q1q2*q4J2*epsrrrr(mu1,mu3,mu4,q3)-q1sq*q4J2*epsrrrr(mu1,mu3,mu4,q3)-
     &   2*q1J2*q1q2*epsrrrr(mu1,mu3,mu4,q4)-
     &   2*q1J2*q1q3*epsrrrr(mu1,mu3,mu4,q4)-q1J2*q1sq*epsrrrr(mu1,mu3,mu4,q4)-
     &   q1J2*q2q3*epsrrrr(mu1,mu3,mu4,q4)-q1J2*q2sq*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1q2*q3J2*epsrrrr(mu1,mu3,mu4,q4)+q1sq*q3J2*epsrrrr(mu1,mu3,mu4,q4)-
     &   2*J2J4*q1q2*epsrrrr(mu1,mu3,q1,q2)-2*J2J4*q1q3*epsrrrr(mu1,mu3,q1,q2)-
     &   J2J4*q1sq*epsrrrr(mu1,mu3,q1,q2)-2*J2J4*q2q3*epsrrrr(mu1,mu3,q1,q2)-
     &   J2J4*q2sq*epsrrrr(mu1,mu3,q1,q2)+2*q1J4*q3J2*epsrrrr(mu1,mu3,q1,q2)+
     &   2*q2J4*q3J2*epsrrrr(mu1,mu3,q1,q2)+2*q3J2*q3J4*epsrrrr(mu1,mu3,q1,q2)-
     &   J2J4*q3q4*epsrrrr(mu1,mu3,q1,q2)-J2J4*q3sq*epsrrrr(mu1,mu3,q1,q2)+
     &   q3J4*q4J2*epsrrrr(mu1,mu3,q1,q2)-4*q1J2*q1J4*epsrrrr(mu1,mu3,q1,q3)-
     &   4*q1J2*q2J4*epsrrrr(mu1,mu3,q1,q3)+J2J4*q2q4*epsrrrr(mu1,mu3,q1,q3)-
     &   4*q1J2*q3J4*epsrrrr(mu1,mu3,q1,q3)-q2J4*q4J2*epsrrrr(mu1,mu3,q1,q3)-
     &   2*J2J4*q1q2*epsrrrr(mu1,mu3,q1,q4)-J2J4*q1sq*epsrrrr(mu1,mu3,q1,q4)-
     &   J2J4*q2q3*epsrrrr(mu1,mu3,q1,q4)-J2J4*q2sq*epsrrrr(mu1,mu3,q1,q4)+
     &   q2J4*q3J2*epsrrrr(mu1,mu3,q1,q4)-2*q1J2*q3J4*epsrrrr(mu1,mu3,q1,q4)-
     &   2*q1J2*q1J4*epsrrrr(mu1,mu3,q2,q3)-J2J4*q1q4*epsrrrr(mu1,mu3,q2,q3)-
     &   2*q1J2*q2J4*epsrrrr(mu1,mu3,q2,q3)-2*q1J2*q3J4*epsrrrr(mu1,mu3,q2,q3)+
     &   q1J4*q4J2*epsrrrr(mu1,mu3,q2,q3)+J2J4*q1q3*epsrrrr(mu1,mu3,q2,q4)-
     &   q1J4*q3J2*epsrrrr(mu1,mu3,q2,q4)-q1J2*q3J4*epsrrrr(mu1,mu3,q2,q4)+
     &   2*q1J2*q1J4*epsrrrr(mu1,mu3,q3,q4)-J2J4*q1q2*epsrrrr(mu1,mu3,q3,q4)-
     &   J2J4*q1sq*epsrrrr(mu1,mu3,q3,q4)+q1J2*q2J4*epsrrrr(mu1,mu3,q3,q4)+
     &   2*J2J3*q1q2*epsrrrr(mu1,mu4,q1,q2)+2*J2J3*q1q3*epsrrrr(mu1,mu4,q1,q2)+
     &   J2J3*q1sq*epsrrrr(mu1,mu4,q1,q2)+2*J2J3*q2q3*epsrrrr(mu1,mu4,q1,q2)+
     &   J2J3*q2sq*epsrrrr(mu1,mu4,q1,q2)+J2J3*q3q4*epsrrrr(mu1,mu4,q1,q2)+
     &   J2J3*q3sq*epsrrrr(mu1,mu4,q1,q2)+2*q1J3*q4J2*epsrrrr(mu1,mu4,q1,q2)+
     &   2*q2J3*q4J2*epsrrrr(mu1,mu4,q1,q2)-q3J2*q4J3*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q2q4*epsrrrr(mu1,mu4,q1,q3)+q2J3*q4J2*epsrrrr(mu1,mu4,q1,q3)+
     &   2*q1J2*q4J3*epsrrrr(mu1,mu4,q1,q3)-4*q1J2*q1J3*epsrrrr(mu1,mu4,q1,q4)+
     &   2*J2J3*q1q2*epsrrrr(mu1,mu4,q1,q4)+J2J3*q1sq*epsrrrr(mu1,mu4,q1,q4)-
     &   4*q1J2*q2J3*epsrrrr(mu1,mu4,q1,q4)+J2J3*q2q3*epsrrrr(mu1,mu4,q1,q4)+
     &   J2J3*q2sq*epsrrrr(mu1,mu4,q1,q4)-q2J3*q3J2*epsrrrr(mu1,mu4,q1,q4)+
     &   J2J3*q1q4*epsrrrr(mu1,mu4,q2,q3)-q1J3*q4J2*epsrrrr(mu1,mu4,q2,q3)+
     &   q1J2*q4J3*epsrrrr(mu1,mu4,q2,q3)-2*q1J2*q1J3*epsrrrr(mu1,mu4,q2,q4)-
     &   J2J3*q1q3*epsrrrr(mu1,mu4,q2,q4)-2*q1J2*q2J3*epsrrrr(mu1,mu4,q2,q4)+
     &   q1J3*q3J2*epsrrrr(mu1,mu4,q2,q4)-2*q1J2*q1J3*epsrrrr(mu1,mu4,q3,q4)+
     &   J2J3*q1q2*epsrrrr(mu1,mu4,q3,q4)+J2J3*q1sq*epsrrrr(mu1,mu4,q3,q4)-
     &   q1J2*q2J3*epsrrrr(mu1,mu4,q3,q4)-2*J2J3*q1J4*epsrrrr(mu1,q1,q2,q3)-
     &   2*J2J3*q2J4*epsrrrr(mu1,q1,q2,q3)-2*J2J3*q3J4*epsrrrr(mu1,q1,q2,q3)-
     &   J3J4*q4J2*epsrrrr(mu1,q1,q2,q3)+J2J4*q4J3*epsrrrr(mu1,q1,q2,q3)-
     &   2*J2J4*q1J3*epsrrrr(mu1,q1,q2,q4)-2*J2J4*q2J3*epsrrrr(mu1,q1,q2,q4)+
     &   J3J4*q3J2*epsrrrr(mu1,q1,q2,q4)-J2J3*q3J4*epsrrrr(mu1,q1,q2,q4)-
     &   2*J3J4*q1J2*epsrrrr(mu1,q1,q3,q4)-J2J4*q2J3*epsrrrr(mu1,q1,q3,q4)+
     &   J2J3*q2J4*epsrrrr(mu1,q1,q3,q4)-J3J4*q1J2*epsrrrr(mu1,q2,q3,q4)+
     &   J2J4*q1J3*epsrrrr(mu1,q2,q3,q4)-J2J3*q1J4*epsrrrr(mu1,q2,q3,q4)+
     &   2*q1q2*q2J1*epsrrrr(mu2,mu3,mu4,q1)+
     &   2*q1q3*q2J1*epsrrrr(mu2,mu3,mu4,q1)+q1sq*q2J1*epsrrrr(mu2,mu3,mu4,q1)+
     &   2*q2J1*q2q3*epsrrrr(mu2,mu3,mu4,q1)+q2J1*q2sq*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2q4*q3J1*epsrrrr(mu2,mu3,mu4,q1)+q2J1*q3q4*epsrrrr(mu2,mu3,mu4,q1)+
     &   q2J1*q3sq*epsrrrr(mu2,mu3,mu4,q1)+2*q1q2*q4J1*epsrrrr(mu2,mu3,mu4,q1)+
     &   q1sq*q4J1*epsrrrr(mu2,mu3,mu4,q1)+q2q3*q4J1*epsrrrr(mu2,mu3,mu4,q1)+
     &   q2sq*q4J1*epsrrrr(mu2,mu3,mu4,q1)+q1q4*q3J1*epsrrrr(mu2,mu3,mu4,q2)-
     &   q1q3*q4J1*epsrrrr(mu2,mu3,mu4,q2)-q1q4*q2J1*epsrrrr(mu2,mu3,mu4,q3)+
     &   q1q2*q4J1*epsrrrr(mu2,mu3,mu4,q3)+q1sq*q4J1*epsrrrr(mu2,mu3,mu4,q3)+
     &   q1q3*q2J1*epsrrrr(mu2,mu3,mu4,q4)-q1q2*q3J1*epsrrrr(mu2,mu3,mu4,q4)-
     &   q1sq*q3J1*epsrrrr(mu2,mu3,mu4,q4)+2*J1J4*q1q2*epsrrrr(mu2,mu3,q1,q2)+
     &   2*J1J4*q1q3*epsrrrr(mu2,mu3,q1,q2)+J1J4*q1sq*epsrrrr(mu2,mu3,q1,q2)+
     &   2*J1J4*q2q3*epsrrrr(mu2,mu3,q1,q2)+J1J4*q2sq*epsrrrr(mu2,mu3,q1,q2)-
     &   2*q1J4*q3J1*epsrrrr(mu2,mu3,q1,q2)-2*q2J4*q3J1*epsrrrr(mu2,mu3,q1,q2)-
     &   2*q3J1*q3J4*epsrrrr(mu2,mu3,q1,q2)+J1J4*q3q4*epsrrrr(mu2,mu3,q1,q2)+
     &   J1J4*q3sq*epsrrrr(mu2,mu3,q1,q2)-q3J4*q4J1*epsrrrr(mu2,mu3,q1,q2)+
     &   2*q1J4*q2J1*epsrrrr(mu2,mu3,q1,q3)+2*q2J1*q2J4*epsrrrr(mu2,mu3,q1,q3)-
     &   J1J4*q2q4*epsrrrr(mu2,mu3,q1,q3)+2*q2J1*q3J4*epsrrrr(mu2,mu3,q1,q3)+
     &   q2J4*q4J1*epsrrrr(mu2,mu3,q1,q3)+2*J1J4*q1q2*epsrrrr(mu2,mu3,q1,q4)+
     &   J1J4*q1sq*epsrrrr(mu2,mu3,q1,q4)+J1J4*q2q3*epsrrrr(mu2,mu3,q1,q4)+
     &   J1J4*q2sq*epsrrrr(mu2,mu3,q1,q4)-q2J4*q3J1*epsrrrr(mu2,mu3,q1,q4)+
     &   q2J1*q3J4*epsrrrr(mu2,mu3,q1,q4)+J1J4*q1q4*epsrrrr(mu2,mu3,q2,q3)-
     &   q1J4*q4J1*epsrrrr(mu2,mu3,q2,q3)-J1J4*q1q3*epsrrrr(mu2,mu3,q2,q4)+
     &   q1J4*q3J1*epsrrrr(mu2,mu3,q2,q4)+J1J4*q1q2*epsrrrr(mu2,mu3,q3,q4)+
     &   J1J4*q1sq*epsrrrr(mu2,mu3,q3,q4)-q1J4*q2J1*epsrrrr(mu2,mu3,q3,q4)-
     &   2*J1J3*q1q2*epsrrrr(mu2,mu4,q1,q2)-2*J1J3*q1q3*epsrrrr(mu2,mu4,q1,q2)-
     &   J1J3*q1sq*epsrrrr(mu2,mu4,q1,q2)-2*J1J3*q2q3*epsrrrr(mu2,mu4,q1,q2)-
     &   J1J3*q2sq*epsrrrr(mu2,mu4,q1,q2)-J1J3*q3q4*epsrrrr(mu2,mu4,q1,q2)-
     &   J1J3*q3sq*epsrrrr(mu2,mu4,q1,q2)-2*q1J3*q4J1*epsrrrr(mu2,mu4,q1,q2)-
     &   2*q2J3*q4J1*epsrrrr(mu2,mu4,q1,q2)+q3J1*q4J3*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q2q4*epsrrrr(mu2,mu4,q1,q3)-q2J3*q4J1*epsrrrr(mu2,mu4,q1,q3)-
     &   q2J1*q4J3*epsrrrr(mu2,mu4,q1,q3)-2*J1J3*q1q2*epsrrrr(mu2,mu4,q1,q4)-
     &   J1J3*q1sq*epsrrrr(mu2,mu4,q1,q4)+2*q1J3*q2J1*epsrrrr(mu2,mu4,q1,q4)+
     &   2*q2J1*q2J3*epsrrrr(mu2,mu4,q1,q4)-J1J3*q2q3*epsrrrr(mu2,mu4,q1,q4)-
     &   J1J3*q2sq*epsrrrr(mu2,mu4,q1,q4)+q2J3*q3J1*epsrrrr(mu2,mu4,q1,q4)-
     &   J1J3*q1q4*epsrrrr(mu2,mu4,q2,q3)+q1J3*q4J1*epsrrrr(mu2,mu4,q2,q3)+
     &   J1J3*q1q3*epsrrrr(mu2,mu4,q2,q4)-q1J3*q3J1*epsrrrr(mu2,mu4,q2,q4)-
     &   J1J3*q1q2*epsrrrr(mu2,mu4,q3,q4)-J1J3*q1sq*epsrrrr(mu2,mu4,q3,q4)+
     &   q1J3*q2J1*epsrrrr(mu2,mu4,q3,q4)+2*J1J3*q1J4*epsrrrr(mu2,q1,q2,q3)+
     &   2*J1J3*q2J4*epsrrrr(mu2,q1,q2,q3)+2*J1J3*q3J4*epsrrrr(mu2,q1,q2,q3)+
     &   J3J4*q4J1*epsrrrr(mu2,q1,q2,q3)-J1J4*q4J3*epsrrrr(mu2,q1,q2,q3)+
     &   2*J1J4*q1J3*epsrrrr(mu2,q1,q2,q4)+2*J1J4*q2J3*epsrrrr(mu2,q1,q2,q4)-
     &   J3J4*q3J1*epsrrrr(mu2,q1,q2,q4)+J1J3*q3J4*epsrrrr(mu2,q1,q2,q4)+
     &   J3J4*q2J1*epsrrrr(mu2,q1,q3,q4)+J1J4*q2J3*epsrrrr(mu2,q1,q3,q4)-
     &   J1J3*q2J4*epsrrrr(mu2,q1,q3,q4)-J1J4*q1J3*epsrrrr(mu2,q2,q3,q4)+
     &   J1J3*q1J4*epsrrrr(mu2,q2,q3,q4)+2*J1J2*q1q2*epsrrrr(mu3,mu4,q1,q2)+
     &   2*J1J2*q1q3*epsrrrr(mu3,mu4,q1,q2)+J1J2*q1sq*epsrrrr(mu3,mu4,q1,q2)+
     &   2*J1J2*q2q3*epsrrrr(mu3,mu4,q1,q2)+J1J2*q2sq*epsrrrr(mu3,mu4,q1,q2)+
     &   J1J2*q3q4*epsrrrr(mu3,mu4,q1,q2)+J1J2*q3sq*epsrrrr(mu3,mu4,q1,q2)+
     &   q3J2*q4J1*epsrrrr(mu3,mu4,q1,q2)-q3J1*q4J2*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q2q4*epsrrrr(mu3,mu4,q1,q3)-2*q1J2*q4J1*epsrrrr(mu3,mu4,q1,q3)+
     &   q2J1*q4J2*epsrrrr(mu3,mu4,q1,q3)+2*J1J2*q1q2*epsrrrr(mu3,mu4,q1,q4)+
     &   J1J2*q1sq*epsrrrr(mu3,mu4,q1,q4)+J1J2*q2q3*epsrrrr(mu3,mu4,q1,q4)+
     &   J1J2*q2sq*epsrrrr(mu3,mu4,q1,q4)+2*q1J2*q3J1*epsrrrr(mu3,mu4,q1,q4)-
     &   q2J1*q3J2*epsrrrr(mu3,mu4,q1,q4)+J1J2*q1q4*epsrrrr(mu3,mu4,q2,q3)-
     &   q1J2*q4J1*epsrrrr(mu3,mu4,q2,q3)-J1J2*q1q3*epsrrrr(mu3,mu4,q2,q4)+
     &   q1J2*q3J1*epsrrrr(mu3,mu4,q2,q4)+J1J2*q1q2*epsrrrr(mu3,mu4,q3,q4)+
     &   J1J2*q1sq*epsrrrr(mu3,mu4,q3,q4)-q1J2*q2J1*epsrrrr(mu3,mu4,q3,q4)-
     &   2*J1J2*q1J4*epsrrrr(mu3,q1,q2,q3)-2*J1J2*q2J4*epsrrrr(mu3,q1,q2,q3)-
     &   2*J1J2*q3J4*epsrrrr(mu3,q1,q2,q3)-J2J4*q4J1*epsrrrr(mu3,q1,q2,q3)+
     &   J1J4*q4J2*epsrrrr(mu3,q1,q2,q3)+J2J4*q3J1*epsrrrr(mu3,q1,q2,q4)-
     &   J1J4*q3J2*epsrrrr(mu3,q1,q2,q4)-J1J2*q3J4*epsrrrr(mu3,q1,q2,q4)+
     &   2*J1J4*q1J2*epsrrrr(mu3,q1,q3,q4)-J2J4*q2J1*epsrrrr(mu3,q1,q3,q4)+
     &   J1J2*q2J4*epsrrrr(mu3,q1,q3,q4)+J1J4*q1J2*epsrrrr(mu3,q2,q3,q4)-
     &   J1J2*q1J4*epsrrrr(mu3,q2,q3,q4)+J2J3*q4J1*epsrrrr(mu4,q1,q2,q3)-
     &   J1J3*q4J2*epsrrrr(mu4,q1,q2,q3)+J1J2*q4J3*epsrrrr(mu4,q1,q2,q3)-
     &   2*J1J2*q1J3*epsrrrr(mu4,q1,q2,q4)-2*J1J2*q2J3*epsrrrr(mu4,q1,q2,q4)-
     &   J2J3*q3J1*epsrrrr(mu4,q1,q2,q4)+J1J3*q3J2*epsrrrr(mu4,q1,q2,q4)-
     &   2*J1J3*q1J2*epsrrrr(mu4,q1,q3,q4)+J2J3*q2J1*epsrrrr(mu4,q1,q3,q4)-
     &   J1J2*q2J3*epsrrrr(mu4,q1,q3,q4)-J1J3*q1J2*epsrrrr(mu4,q2,q3,q4)+
     &   J1J2*q1J3*epsrrrr(mu4,q2,q3,q4)+
     &   (J1J4*J2J3-J1J3*J2J4+J1J2*J3J4)*epsrrrr(q1,q2,q3,q4))

      coeffEmu=8*m*(-(q1q2**2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4))-
     &   q1q2*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q1q4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   4*q1q2*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   3*q1q3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1sq**2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1sq*q2q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3q4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   3*q1q2**2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q1q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q1q4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q1sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   3*q1q2*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3q4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   3*q1q2*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3**2*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q1q4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*q1sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q2sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q3q4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   3*q1q2*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4**2*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q4*q1sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q4*q2q3*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q3sq*Eij(1,4)*epsrrrr(mu1,mu2,mu3,mu4)+
     &   4*q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   3*q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q2*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)-
     &   q1q3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1sq*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q2q3*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q2sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q2*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   3*q1q2*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1sq*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q2sq*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q3*q2J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q2*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   3*q2q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)-
     &   q2J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q4*q2J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   2*q2q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)-
     &   q2J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)-
     &   q1q3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)-
     &   q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q2)+
     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   3*q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1sq*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   3*q1q2*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1sq*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1q3*q2J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q2*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1sq*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1q4*q2J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
     &   2*q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4))

      coeffEmu=coeffEmu + 8*m*(q1sq*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1q2*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1q3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1q4*q3J4*Eij(1,4)*epsrrrr(mu1,mu2,mu3,q4)+
     &   q1q2*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1q3*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1sq*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1J3*q2q4*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1q2*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1sq*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1q2*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q2q4*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1q2*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q2sq*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1q3*q2J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q2q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q1q4*q2J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)+
     &   q2J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q1)-
     &   q1J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1q4*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1q4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q2)-
     &   q1sq*q4J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q3)-
     &   q1q2*q4J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q3)-
     &   q1q3*q4J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q3)-
     &   q1q4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q3)+
     &   2*q1J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
     &   2*q1sq*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)+
     &   3*q1q2*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4))

      coeffEmu=coeffEmu + 8*m*(q1sq*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)-
     &   q1J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
     &   2*q1q3*q2J3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)-
     &   q1J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
     &   2*q1q4*q2J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)-
     &   q1J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q2*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1sq*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,mu4,q4)+
     &   4*q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1q3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J4*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J3*q3J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
     &   4*q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q2J3*q3J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q2J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)-
     &   J3J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q1J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q2J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
     &   2*q3J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J4*q2J3*Eij(1,1)*epsrrrr(mu1,mu2,q1,q3)+
     &   q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q3)+
     &   q2J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q1,q3)-
     &   J3J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)+
     &   q1J3*q2J4*Eij(1,1)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)+
     &   q2J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)-
     &   J3J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)+
     &   q2J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q1,q4)-
     &   q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q2,q3)-
     &   q1J3*q2J4*Eij(1,2)*epsrrrr(mu1,mu2,q2,q3)-
     &   q1J3*q3J4*Eij(1,3)*epsrrrr(mu1,mu2,q2,q3)-
     &   q1J3*q1J4*Eij(1,1)*epsrrrr(mu1,mu2,q2,q4)-
     &   q1J4*q2J3*Eij(1,2)*epsrrrr(mu1,mu2,q2,q4)-
     &   q1J4*q4J3*Eij(1,4)*epsrrrr(mu1,mu2,q2,q4)-
     &   J3J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu2,q3,q4)-
     &   J3J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu2,q3,q4)-
     &   J3J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu2,q3,q4)-
     &   J3J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu2,q3,q4)+
     &   4*q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   4*q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   3*q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q2q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q2sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q3q4*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q3sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1sq*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2sq*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1q2*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1q3*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1sq*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q2q3*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2sq*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3J2*q3q4*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1))

      coeffEmu=coeffEmu + 8*m*(q3J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1J2*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1q2*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q1q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   2*q2q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3q4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3sq*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)+
     &   q1J2*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q2)-
     &   q1J2*q1q4*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q3)+
     &   q1sq*q4J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q3)+
     &   q1q2*q4J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q3)-
     &   q1q4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q3)+
     &   q1q3*q4J2*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q3)+
     &   q1J2*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)-
     &   q1sq*q3J2*Eij(1,1)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)-
     &   q1q2*q3J2*Eij(1,2)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1J2*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)-
     &   q1q4*q3J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
     &   q1q3*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,mu4,q4)+
     &   J2J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q1q3*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)-
     &   q1J4*q3J2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q2q3*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)-
     &   q2J4*q3J2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)-
     &   q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q3sq*Eij(1,3)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
     &   J2J4*q3q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q2)+
     &   4*q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J2*q2J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J2*q3J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J2*q2J4*Eij(1,2)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q2J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q1J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q2J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
     &   2*q3J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q3)+
     &   J2J4*q1q2*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
     &   q1J2*q3J4*Eij(1,1)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q2sq*Eij(1,2)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q2q3*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
     &   q3J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q1,q4)+
     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4))

      coeffEmu=coeffEmu + 8*m*(J2J4*q2q4*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4)+
     &   q3J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q1,q4)+
     &   q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q2,q3)+
     &   q1J2*q2J4*Eij(1,2)*epsrrrr(mu1,mu3,q2,q3)+
     &   q1J2*q3J4*Eij(1,3)*epsrrrr(mu1,mu3,q2,q3)-
     &   q1J2*q1J4*Eij(1,1)*epsrrrr(mu1,mu3,q3,q4)+
     &   J2J4*q1sq*Eij(1,1)*epsrrrr(mu1,mu3,q3,q4)+
     &   J2J4*q1q2*Eij(1,2)*epsrrrr(mu1,mu3,q3,q4)+
     &   J2J4*q1q3*Eij(1,3)*epsrrrr(mu1,mu3,q3,q4)-
     &   q1J4*q3J2*Eij(1,3)*epsrrrr(mu1,mu3,q3,q4)+
     &   J2J4*q1q4*Eij(1,4)*epsrrrr(mu1,mu3,q3,q4)-
     &   q1J4*q4J2*Eij(1,4)*epsrrrr(mu1,mu3,q3,q4)-
     &   J2J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q1q3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
     &   q1J3*q4J2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q2q3*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
     &   q2J3*q4J2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q3sq*Eij(1,3)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
     &   J2J3*q3q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
     &   q4J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q2)-
     &   q1J2*q4J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q3)-
     &   q3J2*q4J3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q3)-
     &   q4J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q3)+
     &   4*q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q1q2*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q1J2*q2J3*Eij(1,1)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q1J2*q2J3*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q2sq*Eij(1,2)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q2q3*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q1J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q2J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)-
     &   J2J3*q2q4*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q1J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q2J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
     &   2*q1J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q1,q4)+
     &   q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q2,q4)+
     &   q1J2*q2J3*Eij(1,2)*epsrrrr(mu1,mu4,q2,q4)+
     &   q1J2*q4J3*Eij(1,4)*epsrrrr(mu1,mu4,q2,q4)+
     &   q1J2*q1J3*Eij(1,1)*epsrrrr(mu1,mu4,q3,q4)-
     &   J2J3*q1sq*Eij(1,1)*epsrrrr(mu1,mu4,q3,q4)-
     &   J2J3*q1q2*Eij(1,2)*epsrrrr(mu1,mu4,q3,q4)-
     &   J2J3*q1q3*Eij(1,3)*epsrrrr(mu1,mu4,q3,q4)+
     &   q1J3*q3J2*Eij(1,3)*epsrrrr(mu1,mu4,q3,q4)-
     &   J2J3*q1q4*Eij(1,4)*epsrrrr(mu1,mu4,q3,q4)+
     &   q1J3*q4J2*Eij(1,4)*epsrrrr(mu1,mu4,q3,q4)+
     &   J2J3*q1J4*Eij(1,1)*epsrrrr(mu1,q1,q2,q3)+
     &   J2J3*q2J4*Eij(1,2)*epsrrrr(mu1,q1,q2,q3)+
     &   J2J3*q3J4*Eij(1,3)*epsrrrr(mu1,q1,q2,q3)+
     &   J2J4*q1J3*Eij(1,1)*epsrrrr(mu1,q1,q2,q4)+
     &   J2J4*q2J3*Eij(1,2)*epsrrrr(mu1,q1,q2,q4)+
     &   J2J4*q4J3*Eij(1,4)*epsrrrr(mu1,q1,q2,q4)+
     &   J3J4*q1J2*Eij(1,1)*epsrrrr(mu1,q1,q3,q4)+
     &   J3J4*q3J2*Eij(1,3)*epsrrrr(mu1,q1,q3,q4)+
     &   J3J4*q4J2*Eij(1,4)*epsrrrr(mu1,q1,q3,q4)-
     &   q1q2*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q3*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1sq*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q2*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1sq*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q2*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2sq*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q3*q2J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q4*q2J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2J1*q3q4*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q1q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)-
     &   q2q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q1)+
     &   2*q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q1q3*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q3q4*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q3sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q1q2*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q1q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q2q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3J1*q3q4*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3J1*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q1q2*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q1q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   2*q2q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q2)-
     &   q1sq*q4J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q2J1*q2q4*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q1q2*q4J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q2q4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q1q3*q4J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q1q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q3)-
     &   q2q4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q3)+
     &   q1sq*q3J1*Eij(1,1)*epsrrrr(mu2,mu3,mu4,q4)+
     &   2*q1q2*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4))

      coeffEmu=coeffEmu + 8*m*(q1sq*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q3J1*Eij(1,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   2*q1q2*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q3*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2sq*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q4*q3J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   2*q1q2*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q3*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2sq*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,mu4,q4)-
     &   J1J4*q1q2*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q1q3*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)+
     &   q1J4*q3J1*Eij(1,1)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q2q3*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)+
     &   q2J4*q3J1*Eij(1,2)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)+
     &   q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q3sq*Eij(1,3)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
     &   J1J4*q3q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q2)-
     &   q1J4*q2J1*Eij(1,1)*epsrrrr(mu2,mu3,q1,q3)-
     &   q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q1,q3)-
     &   q2J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q1,q3)-
     &   J1J4*q1q2*Eij(1,1)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q2sq*Eij(1,2)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q2q3*Eij(1,3)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q4)-
     &   J1J4*q2q4*Eij(1,4)*epsrrrr(mu2,mu3,q1,q4)+
     &   2*q1J4*q2J1*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q2J1*q3J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q1J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q2J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q1J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q2J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
     &   2*q3J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q3)+
     &   q2J1*q3J4*Eij(1,2)*epsrrrr(mu2,mu3,q2,q4)+
     &   q3J1*q3J4*Eij(1,3)*epsrrrr(mu2,mu3,q2,q4)+
     &   q3J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q2,q4)-
     &   J1J4*q1sq*Eij(1,1)*epsrrrr(mu2,mu3,q3,q4)-
     &   J1J4*q1q2*Eij(1,2)*epsrrrr(mu2,mu3,q3,q4)-
     &   q2J1*q2J4*Eij(1,2)*epsrrrr(mu2,mu3,q3,q4)-
     &   J1J4*q1q3*Eij(1,3)*epsrrrr(mu2,mu3,q3,q4)-
     &   q2J4*q3J1*Eij(1,3)*epsrrrr(mu2,mu3,q3,q4)-
     &   J1J4*q1q4*Eij(1,4)*epsrrrr(mu2,mu3,q3,q4)-
     &   q2J4*q4J1*Eij(1,4)*epsrrrr(mu2,mu3,q3,q4)+
     &   J1J3*q1q2*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1q3*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
     &   q1J3*q4J1*Eij(1,1)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q2q3*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q2sq*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
     &   q2J3*q4J1*Eij(1,2)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2))

      coeffEmu=coeffEmu + 8*m*(J1J3*q2q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q3sq*Eij(1,3)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q2q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q3q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
     &   q4J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q1,q2)+
     &   J1J3*q1q2*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)-
     &   q1J3*q2J1*Eij(1,1)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)-
     &   q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q2sq*Eij(1,2)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q2q3*Eij(1,3)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)+
     &   J1J3*q2q4*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)-
     &   q2J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q1,q4)-
     &   q2J1*q4J3*Eij(1,2)*epsrrrr(mu2,mu4,q2,q3)-
     &   q3J1*q4J3*Eij(1,3)*epsrrrr(mu2,mu4,q2,q3)-
     &   q4J1*q4J3*Eij(1,4)*epsrrrr(mu2,mu4,q2,q3)+
     &   2*q1J3*q2J1*Eij(1,2)*epsrrrr(mu2,mu4,q2,q4)+
     &   2*q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q2,q4)+
     &   2*q1J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q2,q4)+
     &   2*q2J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q2,q4)+
     &   2*q1J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q2,q4)+
     &   2*q2J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q2,q4)+
     &   J1J3*q1sq*Eij(1,1)*epsrrrr(mu2,mu4,q3,q4)+
     &   J1J3*q1q2*Eij(1,2)*epsrrrr(mu2,mu4,q3,q4)+
     &   q2J1*q2J3*Eij(1,2)*epsrrrr(mu2,mu4,q3,q4)+
     &   J1J3*q1q3*Eij(1,3)*epsrrrr(mu2,mu4,q3,q4)+
     &   q2J3*q3J1*Eij(1,3)*epsrrrr(mu2,mu4,q3,q4)+
     &   J1J3*q1q4*Eij(1,4)*epsrrrr(mu2,mu4,q3,q4)+
     &   q2J3*q4J1*Eij(1,4)*epsrrrr(mu2,mu4,q3,q4)-
     &   J1J3*q1J4*Eij(1,1)*epsrrrr(mu2,q1,q2,q3)-
     &   J1J3*q2J4*Eij(1,2)*epsrrrr(mu2,q1,q2,q3)-
     &   J1J3*q3J4*Eij(1,3)*epsrrrr(mu2,q1,q2,q3)-
     &   J1J4*q1J3*Eij(1,1)*epsrrrr(mu2,q1,q2,q4)-
     &   J1J4*q2J3*Eij(1,2)*epsrrrr(mu2,q1,q2,q4)-
     &   J1J4*q4J3*Eij(1,4)*epsrrrr(mu2,q1,q2,q4)+
     &   J3J4*q2J1*Eij(1,2)*epsrrrr(mu2,q2,q3,q4)+
     &   J3J4*q3J1*Eij(1,3)*epsrrrr(mu2,q2,q3,q4)+
     &   J3J4*q4J1*Eij(1,4)*epsrrrr(mu2,q2,q3,q4)-
     &   J1J2*q1q2*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q1q3*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q2q3*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q2sq*Eij(1,2)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q2q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q3sq*Eij(1,3)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q2q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)-
     &   J1J2*q3q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q2)+
     &   q1J2*q4J1*Eij(1,1)*epsrrrr(mu3,mu4,q1,q3)+
     &   q3J2*q4J1*Eij(1,3)*epsrrrr(mu3,mu4,q1,q3)+
     &   q4J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q1,q3)-
     &   J1J2*q1q2*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
     &   q1J2*q3J1*Eij(1,1)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q2sq*Eij(1,2)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q2q3*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
     &   q3J1*q3J2*Eij(1,3)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)-
     &   J1J2*q2q4*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)-
     &   q3J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q1,q4)+
     &   q2J1*q4J2*Eij(1,2)*epsrrrr(mu3,mu4,q2,q3)+
     &   q3J1*q4J2*Eij(1,3)*epsrrrr(mu3,mu4,q2,q3)+
     &   q4J1*q4J2*Eij(1,4)*epsrrrr(mu3,mu4,q2,q3)-
     &   q2J1*q3J2*Eij(1,2)*epsrrrr(mu3,mu4,q2,q4)-
     &   q3J1*q3J2*Eij(1,3)*epsrrrr(mu3,mu4,q2,q4)-
     &   q3J2*q4J1*Eij(1,4)*epsrrrr(mu3,mu4,q2,q4)-
     &   J1J2*q1sq*Eij(1,1)*epsrrrr(mu3,mu4,q3,q4)-
     &   J1J2*q1q2*Eij(1,2)*epsrrrr(mu3,mu4,q3,q4)+
     &   2*q1J2*q2J1*Eij(1,2)*epsrrrr(mu3,mu4,q3,q4)-
     &   J1J2*q1q3*Eij(1,3)*epsrrrr(mu3,mu4,q3,q4))

      coeffEmu=coeffEmu + 8*m*(2*q1J2*q3J1*Eij(1,3)*epsrrrr(mu3,mu4,q3,q4)-
     &   J1J2*q1q4*Eij(1,4)*epsrrrr(mu3,mu4,q3,q4)+
     &   2*q1J2*q4J1*Eij(1,4)*epsrrrr(mu3,mu4,q3,q4)+
     &   J1J2*q1J4*Eij(1,1)*epsrrrr(mu3,q1,q2,q3)+
     &   J1J2*q2J4*Eij(1,2)*epsrrrr(mu3,q1,q2,q3)+
     &   J1J2*q3J4*Eij(1,3)*epsrrrr(mu3,q1,q2,q3)-
     &   J1J4*(q1J2*Eij(1,1)+q3J2*Eij(1,3)+q4J2*Eij(1,4))*
     &    epsrrrr(mu3,q1,q3,q4)-J2J4*q2J1*Eij(1,2)*epsrrrr(mu3,q2,q3,q4)-
     &   J2J4*q3J1*Eij(1,3)*epsrrrr(mu3,q2,q3,q4)-
     &   J2J4*q4J1*Eij(1,4)*epsrrrr(mu3,q2,q3,q4)+
     &   J1J2*q1J3*Eij(1,1)*epsrrrr(mu4,q1,q2,q4)+
     &   J1J2*q2J3*Eij(1,2)*epsrrrr(mu4,q1,q2,q4)+
     &   J1J2*q4J3*Eij(1,4)*epsrrrr(mu4,q1,q2,q4)+
     &   J1J3*(q1J2*Eij(1,1)+q3J2*Eij(1,3)+q4J2*Eij(1,4))*
     &    epsrrrr(mu4,q1,q3,q4)+J2J3*q2J1*Eij(1,2)*epsrrrr(mu4,q2,q3,q4)+
     &   J2J3*q3J1*Eij(1,3)*epsrrrr(mu4,q2,q3,q4)+
     &   J2J3*q4J1*Eij(1,4)*epsrrrr(mu4,q2,q3,q4))

      coeffEmunu=16*m*(-(q1q2*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4))-
     &   q1q3*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq**2*Eij(2,1)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2**2*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2q3*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2sq*Eij(2,2)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3**2*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q2q3*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q3sq*Eij(2,3)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4**2*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2q4*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q3q4*Eij(2,4)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2**2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q1q3*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q1sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q3*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3**2*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*q1sq*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3sq*Eij(2,6)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q4*q1sq*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q2q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*q3q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q1q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q2sq*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3sq*Eij(2,8)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q2*q1q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2q3*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q2q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2sq*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*q3q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   2*q1q3*q1q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q2q3*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q2q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*q3q4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q4*q3sq*Eij(2,10)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q2*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1q3*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)-
     &   q1sq*Eij(2,11)*epsrrrr(mu1,mu2,mu3,mu4)+
     &   q1J4*q1q2*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q2*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q2sq*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q2*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1sq*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q2sq*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q2q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q2*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1sq*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q2q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q3*q2J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q2q3*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q2*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2sq*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q4*q2J4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*q2q4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q2J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q1)+
     &   q1J4*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q2*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q3*q3J4*Eij(2,3)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J4*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1sq*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J4*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1sq*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J4*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q3*q2J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q2*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q4*q2J4*Eij(2,9)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1q4*q3J4*Eij(2,10)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J4*Eij(2,11)*epsrrrr(mu1,mu2,mu3,q3)+
     &   q1J3*q1sq*Eij(2,1)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q2*q2J3*Eij(2,2)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q4*q4J3*Eij(2,4)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q2*Eij(2,5)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1sq*q2J3*Eij(2,5)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q3*Eij(2,6)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1q4*Eij(2,7)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1sq*q4J3*Eij(2,7)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q3*q2J3*Eij(2,8)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q4*q2J3*Eij(2,9)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q2*q4J3*Eij(2,9)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1q3*q4J3*Eij(2,10)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*Eij(2,11)*epsrrrr(mu1,mu2,mu4,q4)+
     &   q1J3*q1J4*Eij(2,1)*epsrrrr(mu1,mu2,q1,q2)+
     &   q2J3*q2J4*Eij(2,2)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J4*q2J3*Eij(2,5)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J3*q2J4*Eij(2,5)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J3*q3J4*Eij(2,6)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J4*q4J3*Eij(2,7)*epsrrrr(mu1,mu2,q1,q2)+
     &   q2J3*q3J4*Eij(2,8)*epsrrrr(mu1,mu2,q1,q2)+
     &   q2J4*q4J3*Eij(2,9)*epsrrrr(mu1,mu2,q1,q2)+
     &   q3J4*q4J3*Eij(2,10)*epsrrrr(mu1,mu2,q1,q2)+
     &   J3J4*Eij(2,11)*epsrrrr(mu1,mu2,q1,q2)+
     &   q1J2*q1q2*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1q3*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1sq*Eij(2,1)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q3*q3J2*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q3*q3J2*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3J2*q3sq*Eij(2,3)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3q4*q4J2*Eij(2,4)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1q2*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q2q3*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q2sq*Eij(2,5)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1q3*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q2q3*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q3*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1sq*q3J2*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q3sq*Eij(2,6)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q1q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q2q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*q3q4*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q3*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1sq*q4J2*Eij(2,7)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q3*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2sq*q3J2*Eij(2,8)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q2*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q3*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2sq*q4J2*Eij(2,9)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q4*q3J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q4*q3J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3J2*q3q4*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1q3*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q2q3*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3sq*q4J2*Eij(2,10)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q1J2*Eij(2,11)*epsrrrr(mu1,mu3,mu4,q1)+
     &   q3J2*Eij(2,11)*epsrrrr(mu1,mu3,mu4,q1)+
     &   (q3J2*q3J4*Eij(2,3)+q1J4*q3J2*Eij(2,6)+
     &      q1J2*(q1J4*Eij(2,1)+q2J4*Eij(2,5)+q3J4*Eij(2,6))+
     &      q2J4*q3J2*Eij(2,8)+q4J2*
     &       (q1J4*Eij(2,7)+q2J4*Eij(2,9)+q3J4*Eij(2,10))+J2J4*Eij(2,11))*
     &    epsrrrr(mu1,mu3,q1,q3)+
     &   (q4J2*q4J3*Eij(2,4)+q1J3*q4J2*Eij(2,7)+
     &      q1J2*(q1J3*Eij(2,1)+q2J3*Eij(2,5)+q4J3*Eij(2,7))+
     &      q2J3*q4J2*Eij(2,9)+q3J2*
     &       (q1J3*Eij(2,6)+q2J3*Eij(2,8)+q4J3*Eij(2,10))+J2J3*Eij(2,11))*
     &    epsrrrr(mu1,mu4,q1,q4)+q1q2*q2J1*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q2q3*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q2sq*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3J1*q3sq*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1sq*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q2J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q2q3*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q3*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2sq*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q3sq*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q4*q2J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q2q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*q3q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q3*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2sq*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3J1*q3q4*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3sq*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q2J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q3J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q2)+
     &   q1q2*q2J1*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2sq*Eij(2,2)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q3*q3J1*Eij(2,3)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q4*q4J1*Eij(2,4)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1sq*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1sq*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1sq*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q3*q2J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2q3*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2sq*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q4*q2J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2q4*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q2*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2sq*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q4*q3J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q1q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2q3*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*Eij(2,11)*epsrrrr(mu2,mu3,mu4,q4)+
     &   q2J1*q2J4*Eij(2,2)*epsrrrr(mu2,mu3,q2,q3)+
     &   q3J1*q3J4*Eij(2,3)*epsrrrr(mu2,mu3,q2,q3)+
     &   q1J4*q2J1*Eij(2,5)*epsrrrr(mu2,mu3,q2,q3)+
     &   q1J4*q3J1*Eij(2,6)*epsrrrr(mu2,mu3,q2,q3)+
     &   q1J4*q4J1*Eij(2,7)*epsrrrr(mu2,mu3,q2,q3)+
     &   q2J4*q3J1*Eij(2,8)*epsrrrr(mu2,mu3,q2,q3)+
     &   q2J1*q3J4*Eij(2,8)*epsrrrr(mu2,mu3,q2,q3)+
     &   q2J4*q4J1*Eij(2,9)*epsrrrr(mu2,mu3,q2,q3)+
     &   q3J4*q4J1*Eij(2,10)*epsrrrr(mu2,mu3,q2,q3)+
     &   J1J4*Eij(2,11)*epsrrrr(mu2,mu3,q2,q3)+
     &   q2J1*q2J3*Eij(2,2)*epsrrrr(mu2,mu4,q2,q4)+
     &   q4J1*q4J3*Eij(2,4)*epsrrrr(mu2,mu4,q2,q4)+
     &   q1J3*q2J1*Eij(2,5)*epsrrrr(mu2,mu4,q2,q4)+
     &   q1J3*q3J1*Eij(2,6)*epsrrrr(mu2,mu4,q2,q4)+
     &   q1J3*q4J1*Eij(2,7)*epsrrrr(mu2,mu4,q2,q4)+
     &   q2J3*q3J1*Eij(2,8)*epsrrrr(mu2,mu4,q2,q4)+
     &   q2J3*q4J1*Eij(2,9)*epsrrrr(mu2,mu4,q2,q4)+
     &   q2J1*q4J3*Eij(2,9)*epsrrrr(mu2,mu4,q2,q4)+
     &   q3J1*q4J3*Eij(2,10)*epsrrrr(mu2,mu4,q2,q4)+
     &   J1J3*Eij(2,11)*epsrrrr(mu2,mu4,q2,q4)+
     &   q3J1*q3J2*Eij(2,3)*epsrrrr(mu3,mu4,q3,q4)+
     &   q4J1*q4J2*Eij(2,4)*epsrrrr(mu3,mu4,q3,q4)+
     &   q1J2*q2J1*Eij(2,5)*epsrrrr(mu3,mu4,q3,q4)+
     &   q1J2*q3J1*Eij(2,6)*epsrrrr(mu3,mu4,q3,q4)+
     &   q1J2*q4J1*Eij(2,7)*epsrrrr(mu3,mu4,q3,q4)+
     &   q2J1*q3J2*Eij(2,8)*epsrrrr(mu3,mu4,q3,q4)+
     &   q2J1*q4J2*Eij(2,9)*epsrrrr(mu3,mu4,q3,q4)+
     &   q3J2*q4J1*Eij(2,10)*epsrrrr(mu3,mu4,q3,q4)+
     &   q3J1*q4J2*Eij(2,10)*epsrrrr(mu3,mu4,q3,q4)+
     &   J1J2*Eij(2,11)*epsrrrr(mu3,mu4,q3,q4))

      A4G_LoopDD=(coeffD0+coeffDmu+coeffE0C0+coeffE0+coeffEmu
     &	         +coeffEmunu)/(4d0*m)

      RETURN
      END
