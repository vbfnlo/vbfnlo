***********************************************************************
      SUBROUTINE H4G_LoopDD(i1,i2,i3,i4,m,pent)

      IMPLICIT NONE
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

      q1sq=q_dot_q(i1,i1)
      q2sq=q_dot_q(i2,i2)
      q3sq=q_dot_q(i3,i3)
      q4sq=q_dot_q(i4,i4)

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

      C034=-C0_1M(m,q3sq,q4sq,q3sq+q4sq+2d0*q3q4) 

      D0234=D0_1M(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4)
      CALL D_1M_ALL2(m,q2sq,q3sq,q4sq,q2q3,q2q4,q3q4,Dij)
      Dij(2,7) =-Dij(2,7)
      Dij(3,11)=-Dij(3,11)
      Dij(3,12)=-Dij(3,12)
      Dij(3,13)=-Dij(3,13)

      a = H4G_LoopDD2a(q1sq,q2sq,q3sq,q4sq,
     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
     &                 C034,D0234,Dij,E0L,Eij,m)

      b = H4G_LoopDD2b(q1sq,q2sq,q3sq,q4sq,
     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
     &                 C034,D0234,Dij,E0L,Eij,m)

      c = H4G_LoopDD2c(q1sq,q2sq,q3sq,q4sq,
     &             q1q2,q1q3,q1q4,q2q3,q2q4,q3q4,
     &   q1J1,q2J1,q3J1,q4J1,q1J2,q2J2,q3J2,q4J2,
     &   q1J3,q2J3,q3J3,q4J3,q1J4,q2J4,q3J4,q4J4,
     &             J1J2,J1J3,J1J4,J2J3,J2J4,J3J4,
     &                 C034,D0234,Dij,E0L,Eij,m)

      pent=1d0/(4d0*m)*(a+b+c)

      END
