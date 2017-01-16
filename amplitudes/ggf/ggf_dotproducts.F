***********************************************************************
      SUBROUTINE auxiliary_calculation(epvec,qvec)
*----------------------------------------------------------------------
c  assembled and modified: Michael Kubocz
c  date                  : 14.08.06    
c  modified, date        : 30.01.11
***********************************************************************

      IMPLICIT NONE
      
#include "ggf_process.inc"
#include "ggf_dotproducts.inc"

      real*8 eps,qi
      common /bh4g_mom_eps/ qi(0:3,4),eps(0:3,4)

      complex*16 epvec(0:3,4)
      real*8 qvec(0:4,4)
      integer i,j,k
      real*8 tmp

      real*8 dotrr
      complex*16 dotrc,dotcc

      EXTERNAL dotrr,dotrc,dotcc
           
      do i=1,4
         do j=i,4
            if (i.eq.j) then
               q_dot_q(i,j)=0d0
            else
               q_dot_q(i,j)=dotrr(qvec(0,i),qvec(0,j))
               q_dot_q(j,i) = q_dot_q(i,j)
            endif
         enddo
      enddo
      
      CALL Xij_fill(q_dot_q,Xi,det)   
      
      do i=1,4
         do j=i,4
            if (i.eq.j) then
               q_dot_ep(i,j) =DCMPLX(0d0,0d0)
               ep_dot_ep(i,j)=DCMPLX(-1d0,0d0)
            else
               q_dot_ep(i,j) =dotrc(qvec(0,i),epvec(0,j))
               ep_dot_ep(i,j)=dotcc(epvec(0,i),epvec(0,j))
            endif
            q_dot_ep(j,i) = dotrc(qvec(0,j),epvec(0,i))
            ep_dot_ep(j,i)= ep_dot_ep(i,j)
         enddo
      enddo

      if (det.eq.0) then
         do i=1,4
            do j=1,4            
               ep(i,j) = 0d0
            enddo
         enddo
         return
      endif
         
      do i=1,4
         do j=1,4            
            tmp = 0
            do k=1,4
               tmp = tmp - Xi(j,k)*dreal(q_dot_ep(k,i))
            enddo
            ep(i,j) = tmp/det
         enddo
      enddo

      END


      SUBROUTINE calculate_E_functions()
*----------------------------------------------------------------------
c... author    : Michael Kubocz
c... date      : 20.07.08
c... mopdified : 30.01.11
***********************************************************************

      IMPLICIT NONE

#include "ggf_process.inc"
#include "ggf_dotproducts.inc"

      if (HloopTop.or.AloopTop.or.sHloopTop.or.bHloopTop) then
         CALL E01M_off(4,2,1,3,m_top,E0t(1))
         CALL E01M_off(2,1,3,4,m_top,E0t(2))
         CALL E01M_off(3,2,1,4,m_top,E0t(3))
         CALL E01M_off(2,1,4,3,m_top,E0t(4))
         CALL E01M_off(4,3,2,1,m_top,E0t(5))
         CALL E01M_off(3,4,2,1,m_top,E0t(6))
         CALL E01M_off(2,3,1,4,m_top,E0t(7))
         CALL E01M_off(3,1,4,2,m_top,E0t(8))
         CALL E01M_off(4,2,3,1,m_top,E0t(9))
         CALL E01M_off(2,4,3,1,m_top,E0t(10))
         CALL E01M_off(3,2,4,1,m_top,E0t(11))
         CALL E01M_off(2,3,4,1,m_top,E0t(12))

         CALL DDtens51m2(4,2,1,3,m_top,DDEijtop(1,1,1))
         CALL DDtens51m2(2,1,3,4,m_top,DDEijtop(1,1,2)) 
         CALL DDtens51m2(3,2,1,4,m_top,DDEijtop(1,1,3))
         CALL DDtens51m2(2,1,4,3,m_top,DDEijtop(1,1,4))
         CALL DDtens51m2(4,3,2,1,m_top,DDEijtop(1,1,5))
         CALL DDtens51m2(3,4,2,1,m_top,DDEijtop(1,1,6))
         CALL DDtens51m2(2,3,1,4,m_top,DDEijtop(1,1,7))
         CALL DDtens51m2(3,1,4,2,m_top,DDEijtop(1,1,8))
         CALL DDtens51m2(4,2,3,1,m_top,DDEijtop(1,1,9))
         CALL DDtens51m2(2,4,3,1,m_top,DDEijtop(1,1,10))
         CALL DDtens51m2(3,2,4,1,m_top,DDEijtop(1,1,11))
         CALL DDtens51m2(2,3,4,1,m_top,DDEijtop(1,1,12))
      endif

      if (HloopBottom.or.AloopBottom.or.sHloopBottom.or.bHloopBottom) then
         CALL E01M_off(4,2,1,3,m_bottom,E0b(1))
         CALL E01M_off(2,1,3,4,m_bottom,E0b(2))
         CALL E01M_off(3,2,1,4,m_bottom,E0b(3))
         CALL E01M_off(2,1,4,3,m_bottom,E0b(4))
         CALL E01M_off(4,3,2,1,m_bottom,E0b(5))
         CALL E01M_off(3,4,2,1,m_bottom,E0b(6)) 
         CALL E01M_off(2,3,1,4,m_bottom,E0b(7))
         CALL E01M_off(3,1,4,2,m_bottom,E0b(8))
         CALL E01M_off(4,2,3,1,m_bottom,E0b(9))
         CALL E01M_off(2,4,3,1,m_bottom,E0b(10))
         CALL E01M_off(3,2,4,1,m_bottom,E0b(11))
         CALL E01M_off(2,3,4,1,m_bottom,E0b(12))

         CALL DDtens51m2(4,2,1,3,m_bottom,DDEijbottom(1,1,1))
         CALL DDtens51m2(2,1,3,4,m_bottom,DDEijbottom(1,1,2)) 
         CALL DDtens51m2(3,2,1,4,m_bottom,DDEijbottom(1,1,3))
         CALL DDtens51m2(2,1,4,3,m_bottom,DDEijbottom(1,1,4))
         CALL DDtens51m2(4,3,2,1,m_bottom,DDEijbottom(1,1,5))
         CALL DDtens51m2(3,4,2,1,m_bottom,DDEijbottom(1,1,6))
         CALL DDtens51m2(2,3,1,4,m_bottom,DDEijbottom(1,1,7))
         CALL DDtens51m2(3,1,4,2,m_bottom,DDEijbottom(1,1,8))
         CALL DDtens51m2(4,2,3,1,m_bottom,DDEijbottom(1,1,9))
         CALL DDtens51m2(2,4,3,1,m_bottom,DDEijbottom(1,1,10))
         CALL DDtens51m2(3,2,4,1,m_bottom,DDEijbottom(1,1,11))
         CALL DDtens51m2(2,3,4,1,m_bottom,DDEijbottom(1,1,12))
      endif

c$$$      if (sHloopStop1.or.bHloopStop1) then
c$$$         CALL DDtens51m(4,2,1,3,m_stop(1),DDEijstop1(1,1,1))
c$$$         CALL DDtens51m(2,1,3,4,m_stop(1),DDEijstop1(1,1,2)) 
c$$$         CALL DDtens51m(3,2,1,4,m_stop(1),DDEijstop1(1,1,3))
c$$$         CALL DDtens51m(2,1,4,3,m_stop(1),DDEijstop1(1,1,4))
c$$$         CALL DDtens51m(4,3,2,1,m_stop(1),DDEijstop1(1,1,5))
c$$$         CALL DDtens51m(3,4,2,1,m_stop(1),DDEijstop1(1,1,6))
c$$$         CALL DDtens51m(2,3,1,4,m_stop(1),DDEijstop1(1,1,7))
c$$$         CALL DDtens51m(3,1,4,2,m_stop(1),DDEijstop1(1,1,8))
c$$$         CALL DDtens51m(4,2,3,1,m_stop(1),DDEijstop1(1,1,9))
c$$$         CALL DDtens51m(2,4,3,1,m_stop(1),DDEijstop1(1,1,10))
c$$$         CALL DDtens51m(3,2,4,1,m_stop(1),DDEijstop1(1,1,11))
c$$$         CALL DDtens51m(2,3,4,1,m_stop(1),DDEijstop1(1,1,12))
c$$$      endif
c$$$      if (sHloopStop2.or.bHloopStop2) then
c$$$         CALL DDtens51m(4,2,1,3,m_stop(2),DDEijstop2(1,1,1))
c$$$         CALL DDtens51m(2,1,3,4,m_stop(2),DDEijstop2(1,1,2)) 
c$$$         CALL DDtens51m(3,2,1,4,m_stop(2),DDEijstop2(1,1,3))
c$$$         CALL DDtens51m(2,1,4,3,m_stop(2),DDEijstop2(1,1,4))
c$$$         CALL DDtens51m(4,3,2,1,m_stop(2),DDEijstop2(1,1,5))
c$$$         CALL DDtens51m(3,4,2,1,m_stop(2),DDEijstop2(1,1,6))
c$$$         CALL DDtens51m(2,3,1,4,m_stop(2),DDEijstop2(1,1,7))
c$$$         CALL DDtens51m(3,1,4,2,m_stop(2),DDEijstop2(1,1,8))
c$$$         CALL DDtens51m(4,2,3,1,m_stop(2),DDEijstop2(1,1,9))
c$$$         CALL DDtens51m(2,4,3,1,m_stop(2),DDEijstop2(1,1,10))
c$$$         CALL DDtens51m(3,2,4,1,m_stop(2),DDEijstop2(1,1,11))
c$$$         CALL DDtens51m(2,3,4,1,m_stop(2),DDEijstop2(1,1,12))
c$$$      endif
c$$$      if (sHloopSbottom1.or.bHloopSbottom1) then
c$$$         CALL DDtens51m(4,2,1,3,m_sbottom(1),DDEijsbottom1(1,1,1))
c$$$         CALL DDtens51m(2,1,3,4,m_sbottom(1),DDEijsbottom1(1,1,2)) 
c$$$         CALL DDtens51m(3,2,1,4,m_sbottom(1),DDEijsbottom1(1,1,3))
c$$$         CALL DDtens51m(2,1,4,3,m_sbottom(1),DDEijsbottom1(1,1,4))
c$$$         CALL DDtens51m(4,3,2,1,m_sbottom(1),DDEijsbottom1(1,1,5))
c$$$         CALL DDtens51m(3,4,2,1,m_sbottom(1),DDEijsbottom1(1,1,6))
c$$$         CALL DDtens51m(2,3,1,4,m_sbottom(1),DDEijsbottom1(1,1,7))
c$$$         CALL DDtens51m(3,1,4,2,m_sbottom(1),DDEijsbottom1(1,1,8))
c$$$         CALL DDtens51m(4,2,3,1,m_sbottom(1),DDEijsbottom1(1,1,9))
c$$$         CALL DDtens51m(2,4,3,1,m_sbottom(1),DDEijsbottom1(1,1,10))
c$$$         CALL DDtens51m(3,2,4,1,m_sbottom(1),DDEijsbottom1(1,1,11))
c$$$         CALL DDtens51m(2,3,4,1,m_sbottom(1),DDEijsbottom1(1,1,12))
c$$$      endif
c$$$      if (sHloopSbottom2.or.bHloopSbottom2) then
c$$$         CALL DDtens51m(4,2,1,3,m_sbottom(2),DDEijsbottom2(1,1,1))
c$$$         CALL DDtens51m(2,1,3,4,m_sbottom(2),DDEijsbottom2(1,1,2)) 
c$$$         CALL DDtens51m(3,2,1,4,m_sbottom(2),DDEijsbottom2(1,1,3))
c$$$         CALL DDtens51m(2,1,4,3,m_sbottom(2),DDEijsbottom2(1,1,4))
c$$$         CALL DDtens51m(4,3,2,1,m_sbottom(2),DDEijsbottom2(1,1,5))
c$$$         CALL DDtens51m(3,4,2,1,m_sbottom(2),DDEijsbottom2(1,1,6))
c$$$         CALL DDtens51m(2,3,1,4,m_sbottom(2),DDEijsbottom2(1,1,7))
c$$$         CALL DDtens51m(3,1,4,2,m_sbottom(2),DDEijsbottom2(1,1,8))
c$$$         CALL DDtens51m(4,2,3,1,m_sbottom(2),DDEijsbottom2(1,1,9))
c$$$         CALL DDtens51m(2,4,3,1,m_sbottom(2),DDEijsbottom2(1,1,10))
c$$$         CALL DDtens51m(3,2,4,1,m_sbottom(2),DDEijsbottom2(1,1,11))
c$$$         CALL DDtens51m(2,3,4,1,m_sbottom(2),DDEijsbottom2(1,1,12))
c$$$      endif

      END

***********************************************************************      
      SUBROUTINE Xij_fill(q_dot_q,Xi,detx) 
*----------------------------------------------------------------------      
c   Xij_fill calculates the inverse of the 4 x 4 matrix X_ij = -q_i.q_j
c   and its determinant

c   input: q1sq =q_1.q_1 ... q3q4 = q_3.q_4     entry Xij_fill
c          q_dot_q(i,j) = q_i.q_j                   entry Xij_fillP

c   output: detx = det(X)
c           Xi = (inverse of X)*detx
***********************************************************************

      IMPLICIT NONE
      
      real*8 p1sq,p2sq,p3sq,p4sq, p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      real*8 q_dot_q(4,4), X(4,4), Xi(4,4), detx
      integer i,j
      logical ldebug
      parameter (ldebug=.false.)

      p1sq = q_dot_q(1,1)
      p2sq = q_dot_q(2,2)
      p3sq = q_dot_q(3,3)
      p4sq = q_dot_q(4,4)
      p1p2 = q_dot_q(1,2)
      p1p3 = q_dot_q(1,3)
      p1p4 = q_dot_q(1,4)
      p2p3 = q_dot_q(2,3)
      p2p4 = q_dot_q(2,4)
      p3p4 = q_dot_q(3,4)

      if (p1sq.eq.0d0 .and. p2sq.eq.0d0 .and. p3sq.eq.0d0 
     &    .and. p4sq.eq.0d0 ) then
         xi(1,1)=-2d0*p2p3*p2p4*p3p4
         xi(1,2)=p3p4*(-p1p2*p3p4+p1p3*p2p4+p1p4*p2p3)
         xi(1,3)=p2p4*(p1p2*p3p4-p1p3*p2p4+p1p4*p2p3)
         xi(1,4)=p2p3*(p1p2*p3p4+p1p3*p2p4-p1p4*p2p3)
         xi(2,1)=xi(1,2)
         xi(2,2)=-2d0*p1p3*p1p4*p3p4
         xi(2,3)=p1p4*(p1p2*p3p4+p1p3*p2p4-p1p4*p2p3)
         xi(2,4)=p1p3*(p1p2*p3p4-p1p3*p2p4+p1p4*p2p3)
         xi(3,1)=xi(1,3)
         xi(3,2)=xi(2,3)
         xi(3,3)=-2d0*p1p2*p1p4*p2p4
         xi(3,4)=p1p2*(-p1p2*p3p4+p1p3*p2p4+p1p4*p2p3)
         xi(4,1)=xi(1,4)
         xi(4,2)=xi(2,4)
         xi(4,3)=xi(3,4)
         xi(4,4)=-2d0*p1p2*p1p3*p2p3

         detx=p1p2**2*p3p4**2-2d0*p1p2*(p1p3*p2p4*p3p4+p1p4*p2p3*p3p4)
     &       +p1p3**2*p2p4**2-2d0*p1p3*p1p4*p2p3*p2p4+p1p4**2*p2p3**2
c         print*,'detX 1 = ',detx
      else
         detx = p1p2**2*p3p4**2-p1p2**2*p3sq*p4sq+
     &      2d0*p1p2*p1p3*p2p3*p4sq-2d0*p1p2*p1p3*p2p4*p3p4-
     &      2d0*p1p2*p1p4*p2p3*p3p4+2d0*p1p2*p1p4*p2p4*p3sq+
     &      p1p3**2*p2p4**2-p1p3**2*p2sq*p4sq-2d0*p1p3*p1p4*p2p3*p2p4+
     &      2d0*p1p3*p1p4*p2sq*p3p4+p1p4**2*p2p3**2-p1p4**2*p2sq*p3sq-
     &      p1sq*p2p3**2*p4sq+2d0*p1sq*p2p3*p2p4*p3p4-p1sq*p2p4**2*p3sq-
     &      p1sq*p2sq*p3p4**2+p1sq*p2sq*p3sq*p4sq

c         print*,'detX 2 = ',detx

         xi(1,1)=p2p3**2*p4sq-2d0*p2p3*p2p4*p3p4+p2p4**2*p3sq+
     &           p2sq*p3p4**2-p2sq*p3sq*p4sq
         xi(1,2)=
     &        -p1p2*p3p4**2+p1p2*p3sq*p4sq-p1p3*p2p3*p4sq+p1p3*p2p4*
     &        p3p4+p1p4*p2p3*p3p4-p1p4*p2p4*p3sq
         xi(1,3)=
     &        -p1p2*p2p3*p4sq+p1p2*p2p4*p3p4-p1p3*p2p4**2+p1p3*p2sq*
     &        p4sq+p1p4*p2p3*p2p4-p1p4*p2sq*p3p4
         xi(1,4)=
     &        p1p2*p2p3*p3p4-p1p2*p2p4*p3sq+p1p3*p2p3*p2p4-p1p3*p2sq
     &        *p3p4-p1p4*p2p3**2+p1p4*p2sq*p3sq
         xi(2,2)=
     &        p1p3**2*p4sq-2d0*p1p3*p1p4*p3p4+p1p4**2*p3sq+
     &        p1sq*p3p4**2-p1sq*p3sq*p4sq
         xi(2,3)=
     &        -p1p2*p1p3*p4sq+p1p2*p1p4*p3p4+p1p3*p1p4*p2p4-p1p4**2*
     &        p2p3+p1sq*p2p3*p4sq-p1sq*p2p4*p3p4
         xi(2,4)=
     &        p1p2*p1p3*p3p4-p1p2*p1p4*p3sq-p1p3**2*p2p4+p1p3*p1p4*
     &        p2p3-p1sq*p2p3*p3p4+p1sq*p2p4*p3sq
         xi(3,3)=
     &        p1p2**2*p4sq-2d0*p1p2*p1p4*p2p4+p1p4**2*p2sq+
     &        p1sq*p2p4**2-p1sq*p2sq*p4sq
         xi(3,4)=
     &        -p1p2**2*p3p4+p1p2*p1p3*p2p4+p1p2*p1p4*p2p3-p1p3*p1p4*
     &        p2sq-p1sq*p2p3*p2p4+p1sq*p2sq*p3p4
         xi(4,4)=
     &        p1p2**2*p3sq-2d0*p1p2*p1p3*p2p3+p1p3**2*p2sq+
     &        p1sq*p2p3**2-p1sq*p2sq*p3sq

         do i = 2,4
            do j = 1,i-1
               xi(i,j) = xi(j,i)
            enddo
         enddo
      endif

      if (ldebug) then
         print*,' detx = ',detx
         print*,' check inversion of X matrix for n = 4 '
         X(1,1) = -p1sq
         X(1,2) = -p1p2
         X(1,3) = -p1p3
         X(1,4) = -p1p4
         X(2,1) = -p1p2
         X(2,2) = -p2sq
         X(2,3) = -p2p3
         X(2,4) = -p2p4
         X(3,1) = -p1p3
         X(3,2) = -p2p3
         X(3,3) = -p3sq
         X(3,4) = -p3p4
         X(4,1) = -p1p4
         X(4,2) = -p2p4
         X(4,3) = -p3p4
         X(4,4) = -p4sq

         do i = 1,4
            do j = 1,4
               print*,'X*X^-1(',i,j,') = ', 
     &              (X(i,1)*Xi(1,j)+X(i,2)*Xi(2,j)+X(i,3)*Xi(3,j)+
     &               X(i,4)*Xi(4,j))/detx
            enddo
         enddo

         print*,' Xi = '

         do i = 1,4
            print 13,(Xi(i,j)/detx,j=1,4)
         enddo

         print*,' X = '

         do i = 1,4
            print 13,(X(i,j),j=1,4)
         enddo
 13      format (4g12.2)
      endif

      RETURN
	
      END
