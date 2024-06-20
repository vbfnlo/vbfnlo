***********************************************************************
      SUBROUTINE DDtens51m2(i1,i2,i3,i4,m,Eij)
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... modified : 30.31.11
c... prepares all Eij coefficients for pentagon topologies 
*----------------------------------------------------------------------

      IMPLICIT NONE

#include "d_coefficients.inc"

      integer i1,i2,i3,i4,i,j

      double precision q_dot_q
      dimension q_dot_q(4,4)
      complex*16 q_dot_ep
      dimension q_dot_ep(4,4)
      complex*16 ep_dot_ep
      dimension ep_dot_ep(4,4)
      double precision ep
      dimension ep(4,4)
      double precision det,Xi(4,4)
      COMMON/auxiliary_calculation_GG/q_dot_ep,ep_dot_ep,ep,q_dot_q,
     &                                det,Xi

      real*8 p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,m
      complex*16 Eij(4,46)

      real*8 psq,qsq,lsq,tsq,pq,pl,pt,ql,qt,lt
      real*8 p_q,q_l,l_t,p_ql,p_qt,q_lp,l_tp,q_lt,l_tq

      complex*16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234
      complex*16 Dij_2345(13,3),Dij_1345(13,3),Dij_1245(13,3),
     &           Dij_1235(13,3),Dij_1234(13,3)
      real*8 DR2345,DR1345,DR1245,DR1235,DR1234,DRR1234(13,3),E(46,4)
      real*8 DRR2345(13,3),DRR1345(13,3),DRR1245(13,3),DRR1235(13,3)
      real*8 DI2345,DI1345,DI1245,DI1235,DI1234,DII1234(13,3),C(46,4)
      real*8 DII2345(13,3),DII1345(13,3),DII1245(13,3),DII1235(13,3)

      complex*16 D0_1M
      EXTERNAL D0_1M

      psq=q_dot_q(i1,i1)
      qsq=q_dot_q(i2,i2)
      lsq=q_dot_q(i3,i3)
      tsq=q_dot_q(i4,i4)

      pq=q_dot_q(i1,i2)
      pl=q_dot_q(i1,i3)
      pt=q_dot_q(i1,i4)
      ql=q_dot_q(i2,i3)
      qt=q_dot_q(i2,i4)
      lt=q_dot_q(i3,i4)

      p_q=(psq+qsq+2d0*pq)
      q_l=(qsq+lsq+2d0*ql)
      l_t=(lsq+tsq+2d0*lt)
      p_ql=(pl+ql)
      p_qt=(pt+qt)
      q_lp=(pq+pl)
      q_lt=(qt+lt)
      l_tp=(pl+pt)
      l_tq=(ql+qt)

c... calculate required D0 and Dij functions
      D0_1345=D0_1M (m,p_q,lsq,tsq,p_ql,p_qt,lt)
      CALL D_1M_ALL3(m,p_q,lsq,tsq,p_ql,p_qt,lt  ,Dij_1345) 
      D0_1245=D0_1M (m,psq,q_l,tsq,q_lp,pt  ,q_lt)
      CALL D_1M_ALL3(m,psq,q_l,tsq,q_lp,pt  ,q_lt,Dij_1245)
      D0_1234=D0_1M (m,psq,qsq,lsq,pq  ,pl  ,ql)
      CALL D_1M_ALL3(m,psq,qsq,lsq,pq  ,pl  ,ql  ,Dij_1234)
      D0_1235=D0_1M (m,psq,qsq,l_t,pq  ,l_tp,l_tq)
      CALL D_1M_ALL3(m,psq,qsq,l_t,pq  ,l_tp,l_tq,Dij_1235) 
      D0_2345=D0_1M (m,qsq,lsq,tsq,ql  ,qt  ,lt)
      CALL D_1M_ALL3(m,qsq,lsq,tsq,ql  ,qt  ,lt  ,Dij_2345) 

      DR2345=Dble(D0_2345)
      DR1345=Dble(D0_1345)
      DR1245=Dble(D0_1245)
      DR1235=Dble(D0_1235)
      DR1234=Dble(D0_1234)

      DI2345=DImag(D0_2345)
      DI1345=DImag(D0_1345)
      DI1245=DImag(D0_1245)
      DI1235=DImag(D0_1235)
      DI1234=DImag(D0_1234)

      do i=1,13
         do j=1,3
            DRR2345(i,j)=Dble(Dij_2345(i,j))
            DRR1345(i,j)=Dble(Dij_1345(i,j))
            DRR1245(i,j)=Dble(Dij_1245(i,j))
            DRR1235(i,j)=Dble(Dij_1235(i,j))
            DRR1234(i,j)=Dble(Dij_1234(i,j))
            DII2345(i,j)=DImag(Dij_2345(i,j))
            DII1345(i,j)=DImag(Dij_1345(i,j))
            DII1245(i,j)=DImag(Dij_1245(i,j))
            DII1235(i,j)=DImag(Dij_1235(i,j))
            DII1234(i,j)=DImag(Dij_1234(i,j))           
         enddo
      enddo

      p1=q_dot_q(i1,i1)
      p2=q_dot_q(i2,i2)
      p3=q_dot_q(i3,i3)
      p4=q_dot_q(i4,i4)
      p1p2=q_dot_q(i1,i2)
      p1p3=q_dot_q(i1,i3)
      p1p4=q_dot_q(i1,i4)
      p2p3=q_dot_q(i2,i3)
      p2p4=q_dot_q(i2,i4)
      p3p4=q_dot_q(i3,i4)

      CALL tens_red5_new_Re_Com_1M(m,p1, p2, p3, p4,
     &                     p1p2, p1p3,p1p4, p2p3,p2p4,p3p4, 
     &                     DR2345,DR1345,DR1245,DR1235,DR1234,
     &                     DRR2345,DRR1345,DRR1245,DRR1235,DRR1234,
     &                     DI2345,DI1345,DI1245,DI1235,DI1234,
     &                     DII2345,DII1345,DII1245,DII1235,DII1234,
     &                     E,C)

      Eij(1,1)=DCMPLX(E(1,1),C(1,1))
      Eij(1,2)=DCMPLX(E(2,1),C(2,1))
      Eij(1,3)=DCMPLX(E(3,1),C(3,1))    
      Eij(1,4)=DCMPLX(E(4,1),C(4,1))
c      write(*,*)  '************ 2 ***********'
c      write(*,*)  '************ 2 ***********'
      Eij(2,1)=DCMPLX(E(1,2),C(1,2))
      Eij(2,2)=DCMPLX(E(2,2),C(2,2))
      Eij(2,3)=DCMPLX(E(3,2),C(3,2))
      Eij(2,4)=DCMPLX(E(4,2),C(4,2))
      Eij(2,5)=DCMPLX(E(5,2),C(5,2))
      Eij(2,6)=DCMPLX(E(6,2),C(6,2))
      Eij(2,7)=DCMPLX(E(7,2),C(7,2) )  
      Eij(2,8)=DCMPLX(E(8,2),C(8,2))
      Eij(2,9)=DCMPLX(E(9,2),C(9,2))
      Eij(2,10)=DCMPLX(E(10,2),C(10,2))
      Eij(2,11)=DCMPLX(E(11,2),C(11,2))
c      write(*,*)  '************ 3 ***********'
c      write(*,*)  '************ 3 ***********'
      Eij(3,1)=DCMPLX(E(1,3), C(1,3))
      Eij(3,2)=DCMPLX(E(2,3), C(2,3))
      Eij(3,3)=DCMPLX(E(3,3), C(3,3))
      Eij(3,4)=DCMPLX(E(4,3), C(4,3))
      Eij(3,5)=DCMPLX(E(5,3), C(5,3))
      Eij(3,6)=DCMPLX(E(6,3), C(6,3))
      Eij(3,7)=DCMPLX(E(7,3), C(7,3) )   
      Eij(3,8)=DCMPLX(E(8,3), C(8,3))
      Eij(3,9)=DCMPLX(E(9,3), C(9,3))
      Eij(3,10)=DCMPLX(E(10,3), C(10,3))
      Eij(3,11)=DCMPLX(E(11,3), C(11,3))
      Eij(3,12)=DCMPLX(E(12,3), C(12,3))
      Eij(3,13)=DCMPLX(E(13,3), C(13,3))
      Eij(3,14)=DCMPLX(E(14,3), C(14,3))
      Eij(3,15)=DCMPLX(E(15,3), C(15,3))
      Eij(3,16)=DCMPLX(E(16,3), C(16,3))
      Eij(3,17)=DCMPLX(E(17,3), C(17,3) )   
      Eij(3,18)=DCMPLX(E(18,3), C(18,3))
      Eij(3,19)=DCMPLX(E(19,3), C(19,3))
      Eij(3,20)=DCMPLX(E(20,3), C(20,3))
      Eij(3,21)=DCMPLX(E(21,3), C(21,3))
      Eij(3,22)=DCMPLX(E(22,3), C(22,3))
      Eij(3,23)=DCMPLX(E(23,3), C(23,3))
      Eij(3,24)=DCMPLX(E(24,3), C(24,3))
c      write(*,*)  '************ 4 ***********'
c      write(*,*)  '************ 4 **0*********'
      Eij(4,1)=DCMPLX(E(1,4), C(1,4))
      Eij(4,2)=DCMPLX(E(2,4), C(2,4))
      Eij(4,3)=DCMPLX(E(3,4), C(3,4))
      Eij(4,4)=DCMPLX(E(4,4), C(4,4))
      Eij(4,5)=DCMPLX(E(5,4), C(5,4))
      Eij(4,6)=DCMPLX(E(6,4), C(6,4))
      Eij(4,7)=DCMPLX(E(7,4), C(7,4))    
      Eij(4,8)=DCMPLX(E(8,4), C(8,4))
      Eij(4,9)=DCMPLX(E(9,4), C(9,4))
      Eij(4,10)=DCMPLX(E(10,4), C(10,4))
      Eij(4,11)=DCMPLX(E(11,4), C(11,4))
      Eij(4,12)=DCMPLX(E(12,4), C(12,4))
      Eij(4,13)=DCMPLX(E(13,4), C(13,4))
      Eij(4,14)=DCMPLX(E(14,4), C(14,4))
      Eij(4,15)=DCMPLX(E(15,4), C(15,4))
      Eij(4,16)=DCMPLX(E(16,4), C(16,4))
      Eij(4,17)=DCMPLX(E(17,4), C(17,4))
      Eij(4,18)=DCMPLX(E(18,4), C(18,4))
      Eij(4,19)=DCMPLX(E(19,4), C(19,4))
      Eij(4,20)=DCMPLX(E(20,4), C(20,4))
      Eij(4,21)=DCMPLX(E(21,4), C(21,4))
      Eij(4,22)=DCMPLX(E(22,4), C(22,4))
      Eij(4,23)=DCMPLX(E(23,4), C(23,4))
      Eij(4,24)=DCMPLX(E(24,4), C(24,4))
      Eij(4,25)=DCMPLX(E(25,4), C(25,4))
      Eij(4,26)=DCMPLX(E(26,4), C(26,4))
      Eij(4,27)=DCMPLX(E(27,4), C(27,4))   
      Eij(4,28)=DCMPLX(E(28,4), C(28,4))
      Eij(4,29)=DCMPLX(E(29,4), C(29,4))
      Eij(4,30)=DCMPLX(E(30,4), C(30,4))
      Eij(4,31)=DCMPLX(E(31,4), C(31,4))
      Eij(4,32)=DCMPLX(E(32,4), C(32,4))
      Eij(4,33)=DCMPLX(E(33,4), C(33,4))
      Eij(4,34)=DCMPLX(E(34,4), C(34,4))
      Eij(4,35)=DCMPLX(E(35,4), C(35,4))
      Eij(4,36)=DCMPLX(E(36,4), C(36,4))
      Eij(4,37)=DCMPLX(E(37,4), C(37,4))    
      Eij(4,38)=DCMPLX(E(38,4), C(38,4))
      Eij(4,39)=DCMPLX(E(39,4), C(39,4))
      Eij(4,40)=DCMPLX(E(40,4), C(40,4))
      Eij(4,41)=DCMPLX(E(41,4), C(41,4))
      Eij(4,42)=DCMPLX(E(42,4), C(42,4))
      Eij(4,43)=DCMPLX(E(43,4), C(43,4))
      Eij(4,44)=DCMPLX(E(44,4), C(44,4))
      Eij(4,45)=DCMPLX(E(45,4), C(45,4))
      Eij(4,46)=DCMPLX(E(46,4), C(46,4))

      END                                  


c
c--------------  F1A4 -------------------------------------------------
c
      subroutine f1a4( mt, q1sq, q2sq, psq, f1, a4)
      implicit none
      double precision mt, q1sq, q2sq, psq
      complex*16 f1, a4
      complex*16 B0(3), C0
c
c determine the scalar integrals needed for the tensor T^munu describing
c g(q1)+g(q2)--> H(P) with virtual gluons and on- or offshell Higgs
c
c  in:    mt      quark mass in loop
c         q1sq    virtuality q1.q1 of first gluon
c         q2sq    virtuality q2.q2 of second gluon
c         pq2     momentum square P.P of Higgs 
c
c  out:   f1      coefficient of T_L^munu
c         a4      coefficient of T_4^munu
c
c	Dieter Zeppenfeld, <dieter@pheno.physics.wisc.edu>
c	Initial version:  2000 April 7
c	Last modified:    2000 November 12
c
c This version has been extensively checked against an older one based
c on the Taylor expansion in qi^2/mt^2, P^2/mt^2  

c  local variables
      double precision lambda, q1q2
      complex*16 J, L, Delta, b0tq1, b0tq2, b0tp, I3
c      real*8 a1, b1, c1,  den, q1sqn, q2sqn, qisq, qisqn, psqn, 
c     1       phi, beta, srt, lnfac, re, im
c      complex*16  i1q, i1p, b0tq
      complex*16 B0tM, I3point
      external B0tM, I3point

      b0tq1 = B0tM(mt,q1sq) 
      b0tq2 = B0tM(mt,q2sq)
      b0tp  = B0tM(mt,psq)
      I3 = I3point(mt,q1sq,q2sq,psq)
c      PRint*,' i3*2mt^2 =',i3*2*mt**2/(1 + (q1sq+q2sq+psq)/(12*mt**2) )
      goto 10

      entry F1A4IN( mt, q1sq, q2sq, psq, B0, C0, f1, a4)
c
c  entry F1A4IN allows to tranfer the values of the B and C function 
c  which are needed to evaluate F1 and A4
c
c  B0(1) = B0(q1sq)
c  B0(2) = B0(q2sq)
c  B0(3) = B0(psq)
c  C0 = C0(q1,q2) = -I3
c
      b0tq1 = B0(1)
      b0tq2 = B0(2)
      b0tp  = B0(3)
      I3 = -C0

 10   continue
      q1q2 = 0.5d0*(psq-q1sq-q2sq)
      lambda = q1q2**2-q1sq*q2sq
      if ( lambda.le.0d0 ) then 
         print*,' singular lambda in F1A4 ',lambda
         stop
      endif

      L = -2*q1q2*I3 - 2*b0tp + b0tq1 + b0tq2 
      Delta = (q1sq - q2sq)*(b0tq1 - b0tq2) 
      J = (4*mt**2-psq) * I3

      F1 = 2*lambda*(2+J) + L*(2*lambda - 3*q1q2*psq) - 3*q1q2*Delta
      F1 = -f1/(8*lambda**2)

      A4 = -2*q1q2*f1 +(psq*L + Delta)/(2*lambda)

      return
      end


***********************************************************************
      SUBROUTINE E01M_off(i1,i2,i3,i4,m,E051m)
*----------------------------------------------------------------------
c... author   : Michael Kubocz
c... date     : 26.08.08
c... modified : 30.01.11
c... calculates the scalar E-function for gluon-gluon scattering in
c... Higgs + 2 jet process
*----------------------------------------------------------------------
      IMPLICIT NONE

#include "ggf_dotproducts.inc"

      integer i1,i2,i3,i4
      real*8 m
      complex*16 E051m
      real * 8 p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      complex * 16 D0_2345,D0_1345,D0_1245,D0_1235,D0_1234
      
      complex*16 E01M
      complex*16 D0_1M
      EXTERNAL D0_1M,E01M

      p1=q_dot_q(i1,i1)
      p2=q_dot_q(i2,i2)
      p3=q_dot_q(i3,i3)
      p4=q_dot_q(i4,i4)
      p1p2=q_dot_q(i1,i2)
      p1p3=q_dot_q(i1,i3)
      p1p4=q_dot_q(i1,i4)
      p2p3=q_dot_q(i2,i3)
      p2p4=q_dot_q(i2,i4)
      p3p4=q_dot_q(i3,i4)

      D0_2345=D0_1M (m,p2,p3,p4,p2p3,p2p4,p3p4)
      D0_1345=D0_1M (m,p1+p2+2d0*p1p2,p3,p4,p1p3+p2p3,p1p4+p2p4,p3p4)
      D0_1245=D0_1M (m,p1,p2+p3+2d0*p2p3,p4,p1p2+p1p3,p1p4,p2p4+p3p4)
      D0_1235=D0_1M (m,p1,p2,p3+p4+2d0*p3p4,p1p2,p1p3+p1p4,p2p3+p2p4)
      D0_1234=D0_1M (m,p1,p2,p3,p1p2,p1p3,p2p3)

      E051m=-E01M(m,p1,p2,p3,p4,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,
     &            D0_2345,D0_1345,D0_1245,D0_1235,D0_1234)

      END


