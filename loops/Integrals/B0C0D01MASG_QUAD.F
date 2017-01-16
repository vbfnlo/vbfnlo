      function A0finG_QUAD(m,musq)
      IMPLICIT NONE
      complex*32 A0finG_QUAD
      complex*32 A0finG1_QUAD,A0finG2_QUAD
      External A0finG1_QUAD,A0finG2_QUAD
      real*16 eps
      parameter (eps=1d-7)
      real*16 m,musq
         if(m.ge.1d-7) then
            A0finG_QUAD=A0finG2_QUAD(m,musq)
         else
            A0finG_QUAD=A0finG1_QUAD(musq)
         endif
      return
      end

c   A0finG1=A0finG(0,musq)                                 
c --------------------------------------
      complex*32 function A0finG1_QUAD(musq)
      implicit none
      real*16 musq
      A0finG1_QUAD = (0.q0,0.q0)
      end

c --------------------------------------
c   A0finG2=A0finG(M,musq)                                 
c --------------------------------------
      complex*32 function A0finG2_QUAD(M,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),Ieps=(0q0,1q-78))
      real*16 M,musq
      A0finG2_QUAD = M*M*(-Log(M*M/musq)+1.q0)
      end


c --------------------------------------
c   B0finG1=B0finG(0,0,s,musq)  s<=>0                               
c --------------------------------------
      complex*32 function B0finG1_QUAD(s,musq)
c
c    0,0,s
c    0,0,0
c
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),Ieps=(0q0,1q-78))
      real*16 s,musq

      if(abs(s).lt.1d-7) then
       B0finG1_QUAD=0q0
       return
       elseif (s.gt.0.q0) then 
         B0finG1_QUAD = -(Log(s/musq)-Ipi)+2.q0
      else
         B0finG1_QUAD = -Log(-s/musq)+2.q0
      endif
      end

c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)                                 
c --------------------------------------
c     m1,0,s
c     m1,0,0
c     m1,0,m1   -> B0g3
c     0,0,s    ->  B0G1
      complex*32 function B0finG2_QUAD(M1,s,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),
     &  Ieps=(0q0,1q-78))
      real*16 M1,s,musq,msq
      complex*32 B0finG3_QUAD,B0finG1_QUAD
      external B0finG3_QUAD,B0finG1_QUAD
       msq = M1*M1
      IF(MSQ.GT.1D-7) THEN
      if (abs(s).lt.1d-10) then
        B0finG2_QUAD = -Log(msq/musq)+ 1.q0
      elseif(abs(s-msq).lt.1d-12) then
c         print*, "HEREB"
        B0finG2_QUAD = B0finG3_QUAD(M1,musq)
      elseif (s.gt.msq) then 
         B0finG2_QUAD = -Log(msq/musq)
     &          +(msq-s)/s*(Log((s-msq)/msq)-Ipi)+2.q0
      else
         B0finG2_QUAD = -Log(msq/musq)
     &          +(msq-s)/s*Log(-(s-msq)/msq)+2.q0
      endif 
      ELSE
         B0finG2_QUAD = B0finG1_QUAD(S,musq)
      ENDIF   
      end


c --------------------------------------
c   B0finG3=B0finG(M1,0,s=M^2,musq)          
c   one leg on-shell                       
c --------------------------------------
c      M1 0,M1
      complex*32 function B0finG3_QUAD(M1,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),Ieps=(0q0,1q-78))
      real*16 M1,musq,msq

      msq = M1*M1
      B0finG3_QUAD = -Log(msq/musq)+2.q0
      end


c --------------------------------------
c   B0finG3=B0finG(M1,M1,s,musq)                                 
c --------------------------------------
      complex*32 function B0finG4_QUAD(M1,s,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),Ieps=(0q0,1q-78))
      real*16 M1,s,musq,msq
      complex*32 b,g1,g2

      msq = M1*M1
      b = Sqrt(s*s-4.q0*s*(msq-Ieps))
      g1 = .5q0*(s+b)/s
      g2 = .5q0*(s-b)/s

      B0finG4_QUAD = -log((s-Ieps)/musq) 
     &        + (g1*log((g1-1.q0)/g1)-log(g1-1))
     &        + (g2*log((g2-1.q0)/g2)-log(g2-1)) + 2.q0 
      end


c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
      complex*32 function B0finG5_QUAD(M1,M2,s,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     & Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0),Ieps=(0q0,1q-78))
      real*16 M1,M2,s,musq,m1sq,m2sq
      complex*32 b,g1,g2
      complex*32 B0finG6_QUAD
      external B0finG6_QUAD

      if(abs(s).gt.1d-7) then
      m1sq = M1*M1
      m2sq = M2*M2
      b = Sqrt((s-m2sq+m1sq)**2-4.q0*s*(m1sq-Ieps))
      g1 = .5q0*(s-m2sq+m1sq+b)/s
      g2 = .5q0*(s-m2sq+m1sq-b)/s

      B0finG5_QUAD = -log((s-Ieps)/musq) 
     &        + (g1*log((g1-1.q0)/g1)-log(g1-1))
     &        + (g2*log((g2-1.q0)/g2)-log(g2-1)) + 2.q0 
      else
       if(m1.eq.m2) then
       B0finG5_QUAD = B0finG6_QUAD(M1,musq)
       else
       print*,"B:Not program yet" 
       endif
       endif

      end


c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
      complex*32 function B0finG6_QUAD(M1,musq)
c      include 'QUADmuonc0d0.inc'
      implicit none
      real*16 pi
      complex*32 Ipi, Ieps 
      parameter(pi=3.141 592 653 589 793 238 462 643 383 279 502 88q0,
     &  Ipi=(0q0,3.141 592 653 589 793 238 462 643 383 279 502 88q0)
     & ,Ieps=(0q0,1q-78))
      real*16 M1,musq,msq
      if(abs(m1).gt.1d-7) then
      msq = M1*M1
      B0finG6_QUAD = -Log(msq/musq)
      else
      B0finG6_QUAD=0q0
      endif
      end


      function B0finG_QUAD(m0,m1,p1sq,musq)
      implicit none
      complex*32 B0finG_QUAD,B0fin_QUAD
      complex*32 B0finG1_QUAD,B0finG2_QUAD,B0finG3_QUAD
      complex*32 B0finG4_QUAD,B0finG5_QUAD,B0finG6_QUAD
      complex*32 B2i1e_QUAD,B1i1e_QUAD,B1i0e_QUAD,B0t1_QUAD
      External B0finG1_QUAD,B0finG2_QUAD,B0finG3_QUAD,B0fin_QUAD
      External B0finG4_QUAD,B0finG5_QUAD,B0finG6_QUAD
      External B2i1e_QUAD,B1i1e_QUAD,B1i0e_QUAD,B0t1_QUAD
      real*16 m0,m1,p1sq,musq
      real*16 eps
c      complex*32 rslt(0:2)
      parameter (eps=1d-7)
      logical NotAVH
      parameter (NotAVH=.true.)
c      print*,"NotAVH",NotAVH
c --------------------------------------
c   B0finG2=B0finG(M1,0,s,musq)  s<=>0                               
c --------------------------------------
c   M1,0,s, 
c   M1,0,0
c   0,M,s                ?     = ?B1i1e(s,mis,mus,id)
c   0,M,0
c   B0finG2(m1,0d0,musq)     =   B1i0e(m02,musq,0)

c --------------------------------------
c   B0finG5=B0finG(M1,M2,s,musq)                                 
c --------------------------------------
c  M1,M2,s
c  M,M,s, (but there is another one for this case)
c
c   B0finG5(m0,m1,p1sq,musq) =   B2i1e(p1sq,m02,m12,musq,0)

c --------------------------------------
c   B0finG6=B0finG(M1,M1,0,musq)          
c --------------------------------------
c M,M, 0

c   B0finG6(m1,musq)         =   B0t1(m1,0d0,musq)     
   
c   M1,M2,S    v
c   M1,M2,0    not yet
c   M,M,S      v
c   M,M,0      v
c   M,0,S      v
c   M,0,0      v
c   0,M,S      v
c   0,M,0      v
C   0,0,0      v
c   0,0,s      v
  
      if(NotAVH) then
      if (m0.ge.eps) then
         if (m1.ge.eps) then
            if (abs(p1sq).ge.eps) then
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
               B0finG_QUAD= B0finG5_QUAD(m0,m1,p1sq,musq)
               return
               else
c B0finG(m0,m0,p1sq,musq)                  
               B0finG_QUAD= B0finG4_QUAD(m0,p1sq,musq)
               return
               endif
            else   ! p1s =0
c B0finG(m0,m1,0,musq)'
               if(abs((m0-m1)).gt.eps) then
c B0finG(m0,m1,p1sq,musq)
                Print*, 'No program yet'
                stop
               return
               else
c B0finG(m0,m0,0,musq)                  
               B0finG_QUAD= B0finG6_QUAD(m0,musq)
               return
               endif
               return
            endif
        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(m0,0,p1sq,musq)
c              B0finG=B0finG2(m0,p1sq,musq)
c               return
c             else
c B0finG(m0,0,0,musq)
             B0finG_QUAD= B0finG2_QUAD(m0,p1sq,musq)! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
        endif  ! m1s=0
c m0=0       
      else  
c        if (m1.ge.eps) then
c            if (abs(p1sq).ge.eps) then
c B0finG(0,m1,p1sq,musq)                  
c               B0finG= B0finG2(m1,p1sq,musq)
c               return
c               else   ! p1s =0
c B0finG(0,m1,0,musq)'
c           print*, "HERE, HERE"
               B0finG_QUAD= B0finG2_QUAD(m1,p1sq,musq) ! B0finG2(m1,0d0,musq)= B1i0e(m02,musq,0)
               return
c            endif
c        else   ! m1s=0
c            if (abs(p1sq).ge.eps) then
c B0finG(0,0,p1sq,musq)
c              B0finG=B0fin(p1sq,musq)
c             return
c             else
c B0finG(0,0,0,musq)
c             B0finG= 0d0
c             return
c            endif
c        endif  ! m1s=0
       endif   ! m0s =0 
       return
       else
         print*, "B: wrong switch"
         stop
c         call avh_olo_mu_set(sqrt(musq))
c         call avh_olo_onshell(eps)
c         call avh_olo_b0m(rslt,p1sq,m0**2,m1**2)
c         B0finG = rslt(0)
c         return
      endif
    
      end

