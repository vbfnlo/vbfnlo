c
c---------------- BCD_fill_v(k1,k2,q1,q2) -------------------------------
c
      subroutine BCD_fill_v(k1,k2,q1,q2)
      implicit none
      double precision k1(0:3),k2(0:3),q1(0:4),q2(0:4)
c
c     q1s = 0 for  the massless gluon
c
c This subroutine does the same as BCD_fill for a single momentum
c assignment as in the diagram
c
c     ID = 1  D0(k2,q2,q1),Dij(k2,q2,q1)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q1      $  q2
c                 $         $ 
c               gluon     vector       
c 
c     ID = 2  D0(k2,q1,q2),Dij(k2,q1,q2)
c     k1  -->---------->------------- k2
c                 $         $
c                 $         $
c                 $ q2      $  q1
c                 $         $ 
c               vector     gluon      
c
c     ID = 3 D0(q1,k2,q2),Dij(q1,k2,q2)
c
c
c results are stored in Dijv(*,*,ID) etc. i.e.
c
c     renormalization scale
c
#include "VBFNLO/utilities/scales.inc"
c
Common block for the output
      complex*16 D0v(3), Dijv(3,13,3), Teps1(2), Teps2(2), Tborn(2,2),
     $     Tg1(2),Tg2(2),Ceps(2)
      common /bcd_qqvg/ D0v,Dijv,Teps1,Teps2,Tborn,Tg1,Tg2,Ceps
      double precision s, t, q1sq, q2sq,u
c
c local variables, external functions etc.
      double precision p1p2, p1p3, p2p3
      double precision musq     !renormalization scale
      complex*16 C0123(3),  C0124(3),  C0134(3),  C0234(3)
      complex*16 Cij123(2,4), Cij124(2,4), 
     1           Cij134(2,4), Cij234(2,4)
      complex*16 B0a(3), B0aq1, B0aq2, B0at, B0as,B0au
      complex*16 B0t, C0t, D01m_fin
      external B0t, C0t, D01m_fin
      double precision tmq2i,umq2i
      double precision dot0p, psumsq
      external dot0p, psumsq
      double precision pi,pi2o3
      parameter (pi=3.14159 26535 89793d0,pi2o3=pi**2/3d0)

c
      s = -2*dot0p(k1,k2)
      t = psumsq(k2,q2)
      u = psumsq(k2,q1)
      q1sq = dot0p(q1,q1)        !q1s = 0 
c      print*,'q1sq=',q1sq
      q2sq = dot0p(q2,q2)
c      print*,'q2sq=',q2sq
      musq = mursq(1,1)        !set in scales.inc
c
c now determine the required B0, C0 and D0 functions
c D0v(1) = D0t(s,t,0,q2sq,musq) = D0t(s,t,q2sq,0,musq)
c D0v(2) = D0t(s,u,0,q2sq,musq) = D0t(s,u,q2sq,0,musq)
      D0v(1) = D01m_fin(s,t,q2sq,musq)
      D0v(2) = D01m_fin(s,u,q2sq,musq)
      D0v(3) = D01m_fin(t,u,q2sq,musq)
c ID=1
c C0123 = C0(0,q2sq,t), C0124 = C0(0,s,0), C0134 = C0(t,0,0),
c C0234 = C0(0,q2sq,s)
c ID=2
c C0123 = C0(0,0,u), C0124 =C0(0,s,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,s)
c ID=3 
c  C0123 = C0(0,0,u), C0124 =C0(0,t,0) , C0134 = C0(u,q2s,0),
c C0234 = C0(0,q2s,t)
      
      C0124(1) = C0t(0d0,s,0d0,musq)
      C0234(1) = C0t(q2sq,0d0,s,musq)
      C0123(1) = C0t(0d0,q2sq,t,musq)
      C0134(1) = C0t(t,0d0,0d0,musq) !
c
      C0124(2) = C0t(0d0,s,0d0,musq) !=C0124(1)
      C0234(2) = C0t(0d0,q2sq,s,musq) !=C0234(1)
      C0123(2) = C0t(0d0,0d0,u,musq)
      C0134(2) = C0t(u,q2sq,0d0,musq)
c     
      C0124(3) = C0t(0d0,t,0d0,musq) !=C0134(1) 
      C0234(3) = C0t(0d0,q2sq,t,musq) !=C0123(1) 
      C0123(3) = C0t(0d0,0d0,u,musq) !=C0123(2)
      C0134(3) = C0t(u,q2sq,0d0,musq) !=C0134(2)    
c
c     and now the B0 functions: B0t(q2sq), B0t(t)
c     Note: B0t(s,musq=s) = 2
      B0aq1 = 0.0d0
      B0aq2 = B0t(q2sq,musq)
      B0at  = B0t(t,musq)
      B0as  = B0t(s,musq)
      B0au  = B0t(u,musq)
c
c Now everything is ready to calculate the C_ij and D_ij functions
c     ID = 1
c     p1 = k2,p2 = q2,p3 = q1
      B0a(1) = 0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0123(1),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(1),Cij124(1,1))
      B0a(1) = B0at
      B0a(2) = 0d0
      B0a(3) = 0d0
      call tens3(0d0,t,0d0,0d0,B0a,C0134(1),Cij134(1,1))
      B0a(1) = b0aq2
      B0a(2) = 0d0
      B0a(3) = B0as
      call tens3(0d0,q2sq,0d0,s,B0a,C0234(1),Cij234(1,1))
      p1p2 = dot0p(k2,q2)
      p1p3 = dot0p(k2,q1)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, q2sq, 0d0, p1p2, p1p3, p2p3, 
     &            C0123(1),  C0124(1),  C0134(1),  C0234(1),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(1), Dijv(1,1,1) )
c
c     ID = 2
c     p1 = k2,p2=q1,p3=q2 
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(2),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0as
      B0a(3) = 0d0
      call tens3(0d0,0d0,s,0d0,B0a,C0124(2),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(2),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0as
      call tens3(0d0,0d0,q2sq,s,B0a,C0234(2),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(k2,q2)
      p2p3 = dot0p(q2,q1)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(2),  C0124(2),  C0134(2),  C0234(2),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(2), Dijv(1,1,2) )
c
c      ID = 3
c     p1 = q1,p2=k2,p3=q2
c
      B0a(1) = 0
      B0a(2) = 0
      B0a(3) = B0au
      call tens3(0d0,0d0,0d0,u,B0a,C0123(3),Cij123(1,1))
      B0a(1) = 0d0
      B0a(2) = B0at
      B0a(3) = 0d0
      call tens3(0d0,0d0,t,0d0,B0a,C0124(3),Cij124(1,1))
      B0a(1) = B0au
      B0a(2) = B0aq2
      B0a(3) = 0d0
      call tens3(0d0,u,q2sq,0d0,B0a,C0134(3),Cij134(1,1))
      B0a(1) = 0d0
      B0a(2) = B0aq2
      B0a(3) = B0at
      call tens3(0d0,0d0,q2sq,t,B0a,C0234(3),Cij234(1,1))
      p1p2 = dot0p(k2,q1)
      p1p3 = dot0p(q2,q1)
      p2p3 = dot0p(k2,q2)
      call tens4( 0d0, 0d0, 0d0, q2sq, p1p2, p1p3, p2p3, 
     &            C0123(3),  C0124(3),  C0134(3),  C0234(3),
     &            Cij123(1,1), Cij124(1,1), 
     &            Cij134(1,1), Cij234(1,1), 
     &            D0v(3), Dijv(1,1,3) )
c     
c     
c and the coefficients of the tensors in M_6:
c     
c     check these formulas
      tmq2i = 1d0/(t-q2sq)
      umq2i = 1d0/(u-q2sq)
      
      Teps1(1) = (2*B0at+1)/t
      Teps2(1) = ( (B0at-b0aq2)*(2*t+3*q2sq)*tmq2i
     &     +2*b0aq2+1-2*q2sq*C0123(1) ) * tmq2i
      Tg1(1) = B0at/t
      Tg2(1) = (B0at-B0aq2)*tmq2i

      Ceps(1) = 2*(t*C0134(1)+1d0)/t
      Ceps(2) = 2*(u*C0123(2)+1d0)/u
c
      Tborn(1,1) =  ( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     +2*B0at -5d0 +pi2o3 !-5d0 + pi2o3  
c
c      print*,'C0123-c0123=',C0123(1)-C0123(1)
      Teps1(2) = (2*B0au+1)/u
      Teps2(2) = ( (B0au-b0aq2)*(2*u+3*q2sq)*umq2i
     &     +2*b0aq2+1-2*q2sq*C0134(2) ) * umq2i
      Tg1(2) = B0au/u
      Tg2(2) = (B0au-B0aq2)*umq2i
      Tborn(2,1) =  ( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     +2*B0au -5d0 + pi2o3 !-5d0 + pi2o3 
c 
c     5-pi2o3 + cvirt = -3 see expression for div piece  
c     
      Tborn(1,2) = -2*(t*C0134(1)+1) +( 2*q2sq*(B0at-b0aq2) + t*B0at
     &     -q2sq*b0aq2 )*tmq2i - 2*q2sq*C0123(1)
     &     + B0at -5d0 + pi2o3 !-5d0 + pi2o3 
c
      Tborn(2,2) = -2*(u*C0123(2)+1) +( 2*q2sq*(B0au-b0aq2) + u*B0au
     &     -q2sq*b0aq2 )*umq2i - 2*q2sq*C0134(2)
     &     + B0au -5d0 + pi2o3 !-5d0 + pi2o3
c      print*,'bcd_fill'
c      print*,'Teps1=',Teps1
c      print*,'Teps2=',Teps2
c      print*,'Tg1=',Tg1
c      print*,'Tg2=',Tg2
c      print*,'Tborn(1,1)=',Tborn(1,1)
c      print*,'Tborn(2,1)=',Tborn(2,1)
c      print*,'Tborn(1,2)=',Tborn(1,2)
c      print*,'Tborn(2,2)=',Tborn(2,2)
      return
      end
