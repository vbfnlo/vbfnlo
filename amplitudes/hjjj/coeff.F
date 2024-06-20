      SUBROUTINE boxline_vg(psi1,psi2,
     $     k1,k2,isig,eps1,eps2,
     $     q1,q2,mborn,mvirt)
c     So far the syntax in correct.
      implicit none
c     arguments
      integer isig
      real*8 k1(0:3),k2(0:3),q1(0:4),q2(0:4)
      complex*16 psi1(2),psi2(2),eps1(0:3),eps2(0:3)
      complex*16 mborn(2),mvirt(4)
c     local variables
      integer mu,k
      real*8 mom(0:3)
      real*8 s,t,u,q2s,q1s
      complex*16 ce1(3),ce2(3),cq(3),cborn(4)
      complex*16 d0,d11,d12,d13,d21,d22,d23,d24,d25,d26,d27
      complex*16 d31,d32,d33,d34,d35,d36,d37,d38,d39,d310
      complex*16 d311,d312,d313
      complex*16 eps1k2,eps1q2,eps1q1,eps2k2
      complex*16 eps2q2,eps2q1,eps1eps2,eps2k1,eps1k1
      complex*16 ucborn,tcborn
      complex*16 me1,me2,mq,s1c,s1r,sc3,z
      external s1c,s1r,sc3
c     Common block for the output
      complex*16 D0v(3), Dijv(3,13,3), Teps1(2), Teps2(2), Tborn(2,2)
      complex*16 Tg1(2),Tg2(2),Ceps(2)
      common /bcd_qqvg/ D0v,Dijv,Teps1,Teps2,Tborn,Tg1,Tg2,Ceps
c     
      complex*16 dotcr,dotcc
      double precision dot0p,psumsq
      external dot0p,psumsq,dotcr,dotcc
      logical ldebug,lskip,lgtree,ltriangles
      parameter(ltriangles =.false.)
      parameter(ldebug = .false.,lskip=.false.,lgtree=.false.)
      double precision colfac(2),ncol
      parameter(ncol = 3d0)	!number of colors
      complex*16 czero
      parameter(czero = (0d0,0d0))
c     
c     
c     print*,'inside coeff'
c     
      colfac(1) = -1d0/(2d0*ncol)
      colfac(2) = ncol/2d0
c     
c     
c     born processes
c     
c     k1---->--->----->--k2
c            g       v
c            g       v
c            g       v
c          q1       q2
c       the amplitude for this is mborn(1)=mborn1,
c       k1---->--->----->--k2
c            v       g
c            v       g
c            v       g
c          q2        q1
c       the amplitude for this is mborn(2)=mborn2,
c
c       define coefficeints as an array to innumerate the 
c       three box graphs
c       
c       box graph 1:   graphid =1
c       k2 ---<------vvvvvvvv q2
c               g    |
c               g    |
c       k1 ->--------gggggggg q1
c       box graph 2:  graphid =2
c       k2 ---<------gggggggg q1
c              g    |
c              g    |
c       k1 ->--------vvvvvvv q2 
c       
c       box 1 and box 2 have color factor (CF-1/2*CA)
c       m_virt(box 1 + box2) = mvirt(1)+mvirt(2)
c       box graph 3: graphid = 3
c       
c       q2 vvvvvvvv-------->--- k2
c                  |    g
c                  |    g
c       k1 --->----gggggggggg q1
c       
c       box 3 has color factor 1/2*CA 
c       m_virt(box 3) = mvirt(3)
c       For the sum of all graphs: graphid = 4
c       m_virt(total) = mvirt(4) =(CF-1/2CA)mvirt(1+2)+1/2CA mvirt(3)
c
      me1=czero
      me2=czero
      mq = czero
      eps1k2=czero
      eps1q2=czero
      eps1q1=czero
      eps2k2=czero
      eps2q2=czero
      eps2q1=czero
      eps1eps2=czero
      eps2k1=czero
      eps1k1=czero
      do k=1,3
         ce1(k)=czero
         ce2(k)=czero
         cq(k)=czero
         cborn(k)=czero
      enddo
      cborn(4)=czero
c     One should also include triagle graphs and the 
c     propagator correction to the fermion.
c     get matrix elements
c       
c     me1 = psibar(2) eps1slash psi(1)
c     me2 = psibar(2) eps2slash psi(1)
c     mq =  psibar(2) (q1-q2)slash  psi(1)
      me1 = s1c(psi2,eps1,.true.,isig,psi1)
      me2 = s1c(psi2,eps2,.true.,isig,psi1)
      mq  = 2*s1r(psi2,q1,.true.,isig,psi1)
c     born-like ampitudes 
c     mborn(1) = psibar(2) eps2slash (k2+q2)slash eps1slash psi(1)/t
c     mborn(2) = psibar(2) eps1slash (k2+q1)slash eps2slash psi(1)/u
      
c     
      if (ldebug) then
         do mu =0,3
            print*,'k1-k2-q1-q2,mu=',mu,k1(mu)-k2(mu)-q1(mu)-q2(mu)
         enddo
         print*,'check M_q:',mq/s1r(psi2,q2,.true.,isig,psi1),' = 2?'
         do mu = 0,3
            mom(mu)=k2(mu)+q2(mu)
         enddo
         t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
         z = sc3(psi2,eps2,mom,eps1,psi1,isig)/t
         print*,'M_born 1:    ',mborn(1)/z,'=1? (ratio)=',(mborn(1)/z)
         do mu = 0,3
            mom(mu)=k2(mu)+q1(mu)
         enddo
         u = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
         z = sc3(psi2,eps1,mom,eps2,psi1,isig)/u
         print*,'M_born 2:    ',mborn(2)/z,'=1? (ratio)=',(mborn(2)/z)
c     print*,'mborn1/mborn2=',mborn(1)/mborn(2)
      endif
      
      
c     prepare coefficients
c     need to compute s,t,u
      s = -2*dot0p(k1,k2)
      t = psumsq(k2,q2)
      u = psumsq(k2,q1)
      q1s = dot0p(q1,q1)
      q2s = dot0p(q2,q2)
      if(ldebug) then
         print*,'(k2+q1)^2=',u,'=?',q1s+q2s-s-t
         print*,'(k2+q2)^2=',t
         print*,'k1^2=',dot0p(k1,k1)
         print*,'k2^2=',dot0p(k2,k2)
         print*,'s=',s
c     if(s.ge.0d0) pause
      endif
      
      eps1k2 = dotcr(eps1,k2)
      eps1k1 = dotcr(eps1,k1)
      eps1q1 = dotcr(eps1,q1)
      eps1q2 = dotcr(eps1,q2)
      eps2k2 = dotcr(eps2,k2)
      eps2k1 = dotcr(eps2,k1)
      eps2q1 = dotcr(eps2,q1)
      eps2q2 = dotcr(eps2,q2)
      eps1eps2 = dotcc(eps1,eps2)
c     
c     
      if(ldebug) then
         print*,'eps1q1=',eps1q1
      endif
c     compute dij, this is called once
c     call BCD_fill_v(k1,k2,q1,q2)
c     
      d0 = D0v(1)
c     
      d11 = Dijv(1,1,1)
      d12 = Dijv(1,2,1)
      d13 = Dijv(1,3,1)
c     
      d21 = Dijv(2,1,1)
      d22 = Dijv(2,2,1)
      d23 = Dijv(2,3,1)
      d24 = Dijv(2,4,1)
      d25 = Dijv(2,5,1)
      d26 = Dijv(2,6,1)
      d27 = Dijv(2,7,1)
c     
      d31 = Dijv(3,1,1)
      d32 = Dijv(3,2,1)
      d33 = Dijv(3,3,1)
      d34 = Dijv(3,4,1)
      d35 = Dijv(3,5,1)
      d36 = Dijv(3,6,1)
      d37 = Dijv(3,7,1)
      d38 = Dijv(3,8,1)
      d39 = Dijv(3,9,1)
      d310 = Dijv(3,10,1)
      d311 = Dijv(3,11,1)
      d312 = Dijv(3,12,1)
      d313 = Dijv(3,13,1)
c
c
c     box graph 1
      
      ce1(1) = -(eps2q2*(-8*d27-8*d312-(d11-d12+d13-4*d22+4*d24)*q2s+d11*s- 
     $     d12*s+d13*s-4*d22*s+4*d24*s+d11*t-d12*t+d13*t-4*d22*t+ 
     $     4*d24*t))+eps2q1*(8*d27+8*d313+(d11-d12+d13-4*d22+4*d24)*q2s-
     $     d11*t-3*d12*t+3*d13*t-4*d24*t+4*d26*t)-
     $     eps2k2*(16*d311-24*d312-(d11-d12+d13+4*d25-4*d26-8*d310-
     $     4*d32-4*d34+4*d35+8*d36+4*d38)*q2s+5*d11*s-5*d12*s+d13*s+
     $     4*d21*s-4*d24*s+4*d25*s-4*d26*s-8*d310*s+4*d35*s+4*d38*s-
     $     4*d12*t+4*d13*t+4*d22*t-8*d24*t+8*d25*t-4*d26*t-4*d310*t-
     $     4*d34*t+4*d35*t+4*d36*t)
c     
      ce2(1) =eps1k2*(8*d311-24*d313-(d11+3*d12-3*d13-4*d22+8*d24-4*d25+
     $     4*d310-4*d37-4*d38+4*d39)*q2s-d11*s+d12*s-5*d13*s-
     $     8*d25*s+4*d26*s-4*d37*s+4*d39*s+4*d12*t-4*d13*t-4*d23*t+
     $     4*d24*t-4*d25*t+4*d26*t+4*d310*t-4*d37*t)-
     $     eps1q2*(8*d27-8*d312+24*d313+(d11+3*d12-3*d13-4*d22+8*d24-
     $     4*d25+4*d310-4*d37-4*d38+4*d39)*q2s+d11*s-d12*s+5*d13*s+
     $     8*d25*s-4*d26*s+4*d37*s-4*d39*s-d11*t-7*d12*t+7*d13*t+
     $     4*d23*t-8*d24*t+4*d25*t-4*d310*t+4*d37*t)
c     cq = 1/2*cq1 in mathematica code
      cq(1) =8*d312*eps1eps2-8*d313*eps1eps2+8*d12*eps1k2*eps2k2-
     $     8*d13*eps1k2*eps2k2+12*d24*eps1k2*eps2k2-12*d25*eps1k2*eps2k2+
     $     4*d34*eps1k2*eps2k2-4*d35*eps1k2*eps2k2+d11*eps1q2*eps2k2+
     $     7*d12*eps1q2*eps2k2-7*d13*eps1q2*eps2k2+4*d22*eps1q2*eps2k2+
     $     8*d24*eps1q2*eps2k2-4*d25*eps1q2*eps2k2-8*d26*eps1q2*eps2k2-
     $     4*d310*eps1q2*eps2k2+4*d36*eps1q2*eps2k2-d11*eps1k2*eps2q1+
     $     d12*eps1k2*eps2q1-d13*eps1k2*eps2q1-4*d23*eps1k2*eps2q1-
     $     4*d25*eps1k2*eps2q1+8*d26*eps1k2*eps2q1+4*d310*eps1k2*eps2q1-
     $     4*d37*eps1k2*eps2q1-4*d23*eps1q2*eps2q1+4*d26*eps1q2*eps2q1+
     $     4*d38*eps1q2*eps2q1-4*d39*eps1q2*eps2q1-d11*eps1k2*eps2q2+
     $     5*d12*eps1k2*eps2q2-5*d13*eps1k2*eps2q2+8*d22*eps1k2*eps2q2-
     $     4*d25*eps1k2*eps2q2-4*d26*eps1k2*eps2q2-4*d310*eps1k2*eps2q2+
     $     4*d36*eps1k2*eps2q2+4*d12*eps1q2*eps2q2-4*d13*eps1q2*eps2q2+
     $     8*d22*eps1q2*eps2q2-8*d26*eps1q2*eps2q2+4*d32*eps1q2*eps2q2-
     $     4*d38*eps1q2*eps2q2-4*d310*eps1eps2*q2s-2*d32*eps1eps2*q2s+
     $       2*d36*eps1eps2*q2s+2*d37*eps1eps2*q2s+4*d38*eps1eps2*q2s-
     $     2*d39*eps1eps2*q2s-(d11*eps1eps2*s)/2+(d12*eps1eps2*s)/2-
     $     (d13*eps1eps2*s)/2-2*d25*eps1eps2*s+2*d26*eps1eps2*s+
     $     2*d310*eps1eps2*s-2*d37*eps1eps2*s-2*d38*eps1eps2*s+2*d39*eps1eps2*s-
     $     2*d22*eps1eps2*t-2*d23*eps1eps2*t+4*d26*eps1eps2*t+4*d310*eps1eps2*t-
     $     2*d36*eps1eps2*t-2*d37*eps1eps2*t
c     
      cborn(1) =-4*d27-12*d312+12*d313+4*d310*q2s+2*d32*q2s-2*d36*q2s-2*d37*q2s-
     $     4*d38*q2s+2*d39*q2s-2*d0*s-d11*s-d12*s+d13*s+2*d25*s-2*d26*s-
     $     2*d310*s+2*d37*s+2*d38*s-2*d39*s+2*d22*t+2*d23*t-4*d26*t-
     $     4*d310*t+2*d36*t+2*d37*t
c     
c     only triagle graphs
      if(ltriangles) then
         ce1(1)=czero
         ce2(1)=czero
         cq(1) = czero
         cborn(1)=czero
      endif
c     include vertex and propagator corrections
c     
      ce1(1) = ce1(1) + 2*eps2k2*Teps2(1)-2*eps2q2*Tg2(1)
      ce2(1) = ce2(1) + 2*eps1k1*Teps1(1)
      
      if(ltriangles) then
         cborn(1)=-Tborn(1,1)
      else
         cborn(1) = cborn(1)*t - Tborn(1,1)
      endif
c     
c     u-channel diagram
      
      d0 = D0v(2)
c     
      d11 = Dijv(1,1,2)
      d12 = Dijv(1,2,2)
      d13 = Dijv(1,3,2)
c     
      d21 = Dijv(2,1,2)
      d22 = Dijv(2,2,2)
      d23 = Dijv(2,3,2)
      d24 = Dijv(2,4,2)
      d25 = Dijv(2,5,2)
      d26 = Dijv(2,6,2)
      d27 = Dijv(2,7,2)
c     
      d31 = Dijv(3,1,2)
      d32 = Dijv(3,2,2)
      d33 = Dijv(3,3,2)
      d34 = Dijv(3,4,2)
      d35 = Dijv(3,5,2)
      d36 = Dijv(3,6,2)
      d37 = Dijv(3,7,2)
      d38 = Dijv(3,8,2)
      d39 = Dijv(3,9,2)
      d310 = Dijv(3,10,2)
      d311 = Dijv(3,11,2)
      d312 = Dijv(3,12,2)
      d313 = Dijv(3,13,2)
c     
c     box graph 2
      ce1(2)=-(eps2q1*(8*d27-8*d312+24*d313-4*(d23-d26+d33-d39)*q2s+d11*s-
     $     d12*s+5*d13*s+8*d25*s-4*d26*s+4*d37*s-4*d39*s-d11*u-
     $     7*d12*u+7*d13*u+4*d23*u-8*d24*u+4*d25*u-4*d310*u+4*d37*u))-
     $     eps2q2*(8*d27+16*d313-4*(d23-d26+d33-d39)*q2s+4*d13*s+
     $     4*d23*s+4*d25*s-4*d26*s+4*d37*s-4*d39*s-d11*u-3*d12*u+
     $     3*d13*u+8*d23*u-4*d24*u-4*d26*u-4*d310*u+4*d37*u)-
     $     eps2k2*(-8*d311+24*d313-(d11-d12+d13+4*d25-4*d26+4*d33-
     $       4*d39)*q2s+d11*s-d12*s+5*d13*s+8*d25*s-4*d26*s+4*d37*s-
     $     4*d39*s-4*d12*u+4*d13*u+4*d23*u-4*d24*u+4*d25*u-4*d26*u-
     $     4*d310*u+4*d37*u)
c     
      ce2(2)=-(eps1q2*(-8*d27-8*d313+(d11+3*d12-3*d13+4*d24-4*d26)*u))-
     $     eps1k2*(16*d311-24*d312+(d11+3*d12-3*d13-4*d23+4*d24+4*d310-
     $     4*d37-4*d38+4*d39)*q2s+5*d11*s-5*d12*s+d13*s+4*d21*s-
     $     4*d24*s+4*d25*s-4*d26*s-8*d310*s+4*d35*s+4*d38*s-4*d12*u+
     $     4*d13*u+4*d22*u-8*d24*u+8*d25*u-4*d26*u-4*d310*u-4*d34*u+
     $     4*d35*u+4*d36*u)
c     
      cq(2)=-8*d312*eps1eps2+8*d313*eps1eps2-8*d12*eps1k2*eps2k2+
     $       8*d13*eps1k2*eps2k2-12*d24*eps1k2*eps2k2+12*d25*eps1k2*eps2k2-
     $     4*d34*eps1k2*eps2k2+4*d35*eps1k2*eps2k2+d11*eps1q2*eps2k2-
     $     d12*eps1q2*eps2k2+d13*eps1q2*eps2k2+4*d23*eps1q2*eps2k2+
     $     4*d25*eps1q2*eps2k2-8*d26*eps1q2*eps2k2-4*d310*eps1q2*eps2k2+
     $     4*d37*eps1q2*eps2k2-d11*eps1k2*eps2q1-7*d12*eps1k2*eps2q1+
     $     7*d13*eps1k2*eps2q1-4*d22*eps1k2*eps2q1-8*d24*eps1k2*eps2q1+
     $     4*d25*eps1k2*eps2q1+8*d26*eps1k2*eps2q1+4*d310*eps1k2*eps2q1-
     $     4*d36*eps1k2*eps2q1+4*d23*eps1q2*eps2q1-4*d26*eps1q2*eps2q1-
     $     4*d38*eps1q2*eps2q1+4*d39*eps1q2*eps2q1-d11*eps1k2*eps2q2-
     $     3*d12*eps1k2*eps2q2+3*d13*eps1k2*eps2q2+8*d23*eps1k2*eps2q2-
     $     4*d24*eps1k2*eps2q2-4*d26*eps1k2*eps2q2-4*d310*eps1k2*eps2q2+
     $     4*d37*eps1k2*eps2q2+4*d23*eps1q2*eps2q2-4*d26*eps1q2*eps2q2+
     $     4*d33*eps1q2*eps2q2-4*d39*eps1q2*eps2q2-2*d33*eps1eps2*q2s-
     $     2*d38*eps1eps2*q2s+4*d39*eps1eps2*q2s+(d11*eps1eps2*s)/2-
     $     (d12*eps1eps2*s)/2+(d13*eps1eps2*s)/2+2*d25*eps1eps2*s-
     $     2*d26*eps1eps2*s-2*d310*eps1eps2*s+2*d37*eps1eps2*s+2*d38*eps1eps2*s-
     $     2*d39*eps1eps2*s+2*d22*eps1eps2*u+2*d23*eps1eps2*u-4*d26*eps1eps2*u-
     $     4*d310*eps1eps2*u+2*d36*eps1eps2*u+2*d37*eps1eps2*u
c     
      cborn(2) =-4*d27-12*d312+12*d313-2*d33*q2s-2*d38*q2s+4*d39*q2s-2*d0*s-
     $     d11*s-d12*s+d13*s+2*d25*s-2*d26*s-2*d310*s+2*d37*s+2*d38*s-
     $     2*d39*s+2*d22*u+2*d23*u-4*d26*u-4*d310*u+2*d36*u+2*d37*u
c     only triagle graphs
      if(ltriangles) then
         ce1(2)=czero
         ce2(2)=czero
         cq(2) = czero
         cborn(2)=czero
      endif
c     include vertex and propagator corrections
c     
      ce1(2) = ce1(2) + 2*eps2k1*Teps2(2) + 2*eps2q2*Tg2(2)
      ce2(2) = ce2(2) + 2*eps1k2*Teps1(2)
      if(ltriangles) then
         cborn(2)=-Tborn(2,1)
      else
         cborn(2) = cborn(2)*u - Tborn(2,1)
      endif
c     
c     
      if(lskip) goto 20
c     
      d0 = d0v(3)
      d11 = Dijv(1,1,3)
      d12 = Dijv(1,2,3)
      d13 = Dijv(1,3,3)
c     
      d21 = Dijv(2,1,3)
      d22 = Dijv(2,2,3)
      d23 = Dijv(2,3,3)
      d24 = Dijv(2,4,3)
      d25 = Dijv(2,5,3)
      d26 = Dijv(2,6,3)
      d27 = Dijv(2,7,3)
c     
      d31 = Dijv(3,1,3)
      d32 = Dijv(3,2,3)
      d33 = Dijv(3,3,3)
      d34 = Dijv(3,4,3)
      d35 = Dijv(3,5,3)
      d36 = Dijv(3,6,3)
      d37 = Dijv(3,7,3)
      d38 = Dijv(3,8,3)
      d39 = Dijv(3,9,3)
      d310 = Dijv(3,10,3)
      d311 = Dijv(3,11,3)
      d312 = Dijv(3,12,3)
      d313 = Dijv(3,13,3)
c     
c     box graph 3
      ce1(3)=24*d27*eps2k2+20*d312*eps2k2+22*d27*eps2q1+20*d311*eps2q1+
     $     12*d27*eps2q2+20*d313*eps2q2-4*d23*eps2k2*q2s+4*d26*eps2k2*q2s+
     $     2*d38*eps2k2*q2s-2*d39*eps2k2*q2s+d12*eps2q1*q2s-d13*eps2q1*q2s-
     $     3*d23*eps2q1*q2s+d24*eps2q1*q2s+2*d26*eps2q1*q2s+2*d310*eps2q1*q2s-
     $     2*d37*eps2q1*q2s-2*d23*eps2q2*q2s+2*d26*eps2q2*q2s-2*d33*eps2q2*q2s+
     $     2*d39*eps2q2*q2s+(3*d0*eps2k2*t)/2+(3*d12*eps2k2*t)/2+4*d13*eps2k2*t+
     $     2*d25*eps2k2*t+2*d26*eps2k2*t+2*d310*eps2k2*t-2*d38*eps2k2*t+
     $     2*d0*eps2q1*t+3*d11*eps2q1*t-d12*eps2q1*t+4*d13*eps2q1*t+
     $     d21*eps2q1*t-d24*eps2q1*t+6*d25*eps2q1*t-2*d26*eps2q1*t-
     $     2*d310*eps2q1*t+2*d35*eps2q1*t+4*d13*eps2q2*t+4*d23*eps2q2*t+
     $     2*d25*eps2q2*t-2*d26*eps2q2*t+2*d37*eps2q2*t-2*d39*eps2q2*t+
     $     (3*d0*eps2k2*u)/2+(7*d12*eps2k2*u)/2-2*d13*eps2k2*u+2*d22*eps2k2*u-
     $     2*d24*eps2k2*u+2*d25*eps2k2*u-2*d26*eps2k2*u+2*d310*eps2k2*u-
     $     2*d36*eps2k2*u+(3*d0*eps2q1*u)/2+d11*eps2q1*u+(5*d12*eps2q1*u)/2-
     $     2*d13*eps2q1*u-d21*eps2q1*u+d24*eps2q1*u-2*d34*eps2q1*u+
     $     2*d35*eps2q1*u+(3*d0*eps2q2*u)/2+(3*d12*eps2q2*u)/2-2*d23*eps2q2*u+
     $     2*d26*eps2q2*u-2*d310*eps2q2*u+2*d37*eps2q2*u
c     
      ce2(3)=-12*d27*eps1k2-4*d312*eps1k2-6*d27*eps1q2-4*d313*eps1q2-
     $     3*d0*eps1k2*q2s-7*d12*eps1k2*q2s+2*d13*eps1k2*q2s-2*d22*eps1k2*q2s+
     $     6*d23*eps1k2*q2s-2*d24*eps1k2*q2s+4*d25*eps1k2*q2s-8*d26*eps1k2*q2s-
     $     2*d38*eps1k2*q2s+2*d39*eps1k2*q2s-(3*d0*eps1q2*q2s)/2-
     $     (5*d12*eps1q2*q2s)/2-d13*eps1q2*q2s+3*d23*eps1q2*q2s-d24*eps1q2*q2s+
     $     2*d25*eps1q2*q2s-6*d26*eps1q2*q2s+2*d33*eps1q2*q2s-2*d39*eps1q2*q2s+
     $     (3*d0*eps1k2*t)/2-4*d11*eps1k2*t+(11*d12*eps1k2*t)/2-2*d21*eps1k2*t+
     $     2*d22*eps1k2*t-6*d25*eps1k2*t+6*d26*eps1k2*t-2*d310*eps1k2*t+
     $     2*d38*eps1k2*t-(d0*eps1q2*t)/2-3*d11*eps1q2*t+(5*d12*eps1q2*t)/2-
     $     d21*eps1q2*t+d24*eps1q2*t-6*d25*eps1q2*t+6*d26*eps1q2*t-
     $     2*d37*eps1q2*t+2*d39*eps1q2*t+(7*d0*eps1k2*u)/2+(23*d12*eps1k2*u)/2-
     $     6*d13*eps1k2*u-2*d21*eps1k2*u+2*d22*eps1k2*u+10*d24*eps1k2*u-
     $     6*d25*eps1k2*u-2*d26*eps1k2*u-2*d310*eps1k2*u+2*d36*eps1k2*u+
     $     2*d0*eps1q2*u-d11*eps1q2*u+5*d12*eps1q2*u-d21*eps1q2*u-
     $     2*d23*eps1q2*u+3*d24*eps1q2*u+2*d26*eps1q2*u+2*d310*eps1q2*u-
     $     2*d37*eps1q2*u
c     
      cq(3)=-(d27*eps1eps2)-2*d311*eps1eps2+2*d313*eps1eps2-8*d12*eps1k2*eps2k2+
     $     8*d13*eps1k2*eps2k2-4*d22*eps1k2*eps2k2-8*d24*eps1k2*eps2k2+
     $     12*d26*eps1k2*eps2k2-4*d36*eps1k2*eps2k2+4*d38*eps1k2*eps2k2+
     $     (3*d0*eps1q2*eps2k2)/2+(3*d12*eps1q2*eps2k2)/2+4*d23*eps1q2*eps2k2-
     $     8*d25*eps1q2*eps2k2+4*d26*eps1q2*eps2k2-4*d310*eps1q2*eps2k2+
     $     4*d39*eps1q2*eps2k2-(3*d0*eps1k2*eps2q1)/2-(19*d12*eps1k2*eps2q1)/2+
     $     8*d13*eps1k2*eps2q1-12*d24*eps1k2*eps2q1+8*d25*eps1k2*eps2q1+
     $     4*d26*eps1k2*eps2q1+4*d310*eps1k2*eps2q1-4*d34*eps1k2*eps2q1+
     $     4*d23*eps1q2*eps2q1-4*d25*eps1q2*eps2q1-4*d35*eps1q2*eps2q1+
     $     4*d37*eps1q2*eps2q1-(3*d0*eps1k2*eps2q2)/2-(11*d12*eps1k2*eps2q2)/2+
     $     4*d13*eps1k2*eps2q2+8*d23*eps1k2*eps2q2-4*d24*eps1k2*eps2q2-
     $     4*d26*eps1k2*eps2q2-4*d310*eps1k2*eps2q2+4*d39*eps1k2*eps2q2+
     $     4*d23*eps1q2*eps2q2-4*d25*eps1q2*eps2q2+4*d33*eps1q2*eps2q2-
     $     4*d37*eps1q2*eps2q2-(d12*eps1eps2*q2s)/2+(d13*eps1eps2*q2s)/2+
     $     (d23*eps1eps2*q2s)/2-(d24*eps1eps2*q2s)/2-d310*eps1eps2*q2s-
     $     d33*eps1eps2*q2s+d37*eps1eps2*q2s+d39*eps1eps2*q2s-d0*eps1eps2*t-
     $     (3*d11*eps1eps2*t)/2+(d12*eps1eps2*t)/2-(d21*eps1eps2*t)/2+
     $     (d24*eps1eps2*t)/2+d310*eps1eps2*t-d35*eps1eps2*t+d37*eps1eps2*t-
     $     d39*eps1eps2*t-(d11*eps1eps2*u)/2+(3*d12*eps1eps2*u)/2-
     $     d13*eps1eps2*u+(d21*eps1eps2*u)/2+d23*eps1eps2*u+
     $     (3*d24*eps1eps2*u)/2-2*d25*eps1eps2*u-d26*eps1eps2*u-d310*eps1eps2*u+
     $     d34*eps1eps2*u-d35*eps1eps2*u+d37*eps1eps2*u
c     
      cborn(3)=(6*d27+(3*d0*q2s)/2+(5*d12*q2s)/2-d13*q2s-3*d23*q2s+d24*q2s-
     $     2*d25*q2s+4*d26*q2s+(d0*t)/2+3*d11*t-(5*d12*t)/2+d21*t-d24*t+
     $     4*d25*t-4*d26*t-(3*d0*u)/2-d11*u-(9*d12*u)/2+4*d13*u+d21*u-
     $     5*d24*u+4*d25*u)/2
c     vertex and propagator corrections
c       Teps2(1) = Teps(q2sq,t),     Tg2(1) = Tg(q2sq,t)
c       Teps2(2) = Teps(q2sq,u),     Tg2(2) = Tg(q2sq,u)
c       Ceps(1) = 2*(t*C0t(t)+1)/t
c       Ceps(2) = 2*(t*C0t(u)+1)/u
c       Tborn(1,2) = Tborn(q2sq,t) = B0t(t)-t*Ceps(1) + T_b(q2sq,t)
c       Tborn(2,2) = Tborn(q2sq,u) = B0t(u)-t*Ceps(2) + T_b(q2sq,u)
c
      ce1(3) = ce1(3)-2*eps2k2*Teps2(1) + 2*eps2q2 * Tg2(1)
     $     -2*eps2k1*Teps2(2) - 2*eps2q2 * Tg2(2)
c     
      ce2(3) = ce2(3) + Ceps(1)*eps1k1 + Ceps(2)*eps1k2
      tcborn = cborn(3)*t +Tborn(1,2)
      ucborn = cborn(3)*u +Tborn(2,2)
c     
      cborn(3) = tcborn
      cborn(4) = ucborn
      
 20   continue
c     now add the terms from the vertex and propagator corrections
c     
c     to check tree level gauge invariance
      if(lgtree) then
         cborn(1)=czero
         cborn(2)=czero
         cborn(3)=czero
         cborn(4)=czero
         do k=1,3
            ce1(k)=czero
            ce2(k)=czero
            cq(k)=czero
         enddo
      endif
c     and the final result
c     box 1,box 2
      do k=1,2
         mvirt(k) = me1*ce1(k) + me2*ce2(k) +
     $        mq*cq(k) + cborn(k)*mborn(k)
      enddo
c     
c     box 3
      if(lskip) then
         mvirt(3) = czero
      else
         mvirt(3) = -(me1*ce1(3)+me2*ce2(3) + mq*cq(3) +
     $          cborn(3)*mborn(1) + cborn(4)*mborn(2))
      endif
c     Total result with color factors
c     colfac(1) = CF-1/2CA  
c     colfac(2) = 1/2 CA  
      
      mvirt(4) = colfac(1)*(mvirt(1)+mvirt(2))+colfac(2)*mvirt(3)
      
c     
      if(ldebug) print*,'here'
      return
      end 
