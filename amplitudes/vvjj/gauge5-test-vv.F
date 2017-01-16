      subroutine qqwwqq_comp(pbar,sign, nlo, 
     1                  uucc,uuss,ddcc,ddss,udsc,ducs)

      implicit none                

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
c
c and which are calculated from output of KOPPLN
c
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
#include "VBFNLO/utilities/koppln.inc"
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
c
c variables for the main part of the program
c
      double precision  pbar(0:3,4+max_v),  pbs(0:3,4+max_v), ans
      double precision uucc,uuss,ddcc,ddss,udsc,ducs

      integer  sign(4+max_v), fs(4+max_v), nlo, mu,i
      double complex zz
c assign momenta for switched quark lines
      do mu = 0,3
         pbs(mu,1) = pbar(mu,3)
         pbs(mu,2) = pbar(mu,4)
         pbs(mu,3) = pbar(mu,1)
         pbs(mu,4) = pbar(mu,2)
         do i = 5,8
            pbs(mu,i) = pbar(mu,i)
         enddo
      enddo
      fs(1) = sign(3)
      fs(2) = sign(4)
      fs(3) = sign(1)
      fs(4) = sign(2)
      do i = 5,8
         fs(i) = sign(i)
      enddo
         
#ifdef WITH_MADGRAPH
      if (.false.) then
         call sucuc(pbar,sign,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(uucc/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," uucc: ",uucc," mg2 ",ans," ratio = ",uucc/ans
c         if (abs(uucc/ans-1).gt.3d-3) read(*,"(a)") text
c      endif
         call susus(pbar,sign,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(uuss/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," uuss: ",uuss," mg2 ",ans," ratio = ",uuss/ans
c         if (abs(uuss/ans-1).gt.3d-3) read(*,"(a)") text
c      endif
         call susus(pbs,fs,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(ddcc/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," ddcc: ",ddcc," mg2 ",ans," ratio = ",ddcc/ans
c         if (abs(ddcc/ans-1).gt.3d-3) read(*,"(a)") text
c      endif
         call sdsds(pbar,sign,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(ddss/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," ddss: ",ddss," mg2 ",ans," ratio = ",ddss/ans
c         if (abs(ddss/ans-1).gt.1d-3) read(*,"(a)") text
c      endif
         call susdc(pbar,sign,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(udsc/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," udsc: ",ddss," mg2 ",ans," ratio = ",udsc/ans
!         if (abs(udsc/ans-1).gt.2d-3) read(*,"(a)") text
c      endif
         call susdc(pbs,fs,ans)
c      ans = ans*36     ! remove color/spin averaging to color/spin sum
c      if (abs(ducs/ans-1).gt.1d-3) then 
         print*," sign = ",sign
         print*," ducs: ",ddcc," mg2 ",ans," ratio = ",ducs/ans
!         if (abs(ducs/ans-1).gt.2d-3) read(*,"(a)") text
c      endif
      endif

      ip = 1
      call sucuc(pbar,sign,uucc)
      ip = 2
      call susus(pbar,sign,uuss)
      ip = 3
      call susus(pbs ,fs  ,ddcc)
      do i = 1,200
         zz = tamp(i,-1,1,ip)
         tamp(i,-1,1,ip) = tamp(i,1,-1,ip)
         tamp(i,1,-1,ip) = zz
      enddo
      ip = 4
      call sdsds(pbar,sign,ddss)
      ip = 5
      call susdc(pbar,sign,udsc)
      ip = 6
      call susdc(pbs ,fs  ,ducs)
#endif


      end


c
c----------  boxline6(id, ... ) ----------------------------------------
c
      subroutine boxline6(id, psi1,psi2,k1,k2, isig,
     1                    eps1,eps2,
     2                    mborn,    mvirt)
      implicit none
      integer id, isig
      double precision k1(0:3),k2(0:3), q1(0:4),q2(0:4)
      double complex psi1(2),psi2(2),eps1(0:5),eps2(0:5),mborn,mvirt
c 
c  determine the finite virtual corrections to the fundamental box-type
c  graph, i.e. the sum of all virtual corrections to the Born diagram
c
c    psi1(k1)    ---->---------->-------------->--  psi2(k2)
c                         S              S
c                         S              S
c                         S              S
c                     q1  V           q2 V
c                         S              S
c                         S              S
c
c                      eps1(q1)       eps2(q2)
c
c   note: k1 = k2+q1+q2, k1.k1 = 0 = k2.k2, q1.q1,q2.q2 arbitrary
c
c  INPUT:      id = 1...8      identifier for graph to determine 
c                              correct D_ij function etc., see BCD_fill
c              psi1(2)         the incoming Weyl spinor
c              psi2(2)         the outgoing Weyl spinor, psi^+(k2)
c       k1(0:3),k2(0:3)        the corresponding momenta in direction of
c                              the fermion line
c              isig            chirality of Weyl fermion (=tau in HZ notation)
c       eps1,eps2(0:3)         the "polarization vectors/currents" 
c                              corresponding to the 2 external EW bosons
c       eps1,eps2(4:5)         outflowing momenta of EW bosons =
c                              q1(0:3),q2(0:3) in DHELAS notation
c              mborn           the amplitude for the corresponding Born diagram
c
c  OUTPUT:     mvirt           the amplitude for the finite part of all 
c                              virtual 1-loop QCD corrections
c
Common block for D_ij etc. information, calculated in BCD_fill
      double complex Teps1(8), Teps2(8), Tborn(8),
     1               Tg1(8), Tg2(8)
      double precision ss(8), tt(8), q1sq(8), q2sq(8)
      double precision aq1(0:3,8), aq2(0:3,8), ak1(0:3,8), ak2(0:3,8)
c      common /bcd_qqv/ ss,tt,q1sq,q2sq, D0v,Dijv, Teps1,Teps2,Tborn,
c     1                                            Tg1,Tg2
      common /bcd_qqv_aux/ aq1,aq2,ak1,ak2
Common block bcd_qqv_aux is used for debugging purposes only
      logical ldebug, gaugecheck, loutput, lwarn
      parameter (ldebug=.false., loutput = .false., lwarn = .true.)
      character text
common block for debugging via bcd_fill_c
      double complex D0c, Dijc(3,13), Teps1c, Teps2c, Tbornc,Tg1c,Tg2c
      double precision sc, tc, q1sqc, q2sqc
      common /bcd_qqv_c/ D0c,Dijc,Teps1c,Teps2c,Tbornc,Tg1c,Tg2c, 
     1                   sc,tc,q1sqc,q2sqc
C   local variables
      integer mu
      double precision s,t,mom(0:3),dotrr
      external dotrr
      double complex me1, me2, mq, d0t, s1c, s1r, sc3, z
      external d0t, s1c, s1r, sc3
      double complex ce1,ce2,cq,cborn
      double complex d0, d11,d12,d13, d21,d22,d23,d24,d25,d26,d27
      double complex d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,
     1               d311,d312,d313
      double complex eps1k1, eps1k2, eps1q2, eps1eps2, 
     1               eps2k1, eps2q1, eps2k2, eps1q1, eps2q2,
     2               gce1, gce2, gcq

c transfer momentum information of EW bosons to q1,q2

      q1(0) = dreal(eps1(4))
      q1(1) = dreal(eps1(5))
      q1(2) = dimag(eps1(5))
      q1(3) = dimag(eps1(4))
      q1(4) = q1(0)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2(0) = dreal(eps2(4))
      q2(1) = dreal(eps2(5))
      q2(2) = dimag(eps2(5))
      q2(3) = dimag(eps2(4))
      q2(4) = q2(0)**2-q2(1)**2-q2(2)**2-q2(3)**2
      if (id.eq.-1) then
         do mu = 0,3
            q1(mu) = - q1(mu)
            if(ldebug) print*," mu = ",mu,(k2(mu)+q1(mu)+q2(mu)),(k1(mu)+1d-33)
         enddo
      elseif (id.eq.-2) then
         do mu = 0,3
            q2(mu) = - q2(mu)
            if(ldebug) print*," mu = ",mu,(k2(mu)+q1(mu)+q2(mu)),(k1(mu)+1d-33)
         enddo
      endif
      do mu = 0,3
         if (abs(k2(mu)+q1(mu)+q2(mu)-k1(mu)).gt.1d-5) then
            print*," mu = ",mu,":",k2(mu)+q1(mu)+q2(mu)-k1(mu)
            print*,k1(mu),k2(mu),q1(mu),q2(mu)
         endif
      enddo
      id = 1

      call BCD_fill_c(k1,k2,q1,q2)
      ss(id) = sc
      tt(id) = tc
      q1sq(id) = q1sqc
      q2sq(id) = q2sqc
      Teps1(id) =Teps1c
      Teps2(id) =Teps2c
      Tborn(id) =Tbornc
      Tg1(id) =  Tg1c 
      Tg2(id) =  Tg2c 

Check consistency of kinematics
      if (ldebug) then
         print*," "
         print*," comparison of input to BOXLINE for id = ",id
         print("(a,3f16.4)")," q1sq = ",q1sq(id),q1(4),q1sq(id)/q1(4)
         print("(a,3f16.4)")," q2sq = ",q2sq(id),q2(4),q2sq(id)/q2(4)
         s = -2*dotrr(k1,k2)
         t =  2*dotrr(k2,q2)+q2(4)
         print("(a,3f16.4)")," s = ",ss(id),s,ss(id)/s
         print("(a,3f16.4)")," t = ",tt(id),t,tt(id)/t
c         do mu = 0,3
c            if (abs(k1(mu)-ak1(mu,id)).gt.(1d-10*k1(0)) .or.
c     1          abs(k2(mu)-ak2(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q1(mu)-aq1(mu,id)).gt.(1d-10*k2(0)) .or.
c     1          abs(q2(mu)-aq2(mu,id)).gt.(1d-10*k2(0)) ) then
c               print*," WARNING: momentum mismatch in boxline, mu =",mu
c               print*," k1: ",k1(mu)/(ak1(mu,id)+1d-80)
c               print*," k2: ",k2(mu)/(ak2(mu,id)+1d-80)
c               print*," q1: ",q1(mu)/(aq1(mu,id)+1d-80)
c               print*," q2: ",q2(mu)/(aq2(mu,id)+1d-80)
c            endif
c         enddo
         print*," input psi1 ",psi1
         print*," input psi2 ",psi2
         print*," input eps1 ",eps1
         print*," input eps2 ",eps2
      endif

c  get matrix elements 
c
c  me1 = psibar(2) eps1slash psi(1)
c  me2 = psibar(2) eps2slash psi(1)
c  mq  = psibar(2) (q1-q2)slash psi(1)
      me1 = -s1c(psi2,eps1,.true.,isig,psi1)
      me2 = -s1c(psi2,eps2,.true.,isig,psi1)
      mq  = -2*s1r(psi2,q1,.true.,isig,psi1)
      if (ldebug) then
         print*," check M_q: ",mq/s1r(psi2,q2,.true.,isig,psi1),
     1          " = 2?"
         do mu = 0,3
            mom(mu) = k2(mu)+q2(mu)
         enddo
         t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
         z = -sc3(psi2,eps2,mom,eps1,psi1,isig)/t
c         print*," M_born: 1=?   ",mborn/z," Re(ratio) = ",
c     &                     dreal(mborn/z)
         mborn = z
      endif
      do mu = 0,3
         mom(mu) = k2(mu)+q2(mu)
      enddo
      t = mom(0)**2-mom(1)**2-mom(2)**2-mom(3)**2
      mborn = -sc3(psi2,eps2,mom,eps1,psi1,isig)/t
c
c prepare the coefficients
      s = ss(id)
      t = tt(id)

      eps1k1 = eps1(0)*k1(0)-eps1(1)*k1(1)-eps1(2)*k1(2)-eps1(3)*k1(3)
      eps1k2 = eps1(0)*k2(0)-eps1(1)*k2(1)-eps1(2)*k2(2)-eps1(3)*k2(3)
      eps1q2 = eps1(0)*q2(0)-eps1(1)*q2(1)-eps1(2)*q2(2)-eps1(3)*q2(3)
      eps1eps2 = eps1(0)*eps2(0)-eps1(1)*eps2(1)-
     1           eps1(2)*eps2(2)-eps1(3)*eps2(3)
      eps2k1 = eps2(0)*k1(0)-eps2(1)*k1(1)-eps2(2)*k1(2)-eps2(3)*k1(3)
      eps2k2 = eps2(0)*k2(0)-eps2(1)*k2(1)-eps2(2)*k2(2)-eps2(3)*k2(3)
      eps2q1 = eps2(0)*q1(0)-eps2(1)*q1(1)-eps2(2)*q1(2)-eps2(3)*q1(3)

      eps1q1 = eps1(0)*q1(0)-eps1(1)*q1(1)-eps1(2)*q1(2)-eps1(3)*q1(3)
      eps2q2 = eps2(0)*q2(0)-eps2(1)*q2(1)-eps2(2)*q2(2)-eps2(3)*q2(3)

      gaugecheck = abs(eps1q1).gt.1d-10*abs(eps1q2) .or.
     1             abs(eps2q2).gt.1d-10*abs(eps2q1)
      if (ldebug) print*," gaugecheck = ",gaugecheck
      gaugecheck = .true. !gaugecheck .or. ldebug
      
      d0 = d0c
      
      d11 = dijc(1,1)
      d12 = dijc(1,2)
      d13 = dijc(1,3)
         
      d21 = dijc(2,1)
      d22 = dijc(2,2)
      d23 = dijc(2,3)
      d24 = dijc(2,4)
      d25 = dijc(2,5)
      d26 = dijc(2,6)
      d27 = dijc(2,7)

      d31 = dijc(3,1)
      d32 = dijc(3,2)
      d33 = dijc(3,3)
      d34 = dijc(3,4)
      d35 = dijc(3,5)
      d36 = dijc(3,6)
      d37 = dijc(3,7)
      d38 = dijc(3,8)
      d39 = dijc(3,9)
      d310 = dijc(3,10)
      d311 = dijc(3,11)
      d312 = dijc(3,12)
      d313 = dijc(3,13)

      if (ldebug) then
         call BCD_fill_c(k1,k2,q1,q2)
         print*," D0: ",abs(d0/d0c-1)
         print("(a,3g10.2)")," D1j ",abs(d11/dijc(1,1)-1),
     1                               abs(d12/dijc(1,2)-1),
     2                               abs(d13/dijc(1,3)-1)
         print("(a,7g10.2)")," D2j ",abs(d21/dijc(2,1)-1),
     1                               abs(d22/dijc(2,2)-1),
     1                               abs(d23/dijc(2,3)-1),
     1                               abs(d24/dijc(2,4)-1),
     1                               abs(d25/dijc(2,5)-1),
     1                               abs(d26/dijc(2,6)-1),
     1                               abs(d27/dijc(2,7)-1)
         print("(a,7g10.2)")," D3j ",abs(d31/dijc(3,1)-1),
     1                               abs(d32/dijc(3,2)-1),
     1                               abs(d33/dijc(3,3)-1),
     1                               abs(d34/dijc(3,4)-1),
     1                               abs(d35/dijc(3,5)-1),
     1                               abs(d36/dijc(3,6)-1),
     1                               abs(d37/dijc(3,7)-1)
         print("(a,7g10.2)"),"     ",abs(d38/dijc(3,8)-1),
     1                               abs(d39/dijc(3,9)-1),
     1                               abs(d310/dijc(3,10)-1),
     1                               abs(d311/dijc(3,11)-1),
     1                               abs(d312/dijc(3,12)-1),
     1                               abs(d313/dijc(3,13)-1)
      endif
c
      ce1 = -2*(d12+d24)*eps2k2*(q1sq(id)+q2sq(id)-3*s-4*t) - 
     1      2*(d12-d24)*eps2q1*(q2sq(id)-t) 
     2  +4*(- d11*eps2k2*s - d12*eps2k1*t + d13*eps2k2*
     3  (q2sq(id)-s-t) + d13*eps2q1*q2sq(id) - d21*eps2k2*s 
     4 - d22*eps2k2*t - d22*eps2q1*q2sq(id) + d23*eps2k2*q1sq(id)
     5 + d25*eps2k2*(q2sq(id)-s-2*t) 
     6 - d26*eps2k2*(q2sq(id)-s-t)+d26*eps2q1*t+2*d27*eps2q1
     7 - d32*eps2k2*q2sq(id) - d34*eps2k2*(q2sq(id)-t) 
     8 + d35*eps2k2*(q2sq(id)-s-t)+d36*eps2k2*(2*q2sq(id)-t) 
     9 + d37*eps2k2*q1sq(id) + d38*eps2k2*(q1sq(id)+q2sq(id)-s)
     & - d39*eps2k2*q1sq(id) 
     1 - d310*eps2k2*(q1sq(id)+2*q2sq(id)-2*s-t)
     2 - 4*d311*eps2k2 + 6*d312*eps2k2 + 2*d313*eps2q1)

      ce2= -2*(d12+d24)*(eps1k2*(q1sq(id)+q2sq(id)-s-2*t)
     1  +eps1q2*(q2sq(id)-s-3*t)) + 
     2  4*(d13*eps1k2*(q1sq(id)-2*s-t) - 
     3  d13*eps1q2*(2*s+t) + d22*eps1k2*q2sq(id)
     4  + d22*eps1q2*q2sq(id) - d23*eps1k2*t + 
     5    d23*eps1q2*(q1sq(id)-t) - d24*eps1k1*q2sq(id) 
     6  + d25*eps1k2*q1sq(id) + d25*eps1k1*(q2sq(id)-2*s-t) 
     7  + d26*eps1k2*t - d26*eps1k1*(q1sq(id)-s) - 
     8  2*d27*eps1q2 + d33*eps1k2*q1sq(id) + d33*eps1q2*q1sq(id)
     9  + d37*eps1k1*(q2sq(id)-s-t) + d38*eps1k1*q2sq(id)
     &  - d39*eps1k1*(q1sq(id)+q2sq(id)-s) - d310*eps1k1*(q2sq(id)-t)
     1  +2*d311*eps1k2+2*d312*eps1q2-6*d313*eps1k1)

C NOTE: cq = 1/2 * cq1 in my reduce code
      cq = (d12+d24)*eps1eps2*s + 2*(4*d12*eps2k2*eps1k2 + 
     1 3*d12*eps2k2*eps1q2+d12*eps2q1*eps1k2-4*d13*eps2k2*eps1k2-
     2 2*d13*eps2k2*eps1q2-2*d13*eps2q1*eps1k2-d13*eps1eps2*s+
     3 2*d22*eps2k2*eps1q2-d22*eps1eps2*t-2*d23*eps2q1*eps1k2-
     4 2*d23*eps2q1*eps1q2-d23*eps1eps2*t+6*d24*eps2k2*eps1k2+
     5 3*d24*eps2k2*eps1q2+d24*eps2q1*eps1k2-6*d25*eps2k2*eps1k2-
     6 2*d25*eps2k2*eps1q2-2*d25*eps2q1*eps1k2-d25*eps1eps2*s-
     7 4*d26*eps2k2*eps1q2+4*d26*eps2q1*eps1k2+2*d26*eps2q1*eps1q2+ 
     8 d26*eps1eps2*(s+2*t)-d32*eps1eps2*q2sq(id)+
     9 d33*eps1eps2*q1sq(id)+2*d34*eps2k2*eps1k2-2*d35*eps2k2*eps1k2+
     & 2*d36*eps2k2*eps1q2+d36*eps1eps2*(q2sq(id)-t)-
     1 2*d37*eps2q1*eps1k2+d37*eps1eps2*(q2sq(id)-s-t)+
     2 2*d38*eps2q1*eps1q2+d38*eps1eps2*(q1sq(id)+2*q2sq(id)-s) -
     3 2*d39*eps2q1*eps1q2-d39*eps1eps2*(2*q1sq(id)+q2sq(id)-s) -
     4 2*d310*eps2k2*eps1q2+2*d310*eps2q1*eps1k2-d310*eps1eps2*
     5 (2*q2sq(id)-s-2*t)+4*d312*eps1eps2-4*d313*eps1eps2)

      cborn = -2*((d36+d37-2*d310)*(q2sq(id)-t) + 
     1    d38*(q1sq(id)+2*q2sq(id)) - d39*(2*q1sq(id)+q2sq(id)) )
     2    -2*(d0+d11+d12-2*d13+d24-d25+d26-d37-d38+d39+d310)*s
     3    -2*(-(d22+d23-2*d26)*t + 2*d27 - d32*q2sq(id)
     4    +d33*q1sq(id)+ 6*(d312-d313))

      if (gaugecheck) then
         gce1 = -2*eps2q2*( d12*(2*q1sq(id)+q2sq(id)-s-3*t) - 
     1       2*(d13-d22)*(q1sq(id)+q2sq(id)-s-t)
     2     - d24*(q2sq(id)-s-t)-2*d26*q1sq(id)-4*d27-4*d312)

         gce2 = 2*eps1q1*(d12*(q2sq(id)+t)-2*d13*(q2sq(id)+s) + 
     1      2*d22*q2sq(id) + 2*d23*(q1sq(id)+q2sq(id)-s-2*t)
     2     +d24*(q2sq(id)+t)+2*d25*(-q2sq(id)+s+t)-2*d26*(q2sq(id)-t)
     3     -4*d27+2*d33*q1sq(id) + 4*d313)

         gcq = 4*eps1q1*eps2q2*(d12-d13+d22-d23+d38-d39) +
     1        2*eps1q1*eps2k2*(d12-4*d23+d24+2*d26-2*d37+2*d310) -
     2        4*eps1q1*eps2q1*(d23-d26+d33-d39) +
     3        2*eps2q2*eps1k2*(3*d12-4*d13+4*d22+d24-2*d25-
     4                                     2*d26+2*d36-2*d310) +
     5        4*eps2q2*eps1q2*(d12-d13+2*d22-2*d26+d32-d38)
         ce1 = ce1 + gce1
         ce2 = ce2 + gce2
         cq = cq + gcq
      endif
c
c  now add the terms from the vertex and propagator corrections
      ce1 = ce1 + 2*eps2k2*Teps2(id)
      ce2 = ce2 + 2*eps1k1*Teps1(id)
      if (gaugecheck) then
         ce1 = ce1 + 2*eps2q2*Tg2(id)
         ce2 = ce2 - 2*eps1q1*Tg1(id)
      endif
      cborn = cborn*t - Tborn(id)
c  and the final result
      mvirt = me1*ce1 + me2*ce2 + mq*cq + mborn*cborn
      if (ldebug) then
         print*," mvirt/mborn ",mvirt/mborn
         if (abs(eps1(0)/q1(0)-1).lt.1d-10) then
            print*," mborn type ",(mborn*cborn+ me2*ce2)/mborn
            print*," mq type    ",(me1*ce1+ mq*cq)/mborn
         else
            print*," mborn type ",(mborn*cborn+ me1*ce1)/mborn
            print*," mq type    ",(me2*ce2+ mq*cq)/mborn
         endif
         read(*,"(a)") text
      endif

      return
      end

c---------------------------------------------------------------------

      subroutine gauge_check_5(j5,qin,qsign,psi,p,nq1,nq2,
     1                         jal,jbe,cvirtc)
      implicit none
      integer qsign,nq1,nq2
      double complex j5(0:3,3), psi(2,-1:1,4), jal(0:5),jbe(0:5)
      double precision q(0:3), qal(0:3), qbe(0:3), qin(0:3), p(0:3,4)
c
      double complex mv5(3), mvc(3), mb(3), mbc(3), braket(2,2), jx(0:5)
      integer mu,j,id
      double precision momsum, dummy(0:4), cvirtc
      double complex s1c,dotrc
      external s1c, dotrc
c
	logical bad_gauge,bad_gauge_sin,bad_gauge_isig
	common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
	
	bad_gauge_sin = .false. ! set F at beginning of each run
	
      do mu = 0,3
         q(mu) = qsign*qin(mu)
      enddo
c momenta are outgoing
      qal(0) = dreal(jal(4))
      qal(1) = dreal(jal(5))
      qal(2) = dimag(jal(5))
      qal(3) = dimag(jal(4))
      qbe(0) = dreal(jbe(4))
      qbe(1) = dreal(jbe(5))
      qbe(2) = dimag(jbe(5))
      qbe(3) = dimag(jbe(4))
c compute Born amplitudes contracted with momentum instead of other quark line      
      call ket2c(psi(1,-1,nq1),.true.,p(0,nq1),-1,
     1        qal,jal,braket(1,1),dummy)
      call bra2c(psi(1,-1,nq2),.true.,p(0,nq2),-1,
     1        qbe,jbe,braket(1,2),dummy)
      mb(3) = s1c(psi(1,-1,nq2),jbe,.true.,-1,braket(1,1))
      mb(1) = -s1c(braket(1,2),jal,.true.,-1,psi(1,-1,nq1))
      mb(2) = -mb(1)-mb(3)
c
      do mu = 0,3
         jx(mu) = jal(mu)
      enddo
      jx(4) = jal(4)+dcmplx(q(0),q(3))
      jx(5) = jal(5)+dcmplx(q(1),q(2))
      id = 1
      call boxline6(id, 
     1     psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),-1,
     2     jx, jbe,
     3     mbc(1),   mvc(1))
c      mvc(1) = - mvc(1)     not needed: boxline6 alread return -amplitude
c      mbc(1) = - mbc(1)

      do mu = 0,3
         jx(mu) = jbe(mu)
      enddo
      jx(4) = jbe(4)+dcmplx(q(0),q(3))
      jx(5) = jbe(5)+dcmplx(q(1),q(2))
      id = 1
      call boxline6(id, 
     1     psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),-1,
     2     jal, jx,
     3     mbc(3),   mvc(3))
      mvc(3) = - mvc(3)
      mbc(3) = - mbc(3)

      mvc(2) = - mvc(1) - mvc(3)
      mbc(2) = - mbc(1) - mbc(3)

      if (.false.) then
         print*," momentum check in gauge_check_5 "
         do mu = 0,3
            momsum = p(mu,nq1)-qal(mu)-qbe(mu)-q(mu)-p(mu,nq2)
            print*," mu = ",mu," momsum = ",momsum
         enddo
         do j = 1,3
            print*," j = ",j," tree ratio = ",mbc(j)/mb(j)
         enddo
      endif

c      print*," gauge check ratios "
      do j = 1,3
         mv5(j) = dotrc(q,j5(0,j))-cvirtc*mb(j)

c         print*," j = ",j," virtual ward ratio = ",mv5(j)/mvc(j)

	 if( abs((abs(mv5(j)/mvc(j))-1d0)).gt.1e-1 ) then
		bad_gauge_sin = .true.
	 endif
         if ((mv5(j).ne.mv5(j)).or.
     &       (mvc(j).ne.mvc(j))) bad_gauge_sin = .true.
      enddo
c      print*," virtual sum ",mv5(1)+mv5(2)+mv5(3),mvc(1),mvc(2),mvc(3)

      end
c---------------------------------------------------------------------

      subroutine gauge_check_isig_5(j5,qin,qsign,psi,p,nq1,nq2,
     1                         jal,jbe,cvirtc)
      implicit none
      integer qsign,nq1,nq2
      double complex j5(0:3,-1:1,3), psi(2,-1:1,4), jal(0:5),jbe(0:5)
      double precision q(0:3), qal(0:3), qbe(0:3), qin(0:3), p(0:3,4)
c
      double complex mv5(3,-1:1),mvc(3,-1:1),mb(3,-1:1),mbc(3,-1:1),
     1 		     braket(2,-1:1,2),jx(0:5)
      integer mu,j,id,isig
      double precision momsum, dummy(0:4), cvirtc
      double complex s1c,dotrc
      external s1c, dotrc

	logical bad_gauge,bad_gauge_sin,bad_gauge_isig
	common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
	
	bad_gauge_isig = .false. ! set F at beginning of each run
	
      do mu = 0,3
         q(mu) = qsign*qin(mu)
      enddo
c momenta are outgoing
      qal(0) = dreal(jal(4))
      qal(1) = dreal(jal(5))
      qal(2) = dimag(jal(5))
      qal(3) = dimag(jal(4))
      qbe(0) = dreal(jbe(4))
      qbe(1) = dreal(jbe(5))
      qbe(2) = dimag(jbe(5))
      qbe(3) = dimag(jbe(4))
c compute Born amplitudes contracted with momentum instead of other quark line      
      
      do isig = -1,1,2
      
      call ket2c(psi(1,isig,nq1),.true.,p(0,nq1),isig,
     1        qal,jal,braket(1,isig,1),dummy)
      call bra2c(psi(1,isig,nq2),.true.,p(0,nq2),isig,
     1        qbe,jbe,braket(1,isig,2),dummy)
      mb(3,isig) =  s1c(psi(1,isig,nq2),jbe,.true.,isig,braket(1,isig,1))
      mb(1,isig) = -s1c(braket(1,isig,2),jal,.true.,isig,psi(1,isig,nq1))
      mb(2,isig) = -mb(1,isig)-mb(3,isig)
      
      enddo !isig
c
      do mu = 0,3
         jx(mu) = jal(mu)
      enddo
      jx(4) = jal(4)+dcmplx(q(0),q(3))
      jx(5) = jal(5)+dcmplx(q(1),q(2))
      id = 1
      do isig = -1,1,2
      call boxline6(id, 
     1     psi(1,isig,nq1),psi(1,isig,nq2),p(0,nq1),p(0,nq2),isig,
     2     jx, jbe,
     3     mbc(1,isig),   mvc(1,isig))
c      mvc(1) = - mvc(1)     not needed: boxline6 alread return -amplitude
c      mbc(1) = - mbc(1)
      enddo !isig

      do mu = 0,3
         jx(mu) = jbe(mu)
      enddo
      jx(4) = jbe(4)+dcmplx(q(0),q(3))
      jx(5) = jbe(5)+dcmplx(q(1),q(2))
      id = 1
      do isig = -1,1,2
      call boxline6(id, 
     1     psi(1,isig,nq1),psi(1,isig,nq2),p(0,nq1),p(0,nq2),isig,
     2     jal, jx,
     3     mbc(3,isig),   mvc(3,isig))
      mvc(3,isig) = - mvc(3,isig)
      mbc(3,isig) = - mbc(3,isig)

      mvc(2,isig) = - mvc(1,isig) - mvc(3,isig)
      mbc(2,isig) = - mbc(1,isig) - mbc(3,isig)
      enddo !isig

      if (.false.) then
         print*," momentum check in gauge_check_5 "
         do mu = 0,3
            momsum = p(mu,nq1)-qal(mu)-qbe(mu)-q(mu)-p(mu,nq2)
            print*," mu = ",mu," momsum = ",momsum
         enddo
         do j = 1,3
            print*," j = ",j," tree ratio for isig = -1: ",mbc(j,-1)/mb(j,-1)
         enddo
      endif

c      print*," gauge check ratios "
      do isig = -1,1,2
      do j = 1,3
         mv5(j,isig) = dotrc(q,j5(0,isig,j))-cvirtc*mb(j,isig)
c         print*," j = ",j," isig = ",isig,
c     #		" virtual ward ratio = ",mv5(j,isig)/mvc(j,isig)

 	 if( abs((abs(mv5(j,isig)/mvc(j,isig))-1d0)).gt.1e0 ) then
		bad_gauge_isig = .true.
 	 endif
         if ((mv5(j,isig).ne.mv5(j,isig)).or.
     &       (mvc(j,isig).ne.mvc(j,isig))) bad_gauge_isig = .true.
   
      enddo
      enddo !isig

      end
     

c-------------------------------------------------------------------------------

      subroutine gauge_check0_5(j5,qin,qsign,psi,p,nq1,nq2,
     1                          jal,jbe,cvirtc)
      implicit none
      integer qsign,nq1,nq2
      double complex j5(0:3), psi(2,-1:1,4), jal(0:5),jbe(0:5)
      double precision q(0:3), qal(0:3), qbe(0:3), qin(0:3), p(0:3,4)
c
      double complex mv5(3), mvc(3), mb(3), mbc(3), braket(2,2), jx(0:5)
      integer mu,j,id
      double precision momsum, dummy(0:4), cvirtc
      double complex s1c,dotrc
      external s1c, dotrc
c
      logical bad_gauge,bad_gauge_sin,bad_gauge_isig
      common /gauge / bad_gauge,bad_gauge_sin,bad_gauge_isig
      
      ! test only:


      double complex dotcc
      external dotcc
      

*************************************************************************

      bad_gauge_sin = .false. ! set F at beginning of each run
	
      do mu = 0,3
         q(mu) = qsign*qin(mu)
      enddo
c momenta are outgoing
      qal(0) = dreal(jal(4))
      qal(1) = dreal(jal(5))
      qal(2) = dimag(jal(5))
      qal(3) = dimag(jal(4))
      qbe(0) = dreal(jbe(4))
      qbe(1) = dreal(jbe(5))
      qbe(2) = dimag(jbe(5))
      qbe(3) = dimag(jbe(4))
c compute Born amplitudes contracted with momentum instead of other quark line      
      call ket2c(psi(1,-1,nq1),.true.,p(0,nq1),-1,
     1        qal,jal,braket(1,1),dummy)
      call bra2c(psi(1,-1,nq2),.true.,p(0,nq2),-1,
     1        qbe,jbe,braket(1,2),dummy)
      mb(3) = s1c(psi(1,-1,nq2),jbe,.true.,-1,braket(1,1))
      mb(1) = -s1c(braket(1,2),jal,.true.,-1,psi(1,-1,nq1))
      mb(2) = -mb(1)-mb(3)

      do mu = 0,3
         jx(mu) = jal(mu)
      enddo
      jx(4) = jal(4)+dcmplx(q(0),q(3))
      jx(5) = jal(5)+dcmplx(q(1),q(2))
      id = 1
      call boxline6(id, 
     1     psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),-1,
     2     jx, jbe,
     3     mbc(1),   mvc(1))

   
      do mu = 0,3
         jx(mu) = jbe(mu)
      enddo
      jx(4) = jbe(4)+dcmplx(q(0),q(3))
      jx(5) = jbe(5)+dcmplx(q(1),q(2))
      id = 1
      call boxline6(id, 
     1     psi(1,-1,nq1),psi(1,-1,nq2),p(0,nq1),p(0,nq2),-1,
     2     jal, jx,
     3     mbc(3),   mvc(3))
 
      mvc(3) = - mvc(3)
      mbc(3) = - mbc(3)

      mvc(2) = - mvc(1) - mvc(3)
      mbc(2) = - mbc(1) - mbc(3)
      
      if (.false.) then
         print*,' momentum check in gauge_check_5 '
         do mu = 0,3
            momsum = p(mu,nq1)-qal(mu)-qbe(mu)-q(mu)-p(mu,nq2)
            print*,' mu = ',mu,' momsum = ',momsum
         enddo
         j=2
         print*,' j = ',j,' tree ratio = ',mbc(j)/mb(j)
      endif

         mv5(2) = dotrc(q,j5(0))-cvirtc*mb(2)
c	 print*
	 if( (abs((abs(mv5(2)/mvc(2))-1d0)).gt.1e-1).or.
     &	     ((abs((abs(dreal(mv5(2))/dreal(mvc(2))))-1d0)).gt.1e-1).or.
     &	     ((abs((abs(dimag(mv5(2))/dimag(mvc(2))))-1d0)).gt.1e-1)) then
     
		bad_gauge_sin = .true.
c		print*,'mv5,mvc:',mv5(2),mvc(2)
	 endif
         if ((mv5(2).ne.mv5(2)).or.
     &       (mvc(2).ne.mvc(2))) bad_gauge_sin = .true.
c	 print*,'bad_gauge=',bad_gauge_sin
c	 print*,'re-im:',dreal(mv5(2))/dreal(mvc(2)),dimag(mv5(2))/dimag(mvc(2))
c	 print*,'rat-abs:',abs(mv5(2))/abs(mvc(2))
c	 print*,'mv5/mvc:',mv5(2)/mvc(2),abs(mv5(2)/mvc(2))
c	 print*
	
	 
      end
