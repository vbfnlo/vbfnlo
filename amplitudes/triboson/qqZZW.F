c****************subroutine QQZZW *************************************
c
c	Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008, January
c	Last modified: 2008, February
C
C  QQZZW calculates the matrix elements**2 for 
C        q1 q2    ---->   Z Z W,   
C         Z ---> f1-bar f2, Z ---> f3-bar f4, W ---> f5-bar f6
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
C  This code is modified to allow for virtual corrections, more precisely
C  the interference of Born with the finite part of virtual diagrams
C  for 
c
c  INPUT:  NLO = 1       return uuww = |M_born|^2 + 2Re(M_Born^* M_virt)
c          NLO = 0       return uuww = |M_born|^2   etc.
c  for steering the calculation of finite box and pentagon contributions
c  the following additional options are implemented
c          NLO = +4      set all finite box and pentagon contributions to 0
c                        but include cvirt*M_Born contribution
c          NLO = +5      true finite pentagons = 0, all other virtual
c                        contributions included, including leading q^mu/m_V
c                        terms of pentagons which can be related to boxes
c       remaining cases:   return uuww = 2Re(M_born^* M_virt). more precisely
c       the piece which is missing for the full NLO=1 when calling +|NLO|     
c          NLO = -4      finite Box + q^mu/m_V terms of Pentagon contributions
c          NLO = -5      true finite Pentagon contributions only
c
c          N = 1         Momentum configuration 1 (LO)
c          N = 2         Momentum configuration 2: emission of parton 1
c          N = 3         Momentum configuration 3: emission of parton 2
c
c*********************************************************************

      subroutine qqZZW(pbar,sign, nlo, N, qqWZZ, bos, NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"


c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

      integer bos
      double precision  pbar(0:3,2+n_v)
      double precision  qqWZZ, res(2),resv(2)
      double precision  p(0:3,2+n_v), pwz1(0:4), pwz2(0:4),
     1                  fw(0:4,2),fqze(0:4,2), fqzmu(0:4,2)
      integer  sign(2+n_v), nlo, mu, i, k, isig, N
      integer  u(2), d(2)
      double complex psize(2,-1:1,2), psiae(2,-1:1,2),
     1               psizmu(2,-1:1,2), psiamu(2,-1:1,2),
     2               psiw(2,2)
      double complex jae(0:5,-1:1,2), jam(0:5,-1:1,2),
     1               jze(0:5,-1:1,2),jzm(0:5,-1:1,2),
     2               jw(0:5,-1:1,2) 

      double complex temp, temp1
      double complex mat(2,-1:1,10), matv(2,-1:1,10), mborn
      double complex madgr(10,-1:1,2)
      double complex mm(2), mv(2), zero
      double complex m1, m2, m3, m4, m5, m6, m7, m8
      double complex mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8
      double complex psi(2,-1:1,2), jqq(0:5,-1:1)
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save u, k, d
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons
      double complex zm2i(2:3)
      save zm2i

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebug, ldebugm, linit
      integer wardidtest1, wardidtest2, wardidtest3
      integer wardidtest4, wardidtest5, wardidtest6
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)
      integer pentcount, pentcount2
      data pentcount /0/
      data pentcount2 /0/
      save pentcount, pentcount2
      integer NW !NW = 1 means Narrow width approximaion
                
c variables for pacos implementation of the pentagon routine
      double complex j5z1z2(3), j5z2z1(3), j5a1a2(3), j5a2a1(3)
      double complex j5a1z2(3), j5a2z1(3), j5z1a2(3), j5z2a1(3)
      double complex jbz1z2(3), jbz2z1(3), jba1a2(3), jba2a1(3)
      double complex jba1z2(3), jba2z1(3), jbz1a2(3), jbz2a1(3)
      double complex test(3), test2(3)
      double precision minus_qz1(0:3),minus_qz2(0:3),minus_qw(0:3)
      double precision minus_qzz(0:3),minus_qz1w(0:3),minus_qz2w(0:3)
      double precision p1(0:3), p2(0:3) 
      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)

c variables for pacos implementation of the box routine
      double complex gaugebox1(2), gaugebox1_born(2), boxtest(2), treetest(2)
      double complex gaugebox2(2), gaugebox2_born(2)
      double complex gaugebox3(2), gaugebox3_born(2)
      double complex gaugebox4(2), gaugebox4_born(2)
      double complex result(-1:1,8),result_born(-1:1,8)
      double complex cte

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt .or. lpq

c define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
         if (bos.eq.223) then
            k = 1
            u(1) = 3
            d(1) = 4
         else if (bos.eq.224) then
            k = 2
            u(2) = 4
            d(2) = 3
         endif
      endif

      do i = 1,8
         mat(k,-1,i) = (0.0d0,0.0d0)
         matv(k,-1,i) = (0.0d0,0.0d0)
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwz1(mu) = qwbos(mu,N) + qze(mu,N)
         pwz2(mu) = qwbos(mu,N) + qzmu(mu,N)
      enddo
      pwz1(4) = pwz1(0)**2-pwz1(1)**2-pwz1(2)**2-pwz1(3)**2
      pwz2(4) = pwz2(0)**2-pwz2(1)**2-pwz2(2)**2-pwz2(3)**2

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qwbos(0,N),wbos(1,N),
     &              psiw(1,2),fw(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qwbos(0,N),wbos(1,N),
     &              psiw(1,1),fw(0,1))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qze(0,N),
     &              ze(1,N),psize(1,-1,2),fqze(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qzmu(0,N),
     &              zmu(1,N),psizmu(1,-1,2),fqzmu(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qze(0,N),
     &              ae(1,N),psiae(1,-1,2),fqze(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qzmu(0,N),
     &              amu(1,N),psiamu(1,-1,2),fqzmu(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qze(0,N),
     &              ze(1,N), psize(1,-1,1),fqze(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qzmu(0,N),
     &              zmu(1,N), psizmu(1,-1,1),fqzmu(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qze(0,N),
     &              ae(1,N), psiae(1,-1,1),fqze(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qzmu(0,N),
     &              amu(1,N), psiamu(1,-1,1),fqzmu(0,1))

**************************************************************
*
* Calculation of the born matrix element
*
**************************************************************

c Calculation of processes qq -> Z Z W **********************

c   q qbar -> Z Z W
      m1=-s1c(psize(1,-1,2),wbos(1,N),.true.,-1,psizmu(1,-1,1))
     1    *clr(d(k),2,-1)*clr(u(k),2,-1)
     2   -s1c(psizmu(1,-1,2),wbos(1,N),.true.,-1,psize(1,-1,1))
     2    *clr(d(k),2,-1)*clr(u(k),2,-1)
     3   -s1c(psiw(1,2),ze(1,N),.true.,-1,psizmu(1,-1,1))
     3    *clr(u(k),2,-1)**2
     4   -s1c(psiw(1,2),zmu(1,N),.true.,-1,psize(1,-1,1))
     4    *clr(u(k),2,-1)**2
     5   -s1c(psize(1,-1,2),zmu(1,N),.true.,-1,psiw(1,1))
     5    *clr(d(k),2,-1)**2
     6   -s1c(psizmu(1,-1,2),ze(1,N),.true.,-1,psiw(1,1))
     6    *clr(d(k),2,-1)**2

      if (NW.eq.0) then

      m2=-s1c(psiae(1,-1,2),wbos(1,N),.true.,-1,psizmu(1,-1,1))
     1    *clr(d(k),1,-1)*clr(u(k),2,-1)
     2   -s1c(psizmu(1,-1,2),wbos(1,N),.true.,-1,psiae(1,-1,1))
     2    *clr(d(k),2,-1)*clr(u(k),1,-1)
     3   -s1c(psiw(1,2),ae(1,N),.true.,-1,psizmu(1,-1,1))
     3    *clr(u(k),2,-1)*clr(u(k),1,-1)
     4   -s1c(psiw(1,2),zmu(1,N),.true.,-1,psiae(1,-1,1))
     4    *clr(u(k),2,-1)*clr(u(k),1,-1)
     5   -s1c(psiae(1,-1,2),zmu(1,N),.true.,-1,psiw(1,1))
     5    *clr(d(k),2,-1)*clr(d(k),1,-1)
     6   -s1c(psizmu(1,-1,2),ae(1,N),.true.,-1,psiw(1,1))
     6    *clr(d(k),2,-1)*clr(d(k),1,-1)

      m3=-s1c(psize(1,-1,2),wbos(1,N),.true.,-1,psiamu(1,-1,1))
     1    *clr(u(k),1,-1)*clr(d(k),2,-1)
     2   -s1c(psiamu(1,-1,2),wbos(1,N),.true.,-1,psize(1,-1,1))
     2    *clr(d(k),1,-1)*clr(u(k),2,-1)
     3   -s1c(psiw(1,2),ze(1,N),.true.,-1,psiamu(1,-1,1))
     3    *clr(u(k),1,-1)*clr(u(k),2,-1)
     4   -s1c(psiw(1,2),amu(1,N),.true.,-1,psize(1,-1,1))
     4    *clr(u(k),1,-1)*clr(u(k),2,-1)
     5   -s1c(psize(1,-1,2),amu(1,N),.true.,-1,psiw(1,1))
     5    *clr(d(k),1,-1)*clr(d(k),2,-1)
     6   -s1c(psiamu(1,-1,2),ze(1,N),.true.,-1,psiw(1,1))
     6    *clr(d(k),1,-1)*clr(d(k),2,-1)

      m4=-s1c(psiae(1,-1,2),wbos(1,N),.true.,-1,psiamu(1,-1,1))
     1    *clr(u(k),1,-1)*clr(d(k),1,-1)
     2   -s1c(psiamu(1,-1,2),wbos(1,N),.true.,-1,psiae(1,-1,1))
     2    *clr(d(k),1,-1)*clr(u(k),1,-1)
     3   -s1c(psiw(1,2),ae(1,N),.true.,-1,psiamu(1,-1,1))
     3    *clr(u(k),1,-1)**2
     4   -s1c(psiw(1,2),amu(1,N),.true.,-1,psiae(1,-1,1))
     4    *clr(u(k),1,-1)**2
     5   -s1c(psiae(1,-1,2),amu(1,N),.true.,-1,psiw(1,1))
     5    *clr(d(k),1,-1)**2
     6   -s1c(psiamu(1,-1,2),ae(1,N),.true.,-1,psiw(1,1))
     6    *clr(d(k),1,-1)**2

      else if (NW.eq.1) then

         m2 = (0.0d0, 0.0d0)
         m3 = (0.0d0, 0.0d0)
         m4 = (0.0d0, 0.0d0)

      endif ! NW

      mat(k,-1,1) = (m1 + m2 + m3 + m4) * clr(3,3,-1)


c Calculation of processes qq -> V V -> 2l + 4l ****************

c calculate currents needed for W A/Ztozz diagrams
      call curr6(-1, psi(1,-1,2), p(0,2),
     &               psiw(1,1), fw(0,1), jw(0,-1,1) )
      call curr6(-1, psiw(1,2), fw(0,2),
     &               psi(1,-1,1), p(0,1), jw(0,-1,2) )

c qq -> W A/ZtoZZ - contract with atozz/ztozz tensor
      m1 = -dotcc(azezm(0,N),jw(0,-1,1))
      m2 = -dotcc(zzezm(0,N),jw(0,-1,1))
      m3 = -dotcc(azezm(0,N),jw(0,-1,2))
      m4 = -dotcc(zzezm(0,N),jw(0,-1,2))

      mat(k,-1,2) = ( m1*clr(d(k),1,-1) + m2*clr(d(k),2,-1) )
     1                               *clr(3,3,-1)!qq -> W A/ZtoZZ
      mat(k,-1,3) = ( m3*clr(u(k),1,-1) + m4*clr(u(k),2,-1) )
     1                               *clr(3,3,-1)!qq -> A/ZtoZZ W

c calculate currents needed for A/Z WtoZW diagrams
      if (NW.eq.0) then
      call curr6(-1,psiae(1,-1,2),fqze(0,2),
     &           psi(1,-1,1),p(0,1), jae(0,-1,2) )
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiae(1,-1,1),fqze(0,1), jae(0,-1,1) )
      call curr6(-1,psiamu(1,-1,2),fqzmu(0,2),
     &           psi(1,-1,1),p(0,1), jam(0,-1,2) )
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiamu(1,-1,1),fqzmu(0,1), jam(0,-1,1) )
      endif

      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psizmu(1,-1,1),fqzmu(0,1), jzm(0,-1,1) )
      call curr6(-1,psize(1,-1,2),fqze(0,2),
     &           psi(1,-1,1),p(0,1), jze(0,-1,2) )
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psize(1,-1,1),fqze(0,1), jze(0,-1,1) )
      call curr6(-1,psizmu(1,-1,2),fqzmu(0,2),
     &           psi(1,-1,1),p(0,1), jzm(0,-1,2) )

c contract with wtozw tensors
      if (NW.eq.0) then
         m1 = -dotcc(wzw2(0,N),jae(0,-1,2))
         m3 = -dotcc(wzw2(0,N),jae(0,-1,1))
         m5 = -dotcc(wzw1(0,N),jam(0,-1,2))
         m7 = -dotcc(wzw1(0,N),jam(0,-1,1))
c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m1 = m1+zm2i(3)*dotrc(pwz2(0),wzw2(0,N))*dotrc(pwz2(0),jae(0,-1,2))
      m3 = m3+zm2i(3)*dotrc(pwz2(0),wzw2(0,N))*dotrc(pwz2(0),jae(0,-1,1))
      m5 = m5+zm2i(3)*dotrc(pwz1(0),wzw1(0,N))*dotrc(pwz1(0),jam(0,-1,2))
      m7 = m7+zm2i(3)*dotrc(pwz1(0),wzw1(0,N))*dotrc(pwz1(0),jam(0,-1,1))
      else if (NW.eq.1) then
         m1 = (0.0d0,0.0d0)
         m3 = (0.0d0,0.0d0)
         m5 = (0.0d0,0.0d0)
         m7 = (0.0d0,0.0d0)
      endif

      m2 = -dotcc(wzw2(0,N),jze(0,-1,2))
      m4 = -dotcc(wzw2(0,N),jze(0,-1,1))
      m6 = -dotcc(wzw1(0,N),jzm(0,-1,2))
      m8 = -dotcc(wzw1(0,N),jzm(0,-1,1))
c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      m2 = m2+zm2i(3)*dotrc(pwz2(0),wzw2(0,N))*dotrc(pwz2(0),jze(0,-1,2))
      m4 = m4+zm2i(3)*dotrc(pwz2(0),wzw2(0,N))*dotrc(pwz2(0),jze(0,-1,1))
      m6 = m6+zm2i(3)*dotrc(pwz1(0),wzw1(0,N))*dotrc(pwz1(0),jzm(0,-1,2))
      m8 = m8+zm2i(3)*dotrc(pwz1(0),wzw1(0,N))*dotrc(pwz1(0),jzm(0,-1,1))

      mat(k,-1,4)=(m1*clr(d(k),1,-1)+m2*clr(d(k),2,-1))*clr(3,3,-1)!WtoWZ2 A/Z
      mat(k,-1,5)=(m3*clr(u(k),1,-1)+m4*clr(u(k),2,-1))*clr(3,3,-1)!A/Z WtoWZ2
      mat(k,-1,6)=(m5*clr(d(k),1,-1)+m6*clr(d(k),2,-1))*clr(3,3,-1)!WtoWZ1 A/Z
      mat(k,-1,7)=(m7*clr(u(k),1,-1)+m8*clr(u(k),2,-1))*clr(3,3,-1)!A/Z WtoWZ1


c Calculation of processes qq -> W -> Z Z W ********************

      call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with wtozzw tensor
      m1 = dotcc(wzzw(0,N),jqq(0,-1))

      mat(k,-1,8) = m1*clr(3,3,-1)


*************************************************************
*
* Calculation of the virtual contributions
*
*************************************************************

c Calculation of processes qq -> Z Z W **********************

#ifdef WITH_NLO
      if (lpent) then
c calculate pentagon contributions
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(k,-1,1)+mat(k,-1,2)+mat(k,-1,3)+mat(k,-1,4)
     1        +mat(k,-1,5)+mat(k,-1,6)+mat(k,-1,7)+mat(k,-1,8)

            do mu = 0,3
               p1(mu) = pbar(mu,1)
               p2(mu) = pbar(mu,2)
               minus_qz1(mu) = -qze(mu,1)
               minus_qz2(mu) = -qzmu(mu,1)
               minus_qw(mu) = -p1(mu)-p2(mu)+qze(mu,1)+qzmu(mu,1)
            enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         if (NLO.eq.1) then ! complete NLO

            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),zmu(1,1),wbos(1,1),
     2         -1,scale,4,1,pentcontr1,j5z1z2(3),test,jbz1z2(3))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),zmu(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(3),test,jba1z2(3))
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),amu(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(3),test,jbz1a2(3))
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),amu(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(3),test,jba1a2(3))
            endif

            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),wbos(1,1),zmu(1,1),
     2         -1,scale,4,1,pentcontr2,j5z1z2(2),test,jbz1z2(2))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),wbos(1,1),zmu(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(2),test,jba1z2(2))
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),wbos(1,1),amu(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(2),test,jbz1a2(2))
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),wbos(1,1),amu(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(2),test,jba1a2(2))
            endif

            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),ze(1,1),zmu(1,1),
     2         -1,scale,4,1,pentcontr3,j5z1z2(1),test,jbz1z2(1))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),ae(1,1),zmu(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(1),test,jba1z2(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),ze(1,1),amu(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(1),test,jbz1a2(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),ae(1,1),amu(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(1),test,jba1a2(1))
            endif

            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),ze(1,1),wbos(1,1),
     2         -1,scale,4,1,pentcontr4,j5z2z1(3),test,jbz2z1(3))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),ze(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(3),test,jba2z1(3))
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),ae(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(3),test,jbz2a1(3))
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),ae(1,1),wbos(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(3),test,jba2a1(3))
            endif

            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),wbos(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr5,j5z2z1(2),test,jbz2z1(2))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),wbos(1,1),ze(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(2),test,jba2z1(2))
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),wbos(1,1),ae(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(2),test,jbz2a1(2))
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),wbos(1,1),ae(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(2),test,jba2a1(2))
            endif

            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),zmu(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr6,j5z2z1(1),test,jbz2z1(1))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),amu(1,1),ze(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(1),test,jba2z1(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),zmu(1,1),ae(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(1),test,jbz2a1(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbos(1,1),amu(1,1),ae(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(1),test,jba2a1(1))
            endif

            call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qz2, minus_qw,
     1               ze(1,1), zmu(1,1), wbos(1,1), scale, als(1,1),
     2               pentcontr1, mborn, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qw, minus_qz2,
     1               ze(1,1), wbos(1,1), zmu(1,1), scale, als(1,1),
     2               pentcontr2, mborn, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qw, minus_qz1, minus_qz2,
     1               wbos(1,1), ze(1,1), zmu(1,1), scale, als(1,1),
     2               pentcontr3, mborn, wardidtest3)
            call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qz1, minus_qw,
     1               zmu(1,1), ze(1,1), wbos(1,1), scale, als(1,1),
     2               pentcontr4, mborn, wardidtest4)
            call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qw, minus_qz1,
     1               zmu(1,1), wbos(1,1), ze(1,1), scale, als(1,1), 
     2               pentcontr5, mborn, wardidtest5)
            call Wardtest_Pent(psi, p1, p2, minus_qw, minus_qz2, minus_qz1,
     1               wbos(1,1), zmu(1,1), ze(1,1), scale, als(1,1),
     2               pentcontr6, mborn, wardidtest6)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1          .or. (wardidtest3.eq.1) .or. (wardidtest4.eq.1) 
     2          .or. (wardidtest5.eq.1) .or. (wardidtest6.eq.1) ) then
               do i = 1,3
                  j5z1z2(i) = (0.0d0, 0.0d0)
                  j5z2z1(i) = (0.0d0, 0.0d0)
                  if(NW.eq.0) then
                     j5a1z2(i) = (0.0d0, 0.0d0)
                     j5a2z1(i) = (0.0d0, 0.0d0)
                     j5z1a2(i) = (0.0d0, 0.0d0)
                     j5z2a1(i) = (0.0d0, 0.0d0)
                     j5a1a2(i) = (0.0d0, 0.0d0)
                     j5a2a1(i) = (0.0d0, 0.0d0)
                  endif
                enddo
c               pentcount = pentcount+1
            else
             do i = 1,3
                j5z1z2(i)=j5z1z2(i)+cte*jbz1z2(i)
                j5z2z1(i)=j5z2z1(i)+cte*jbz2z1(i)
                if(NW.eq.0) then
                j5a1z2(i)=j5a1z2(i)+cte*jba1z2(i)
                j5a2z1(i)=j5a2z1(i)+cte*jba2z1(i)
                j5z1a2(i)=j5z1a2(i)+cte*jbz1a2(i)
                j5z2a1(i)=j5z2a1(i)+cte*jbz2a1(i)
                j5a1a2(i)=j5a1a2(i)+cte*jba1a2(i)
                j5a2a1(i)=j5a2a1(i)+cte*jba2a1(i)
                endif
             enddo
c               pentcount2 = pentcount2+1
            endif

         elseif (NLO.eq.-5) then ! only true pentagon
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),zmup(1,1),wbosp(1,1),
     2         -1,scale,4,1,pentcontr1,j5z1z2(3),test,jbz1z2(3))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),zmup(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(3),test,jba1z2(3))
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),amup(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(3),test,jbz1a2(3))
            call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),amup(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(3),test,jba1a2(3))
            endif

            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),wbosp(1,1),zmup(1,1),
     2         -1,scale,4,1,pentcontr2,j5z1z2(2),test,jbz1z2(2))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),wbosp(1,1),zmup(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(2),test,jba1z2(2))
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),wbosp(1,1),amup(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(2),test,jbz1a2(2))
            call penlinemmm(p1(0),minus_qz1(0),minus_qw(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),wbosp(1,1),amup(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(2),test,jba1a2(2))
            endif

            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),zepp(1,1),zmup(1,1),
     2         -1,scale,4,1,pentcontr3,j5z1z2(1),test,jbz1z2(1))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),aep(1,1),zmup(1,1),
     2         -1,scale,9,-1,test2,j5a1z2(1),test,jba1z2(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),zepp(1,1),amup(1,1),
     2         -1,scale,9,-1,test2,j5z1a2(1),test,jbz1a2(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),aep(1,1),amup(1,1),
     2         -1,scale,9,-1,test2,j5a1a2(1),test,jba1a2(1))
            endif

            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmup(1,1),zepp(1,1),wbosp(1,1),
     2         -1,scale,4,1,pentcontr4,j5z2z1(3),test,jbz2z1(3))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amup(1,1),zepp(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(3),test,jba2z1(3))
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmup(1,1),aep(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(3),test,jbz2a1(3))
            call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amup(1,1),aep(1,1),wbosp(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(3),test,jba2a1(3))
            endif

            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmup(1,1),wbosp(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr5,j5z2z1(2),test,jbz2z1(2))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amup(1,1),wbosp(1,1),zepp(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(2),test,jba2z1(2))
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmup(1,1),wbosp(1,1),aep(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(2),test,jbz2a1(2))
            call penlinemmm(p1(0),minus_qz2(0),minus_qw(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amup(1,1),wbosp(1,1),aep(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(2),test,jba2a1(2))
            endif

            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),zmup(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr6,j5z2z1(1),test,jbz2z1(1))
            if (NW.eq.0) then
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),amup(1,1),zepp(1,1),
     2         -1,scale,9,-1,test2,j5a2z1(1),test,jba2z1(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),zmup(1,1),aep(1,1),
     2         -1,scale,9,-1,test2,j5z2a1(1),test,jbz2a1(1))
            call penlinemmm(p1(0),minus_qw(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),wbosp(1,1),amup(1,1),aep(1,1),
     2         -1,scale,9,-1,test2,j5a2a1(1),test,jba2a1(1))
            endif

            call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qz2, minus_qw,
     1               zepp(1,1), zmup(1,1), wbosp(1,1), scale, als(1,1),
     2               pentcontr1, mborn, wardidtest1)
            call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qw, minus_qz2,
     1               zepp(1,1), wbosp(1,1), zmup(1,1), scale, als(1,1),
     2               pentcontr2, mborn, wardidtest2)
            call Wardtest_Pent(psi, p1, p2, minus_qw, minus_qz1, minus_qz2,
     1               wbosp(1,1), zepp(1,1), zmup(1,1), scale, als(1,1),
     2               pentcontr3, mborn, wardidtest3)
            call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qz1, minus_qw,
     1               zmup(1,1), zepp(1,1), wbosp(1,1), scale, als(1,1),
     2               pentcontr4, mborn, wardidtest4)
            call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qw, minus_qz1,
     1               zmup(1,1), wbosp(1,1), zepp(1,1), scale, als(1,1),
     2               pentcontr5, mborn, wardidtest5)
            call Wardtest_Pent(psi, p1, p2, minus_qw, minus_qz2, minus_qz1,
     1               wbosp(1,1), zmup(1,1), zepp(1,1), scale, als(1,1),
     2               pentcontr6, mborn, wardidtest6)

            if ( (wardidtest1.eq.1) .or. (wardidtest2.eq.1) 
     1          .or. (wardidtest3.eq.1) .or. (wardidtest4.eq.1) 
     2          .or. (wardidtest5.eq.1) .or. (wardidtest6.eq.1) ) then
               do i = 1,3
                  j5z1z2(i) = (0.0d0, 0.0d0)
                  j5z2z1(i) = (0.0d0, 0.0d0)
                  if(NW.eq.0) then
                     j5a1z2(i) = (0.0d0, 0.0d0)
                     j5a2z1(i) = (0.0d0, 0.0d0)
                     j5z1a2(i) = (0.0d0, 0.0d0)
                     j5z2a1(i) = (0.0d0, 0.0d0)
                     j5a1a2(i) = (0.0d0, 0.0d0)
                     j5a2a1(i) = (0.0d0, 0.0d0)
                  endif
                enddo
c               pentcount = pentcount+1
            else
             do i = 1,3
                j5z1z2(i)=j5z1z2(i)+cte*jbz1z2(i)
                j5z2z1(i)=j5z2z1(i)+cte*jbz2z1(i)
                if(NW.eq.0) then
                j5a1z2(i)=j5a1z2(i)+cte*jba1z2(i)
                j5a2z1(i)=j5a2z1(i)+cte*jba2z1(i)
                j5z1a2(i)=j5z1a2(i)+cte*jbz1a2(i)
                j5z2a1(i)=j5z2a1(i)+cte*jbz2a1(i)
                j5a1a2(i)=j5a1a2(i)+cte*jba1a2(i)
                j5a2a1(i)=j5a2a1(i)+cte*jba2a1(i)
                endif
             enddo
c               pentcount2 = pentcount2+1
            endif

c      do i = 10,10000000,1000
c      if (pentcount2.eq.1d2*i) print*,"ratio, pent= ",pentcount*1.0d0/
c     1                        (pentcount+pentcount2)
c      enddo

         else ! pentagon parts ~ boxes

            call pent_to_box(p1, minus_qz1, minus_qz2, minus_qw, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), zmu(1,1),wbos(1,1),
     2         xze, xzmu, xw, -1, scale, j5z1z2(3), jbz1z2(3))
            call pent_to_box(p1, minus_qz1, minus_qw, minus_qz2, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), wbos(1,1),zmu(1,1),
     2         xze, xw, xzmu, -1, scale, j5z1z2(2), jbz1z2(2))
            call pent_to_box(p1, minus_qw, minus_qz1, minus_qz2, p2,
     1         psi(1,-1,2), psi(1,-1,1), wbos(1,1), ze(1,1),zmu(1,1),
     2         xw, xze, xzmu, -1, scale, j5z1z2(1), jbz1z2(1))

            call pent_to_box(p1, minus_qz2, minus_qz1, minus_qw, p2,
     1         psi(1,-1,2), psi(1,-1,1), zmu(1,1), ze(1,1),wbos(1,1),
     2         xzmu, xze, xw, -1, scale, j5z2z1(3), jbz2z1(3))
            call pent_to_box(p1, minus_qz2, minus_qw, minus_qz1, p2,
     1         psi(1,-1,2), psi(1,-1,1), zmu(1,1), wbos(1,1),ze(1,1),
     2         xzmu, xw, xze, -1, scale, j5z2z1(2), jbz2z1(2))
            call pent_to_box(p1, minus_qw, minus_qz2, minus_qz1, p2,
     1         psi(1,-1,2), psi(1,-1,1), wbos(1,1), zmu(1,1),ze(1,1),
     2         xw, xzmu, xze, -1, scale, j5z2z1(1), jbz2z1(1))

            if (NW.eq.0) then
               call pent_to_box(p1, minus_qz1, minus_qz2, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), ae(1,1), zmu(1,1),wbos(1,1),
     2            xae, xzmu, xw, -1, scale, j5a1z2(3), jba1z2(3))
               call pent_to_box(p1, minus_qz1, minus_qw, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), ae(1,1), wbos(1,1),zmu(1,1),
     2            xae, xw, xzmu, -1, scale, j5a1z2(2), jba1z2(2))
               call pent_to_box(p1, minus_qw, minus_qz1, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), ae(1,1),zmu(1,1),
     2            xw, xae, xzmu, -1, scale, j5a1z2(1), jba1z2(1))

               call pent_to_box(p1, minus_qz2, minus_qz1, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), zmu(1,1), ae(1,1),wbos(1,1),
     2            xzmu, xae, xw, -1, scale, j5z2a1(3), jbz2a1(3))
               call pent_to_box(p1, minus_qz2, minus_qw, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), zmu(1,1), wbos(1,1),ae(1,1),
     2            xzmu, xw, xae, -1, scale, j5z2a1(2), jbz2a1(2))
               call pent_to_box(p1, minus_qw, minus_qz2, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), zmu(1,1),ae(1,1),
     2            xw, xzmu, xae, -1, scale, j5z2a1(1), jbz2a1(1))

               call pent_to_box(p1, minus_qz1, minus_qz2, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), ze(1,1), amu(1,1),wbos(1,1),
     2            xze, xamu, xw, -1, scale, j5z1a2(3), jbz1a2(3))
               call pent_to_box(p1, minus_qz1, minus_qw, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), ze(1,1), wbos(1,1),amu(1,1),
     2            xze, xw, xamu, -1, scale, j5z1a2(2), jbz1a2(2))
               call pent_to_box(p1, minus_qw, minus_qz1, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), ze(1,1),amu(1,1),
     2            xw, xze, xamu, -1, scale, j5z1a2(1), jbz1a2(1))

               call pent_to_box(p1, minus_qz2, minus_qz1, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), amu(1,1), ze(1,1),wbos(1,1),
     2            xamu, xze, xw, -1, scale, j5a2z1(3), jba2z1(3))
               call pent_to_box(p1, minus_qz2, minus_qw, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), amu(1,1), wbos(1,1),ze(1,1),
     2            xamu, xw, xze, -1, scale, j5a2z1(2), jba2z1(2))
               call pent_to_box(p1, minus_qw, minus_qz2, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), amu(1,1),ze(1,1),
     2            xw, xamu, xze, -1, scale, j5a2z1(1), jba2z1(1))

               call pent_to_box(p1, minus_qz1, minus_qz2, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), ae(1,1), amu(1,1),wbos(1,1),
     2            xae, xamu, xw, -1, scale, j5a1a2(3), jba1a2(3))
               call pent_to_box(p1, minus_qz1, minus_qw, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), ae(1,1), wbos(1,1),amu(1,1),
     2            xae, xw, xamu, -1, scale, j5a1a2(2), jba1a2(2))
               call pent_to_box(p1, minus_qw, minus_qz1, minus_qz2, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), ae(1,1),amu(1,1),
     2            xw, xae, xamu, -1, scale, j5a1a2(1), jba1a2(1))

               call pent_to_box(p1, minus_qz2, minus_qz1, minus_qw, p2,
     1            psi(1,-1,2), psi(1,-1,1), amu(1,1), ae(1,1),wbos(1,1),
     2            xamu, xae, xw, -1, scale, j5a2a1(3), jba2a1(3))
               call pent_to_box(p1, minus_qz2, minus_qw, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), amu(1,1), wbos(1,1),ae(1,1),
     2            xamu, xw, xae, -1, scale, j5a2a1(2), jba2a1(2))
               call pent_to_box(p1, minus_qw, minus_qz2, minus_qz1, p2,
     1            psi(1,-1,2), psi(1,-1,1), wbos(1,1), amu(1,1),ae(1,1),
     2            xw, xamu, xae, -1, scale, j5a2a1(1), jba2a1(1))
            endif

            j5z1z2(1) = j5z1z2(1)+cte*jbz1z2(1)
            j5z1z2(2) = j5z1z2(2)+cte*jbz1z2(2)
            j5z1z2(3) = j5z1z2(3)+cte*jbz1z2(3)
            j5z2z1(1) = j5z2z1(1)+cte*jbz2z1(1)
            j5z2z1(2) = j5z2z1(2)+cte*jbz2z1(2)
            j5z2z1(3) = j5z2z1(3)+cte*jbz2z1(3)
       
            if (NW.eq.0) then     
               j5a1z2(1) = j5a1z2(1)+cte*jba1z2(1)
               j5a1z2(2) = j5a1z2(2)+cte*jba1z2(2)
               j5a1z2(3) = j5a1z2(3)+cte*jba1z2(3)
               j5z2a1(1) = j5z2a1(1)+cte*jbz2a1(1)
               j5z2a1(2) = j5z2a1(2)+cte*jbz2a1(2)
               j5z2a1(3) = j5z2a1(3)+cte*jbz2a1(3)
               j5z1a2(1) = j5z1a2(1)+cte*jbz1a2(1)
               j5z1a2(2) = j5z1a2(2)+cte*jbz1a2(2)
               j5z1a2(3) = j5z1a2(3)+cte*jbz1a2(3)
               j5a2z1(1) = j5a2z1(1)+cte*jba2z1(1)
               j5a2z1(2) = j5a2z1(2)+cte*jba2z1(2)
               j5a2z1(3) = j5a2z1(3)+cte*jba2z1(3)
               j5a1a2(1) = j5a1a2(1)+cte*jba1a2(1)
               j5a1a2(2) = j5a1a2(2)+cte*jba1a2(2)
               j5a1a2(3) = j5a1a2(3)+cte*jba1a2(3)
               j5a2a1(1) = j5a2a1(1)+cte*jba2a1(1)
               j5a2a1(2) = j5a2a1(2)+cte*jba2a1(2)
               j5a2a1(3) = j5a2a1(3)+cte*jba2a1(3)
            endif

         endif !NLO

c   q qbar -> Z Z W
         mv1 = j5z1z2(3)*clr(u(k),2,-1)**2
     2       + j5z1z2(2)*clr(u(k),2,-1)*clr(d(k),2,-1)
     3       + j5z1z2(1)*clr(d(k),2,-1)**2
     4       + j5z2z1(3)*clr(u(k),2,-1)**2
     5       + j5z2z1(2)*clr(u(k),2,-1)*clr(d(k),2,-1)
     6       + j5z2z1(1)*clr(d(k),2,-1)**2

         if (NW.eq.0) then
            mv2 = j5a1z2(3)*clr(u(k),1,-1)*clr(u(k),2,-1)
     2          + j5a1z2(2)*clr(u(k),1,-1)*clr(d(k),2,-1)
     3          + j5a1z2(1)*clr(d(k),1,-1)*clr(d(k),2,-1)
     4          + j5z2a1(3)*clr(u(k),1,-1)*clr(u(k),2,-1)
     5          + j5z2a1(2)*clr(u(k),2,-1)*clr(d(k),1,-1)
     6          + j5z2a1(1)*clr(d(k),1,-1)*clr(d(k),2,-1)

            mv3 = j5z1a2(3)*clr(u(k),1,-1)*clr(u(k),2,-1)
     2          + j5z1a2(2)*clr(u(k),2,-1)*clr(d(k),1,-1)
     3          + j5z1a2(1)*clr(d(k),1,-1)*clr(d(k),2,-1)
     4          + j5a2z1(3)*clr(u(k),1,-1)*clr(u(k),2,-1)
     5          + j5a2z1(2)*clr(u(k),1,-1)*clr(d(k),2,-1)
     6          + j5a2z1(1)*clr(d(k),1,-1)*clr(d(k),2,-1)

            mv4 = j5a1a2(3)*clr(u(k),1,-1)**2
     2          + j5a1a2(2)*clr(u(k),1,-1)*clr(d(k),1,-1)
     3          + j5a1a2(1)*clr(d(k),1,-1)**2
     4          + j5a2a1(3)*clr(u(k),1,-1)**2
     5          + j5a2a1(2)*clr(u(k),1,-1)*clr(d(k),1,-1)
     6          + j5a2a1(1)*clr(d(k),1,-1)**2
         else if (NW.eq.1) then
            mv2 = (0.0d0, 0.0d0)
            mv3 = (0.0d0, 0.0d0)
            mv4 = (0.0d0, 0.0d0)
         endif


         matv(k,-1,1) = (mv1 + mv2 + mv3 + mv4)*clr(3,3,-1)

      endif  !lpent
#endif

c Calculation of processes qq -> V V -> 2l + 4l ****************
#ifdef WITH_NLO
      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge check:
         mborn=mat(k,-1,1)+mat(k,-1,2)+mat(k,-1,3)+mat(k,-1,4)
     1        +mat(k,-1,5)+mat(k,-1,6)+mat(k,-1,7)+mat(k,-1,8)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qzz(mu) = -qze(mu,1)-qzmu(mu,1)
            minus_qz1(mu) = -qze(mu,1)
            minus_qz2(mu) = -qzmu(mu,1)
            minus_qw(mu) =  -p1(mu)-p2(mu)+qze(mu,1)+qzmu(mu,1)
            minus_qz1w(mu) = -p1(mu)-p2(mu)+qzmu(mu,1)
            minus_qz2w(mu) = -p1(mu)-p2(mu)+qze(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call boxlinemm(p1,minus_qw,minus_qzz,p2,
     1      psi(1,-1,2),psi(1,-1,1),wbos(1,1),azezm(0,1),
     2      -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3      result_born(-1,1))

         call boxlinemm(p1,minus_qw,minus_qzz,p2,
     1      psi(1,-1,2),psi(1,-1,1),wbos(1,1),zzezm(0,1),
     2      -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3      result_born(-1,2))

         call boxlinemm(p1,minus_qzz,minus_qw,p2,
     1      psi(1,-1,2),psi(1,-1,1),azezm(0,1),wbos(1,1),
     2      -1,scale,3,1,gaugebox2,result(-1,3),gaugebox2_born,
     3      result_born(-1,3))

         call boxlinemm(p1,minus_qzz,minus_qw,p2,
     1      psi(1,-1,2),psi(1,-1,1),zzezm(0,1),wbos(1,1),
     2      -1,scale,9,-1,boxtest,result(-1,4),treetest,
     3      result_born(-1,4))

         call Wardtest_Box(wbos(1,1), qwbos(0,1), azezm(0,1), qzz(0,1),
     1            gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         call Wardtest_Box(azezm(0,1), qzz(0,1), wbos(1,1), qwbos(0,1),
     1            gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         if (wardidtest1.eq.1) then
            matv(k,-1,2) = (0.0d0, 0.0d0)
         else
            mv1=-(result(-1,1)+cte*result_born(-1,1))
            mv2=-(result(-1,2)+cte*result_born(-1,2))
            matv(k,-1,2) = ( mv1*clr(d(k),1,-1) + mv2*clr(d(k),2,-1) )
     1                               *clr(3,3,-1)!qq -> W A/ZtoZZ
         endif

         if (wardidtest2.eq.1) then
            matv(k,-1,3) = (0.0d0, 0.0d0)
         else
            mv3=-(result(-1,3)+cte*result_born(-1,3))
            mv4=-(result(-1,4)+cte*result_born(-1,4))
            matv(k,-1,3) = ( mv3*clr(u(k),1,-1) + mv4*clr(u(k),2,-1) )
     1                               *clr(3,3,-1)!qq -> A/ZtoZZ W
         endif

         call boxlinemm(p1,minus_qz2w,minus_qz1,p2,
     1      psi(1,-1,2),psi(1,-1,1),wzw2(0,1),ze(1,1),
     2      -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3      result_born(-1,1))
         if (NW.eq.0) then
            call boxlinemm(p1,minus_qz2w,minus_qz1,p2,
     1         psi(1,-1,2),psi(1,-1,1),wzw2(0,1),ae(1,1),
     2         -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3         result_born(-1,2))
         endif

         call boxlinemm(p1,minus_qz1,minus_qz2w,p2,
     1      psi(1,-1,2),psi(1,-1,1),ze(1,1),wzw2(0,1),
     2      -1,scale,3,1,gaugebox2,result(-1,3),gaugebox2_born,
     3      result_born(-1,3))
         if (NW.eq.0) then
            call boxlinemm(p1,minus_qz1,minus_qz2w,p2,
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),wzw2(0,1),
     2         -1,scale,9,-1,boxtest,result(-1,4),treetest,
     3         result_born(-1,4))
         endif

         call boxlinemm(p1,minus_qz1w,minus_qz2,p2,
     1      psi(1,-1,2),psi(1,-1,1),wzw1(0,1),zmu(1,1),
     2      -1,scale,3,1,gaugebox3,result(-1,5),gaugebox3_born,
     3      result_born(-1,5))
         if (NW.eq.0) then
            call boxlinemm(p1,minus_qz1w,minus_qz2,p2,
     1         psi(1,-1,2),psi(1,-1,1),wzw1(0,1),amu(1,1),
     2         -1,scale,9,-1,boxtest,result(-1,6),treetest,
     3         result_born(-1,6))
         endif

         call boxlinemm(p1,minus_qz2,minus_qz1w,p2,
     1      psi(1,-1,2),psi(1,-1,1),zmu(1,1),wzw1(0,1),
     2      -1,scale,3,1,gaugebox4,result(-1,7),gaugebox4_born,
     3      result_born(-1,7))
         if (NW.eq.0) then
            call boxlinemm(p1,minus_qz2,minus_qz1w,p2,
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),wzw1(0,1),
     2         -1,scale,9,-1,boxtest,result(-1,8),treetest,
     3         result_born(-1,8))
         endif

         call Wardtest_Box(wzw2(0,1), pwz2(0), ze(1,1), qze(0,1),
     1            gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         call Wardtest_Box(ze(1,1), qze(0,1), wzw2(0,1), pwz2(0),
     1            gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         call Wardtest_Box(wzw1(0,1), pwz1(0), zmu(1,1), qzmu(0,1),
     1            gaugebox3, gaugebox3_born, mborn, wardidtest3, als(1,1))

         call Wardtest_Box(zmu(1,1), qzmu(0,1), wzw1(0,1), pwz1(0),
     1            gaugebox4, gaugebox4_born, mborn, wardidtest4, als(1,1))

         if (wardidtest1.eq.1) then
            matv(k,-1,4) = (0.0d0, 0.0d0)
         else
            if (NW.eq.0) then
               mv1 = -(result(-1,2)+cte*result_born(-1,2))
            else
               mv1 = (0.0d0,0.0d0)
            endif
            mv2=-(result(-1,1)+cte*result_born(-1,1))

            matv(k,-1,4)=(mv1*clr(d(k),1,-1)+mv2*clr(d(k),2,-1))
     &                                                *clr(3,3,-1)
         endif

         if (wardidtest2.eq.1) then
            matv(k,-1,5) = (0.0d0, 0.0d0)
         else
            if (NW.eq.0) then
               mv3 = -(result(-1,4)+cte*result_born(-1,4))
            else
               mv3 = (0.0d0,0.0d0)
            endif
            mv4=-(result(-1,3)+cte*result_born(-1,3))

            matv(k,-1,5)=(mv3*clr(u(k),1,-1)+mv4*clr(u(k),2,-1))
     &                                                *clr(3,3,-1)
         endif

         if (wardidtest3.eq.1) then
            matv(k,-1,6) = (0.0d0, 0.0d0)
         else
            if (NW.eq.0) then
               mv5 = -(result(-1,6)+cte*result_born(-1,6))
            else
               mv5 = (0.0d0,0.0d0)
            endif
            mv6=-(result(-1,5)+cte*result_born(-1,5))

            matv(k,-1,6)=(mv5*clr(d(k),1,-1)+mv6*clr(d(k),2,-1))
     &                                                *clr(3,3,-1)
         endif

         if (wardidtest4.eq.1) then
            matv(k,-1,7) = (0.0d0, 0.0d0)
         else
            if (NW.eq.0) then
               mv7 = -(result(-1,8)+cte*result_born(-1,8))
            else
               mv7 = (0.0d0,0.0d0)
            endif
            mv8=-(result(-1,7)+cte*result_born(-1,7))

            matv(k,-1,7)=(mv7*clr(u(k),1,-1)+mv8*clr(u(k),2,-1))
     &                                                *clr(3,3,-1)
         endif

      endif !lbox
#endif

c Calculation of processes qq -> W -> Z Z W ********************

      matv(k,-1,8) = (0.0d0, 0.0d0)

c**************************************************************
c
c comparison with madgraph ************************************
c
c**************************************************************
#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,8
               madgr(i,isig,k) = 0d0
            enddo
         enddo

         do i =1,8
            madgr(i,-1,k) = tamp(i,-1,k)
         enddo

c print results from checks:
         do i=1,8
    	    print*,"ratio for k = ",k," isig = ",-1," i = ",i,":",
     &	    	    mat(k,-1,i)/madgr(i,-1,k)
            print*,"mat(",k,",-1,",i,") = ",mat(k,-1,i)
            temp = dreal(mat(k,-1,i))**2+dimag(mat(k,-1,i))**2
            print*,temp
            print*,"madgr(",i,",-1,",k,") = ",madgr(i,-1,k)
            temp1 = dreal(madgr(i,-1,k))**2+dimag(madgr(i,-1,k))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "
         enddo

      endif ! end of debugging block
#endif

c**************************************************************
c
c sum the graphs, square them and map them onto qqwzz *********
c
c**************************************************************
      res(k) = 0
      resv(k) = 0

      mm(k) = 0
      do i = 1,8
         mm(k) = mm(k) +  mat(k,-1,i)
      enddo
      res(k)= res(k) + dreal(mm(k))**2 + dimag(mm(k))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqZZWj.f)

#ifdef WITH_NLO
      if (lnlo) then
         mv(k)=matv(k,-1,1)+matv(k,-1,2)+matv(k,-1,3)+matv(k,-1,4)
     1        +matv(k,-1,5)+matv(k,-1,6)+matv(k,-1,7)+matv(k,-1,8)

      endif
#endif

      if (nlo.gt.0) then
#ifdef WITH_NLO
         mv(k) = als(1,1)*c2o4pi*(mv(k) + mm(k)* cvirtVV)
#endif 
      else if(nlo.lt.0) then
#ifdef WITH_NLO
         mv(k) = als(1,1)*c2o4pi*mv(k)
#endif
      else
         mv(k) = (0d0, 0d0)
      endif

      resv(k) = resv(k) + 2d0*dreal( mm(k)*conjg( mv(k) ) )

      if (nlo.ge.0) then
         res(k) = (res(k)+resv(k))*3d0      ! 3 is the color sum factor
      else
         res(k) = resv(k)*3d0               ! 3 is the color sum factor
      endif  

      qqWZZ = res(k)

      return
      end


c****************************************************************************
