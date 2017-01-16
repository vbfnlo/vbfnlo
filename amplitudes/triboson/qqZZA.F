c****************subroutine QQZZA *************************************
c
c	Giuseppe Bozzi, <giuseppe@particle.uni-karlsruhe.de>
c	Modified by Francisco Campanario <francam@particle.uni-karlsruhe.de>
c	Last modified: Jan. 2009 by V. Hankele <vera@particle.uni-karlsruhe.de>
C
C  QQZZA calculates the matrix elements**2 for 
C        q1 q1bar    ---->   Z Z A,   
C         Z ---> f1-bar f1, Z ---> f2-bar f2
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
c       no reduction of pentagon contribution as in the other tri-boson
c       codes is performed because it drops out in the end due to gauge
c       invariance.
c
c*********************************************************************

      subroutine qqZZA(pbar,fsign, nlo, N, uuzza, ddzza)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,3+n_v)
      double precision  uuzza, ddzza, res(2),resv(2)
      double precision  p(0:3,3+n_v)
      double precision  fqa(0:4,2), fqze(0:4,2), fqzmu(0:4,2)
      integer  fsign(2+n_v), sign(3+n_v), nlo, mu, i, k, isig, N
      integer  ifl(2)
      double complex psize(2,-1:1,2), psiae(2,-1:1,2), psia(2,-1:1,2),
     1               psizmu(2,-1:1,2), psiamu(2,-1:1,2)
      double complex jae(0:5,-1:1,2), jam(0:5,-1:1,2), ja(0:5,-1:1,2),
     1               jze(0:5,-1:1,2),jzm(0:5,-1:1,2)


      double complex temp, temp1
      double complex mat(2,-1:1,7), matv(2,-1:1,7)
      double complex madgr(7,-1:1,2)
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex ma, mz, m1, m2, m3, m4
      double complex mv1(3), mv2(3), mv3(3), mv4(3)
      double complex psi(2,-1:1,2), jqq(0:5,-1:1)
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (zero = (0d0,0d0) )

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      double precision scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c variables for new implementation of the pentagon routine
      double complex j55z1z2(-1:1,3), j55a1z2(-1:1,3),j55z1a2(-1:1,3)
     - ,j55a1a2(-1:1,3)

      double complex j55z2z1(-1:1,1:3), j55a2z1(-1:1,1:3),j55z2a1(-1:1,1:3)
     -  ,j55a2a1(-1:1,1:3)
      double complex jbz1z2(-1:1,3), jba1z2(-1:1,3),jbz1a2(-1:1,3)
     - ,jba1a2(-1:1,3)
      double complex jbz2z1(-1:1,1:3), jba2z1(-1:1,1:3),jbz2a1(-1:1,1:3)
     - ,jba2a1(-1:1,1:3)

     - 

     -  

     - 

     - 
      double complex borntest(3), penttest(3)
      double precision minus_qz1(0:3),minus_qz2(0:3),minus_qa(0:3)

      double precision p1(0:3), p2(0:3) 
      double precision qz1a(0:3), qz2a(0:3)

      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)


      double complex mborn, test, cte
      integer wardidtest1,wardidtest2,wardidtest3
      integer wardidtest4,wardidtest5,wardidtest6     


c variables for new implementation of boxes
      double complex gaugebox1(2),gaugebox1_born(2)
      double complex gaugebox2(2),gaugebox2_born(2)
      double complex result(-1:1,8),result_born(-1:1,8)
      double precision  minus_qz1z2(0:3),minus_qaz1(0:3),minus_qaz2(0:3)

      double complex boxtest(2), treetest(2)

c variables for test purposes


c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0           ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5) ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5) ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4) ! q^mu/m_V terms of Pentagon
      lpq=.false.
      lpent = lpt .or. lpq

c     define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzza: u = 3
         ifl(k) = 3
         k = 2                  ! ddzza: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,7
                  mat(k,isig,i) = (0.0d0,0.0d0)
                  matv(k,isig,i) = (0.0d0,0.0d0)
            enddo
         enddo
      enddo

c define the internal momenta
      do i = 1,2+n_v
         sign(i) = fsign(i)
      enddo
      sign(8) = -1 
      do mu = 0,3
         do i = 1,2+n_v+1
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

      do mu = 0,3
         qz1a(mu) = qz(mu,1) + qze(mu,1)
         qz2a(mu) = qz(mu,1) + qzmu(mu,1)
      enddo

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,N),
     &              ze(1,N),psize(1,isig,2),fqze(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qzmu(0,N),
     &              zmu(1,N),psizmu(1,isig,2),fqzmu(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,N),
     &              ae(1,N),psiae(1,isig,2),fqze(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qzmu(0,N),
     &              amu(1,N),psiamu(1,isig,2),fqzmu(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,N),
     &              atau(1,N), psia(1,isig,2),fqa(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,N),
     &              ze(1,N), psize(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qzmu(0,N),
     &              zmu(1,N), psizmu(1,isig,1),fqzmu(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,N),
     &              ae(1,N), psiae(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qzmu(0,N),
     &              amu(1,N), psiamu(1,isig,1),fqzmu(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,N),
     &              atau(1,N), psia(1,isig,1),fqa(0,1))

      enddo

******************************************************
*
C Calculation of born matrix elements
*
**************************************************************

c Calculation of processes qq -> Z Z Z ---------------------
     
c q q -> Z Z A
      do isig = -1,1,2
         m1=-s1c(psize(1,isig,2),atau(1,N),.true.,isig,psizmu(1,isig,1))
     1      -s1c(psize(1,isig,2),zmu(1,N),.true.,isig,psia(1,isig,1))
     2      -s1c(psizmu(1,isig,2),atau(1,N),.true.,isig,psize(1,isig,1))
     3      -s1c(psizmu(1,isig,2),ze(1,N),.true.,isig,psia(1,isig,1))
     4      -s1c(psia(1,isig,2),zmu(1,N),.true.,isig,psize(1,isig,1))
     5      -s1c(psia(1,isig,2),ze(1,N),.true.,isig,psizmu(1,isig,1))

         m2=-s1c(psiae(1,isig,2),atau(1,N),.true.,isig,psiamu(1,isig,1))
     1      -s1c(psiae(1,isig,2),amu(1,N),.true.,isig,psia(1,isig,1))
     2      -s1c(psiamu(1,isig,2),atau(1,N),.true.,isig,psiae(1,isig,1))
     3      -s1c(psiamu(1,isig,2),ae(1,N),.true.,isig,psia(1,isig,1))
     4      -s1c(psia(1,isig,2),amu(1,N),.true.,isig,psiae(1,isig,1))
     5      -s1c(psia(1,isig,2),ae(1,N),.true.,isig,psiamu(1,isig,1))

         m3=-s1c(psiae(1,isig,2),atau(1,N),.true.,isig,psizmu(1,isig,1))
     1      -s1c(psiae(1,isig,2),zmu(1,N),.true.,isig,psia(1,isig,1))
     2      -s1c(psiamu(1,isig,2),atau(1,N),.true.,isig,psize(1,isig,1))
     3      -s1c(psiamu(1,isig,2),ze(1,N),.true.,isig,psia(1,isig,1))
     4      -s1c(psia(1,isig,2),amu(1,N),.true.,isig,psize(1,isig,1))
     5      -s1c(psia(1,isig,2),ae(1,N),.true.,isig,psizmu(1,isig,1))
     6      -s1c(psia(1,isig,2),zmu(1,N),.true.,isig,psiae(1,isig,1))
     7      -s1c(psia(1,isig,2),ze(1,N),.true.,isig,psiamu(1,isig,1))
     8      -s1c(psize(1,isig,2),atau(1,N),.true.,isig,psiamu(1,isig,1))
     9      -s1c(psize(1,isig,2),amu(1,N),.true.,isig,psia(1,isig,1))
     &      -s1c(psizmu(1,isig,2),atau(1,N),.true.,isig,psiae(1,isig,1))
     &      -s1c(psizmu(1,isig,2),ae(1,N),.true.,isig,psia(1,isig,1))

         mat(1,isig,1) =  m1 * clr(3,1,isig) * clr(3,2,isig)**2 ! uu -> ZZA
     1                  + m2 * clr(3,1,isig)**3                 ! uu -> AAA
     3                  + m3 * clr(3,2,isig) * clr(3,1,isig)**2 ! uu -> ZAA

         mat(2,isig,1) =  m1 * clr(4,1,isig) * clr(4,2,isig)**2 ! dd -> ZZA
     1                  + m2 * clr(4,1,isig)**3                 ! dd -> AAA
     3                  + m3 * clr(4,2,isig) * clr(4,1,isig)**2 ! dd -> ZAA

      enddo
  
c**********************************************************************

c Calculation of processes qq -> V V -> Z Z A -------------

c qq -> A A/Z (-> Z Z)
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psia(1,-1,1),fqa(0,1), ja(0,-1,1) )

c contract with atozz/ztozz tensor
      do isig = -1,1,2
         m1 = -dotcc(azz(0,N),ja(0,isig,1))
         m2 = -dotcc(zzztens(0,N),ja(0,isig,1))
         do k = 1,2
            mat(k,isig,2) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
         enddo
      enddo

c qq -> A/Z(-> Z Z) A
      call curr6(1,psia(1,-1,2),fqa(0,2),
     &           psi(1,-1,1),p(0,1), ja(0,-1,2) )

c contract with atozz/ztozz tensor
      do isig = -1,1,2
         m1 = -dotcc(azz(0,N),ja(0,isig,2))
         m2 = -dotcc(zzztens(0,N),ja(0,isig,2))
         do k = 1,2
            mat(k,isig,3) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
         enddo
      enddo

c qq -> Z A/Z (-> Z A)
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psize(1,-1,1),fqze(0,1), jze(0,-1,1) )
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psizmu(1,-1,1),fqzmu(0,1), jzm(0,-1,1) )
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psiae(1,-1,1),fqze(0,1), jae(0,-1,1) )
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psiamu(1,-1,1),fqzmu(0,1), jam(0,-1,1) )

c contract with atozz/ztozz tensor
      do isig = -1,1,2
         m1 = +dotcc(aazm(0,N),jze(0,isig,1))
     1        +dotcc(aaze(0,N),jzm(0,isig,1))
         m2 = +dotcc(aazm(0,N),jae(0,isig,1))
     1        +dotcc(aaze(0,N),jam(0,isig,1))
         m3 = +dotcc(zzam(0,N),jze(0,isig,1))
     1        +dotcc(zzae(0,N),jzm(0,isig,1))
         m4 = +dotcc(zzam(0,N),jae(0,isig,1))
     1        +dotcc(zzae(0,N),jam(0,isig,1))
         do k = 1,2
            mat(k,isig,4) =  m1*clr(ifl(k),1,isig)*clr(ifl(k),2,isig)
     1                     + m2*clr(ifl(k),1,isig)**2
     1                     + m3*clr(ifl(k),2,isig)**2
     1                     + m4*clr(ifl(k),2,isig)*clr(ifl(k),1,isig)
         enddo
      enddo

c qq -> A/Z(-> Z A) Z
      call curr6(1,psize(1,-1,2),fqze(0,2),
     &           psi(1,-1,1),p(0,1), jze(0,-1,2) )
      call curr6(1,psizmu(1,-1,2),fqzmu(0,2),
     &           psi(1,-1,1),p(0,1), jzm(0,-1,2) )
      call curr6(1,psiae(1,-1,2),fqze(0,2),
     &           psi(1,-1,1),p(0,1), jae(0,-1,2) )
      call curr6(1,psiamu(1,-1,2),fqzmu(0,2),
     &           psi(1,-1,1),p(0,1), jam(0,-1,2) )

c contract with atozz/ztozz tensor
      do isig = -1,1,2
         m1 = +dotcc(aazm(0,N),jze(0,isig,2))
     1        +dotcc(aaze(0,N),jzm(0,isig,2))
         m2 = +dotcc(aazm(0,N),jae(0,isig,2))
     1        +dotcc(aaze(0,N),jam(0,isig,2))
         m3 = +dotcc(zzam(0,N),jze(0,isig,2))
     1        +dotcc(zzae(0,N),jzm(0,isig,2))
         m4 = +dotcc(zzam(0,N),jae(0,isig,2))
     1        +dotcc(zzae(0,N),jam(0,isig,2))
         do k = 1,2
            mat(k,isig,5) =  m1*clr(ifl(k),1,isig)*clr(ifl(k),2,isig)
     1                     + m2*clr(ifl(k),1,isig)**2
     1                     + m3*clr(ifl(k),2,isig)**2
     1                     + m4*clr(ifl(k),2,isig)*clr(ifl(k),1,isig)
         enddo
      enddo

c  Calculation of processes qq -> A/Z -> Z Z Z ---------------------
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atowwz/ztowwz tensor
      do isig = -1,1,2
         ma = -dotcc(azza(0,N),jqq(0,isig))
         mz = -dotcc(zzza(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,6)  = ma*clr(ifl(k),1,isig)
            mat(k,isig,7)  = mz*clr(ifl(k),2,isig)
         enddo
      enddo

******************************************************
*
* Calculation of virtual contribution matrix elements
*
**************************************************************

#ifdef WITH_NLO

c Calculation of processes qq -> Z Z A ---------------------

c**********************************************************************
c calculate pentagon contributions

      if (lpent) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge test
         mborn=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)+mat(1,-1,4)+mat(1,-1,5)
     1     +mat(1,-1,6)+mat(1,-1,7)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qz1(mu) = -qze(mu,1)
            minus_qa(mu)  = -qz(mu,1)
            minus_qz2(mu) = -p1(mu)-p2(mu)+qze(mu,1)+qz(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),zmu(1,1),atau(1,1),
     2         -1,scale,4,1,pentcontr1,j55z1z2(-1,3),borntest,jbz1z2(-1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),zmu(1,1),atau(1,1),
     2         1,scale,4,1,penttest,j55z1z2(1,3),borntest,jbz1z2(1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),amu(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55z1a2(-1,3),test,jbz1a2(-1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),amu(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55z1a2(1,3),borntest,jbz1a2(1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),zmu(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55a1z2(-1,3),borntest,jba1z2(-1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),zmu(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55a1z2(1,3),borntest,jba1z2(1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),amu(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55a1a2(-1,3),borntest,jba1a2(-1,3))
         call penlinemmm(p1(0),minus_qz1(0),minus_qz2(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),amu(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55a1a2(1,3),borntest,jba1a2(1,3))

         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau(1,1),amu(1,1),
     2         -1,scale,4,1,pentcontr2,j55z1a2(-1,2),borntest,jbz1a2(-1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau(1,1),amu(1,1),
     2         1,scale,9,-1,penttest,j55z1a2(1,2),borntest,jbz1a2(1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau(1,1),zmu(1,1),
     2         -1,scale,9,-1,penttest,j55z1z2(-1,2),borntest,jbz1z2(-1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau(1,1),zmu(1,1),
     2         1,scale,9,-1,penttest,j55z1z2(1,2),borntest,jbz1z2(1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),atau(1,1),zmu(1,1),
     2         -1,scale,9,-1,penttest,j55a1z2(-1,2),borntest,jba1z2(-1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),atau(1,1),zmu(1,1),
     2         1,scale,9,-1,penttest,j55a1z2(1,2),borntest,jba1z2(1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),atau(1,1),amu(1,1),
     2         -1,scale,9,-1,penttest,j55a1a2(-1,2),borntest,jba1a2(-1,2))
         call penlinemmm(p1(0),minus_qz1(0),minus_qa(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),atau(1,1),amu(1,1),
     2         1,scale,9,-1,penttest,j55a1a2(1,2),borntest,jba1a2(1,2))

         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),ae(1,1),zmu(1,1),
     2         -1,scale,4,1,pentcontr3,j55a1z2(-1,1),borntest,jba1z2(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),ae(1,1),zmu(1,1),
     2         1,scale,9,-1,penttest,j55a1z2(1,1),borntest,jba1z2(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),ze(1,1),zmu(1,1),
     2         -1,scale,9,-1,penttest,j55z1z2(-1,1),borntest,jbz1z2(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),ze(1,1),zmu(1,1),
     2         1,scale,9,-1,penttest,j55z1z2(1,1),borntest,jbz1z2(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),ze(1,1),amu(1,1),
     2         -1,scale,9,-1,penttest,j55z1a2(-1,1),borntest,jbz1a2(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),ze(1,1),amu(1,1),
     2         1,scale,9,-1,penttest,j55z1a2(1,1),borntest,jbz1a2(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),ae(1,1),amu(1,1),
     2         -1,scale,9,-1,penttest,j55a1a2(-1,1),borntest,jba1a2(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz1(0),minus_qz2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),ae(1,1),amu(1,1),
     2         1,scale,9,-1,penttest,j55a1a2(1,1),borntest,jba1a2(1,1))

         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),ae(1,1),atau(1,1),
     2         -1,scale,4,1,pentcontr4,j55a2a1(-1,3),borntest,jba2a1(-1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),amu(1,1),ae(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55a2a1(1,3),borntest,jba2a1(1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),ze(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55z2z1(-1,3),borntest,jbz2z1(-1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zmu(1,1),ze(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55z2z1(1,3),borntest,jbz2z1(1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),ae(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55z2a1(-1,3),borntest,jbz2a1(-1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zmu(1,1),ae(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55z2a1(1,3),borntest,jbz2a1(1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),ze(1,1),atau(1,1),
     2         -1,scale,9,-1,penttest,j55a2z1(-1,3),borntest,jba2z1(-1,3))
         call penlinemmm(p1(0),minus_qz2(0),minus_qz1(0),minus_qa(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),amu(1,1),ze(1,1),atau(1,1),
     2         1,scale,9,-1,penttest,j55a2z1(1,3),borntest,jba2z1(1,3))

         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),atau(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr5,j55z2z1(-1,2),borntest,jbz2z1(-1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zmu(1,1),atau(1,1),ze(1,1),
     2         1,scale,9,-1,penttest,j55z2z1(1,2),borntest,jbz2z1(1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zmu(1,1),atau(1,1),ae(1,1),
     2         -1,scale,9,-1,penttest,j55z2a1(-1,2),borntest,jbz2a1(-1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zmu(1,1),atau(1,1),ae(1,1),
     2         1,scale,9,-1,penttest,j55z2a1(1,2),borntest,jbz2a1(1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),atau(1,1),ze(1,1),
     2         -1,scale,9,-1,penttest,j55a2z1(-1,2),borntest,jba2z1(-1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),amu(1,1),atau(1,1),ze(1,1),
     2         1,scale,9,-1,penttest,j55a2z1(1,2),borntest,jba2z1(1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),amu(1,1),atau(1,1),ae(1,1),
     2         -1,scale,9,-1,penttest,j55a2a1(-1,2),borntest,jba2a1(-1,2))
         call penlinemmm(p1(0),minus_qz2(0),minus_qa(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),amu(1,1),atau(1,1),ae(1,1),
     2         1,scale,9,-1,penttest,j55a2a1(1,2),borntest,jba2a1(1,2))

         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),zmu(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr6,j55z2z1(-1,1),borntest,jbz2z1(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),zmu(1,1),ze(1,1),
     2         1,scale,9,-1,penttest,j55z2z1(1,1),borntest,jbz2z1(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),zmu(1,1),ae(1,1),
     2         -1,scale,9,-1,penttest,j55z2a1(-1,1),borntest,jbz2a1(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),zmu(1,1),ae(1,1),
     2         1,scale,9,-1,penttest,j55z2a1(1,1),borntest,jbz2a1(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),amu(1,1),ze(1,1),
     2         -1,scale,9,-1,penttest,j55a2z1(-1,1),borntest,jba2z1(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),amu(1,1),ze(1,1),
     2         1,scale,9,-1,penttest,j55a2z1(1,1),borntest,jba2z1(1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau(1,1),amu(1,1),ae(1,1),
     2         -1,scale,9,-1,penttest,j55a2a1(-1,1),borntest,jba2a1(-1,1))
         call penlinemmm(p1(0),minus_qa(0),minus_qz2(0),minus_qz1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau(1,1),amu(1,1),ae(1,1),
     2         1,scale,9,-1,penttest,j55a2a1(1,1),borntest,jba2a1(1,1))

         call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qz2, minus_qa,
     1               ze(1,1), zmu(1,1), atau(1,1), scale, als(1,1),
     2               pentcontr1, mborn, wardidtest1)
         call Wardtest_Pent(psi, p1, p2, minus_qz1, minus_qa, minus_qz2,
     1               ze(1,1), atau(1,1), amu(1,1), scale, als(1,1),
     2               pentcontr2, mborn, wardidtest2)
         call Wardtest_Pent(psi, p1, p2, minus_qa, minus_qz1, minus_qz2,
     1               atau(1,1), ae(1,1), zmu(1,1), scale, als(1,1),
     2               pentcontr3, mborn, wardidtest3)
         call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qz1, minus_qa,
     1               amu(1,1), ae(1,1), atau(1,1), scale, als(1,1),
     2               pentcontr4, mborn, wardidtest4)
         call Wardtest_Pent(psi, p1, p2, minus_qz2, minus_qa, minus_qz1,
     1               zmu(1,1), atau(1,1), ze(1,1), scale, als(1,1),
     2               pentcontr5, mborn, wardidtest5)
         call Wardtest_Pent(psi, p1, p2, minus_qa, minus_qz2, minus_qz1,
     1               atau(1,1), zmu(1,1), ze(1,1), scale, als(1,1),
     2               pentcontr6, mborn, wardidtest6)

         if (wardidtest1.eq.1 .or. wardidtest2.eq.1 .or. 
     1               wardidtest3.eq.1 .or. wardidtest4.eq.1 .or. 
     2               wardidtest5.eq.1 .or. wardidtest6.eq.1) then 
            do isig=-1,1,2  
               do i = 1,3
                  j55z1z2(isig,i) = (0.0d0, 0.0d0)
                  j55z2z1(isig,i) = (0.0d0, 0.0d0)
                  j55a1z2(isig,i) = (0.0d0, 0.0d0)
                  j55a2z1(isig,i) = (0.0d0, 0.0d0)
                  j55z1a2(isig,i) = (0.0d0, 0.0d0)
                  j55z2a1(isig,i) = (0.0d0, 0.0d0)
                  j55a1a2(isig,i) = (0.0d0, 0.0d0)
                  j55a2a1(isig,i) = (0.0d0, 0.0d0)
               enddo
            enddo
         else
            do isig=-1,1,2  
               do i = 1,3
                   j55z1z2(isig,i)=j55z1z2(isig,i)+cte*jbz1z2(isig,i)
                   j55z2z1(isig,i)=j55z2z1(isig,i)+cte*jbz2z1(isig,i)
                   j55a1z2(isig,i)=j55a1z2(isig,i)+cte*jba1z2(isig,i)
                   j55a2z1(isig,i)=j55a2z1(isig,i)+cte*jba2z1(isig,i)
                   j55z1a2(isig,i)=j55z1a2(isig,i)+cte*jbz1a2(isig,i)
                   j55z2a1(isig,i)=j55z2a1(isig,i)+cte*jbz2a1(isig,i)
                   j55a1a2(isig,i)=j55a1a2(isig,i)+cte*jba1a2(isig,i)
                   j55a2a1(isig,i)=j55a2a1(isig,i)+cte*jba2a1(isig,i)
               enddo
            enddo
         endif

         do isig = -1,1,2
            do i = 1,3
               mv1(i) = j55a1a2(isig,i)*clr(3,1,isig)**3+
     1                  j55a2a1(isig,i)*clr(3,1,isig)**3+
     2                  j55z1z2(isig,i)*clr(3,1,isig)*clr(3,2,isig)**2+
     3                  j55z2z1(isig,i)*clr(3,1,isig)*clr(3,2,isig)**2+
     4                  j55a1z2(isig,i)*clr(3,2,isig)*clr(3,1,isig)**2+
     5                  j55a2z1(isig,i)*clr(3,2,isig)*clr(3,1,isig)**2+
     6                  j55z1a2(isig,i)*clr(3,2,isig)*clr(3,1,isig)**2+
     7                  j55z2a1(isig,i)*clr(3,2,isig)*clr(3,1,isig)**2
             
               mv2(i) = j55a1a2(isig,i)*clr(4,1,isig)**3+
     1                  j55a2a1(isig,i)*clr(4,1,isig)**3+
     2                  j55z1z2(isig,i)*clr(4,1,isig)*clr(4,2,isig)**2+
     3                  j55z2z1(isig,i)*clr(4,1,isig)*clr(4,2,isig)**2+
     4                  j55a1z2(isig,i)*clr(4,2,isig)*clr(4,1,isig)**2+
     5                  j55a2z1(isig,i)*clr(4,2,isig)*clr(4,1,isig)**2+
     6                  j55z1a2(isig,i)*clr(4,2,isig)*clr(4,1,isig)**2+
     7                  j55z2a1(isig,i)*clr(4,2,isig)*clr(4,1,isig)**2
            enddo
            matv(1,isig,1) = mv1(1)+mv1(2)+mv1(3)
            matv(2,isig,1)=  mv2(1)+mv2(2)+mv2(3)
         enddo
      endif                     !lpent

c**********************************************************************
c Calculation of processes qq -> V V -> Z Z A -------------

      if (lbox) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

c calculate born matrix element and store it in mborn for gauge test
         mborn=mat(1,-1,1)+mat(1,-1,2)+mat(1,-1,3)+mat(1,-1,4)+mat(1,-1,5)
     1     +mat(1,-1,6)+mat(1,-1,7)

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qz1(mu) = -qze(mu,1)
            minus_qa(mu)  = -qz(mu,1)
            minus_qz2(mu) = -p1(mu)-p2(mu)-minus_qz1(mu)-minus_qa(mu)
            minus_qz1z2(mu) =minus_qz2(mu)+minus_qz1(mu) 
            minus_qaz1(mu) =minus_qa(mu)+minus_qz1(mu) 
            minus_qaz2(mu) =minus_qa(mu)+minus_qz2(mu) 
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi**2/3.0d0+7.0d0)

c qq -> A A/Z (-> Z Z)
         call boxlinemm(p1,minus_qa,minus_qz1z2,p2,
     1    psi(1,-1,2),psi(1,-1,1),atau(1,1),azz(0,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qa,minus_qz1z2,p2,
     1    psi(1,1,2),psi(1,1,1),atau(1,1),azz(0,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))
         call boxlinemm(p1,minus_qa,minus_qz1z2,p2,
     1    psi(1,-1,2),psi(1,-1,1),atau(1,1),zzztens(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qa,minus_qz1z2,p2,
     1    psi(1,1,2),psi(1,1,1),atau(1,1),zzztens(0,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))

         call Wardtest_Box(atau(1,1), qz(0,1), zzztens(0,1), qzz(0,1),
     1         gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         if (wardidtest1.eq.1) then
            do isig=-1,1,2
               do k = 1,2
                  matv(k,isig,2) = (0.0d0, 0.0d0)
               enddo
            enddo   
         else
            do isig=-1,1,2
               mv1(1)=-(result(isig,1)+cte*result_born(isig,1))
               mv2(1)=-(result(isig,2)+cte*result_born(isig,2))
               do k = 1,2
                  matv(k,isig,2) = (mv1(1)*clr(ifl(k),1,isig)
     1                           +  mv2(1)*clr(ifl(k),2,isig))
     2                                 *clr(ifl(k),1,isig)
               enddo
            enddo   
         endif

c qq -> A/Z(-> Z Z) A
         call boxlinemm(p1,minus_qz1z2,minus_qa,p2,
     1    psi(1,-1,2),psi(1,-1,1),azz(0,1),atau(1,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qz1z2,minus_qa,p2,
     1    psi(1,1,2),psi(1,1,1),azz(0,1),atau(1,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))
         call boxlinemm(p1,minus_qz1z2,minus_qa,p2,
     1    psi(1,-1,2),psi(1,-1,1),zzztens(0,1),atau(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qz1z2,minus_qa,p2,
     1    psi(1,1,2),psi(1,1,1),zzztens(0,1),atau(1,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))

         call Wardtest_Box(azz(0,1), qzz(0,1), atau(1,1), qz(0,1),
     1         gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))

         if (wardidtest1.eq.1) then
            do isig=-1,1,2
               do k = 1,2
                  matv(k,isig,3) = (0.0d0, 0.0d0)
               enddo
            enddo   
         else
            do isig=-1,1,2
               mv1(1)=-(result(isig,1)+cte*result_born(isig,1))
               mv2(1)=-(result(isig,2)+cte*result_born(isig,2))
               do k = 1,2
                  matv(k,isig,3) = (mv1(1)*clr(ifl(k),1,isig)
     1                           +  mv2(1)*clr(ifl(k),2,isig))
     2                                 *clr(ifl(k),1,isig)
               enddo
            enddo
         endif

c qq -> Z A/Z (-> Z A)
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),ze(1,1),aazm(0,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,1,2),psi(1,1,1),ze(1,1),aazm(0,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),ze(1,1),zzam(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,1,2),psi(1,1,1),ze(1,1),zzam(0,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),ae(1,1),aazm(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,3),treetest,
     3    result_born(-1,3))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,1,2),psi(1,1,1),ae(1,1),aazm(0,1),
     2    1,scale,9,-1,boxtest,result(1,3),treetest,
     3    result_born(1,3))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),ae(1,1),zzam(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,4),treetest,
     3    result_born(-1,4))
         call boxlinemm(p1,minus_qz1,minus_qaz2,p2,
     1    psi(1,1,2),psi(1,1,1),ae(1,1),zzam(0,1),
     2    1,scale,9,-1,boxtest,result(1,4),treetest,
     3    result_born(1,4))

         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),zmu(1,1),aaze(0,1),
     2    -1,scale,3,1,gaugebox2,result(-1,5),gaugebox2_born,
     3    result_born(-1,5))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,1,2),psi(1,1,1),zmu(1,1),aaze(0,1),
     2    1,scale,9,-1,boxtest,result(1,5),treetest,
     3    result_born(1,5))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),zmu(1,1),zzae(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,6),treetest,
     3    result_born(-1,6))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,1,2),psi(1,1,1),zmu(1,1),zzae(0,1),
     2    1,scale,9,-1,boxtest,result(1,6),treetest,
     3    result_born(1,6))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),amu(1,1),aaze(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,7),treetest,
     3    result_born(-1,7))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,1,2),psi(1,1,1),amu(1,1),aaze(0,1),
     2    1,scale,9,-1,boxtest,result(1,7),treetest,
     3    result_born(1,7))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),amu(1,1),zzae(0,1),
     2    -1,scale,9,-1,boxtest,result(-1,8),treetest,
     3    result_born(-1,8))
         call boxlinemm(p1,minus_qz2,minus_qaz1,p2,
     1    psi(1,1,2),psi(1,1,1),amu(1,1),zzae(0,1),
     2    1,scale,9,-1,boxtest,result(1,8),treetest,
     3    result_born(1,8))

         call Wardtest_Box(ze(1,1), qze(0,1), aazm(0,1), qz2a(0),
     1         gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))
         call Wardtest_Box(zmu(1,1), qzmu(0,1), aaze(0,1), qz1a(0),
     1         gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
            do isig= -1,1,2
               do k = 1,2
                  matv(k,isig,4) = (0.0d0, 0.0d0)
               enddo
            enddo
         else
            do isig=-1,1,2
               do i=1,8
                  result(isig,i)=(result(isig,i)+cte*result_born(isig,i))
               enddo
            enddo   

            do isig= -1,1,2
               mv1(1)=result(isig,1)+result(isig,5)   
               mv2(1)=result(isig,3)+result(isig,7)
               mv3(1)=result(isig,2)+result(isig,6)
               mv4(1)=result(isig,4)+result(isig,8)   
               do k = 1,2
                  matv(k,isig,4) = mv1(1)*clr(ifl(k),1,isig)*clr(ifl(k),2,isig)
     1                        +  mv2(1)*clr(ifl(k),1,isig)**2
     1                        +  mv3(1)*clr(ifl(k),2,isig)**2
     1                        +  mv4(1)*clr(ifl(k),2,isig)*clr(ifl(k),1,isig)
               enddo   
            enddo
         endif

c qq -> A/Z(-> Z A) Z
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),aazm(0,1),ze(1,1),
     2    -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3    result_born(-1,1))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,1,2),psi(1,1,1),aazm(0,1),ze(1,1),
     2    1,scale,9,-1,boxtest,result(1,1),treetest,
     3    result_born(1,1))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),zzam(0,1),ze(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,2),treetest,
     3    result_born(-1,2))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,1,2),psi(1,1,1),zzam(0,1),ze(1,1),
     2    1,scale,9,-1,boxtest,result(1,2),treetest,
     3    result_born(1,2))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),aazm(0,1),ae(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,3),treetest,
     3    result_born(-1,3))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,1,2),psi(1,1,1),aazm(0,1),ae(1,1),
     2    1,scale,9,-1,boxtest,result(1,3),treetest,
     3    result_born(1,3))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,-1,2),psi(1,-1,1),zzam(0,1),ae(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,4),treetest,
     3    result_born(-1,4))
         call boxlinemm(p1,minus_qaz2,minus_qz1,p2,
     1    psi(1,1,2),psi(1,1,1),zzam(0,1),ae(1,1),
     2    1,scale,9,-1,boxtest,result(1,4),treetest,
     3    result_born(1,4))

         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),aaze(0,1),zmu(1,1),
     2    -1,scale,3,1,gaugebox2,result(-1,5),gaugebox2_born,
     3    result_born(-1,5))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,1,2),psi(1,1,1),aaze(0,1),zmu(1,1),
     2    1,scale,9,-1,boxtest,result(1,5),treetest,
     3    result_born(1,5))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),zzae(0,1),zmu(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,6),treetest,
     3    result_born(-1,6))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,1,2),psi(1,1,1),zzae(0,1),zmu(1,1),
     2    1,scale,9,-1,boxtest,result(1,6),treetest,
     3    result_born(1,6))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),aaze(0,1),amu(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,7),treetest,
     3    result_born(-1,7))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,1,2),psi(1,1,1),aaze(0,1),amu(1,1),
     2    1,scale,9,-1,boxtest,result(1,7),treetest,
     3    result_born(1,7))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,-1,2),psi(1,-1,1),zzae(0,1),amu(1,1),
     2    -1,scale,9,-1,boxtest,result(-1,8),treetest,
     3    result_born(-1,8))
         call boxlinemm(p1,minus_qaz1,minus_qz2,p2,
     1    psi(1,1,2),psi(1,1,1),zzae(0,1),amu(1,1),
     2    1,scale,9,-1,boxtest,result(1,8),treetest,
     3    result_born(1,8))

         call Wardtest_Box(aazm(0,1), qz2a(0), ze(1,1), qze(0,1),
     1         gaugebox1, gaugebox1_born, mborn, wardidtest1, als(1,1))
         call Wardtest_Box(aaze(0,1), qz1a(0), zmu(1,1), qzmu(0,1),
     1         gaugebox2, gaugebox2_born, mborn, wardidtest2, als(1,1))

         if (wardidtest1.eq.1 .or. wardidtest2.eq.1) then
            do isig= -1,1,2
               do k = 1,2
                  matv(k,isig,5) = (0.0d0, 0.0d0)
               enddo
            enddo
         else
            do isig=-1,1,2
               do i=1,8
                  result(isig,i)=(result(isig,i)+cte*result_born(isig,i))
               enddo
            enddo   
            do isig= -1,1,2
               mv1(1)=result(isig,1)+result(isig,5)   
               mv2(1)=result(isig,3)+result(isig,7)
               mv3(1)=result(isig,2)+result(isig,6)
               mv4(1)=result(isig,4)+result(isig,8)   
               do k = 1,2
                  matv(k,isig,5) = mv1(1)*clr(ifl(k),1,isig)*clr(ifl(k),2,isig)
     1                        +  mv2(1)*clr(ifl(k),1,isig)**2
     1                        +  mv3(1)*clr(ifl(k),2,isig)**2
     1                        +  mv4(1)*clr(ifl(k),2,isig)*clr(ifl(k),1,isig)
               enddo  
            enddo
         endif

      endif !lbox

c  Calculation of processes qq -> A/Z -> Z Z A ---------------------
      do isig = -1,1,2
         do k = 1,2
            matv(k,isig,6) = (0.0d0, 0.0d0)
            matv(k,isig,7) = (0.0d0, 0.0d0)
         enddo
      enddo


#endif /* WITH_NLO */
c----------------------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then
         do isig = -1,1,2
            do i = 1,7
               do k = 1,2
                  madgr(i,isig,k) = 0d0
               enddo
            enddo
         enddo

         do isig = -1,1,2
            do k = 1,2
               do i =1,7
                  madgr(i,isig,k) = tamp(i,isig,k)
               enddo
            enddo
         enddo

c print results from checks:

         do i=1,7
            do k = 1,2
               do isig = -1,1,2
               temp = dreal(mat(k,isig,i))**2+dimag(mat(k,isig,i))**2
               temp1 = dreal(madgr(i,isig,k))**2+dimag(madgr(i,isig,k))**2
               print*,'vbfnlo/madgraph(',i,k,isig,') = ',mat(k,isig,i)/madgr(i,isig,k)
               print*,'vbfnlo2/madgraph2(',i,k,isig,') = ',temp/temp1
               enddo
            enddo
         enddo
         
      endif                     ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto uuww and ddww

      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig = -1,1,2
            mm(k,isig) = 0
            do i = 1,7
               mm(k,isig) = mm(k,isig) +  mat(k,isig,i)
            enddo
            res(k)= res(k) + dreal(mm(k,isig))**2
     &                     + dimag(mm(k,isig))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

            if (lnlo) then
               mv(k,isig)=matv(k,isig,1)+matv(k,isig,2)+matv(k,isig,3)
     1                   +matv(k,isig,4)+matv(k,isig,5)+matv(k,isig,6)
     2                   +matv(k,isig,7)
            endif
            if (nlo.gt.0) then
               mv(k,isig) = als(1,1)*c2o4pi*
     1                  (mv(k,isig) + mm(k,isig)* cvirtVV)
            else if(nlo.lt.0) then
               mv(k,isig) = als(1,1)*c2o4pi*mv(k,isig)
            else
                mv(k,isig) = (0d0, 0d0)
            endif

           resv(k) = resv(k) + 2d0*dreal( mm(k,isig)   *
     1              conjg( mv(k,isig) ) )

         enddo

         if (nlo.ge.0) then
            res(k) = (res(k)+resv(k))*3d0 ! 3 is the color sum factor
         else
            res(k) = resv(k)*3d0 ! 3 is the color sum factor
         endif  
      enddo
      
      uuzza = res(1)
      ddzza = res(2)

      return
      end

