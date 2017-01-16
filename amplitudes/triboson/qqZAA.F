c****************subroutine QQZAA *************************************
c
c	Giuseppe Bozzi, <giuseppe@particle.uni-karlsruhe.de>
c	Last modified: Nov. 2010 (giuseppe)
C
C       Francisco Campanario, francam@particle.uni-karlsruhe.de
C       Include Virtual corrections
C       Last modified: Nov 2010 (Paco)
C
C 
C  QQZAA calculates the matrix elements**2 for 
C        q1 q1bar    ---->   Z A A  
C        Z ---->   f1-bar f1
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

      subroutine qqZAA(bos,pbar,fsign,nlo,N,uuzaa,ddzaa,ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      double precision  pbar(0:3,2+n_v)
      double precision  uuzaa, ddzaa, res(2),resv(2)
      double precision  p(0:3,3+n_v)
      double precision  fqa1(0:4,2), fqa2(0:4,2), fqze(0:4,2)
      integer bos,fsign(2+n_v),sign(3+n_v),nlo,mu,i,k,isig,N,ps_number
      integer ifl(2)
      double complex psize(2,-1:1,2), psiae(2,-1:1,2),
     1               psia1(2,-1:1,2),psia2(2,-1:1,2)
      double complex ja1(0:5,-1:1,2),ja2(0:5,-1:1,2)

      double complex temp, temp1
      double complex mat(2,-1:1,5)
      double complex madgr(5,-1:1,2)
      double complex mm(2,-1:1), mv(2,-1:1), zero
      double complex ma, mz, m1, m2
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

c   pentagons 
 
      integer i1,i2
      double complex cte
      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)
      double complex pentcontr1b(3), pentcontr2b(3), pentcontr3b(3)
      double complex pentcontr4b(3), pentcontr5b(3), pentcontr6b(3)
      double complex testg(3)
      Logical wardidtest1,wardidtest2,wardidtest3
      Logical wardidtest4,wardidtest5,wardidtest6


      real*8 p1(0:3),p2(0:3),minus_qa1(0:3),minus_qa2(0:3),minus_qz(0:3)
      real*8 minus_qa1z(0:3),minus_qa2z(0:3)


      double complex j5a1a2(3,-1:1,2),jba1a2(3,-1:1,2)
      double complex j5a2a1(3,-1:1,2),jba2a1(3,-1:1,2)
      double complex j5pba1a2(-1:1,3,2),jbpba1a2(-1:1,3,2)
      double complex j5pba2a1(-1:1,3,2),jbpba2a1(-1:1,3,2)
      double complex result(-1:1,4,2),result_born(-1:1,4,2)
      double complex mv1, mv2
      double complex  mv1m,mv2m

      double complex  matv(2,-1:1,5),matm(2,-1:1,5)
      double complex  matvbf(2,-1:1,3)


      double complex gaugebox1(2), gaugebox2(2), gaugebox3(2)
      double complex gaugebox4(2)
      
      double complex gaugebox1_born(2), gaugebox2_born(2), gaugebox3_born(2)
      double complex gaugebox4_born(2)


      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0           ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5) ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5) ! true pentagon contribution
!      lpt = (NLO.eq.-1) .or. (NLO.eq.-5) ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4) ! q^mu/m_V terms of Pentagon
      lpq=.false.
      lpent = lpt .or. lpq
c   Here the Full pentagon is the total contribution. The q terms sum to zero and
c in the shifted one is equal to the Full but with additional cancellations.
C Therefore, i change:
c     if(nlo.eq.-5) to if(nlo.eq.-8) 
c     if(nlo.eq.1) to if(nlo.eq.-5) 

cccccc
c    FOR CHEQUING
!      lnlo=.true.
!      lpent=.false.
!      NLO=1
 
cccccc
    

c     define flavors of external quarks for the 2 subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzaa: u = 3
         ifl(k) = 3
         k = 2                  ! ddzaa: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,5
                  mat(k,isig,i) = (0.0d0,0.0d0)
                  matv(k,isig,i) = (0.0d0,0.0d0)
            enddo
         enddo
      enddo

c define the internal momenta
      do i = 1,2+n_v
         sign(i) = fsign(i)
      enddo
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,N),
     &              ze(1,N),psize(1,isig,2),fqze(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa1(0,N),
     &              atau1(1,N), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa2(0,N),
     &              atau2(1,N), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,N),
     &              ze(1,N), psize(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa1(0,N),
     &              atau1(1,N), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa2(0,N),
     &              atau2(1,N), psia2(1,isig,1),fqa2(0,1))
         if (bos.eq.211) then ! only for bos=211 (e-e+aa)
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,N),
     &           ae(1,N),psiae(1,isig,2),fqze(0,2))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,N),
     &           ae(1,N), psiae(1,isig,1),fqze(0,1))
         endif
      enddo
     
*******************************************************
C Calculation of born matrix elements
*******************************************************

c *****************************************************
      if (bos.eq.211) then
c *****************************************************

c q qbar -> Z/A A A
      do isig = -1,1,2
         m1=-s1c(psiae(1,isig,2),atau1(1,N),.true.,isig,psia2(1,isig,1))
     1      -s1c(psiae(1,isig,2),atau2(1,N),.true.,isig,psia1(1,isig,1))
     2      -s1c(psia1(1,isig,2),ae(1,N),.true.,isig,psia2(1,isig,1))
     3      -s1c(psia1(1,isig,2),atau2(1,N),.true.,isig,psiae(1,isig,1))
     4      -s1c(psia2(1,isig,2),ae(1,N),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),atau1(1,N),.true.,isig,psiae(1,isig,1))

         m2=-s1c(psize(1,isig,2),atau1(1,N),.true.,isig,psia2(1,isig,1))
     1      -s1c(psize(1,isig,2),atau2(1,N),.true.,isig,psia1(1,isig,1))
     2      -s1c(psia1(1,isig,2),ze(1,N),.true.,isig,psia2(1,isig,1))
     3      -s1c(psia1(1,isig,2),atau2(1,N),.true.,isig,psize(1,isig,1))
     4      -s1c(psia2(1,isig,2),ze(1,N),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),atau1(1,N),.true.,isig,psize(1,isig,1))

         mat(1,isig,1) =  m1 * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + m2 * clr(3,2,isig) * clr(3,1,isig)**2 ! uu -> ZAA

         mat(2,isig,1) =  m1 * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + m2 * clr(4,2,isig) * clr(4,1,isig)**2 ! dd -> ZAA
         
!         print*, "mat(1,isig,1)",mat(1,isig,1)
!         print*, "mat(2,isig,1)",mat(2,isig,1)
!         stop
      enddo

c qqbar -> A A/Z (-> Z A)
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psia1(1,-1,1),fqa1(0,1), ja1(0,-1,1) )
      call curr6(1,psi(1,-1,2),p(0,2),
     &           psia2(1,-1,1),fqa2(0,1), ja2(0,-1,1) )

c contract with atoza/ztoza tensor
      do isig = -1,1,2
         m1 = -dotcc(aaze2(0,N),ja1(0,isig,1))
     1        -dotcc(aaze1(0,N),ja2(0,isig,1))
         m2 = -dotcc(zzae2(0,N),ja1(0,isig,1))
     1        -dotcc(zzae1(0,N),ja2(0,isig,1))
         do k = 1,2
            mat(k,isig,2) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
         enddo
      enddo

c qqbar -> A/Z(-> Z A) A
      call curr6(1,psia1(1,-1,2),fqa1(0,2),
     &           psi(1,-1,1),p(0,1), ja1(0,-1,2) )
      call curr6(1,psia2(1,-1,2),fqa2(0,2),
     &           psi(1,-1,1),p(0,1), ja2(0,-1,2) )

c contract with atoza/ztoza tensor
      do isig = -1,1,2
         m1 = -dotcc(aaze2(0,N),ja1(0,isig,2))
     1        -dotcc(aaze1(0,N),ja2(0,isig,2))
         m2 = -dotcc(zzae2(0,N),ja1(0,isig,2))
     1        -dotcc(zzae1(0,N),ja2(0,isig,2))
         do k = 1,2
            mat(k,isig,3) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
         enddo
      enddo

c qqbar -> A/Z -> Z A A 
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atozaa/ztozaa tensor
      do isig = -1,1,2
         ma = -dotcc(azaa(0,N),jqq(0,isig))
         mz = -dotcc(zzaa(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,4)  = ma*clr(ifl(k),1,isig)
            mat(k,isig,5)  = mz*clr(ifl(k),2,isig)
         enddo
      enddo

c *****************************************************
      elseif(bos.eq.511) then
c *****************************************************

c q qbar -> Z A A
      do isig = -1,1,2

         m1=-s1c(psize(1,isig,2),atau1(1,N),.true.,isig,psia2(1,isig,1))
     1      -s1c(psize(1,isig,2),atau2(1,N),.true.,isig,psia1(1,isig,1))
     2      -s1c(psia1(1,isig,2),ze(1,N),.true.,isig,psia2(1,isig,1))
     3      -s1c(psia1(1,isig,2),atau2(1,N),.true.,isig,psize(1,isig,1))
     4      -s1c(psia2(1,isig,2),ze(1,N),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),atau1(1,N),.true.,isig,psize(1,isig,1))

         mat(1,isig,1) = m1 * clr(3,2,isig) * clr(3,1,isig)**2 ! uu -> ZAA
         mat(2,isig,1) = m1 * clr(4,2,isig) * clr(4,1,isig)**2 ! dd -> ZAA
         
      enddo

c qqbar -> A/Z -> Z A A 
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with atozaa/ztozaa tensor
      do isig = -1,1,2
         ma = -dotcc(azaa(0,N),jqq(0,isig))
         mz = -dotcc(zzaa(0,N),jqq(0,isig))
         do k = 1,2
            mat(k,isig,4)  = ma*clr(ifl(k),1,isig)
            mat(k,isig,5)  = mz*clr(ifl(k),2,isig)
         enddo
      enddo

      endif

*******************************************************
*******************************************************
*
*
* Calculation of virtual contribution matrix elements
*
*
*******************************************************
*******************************************************
   
       cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)
!       cte=-cvirtVV -1d0 + 3*pi*(0d0,1d0)


c *****************************************************
      if (bos.eq.211) then
c *****************************************************
#ifdef WITH_NLO
      if (lpent) then
c Calculation of processes qq -> Z/A A A ---------------------
      scale=-2d0*dotrr(p(0,1),p(0,2))
!        Print*, "scale1",scale

            do mu = 0,3
               p1(mu) = pbar(mu,1)
               p2(mu) = pbar(mu,2)
               minus_qa1(mu) = -qa1(mu,1)
               minus_qa2(mu) = -qa2(mu,1)
               minus_qz(mu) = -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            enddo
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE for testing shift invariance change to -5
        if (NLO.eq.-8) then ! only true pentagon 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3,-1,1),pentcontr1b,jba1a2(3,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,1),testg,jba1a2(3,1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),aep(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(3,-1,2),testg,jba1a2(3,-1,2))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),atau2p(1,1),aep(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,2),testg,jba1a2(3,1,2))
c A1 Z A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2,-1,1),pentcontr2b,jba1a2(2,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,1),testg,jba1a2(2,1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),aep(1,1),atau2p(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(2,-1,2),testg,jba1a2(2,-1,2))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),aep(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,2),testg,jba1a2(2,1,2))
c Z A1 A2
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1,-1,1),pentcontr3b,jba1a2(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,1),testg,jba1a2(1,1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(1,-1,2),testg,jba1a2(1,-1,2))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),aep(1,1),atau1p(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,2),testg,jba1a2(1,1,2))
c A2 A1 Z
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3,-1,1),pentcontr4b,jba2a1(3,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,1),testg,jba2a1(3,1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),aep(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(3,-1,2),testg,jba2a1(3,-1,2))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),atau1p(1,1),aep(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,2),testg,jba2a1(3,1,2))
c A2 Z A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2,-1,1),pentcontr5b,jba2a1(2,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,1),testg,jba2a1(2,1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),aep(1,1),atau1p(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(2,-1,2),testg,jba2a1(2,-1,2))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),aep(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,2),testg,jba2a1(2,1,2))
c Z A2 A1
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1,-1,1),pentcontr6b,jba2a1(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,1),testg,jba2a1(1,1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),aep(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(1,-1,2),testg,jba2a1(1,-1,2))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),aep(1,1),atau2p(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,2),testg,jba2a1(1,1,2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 Z A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c Z A1 A2 
            call pent_gauge(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 Z
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 Z A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c Z A2 A1
            call pent_gauge(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)
	
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6

            if (.not.wardidtest1) then
                do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                       j5a1a2(i,i1,i2) = (0.0d0, 0.0d0)
                       j5a2a1(i,i1,i2) = (0.0d0, 0.0d0)
                     enddo
                  enddo
                 enddo
             pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1		 
c               pentcount = pentcount+1
            else
             do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                        j5a1a2(i,i1,i2)=j5a1a2(i,i1,i2)+cte*jba1a2(i,i1,i2)
                        j5a2a1(i,i1,i2)=j5a2a1(i,i1,i2)+cte*jba2a1(i,i1,i2)
                     enddo
                 enddo
              enddo
             pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1      
c               pentcount2 = pentcount2+1
            endif

c   q qbar -> A A Z

            do isig=-1,1,2

            mv1=0d0
            mv2=0d0
            mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
             mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
             mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)

           enddo

         matv(1,isig,1) =  mv1 * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matv(2,isig,1) =  mv1 * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

        matm(1,isig,1) =  mv1m * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matm(2,isig,1) =  mv1m * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

 
         enddo
         
            
!       Print*, 'True Pentagon matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE Should be -4
         elseif(NLO.eq.-4) then ! pentagon parts ~ boxes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
            call pent_to_box_2hel(p1, minus_qa1, minus_qa2, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), atau2(1,1),ze(1,1),
     2         xatau1, xatau2, xze, scale, j5pba1a2(-1,3,1), jbpba1a2(-1,3,1))
            call pent_to_box_2hel(p1, minus_qa1, minus_qa2, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), atau2(1,1),ae(1,1),
     2         xatau1, xatau2, xae, scale, j5pba1a2(-1,3,2), jbpba1a2(-1,3,2))
c A1 Z A2
            call pent_to_box_2hel(p1, minus_qa1, minus_qz, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), ze(1,1),atau2(1,1),
     2         xatau1, xze, xatau2, scale, j5pba1a2(-1,2,1), jbpba1a2(-1,2,1))
            call pent_to_box_2hel(p1, minus_qa1, minus_qz, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), ae(1,1),atau2(1,1),
     2         xatau1, xae, xatau2, scale, j5pba1a2(-1,2,2), jbpba1a2(-1,2,2))
c Z A1 A2     
            call pent_to_box_2hel(p1, minus_qz, minus_qa1, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau1(1,1),atau2(1,1),
     2         xze, xatau1, xatau2, scale, j5pba1a2(-1,1,1), jbpba1a2(-1,1,1))
            call pent_to_box_2hel(p1, minus_qz, minus_qa1, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau1(1,1),atau2(1,1),
     2         xae, xatau1, xatau2, scale, j5pba1a2(-1,1,2), jbpba1a2(-1,1,2))
c A2 A1 Z
            call pent_to_box_2hel(p1, minus_qa2, minus_qa1, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), atau1(1,1),ze(1,1),
     2         xatau2, xatau1, xze, scale, j5pba2a1(-1,3,1), jbpba2a1(-1,3,1))
            call pent_to_box_2hel(p1, minus_qa2, minus_qa1, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), atau1(1,1),ze(1,1),
     2         xatau2, xatau1, xae, scale, j5pba2a1(-1,3,2), jbpba2a1(-1,3,2))
c A2 Z A1     
            call pent_to_box_2hel(p1, minus_qa2, minus_qz, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), ze(1,1),atau1(1,1),
     2         xatau2, xze, xatau1, scale, j5pba2a1(-1,2,1), jbpba2a1(-1,2,1))
           call pent_to_box_2hel(p1, minus_qa2, minus_qz, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), ze(1,1),atau1(1,1),
     2         xatau2, xae, xatau1, scale, j5pba2a1(-1,2,2), jbpba2a1(-1,2,2))
c Z A2 A1     
            call pent_to_box_2hel(p1, minus_qz, minus_qa2, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau2(1,1),atau1(1,1),
     2         xze, xatau2, xatau1, scale, j5pba2a1(-1,1,1), jbpba2a1(-1,1,1))
           call pent_to_box_2hel(p1, minus_qz, minus_qa2, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau2(1,1),atau1(1,1),
     2         xae, xatau2, xatau1, scale, j5pba2a1(-1,1,2), jbpba2a1(-1,1,2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           do i1=-1,1,2
              do i2=1,2
                 j5a1a2(1,i1,i2) = j5pba1a2(i1,1,i2)+cte*jbpba1a2(i1,1,i2)
                 j5a1a2(2,i1,i2) = j5pba1a2(i1,2,i2)+cte*jbpba1a2(i1,2,i2)
                 j5a1a2(3,i1,i2) = j5pba1a2(i1,3,i2)+cte*jbpba1a2(i1,3,i2)
                 j5a2a1(1,i1,i2) = j5pba2a1(i1,1,i2)+cte*jbpba2a1(i1,1,i2)
                 j5a2a1(2,i1,i2) = j5pba2a1(i1,2,i2)+cte*jbpba2a1(i1,2,i2)
                 j5a2a1(3,i1,i2) = j5pba2a1(i1,3,i2)+cte*jbpba2a1(i1,3,i2)
               enddo 
            enddo
c   q qbar -> A A W

            do isig=-1,1,2

            mv1=0d0
            mv2=0d0
             mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
             mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
             mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)

           enddo


c$$$            Print*, 'j5a1a2(i1,-1,2)',(j5a1a2(i1,-1,2),i1=1,3)
c$$$            Print*, 'j5a1a2(i1,1,2)',(j5a1a2(i1,1,2),i1=1,3)
c$$$             Print*, 'j5a1a2(i1,-1,1)',(j5a1a2(i1,-1,1),i1=1,3)
c$$$             Print*, 'j5a1a2(i1,1,1)',(j5a1a2(i1,1,1),i1=1,2)
c$$$
c$$$           Print*, ' mv1', mv1
c$$$           Print*, ' mv2', mv2
c$$$           stop

         matv(1,isig,1) =  mv1 * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matv(2,isig,1) =  mv1 * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

         matm(1,isig,1) =  mv1m * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matm(2,isig,1) =  mv1m * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

 
            enddo

!     7    +mv1  ! This line when checking sum and remove if(nlo)
!       Print*, 'Full SUM matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!CHANGE for testing shift invariance change to 1
	  elseif(NLO.eq.-5) then ! FULL pentagon contributions NLO 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3,-1,1),pentcontr1b,jba1a2(3,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,1),testg,jba1a2(3,1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),ae(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(3,-1,2),testg,jba1a2(3,-1,2))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),atau2(1,1),ae(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,2),testg,jba1a2(3,1,2))
c A1 Z A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2,-1,1),pentcontr2b,jba1a2(2,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,1),testg,jba1a2(2,1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),ae(1,1),atau2(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(2,-1,2),testg,jba1a2(2,-1,2))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),ae(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,2),testg,jba1a2(2,1,2))
c Z A1 A2
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1,-1,1),pentcontr3b,jba1a2(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,1),testg,jba1a2(1,1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,9,-1,testg,j5a1a2(1,-1,2),testg,jba1a2(1,-1,2))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),atau1(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,2),testg,jba1a2(1,1,2))
c A2 A1 Z
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3,-1,1),pentcontr4b,jba2a1(3,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,1),testg,jba2a1(3,1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),ae(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(3,-1,2),testg,jba2a1(3,-1,2))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),atau1(1,1),ae(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,2),testg,jba2a1(3,1,2))
c A2 Z A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2,-1,1),pentcontr5b,jba2a1(2,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,1),testg,jba2a1(2,1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),ae(1,1),atau1(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(2,-1,2),testg,jba2a1(2,-1,2))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),ae(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,2),testg,jba2a1(2,1,2))
c Z A2 A1
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1,-1,1),pentcontr6b,jba2a1(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,1),testg,jba2a1(1,1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ae(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,9,-1,testg,j5a2a1(1,-1,2),testg,jba2a1(1,-1,2))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ae(1,1),atau2(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,2),testg,jba2a1(1,1,2))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 Z A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c Z A1 A2 
            call pent_gauge(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 Z
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 Z A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c Z A2 A1
            call pent_gauge(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)

  
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1             wardidtest4.and.wardidtest5.and.wardidtest6

 
            if (.not.wardidtest1) then
                do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                       j5a1a2(i,i1,i2) = (0.0d0, 0.0d0)
                       j5a2a1(i,i1,i2) = (0.0d0, 0.0d0)
                     enddo
                  enddo
                 enddo
	     pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
c               pentcount = pentcount+1
            else
             do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                        j5a1a2(i,i1,i2)=j5a1a2(i,i1,i2)+cte*jba1a2(i,i1,i2)
                        j5a2a1(i,i1,i2)=j5a2a1(i,i1,i2)+cte*jba2a1(i,i1,i2)
                     enddo
                 enddo
              enddo
	     pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
c               pentcount2 = pentcount2+1
            endif
c   q qbar -> A A W


            do isig=-1,1,2
            mv1=0d0
            mv2=0d0
            mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
             mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
             mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)
           enddo

         matv(1,isig,1) =  mv1 * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matv(2,isig,1) =  mv1 * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA


        matm(1,isig,1) =  mv1m * clr(3,1,isig)**3                 ! uu -> AAA
     1                  + mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA

         matm(2,isig,1) =  mv1m * clr(4,1,isig)**3                 ! dd -> AAA
     1                  + mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA


c$$$         print*, "mv1m",mv1m
c$$$         print*, "mv2m",mv2m
c$$$!         print*, "jba1a2",jba1a2
c$$$!         print*, "jba2a1",jba2a1
c$$$         print*, "mat(1,isig,1)",mat(1,isig,1)
c$$$         print*, "mat(2,isig,1)",mat(2,isig,1)
c$$$         print*, "matm(1,isig,1)",matm(1,isig,1)
c$$$         print*, "matm(2,isig,1)",matm(2,isig,1)
c$$$         print*, "ratio",matm(1,isig,1)/mat(1,isig,1)
c$$$         print*, "ratio",matm(2,isig,1)/mat(2,isig,1)
c$$$         print*, "matv(1,isig,1)",matv(1,isig,1)
c$$$         print*, "matv(2,isig,1)",matv(2,isig,1)
            enddo
c$$$            stop

!         Print*, 'Full Pen matv(k,-1,1)',matv(k,-1,1)
!         Print*, ''
         endif !NLO

      endif  !lpent

#endif /* WITH_NLO */
      
c Calculation of processes qq -> A V -> A A + l+neutriono ****************
#ifdef WITH_NLO
      

      if (lbox) then


!      if (.false.) then
         scale = -2d0*dotrr(p(0,1),p(0,2))

!	Print*, 'scale2',scale
c calculate born matrix element and store it in mborn for gauge check:
   

         do mu = 0,3
            p1(mu) = pbar(mu,1)
            p2(mu) = pbar(mu,2)
            minus_qa1(mu) = -qa1(mu,1)
            minus_qa2(mu) = -qa2(mu,1)
            minus_qz(mu) =  -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            minus_qa1z(mu) = -p1(mu)-p2(mu)+qa2(mu,1)
            minus_qa2z(mu) = -p1(mu)-p2(mu)+qa1(mu,1)
         enddo

!         cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)
c ZA2 A1
         call boxlinemm(p1,minus_qa2z,minus_qa1,p2,
     1      psi(1,-1,2),psi(1,-1,1),zzae2(0,1),atau1(1,1),
     2      -1,scale,3,1,gaugebox1,result(-1,1,1),gaugebox1_born,
     3      result_born(-1,1,1))
         call boxlinemm(p1,minus_qa2z,minus_qa1,p2,
     1      psi(1,1,2),psi(1,1,1),zzae2(0,1),atau1(1,1),
     2      1,scale,9,-1,testg,result(1,1,1),testg,
     3      result_born(1,1,1))
         call boxlinemm(p1,minus_qa2z,minus_qa1,p2,
     1      psi(1,-1,2),psi(1,-1,1),aaze2(0,1),atau1(1,1),
     2      -1,scale,9,-1,testg,result(-1,1,2),testg,
     3      result_born(-1,1,2))
         call boxlinemm(p1,minus_qa2z,minus_qa1,p2,
     1      psi(1,1,2),psi(1,1,1),aaze2(0,1),atau1(1,1),
     2      1,scale,9,-1,testg,result(1,1,2),testg,
     3      result_born(1,1,2))
c A1 ZA2 
         call boxlinemm(p1,minus_qa1,minus_qa2z,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau1(1,1),zzae2(0,1),
     2      -1,scale,3,1,gaugebox2,result(-1,2,1),gaugebox2_born,
     3      result_born(-1,2,1))
        call boxlinemm(p1,minus_qa1,minus_qa2z,p2,
     1      psi(1,1,2),psi(1,1,1),atau1(1,1),zzae2(0,1),
     2      1,scale,9,-1,testg,result(1,2,1),testg,
     3      result_born(1,2,1))
        call boxlinemm(p1,minus_qa1,minus_qa2z,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau1(1,1),aaze2(0,1),
     2      -1,scale,9,-1,testg,result(-1,2,2),testg,
     3      result_born(-1,2,2))
        call boxlinemm(p1,minus_qa1,minus_qa2z,p2,
     1      psi(1,1,2),psi(1,1,1),atau1(1,1),aaze2(0,1),
     2      1,scale,9,-1,testg,result(1,2,2),testg,
     3      result_born(1,2,2))
c ZA1 A2
         call boxlinemm(p1,minus_qa1z,minus_qa2,p2,
     1      psi(1,-1,2),psi(1,-1,1),zzae1(0,1),atau2(1,1),
     2      -1,scale,3,1,gaugebox3,result(-1,3,1),gaugebox3_born,
     3      result_born(-1,3,1))
         call boxlinemm(p1,minus_qa1z,minus_qa2,p2,
     1      psi(1,1,2),psi(1,1,1),zzae1(0,1),atau2(1,1),
     2      1,scale,9,-1,testg,result(1,3,1),testg,
     3      result_born(1,3,1))
         call boxlinemm(p1,minus_qa1z,minus_qa2,p2,
     1      psi(1,-1,2),psi(1,-1,1),aaze1(0,1),atau2(1,1),
     2      -1,scale,9,-1,testg,result(-1,3,2),testg,
     3      result_born(-1,3,2))
         call boxlinemm(p1,minus_qa1z,minus_qa2,p2,
     1      psi(1,1,2),psi(1,1,1),aaze1(0,1),atau2(1,1),
     2      1,scale,9,-1,testg,result(1,3,2),testg,
     3      result_born(1,3,2))
c A2 ZA1 
         call boxlinemm(p1,minus_qa2,minus_qa1z,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau2(1,1),zzae1(0,1),
     2      -1,scale,3,1,gaugebox4,result(-1,4,1),gaugebox4_born,
     3      result_born(-1,4,1))
         call boxlinemm(p1,minus_qa2,minus_qa1z,p2,
     1      psi(1,1,2),psi(1,1,1),atau2(1,1),zzae1(0,1),
     2      1,scale,9,-1,testg,result(1,4,1),testg,
     3      result_born(1,4,1))
         call boxlinemm(p1,minus_qa2,minus_qa1z,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau2(1,1),aaze1(0,1),
     2      -1,scale,9,-1,testg,result(-1,4,2),testg,
     3      result_born(-1,4,2))
         call boxlinemm(p1,minus_qa2,minus_qa1z,p2,
     1      psi(1,1,2),psi(1,1,1),atau2(1,1),aaze1(0,1),
     2      1,scale,9,-1,testg,result(1,4,2),testg,
     3      result_born(1,4,2))
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c WA2 A1
        call box_gauge(p1,minus_qa2z,minus_qa1,p2,
     1  scale,gaugebox1,gaugebox1_born,wardidtest1)
c A1 WA2 
         call box_gauge(p1,minus_qa1,minus_qa2z,p2,
     1   scale,gaugebox2,gaugebox2_born,wardidtest2)
c WA1 A2
         call box_gauge(p1,minus_qa1z,minus_qa2,p2,
     1      scale,gaugebox3,gaugebox3_born,wardidtest3)
c A2 WA1 
         call box_gauge(p1,minus_qa2,minus_qa1z,p2,
     2      scale,gaugebox4,gaugebox4_born,wardidtest4)
	 
         wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3
     1               .and.wardidtest4

     
         if (.NOT.wardidtest1) then

           do k = 1,2
              do isig=-1,1
                 matv(k,isig,2) = (0.0d0, 0.0d0)
                 matv(k,isig,3) = (0.0d0, 0.0d0)
              enddo
           enddo

         else
  
           do isig=-1,1,2
              do i1=1,4
                 do i2=1,2
                  result(isig,i1,i2)= -(result(isig,i1,i2)+cte*result_born(isig,i1,i2))
                 enddo
              enddo
!          A   Z*/A*  
               mv1=result(isig,2,2)+result(isig,4,2)  ! A*
               mv2=result(isig,2,1)+result(isig,4,1)  ! Z*
               mv1m=result_born(isig,2,2)+result_born(isig,4,2) ! A*
               mv2m=result_born(isig,2,1)+result_born(isig,4,1)  ! Z*

                  do k = 1,2 
                   matv(k,isig,2)=(mv1*clr(ifl(k),1,isig)
     1                           +  mv2*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)
           
                  matm(k,isig,2)=(mv1m*clr(ifl(k),1,isig)
     1                           +  mv2m*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)
                   enddo

!          Z*/A* A                   
               mv1=result(isig,1,2)+result(isig,3,2)
               mv2=result(isig,1,1)+result(isig,3,1)
               mv1m=result_born(isig,1,2)+result_born(isig,3,2)
               mv2m=result_born(isig,1,1)+result_born(isig,3,1)
                  do k = 1,2 
                     matv(k,isig,3)=(mv1*clr(ifl(k),1,isig)
     1                           +  mv2*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)
                    matm(k,isig,3)=(mv1m*clr(ifl(k),1,isig)
     1                           +  mv2m*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)
                   enddo


c$$$         print*, "mat(1,isig,2)",mat(1,isig,2)
c$$$         print*, "mat(2,isig,2)",mat(2,isig,2)
c$$$         print*, "matm(1,isig,2)",matm(1,isig,2)
c$$$         print*, "matm(2,isig,2)",matm(2,isig,2)
c$$$         print*, "ratio",-matm(1,isig,2)/mat(1,isig,2)
c$$$         print*, "ratio",-matm(2,isig,2)/mat(2,isig,2)
c$$$         print*, "mat(1,isig,3)",mat(1,isig,3)
c$$$         print*, "mat(2,isig,3)",mat(2,isig,3)
c$$$         print*, "matm(1,isig,3)",-matm(1,isig,3)
c$$$         print*, "matm(2,isig,3)",-matm(2,isig,3)
c$$$         print*, "ratio",-matm(1,isig,3)/mat(1,isig,3)
c$$$         print*, "ratio",-matm(2,isig,3)/mat(2,isig,3)
c$$$         print*, "matv(1,isig,2)",matv(1,isig,2)
c$$$         print*, "matv(2,isig,2)",matv(2,isig,2)
c$$$         print*, "matv(1,isig,3)",matv(1,isig,3)
c$$$         print*, "matv(2,isig,3)",matv(2,isig,3)
        enddo
c$$$        stop           
        
	 endif
      endif !lbox
#endif

      if (nlo.eq.4) then
c Calculation of processes qq -> Z/A -> A A Z ********************
C Factorize to the born
        do k = 1,2
        do isig = -1,1,2
           matv(k,isig,4) = (0.0d0, 0.0d0)     
           matv(k,isig,5) = (0.0d0, 0.0d0) 
c          matv(k,isig,4)=mat(k,isig,4)*(-8d0+4.0d0*pi*pi/3.0d0-cvirtVV)
c          matv(k,isig,5)=mat(k,isig,5)*(-8d0+4.0d0*pi*pi/3.0d0-cvirtVV)
        enddo
        enddo
        endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
CCCC For neutrino 
C
C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c *****************************************************
      elseif(bos.eq.511) then
c *****************************************************
#ifdef WITH_NLO
      if (lpent) then
c Calculation of processes qq -> Z/A A A ---------------------
      scale=-2d0*dotrr(p(0,1),p(0,2))
!        Print*, "scale1",scale

            do mu = 0,3
               p1(mu) = pbar(mu,1)
               p2(mu) = pbar(mu,2)
               minus_qa1(mu) = -qa1(mu,1)
               minus_qa2(mu) = -qa2(mu,1)
               minus_qz(mu) = -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            enddo
      
!      cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE for testing shift invariance change to -5
        if (NLO.eq.-8) then ! only true pentagon 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3,-1,1),pentcontr1b,jba1a2(3,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,1),testg,jba1a2(3,1,1))
c A1 Z A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2,-1,1),pentcontr2b,jba1a2(2,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,1),testg,jba1a2(2,1,1))
c Z A1 A2
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1,-1,1),pentcontr3b,jba1a2(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,1),testg,jba1a2(1,1,1))
c A2 A1 Z
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3,-1,1),pentcontr4b,jba2a1(3,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,1),testg,jba2a1(3,1,1))
c A2 Z A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2,-1,1),pentcontr5b,jba2a1(2,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,1),testg,jba2a1(2,1,1))
c Z A2 A1
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1,-1,1),pentcontr6b,jba2a1(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,1),testg,jba2a1(1,1,1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),zepp(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 Z A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),zepp(1,1),atau2p(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c Z A1 A2 
            call pent_gauge(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 Z
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),zepp(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 Z A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),zepp(1,1),atau1p(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c Z A2 A1
            call pent_gauge(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),zepp(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)
	
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6

            if (.not.wardidtest1) then
                do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                       j5a1a2(i,i1,i2) = (0.0d0, 0.0d0)
                       j5a2a1(i,i1,i2) = (0.0d0, 0.0d0)
                     enddo
                  enddo
                 enddo
             pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1		 
c               pentcount = pentcount+1
            else
             do i = 1,3
                  do i1=-1,1,2
                     do i2=1,1
                        j5a1a2(i,i1,i2)=j5a1a2(i,i1,i2)+cte*jba1a2(i,i1,i2)
                        j5a2a1(i,i1,i2)=j5a2a1(i,i1,i2)+cte*jba2a1(i,i1,i2)
                     enddo
                 enddo
              enddo
             pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1      
c               pentcount2 = pentcount2+1
            endif

c   q qbar -> A A Z

            do isig=-1,1,2

            mv1=0d0
            mv2=0d0
            mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
c            mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
c            mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)

           enddo

         matv(1,isig,1) =  mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matv(2,isig,1) =  mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA
         matm(1,isig,1) =  mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matm(2,isig,1) = mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

 
         enddo
         
            
!       Print*, 'True Pentagon matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE Should be -4
         elseif(NLO.eq.-4) then ! pentagon parts ~ boxes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
            call pent_to_box_2hel(p1, minus_qa1, minus_qa2, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), atau2(1,1),ze(1,1),
     2         xatau1, xatau2, xze, scale, j5pba1a2(-1,3,1), jbpba1a2(-1,3,1))
c A1 Z A2
            call pent_to_box_2hel(p1, minus_qa1, minus_qz, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), ze(1,1),atau2(1,1),
     2         xatau1, xze, xatau2, scale, j5pba1a2(-1,2,1), jbpba1a2(-1,2,1))
c Z A1 A2     
            call pent_to_box_2hel(p1, minus_qz, minus_qa1, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau1(1,1),atau2(1,1),
     2         xze, xatau1, xatau2, scale, j5pba1a2(-1,1,1), jbpba1a2(-1,1,1))
c A2 A1 Z
            call pent_to_box_2hel(p1, minus_qa2, minus_qa1, minus_qz, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), atau1(1,1),ze(1,1),
     2         xatau2, xatau1, xze, scale, j5pba2a1(-1,3,1), jbpba2a1(-1,3,1))
c A2 Z A1     
            call pent_to_box_2hel(p1, minus_qa2, minus_qz, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), ze(1,1),atau1(1,1),
     2         xatau2, xze, xatau1, scale, j5pba2a1(-1,2,1), jbpba2a1(-1,2,1))
c Z A2 A1     
            call pent_to_box_2hel(p1, minus_qz, minus_qa2, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), ze(1,1), atau2(1,1),atau1(1,1),
     2         xze, xatau2, xatau1, scale, j5pba2a1(-1,1,1), jbpba2a1(-1,1,1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           do i1=-1,1,2
              do i2=1,2
                 j5a1a2(1,i1,i2) = j5pba1a2(i1,1,i2)+cte*jbpba1a2(i1,1,i2)
                 j5a1a2(2,i1,i2) = j5pba1a2(i1,2,i2)+cte*jbpba1a2(i1,2,i2)
                 j5a1a2(3,i1,i2) = j5pba1a2(i1,3,i2)+cte*jbpba1a2(i1,3,i2)
                 j5a2a1(1,i1,i2) = j5pba2a1(i1,1,i2)+cte*jbpba2a1(i1,1,i2)
                 j5a2a1(2,i1,i2) = j5pba2a1(i1,2,i2)+cte*jbpba2a1(i1,2,i2)
                 j5a2a1(3,i1,i2) = j5pba2a1(i1,3,i2)+cte*jbpba2a1(i1,3,i2)
               enddo 
            enddo
c   q qbar -> A A W

            do isig=-1,1,2

            mv1=0d0
            mv2=0d0
             mv1m=0d0
            mv2m=0d0


           do i1=1,3 
c             mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
c             mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)

           enddo
           
c$$$
c$$$            Print*, 'nEUTRINO'
c$$$
c$$$            Print*, 'j5a1a2(i1,-1,2)',(j5a1a2(i1,-1,2),i1=1,3)
c$$$            Print*, 'j5a1a2(i1,1,2)',(j5a1a2(i1,1,2),i1=1,3)
c$$$            Print*, 'j5a1a2(i1,-1,1)',(j5a1a2(i1,-1,1),i1=1,3)
c$$$            Print*, 'j5a1a2(i1,1,1)',(j5a1a2(i1,1,1),i1=1,3)
c$$$
c$$$           Print*, ' mv1', mv1
c$$$           Print*, ' mv2', mv2
c$$$           stop


         matv(1,isig,1) = mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matv(2,isig,1) = mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA
         matm(1,isig,1) = mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matm(2,isig,1) = mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA

 
            enddo

!     7    +mv1  ! This line when checking sum and remove if(nlo)
!       Print*, 'Full SUM matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!CHANGE for testing shift invariance change to 1
	  elseif(NLO.eq.-5) then ! FULL pentagon contributions NLO 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3,-1,1),pentcontr1b,jba1a2(3,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         1,scale,9,-1,testg,j5a1a2(3,1,1),testg,jba1a2(3,1,1))
c A1 Z A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2,-1,1),pentcontr2b,jba1a2(2,-1,1))
            call penlinemmm(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(2,1,1),testg,jba1a2(2,1,1))
c Z A1 A2
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1,-1,1),pentcontr3b,jba1a2(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         1,scale,9,-1,testg,j5a1a2(1,1,1),testg,jba1a2(1,1,1))
c A2 A1 Z
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3,-1,1),pentcontr4b,jba2a1(3,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         1,scale,9,-1,testg,j5a2a1(3,1,1),testg,jba2a1(3,1,1))
c A2 Z A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2,-1,1),pentcontr5b,jba2a1(2,-1,1))
            call penlinemmm(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(2,1,1),testg,jba2a1(2,1,1))
c Z A2 A1
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1,-1,1),pentcontr6b,jba2a1(1,-1,1))
            call penlinemmm(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,1,2),psi(1,1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         1,scale,9,-1,testg,j5a2a1(1,1,1),testg,jba2a1(1,1,1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 Z
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),ze(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 Z A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qz(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),ze(1,1),atau2(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c Z A1 A2 
            call pent_gauge(p1(0),minus_qz(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 Z
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qz(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),ze(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 Z A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qz(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),ze(1,1),atau1(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c Z A2 A1
            call pent_gauge(p1(0),minus_qz(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),ze(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)

  
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1             wardidtest4.and.wardidtest5.and.wardidtest6

 
            if (.not.wardidtest1) then
                do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                       j5a1a2(i,i1,i2) = (0.0d0, 0.0d0)
                       j5a2a1(i,i1,i2) = (0.0d0, 0.0d0)
                     enddo
                  enddo
                 enddo
	     pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
c               pentcount = pentcount+1
            else
             do i = 1,3
                  do i1=-1,1,2
                     do i2=1,2
                        j5a1a2(i,i1,i2)=j5a1a2(i,i1,i2)+cte*jba1a2(i,i1,i2)
                        j5a2a1(i,i1,i2)=j5a2a1(i,i1,i2)+cte*jba2a1(i,i1,i2)
                     enddo
                 enddo
              enddo
	     pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
c               pentcount2 = pentcount2+1
            endif
c   q qbar -> A A W


            do isig=-1,1,2
            mv1=0d0
            mv2=0d0
            mv1m=0d0
            mv2m=0d0
           
           do i1=1,3 
c             mv1=mv1+ j5a1a2(i1,isig,2)+ j5a2a1(i1,isig,2)
             mv2=mv2+ j5a1a2(i1,isig,1)+ j5a2a1(i1,isig,1)
c             mv1m=mv1m+ jba1a2(i1,isig,2)+ jba2a1(i1,isig,2)
             mv2m=mv2m+ jba1a2(i1,isig,1)+ jba2a1(i1,isig,1)
           enddo

         matv(1,isig,1) = mv2 * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matv(2,isig,1) =  mv2 * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA
         matm(1,isig,1) =  mv2m * clr(3,2,isig) * clr(3,1,isig)**2  ! uu -> ZAA
         matm(2,isig,1) =  mv2m * clr(4,2,isig) * clr(4,1,isig)**2  ! dd -> ZAA


c$$$         print*, "mv1m",mv1m
c$$$         print*, "mv2m",mv2m
c$$$!         print*, "jba1a2",jba1a2
c$$$!         print*, "jba2a1",jba2a1
c$$$         print*, "mat(1,isig,1)",mat(1,isig,1)
c$$$         print*, "mat(2,isig,1)",mat(2,isig,1)
c$$$         print*, "matm(1,isig,1)",matm(1,isig,1)
c$$$         print*, "matm(2,isig,1)",matm(2,isig,1)
c$$$         print*, "ratio",matm(1,isig,1)/mat(1,isig,1)
c$$$         print*, "ratio",matm(2,isig,1)/mat(2,isig,1)
c$$$         print*, "matv(1,isig,1)",matv(1,isig,1)
c$$$         print*, "matv(2,isig,1)",matv(2,isig,1)
            enddo
c$$$            stop

!         Print*, 'Full Pen matv(k,-1,1)',matv(k,-1,1)
!         Print*, ''
         endif !NLO

      endif  !lpent

#endif /* WITH_NLO */
      
c Calculation of processes qq -> A V -> A A + l+neutriono ****************
#ifdef WITH_NLO
      

      if (lbox) then
             do k = 1,2
              do isig=-1,1
                 matv(k,isig,2) = (0.0d0, 0.0d0)
                 matv(k,isig,3) = (0.0d0, 0.0d0)
              enddo
           enddo
       endif !lbox
#endif


c Calculation of processes qq -> Z/A -> A A Z ********************
C Factorize to the born
        do k = 1,2
        do isig = -1,1,2
           matv(k,isig,4) = (0.0d0, 0.0d0)     
           matv(k,isig,5) = (0.0d0, 0.0d0) 
        enddo
        enddo


c *****************************************************
        endif   !!!bos 211 or 511
c *****************************************************

c#endif /* WITH_NLO */

c**************************************************************
c------------------------------------------------------------
c
c
c comparison with madgraph
c
c
c------------------------------------------------------------
c**************************************************************
c-----------------------------------------------------
c------------------------------------------------------
      if (ldebugm) then
         do isig = -1,1,2
            do k = 1,2
               do i =1,3
                  madgr(i,isig,k) = tamp(i,isig,k)
               enddo
               matvbf(k,isig,1) = mat(k,isig,1)
               matvbf(k,isig,2) = mat(k,isig,2)+mat(k,isig,3)
               matvbf(k,isig,3) = mat(k,isig,4)+mat(k,isig,5)
            enddo
         enddo

c print results from checks:

         do k = 1,2
            do isig = -1,1,2
               do i = 1,3
                 print*,'mat(',k,',',isig,',',i,') = ',matvbf(k,isig,i)
                 temp = dreal(matvbf(k,isig,i))**2+dimag(matvbf(k,isig,i))**2
                 print*,'madgr(',i,',',isig,',',k,') = ',madgr(i,isig,k)
                 temp1 = dreal(madgr(i,isig,k))**2+dimag(madgr(i,isig,k))**2
                 print*,'temp/temp1 = ',temp/temp1
               enddo
            enddo
         enddo
         
      endif                     ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto uuzaa and ddzaa

      do k = 1,2
         res(k) = 0
         resv(k) = 0
         do isig = -1,1,2
            mm(k,isig) = 0
            do i = 1,5
               mm(k,isig) = mm(k,isig) +  mat(k,isig,i)
            enddo
            res(k)= res(k) + dreal(mm(k,isig))**2
     &                     + dimag(mm(k,isig))**2

c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

            if (lnlo) then
               mv(k,isig)=matv(k,isig,1)+matv(k,isig,2)+matv(k,isig,3)!+matv(k,isig,4)+matv(k,isig,5)
  
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
      
      uuzaa = res(1)
      ddzaa = res(2)

      return
      end

