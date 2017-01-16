c****************subroutine QQAAW *************************************
c
c	giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c	Last modified: December 2009 
C       Franciso Campanario, francam@particle.uni-karlsruhe.de
C       Include Virtual corrections
C       Last modified: September 2010
C
C  QQVAA calculates the matrix elements**2 for 
C     q1 q2    ---->   W+- A A   
C     W+- ---> f1-bar f2
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

      SUBROUTINE qqAAW(pbar,sign, nlo, N, final, bos,ps_number)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"

c     electroweak couplings are taken from KOPPLN

      DOUBLE PRECISION  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)
      DOUBLE PRECISION  pbar(0:3,2+n_v)
      DOUBLE PRECISION  res(2),resv(2),final
      DOUBLE PRECISION  p(0:3,2+n_v),
     2                  fqa1(0:4,2), fqa2(0:4,2),
     3                  pwa1(0:4),pwa2(0:4),
     4                  fq0(0:4,2)

      INTEGER sign(2+n_v), nlo, mu, i, k, isig, N,ps_number
      INTEGER  bos
      DOUBLE COMPLEX ja1(0:5,-1:1,2), ja2(0:5,-1:1,2),
     4               psiw0(2,-1:1,2),psia1(2,-1:1,2),psia2(2,-1:1,2)

      DOUBLE COMPLEX j5a1a2(3),j5a2a1(3), jba1a2(3), jba2a1(3)
      DOUBLE COMPLEX temp, temp1
      DOUBLE COMPLEX mat(2,-1:1,3), matv(2,-1:1,6)
      DOUBLE COMPLEX madgr(10,-1:1,2)
      DOUBLE COMPLEX mm(2), mv(2), zero
      DOUBLE COMPLEX ma, m1, m2, m3, m4
      DOUBLE COMPLEX mv1, mv2, mv3, mv4
      DOUBLE COMPLEX psi(2,-1:1,2), jqq(0:5,-1:1)
      DOUBLE COMPLEX dotcc, s1c, dotrc
      DOUBLE PRECISION dotrr, psumsq
      EXTERNAL dotcc, s1c, dotrr, dotrc, psumsq
      integer u(2),d(2)
      SAVE u,k,d
      PARAMETER (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      DOUBLE COMPLEX zm2i(2:3)
      SAVE zm2i

c  helicity selection

      INTEGER h
      COMMON /hcount / h

c variables for virtual corrections
      DOUBLE PRECISION scale
      double precision c2,c2o4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi)
      logical lnlo, lbox, lpq, lpt, lpent, ldebug, ldebugm, linit
      logical wardidtest1, wardidtest2, wardidtest3
      logical wardidtest4, wardidtest5, wardidtest6

      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)

c variables for implementation of the pentagon routine

      double precision minus_qa1(0:3),minus_qa2(0:3),minus_qw(0:3)     

      double precision minus_qa1w(0:3), minus_qa2w(0:3)
      double precision p1(0:3), p2(0:3) 

      double complex pentcontr1(3), pentcontr2(3), pentcontr3(3)
      double complex pentcontr4(3), pentcontr5(3), pentcontr6(3)
      double complex pentcontr1b(3), pentcontr2b(3), pentcontr3b(3)
      double complex pentcontr4b(3), pentcontr5b(3), pentcontr6b(3)


c variables for implementation of the box routine

      double complex result(-1:1,4),result_born(-1:1,4)
      double complex cte
      double complex gaugebox1(2), gaugebox2(2), gaugebox3(2)
      double complex gaugebox4(2)
      
      double complex gaugebox1_born(2), gaugebox2_born(2), gaugebox3_born(2)
      double complex gaugebox4_born(2)

c variables for test purposes


c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4) .or. (NLO.eq.1) .or. (NLO.eq.5)    ! call boxline if T
      lpt = (NLO.eq.1) .or. (NLO.eq.-5)   ! true pentagon contribution
      lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
C### CHANGE 
!      lpq=.false.
C### up to here
      lpent = lpt .or. lpq
c define flavors of external quarks for the 2 subprocesses

      if (linit) then
         linit = .false.
         if (bos.eq.311) then
            k = 1               ! udwaa
            u(1) = 3
            d(1) = 4
         elseif (bos.eq.411) then
            k = 2               ! duwaa
            u(2) = 4
            d(2) = 3
         endif
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif
c only left-handed quarks can appear     
      isig=-1
c      do isig = -1,1,2
         do i = 1,3
               mat(k,isig,i) = (0.0d0,0.0d0)
         enddo
         do i = 1,6
               matv(k,isig,i) = (0.0d0,0.0d0)
         enddo
c      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwa1(mu) = qw0(mu,N) + qa1(mu,N)
         pwa2(mu) = qw0(mu,N) + qa2(mu,N)
      enddo
      pwa1(4) = pwa1(0)**2-pwa1(1)**2-pwa1(2)**2-pwa1(3)**2
      pwa2(4) = pwa2(0)**2-pwa2(1)**2-pwa2(2)**2-pwa2(3)**2

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qw0(0,N),w0(1,N),
     &           psiw0(1,-1,2),fq0(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qw0(0,N),w0(1,N),
     &           psiw0(1,-1,1),fq0(0,1))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,2),fqa1(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,2),fqa2(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,1),fqa1(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,1),fqa2(0,1))
      
c***************************************************************
c
c        Tree-level contributions
c
c***************************************************************

c Calculation of processes qq -> W A A ---------------------

c   q qbar -> W A A
      m1= -s1c(psiw0(1,-1,2),atau1(1,N),.true.,-1,psia2(1,-1,1))
     &          *clr(u(k),1,-1)**2
     &    -s1c(psiw0(1,-1,2),atau2(1,N),.true.,-1,psia1(1,-1,1))
     &          *clr(u(k),1,-1)**2
     &    -s1c(psia1(1,-1,2),w0(1,N),.true.,-1,psia2(1,-1,1))
     &          *clr(u(k),1,-1)*clr(d(k),1,-1)
     &    -s1c(psia2(1,-1,2),w0(1,N),.true.,-1,psia1(1,-1,1))
     &          *clr(u(k),1,-1)*clr(d(k),1,-1)
     &    -s1c(psia1(1,-1,2),atau2(1,N),.true.,-1,psiw0(1,-1,1))
     &          *clr(d(k),1,-1)**2
     &    -s1c(psia2(1,-1,2),atau1(1,N),.true.,-1,psiw0(1,-1,1))
     &          *clr(d(k),1,-1)**2

      mat(k,-1,1) = m1 * clr(3,3,-1)

c Calculation of processes qq -> W A -> W A A -------------

c qq -> A W(-> W A)
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psia1(1,-1,1),fqa1(0,1), ja1(0,-1,1) )
      call curr6(-1,psi(1,-1,2),p(0,2),
     &           psia2(1,-1,1),fqa2(0,1), ja2(0,-1,1) )
      call curr6(-1,psia1(1,-1,2),fqa1(0,2),
     &           psi(1,-1,1),p(0,1), ja1(0,-1,2) )
      call curr6(-1,psia2(1,-1,2),fqa2(0,2),
     &           psi(1,-1,1),p(0,1), ja2(0,-1,2) )

c contract with atoww/ztoww tensor and
c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
      if (bos.eq.311) then
         m1 = -dotcc(wpwa2(0,N),ja1(0,-1,1))
         m2 = -dotcc(wpwa2(0,N),ja1(0,-1,2))
         m3 = -dotcc(wpwa1(0,N),ja2(0,-1,1))
         m4 = -dotcc(wpwa1(0,N),ja2(0,-1,2))
         m1 = m1+zm2i(3)*dotrc(pwa2(0),wpwa2(0,N))
     1        *dotrc(pwa2(0),ja1(0,-1,1))
         m2 = m2+zm2i(3)*dotrc(pwa2(0),wpwa2(0,N))
     1        *dotrc(pwa2(0),ja1(0,-1,2))
         m3 = m3+zm2i(3)*dotrc(pwa1(0),wpwa1(0,N))
     1        *dotrc(pwa1(0),ja2(0,-1,1))
         m4 = m4+zm2i(3)*dotrc(pwa1(0),wpwa1(0,N))
     1        *dotrc(pwa1(0),ja2(0,-1,2))
      elseif (bos.eq.411) then
         m1 = -dotcc(wmwa2(0,N),ja1(0,-1,1))
         m2 = -dotcc(wmwa2(0,N),ja1(0,-1,2))
         m3 = -dotcc(wmwa1(0,N),ja2(0,-1,1))
         m4 = -dotcc(wmwa1(0,N),ja2(0,-1,2))
         m1 = m1+zm2i(3)*dotrc(pwa2(0),wmwa2(0,N))
     1        *dotrc(pwa2(0),ja1(0,-1,1))
         m2 = m2+zm2i(3)*dotrc(pwa2(0),wmwa2(0,N))
     1        *dotrc(pwa2(0),ja1(0,-1,2))
         m3 = m3+zm2i(3)*dotrc(pwa1(0),wmwa1(0,N))
     1        *dotrc(pwa1(0),ja2(0,-1,1))
         m4 = m4+zm2i(3)*dotrc(pwa1(0),wmwa1(0,N))
     1        *dotrc(pwa1(0),ja2(0,-1,2))
      endif

      mat(k,-1,2) = ((m1+m3)*clr(u(k),1,-1)
     1     +(m2+m4)*clr(d(k),1,-1))* clr(3,3,-1) 

c Calculation of processes qq -> W -> W A A ---------------------
      call curr6(-1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),jqq(0,-1))

c contract with wtowaa tensor
      if (bos.eq.311) then
         ma = -dotcc(wpwaa(0,N),jqq(0,-1))
      elseif (bos.eq.411) then
         ma = -dotcc(wmwaa(0,N),jqq(0,-1))
      endif
 
       mat(k,-1,3)  = ma*clr(3,3,-1)
       
c**************************************************************
c
c    Virtual contributions
c
c**************************************************************


#ifdef WITH_NLO
      if (lpent) then
c Calculation of processes qq -> W A A ---------------------
      scale=-2d0*dotrr(p(0,1),p(0,2))
!        Print*, "scale1",scale

            do mu = 0,3
               p1(mu) = pbar(mu,1)
               p2(mu) = pbar(mu,2)
               minus_qa1(mu) = -qa1(mu,1)
               minus_qa2(mu) = -qa2(mu,1)
               minus_qw(mu) = -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            enddo
      
      cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE Should be -5
        if (NLO.eq.-5) then ! only true pentagon 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W	 
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),w0p(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3),pentcontr1b,jba1a2(3))
c A1 W A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qw(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),w0p(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2),pentcontr2b,jba1a2(2))
c W A1 A2
            call penlinemmm(p1(0),minus_qw(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0p(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1),pentcontr3b,jba1a2(1))
c A2 A1 W
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),w0p(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3),pentcontr4b,jba2a1(3))
c A2 W A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qw(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),w0p(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2),pentcontr5b,jba2a1(2))
c W A2 A1
            call penlinemmm(p1(0),minus_qw(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0p(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1),pentcontr6b,jba2a1(1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),atau2p(1,1),w0p(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 W A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qw(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1p(1,1),w0p(1,1),atau2p(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c W A1 A2 
            call pent_gauge(p1(0),minus_qw(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0p(1,1),atau1p(1,1),atau2p(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 W 
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),atau1p(1,1),w0p(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 W A1 
            call pent_gauge(p1(0),minus_qa2(0),minus_qw(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2p(1,1),w0p(1,1),atau1p(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c W A2 A1
            call pent_gauge(p1(0),minus_qw(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0p(1,1),atau2p(1,1),atau1p(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)
	
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6

            if (.not.wardidtest1) then
 
               do i = 1,3
                  j5a1a2(i) = (0.0d0, 0.0d0)
                  j5a2a1(i) = (0.0d0, 0.0d0)
                 enddo
             pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1		 
c               pentcount = pentcount+1
            else
             do i = 1,3
                j5a1a2(i)=j5a1a2(i)+cte*jba1a2(i)
                j5a2a1(i)=j5a2a1(i)+cte*jba2a1(i)
              enddo
             pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1      
c               pentcount2 = pentcount2+1
            endif

c   q qbar -> A A W

          mv1 = j5a1a2(3)*clr(u(k),1,-1)**2 !A1A2W
     2          + j5a1a2(2)*clr(u(k),1,-1)*clr(d(k),1,-1)!A1WA2
     3          + j5a1a2(1)*clr(d(k),1,-1)**2!WA1A2
     4          + j5a2a1(3)*clr(u(k),1,-1)**2 !A2A1W
     5          + j5a2a1(2)*clr(u(k),1,-1)*clr(d(k),1,-1) !A2WA1
     6          + j5a2a1(1)*clr(d(k),1,-1)**2 !WA2A1

 
         matv(k,-1,1) = (mv1)*clr(3,3,-1)
!       Print*, 'True Pentagon matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CHANGE Should be -4
         elseif(NLO.eq.-4) then ! pentagon parts ~ boxes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W
            call pent_to_box(p1, minus_qa1, minus_qa2, minus_qw, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), atau2(1,1),w0(1,1),
     2         xatau1, xatau2, xw, -1, scale, j5a1a2(3), jba1a2(3))

c A1 W A2
            call pent_to_box(p1, minus_qa1, minus_qw, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau1(1,1), w0(1,1),atau2(1,1),
     2         xatau1, xw, xatau2, -1, scale, j5a1a2(2), jba1a2(2))
c W A1 A2     
            call pent_to_box(p1, minus_qw, minus_qa1, minus_qa2, p2,
     1         psi(1,-1,2), psi(1,-1,1), w0(1,1), atau1(1,1),atau2(1,1),
     2         xw, xatau1, xatau2, -1, scale, j5a1a2(1), jba1a2(1))

c A2 A1 W
            call pent_to_box(p1, minus_qa2, minus_qa1, minus_qw, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), atau1(1,1),w0(1,1),
     2         xatau2, xatau1, xw, -1, scale, j5a2a1(3), jba2a1(3))
c A2 W A1     
            call pent_to_box(p1, minus_qa2, minus_qw, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), atau2(1,1), w0(1,1),atau1(1,1),
     2         xatau2, xw, xatau1, -1, scale, j5a2a1(2), jba2a1(2))
c W A2 A1     
            call pent_to_box(p1, minus_qw, minus_qa2, minus_qa1, p2,
     1         psi(1,-1,2), psi(1,-1,1), w0(1,1), atau2(1,1),atau1(1,1),
     2         xw, xatau2, xatau1, -1, scale, j5a2a1(1), jba2a1(1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            j5a1a2(1) = j5a1a2(1)+cte*jba1a2(1)
            j5a1a2(2) = j5a1a2(2)+cte*jba1a2(2)
            j5a1a2(3) = j5a1a2(3)+cte*jba1a2(3)
            j5a2a1(1) = j5a2a1(1)+cte*jba2a1(1)
            j5a2a1(2) = j5a2a1(2)+cte*jba2a1(2)
            j5a2a1(3) = j5a2a1(3)+cte*jba2a1(3)
c   q qbar -> A A W

          mv1 =  j5a1a2(3)*clr(u(k),1,-1)**2 !A1A2W
     2          + j5a1a2(2)*clr(u(k),1,-1)*clr(d(k),1,-1)!A1WA2
     3          + j5a1a2(1)*clr(d(k),1,-1)**2!WA1A2
     4          + j5a2a1(3)*clr(u(k),1,-1)**2 !A2A1W
     5          + j5a2a1(2)*clr(u(k),1,-1)*clr(d(k),1,-1) !A2WA1
     6          + j5a2a1(1)*clr(d(k),1,-1)**2 !WA2A1
!     7    +mv1  ! This line when checking sum and remove if(nlo)
  
         matv(k,-1,1) = (mv1)*clr(3,3,-1)
!       Print*, 'Full SUM matv(k,-1,1)',matv(k,-1,1)
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
!CHANGE Should be 1
	  elseif(NLO.eq.1) then ! FULL pentagon contributions NLO 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W
            call penlinemmm(p1(0),minus_qa1(0),minus_qa2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),w0(1,1),
     2         -1,scale,4,1,pentcontr1,j5a1a2(3),pentcontr1b,jba1a2(3))
c A1 W A2
            call penlinemmm(p1(0),minus_qa1(0),minus_qw(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),w0(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr2,j5a1a2(2),pentcontr2b,jba1a2(2))
c W A1 A2 
            call penlinemmm(p1(0),minus_qw(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,4,1,pentcontr3,j5a1a2(1),pentcontr3b,jba1a2(1))
c A2 A1 W
            call penlinemmm(p1(0),minus_qa2(0),minus_qa1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),w0(1,1),
     2         -1,scale,4,1,pentcontr4,j5a2a1(3),pentcontr4b,jba2a1(3))
c A2 W A1
            call penlinemmm(p1(0),minus_qa2(0),minus_qw(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),w0(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr5,j5a2a1(2),pentcontr5b,jba2a1(2))
c W A2 A1
            call penlinemmm(p1(0),minus_qw(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,4,1,pentcontr6,j5a2a1(1),pentcontr6b,jba2a1(1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W
           call pent_gauge(p1(0),minus_qa1(0),minus_qa2(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),atau2(1,1),w0(1,1),
     2         -1,scale,pentcontr1,pentcontr1b,wardidtest1)
c A1 W A2
            call pent_gauge(p1(0),minus_qa1(0),minus_qw(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau1(1,1),w0(1,1),atau2(1,1),
     2         -1,scale,pentcontr2,pentcontr2b,wardidtest2)
c W A1 A2
            call pent_gauge(p1(0),minus_qw(0),minus_qa1(0),minus_qa2(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0(1,1),atau1(1,1),atau2(1,1),
     2         -1,scale,pentcontr3,pentcontr3b,wardidtest3)
c A2 A1 W 
            call pent_gauge(p1(0),minus_qa2(0),minus_qa1(0),minus_qw(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),atau1(1,1),w0(1,1),
     2         -1,scale,pentcontr4,pentcontr4b,wardidtest4)
c A2 W A1
            call pent_gauge(p1(0),minus_qa2(0),minus_qw(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),atau2(1,1),w0(1,1),atau1(1,1),
     2         -1,scale,pentcontr5,pentcontr5b,wardidtest5)
c W A2 A1
            call pent_gauge(p1(0),minus_qw(0),minus_qa2(0),minus_qa1(0),p2(0),
     1         psi(1,-1,2),psi(1,-1,1),w0(1,1),atau2(1,1),atau1(1,1),
     2         -1,scale,pentcontr6,pentcontr6b,wardidtest6)
 
  
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1             wardidtest4.and.wardidtest5.and.wardidtest6

 
            if (.not.wardidtest1) then
               do i = 1,3
                  j5a1a2(i) = (0.0d0, 0.0d0)
                  j5a2a1(i) = (0.0d0, 0.0d0)
               enddo
	     pentcountm = pentcountm+1
             pentcountps(ps_number) = pentcountps(ps_number)+1
c               pentcount = pentcount+1
            else
             do i = 1,3
                j5a1a2(i)=j5a1a2(i)+cte*jba1a2(i)
                j5a2a1(i)=j5a2a1(i)+cte*jba2a1(i)
               enddo
	     pentcount2m = pentcount2m+1
             pentcount2ps(ps_number) = pentcount2ps(ps_number)+1
c               pentcount2 = pentcount2+1
            endif
c   q qbar -> A A W

          mv1 = j5a1a2(3)*clr(u(k),1,-1)**2 !A1A2W
     2          + j5a1a2(2)*clr(u(k),1,-1)*clr(d(k),1,-1)!A1WA2
     3          + j5a1a2(1)*clr(d(k),1,-1)**2!WA1A2
     4          + j5a2a1(3)*clr(u(k),1,-1)**2 !A2A1W
     5          + j5a2a1(2)*clr(u(k),1,-1)*clr(d(k),1,-1) !A2WA1
     6          + j5a2a1(1)*clr(d(k),1,-1)**2 !WA2A1

         matv(k,-1,1) = (mv1)*clr(3,3,-1)

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
            minus_qw(mu) =  -p1(mu)-p2(mu)+qa1(mu,1)+qa2(mu,1)
            minus_qa1w(mu) = -p1(mu)-p2(mu)+qa2(mu,1)
            minus_qa2w(mu) = -p1(mu)-p2(mu)+qa1(mu,1)
         enddo

         cte=(3*pi*(0d0,1d0)-4.0d0*pi*pi/3.0d0+7.0d0)


	if(bos.eq.311) then
c WA2 A1
         call boxlinemm(p1,minus_qa2w,minus_qa1,p2,
     1      psi(1,-1,2),psi(1,-1,1),wpwa2(0,1),atau1(1,1),
     2      -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3      result_born(-1,1))
c A1 WA2 
         call boxlinemm(p1,minus_qa1,minus_qa2w,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau1(1,1),wpwa2(0,1),
     2      -1,scale,3,1,gaugebox2,result(-1,2),gaugebox2_born,
     3      result_born(-1,2))
c WA1 A2
         call boxlinemm(p1,minus_qa1w,minus_qa2,p2,
     1      psi(1,-1,2),psi(1,-1,1),wpwa1(0,1),atau2(1,1),
     2      -1,scale,3,1,gaugebox3,result(-1,3),gaugebox3_born,
     3      result_born(-1,3))
c A2 WA1 
         call boxlinemm(p1,minus_qa2,minus_qa1w,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau2(1,1),wpwa1(0,1),
     2      -1,scale,3,1,gaugebox4,result(-1,4),gaugebox4_born,
     3      result_born(-1,4))

       elseif(bos.eq.411) then
c WA2 A1       
           call boxlinemm(p1,minus_qa2w,minus_qa1,p2,
     1      psi(1,-1,2),psi(1,-1,1),wmwa2(0,1),atau1(1,1),
     2      -1,scale,3,1,gaugebox1,result(-1,1),gaugebox1_born,
     3      result_born(-1,1))
c A1 WA2 
         call boxlinemm(p1,minus_qa1,minus_qa2w,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau1(1,1),wmwa2(0,1),
     2      -1,scale,3,1,gaugebox2,result(-1,2),gaugebox2_born,
     3      result_born(-1,2))
c WA1 A2
         call boxlinemm(p1,minus_qa1w,minus_qa2,p2,
     1      psi(1,-1,2),psi(1,-1,1),wmwa1(0,1),atau2(1,1),
     2      -1,scale,3,1,gaugebox3,result(-1,3),gaugebox3_born,
     3      result_born(-1,3))
c A2 WA1 
         call boxlinemm(p1,minus_qa2,minus_qa1w,p2,
     1      psi(1,-1,2),psi(1,-1,1),atau2(1,1),wmwa1(0,1),
     2      -1,scale,3,1,gaugebox4,result(-1,4),gaugebox4_born,
     3      result_born(-1,4))
  
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c WA2 A1
        call box_gauge(p1,minus_qa2w,minus_qa1,p2,
     1  scale,gaugebox1,gaugebox1_born,wardidtest1)
c A1 WA2 
         call box_gauge(p1,minus_qa1,minus_qa2w,p2,
     1   scale,gaugebox2,gaugebox2_born,wardidtest2)
c WA1 A2
         call box_gauge(p1,minus_qa1w,minus_qa2,p2,
     1      scale,gaugebox3,gaugebox3_born,wardidtest3)
c A2 WA1 
         call box_gauge(p1,minus_qa2,minus_qa1w,p2,
     2      scale,gaugebox4,gaugebox4_born,wardidtest4)
     
	
c$$$	if (.NOT.wardidtest1) then
c$$$	  call dt4M(p1,minus_qa2w,minus_qa1,p2,ratioD)
c$$$	 Print*, 'RATIOD1',ratioD
c$$$	 Print*, 'gaugebox1',gaugebox1
c$$$	 Print*, 'gaugebox1_born',gaugebox1_born
c$$$	 Print*, 'result(-1,1)',result(-1,1)
c$$$	 Print*, 'result_born(-1,1)',result_born(-1,1)
c$$$         Print*, ''	 
c$$$	 endif
c$$$	 
c$$$     	if (.NOT.wardidtest2) then
c$$$	Print*, 'WERDIDTEST2',WaRDIDTEST2
c$$$	  call dt4M(p1,minus_qa1,minus_qa2w,p2,ratioD)
c$$$		 Print*, 'RATIOD2',ratioD
c$$$	 Print*, 'gaugebox2',gaugebox2
c$$$	  Print*, 'gaugebox2_born',gaugebox2_born
c$$$	  Print*, 'result(-1,2)',result(-1,2)
c$$$	 Print*, 'result_born(-1,2)',result_born(-1,2)
c$$$	  Print*, ''	
c$$$	 endif
c$$$
c$$$     	if (.NOT.wardidtest3) then
c$$$	Print*, 'WERDIDTEST3',WaRDIDTEST3
c$$$	    call dt4M(p1,minus_qa1w,minus_qa2,p2,ratioD)
c$$$	 Print*, 'RATIOD3',ratioD
c$$$	 Print*, 'gaugebox3',gaugebox3
c$$$	 Print*, 'gaugebox3_born',gaugebox3_born
c$$$	 Print*, 'result(-1,3)',result(-1,3)
c$$$	 Print*, 'result_born(-1,3)',result_born(-1,3)
c$$$         Print*, ''	 	 
c$$$	 endif
c$$$	 
c$$$     	if (.NOT.wardidtest4) then
c$$$	Print*, 'WERDIDTEST4',WaRDIDTEST4
c$$$	  call dt4M(p1,minus_qa2,minus_qa1w,p2,ratioD)
c$$$	 Print*, 'RATIOD4',ratioD
c$$$	 Print*, 'gaugebox4',gaugebox4
c$$$	 Print*, 'gaugebox4_born',gaugebox4_born
c$$$	 Print*, 'result(-1,4)',result(-1,4)
c$$$	 Print*, 'result_born(-1,4)',result_born(-1,4)
c$$$	  Print*, ''
c$$$	 endif
	 
         wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3
     1               .and.wardidtest4

     
         if (.NOT.wardidtest1) then
	    matv(k,-1,2) = (0.0d0, 0.0d0)
	    matv(k,-1,3) = (0.0d0, 0.0d0)
	    matv(k,-1,4) = (0.0d0, 0.0d0)
	    matv(k,-1,5) = (0.0d0, 0.0d0)
         else
          mv1=-(result(-1,1)+cte*result_born(-1,1))
	  mv2=-(result(-1,2)+cte*result_born(-1,2))
          mv3=-(result(-1,3)+cte*result_born(-1,3))
	  mv4=-(result(-1,4)+cte*result_born(-1,4))
	  
	  matv(k,-1,2)=(mv1*clr(d(k),1,-1))*clr(3,3,-1)
          matv(k,-1,3)=(mv2*clr(u(k),1,-1))*clr(3,3,-1)  
          matv(k,-1,4)=(mv3*clr(d(k),1,-1))*clr(3,3,-1)    
          matv(k,-1,5)=(mv4*clr(u(k),1,-1))*clr(3,3,-1)
  
	 endif
      endif !lbox
#endif
  
c Calculation of processes qq -> W -> A A W ********************

      matv(k,-1,6) = (0.0d0, 0.0d0)      
c**************************************************************
c
c comparison with madgraph ************************************
c
c**************************************************************
c-----------------------------------------------------
c------------------------------------------------------
#ifdef WITH_MADGRAPH
c comparison with madgraph
      if (ldebugm) then

c         do isig = -1,1,2
            do i = 1,3

               madgr(i,isig,k) = 0d0

            enddo
c         enddo

c         do isig = -1,1,2
            do i =1,3

               madgr(i,isig,k) = tamp(i,isig,k)

            enddo
c         enddo

c print results from checks:

         do i=1,3
c            do isig = -1,-1,2
               print*,'mat(',k,',',isig,',',i,') = ',mat(k,isig,i)
               temp = dreal(mat(k,isig,i))**2+dimag(mat(k,isig,i))**2
               print*,'madgr(',i,',',isig,',',k,') = ',madgr(i,isig,k)
               temp1 = dreal(madgr(i,isig,k))**2+dimag(madgr(i,isig,k))**2
               print*,'temp/temp1 = ',temp/temp1
c            enddo
         enddo

      endif ! end of debugging block
#endif
c------------------------------------------------------
c**************************************************************
c
c sum the graphs, square them and map them onto qqwaa *********
c
c**************************************************************
      res(k) = 0
      resv(k) = 0

      mm(k) = 0
      do i = 1,3
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
     1        +matv(k,-1,5)+matv(k,-1,6)

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
      final = res(k)
  
      return
      end


c**************************************************************************** 

