c****************subroutine QQZZJ *************************************
c
c	Michael Rauch, <rauch@particle.physik.uni-karlsruhe.de>
c	Last modified: June 2011
c  
c      Adding Virtual contributions:     
C     Francisco Campanario, francam@particle.uni-physik.uni-karlsruhe.de
c      Last modified:  April 2012
c
c
C  QQZZJ calculates the matrix elements**2 for 
C        q1 q1bar    ---->   Z Z g,   
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis
c  (comparison with madgraph amplitudes tampr requires to express these in
c   kart. basis also, even though they are given in the helicity basis
c
c  This code includes only real emission contributions, i.e.
c
c      return uuzz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuzzj,2:ddzzj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqZZj(pbar,sign, qbar, gsign, uuzzj, ddzzj,nlo)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
c#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/koppln.inc"

c electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program
      double precision fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  uuzzj, ddzzj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, j, k, isig
      integer  ifl(2), lh, l

      double precision fqz1(0:4,2),fqz2(0:4,2)
      double precision fqa1(0:4,2),fqa2(0:4,2)
      double precision pq(0:4,2)
      double precision pgZ1(0:4,2),pZ1g(0:4,2)
      double precision pgZ2(0:4,2),pZ2g(0:4,2)
      double precision pgA1(0:4,2),pA1g(0:4,2)
      double precision pgA2(0:4,2),pA2g(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psiz1(2,-1:1,2),psiz2(2,-1:1,2)
      double complex psia1(2,-1:1,2),psia2(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgZ1(2,-1:1,2,2)
      double complex braketgZ2(2,-1:1,2,2)
      double complex braketgA1(2,-1:1,2,2)
      double complex braketgA2(2,-1:1,2,2)
      double complex psiZ1g(2,-1:1,2,2)
      double complex psiZ2g(2,-1:1,2,2)
      double complex psiA1g(2,-1:1,2,2)
      double complex psiA2g(2,-1:1,2,2)
      double complex ja1(0:5,-1:1,4,2)
      double complex ja2(0:5,-1:1,4,2)
      double complex mat(2,-1:1,2,2), madgr(-1:1,-1:1,2)
      double complex matt(2,-1:1,2)
      double complex m1kb(2,-1:1), m2kb(2,-1:1), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex ma, mz, m1, m2, m3, m4
      double complex psi(2,-1:1,2), jqq(0:5,-1:1,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2)
      double complex jqqa(0:5,-1:1,2)

      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

      double complex matb(2,-1:1,2,4)
      double complex matv(2,-1:1,2,4)

c  helicity selection
      INTEGER h
      COMMON /hcount / h

c      logical ldebugm, linit
c      data linit /.true./
c      save linit
c      parameter(ldebugm = .false.)


      logical ldebug, ldebugm, linit, christophmods
      logical lnlo, lbox, lpq, lpt, lpent,lhex,lbox2,lhiggs
      logical addkmuknu
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false., christophmods=.false.)

c variables for virtual corrections
      DOUBLE PRECISION scale
      double precision c2,c2o4pi,oneo4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi, oneo4pi=1d0/4d0/pi)
c Variables for boxes
      real*8 mqza(0:3)
      real*8 td
      complex*16 emqza(0:3),empg(0:3)
      complex*16 boxcontr1(3),boxcontr1b(3)
      complex*16 boxcontr2(3),boxcontr2b(3)
      complex*16 boxcontrNoAbe1(3)
      complex*16 nothing(4)
      
      complex*16 j5g(5,2,2,-1:1),jbg(2,2,-1:1)
      complex*16 j5gg(5,2,2,-1:1),jbgg(2,2,-1:1)
      complex*16 j5gg2(5,2,2,-1:1),jbgg2(2,2,-1:1)
      complex*16 j5g1(5,2,2,-1:1),jbg1(2,2,-1:1)
      complex*16 jg5n(2,-1:1)
      complex*16 jg5ng(2,-1:1)
      complex*16 jg5ng2(2,-1:1)
      complex*16 jg5n1(2,-1:1)

      
c Variables for pentagons      
      real*8 mqa2(0:3)
      real*8 mqz(0:3)
      complex*16 cmqz(0:3)
      real*8 p1(0:3),p2(0:3)
      real*8 mpg(0:3),mpz(0:3)
      complex*16 pentcontr1(3),pentcontr1b(3)
      complex*16 pentcontr2(3),pentcontr2b(3)
      complex*16 pentcontr3(3),pentcontr3b(3)
      complex*16 pentcontr4(3),pentcontr4b(3)
      complex*16 pentcontr5(3),pentcontr5b(3)
      complex*16 pentcontr6(3),pentcontr6b(3)
      complex*16 pentcontr7(3),pentcontr7b(3)
      complex*16 pentcontr8(3),pentcontr8b(3)
c
      complex*16 j5a2g(8,4,2,-1:1),jba2g(4,2,-1:1)
      complex*16 j5ga2(8,4,2,-1:1),jbga2(4,2,-1:1)
      complex*16 j5a2g1(8,4,2,-1:1),jba2g1(4,2,-1:1)
      complex*16 j5ga21(8,4,2,-1:1),jbga21(4,2,-1:1)

      complex*16 j5zg(8,4,2,-1:1),jbzg(4,2,-1:1)
      complex*16 j5gz(8,4,2,-1:1),jbgz(4,2,-1:1)
      complex*16 j5zg1(8,4,2,-1:1),jbzg1(4,2,-1:1)
      complex*16 j5gz1(8,4,2,-1:1),jbgz1(4,2,-1:1)
      complex*16 j5zg1T(8,4,2,-1:1),jbzg1T(4,2,-1:1)
      complex*16 j5gz1T(8,4,2,-1:1),jbgz1T(4,2,-1:1)
      complex*16 j5zg1T2(8,4,2,-1:1),jbzg1T2(4,2,-1:1)
      complex*16 j5gz1T2(8,4,2,-1:1),jbgz1T2(4,2,-1:1)
c

c


c      complex*16 j5ma2g(8,4,2,-1:1),jbma2g(4,2,-1:1)
c      complex*16 j5mga2(8,4,2,-1:1),jbmga2(4,2,-1:1)
c
c      complex*16 j5ma2g1(8,4,2,-1:1),jbma2g1(4,2,-1:1)
c      complex*16 j5mga21(8,4,2,-1:1),jbmga21(4,2,-1:1)

c NoAbe
      complex*16 jg5a2(2,2,-1:1)
      complex*16 jg5z(2,2,-1:1)
      complex*16 jg5a21(2,2,-1:1)
      complex*16 jg5z1(2,2,-1:1)
      complex*16 jg5z1T(2,2,-1:1)
      complex*16 jg5z1T2(2,2,-1:1)

c      complex*16 jg5ma2(2,2,-1:1)
c      complex*16 jg5ma21(2,2,-1:1)
c      complex*16 jg5ma2(2,2,-1:1)
      complex*16 pentcontrNoAbe1(3),pentcontrNoAbe2(3)
      complex*16 pentemp(3)
      real*8      accuracyres


c  Colors Parameters       
      real*8 NF,inv2,inv3,CA,TR,tf,CF,eulergamma
      parameter (NF=5d0,inv2=1d0/2d0,inv3=1d0/3d0)
      parameter (CA=3d0)
      parameter (TR=1d0/2d0,tf=1d0/2d0)
      parameter (CF=4d0/3d0)
      parameter (eulergamma=5.772156649015328D-1)
      complex*16 cfc,cac,cfcac,pic2,cacgraph,convfact
      real*8 mans,mant,manu,musq
c     Variables for virtual
      complex*16 mv1,mv2,mv3,mv4,mv5,mv6
      complex*16 mmv2,mmv3,ratioFACT,mmv2b,mmv3b
      Complex*16 mmv2CA,mmv3CA,mmv4
      complex*16 mmv1,mmv1b,mmv1CA
      complex*16 Myeps(0:3,2) ! g in kartesian basis
c  Compare born from penlines, hexline,boxline 
C against NLO_real
       real*8 ratioComp
       complex*16 K_cont,K_contcheck,K_contcheck1,K_contcheck2,K_contcheck3,K_contcheck4
       integer Reno_inv
C  For adding Christoph factor from cancellation of 1/eps  divergences
       integer physToDiag1(3)
       common/Diag/physToDiag1
       complex*16 finiteubd,factortoadd,finiteubgcr,finiteubg
       real*8 s,t,mmu,mursq1
       complex*16 check(2,-1:1,2),check1(2,-1:1,2)
       integer GaugeAbe,GaugeNoAbe
c
       integer isig2,comp
       double complex cnf
      Integer Div,Gauge,PerformGauge,mutemp
      complex*16 mv1CA,mv2CA,mv3CA,mv1b,mv2b,mv3b
      complex*16 rvZERO(4,4),rbZERO(4,4)
      complex*16 lns,lnt,lnu,cteEps1

      Logical wardidtest1,wardidtest2,wardidtest3,wardidtest4
      Logical wardidtest5,wardidtest6,wardidtest7,wardidtest8
      integer N
      real*8 theta
      external theta

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i
      integer nlo

c Variables for fermion loops
      integer flgauge, flcomp
      double precision mf(4)
      integer cpl(4)
      double precision vv, temp_m2s
      double complex rfbox(6,4,4,2,2)
      double precision rgfbox(6,4,4,2,2) 
      double complex rfvert(2,3,4,2,2)
      double precision rgfvert(2,3,4,2,2) 
      double precision fcresult, fctreeresult
      integer*8 fchelicities, thirtytwo
      integer fcflags
      double complex fcamp(-2:2,-2:2), fckamp(2,2)
      parameter(thirtytwo=32)

      real*8 p12(0:3),s12
      real*8 minus_qa1z(0:3),minus_qa2z(0:3),minus_qaa(0:3)
      double complex result(-1:1,4,2),result_born(-1:1,4,2)
      integer poln
      external deltakro
      real*8 deltakro
      complex*16 v1(0:3), curr1(0:3,2),curr2(0:5,-1:1)
      double precision p12c(0:4,-1:1)
      integer l1,l2
      external sc1c
      complex*16 sc1c
!      parameter (lhiggs=.true.)
      double precision resv(2)
      complex*16 mmb(2,-1:1,2),mv(2,-1:1,2)

      complex*16 eff_curr1(0:3,-1:1),eff_curr2(0:3,-1:1)
      complex*16 eff_curr3(0:3,-1:1),eff_curr4(0:3,-1:1)



c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:


      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4)!.or.(NLO.eq.1).or.(NLO.eq.5).or.(NLO.eq.7).or.(NLO.eq.-7)  ! call boxline if T
      lpt = (NLO.eq.-5)!.or.(NLO.eq.1).or.(NLO.eq.7).or.(NLO.eq.-7)    ! true pentagon contribution
      !lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt !.or. lpq
c      Fermionloops
       lbox2=(NLO.eq.-8)
       lhiggs=lbox2

cFC      lnlo=.true.
cFC      lpent=.true.
cFC      lbox=.false.
cFC      lbox2=.false.

c       print*, "lbox",lbox
c       print*, "lpent",lpent
c       print*, "lbox2",lbox2

c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzzj: u = 3
         ifl(k) = 3
         k = 2                  ! ddzzj: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do l = 1,2
               do j = 1,2
                  mat(k,isig,l,j) = (0.0d0,0.0d0)
                  matb(k,isig,l,j) = (0.0d0,0.0d0)
                  matv(k,isig,l,j) = (0.0d0,0.0d0)
                  matv(k,isig,l,j+2) = (0.0d0,0.0d0)
               enddo
            enddo
         enddo
      enddo


      m1=0d0
      m2=0d0
      m3=0d0
      m4=0d0
      mv1=0d0
      mv2=0d0
      mv3=0d0
      mv4=0d0
      mv5=0d0
      mv6=0d0
      mv1b=0d0
      mv2b=0d0
      mv3b=0d0





c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
      enddo

c------------------------------------------------
c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      do isig =-1,1,2
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qp(0,1),
     &              wp(1,1), psiz1(1,isig,2),fqz1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qp(0,1),
     &              atau1(1,1), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qm(0,1),
     &              wm(1,1), psiz2(1,isig,2),fqz2(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qm(0,1),
     &              atau2(1,1), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,1),
     &              wp(1,1), psiz1(1,isig,1),fqz1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,1),
     &              atau1(1,1), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qm(0,1),
     &              wm(1,1), psiz2(1,isig,1),fqz2(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qm(0,1),
     &              atau2(1,1), psia2(1,isig,1),fqa2(0,1))
      enddo

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2        ! 2 gluon polarizations
         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors
c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

         do mu = 0,3                   !QCD gaugecheck
c            eps(mu,l) = qbar(mu)
           Myeps(mu,l)=eps(mu,l)
         enddo

c--- End of QCD gaugecheck -------------------------------------------
         do isig = -1,1,2       ! fermion helicity
c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     &           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig
            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     &           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash.
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qp(0,1),wp(1,1),braketgZ1(1,isig,1,l),pgZ1(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qp(0,1),wp(1,1),braketgZ1(1,isig,2,l),pgZ1(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qp(0,1),atau1(1,1),braketgA1(1,isig,1,l),pgA1(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qp(0,1),atau1(1,1),braketgA1(1,isig,2,l),pgA1(0,2))

            call ket2r(psiz1(1,isig,1),.false.,fqz1(0,1),isig,q,eps(0,l),
     &           psiZ1g(1,isig,1,l),pZ1g(0,1))
            call bra2r(psiz1(1,isig,2),.false.,fqz1(0,2),isig,q,eps(0,l),
     &           psiZ1g(1,isig,2,l),pZ1g(0,2))
            call ket2r(psia1(1,isig,1),.false.,fqa1(0,1),isig,q,eps(0,l),
     &           psiA1g(1,isig,1,l),pA1g(0,1))
            call bra2r(psia1(1,isig,2),.false.,fqa1(0,2),isig,q,eps(0,l),
     &           psiA1g(1,isig,2,l),pA1g(0,2))

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qm(0,1),wm(1,1),braketgZ2(1,isig,1,l),pgZ2(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qm(0,1),wm(1,1),braketgZ2(1,isig,2,l),pgZ2(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qm(0,1),atau2(1,1),braketgA2(1,isig,1,l),pgA2(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qm(0,1),atau2(1,1),braketgA2(1,isig,2,l),pgA2(0,2))

            call ket2r(psiz2(1,isig,1),.false.,fqz2(0,1),isig,q,eps(0,l),
     &           psiZ2g(1,isig,1,l),pZ2g(0,1))
            call bra2r(psiz2(1,isig,2),.false.,fqz2(0,2),isig,q,eps(0,l),
     &           psiZ2g(1,isig,2,l),pZ2g(0,2))
            call ket2r(psia2(1,isig,1),.false.,fqa2(0,1),isig,q,eps(0,l),
     &           psiA2g(1,isig,1,l),pA2g(0,1))
            call bra2r(psia2(1,isig,2),.false.,fqa2(0,2),isig,q,eps(0,l),
     &           psiA2g(1,isig,2,l),pA2g(0,2))
         enddo                   !fermion helicity
      enddo                      !gluon polarization

c qqbar -> Z Z g
      do l=1,2
         do isig = -1,1,2
cc ZZ
           m1=
     1      +s1c(braketg(1,isig,2,l),wp(1,1),.true.,isig,psiz2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),wm(1,1),.true.,isig,psiz1(1,isig,1))
     2      +s1c(psiz1(1,isig,2),wm(1,1),.true.,isig,braketg(1,isig,1,l))
     3      +s1c(psiz2(1,isig,2),wp(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psiz2g(1,isig,2,l),wp(1,1),.true.,isig,psi(1,isig,1))
     5      +s1c(psiz1g(1,isig,2,l),wm(1,1),.true.,isig,psi(1,isig,1))
cc ZA+AZ
           m2=
     1      +s1c(braketg(1,isig,2,l),wp(1,1),.true.,isig,psia2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),wm(1,1),.true.,isig,psia1(1,isig,1))
     2      +s1c(psiz1(1,isig,2),atau2(1,1),.true.,isig,braketg(1,isig,1,l))
     3      +s1c(psiz2(1,isig,2),atau1(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psiz2g(1,isig,2,l),atau1(1,1),.true.,isig,psi(1,isig,1))
     5      +s1c(psiz1g(1,isig,2,l),atau2(1,1),.true.,isig,psi(1,isig,1))
     1      +s1c(braketg(1,isig,2,l),atau1(1,1),.true.,isig,psiz2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),atau2(1,1),.true.,isig,psiz1(1,isig,1))
     2      +s1c(psia1(1,isig,2),wm(1,1),.true.,isig,braketg(1,isig,1,l))
     3      +s1c(psia2(1,isig,2),wp(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psia2g(1,isig,2,l),wp(1,1),.true.,isig,psi(1,isig,1))
     5      +s1c(psia1g(1,isig,2,l),wm(1,1),.true.,isig,psi(1,isig,1))
cc AA
           m3=
     1      +s1c(braketg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     2      +s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketg(1,isig,1,l))
     3      +s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psi(1,isig,1))
     5      +s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psi(1,isig,1))
           
           mat(1,isig,l,1) =  m1 * clr(3,2,isig)**2                 ! uu -> ZZ
     1                      + m2 * clr(3,2,isig)*clr(3,1,isig)
     1                      + m3 * clr(3,1,isig)**2
           mat(2,isig,l,1) =  m1 * clr(4,2,isig)**2                 ! dd -> ZZ
     1                      + m2 * clr(4,2,isig)*clr(4,1,isig)
     1                      + m3 * clr(4,1,isig)**2
        enddo
      enddo

c q qbar -> Z/A 

      do l = 1,2        ! 2 gluon polarizations
        call curr6(1,braketg(1,-1,2,l),pq(0,2),psi(1,-1,1),p(0,1),jqq(0,-1,1))
        call curr6(1,psi(1,-1,2),p(0,2),braketg(1,-1,1,l),pq(0,1),jqq(0,-1,2))

c contract with wz/wa tensor
        do isig = -1,1,2
          mat(1,isig,l,2) =  -dotcc(wz(1,1),jqq(0,isig,1))*clr(3,2,isig)
     &                       -dotcc(wz(1,1),jqq(0,isig,2))*clr(3,2,isig)
     &                       -dotcc(wa(1,1),jqq(0,isig,1))*clr(3,1,isig)
     &                       -dotcc(wa(1,1),jqq(0,isig,2))*clr(3,1,isig)
          mat(2,isig,l,2) =  -dotcc(wz(1,1),jqq(0,isig,1))*clr(4,2,isig)
     &                       -dotcc(wz(1,1),jqq(0,isig,2))*clr(4,2,isig)
     &                       -dotcc(wa(1,1),jqq(0,isig,1))*clr(4,1,isig)
     &                       -dotcc(wa(1,1),jqq(0,isig,2))*clr(4,1,isig)
        enddo
      enddo

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph
      if (ldebugm) then
         do j=1,2
            do isig = -1,1,2
               do k = 1,2
                  m1kb(k,isig) = 0d0 ! kartesian basis l = 1
                  m2kb(k,isig) = 0d0 ! kartesian basis l = 2
                  do lh = -1,1,2
                     madgr(isig,lh,k) = tampr(j,isig,lh,k)
                  enddo
               enddo
            enddo

            do isig = -1,1,2
               do k=1,2
                  m1kb(k,isig) = (madgr(isig,-1,k) -
     &                  madgr(isig,1,k))/sqrt(2d0)
                  m2kb(k,isig) = -im*(madgr(isig,-1,k) +
     &                  madgr(isig,1,k))/sqrt(2d0)
                  enddo
            enddo

c print results from checks:
            do k = 1,2
               do isig = -1,1,2
                  print*,'mat(',k,',',isig,',1,',j,') = ',mat(k,isig,1,j)
                  temp = dreal(mat(k,isig,1,j))**2+dimag(mat(k,isig,1,j))**2
                  print*,'m1kb(',k,',',isig,') = ',m1kb(k,isig)
                  temp1 = dreal(m1kb(k,isig))**2+dimag(m1kb(k,isig))**2
                  print*,'temp/temp1 = ',temp/temp1
                  print*,'mat(',k,',',isig,',2,',j,') = ',mat(k,isig,2,j)
                  temp = dreal(mat(k,isig,2,j))**2+dimag(mat(k,isig,2,j))**2
                  print*,'m2kb(',k,',',isig,') = ',m2kb(k,isig)
                  temp1 = dreal(m2kb(k,isig))**2+dimag(m2kb(k,isig))**2
                  print*,'temp/temp1 = ',temp/temp1
               enddo
            enddo
         enddo

      endif ! end of debugging block

c-----------------------------------------------------
c------------------------------------------------------
c sum the graphs, square them and map them onto uuzaaj and ddzaaj
      do k = 1,2
         res(k) = 0
         resv(k)=0
         do isig = -1,1,2
            do l = 1,2
               mm(k,isig,l) = 0
               mmb(k,isig,l) = 0
               do j = 1,2
                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,j)
                  mmb(k,isig,l) = mmb(k,isig,l) +  matb(k,isig,l,j)
               enddo
               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2
            enddo
         enddo
!         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor
      enddo
c      uuzzj = res(1)
c      ddzzj = res(2)

c      if (lnlo) then
         do l=1,2
            do isig=-1,1,2
               do isig2=1,2
                 mv(l,isig,isig2)=0d0
               enddo
            enddo

            do isig=-1,1,2
               do isig2=1,2 
                  do k=1,4
                   mv(l,isig,isig2)=mv(l,isig,isig2)+matv(l,isig,isig2,k)
                   enddo
                enddo
             enddo   
           enddo

       do l=1,2
           do  isig= -1,1,2
              do isig2=1,2 
                 if (nlo.eq.4) then
                   mv(l,isig,isig2) = als(1,1)*oneo4pi*
     1                (mv(1,isig,isig2) + mm(l,isig,isig2)*K_cont) 

                 elseif(nlo.eq.7) then
                  mv(l,isig,isig2) = als(1,1)*oneo4pi*
     1                (mv(1,isig,isig2) + mm(l,isig,isig2)*K_cont) 
     2                 +factortoadd*mmb(l,isig,isig2)

                 else if(nlo.lt.0) then
                  mv(l,isig,isig2) = als(1,1)*oneo4pi*
     1                (mv(1,isig,isig2)) 
     2                 +factortoadd*mmb(l,isig,isig2)

                  else
                     mv(l,isig,isig2) = (0d0, 0d0)
                endif

        resv(l) = resv(l) + 2d0*dreal(mm(l,isig,isig2)*conjg(mv(l,isig,isig2)))
        
        enddo
        enddo


           if (nlo.ge.0) then
c              print*, "HERE 2", nlo, res(l)

             res(l) = (res(l)+resv(l))*4d0*fpials(2)  ! C_2*3 is the color factor
            else
             res(l) = resv(l)*4d0*fpials(2)           ! C_2*3 is the color factor
           endif  

        enddo
c      endif

      uuzzj = res(1)
      ddzzj = res(2)

c      print*, "HERE 3", nlo, uuzaj,ddzaj

      return
      end
