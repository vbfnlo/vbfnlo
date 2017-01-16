c****************subroutine QQZAJ *************************************
c
c	Michael Rauch, <rauch@particle.physik.uni-karlsruhe.de>
c	Last modified: June 2011
c  
c      Adding Virtual contributions:     
C     Francisco Campanario, francam@particle.uni-physik.uni-karlsruhe.de
c      Last modified:  April 2012
c            
c
c
C  QQZAJ calculates the matrix elements**2 for 
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
c      return uuza = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuzaj,2:ddzaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqZAj(pbar,sign, qbar, gsign, nlo, N, uuzaj, ddzaj,
     &                 ps_number)

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
      double precision  uuzaj, ddzaj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, j, k, isig
      integer  ifl(2), lh, l

      double precision fqz1(0:4,2)
      double precision fqa1(0:4,2),fqa2(0:4,2)
      double precision pq(0:4,2)
      double precision pgZ1(0:4,2),pZ1g(0:4,2)

      double precision pgA1(0:4,2),pA1g(0:4,2)
      double precision pgA2(0:4,2),pA2g(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psiz1(2,-1:1,2)
      double complex psia1(2,-1:1,2),psia2(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgZ1(2,-1:1,2,2)

      double complex braketgA1(2,-1:1,2,2)
      double complex braketgA2(2,-1:1,2,2)
      double complex psiZ1g(2,-1:1,2,2)

      double complex psiA1g(2,-1:1,2,2)
      double complex psiA2g(2,-1:1,2,2)


      double complex mat(2,-1:1,2,2), madgr(-1:1,-1:1,2)

      double complex m1kb(2,-1:1), m2kb(2,-1:1), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex m1, m2, m3, m4
      double complex psi(2,-1:1,2), jqq(0:5,-1:1,2)


      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

      double complex matb(2,-1:1,2,4)
      double complex matv(2,-1:1,2,5)
c  helicity selection
      INTEGER h
      COMMON /hcount / h

      logical ldebug, ldebugm, linit, christophmods
      logical lnlo, lbox, lpt, lpent,lbox2,lhiggs,lvert

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

      real*8 p1(0:3),p2(0:3)
      real*8 mpg(0:3)
      complex*16 pentcontr1(3),pentcontr1b(3)
      complex*16 pentcontr2(3),pentcontr2b(3)
      complex*16 pentcontr3(3),pentcontr3b(3)
      complex*16 pentcontr4(3),pentcontr4b(3)
      complex*16 pentcontr5(3),pentcontr5b(3)
      complex*16 pentcontr6(3),pentcontr6b(3)
c
      complex*16 j5a2g(8,4,2,-1:1),jba2g(4,2,-1:1)
      complex*16 j5ga2(8,4,2,-1:1),jbga2(4,2,-1:1)
      complex*16 j5a2g1(8,4,2,-1:1),jba2g1(4,2,-1:1)
      complex*16 j5ga21(8,4,2,-1:1),jbga21(4,2,-1:1)
c


c


c NoAbe
      complex*16 jg5a2(2,2,-1:1)

      complex*16 jg5a21(2,2,-1:1)

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
      complex*16 mmv2,mmv3,ratioFACT,mmv2b
      Complex*16 mmv2CA,mmv4
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

c
       integer isig2,comp
       double complex cnf
      Integer Div,Performgauge,mutemp
      complex*16 mv1CA,mv2CA,mv3CA,mv1b,mv2b,mv3b

      complex*16 lns,lnt,lnu,cteEps1

      Logical wardidtest1,wardidtest2,wardidtest3,wardidtest4
      Logical wardidtest5,wardidtest6,wardidtest7,wardidtest8
      real*8 theta
      external theta

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i
      integer nlo, N, ps_number

c Variables for fermion loops

      double precision mf(4)
      integer cpl(4)
      double precision vv
      double complex rfbox(6,2,4,2,2)
      double precision rgfbox(6,2,4,2,2) 
      double complex rfvert(2,3,4,2,2)
      double precision rgfvert(2,3,4,2,2) 

      integer*8 thirtytwo


      parameter(thirtytwo=32)

      real*8 p12(0:3)
      real*8 minus_qa1z(0:3),minus_qa2z(0:3),minus_qaa(0:3)
      integer poln
      external deltakro
      real*8 deltakro
      complex*16 curr1(0:3,2),curr2(0:5,-1:1)
      double precision p12c(0:4,-1:1)
      integer l1,l2
      external sc1c
      complex*16 sc1c
!      parameter (lhiggs=.true.)
      double precision resv(2)
      complex*16 mmb(2,-1:1,2),mv(2,-1:1,2)
c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:


      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4)!.or.(NLO.eq.1).or.(NLO.eq.5).or.(NLO.eq.7).or.(NLO.eq.-7)  ! call boxline if T
      lpt = (NLO.eq.-5)!.or.(NLO.eq.1).or.(NLO.eq.7).or.(NLO.eq.-7)    ! true pentagon contribution
      !lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt !.or. lpq
c      Fermionloops
       lbox2=(NLO.eq.-6)
       lhiggs=lbox2
       lvert=lbox2

!      lnlo=.true.
!      lpent=.false.
!      lbox=.false.
!      lbox2=.true.

c       print*, "lbox",lbox
c       print*, "lpent",lpent
c       print*, "lbox2",lbox2


c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzaj: u = 3
         ifl(k) = 3
         k = 2                  ! ddzaj: d = 4
         ifl(k) = 4
      endif
cFC To check this

         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
       
cFC

      do k = 1,2
         do isig = -1,1,2
            do l = 1,2
               do j = 1,2
                  mat(k,isig,l,j) = (0.0d0,0.0d0)
                  matb(k,isig,l,j) = (0.0d0,0.0d0)
                  matv(k,isig,l,j) = (0.0d0,0.0d0)
                  matv(k,isig,l,j+2) = (0.0d0,0.0d0)
                  matv(k,isig,l,j+3) = (0.0d0,0.0d0)
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


c define the internal momenta(quarks)
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo

c define momenta of the gluon
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
     &              atau2(1,1), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,1),
     &              wp(1,1), psiz1(1,isig,1),fqz1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qp(0,1),
     &              atau1(1,1), psia1(1,isig,1),fqa1(0,1))
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
     &           qm(0,1),atau2(1,1),braketgA2(1,isig,1,l),pgA2(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qm(0,1),atau2(1,1),braketgA2(1,isig,2,l),pgA2(0,2))

            call ket2r(psia2(1,isig,1),.false.,fqa2(0,1),isig,q,eps(0,l),
     &           psiA2g(1,isig,1,l),pA2g(0,1))
            call bra2r(psia2(1,isig,2),.false.,fqa2(0,2),isig,q,eps(0,l),
     &           psiA2g(1,isig,2,l),pA2g(0,2))
         enddo                   !fermion helicity
      enddo                      !gluon polarization

c qqbar -> Z Z g
      do l=1,2
         do isig = -1,1,2
cc ZA+AZ
           m2=
     1      +s1c(braketg(1,isig,2,l),wp(1,1),.true.,isig,psia2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),atau2(1,1),.true.,isig,psiz1(1,isig,1))
     3      +s1c(psiz1(1,isig,2),atau2(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psia2(1,isig,2),wp(1,1),.true.,isig,braketg(1,isig,1,l))
     5      +s1c(psia2g(1,isig,2,l),wp(1,1),.true.,isig,psi(1,isig,1))
     6      +s1c(psiz1g(1,isig,2,l),atau2(1,1),.true.,isig,psi(1,isig,1))
cc AA
           m3=
     1      +s1c(braketg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     2      +s1c(braketg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     2      +s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketg(1,isig,1,l))
     3      +s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketg(1,isig,1,l))
     4      +s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psi(1,isig,1))
     5      +s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psi(1,isig,1))
           
           mat(1,isig,l,1) =                                        ! uu -> ZA
     1                      + m2 * clr(3,2,isig)*clr(3,1,isig)
     1                      + m3 * clr(3,1,isig)**2
           mat(2,isig,l,1) =                                        ! dd -> ZA
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


c**************************************************************
c
c    Virtual contributions
c
c**************************************************************
#ifdef WITH_NLO
      if (lnlo) then
C Som variables to perform checks of the virtual
C contributions
C To checl Factorization of Divergences
C  Div=1 -> 1/Eps;   Div=2 ->1/Eps^2
!         do Div=2,0,-1
        Div=0
	PerformGauge=0
	Reno_inv=0
C The only piece that is not checked with this variables 
C is the finite terms from the cancellation of the poles
C  against the Real part. Ask Christoph for this

c Some factors
      
      cfc=dcmplx(cf)
      cnf=dcmplx(nf)
      cac=dcmplx(ca)
      cacgraph=cac
      cfcac=(cfc-cac/2d0)
      pic2=dcmplx(pi**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Finite Terms from Pacos expansion     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do mu = 0,3
            p1(mu) = p(mu,1)
            p2(mu) = -p(mu,2)
            mpg(mu) = -q(mu)
         enddo
         
      mans=2d0*dotrr(p1,p2) !assuming massless (p1)**2,(p2)**2
      manu=2d0*dotrr(p2,mpg)!assuming massless (p2)**2,(mpg)**2
      mant=2d0*dotrr(p1,mpg)!assuming massless (p1)**2,(mpg)**2
      musq=mursq(1,1) ! renormalization scale
                     
                     
 100    lns=log(Abs(mans))-(0,1)*pi*theta(mans)-log(musq)
	lnt=log(Abs(mant))-(0,1)*pi*theta(mant)-log(musq)
	lnu=log(Abs(manu))-(0,1)*pi*theta(manu)-log(musq)

! This terms is a guess to make a substraction of each of the contributions to make
! them small, so that, we can decrease statistics in boxline,penlines, hexlines.
! So, it will give the M_v_born: K_cont*M_born
!       K_cont=(-2d0*pi*pi/3d0+50d0/9d0)*CAC-16d0/9d0*TF*cnf
!     1 	      +2d0*CFC*(-pi*pi+3d0)	
       K_cont=0d0
! This is the finite piece that has to be added to the finite contributions.
! It includes also the counterterms. We have factor out musq^-eps which goes to the 
! running of the alpha_s
     
      convfact=-2*EulerGamma*Inv3*cnf*TF + 2*Inv3*cnf*TF*Log(4*Pi) + CFC*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 
     &   + Log(musq)*(2*EulerGamma - 2*Log(4*Pi)) - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(mans)) + (0,2)*EulerGamma*Pi*theta(mans) + Log(4*Pi)
     & *(-3 + 2*EulerGamma + 2*Log(Abs(mans)) - (0,2)*Pi*theta(mans))) + 
     -     CAC*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*Pi**2 + Log(musq)*(EulerGamma - Log(4*Pi)) - 
     &  Inv2*Log(4*Pi)**2 + EulerGamma*Log(Abs(mans)) - 
     -     EulerGamma*Log(Abs(mant)) - EulerGamma*Log(Abs(manu)) + (0,1)*EulerGamma*Pi*(-theta(mans)
     &    + theta(mant) + theta(manu)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(mans)) + Log(Abs(mant)) + Log(Abs(manu))
     &  - (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))))

! The K_contcheck is used to check the renormalization invariance of the hexline
! hexlinenoAbe, penline, penlinenoAbe, boxline, boxlineNoABe
! Mv_(musq=s) + (convfact+ Kcontcheck= Kcontcheck4-K_cont)*mborn = Mv_(musq=s^prime) + (convfact+ Kcontcheck+ Kcontcheck4-K_cont)*mborn
! Basically, K_contcheck, include additional terms from the musq^-eps that is factor 
! out to define alpha_s(musq). Still, we have to subtract the counterterms to get
! the finite part of the hexlines,penline's independent of musq.
!
!
	If(Reno_inv.gt.0) then
      K_contcheck=-4*Inv2*Inv3*cnf*TF*Log(musq) + 
     1              CFC*Log(musq)*(3 - 2*EulerGamma + Log(musq) + 2*Log(4*Pi) - 
     2               2*Log(Abs(mans)) + (0,2)*Pi*theta(mans)) + 
     3              CAC*Log(musq)*(-EulerGamma + 11*Inv2*Inv3 + Inv2*Log(musq) + 
     4         Log((4*Pi*Abs(mans))/(Abs(mant)*Abs(manu))) - (0,1)*Pi*theta(mans) + 
     5            (0,1)*Pi*theta(mant) + (0,1)*Pi*theta(manu))
! This term cancel the musq dependence of the counterterms.
!  It should be added
! to prove renormalization invariance of the finite contribution.
        K_contcheck4 =4*Inv2*Inv3*cnf*TF*(EulerGamma + Log(musq) - Log(4*Pi)) + 
     -  Inv2*Inv3*(-11*EulerGamma + Log(4194304d0) - 11*Log(musq) + 11*Log(Pi))*CAC
!    This terms provide the full finite part that have to be added to the result.
!  Nevertheless, it is not the term that we have to add since musq^-eps must be 
!  factorize. It is only to check that convfact is right: 
!   K_contcheck1=convfact + K_contcheck
        K_contcheck1=  -2*EulerGamma*Inv3*cnf*TF - 2*Inv3*cnf*TF*Log(musq) + 2*Inv3*cnf*TF*Log(4*Pi) + 
     -  CFC*(-1 + 3*EulerGamma - EulerGamma**2 - Inv2*Inv3*Pi**2 + Log(musq)**2 - 
     -     Log(4*Pi)**2 - 2*EulerGamma*Log(Abs(mans)) + (0,2)*EulerGamma*Pi*theta(mans) + 
     -     Log(4*Pi)*(-3 + 2*EulerGamma + 2*Log(Abs(mans)) - (0,2)*Pi*theta(mans)) + 
     -     Log(musq)*(3 - 2*Log(Abs(mans)) + (0,2)*Pi*theta(mans))) + 
     -  CAC*((11 - 3*EulerGamma)*EulerGamma*Inv2*Inv3 - Inv2**2*Inv3*Pi**2 + 
     -     Inv2*Log(musq)**2 - Inv2*Log(4*Pi)**2 + EulerGamma*Log(Abs(mans)) - 
     -     EulerGamma*Log(Abs(mant)) - EulerGamma*Log(Abs(manu)) + 
     -     (0,1)*EulerGamma*Pi*(-theta(mans) + theta(mant) + theta(manu)) + 
     -     Log(4*Pi)*(EulerGamma - 11*Inv2*Inv3 - Log(Abs(mans)) + Log(Abs(mant)) + Log(Abs(manu)) - 
     -        (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))) + 
     -     Log(musq)*(11*Inv2*Inv3 + Log(Abs(mans)) - Log(Abs(mant)) - Log(Abs(manu)) + 
     -        (0,1)*Pi*(-theta(mans) + theta(mant) + theta(manu))))
      print*, ''
      print*, 'convfact',convfact
      Print*, 'K_contcheck',K_contcheck
      Print*, 'K_contcheck+convfact',K_contcheck+convfact
      Print*, 'K_contcheck1',K_contcheck1
      Print*, 'K_contcheck4',K_contcheck4
      Print*, 'K_contcheck+convfact-K_contcheck4',K_contcheck+convfact-K_contcheck4
!   This part proof renormalization invariance following Eq.2.38 of my paper 
        K_contcheck2=CAC*(1/2d0*(lns*lns-lnt*lnt-lnu*lnu)+
     &  3D0*(-2*lns+lnt+lnu))
	K_contcheck3=-CFC*(lns-3d0)*lns
!       Print*, 'K_contcheck2',K_contcheck2
        endif ! Finish the scale renormalization check

! This is to add the finite terms coming from the cancelation of the poles 
!  against the Real part. It comes from Christoph's code of WAj, WZj
      if(gsign.eq.1) then ! q q~

      s= 2d0*dotrr(pbar(0,1),pbar(0,2))
      t=-2d0*dotrr(pbar(0,1),qbar(0)) 
      mmu=-2d0*dotrr(pbar(0,2),qbar(0)) 
      mursq1=4.d0*pi*mursq(1,1)


      finiteubd=(25*als(1,N)*CA)/(9D0*Pi) + (5*als(1,N)*CF)/Pi - 
     &  (11*als(1,N)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,N)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,N)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,N)*CA*Pi)/24d0 - 
     &  (7*als(1,N)*CF*Pi)/12d0 - (8*als(1,N)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,N)*EulerGamma*Nf*TR)/(3D0*Pi) - 
     -  (3*als(1,N)*CA*Log(mursq1/S))/(4D0*Pi) + 
     &  (3*als(1,N)*CF*Log(mursq1/S))/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     -  (als(1,N)*CF*EulerGamma*Log(mursq1/S))/Pi - 
     &  (als(1,N)*CA*Log(mursq1/S)**2)/(4D0*Pi) + 
     &  (als(1,N)*CF*Log(mursq1/S)**2)/(2D0*Pi) + 
     -  (5*als(1,N)*CA*Log(-(mursq1/T)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,N)*Nf*TR*Log(-(mursq1/T)))/(6D0*Pi) + 
     -  (als(1,N)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     &  (5*als(1,N)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,N)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,N)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

       factortoadd=finiteubd

      elseif (gsign.eq.-1) then ! gluon iniciated
        
         if(phystoDiag1(1).eq.3) then ! g q / g q~
            if(phystoDiag1(2).eq.1) then
      s= 2d0*dotrr(qbar(0),pbar(0,1))
      t=-2d0*dotrr(pbar(0,1),pbar(0,2)) 
      mmU=-2d0*dotrr(qbar(0),pbar(0,2))  
      mursq1=4.d0*pi*mursq(1,1)
             else
      s= 2d0*dotrr(qbar(0),pbar(0,2))
      t=-2d0*dotrr(pbar(0,2),pbar(0,1)) 
      mmU=-2d0*dotrr(qbar(0),pbar(0,1))   
      mursq1=4.d0*pi*mursq(1,1)
            endif

      finiteubgcr=(25*als(1,N)*CA)/(9D0*Pi) + (5*als(1,N)*CF)/Pi - 
     &  (11*als(1,N)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,N)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,N)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,N)*CA*Pi)/24d0 - (7*als(1,N)*CF*Pi)/12d0 - 
     &  (8*als(1,N)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,N)*EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*als(1,N)*CA*Log(mursq1/S))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,N)*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (als(1,N)*CA*Log(mursq1/S)**2)/(4D0*Pi) - 
     &  (3*als(1,N)*CA*Log(-(mursq1/T)))/(4D0*Pi) + 
     &  (3*als(1,N)*CF*Log(-(mursq1/T)))/(2D0*Pi) + 
     -  (als(1,N)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,N)*CF*EulerGamma*Log(-(mursq1/T)))/Pi - 
     &  (als(1,N)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (als(1,N)*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + 
     &  (5*als(1,N)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,N)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,N)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

       factortoadd=finiteubgcr
            
         elseif(phystoDiag1(2).eq.3) then ! q g/ q~ g
            if(phystoDiag1(1).eq.1) then
      s=2d0*dotrr(pbar(0,1),qbar(0))      
      t=-2d0*dotrr(pbar(0,1),pbar(0,2))
      mmU=-2d0*dotrr(qbar(0),pbar(0,2)) 
      mursq1=4.d0*pi*mursq(1,1)
             else
      s= 2d0*dotrr(pbar(0,2),qbar(0))
      t=-2d0*dotrr(pbar(0,2),pbar(0,1))   
      mmU=-2d0*dotrr(qbar(0),pbar(0,1)) 
      mursq1=4.d0*pi*mursq(1,1)
            endif

      finiteubg=(25*als(1,N)*CA)/(9D0*Pi) + 
     &  (5*als(1,N)*CF)/Pi - (11*als(1,N)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,N)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,N)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,N)*CA*Pi)/24d0 - (7*als(1,N)*CF*Pi)/12d0 - 
     &  (8*als(1,N)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,N)*EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*als(1,N)*CA*Log(mursq1/S))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,N)*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (als(1,N)*CA*Log(mursq1/S)**2)/(4D0*Pi) - 
     &  (3*als(1,N)*CA*Log(-(mursq1/T)))/(4D0*Pi) + 
     &  (3*als(1,N)*CF*Log(-(mursq1/T)))/(2D0*Pi) + 
     -  (als(1,N)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,N)*CF*EulerGamma*Log(-(mursq1/T)))/Pi - 
     &  (als(1,N)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (als(1,N)*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + 
     &  (5*als(1,N)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,N)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,N)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

       factortoadd=finiteubg

          endif

      endif     ! Finish gsign

      factortoadd=factortoadd/2d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                 Here start the virtual staff
C
!	 LPENT=.FALSE.
!	 LBOX=.FALSE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc


CCCCCCCCCCCCCCCCCCCCCCCCC
C Pentagons
CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC
C the results are store in rvga2(x,y,l)
C The name indicate the order of gluon and the photons
C x: Number of graphs, 
C x=1 collect results proportional CF
C x=2 collect results proportional to CF-1/2CA
C y: Indicates the position of the Z. This will affect 
C    the couplings
C l: Polarization of the gluonl l=1,l=2
CCCCCCCCCCCCCCCCCCCCCCCCC
C qq-> A  Z g
CCCCCCCCCCCCCCCCCCCCCCCCC
       If(lpent) then

cFC
c  Reno_inv=1 Checked!!! CA,CF color structure checked
c  Div=1,2  Checked!! Factorization of divergences,CA,CF color checked
c  Performed_gauge=1,2,3 ! Checked
c
c
cFC 


      scale= musq!mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
      
            do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
	       mpg(mu) = -q(mu)
c photon
               mqa2(mu) = -qm(mu,1)
c Z,A*
               mqz(mu) = -qp(mu,1)
	   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
! Check gauge invariance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If(PerformGauge.eq.1) then
         do mu = 0,3
	    Myeps(mu,1)=mpg(mu)
	    Myeps(mu,2)=mpg(mu)
	  enddo
      elseif(PerformGauge.eq.2) then
         do mu = 0,3
	     mutemp=mu+1
c Polarization of Z,A*
	    wp(mutemp,1) = mqz(mu)
	    atau1(mutemp,1) = mqz(mu)
	 enddo
      elseif(PerformGauge.eq.3) then
         do mu = 0,3
	    mutemp=mu+1
c Polariazation A
	    atau2(mutemp,1) = mqa2(mu)
	 enddo
       endif
CCCCCCCCCCCCCCCCCCCCCCCCC
C  Abelian PART
CCCCCCCCCCCCCCCCCCCCCCCC
C
C A2 g Z
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2
       call penlineABETotal(p1,mqa2,mpg,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),Myeps(0,isig2),wp(1,1),isig,scale,1,3,
     2    comp*4,comp,pentcontr1,j5a2g(1,3,isig2,isig),pentcontr1b
     3    ,jba2g(3,isig2,isig),Div)

        comp=-1

       call penlineABETotal(p1,mqa2,mpg,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),Myeps(0,isig2),atau1(1,1),isig,scale,1,3,
     2    comp*4,comp,pentemp,j5a2g1(1,3,isig2,isig),pentemp
     3    ,jba2g1(3,isig2,isig),Div)

         enddo
       enddo   

C
C g A2 Z
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

       call penlineABETotal(p1,mpg,mqa2,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),atau2(1,1),wp(1,1),isig,scale,1,2,
     2    comp*4,comp,pentcontr2,j5ga2(1,3,isig2,isig),pentcontr2b
     3    ,jbga2(3,isig2,isig),Div)
 
       comp=-1

       call penlineABETotal(p1,mpg,mqa2,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),atau2(1,1),atau1(1,1),isig,scale,1,2,
     2    comp*4,comp,pentemp,j5ga21(1,3,isig2,isig),pentemp
     3    ,jbga21(3,isig2,isig),Div)

         enddo
       enddo   
C g A2 Z

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C A2 Z g
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

       call penlineABETotal(p1,mqa2,mqz,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),wp(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentcontr3,j5a2g(1,2,isig2,isig),pentcontr3b
     3    ,jba2g(2,isig2,isig),Div)
 
        comp=-1

       call penlineABETotal(p1,mqa2,mqz,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),atau1(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentemp,j5a2g1(1,2,isig2,isig),pentemp
     3    ,jba2g1(2,isig2,isig),Div)
 
         enddo
       enddo   
         
C
C g Z A2
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

       call penlineABETotal(p1,mpg,mqz,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),wp(1,1),atau2(1,1),isig,scale,1,2,
     2    comp*4,comp,pentcontr4,j5ga2(1,2,isig2,isig),pentcontr4b
     3    ,jbga2(2,isig2,isig),Div)
 
       comp=-1

       call penlineABETotal(p1,mpg,mqz,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),atau1(1,1),atau2(1,1),isig,scale,1,2,
     2    comp*4,comp,pentemp,j5ga21(1,2,isig2,isig),pentemp
     3    ,jbga21(2,isig2,isig),Div)

         enddo
       enddo   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Z A2 g
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

       call penlineABETotal(p1,mqz,mqa2,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),atau2(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentcontr5,j5a2g(1,1,isig2,isig),pentcontr5b
     3    ,jba2g(1,isig2,isig),Div)
 
        comp=-1

       call penlineABETotal(p1,mqz,mqa2,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau1(1,1),atau2(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentemp,j5a2g1(1,1,isig2,isig),pentemp
     3    ,jba2g1(1,isig2,isig),Div)

         enddo
       enddo   

C
C Z g A2
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

       call penlineABETotal(p1,mqz,mpg,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),Myeps(0,isig2),atau2(1,1),isig,scale,1,3,
     2    comp*4,comp,pentcontr6,j5ga2(1,1,isig2,isig),pentcontr6b
     3    ,jbga2(1,isig2,isig),Div)
 
       comp=-1

       call penlineABETotal(p1,mqz,mpg,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    atau1(1,1),Myeps(0,isig2),atau2(1,1),isig,scale,1,3,
     2    comp*4,comp,pentemp,j5ga21(1,1,isig2,isig),pentemp
     3    ,jbga21(1,isig2,isig),Div)

         enddo
       enddo   

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C  NO ABELIAN PART
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C A2 W
C
       comp=1
      do isig= -1,1,2 ! spinors
        do isig2 =1,2 ! Gluon

      call penlineNoAbeTotal(p1,mqa2,mqz,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),wp(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentcontrNoAbe1,jg5a2(2,isig2,isig),Div)

      comp=-1

      call penlineNoAbeTotal(p1,mqa2,mqz,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    atau2(1,1),atau1(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentemp,jg5a21(2,isig2,isig),Div)

         enddo
       enddo   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A2
C
       comp=1
      do isig= -1,1,2 ! spinors
        do isig2 =1,2 ! Gluon
      call penlineNoAbeTotal(p1,mqz,mqa2,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),atau2(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentcontrNoAbe2,jg5a2(1,isig2,isig),Div)

       comp=-1

      call penlineNoAbeTotal(p1,mqz,mqa2,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    atau1(1,1),atau2(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentemp,jg5a21(1,isig2,isig),Div)

         enddo
       enddo   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Gauge Test
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 g Z
      call pent_gauge(p1,mqa2,mpg,mqz,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,1),wp(1,1),-1,scale,
     2    pentcontr1,pentcontr1b,wardidtest1)
C g A2 Z
       call pent_gauge(p1,mpg,mqa2,mqz,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau2(1,1),wp(1,1),-1,scale,
     2   pentcontr2,pentcontr2b,wardidtest2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 Z g
       call pent_gauge(p1,mqa2,mqz,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wp(1,1),Myeps(0,1),-1,scale,
     2    pentcontr3,pentcontr3b,wardidtest3)
C g Z A2
      call pent_gauge(p1,mpg,mqz,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wp(1,1),atau2(1,1),-1,scale,
     2    pentcontr4,pentcontr4b,wardidtest4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Z A2 g
      call pent_gauge(p1,mqz,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),atau2(1,1),Myeps(0,1),-1,scale,
     2    pentcontr5,pentcontr5b,wardidtest5)
C Z g A2
      call pent_gauge(p1,mqz,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),Myeps(0,1),atau2(1,1),-1,scale,
     2    pentcontr6,pentcontr6b,wardidtest6)
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  NO ABELIAN PART
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 W
      call pentNoAbe_gauge(1d-1,p1,mqa2,mqz,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wp(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe1,accuracyres,wardidtest7)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A2
      call pentNoAbe_gauge(1d-1,p1,mqz,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),atau2(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe2,accuracyres,wardidtest8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!

       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6
       wardidtest1=wardidtest1.and.wardidtest7.and.wardidtest8


        If(Div.ne.0) then
! Check factorization of Divergences. Set the color equal.
! So that I can check separately the CF and CA factorization
!	Print*,'CFCFCF',cfc
	cfc=1d0
	cacgraph=1d0
!	Print*,'CFCFCF',cfc
	cfcac=cfc
	wardidtest1=.True.
	endif


         if (.not.wardidtest1) then
          do isig=-1,1,2  
	    do  isig2=1,2
                do i = 1,3
		   do mu=1,2

                     j5a2g(mu,i,isig2,isig) = (0.0d0, 0.0d0)
 		     j5ga2(mu,i,isig2,isig) = (0.0d0, 0.0d0)

                     j5a2g1(mu,i,isig2,isig) = (0.0d0, 0.0d0)
 		     j5ga21(mu,i,isig2,isig) = (0.0d0, 0.0d0)

 	          jg5a2(mu,isig2,isig)=(0.0d0, 0.0d0)
 	          jg5a21(mu,isig2,isig)=(0.0d0, 0.0d0)


		   enddo
                enddo
              enddo
           enddo 
              pentcountm = pentcountm+1
              pentcountps(ps_number) = pentcountps(ps_number)+1 

c           print*, "NO PEN GAUGE", ps_number,pentcountm,pentcountps(ps_number)
         else

          do isig=-1,1,2  
	    do isig2=1,2
              do i = 1,3
           j5a2g(1,i,isig2,isig)=(cfc*j5a2g(1,i,isig2,isig)
     &                         +cfcac*j5a2g(2,i,isig2,isig))
           j5ga2(1,i,isig2,isig)=(cfc*j5ga2(1,i,isig2,isig)
     &                          +cfcac*j5ga2(2,i,isig2,isig))


           j5a2g1(1,i,isig2,isig)=(cfc*j5a2g1(1,i,isig2,isig)
     &                          +cfcac*j5a2g1(2,i,isig2,isig))
           j5ga21(1,i,isig2,isig)=(cfc*j5ga21(1,i,isig2,isig)
     &                          +cfcac*j5ga21(2,i,isig2,isig))

               enddo
             do i = 1,2
        	 jg5a2(i,isig2,isig) = cacgraph*jg5a2(i,isig2,isig)
cccccc		 
        	 jg5a21(i,isig2,isig) = cacgraph*jg5a21(i,isig2,isig)
    	      enddo
	    enddo  
	  enddo  
              pentcount2m = pentcount2m+1
              pentcount2ps(ps_number) = pentcount2ps(ps_number)+1      
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   q qbar -> A A W
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do isig= -1,1,2
         do isig2=1,2
!  A2gZ  gA2Z  A2Zg  + 3 other
        mv2b=( 
     1   jba2g(1,isig2,isig)+ jba2g(2,isig2,isig)+jba2g(3,isig2,isig)
     2  +jbga2(1,isig2,isig)+jbga2(2,isig2,isig)+jbga2(3,isig2,isig))

!  A2gA* gA2A* A2A*g + 3 other
        mv3b=( 
     1   jba2g1(1,isig2,isig)+ jba2g1(2,isig2,isig)+jba2g1(3,isig2,isig)
     2  +jbga21(1,isig2,isig)+jbga21(2,isig2,isig)+jbga21(3,isig2,isig))

c      Print*, 'MV22B',jba2g(1,isig2,isig)+ jba2g(2,isig2,isig)+jba2g(3,isig2,isig)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

       matb(1,isig,isig2,1)=
     1     (mv2b)*clr(3,2,isig)*clr(3,1,isig)
     2    +(mv3b)*clr(3,1,isig)**2
    
       matb(2,isig,isig2,1)=
     1     (mv2b)*clr(4,2,isig)*clr(4,1,isig)
     2    +(mv3b)* clr(4,1,isig)**2

c Change of sign to agree with MadGraph

      matb(1,isig,isig2,1)=- matb(1,isig,isig2,1)
      matb(2,isig,isig2,1)=- matb(2,isig,isig2,1)
     
      ratioComp=abs(mat(1,isig,isig2,1)/matb(1,isig,isig2,1)-1d0)
      

      if(ratioComp.gt.1d-6) then

      print*
      print*
      print*, "ISIG",ISIG, "ISIG2",ISIG2
      print*
      print*
      PRINT*, 'Madgr/penline 1', ratioComp


      ratioComp=abs(mat(2,isig,isig2,1)/matb(2,isig,isig2,1)-1d0)

      PRINT*, 'Madgr/penline 2', ratioComp

      Print*, 'wardidtest1_PEN', wardidtest1
      Print*,"Born/b_pen-1d0",ratioComp
      Print*, "Born_uu",mat(1,isig,isig2,1)
      Print*, 'Sum_uu',matb(1,isig,isig2,1)
      Print*, 'ratio_uu',mat(1,isig,isig2,1)/matb(1,isig,isig2,1)

      Print*, "Born_dd",mat(2,isig,isig2,1)
      Print*, 'Sum_dd',matb(2,isig,isig2,1)
      PRINT*,  'RATION_DD', mat(2,isig,isig2,1)/matb(2,isig,isig2,1)

c      Print*, 'mv2b',mv2b
c      Print*, 'mv3b',mv3b

c      stop
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                       VIRTUAL
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Add here the Virtual part. Same as born, but replacing 
C  in the above eq. jb -> j5.
C  In addition the non-Abe terms have to be added.
C No way to check the part.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   q qbar -> A A W
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  A2gZ  gA2Z  A2Zg  + 3 other
        mv2=( 
     1   j5a2g(1,1,isig2,isig)+j5a2g(1,2,isig2,isig)+j5a2g(1,3,isig2,isig)
     2  +j5ga2(1,1,isig2,isig)+j5ga2(1,2,isig2,isig)+j5ga2(1,3,isig2,isig))

!  A2gA* gA2A* A2A*g + 3 other
        mv3=( 
     1   j5a2g1(1,1,isig2,isig)+j5a2g1(1,2,isig2,isig)+j5a2g1(1,3,isig2,isig)
     2  +j5ga21(1,1,isig2,isig)+j5ga21(1,2,isig2,isig)+j5ga21(1,3,isig2,isig))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non-Abelian graphs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ZA2 + A2Z
	mv4=(jg5a2(1,isig2,isig)+jg5a2(2,isig2,isig))
! A*A2 + A2A* 
	mv5=(jg5a21(1,isig2,isig)+jg5a21(2,isig2,isig))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	IF (wardidtest1 ) THEN
	
	IF (PerformGauge.eq.1) then
! gauge test for the gluon. The CF part of ,mv2,mv3,
! is already zero. mv1,mv2,mv3, then, only contains CA part which
! should cancel agains NOABE, mv4,mv5
	Print*, ''
	Print*, 'Div',Div
	Print*, 'PerformGauge',PerformGauge
	
        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	Print*, 'mv2b',mv2b
	Print*, 'mv2',mv2
	Print*, 'mv4',mv4
	Print*, 'Gauge',-mv2/mv4-1d0

	
        Print*, 'mv3',mv3
	Print*, 'mv5',mv5
	Print*, 'Gauge',-mv3/mv5-1d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif(PerformGauge.gt.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!   Here only the sum of all the contributions is gauge
!   invariant, but, with coupling cte equal=1. 
!    Otherwise, there are reamining terms that cancel against
!    penline's and so on due to the difference of clr(u(k),1,-1)
!    and clr(d(k),1,-1)
	Print*, ''
	Print*, 'Div',Div
	Print*, 'PerformGauge',PerformGauge

        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	
	Print*, ''
	Print*, 'mv2b',mv2b
        print*, 'j5a2g(1,1,isig2,isig)',j5a2g(1,1,isig2,isig)
	Print*, 'Gauge mv2',mv2/j5a2g(1,1,isig2,isig)
        print*, 'jg5a2(1,isig2,isig)',jg5a2(1,isig2,isig)
	Print*, 'GAuge mv4',mv4/jg5a2(1,isig2,isig)
c	Print*, 'Gauge',-mv2/mv4-1d0

        print*, 'j5a2g1(1,1,isig2,isig)',j5a2g1(1,1,isig2,isig)
        Print*, 'Gauge mv3',mv3/j5a2g1(1,1,isig2,isig)
        print*, 'jg5a21(1,isig2,isig)',jg5a21(1,isig2,isig)
	Print*, 'Gauge mv5',mv5/jg5a21(1,isig2,isig)
c	Print*, 'Gauge',-mv3/mv5-1d0


	endif
	endif ! Warditest


	If(wardidtest1) then
	
	If(Div.ne.0) then
!    Check Factorization of CF and CA part separately	
!    Mind that the born part is already checked against Madgraph.
!    So this test is quite good!!
!    At this point mv1,mv2,mv3, only contains the CF part.
!    The color factor of result(2) has been set to CF above.
!    I have to duplicate the virtual part now but with the 
!    CA part: Copy and paste of above expression. and replace
!    in the first argument 1->2       

!  A2gZ  gA2Z  A2Zg  + 3 other
        mv2CA=( 
     1   j5a2g(2,1,isig2,isig)+j5a2g(2,2,isig2,isig)+j5a2g(2,3,isig2,isig)
     2  +j5ga2(2,1,isig2,isig)+j5ga2(2,2,isig2,isig)+j5ga2(2,3,isig2,isig))

!  A2gA* gA2A* A2A*g + 3 other
        mv3CA=( 
     1   j5a2g1(2,1,isig2,isig)+j5a2g1(2,2,isig2,isig)+j5a2g1(2,3,isig2,isig)
     2  +j5ga21(2,1,isig2,isig)+j5ga21(2,2,isig2,isig)+j5ga21(2,3,isig2,isig))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	If(Div.eq.2) then
! CF part
	 Print*, ""
	 Print*, "CF PART"
	 Print*, ""	
	 ratioFact=mv2/(-2d0*mv2b)
	 Print*, 'mv2',mv2
	 Print*, 'mv2b',mv2b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 
	 ratioFact=mv3/(-2d0*mv3b)
	 Print*, 'mv3',mv3
	 Print*, 'mv3b',mv3b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 Print*, ""
	 Print*, 'SUM of born', mv2b+mv3b
! CA part
	 Print*, ""
	 Print*, ""
	 Print*, "CA PART"
	 Print*, ""
	 Print*, ""
	 Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv1b',mv2b
	 ratioFact=(mv2CA*(-1d0/2d0)+mv4)/(-1d0*mv2b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""
	 
         Print*, 'mv3CA',mv3CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv5
         Print*, 'mv3b',mv3b
	 ratioFact=(mv3CA*(-1d0/2d0)+mv5)/(-1d0*mv3b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, "" 
	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	else   !  Div 1/Eps
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
	lns=log(Abs(mans))-(0,1)*pi*theta(mans)-log(musq)
	lnt=log(Abs(mant))-(0,1)*pi*theta(mant)-log(musq)
	lnu=log(Abs(manu))-(0,1)*pi*theta(manu)-log(musq)

! CF PART
	 cteEps1=(-3d0+2d0*lns)
	 Print*, ""
	 Print*, "CF PART"
	 Print*, ""
	 Print*, 'CF cte',cteEps1
	 Print*, ""
	 
	 ratioFact=mv2/(mv2b*cteEps1)
	 Print*, 'mv2',mv2
	 Print*, 'mv2b*cte',mv2b*cteEps1
	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""

	 ratioFact=mv3/(mv3b*cteEps1)
	 Print*, 'mv3',mv3
	 Print*, 'mv3b*cte',mv3b*cteEps1
	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""
	 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CA part
	 cteEps1=(-lns+lnt+lnu)
	 
	 Print*, "CA PART"
	 Print*, ""
	 Print*, 'cteEps1',cteEps1
	 Print*, ""
	 Print*, 'mv1CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv2b*cte',mv2b*cteEps1
	 ratioFact=(mv2CA*(-1d0/2d0)+mv4)/(mv2b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	 Print*, 'mv3CA',mv3CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv5
         Print*, 'mv3b*cte',mv3b*cteEps1
	 ratioFact=(mv3CA*(-1d0/2d0)+mv5)/(mv3b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	endif    !Div
	
	endif    ! Div, Compute CA part of ABe contributions
	endif   ! Warditest set to true for Div checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   SUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

       matv(1,isig,isig2,1)=
     1     (mv2+mv4)*clr(3,2,isig)*clr(3,1,isig)
     2    +(mv3+mv5)*clr(3,1,isig)**2
    
       matv(2,isig,isig2,1)=
     1     (mv2+mv4)*clr(4,2,isig)*clr(4,1,isig)
     2    +(mv3+mv5)* clr(4,1,isig)**2

c Change of sign to agree with MadGraph

      matv(1,isig,isig2,1)=- matv(1,isig,isig2,1)
      matv(2,isig,isig2,1)=- matv(2,isig,isig2,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PerformGauge.eq.1) then
! then, matv(k,-1,l,1) should be zeror
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(1,isg,isg2,1)', matv(1,isig,isig2,1)
      Print*,' matv(1,isg,isg2,1)', matv(2,isig,isig2,1)
      Print*, 'mv1*clr(3,3,-1)',mv2*clr(3,2,isig)*clr(3,1,isig)

      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(k,-1,l,1)', matv(1,isig,isig2,1)
     &                /(jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig)  )
      PRINT*, 'jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig) ',jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig) 
c      cteEps1=1!clr(u(k),1,-1)/clr(d(k),1,-1) 
c      matv(k,-1,l,2) = (mv2+mv4+(mv3+mv5)*cteEps1)!*clr(3,3,-1)
c      Print*,'matv(k,-1,l,2) with corrected couplings'
c      Print*,'matv(k,-1,l,2)', matv(k,-1,l,2)
      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,1)',matv(1,isig,isig2,1)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(1,isig,isig2,1)
      Print*, 'ratio',matv(1,isig,isig2,1)/(-4*matb(1,isig,isig2,1))
      print*
      Print*,'matv(k,-1,l,1)',matv(2,isig,isig2,1)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(2,isig,isig2,1)
      Print*, 'ratio',matv(2,isig,isig2,1)/(-4*matb(2,isig,isig2,1))


      endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!              Add and substract cte's
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      	  matv(1,isig,isig2,1)=matv(1,isig,isig2,1)+matb(1,isig,isig2,1)*convfact
	  matv(1,isig,isig2,1)=matv(1,isig,isig2,1)-matb(1,isig,isig2,1)*(K_cont)

     	  matv(2,isig,isig2,1)=matv(2,isig,isig2,1)+matb(2,isig,isig2,1)*convfact
	  matv(2,isig,isig2,1)=matv(2,isig,isig2,1)-matb(2,isig,isig2,1)*(K_cont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	enddo ! end isig2 loop
	enddo ! end isig loop


!  Check renormalization invariance using K_contcheck      
      if(Reno_inv.eq.1) then
       do k=1,2
         do isig=-1,1,2
           do isig2=1,2
      check(k,isig,isig2)=matv(k,isig,isig2,1)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,1)

      check1(k,isig,isig2)=matv(k,isig,isig2,1)
            enddo
          enddo
       enddo   
      Reno_inv=2
      musq=1d0
      goto 100
      endif

      If (Reno_inv.gt.0) then
       do k=1,2
         do isig=-1,1,2
           do isig2=1,2
      print*
      print*
      print*, "K",k,'isig',isig,'isig2',isig2
      print*
       Print*, ' '
       Print*, 'chek1 hel 1 ',check1(k,isig,isig2)
       Print*, 'chek1 hel 1 ',check(k,isig,isig2)

       Print*, ''
       Print*, 'check1 hel 1 ',matv(k,isig,isig2,1)
       Print*, 'check1 hel 1 ',(matv(k,isig,isig2,1)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,1))
      
       Print*, 'ratio1',check(k,isig,isig2)/((matv(k,isig,isig2,1)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,1)))-1d0
       

!       stop
         enddo
         enddo
       enddo

      endif	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c FC
c      Print*, "Here"
c      stop
c FC
      endif ! end lpent  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Finish Pentagons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Start BOXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               qq -> Z g ->Z A g -----
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(lbox) then

cFC
c  Reno_inv=1 Checked!!! CA,CF color structure checked
c  Div=1,2  Checked!! Factorization of divergences,CA,CF color checked
c  Performed_gauge=1,2 ! Checked
c
c
cFC


      scale=musq!mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
            
            do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
	       mpg(mu) = -q(mu)
               mqza(mu) = -p1(mu)-p2(mu)+q(mu)
c for gauge test	       
	       emqza(mu)=mqza(mu)
	       empg(mu)=mpg(mu)
	   enddo
!!! Check Gauge invariance
      If(PerformGauge.eq.1) then
         do mu = 0,3
	    Myeps(mu,1)=mpg(mu)
	    Myeps(mu,2)=mpg(mu)
	  enddo
      elseif(PerformGauge.eq.2) then
         do mu = 0,3
	     mutemp=mu+1
	    wz(mutemp,1) = mqza(mu)
      	    wa(mutemp,1) = mqza(mu)
	 enddo
       endif    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Z g 
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wz(1,1),Myeps(0,isig2),isig,scale,1,3,
     2    7*comp,comp,boxcontr1,j5g(1,1,isig2,isig),boxcontr1b
     3    ,jbg(1,isig2,isig),Div)
       comp=-1
c KK Terms and gauge test
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    emqza,Myeps(0,isig2),isig,scale,1,3,
     2    3*comp,comp,nothing,j5gg(1,1,isig2,isig),nothing
     3    ,jbgg(1,isig2,isig),Div)
 
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wz(1,1),empg,isig,scale,1,3,
     2    3*comp,comp,nothing,j5gg2(1,1,isig2,isig),nothing
     3    ,jbgg2(1,isig2,isig),Div)

       if(isig*isig2.eq.-1) then
       boxcontr1(1)=j5gg(1,1,isig2,isig)+j5gg(2,1,isig2,isig)
       boxcontr1b(1)=jbgg(1,isig2,isig)

       boxcontr1(2)=j5gg2(1,1,isig2,isig)+j5gg2(2,1,isig2,isig)
       boxcontr1b(2)=jbgg2(1,isig2,isig)
       endif

c  A*
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wa(1,1),Myeps(0,isig2),isig,scale,1,3,
     2    3*comp,comp,nothing,j5g1(1,1,isig2,isig),nothing
     3    ,jbg1(1,isig2,isig),Div)
c KK is the same as for Z, at the end they should cancel.
           enddo
       enddo    
C
C g Z
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2

 
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),wz(1,1),isig,scale,1,2,
     2    7*comp,comp,boxcontr2,j5g(1,2,isig2,isig),boxcontr2b
     3    ,jbg(2,isig2,isig),Div)
       comp=-1
c KK Terms and gauge test
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),emqza,isig,scale,1,2,
     2    3*comp,comp,nothing,j5gg(1,2,isig2,isig),nothing
     3    ,jbgg(2,isig2,isig),Div)
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    empg,wz(1,1),isig,scale,1,2,
     2    3*comp,comp,nothing,j5gg2(1,2,isig2,isig),nothing
     3    ,jbgg2(2,isig2,isig),Div)

       if(isig*isig2.eq.-1) then
       boxcontr2(1)=j5gg2(1,2,isig2,isig)+j5gg2(2,2,isig2,isig)
       boxcontr2b(1)=jbgg2(2,isig2,isig)

       boxcontr2(2)=j5gg(1,2,isig2,isig)+j5gg(2,2,isig2,isig)
       boxcontr2b(2)=jbgg(2,isig2,isig)
       endif

c  A*
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),wa(1,1),isig,scale,1,2,
     2    3*comp,comp,nothing,j5g1(1,2,isig2,isig),nothing
     3    ,jbg1(2,isig2,isig),Div)
c KK is the same as for Z, at the end they should cancel.
           enddo
       enddo    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Z
C
       comp=1
      do isig= -1,1,2
        do isig2 =1,2
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wz(1,1),Myeps(0,isig2),isig,scale,20*comp,comp,
     2    boxcontrNoAbe1,jg5n(isig2,isig),Div)
      comp=-1
c For non-abelian gauge and kk/M terms          
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    emqza,Myeps(0,isig2),isig,scale,2*comp,comp,
     2    nothing,jg5ng(isig2,isig),Div)
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wz(1,1),empg,isig,scale,2*comp,comp,
     2    nothing,jg5ng2(isig2,isig),Div)

       if(isig*isig2.eq.-1) then
c          print*,boxcontrNoAbe1(1)
c          print*,boxcontrNoAbe1(2)
        boxcontrNoAbe1(1)=jg5ng2(isig2,isig)
        boxcontrNoAbe1(2)=jg5ng(isig2,isig)
c         print*,boxcontrNoAbe1(1)
c         print*,boxcontrNoAbe1(2)
       endif
c AA 
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wa(1,1),Myeps(0,isig2),isig,scale,2*comp,comp,
     2    nothing,jg5n1(isig2,isig),Div)
        enddo
      enddo  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Gauge Boxes Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C W g 
      call box_gauge(p1,mqza,mpg,p2,scale,
     2    boxcontr1,boxcontr1b,wardidtest1)
C g W
      call box_gauge(p1,mpg,mqza,p2,scale,
     2    boxcontr2,boxcontr2b,wardidtest2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
      wardidtest3=.false.
      if (abs(jg5ng(1,-1)).lt.1d-3) then
         td=abs(-(j5gg(2,2,1,-1)+j5gg(2,1,1,-1))/2d0-jg5ng(1,-1))
      else
         td=abs(-(j5gg(2,2,1,-1)+j5gg(2,1,1,-1))/2d0/jg5ng(1,-1)+1d0)
      endif
      IF(td.le.1.5d-1) then
           if (abs(jg5ng2(1,-1)).lt.1d-3) then
            td=abs(-(j5gg2(2,2,1,-1)+j5gg2(2,1,1,-1))/2d0-jg5ng2(1,-1))
           else
             td=abs(-(j5gg2(2,2,1,-1)+j5gg2(2,1,1,-1))/2d0/jg5ng2(1,-1)+1d0)
           endif
         IF(td.le.1.5d-1) wardidtest3=.true.
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
cFC c	Print*, 'wardidtest1',wardidtest1
cFC c	Print*, 'wardidtest2',wardidtest2
cFC c	Print*, 'wardidtest3',wardidtest3
c      STOP
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
        If(Div.ne.0) then
! Check factorization of Divergences. Set the color equal.
! So that I can check separately the CF and CA factorization
!	Print*,'CFCFCF',cfc
	cfc=1d0
	cacgraph=1d0
!	Print*,'CFCFCF',cfc
	cfcac=cfc
	wardidtest1=.True.
	endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     

	    if (.not.wardidtest1) then
          do isig=-1,1,2  
	    do  isig2=1,2
                do i = 1,2
		  do mu=1,2
                   j5g(mu,i,isig2,isig) = (0.0d0, 0.0d0)
		   j5gg(mu,i,isig2,isig) = (0.0d0, 0.0d0)
                   j5gg2(mu,i,isig2,isig) = (0.0d0, 0.0d0)
		   j5g1(mu,i,isig2,isig) = (0.0d0, 0.0d0)

             	
 	           jg5n(isig2,isig) =     (0.0d0, 0.0d0)
		   jg5ng(isig2,isig)=     (0.0d0, 0.0d0)

 	           jg5ng2(isig2,isig) =     (0.0d0, 0.0d0)
		   jg5n1(isig2,isig)=     (0.0d0, 0.0d0)


                   enddo
                  enddo
	        enddo
	       enddo	

 	      boxcountm = boxcountm+1
              boxcountps(ps_number) = boxcountps(ps_number)+1 
c               print*, "NO BOX GAUGE",ps_number,boxcountm,boxcountps(ps_number) 
             else
         do isig=-1,1,2  
	    do isig2=1,2
                  do i = 1,2

		     j5g(1,i,isig2,isig)=(cfc*j5g(1,i,isig2,isig)+cfcac*j5g(2,i,isig2,isig))
		     j5gg(1,i,isig2,isig)=(cfc*j5gg(1,i,isig2,isig)+cfcac*j5gg(2,i,isig2,isig))
		     j5gg2(1,i,isig2,isig)=(cfc*j5gg2(1,i,isig2,isig)+cfcac*j5gg2(2,i,isig2,isig))
		     j5g1(1,i,isig2,isig)=(cfc*j5g1(1,i,isig2,isig)+cfcac*j5g1(2,i,isig2,isig))

                  enddo
	             jg5n(isig2,isig)=cacgraph*jg5n(isig2,isig)
		     jg5ng(isig2,isig)=cacgraph*jg5ng(isig2,isig)
	             jg5ng2(isig2,isig)=cacgraph*jg5ng2(isig2,isig)
		     jg5n1(isig2,isig)=cacgraph*jg5n1(isig2,isig)

           enddo          
	      enddo 
              boxcount2m = boxcount2m+1
              boxcount2ps(ps_number) = boxcount2ps(ps_number)+1      
             endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c   q qbar -> W^* g-> WAAg
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      do isig= -1,1,2
         do isig2=1,2
        mv1b =(jbg(2,isig2,isig)+  !g Z*
     2         jbg(1,isig2,isig))   !Z* g
c
        mv2b =(jbg1(2,isig2,isig)+  !g A*
     2         jbg1(1,isig2,isig))   !A* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1b =((jbgg(2,isig2,isig))*dotrc(mqza,wz(1,1))+  !g W^*
     2        (jbgg(1,isig2,isig))*dotrc(mqza,wz(1,1)))*zm2i(2) !W^* g

c Check here the zm2i
       mmv2b =0d0 !((jbgg(2,isig2,isig))*dotrc(mqza,wa(1,1))+  !g W^*
!     2       (jbgg(1,isig2,isig))*dotrc(mqza,wa(1,1)))*zm2i(3) !W^* g

C Change of sign to agree with Madgraph
        mv1b =-mv1b
	mmv1b =0!-mmv1b
        mv2b =-mv2b
	mmv2b =0!-mmv2b

       matb(1,isig,isig2,2)=(mv1b+mmv1b)*clr(3,2,isig)
     1                      +(mv2b+mmv2b)*clr(3,1,isig)
       matb(2,isig,isig2,2)=(mv1b+mmv1b)*clr(4,2,isig)
     1                     +(mv2b+mmv2b)*clr(4,1,isig)

      
       ratioComp=abs(mat(1,isig,isig2,2)/matb(1,isig,isig2,2)  -1d0)
	
      if(ratioComp.gt.1d-6) then

      print*
      print*
      print*, "ISIG",ISIG, "ISIG2",ISIG2
      print*
      print*
 
      Print*,"1: Born/b_box-1d0",ratioComp
      ratioComp=abs(mat(2,isig,isig2,2)/matb(2,isig,isig2,2)  -1d0)
      Print*,"2: Born/b_box-1d0",ratioComp
      Print*, 'wardidtest1_BOX', wardidtest1

      Print*, "Born_uu",mat(1,isig,isig2,2)
      Print*, "Sum_uu",matb(1,isig,isig2,2)
      print*
      Print*, "Born_dd",mat(2,isig,isig2,2)
      Print*, "Sum_dd",matb(2,isig,isig2,2)

      Print*, 'mv1',mv1*clr(3,3,-1)
      Print*, 'mmv1',mmv1
c      stop
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        VIRTUAL
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Add here the Virtual part. Same as born, but replacing 
C  in the above eq. jb -> j5.
C  In addition the non-Abe terms have to be added.
C No way to check the part.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        mv1 =(j5g(1,2,isig2,isig)+  !g Z*
     2        j5g(1,1,isig2,isig))   !Z* g
c
        mv2 =(j5g1(1,2,isig2,isig)+  !g A*
     2         j5g1(1,1,isig2,isig))   !A* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

       mmv1 =((j5gg(1,2,isig2,isig))*dotrc(mqza,wz(1,1))+  !g W^*
     2        (j5gg(1,1,isig2,isig))*dotrc(mqza,wz(1,1)))*zm2i(2) !W^* g

c Check here the zm2i
       mmv2 =0d0 !((jbgg(2,isig2,isig))*dotrc(mqza,wa(1,1))+  !g W^*
!     2       (jbgg(1,isig2,isig))*dotrc(mqza,wa(1,1)))*zm2i(3) !W^* g


       mmv1 =0d0
       mmv2 =0d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
C  No abelian part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
	mv3=jg5n(isig2,isig)
        mv4=jg5n1(isig2,isig)

	mmv3=0d0
        mmv4=0d0

C Change of sign to agree with Madgraph
        mv1=-mv1
        mv2=-mv2
        mv3 =-mv3
        mv4 =-mv4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	IF (wardidtest1 ) THEN
	
	IF (PerformGauge.eq.1) then
! gauge test for the gluon. The CF part of mv1,mv2,mv3,
! is already zero. mv1,mv2,mv3, then, only contains CA part which
! should cancel agains NOABE, mv4,mv5,mv6
	Print*, ''
	Print*, 'Div',Div
	Print*, 'PerformGauge',PerformGauge
	
        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	Print*, ''
	Print*, 'mv1b',mv1b
	Print*, 'mv1',mv1
	Print*, 'mv3',mv3
	Print*, 'Gauge',-mv1/mv3-1d0

	
        Print*, 'mv2',mv2
	Print*, 'mv4',mv4
	Print*, 'Gauge',-mv2/mv4-1d0
	
!         stop
        elseif(PerformGauge.gt.1) then
!   Here only the sum of all the contributions is gauge
!   invariant, but, with coupling cte equal=1. 
!    Otherwise, there are reamining terms that cancel against
!    penline's and so on due to the difference of clr(u(k),1,-1)
!    and clr(d(k),1,-1)

        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	Print*, ''
	Print*, 'mv1b',mv1b
        print*, 'j5g(1,1,isig2,isig)',j5g(1,1,isig2,isig)
        print*, 'j5g(1,1,isig2,isig)',j5g(1,2,isig2,isig)
	Print*, 'Gauge mv1',mv1!/j5g(1,1,isig2,isig)

        print*, 'jg5a2(1,isig2,isig)',jg5n(isig2,isig)
	Print*, 'GAuge mv3',mv3
c	Print*, 'Gauge',-mv2/mv4-1d0

        print*, 'j5a2g1(1,1,isig2,isig)',j5g1(1,1,isig2,isig)
        Print*, 'Gauge mv2',mv2!/j5a2g1(1,1,isig2,isig)

        print*, 'jg5a21(1,isig2,isig)',jg5n1(isig2,isig)
	Print*, 'Gauge mv5',mv4!/jg5a21(1,isig2,isig)
c	Print*, 'Gauge',-mv3/mv5-1d0
	
!        stop
	endif
	endif ! Warditest


	If(wardidtest1) then
	
	If(Div.ne.0) then
!    Check Factorization of CF and CA part separately	
!    Mind that the born part is already checked against Madgraph.
!    So this test is quite good!!
!    At this point mv1,mv2,mv3, only contains the CF part.
!    The color factor of result(2) has been set to CF above.
!    I have to duplicate the virtual part now but with the 
!    CA part: Copy and paste of above expression. and replace
!    in the first argument 1->2       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        mv1CA =(j5g(2,2,isig2,isig)+  !g Z*
     2        j5g(2,1,isig2,isig))   !Z* g
c
        mv2CA =(j5g1(2,2,isig2,isig)+  !g A*
     2         j5g1(2,1,isig2,isig))   !A* g

        mv1CA=-mv1CA
        mv2CA=-mv2CA

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1CA =0d0
       mmv2CA =0d0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	If(Div.eq.2) then
! CF part
	 Print*, ""
	 Print*, "CF PART"
	 Print*, ""	
	 ratioFact=mv1/(-2d0*mv1b)
	 Print*, 'mv1',mv1
	 Print*, 'mv1b',mv1b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 
	 ratioFact=mv2/(-2d0*mv2b)
	 Print*, 'mv2',mv2
	 Print*, 'mv2b',mv2b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 	 
	 Print*, ""
	 !Print*, 'SUM of born', mmv1b+mmv2b	 
	 
! CA part
	 Print*, ""
	 Print*, ""
	 Print*, "CA PART"
	 Print*, ""
	 Print*, ""
	 Print*, 'mv1CA',mv1CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv3
         Print*, 'mv1b',mv1b
	 ratioFact=(mv1CA*(-1d0/2d0)+mv3)/(-1d0*mv1b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""

 	 Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv2b',mv2b
	 ratioFact=(mv2CA*(-1d0/2d0)+mv4)/(-1d0*mv2b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""

	 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	else   !  Div 1/Eps
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
	lns=log(Abs(mans))-(0,1)*pi*theta(mans)-log(musq)
	lnt=log(Abs(mant))-(0,1)*pi*theta(mant)-log(musq)
	lnu=log(Abs(manu))-(0,1)*pi*theta(manu)-log(musq)

! CF PART
	 cteEps1=(-3d0+2d0*lns)
	 Print*, ""
	 Print*, "CF PART"
	 Print*, ""
	 Print*, 'CF cte',cteEps1
	 Print*, ""
	 
	 ratioFact=mv1/(mv1b*cteEps1)
	 Print*, 'mv1',mv1
	 Print*, 'mv1b*cte',mv1b*cteEps1
	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""

	 ratioFact=mv2/(mv2b*cteEps1)
	 Print*, 'mv2',mv2
	 Print*, 'mv2b*cte',mv2b*cteEps1
	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CA part
	 cteEps1=(-lns+lnt+lnu)
	 
	 Print*, "CA PART"
	 Print*, ""
	 Print*, 'cteEps1',cteEps1
	 Print*, ""
	 Print*, 'mv1CA',mv1CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv3
         Print*, 'mv1b*cte',mv1b*cteEps1
	 ratioFact=(mv1CA*(-1d0/2d0)+mv3)/(mv1b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""
	 Print*, ""
	 Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv2b*cte',mv2b*cteEps1
	 ratioFact=(mv2CA*(-1d0/2d0)+mv4)/(mv2b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	endif    !Div
	
	endif    ! Div, Compute CA part of ABe contributions
	endif   ! Warditest set to true for Div checks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   SUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       matv(1,isig,isig2,2)=(mv1+mv3)*clr(3,2,isig)
     1                     +(mv2+mv4)*clr(3,1,isig)
       matv(2,isig,isig2,2)=(mv1+mv3)*clr(4,2,isig)
     1                     +(mv2+mv4)*clr(4,1,isig)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (PerformGauge.ge.1) then
! then, matv(k,-1,l,1) should be zero
 
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(1,isg,isg2,1)', matv(1,isig,isig2,2)
      Print*,' matv(1,isg,isg2,1)', matv(2,isig,isig2,2)
      Print*, 'mv1*clr(3,3,-1)',mv2*clr(3,1,isig)
      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(k,-1,l,1)', matv(1,isig,isig2,2)
      Print*,'j5', (jg5n(isig2,isig)*clr(3,2,isig))
      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,2)',matv(1,isig,isig2,2)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(1,isig,isig2,2)
      Print*, 'ratio',matv(1,isig,isig2,2)/(-4*matb(1,isig,isig2,2))
      print*
      Print*,'matv(k,-1,l,1)',matv(2,isig,isig2,2)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(2,isig,isig2,2)
      Print*, 'ratio',matv(2,isig,isig2,2)/(-4*matb(2,isig,isig2,2))
      endif  

c      stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Add and subtract cte's 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      matv(1,isig,isig2,2)=matv(1,isig,isig2,2)+matb(1,isig,isig2,2)*convfact
      matv(1,isig,isig2,2)=matv(1,isig,isig2,2)-matb(1,isig,isig2,2)*(K_cont)

      matv(2,isig,isig2,2)=matv(2,isig,isig2,2)+matb(2,isig,isig2,2)*convfact
      matv(2,isig,isig2,2)=matv(2,isig,isig2,2)-matb(2,isig,isig2,2)*(K_cont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo  !  isig2
        enddo ! isig
        
!  Check renormalization invariance using K_contcheck      
      if(Reno_inv.eq.1) then
       do k=1,2
         do isig=-1,1,2
           do isig2=1,2
      check(k,isig,isig2)=matv(k,isig,isig2,2)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,2)

      check1(k,isig,isig2)=matv(k,isig,isig2,2)
            enddo
          enddo
       enddo   
      Reno_inv=2
      musq=1d0
      goto 100
      endif

      If (Reno_inv.gt.0) then
       do k=1,2
         do isig=-1,1,2
           do isig2=1,2
      print*
      print*
      print*, "K",k,'isig',isig,'isig2',isig2
      print*
       Print*, ' '
       Print*, 'chek1 hel 1 ',check1(k,isig,isig2)
       Print*, 'chek1 hel 1 ',check(k,isig,isig2)

       Print*, ''
       Print*, 'check1 hel 1 ',matv(k,isig,isig2,2)
       Print*, 'check1 hel 1 ',(matv(k,isig,isig2,2)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,2))
      
       Print*, 'ratio1',check(k,isig,isig2)/((matv(k,isig,isig2,2)+(K_contcheck+K_contcheck4)*matb(k,isig,isig2,2)))-1d0
       

!       stop
         enddo
         enddo
       enddo

      endif

c FC
c      Print*, "End BOX Here"
c      stop
c FC     
	
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif ! end lbox
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C End Bosonic Virtual Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (lbox2) then
C Begin Fermionic Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c * W always attached to the fermion line (charge conservation)
c * gluon always attached to the loop (colour trace)
c * ggA-triangles are zero (Furry's theorem)
c * ggAA-boxes remaining
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       scale=musq  !mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
       Div=0

      mf(1) = 0d0 ! light up quarks
      mf(2) = 0d0 ! light down quarks
      mf(3) = xmt ! top
      mf(4) = bmass_h ! bottom

      vv = 2.d0*sqrt(xm2(3)*(1-xm2(3)/xm2(2)))/(-clr(2,1,-1)) ! 2*MW*SW/e

      cpl(1) = 3 ! light up quarks
      cpl(2) = 4 ! light down quarks
      cpl(3) = 3 ! top
      cpl(4) = 4 ! bottom
     
C transfer to M.Rauch notation of gg->ZA
c p(0,2) amd pbar(0,1) are change from his original code
     
            do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
c propapator momenta
               p12(mu)= p1(mu)+p2(mu)
c On-shell gluon
           mpg(mu) = -q(mu)
c  Momenta photon 
        minus_qa2z(mu)=  -qm(mu,1)
c Momenta  Z
        minus_qa1z(mu)= -qp(mu,1)
c Photon +Z
        minus_qaa(mu)=minus_qa2z(mu)+minus_qa1z(mu)
	   enddo

cFCc Propagator of virtual gluon:
cFC          s12=2d0*dotrr(p1,p2)+dotrr(p1,p1)+dotrr(p2,p2)
cFCc prepare gluon current
cFC      do isig=-1,1,2
cFC         if(isig.eq.-1) then
cFC            poln=1
cFC         else
cFC            poln=2
cFC          endif
cFC
cFC        do i=0,3
cFC          v1(0)=deltakro(i,0)
cFC          v1(1)=deltakro(i,1)
cFC          v1(2)=deltakro(i,2)
cFC          v1(3)=deltakro(i,3)
cFC          curr1(i,poln)=
cFC     2             (1d0)/s12*   ! From the propagator
cFC     1            sc1c(psi(1,isig,2),v1,psi(1,isig,1),isig)
cFC         enddo
cFC      enddo
cFC
    
      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),curr2(0,-1))
c      print*,"p12",(p12(i), i=0,3)
       call propagate(1,curr2(0,-1),p12c)
c        print*,"p12q",(p12c(i,-1), i=0,3)
c        print*,"p12q",(p12c(i,1), i=0,3)
! Checked against curr1 above: Now commented out

      do isig=-1,1,2
         if(isig.eq.-1) then
            poln=1
         else
            poln=2
          endif
        do i=0,3
          curr1(i,poln)=curr2(i,isig)
       enddo
      enddo 
          
cFC      print*, "curr1   1", (curr1(i,1),i=0,3)
cFC      print*, "curr2  -1", (curr2(i,-1),i=0,3)
cFC      print*,
cFC      print*, "curr1 2", (curr1(i,2),i=0,3)
cFC      print*, "curr2 1", (curr2(i,1),i=0,3)
cFC      print*,"Conser",(p12(i)+mpg(i)+minus_qa1z(i)+minus_qa2z(i) ,i=0,3)
cFC
cFC
cFC      stop
c Kinematics
cccccccccccccccccccccccccccccccccccccc
      if (N.ne.1) write(*,*) "Error in qqZZj.F: wrong kinematics: N = ", N, " =/= 1"
ccccccccccccccccccccccccccccccccc
c box loops
ccc first anchor point is always first gluon
ccc permutate second g and the two Z/photons over the other ones
ccc strong coupling is added later
cc g-1-2
       do j=1,4
          if(j.ne.2) comp=1
         do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),2,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),clr(cpl(j),1,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,mpg,minus_qa1z(0),minus_qa2z(0),
     &     curr1(0,l1),MyEps(0,l2),wp(1,N),atau2(1,N),
     &     scale,comp,comp,rfbox(1,1,j,l1,l2),rgfbox(1,1,j,l1,l2),
     &     Div)
         comp=-1
         
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,mpg,minus_qa1z(0),minus_qa2z(0),
     &     curr1(0,l1),MyEps(0,l2),atau1(1,N),atau2(1,N),
     &     scale,comp,comp,rfbox(1,2,j,l1,l2),rgfbox(1,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc g-2-1
       do j=1,4
          if(j.ne.2) comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),2,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,mpg,minus_qa2z(0),minus_qa1z(0),
     &     curr1(0,l1),MyEps(0,l2),atau2(1,N),wp(1,N),
     &     scale,comp,comp,rfbox(2,1,j,l1,l2),rgfbox(2,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,mpg,minus_qa2z(0),minus_qa1z(0),
     &     curr1(0,l1),MyEps(0,l2),atau2(1,N),atau1(1,N),
     &     scale,comp,comp,rfbox(2,2,j,l1,l2),rgfbox(2,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-g-2
       do j=1,4
          if(j.ne.2) comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),2,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),2,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j), mf(j), mf(j), mf(j),
     &     p12,minus_qa1z(0),mpg,minus_qa2z(0),
     &     curr1(0,l1),wp(1,N),MyEps(0,l2),atau2(1,N),
     &     scale,comp,comp,rfbox(3,1,j,l1,l2),rgfbox(3,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j), mf(j), mf(j), mf(j),
     &     p12,minus_qa1z(0),mpg,minus_qa2z(0),
     &     curr1(0,l1),atau1(1,N),MyEps(0,l2),atau2(1,N),
     &     scale,comp,comp,rfbox(3,2,j,l1,l2),rgfbox(3,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 2-g-1
       do j=1,4
       if(j.ne.2) comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),2,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),2,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa2z(0),mpg,minus_qa1z(0),
     &     curr1(0,l1),atau2(1,N),MyEps(0,l2),wp(1,N),
     &     scale,comp,comp,rfbox(4,1,j,l1,l2),rgfbox(4,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),1d0,clr(cpl(j),1,+1),
     &     1d0,clr(cpl(j),1,-1),1d0,clr(cpl(j),1,-1),
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa2z(0),mpg,minus_qa1z(0),
     &     curr1(0,l1),atau2(1,N),MyEps(0,l2),atau1(1,N),
     &     scale,comp,comp,rfbox(4,2,j,l1,l2),rgfbox(4,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo
cc 1-2-g
       do j=1,4
       if(j.ne.2) comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),2,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),2,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa1z(0),minus_qa2z(0),mpg,
     &     curr1(0,l1),wp(1,N),atau2(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(5,1,j,l1,l2),rgfbox(5,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa1z(0),minus_qa2z(0),mpg,
     &     curr1(0,l1),atau1(1,N),atau2(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(5,2,j,l1,l2),rgfbox(5,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo

cc 2-1-g
       do j=1,4
       if(j.ne.2) comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),2,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),2,-1),1d0,
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa2z(0),minus_qa1z(0),mpg,
     &     curr1(0,l1),atau2(1,N),wp(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(6,1,j,l1,l2),rgfbox(6,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),1,+1),clr(cpl(j),1,+1),1d0,
     &     1d0,clr(cpl(j),1,-1),clr(cpl(j),1,-1),1d0,
     &     mf(j),mf(j),mf(j),mf(j),
     &     p12,minus_qa2z(0),minus_qa1z(0),mpg,
     &     curr1(0,l1),atau2(1,N),atau1(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(6,2,j,l1,l2),rgfbox(6,2,j,l1,l2),
     &     Div)
       enddo
       enddo
       enddo

       wardidtest1 = .true.
       do i=1,6
           if (rgfbox(i,1,1,1,1).gt.1d-1) then 
c            if (wardidtest1) then
c           write(*,*) "MARKER1: gauge test failed", i,rgfbox(i,1,1,1,1)
c           write(*,*) "MARKER2: gauge test failed", j,i,rgfbox(i,2,1,1,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER4: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            else
c              write(*,*) "MARKER2: gauge test failed", j,i,rgfbox(i,j,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER5: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            endif
             wardidtest1=.false.
           endif
       enddo

c       write(*,*) "MARKER1: ", wardidtest1
c add everything together
       if (wardidtest1) then
             box2count2m = box2count2m+1
             box2count2ps(ps_number) = box2count2ps(ps_number)+1   
ccc         write(*,*) "MARKER3: gauge test ok"
         do l1=1,2
         do l2=1,2
           m1 = 0d0
           m2 = 0d0
           do i=1,6
           do j=1,2
             m1 = m1 + 2*rfbox(i,j,1,l1,l2) ! u,c
     &               + 1*rfbox(i,j,3,l1,l2) ! t
             m2 = m2 + 2*rfbox(i,j,2,l1,l2) ! d,s
     &               + 1*rfbox(i,j,4,l1,l2) ! b
           enddo
           enddo
        if(l1.eq.1) then
            poln=-1
         else
            poln=1
          endif
c and multiply with color factor TR

          matv(1,poln,l2,3) =0d0
          matv(2,poln,l2,3) =0d0

          matv(1,poln,l2,3) = (m1 + m2) * TR 
          matv(2,poln,l2,3) = (m1 + m2) * TR 

c          print*,"matv(1,poln,l2,3)",matv(1,poln,l2,3)

         enddo
         enddo
c         stop
        else !wardidtest
 	      box2countm = box2countm+1
              box2countps(ps_number) = box2countps(ps_number)+1 
c              print*, "NO GUAGE BOX2",ps_number, box2countm,box2countps(ps_number)
       endif

       if (lhiggs) then
c vertex diagrams
ccc strong coupling is added later
ccc triangles with Z/photon vanish, for Higgs only top relevant
ccc Both fermion-line directions identical - factor 2 added later
         do j=3,3 ! no contribution from ~ massless loops
         do l1=1,2 
         l2=l1 ! 2-1 and 1-2 vanish
           call gluonverHiggsEvenT(
     &       mf(j),
     &       p12,mpg,minus_qaa,
     &       curr1(0,l1),MyEps(0,l2),
     &       scale,1,rfvert(1,3,j,l1,l2),
     &       Div)
         enddo
         rfvert(1,3,j,1,2) = 0d0
         rfvert(1,3,j,2,1) = 0d0
         enddo

c add everything together
         do l1=1,2
         do l2=1,2

        if(l1.eq.1) then
            poln=-1
         else
            poln=1
          endif

           matv(1,poln,l2,4) = rfvert(1,3,3,l1,l2) ! *2*TR ! amplitude (factor 2 for both fermion directions) * colour factor
     &           * mf(3)/vv    ! top Yukawa
     &           * w4l(h,1)/dcmplx(qaa(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator


           matv(2,poln,l2,4) = matv(1,poln,l2,4) 

         enddo
         enddo
       endif

c       print*, "MH",qaa(4,1)
c       print*, "psq", dotrr(minus_qaa,minus_qaa)


c      stop


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(lvert) then

        do l1=1,2
         do l2=1,2
             do i=1,2
           do j=1,2
             rfvert(i,j,1,l1,l2)=0 ! t
             rfvert(i,j,2,l1,l2)=0 ! b
             rfvert(i,j,3,l1,l2)=0 ! t
             rfvert(i,j,4,l1,l2)=0 ! b
           enddo
           enddo
           enddo
           enddo


cc 1-2
       do j=3,4 ! no contribution from massless loops
          comp=1
       do l1=1,2
       do l2=1,2
 
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     p12,mpg,minus_qaa(0),
     &     curr1(0,l1),MyEps(0,l2),wz(1,1),
     &     scale,comp,comp,rfvert(1,1,j,l1,l2),rgfvert(1,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     p12,mpg,minus_qaa(0),
     &     curr1(0,l1),MyEps(0,l2),wa(1,1),
     &     scale,comp,comp,rfvert(1,2,j,l1,l2),rgfvert(1,2,j,l1,l2),
     &     Div)

       enddo
       enddo
       enddo
cc 2-1
       do j=3,4 ! no contribution from massless loops
          comp=1
       do l1=1,2
       do l2=1,2

         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     mpg,p12,minus_qaa(0),
     &     MyEps(0,l2),curr1(0,l1),wz(1,1),
     &     scale,comp,comp,rfvert(2,1,j,l1,l2),rgfvert(2,1,j,l1,l2),
     &     Div)
         comp=-1
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),1,+1),
     &     1d0,1d0,clr(cpl(j),1,-1),
     &     mf(j),
     &     mpg,p12,minus_qaa(0),
     &     MyEps(0,l2),curr1(0,l1),wa(1,1),
     &     scale,comp,comp,rfvert(2,2,j,l1,l2),rgfvert(2,2,j,l1,l2),
     &     Div)

       enddo
       enddo
       enddo

       wardidtest1 = .true.
       do i=1,2
           if (rgfvert(i,1,1,1,1).gt.1d-1) then 
c            if (wardidtest1) then
c              write(*,*) "MARKER1: gauge test failed", j,i,rgfbox(i,j,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER4: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            else
c              write(*,*) "MARKER2: gauge test failed", j,i,rgfbox(i,j,1)
c              if (Div.ne.0) then
c              write(*,*) "MARKER5: gauge test failed", j,i,rfbox(i,j,1)
c              endif
c            endif
             wardidtest1=.false.
           endif
       enddo
c add everything together
       if (wardidtest1) then
         do l1=1,2
         do l2=1,2
           mv1 = 0d0
           mv2 = 0d0
           mv3 = 0d0
           do i=1,2
           do j=1,2
             mv1 = mv1 + 2*rfvert(i,j,1,l1,l2) ! u,c
     &                  + rfvert(i,j,3,l1,l2) ! t
             mv2 = mv2 + 2*rfvert(i,j,2,l1,l2) ! d,s
     &                  + rfvert(i,j,4,l1,l2) ! b
           enddo
           enddo

        if(l1.eq.1) then
            poln=-1
         else
            poln=1
          endif
c and multiply with color factor TR

          matv(1,poln,l2,5) =0d0
          matv(2,poln,l2,5) =0d0
cFC
          matv(1,poln,l2,5) = (mv1 + mv2) * TR 
          matv(2,poln,l2,5) = (mv1 + mv2) * TR 
c        print*,"matv(1,poln,l2,5)",abs(matv(1,poln,l2,5)),mv1,mv2
         enddo
         enddo
 
        else !wardidtest
 	      box2countm = box2countm+1
              box2countps(ps_number) = box2countps(ps_number)+1 
c              print*, "NO GUAGE BOX2",ps_number, box2countm,box2countps(ps_number)
       endif
         

         endif ! lvert


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif !lbox

      end if ! end lnlo


#endif   ! with NLO


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
c      uuzaj = res(1)
c      ddzaj = res(2)


c      if (lnlo) then
         do l=1,2
            do isig=-1,1,2
               do isig2=1,2
                 mv(l,isig,isig2)=0d0
               enddo
            enddo

            do isig=-1,1,2
               do isig2=1,2 
                  do k=1,5
                   mv(l,isig,isig2)=mv(l,isig,isig2)+matv(l,isig,isig2,k)
                   enddo
                enddo
             enddo   
           enddo

       do l=1,2
           do  isig= -1,1,2
              do isig2=1,2 
                 if (nlo.eq.4) then
                   mv(l,isig,isig2) = als(1,N)*oneo4pi*
     1                (mv(1,isig,isig2) + mm(l,isig,isig2)*K_cont) 

                 elseif(nlo.eq.7) then
                  mv(l,isig,isig2) = als(1,N)*oneo4pi*
     1                (mv(1,isig,isig2) + mm(l,isig,isig2)*K_cont) 
     2                 +factortoadd*mmb(l,isig,isig2)

                 else if(nlo.lt.0) then
                  mv(l,isig,isig2) = als(1,N)*oneo4pi*
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

      uuzaj = res(1)
      ddzaj = res(2)

c      print*, "HERE 3", nlo, uuzaj,ddzaj
         

      return


      end


      function deltakro(x,y)
      implicit none
      real * 8 deltakro
      integer x,y
      if (x.eq.y) then
         deltakro = 1d0
      else
         deltakro = 0d0
      endif
      end

