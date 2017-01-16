c****************subroutine QQAAWJ *************************************
c
c	giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c	Last modified: December 2009
C
c       Francisco Campanario, <francam@particle.uni-karlsruhe.de>
c       Last modified: Jul 2010
C       Mod: including NLO: hexagons,pentagons, boxes and different 
C      checks, including gauge and factorization of the IR 
C      divergences.....
C
C  QQVAAj calculates the matrix elements**2 for 
C     q1 q2    ---->   W+- A A g,   
C     W+- ---> f1-bar f2
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis
c  (comparison with madgraph amplitudes tampr requires to express these in
c   kart. basis also, even though they are given in the helicity basis
c   a priori).
c
c  This code includes only real emission contributions, i.e.
c
c      return udwwz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuwwaj,2:ddwwaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqAAWj(pbar,sign, qbar, gsign, nlo, N, final, sbos
     & ,ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

      integer ps_number

c     electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program      

      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      double precision  pbar(0:3,n_p+n_v), qbar(0:4), q(0:3)
      double precision  res(2), resv(2), final
      double precision  p(0:3,n_p+n_v),fq0(0:4,2),fqa1(0:4,2),fqa2(0:4,2)

      integer  sign(n_p+n_v), gsign, mu, i, j, k, isig, nlo, N
      integer  lh, l, sbos,bos


      double precision pq(0:4,2)
      double precision pgA1(0:4,2), pA1g(0:4,2)
      double precision pgA2(0:4,2), pA2g(0:4,2)
      double precision pgw0(0:4,2), pw0g(0:4,2)


      double precision pwa1(0:4), pwa2(0:4) 
      double precision pqqw(0:4), minus_pqqw(0:4)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psi(2,-1:1,2)

      double complex psia1g(2,-1:1,2,2), psia2g(2,-1:1,2,2)
      double complex psiw0g(2,-1:1,2,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgw0(2,-1:1,2,2)


      double complex braketgA1(2,-1:1,2,2), braketgA2(2,-1:1,2,2)

      double complex psiw0(2,-1:1, 2),psia1(2,-1:1, 2),
     1               psia2(2,-1:1, 2)

      double complex mat(2,-1:1,2,3), madgr(3,-1:1,-1:1,2),matb(2,-1:1,2,3)
      double complex matv(2,-1:1,2,3)
      double complex mm(2, -1:1,2), mv(2, -1:1,2)
      double complex mmb(2, -1:1,2)
      double complex m1kb(2,-1:1,3), m2kb(2,-1:1,3), temp, temp1
      double complex m1(3,2), m2(3,2)
      double complex ma
      double complex janew(0:5,-1:1,14,2)


      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2)
      double complex jqqa(0:5,-1:1,2)
      double complex jqqw(0:5,-1:1)
   
      real*8 dotrr,theta
      double complex im
      double complex dotcc, s1c, dotrc,cnf
      double complex zero 
      external dotcc, s1c, dotrc,dotrr,theta
      integer u(2),d(2)
      save u,k,d
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i

c  helicity selection

      INTEGER h
      COMMON /hcount / h
      COMMON /icount / ie,iu,it

      logical ldebug, ldebugm, linit, christophmods
      logical lnlo, lbox, lpt, lpent,lhex,lbox2
      logical addkmuknu
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false., christophmods=.false.)

c variables for virtual corrections
      DOUBLE PRECISION scale
      double precision c2,c2o4pi,oneo4pi
      parameter (c2=4d0/3d0, c2o4pi=c2/4d0/pi, oneo4pi=1d0/4d0/pi)
c Variables for boxes
      complex*16 j5g(5,2,2),jbg(2,2),jg5n(2)
      complex*16 j5gg(5,2,2),jbgg(2,2),jg5ng(2)
      complex*16 j5gg2(5,2,2),jbgg2(2,2),jg5ng2(2)
      complex*16 nothing(2)
      complex*16 boxcontr1(2),boxcontr1b(2)
      complex*16 boxcontr2(2),boxcontr2b(2)
      complex*16 boxcontrNoAbe3(2)
      real*8 mqwa12(0:3),td
      complex*16 wwa12(0:3,3)
      complex*16 emqwa12(0:3),empg(0:3)
c Variables for pentagons      
      real*8 mqa1(0:3),mqa2(0:3)
      real*8 mqwa1(0:3),mqwa2(0:3)
      complex*16 cmqwa1(0:3),cmqwa2(0:3)
      real*8 p1(0:3),p2(0:3)
      real*8 mpg(0:3),mpw(0:3)
      complex*16 pentcontr1(3),pentcontr1b(3)
      complex*16 pentcontr2(3),pentcontr2b(3)
      complex*16 pentcontr3(3),pentcontr3b(3)
      complex*16 pentcontr4(3),pentcontr4b(3)
      complex*16 pentcontr5(3),pentcontr5b(3)
      complex*16 pentcontr6(3),pentcontr6b(3)
      complex*16 pentcontr7(3),pentcontr7b(3)
      complex*16 pentcontr8(3),pentcontr8b(3)
      complex*16 pentcontr9(3),pentcontr9b(3)
      complex*16 pentcontr10(3),pentcontr10b(3)
      complex*16 pentcontr11(3),pentcontr11b(3)
      complex*16 pentcontr12(3),pentcontr12b(3)
      complex*16 j5a1g(8,4,2),jba1g(4,2)
      complex*16 j5ga1(8,4,2),jbga1(4,2)
      complex*16 j5a2g(8,4,2),jba2g(4,2)
      complex*16 j5ga2(8,4,2),jbga2(4,2)
      complex*16 wwa2(0:3,3),wwa1(0:3,3)
      complex*16 j5ma1g(8,4,2),jbma1g(4,2)
      complex*16 j5mga1(8,4,2),jbmga1(4,2)
      complex*16 j5ma2g(8,4,2),jbma2g(4,2)
      complex*16 j5mga2(8,4,2),jbmga2(4,2)
c NoAbe
      complex*16 jg5a1(2,2),jg5a2(2,2)
      complex*16 jg5ma1(2,2),jg5ma2(2,2)
      complex*16 pentcontrNoAbe1(3),pentcontrNoAbe2(3)
      complex*16 pentcontrNoAbe3(3),pentcontrNoAbe4(3),pentemp(3)
      real*8      accuracyres
c Variables for hexagons
      complex*16 rvga1a2(13,4,2),rbga1a2(4,2)
      complex*16 rgvga1a2(4,4),rgbga1a2(4,4)
      complex*16 rvga2a1(13,4,2),rbga2a1(4,2)
      complex*16 rgvga2a1(4,4),rgbga2a1(4,4)
      complex*16 rva1ga2(13,4,2),rba1ga2(4,2)
      complex*16 rgva1ga2(4,4),rgba1ga2(4,4)
      complex*16 rva2ga1(13,4,2),rba2ga1(4,2)
      complex*16 rgva2ga1(4,4),rgba2ga1(4,4)
      complex*16 rva1a2g(13,4,2),rba1a2g(4,2)
      complex*16 rgva1a2g(4,4),rgba1a2g(4,4)
      complex*16 rva2a1g(13,4,2),rba2a1g(4,2)
      complex*16 rgva2a1g(4,4),rgba2a1g(4,4)
      Integer Div,PerformGauge,mutemp
      complex*16 mv1CA,mv2CA,mv3CA,mv1b,mv2b,mv3b
      complex*16 rvZERO(4,4),rbZERO(4,4)
      complex*16 lns,lnt,lnu,cteEps1
c No abe
      complex*16 rva1a2n(10,3,2),rgva1a2n(3,3)
      complex*16 rva2a1n(10,3,2),rgva2a1n(3,3)
c Variables for fermion loops
      integer flgauge, flcomp
      double precision mf(3)
      double complex rfbox(6,3,2)
      double precision rgfbox(6,3,2) 

      integer*8 thirtytwo


      parameter(thirtytwo=32)
c      
      Logical wardidtest1,wardidtest2,wardidtest3,wardidtest4
      Logical wardidtest5,wardidtest6,wardidtest7,wardidtest8
      Logical  wardidtest9,wardidtest10,wardidtest11,wardidtest12    
      Logical wardidtest13,wardidtest14,wardidtest15,wardidtest16
      Logical wardidtest17,wardidtest18,wardidtest19,wardidtest20
      Logical wardidtest21,wardidtest22,wardidtest23,wardidtest24
      Logical wardidtest25,wardidtest26,wardidtest27,wardidtest28
      Logical  wardidtest29,wardidtest30
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
      complex*16 mmv1,mmv2,ratioFACT,mmv1b,mmv2b
      Complex*16 mmv1CA,mmv2CA,mmv3,mmv4
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
       complex*16 check(2),check1(2)
       integer GaugeAbe,GaugeNoAbe

       double precision RandomNumber
       external RandomNumber

              
c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c      fpials(2) = fpi*alfas

c---------------------------- BEGIN CODE -----------------------

      lnlo = NLO.ne.0       ! include some virtual stuff if T
      lbox = (NLO.eq.-4).or.(NLO.eq.1).or.(NLO.eq.5).or.(NLO.eq.7).or.(NLO.eq.-7)  ! call boxline if T
      lpt = (NLO.eq.1).or.(NLO.eq.-5).or.(NLO.eq.7).or.(NLO.eq.-7)    ! true pentagon contribution
      !lpq = (NLO.eq.5) .or. (NLO.eq.-4)   ! q^mu/m_V terms of Pentagon
      lpent = lpt !.or. lpq
c      Hexagons
      lhex= (NLO.eq.-6) .or. (NLO.eq.7).or. (NLO.eq.-7)      
      lbox2=(NLO.eq.-8)

      addkmuknu=.true.
      if (sbos.lt.0) addkmuknu=.false.
      bos=abs(sbos)

c define flavors of external quarks for the 2 NC subprocesses

c      if (linit) then
c         linit = .false.
         if (bos.eq.311) then
            k = 1               ! udwaaj
            u(1) = 3
            d(1) = 4 
         elseif (bos.eq.411) then
            k = 2               ! duwwaj
            u(2) = 4
            d(2) = 3            
         endif
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
c      endif

      do isig = -1,1,2
         do l = 1,2
            do i = 1,3
               mat(k,isig,l,i) = (0.0d0, 0.0d0)
               matv(k,isig,l,i) = (0.0d0, 0.0d0)
               matb(k,isig,l,i) = (0.0d0, 0.0d0)
            enddo
         enddo
      enddo

c define the internal momenta
 
      
      do mu = 0,3
         do i = 1,n_p+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
         pwa1(mu) = qw0(mu,N) + qa1(mu,N)
         pwa2(mu) = qw0(mu,N) + qa2(mu,N)
      enddo
      pwa1(4) = pwa1(0)**2-pwa1(1)**2-pwa1(2)**2-pwa1(3)**2
      pwa2(4) = pwa2(0)**2-pwa2(1)**2-pwa2(2)**2-pwa2(3)**2

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qw0(0,N),w0(1,N),
     &              psiw0(1,-1,2),fq0(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qw0(0,N),w0(1,N),
     &              psiw0(1,-1,1),fq0(0,1))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,2),fqa1(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,1),fqa1(0,1))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,2),fqa2(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,1),fqa2(0,1))
      
c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2                ! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l)) ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

         do mu = 0,3                   !QCD gaugecheck      
!            eps(mu,l) = qbar(mu)
	    Myeps(mu,l)=eps(mu,l)
         enddo

c--- End of QCD gaugecheck -------------------------------------------    

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
         call ket2r(psi(1,-1,1),.true.,p(0,1),-1,q,eps(0,l),
     &        braketg(1,-1,1,l),pq(0,1)) ! |q,1>_l,isig
         
         call bra2r(psi(1,-1,2),.true.,p(0,2),-1,q,eps(0,l),
     &        braketg(1,-1,2,l),pq(0,2)) ! <2,q|_l,isig2
         
c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)

         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qa1(0,N),atau1(1,N),braketgA1(1,-1,1,l),pgA1(0,1))
         
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qa1(0,N),atau1(1,N),braketgA1(1,-1,2,l),pgA1(0,2))
         
         call ket2r(psia1(1,-1,1),.false.,fqa1(0,1),-1,q,eps(0,l),
     &        psiA1g(1,-1,1,l),pA1g(0,1))
         
         call bra2r(psia1(1,-1,2),.false.,fqa1(0,2),-1,q,eps(0,l),
     &        psiA1g(1,-1,2,l),pA1g(0,2))

         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qa2(0,N),atau2(1,N),braketgA2(1,-1,1,l),pgA2(0,1))
         
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qa2(0,N),atau2(1,N),braketgA2(1,-1,2,l),pgA2(0,2))
         
         call ket2r(psia2(1,-1,1),.false.,fqa2(0,1),-1,q,eps(0,l),
     &        psiA2g(1,-1,1,l),pA2g(0,1))
         
         call bra2r(psia2(1,-1,2),.false.,fqa2(0,2),-1,q,eps(0,l),
     &        psiA2g(1,-1,2,l),pA2g(0,2))
         
         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qw0(0,N),w0(1,N),braketgw0(1,-1,1,l),pgw0(0,1)) 

         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qw0(0,N),w0(1,N),braketgw0(1,-1,2,l),pgw0(0,2))

         call ket2r(psiw0(1,-1,1),.false.,fq0(0,1),-1,
     &        q,eps(0,l),psiw0g(1,-1,1,l),pw0g(0,1))

         call bra2r(psiw0(1,-1,2),.false.,fq0(0,2),-1,
     &        q,eps(0,l),psiw0g(1,-1,2,l),pw0g(0,2))
         
      enddo    ! 2 gluon polarizations

c Calculation of processes qq -> W A A g ---------------------

c   u d -> W A A g

      do l=1,2

         m1(1,l) =
     &      -s1c(psiA1(1,-1,2),atau2(1,N),.true.,-1,braketgw0(1,-1,1,l))
     &      -s1c(psiA2(1,-1,2),atau1(1,N),.true.,-1,braketgw0(1,-1,1,l))
     &      -s1c(psiA1(1,-1,2),atau2(1,N),.true.,-1,psiw0g(1,-1,1,l))
     &      -s1c(psiA2(1,-1,2),atau1(1,N),.true.,-1,psiw0g(1,-1,1,l))
     &      -s1c(braketgA1(1,-1,2,l),atau2(1,N),.true.,-1,psiw0(1,-1,1))
     &      -s1c(braketgA2(1,-1,2,l),atau1(1,N),.true.,-1,psiw0(1,-1,1))
     &      -s1c(psiA1g(1,-1,2,l),atau2(1,N),.true.,-1,psiw0(1,-1,1))
     &      -s1c(psiA2g(1,-1,2,l),atau1(1,N),.true.,-1,psiw0(1,-1,1))

         m1(2,l) =
     &      -s1c(psiw0(1,-1,2),atau1(1,N),.true.,-1,braketgA2(1,-1,1,l))
     &      -s1c(psiw0(1,-1,2),atau2(1,N),.true.,-1,braketgA1(1,-1,1,l))
     &      -s1c(psiw0(1,-1,2),atau1(1,N),.true.,-1,psiA2g(1,-1,1,l))
     &      -s1c(psiw0(1,-1,2),atau2(1,N),.true.,-1,psiA1g(1,-1,1,l))
     &      -s1c(braketgw0(1,-1,2,l),atau1(1,N),.true.,-1,psiA2(1,-1,1))
     &      -s1c(braketgw0(1,-1,2,l),atau2(1,N),.true.,-1,psiA1(1,-1,1))
     &      -s1c(psiw0g(1,-1,2,l),atau1(1,N),.true.,-1,psiA2(1,-1,1))
     &      -s1c(psiw0g(1,-1,2,l),atau2(1,N),.true.,-1,psiA1(1,-1,1))

         m1(3,l) =
     &      -s1c(psiA1(1,-1,2),w0(1,N),.true.,-1,braketgA2(1,-1,1,l))
     &      -s1c(psiA2(1,-1,2),w0(1,N),.true.,-1,braketgA1(1,-1,1,l))
     &      -s1c(psiA1(1,-1,2),w0(1,N),.true.,-1,psiA2g(1,-1,1,l))
     &      -s1c(psiA2(1,-1,2),w0(1,N),.true.,-1,psiA1g(1,-1,1,l))
     &      -s1c(braketgA1(1,-1,2,l),w0(1,N),.true.,-1,psiA2(1,-1,1))
     &      -s1c(braketgA2(1,-1,2,l),w0(1,N),.true.,-1,psiA1(1,-1,1))
     &      -s1c(psiA1g(1,-1,2,l),w0(1,N),.true.,-1,psiA2(1,-1,1))
     &      -s1c(psiA2g(1,-1,2,l),w0(1,N),.true.,-1,psiA1(1,-1,1))

         mat(k,-1,l,1) = (m1(1,l) * clr(d(k),1,-1)**2 
     &                  + m1(2,l) * clr(u(k),1,-1)**2
     &                  + m1(3,l) * clr(u(k),1,-1)*clr(d(k),1,-1)) 
     &                  * clr(3,3,-1) 

      enddo


c Calculation of processes qq -> W A g -> W A A g -------------

c qq -> A W(-> W A) g

      do l=1,2

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgA1(1,-1,1,l),pgA1(0,1),janew(0,-1,1,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgA2(1,-1,1,l),pgA2(0,1),janew(0,-1,2,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiA1g(1,-1,1,l),pA1g(0,1),janew(0,-1,3,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiA2g(1,-1,1,l),pA2g(0,1),janew(0,-1,4,l) )

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiA1(1,-1,1),fqA1(0,1),janew(0,-1,5,l) )

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiA2(1,-1,1),fqA2(0,1),janew(0,-1,6,l) )


         call curr6(-1,psiA1(1,-1,2),fqA1(0,2),
     &              braketg(1,-1,1,l),pq(0,1),janew(0,-1,7,l) )

         call curr6(-1,psiA2(1,-1,2),fqA2(0,2),
     &              braketg(1,-1,1,l),pq(0,1),janew(0,-1,8,l) )

         call curr6(-1,psiA1g(1,-1,2,l),pA1g(0,2),
     &              psi(1,-1,1),p(0,1),janew(0,-1,9,l) )

         call curr6(-1,psiA2g(1,-1,2,l),pA2g(0,2),
     &              psi(1,-1,1),p(0,1),janew(0,-1,10,l) )

         call curr6(-1,braketgA1(1,-1,2,l),pgA1(0,2),
     &              psi(1,-1,1),p(0,1),janew(0,-1,11,l) )

         call curr6(-1,braketgA2(1,-1,2,l),pgA2(0,2),
     &              psi(1,-1,1),p(0,1),janew(0,-1,12,l) )

         do mu = 0,5
            janew(mu,-1,13,l) = janew(mu,-1,1,l) *clr(u(k),1,-1)
     &                        + janew(mu,-1,3,l) *clr(u(k),1,-1)
     &                        + janew(mu,-1,5,l) *clr(u(k),1,-1)
     &                        + janew(mu,-1,7,l) *clr(d(k),1,-1)
     &                        + janew(mu,-1,9,l) *clr(d(k),1,-1)
     &                        + janew(mu,-1,11,l)*clr(d(k),1,-1)
            janew(mu,-1,14,l) = janew(mu,-1,2,l) *clr(u(k),1,-1) 
     &                        + janew(mu,-1,4,l) *clr(u(k),1,-1)
     &                        + janew(mu,-1,6,l) *clr(u(k),1,-1)
     &                        + janew(mu,-1,8,l) *clr(d(k),1,-1)
     &                        + janew(mu,-1,10,l)*clr(d(k),1,-1)
     &                        + janew(mu,-1,12,l)*clr(d(k),1,-1)
         enddo

c contract with atoww/ztoww tensor and
c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         if (bos.eq.311) then
            m1(1,l) = -dotcc(wpwa2(0,N),janew(0,-1,13,l))
            m2(1,l) = -dotcc(wpwa1(0,N),janew(0,-1,14,l))
          if (addkmuknu) then
            m1(1,l) = m1(1,l)+zm2i(3)*dotrc(pwa2(0),wpwa2(0,N))
     1           *dotrc(pwa2(0),janew(0,-1,13,l))
            m2(1,l) = m2(1,l)+zm2i(3)*dotrc(pwa1(0),wpwa1(0,N))
     1           *dotrc(pwa1(0),janew(0,-1,14,l))
          endif
         elseif (bos.eq.411) then
            m1(1,l) = -dotcc(wmwa2(0,N),janew(0,-1,13,l))
            m2(1,l) = -dotcc(wmwa1(0,N),janew(0,-1,14,l))
          if (addkmuknu) then
            m1(1,l) = m1(1,l)+zm2i(3)*dotrc(pwa2(0),wmwa2(0,N))
     1           *dotrc(pwa2(0),janew(0,-1,13,l))
            m2(1,l) = m2(1,l)+zm2i(3)*dotrc(pwa1(0),wmwa1(0,N))
     1           *dotrc(pwa1(0),janew(0,-1,14,l))
          endif
         endif
         
         mat(k,-1,l,2) = (m1(1,l)+m2(1,l))*clr(3,3,-1)
         
      enddo                     ! gluon polarization
      

c Calculation of processes qq -> A/Z g -> W W A g ---------------------

      do l =1,2
         
         call curr6(-1,psi(1,-1,2),p(0,2),
     &        braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))
         
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &        psi(1,-1,1),p(0,1),jqagq(0,-1,l))
         
         do mu = 0,5
            jqqa(mu,-1,l) = jqgaq(mu,-1,l) + jqagq(mu,-1,l)
         enddo
         
c     contract with atowwa/ztowwa tensor
         if (bos.eq.311) then
            ma = -dotcc(wpwaa(0,N),jqqa(0,-1,l))
         elseif (bos.eq.411) then
            ma = -dotcc(wmwaa(0,N),jqqa(0,-1,l))
         endif

         mat(k,-1,l,3) = ma*clr(3,3,-1)
         
      enddo                     ! gluon polarization
      
c**************************************************************
c
c    Virtual contributions
c
c**************************************************************
#ifdef WITH_TRIBOSONJET
#ifdef WITH_NLO
      if (lnlo) then
C Som variables to perform checks of the virtual
C contributions
C To checl Factorization of Divergences
C  Div=1 -> 1/Eps;   Div=2 ->1/Eps^2
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


      finiteubd=(25*als(1,1)*CA)/(9D0*Pi) + (5*als(1,1)*CF)/Pi - 
     &  (11*als(1,1)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,1)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,1)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,1)*CA*Pi)/24d0 - 
     &  (7*als(1,1)*CF*Pi)/12d0 - (8*als(1,1)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,1)*EulerGamma*Nf*TR)/(3D0*Pi) - 
     -  (3*als(1,1)*CA*Log(mursq1/S))/(4D0*Pi) + 
     &  (3*als(1,1)*CF*Log(mursq1/S))/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     -  (als(1,1)*CF*EulerGamma*Log(mursq1/S))/Pi - 
     &  (als(1,1)*CA*Log(mursq1/S)**2)/(4D0*Pi) + 
     &  (als(1,1)*CF*Log(mursq1/S)**2)/(2D0*Pi) + 
     -  (5*als(1,1)*CA*Log(-(mursq1/T)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,1)*Nf*TR*Log(-(mursq1/T)))/(6D0*Pi) + 
     -  (als(1,1)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     &  (5*als(1,1)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,1)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,1)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

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

      finiteubgcr=(25*als(1,1)*CA)/(9D0*Pi) + (5*als(1,1)*CF)/Pi - 
     &  (11*als(1,1)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,1)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,1)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,1)*CA*Pi)/24d0 - (7*als(1,1)*CF*Pi)/12d0 - 
     &  (8*als(1,1)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,1)*EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*als(1,1)*CA*Log(mursq1/S))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,1)*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (als(1,1)*CA*Log(mursq1/S)**2)/(4D0*Pi) - 
     &  (3*als(1,1)*CA*Log(-(mursq1/T)))/(4D0*Pi) + 
     &  (3*als(1,1)*CF*Log(-(mursq1/T)))/(2D0*Pi) + 
     -  (als(1,1)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,1)*CF*EulerGamma*Log(-(mursq1/T)))/Pi - 
     &  (als(1,1)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (als(1,1)*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + 
     &  (5*als(1,1)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,1)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,1)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

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

      finiteubg=(25*als(1,1)*CA)/(9D0*Pi) + 
     &  (5*als(1,1)*CF)/Pi - (11*als(1,1)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,1)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,1)*CA*EulerGamma**2)/(4D0*Pi) + 
     -  (als(1,1)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,1)*CA*Pi)/24d0 - (7*als(1,1)*CF*Pi)/12d0 - 
     &  (8*als(1,1)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,1)*EulerGamma*Nf*TR)/(3D0*Pi) + 
     -  (5*als(1,1)*CA*Log(mursq1/S))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,1)*Nf*TR*Log(mursq1/S))/(6D0*Pi) + 
     -  (als(1,1)*CA*Log(mursq1/S)**2)/(4D0*Pi) - 
     &  (3*als(1,1)*CA*Log(-(mursq1/T)))/(4D0*Pi) + 
     &  (3*als(1,1)*CF*Log(-(mursq1/T)))/(2D0*Pi) + 
     -  (als(1,1)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,1)*CF*EulerGamma*Log(-(mursq1/T)))/Pi - 
     &  (als(1,1)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     -  (als(1,1)*CF*Log(-(mursq1/T))**2)/(2D0*Pi) + 
     &  (5*als(1,1)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,1)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     -  (als(1,1)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,1)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

       factortoadd=finiteubg

          endif

      endif     ! Finish gsign

      factortoadd=factortoadd/2d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                 Here start the virtual staff
C
!         Lhex=.false.
!	 LPENT=.FALSE.
!	 LBOX=.FALSE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
c Calculation of processes qq -> W A A g ---------------------
c   u d -> W A A g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	if(lhex) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        scale=musq!mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))

	  do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
	       mpg(mu) = -q(mu)
               mqa1(mu) = -qa1(mu,1)
               mqa2(mu) = -qa2(mu,1)
	       mpw(mu)= -qw0(mu,1)
 	   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  Here I do Ward identity test for the different momenta:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If(PerformGauge.eq.1) then
         do mu = 0,3
	    Myeps(mu,1)=mpg(mu)
	    Myeps(mu,2)=mpg(mu)
	  enddo
      elseif(PerformGauge.eq.2) then
         do mu = 0,3
	     mutemp=mu+1
	    atau1(mutemp,N) = mqa1(mu)
	 enddo
      elseif(PerformGauge.eq.3) then
         do mu = 0,3
	    mutemp=mu+1
	    atau2(mutemp,N) = mqa2(mu)
	 enddo
        elseif(PerformGauge.eq.4) then
         do mu = 0,3
	    mutemp=mu+1
	    w0(mutemp,N) = mpw(mu)
	 enddo
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C the results are store in rvga1a2(x,y,l)
C The name indicate the order of gluon and the photons
C x: Number of graphs, 
C x=1 collect results proportional CF
C x=2 collect results proportional to CF-1/2CA
C y: Indicates the position of the W. This will affect 
C    the couplings
C l: Polarization of the gluonl l=1,l=2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Hexagon Abelian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       gaugeAbe=4*RandomNumber()+1
!       Print*,'gaugeAbe',gaugeAbe
       if(gaugeAbe.eq.5) then
          gaugeAbe=4
       endif
C g A1 A2 W
       call hexlineABETotal(p1,mpg,mqa1,mqa2,mpw,p2,
     1      psi(1,-1,2),psi(1,-1,1),
     2      Myeps(0,1),atau1(1,N),atau2(1,N),w0(1,N),-1,scale,
     3      1,2,1,
     4      rvga1a2(1,4,1),rbga1a2(4,1),rgvga1a2(1,4),rgbga1a2(1,4),
     5      Div,GaugeAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time   
      call hexlineABETotal(p1,mpg,mqa1,mqa2,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,2),atau1(1,N),atau2(1,N),w0(1,N),-1,scale,
     3     1,2,-1,rvga1a2(1,4,2),rbga1a2(4,2),rvZERO(1,4),rbZERO(1,4),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C g A2 A1 W
        call hexlineABETotal(p1,mpg,mqa2,mqa1,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau2(1,N),atau1(1,N),w0(1,N),-1,scale,
     3     1,2,1,rvga2a1(1,4,1),rbga2a1(4,1),rgvga2a1(1,4),rgbga2a1(1,4),Div,GaugeAbe)
c  gluon helicity 2    
        call hexlineABETotal(p1,mpg,mqa2,mqa1,mpw,p2,
     1       psi(1,-1,2),psi(1,-1,1),
     2       Myeps(0,2),atau2(1,N),atau1(1,N),w0(1,N),-1,scale,
     3       1,2,-1,
     4       rvga2a1(1,4,2),rbga2a1(4,2),rvZERO(1,4),rbZERO(1,4),
     5       Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C A1 g A2 W 
        call hexlineABETotal(p1,mqa1,mpg,mqa2,mpw,p2,
     1       psi(1,-1,2),psi(1,-1,1),
     2       atau1(1,N),Myeps(0,1),atau2(1,N),w0(1,N),-1,scale,
     3       1,3,1,
     4       rva1ga2(1,4,1),rba1ga2(4,1),rgva1ga2(1,4),rgba1ga2(1,4),
     5       Div,GaugeAbe)
c  gluon helicity 2   
        call hexlineABETotal(p1,mqa1,mpg,mqa2,mpw,p2,
     1       psi(1,-1,2),psi(1,-1,1),
     2       atau1(1,N),Myeps(0,2),atau2(1,N),w0(1,N),-1,scale,
     3       1,3,-1,
     4       rva1ga2(1,4,2),rba1ga2(4,2),rvZERO(1,4),rbZERO(1,4),
     5       Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 g A1 W
         call hexlineABETotal(p1,mqa2,mpg,mqa1,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,1),atau1(1,N),w0(1,N),-1,scale,
     3     1,3,1,rva2ga1(1,4,1),rba2ga1(4,1),rgva2ga1(1,4),rgba2ga1(1,4),Div,GaugeAbe)
c  gluon helicity 2     
         call hexlineABETotal(p1,mqa2,mpg,mqa1,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,2),atau1(1,N),w0(1,N),-1,scale,
     3     1,3,-1,rva2ga1(1,4,2),rba2ga1(4,2),rvZERO(1,4),rbZERO(1,4),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 g W
        call hexlineABETotal(p1,mqa1,mqa2,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),Myeps(0,1),w0(1,N),-1,scale,
     3     1,4,1,rva1a2g(1,4,1),rba1a2g(4,1),rgva1a2g(1,4),rgba1a2g(1,4),Div,GaugeAbe)
c  gluon helicity 2    
        call hexlineABETotal(p1,mqa1,mqa2,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),Myeps(0,2),w0(1,N),-1,scale,
     3     1,4,-1,rva1a2g(1,4,2),rba1a2g(4,2),rvZERO(1,4),rbZERO(1,4),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 A1 g W
 	call hexlineABETotal(p1,mqa2,mqa1,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),Myeps(0,1),w0(1,N),-1,scale,
     3     1,4,1,rva2a1g(1,4,1),rba2a1g(4,1),rgva2a1g(1,4),rgba2a1g(1,4),Div,GaugeAbe)
c  gluon helicity 2          
 	call hexlineABETotal(p1,mqa2,mqa1,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),Myeps(0,2),w0(1,N),-1,scale,
     3     1,4,-1,rva2a1g(1,4,2),rba2a1g(4,2),rvZERO(1,4),rbZERO(1,4),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A1 A2 W g
       call hexlineABETotal(p1,mqa1,mqa2,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva1a2g(1,3,1),rba1a2g(3,1),rgva1a2g(1,3),rgba1a2g(1,3),Div,GaugeAbe)
c  gluon helicity 2    
       call hexlineABETotal(p1,mqa1,mqa2,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva1a2g(1,3,2),rba1a2g(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c A2 A1 W g
       call hexlineABETotal(p1,mqa2,mqa1,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva2a1g(1,3,1),rba2a1g(3,1),rgva2a1g(1,3),rgba2a1g(1,3),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mqa2,mqa1,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva2a1g(1,3,2),rba2a1g(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c g W A1 A2
      call hexlineABETotal(p1,mpg,mpw,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),w0(1,N),atau1(1,N),atau2(1,N),-1,scale,
     3     1,2,1,rvga1a2(1,2,1),rbga1a2(2,1),rgvga1a2(1,2),rgbga1a2(1,2),Div,GaugeAbe)
c  gluon helicity 2 
      call hexlineABETotal(p1,mpg,mpw,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,2),w0(1,N),atau1(1,N),atau2(1,N),-1,scale,
     3     1,2,-1,rvga1a2(1,2,2),rbga1a2(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C g W A2 A1
       call hexlineABETotal(p1,mpg,mpw,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),w0(1,N),atau2(1,N),atau1(1,N),-1,scale,
     3     1,2,1,rvga2a1(1,2,1),rbga2a1(2,1),rgvga2a1(1,2),rgbga2a1(1,2),Div,GaugeAbe)
c  gluon helicity 2      
       call hexlineABETotal(p1,mpg,mpw,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,2),w0(1,N),atau2(1,N),atau1(1,N),-1,scale,
     3     1,2,-1,rvga2a1(1,2,2),rbga2a1(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W g A1 A2
       call hexlineABETotal(p1,mpw,mpg,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,1),atau1(1,N),atau2(1,N),-1,scale,
     3     1,3,1,rvga1a2(1,1,1),rbga1a2(1,1),rgvga1a2(1,1),rgbga1a2(1,1),Div,GaugeAbe)
c  gluon helicity 2   
       call hexlineABETotal(p1,mpw,mpg,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,2),atau1(1,N),atau2(1,N),-1,scale,
     3     1,3,-1,rvga1a2(1,1,2),rbga1a2(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W g A2 A1
      call hexlineABETotal(p1,mpw,mpg,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,1),atau2(1,N),atau1(1,N),-1,scale,
     3     1,3,1,rvga2a1(1,1,1),rbga2a1(1,1),rgvga2a1(1,1),rgbga2a1(1,1),Div,GaugeAbe)
c  gluon helicity 2        
      call hexlineABETotal(p1,mpw,mpg,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,2),atau2(1,N),atau1(1,N),-1,scale,
     3     1,3,-1,rvga2a1(1,1,2),rbga2a1(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A1 g A2
      call hexlineABETotal(p1,mpw,mqa1,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),Myeps(0,1),atau2(1,N),-1,scale,
     3     1,4,1,rva1ga2(1,1,1),rba1ga2(1,1),rgva1ga2(1,1),rgba1ga2(1,1),Div,GaugeAbe)
c  gluon helicity 2          
      call hexlineABETotal(p1,mpw,mqa1,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),Myeps(0,2),atau2(1,N),-1,scale,
     3     1,4,-1,rva1ga2(1,1,2),rba1ga2(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A2 g A1
        call hexlineABETotal(p1,mpw,mqa2,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),Myeps(0,1),atau1(1,N),-1,scale,
     3     1,4,1,rva2ga1(1,1,1),rba2ga1(1,1),rgva2ga1(1,1),rgba2ga1(1,1),Div,GaugeAbe)
c  gluon helicity 2   
        call hexlineABETotal(p1,mpw,mqa2,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),Myeps(0,2),atau1(1,N),-1,scale,
     3     1,4,-1,rva2ga1(1,1,2),rba2ga1(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A1 A2 g
       call hexlineABETotal(p1,mpw,mqa1,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva1a2g(1,1,1),rba1a2g(1,1),rgva1a2g(1,1),rgba1a2g(1,1),Div,GaugeAbe)
c  gluon helicity 2          
       call hexlineABETotal(p1,mpw,mqa1,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva1a2g(1,1,2),rba1a2g(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A2 A1 g
       call hexlineABETotal(p1,mpw,mqa2,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva2a1g(1,1,1),rba2a1g(1,1),rgva2a1g(1,1),rgba2a1g(1,1),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mpw,mqa2,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva2a1g(1,1,2),rba2a1g(1,2),rvZERO(1,1),rbZERO(1,1),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C g A1 W A2 
      call hexlineABETotal(p1,mpg,mqa1,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau1(1,N),w0(1,N),atau2(1,N),-1,scale,
     3     1,2,1,rvga1a2(1,3,1),rbga1a2(3,1),rgvga1a2(1,3),rgbga1a2(1,3),Div,GaugeAbe)
c  gluon helicity 2   
      call hexlineABETotal(p1,mpg,mqa1,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,2),atau1(1,N),w0(1,N),atau2(1,N),-1,scale,
     3     1,2,-1,rvga1a2(1,3,2),rbga1a2(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
Cg A2 W A1
       call hexlineABETotal(p1,mpg,mqa2,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau2(1,N),w0(1,N),atau1(1,N),-1,scale,
     3     1,2,1,rvga2a1(1,3,1),rbga2a1(3,1),rgvga2a1(1,3),rgbga2a1(1,3),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mpg,mqa2,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,2),atau2(1,N),w0(1,N),atau1(1,N),-1,scale,
     3     1,2,-1,rvga2a1(1,3,2),rbga2a1(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A1 g W A2
       call hexlineABETotal(p1,mqa1,mpg,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),Myeps(0,1),w0(1,N),atau2(1,N),-1,scale,
     3     1,3,1,rva1ga2(1,3,1),rba1ga2(3,1),rgva1ga2(1,3),rgba1ga2(1,3),Div,GaugeAbe)
c  gluon helicity 2    
       call hexlineABETotal(p1,mqa1,mpg,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),Myeps(0,2),w0(1,N),atau2(1,N),-1,scale,
     3     1,3,-1,rva1ga2(1,3,2),rba1ga2(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 g W A1
       call hexlineABETotal(p1,mqa2,mpg,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,1),w0(1,N),atau1(1,N),-1,scale,
     3     1,3,1,rva2ga1(1,3,1),rba2ga1(3,1),rgva2ga1(1,3),rgba2ga1(1,3),Div,GaugeAbe)
c  gluon helicity 2  
       call hexlineABETotal(p1,mqa2,mpg,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,2),w0(1,N),atau1(1,N),-1,scale,
     3     1,3,-1,rva2ga1(1,3,2),rba2ga1(3,2),rvZERO(1,3),rbZERO(1,3),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A1 W g A2
       call hexlineABETotal(p1,mqa1,mpw,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),Myeps(0,1),atau2(1,N),-1,scale,
     3     1,4,1,rva1ga2(1,2,1),rba1ga2(2,1),rgva1ga2(1,2),rgba1ga2(1,2),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mqa1,mpw,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),Myeps(0,2),atau2(1,N),-1,scale,
     3     1,4,-1,rva1ga2(1,2,2),rba1ga2(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 W g A1
        call hexlineABETotal(p1,mqa2,mpw,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),Myeps(0,1),atau1(1,N),-1,scale,
     3     1,4,1,rva2ga1(1,2,1),rba2ga1(2,1),rgva2ga1(1,2),rgba2ga1(1,2),Div,GaugeAbe)
c  gluon helicity 2     
        call hexlineABETotal(p1,mqa2,mpw,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),Myeps(0,2),atau1(1,N),-1,scale,
     3     1,4,-1,rva2ga1(1,2,2),rba2ga1(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A1 W A2 g
      call hexlineABETotal(p1,mqa1,mpw,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva1a2g(1,2,1),rba1a2g(2,1),rgva1a2g(1,2),rgba1a2g(1,2),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mqa1,mpw,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva1a2g(1,2,2),rba1a2g(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A2 W A1 g
       call hexlineABETotal(p1,mqa2,mpw,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     1,5,1,rva2a1g(1,2,1),rba2a1g(2,1),rgva2a1g(1,2),rgba2a1g(1,2),Div,GaugeAbe)
c  gluon helicity 2     
       call hexlineABETotal(p1,mqa2,mpw,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,2),-1,scale,
     3     1,5,-1,rva2a1g(1,2,2),rba2a1g(2,2),rvZERO(1,2),rbZERO(1,2),Div,-1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Hexagon NoAbelian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       gaugeNoAbe=3*RandomNumber()+1
       if(gaugeNoAbe.eq.4) then
          gaugeNoAbe=3
       endif
!       Print*,'gaugeNoAbe',gaugeNoAbe
C A1 A2 W
	call HexlineNoAbe1gTotal(p1,mqa1,mqa2,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva1a2n(1,3,1),rgva1a2n(1,3),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
	call test_HexlineNoAbe1g(p1,mqa1,mqa2,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva1a2n(1,3,2),temp)
C A2 A1 W
 	call HexlineNoAbe1gTotal(p1,mqa2,mqa1,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva2a1n(1,3,1),rgva2a1n(1,3),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
 	call test_HexlineNoAbe1g(p1,mqa2,mqa1,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva2a1n(1,3,2),temp)
CCCCCCCCCCCCCCCCCC
C A1 W A2
	call HexlineNoAbe1gTotal(p1,mqa1,mpw,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva1a2n(1,2,1),rgva1a2n(1,2),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time  
	call test_HexlineNoAbe1g(p1,mqa1,mpw,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva1a2n(1,2,2),temp)
C A2 W A1
	call HexlineNoAbe1gTotal(p1,mqa2,mpw,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva2a1n(1,2,1),rgva2a1n(1,2),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
	call test_HexlineNoAbe1g(p1,mqa2,mpw,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva2a1n(1,2,2),temp)
CCCCCCCCCCCCCCC
C W A1 A2
	call HexlineNoAbe1gTotal(p1,mpw,mqa1,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva1a2n(1,1,1),rgva1a2n(1,1),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
	call test_HexlineNoAbe1g(p1,mpw,mqa1,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva1a2n(1,1,2),rgva1a2n(1,1))
C W A2 A1 
	call HexlineNoAbe1gTotal(p1,mpw,mqa2,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     1,1,1,rva2a1n(1,1,1),rgva2a1n(1,1),Div,GaugeNoAbe)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time   
	call test_HexlineNoAbe1g(p1,mpw,mqa2,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,2),-1,scale,
     3     1,1,-1,rva2a1n(1,1,2),rgva2a1n(1,1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC Call Gauge for Hexlines     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Hexagon Abelian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c g A1 A2 W
        call HexGaugeABE(p1,mpg,mqa1,mqa2,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau1(1,N),atau2(1,N),w0(1,N),-1,scale,
     3     rgvga1a2(1,4),rgbga1a2(1,4), wardidtest1,GaugeAbe)
c g A2 A1 W
        call HexGaugeABE(p1,mpg,mqa2,mqa1,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau2(1,N),atau1(1,N),w0(1,N),-1,scale,
     3     rgvga2a1(1,4),rgbga2a1(1,4), wardidtest2,GaugeAbe)
c A1 g A2 W 
       call HexGaugeABE(p1,mqa1,mpg,mqa2,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),Myeps(0,1),atau2(1,N),w0(1,N),-1,scale,
     3     rgva1ga2(1,4),rgba1ga2(1,4), wardidtest3,GaugeAbe)
c A2 g A1 W
        call HexGaugeABE(p1,mqa2,mpg,mqa1,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,1),atau1(1,N),w0(1,N),-1,scale,
     3     rgva2ga1(1,4),rgba2ga1(1,4), wardidtest4,GaugeAbe)
c A1 A2 g W
       call HexGaugeABE(p1,mqa1,mqa2,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),Myeps(0,1),w0(1,N),-1,scale,
     3     rgva1a2g(1,4),rgba1a2g(1,4), wardidtest5,GaugeAbe)
c A2 A1 g W
	call HexGaugeABE(p1,mqa2,mqa1,mpg,mpw,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),Myeps(0,1),w0(1,N),-1,scale,
     3     rgva2a1g(1,4),rgba2a1g(1,4), wardidtest6,GaugeAbe)
c A1 A2 W g
       call HexGaugeABE(p1,mqa1,mqa2,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2g(1,3),rgba1a2g(1,3), wardidtest7,GaugeAbe)
c A2 A1 W g
       call HexGaugeABE(p1,mqa2,mqa1,mpw,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1g(1,3),rgba2a1g(1,3), wardidtest8,GaugeAbe)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c g W A1 A2
      call HexGaugeABE(p1,mpg,mpw,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),w0(1,N),atau1(1,N),atau2(1,N),-1,scale,
     3     rgvga1a2(1,2),rgbga1a2(1,2), wardidtest9,GaugeAbe)
c g W A2 A1
       call HexGaugeABE(p1,mpg,mpw,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),w0(1,N),atau2(1,N),atau1(1,N),-1,scale,
     3     rgvga2a1(1,2),rgbga2a1(1,2), wardidtest10,GaugeAbe)
c W g A1 A2
       call HexGaugeABE(p1,mpw,mpg,mqa1,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,1),atau1(1,N),atau2(1,N),-1,scale,
     3     rgvga1a2(1,1),rgbga1a2(1,1), wardidtest11,GaugeAbe)
c W g A2 A1
      call HexGaugeABE(p1,mpw,mpg,mqa2,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),Myeps(0,1),atau2(1,N),atau1(1,N),-1,scale,
     3     rgvga2a1(1,1),rgbga2a1(1,1), wardidtest12,GaugeAbe)
c W A1 g A2
      call HexGaugeABE(p1,mpw,mqa1,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),Myeps(0,1),atau2(1,N),-1,scale,
     3     rgva1ga2(1,1),rgba1ga2(1,1), wardidtest13,GaugeAbe)
c W A2 g A1
        call HexGaugeABE(p1,mpw,mqa2,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),Myeps(0,1),atau1(1,N),-1,scale,
     3     rgva2ga1(1,1),rgba2ga1(1,1), wardidtest14,GaugeAbe)

c W A1 A2 g
       call HexGaugeABE(p1,mpw,mqa1,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2g(1,1),rgba1a2g(1,1), wardidtest15,GaugeAbe)
c W A2 A1 g
       call HexGaugeABE(p1,mpw,mqa2,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1g(1,1),rgba2a1g(1,1), wardidtest16,GaugeAbe)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c g A1 W A2 
      call HexGaugeABE(p1,mpg,mqa1,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau1(1,N),w0(1,N),atau2(1,N),-1,scale,
     3     rgvga1a2(1,3),rgbga1a2(1,3), wardidtest17,GaugeAbe)
c g A2 W A1
       call HexGaugeABE(p1,mpg,mqa2,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     Myeps(0,1),atau2(1,N),w0(1,N),atau1(1,N),-1,scale,
     3     rgvga2a1(1,3),rgbga2a1(1,3), wardidtest18,GaugeAbe)
c A1 g W A2
       call HexGaugeABE(p1,mqa1,mpg,mpw,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),Myeps(0,1),w0(1,N),atau2(1,N),-1,scale,
     3     rgva1ga2(1,3),rgba1ga2(1,3), wardidtest19,GaugeAbe)
c A2 g W A1
       call HexGaugeABE(p1,mqa2,mpg,mpw,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),Myeps(0,1),w0(1,N),atau1(1,N),-1,scale,
     3     rgva2ga1(1,3),rgba2ga1(1,3), wardidtest20,GaugeAbe)
c A1 W g A2
       call HexGaugeABE(p1,mqa1,mpw,mpg,mqa2,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),Myeps(0,1),atau2(1,N),-1,scale,
     3     rgva1ga2(1,2),rgba1ga2(1,2), wardidtest21,GaugeAbe)
c A2 W g A1
        call HexGaugeABE(p1,mqa2,mpw,mpg,mqa1,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),Myeps(0,1),atau1(1,N),-1,scale,
     3     rgva2ga1(1,2),rgba2ga1(1,2), wardidtest22,GaugeAbe)

c A1 W A2 g
      call HexGaugeABE(p1,mqa1,mpw,mqa2,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2g(1,2),rgba1a2g(1,2), wardidtest23,GaugeAbe)

c A2 W A1 g
       call HexGaugeABE(p1,mqa2,mpw,mqa1,mpg,p2,
     1     psi(1,-1,2),psi(1,-1,1),
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1g(1,2),rgba2a1g(1,2), wardidtest24,GaugeAbe)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Hexagon NoAbelian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Accuracyres=0.1d0
C A1 A2 W
	call HexGauge1gNoAbe(accuracyres,p1,mqa1,mqa2,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),atau2(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2n(1,3),wardidtest25,GaugeNoAbe)
C A2 A1 W
	call HexGauge1gNoAbe(accuracyres,p1,mqa2,mqa1,mpw,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),atau1(1,N),w0(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1n(1,3),wardidtest26,GaugeNoAbe)
CCCCCCCCCCCCCCCCCC
C A1 W A2
	call HexGauge1gNoAbe(accuracyres,p1,mqa1,mpw,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau1(1,N),w0(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2n(1,2),wardidtest27,GaugeNoAbe)
C A2 W A1
	call HexGauge1gNoAbe(accuracyres,p1,mqa2,mpw,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     atau2(1,N),w0(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1n(1,2),wardidtest28,GaugeNoAbe)
CCCCCCCCCCCCCCC
C W A1 A2
	call HexGauge1gNoAbe(accuracyres,p1,mpw,mqa1,mqa2,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau1(1,N),atau2(1,N),Myeps(0,1),-1,scale,
     3     rgva1a2n(1,1),wardidtest29,GaugeNoAbe)
C W A2 A1 
	call HexGauge1gNoAbe(accuracyres,p1,mpw,mqa2,mqa1,p2,mpg,
     1     psi(1,-1,2),psi(1,-1,1),	
     2     w0(1,N),atau2(1,N),atau1(1,N),Myeps(0,1),-1,scale,
     3     rgva2a1n(1,1),wardidtest30,GaugeNoAbe)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C***
C*****************************************************************
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6
!      Print*, 'wardidtest1',wardidtest1
       wardidtest1=wardidtest1.and.wardidtest7.and.wardidtest8.and.
     1 wardidtest9.and.wardidtest10.and.wardidtest11
!       Print*, 'wardidtest2',wardidtest1
        wardidtest1=wardidtest1.and.wardidtest12.and.wardidtest13.and.
     1 wardidtest14.and.wardidtest15.and.wardidtest16
!        Print*, 'wardidtest3',wardidtest1
       wardidtest1=wardidtest1.and.wardidtest17.and.wardidtest18.and.
     1 wardidtest19.and.wardidtest20.and.wardidtest21
!        Print*, 'wardidtest4',wardidtest1     
       wardidtest1=wardidtest1.and.wardidtest22.and.wardidtest23.and.
     1 wardidtest24
!        Print*, 'wardidtest5',wardidtest1    
        wardidtest1=wardidtest1.and.wardidtest25.and.wardidtest26.and.
     1 wardidtest27
        wardidtest1=wardidtest1.and.wardidtest28.and.wardidtest29.and.
     1 wardidtest30
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC1CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!	Print*, 'TOTAL wardidtest1',wardidtest1 
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	    if (.not.wardidtest1) then
	      do l=1,2
               do i = 1,4
 	         do mu=1,2
                    rva1a2g(mu,i,l) = (0.0d0, 0.0d0)
                    rva2a1g(mu,i,l) = (0.0d0, 0.0d0)
         	    rva1ga2(mu,i,l) = (0.0d0, 0.0d0)
                    rva2ga1(mu,i,l) = (0.0d0, 0.0d0)
  		    rvga1a2(mu,i,l) = (0.0d0, 0.0d0)
                    rvga2a1(mu,i,l) = (0.0d0, 0.0d0)
 		  enddo   
               enddo
	      enddo
	     do l=1,2 
  	      do i = 1,3
	        do mu=1,2
  	           rva1a2n(mu,i,l) = (0.0d0, 0.0d0)
  	           rva2a1n(mu,i,l) = (0.0d0, 0.0d0)
	         enddo
  	      enddo
	     enddo 
              hexcountm = hexcountm+1
!              Print*,'BAD hexcountm',hexcountm
              hexcountps(ps_number) = hexcountps(ps_number)+1 
!              Print*,'hexcountps(ps_number) ',ps_number,hexcountps(ps_number)
             else
	     do l=1,2
              do i = 1,4
	        rva1a2g(1,i,l) = (cfc*rva1a2g(1,i,l)+cfcac*rva1a2g(2,i,l))
                rva2a1g(1,i,l) = (cfc*rva2a1g(1,i,l)+cfcac*rva2a1g(2,i,l))
 		rva1ga2(1,i,l) = (cfc*rva1ga2(1,i,l)+cfcac*rva1ga2(2,i,l))
                rva2ga1(1,i,l) = (cfc*rva2ga1(1,i,l)+cfcac*rva2ga1(2,i,l))
 		rvga1a2(1,i,l) = (cfc*rvga1a2(1,i,l)+cfcac*rvga1a2(2,i,l))
                rvga2a1(1,i,l) = (cfc*rvga2a1(1,i,l)+cfcac*rvga2a1(2,i,l))
               enddo
	      enddo
	     do l=1,2 
 	      do i = 1,3
 	      rva1a2n(1,i,l) = cacgraph*rva1a2n(1,i,l)
 	      rva2a1n(1,i,l) = cacgraph*rva2a1n(1,i,l)
 	      enddo
	     enddo 
                hexcount2m = hexcount2m+1
              hexcount2ps(ps_number) = hexcount2ps(ps_number)+1   
!              Print*,'GOOD hexcount2m',hexcount2m
!              Print*,'hexcount2ps(ps_number) ',ps_number,hexcount2ps(ps_number)
             endif

      !Print*, 'rvga1a2(1,1,1)',rvga1a2(1,1,1)

        do l=1,2
c   q qbar -> A A W
!Check Born
!gA1A2W+gA2A1W+A1gA2W+A2gA1W+A1A2gW+A2A1gW+A1A2Wg+A2A1Wg
        mv1b =clr(u(k),1,-1)**2*(
     1  rbga1a2(4,l)+rbga2a1(4,l)+rba1ga2(4,l)+rba2ga1(4,l)+
     2  rba1a2g(4,l)+rba2a1g(4,l)+rba1a2g(3,l)+rba2a1g(3,l)) 
!gA1WA2+gA2WA1+A1gWA2+A2gWA1+A1WgA2+A2WgA1+A1WA2g+A2WA1g    
	mv2b=clr(u(k),1,-1)*clr(d(k),1,-1)*(
     1  rbga1a2(3,l)+rbga2a1(3,l)+rba1ga2(3,l)+rba2ga1(3,l)+
     2  rba1ga2(2,l)+rba2ga1(2,l)+rba1a2g(2,l)+rba2a1g(2,l) )
!gWA1A2+gWA2A1+WgA1A2+WgA2A1+WA1gA2+WA2gA1+WA1A2g+WA2A1g
	mv3b= clr(d(k),1,-1)**2*(
     1  rbga1a2(2,l)+rbga2a1(2,l)+rbga1a2(1,l)+rbga2a1(1,l)+
     2  rba1ga2(1,l)+rba2ga1(1,l)+rba1a2g(1,l)+rba2a1g(1,l))
 
        matb(k,-1,l,1)=(mv1b+mv2b+mv3b)* clr(3,3,-1)
 
       ratioComp=abs(mat(k,-1,l,1)/matb(k,-1,l,1)-1d0)
      
ccc       if(PerformGauge.eq.0) then
ccc       if(ratioComp.gt.1d-6) then
ccc       Print*, ""
ccc       Print*,'l',l
ccc       Print*, 'wardidtest1',wardidtest1
ccc       Print*, "Born/bor_hex",	mat(k,-1,l,1)/(matb(k,-1,l,1))
ccc       Print*, "Born",mat(k,-1,l,1)
ccc       !stop
ccc       endif 
ccc       endif
C****************************************************************      
C*****************************************************************      
c   q qbar -> A A W
! Virtual: Here I just do the replacement rbg...->rvg..
! From the checked mv1 against born.
C*****************************************************************
C*****************************************************************
!gA1A2W+gA2A1W+A1gA2W+A2gA1W+A1A2gW+A2A1gW+A1A2Wg+A2A1Wg
        mv1 =clr(u(k),1,-1)**2*(
     1  rvga1a2(1,4,l)+rvga2a1(1,4,l)+rva1ga2(1,4,l)+rva2ga1(1,4,l)+
     2  rva1a2g(1,4,l)+rva2a1g(1,4,l)+rva1a2g(1,3,l)+rva2a1g(1,3,l)) 
!gA1WA2+gA2WA1+A1gWA2+A2gWA1+A1WgA2+A2WgA1+A1WA2g+A2WA1g    
	mv2=clr(u(k),1,-1)*clr(d(k),1,-1)*(
     1  rvga1a2(1,3,l)+rvga2a1(1,3,l)+rva1ga2(1,3,l)+rva2ga1(1,3,l)+
     2  rva1ga2(1,2,l)+rva2ga1(1,2,l)+rva1a2g(1,2,l)+rva2a1g(1,2,l) )
!gWA1A2+gWA2A1+WgA1A2+WgA2A1+WA1gA2+WA2gA1+WA1A2g+WA2A1g
	mv3= clr(d(k),1,-1)**2*(
     1  rvga1a2(1,2,l)+rvga2a1(1,2,l)+rvga1a2(1,1,l)+rvga2a1(1,1,l)+
     2  rva1ga2(1,1,l)+rva2ga1(1,1,l)+rva1a2g(1,1,l)+rva2a1g(1,1,l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non-Abelian graphs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a1a2W+a2a1W	  
	  mv4=clr(u(k),1,-1)**2*(rva1a2n(1,3,l)+rva2a1n(1,3,l))
! a1Wa2 +a2Wa1
          mv5=clr(u(k),1,-1)*clr(d(k),1,-1)*(rva1a2n(1,2,l)+rva2a1n(1,2,l))
! Wa1a2 + Wa2a1
	  mv6=clr(d(k),1,-1)**2*(rva1a2n(1,1,l)+rva2a1n(1,1,l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
	IF (wardidtest1 ) THEN
	
	IF (PerformGauge.eq.1) then
! gauge test for the gluon. The CF part of mv1,mv2,mv3,
! is already zero. mv1,mv2,mv3  only contains CA part which
! should cancel agains NOABE, mv4,mv5,mv6

!  This test checks the couplings and the NOABe which can 
!  not be check against born
	
        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	Print*, 'mv1b',mv1b
	Print*, 'mv1',mv1
	Print*, 'mv4',mv4
	Print*, 'Gauge',-mv4/mv1-1d0
!	Print*, 'rva1a2n(1,3,l)',rva1a2n(1,3,1)
	
        Print*, 'mv2',mv2
	Print*, 'mv5',mv5
	Print*, 'Gauge',-mv2/mv5-1d0
	
	Print*, 'mv3',mv3
	Print*, 'mv6',mv6
	Print*, 'Gauge',-mv3/mv6-1d0
	!stop
        elseif(PerformGauge.gt.1) then
!   Here only the sum of all the contributions is gauge
!   invariant, but, with coupling cte equal=1. 
!    Otherwise, there are reamining terms that cancel against
!    penline's and so on due to the difference of clr(u(k),1,-1)
!    and clr(d(k),1,-1)
       Print*,'Here only the sum of all the contributions is gauge'
       Print*,'invariant, but, with coupling cte equal=1.'
       Print*,'This is corrected with the fac cte '
       Print*,'Otherwise, there are reamining terms that cancel against'
       Print*,'penlines and so on due to the difference in
     -    clr(u(k),1,-1)'
       Print*,'and clr(d(k),1,-1)'


        Print*, ''
	Print*, 'wardidtest1',wardidtest1 
	Print*, 'cacgraph',cacgraph
	Print*, 'cfc',cfc
	Print*, 'cfca',cfcac
	Print*, 'cfca-CF',cfcac-cfc
	Print*, ''
	
	Print*, 'mv1b',mv1b
	Print*, 'mv2b',mv2b
	Print*, 'mv3b',mv3b
	Print*, 'mv1b+mv2b+mv3b',mv1b+mv2b+mv3b
		
	Print*, ''
	
	Print*, 'mv1               ',mv1
	Print*, 'fac*mv2+fac**2*mv3',clr(u(k),1,-1)/clr(d(k),1,-1)*(mv2+ 
     2   clr(u(k),1,-1)/clr(d(k),1,-1)*mv3)
	Print*, 'ratio             ', -mv1/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mv2+ 
     2   clr(u(k),1,-1)/clr(d(k),1,-1)*mv3))
	
	Print*, 'mv4               ',mv4
	
	Print*, 'fac*mv5+fac**2*mv6',clr(u(k),1,-1)/clr(d(k),1,-1)*(mv5+
     2	clr(u(k),1,-1)/clr(d(k),1,-1)*mv6)
	Print*, 'ratio             ',-mv4/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mv5+
     2	clr(u(k),1,-1)/clr(d(k),1,-1)*mv6))


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
        mv1CA =clr(u(k),1,-1)**2*(
     1  rvga1a2(2,4,l)+rvga2a1(2,4,l)+rva1ga2(2,4,l)+rva2ga1(2,4,l)+
     2  rva1a2g(2,4,l)+rva2a1g(2,4,l)+rva1a2g(2,3,l)+rva2a1g(2,3,l)) 
!gA1WA2+gA2WA1+A1gWA2+A2gWA1+A1WgA2+A2WgA1+A1WA2g+A2WA1g    
	mv2CA=clr(u(k),1,-1)*clr(d(k),1,-1)*(
     1  rvga1a2(2,3,l)+rvga2a1(2,3,l)+rva1ga2(2,3,l)+rva2ga1(2,3,l)+
     2  rva1ga2(2,2,l)+rva2ga1(2,2,l)+rva1a2g(2,2,l)+rva2a1g(2,2,l) )
!gWA1A2+gWA2A1+WgA1A2+WgA2A1+WA1gA2+WA2gA1+WA1A2g+WA2A1g
	mv3CA= clr(d(k),1,-1)**2*(
     1  rvga1a2(2,2,l)+rvga2a1(2,2,l)+rvga1a2(2,1,l)+rvga2a1(2,1,l)+
     2  rva1ga2(2,1,l)+rva2ga1(2,1,l)+rva1a2g(2,1,l)+rva2a1g(2,1,l))

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
	 
	 ratioFact=mv3/(-2d0*mv3b)
	 Print*, 'mv3',mv3
	 Print*, 'mv3b',mv3b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 Print*, 'SUM of born', mv1b+mv2b+mv3b
! CA part
	 Print*, ""
	 Print*, ""
	 Print*, "CA PART"
	 Print*, ""
	 Print*, ""
	 Print*, 'mv1CA',mv1CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv1b',mv1b
	 ratioFact=(mv1CA*(-1d0/2d0)+mv4)/(-1d0*mv1b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""
	 
         Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv5
         Print*, 'mv2b',mv2b
	 ratioFact=(mv2CA*(-1d0/2d0)+mv5)/(-1d0*mv2b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, "" 
	 
         Print*, 'mv3CA',mv3CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv6
         Print*, 'mv3b',mv3b
	 ratioFact=(mv3CA*(-1d0/2d0)+mv6)/(-1d0*mv3b)
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
	 Print*, 'mv1CA',mv1CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv1b*cte',mv1b*cteEps1
	 ratioFact=(mv1CA*(-1d0/2d0)+mv4)/(mv1b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	 Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv5
         Print*, 'mv2b*cte',mv2b*cteEps1
	 ratioFact=(mv2CA*(-1d0/2d0)+mv5)/(mv2b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""
	  
	 Print*, 'mv3CA',mv3CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv6
         Print*, 'mv3b*cte',mv3b*cteEps1
	 ratioFact=(mv3CA*(-1d0/2d0)+mv6)/(mv3b*cteEps1)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C**************************************************************
C   SUM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      matv(k,-1,l,1) = (mv1+mv2+mv3+mv4+mv5+mv6)*clr(3,3,-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PerformGauge.eq.1) then
! then, matv(k,-1,l,1) should be zero
      Print*,''
      Print*,'matv(k,-1,l,1) should be zero'
      Print*,' matv(k,-1,l,1)', matv(k,-1,l,1)
      Print*, 'mv1*clr(3,3,-1)',mv1*clr(3,3,-1)
      Print*,''
      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
      Print*,''
      Print*,'matv(k,-1,l,1) with couplings corrected should be zero'
      Print*,' matv(k,-1,l,1)', matv(k,-1,l,1)
      cteEps1=clr(u(k),1,-1)/clr(d(k),1,-1) 
      matv(k,-1,l,1) = (mv1+mv4+(mv2+mv5)*cteEps1+(mv3+mv6)*cteEps1**2)*clr(3,3,-1)
      Print*,'matv(k,-1,l,1) with corrected couplings'
      Print*,'matv(k,-1,l,1)', matv(k,-1,l,1)
      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,1)',matv(k,-1,l,1)
      Print*, '-4*matvb(k,-1,l,1)',-4*matb(k,-1,l,1)
      Print*, 'ratio',matv(k,-1,l,1)/(-4*matb(k,-1,l,1))
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
!!!! Here I add the convfact terms. It comes from the expansion of the 1/eps terms in 
!  virtual contributions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      matv(k,-1,l,1)=matv(k,-1,l,1)+matb(k,-1,l,1)*convfact
! Here I substract k_cont. It is a guess to make the contribution small. 
! Not clear if it works
      matv(k,-1,l,1)=matv(k,-1,l,1) -matb(k,-1,l,1)*(K_cont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo ! l gluon polarization
      !If(wardidtest1) stop
!  Check renormalization invariance using K_contcheck      
! This test actually check the finite result of the hexagon contribution!!!!
! The finite terms are musq dependent due to the definition of C0,D0,B0.
! however, I know that this depencence is a term proportional to the born.
! Therefore, we can chech that M_v(musq) + f(musq)= M_v(musq1)+f(musq1)
      if(Reno_inv.eq.1) then
      check(1)=matv(k,-1,1,1)+(K_contcheck+K_contcheck4)*matb(k,-1,1,1)
      check(2)=matv(k,-1,2,1)+(K_contcheck+K_contcheck4)*matb(k,-1,2,1)
      check1(1)=matv(k,-1,1,1)
      check1(2)=matv(k,-1,2,1)
      Reno_inv=2
      musq=1d0
      goto 100
      endif
      If (Reno_inv.gt.0) then
       Print*, ' '
       Print*, 'chek1 hel 1 ',check1(1)
       Print*, 'check2 hel -1',check1(2)
       Print*, 'chek1 hel 1 ',check(1)
       Print*, 'check2 hel -1',check(2)
       Print*, ''
       Print*, 'check1 hel 1 ',matv(k,-1,1,1)
       Print*, 'check2 hel -1',matv(k,-1,2,1)
       Print*, 'check1 hel 1 ',(matv(k,-1,1,1)+(K_contcheck+K_contcheck4)*matb(k,-1,1,1))
       Print*, 'check2 hel -1',(matv(k,-1,2,1)+(K_contcheck+K_contcheck4)*matb(k,-1,2,1))
      
       Print*, 'ratio1',check(1)/((matv(k,-1,1,1)+(K_contcheck+K_contcheck4)*matb(k,-1,1,1)))-1d0
       Print*, 'ratio2',check(2)/((matv(k,-1,2,1)+(K_contcheck+K_contcheck4)*matb(k,-1,2,1)))-1d0
       stop
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
      endif ! Finish lhex
CCCCCCCCCCCCCCCCCCCCCCCCC
C Finish with Hexagons
CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC
C Pentagons
CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC
C the results are store in rvga2(x,y,l)
C The name indicate the order of gluon and the photons
C x: Number of graphs, 
C x=1 collect results proportional CF
C x=2 collect results proportional to CF-1/2CA
C y: Indicates the position of the W. This will affect 
C    the couplings
C l: Polarization of the gluonl l=1,l=2
CCCCCCCCCCCCCCCCCCCCCCCCC
C qq-> A W(W->WA) g
CCCCCCCCCCCCCCCCCCCCCCCCC
       If(lpent) then
      scale= musq!mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))
      
            do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
	       mpg(mu) = -q(mu)
               mqa1(mu) = -qa1(mu,1)
               mqa2(mu) = -qa2(mu,1)
               mqwa2(mu) = -p1(mu)-p2(mu)+qa1(mu,1)+q(mu)
	       mqwa1(mu) = -p1(mu)-p2(mu)+qa2(mu,1)+q(mu)
	       cmqwa2(mu) = mqwa2(mu)
	       cmqwa1(mu) = mqwa1(mu)	       
	   enddo
	      
      if(bos.eq.311) then
      do mu=0,3
      wwa2(mu,1)=wpwa2(mu,1)
      wwa1(mu,1)=wpwa1(mu,1)
      enddo
      elseif(bos.eq.411) then
      do mu=0,3
      wwa2(mu,1)=wmwa2(mu,1)
      wwa1(mu,1)=wmwa1(mu,1)
      enddo
      endif
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
	     mutemp=mu
	    wwa1(mutemp,N) = mqwa1(mu)
	    wwa2(mutemp,N) = mqwa2(mu)
	 enddo
      elseif(PerformGauge.eq.3) then
         do mu = 0,3
	    mutemp=mu+1
	    atau2(mutemp,N) = mqa2(mu)
	    atau1(mutemp,N) = mqa1(mu)
	 enddo
       endif
CCCCCCCCCCCCCCCCCCCCCCCCC
C  Abelian PART
CCCCCCCCCCCCCCCCCCCCCCCC
C
C  A1 g W
C
       call penlineABETotal(p1,mqa1,mpg,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1   atau1(1,1),Myeps(0,1),wwa2(0,1),-1,scale,1,3,
     2   4,1,pentcontr1,j5a1g(1,3,1),pentcontr1b,jba1g(3,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mqa1,mpg,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1   atau1(1,1),Myeps(0,1),cmqwa2,-1,scale,1,3,
     2   9,-1,pentemp,j5ma1g(1,3,1),pentemp,jbma1g(3,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time   
       call penlineABETotal(p1,mqa1,mpg,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1   atau1(1,1),Myeps(0,2),wwa2(0,1),-1,scale,1,3,
     2   9,-1,pentemp,j5a1g(1,3,2),pentemp,jba1g(3,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mqa1,mpg,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1   atau1(1,1),Myeps(0,2),cmqwa2,-1,scale,1,3,
     2   9,-1,pentemp,j5ma1g(1,3,2),pentemp,jbma1g(3,2),Div)
C
C  g A1 W
C
       call penlineABETotal(p1,mpg,mqa1,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau1(1,1),wwa2(0,1),-1,scale,1,2,
     2    4,1,pentcontr2,j5ga1(1,3,1),pentcontr2b,jbga1(3,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqa1,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau1(1,1),cmqwa2,-1,scale,1,2,
     2    9,-1,pentemp,j5mga1(1,3,1),pentemp,jbmga1(3,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
       call penlineABETotal(p1,mpg,mqa1,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),atau1(1,1),wwa2(0,1),-1,scale,1,2,
     2    9,-1,pentemp,j5ga1(1,3,2),pentemp,jbga1(3,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqa1,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),atau1(1,1),cmqwa2,-1,scale,1,2,
     2    9,-1,pentemp,j5mga1(1,3,2),pentemp,jbmga1(3,2),Div)     
C
C A2 g W
C
      call penlineABETotal(p1,mqa2,mpg,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,1),wwa1(0,1),-1,scale,1,3,
     2    4,1,pentcontr3,j5a2g(1,3,1),pentcontr3b,jba2g(3,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqa2,mpg,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,1),cmqwa1,-1,scale,1,3,
     2    9,-1,pentemp,j5ma2g(1,3,1),pentemp,jbma2g(3,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineABETotal(p1,mqa2,mpg,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,2),wwa1(0,1),-1,scale,1,3,
     2    9,-1,pentemp,j5a2g(1,3,2),pentemp,jba2g(3,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqa2,mpg,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,2),cmqwa1,-1,scale,1,3,
     2    9,-1,pentemp,j5ma2g(1,3,2),pentemp,jbma2g(3,2),Div)          
C
C g A2 W
C
       call penlineABETotal(p1,mpg,mqa2,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau2(1,1),wwa1(0,1),-1,scale,1,2,
     2    4,1,pentcontr4,j5ga2(1,3,1),pentcontr4b,jbga2(3,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqa2,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau2(1,1),cmqwa1,-1,scale,1,2,
     2    9,-1,pentemp,j5mga2(1,3,1),pentemp,jbmga2(3,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
       call penlineABETotal(p1,mpg,mqa2,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),atau2(1,1),wwa1(0,1),-1,scale,1,2,
     2    9,-1,pentemp,j5ga2(1,3,2),pentemp,jbga2(3,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqa2,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),atau2(1,1),cmqwa1,-1,scale,1,2,
     2    9,-1,pentemp,j5mga2(1,3,2),pentemp,jbmga2(3,2),Div)     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C A1 W g
C
       call penlineABETotal(p1,mqa1,mqwa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,1),-1,scale,1,4,
     2    4,1,pentcontr5,j5a1g(1,2,1),pentcontr5b,jba1g(2,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mqa1,mqwa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),cmqwa2,Myeps(0,1),-1,scale,1,4,
     2    9,-1,pentemp,j5ma1g(1,2,1),pentemp,jbma1g(2,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
       call penlineABETotal(p1,mqa1,mqwa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5a1g(1,2,2),pentemp,jba1g(2,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mqa1,mqwa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),cmqwa2,Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5ma1g(1,2,2),pentemp,jbma1g(2,2),Div)          
C
C g W A1
C
       call penlineABETotal(p1,mpg,mqwa2,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwa2(0,1),atau1(1,1),-1,scale,1,2,
     2    4,1,pentcontr6,j5ga1(1,2,1),pentcontr6b,jbga1(2,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqwa2,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),cmqwa2,atau1(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5mga1(1,2,1),pentemp,jbmga1(2,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
       call penlineABETotal(p1,mpg,mqwa2,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),wwa2(0,1),atau1(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5ga1(1,2,2),pentemp,jbga1(2,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqwa2,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),cmqwa2,atau1(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5mga1(1,2,2),pentemp,jbmga1(2,2),Div)          
C
C A2 W g
C
      call penlineABETotal(p1,mqa2,mqwa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,1),-1,scale,1,4,
     2    4,1,pentcontr7,j5a2g(1,2,1),pentcontr7b,jba2g(2,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqa2,mqwa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),cmqwa1,Myeps(0,1),-1,scale,1,4,
     2    9,-1,pentemp,j5ma2g(1,2,1),pentemp,jbma2g(2,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time   
      call penlineABETotal(p1,mqa2,mqwa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5a2g(1,2,2),pentemp,jba2g(2,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqa2,mqwa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),cmqwa1,Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5ma2g(1,2,2),pentemp,jbma2g(2,2),Div)          
C
C g W A2
C
       call penlineABETotal(p1,mpg,mqwa1,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwa1(0,1),atau2(1,1),-1,scale,1,2,
     2    4,1,pentcontr8,j5ga2(1,2,1),pentcontr8b,jbga2(2,1),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqwa1,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),cmqwa1,atau2(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5mga2(1,2,1),pentemp,jbmga2(2,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
       call penlineABETotal(p1,mpg,mqwa1,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),wwa1(0,1),atau2(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5ga2(1,2,2),pentemp,jbga2(2,2),Div)
c Additional contributions from the propatos kk/M     
       call penlineABETotal(p1,mpg,mqwa1,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),cmqwa1,atau2(1,1),-1,scale,1,2,
     2    9,-1,pentemp,j5mga2(1,2,2),pentemp,jbmga2(2,2),Div)          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A1 g
C
      call penlineABETotal(p1,mqwa2,mqa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,1),-1,scale,1,4,
     2    4,1,pentcontr9,j5a1g(1,1,1),pentcontr9b,jba1g(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa2,mqa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,atau1(1,1),Myeps(0,1),-1,scale,1,4,
     2    9,-1,pentemp,j5ma1g(1,1,1),pentemp,jbma1g(1,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineABETotal(p1,mqwa2,mqa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5a1g(1,1,2),pentemp,jba1g(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa2,mqa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,atau1(1,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5ma1g(1,1,2),pentemp,jbma1g(1,2),Div)          
C
C W g A1
C
      call penlineABETotal(p1,mqwa2,mpg,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),Myeps(0,1),atau1(1,1),-1,scale,1,3,
     2    4,1,pentcontr10,j5ga1(1,1,1),pentcontr10b,jbga1(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa2,mpg,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,Myeps(0,1),atau1(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5mga1(1,1,1),pentemp,jbmga1(1,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineABETotal(p1,mqwa2,mpg,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),Myeps(0,2),atau1(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5ga1(1,1,2),pentemp,jbga1(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa2,mpg,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,Myeps(0,2),atau1(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5mga1(1,1,2),pentemp,jbmga1(1,2),Div)          
C
C W A2 g
C
      call penlineABETotal(p1,mqwa1,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,1),-1,scale,1,4,
     2    4,1,pentcontr11,j5a2g(1,1,1),pentcontr11b,jba2g(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa1,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,atau2(1,1),Myeps(0,1),-1,scale,1,4,
     2    9,-1,pentemp,j5ma2g(1,1,1),pentemp,jbma2g(1,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineABETotal(p1,mqwa1,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5a2g(1,1,2),pentemp,jba2g(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa1,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,atau2(1,1),Myeps(0,2),-1,scale,1,4,
     2    9,-1,pentemp,j5ma2g(1,1,2),pentemp,jbma2g(1,2),Div)          
C
C W g A2
C
      call penlineABETotal(p1,mqwa1,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),Myeps(0,1),atau2(1,1),-1,scale,1,3,
     2    4,1,pentcontr12,j5ga2(1,1,1),pentcontr12b,jbga2(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa1,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,Myeps(0,1),atau2(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5mga2(1,1,1),pentemp,jbmga2(1,1),Div)
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineABETotal(p1,mqwa1,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),Myeps(0,2),atau2(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5ga2(1,1,2),pentemp,jbga2(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineABETotal(p1,mqwa1,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,Myeps(0,2),atau2(1,1),-1,scale,1,3,
     2    9,-1,pentemp,j5mga2(1,1,2),pentemp,jbmga2(1,2),Div)          
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C  NO ABELIAN PART
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C A1 W
C
      call penlineNoAbeTotal(p1,mqa1,mqwa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,1),-1,scale,4,1,
     2    pentcontrNoAbe1,jg5a1(2,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqa1,mqwa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),cmqwa2,Myeps(0,1),-1,scale,9,-1,
     2    pentemp,jg5ma1(2,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineNoAbeTotal(p1,mqa1,mqwa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5a1(2,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqa1,mqwa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),cmqwa2,Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5ma1(2,2),Div)          
C
C A2 W
C
      call penlineNoAbeTotal(p1,mqa2,mqwa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,1),-1,scale,4,1,
     2    pentcontrNoAbe2,jg5a2(2,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqa2,mqwa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),cmqwa1,Myeps(0,1),-1,scale,9,-1,
     2    pentemp,jg5ma2(2,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineNoAbeTotal(p1,mqa2,mqwa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5a2(2,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqa2,mqwa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),cmqwa1,Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5ma2(2,2),Div)          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A1
C
      call penlineNoAbeTotal(p1,mqwa2,mqa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,1),-1,scale,4,1,
     2    pentcontrNoAbe3,jg5a1(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqwa2,mqa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,atau1(1,1),Myeps(0,1),-1,scale,9,-1,
     2    pentemp,jg5ma1(1,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineNoAbeTotal(p1,mqwa2,mqa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5a1(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqwa2,mqa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa2,atau1(1,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5ma1(1,2),Div)          
C
C W A2
C
      call penlineNoAbeTotal(p1,mqwa1,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,1),-1,scale,4,1,
     2    pentcontrNoAbe4,jg5a2(1,1),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqwa1,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,atau2(1,1),Myeps(0,1),-1,scale,9,-1,
     2    pentemp,jg5ma2(1,1),Div)     
c  gluon helicity 2 . Comp=-1. Reduce 1/2 the CPU time        
      call penlineNoAbeTotal(p1,mqwa1,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5a2(1,2),Div)
c Additional contributions from the propatos kk/M     
      call penlineNoAbeTotal(p1,mqwa1,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    cmqwa1,atau2(1,1),Myeps(0,2),-1,scale,9,-1,
     2    pentemp,jg5ma2(1,2),Div)          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Gauge Test
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  A1 g W
       call pent_gauge(p1,mqa1,mpg,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),Myeps(0,1),wwa2(0,1),-1,scale,
     2    pentcontr1,pentcontr1b,wardidtest1)
C  g A1 W
       call pent_gauge(p1,mpg,mqa1,mqwa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau1(1,1),wwa2(0,1),-1,scale,
     2    pentcontr2,pentcontr2b,wardidtest2)
C A2 g W
      call pent_gauge(p1,mqa2,mpg,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),Myeps(0,1),wwa1(0,1),-1,scale,
     2    pentcontr3,pentcontr3b,wardidtest3)
C g A2 W
       call pent_gauge(p1,mpg,mqa2,mqwa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),atau2(1,1),wwa1(0,1),-1,scale,
     2   pentcontr4,pentcontr4b,wardidtest4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A1 W g
       call pent_gauge(p1,mqa1,mqwa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,1),-1,scale,
     2   pentcontr5,pentcontr5b,wardidtest5)
C g W A1
       call pent_gauge(p1,mpg,mqwa2,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwa2(0,1),atau1(1,1),-1,scale,
     2   pentcontr6,pentcontr6b,wardidtest6)
C A2 W g
       call pent_gauge(p1,mqa2,mqwa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,1),-1,scale,
     2    pentcontr7,pentcontr7b,wardidtest7)
C g W A2
      call pent_gauge(p1,mpg,mqwa1,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwa1(0,1),atau2(1,1),-1,scale,
     2    pentcontr8,pentcontr8b,wardidtest8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A1 g
      call pent_gauge(p1,mqwa2,mqa1,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,1),-1,scale,
     2    pentcontr9,pentcontr9b,wardidtest9)
C W g A1
      call pent_gauge(p1,mqwa2,mpg,mqa1,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),Myeps(0,1),atau1(1,1),-1,scale,
     2    pentcontr10,pentcontr10b,wardidtest10)
C W A2 g
      call pent_gauge(p1,mqwa1,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,1),-1,scale,
     2    pentcontr11,pentcontr11b,wardidtest11)
C W g A2
      call pent_gauge(p1,mqwa1,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),Myeps(0,1),atau2(1,1),-1,scale,
     2    pentcontr12,pentcontr12b,wardidtest12)
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  NO ABELIAN PART
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C A1 W
      call pentNoAbe_gauge(1d-1,p1,mqa1,mqwa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau1(1,1),wwa2(0,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe1,accuracyres,wardidtest13)
C A2 W
      call pentNoAbe_gauge(1d-1,p1,mqa2,mqwa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    atau2(1,1),wwa1(0,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe2,accuracyres,wardidtest14)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W A1
       call pentNoAbe_gauge(1d-1,p1,mqwa2,mqa1,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa2(0,1),atau1(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe3,accuracyres,wardidtest15)
C W A2
      call pentNoAbe_gauge(1d-1,p1,mqwa1,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa1(0,1),atau2(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe4,accuracyres,wardidtest16)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6
       wardidtest1=wardidtest1.and.wardidtest7.and.wardidtest8.and.
     1 wardidtest9.and.wardidtest10.and.wardidtest11
       wardidtest1=wardidtest1.and.wardidtest12.and.wardidtest13.and.
     1 wardidtest14.and.wardidtest15.and.wardidtest16

     
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
	    do  l=1,2
                do i = 1,3
		   do mu=1,2
                     j5a1g(mu,i,l) = (0.0d0, 0.0d0)
                     j5a2g(mu,i,l) = (0.0d0, 0.0d0)
 		     j5ga2(mu,i,l) = (0.0d0, 0.0d0)
                     j5ga1(mu,i,l) = (0.0d0, 0.0d0)
		     
		     j5ma1g(mu,i,l) = (0.0d0, 0.0d0)
                     j5ma2g(mu,i,l) = (0.0d0, 0.0d0)
 		     j5mga2(mu,i,l) = (0.0d0, 0.0d0)
                     j5mga1(mu,i,l) = (0.0d0, 0.0d0)
		   enddo
                enddo
	     	
 	          jg5a2(1,l)=(0.0d0, 0.0d0)
 		  jg5a1(1,l)=(0.0d0, 0.0d0)
 	          jg5a2(2,l)=(0.0d0, 0.0d0)
 		  jg5a1(2,l)=(0.0d0, 0.0d0)
		  
		  jg5ma2(1,l)=(0.0d0, 0.0d0)
 		  jg5ma1(1,l)=(0.0d0, 0.0d0)
 	          jg5ma2(2,l)=(0.0d0, 0.0d0)
 		  jg5ma1(2,l)=(0.0d0, 0.0d0)
		  
	    enddo	  
              pentcountm = pentcountm+1
              pentcountps(ps_number) = pentcountps(ps_number)+1 
         else
	    do l=1,2
              do i = 1,3
                 j5a1g(1,i,l)=(cfc*j5a1g(1,i,l)+cfcac*j5a1g(2,i,l))
                 j5a2g(1,i,l)=(cfc*j5a2g(1,i,l)+cfcac*j5a2g(2,i,l))
		 j5ga1(1,i,l)=(cfc*j5ga1(1,i,l)+cfcac*j5ga1(2,i,l))
                 j5ga2(1,i,l)=(cfc*j5ga2(1,i,l)+cfcac*j5ga2(2,i,l))
ccccc
		 j5ma1g(1,i,l)=(cfc*j5ma1g(1,i,l)+cfcac*j5ma1g(2,i,l))
                 j5ma2g(1,i,l)=(cfc*j5ma2g(1,i,l)+cfcac*j5ma2g(2,i,l))
		 j5mga1(1,i,l)=(cfc*j5mga1(1,i,l)+cfcac*j5mga1(2,i,l))
                 j5mga2(1,i,l)=(cfc*j5mga2(1,i,l)+cfcac*j5mga2(2,i,l))
	      enddo
	      do i=1,2
		 jg5a1(i,l) = cacgraph*jg5a1(i,l)
		 jg5a2(i,l) = cacgraph*jg5a2(i,l)
cccccc		 
		 jg5ma1(i,l) = cacgraph*jg5ma1(i,l)
		 jg5ma2(i,l) = cacgraph*jg5ma2(i,l)
	      enddo  
	    enddo  
	    
              pentcount2m = pentcount2m+1
              pentcount2ps(ps_number) = pentcount2ps(ps_number)+1      
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   q qbar -> A A W
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do l=1,2

! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
        mv1b=clr(u(k),1,-1)*(jba1g(3,l)+jbga1(3,l)+jba2g(3,l)
     1         +jbga2(3,l)+jba1g(2,l)+jba2g(2,l))
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
        mv2b=clr(d(k),1,-1)*(jba1g(1,l)+jbga1(1,l)+jba2g(1,l)
     1         +jbga2(1,l)+jbga1(2,l)+jbga2(2,l))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
        mmv1b=-clr(u(k),1,-1)*(
     1  (jbma1g(3,l)+jbmga1(3,l)+jbma1g(2,l))*dotrc(mqwa2,wwa2(0,1))+
     2  (jbma2g(3,l)+jbmga2(3,l)+jbma2g(2,l))*dotrc(mqwa1,wwa1(0,1))
     3         )*zm2i(3)
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
      mmv2b=-clr(d(k),1,-1)*(
     1 (jbma1g(1,l)+jbmga1(1,l)+jbmga1(2,l))*dotrc(mqwa2,wwa2(0,1))+
     2 (jbma2g(1,l)+jbmga2(1,l)+jbmga2(2,l))*dotrc(mqwa1,wwa1(0,1))
     3       )*zm2i(3)
    
       matb(k,-1,l,2)=(mv1b+mv2b+mmv1b+mmv2b)* clr(3,3,-1)
    
      ratioComp=abs(mat(k,-1,l,2)/matb(k,-1,l,2)-1d0)
      
      if(ratioComp.gt.1d-6) then
      Print*, 'wardidtest1_PEN', wardidtest1
      Print*,"Born/b_pen-1d0",ratioComp
      Print*, "Born",mat(k,-1,l,2)
      Print*, 'Sum',matb(k,-1,l,2)
      Print*, 'mv1',mv1b*clr(3,3,-1)
      Print*, 'mmv1',mmv1b*clr(3,3,-1)
      Print*, 'mv2',mv2b*clr(3,3,-1)
      Print*, 'mmv2',mmv2b*   clr(3,3,-1)   
      !stop
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

! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
        mv1=clr(u(k),1,-1)*(j5a1g(1,3,l)+j5ga1(1,3,l)+j5a2g(1,3,l)
     1         +j5ga2(1,3,l)+j5a1g(1,2,l)+j5a2g(1,2,l))
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
        mv2=clr(d(k),1,-1)*(j5a1g(1,1,l)+j5ga1(1,1,l)+j5a2g(1,1,l)
     1         +j5ga2(1,1,l)+j5ga1(1,2,l)+j5ga2(1,2,l))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
       mmv1=-clr(u(k),1,-1)*(
     1  (j5ma1g(1,3,l)+j5mga1(1,3,l)+j5ma1g(1,2,l))*dotrc(mqwa2,wwa2(0,1))+
     2  (j5ma2g(1,3,l)+j5mga2(1,3,l)+j5ma2g(1,2,l))*dotrc(mqwa1,wwa1(0,1))
     3         )*zm2i(3)
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
      mmv2=-clr(d(k),1,-1)*(
     1 (j5ma1g(1,1,l)+j5mga1(1,1,l)+j5mga1(1,2,l))*dotrc(mqwa2,wwa2(0,1))+
     2 (j5ma2g(1,1,l)+j5mga2(1,1,l)+j5mga2(1,2,l))*dotrc(mqwa1,wwa1(0,1))
     3       )*zm2i(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non-Abelian graphs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A1W + A2W
	mv3=clr(u(k),1,-1)*(jg5a1(2,l)+jg5a2(2,l))
! WA1 + WA2 
	mv4=clr(d(k),1,-1)*(jg5a1(1,l)+jg5a2(1,l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional KK terms	
! A1W + A2W
	mmv3=-clr(u(k),1,-1)*(jg5ma1(2,l)*dotrc(mqwa2,wwa2(0,1))
     1	                    +jg5ma2(2,l)*dotrc(mqwa1,wwa1(0,1)))*zm2i(3)
! WA1 + WA2 
	mmv4=-clr(d(k),1,-1)*(jg5ma1(1,l)*dotrc(mqwa2,wwa2(0,1))
     1	          +jg5ma2(1,l)*dotrc(mqwa1,wwa1(0,1)))*zm2i(3)
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
	Print*, 'mv1b',mv1b
	Print*, 'mv1',mv1
	Print*, 'mv3',mv3
	Print*, 'Gauge',-mv3/mv1-1d0
	Print*, 'mmv1',mmv1
	Print*, 'mmv3',mmv3
	Print*, 'Gauge',-mmv3/mmv1-1d0

	
        Print*, 'mv2',mv2
	Print*, 'mv4',mv4
	Print*, 'Gauge',-mv2/mv4-1d0
        Print*, 'mmv2',mmv2
	Print*, 'mmv4',mmv4
	Print*, 'Gauge',-mmv2/mmv4-1d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif(PerformGauge.gt.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
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
	
	Print*, 'mv1b',mv1b
	Print*, 'mv2b',mv2b
	
	Print*, 'mv1b+fac*mv2b',mv1b+clr(u(k),1,-1)/clr(d(k),1,-1)*mv2b
		
	Print*, ''
	
	Print*, 'mv1    ',mv1
	Print*, 'fac*mv2',clr(u(k),1,-1)/clr(d(k),1,-1)*(mv2)
	Print*, 'ratio  ', -mv1/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mv2))
	Print*, ''
	
	Print*, 'mv3    ',mv3
	Print*, 'fac*mv4',clr(u(k),1,-1)/clr(d(k),1,-1)*(mv4)
	Print*, 'ratio  ',-mv3/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mv4))	
	Print*, ''
	Print*, 'mmv1    ',mmv1
	Print*, 'fac*mmv2',clr(u(k),1,-1)/clr(d(k),1,-1)*(mmv2)
	Print*, 'ratio   ', -mmv1/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mmv2))
	Print*, ''
	Print*, 'mmv3    ',mmv3
	
	Print*, 'fac*mmv4',clr(u(k),1,-1)/clr(d(k),1,-1)*(mmv4)
	Print*, 'ratio   ',-mmv3/(clr(u(k),1,-1)/clr(d(k),1,-1)*(mmv4))	
	Print*, ''
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
! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
        mv1CA=clr(u(k),1,-1)*(j5a1g(2,3,l)+j5ga1(2,3,l)+j5a2g(2,3,l)
     1         +j5ga2(2,3,l)+j5a1g(2,2,l)+j5a2g(2,2,l))
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
        mv2CA=clr(d(k),1,-1)*(j5a1g(2,1,l)+j5ga1(2,1,l)+j5a2g(2,1,l)
     1         +j5ga2(2,1,l)+j5ga1(2,2,l)+j5ga2(2,2,l))

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! A1gW+ gA1W+ A2gW+ gA2W+ A1Wg+ A2Wg
        mmv1CA=-clr(u(k),1,-1)*(
     1  (j5ma1g(2,3,l)+j5mga1(2,3,l)+j5ma1g(2,2,l))*dotrc(mqwa2,wwa2(0,1))+
     2  (j5ma2g(2,3,l)+j5mga2(2,3,l)+j5ma2g(2,2,l))*dotrc(mqwa1,wwa1(0,1))
     3         )*zm2i(3)
! WA1g+ WgA1+ WA2g+ WgA2+ gWA1+ gWA2 
      mmv2CA=-clr(d(k),1,-1)*(
     1 (j5ma1g(2,1,l)+j5mga1(2,1,l)+j5mga1(2,2,l))*dotrc(mqwa2,wwa2(0,1))+
     2 (j5ma2g(2,1,l)+j5mga2(2,1,l)+j5mga2(2,2,l))*dotrc(mqwa1,wwa1(0,1))
     3       )*zm2i(3)
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
	 Print*, 'SUM of born', mv1b+mv2b

         ratioFact=mmv1/(-2d0*mmv1b)
	 Print*, 'mmv1',mmv1
	 Print*, 'mmv1b',mmv1b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 
	 ratioFact=mmv2/(-2d0*mmv2b)
	 Print*, 'mmv2',mmv2
	 Print*, 'mmv2b',mmv2b
	 Print*, "CF Div/Eps^2", ratioFact
	 Print*, ""
	 Print*, ""
	 Print*, 'SUM of born', mmv1b+mmv2b	 
	 
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
	 
	 Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv3
         Print*, 'mmv1b',mmv1b
	 ratioFact=(mmv1CA*(-1d0/2d0)+mmv3)/(-1d0*mmv1b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""
	 
         Print*, 'mmv2CA',mmv2CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv4
         Print*, 'mmv2b',mmv2b
	 ratioFact=(mmv2CA*(-1d0/2d0)+mmv4)/(-1d0*mmv2b)
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
	 
	 ratioFact=mmv1/(mmv1b*cteEps1)
	 Print*, 'mmv1',mmv1
	 Print*, 'mmv1b*cte',mmv1b*cteEps1
	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""

	 ratioFact=mmv2/(mmv2b*cteEps1)
	 Print*, 'mmv2',mmv2
	 Print*, 'mmv2b*cte',mmv2b*cteEps1
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

	 Print*, 'mv2CA',mv2CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv4
         Print*, 'mv2b*cte',mv2b*cteEps1
	 ratioFact=(mv2CA*(-1d0/2d0)+mv4)/(mv2b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	 Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv3
         Print*, 'mmv1b*cte',mmv1b*cteEps1
	 ratioFact=(mmv1CA*(-1d0/2d0)+mmv3)/(mmv1b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	 Print*, 'mmv2CA',mmv2CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv4
         Print*, 'mmv2b*cte',mmv2b*cteEps1
	 ratioFact=(mmv2CA*(-1d0/2d0)+mmv4)/(mmv2b*cteEps1)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          matv(k,-1,l,2) = (mv1+mv2+mmv1+mmv2+mv3+mv4+mmv3+mmv4)*clr(3,3,-1)
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
! then, matv(k,-1,l,1) should be zero
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(k,-1,l,2)', matv(k,-1,l,2)
      Print*, 'mv1*clr(3,3,-1)',mv1*clr(3,3,-1)
      Print*,''
      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
      Print*,''
      Print*,'matv(k,-1,l,2) with couplings corrected should be zero'
      Print*,' matv(k,-1,l,1)', matv(k,-1,l,2)
      cteEps1=clr(u(k),1,-1)/clr(d(k),1,-1) 
      matv(k,-1,l,2) = (mv1+mv3+mmv1+mmv3+(mv2+mv4+mmv2+mmv4)*cteEps1)*clr(3,3,-1)
      Print*,'matv(k,-1,l,2) with corrected couplings'
      Print*,'matv(k,-1,l,2)', matv(k,-1,l,2)
      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,1)',matv(k,-1,l,2)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(k,-1,l,2)
      Print*, 'ratio',matv(k,-1,l,2)/(-4*matb(k,-1,l,2))
      endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!              Add and substract cte's
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      	  matv(k,-1,l,2)=matv(k,-1,l,2)+matb(k,-1,l,2)*convfact
	  matv(k,-1,l,2)=matv(k,-1,l,2)-matb(k,-1,l,2)*(K_cont)
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
	enddo ! end l loop
	
!  Check renormalization invariance using K_contcheck      
      if(Reno_inv.eq.1) then
      check(1)=matv(k,-1,1,2)+(K_contcheck+K_contcheck4)*matb(k,-1,1,2)
      check(2)=matv(k,-1,2,2)+(K_contcheck+K_contcheck4)*matb(k,-1,2,2)
      check1(1)=matv(k,-1,1,2)
      check1(2)=matv(k,-1,2,2)
      Reno_inv=2
      musq=1d0
      goto 100
      endif
      If (Reno_inv.gt.0) then
       Print*, ' '
       Print*, 'chek1 hel 1 ',check1(1)
       Print*, 'check2 hel -1',check1(2)
       Print*, 'chek1 hel 1 ',check(1)
       Print*, 'check2 hel -1',check(2)
       Print*, ''
       Print*, 'check1 hel 1 ',matv(k,-1,1,2)
       Print*, 'check2 hel -1',matv(k,-1,2,2)
       Print*, 'check1 hel 1 ',(matv(k,-1,1,2)+(K_contcheck+K_contcheck4)*matb(k,-1,1,2))
       Print*, 'check2 hel -1',(matv(k,-1,2,2)+(K_contcheck+K_contcheck4)*matb(k,-1,2,2))
      
       Print*, 'ratio1',check(1)/((matv(k,-1,1,2)+(K_contcheck+K_contcheck4)*matb(k,-1,1,2)))-1d0
       Print*, 'ratio2',check(2)/((matv(k,-1,2,2)+(K_contcheck+K_contcheck4)*matb(k,-1,2,2)))-1d0
!       stop
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
      endif ! end lpent  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Finish Pentagons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Start BOXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               qq -> W g -> W A A g -----
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(lbox) then

      scale=musq!mursq(1,1)!-2d0*dotrr(p(0,1),p(0,2))

            
            do mu = 0,3
               p1(mu) = p(mu,1)
               p2(mu) = -p(mu,2)
	       mpg(mu) = -q(mu)
               mqwa12(mu) = -p1(mu)-p2(mu)+q(mu)
c for gauge test	       
	       emqwa12(mu)=mqwa12(mu)
	       empg(mu)=mpg(mu)
	   enddo
    
      
      if(bos.eq.311) then
      do mu=0,3
      wwa12(mu,1)=wpwaa(mu,1)
      enddo
      elseif(bos.eq.411) then
      do mu=0,3
      wwa12(mu,1)=wmwaa(mu,1)
      enddo
      endif
!!! Check Gauge invariance
      If(PerformGauge.eq.1) then
         do mu = 0,3
	    Myeps(mu,1)=mpg(mu)
	    Myeps(mu,2)=mpg(mu)
	  enddo
      elseif(PerformGauge.eq.2) then
         do mu = 0,3
	     mutemp=mu
	    wwa12(mutemp,N) = mqwa12(mu)
	 enddo
       endif    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W g 
C
      call boxlineABETotal(p1,mqwa12,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),Myeps(0,1),-1,scale,1,3,
     2    3,1,boxcontr1,j5g(1,1,1),boxcontr1b,jbg(1,1),Div)
c For non-abelian gauge and kk/M terms
      call boxlineABETotal(p1,mqwa12,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    emqwa12,Myeps(0,1),-1,scale,1,3,
     2    9,-1,nothing,j5gg(1,1,1),nothing,jbgg(1,1),Div)
      call boxlineABETotal(p1,mqwa12,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),empg,-1,scale,1,3,
     2    9,-1,nothing,j5gg2(1,1,1),nothing,jbgg2(1,1),Div)
c Gluon Helicty
       call boxlineABETotal(p1,mqwa12,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),Myeps(0,2),-1,scale,1,3,
     2    9,-1,boxcontr1,j5g(1,1,2),boxcontr1b,jbg(1,2),Div)
      call boxlineABETotal(p1,mqwa12,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    emqwa12,Myeps(0,2),-1,scale,1,3,
     2    9,-1,nothing,j5gg(1,1,2),nothing,jbgg(1,2),Div)
C
C g W
C
      call boxlineABETotal(p1,mpg,mqwa12,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wwa12(0,1),-1,scale,1,2,
     2    3,1,boxcontr2,j5g(1,2,1),boxcontr2b,jbg(2,1),Div)
c For non-abelian gauge and kk/M terms     
      call boxlineABETotal(p1,mpg,mqwa12,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),emqwa12,-1,scale,1,2,
     2    9,-1,nothing,j5gg(1,2,1),nothing,jbgg(2,1),Div)
      call boxlineABETotal(p1,mpg,mqwa12,p2,psi(1,-1,2),psi(1,-1,1),
     1    empg,wwa12(0,1),-1,scale,1,2,
     2    9,-1,nothing,j5gg2(1,2,1),nothing,jbgg2(2,1),Div)
c Gluon Helicty
      call boxlineABETotal(p1,mpg,mqwa12,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),wwa12(0,1),-1,scale,1,2,
     2    9,-1,boxcontr2,j5g(1,2,2),boxcontr2b,jbg(2,2),Div)
      call boxlineABETotal(p1,mpg,mqwa12,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,2),emqwa12,-1,scale,1,2,
     2    9,-1,nothing,j5gg(1,2,2),nothing,jbgg(2,2),Div)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W
C
      call boxlineNoAbeTotal(p1,mqwa12,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),Myeps(0,1),-1,scale,2,1,
     2    boxcontrNoAbe3,jg5n(1),Div)
c For non-abelian gauge and kk/M terms          
      call boxlineNoAbeTotal(p1,mqwa12,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    emqwa12,Myeps(0,1),-1,scale,9,-1,
     2    nothing,jg5ng(1),Div)
      call boxlineNoAbeTotal(p1,mqwa12,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),empg,-1,scale,9,-1,
     2    nothing,jg5ng2(1),Div)
c Gluon Helicty
      call boxlineNoAbeTotal(p1,mqwa12,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wwa12(0,1),Myeps(0,2),-1,scale,9,-1,
     2    boxcontrNoAbe3,jg5n(2),Div)
      call boxlineNoAbeTotal(p1,mqwa12,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    emqwa12,Myeps(0,2),-1,scale,9,-1,
     2    nothing,jg5ng(2),Div)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Gauge Boxes Abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C W g 
      call box_gauge(p1,mqwa12,mpg,p2,scale,
     2    boxcontr1,boxcontr1b,wardidtest1)
C g W
      call box_gauge(p1,mpg,mqwa12,p2,scale,
     2    boxcontr2,boxcontr2b,wardidtest2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Non-abelian Part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W
      wardidtest3=.false.
      if (abs(jg5ng(1)).lt.1d-3) then
           td=abs(-(j5gg(2,2,1)+j5gg(2,1,1))/2d0-jg5ng(1))
      else
            td=abs(-(j5gg(2,2,1)+j5gg(2,1,1))/2d0/jg5ng(1)+1d0)
      endif
      IF(td.le.1.5d-1) then
           if (abs(jg5ng2(1)).lt.1d-3) then
               td=abs(-(j5gg2(2,2,1)+j5gg2(2,1,1))/2d0-jg5ng2(1))
           else
                td=abs(-(j5gg2(2,2,1)+j5gg2(2,1,1))/2d0/jg5ng2(1)+1d0)
           endif
         IF(td.le.1.5d-1) wardidtest3=.true.
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
!	Print*, 'wardidtest1',wardidtest1
!	Print*, 'wardidtest2',wardidtest2
!	Print*, 'wardidtest3',wardidtest3
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
	      do l =1,2
                do i = 1,2
		  do mu=1,2
                   j5g(mu,i,l) = (0.0d0, 0.0d0)
		   
		   j5gg(mu,i,l) = (0.0d0, 0.0d0)
                  enddo
	        enddo	
 	           jg5n(l) =     (0.0d0, 0.0d0)
		   jg5ng(l)=     (0.0d0, 0.0d0)
	       enddo		   
 	      boxcountm = boxcountm+1
              boxcountps(ps_number) = boxcountps(ps_number)+1 
             else
	      do l=1,2
                  do i = 1,2
		     j5g(1,i,l)=(cfc*j5g(1,i,l)+cfcac*j5g(2,i,l))
		     j5gg(1,i,l)=(cfc*j5gg(1,i,l)+cfcac*j5gg(2,i,l))
                  enddo
	             jg5n(l)=cacgraph*jg5n(l)
		     jg5ng(l)=cacgraph*jg5ng(l)
	      enddo 
              boxcount2m = boxcount2m+1
              boxcount2ps(ps_number) = boxcount2ps(ps_number)+1      
             endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c   q qbar -> W^* g-> WAAg
c CHECK BORN 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       do l=1,2

        mv1b =(jbg(2,l)+  !g W^*
     2       jbg(1,l))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1b =((jbgg(2,l))*dotrc(mqwa12,wwa12(0,1))+  !g W^*
     2       (jbgg(1,l))*dotrc(mqwa12,wwa12(0,1)))*zm2i(3) !W^* g
C Change of sign to agree with Madgraph
        mv1b =-mv1b
	mmv1b =-mmv1b

       matb(k,-1,l,3)=(mv1b+mmv1b)* clr(3,3,-1)
      
       
       ratioComp=abs(mat(k,-1,l,3)/matb(k,-1,l,3)  -1d0)
	
      if(ratioComp.gt.1d-6) then	
      Print*,"Born/b_box-1d0",ratioComp
      Print*, "Born",mat(k,-1,l,3)
      Print*, 'mv1',mv1*clr(3,3,-1)
      Print*, 'mmv1',mmv1
      !stop
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
       mv1 =(j5g(1,2,l)+  !g W^*
     2       j5g(1,1,l))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1 =((j5gg(1,2,l))*dotrc(mqwa12,wwa12(0,1))+  !g W^*
     2       (j5gg(1,1,l))*dotrc(mqwa12,wwa12(0,1)))*zm2i(3)  !W^* g
C Change of sign to agree with Madgraph
        mv1 =-mv1
	mmv1 =-mmv1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
C  No abelian part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
	mv2=jg5n(l)
	mmv2=(jg5ng(l)*dotrc(mqwa12,wwa12(0,1)))*zm2i(3) 
        mv2 =-mv2
	mmv2 =-mmv2
	
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
	Print*, 'mv1b',mv1b
	Print*, 'mv1',mv1
	Print*, 'mv2',mv2
	Print*, 'Gauge',-mv2/mv1-1d0
	Print*, 'mmv1',mmv1
	Print*, 'mmv2',mmv2
	
	!stop
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
	
	Print*, 'mv1b',mv1b
	Print*, ''
	
	Print*, 'mv1    ',mv1
	Print*, ''
	Print*, 'mv2    ',mv2
	Print*, ''
	Print*, 'mmv1    ',mmv1
	Print*, ''
	Print*, 'mmv2    ',mmv2
	Print*, ''
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
       mv1CA =(j5g(2,2,l)+  !g W^*
     2       j5g(2,1,l))   !W^* g
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
        mmv1CA =(j5gg(2,2,l))*dotrc(mqwa12,wwa12(0,1))+  !g W^*
     2       (j5gg(2,1,l))*dotrc(mqwa12,wwa12(0,1))  !W^* g
C Change of sign to agree with Madgraph
        mv1CA =-mv1CA
	mmv1CA =-mmv1CA

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
	 
         ratioFact=mmv1/(-2d0*mmv1b)
	 Print*, 'mmv1',mmv1
	 Print*, 'mmv1b',mmv1b
!	 Print*, "CF Div/Eps^2", ratioFact
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
	 Print*, 'mvNoAbe',mv2
         Print*, 'mv1b',mv1b
	 ratioFact=(mv1CA*(-1d0/2d0)+mv2)/(-1d0*mv1b)
	 Print*, "CA Div/Eps^2", ratioFact
	 Print*, ""
 
	 Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv2
         Print*, 'mmv1b',mmv1b
!	 ratioFact=(mmv1CA*(-1d0/2d0)+mmv2)/(-1d0*mmv1b)
!	 Print*, "CA Div/Eps^2", ratioFact
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
	 
	 ratioFact=mmv1/(mmv1b*cteEps1)
	 Print*, 'mmv1',mmv1
	 Print*, 'mmv1b*cte',mmv1b*cteEps1
!	 Print*, "CF Div/Eps", ratioFact
	 Print*, ""
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CA part
	 cteEps1=(-lns+lnt+lnu)
	 
	 Print*, "CA PART"
	 Print*, ""
	 Print*, 'cteEps1',cteEps1
	 Print*, ""
	 Print*, 'mv1CA',mv1CA*(-1d0/2d0)
	 Print*, 'mvNoAbe',mv2
         Print*, 'mv1b*cte',mv1b*cteEps1
	 ratioFact=(mv1CA*(-1d0/2d0)+mv2)/(mv1b*cteEps1)
	 Print*, "CA Div/Eps", ratioFact
	 Print*, ""

	 Print*, 'mmv1CA',mmv1CA*(-1d0/2d0)
	 Print*, 'mmvNoAbe',mmv2
         Print*, 'mmv1b*cte',mmv1b*cteEps1
!	 ratioFact=(mmv1CA*(-1d0/2d0)+mmv2)/(mmv1b*cteEps1)
!	 Print*, "CA Div/Eps", ratioFact
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
        matv(k,-1,l,3) = (mv1+mv2+mmv1+mmv2)*clr(3,3,-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (PerformGauge.ge.1) then
! then, matv(k,-1,l,1) should be zero
      Print*,''
      Print*,'matv(k,-1,l,3) should be zero'
      Print*,' matv(k,-1,l,3)', matv(k,-1,l,3)
      Print*, 'mv1*clr(3,3,-1)',mv1*clr(3,3,-1)
      Print*,''
!      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
!      Print*,''
!      Print*,'then, matv(k,-1,l,3) with couplings corrected should be zero'
!      Print*,' matv(k,-1,l,1)', matv(k,-1,l,2)
!      cteEps1=clr(u(k),1,-1)/clr(d(k),1,-1) 
!      matv(k,-1,l,2) = (mv1+mv2+mmv1+mmv2)*clr(3,3,-1)
!      Print*,'matv(k,-1,l,2) with corrected couplings'
!      Print*,'matv(k,-1,l,2)', matv(k,-1,l,2)
!      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,3)',matv(k,-1,l,3)
      Print*, '-4*matvb(k,-1,l,3)',-4*matb(k,-1,l,3)
      Print*, 'ratio',matv(k,-1,l,3)/(-4*matb(k,-1,l,3))
      endif  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Add and subtract cte's 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	matv(k,-1,l,3)=matv(k,-1,l,3)+matb(k,-1,l,3)*convfact
	matv(k,-1,l,3)=matv(k,-1,l,3)-matb(k,-1,l,3)*(K_cont)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo  !  l helicity
        
!  Check renormalization invariance using K_contcheck      
      if(Reno_inv.eq.1) then
      check(1)=matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3)
      check(2)=matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3)
      check1(1)=matv(k,-1,1,3)
      check1(2)=matv(k,-1,2,3)
      Reno_inv=2
      musq=1d0
      goto 100
      endif
      If (Reno_inv.gt.0) then
       Print*, ' '
       Print*, 'chek1 hel 1 ',check1(1)
       Print*, 'check2 hel -1',check1(2)
       Print*, 'chek1 hel 1 ',check(1)
       Print*, 'check2 hel -1',check(2)
       Print*, ''
       Print*, 'check1 hel 1 ',matv(k,-1,1,3)
       Print*, 'check2 hel -1',matv(k,-1,2,3)
       Print*, 'check1 hel 1 ',(matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3))
       Print*, 'check2 hel -1',(matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3))
      
       Print*, 'ratio1',check(1)/((matv(k,-1,1,3)+(K_contcheck+K_contcheck4)*matb(k,-1,1,3)))-1d0
       Print*, 'ratio2',check(2)/((matv(k,-1,2,3)+(K_contcheck+K_contcheck4)*matb(k,-1,2,3)))-1d0
!       stop
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

       mf(1) = 0d0 ! light quarks
       mf(2) = tmass ! top
       mf(3) = bmass ! bottom
c prepare gluon current
cc <2W|1>
       call curr6(-1,psiw0(1,-1,2),fq0(0,2),
     &            psi(1,-1,1),p(0,1),jqqw(0,-1))
cc <2|W1>
       call curr6add(-1,psi(1,-1,2),p(0,2),
     &            psiw0(1,-1,1),fq0(0,1),jqqw(0,-1))
cc g propagator
       call propagate(-1,jqqw(0,-1),pqqw(0))
       do mu=0,3
         minus_pqqw(mu) = -pqqw(mu)
       enddo

c contract with box loops
ccc first anchor point is always internal propagator
ccc permutate g and the two photons over the other ones
cc g-1-2
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),q(0),qa1(0,N),qa2(0,N),
     &     jqqw(0,-1),myeps(0,l),atau1(1,N),atau2(1,N),
     &     scale,flgauge,flcomp,rfbox(1,j,l),rgfbox(1,j,l),Div)
       enddo
       enddo
cc g-2-1
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),q(0),qa2(0,N),qa1(0,N),
     &     jqqw(0,-1),myeps(0,l),atau2(1,N),atau1(1,N),
     &     scale,flgauge,flcomp,rfbox(2,j,l),rgfbox(2,j,l),Div)
       enddo
       enddo
cc 1-g-2
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),qa1(0,N),q(0),qa2(0,N),
     &     jqqw(0,-1),atau1(1,N),myeps(0,l),atau2(1,N),
     &     scale,flgauge,flcomp,rfbox(3,j,l),rgfbox(3,j,l),Div)
       enddo
       enddo
cc 2-g-1
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),qa2(0,N),q(0),qa1(0,N),
     &     jqqw(0,-1),atau2(1,N),myeps(0,l),atau1(1,N),
     &     scale,flgauge,flcomp,rfbox(4,j,l),rgfbox(4,j,l),Div)
       enddo
       enddo
cc 1-2-g
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),qa1(0,N),qa2(0,N),q(0),
     &     jqqw(0,-1),atau1(1,N),atau2(1,N),myeps(0,l),
     &     scale,flgauge,flcomp,rfbox(5,j,l),rgfbox(5,j,l),Div)
       enddo
       enddo
cc 2-1-g
       do j=1,3
       do l=1,2
         if (l.eq.1) then
           flgauge=1
         else
           flgauge=-1
         endif
         flcomp=flgauge
         call fermionboxT(mf(j),
     &     minus_pqqw(0),qa2(0,N),qa1(0,N),q(0),
     &     jqqw(0,-1),atau2(1,N),atau1(1,N),myeps(0,l),
     &     scale,flgauge,flcomp,rfbox(6,j,l),rgfbox(6,j,l),Div)
       enddo
       enddo

       wardidtest1 = .true.
       do j=1,3
         do i=1,6
           if (rgfbox(i,j,1).gt.1d-1) then 
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
       enddo
c add everything together
       if (wardidtest1) then
             box2count2m = box2count2m+1
             box2count2ps(ps_number) = box2count2ps(ps_number)+1   
ccc         write(*,*) "MARKER3: gauge test ok"
         do l=1,2
           mv1 = 0d0
           mv2 = 0d0
           do i=1,6
             mv1 = mv1 + 2*rfbox(i,1,l) ! d,s
     &                 +   rfbox(i,3,l) ! b
             mv2 = mv2 + 2*rfbox(i,1,l) ! u,c
     &                 +   rfbox(i,2,l) ! t
           enddo
c and multiply with color factor TR
           matv(k,-1,l,3)=0d0
           matv(k,-1,l,3) = matv(k,-1,l,3) +
     &                      (mv1*clr(d(k),1,-1)**2 +
     &                       mv2*clr(u(k),1,-1)**2)
     &                      * clr(3,3,-1) * TR 
ccc         mm(k,-1,l) = 0
ccc         do i = 1,3
ccc            mm(k,-1,l) = mm(k,-1,l) +  mat(k,-1,l,i)
ccc         enddo
ccc	write(*,*) "TTT1",l, mm(k,-1,l), abs(mm(k,-1,l)),
ccc     &                      ((mv1*clr(d(k),1,-1)**2 +
ccc     &                     mv2*clr(u(k),1,-1)**2)
ccc     &                    * clr(3,3,-1) * TR)*dconjg(mm(k,-1,l))
ccc	write(*,*) "RRR1",l,(mv1*clr(d(k),1,-1)**2 +
ccc     &                     mv2*clr(u(k),1,-1)**2)
ccc     &                    * clr(3,3,-1) * TR, abs(
ccc     &                      (mv1*clr(d(k),1,-1)**2 +
ccc     &                     mv2*clr(u(k),1,-1)**2)
ccc     &                    * clr(3,3,-1) * TR) 
         enddo
ccc	read(*,*)
        else !wardidtest
 	      box2countm = box2countm+1
              box2countps(ps_number) = box2countps(ps_number)+1 
ccc           if(Div.ne.0) then
ccc        write(*,*) "MARKER6: gauge test failed"
ccc       do l=1,2
ccc         mv1 = 0d0
ccc         mv2 = 0d0
ccc         do i=1,6
ccc           mv1 = mv1 + 2*rfbox(i,1,l) ! d,s
ccc     &               +   rfbox(i,3,l) ! b
ccc           mv2 = mv2 + 2*rfbox(i,1,l) ! u,c
ccc     &               +   rfbox(i,2,l) ! t
ccc         enddo
ccc	write(*,*) "failedRRR1",l,(mv1*clr(d(k),1,-1)**2 +
ccc     &                     mv2*clr(u(k),1,-1)**2)
ccc     &                    * clr(3,3,-1) * TR
ccc       enddo
ccc       read(*,*)
ccc       endif
       endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif !lbox
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif !lnlo
#endif 
#endif
! end of with_nlo and tribosonjet

#ifdef WITH_MADGRAPH
c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then
         
         do isig = -1,1,2
            do i = 1,3
                  
                  m1kb(1,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(1,isig,i) = 0d0 ! kartesian basis l = 2
                  
                  do lh = -1,1,2
                     madgr(i,isig,lh,1) = 0d0
                  enddo
                  
            enddo
         enddo
         
         do isig = -1,1,2
            do lh = -1,1,2
                  do i =1,3
                     
                     madgr(i,isig,lh,1) = tampr(i,isig,lh,1)
                     
                  enddo
            enddo
         enddo
         
         do isig = -1,1,2
            do i=1,3
                  
                  m1kb(1,isig,i) = (madgr(i,isig,-1,1) - 
     &                 madgr(i,isig,1,1))/sqrt(2d0)       
                  m2kb(1,isig,i) = -im*(madgr(i,isig,-1,1) + 
     &                 madgr(i,isig,1,1))/sqrt(2d0)
                  
            enddo
         enddo
         
c print results from checks:

         do i=1,3
               do isig = -1,-1,2
                  
                  print*,"mat(",k,",",isig,",1,",i,") = ",mat(k,isig,1,i)
                  temp = dreal(mat(k,isig,1,i))**2+dimag(mat(k,isig,1,i))**2
                  print*,"m1kb(",1,",",isig,",",i,") = ",m1kb(1,isig,i)
                  temp1 = dreal(m1kb(1,isig,i))**2+dimag(m1kb(1,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
                  
                  print*,"mat(",k,",",isig,",2,",i,") = ",mat(k,isig,2,i)
                  temp = dreal(mat(k,isig,2,i))**2+dimag(mat(k,isig,2,i))**2
                  print*,"m2kb(",1,",",isig,",",i,") = ",m2kb(1,isig,i)
                  temp1 = dreal(m2kb(1,isig,i))**2+dimag(m2kb(1,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
               enddo
         enddo
         

      endif ! end of debugging block
c------------------------------------------------------
#endif


c sum the graphs, square them and map them onto uuwwaj and ddwwaj

      res(k) = 0
      resv(k) = 0
      do l = 1,2
         mm(k,-1,l) = 0
         mmb(k,-1,l) = 0
         do i = 1,3
            mm(k,-1,l) = mm(k,-1,l) +  mat(k,-1,l,i)
            mmb(k,-1,l) = mmb(k,-1,l) +  matb(k,-1,l,i)
         enddo
         
         res(k)= res(k) + dreal(mm(k,-1,l))**2
     &        + dimag(mm(k,-1,l))**2*(1d0)
         
      enddo
c  add Born type term and multiply by F_q = alphas*C_2/4pi
c  the additional factor for the born term is after adding the 
c  subtraction term and the counter term for the renormalization
c  of the pdfs (see also D(x) in m2s_qqWWZj.f)

#ifdef WITH_NLO
      if (lnlo) then
         do l = 1,2
            mv(k,-1,l) = 0
            do i = 1,3
               mv(k,-1,l) = mv(k,-1,l) +  matv(k,-1,l,i)
            enddo
         enddo
      endif
#endif
      do l = 1,2
         if (nlo.eq.4) then
#ifdef WITH_NLO
            mv(k,-1,l) = als(1,1)*oneo4pi*
     1         (mv(k,-1,l) + mm(k,-1,l)*K_cont) 
#endif
         elseif(nlo.eq.7) then
#ifdef WITH_NLO
            mv(k,-1,l) = als(1,1)*oneo4pi*
     1         (mv(k,-1,l) + mm(k,-1,l)*K_cont) +factortoadd*mmb(k,-1,l)
#endif
         else if(nlo.lt.0) then
#ifdef WITH_NLO
            mv(k,-1,l) = als(1,1)*oneo4pi*(mv(k,-1,l))+factortoadd*mmb(k,-1,l)
#endif
         else
            mv(k,-1,l) = (0d0, 0d0)
         endif

         resv(k) = resv(k) + 2d0*dreal(mm(k,-1,l)*conjg(mv(k,-1,l)))

      enddo

      if (nlo.ge.0) then
          res(k) = (res(k)+resv(k))*4d0*fpials(2)  ! C_2*3 is the color factor
      else
         res(k) = resv(k)*4d0*fpials(2)           ! C_2*3 is the color factor
      endif  
      
      final = res(k)
      
      return
      end


