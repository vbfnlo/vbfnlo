****************subroutine QQVVj *************************************
c
c        Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c        Initial version:  2006 April
c        Last modified: 
C     Francisco Campanario, francam@particle.uni-physik.uni-karlsruhe.de
c      Last modified:  April 2012
c
c  QQWWj calculates the matrix elements**2 for qq -> WW g -> ve e+ mu- vm~ g
c
c        q1 q2 ----> W+ W- g,   W+ ---> f4-bar f3, W- ---> f6-bar f5
c
c  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c  Gluon polarization vectors are given in the kartesian basis 
c  (comparison with madgraph amplitudes tampr requires to express these in 
c   kart. basis also, even though they are given in the helicity basis 
c   a priori). 
c
c  This code includes only real emission contributions, i.e.
c
c      return uucc = |M_real|^2   etc.
c
c        fpials is attached only in the end of the code
c
c index j = 2:3 indicates, whether g is emitted from 
c                upper line (j=2) or lower line (j=3)
c        l is the gluon polarization in the kartesian basis (l=1,2)
c                l=0 stands for building blocks without gluon emission
c        k is the process ID (1:uuww,2:ddww)
c        isig1/isig3 are the helicities of parton1/2 
c
c*********************************************************************

      subroutine qqVVj(pbar,sign, qbar,gsign, nlo, N, uuww, ddww, 
     &                 ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
C#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/koppln.inc"

c electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program      
      real*8 fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), musq, q(0:3)
      double precision  uuww, ddww, res(2),resv(2)
      double precision  p(0:3,2+n_v), pq(0:4,2)
      double precision pgWP(0:4,2), pgWM(0:4,2), pWP(0:4,2), pWM(0:4,2)

      integer  sign(2+n_v), gsign, mu, i, k, isig, nlo, N
      integer  ifl(2,2), ps_number
      integer  l,j   ! gluon polariz. (l=0:no g, l=1,2:g in kartesian basis)


      double complex mat(2,-1:1,0:2,9),temp, temp1
      double complex mata(2,-1:1,0:2,9),
     1                      matz(2,-1:1,0:2,9)
      double complex mm(2,-1:1,2)
      double complex  m1(-1:1,2,3), m2(-1:1,2,3)
      double complex  raz! for checks only
      double precision eps(0:3,2) ! g in kartesian basis 
      double complex psi(2,-1:1,4), psiWM(2,-1:1,2), psiWP(2,-1:1,2)
      double complex psiWMg(2,-1:1, 2,2), psiWPg(2,-1:1, 2,2)
      double complex braketg(2,-1:1,2,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2), 
     &               jqqa(0:5,-1:1,2)
      double complex jqqgWP(0:5,-1:1,2,0:3), jqqgWM(0:5,-1:1,2,0:3)

      double complex braketgWP(2,-1:1,2,2),braketgWM(2,-1:1,2,2)         
      double precision  pwpg(0:4,4),pwmg(0:4,4)
      double complex ma(-1:1,2,0:2), mz(-1:1,2,0:2)
      double complex  zm2i(2:4)
      double complex contract_Tjj,! contract_T1j,contract_T2j,
     1                      dotcc, dotrc, dotqj, s1c
      double precision dotrr
      external contract_Tjj, !contract_T1j,contract_T2j,
     1                dotcc, dotrc, dotqj, s1c,dotrr
      integer lh
      double complex madgr(9,-1:1,-1:1,2) !checks

      double complex m1kb(2,-1:1,9),m2kb(2,-1:1,9) !checks

      double complex im
      parameter (im = (0d0,1d0))
      logical lerror
     
      save ifl, zm2i!, linit
      data  lerror /.true./

      logical ldebug, ldebugm, ldebuglt, linit, christophmods
      logical lnlo, lbox, lpt, lpent,lbox2,lhiggs,lvert

      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false., ldebuglt=.false.)
      parameter(christophmods=.false.)

c variables for virtual corrections
      double complex matb(2,-1:1,0:2,9)
      double complex matv(2,-1:1,0:2,9)
      double complex tmpv(2,-1:1,0:2,9)
      double complex tmpb(2,-1:1,0:2,9)
      double complex    m3, m4
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
!      complex*16 j5a2g1(8,4,2,-1:1),jba2g1(4,2,-1:1)
!      complex*16 j5ga21(8,4,2,-1:1),jbga21(4,2,-1:1)
c


c
!      complex*16 j5ma2g1(8,4,2,-1:1),jbma2g1(4,2,-1:1)
!      complex*16 j5mga21(8,4,2,-1:1),jbmga21(4,2,-1:1)

c NoAbe
      complex*16 jg5a2(2,2,-1:1)

!      complex*16 jg5a21(2,2,-1:1)
!      complex*16 jg5ma21(2,2,-1:1)
c      complex*16 jg5ma2(2,2,-1:1)
      complex*16 pentcontrNoAbe1(3),pentcontrNoAbe2(3)

      real*8      accuracyres


c  Colors Parameters       
      real*8 NF,inv2,inv3,CA,TR,tf,CF,eulergamma
      parameter (NF=5d0,inv2=1d0/2d0,inv3=1d0/3d0)
      parameter (CA=3d0)
      parameter (TR=1d0/2d0,tf=1d0/2d0)
      parameter (CF=4d0/3d0)
      parameter (eulergamma=5.772156649015328D-1)
      complex*16 cfc,cac,cfcac,pic2,cacgraph,convfact
      real*8 mans,mant,manu!,musq
c     Variables for virtual
      complex*16 mv1,mv2,mv3,mv4,mv5,mv6
      complex*16 mmv2,mmv3,ratioFACT,mmv2b,mmv3b,mmv5,mmv6
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
       complex*16 factortoadd
       real*8 s,t,mmu,mursq1
       complex*16 check(2,-1:1,2),check1(2,-1:1,2)
       complex*16 check2(2,-1:1,2,7),check3(2,-1:1,2,7)

c
       integer isig2,comp
       double complex cnf
      Integer Div,PerformGauge,mutemp
      complex*16 mv1CA,mv2CA,mv3CA,mv1b,mv2b,mv3b,mv4b
      complex*16 mmv4b,mv4CA

      complex*16 lns,lnt,lnu,cteEps1

      Logical wardidtest1,wardidtest2,wardidtest3,wardidtest4
      Logical wardidtest5,wardidtest6,wardidtest7,wardidtest8
!      integer N
      real*8 theta
      external theta

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

!      double complex zm2i(2:3)
!      save zm2i
!      integer nlo

c Variables for fermion loops

      double precision mf(4)
      integer cpl(4)
      double precision vv
      double complex rfbox(6,2,4,2,2)
      double precision rgfbox(6,2,4,2,2) 
      double complex rfvert(2,3,4,2,2)
      double precision rgfvert(2,3,4,2,2) 

      double precision fcresult
      integer*8 fchelicities, thirtytwo
      integer fcflags
      double complex fcamp(0:1,-2:2,-2:2), fckamp(0:1,2)
      parameter(thirtytwo=32)

      real*8 minus_qp(0:3),minus_qm(0:3),minus_qww(0:3)
      integer poln
      external deltakro
      real*8 deltakro
      complex*16 curr1(0:3,2),curr2(0:5,-1:1)
      double precision p12(0:4,-1:1)
      integer l1,l2
      external sc1c
      complex*16 sc1c
!      parameter (lhiggs=.true.)
!      double precision resv(2)
      complex*16 mmb(2,-1:1,2),mv(2,-1:1,2)
c  helicity selection
      INTEGER h
      COMMON /hcount / h

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
!       print*, "lnlo",lnlo

!      lnlo=.true.
!     lpent=.true.
!       lbox=.FALSE.
!      lbox2=.true.

c       print*, "lbox",lbox
c       print*, "lpent",lpent
c       print*, "lbox2",lbox2


c  ---------------------------------------------------------------------
 
c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the 2 quark lines:

      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c       fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.

         k = 1                  ! uuww
         ifl(1,k) = 3
         ifl(2,k) = 3
         k = 2                  ! ddww
         ifl(1,k) = 4
         ifl(2,k) = 4

         zm2i(2) = 1/dcmplx(xm2(2),-xmg(2))
c         zm2i(3) = 1/dcmplx(xm2(3),-xmg(3))
c         zm2i(4) = 1/dcmplx(xm2(4),-xmg(4))

      endif ! linit

      do k = 1,2
         do isig = -1,1,2
            do l = 0,2
               do i = 1,9
                  mat(k,isig,l,i)  = 0
                  mata(k,isig,l,i) = 0
                  matz(k,isig,l,i) = 0
                  matb(k,isig,l,i) =0 
                  matv(k,isig,l,i) =0
                  tmpb(k,isig,l,i) =0 
                  tmpv(k,isig,l,i) =0
               enddo
            enddo
         enddo
      enddo

!      m1=0d0
!      m2=0d0
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

c---------calculation of matrix elements--------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c     NOTE: psi(2,-1:1,2) is a two component spinor with helicity -1
c     or 1.  The last entry identifies the fermion.  If this entry is
c     odd psi is a ket, if even psi is a bra.
c     psi(1,isig1,1) = |1>_isig1
c     psi(1,isig1,2) = <2|_isig

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2        ! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

         do mu = 0,3                   !QCD gaugecheck      
c            eps(mu,l) = qbar(mu)
            Myeps(mu,l)=eps(mu,l)
         enddo

c--- End of QCD gaugecheck -------------------------------------------    

         do isig = -1,1,2        ! fermion helicity 

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     $           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig

            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     $           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2

         enddo
      enddo

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)

c------- T-CHANNEL -----------------

c possibility 1 - emission from upper line

      do l=1,2

      call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qp(0,N),wp(1,N),braketgWP(1,-1,1,l),pgWP(0,1)) !|W+,g,1>_l,-1

      call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qm(0,N),wm(1,N),braketgWM(1,-1,1,l),pgWM(0,1)) !|W-,g,1>_l,-1

      call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgWP(1,-1,1,l),pgWP(0,1),jqqgWP(0,-1,l,1))

      call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgWM(1,-1,1,l),pgWM(0,1),jqqgWM(0,-1,l,1))

      m1(-1,l,1) = dotcc(jqqgWP(0,-1,l,1),wm(1,N))
      m2(-1,l,1) = dotcc(jqqgWM(0,-1,l,1),wp(1,N))

      mat(1,-1,l,1) = clr(3,3,-1)**2 * m1(-1,l,1)
      mat(2,-1,l,1) = clr(3,3,-1)**2 * m2(-1,l,1)

      enddo

c             - emission from lower line

      do l=1,2

      call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &       qm(0,N),wm(1,N),braketgWM(1,-1,2,l),pgWM(0,2)) !<2,g,W-|_l,-1

      call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &       qp(0,N),wp(1,N),braketgWP(1,-1,2,l),pgWP(0,2)) !<2,g,W+|_l,-1

      call curr6(1,braketgWM(1,-1,2,l),pgWM(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWM(0,-1,l,2))

      call curr6(1,braketgWP(1,-1,2,l),pgWP(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWP(0,-1,l,2))

      m1(-1,l,2) = dotcc(jqqgWM(0,-1,l,2),wp(1,N))
      m2(-1,l,2) = dotcc(jqqgWP(0,-1,l,2),wm(1,N))

      mat(1,-1,l,2) = clr(3,3,-1)**2 * m1(-1,l,2)
      mat(2,-1,l,2) = clr(3,3,-1)**2 * m2(-1,l,2)

      enddo

c             - emission from intermediate quark line

      do l=1,2

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     &       qm(0,N),wm(1,N),psiWM(1,-1,2),pWM(0,2)) !<2,W-,g|_-1

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,
     &       qp(0,N),wp(1,N),psiWP(1,-1,2),pWP(0,2)) !<2,W+,g|_-1

      call bra2r(psiWM(1,-1,2),.false.,pWM(0,2),-1,
     &            q,eps(0,l),psiWMg(1,-1,2,l),pWMg(0,2)) !<2,W-,g|_l,-1

      call bra2r(psiWP(1,-1,2),.false.,pWP(0,2),-1,
     &            q,eps(0,l),psiWPg(1,-1,2,l),pWPg(0,2)) !<2,W+,g|_l,-1

      call curr6(1,psiWMg(1,-1,2,l),pWMg(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWM(0,-1,l,3))

      call curr6(1,psiWPg(1,-1,2,l),pWPg(0,2),
     &            psi(1,-1,1),p(0,1),jqqgWP(0,-1,l,3))

      m1(-1,l,3) = dotcc(jqqgWM(0,-1,l,3),wp(1,N))
      m2(-1,l,3) = dotcc(jqqgWP(0,-1,l,3),wm(1,N))

      mat(1,-1,l,3) = clr(3,3,-1)**2 * m1(-1,l,3)
      mat(2,-1,l,3) = clr(3,3,-1)**2 * m2(-1,l,3)

      enddo

      do l=1,2

            mat(1,-1,l,8) = mat(1,-1,l,1) + mat(1,-1,l,2)
     &                      + mat(1,-1,l,3)

            mat(2,-1,l,8) = mat(2,-1,l,1) + mat(2,-1,l,2)
     &                      + mat(2,-1,l,3)

      enddo

c------- end of T-CHANNEL computation

c-------  S-CHANNEL ----------------

c     Get the f-fbar currents with one gluon radiated from the
c     current line.  There are two terms, one for gluon emission to
c     either side of the ffV vertex:

      do l=1,2
         do isig = -1,1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))

         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqagq(0,-1,l))

            do mu = 0,5
               jqqa(mu,isig,l) = jqgaq(mu,isig,l) + jqagq(mu,isig,l)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig)
            enddo
         enddo
      enddo

      do l=1,2
         do isig=-1,1,2

            ma(isig,l,1) = -dotcc(aww(0,N),jqgaq(0,isig,l))
            ma(isig,l,2) = -dotcc(aww(0,N),jqagq(0,isig,l))
            ma(isig,l,0) = -dotcc(aww(0,N),jqqa(0,isig,l))

            mz(isig,l,1) = -dotcc(zww(0,N),jqgaq(0,isig,l))
            mz(isig,l,2) = -dotcc(zww(0,N),jqagq(0,isig,l))
            mz(isig,l,0) = -dotcc(zww(0,N),jqqa(0,isig,l))

            do k = 1,2          !k=1: uuww, k=2: ddww
               mat(k,isig,l,4) = ma(isig,l,1)*clr(ifl(1,k),1,isig)
               mat(k,isig,l,5) = ma(isig,l,2)*clr(ifl(1,k),1,isig)
               mat(k,isig,l,6) = mz(isig,l,1)*clr(ifl(1,k),2,isig)
               mat(k,isig,l,7) = mz(isig,l,2)*clr(ifl(1,k),2,isig)
               mat(k,isig,l,9) = mat(k,isig,l,4) + mat(k,isig,l,5)
     &                         + mat(k,isig,l,6) + mat(k,isig,l,7)
            enddo
         enddo
      enddo   !end gluon polarization l

c-------- end of S-CHANNEL computation ---------------------------
#ifdef WITH_MADGRAPH
c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,9
               do k = 1,2

                  m1kb(k,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(k,isig,i) = 0d0 ! kartesian basis l = 2

                  do lh = -1,1,2
                     madgr(i,isig,lh,k) = 0d0
                  enddo

               enddo
            enddo
         enddo

c compute m1/m2 in kartesian basis:

         do isig = -1,1,2
            do lh = -1,1,2
               do k = 1,2

         madgr(1,isig,lh,k) = tampr(1,isig,lh,k)
         madgr(2,isig,lh,k) = tampr(2,isig,lh,k)
         madgr(3,isig,lh,k) = tampr(3,isig,lh,k)
         madgr(4,isig,lh,k) = tampr(4,isig,lh,k)
     &                 + tampr(10,isig,lh,k) + tampr(14,isig,lh,k)
         madgr(5,isig,lh,k) = tampr(6,isig,lh,k)
     &                 + tampr(12,isig,lh,k) + tampr(16,isig,lh,k)
         madgr(6,isig,lh,k) = tampr(5,isig,lh,k)
     &                 + tampr(8,isig,lh,k) + tampr(11,isig,lh,k)
     &                 + tampr(15,isig,lh,k) + tampr(18,isig,lh,k)
         madgr(7,isig,lh,k) = tampr(7,isig,lh,k)
     &                 + tampr(9,isig,lh,k) + tampr(13,isig,lh,k)
     &                 + tampr(17,isig,lh,k) + tampr(19,isig,lh,k)
         madgr(8,isig,lh,k) = madgr(1,isig,lh,k)
     &                 + madgr(2,isig,lh,k) + madgr(3,isig,lh,k)
         madgr(9,isig,lh,k) = madgr(4,isig,lh,k) + madgr(5,isig,lh,k)
     &                 + madgr(6,isig,lh,k) + madgr(7,isig,lh,k)

               enddo
            enddo
         enddo


         do isig = -1,1,2
            do i=1,9
               do k=1,2
               
         m1kb(k,isig,i) = (madgr(i,isig,-1,k) - 
     &            madgr(i,isig,1,k))/sqrt(2d0)       
         m2kb(k,isig,i) = -im*(madgr(i,isig,-1,k) + 
     &            madgr(i,isig,1,k))/sqrt(2d0)

               enddo
            enddo
         enddo


c print results from checks:
         do isig = -1,1,2
            do k = 2,2
               do i = 4,7

                print*,"ratio for l=1, isig = ",isig,"i = ",i,":",
     &                    mat(k,isig,1,i)/m1kb(k,isig,i)
            print*,"mat(",k,",",isig,",1,",i,") = ",mat(k,isig,1,i)
            temp = dreal(mat(k,isig,1,i))**2+dimag(mat(k,isig,1,i))**2
            print*,temp
            print*,"m1kb(",k,",",isig,",",i,") = ",m1kb(k,isig,i)
            temp1 = dreal(m1kb(k,isig,i))**2+dimag(m1kb(k,isig,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

                print*,"ratio for l=2, isig = ",isig,"i = ",i,":",
     &                    mat(k,isig,2,i)/m2kb(k,isig,i)
            print*,"mat(",k,",",isig,",2,",i,") = ",mat(k,isig,2,i)
            temp = dreal(mat(k,isig,2,i))**2+dimag(mat(k,isig,2,i))**2
            print*,temp
            print*,"m2kb(",k,",",isig,",",i,") = ",m2kb(k,isig,i)
            temp1 = dreal(m2kb(k,isig,i))**2+dimag(m2kb(k,isig,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

               enddo
            enddo
         enddo

      endif ! end of debugging block
c------------------------------------------------------
#endif /* MADGRAPH */
c-----------------END of calculation of matrix elements------------------


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
!     1               +2d0*CFC*(-pi*pi+3d0)        
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
! whether external particles are initial/final does not matter, so we use 
! the kinematics of the internal ME ordering

      s= 2d0*dotrr(pbar(0,1),pbar(0,2))
      t=-2d0*dotrr(pbar(0,1),qbar(0)) 
      mmu=-2d0*dotrr(pbar(0,2),qbar(0)) 
      mursq1=4.d0*pi*mursq(1,1)

      factortoadd=
     &  (25*als(1,N)*CA)/(9D0*Pi) + (5*als(1,N)*CF)/Pi - 
     &  (11*als(1,N)*CA*EulerGamma)/(12d0*Pi) - 
     &  (3*als(1,N)*CF*EulerGamma)/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma**2)/(4D0*Pi) + 
     &  (als(1,N)*CF*EulerGamma**2)/(2D0*Pi) - 
     &  (7*als(1,N)*CA*Pi)/24d0 - 
     &  (7*als(1,N)*CF*Pi)/12d0 - (8*als(1,N)*Nf*TR)/(9D0*Pi) + 
     &  (als(1,N)*EulerGamma*Nf*TR)/(3D0*Pi) - 
     &  (3*als(1,N)*CA*Log(mursq1/S))/(4D0*Pi) + 
     &  (3*als(1,N)*CF*Log(mursq1/S))/(2D0*Pi) + 
     &  (als(1,N)*CA*EulerGamma*Log(mursq1/S))/(2D0*Pi) - 
     &  (als(1,N)*CF*EulerGamma*Log(mursq1/S))/Pi - 
     &  (als(1,N)*CA*Log(mursq1/S)**2)/(4D0*Pi) + 
     &  (als(1,N)*CF*Log(mursq1/S)**2)/(2D0*Pi) + 
     &  (5*als(1,N)*CA*Log(-(mursq1/T)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/T)))/(2D0*Pi) - 
     &  (als(1,N)*Nf*TR*Log(-(mursq1/T)))/(6D0*Pi) + 
     &  (als(1,N)*CA*Log(-(mursq1/T))**2)/(4D0*Pi) + 
     &  (5*als(1,N)*CA*Log(-(mursq1/mmU)))/(6D0*Pi) - 
     &  (als(1,N)*CA*EulerGamma*Log(-(mursq1/mmU)))/(2D0*Pi) - 
     &  (als(1,N)*Nf*TR*Log(-(mursq1/mmU)))/(6D0*Pi) + 
     &  (als(1,N)*CA*Log(-(mursq1/mmU))**2)/(4D0*Pi)

      factortoadd=factortoadd/2d0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                 Here start the virtual staff
C
!         LPENT=.FALSE.
!         LBOX=.FALSE.
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
c W-
               mqa2(mu) = -qm(mu,1)
c W+
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
c Polarization of W+
            wp(mutemp,1) = mqz(mu)
c            atau1(mutemp,1) = mqz(mu)
         enddo
      elseif(PerformGauge.eq.3) then
         do mu = 0,3
            mutemp=mu+1
c Polariazation A
            wm(mutemp,1) = mqa2(mu)
         enddo
       endif
CCCCCCCCCCCCCCCCCCCCCCCCC
C  Abelian PART
CCCCCCCCCCCCCCCCCCCCCCCC
C
C W- g W+
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2
       call penlineABETotal(p1,mqa2,mpg,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    wm(1,1),Myeps(0,isig2),wp(1,1),isig,scale,1,3,
     2    comp*4,comp,pentcontr1,j5a2g(1,3,isig2,isig),pentcontr1b
     3    ,jba2g(3,isig2,isig),Div)

        comp=-1

         enddo
       enddo   
C
C g W- W+
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2

       call penlineABETotal(p1,mpg,mqa2,mqz,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),wm(1,1),wp(1,1),isig,scale,1,2,
     2    comp*4,comp,pentcontr2,j5ga2(1,3,isig2,isig),pentcontr2b
     3    ,jbga2(3,isig2,isig),Div)

       comp=-1
         enddo
       enddo   
C g W- W+

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W- W+ g
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2
       call penlineABETotal(p1,mqa2,mqz,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wm(1,1),wp(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentcontr3,j5a2g(1,2,isig2,isig),pentcontr3b
     3    ,jba2g(2,isig2,isig),Div)
        comp=-1
         enddo
       enddo   
         
C
C g W+ W-
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2
       call penlineABETotal(p1,mpg,mqz,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),wp(1,1),wm(1,1),isig,scale,1,2,
     2    comp*4,comp,pentcontr4,j5ga2(1,2,isig2,isig),pentcontr4b
     3    ,jbga2(2,isig2,isig),Div)
       comp=-1
         enddo
       enddo   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W+ W- g
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2
       call penlineABETotal(p1,mqz,mqa2,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),wm(1,1),Myeps(0,isig2),isig,scale,1,4,
     2    comp*4,comp,pentcontr5,j5a2g(1,1,isig2,isig),pentcontr5b
     3    ,jba2g(1,isig2,isig),Div)
        comp=-1
         enddo
       enddo   

C
C W+ g W-
C
       comp=1
      do isig= -1,-1,2
        do isig2 =1,2
       call penlineABETotal(p1,mqz,mpg,mqa2,p2,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),Myeps(0,isig2),wm(1,1),isig,scale,1,3,
     2    comp*4,comp,pentcontr6,j5ga2(1,1,isig2,isig),pentcontr6b
     3    ,jbga2(1,isig2,isig),Div)
       comp=-1
         enddo
       enddo   

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C  NO ABELIAN PART
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W- W
C
       comp=1
      do isig= -1,-1,2 ! spinors
        do isig2 =1,2 ! Gluon
      call penlineNoAbeTotal(p1,mqa2,mqz,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wm(1,1),wp(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentcontrNoAbe1,jg5a2(2,isig2,isig),Div)
      comp=-1
         enddo
       enddo   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W W-
C
       comp=1
      do isig= -1,-1,2 ! spinors
        do isig2 =1,2 ! Gluon
      call penlineNoAbeTotal(p1,mqz,mqa2,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    wp(1,1),wm(1,1),Myeps(0,isig2),isig,scale,4*comp,comp,
     2    pentcontrNoAbe2,jg5a2(1,isig2,isig),Div)
       comp=-1
         enddo
       enddo   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Gauge Test
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W- g W+
      call pent_gauge(p1,mqa2,mpg,mqz,p2,psi(1,-1,2),psi(1,-1,1),
     1    wm(1,1),Myeps(0,1),wp(1,1),-1,scale,
     2    pentcontr1,pentcontr1b,wardidtest1)
C g W- W+
       call pent_gauge(p1,mpg,mqa2,mqz,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wm(1,1),wp(1,1),-1,scale,
     2   pentcontr2,pentcontr2b,wardidtest2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W- W+ g
       call pent_gauge(p1,mqa2,mqz,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wm(1,1),wp(1,1),Myeps(0,1),-1,scale,
     2    pentcontr3,pentcontr3b,wardidtest3)
C g W+ W-
      call pent_gauge(p1,mpg,mqz,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    Myeps(0,1),wp(1,1),wm(1,1),-1,scale,
     2    pentcontr4,pentcontr4b,wardidtest4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W+ W- g
      call pent_gauge(p1,mqz,mqa2,mpg,p2,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),wm(1,1),Myeps(0,1),-1,scale,
     2    pentcontr5,pentcontr5b,wardidtest5)
C W+ g W-
      call pent_gauge(p1,mqz,mpg,mqa2,p2,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),Myeps(0,1),wm(1,1),-1,scale,
     2    pentcontr6,pentcontr6b,wardidtest6)
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  NO ABELIAN PART
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W- W
      call pentNoAbe_gauge(1d-1,p1,mqa2,mqz,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wm(1,1),wp(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe1,accuracyres,wardidtest7)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C W W-
      call pentNoAbe_gauge(1d-1,p1,mqz,mqa2,p2,mpg,psi(1,-1,2),psi(1,-1,1),
     1    wp(1,1),wm(1,1),Myeps(0,1),-1,scale,
     2    pentcontrNoAbe2,accuracyres,wardidtest8)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!

       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3.and.
     1 wardidtest4.and.wardidtest5.and.wardidtest6
       wardidtest1=wardidtest1.and.wardidtest7.and.wardidtest8


        If(Div.ne.0) then
! Check factorization of Divergences. Set the color equal.
! So that I can check separately the CF and CA factorization
!        Print*,'CFCFCF',cfc
        cfc=1d0
        cacgraph=1d0
!        Print*,'CFCFCF',cfc
        cfcac=cfc
        wardidtest1=.True.
        endif
        
        
         if (.not.wardidtest1) then
          do isig=-1,-1,2  
            do  isig2=1,2
                do i = 1,3
                   do mu=1,2
                     j5a2g(mu,i,isig2,isig) = (0.0d0, 0.0d0)
                      j5ga2(mu,i,isig2,isig) = (0.0d0, 0.0d0)
                   jg5a2(mu,isig2,isig)=(0.0d0, 0.0d0)
                   enddo
                enddo
              enddo
           enddo 
              pentcountm = pentcountm+1
              pentcountps(ps_number) = pentcountps(ps_number)+1 

c           print*, "NO PEN GAUGE", ps_number,pentcountm,pentcountps(ps_number)
         else

          do isig=-1,-1,2  
            do isig2=1,2
              do i = 1,3
           j5a2g(1,i,isig2,isig)=(cfc*j5a2g(1,i,isig2,isig)
     &                         +cfcac*j5a2g(2,i,isig2,isig))
           j5ga2(1,i,isig2,isig)=(cfc*j5ga2(1,i,isig2,isig)
     &                          +cfcac*j5ga2(2,i,isig2,isig))


               enddo
             do i = 1,2
                 jg5a2(i,isig2,isig) = cacgraph*jg5a2(i,isig2,isig)

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
      do isig= -1,-1,2
         do isig2=1,2

!  gW+W-  W+gW-  W+W-g
        mv1b= jbga2(1,isig2,isig)+ jbga2(2,isig2,isig)+ jba2g(1,isig2,isig)
!  gW-W+  W-gW+  W-W+g
        mv2b= jba2g(2,isig2,isig)+ jba2g(3,isig2,isig)+ jbga2(3,isig2,isig)
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   

       matb(1,isig,isig2,1)= mv1b*clr(3,3,-1)**2
       matb(2,isig,isig2,1)= mv2b*clr(3,3,-1)**2
     
c Change of sign to agree with MadGraph

      matb(1,isig,isig2,1)=- matb(1,isig,isig2,1)
      matb(2,isig,isig2,1)=- matb(2,isig,isig2,1)


       ratioComp=
     &    abs(mat(1,isig,isig2,8)/matb(1,isig,isig2,1) -1d0)
     &   +abs(mat(2,isig,isig2,8)/matb(2,isig,isig2,1) -1d0)
      
      if(ratioComp.gt.1d-6) then

      print*
      print*
      print*, "ISIG",ISIG, "ISIG2",ISIG2
      print*
      print*
      PRINT*, 'Madgr/penline 1', ratioComp

      Print*, 'wardidtest1_PEN', wardidtest1
      Print*,"Born/b_pen-1d0",ratioComp
      Print*, "Born_uu",mat(1,isig,isig2,8),mat(2,isig,isig2,8)
      Print*, 'Sum_uu',matb(1,isig,isig2,1),matb(2,isig,isig2,1)
      Print*, 'ratio_uu',mat(1,isig,isig2,8)/matb(1,isig,isig2,1)
      Print*, 'ratio_uu',mat(2,isig,isig2,8)/matb(2,isig,isig2,1)

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

!    W+gW-  gW+W-  W+W-g
        mv1= j5ga2(1,1,isig2,isig)+ j5ga2(1,2,isig2,isig)+ j5a2g(1,1,isig2,isig)
!  gW-W+  W-gW+  W-W+g
        mv2= j5a2g(1,2,isig2,isig)+ j5a2g(1,3,isig2,isig)+ j5ga2(1,3,isig2,isig)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Non-Abelian graphs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W+W-
        mv3=(jg5a2(1,isig2,isig))
! W-W+
        mv4=(jg5a2(2,isig2,isig))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   SUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

       matv(1,isig,isig2,1)=
     1     (mv1+mv3)*clr(3,3,-1)**2

       matv(2,isig,isig2,1)=
     1     (mv2+mv4)*clr(3,3,-1)**2

c Change of sign to agree with MadGraph

      matv(1,isig,isig2,1)=- matv(1,isig,isig2,1)
      matv(2,isig,isig2,1)=- matv(2,isig,isig2,1)

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


        Print*, ''
        Print*, 'mv1b',mv1b
        Print*, 'mv1',mv1
        Print*, 'mv3',mv3
        Print*, 'Gauge',-mv1/mv3-1d0


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
        Print*, 'matb1', matb(1,isig,isig2,1)
        Print*, 'matb2', matb(2,isig,isig2,1)
        print*,"-m1/m2-1d0 ",-matb(1,isig,isig2,1)/matb(2,isig,isig2,1)-1d0

        Print*, 'matb1', matv(1,isig,isig2,1)
        Print*, 'matb2', matv(2,isig,isig2,1)
        print*, "-m1/m2-1d0",-matv(1,isig,isig2,1)/matv(2,isig,isig2,1)


c        Print*, 'Gauge',-mv2/mv4-1d0
 
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

!    W+gW-  gW+W-  W+W-g
        mv1CA= j5ga2(2,1,isig2,isig)+ j5ga2(2,2,isig2,isig)+ j5a2g(2,1,isig2,isig)
!  gW-W+  W-gW+  W-W+g
        mv2CA= j5a2g(2,2,isig2,isig)+ j5a2g(2,3,isig2,isig)+ j5ga2(2,3,isig2,isig)

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
         ratioFact=mv2/(-2d0*mv2b)
         Print*, 'mv1',mv1
         Print*, 'mv1b',mv1b
         Print*, "CF Div/Eps^2", ratioFact
         Print*, ""
         

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
         Print*, ""
         Print*, ""
         Print*, 'mv1CA',mv1CA*(-1d0/2d0)
         Print*, 'mvNoAbe',mv3
         Print*, 'mv1b',mv1b
         ratioFact=(mv1CA*(-1d0/2d0)+mv3)/(-1d0*mv1b)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (PerformGauge.eq.1) then
! then, matv(k,-1,l,1) should be zeror
      Print*,''
      Print*,' should be zero'
      Print*,' matv(1,isg,isg2,1)', matv(1,isig,isig2,1)
      Print*,' matv(2,isg,isg2,1)', matv(2,isig,isig2,1)
      Print*, 'mv1*clr(3,3,-1)',mv2*clr(3,3,isig)*clr(3,3,isig)
c      Print*,' matv(2,isg,isg2,1)', -matv(2,isig,isig2,1)/ matv(1,isig,isig2,1) 

      elseif(PerformGauge.gt.1) then
cFC! then, matv(k,-1,l,1) with couplings corrected should be zero
cFC      Print*,''
cFC      Print*,'matv(k,-1,l,2) should be zero'
cFC      Print*,' matv(k,-1,l,1)', matv(1,isig,isig2,1)
cFC     &                /(jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig)  )
cFC      PRINT*, 'jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig) ',jg5a2(1,isig2,isig)*clr(3,2,isig)*clr(3,1,isig) 
cFCc      cteEps1=1!clr(u(k),1,-1)/clr(d(k),1,-1) 
cFCc      matv(k,-1,l,2) = (mv2+mv4+(mv3+mv5)*cteEps1)!*clr(3,3,-1)
cFCc      Print*,'matv(k,-1,l,2) with corrected couplings'
cFCc      Print*,'matv(k,-1,l,2)', matv(k,-1,l,2)
cFC      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(1,-1,l,1)',matv(1,isig,isig2,1)
      Print*, '-4*matvb(k,-1,l,2)',-4*matb(1,isig,isig2,1)
      Print*, 'ratio',matv(1,isig,isig2,1)/(-4*matb(1,isig,isig2,1))
      print*
        Print*,'matv(2,-1,l,1)',matv(2,isig,isig2,1)
      Print*, '-4*matvb(2,-1,l,2)',-4*matb(2,isig,isig2,1)
      Print*, 'ratio',matv(2,isig,isig2,1)/(-4*matb(2,isig,isig2,1))
      print*

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
         do isig=-1,-1,2
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
         do isig=-1,-1,2
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
C      stop
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
             mutemp=mu!+1
            zww(mutemp,1) = mqza(mu)
                  aww(mutemp,1) = mqza(mu)
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
     1    zww(0,1),Myeps(0,isig2),isig,scale,1,3,
     2    7*comp,comp,boxcontr1,j5g(1,1,isig2,isig),boxcontr1b
     3    ,jbg(1,isig2,isig),Div)
       comp=-1
c KK Terms and gauge test
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    emqza,Myeps(0,isig2),isig,scale,1,3,
     2    3*comp,comp,nothing,j5gg(1,1,isig2,isig),nothing
     3    ,jbgg(1,isig2,isig),Div)
 
       call boxlineABETotal(p1,mqza,mpg,p2,psi(1,isig,2),psi(1,isig,1),
     1    zww(0,1),empg,isig,scale,1,3,
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
     1    aww(0,1),Myeps(0,isig2),isig,scale,1,3,
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
     1    Myeps(0,isig2),zww(0,1),isig,scale,1,2,
     2    7*comp,comp,boxcontr2,j5g(1,2,isig2,isig),boxcontr2b
     3    ,jbg(2,isig2,isig),Div)
       comp=-1
c KK Terms and gauge test
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    Myeps(0,isig2),emqza,isig,scale,1,2,
     2    3*comp,comp,nothing,j5gg(1,2,isig2,isig),nothing
     3    ,jbgg(2,isig2,isig),Div)
       call boxlineABETotal(p1,mpg,mqza,p2,psi(1,isig,2),psi(1,isig,1),
     1    empg,zww(0,1),isig,scale,1,2,
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
     1    Myeps(0,isig2),aww(0,1),isig,scale,1,2,
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
     1    zww(0,1),Myeps(0,isig2),isig,scale,20*comp,comp,
     2    boxcontrNoAbe1,jg5n(isig2,isig),Div)
      comp=-1
c For non-abelian gauge and kk/M terms          
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    emqza,Myeps(0,isig2),isig,scale,2*comp,comp,
     2    nothing,jg5ng(isig2,isig),Div)
      call boxlineNoAbeTotal(p1,mqza,p2,mpg,psi(1,isig,2),psi(1,isig,1),
     1    zww(0,1),empg,isig,scale,2*comp,comp,
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
     1    aww(0,1),Myeps(0,isig2),isig,scale,2*comp,comp,
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
cFC c        Print*, 'wardidtest1',wardidtest1
cFC c        Print*, 'wardidtest2',wardidtest2
cFC c        Print*, 'wardidtest3',wardidtest3
c      STOP
       wardidtest1=wardidtest1.and.wardidtest2.and.wardidtest3
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
        If(Div.ne.0) then
! Check factorization of Divergences. Set the color equal.
! So that I can check separately the CF and CA factorization
!        Print*,'CFCFCF',cfc
        cfc=1d0
        cacgraph=1d0
!        Print*,'CFCFCF',cfc
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
        mv1b =jbg(2,isig2,isig)  !g Z*
        mv2b =jbg(1,isig2,isig)   !Z* g
c
        mv3b =jbg1(2,isig2,isig)  !g A*
        mv4b =jbg1(1,isig2,isig)    !A* g

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1b =(jbgg(2,isig2,isig))*dotrc(mqza,zww(0,1))*zm2i(2) !g W^*
       mmv2b =(jbgg(1,isig2,isig))*dotrc(mqza,zww(0,1))*zm2i(2) !W^* g

c I set to zero these terms. There contribution start at the 13 digit when 
c comparing the the previous born which is checked against madgraph
       mmv1b=0d0
       mmv2b=0d0

       mmv3b =0d0 !g A^*
       mmv4b =0d0 !A^* g

C Change of sign to agree with Madgraph
        mv1b =-mv1b
        mmv1b =-mmv1b
        mv2b =-mv2b
        mmv2b =-mmv2b

        mv3b =-mv3b
        mmv3b =-mmv3b
        mv4b =-mv4b
        mmv4b =-mmv4b

       do k=1,2
       tmpb(k,isig,isig2,2)=(mv3b+mmv3b)*clr(ifl(1,k),1,isig)
       tmpb(k,isig,isig2,3)=(mv4b+mmv4b)*clr(ifl(1,k),1,isig)
       tmpb(k,isig,isig2,4)=(mv1b+mmv1b)*clr(ifl(1,k),2,isig)
       tmpb(k,isig,isig2,5)=(mv2b+mmv2b)*clr(ifl(1,k),2,isig)

       tmpb(k,isig,isig2,6)=tmpb(k,isig,isig2,2)+tmpb(k,isig,isig2,3)
       tmpb(k,isig,isig2,7)=tmpb(k,isig,isig2,4)+tmpb(k,isig,isig2,5)


       enddo
      
       ratioComp=
     &    abs(mat(1,isig,isig2,4)/tmpb(1,isig,isig2,2) -1d0)
     &   +abs(mat(1,isig,isig2,5)/tmpb(1,isig,isig2,3) -1d0)
     &   +abs(mat(1,isig,isig2,6)/tmpb(1,isig,isig2,4) -1d0)
     &   +abs(mat(1,isig,isig2,7)/tmpb(1,isig,isig2,5) -1d0)

      if(ratioComp.gt.1d-6) then

      print*
      print*
      print*, "ISIG",ISIG, "ISIG2",ISIG2
      print*
      print*
      Print*,"2: Born/b_box-1d0",ratioComp
      do k=1,2
      print*,"matb(1,isig,isig2,2)", tmpb(k,isig,isig2,2)
      print*,"mat(1,isig,isig2,4)", mat(k,isig,isig2,4)
      print*, tmpb(k,isig,isig2,2)/mat(k,isig,isig2,4)
      print* 
      print*,"matb(1,isig,isig2,3)", tmpb(k,isig,isig2,3)
      print*,"mat(1,isig,isig2,5)", mat(k,isig,isig2,5)
      print*,tmpb(k,isig,isig2,3)/mat(k,isig,isig2,5)
      print*
      print*,"matb(1,isig,isig2,4)", tmpb(k,isig,isig2,4)
      print*,"mat(1,isig,isig2,6)", mat(k,isig,isig2,6)
      print*, tmpb(k,isig,isig2,4)/mat(k,isig,isig2,6)
      print* 
      print*,"matb(1,isig,isig2,5)", tmpb(k,isig,isig2,5)
      print*,"mat(1,isig,isig2,7)", mat(k,isig,isig2,7)
      print*,tmpb(k,isig,isig2,5)/mat(k,isig,isig2,7)
      print*
      print*,"matb(1,isig,isig2,6)", tmpb(k,isig,isig2,6)
      print*,"mat(1,isig,isig2,8)", mat(k,isig,isig2,4)+mat(k,isig,isig2,5)
      print*, tmpb(k,isig,isig2,6)/(mat(k,isig,isig2,4)+mat(k,isig,isig2,5))
      print* 
      print*,"matb(1,isig,isig2,7)", tmpb(k,isig,isig2,7)
      print*,"mat(1,isig,isig2,8)", (mat(k,isig,isig2,6)+mat(k,isig,isig2,7))
      print*,tmpb(k,isig,isig2,7)/(mat(k,isig,isig2,6)+mat(k,isig,isig2,7))
      print*


      enddo
 
!      stop
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
        mv1 =j5g(1,2,isig2,isig)  !g Z*
        mv2 =j5g(1,1,isig2,isig)   !Z* g
c
        mv3 =j5g1(1,2,isig2,isig)  !g A*
        mv4 =j5g1(1,1,isig2,isig)    !A* g

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Additional KK terms
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
       mmv1 =(j5gg(1,2,isig2,isig))*dotrc(mqza,zww(0,1))*zm2i(2) !g W^*
       mmv2 =(j5gg(1,1,isig2,isig))*dotrc(mqza,zww(0,1))*zm2i(2) !W^* g
!       print*, "mmv1" ,mmv1,mv1 ,mmv1/mv1
!       print*, "mmv2" ,mmv2,mv2 ,mmv2/mv2

c I set to zero these terms. There contribution start at the 13 digit when 
c comparing the the previous born which is checked against madgraph for the 
c individual terms. When added, only affects the 15-16 digit.

       mmv1=0d0
       mmv2=0d0
c
       mmv3 =0d0 !g A^*
       mmv4 =0d0 !A^* g

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
C  No abelian part
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
        mv5=jg5n(isig2,isig)
        mv6=jg5n1(isig2,isig)

        mmv5=0d0
        mmv6=0d0

C Change of sign to agree with Madgraph
        mv1 =-mv1
        mmv1 =-mmv1
        mv2 =-mv2
        mmv2 =-mmv2
        mv3 =-mv3
        mmv3 =-mmv3
        mv4 =-mv4
        mmv4 =-mmv4
        mv5 =-mv5
        mmv5 =-mmv5
        mv6 =-mv6
        mmv6 =-mmv6

       do k=1,2
       tmpv(k,isig,isig2,2)=(mv3+mmv3+mv6+mmv6)*clr(ifl(1,k),1,isig)
       tmpv(k,isig,isig2,3)=(mv4+mmv4)*clr(ifl(1,k),1,isig)

       tmpv(k,isig,isig2,4)=(mv1+mmv1+mv5+mmv5)*clr(ifl(1,k),2,isig)
       tmpv(k,isig,isig2,5)=(mv2+mmv2)*clr(ifl(1,k),2,isig)

       tmpv(k,isig,isig2,6)=tmpv(k,isig,isig2,2)+tmpv(k,isig,isig2,3)
       tmpv(k,isig,isig2,7)=tmpv(k,isig,isig2,4)+tmpv(k,isig,isig2,5)

       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
        IF (wardidtest1 ) THEN
        
        IF (PerformGauge.gt.0) then
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

        do k=1,2
        Print*, 'mv2b',tmpb(k,isig,isig2,2)
        Print*, 'mv3b',tmpb(k,isig,isig2,3)
        Print*, 'Gauge',-tmpb(k,isig,isig2,2)/tmpb(k,isig,isig2,3)

        Print*, 'mv2',tmpv(k,isig,isig2,2)
        Print*, 'mv3',tmpv(k,isig,isig2,3)
        Print*, 'Gauge',-tmpv(k,isig,isig2,2)/tmpv(k,isig,isig2,3)

        Print*, 'mv4b',tmpb(k,isig,isig2,4)
        Print*, 'mv5b',tmpb(k,isig,isig2,5)
        Print*, 'Gauge',-tmpb(k,isig,isig2,4)/tmpb(k,isig,isig2,5)

        Print*, 'mv4',tmpv(k,isig,isig2,4)
        Print*, 'mv5',tmpv(k,isig,isig2,5)
        Print*, 'Gauge',-tmpv(k,isig,isig2,4)/tmpv(k,isig,isig2,5)

        enddo
        
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
c        Print*, 'Gauge',-mv2/mv4-1d0

        print*, 'j5a2g1(1,1,isig2,isig)',j5g1(1,1,isig2,isig)
        Print*, 'Gauge mv2',mv2!/j5a2g1(1,1,isig2,isig)

        print*, 'jg5a21(1,isig2,isig)',jg5n1(isig2,isig)
        Print*, 'Gauge mv5',mv4!/jg5a21(1,isig2,isig)
c        Print*, 'Gauge',-mv3/mv5-1d0
        
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

        mv1CA =j5g(2,2,isig2,isig)  !g Z*
        mv2CA =j5g(2,1,isig2,isig)   !Z* g
c
        mv3CA =j5g1(2,2,isig2,isig)  !g A*
        mv4CA =j5g1(2,1,isig2,isig)    !A* g

        mv1CA=-mv1CA
        mv2CA=-mv2CA
        mv3CA=-mv3CA
        mv4CA=-mv4CA
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

         do k=1,2
         ratioFact=((mv1+mv2))/(-2d0*(mv1b+mv2b))!tmpv(k,isig,isig2,6)/(-2d0*tmpb(k,isig,isig2,6))

         Print*, 'mv1',mv1+mv2!tmpv(k,isig,isig2,6)
         Print*, 'mv1b',mv1b+mv2b!tmpb(k,isig,isig2,6)
         Print*, "CF Div/Eps^2", ratioFact
         Print*, ""
         
         ratioFact=(mv3+mv4)/(-2d0*(mv3b+mv4b))!tmpv(k,isig,isig2,7)/(-2d0*tmpb(k,isig,isig2,7))

         Print*, 'mv1',mv3+mv4!tmpv(k,isig,isig2,7)
         Print*, 'mv1b',mv3B+mv4B!tmpb(k,isig,isig2,7)
         Print*, "CF Div/Eps^2", ratioFact
         Print*, ""
         enddo
         

         Print*, ""
         !Print*, 'SUM of born', mmv1b+mmv2b         
         
! CA part
         Print*, ""
         Print*, ""
         Print*, "CA PART"
         Print*, ""
         Print*, ""

         do k=1,2

         Print*, 'mv1CA',mv1CA*(-1d0/2d0)
         Print*, 'mvNoAbe',mv3
         Print*, 'mv1b',mv1b
         ratioFact=((mv1CA+mv2CA)*(-1d0/2d0)+mv5)/(-1d0*(mv1b+mv2b))
         Print*, "CA Div/Eps^2", ratioFact
         Print*, ""

          Print*, 'mv2CA',mv2CA*(-1d0/2d0)
         Print*, 'mvNoAbe',mv4
         Print*, 'mv2b',mv2b
         ratioFact=((mv3CA+mv4CA)*(-1d0/2d0)+mv6)/(-1d0*(mv3b+mv4b))
         Print*, "CA Div/Eps^2", ratioFact
         Print*, ""

         enddo
         
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
         
         ratioFact=(mv1+mv2)/((mv1b+mv2B)*cteEps1)
         Print*, 'mv1',mv1+MV2
         Print*, 'mv1b*cte',(mv1b+mv2b)*cteEps1
         Print*, "CF Div/Eps", ratioFact
         Print*, ""

         ratioFact=(mv3+mv4)/((mv3b+mv4b)*cteEps1)
         Print*, 'mv2',(mv3+mv4)
         Print*, 'mv2b*cte',(mv3b+mv4b)*cteEps1
         Print*, "CF Div/Eps", ratioFact
         Print*, ""

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CA part
         cteEps1=(-lns+lnt+lnu)
         
         Print*, "CA PART"
         Print*, ""
         Print*, 'cteEps1',cteEps1
         Print*, ""
         Print*, 'mv1CA',(mv1CA+mv2CA)*(-1d0/2d0)
         Print*, 'mvNoAbe',mv5
         Print*, 'mv1b*cte',(mv1b+mv2b)*cteEps1
         ratioFact=((mv1CA+mv2CA)*(-1d0/2d0)+mv5)/((mv1b+mv2b)*cteEps1)
         Print*, "CA Div/Eps", ratioFact
         Print*, ""
         Print*, ""
         Print*, 'mv2CA',(mv3CA+mv4CA)*(-1d0/2d0)
         Print*, 'mvNoAbe',(mv6)
         Print*, 'mv2b*cte',(mv3b+mv4b)*cteEps1
         ratioFact=((mv3CA+mv4CA)*(-1d0/2d0)+mv6)/((mv3b+mv4b)*cteEps1)
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

      IF (PerformGauge.ge.1) then
! then, tmpv(k,-1,l,1) should be zero
 
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(1,isg,isg2,6)', tmpv(1,isig,isig2,6)
      Print*,' matv(2,isg,isg2,6)', tmpv(2,isig,isig2,6)
      Print*,' matv(1,isg,isg2,7)', tmpv(1,isig,isig2,7)
      Print*,' matv(2,isg,isg2,7)', tmpv(2,isig,isig2,7)
      Print*, "matv", tmpv(2,isig,isig2,5)

      elseif(PerformGauge.gt.1) then
! then, matv(k,-1,l,1) with couplings corrected should be zero
      Print*,''
      Print*,'matv(k,-1,l,2) should be zero'
      Print*,' matv(k,-1,l,1)', tmpv(1,isig,isig2,2)
      Print*,'j5', (jg5n(isig2,isig)*clr(3,2,isig))
      Print*,''
      endif
      
      if (Div.eq.2) then
      Print*,''
      Print*,'The total contributions 1/eps^2'
      Print*, 'shoulb be -4*born since:'
      Print*, 'set CF=1 and CA=0 in ABE contribu'
      Print*,'matv(k,-1,l,6)',tmpv(1,isig,isig2,6)
      Print*, '-4*matvb(k,-1,l,6)',-4*tmpb(1,isig,isig2,6)
      Print*, 'ratio',tmpv(1,isig,isig2,6)/(-4*tmpb(1,isig,isig2,6))
      print*

      Print*,'matv(2,-1,l,6)',tmpv(2,isig,isig2,6)
      Print*, '-4*matvb(k,-1,l,6)',-4*tmpb(2,isig,isig2,6)
      Print*, 'ratio',tmpv(2,isig,isig2,6)/(-4*tmpb(2,isig,isig2,6))
      print*


      Print*,'matv(k,-1,l,7)',tmpv(1,isig,isig2,7)
      Print*, '-4*matvb(k,-1,l,7)',-4*tmpb(1,isig,isig2,7)
      Print*, 'ratio',tmpv(1,isig,isig2,7)/(-4*tmpb(1,isig,isig2,7))
      print*

      Print*,'matv(2,-1,l,7)',tmpv(2,isig,isig2,7)
      Print*, '-4*matvb(k,-1,l,7)',-4*tmpb(2,isig,isig2,7)
      Print*, 'ratio',tmpv(2,isig,isig2,7)/(-4*tmpb(2,isig,isig2,7))
      print*

      ENDIF

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
      tmpv(1,isig,isig2,6)=tmpv(1,isig,isig2,6)+tmpb(1,isig,isig2,6)*convfact
      tmpv(1,isig,isig2,6)=tmpv(1,isig,isig2,6)-tmpb(1,isig,isig2,6)*(K_cont)
      tmpv(2,isig,isig2,6)=tmpv(2,isig,isig2,6)+tmpb(2,isig,isig2,6)*convfact
      tmpv(2,isig,isig2,6)=tmpv(2,isig,isig2,6)-tmpb(2,isig,isig2,6)*(K_cont)

      tmpv(1,isig,isig2,7)=tmpv(1,isig,isig2,7)+tmpb(1,isig,isig2,7)*convfact
      tmpv(1,isig,isig2,7)=tmpv(1,isig,isig2,7)-tmpb(1,isig,isig2,7)*(K_cont)
      tmpv(2,isig,isig2,7)=tmpv(2,isig,isig2,7)+tmpb(2,isig,isig2,7)*convfact
      tmpv(2,isig,isig2,7)=tmpv(2,isig,isig2,7)-tmpb(2,isig,isig2,7)*(K_cont)

c copy to matv / matb
      matv(1,isig,isig2,2) = tmpv(1,isig,isig2,6)
      matv(2,isig,isig2,2) = tmpv(2,isig,isig2,6)
      matv(1,isig,isig2,3) = tmpv(1,isig,isig2,7)
      matv(2,isig,isig2,3) = tmpv(2,isig,isig2,7)
      matb(1,isig,isig2,2) = tmpb(1,isig,isig2,6)
      matb(2,isig,isig2,2) = tmpb(2,isig,isig2,6)
      matb(1,isig,isig2,3) = tmpb(1,isig,isig2,7)
      matb(2,isig,isig2,3) = tmpb(2,isig,isig2,7)


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
        do l1=6,7
       do k=1,2
         do isig=-1,1,2
           do isig2=1,2
      check2(k,isig,isig2,l1)=tmpv(k,isig,isig2,l1)+(K_contcheck+K_contcheck4)*tmpb(k,isig,isig2,l1)

      check3(k,isig,isig2,l1)=tmpv(k,isig,isig2,l1)
            enddo
          enddo
       enddo   
      enddo
      Reno_inv=2
      musq=1d0
      goto 100
      endif

      If (Reno_inv.gt.0) then
        do l1=6,7
         do k=1,2
         do isig=-1,1,2
           do isig2=1,2

      print*
      print*
      print*, "l1",l1,"K",k,'isig',isig,'isig2',isig2
      print*
       Print*, ' '
       Print*, 'chek1 hel 1 ',check3(k,isig,isig2,l1)
       Print*, 'chek1 hel 1 ',check2(k,isig,isig2,l1)

       Print*, ''
       Print*, 'check1 hel 1 ',tmpv(k,isig,isig2,l1)
       Print*, 'check1 hel 1 ',(tmpv(k,isig,isig2,l1)+(K_contcheck+K_contcheck4)*tmpb(k,isig,isig2,l1))
      
       Print*, 'ratio1',check2(k,isig,isig2,l1)/((tmpv(k,isig,isig2,l1)+(K_contcheck+K_contcheck4)*tmpb(k,isig,isig2,l1)))-1d0
       
       enddo
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      c FC
c      Print*, "Here"
c      stop
c FC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif ! end lbox
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C End Bosonic Virtual Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (lbox2) then

C Begin Fermionic Corrections
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c * gluon always attached to the loop (colour trace)
c??? * ggA-triangles are zero (Furry's theorem)
c??? * ggAA-boxes remaining
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
     
      do mu=0,3
c On-shell gluon
        mpg(mu) = -q(mu)
c Momentum W+
        minus_qp(mu) = -qp(mu,N)
c Momentum W- 
        minus_qm(mu) = -qm(mu,N)
c Higgs, Z, A*, momenta for triangles
        minus_qww(mu) = -qww(mu,N)
      enddo

      call curr6(1,psi(1,-1,2),p(0,2),psi(1,-1,1),p(0,1),curr2(0,-1))
      call propagate(1,curr2(0,-1),p12)

      do isig=-1,1,2
        l1 = (isig+3)/2
        do i=0,3
          curr1(i,l1) = curr2(i,isig)
       enddo
      enddo 
  
ccccccccccccccccccccccccccccccccc
c box loops
ccc first anchor point is always first gluon
ccc permutate second g and the two W over the other ones
ccc strong coupling is added later
cc g-p-m

  
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),
     &     1d0,1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),
     &     mf(j),mf(j),mf(j+1),mf(j),
     &     p12,mpg,minus_qp(0),minus_qm(0),
     &     curr1(0,l1),MyEps(0,l2),wp(1,N),wm(1,N),
     &     scale,comp,comp,rfbox(1,1,j,l1,l2),rgfbox(1,1,j,l1,l2),
     &     Div)
         comp=-1
       enddo
       enddo
       enddo
cc g-m-p
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingt1_gm(
     &     1d0,1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),
     &     1d0,1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),
     &     mf(j+1),mf(j+1),mf(j),mf(j+1),
     &     p12,mpg,minus_qm(0),minus_qp(0),
     &     curr1(0,l1),MyEps(0,l2),wm(1,N),wp(1,N),
     &     scale,comp,comp,rfbox(2,1,j,l1,l2),rgfbox(2,1,j,l1,l2),
     &     Div)
         comp=-1
       enddo
       enddo
       enddo
cc p-g-m
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),1d0,clr(cpl(j),3,+1),
     &     1d0,clr(cpl(j),3,-1),1d0,clr(cpl(j),3,-1),
     &     mf(j), mf(j+1), mf(j+1), mf(j),
     &     p12,minus_qp(0),mpg,minus_qm(0),
     &     curr1(0,l1),wp(1,N),MyEps(0,l2),wm(1,N),
     &     scale,comp,comp,rfbox(3,1,j,l1,l2),rgfbox(3,1,j,l1,l2),
     &     Div)
         comp=-1
       enddo
       enddo
       enddo
cc m-g-p
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),1d0,clr(cpl(j),3,+1),
     &     1d0,clr(cpl(j),3,-1),1d0,clr(cpl(j),3,-1),
     &     mf(j+1),mf(j),mf(j),mf(j+1),
     &     p12,minus_qm(0),mpg,minus_qp(0),
     &     curr1(0,l1),wm(1,N),MyEps(0,l2),wp(1,N),
     &     scale,comp,comp,rfbox(4,1,j,l1,l2),rgfbox(4,1,j,l1,l2),
     &     Div)
         comp=-1
       enddo
       enddo
       enddo
cc p-m-g
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),1d0,
     &     1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),1d0,
     &     mf(j),mf(j+1),mf(j),mf(j),
     &     p12,minus_qp(0),minus_qm(0),mpg,
     &     curr1(0,l1),wp(1,N),wm(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(5,1,j,l1,l2),rgfbox(5,1,j,l1,l2),
     &     Div)
         comp=-1
       enddo
       enddo
       enddo
cc m-p-g
       do j=1,3,2
       comp=1
       do l1=1,2
       do l2=1,2
         call fermionboxcouplingT1_GM(
     &     1d0,clr(cpl(j),3,+1),clr(cpl(j),3,+1),1d0,
     &     1d0,clr(cpl(j),3,-1),clr(cpl(j),3,-1),1d0,
     &     mf(j+1),mf(j),mf(j+1),mf(j+1),
     &     p12,minus_qm(0),minus_qp(0),mpg,
     &     curr1(0,l1),wm(1,N),wp(1,N),MyEps(0,l2),
     &     scale,comp,comp,rfbox(6,1,j,l1,l2),rgfbox(6,1,j,l1,l2),
     &     Div)
         comp=-1
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
c add everything together
       if (wardidtest1) then
         box2count2m = box2count2m+1
         box2count2ps(ps_number) = box2count2ps(ps_number)+1   
ccc         write(*,*) "MARKER3: gauge test ok"
         do l1=1,2
         isig = 2*l1-3
         do l2=1,2
           mv1 = 0d0
           do i=1,6
           do j=1,1
             mv1 = mv1 + 2*rfbox(i,j,1,l1,l2) ! u,c
     &                 + 1*rfbox(i,j,3,l1,l2) ! t
           enddo
           enddo
c and multiply with color factor TR
           matv(1,isig,l2,4) = (-mv1) * TR 
           matv(2,isig,l2,4) = matv(1,isig,l2,4)
         enddo
         enddo
       else !wardidtest
         box2countm = box2countm+1
         box2countps(ps_number) = box2countps(ps_number)+1 
       endif
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(lvert) then

cc 1-2
       do j=3,4 ! no contribution from massless loops
       comp=1
       do l1=1,2
       do l2=1,2
ccc Z 
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     p12,mpg,minus_qww(0),
     &     curr1(0,l1),MyEps(0,l2),zww(0,N),
     &     scale,comp,comp,rfvert(1,1,j,l1,l2),rgfvert(1,1,j,l1,l2),
     &     Div)
         comp=+1
cccc A -- no contribution due to Furry's theorem
c         call fermionvertcouplingT(
c     &     1d0,1d0,clr(cpl(j),1,+1),
c     &     1d0,1d0,clr(cpl(j),1,-1),
c     &     mf(j),
c     &     p12,mpg,minus_qww(0),
c     &     curr1(0,l1),MyEps(0,l2),aww(0,N),
c     &     scale,comp,comp,rfvert(1,2,j,l1,l2),rgfvert(1,2,j,l1,l2),
c     &     Div)
       enddo
       enddo
       enddo
cc 2-1
       do j=3,4 ! no contribution from massless loops
       comp=1
       do l1=1,2
       do l2=1,2
ccc Z 
         call fermionvertcouplingT(
     &     1d0,1d0,clr(cpl(j),2,+1),
     &     1d0,1d0,clr(cpl(j),2,-1),
     &     mf(j),
     &     mpg,p12,minus_qww(0),
     &     MyEps(0,l2),curr1(0,l1),zww(0,1),
     &     scale,comp,comp,rfvert(2,1,j,l1,l2),rgfvert(2,1,j,l1,l2),
     &     Div)
         comp=+1
cccc A -- no contribution due to Furry's theorem
c         call fermionvertcouplingT(
c     &     1d0,1d0,clr(cpl(j),1,+1),
c     &     1d0,1d0,clr(cpl(j),1,-1),
c     &     mf(j),
c     &     mpg,p12,minus_qww(0),
c     &     MyEps(0,l2),curr1(0,l1),aww(0,1),
c     &     scale,comp,comp,rfvert(2,2,j,l1,l2),rgfvert(2,2,j,l1,l2),
c     &     Div)
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
         isig = 2*l1-3
         do l2=1,2
           mv1 = 0d0
           mv2 = 0d0
           do i=1,2
           do j=1,1
             mv1 = mv1 + rfvert(i,j,3,l1,l2) ! t
             mv2 = mv2 + rfvert(i,j,4,l1,l2) ! b
           enddo
           enddo
c and multiply with color factor TR
           matv(1,isig,l2,5) = - (mv1 + mv2) * TR 
           matv(2,isig,l2,5) = matv(1,isig,l2,5)
         enddo
         enddo
       else !wardidtest
         box2countm = box2countm+1
         box2countps(ps_number) = box2countps(ps_number)+1 
       endif

       endif ! lvert

       if (lhiggs) then
c vertex diagrams
ccc strong coupling is added later
ccc triangles with Z/photon vanish, for Higgs only top relevant
ccc Both fermion-line directions identical - factor 2 added later
         do j=3,4 ! no contribution from ~ massless loops
         do l1=1,2 
         do l2=1,2 
           call gluonverHiggsEvenT(
     &       mf(j),
     &       p12,mpg,minus_qww,
     &       curr1(0,l1),MyEps(0,l2),
     &       scale,1,rfvert(1,3,j,l1,l2),
     &       Div)
         enddo
         enddo
         enddo

c add everything together
         do l1=1,2
         do l2=1,2

         if (l1.eq.1) then
           poln=-1
         else
           poln=1
         endif

         matv(1,poln,l2,6) = 0d0
         do j=3,4 ! no contribution from ~ massless loops
           matv(1,poln,l2,6) = matv(1,poln,l2,5) 
     &           + rfvert(1,3,j,l1,l2) ! *2*TR ! amplitude (factor 2 for both fermion directions) * colour factor
     &           * mf(j)/vv    ! top Yukawa
     &           * w4l(h,1)/dcmplx(qww(4,1)-xm2(6),xmg(6)) ! decay matrix element * Higgs propagator
         enddo

         matv(2,poln,l2,6) = matv(1,poln,l2,6) 

         enddo
         enddo
       endif

c       print*, "MH",qaa(4,1)
c       print*, "psq", dotrr(minus_qww,minus_qww)


c      stop


#ifdef WITH_LOOPTOOLS
        if (ldebuglt) then
        fcflags = 2**0+2**1
        fchelicities = 0
     &   + (2**3+2**1) * thirtytwo**6 ! q
     &   + (2**3+2**1) * thirtytwo**5 ! q
     &   + (     2**1) * thirtytwo**4 ! nu_e
     &   + (2**3     ) * thirtytwo**3 ! e+
     &   + (     2**1) * thirtytwo**2 ! mu-
     &   + (2**3     ) * thirtytwo**1 ! nu_mu
     &   + (2**3+2**1) * thirtytwo**0 ! g
        call fafcWWj_SquaredME(fcresult, fcamp, fchelicities, fcflags)
c convert to kartesian basis 
        do isig=-1,1,2
        do i=0,1
          fckamp(i,1) = -im*(fcamp(i,isig,-1) + fcamp(i,isig,1))/sqrt(2d0)
          fckamp(i,2) = (fcamp(i,isig,-1) - fcamp(i,isig,1))/sqrt(2d0)
        enddo
c divide out g_s and typical loop factor 16*pi**2 
        do l=1,2
c            write(*,*) "VBFNLO tree", l, isig,
c     &               mat(1,isig,l,8)+mat(1,isig,l,9),
c     &           abs(mat(1,isig,l,8)+mat(1,isig,l,9)), 
c     &              (mat(1,isig,l,8)+mat(1,isig,l,9))*
c     &        dconjg(mat(1,isig,l,8)+mat(1,isig,l,9))
c            write(*,*) "FAFC tree", l, isig,
c     &               fckamp(0,l)/(fpials(2)**(1/2d0)), 
c     &           abs(fckamp(0,l)/(fpials(2)**(1/2d0))), 
c     &              (fckamp(0,l)/(fpials(2)**(1/2d0)))*
c     &        dconjg(fckamp(0,l)/(fpials(2)**(1/2d0)))
            write(*,*) "VBFNLO loop", l, isig,
     &               matv(1,isig,l,4)+matv(1,isig,l,5)+matv(1,isig,l,6),
     &           abs(matv(1,isig,l,4)+matv(1,isig,l,5)+matv(1,isig,l,6)), 
     &         dble((matv(1,isig,l,4)+matv(1,isig,l,5)+matv(1,isig,l,6))*
     &        dconjg(mat(1,isig,l,8)+mat(1,isig,l,9)))
            write(*,*) "FAFC loop", l, isig,
     &               fckamp(1,l)/(fpials(2)**(3/2d0)/(16d0*pi**2)), 
     &           abs(fckamp(1,l)/(fpials(2)**(3/2d0)/(16d0*pi**2))), 
     &         dble((fckamp(1,l)/(fpials(2)**(3/2d0)/(16d0*pi**2)))*
     &        dconjg(fckamp(0,l)/(fpials(2)**(1/2d0))))
            write(*,*) "Ratio tree", l, isig,
     &              (mat(1,isig,l,8)+mat(1,isig,l,9))*
     &        dconjg(mat(1,isig,l,8)+mat(1,isig,l,9))/
     &       (      (fckamp(0,l)/(fpials(2)**(1/2d0)))*
     &        dconjg(fckamp(0,l)/(fpials(2)**(1/2d0))))
            write(*,*) "Ratio loop", l, isig,
     &         dble((matv(1,isig,l,4)+matv(1,isig,l,5)+matv(1,isig,l,6))*
     &        dconjg(mat(1,isig,l,8)+mat(1,isig,l,9)))/
     &       ( dble((fckamp(1,l)/(fpials(2)**(3/2d0)/(16d0*pi**2)))*
     &        dconjg(fckamp(0,l)/(fpials(2)**(1/2d0)))))
          enddo
c          do l=1,2
c            write(*,*) "TTT3", l, fcktreeamp(l)/(fpials(2)**(1/2d0))/mm(k,-1,l),
c     &        abs(fcktreeamp(l)/(fpials(2)**(1/2d0))/mm(k,-1,l)),
c     &        (fckamp(l)/(fpials(2)**(3/2d0))*16d0*pi**2)*
c     &        dconjg(fcktreeamp(l)/(fpials(2)**(1/2d0)))/(
c     &                      matv(k,-1,l,4)*dconjg(mm(k,-1,l)))
c            write(*,*) "RRR3", l, fckamp(l)/(fpials(2)**(3/2d0))*16d0*pi**2/matv(k,-1,l,4), 
c     &        abs(fckamp(l)/(fpials(2)**(3/2d0))*16d0*pi**2/matv(k,-1,l,4))
c          enddo
        enddo
        read(*,*)
        endif
#endif


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      endif !lbox


      end if ! end lnlo


#endif


c sum the graphs, square them and map them onto uuww and ddww
      do k = 1,2
         res(k) = 0
         resv(k)=0
         do isig = -1,1,2
             do l = 1,2
               mm(k,isig,l) = 0
               mmb(k,isig,l) = 0
               do i = 1,7
                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,i)
               enddo
               do i=1,3
                  mmb(k,isig,l) = mmb(k,isig,l) +  matb(k,isig,l,i)
c                  print*, "mmb(here)", matb(k,isig,l,i)
               enddo

               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2

            enddo
         enddo
c         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor     

      enddo
c      stop
c      uuww = res(1)
c      ddww = res(2)


         do l=1,2
            do isig=-1,1,2
               do isig2=1,2
                 mv(l,isig,isig2)=0d0
               enddo
            enddo

            do isig=-1,1,2
               do isig2=1,2 
                  do k=1,6
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
     1                (mv(l,isig,isig2) + mm(l,isig,isig2)*K_cont) 

                 elseif(nlo.eq.7) then
                  mv(l,isig,isig2) = als(1,N)*oneo4pi*
     1                (mv(l,isig,isig2) + mm(l,isig,isig2)*K_cont) 
     2                 +factortoadd*mmb(l,isig,isig2)

                 else if(nlo.lt.0) then
                  mv(l,isig,isig2) = als(1,N)*oneo4pi*
     1                (mv(l,isig,isig2)) 
     2                 +factortoadd*mmb(l,isig,isig2)
c         PRINT*, L, ISIG, ISIG2
c         print*, "mmb()",mmb(l,isig,isig2)
c         print*, "mm()",mm(l,isig,isig2)
c         PRINT*
                  else
                     mv(l,isig,isig2) = (0d0, 0d0)
                endif

        resv(l) = resv(l) + 2d0*dreal(mm(l,isig,isig2)*conjg(mv(l,isig,isig2)))
c experimental: add loop-squared for fermion loop
cc        if (nlo .eq. -6) then
cc          resv(l) = resv(l) + dreal(mv(l,isig,isig2))**2 
cc     &                      + dimag(mv(l,isig,isig2))**2
cc        endif
        
        enddo
        enddo


           if (nlo.ge.0) then
c              print*, "HERE 2", nlo, res(l)

             res(l) = (res(l)+resv(l))*4d0*fpials(2)  ! C_2*3 is the color factor
            else
             res(l) = resv(l)*4d0*fpials(2)           ! C_2*3 is the color factor
           endif  

        enddo


      uuww = res(1)
      ddww = res(2)

      return
      end


