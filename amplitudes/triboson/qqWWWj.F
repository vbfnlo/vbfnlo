


c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c
c        Stefan Prestel : sprestel@particle.uni-karlsruhe.de
c        Initial version: February 2008
c
c        This file contains subroutines to calculate matrix elememts**2 for:
c
c     q1 q2~ --> W+W-W+g, W+(1) ---> f1-bar f2, W-(1) ---> f3-bar f4, W+(2) ---> f5-bar f6
c     subroutine qqWWWpj
c
c     q1 q2~ --> W-W+W-g, W-(1) ---> f1-bar f2, W+(1) ---> f3-bar f4, W-(2) ---> f5-bar f6
c     subroutine qqWWWmj
c
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------


c****************subroutine QQWWWpJ ********************************************
C
C  QQWWWpj calculates the matrix elements**2 for 
C        q1 q2~ ----> W+ W- W+ g,   
C         W+ ---> f1-bar f2, W- ---> f3-bar f4, Z ---> f5-bar f6
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
c      return uuwwz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       helicity of partons 1 and 2 is fixed
c
c*******************************************************************************

      subroutine qqWWWpj(pbar,sign, qbar, gsign, qqbarWWWj,NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c     electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program

      double precision fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision res, qqbarWWWj
      double precision p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i
      integer  NW
      integer  lh, l

c internal variables

c momenta of particles and different internal lines
      double precision fqp1(0:4,2), fqp2(0:4,2)
      double precision pq(0:4,2)
      double precision pgp1(0:4,2),pgp2(0:4,2),pp1g(0:4,2),pp2g(0:4,2)
      double precision pwp1wm1(0:4), pwp2wm1(0:4)

c wavefunctions
      double complex psi(2,-1:1,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psiwp1(2,-1:1,2),psiwp2(2,-1:1,2)
      double complex psig(2,-1:1,2,2)
      double complex psigp1(2,-1:1,2,2),psigp2(2,-1:1,2,2),
     1               psip1g(2,-1:1,2,2),psip2g(2,-1:1,2,2)

c currents
      double complex jgqqbar(0:5,-1:1,2),jqqbarg(0:5,-1:1,2)
      double complex jgp1z(0:5,-1:1,2,2),jgp2z(0:5,-1:1,2,2),
     1               jp1zg(0:5,-1:1,2,2),jp2zg(0:5,-1:1,2,2),
     2               jp1gz(0:5,-1:1,2,2),jp2gz(0:5,-1:1,2,2)

c matrix variables
      double complex mat(2,6), madgr(6,-1:1)
      double complex mm(2)
      double complex m1kb(6), m2kb(6), temp, temp1
      double complex m1(3,2)


      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr
      double complex czero 
      external dotcc, s1c, dotrc, dotrr
      parameter (im = (0d0,1d0))
      parameter (czero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i

      logical ldebugs, ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebugs=.false.)

      integer N
      parameter( N = 1)

c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,1)    !als(1,1) : NLO alpha_s if NLO.ne.0
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas       !alfas : LO alpha_s

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do l = 1,2
         do i = 1,6
            mat(l,i) = czero
         enddo
      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
         pwp1wm1(mu) = qp1(mu,N) + q1m(mu,N)
         pwp2wm1(mu) = qp2(mu,N) + q1m(mu,N)
      enddo
      pwp1wm1(4) = dotrr(pwp1wm1(0),pwp1wm1(0))
      pwp2wm1(4) = dotrr(pwp2wm1(0),pwp2wm1(0))

c------------------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp1(0,N),wp1(1,N),psiwp1(1,-1,2),fqp1(0,2))     ! <q~W+(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp2(0,N),wp2(1,N),psiwp2(1,-1,2),fqp2(0,2))     ! <q~W+(2)|

      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp1(0,N),wp1(1,N),psiwp1(1,-1,1),fqp1(0,1))     ! |W+(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp2(0,N),wp2(1,N),psiwp2(1,-1,1),fqp2(0,1))     ! |W+(2)q>

      do l =1,2                        !two gluon polarisations

         call polvec(qbar(0),l,eps(0,l))

         call bra2r(psi(1,-1,2),.true.,p(0,2),-1,q(0),eps(0,l),psig(1,-1,2,l),pq(0,2))        !<q~g|_l
         call ket2r(psi(1,-1,1),.true.,p(0,1),-1,q(0),eps(0,l),psig(1,-1,1,l),pq(0,1))        !|g q>_l

c get other bra and ket vectors

         call bra2c(psig(1,-1,2,l),.false.,pq(0,2),-1,qp1(0,N),wp1(1,N),psigp1(1,-1,2,l),pgp1(0,2))  ! <q~gW+(1)|
         call bra2c(psig(1,-1,2,l),.false.,pq(0,2),-1,qp2(0,N),wp2(1,N),psigp2(1,-1,2,l),pgp2(0,2))  ! <q~gW+(2)|

         call bra2r(psiwp1(1,-1,2),.false.,fqp1(0,2),-1,q(0),eps(0,l),psip1g(1,-1,2,l),pp1g(0,2))    ! <q~W+(1)g|
         call bra2r(psiwp2(1,-1,2),.false.,fqp2(0,2),-1,q(0),eps(0,l),psip2g(1,-1,2,l),pp2g(0,2))    ! <q~W+(2)g|

         call ket2c(psig(1,-1,1,l),.false.,pq(0,1),-1,qp1(0,N),wp1(1,N),psigp1(1,-1,1,l),pgp1(0,1))  ! |W+(1)gq>
         call ket2c(psig(1,-1,1,l),.false.,pq(0,1),-1,qp2(0,N),wp2(1,N),psigp2(1,-1,1,l),pgp2(0,1))  ! |W+(2)gq>

         call ket2r(psiwp1(1,-1,1),.false.,fqp1(0,1),-1,q(0),eps(0,l),psip1g(1,-1,1,l),pp1g(0,1))    ! |gW+(1)q>
         call ket2r(psiwp2(1,-1,1),.false.,fqp2(0,1),-1,q(0),eps(0,l),psip2g(1,-1,1,l),pp2g(0,1))    ! |gW+(2)q>

c get currents for contraction with W->WWW leptonic tensor

         call curr6(-1,psi(1,-1,2)   ,p(0,2) ,psig(1,-1,1,l),pq(0,1),jgqqbar(0,-1,l))      ! gz : <q~|gamma_mu|g q>_l
         call curr6(-1,psig(1,-1,2,l),pq(0,2),psi(1,-1,1)   ,p(0,1) ,jqqbarg(0,-1,l))      ! zg : <q~g|gamma_mu| q>_l

c get currents for contraction with A/Z->WW leptonic tensor

c currents with W+(1)
         call curr6(-1,psi(1,-1,2),p(0,2),psigp1(1,-1,1,l),pgp1(0,1),jgp1z(0,-1,1,l))      ! gwz: <q~|gamma_mu|W+(1)gq>
         call curr6(-1,psiwp1(1,-1,2),fqp1(0,2),psig(1,-1,1,l),pq(0,1),jgp1z(0,-1,2,l))    ! gzw: <q~W+(1)|gamma_mu|gq>

         call curr6(-1,psig(1,-1,2,l),pq(0,2),psiwp1(1,-1,1),fqp1(0,1),jp1zg(0,-1,1,l))    ! wzg: <q~g|gamma_mu|W+(1)q>
         call curr6(-1,psigp1(1,-1,2,l),pgp1(0,2),psi(1,-1,1),p(0,1),jp1zg(0,-1,2,l))      ! zwg: <q~gW+(1)|gamma_mu|q>

         call curr6(-1,psi(1,-1,2),p(0,2),psip1g(1,-1,1,l),pp1g(0,1),jp1gz(0,-1,1,l))      ! wgz: <q~|gamma_mu|gW+(1)q>
         call curr6(-1,psip1g(1,-1,2,l),pp1g(0,2),psi(1,-1,1),p(0,1),jp1gz(0,-1,2,l))      ! zgw: <q~W+(1)g|gamma_mu|q>

c currents with W+(2)
         call curr6(-1,psi(1,-1,2),p(0,2),psigp2(1,-1,1,l),pgp2(0,1),jgp2z(0,-1,1,l))      ! gwz: <q~|gamma_mu|W+(2)gq>
         call curr6(-1,psiwp2(1,-1,2),fqp2(0,2),psig(1,-1,1,l),pq(0,1),jgp2z(0,-1,2,l))    ! gzw: <q~W+(2)|gamma_mu|gq>

         call curr6(-1,psig(1,-1,2,l),pq(0,2),psiwp2(1,-1,1),fqp2(0,1),jp2zg(0,-1,1,l))    ! wzg: <q~g|gamma_mu|W+(2)q>
         call curr6(-1,psigp2(1,-1,2,l),pgp2(0,2),psi(1,-1,1),p(0,1),jp2zg(0,-1,2,l))      ! zwg: <q~gW+(2)|gamma_mu|q>

         call curr6(-1,psi(1,-1,2),p(0,2),psip2g(1,-1,1,l),pp2g(0,1),jp2gz(0,-1,1,l))      ! wgz: <q~|gamma_mu|gW+(2)q>
         call curr6(-1,psip2g(1,-1,2,l),pp2g(0,2),psi(1,-1,1),p(0,1),jp2gz(0,-1,2,l))      ! zgw: <q~W+(2)g|gamma_mu|q>

      enddo                            !gluon polarisation summation

C------------------------------------------------------------------------------
C  CALCULATION OF UD~ --> W+W-W+g --> LEPTONS+g PROCESSES (TREE LEVEL)

c      U D~ --> W+W-W+g --> LEPTONS+g
      do l=1,2
         m1(1,l) =   !W+(1) W- W+(2) + g
     &      -s1c(psiwp2(1,-1,2),w1m(1,1),.true.,-1,psigp1(1,-1,1,l)) !g W+(1)W-(1)W+(2)
     &      -s1c(psiwp2(1,-1,2),w1m(1,1),.true.,-1,psip1g(1,-1,1,l)) !W+(1) g W-(1)W+(2)
     &      -s1c(psip2g(1,-1,2,l),w1m(1,1),.true.,-1,psiwp1(1,-1,1)) !W+(1)W-(1) g W+(2)
     &      -s1c(psigp2(1,-1,2,l),w1m(1,1),.true.,-1,psiwp1(1,-1,1)) !W+(1)W-(1)W+(2) g

         m1(2,l) =   !W+(2) W- W+(1) + g
     &      -s1c(psiwp1(1,-1,2),w1m(1,1),.true.,-1,psigp2(1,-1,1,l)) !g W+(2)W-(1)W+(1)
     &      -s1c(psiwp1(1,-1,2),w1m(1,1),.true.,-1,psip2g(1,-1,1,l)) !W+(2) g W-(1)W+(1)
     &      -s1c(psip1g(1,-1,2,l),w1m(1,1),.true.,-1,psiwp2(1,-1,1)) !W+(2)W-(1) g W+(1)
     &      -s1c(psigp1(1,-1,2,l),w1m(1,1),.true.,-1,psiwp2(1,-1,1)) !W+(2)W-(1)W+(1) g

         mat(l,1) = m1(1,l)*clr(3,3,-1)**3
     &            + m1(2,l)*clr(3,3,-1)**3         ! !M(UD~->W+W-W+g ->Lg)
         m1(1,l) = czero
         m1(2,l) = czero
      enddo

C CALCULATION OF UD~ --> W+ V g --> LEPTONS+g PROCESSES (TREE LEVEL)

c U D~ --> ZW+g --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_zww(0,1,2),jgp1z(0,-1,2,l))*clr(3,2,-1)     !g Z W+(1)
     &             -dotcc(WWW_zww(0,1,2),jp1gz(0,-1,2,l))*clr(3,2,-1)     !Z g W+(1)
     &             -dotcc(WWW_zww(0,1,2),jp1zg(0,-1,2,l))*clr(3,2,-1)     !Z W+(1) g

         m1(2,l) = -dotcc(WWW_zww(0,1,1),jgp2z(0,-1,2,l))*clr(3,2,-1)     !g Z W+(2)
     &             -dotcc(WWW_zww(0,1,1),jp2gz(0,-1,2,l))*clr(3,2,-1)     !Z g W+(2)
     &             -dotcc(WWW_zww(0,1,1),jp2zg(0,-1,2,l))*clr(3,2,-1)     !Z W+(2) g

c      additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         m1(3,l) = zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jgp1z(0,-1,2,l)) !g Z W+(1)
     &           + zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jp1gz(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jp1zg(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jgp2z(0,-1,2,l)) !g Z W+(2)
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jp2gz(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jp2zg(0,-1,2,l))
         m1(3,l) = m1(3,l)*clr(3,2,-1)

         mat(l,2) = (m1(1,l) + m1(2,l) +m1(3,l))*CLR(3,3,-1)              !M(UD~->ZW+g ->Lg)
         m1(1,l) = czero
         m1(2,l) = czero
         m1(3,l) = czero
      enddo

c U D~ --> AW+g --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_aww(0,1,2),jgp1z(0,-1,2,l))*clr(3,1,-1)     !g A W+(1)
     &             -dotcc(WWW_aww(0,1,2),jp1gz(0,-1,2,l))*clr(3,1,-1)     !A g W+(1)
     &             -dotcc(WWW_aww(0,1,2),jp1zg(0,-1,2,l))*clr(3,1,-1)     !A W+(1) g

         m1(2,l) = -dotcc(WWW_aww(0,1,1),jgp2z(0,-1,2,l))*clr(3,1,-1)     !g A W+(2)
     &             -dotcc(WWW_aww(0,1,1),jp2gz(0,-1,2,l))*clr(3,1,-1)     !A g W+(2)
     &             -dotcc(WWW_aww(0,1,1),jp2zg(0,-1,2,l))*clr(3,1,-1)     !A W+(2) g

         mat(l,3) = (m1(1,l) + m1(2,l))*clr(3,3,-1)              !M(UD~->AW+g ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
      enddo

c U D~ --> W+Zg --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_zww(0,1,2),jgp1z(0,-1,1,l))*clr(4,2,-1)     !g W+(1) Z
     &             -dotcc(WWW_zww(0,1,2),jp1gz(0,-1,1,l))*clr(4,2,-1)     !W+(1) g Z
     &             -dotcc(WWW_zww(0,1,2),jp1zg(0,-1,1,l))*clr(4,2,-1)     !W+(1) Z g

         m1(2,l) = -dotcc(WWW_zww(0,1,1),jgp2z(0,-1,1,l))*clr(4,2,-1)     !g W+(2) Z
     &             -dotcc(WWW_zww(0,1,1),jp2gz(0,-1,1,l))*clr(4,2,-1)     !W+(2) g Z
     &             -dotcc(WWW_zww(0,1,1),jp2zg(0,-1,1,l))*clr(4,2,-1)     !W+(2) Z g

c      additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         m1(3,l) = zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jgp1z(0,-1,1,l)) !g W+(1) Z
     &           + zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jp1gz(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp2wm1(0),WWW_zww(0,1,2))*dotrc(pwp2wm1(0),jp1zg(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jgp2z(0,-1,1,l)) !g W+(2) Z
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jp2gz(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jp2zg(0,-1,1,l))
         m1(3,l) = m1(3,l)*clr(4,2,-1)

         mat(l,4) = (m1(1,l) + m1(2,l) +m1(3,l))*clr(3,3,-1)              !M(UD~->W+Zg ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
         m1(3,l) = czero
      enddo

c U D~ --> W+Ag --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_aww(0,1,2),jgp1z(0,-1,1,l))*clr(4,1,-1)     !g W+(1) A
     &             -dotcc(WWW_aww(0,1,2),jp1gz(0,-1,1,l))*clr(4,1,-1)     !W+(1) g A
     &             -dotcc(WWW_aww(0,1,2),jp1zg(0,-1,1,l))*clr(4,1,-1)     !W+(1) A g

         m1(2,l) = -dotcc(WWW_aww(0,1,1),jgp2z(0,-1,1,l))*clr(4,1,-1)     !g W+(2) A
     &             -dotcc(WWW_aww(0,1,1),jp2gz(0,-1,1,l))*clr(4,1,-1)     !W+(2) g A
     &             -dotcc(WWW_aww(0,1,1),jp2zg(0,-1,1,l))*clr(4,1,-1)     !W+(2) A g

         mat(l,5) = (m1(1,l) + m1(2,l))*clr(3,3,-1)              !M(UD~->W+Ag ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
      enddo

C CALCULATION OF UD~ --> W+W-W+g --> LEPTONS+g PROCESSES (TREE LEVEL)

c U D~ --> W+g --> LEPTONS+g
      do l = 1,2

         m1(1,l) = -dotcc(wptowww(0,1),jgqqbar(0,-1,l))               ! GLUON EMISSION OFF U
     &             -dotcc(wptowww(0,1),jqqbarg(0,-1,l))               ! GLUON EMISSION OFF D~

         mat(l,6) = - m1(1,l)*clr(3,3,-1)
         m1(1,l) = czero
      enddo

C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do i = 1,6
            m1kb(i) = czero ! kartesian basis l = 1
            m2kb(i) = czero ! kartesian basis l = 2
            do lh = -1,1,2
               madgr(i,lh) = czero
            enddo
         enddo

         do lh = -1,1,2
            do i =1,6
               madgr(i,lh) = WWWtampr(i,-1,lh,1)       ! for udbar process
c               madgr(i,lh) = WWWtampr(i,-1,lh,2)       ! for dubar process
            enddo
         enddo

         do i=1,6
            m1kb(i) =     (madgr(i,-1) - madgr(i,1))/sqrt(2d0)
            m2kb(i) =-im*(madgr(i,-1) + madgr(i,1))/sqrt(2d0)
         enddo

c print results from checks:

         do i=1,6
            print*,' '
            print*,'ratio mat/m1kb (for l=1),for i = ',i,':',
     &         '(', dreal(mat(1,i))/dreal(m1kb(i)) ,
     &              dimag(mat(1,i))/dimag(m1kb(i)) , ')'
            print*,'mat(1,',i,') = ',mat(1,i)
            temp = dreal(mat(1,i))**2+dimag(mat(1,i))**2
            print*,'temp = ',temp
            print*,'m1kb(',i,') =  ',m1kb(i)
            temp1 = dreal(m1kb(i))**2+dimag(m1kb(i))**2
            print*,'temp1 = ',temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '

            print*,'ratio mat/m2kb (for l=1),for i = ',i,':',
     &         '(', dreal(mat(2,i))/dreal(m2kb(i)) ,
     &              dimag(mat(2,i))/dimag(m2kb(i)), ')'
            print*,'mat(2,',i,') = ',mat(2,i)
            temp = dreal(mat(2,i))**2+dimag(mat(2,i))**2
            print*,'temp = ',temp
            print*,'m2kb(',i,') =  ',m2kb(i)
            temp1 = dreal(m2kb(i))**2+dimag(m2kb(i))**2
            print*,'temp1 = ',temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '
         enddo

      endif ! end of debugging block
#endif

c------------------------------------------------------------------------------

c sum the graphs, square them and map them onto qqbarwwwj

      res = 0.0d0
      do l = 1,2
         mm(l) = 0
         do i = 1,6
            mm(l) = mm(l) +  mat(l,i)
         enddo
         res = res + dreal(mm(l))**2 + dimag(mm(l))**2
      enddo
      res = res*4d0*fpials(2)  ! C_2*3 is the color factor

      qqbarWWWj = res

      return
      end

c*******************************************************************************
c  end of subroutine qqWWWpj
c*******************************************************************************


c****************subroutine QQWWWmJ ********************************************
C
C  QQWWWpj calculates the matrix elements**2 for 
C        q1 q2~ ----> W- W+ W- g,
C         W- ---> f1-bar f2, W+ ---> f3-bar f4, W- ---> f5-bar f6
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
c      return uuwwz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       helicity of partons 1 and 2 is fixed
c
c*******************************************************************************

      subroutine qqWWWmj(pbar,sign, qbar, gsign, qqbarWWWj,NW)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c electroweak couplings are taken from KOPPLN
      double precision  clr, xm2, xmg, b, v, a
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),
     1                  V(4,5),A(4,5)

c variables for the main part of the program
      double precision fpials(2:3), fpi
      parameter (fpi=4d0*pi)

      double precision pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision res, qqbarWWWj
      double precision p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i
      integer  NW
      integer  lh, l

c internal variables

c momenta of particles and different internal lines
      double precision fqm1(0:4,2), fqm2(0:4,2)
      double precision pq(0:4,2)
      double precision pgm1(0:4,2),pgm2(0:4,2),pm1g(0:4,2),pm2g(0:4,2)
      double precision pwp1wm1(0:4), pwp1wm2(0:4)

c wavefunctions
      double complex psi(2,-1:1,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psiwm1(2,-1:1,2),psiwm2(2,-1:1,2)
      double complex psig(2,-1:1,2,2)
      double complex psigm1(2,-1:1,2,2),psigm2(2,-1:1,2,2),
     1               psim1g(2,-1:1,2,2),psim2g(2,-1:1,2,2)

c currents
      double complex jgqqbar(0:5,-1:1,2),jqqbarg(0:5,-1:1,2)
      double complex jgm1z(0:5,-1:1,2,2),jgm2z(0:5,-1:1,2,2),
     1               jm1zg(0:5,-1:1,2,2),jm2zg(0:5,-1:1,2,2),
     2               jm1gz(0:5,-1:1,2,2),jm2gz(0:5,-1:1,2,2)

c matrix variables
      double complex mat(2,6), madgr(6,-1:1)
      double complex mm(2)
      double complex m1kb(6), m2kb(6), temp, temp1
      double complex m1(3,2)


      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr
      double complex czero 
      external dotcc, s1c, dotrc, dotrr
      parameter (im = (0d0,1d0))
      parameter (czero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i

      logical ldebugs, ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebugs=.false.)

      integer N
      parameter( N = 1)

c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,1)    !als(1,1) : NLO alpha_s
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas       !alfas : LO alpha_s

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do l = 1,2
         do i = 1,6
            mat(l,i) = czero
         enddo
      enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
         pwp1wm1(mu) = q1p(mu,N) + qm1(mu,N)
         pwp1wm2(mu) = q1p(mu,N) + qm2(mu,N)
      enddo
      pwp1wm1(4) = dotrr(pwp1wm1(0),pwp1wm1(0))
      pwp1wm2(4) = dotrr(pwp1wm2(0),pwp1wm2(0))

c------------------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm1(0,N),wm1(1,N),psiwm1(1,-1,2),fqm1(0,2))     ! <q~W-(1)|
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm2(0,N),wm2(1,N),psiwm2(1,-1,2),fqm2(0,2))     ! <q~W-(2)|

      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm1(0,N),wm1(1,N),psiwm1(1,-1,1),fqm1(0,1))     ! |W-(1)q>
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm2(0,N),wm2(1,N),psiwm2(1,-1,1),fqm2(0,1))     ! |W-(2)q>

      do l =1,2                        !two gluon polarisations
         call polvec(qbar(0),l,eps(0,l))

         call bra2r(psi(1,-1,2),.true.,p(0,2),-1,q(0),eps(0,l),psig(1,-1,2,l),pq(0,2))        !<q~g|_l
         call ket2r(psi(1,-1,1),.true.,p(0,1),-1,q(0),eps(0,l),psig(1,-1,1,l),pq(0,1))        !|g q>_l

c get other bra and ket vectors

         call bra2c(psig(1,-1,2,l),.false.,pq(0,2),-1,qm1(0,N),wm1(1,N),psigm1(1,-1,2,l),pgm1(0,2))     ! <q~gW-(1)|
         call bra2c(psig(1,-1,2,l),.false.,pq(0,2),-1,qm2(0,N),wm2(1,N),psigm2(1,-1,2,l),pgm2(0,2))     ! <q~gW-(2)|

         call bra2r(psiwm1(1,-1,2),.false.,fqm1(0,2),-1,q(0),eps(0,l),psim1g(1,-1,2,l),pm1g(0,2))       ! <q~W-(1)g|
         call bra2r(psiwm2(1,-1,2),.false.,fqm2(0,2),-1,q(0),eps(0,l),psim2g(1,-1,2,l),pm2g(0,2))       ! <q~W-(2)g|

         call ket2c(psig(1,-1,1,l),.false.,pq(0,1),-1,qm1(0,N),wm1(1,N),psigm1(1,-1,1,l),pgm1(0,1))     ! |W-(1)gq>
         call ket2c(psig(1,-1,1,l),.false.,pq(0,1),-1,qm2(0,N),wm2(1,N),psigm2(1,-1,1,l),pgm2(0,1))     ! |W-(2)gq>

         call ket2r(psiwm1(1,-1,1),.false.,fqm1(0,1),-1,q(0),eps(0,l),psim1g(1,-1,1,l),pm1g(0,1))       ! |gW-(1)q>
         call ket2r(psiwm2(1,-1,1),.false.,fqm2(0,1),-1,q(0),eps(0,l),psim2g(1,-1,1,l),pm2g(0,1))       ! |gW-(2)q>

c get currents for contraction with W->WWW leptonic tensor

         call curr6(-1,psi(1,-1,2)   ,p(0,2) ,psig(1,-1,1,l),pq(0,1),jgqqbar(0,-1,l))      ! gz : <q~|gamma_mu|g q>_l
         call curr6(-1,psig(1,-1,2,l),pq(0,2),psi(1,-1,1)   ,p(0,1) ,jqqbarg(0,-1,l))      ! zg : <q~g|gamma_mu| q>_l

c get currents for contraction with A/Z->WW leptonic tensor

c currents with W-(1)
         call curr6(-1,psi(1,-1,2),p(0,2),psigm1(1,-1,1,l),pgm1(0,1),jgm1z(0,-1,1,l))      ! gwz: <q~|gamma_mu|W-(1)gq>
         call curr6(-1,psiwm1(1,-1,2),fqm1(0,2),psig(1,-1,1,l),pq(0,1),jgm1z(0,-1,2,l))    ! gzw: <q~W-(1)|gamma_mu|gq>

         call curr6(-1,psig(1,-1,2,l),pq(0,2),psiwm1(1,-1,1),fqm1(0,1),jm1zg(0,-1,1,l))    ! wzg: <q~g|gamma_mu|W-(1)q>
         call curr6(-1,psigm1(1,-1,2,l),pgm1(0,2),psi(1,-1,1),p(0,1),jm1zg(0,-1,2,l))      ! zwg: <q~gW-(1)|gamma_mu|q>

         call curr6(-1,psi(1,-1,2),p(0,2),psim1g(1,-1,1,l),pm1g(0,1),jm1gz(0,-1,1,l))      ! wgz: <q~|gamma_mu|gW-(1)q>
         call curr6(-1,psim1g(1,-1,2,l),pm1g(0,2),psi(1,-1,1),p(0,1),jm1gz(0,-1,2,l))      ! zgw: <q~W-(1)g|gamma_mu|q>

c currents with W-(2)
         call curr6(-1,psi(1,-1,2),p(0,2),psigm2(1,-1,1,l),pgm2(0,1),jgm2z(0,-1,1,l))      ! gwz: <q~|gamma_mu|W-(2)gq>
         call curr6(-1,psiwm2(1,-1,2),fqm2(0,2),psig(1,-1,1,l),pq(0,1),jgm2z(0,-1,2,l))    ! gzw: <q~W-(2)|gamma_mu|gq>

         call curr6(-1,psig(1,-1,2,l),pq(0,2),psiwm2(1,-1,1),fqm2(0,1),jm2zg(0,-1,1,l))    ! wzg: <q~g|gamma_mu|W-(2)q>
         call curr6(-1,psigm2(1,-1,2,l),pgm2(0,2),psi(1,-1,1),p(0,1),jm2zg(0,-1,2,l))      ! zwg: <q~gW-(2)|gamma_mu|q>

         call curr6(-1,psi(1,-1,2),p(0,2),psim2g(1,-1,1,l),pm2g(0,1),jm2gz(0,-1,1,l))      ! wgz: <q~|gamma_mu|gW-(2)q>
         call curr6(-1,psim2g(1,-1,2,l),pm2g(0,2),psi(1,-1,1),p(0,1),jm2gz(0,-1,2,l))      ! zgw: <q~W-(2)g|gamma_mu|q>

      enddo                            !gluon polarisation summation

C------------------------------------------------------------------------------
C  CALCULATION OF DU~ --> W-W+W-g --> LEPTONS+g PROCESSES (TREE LEVEL)

c      D U~ --> W-W+W-g --> LEPTONS+g
      do l=1,2
         m1(1,l) =   !W-(1) W+ W-(2) + g
     &      -s1c(psiwm2(1,-1,2),w1p(1,1),.true.,-1,psigm1(1,-1,1,l)) !g W-(1)W+(1)W-(2)
     &      -s1c(psiwm2(1,-1,2),w1p(1,1),.true.,-1,psim1g(1,-1,1,l)) !W-(1) g W+(1)W-(2)
     &      -s1c(psim2g(1,-1,2,l),w1p(1,1),.true.,-1,psiwm1(1,-1,1)) !W-(1)W+(1) g W-(2)
     &      -s1c(psigm2(1,-1,2,l),w1p(1,1),.true.,-1,psiwm1(1,-1,1)) !W-(1)W+(1)W-(2) g

         m1(2,l) =   !W-(2) W+ W-(1) + g
     &      -s1c(psiwm1(1,-1,2),w1p(1,1),.true.,-1,psigm2(1,-1,1,l)) !g W-(2)W+(1)W-(1)
     &      -s1c(psiwm1(1,-1,2),w1p(1,1),.true.,-1,psim2g(1,-1,1,l)) !W-(2) g W+(1)W-(1)
     &      -s1c(psim1g(1,-1,2,l),w1p(1,1),.true.,-1,psiwm2(1,-1,1)) !W-(2)W+(1) g W-(1)
     &      -s1c(psigm1(1,-1,2,l),w1p(1,1),.true.,-1,psiwm2(1,-1,1)) !W-(2)W+(1)W-(1) g

         mat(l,1) = m1(1,l)*clr(3,3,-1)**3
     &            + m1(2,l)*clr(3,3,-1)**3         ! !M(DU~->W-W+W-g ->Lg)
         m1(1,l) = czero
         m1(2,l) = czero
      enddo

C CALCULATION OF DU~ --> W- V g --> LEPTONS+g PROCESSES (TREE LEVEL)

c D U~ --> ZW-g --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_zww(0,1,3),jgm1z(0,-1,2,l))*clr(4,2,-1)     !g Z W-(1)
     &             -dotcc(WWW_zww(0,1,3),jm1gz(0,-1,2,l))*clr(4,2,-1)     !Z g W-(1)
     &             -dotcc(WWW_zww(0,1,3),jm1zg(0,-1,2,l))*clr(4,2,-1)     !Z W-(1) g

         m1(2,l) = -dotcc(WWW_zww(0,1,1),jgm2z(0,-1,2,l))*clr(4,2,-1)     !g Z W-(2)
     &             -dotcc(WWW_zww(0,1,1),jm2gz(0,-1,2,l))*clr(4,2,-1)     !Z g W-(2)
     &             -dotcc(WWW_zww(0,1,1),jm2zg(0,-1,2,l))*clr(4,2,-1)     !Z W-(2) g

c      additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         m1(3,l) = zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jgm1z(0,-1,2,l)) !g Z W-(1)
     &           + zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jm1gz(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jm1zg(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jgm2z(0,-1,2,l)) !g Z W-(2)
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jm2gz(0,-1,2,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jm2zg(0,-1,2,l))
         m1(3,l) = m1(3,l)*clr(4,2,-1)

         mat(l,2) = (m1(1,l) + m1(2,l) +m1(3,l))*CLR(3,3,-1)              !M(DU~->ZW-g ->Lg)
         m1(1,l) = czero
         m1(2,l) = czero
         m1(3,l) = czero
      enddo

c D U~ --> AW-g --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_aww(0,1,3),jgm1z(0,-1,2,l))*clr(4,1,-1)     !g A W-(1)
     &             -dotcc(WWW_aww(0,1,3),jm1gz(0,-1,2,l))*clr(4,1,-1)     !A g W-(1)
     &             -dotcc(WWW_aww(0,1,3),jm1zg(0,-1,2,l))*clr(4,1,-1)     !A W-(1) g

         m1(2,l) = -dotcc(WWW_aww(0,1,1),jgm2z(0,-1,2,l))*clr(4,1,-1)     !g A W-(2)
     &             -dotcc(WWW_aww(0,1,1),jm2gz(0,-1,2,l))*clr(4,1,-1)     !A g W-(2)
     &             -dotcc(WWW_aww(0,1,1),jm2zg(0,-1,2,l))*clr(4,1,-1)     !A W-(2) g

         mat(l,3) = (m1(1,l) + m1(2,l))*clr(3,3,-1)              !M(DU~->AW-g ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
      enddo

c D U~ --> W-Zg --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_zww(0,1,3),jgm1z(0,-1,1,l))*clr(3,2,-1)     !g W-(1) Z
     &             -dotcc(WWW_zww(0,1,3),jm1gz(0,-1,1,l))*clr(3,2,-1)     !W-(1) g Z
     &             -dotcc(WWW_zww(0,1,3),jm1zg(0,-1,1,l))*clr(3,2,-1)     !W-(1) Z g

         m1(2,l) = -dotcc(WWW_zww(0,1,1),jgm2z(0,-1,1,l))*clr(3,2,-1)     !g W-(2) Z
     &             -dotcc(WWW_zww(0,1,1),jm2gz(0,-1,1,l))*clr(3,2,-1)     !W-(2) g Z
     &             -dotcc(WWW_zww(0,1,1),jm2zg(0,-1,1,l))*clr(3,2,-1)     !W-(2) Z g

c      additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         m1(3,l) = zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jgm1z(0,-1,1,l)) !g W-(1) Z
     &           + zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jm1gz(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp1wm2(0),WWW_zww(0,1,3))*dotrc(pwp1wm2(0),jm1zg(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jgm2z(0,-1,1,l)) !g W-(2) Z
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jm2gz(0,-1,1,l))
     &           + zm2i(2)*dotrc(pwp1wm1(0),WWW_zww(0,1,1))*dotrc(pwp1wm1(0),jm2zg(0,-1,1,l))
         m1(3,l) = m1(3,l)*clr(3,2,-1)

         mat(l,4) = (m1(1,l) + m1(2,l) +m1(3,l))*clr(3,3,-1)              !M(DU~->W-Zg ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
         m1(3,l) = czero
      enddo

c D U~ --> W-Ag --> LEPTONS+g
      do l = 1,2
c      contract current with A/Z->WW tensors
         m1(1,l) = -dotcc(WWW_aww(0,1,3),jgm1z(0,-1,1,l))*clr(3,1,-1)     !g W-(1) A
     &             -dotcc(WWW_aww(0,1,3),jm1gz(0,-1,1,l))*clr(3,1,-1)     !W-(1) g A
     &             -dotcc(WWW_aww(0,1,3),jm1zg(0,-1,1,l))*clr(3,1,-1)     !W-(1) A g

         m1(2,l) = -dotcc(WWW_aww(0,1,1),jgm2z(0,-1,1,l))*clr(3,1,-1)     !g W-(2) A
     &             -dotcc(WWW_aww(0,1,1),jm2gz(0,-1,1,l))*clr(3,1,-1)     !W-(2) g A
     &             -dotcc(WWW_aww(0,1,1),jm2zg(0,-1,1,l))*clr(3,1,-1)     !W-(2) A g

         mat(l,5) = (m1(1,l) + m1(2,l))*clr(3,3,-1)              !M(DU~->W-Ag ->Lg)

         m1(1,l) = czero
         m1(2,l) = czero
      enddo

C CALCULATION OF DU~ --> W-W+W-g --> LEPTONS+g PROCESSES (TREE LEVEL)

c D U~ --> W-g --> LEPTONS+g
      do l = 1,2

         m1(1,l) = -dotcc(wmtowww(0,1),jgqqbar(0,-1,l))               ! GLUON EMISSION OFF D
     &             -dotcc(wmtowww(0,1),jqqbarg(0,-1,l))               ! GLUON EMISSION OFF U~

         mat(l,6) = - m1(1,l)*clr(3,3,-1)
         m1(1,l) = czero
      enddo

C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
c comparison with madgraph

#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do i = 1,6
            m1kb(i) = czero ! kartesian basis l = 1
            m2kb(i) = czero ! kartesian basis l = 2
            do lh = -1,1,2
               madgr(i,lh) = czero
            enddo
         enddo

         do lh = -1,1,2
            do i =1,6
c               madgr(i,lh) = WWWtampr(i,-1,lh,1)       ! for udbar process
               madgr(i,lh) = WWWtampr(i,-1,lh,2)       ! for dubar process
            enddo
         enddo

         do i=1,6
            m1kb(i) =     (madgr(i,-1) - madgr(i,1))/sqrt(2d0)
            m2kb(i) =-im*(madgr(i,-1) + madgr(i,1))/sqrt(2d0)
         enddo

c print results from checks:

         do i=1,6
            print*,' '
            print*,'ratio mat/m1kb (for l=1),for i = ',i,':',
     &         '(', dreal(mat(1,i))/dreal(m1kb(i)) ,
     &              dimag(mat(1,i))/dimag(m1kb(i)) , ')'
            print*,'mat(1,',i,') = ',mat(1,i)
            temp = dreal(mat(1,i))**2+dimag(mat(1,i))**2
            print*,'temp = ',temp
            print*,'m1kb(',i,') =  ',m1kb(i)
            temp1 = dreal(m1kb(i))**2+dimag(m1kb(i))**2
            print*,'temp1 = ',temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '

            print*,'ratio mat/m2kb (for l=1),for i = ',i,':',
     &         '(', dreal(mat(2,i))/dreal(m2kb(i)) ,
     &              dimag(mat(2,i))/dimag(m2kb(i)), ')'
            print*,'mat(2,',i,') = ',mat(2,i)
            temp = dreal(mat(2,i))**2+dimag(mat(2,i))**2
            print*,'temp = ',temp
            print*,'m2kb(',i,') =  ',m2kb(i)
            temp1 = dreal(m2kb(i))**2+dimag(m2kb(i))**2
            print*,'temp1 = ',temp1
            print*,'temp/temp1 = ',temp/temp1
            print*,' '
         enddo

      endif ! end of debugging block
#endif

c------------------------------------------------------------------------------

c sum the graphs, square them and map them onto qqbarwwwj

      res = 0.0d0
      do l = 1,2
         mm(l) = 0
         do i = 1,6
            mm(l) = mm(l) +  mat(l,i)
         enddo
         res = res + dreal(mm(l))**2 + dimag(mm(l))**2
      enddo
      res = res*4d0*fpials(2)  ! C_2*3 is the color factor

      qqbarWWWj = res

      return
      end

C------------------------------------------------------------------------------

c*******************************************************************************
c  end of subroutine qqWWWpj
c*******************************************************************************

