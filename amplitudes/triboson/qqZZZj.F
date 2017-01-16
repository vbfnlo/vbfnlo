c****************subroutine QQZZZJ *************************************
c
c	Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version:  2007, September
c	Last modified: 
C
C  QQZZZ calculates the matrix elements**2 for 
C        q1 q2    ---->   Z Z Z g,   
C         Z ---> f1-bar f2, Z ---> f3-bar f4, Z ---> f5-bar f6
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
c      return uuzzz = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuzzzj,2:ddzzzj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqZZZj(pbar,sign, qbar, gsign, uuzzzj, ddzzzj, NW)

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

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  uuzzzj, ddzzzj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  ifl(2), lh, l

      double precision fqze(0:4,2), fqzmu(0:4,2)
      double precision  fqztau(0:4,2)
      double precision pq(0:4,2), pgZtau(0:4,2), pZtaug(0:4,2)
      double precision pgZmu(0:4,2), pZmug(0:4,2)
      double precision pgZe(0:4,2), pZeg(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psize(2,-1:1,2), psiae(2,-1:1,2),
     1               psizmu(2,-1:1,2), psiamu(2,-1:1,2),
     2               psiztau(2,-1:1,2), psiatau(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgZtau(2,-1:1,2,2), braketgatau(2,-1:1,2,2)
      double complex braketgZmu(2,-1:1,2,2), braketgamu(2,-1:1,2,2)
      double complex braketgZe(2,-1:1,2,2), braketgAe(2,-1:1,2,2)
      double complex psiZtaug(2,-1:1, 2,2), psiAtaug(2,-1:1, 2,2)
      double complex psiZmug(2,-1:1,2,2), psiAmug(2,-1:1,2,2)
      double complex psiZeg(2,-1:1,2,2), psiAeg(2,-1:1,2,2)

      double complex jae(0:5,-1:1,4,2), jam(0:5,-1:1,4,2), jat(0:5,-1:1,4,2),
     1               jze(0:5,-1:1,4,2), jzm(0:5,-1:1,4,2), jzt(0:5,-1:1,4,2)
      double complex mat(2,-1:1,2,10), madgr(10,-1:1,-1:1,2)
      double complex m1kb(2,-1:1,10), m2kb(2,-1:1,10), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex ma, mz, m1, m2, m3, m4
      double complex psi(2,-1:1,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2)
      double complex jqqa(0:5,-1:1,2)

      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      save ifl
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c  helicity selection

      INTEGER h
      COMMON /hcount / h

      logical ldebug, ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false.)
      integer NW !NW = 1 means only the qqbar -> ZZZ diagrams
                 !NW = 0 means all

c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzzzj: u = 3
         ifl(k) = 3
         k = 2                  ! ddzzzj: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,10
               do l = 1,2
                  mat(k,isig,l,i) = (0.0d0,0.0d0)
               enddo
            enddo
         enddo
      enddo

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

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,1),
     &              ze(1,1),psize(1,isig,2),fqze(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qzmu(0,1),
     &              zmu(1,1),psizmu(1,isig,2),fqzmu(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qztau(0,1),
     &              ztau(1,1),psiztau(1,isig,2),fqztau(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,1),
     &              ae(1,1),psiae(1,isig,2),fqze(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qzmu(0,1),
     &              amu(1,1),psiamu(1,isig,2),fqzmu(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qztau(0,1),
     &              atau(1,1),psiatau(1,isig,2),fqztau(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,1),
     &              ze(1,1), psize(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qzmu(0,1),
     &              zmu(1,1), psizmu(1,isig,1),fqzmu(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qztau(0,1),
     &              ztau(1,1), psiztau(1,isig,1),fqztau(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,1),
     &              ae(1,1), psiae(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qzmu(0,1),
     &              amu(1,1), psiamu(1,isig,1),fqzmu(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qztau(0,1),
     &              atau(1,1), psiatau(1,isig,1),fqztau(0,1))

      enddo

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2        ! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck
c            eps(mu,l) = qbar(mu)
c         enddo

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
     &         qztau(0,1),atau(1,1),braketgAtau(1,isig,1,l),pgZtau(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &         qztau(0,1),atau(1,1),braketgAtau(1,isig,2,l),pgZtau(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &         qzmu(0,1),amu(1,1),braketgAmu(1,isig,1,l),pgZmu(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &         qzmu(0,1),amu(1,1),braketgAmu(1,isig,2,l),pgZmu(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &         qze(0,1),ae(1,1),braketgAe(1,isig,1,l),pgZe(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &         qze(0,1),ae(1,1),braketgAe(1,isig,2,l),pgZe(0,2))

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qztau(0,1),ztau(1,1),braketgZtau(1,isig,1,l),pgZtau(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qztau(0,1),ztau(1,1),braketgZtau(1,isig,2,l),pgZtau(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qzmu(0,1),zmu(1,1),braketgZmu(1,isig,1,l),pgZmu(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qzmu(0,1),zmu(1,1),braketgZmu(1,isig,2,l),pgZmu(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qze(0,1),ze(1,1),braketgZe(1,isig,1,l),pgZe(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qze(0,1),ze(1,1),braketgZe(1,isig,2,l),pgZe(0,2))

            call ket2r(psiAtau(1,isig,1),.false.,fqztau(0,1),isig,q,eps(0,l),
     &           psiAtaug(1,isig,1,l),pZtaug(0,1))
            call bra2r(psiAtau(1,isig,2),.false.,fqztau(0,2),isig,q,eps(0,l),
     &           psiAtaug(1,isig,2,l),pZtaug(0,2))
            call ket2r(psiAmu(1,isig,1),.false.,fqzmu(0,1),isig,q,eps(0,l),
     &           psiAmug(1,isig,1,l),pZmug(0,1))
            call bra2r(psiAmu(1,isig,2),.false.,fqzmu(0,2),isig,q,eps(0,l),
     &           psiAmug(1,isig,2,l),pZmug(0,2))
            call ket2r(psiAe(1,isig,1),.false.,fqze(0,1),isig,q,eps(0,l),
     &           psiAeg(1,isig,1,l),pZeg(0,1))
            call bra2r(psiAe(1,isig,2),.false.,fqze(0,2),isig,q,eps(0,l),
     &           psiAeg(1,isig,2,l),pZeg(0,2))

            call ket2r(psiZtau(1,isig,1),.false.,fqztau(0,1),isig,q,eps(0,l),
     &           psiZtaug(1,isig,1,l),pZtaug(0,1))
            call bra2r(psiZtau(1,isig,2),.false.,fqztau(0,2),isig,q,eps(0,l),
     &           psiZtaug(1,isig,2,l),pZtaug(0,2))
            call ket2r(psiZmu(1,isig,1),.false.,fqzmu(0,1),isig,q,eps(0,l),
     &           psiZmug(1,isig,1,l),pZmug(0,1))
            call bra2r(psiZmu(1,isig,2),.false.,fqzmu(0,2),isig,q,eps(0,l),
     &           psiZmug(1,isig,2,l),pZmug(0,2))
            call ket2r(psiZe(1,isig,1),.false.,fqze(0,1),isig,q,eps(0,l),
     &           psiZeg(1,isig,1,l),pZeg(0,1))
            call bra2r(psiZe(1,isig,2),.false.,fqze(0,2),isig,q,eps(0,l),
     &           psiZeg(1,isig,2,l),pZeg(0,2))

         enddo                   !fermion helicity
      enddo                      !gluon polarization


c Calculation of processes qq -> Z Z Z g ---------------------

c   u u -> Z Z Z g

      do l=1,2
         do isig = -1,1,2
           m1=-s1c(psizeg(1,isig,2,l),ztau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(braketgze(1,isig,2,l),ztau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(psize(1,isig,2),ztau(1,1),.true.,isig,psizmug(1,isig,1,l))
     1      -s1c(psize(1,isig,2),ztau(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     2      -s1c(psizeg(1,isig,2,l),zmu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(braketgze(1,isig,2,l),zmu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(psize(1,isig,2),zmu(1,1),.true.,isig,psiztaug(1,isig,1,l))
     2      -s1c(psize(1,isig,2),zmu(1,1),.true.,isig,braketgztau(1,isig,1,l))
     3      -s1c(psizmug(1,isig,2,l),ztau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(braketgzmu(1,isig,2,l),ztau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(psizmu(1,isig,2),ztau(1,1),.true.,isig,psizeg(1,isig,1,l))
     3      -s1c(psizmu(1,isig,2),ztau(1,1),.true.,isig,braketgze(1,isig,1,l))
     4      -s1c(psizmug(1,isig,2,l),ze(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(braketgzmu(1,isig,2,l),ze(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(psizmu(1,isig,2),ze(1,1),.true.,isig,psiztaug(1,isig,1,l))
     4      -s1c(psizmu(1,isig,2),ze(1,1),.true.,isig,braketgztau(1,isig,1,l))
     5      -s1c(psiztaug(1,isig,2,l),zmu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(braketgztau(1,isig,2,l),zmu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(psiztau(1,isig,2),zmu(1,1),.true.,isig,psizeg(1,isig,1,l))
     5      -s1c(psiztau(1,isig,2),zmu(1,1),.true.,isig,braketgze(1,isig,1,l))
     6      -s1c(psiztaug(1,isig,2,l),ze(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(braketgztau(1,isig,2,l),ze(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(psiztau(1,isig,2),ze(1,1),.true.,isig,psizmug(1,isig,1,l))
     6      -s1c(psiztau(1,isig,2),ze(1,1),.true.,isig,braketgzmu(1,isig,1,l))

           if (NW.eq.0) then
           m2=-s1c(psiaeg(1,isig,2,l),atau(1,1),.true.,isig,psiamu(1,isig,1))
     1      -s1c(braketgae(1,isig,2,l),atau(1,1),.true.,isig,psiamu(1,isig,1))
     1      -s1c(psiae(1,isig,2),atau(1,1),.true.,isig,psiamug(1,isig,1,l))
     1      -s1c(psiae(1,isig,2),atau(1,1),.true.,isig,braketgamu(1,isig,1,l))
     2      -s1c(psiaeg(1,isig,2,l),amu(1,1),.true.,isig,psiatau(1,isig,1))
     2      -s1c(braketgae(1,isig,2,l),amu(1,1),.true.,isig,psiatau(1,isig,1))
     2      -s1c(psiae(1,isig,2),amu(1,1),.true.,isig,psiataug(1,isig,1,l))
     2      -s1c(psiae(1,isig,2),amu(1,1),.true.,isig,braketgatau(1,isig,1,l))
     3      -s1c(psiamug(1,isig,2,l),atau(1,1),.true.,isig,psiae(1,isig,1))
     3      -s1c(braketgamu(1,isig,2,l),atau(1,1),.true.,isig,psiae(1,isig,1))
     3      -s1c(psiamu(1,isig,2),atau(1,1),.true.,isig,psiaeg(1,isig,1,l))
     3      -s1c(psiamu(1,isig,2),atau(1,1),.true.,isig,braketgae(1,isig,1,l))
     4      -s1c(psiamug(1,isig,2,l),ae(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(braketgamu(1,isig,2,l),ae(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(psiamu(1,isig,2),ae(1,1),.true.,isig,psiataug(1,isig,1,l))
     4      -s1c(psiamu(1,isig,2),ae(1,1),.true.,isig,braketgatau(1,isig,1,l))
     5      -s1c(psiataug(1,isig,2,l),amu(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(braketgatau(1,isig,2,l),amu(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(psiatau(1,isig,2),amu(1,1),.true.,isig,psiaeg(1,isig,1,l))
     5      -s1c(psiatau(1,isig,2),amu(1,1),.true.,isig,braketgae(1,isig,1,l))
     6      -s1c(psiataug(1,isig,2,l),ae(1,1),.true.,isig,psiamu(1,isig,1))
     6      -s1c(braketgatau(1,isig,2,l),ae(1,1),.true.,isig,psiamu(1,isig,1))
     6      -s1c(psiatau(1,isig,2),ae(1,1),.true.,isig,psiamug(1,isig,1,l))
     6      -s1c(psiatau(1,isig,2),ae(1,1),.true.,isig,braketgamu(1,isig,1,l))

           m3=-s1c(psiaeg(1,isig,2,l),ztau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(braketgae(1,isig,2,l),ztau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(psiae(1,isig,2),ztau(1,1),.true.,isig,psizmug(1,isig,1,l))
     1      -s1c(psiae(1,isig,2),ztau(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     2      -s1c(psiaeg(1,isig,2,l),zmu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(braketgae(1,isig,2,l),zmu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(psiae(1,isig,2),zmu(1,1),.true.,isig,psiztaug(1,isig,1,l))
     2      -s1c(psiae(1,isig,2),zmu(1,1),.true.,isig,braketgztau(1,isig,1,l))
     3      -s1c(psiamug(1,isig,2,l),ztau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(braketgamu(1,isig,2,l),ztau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(psiamu(1,isig,2),ztau(1,1),.true.,isig,psizeg(1,isig,1,l))
     3      -s1c(psiamu(1,isig,2),ztau(1,1),.true.,isig,braketgze(1,isig,1,l))
     4      -s1c(psiamug(1,isig,2,l),ze(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(braketgamu(1,isig,2,l),ze(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(psiamu(1,isig,2),ze(1,1),.true.,isig,psiztaug(1,isig,1,l))
     4      -s1c(psiamu(1,isig,2),ze(1,1),.true.,isig,braketgztau(1,isig,1,l))
     5      -s1c(psiataug(1,isig,2,l),zmu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(braketgatau(1,isig,2,l),zmu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(psiatau(1,isig,2),zmu(1,1),.true.,isig,psizeg(1,isig,1,l))
     5      -s1c(psiatau(1,isig,2),zmu(1,1),.true.,isig,braketgze(1,isig,1,l))
     6      -s1c(psiataug(1,isig,2,l),ze(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(braketgatau(1,isig,2,l),ze(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(psiatau(1,isig,2),ze(1,1),.true.,isig,psizmug(1,isig,1,l))
     6      -s1c(psiatau(1,isig,2),ze(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     7      -s1c(psizeg(1,isig,2,l),atau(1,1),.true.,isig,psizmu(1,isig,1))
     7      -s1c(braketgze(1,isig,2,l),atau(1,1),.true.,isig,psizmu(1,isig,1))
     7      -s1c(psize(1,isig,2),atau(1,1),.true.,isig,psizmug(1,isig,1,l))
     7      -s1c(psize(1,isig,2),atau(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     8      -s1c(psizeg(1,isig,2,l),amu(1,1),.true.,isig,psiztau(1,isig,1))
     8      -s1c(braketgze(1,isig,2,l),amu(1,1),.true.,isig,psiztau(1,isig,1))
     8      -s1c(psize(1,isig,2),amu(1,1),.true.,isig,psiztaug(1,isig,1,l))
     8      -s1c(psize(1,isig,2),amu(1,1),.true.,isig,braketgztau(1,isig,1,l))
     9      -s1c(psizmug(1,isig,2,l),atau(1,1),.true.,isig,psize(1,isig,1))
     9      -s1c(braketgzmu(1,isig,2,l),atau(1,1),.true.,isig,psize(1,isig,1))
     9      -s1c(psizmu(1,isig,2),atau(1,1),.true.,isig,psizeg(1,isig,1,l))
     9      -s1c(psizmu(1,isig,2),atau(1,1),.true.,isig,braketgze(1,isig,1,l))
     &      -s1c(psizmug(1,isig,2,l),ae(1,1),.true.,isig,psiztau(1,isig,1))
     &      -s1c(braketgzmu(1,isig,2,l),ae(1,1),.true.,isig,psiztau(1,isig,1))
     &      -s1c(psizmu(1,isig,2),ae(1,1),.true.,isig,psiztaug(1,isig,1,l))
     &      -s1c(psizmu(1,isig,2),ae(1,1),.true.,isig,braketgztau(1,isig,1,l))
     1      -s1c(psiztaug(1,isig,2,l),amu(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(braketgztau(1,isig,2,l),amu(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(psiztau(1,isig,2),amu(1,1),.true.,isig,psizeg(1,isig,1,l))
     1      -s1c(psiztau(1,isig,2),amu(1,1),.true.,isig,braketgze(1,isig,1,l))
     2      -s1c(psiztaug(1,isig,2,l),ae(1,1),.true.,isig,psizmu(1,isig,1))
     2      -s1c(braketgztau(1,isig,2,l),ae(1,1),.true.,isig,psizmu(1,isig,1))
     2      -s1c(psiztau(1,isig,2),ae(1,1),.true.,isig,psizmug(1,isig,1,l))
     2      -s1c(psiztau(1,isig,2),ae(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     3      -s1c(psizeg(1,isig,2,l),ztau(1,1),.true.,isig,psiamu(1,isig,1))
     3      -s1c(braketgze(1,isig,2,l),ztau(1,1),.true.,isig,psiamu(1,isig,1))
     3      -s1c(psize(1,isig,2),ztau(1,1),.true.,isig,psiamug(1,isig,1,l))
     3      -s1c(psize(1,isig,2),ztau(1,1),.true.,isig,braketgamu(1,isig,1,l))
     4      -s1c(psizeg(1,isig,2,l),zmu(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(braketgze(1,isig,2,l),zmu(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(psize(1,isig,2),zmu(1,1),.true.,isig,psiataug(1,isig,1,l))
     4      -s1c(psize(1,isig,2),zmu(1,1),.true.,isig,braketgatau(1,isig,1,l))
     5      -s1c(psizmug(1,isig,2,l),ztau(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(braketgzmu(1,isig,2,l),ztau(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(psizmu(1,isig,2),ztau(1,1),.true.,isig,psiaeg(1,isig,1,l))
     5      -s1c(psizmu(1,isig,2),ztau(1,1),.true.,isig,braketgae(1,isig,1,l))
     6      -s1c(psizmug(1,isig,2,l),ze(1,1),.true.,isig,psiatau(1,isig,1))
     6      -s1c(braketgzmu(1,isig,2,l),ze(1,1),.true.,isig,psiatau(1,isig,1))
     6      -s1c(psizmu(1,isig,2),ze(1,1),.true.,isig,psiataug(1,isig,1,l))
     6      -s1c(psizmu(1,isig,2),ze(1,1),.true.,isig,braketgatau(1,isig,1,l))
     7      -s1c(psiztaug(1,isig,2,l),zmu(1,1),.true.,isig,psiae(1,isig,1))
     7      -s1c(braketgztau(1,isig,2,l),zmu(1,1),.true.,isig,psiae(1,isig,1))
     7      -s1c(psiztau(1,isig,2),zmu(1,1),.true.,isig,psiaeg(1,isig,1,l))
     7      -s1c(psiztau(1,isig,2),zmu(1,1),.true.,isig,braketgae(1,isig,1,l))
     8      -s1c(psiztaug(1,isig,2,l),ze(1,1),.true.,isig,psiamu(1,isig,1))
     8      -s1c(braketgztau(1,isig,2,l),ze(1,1),.true.,isig,psiamu(1,isig,1))
     8      -s1c(psiztau(1,isig,2),ze(1,1),.true.,isig,psiamug(1,isig,1,l))
     8      -s1c(psiztau(1,isig,2),ze(1,1),.true.,isig,braketgamu(1,isig,1,l))

           m4=-s1c(psiaeg(1,isig,2,l),atau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(braketgae(1,isig,2,l),atau(1,1),.true.,isig,psizmu(1,isig,1))
     1      -s1c(psiae(1,isig,2),atau(1,1),.true.,isig,psizmug(1,isig,1,l))
     1      -s1c(psiae(1,isig,2),atau(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     2      -s1c(psiaeg(1,isig,2,l),amu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(braketgae(1,isig,2,l),amu(1,1),.true.,isig,psiztau(1,isig,1))
     2      -s1c(psiae(1,isig,2),amu(1,1),.true.,isig,psiztaug(1,isig,1,l))
     2      -s1c(psiae(1,isig,2),amu(1,1),.true.,isig,braketgztau(1,isig,1,l))
     3      -s1c(psiamug(1,isig,2,l),atau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(braketgamu(1,isig,2,l),atau(1,1),.true.,isig,psize(1,isig,1))
     3      -s1c(psiamu(1,isig,2),atau(1,1),.true.,isig,psizeg(1,isig,1,l))
     3      -s1c(psiamu(1,isig,2),atau(1,1),.true.,isig,braketgze(1,isig,1,l))
     4      -s1c(psiamug(1,isig,2,l),ae(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(braketgamu(1,isig,2,l),ae(1,1),.true.,isig,psiztau(1,isig,1))
     4      -s1c(psiamu(1,isig,2),ae(1,1),.true.,isig,psiztaug(1,isig,1,l))
     4      -s1c(psiamu(1,isig,2),ae(1,1),.true.,isig,braketgztau(1,isig,1,l))
     5      -s1c(psiataug(1,isig,2,l),amu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(braketgatau(1,isig,2,l),amu(1,1),.true.,isig,psize(1,isig,1))
     5      -s1c(psiatau(1,isig,2),amu(1,1),.true.,isig,psizeg(1,isig,1,l))
     5      -s1c(psiatau(1,isig,2),amu(1,1),.true.,isig,braketgze(1,isig,1,l))
     6      -s1c(psiataug(1,isig,2,l),ae(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(braketgatau(1,isig,2,l),ae(1,1),.true.,isig,psizmu(1,isig,1))
     6      -s1c(psiatau(1,isig,2),ae(1,1),.true.,isig,psizmug(1,isig,1,l))
     6      -s1c(psiatau(1,isig,2),ae(1,1),.true.,isig,braketgzmu(1,isig,1,l))
     7      -s1c(psizeg(1,isig,2,l),atau(1,1),.true.,isig,psiamu(1,isig,1))
     7      -s1c(braketgze(1,isig,2,l),atau(1,1),.true.,isig,psiamu(1,isig,1))
     7      -s1c(psize(1,isig,2),atau(1,1),.true.,isig,psiamug(1,isig,1,l))
     7      -s1c(psize(1,isig,2),atau(1,1),.true.,isig,braketgamu(1,isig,1,l))
     8      -s1c(psizeg(1,isig,2,l),amu(1,1),.true.,isig,psiatau(1,isig,1))
     8      -s1c(braketgze(1,isig,2,l),amu(1,1),.true.,isig,psiatau(1,isig,1))
     8      -s1c(psize(1,isig,2),amu(1,1),.true.,isig,psiataug(1,isig,1,l))
     8      -s1c(psize(1,isig,2),amu(1,1),.true.,isig,braketgatau(1,isig,1,l))
     9      -s1c(psizmug(1,isig,2,l),atau(1,1),.true.,isig,psiae(1,isig,1))
     9      -s1c(braketgzmu(1,isig,2,l),atau(1,1),.true.,isig,psiae(1,isig,1))
     9      -s1c(psizmu(1,isig,2),atau(1,1),.true.,isig,psiaeg(1,isig,1,l))
     9      -s1c(psizmu(1,isig,2),atau(1,1),.true.,isig,braketgae(1,isig,1,l))
     &      -s1c(psizmug(1,isig,2,l),ae(1,1),.true.,isig,psiatau(1,isig,1))
     &      -s1c(braketgzmu(1,isig,2,l),ae(1,1),.true.,isig,psiatau(1,isig,1))
     &      -s1c(psizmu(1,isig,2),ae(1,1),.true.,isig,psiataug(1,isig,1,l))
     &      -s1c(psizmu(1,isig,2),ae(1,1),.true.,isig,braketgatau(1,isig,1,l))
     1      -s1c(psiztaug(1,isig,2,l),amu(1,1),.true.,isig,psiae(1,isig,1))
     1      -s1c(braketgztau(1,isig,2,l),amu(1,1),.true.,isig,psiae(1,isig,1))
     1      -s1c(psiztau(1,isig,2),amu(1,1),.true.,isig,psiaeg(1,isig,1,l))
     1      -s1c(psiztau(1,isig,2),amu(1,1),.true.,isig,braketgae(1,isig,1,l))
     2      -s1c(psiztaug(1,isig,2,l),ae(1,1),.true.,isig,psiamu(1,isig,1))
     2      -s1c(braketgztau(1,isig,2,l),ae(1,1),.true.,isig,psiamu(1,isig,1))
     2      -s1c(psiztau(1,isig,2),ae(1,1),.true.,isig,psiamug(1,isig,1,l))
     2      -s1c(psiztau(1,isig,2),ae(1,1),.true.,isig,braketgamu(1,isig,1,l))
     3      -s1c(psiaeg(1,isig,2,l),ztau(1,1),.true.,isig,psiamu(1,isig,1))
     3      -s1c(braketgae(1,isig,2,l),ztau(1,1),.true.,isig,psiamu(1,isig,1))
     3      -s1c(psiae(1,isig,2),ztau(1,1),.true.,isig,psiamug(1,isig,1,l))
     3      -s1c(psiae(1,isig,2),ztau(1,1),.true.,isig,braketgamu(1,isig,1,l))
     4      -s1c(psiaeg(1,isig,2,l),zmu(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(braketgae(1,isig,2,l),zmu(1,1),.true.,isig,psiatau(1,isig,1))
     4      -s1c(psiae(1,isig,2),zmu(1,1),.true.,isig,psiataug(1,isig,1,l))
     4      -s1c(psiae(1,isig,2),zmu(1,1),.true.,isig,braketgatau(1,isig,1,l))
     5      -s1c(psiamug(1,isig,2,l),ztau(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(braketgamu(1,isig,2,l),ztau(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(psiamu(1,isig,2),ztau(1,1),.true.,isig,psiaeg(1,isig,1,l))
     5      -s1c(psiamu(1,isig,2),ztau(1,1),.true.,isig,braketgae(1,isig,1,l))
     6      -s1c(psiamug(1,isig,2,l),ze(1,1),.true.,isig,psiatau(1,isig,1))
     6      -s1c(braketgamu(1,isig,2,l),ze(1,1),.true.,isig,psiatau(1,isig,1))
     6      -s1c(psiamu(1,isig,2),ze(1,1),.true.,isig,psiataug(1,isig,1,l))
     6      -s1c(psiamu(1,isig,2),ze(1,1),.true.,isig,braketgatau(1,isig,1,l))
     7      -s1c(psiataug(1,isig,2,l),zmu(1,1),.true.,isig,psiae(1,isig,1))
     7      -s1c(braketgatau(1,isig,2,l),zmu(1,1),.true.,isig,psiae(1,isig,1))
     7      -s1c(psiatau(1,isig,2),zmu(1,1),.true.,isig,psiaeg(1,isig,1,l))
     7      -s1c(psiatau(1,isig,2),zmu(1,1),.true.,isig,braketgae(1,isig,1,l))
     8      -s1c(psiataug(1,isig,2,l),ze(1,1),.true.,isig,psiamu(1,isig,1))
     8      -s1c(braketgatau(1,isig,2,l),ze(1,1),.true.,isig,psiamu(1,isig,1))
     8      -s1c(psiatau(1,isig,2),ze(1,1),.true.,isig,psiamug(1,isig,1,l))
     8      -s1c(psiatau(1,isig,2),ze(1,1),.true.,isig,braketgamu(1,isig,1,l))
           endif !NW = 0

           if (NW.eq.0) then
            mat(1,isig,l,1) =  m1 * clr(3,2,isig)**3 ! uu -> ZZZ
     1                  + m2 * clr(3,1,isig)**3 ! uu -> AAA
     2                  + m3 * clr(3,2,isig)**2 * clr(3,1,isig) ! uu -> ZZA
     3                  + m4 * clr(3,2,isig) * clr(3,1,isig)**2 ! uu -> AAZ
            mat(2,isig,l,1) =  m1 * clr(4,2,isig)**3 ! dd -> ZZZ
     1                  + m2 * clr(4,1,isig)**3 ! dd -> AAA
     2                  + m3 * clr(4,2,isig)**2 * clr(4,1,isig) ! dd -> ZZA
     3                  + m4 * clr(4,2,isig) * clr(4,1,isig)**2 ! dd -> AAZ
           else if (NW.eq.1) then
            mat(1,isig,l,1) =  m1 * clr(3,2,isig)**3 ! uu -> ZZZ
            mat(2,isig,l,1) =  m1 * clr(4,2,isig)**3 ! dd -> ZZZ
           endif
         enddo
      enddo

c Calculation of processes qq -> V V g -> Z Z Z g -------------

      if (NW.eq.0) then
c qq -> A A/Z (-> Z Z) g

      do l = 1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgae(1,-1,1,l),pgZe(0,1), jae(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psiaeg(1,-1,1,l),pZeg(0,1), jae(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psiae(1,-1,1),fqze(0,1), jae(0,-1,3,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgamu(1,-1,1,l),pgZmu(0,1), jam(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psiamug(1,-1,1,l),pZmug(0,1), jam(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psiamu(1,-1,1),fqzmu(0,1), jam(0,-1,3,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgatau(1,-1,1,l),pgZtau(0,1), jat(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psiataug(1,-1,1,l),pZtaug(0,1), jat(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psiatau(1,-1,1),fqztau(0,1), jat(0,-1,3,l) )

         do isig = -1,1,2

            do mu = 0,5
               jae(mu,isig,4,l) = jae(mu,isig,1,l) + jae(mu,isig,2,l)
     &                          + jae(mu,isig,3,l)
               jam(mu,isig,4,l) = jam(mu,isig,1,l) + jam(mu,isig,2,l)
     &                          + jam(mu,isig,3,l)
               jat(mu,isig,4,l) = jat(mu,isig,1,l) + jat(mu,isig,2,l)
     &                          + jat(mu,isig,3,l)
            enddo


c contract with atozz/ztozz tensor

            m1 = -dotcc(azmzt(0,1),jae(0,isig,4,l))
     1           -dotcc(azezt(0,1),jam(0,isig,4,l))
     2           -dotcc(azezm(0,1),jat(0,isig,4,l))
            m2 = -dotcc(zzmzt(0,1),jae(0,isig,4,l))
     1           -dotcc(zzezt(0,1),jam(0,isig,4,l))
     2           -dotcc(zzezm(0,1),jat(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,2) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
            enddo

         enddo

      enddo

c qq -> A/Z(-> Z Z) A g

      do l = 1,2

         call curr6(1,psiae(1,-1,2),fqze(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jae(0,-1,1,l) )
         call curr6(1,psiaeg(1,-1,2,l),pZeg(0,2),
     &           psi(1,-1,1),p(0,1), jae(0,-1,2,l) )
         call curr6(1,braketgae(1,-1,2,l),pgZe(0,2),
     &           psi(1,-1,1),p(0,1), jae(0,-1,3,l) )

         call curr6(1,psiamu(1,-1,2),fqzmu(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jam(0,-1,1,l) )
         call curr6(1,psiamug(1,-1,2,l),pZmug(0,2),
     &           psi(1,-1,1),p(0,1), jam(0,-1,2,l) )
         call curr6(1,braketgamu(1,-1,2,l),pgZmu(0,2),
     &           psi(1,-1,1),p(0,1), jam(0,-1,3,l) )

         call curr6(1,psiatau(1,-1,2),fqztau(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jat(0,-1,1,l) )
         call curr6(1,psiataug(1,-1,2,l),pZtaug(0,2),
     &           psi(1,-1,1),p(0,1), jat(0,-1,2,l) )
         call curr6(1,braketgatau(1,-1,2,l),pgZtau(0,2),
     &           psi(1,-1,1),p(0,1), jat(0,-1,3,l) )


         do isig = -1,1,2

            do mu = 0,5
               jae(mu,isig,4,l) = jae(mu,isig,1,l) + jae(mu,isig,2,l)
     &                          + jae(mu,isig,3,l)
               jam(mu,isig,4,l) = jam(mu,isig,1,l) + jam(mu,isig,2,l)
     &                          + jam(mu,isig,3,l)
               jat(mu,isig,4,l) = jat(mu,isig,1,l) + jat(mu,isig,2,l)
     &                          + jat(mu,isig,3,l)
            enddo

c contract with atozz/ztozz tensor

            m1 = -dotcc(azmzt(0,1),jae(0,isig,4,l))
     1           -dotcc(azezt(0,1),jam(0,isig,4,l))
     2           -dotcc(azezm(0,1),jat(0,isig,4,l))
            m2 = -dotcc(zzmzt(0,1),jae(0,isig,4,l))
     1           -dotcc(zzezt(0,1),jam(0,isig,4,l))
     2           -dotcc(zzezm(0,1),jat(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,3) = (m1*clr(ifl(k),1,isig)
     1                     +  m2*clr(ifl(k),2,isig))
     2                       * clr(ifl(k),1,isig)
            enddo

         enddo
      enddo

c qq -> Z A/Z (-> Z Z) g

      do l = 1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgze(1,-1,1,l),pgze(0,1), jze(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psizeg(1,-1,1,l),pzeg(0,1), jze(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psize(1,-1,1),fqze(0,1), jze(0,-1,3,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgzmu(1,-1,1,l),pgzmu(0,1), jzm(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psizmug(1,-1,1,l),pzmug(0,1), jzm(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psizmu(1,-1,1),fqzmu(0,1), jzm(0,-1,3,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketgztau(1,-1,1,l),pgztau(0,1), jzt(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psiztaug(1,-1,1,l),pztaug(0,1), jzt(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psiztau(1,-1,1),fqztau(0,1), jzt(0,-1,3,l) )

         do isig = -1,1,2

            do mu = 0,5
               jze(mu,isig,4,l) = jze(mu,isig,1,l) + jze(mu,isig,2,l)
     &                          + jze(mu,isig,3,l)
               jzm(mu,isig,4,l) = jzm(mu,isig,1,l) + jzm(mu,isig,2,l)
     &                          + jzm(mu,isig,3,l)
               jzt(mu,isig,4,l) = jzt(mu,isig,1,l) + jzt(mu,isig,2,l)
     &                          + jzt(mu,isig,3,l)
            enddo

c contract with atozz/ztozz tensor

            m1 = -dotcc(azmzt(0,1),jze(0,isig,4,l))
     1           -dotcc(azezt(0,1),jzm(0,isig,4,l))
     2           -dotcc(azezm(0,1),jzt(0,isig,4,l))
            m2 = -dotcc(zzmzt(0,1),jze(0,isig,4,l))
     1           -dotcc(zzezt(0,1),jzm(0,isig,4,l))
     2           -dotcc(zzezm(0,1),jzt(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,4) = (m1*clr(ifl(k),1,isig)
     1                         +  m2*clr(ifl(k),2,isig))
     2                              * clr(ifl(k),2,isig)
            enddo

         enddo
      enddo

c qq -> A/Z(-> Z Z) Z g

      do l = 1,2

         call curr6(1,psize(1,-1,2),fqze(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jze(0,-1,1,l) )
         call curr6(1,psizeg(1,-1,2,l),pzeg(0,2),
     &           psi(1,-1,1),p(0,1), jze(0,-1,2,l) )
         call curr6(1,braketgze(1,-1,2,l),pgze(0,2),
     &           psi(1,-1,1),p(0,1), jze(0,-1,3,l) )

         call curr6(1,psizmu(1,-1,2),fqzmu(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jzm(0,-1,1,l) )
         call curr6(1,psizmug(1,-1,2,l),pzmug(0,2),
     &           psi(1,-1,1),p(0,1), jzm(0,-1,2,l) )
         call curr6(1,braketgzmu(1,-1,2,l),pgzmu(0,2),
     &           psi(1,-1,1),p(0,1), jzm(0,-1,3,l) )

         call curr6(1,psiztau(1,-1,2),fqztau(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jzt(0,-1,1,l) )
         call curr6(1,psiztaug(1,-1,2,l),pztaug(0,2),
     &           psi(1,-1,1),p(0,1), jzt(0,-1,2,l) )
         call curr6(1,braketgztau(1,-1,2,l),pgztau(0,2),
     &           psi(1,-1,1),p(0,1), jzt(0,-1,3,l) )

         do isig = -1,1,2

            do mu = 0,5
               jze(mu,isig,4,l) = jze(mu,isig,1,l) + jze(mu,isig,2,l)
     &                          + jze(mu,isig,3,l)
               jzm(mu,isig,4,l) = jzm(mu,isig,1,l) + jzm(mu,isig,2,l)
     &                          + jzm(mu,isig,3,l)
               jzt(mu,isig,4,l) = jzt(mu,isig,1,l) + jzt(mu,isig,2,l)
     &                          + jzt(mu,isig,3,l)
            enddo

c contract with atozz/ztozz tensor

            m1 = -dotcc(azmzt(0,1),jze(0,isig,4,l))
     1           -dotcc(azezt(0,1),jzm(0,isig,4,l))
     2           -dotcc(azezm(0,1),jzt(0,isig,4,l))
            m2 = -dotcc(zzmzt(0,1),jze(0,isig,4,l))
     1           -dotcc(zzezt(0,1),jzm(0,isig,4,l))
     2           -dotcc(zzezm(0,1),jzt(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,5) = (m1*clr(ifl(k),1,isig)
     1                         +  m2*clr(ifl(k),2,isig))
     2                              * clr(ifl(k),2,isig)
            enddo

         enddo
      enddo

      endif !NW = 0

c  Calculation of processes qq -> A/Z g -> Z Z Z g ---------------------

      do l = 1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))

         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqagq(0,-1,l))

         do isig = -1,1,2

            do mu = 0,5
               jqqa(mu,isig,l) = jqgaq(mu,isig,l) + jqagq(mu,isig,l)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig)
            enddo

c contract with atowwz/ztowwz tensor

            if (NW.eq.1) then
               mz = dotcc(zzzz(0,1),jqqa(0,isig,l))
               do k = 1,2
                  mat(k,isig,l,6)  = (0.0d0,0.0d0)
                  mat(k,isig,l,7)  = mz*clr(ifl(k),2,isig)
               enddo
            else !(NW = 0)
               ma = dotcc(azzz(0,1),jqqa(0,isig,l))
               mz = dotcc(zzzz(0,1),jqqa(0,isig,l))
               do k = 1,2
                  mat(k,isig,l,6)  = ma*clr(ifl(k),1,isig)
                  mat(k,isig,l,7)  = mz*clr(ifl(k),2,isig)
               enddo
            endif !NW

         enddo  !isig
      enddo !l = 1,2

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph
#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,7
               do k = 1,2

                  m1kb(k,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(k,isig,i) = 0d0 ! kartesian basis l = 2

                  do lh = -1,1,2
                     madgr(i,isig,lh,k) = 0d0
                  enddo

               enddo
            enddo
         enddo

         do isig = -1,1,2
            do lh = -1,1,2
               do k = 1,2
                  do i =1,7

                     madgr(i,isig,lh,k) = tampr(i,isig,lh,k)

                  enddo
               enddo
            enddo
         enddo

         do isig = -1,1,2
            do i=1,7
               do k=1,2

                  m1kb(k,isig,i) = (madgr(i,isig,-1,k) -
     &                  madgr(i,isig,1,k))/sqrt(2d0)
                  m2kb(k,isig,i) = -im*(madgr(i,isig,-1,k) +
     &                  madgr(i,isig,1,k))/sqrt(2d0)

               enddo
            enddo
         enddo

c print results from checks:

         do i=1,7
            do k = 1,1
               do isig = -1,1,2


            print*,"ratio for l=1, k = ",k,",isig = ",isig,",i = ",i,":",
     &          mat(k,isig,1,i)/m1kb(k,isig,i)
            print*,"mat(",k,",",isig,",1,",i,") = ",mat(k,isig,1,i)
            temp = dreal(mat(k,isig,1,i))**2+dimag(mat(k,isig,1,i))**2
            print*,temp
            print*,"m1kb(",k,",",isig,",",i,") = ",m1kb(k,isig,i)
            temp1 = dreal(m1kb(k,isig,i))**2+dimag(m1kb(k,isig,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

            print*,"ratio for l=2, k = ",k,",isig = ",isig,",i = ",i,":",
     &          mat(k,isig,2,i)/m2kb(k,isig,i)
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
#endif
c------------------------------------------------------
c------------------------------------------------------

c sum the graphs, square them and map them onto uuzzzj and ddzzzj

      do k = 1,2
         res(k) = 0
         do isig = -1,1,2
            do l = 1,2
               mm(k,isig,l) = 0
               do i = 1,7
                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,i)
               enddo

               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2

            enddo
         enddo
         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor
      enddo

      uuzzzj = res(1)
      ddzzzj = res(2)

      return
      end
