c****************subroutine QQZZWJ *************************************
c
c	Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version:  2008, January
c	Last modified:  2008, February
C
C  QQZZWj calculates the matrix elements**2 for 
C        q1 q2    ---->   Z Z W g,   
C         Z ---> f1-bar f2, Z ---> f3-bar f4, W ---> f5-bar f6
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
c      return uuzzw = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuzzzj,2:ddzzzj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqZZWj(pbar,sign, qbar, gsign, qqwzzj, bos, NW)

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

      integer bos
      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  qqWZZj, res(2)
      double precision  p(0:3,2+n_v), pwz1(0:4), pwz2(0:4)

      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  lh, l, u(2), d(2)

      double precision fqze(0:4,2), fqzmu(0:4,2),
     1                 pq(0:4,2), fw(0:4,2)
      double precision pgZmu(0:4,2), pZmug(0:4,2), pWg(0:4,2),
     1                 pgZe(0:4,2), pZeg(0:4,2), pgW(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psize(2,-1:1,2), psiae(2,-1:1,2),
     1               psizmu(2,-1:1,2), psiamu(2,-1:1,2),
     2               psiw(2,2)
      double complex braketg(2,-1:1,2,2), braketgW(2,2,2),
     1               braketgZmu(2,-1:1,2,2), braketgamu(2,-1:1,2,2),
     2               braketgZe(2,-1:1,2,2), braketgAe(2,-1:1,2,2)
      double complex psiZmug(2,-1:1,2,2), psiAmug(2,-1:1,2,2),
     1               psiZeg(2,-1:1,2,2), psiAeg(2,-1:1,2,2),
     2               psiWg(2,2,2)

      double complex jae1(0:5,-1:1,4,2), jam1(0:5,-1:1,4,2),
     1               jze1(0:5,-1:1,4,2), jzm1(0:5,-1:1,4,2),
     2               jae2(0:5,-1:1,4,2), jam2(0:5,-1:1,4,2),
     3               jze2(0:5,-1:1,4,2), jzm2(0:5,-1:1,4,2),
     4               jw1(0:5,-1:1,4,2), jw2(0:5,-1:1,4,2)
      double complex mat(2,-1:1,2,8), madgr(8,-1:1,-1:1,2)
      double complex m1kb(2,-1:1,8), m2kb(2,-1:1,8), temp, temp1
      double complex mm(2,2), zero
      double complex m1, m2, m3, m4, m5, m6, m7, m8
      double complex psi(2,-1:1,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2)
      double complex jqqa(0:5,-1:1,2)

      double complex im
      double complex dotcc, s1c, dotrc
      double precision dotrr, psumsq
      external dotcc, s1c, dotrr, dotrc, psumsq
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons
      double complex zm2i(2:3)
      save zm2i

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      logical ldebug, ldebugm, linit
      data linit /.true./
      save linit,u,k,d
      parameter(ldebugm = .false., ldebug=.false.)
      integer NW !NW = 1 means Narrow width approximaion
                 !NW = 0 means full calculation

c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
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
         do l = 1,2
            mat(k,-1,l,i) = (0.0d0,0.0d0)
         enddo
      enddo

c define the internal momenta
      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwz1(mu) = qwbos(mu,1) + qze(mu,1)
         pwz2(mu) = qwbos(mu,1) + qzmu(mu,1)
         q(mu) = qbar(mu)*gsign
      enddo
      pwz1(4) = pwz1(0)**2-pwz1(1)**2-pwz1(2)**2-pwz1(3)**2
      pwz2(4) = pwz2(0)**2-pwz2(1)**2-pwz2(2)**2-pwz2(3)**2

c------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )
      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qwbos(0,1),wbos(1,1),
     &              psiw(1,2),fw(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qwbos(0,1),wbos(1,1),
     &              psiw(1,1),fw(0,1))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qze(0,1),
     &              ze(1,1),psize(1,-1,2),fqze(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qzmu(0,1),
     &              zmu(1,1),psizmu(1,-1,2),fqzmu(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qze(0,1),
     &              ae(1,1),psiae(1,-1,2),fqze(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qzmu(0,1),
     &              amu(1,1),psiamu(1,-1,2),fqzmu(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qze(0,1),
     &              ze(1,1), psize(1,-1,1),fqze(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qzmu(0,1),
     &              zmu(1,1), psizmu(1,-1,1),fqzmu(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qze(0,1),
     &              ae(1,1), psiae(1,-1,1),fqze(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qzmu(0,1),
     &              amu(1,1), psiamu(1,-1,1),fqzmu(0,1))

c  Get the gluon polarization vector and the gluon emission spinors
      do l = 1,2        ! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c***** QCD gaugecheck **** eps_mu -> k_mu *** k_mu M^mu = 0 ***********
c
c         do mu = 0,3                   !QCD gaugecheck
c            eps(mu,l) = qbar(mu)
c         enddo
c
c***** End of QCD gaugecheck ******************************************

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.

         call ket2r(psi(1,-1,1),.true.,p(0,1),-1,q,eps(0,l),
     &        braketg(1,-1,1,l),pq(0,1))      ! |q,1>_l,isig

         call bra2r(psi(1,-1,2),.true.,p(0,2),-1,q,eps(0,l),
     &        braketg(1,-1,2,l),pq(0,2))      ! <2,q|_l,isig2

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash.
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)


         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &         qwbos(0,1),wbos(1,1),braketgW(1,1,l),pgW(0,1))
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &         qwbos(0,1),wbos(1,1),braketgW(1,2,l),pgW(0,2))
         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &         qzmu(0,1),amu(1,1),braketgAmu(1,-1,1,l),pgZmu(0,1))
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &         qzmu(0,1),amu(1,1),braketgAmu(1,-1,2,l),pgZmu(0,2))
         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &         qze(0,1),ae(1,1),braketgAe(1,-1,1,l),pgZe(0,1))
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &         qze(0,1),ae(1,1),braketgAe(1,-1,2,l),pgZe(0,2))
         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &           qzmu(0,1),zmu(1,1),braketgZmu(1,-1,1,l),pgZmu(0,1))
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &           qzmu(0,1),zmu(1,1),braketgZmu(1,-1,2,l),pgZmu(0,2))
         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &           qze(0,1),ze(1,1),braketgZe(1,-1,1,l),pgZe(0,1))
         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &           qze(0,1),ze(1,1),braketgZe(1,-1,2,l),pgZe(0,2))

         call ket2r(psiW(1,1),.false.,fw(0,1),-1,q,eps(0,l),
     &           psiWg(1,1,l),pWg(0,1))
         call bra2r(psiW(1,2),.false.,fw(0,2),-1,q,eps(0,l),
     &           psiWg(1,2,l),pWg(0,2))
         call ket2r(psiAmu(1,-1,1),.false.,fqzmu(0,1),-1,q,eps(0,l),
     &           psiAmug(1,-1,1,l),pZmug(0,1))
         call bra2r(psiAmu(1,-1,2),.false.,fqzmu(0,2),-1,q,eps(0,l),
     &           psiAmug(1,-1,2,l),pZmug(0,2))
         call ket2r(psiAe(1,-1,1),.false.,fqze(0,1),-1,q,eps(0,l),
     &           psiAeg(1,-1,1,l),pZeg(0,1))
         call bra2r(psiAe(1,-1,2),.false.,fqze(0,2),-1,q,eps(0,l),
     &           psiAeg(1,-1,2,l),pZeg(0,2))
         call ket2r(psiZmu(1,-1,1),.false.,fqzmu(0,1),-1,q,eps(0,l),
     &           psiZmug(1,-1,1,l),pZmug(0,1))
         call bra2r(psiZmu(1,-1,2),.false.,fqzmu(0,2),-1,q,eps(0,l),
     &           psiZmug(1,-1,2,l),pZmug(0,2))
         call ket2r(psiZe(1,-1,1),.false.,fqze(0,1),-1,q,eps(0,l),
     &           psiZeg(1,-1,1,l),pZeg(0,1))
         call bra2r(psiZe(1,-1,2),.false.,fqze(0,2),-1,q,eps(0,l),
     &           psiZeg(1,-1,2,l),pZeg(0,2))

      enddo                      !gluon polarization

c************************************************************
c
c Calculation of processes qq -> Z Z W g ********************
c
c************************************************************

c   q qbar -> Z Z W g
      do l=1,2

        m1=(-s1c(braketgze(1,-1,2,l),wbos(1,1),.true.,-1,psizmu(1,-1,1))
     1      -s1c(psizeg(1,-1,2,l),wbos(1,1),.true.,-1,psizmu(1,-1,1))
     1      -s1c(psize(1,-1,2),wbos(1,1),.true.,-1,psizmug(1,-1,1,l))
     1      -s1c(psize(1,-1,2),wbos(1,1),.true.,-1,braketgzmu(1,-1,1,l)))
     1         *clr(d(k),2,-1)*clr(u(k),2,-1)
     2    +(-s1c(braketgzmu(1,-1,2,l),wbos(1,1),.true.,-1,psize(1,-1,1))
     2      -s1c(psizmug(1,-1,2,l),wbos(1,1),.true.,-1,psize(1,-1,1))
     2      -s1c(psizmu(1,-1,2),wbos(1,1),.true.,-1,psizeg(1,-1,1,l))
     2      -s1c(psizmu(1,-1,2),wbos(1,1),.true.,-1,braketgze(1,-1,1,l)))
     2         *clr(d(k),2,-1)*clr(u(k),2,-1)
     3    +(-s1c(braketgw(1,2,l),ze(1,1),.true.,-1,psizmu(1,-1,1))
     3      -s1c(psiwg(1,2,l),ze(1,1),.true.,-1,psizmu(1,-1,1))
     3      -s1c(psiw(1,2),ze(1,1),.true.,-1,psizmug(1,-1,1,l))
     3      -s1c(psiw(1,2),ze(1,1),.true.,-1,braketgzmu(1,-1,1,l)))
     3         *clr(u(k),2,-1)**2
     4    +(-s1c(braketgw(1,2,l),zmu(1,1),.true.,-1,psize(1,-1,1))
     4      -s1c(psiwg(1,2,l),zmu(1,1),.true.,-1,psize(1,-1,1))
     4      -s1c(psiw(1,2),zmu(1,1),.true.,-1,psizeg(1,-1,1,l))
     4      -s1c(psiw(1,2),zmu(1,1),.true.,-1,braketgze(1,-1,1,l)))
     4         *clr(u(k),2,-1)**2
     5    +(-s1c(braketgze(1,-1,2,l),zmu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psizeg(1,-1,2,l),zmu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psize(1,-1,2),zmu(1,1),.true.,-1,psiwg(1,1,l))
     5      -s1c(psize(1,-1,2),zmu(1,1),.true.,-1,braketgw(1,1,l)))
     5         *clr(d(k),2,-1)**2
     6    +(-s1c(braketgzmu(1,-1,2,l),ze(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psizmug(1,-1,2,l),ze(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psizmu(1,-1,2),ze(1,1),.true.,-1,psiwg(1,1,l))
     6      -s1c(psizmu(1,-1,2),ze(1,1),.true.,-1,braketgw(1,1,l)))
     6         *clr(d(k),2,-1)**2

      if (NW.eq.0) then

        m2=(-s1c(braketgae(1,-1,2,l),wbos(1,1),.true.,-1,psizmu(1,-1,1))
     1      -s1c(psiaeg(1,-1,2,l),wbos(1,1),.true.,-1,psizmu(1,-1,1))
     1      -s1c(psiae(1,-1,2),wbos(1,1),.true.,-1,psizmug(1,-1,1,l))
     1      -s1c(psiae(1,-1,2),wbos(1,1),.true.,-1,braketgzmu(1,-1,1,l)))
     1         *clr(d(k),1,-1)*clr(u(k),2,-1)
     2    +(-s1c(braketgzmu(1,-1,2,l),wbos(1,1),.true.,-1,psiae(1,-1,1))
     2      -s1c(psizmug(1,-1,2,l),wbos(1,1),.true.,-1,psiae(1,-1,1))
     2      -s1c(psizmu(1,-1,2),wbos(1,1),.true.,-1,psiaeg(1,-1,1,l))
     2      -s1c(psizmu(1,-1,2),wbos(1,1),.true.,-1,braketgae(1,-1,1,l)))
     2         *clr(d(k),2,-1)*clr(u(k),1,-1)
     3    +(-s1c(braketgw(1,2,l),ae(1,1),.true.,-1,psizmu(1,-1,1))
     3      -s1c(psiwg(1,2,l),ae(1,1),.true.,-1,psizmu(1,-1,1))
     3      -s1c(psiw(1,2),ae(1,1),.true.,-1,psizmug(1,-1,1,l))
     3      -s1c(psiw(1,2),ae(1,1),.true.,-1,braketgzmu(1,-1,1,l)))
     3         *clr(u(k),2,-1)*clr(u(k),1,-1)
     4    +(-s1c(braketgw(1,2,l),zmu(1,1),.true.,-1,psiae(1,-1,1))
     4      -s1c(psiwg(1,2,l),zmu(1,1),.true.,-1,psiae(1,-1,1))
     4      -s1c(psiw(1,2),zmu(1,1),.true.,-1,psiaeg(1,-1,1,l))
     4      -s1c(psiw(1,2),zmu(1,1),.true.,-1,braketgae(1,-1,1,l)))
     4         *clr(u(k),2,-1)*clr(u(k),1,-1)
     5    +(-s1c(braketgae(1,-1,2,l),zmu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psiaeg(1,-1,2,l),zmu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psiae(1,-1,2),zmu(1,1),.true.,-1,psiwg(1,1,l))
     5      -s1c(psiae(1,-1,2),zmu(1,1),.true.,-1,braketgw(1,1,l)))
     5         *clr(d(k),2,-1)*clr(d(k),1,-1)
     6    +(-s1c(braketgzmu(1,-1,2,l),ae(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psizmug(1,-1,2,l),ae(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psizmu(1,-1,2),ae(1,1),.true.,-1,psiwg(1,1,l))
     6      -s1c(psizmu(1,-1,2),ae(1,1),.true.,-1,braketgw(1,1,l)))
     6         *clr(d(k),2,-1)*clr(d(k),1,-1)

        m3=(-s1c(braketgze(1,-1,2,l),wbos(1,1),.true.,-1,psiamu(1,-1,1))
     1      -s1c(psizeg(1,-1,2,l),wbos(1,1),.true.,-1,psiamu(1,-1,1))
     1      -s1c(psize(1,-1,2),wbos(1,1),.true.,-1,psiamug(1,-1,1,l))
     1      -s1c(psize(1,-1,2),wbos(1,1),.true.,-1,braketgamu(1,-1,1,l)))
     1         *clr(u(k),1,-1)*clr(d(k),2,-1)
     2    +(-s1c(braketgamu(1,-1,2,l),wbos(1,1),.true.,-1,psize(1,-1,1))
     2      -s1c(psiamug(1,-1,2,l),wbos(1,1),.true.,-1,psize(1,-1,1))
     2      -s1c(psiamu(1,-1,2),wbos(1,1),.true.,-1,psizeg(1,-1,1,l))
     2      -s1c(psiamu(1,-1,2),wbos(1,1),.true.,-1,braketgze(1,-1,1,l)))
     2         *clr(d(k),1,-1)*clr(u(k),2,-1)
     3    +(-s1c(braketgw(1,2,l),ze(1,1),.true.,-1,psiamu(1,-1,1))
     3      -s1c(psiwg(1,2,l),ze(1,1),.true.,-1,psiamu(1,-1,1))
     3      -s1c(psiw(1,2),ze(1,1),.true.,-1,psiamug(1,-1,1,l))
     3      -s1c(psiw(1,2),ze(1,1),.true.,-1,braketgamu(1,-1,1,l)))
     3         *clr(u(k),1,-1)*clr(u(k),2,-1)
     4    +(-s1c(braketgw(1,2,l),amu(1,1),.true.,-1,psize(1,-1,1))
     4      -s1c(psiwg(1,2,l),amu(1,1),.true.,-1,psize(1,-1,1))
     4      -s1c(psiw(1,2),amu(1,1),.true.,-1,psizeg(1,-1,1,l))
     4      -s1c(psiw(1,2),amu(1,1),.true.,-1,braketgze(1,-1,1,l)))
     4         *clr(u(k),1,-1)*clr(u(k),2,-1)
     5    +(-s1c(braketgze(1,-1,2,l),amu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psizeg(1,-1,2,l),amu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psize(1,-1,2),amu(1,1),.true.,-1,psiwg(1,1,l))
     5      -s1c(psize(1,-1,2),amu(1,1),.true.,-1,braketgw(1,1,l)))
     5         *clr(d(k),1,-1)*clr(d(k),2,-1)
     6    +(-s1c(braketgamu(1,-1,2,l),ze(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psiamug(1,-1,2,l),ze(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psiamu(1,-1,2),ze(1,1),.true.,-1,psiwg(1,1,l))
     6      -s1c(psiamu(1,-1,2),ze(1,1),.true.,-1,braketgw(1,1,l)))
     6         *clr(d(k),1,-1)*clr(d(k),2,-1)

        m4=(-s1c(braketgae(1,-1,2,l),wbos(1,1),.true.,-1,psiamu(1,-1,1))
     1      -s1c(psiaeg(1,-1,2,l),wbos(1,1),.true.,-1,psiamu(1,-1,1))
     1      -s1c(psiae(1,-1,2),wbos(1,1),.true.,-1,psiamug(1,-1,1,l))
     1      -s1c(psiae(1,-1,2),wbos(1,1),.true.,-1,braketgamu(1,-1,1,l)))
     1         *clr(u(k),1,-1)*clr(d(k),1,-1)
     2    +(-s1c(braketgamu(1,-1,2,l),wbos(1,1),.true.,-1,psiae(1,-1,1))
     2      -s1c(psiamug(1,-1,2,l),wbos(1,1),.true.,-1,psiae(1,-1,1))
     2      -s1c(psiamu(1,-1,2),wbos(1,1),.true.,-1,psiaeg(1,-1,1,l))
     2      -s1c(psiamu(1,-1,2),wbos(1,1),.true.,-1,braketgae(1,-1,1,l)))
     2         *clr(d(k),1,-1)*clr(u(k),1,-1)
     3    +(-s1c(braketgw(1,2,l),ae(1,1),.true.,-1,psiamu(1,-1,1))
     3      -s1c(psiwg(1,2,l),ae(1,1),.true.,-1,psiamu(1,-1,1))
     3      -s1c(psiw(1,2),ae(1,1),.true.,-1,psiamug(1,-1,1,l))
     3      -s1c(psiw(1,2),ae(1,1),.true.,-1,braketgamu(1,-1,1,l)))
     3         *clr(u(k),1,-1)**2
     4    +(-s1c(braketgw(1,2,l),amu(1,1),.true.,-1,psiae(1,-1,1))
     4      -s1c(psiwg(1,2,l),amu(1,1),.true.,-1,psiae(1,-1,1))
     4      -s1c(psiw(1,2),amu(1,1),.true.,-1,psiaeg(1,-1,1,l))
     4      -s1c(psiw(1,2),amu(1,1),.true.,-1,braketgae(1,-1,1,l)))
     4         *clr(u(k),1,-1)**2
     5    +(-s1c(braketgae(1,-1,2,l),amu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psiaeg(1,-1,2,l),amu(1,1),.true.,-1,psiw(1,1))
     5      -s1c(psiae(1,-1,2),amu(1,1),.true.,-1,psiwg(1,1,l))
     5      -s1c(psiae(1,-1,2),amu(1,1),.true.,-1,braketgw(1,1,l)))
     5         *clr(d(k),1,-1)**2
     6    +(-s1c(braketgamu(1,-1,2,l),ae(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psiamug(1,-1,2,l),ae(1,1),.true.,-1,psiw(1,1))
     6      -s1c(psiamu(1,-1,2),ae(1,1),.true.,-1,psiwg(1,1,l))
     6      -s1c(psiamu(1,-1,2),ae(1,1),.true.,-1,braketgw(1,1,l)))
     6         *clr(d(k),1,-1)**2

         endif

         if (NW.eq.0) then
            mat(k,-1,l,1) = (m1 + m2 + m3 + m4) * clr(3,3,-1)
         else if (NW.eq.1) then
            mat(k,-1,l,1) = m1 * clr(3,3,-1)
         endif

      enddo

c***************************************************************
c
c Calculation of processes qq -> V V g -> 2l + 4l g ************
c
c***************************************************************

      do l = 1,2

c calculate currents needed for W A/Ztozz g diagrams
         call curr6(-1, braketg(1,-1,2,l), pq(0,2),
     &               psiw(1,1), fw(0,1), jw1(0,-1,3,l) )
         call curr6(-1, psi(1,-1,2), p(0,2),
     &               braketgW(1,1,l), pgW(0,1), jw1(0,-1,2,l) )
         call curr6(-1, psi(1,-1,2), p(0,2),
     &               psiwg(1,1,l), pWg(0,1), jw1(0,-1,1,l) )

         call curr6(-1, psiw(1,2), fw(0,2),
     &               braketg(1,-1,1,l), pq(0,1), jw2(0,-1,3,l) )
         call curr6(-1, braketgW(1,2,l), pgW(0,2),
     &               psi(1,-1,1), p(0,1), jw2(0,-1,2,l) )
         call curr6(-1, psiwg(1,2,l), pWg(0,2),
     &               psi(1,-1,1), p(0,1), jw2(0,-1,1,l) )

         do mu = 0,5
            jw1(mu,-1,4,l) = jw1(mu,-1,1,l) + jw1(mu,-1,2,l)
     &                     + jw1(mu,-1,3,l)
            jw2(mu,-1,4,l) = jw2(mu,-1,1,l) + jw2(mu,-1,2,l)
     &                     + jw2(mu,-1,3,l)
         enddo

c qq -> W A/ZtoZZ g - contract with atozz/ztozz tensor
         m1 = -dotcc(azezm(0,1),jw1(0,-1,4,l))
         m2 = -dotcc(zzezm(0,1),jw1(0,-1,4,l))
         m3 = -dotcc(azezm(0,1),jw2(0,-1,4,l))
         m4 = -dotcc(zzezm(0,1),jw2(0,-1,4,l))

         mat(k,-1,l,2)=(m1*clr(d(k),1,-1)+m2*clr(d(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,3)=(m3*clr(u(k),1,-1)+m4*clr(u(k),2,-1))*clr(3,3,-1)

      enddo !l = 1,2

      do l = 1,2

c calculate currents needed for A/Z WtoZW g diagrams
         call curr6(-1,psiae(1,-1,2),fqze(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jae2(0,-1,1,l) )
         call curr6(-1,psiaeg(1,-1,2,l),pZeg(0,2),
     &           psi(1,-1,1),p(0,1), jae2(0,-1,2,l) )
         call curr6(-1,braketgae(1,-1,2,l),pgZe(0,2),
     &           psi(1,-1,1),p(0,1), jae2(0,-1,3,l) )

         call curr6(-1,psize(1,-1,2),fqze(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jze2(0,-1,1,l) )
         call curr6(-1,psizeg(1,-1,2,l),pzeg(0,2),
     &           psi(1,-1,1),p(0,1), jze2(0,-1,2,l) )
         call curr6(-1,braketgze(1,-1,2,l),pgze(0,2),
     &           psi(1,-1,1),p(0,1), jze2(0,-1,3,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgae(1,-1,1,l),pgZe(0,1), jae1(0,-1,1,l) )
         call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiaeg(1,-1,1,l),pZeg(0,1), jae1(0,-1,2,l) )
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &           psiae(1,-1,1),fqze(0,1), jae1(0,-1,3,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgze(1,-1,1,l),pgze(0,1), jze1(0,-1,1,l) )
         call curr6(-1,psi(1,-1,2),p(0,2),
     &           psizeg(1,-1,1,l),pzeg(0,1), jze1(0,-1,2,l) )
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &           psize(1,-1,1),fqze(0,1), jze1(0,-1,3,l) )

         call curr6(-1,psiamu(1,-1,2),fqzmu(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jam2(0,-1,1,l) )
         call curr6(-1,psiamug(1,-1,2,l),pZmug(0,2),
     &           psi(1,-1,1),p(0,1), jam2(0,-1,2,l) )
         call curr6(-1,braketgamu(1,-1,2,l),pgZmu(0,2),
     &           psi(1,-1,1),p(0,1), jam2(0,-1,3,l) )

         call curr6(-1,psizmu(1,-1,2),fqzmu(0,2),
     &           braketg(1,-1,1,l),pq(0,1), jzm2(0,-1,1,l) )
         call curr6(-1,psizmug(1,-1,2,l),pzmug(0,2),
     &           psi(1,-1,1),p(0,1), jzm2(0,-1,2,l) )
         call curr6(-1,braketgzmu(1,-1,2,l),pgzmu(0,2),
     &           psi(1,-1,1),p(0,1), jzm2(0,-1,3,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgamu(1,-1,1,l),pgZmu(0,1), jam1(0,-1,1,l) )
         call curr6(-1,psi(1,-1,2),p(0,2),
     &           psiamug(1,-1,1,l),pZmug(0,1), jam1(0,-1,2,l) )
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &           psiamu(1,-1,1),fqzmu(0,1), jam1(0,-1,3,l) )

         call curr6(-1,psi(1,-1,2),p(0,2),
     &           braketgzmu(1,-1,1,l),pgzmu(0,1), jzm1(0,-1,1,l) )
         call curr6(-1,psi(1,-1,2),p(0,2),
     &           psizmug(1,-1,1,l),pzmug(0,1), jzm1(0,-1,2,l) )
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &           psizmu(1,-1,1),fqzmu(0,1), jzm1(0,-1,3,l) )

         do mu = 0,5
            jae1(mu,-1,4,l) = jae1(mu,-1,1,l) + jae1(mu,-1,2,l)
     &                       + jae1(mu,-1,3,l)
            jae2(mu,-1,4,l) = jae2(mu,-1,1,l) + jae2(mu,-1,2,l)
     &                       + jae2(mu,-1,3,l)
            jam1(mu,-1,4,l) = jam1(mu,-1,1,l) + jam1(mu,-1,2,l)
     &                       + jam1(mu,-1,3,l)
            jam2(mu,-1,4,l) = jam2(mu,-1,1,l) + jam2(mu,-1,2,l)
     &                       + jam2(mu,-1,3,l)
            jze1(mu,-1,4,l) = jze1(mu,-1,1,l) + jze1(mu,-1,2,l)
     &                       + jze1(mu,-1,3,l)
            jze2(mu,-1,4,l) = jze2(mu,-1,1,l) + jze2(mu,-1,2,l)
     &                       + jze2(mu,-1,3,l)
            jzm1(mu,-1,4,l) = jzm1(mu,-1,1,l) + jzm1(mu,-1,2,l)
     &                       + jzm1(mu,-1,3,l)
            jzm2(mu,-1,4,l) = jzm2(mu,-1,1,l) + jzm2(mu,-1,2,l)
     &                       + jzm2(mu,-1,3,l)
         enddo

c contract with wtozw1 tensors
         m1 = -dotcc(wzw2(0,1),jae2(0,-1,4,l))
         m2 = -dotcc(wzw2(0,1),jze2(0,-1,4,l))
         m3 = -dotcc(wzw2(0,1),jae1(0,-1,4,l))
         m4 = -dotcc(wzw2(0,1),jze1(0,-1,4,l))
         m5 = -dotcc(wzw1(0,1),jam2(0,-1,4,l))
         m6 = -dotcc(wzw1(0,1),jzm2(0,-1,4,l))
         m7 = -dotcc(wzw1(0,1),jam1(0,-1,4,l))
         m8 = -dotcc(wzw1(0,1),jzm1(0,-1,4,l))

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator
         m1 = m1+zm2i(3)*dotrc(pwz2(0),wzw2(0,1))*dotrc(pwz2(0),jae2(0,-1,4,l))
         m2 = m2+zm2i(3)*dotrc(pwz2(0),wzw2(0,1))*dotrc(pwz2(0),jze2(0,-1,4,l))
         m3 = m3+zm2i(3)*dotrc(pwz2(0),wzw2(0,1))*dotrc(pwz2(0),jae1(0,-1,4,l))
         m4 = m4+zm2i(3)*dotrc(pwz2(0),wzw2(0,1))*dotrc(pwz2(0),jze1(0,-1,4,l))
         m5 = m5+zm2i(3)*dotrc(pwz1(0),wzw1(0,1))*dotrc(pwz1(0),jam2(0,-1,4,l))
         m6 = m6+zm2i(3)*dotrc(pwz1(0),wzw1(0,1))*dotrc(pwz1(0),jzm2(0,-1,4,l))
         m7 = m7+zm2i(3)*dotrc(pwz1(0),wzw1(0,1))*dotrc(pwz1(0),jam1(0,-1,4,l))
         m8 = m8+zm2i(3)*dotrc(pwz1(0),wzw1(0,1))*dotrc(pwz1(0),jzm1(0,-1,4,l))

         if (NW.eq.0) then
         mat(k,-1,l,4)=(m1*clr(d(k),1,-1)+m2*clr(d(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,5)=(m3*clr(u(k),1,-1)+m4*clr(u(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,6)=(m5*clr(d(k),1,-1)+m6*clr(d(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,7)=(m7*clr(u(k),1,-1)+m8*clr(u(k),2,-1))*clr(3,3,-1)
         else if (NW.eq.1) then
         mat(k,-1,l,4)=(m2*clr(d(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,5)=(m4*clr(u(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,6)=(m6*clr(d(k),2,-1))*clr(3,3,-1)
         mat(k,-1,l,7)=(m8*clr(u(k),2,-1))*clr(3,3,-1)
         endif

      enddo ! l = 1,2

c***************************************************************
c
c Calculation of processes qq -> W g -> Z Z W g ****************
c
c***************************************************************
      do l = 1,2

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqagq(0,-1,l))

         do mu = 0,5
            jqqa(mu,-1,l) = jqgaq(mu,-1,l) + jqagq(mu,-1,l)
         enddo

c contract with wtozzw tensor
         m1 = dotcc(wzzw(0,1),jqqa(0,-1,l))

         mat(k,-1,l,8)  = m1*clr(3,3,-1)

      enddo !l = 1,2

c**************************************************************
c
c comparison with madgraph ************************************
c
c**************************************************************
#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do i = 1,8
            m1kb(k,-1,i) = 0d0 ! kartesian basis l = 1
            m2kb(k,-1,i) = 0d0 ! kartesian basis l = 2
            do lh = -1,1,2
               madgr(i,-1,lh,k) = 0d0
            enddo
         enddo

         do lh = -1,1,2
            do i =1,8
               madgr(i,-1,lh,k) = tampr(i,-1,lh,k)
            enddo
         enddo

         do i=1,8
            m1kb(k,-1,i) = (madgr(i,-1,-1,k) -
     &                        madgr(i,-1,1,k))/sqrt(2d0)
            m2kb(k,-1,i) = -im*(madgr(i,-1,-1,k) +
     &                            madgr(i,-1,1,k))/sqrt(2d0)
         enddo

c print results from checks:

         do i=1,8

            print*,"ratio for l=1, k = ",k,",isig = -1, i = ",i,":",
     &          mat(k,-1,1,i)/m1kb(k,-1,i)
            print*,"mat(",k,",",-1,",1,",i,") = ",mat(k,-1,1,i)
            temp = dreal(mat(k,-1,1,i))**2+dimag(mat(k,-1,1,i))**2
            print*,temp
            print*,"m1kb(",k,",",-1,",",i,") = ",m1kb(k,-1,i)
            temp1 = dreal(m1kb(k,-1,i))**2+dimag(m1kb(k,-1,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

            print*,"ratio for l=2, k = ",k,",isig = -1, i = ",i,":",
     &          mat(k,-1,2,i)/m2kb(k,-1,i)
            print*,"mat(",k,",",-1,",2,",i,") = ",mat(k,-1,2,i)
            temp = dreal(mat(k,-1,2,i))**2+dimag(mat(k,-1,2,i))**2
            print*,temp
            print*,"m2kb(",k,",",-1,",",i,") = ",m2kb(k,-1,i)
            temp1 = dreal(m2kb(k,-1,i))**2+dimag(m2kb(k,-1,i))**2
            print*,temp1
            print*,"temp/temp1 = ",temp/temp1
            print*," "

         enddo

      endif ! end of debugging block
#endif
c------------------------------------------------------
c------------------------------------------------------

c sum the graphs, square them and map them onto duzzzj and udzzzj
      res(k) = 0
      do l = 1,2
         mm(k,l) = 0
         do i = 1,8
            mm(k,l) = mm(k,l) +  mat(k,-1,l,i)
         enddo
         res(k)= res(k) + dreal(mm(k,l))**2 + dimag(mm(k,l))**2
      enddo
      res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor

      qqWZZj = res(k)

      return
      end
