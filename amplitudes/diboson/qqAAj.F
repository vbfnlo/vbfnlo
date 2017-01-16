c****************subroutine QQAAJ *************************************
c
c	Michael Rauch, <rauch@particle.physik.uni-karlsruhe.de>
c	Last modified: February 2011
C     Francisco Campanario, francam@particle.uni-physik.uni-karlsruhe.de
c      Last modified:  April 2012
C
C  QQAA calculates the matrix elements**2 for 
C        q1 q1bar    ---->   A A g,   
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
c      return uuaa = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuaaaj,2:ddaaaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqAAj(pbar,sign, qbar, gsign, nlo, N, uuaaj, ddaaj,
     &                 ps_number)

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
      double precision  uuaaj, ddaaj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  ifl(2), lh, l, nlo, N, ps_number

      double precision fqa1(0:4,2),fqa2(0:4,2)
      double precision pq(0:4,2)
      double precision pgA1(0:4,2),pA1g(0:4,2)
      double precision pgA2(0:4,2),pA2g(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psia1(2,-1:1,2),psia2(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgA1(2,-1:1,2,2)
      double complex braketgA2(2,-1:1,2,2)
      double complex psiA1g(2,-1:1,2,2)
      double complex psiA2g(2,-1:1,2,2)


      double complex mat(2,-1:1,2), madgr(-1:1,-1:1,2)

      double complex m1kb(2,-1:1), m2kb(2,-1:1), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex m1
      double complex psi(2,-1:1,2)


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

      logical ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuaaj: u = 3
         ifl(k) = 3
         k = 2                  ! ddaaj: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do l = 1,2
               mat(k,isig,l) = (0.0d0,0.0d0)
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
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa1(0,1),
     &              atau1(1,1), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa2(0,1),
     &              atau2(1,1), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa1(0,1),
     &              atau1(1,1), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa2(0,1),
     &              atau2(1,1), psia2(1,isig,1),fqa2(0,1))
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
     &           qa1(0,1),atau1(1,1),braketgA1(1,isig,1,l),pgA1(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qa1(0,1),atau1(1,1),braketgA1(1,isig,2,l),pgA1(0,2))
            call ket2r(psia1(1,isig,1),.false.,fqa1(0,1),isig,q,eps(0,l),
     &           psiA1g(1,isig,1,l),pA1g(0,1))
            call bra2r(psia1(1,isig,2),.false.,fqa1(0,2),isig,q,eps(0,l),
     &           psiA1g(1,isig,2,l),pA1g(0,2))
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qa2(0,1),atau2(1,1),braketgA2(1,isig,1,l),pgA2(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qa2(0,1),atau2(1,1),braketgA2(1,isig,2,l),pgA2(0,2))
            call ket2r(psia2(1,isig,1),.false.,fqa2(0,1),isig,q,eps(0,l),
     &           psiA2g(1,isig,1,l),pA2g(0,1))
            call bra2r(psia2(1,isig,2),.false.,fqa2(0,2),isig,q,eps(0,l),
     &           psiA2g(1,isig,2,l),pA2g(0,2))
         enddo                   !fermion helicity
      enddo                      !gluon polarization

c qqbar -> A A g
      do l=1,2
         do isig = -1,1,2
           m1=
     1      -s1c(braketg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     2      -s1c(braketg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     2      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketg(1,isig,1,l))
     3      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketg(1,isig,1,l))
     4      -s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psi(1,isig,1))
     5      -s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psi(1,isig,1))
           
           mat(1,isig,l) =  m1 * clr(3,1,isig)**2               !uu->AA
           mat(2,isig,l) =  m1 * clr(4,1,isig)**2               !dd->AA
        enddo
      enddo

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph
      if (ldebugm) then
         do isig = -1,1,2
            do k = 1,2
               m1kb(k,isig) = 0d0 ! kartesian basis l = 1
               m2kb(k,isig) = 0d0 ! kartesian basis l = 2
               do lh = -1,1,2
                  madgr(isig,lh,k) = tampr(1,isig,lh,k)
               enddo
            enddo
         enddo

         do isig = -1,1,2
            do k=1,2
               m1kb(k,isig) = (madgr(isig,-1,k) -
     &               madgr(isig,1,k))/sqrt(2d0)
               m2kb(k,isig) = -im*(madgr(isig,-1,k) +
     &               madgr(isig,1,k))/sqrt(2d0)
               enddo
         enddo

c print results from checks:
         do k = 1,2
            do isig = -1,1,2
               print*,'mat(',k,',',isig,',1,) = ',mat(k,isig,1)
               temp = dreal(mat(k,isig,1))**2+dimag(mat(k,isig,1))**2
               print*,'m1kb(',k,',',isig,') = ',m1kb(k,isig)
               temp1 = dreal(m1kb(k,isig))**2+dimag(m1kb(k,isig))**2
               print*,'temp/temp1 = ',temp/temp1
               print*,'mat(',k,',',isig,',2,) = ',mat(k,isig,2)
               temp = dreal(mat(k,isig,2))**2+dimag(mat(k,isig,2))**2
               print*,'m2kb(',k,',',isig,') = ',m2kb(k,isig)
               temp1 = dreal(m2kb(k,isig))**2+dimag(m2kb(k,isig))**2
               print*,'temp/temp1 = ',temp/temp1
            enddo
         enddo

      endif ! end of debugging block
c-----------------------------------------------------
c------------------------------------------------------
c sum the graphs, square them and map them onto uuzaaj and ddzaaj
      do k = 1,2
         res(k) = 0
         do isig = -1,1,2
            do l = 1,2
               mm(k,isig,l) = 0
               mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l)
               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2
            enddo
         enddo
         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor
      enddo
      uuaaj = res(1)
      ddaaj = res(2)
      return
      end
