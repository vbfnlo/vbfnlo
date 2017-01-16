c****************subroutine QQZZZJ1 *************************************
c
c	Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  April 2009
c	Last modified: June 2009
C
C  QQZZZJ1 calculates the matrix elements**2 for for the Higgs contributions in
C
C        q1 q2    ---->   Z Z Z g,   
C         Z ---> f1-bar f2, Z ---> f3-bar f4, Z ---> f5-bar f6
C
C         i.e. output is |M12|^2 + |M23|^2 + |M13|^2 where
C         M12 = M(q qbar -> H Z3 g -> Z1 Z2 Z3 g).
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

      subroutine qqZZZj1(pbar,sign, qbar, gsign, uuzzzj, ddzzzj, NW)

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

      double precision pq(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex braketg(2,-1:1,2,2)
      double complex mat(2,-1:1,2,3), madgr(3,-1:1,-1:1,2)
      double complex m1kb(2,-1:1,3), m2kb(2,-1:1,3), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex mz1, mz2, mz3
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
            do i = 1,3
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

         enddo                   !fermion helicity
      enddo                      !gluon polarization


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
            mz1 = dotcc(zzzz_higgs1(0,1),jqqa(0,isig,l))
            mz2 = dotcc(zzzz_higgs2(0,1),jqqa(0,isig,l))
            mz3 = dotcc(zzzz_higgs3(0,1),jqqa(0,isig,l))
            do k = 1,2
               mat(k,isig,l,1)  = mz1*clr(ifl(k),2,isig)
               mat(k,isig,l,2)  = mz2*clr(ifl(k),2,isig)
               mat(k,isig,l,3)  = mz3*clr(ifl(k),2,isig)
            enddo
         enddo  !isig
      enddo !l = 1,2

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph
#ifdef WITH_MADGRAPH
      if (ldebugm) then

         do isig = -1,1,2
            do i = 1,3
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
                  do i =1,3

                     madgr(i,isig,lh,k) = tampr(i,isig,lh,k)

                  enddo
               enddo
            enddo
         enddo

         do isig = -1,1,2
            do i=1,3
               do k=1,2

                  m1kb(k,isig,i) = (madgr(i,isig,-1,k) -
     &                  madgr(i,isig,1,k))/sqrt(2d0)
                  m2kb(k,isig,i) = -im*(madgr(i,isig,-1,k) +
     &                  madgr(i,isig,1,k))/sqrt(2d0)

               enddo
            enddo
         enddo

c print results from checks:

         do i=1,3
            do k = 1,2
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
               do i = 1,3
!                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,i)       ! not used: only individual M^2 contributions, as in born!
                  res(k)= res(k) + dreal(mat(k,isig,l,i))**2
     &                           + dimag(mat(k,isig,l,i))**2
               enddo

!               res(k)= res(k) + dreal(mm(k,isig,l))**2
!     &                        + dimag(mm(k,isig,l))**2

            enddo
         enddo
         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor
      enddo

      uuzzzj = res(1)
      ddzzzj = res(2)

      return
      end
