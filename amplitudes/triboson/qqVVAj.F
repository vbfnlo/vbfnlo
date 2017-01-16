c****************subroutine QQVVAJ *************************************
c
c	giuseppe bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c	Initial version:  2008 January
c	Last modified: 2008 June    
C
C  QQVVAj calculates the matrix elements**2 for 
C     q1 q2    ---->   W+ W- A g,   
C     W+ ---> f1-bar f2, W- ---> f3-bar f4,
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
c       k is the process ID (1:uuwwaj,2:ddwwaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqVVAj(pbar,sign, qbar, gsign, uuwwaj, ddwwaj)

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

      DOUBLE PRECISION fpials(2:3), fpi
      PARAMETER (fpi=4d0*pi)

      double precision  pbar(0:3,2+n_v), qbar(0:4), q(0:3)
      double precision  uuwwaj, ddwwaj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  ifl(2), lh, l

      double precision fqm(0:4,2), fqp(0:4,2), fqz(0:4,2)
      double precision pq(0:4,2), pgZ(0:4,2), pZg(0:4,2)
      double precision pgWP(0:4,2), pgWM(0:4,2), pWPg(0:4,2)
      double precision pWMg(0:4,2), pwpz(0:3), pwmz(0:3) 
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psi(2,-1:1,2), psiWM(2,-1:1,2), psiWP(2,-1:1,2)
      double complex psia(2,-1:1,2), psiz(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgWP(2,-1:1,2,2), braketgWM(2,-1:1,2,2)
      double complex braketgA(2,-1:1,2,2), braketgZ(2,-1:1,2,2)
      double complex psiWPg(2,-1:1, 2,2), psiWMg(2,-1:1, 2,2)
      double complex psiAg(2,-1:1,2,2), psiZg(2,-1:1,2,2)
      double complex mat(2,-1:1,2,10), madgr(10,-1:1,-1:1,2)
      double complex matt(2,-1:1,2,10)
      double complex mm(2, -1:1,2)
      double complex m1kb(2,-1:1,10), m2kb(2,-1:1,10), temp, temp1
      double complex m1(3,2), m2(3,2)
      double complex ma, mz
      double complex ja(0:5,-1:1,4,2)
      double complex jwm(0:5,-1:1,4,2), jwp(0:5,-1:1,4,2)
      double complex jqgaq(0:5,-1:1,2), jqagq(0:5,-1:1,2)
      double complex jqqa(0:5,-1:1,2)
   

      double complex im
      double complex dotcc, s1c, dotrc
      double complex zero 
      external dotcc, s1c, dotrc
      save ifl
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
      save linit
      parameter(ldebugm = .false., ldebug=.false.)

c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         k = 1                  ! uuwwaj
         ifl(k) = 3
         k = 2                  ! ddwwaj
         ifl(k) = 4
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

      do k = 1,2
         do isig = -1,1,2
            do l = 1,2
               do i = 1,10
                  mat(k,isig,l,i) = (0.0d0, 0.0d0)
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
         pwpz(mu) = qp(mu,1) + qz(mu,1)
         pwmz(mu) = qm(mu,1) + qz(mu,1)
      enddo

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp(0,1),WP(1,1),
     &              psiWP(1,-1,2),fqp(0,2))
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qm(0,1),WM(1,1),
     &              psiWM(1,-1,2),fqm(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp(0,1),WP(1,1),
     &              psiWP(1,-1,1),fqp(0,1))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qm(0,1),WM(1,1),
     &              psiWM(1,-1,1),fqm(0,1))

      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,1),
     &                 atau(1,1), psia(1,isig,2),fqz(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,1),
     &                 ztau(1,1), psiz(1,isig,2),fqz(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,1),
     &                 atau(1,1), psia(1,isig,1),fqz(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,1),
     &                 ztau(1,1), psiz(1,isig,1),fqz(0,1))
      enddo

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2	! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck      
c            eps(mu,l) = qbar(mu)
c         enddo

c--- End of QCD gaugecheck -------------------------------------------    

         do isig = -1,1,2	! fermion helicity 

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
     &           qz(0,1),atau(1,1),braketgA(1,isig,1,l),pgZ(0,1))

            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qz(0,1),atau(1,1),braketgA(1,isig,2,l),pgZ(0,2))

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qz(0,1),ztau(1,1),braketgZ(1,isig,1,l),pgZ(0,1))

            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qz(0,1),ztau(1,1),braketgZ(1,isig,2,l),pgZ(0,2))

            call ket2r(psiA(1,isig,1),.false.,fqz(0,1),isig,q,eps(0,l),
     &           psiAg(1,isig,1,l),pZg(0,1))

            call bra2r(psiA(1,isig,2),.false.,fqz(0,2),isig,q,eps(0,l),
     &           psiAg(1,isig,2,l),pZg(0,2))

            call ket2r(psiZ(1,isig,1),.false.,fqz(0,1),isig,q,eps(0,l),
     &           psiZg(1,isig,1,l),pZg(0,1))

            call bra2r(psiZ(1,isig,2),.false.,fqz(0,2),isig,q,eps(0,l),
     &           psiZg(1,isig,2,l),pZg(0,2))

         enddo                   !fermion helicity
      enddo                      !gluon polarization

      do l=1,2

         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qp(0,1),WP(1,1),braketgWP(1,-1,1,l),pgWP(0,1)) !|W+,g,1>_l

         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qm(0,1),WM(1,1),braketgWM(1,-1,1,l),pgWM(0,1)) !|W-,g,1>_l

         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qm(0,1),WM(1,1),braketgWM(1,-1,2,l),pgWM(0,2)) !<2,g,W-|_l

         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qp(0,1),WP(1,1),braketgWP(1,-1,2,l),pgWP(0,2)) !<2,g,W+|_l

         call ket2r(psiWP(1,-1,1),.false.,fqp(0,1),-1,
     &        q,eps(0,l),psiWPg(1,-1,1,l),pWPg(0,1))

         call bra2r(psiWP(1,-1,2),.false.,fqp(0,2),-1,
     &        q,eps(0,l),psiWPg(1,-1,2,l),pWPg(0,2)) !<2,W+,g|_l

         call ket2r(psiWM(1,-1,1),.false.,fqm(0,1),-1,
     &        q,eps(0,l),psiWMg(1,-1,1,l),pWMg(0,1))

         call bra2r(psiWM(1,-1,2),.false.,fqm(0,2),-1,
     &        q,eps(0,l),psiWMg(1,-1,2,l),pWMg(0,2)) !<2,W-,g|_l

      enddo

c Calculation of processes qq -> W W A g ---------------------

c   u u -> W W A g

      do l=1,2

         m1(1,l) =
     &      -s1c(psiA(1,-1,2),WM(1,1),.true.,-1,braketgWP(1,-1,1,l))
     &      -s1c(psiA(1,-1,2),WM(1,1),.true.,-1,psiWPg(1,-1,1,l))
     &      -s1c(psiAg(1,-1,2,l),WM(1,1),.true.,-1,psiWP(1,-1,1))
     &      -s1c(braketgA(1,-1,2,l),WM(1,1),.true.,-1,psiWP(1,-1,1))

         m1(2,l) =
     &      -s1c(psiWM(1,-1,2),Atau(1,1),.true.,-1,braketgWP(1,-1,1,l))
     &      -s1c(psiWM(1,-1,2),Atau(1,1),.true.,-1,psiWPg(1,-1,1,l))
     &      -s1c(psiWMg(1,-1,2,l),Atau(1,1),.true.,-1,psiWP(1,-1,1))
     &      -s1c(braketgWM(1,-1,2,l),Atau(1,1),.true.,-1,psiWP(1,-1,1))

         m1(3,l) =
     &      -s1c(psiWM(1,-1,2),WP(1,1),.true.,-1,braketgA(1,-1,1,l))
     &      -s1c(psiWM(1,-1,2),WP(1,1),.true.,-1,psiAg(1,-1,1,l))
     &      -s1c(psiWMg(1,-1,2,l),WP(1,1),.true.,-1,psiA(1,-1,1))
     &      -s1c(braketgWM(1,-1,2,l),WP(1,1),.true.,-1,psiA(1,-1,1))

         mat(1,-1,l,1) = (m1(1,l) * clr(3,1,-1) + m1(2,l) * clr(4,1,-1)
     &                  + m1(3,l) * clr(3,1,-1)) * clr(3,3,-1)**2 ! uu -> WWA

      enddo

c   d d -> W W A g

      do l=1,2

         m2(1,l) =
     &      -s1c(psiA(1,-1,2),WP(1,1),.true.,-1,braketgWM(1,-1,1,l))
     &      -s1c(psiA(1,-1,2),WP(1,1),.true.,-1,psiWMg(1,-1,1,l))
     &      -s1c(psiAg(1,-1,2,l),WP(1,1),.true.,-1,psiWM(1,-1,1))
     &      -s1c(braketgA(1,-1,2,l),WP(1,1),.true.,-1,psiWM(1,-1,1))

         m2(2,l) =
     &      -s1c(psiWP(1,-1,2),Atau(1,1),.true.,-1,braketgWM(1,-1,1,l))
     &      -s1c(psiWP(1,-1,2),Atau(1,1),.true.,-1,psiWMg(1,-1,1,l))
     &      -s1c(psiWPg(1,-1,2,l),Atau(1,1),.true.,-1,psiWM(1,-1,1))
     &      -s1c(braketgWP(1,-1,2,l),Atau(1,1),.true.,-1,psiWM(1,-1,1))

         m2(3,l) =
     &      -s1c(psiWP(1,-1,2),WM(1,1),.true.,-1,braketgA(1,-1,1,l))
     &      -s1c(psiWP(1,-1,2),WM(1,1),.true.,-1,psiAg(1,-1,1,l))
     &      -s1c(psiWPg(1,-1,2,l),WM(1,1),.true.,-1,psiA(1,-1,1))
     &      -s1c(braketgWP(1,-1,2,l),WM(1,1),.true.,-1,psiA(1,-1,1))

         mat(2,-1,l,1) = (m2(1,l) * clr(4,1,-1) + m2(2,l) * clr(3,1,-1)
     &                  + m2(3,l) * clr(4,1,-1)) * clr(3,3,-1)**2 ! dd -> WWA

      enddo

c Calculation of processes qq -> V V g -> W W Z g -------------

c qq -> A A/Z(-> W W) g

      do l=1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketgA(1,-1,1,l),pgZ(0,1),ja(0,-1,1,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &              psiAg(1,-1,1,l),pZg(0,1),ja(0,-1,2,l) )

         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &              psiA(1,-1,1),fqZ(0,1),ja(0,-1,3,l) )

         do isig = -1,1,2

            do mu = 0,5
	       ja(mu,isig,4,l) = ja(mu,isig,1,l) + ja(mu,isig,2,l)
     &                          + ja(mu,isig,3,l)
            enddo

c contract with atoww/ztoww tensor

            m1(1,l) = -dotcc(aww(0,1),ja(0,isig,4,l))
            m2(1,l) = -dotcc(zww(0,1),ja(0,isig,4,l))

            do k = 1,2
               mat(k,isig,l,2) = (m1(1,l)*clr(ifl(k),1,isig)
     1                           + m2(1,l)*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)  ! qq -> A A/Z(->WW)
            enddo

         enddo             ! fermion helicity
      enddo                ! gluon polarization

c qq -> A/Z(-> W W) A g

      do l=1,2

         call curr6(1,psiA(1,-1,2),fqz(0,2),
     &              braketg(1,-1,1,l),pq(0,1),ja(0,-1,1,l) )

         call curr6(1,psiAg(1,-1,2,l),pZg(0,2),
     &              psi(1,-1,1),p(0,1),ja(0,-1,2,l) )

         call curr6(1,braketgA(1,-1,2,l),pgZ(0,2),
     &              psi(1,-1,1),p(0,1),ja(0,-1,3,l) )

         do isig = -1,1,2

            do mu = 0,5
	       ja(mu,isig,4,l) = ja(mu,isig,1,l) + ja(mu,isig,2,l)
     &                          + ja(mu,isig,3,l)
            enddo

c contract with atoww/ztoww tensor

            m1(1,l) = -dotcc(aww(0,1),ja(0,isig,4,l))
            m2(1,l) = -dotcc(zww(0,1),ja(0,isig,4,l))

            do k = 1,2
               mat(k,isig,l,3) = (m1(1,l)*clr(ifl(k),1,isig)
     1                           + m2(1,l)*clr(ifl(k),2,isig))
     2                           * clr(ifl(k),1,isig)  ! qq -> A/Z(->WW) A
            enddo

         enddo             ! fermion helicity
      enddo                ! gluon polarization

c uu -> W+(-> W Z) W- g

      do l=1,2

         call curr6(-1,psiWM(1,-1,2),fqm(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwm(0,-1,1,l) )

         call curr6(-1,psiWMg(1,-1,2,l),pWMg(0,2),
     &              psi(1,-1,1),p(0,1),jwm(0,-1,2,l) )

         call curr6(-1,braketgWM(1,-1,2,l),pgWM(0,2),
     &              psi(1,-1,1),p(0,1),jwm(0,-1,3,l) )

         do mu = 0,5
            jwm(mu,-1,4,l) = jwm(mu,-1,1,l) + jwm(mu,-1,2,l)
     &                          + jwm(mu,-1,3,l)
         enddo

c contract with wptowa tensor

         m1(1,l) = -dotcc(wpwa(0,1),jwm(0,-1,4,l)) ! uu -> W+(->WA) W- g

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator

         m1(2,l) = zm2i(3)*dotrc(pwpz(0),wpwa(0,1))
     &                    *dotrc(pwpz(0),jwm(0,-1,4,l))

         mat(1,-1,l,4) = (m1(1,l) + m1(2,l)) * clr(3,3,-1)**2

      enddo          ! gluon polarization

c dd -> W- W+(-> W Z) g

      do l=1,2

         call curr6(-1, psi(1,-1,2), p(0,2),
     &                 braketgWM(1,-1,1,l), pgWM(0,1), jwm(0,-1,1,l) )

         call curr6(-1, psi(1,-1,2), p(0,2),
     &                 psiWMg(1,-1,1,l), pWMg(0,1), jwm(0,-1,2,l) )

         call curr6(-1, braketg(1,-1,2,l), pq(0,2),
     &                 psiWM(1,-1,1), fqm(0,1), jwm(0,-1,3,l) )

         do mu = 0,5
            jwm(mu,-1,4,l) = jwm(mu,-1,1,l) + jwm(mu,-1,2,l)
     &                          + jwm(mu,-1,3,l)
         enddo

c contract with wptowa tensor

         m1(1,l) = -dotcc(wpwa(0,1),jwm(0,-1,4,l)) ! dd -> W- W+(->WA) g

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator

         m1(2,l) = zm2i(3)*dotrc(pwpz(0),wpwa(0,1))
     &                    *dotrc(pwpz(0),jwm(0,-1,4,l))

         mat(2,-1,l,4) = (m1(1,l) + m1(2,l)) * clr(3,3,-1)**2

      enddo          ! gluon polarization

c uu -> W+ W-(-> W Z) g

      do l=1,2

         call curr6(-1, psi(1,-1,2), p(0,2),
     &                 braketgWP(1,-1,1,l), pgWP(0,1), jwp(0,-1,1,l) )

         call curr6(-1, psi(1,-1,2), p(0,2),
     &                 psiWPg(1,-1,1,l), pWPg(0,1), jwp(0,-1,2,l) )

         call curr6(-1, braketg(1,-1,2,l), pq(0,2),
     &                 psiWP(1,-1,1), fqp(0,1), jwp(0,-1,3,l) )

         do mu = 0,5
            jwp(mu,-1,4,l) = jwp(mu,-1,1,l) + jwp(mu,-1,2,l)
     &                          + jwp(mu,-1,3,l)
         enddo

c contract with wmtowa tensor

         m1(1,l) = -dotcc(wmwa(0,1),jwp(0,-1,4,l)) ! uu -> W+ W-(->WA) g

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator

         m1(2,l) = zm2i(3)*dotrc(pwmz(0),wmwa(0,1))
     &                    *dotrc(pwmz(0),jwp(0,-1,4,l))

         mat(1,-1,l,5) = (m1(1,l) + m1(2,l)) * clr(3,3,-1)**2

      enddo          ! gluon polarization

c dd -> W-(-> W Z) W+ g

      do l=1,2

         call curr6(-1,psiWP(1,-1,2),fqp(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,1,l) )

         call curr6(-1,psiWPg(1,-1,2,l),pWPg(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )

         call curr6(-1,braketgWP(1,-1,2,l),pgWP(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,3,l) )

         do mu = 0,5
            jwp(mu,-1,4,l) = jwp(mu,-1,1,l) + jwp(mu,-1,2,l)
     &                          + jwp(mu,-1,3,l)
         enddo

c contract with wmtowa tensor

         m1(1,l) = -dotcc(wmwa(0,1),jwp(0,-1,4,l)) ! dd -> W-(->WA) W+ g

c additional contributions due to q_mu q_nu/ m^2 terms in the propagator

         m1(2,l) = zm2i(3)*dotrc(pwmz(0),wmwa(0,1))
     &                    *dotrc(pwmz(0),jwp(0,-1,4,l))

         mat(2,-1,l,5) = (m1(1,l) + m1(2,l)) * clr(3,3,-1)**2

      enddo          ! gluon polarization

c Calculation of processes qq -> A/Z g -> W W A g ---------------------

      do l =1,2

         call curr6(1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jqgaq(0,-1,l))

         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jqagq(0,-1,l))

         do isig = -1,1,2

            do mu = 0,5
	       jqqa(mu,isig,l) = jqgaq(mu,isig,l) + jqagq(mu,isig,l)
c                            = (<2|gam.mu|q,1>+<2,q|gam.mu|1>)_l,isig)
            enddo

c contract with atowwa/ztowwa tensor

            ma = -dotcc(awwa(0,1),jqqa(0,isig,l))
            mz = -dotcc(zwwa(0,1),jqqa(0,isig,l))
            
            do k = 1,2          !k=1: uuwwaj, k=2: ddwwaj
               mat(k,isig,l,6)  = ma * clr(ifl(k),1,isig)
               mat(k,isig,l,7) = mz * clr(ifl(k),2,isig)
            enddo

         enddo           ! fermion helicity  
      enddo              ! gluon polarization

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then
         
         do isig = -1,1,2
            do i = 1,10
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
                  do i =1,10
                     
                     madgr(i,isig,lh,k) = tampr(i,isig,lh,k)
                     
                  enddo
               enddo
            enddo
         enddo
         
         do isig = -1,1,2
            do i=1,10
               do k=1,2
                  
	          m1kb(k,isig,i) = (madgr(i,isig,-1,k) - 
     &                 madgr(i,isig,1,k))/sqrt(2d0)       
	          m2kb(k,isig,i) = -im*(madgr(i,isig,-1,k) + 
     &                 madgr(i,isig,1,k))/sqrt(2d0)
                  
               enddo
            enddo
         enddo
         
c print results from checks:

         do i=1,8
            do k = 1,2
               do isig = -1,1,2
                  
c there is a mismatch in labeling diagrams with madgraph:
c for contributions 2 and 3 I will check the sum
                  matt(k,isig,1,i)=mat(k,isig,1,i)
                  matt(k,isig,2,i)=mat(k,isig,2,i)
                  matt(k,isig,1,8)=mat(k,isig,1,2)+mat(k,isig,1,3)
                  matt(k,isig,2,8)=mat(k,isig,2,2)+mat(k,isig,2,3)
                  print*,"mat(",k,",",isig,",1,",i,") = ",matt(k,isig,1,i)
                  temp = dreal(matt(k,isig,1,i))**2+dimag(matt(k,isig,1,i))**2
                  print*,"m1kb(",k,",",isig,",",i,") = ",m1kb(k,isig,i)
                  temp1 = dreal(m1kb(k,isig,i))**2+dimag(m1kb(k,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
                  
                  print*,"mat(",k,",",isig,",2,",i,") = ",matt(k,isig,2,i)
                  temp = dreal(matt(k,isig,2,i))**2+dimag(matt(k,isig,2,i))**2
                  print*,"m2kb(",k,",",isig,",",i,") = ",m2kb(k,isig,i)
                  temp1 = dreal(m2kb(k,isig,i))**2+dimag(m2kb(k,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
               enddo
            enddo
         enddo
         

      endif ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto uuwwaj and ddwwaj

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

      uuwwaj = res(1)
      ddwwaj = res(2)
      
      return
      end


