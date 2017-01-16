c****************subroutine QQWZAJ *************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Initial version:  2009, January
c       Last modified: 2009, January    
c       Adapted from qqVVAj.f by G. Bozzi
C
C  QQWZAj calculates the matrix elements**2 for 
C     q1 q2    ---->   W+ Z A g,   
C     W+ ---> f1-bar f2, Z ---> f3-bar f4,
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
c      return udwza = |M_real|^2   etc.
c
c       fpials is attached only at the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:udwzaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      SUBROUTINE qqWZAj(pbar,sign, qbar, gsign, udwzaj, bos, ps_number)

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
      double precision  udwzaj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, isig, ps_number
      integer  ifl(2), lh, l, bos

      double precision fqa(0:4,2), fqp(0:4,2), fqz(0:4,2)
      double precision pq(0:4,2), pgZ(0:4,2), pZg(0:4,2)
      double precision pgWP(0:4,2), pgA(0:4,2), pWPg(0:4,2)
      double precision pAg(0:4,2), pwpz(0:3), pwpa(0:3), pzpa(0:3)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psi(2,-1:1,2), psiWP(2,-1:1,2)
      double complex psia(2,-1:1,2), psiz(2,-1:1,2), psiaz(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgWP(2,-1:1,2,2), braketgAZ(2,-1:1,2,2)
      double complex braketgA(2,-1:1,2,2), braketgZ(2,-1:1,2,2)
      double complex psiWPg(2,-1:1, 2,2)
      double complex psiAg(2,-1:1,2,2), psiZg(2,-1:1,2,2)
      double complex psiAZg(2,-1:1,2,2)
      double complex mat(2,-1:1,2,10), madgr(10,-1:1,-1:1,2)

      double complex mm(2, -1:1,2)
      double complex m1kb(2,-1:1,10), m2kb(2,-1:1,10), temp, temp1
      double complex m1(3,2), m2(3,2), m3(3,2), m4(3,2)
      double complex mw

      double complex jwp(0:5,-1:1,4,2)


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
         if (bos.eq.321) then
           ifl(1) = 3             ! udwza: u = 3; d = 4
           ifl(2) = 4
         else if (bos.eq.421) then
           ifl(1) = 4             
           ifl(2) = 3
         else 
           print *, "Wrong boson ids:", bos
           stop
         endif
         zm2i(2) = 1d0/dcmplx(xm2(2),-xmg(2))
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

c only left-handed quarks can appear due to W and m_q=0
      isig = -1

c         do isig = -1,1,2
            do l = 1,2
               do i = 1,3
                  mat(1,isig,l,i) = (0.0d0, 0.0d0)
               enddo
            enddo
c         enddo

c define the internal momenta

      do mu = 0,3
         do i = 1,2+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         q(mu) = qbar(mu)*gsign
         pwpz(mu) = qp(mu,1) + qz(mu,1)
         pwpa(mu) = qp(mu,1) + qa(mu,1)
         pzpa(mu) = qz(mu,1) + qa(mu,1)
      enddo

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qp(0,1),
     &              wp(1,1),psiwp(1,-1,2),fqp(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qp(0,1),
     &              wp(1,1),psiwp(1,-1,1),fqp(0,1))

c      do isig =-1,1,2

         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa(0,1),
     &                 wa(1,1), psia(1,isig,2),fqa(0,2)) ! <d|A
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,1),
     &                 ztau(1,1), psiz(1,isig,2),fqz(0,2)) ! <d|Z
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qz(0,1),
     &                 atau(1,1), psiaz(1,isig,2),fqz(0,2)) ! <d|A(Z)
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa(0,1),
     &                 wa(1,1), psia(1,isig,1),fqa(0,1)) ! A|u>
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,1),
     &                 ztau(1,1), psiz(1,isig,1),fqz(0,1)) ! Z|u>
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qz(0,1),
     &                 atau(1,1), psiaz(1,isig,1),fqz(0,1)) ! A(Z)|u>
c      enddo

c  Get the gluon polarization vector and the gluon emission spinors

      do l = 1,2  ! 2 gluon polarizations

         call polvec(qbar,l,eps(0,l))  ! get gluon pol.vectors

c---- QCD gaugecheck --- eps_mu -> k_mu --- k_mu M^mu = 0 ------------

c         do mu = 0,3                   !QCD gaugecheck      
c            eps(mu,l) = qbar(mu)
c         enddo

c--- End of QCD gaugecheck -------------------------------------------    

c         do isig = -1,1,2 ! fermion helicity 

            call bra2r(psi(1,isig,2),.true.,p(0,2),isig,q,eps(0,l),
     &           braketg(1,isig,2,l),pq(0,2))      ! <2,q|_l,isig2

            call ket2r(psi(1,isig,1),.true.,p(0,1),isig,q,eps(0,l),
     &           braketg(1,isig,1,l),pq(0,1))      ! |q,1>_l,isig

c     braketg contains the free quark spinors multiplied by a fermion
c     propagator and a gluon eps_slash. 
c     NOTATION: braketg(2 component spinor, isig =-1 or 1 (fermion hel.),
c     fermion ID = 1:4, gluon polarization l=1:2)

c         enddo

         call bra2c(braketg(1,-1,2,l),.false.,pq(0,2),-1,
     &        qp(0,1),WP(1,1),braketgWP(1,-1,2,l),pgWP(0,2)) !<2,g,W+|_l

         call ket2c(braketg(1,-1,1,l),.false.,pq(0,1),-1,
     &        qp(0,1),WP(1,1),braketgWP(1,-1,1,l),pgWP(0,1)) !|W+,g,1>_l

         call bra2r(psiWP(1,-1,2),.false.,fqp(0,2),-1,
     &        q,eps(0,l),psiWPg(1,-1,2,l),pWPg(0,2)) !<2,W+,g|_l

         call ket2r(psiWP(1,-1,1),.false.,fqp(0,1),-1,
     &        q,eps(0,l),psiWPg(1,-1,1,l),pWPg(0,1))

c         do isig = -1,1,2 ! fermion helicity 

c     NOTES for bras and kets: .true. if psi is a 2-spinor of the chi
c     form as output by psi0m, .false. otherwise.  the last entry is
c     the sum of the two momenta (p plus q) and effectively the
c     momentum of the new spinor.
 
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qa(0,1),wa(1,1),braketgA(1,isig,2,l),pgA(0,2)) ! <d,g|A

            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qz(0,1),ztau(1,1),braketgZ(1,isig,2,l),pgZ(0,2)) ! <d,g|Z

            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qz(0,1),atau(1,1),braketgAZ(1,isig,2,l),pgZ(0,2)) ! <d,g|z

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig, 
     &           qa(0,1),wa(1,1),braketgA(1,isig,1,l),pgA(0,1)) ! A|g,u>

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qz(0,1),ztau(1,1),braketgZ(1,isig,1,l),pgZ(0,1)) ! Z|g,u>

            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qz(0,1),atau(1,1),braketgAZ(1,isig,1,l),pgZ(0,1)) ! z|g,u>

            call bra2r(psiA(1,isig,2),.false.,fqa(0,2),isig,
     &           q,eps(0,l),psiAg(1,isig,2,l),pAg(0,2)) ! <d,A|g

            call bra2r(psiZ(1,isig,2),.false.,fqz(0,2),isig,
     &           q,eps(0,l),psiZg(1,isig,2,l),pZg(0,2)) ! <d,Z|g

            call bra2r(psiAZ(1,isig,2),.false.,fqz(0,2),isig,
     &           q,eps(0,l),psiAZg(1,isig,2,l),pZg(0,2)) ! <d,z|g

            call ket2r(psiA(1,isig,1),.false.,fqa(0,1),isig,
     &           q,eps(0,l),psiAg(1,isig,1,l),pAg(0,1)) ! g|A,u>

            call ket2r(psiZ(1,isig,1),.false.,fqz(0,1),isig,
     &           q,eps(0,l),psiZg(1,isig,1,l),pZg(0,1)) ! g|Z,u>

            call ket2r(psiAZ(1,isig,1),.false.,fqz(0,1),isig,
     &           q,eps(0,l),psiAZg(1,isig,1,l),pZg(0,1)) ! g|Z,u>

c         enddo                   !fermion helicity
      enddo                      !gluon polarization

c Calculation of processes qq -> W Z A g ---------------------

c   u d -> W Z A g

      do l=1,2

         m1(1,l) =
     &     (-s1c(braketgWP(1,-1,2,l),ztau(1,1),.true.,-1,psiA(1,-1,1)) !AZWg
     &      -s1c(braketgWP(1,-1,2,l),wa(1,1),.true.,-1,psiZ(1,-1,1))   !ZAWg
     &      -s1c(psiWPg(1,-1,2,l),ztau(1,1),.true.,-1,psiA(1,-1,1))    !AZgW
     &      -s1c(psiWPg(1,-1,2,l),wa(1,1),.true.,-1,psiZ(1,-1,1))      !ZAgW
     &      -s1c(psiWP(1,-1,2),ztau(1,1),.true.,-1,psiAg(1,-1,1,l))    !AgZW
     &      -s1c(psiWP(1,-1,2),wa(1,1),.true.,-1,psiZg(1,-1,1,l))      !ZAgW
     &      -s1c(psiWP(1,-1,2),ztau(1,1),.true.,-1,braketgA(1,-1,1,l)) !gAZW
     &      -s1c(psiWP(1,-1,2),wa(1,1),.true.,-1,braketgZ(1,-1,1,l))   !gZAW
     A     )*clr(ifl(1),2,-1)*clr(ifl(1),1,-1)+
     &     (-s1c(braketgZ(1,-1,2,l),wp(1,1),.true.,-1,psiA(1,-1,1))    !AWZg
     &      -s1c(psiZg(1,-1,2,l),wp(1,1),.true.,-1,psiA(1,-1,1))       !AWgZ
     &      -s1c(psiZ(1,-1,2),wp(1,1),.true.,-1,psiAg(1,-1,1,l))       !AgWZ
     &      -s1c(psiZ(1,-1,2),wp(1,1),.true.,-1,braketgA(1,-1,1,l))    !gAWZ
     B     )*clr(ifl(2),2,-1)*clr(ifl(1),1,-1)+
     &     (-s1c(braketgA(1,-1,2,l),wp(1,1),.true.,-1,psiZ(1,-1,1))    !ZWAg
     &      -s1c(psiAg(1,-1,2,l),wp(1,1),.true.,-1,psiZ(1,-1,1))       !ZWgA
     &      -s1c(psiA(1,-1,2),wp(1,1),.true.,-1,psiZg(1,-1,1,l))       !ZgWA
     &      -s1c(psiA(1,-1,2),wp(1,1),.true.,-1,braketgZ(1,-1,1,l))    !gZWA
     C     )*clr(ifl(1),2,-1)*clr(ifl(2),1,-1)+
     &     (-s1c(braketgZ(1,-1,2,l),wa(1,1),.true.,-1,psiWP(1,-1,1))   !WAZg
     &      -s1c(braketgA(1,-1,2,l),ztau(1,1),.true.,-1,psiWP(1,-1,1)) !WZAg
     &      -s1c(psiZg(1,-1,2,l),wa(1,1),.true.,-1,psiWP(1,-1,1))      !WAgZ
     &      -s1c(psiAg(1,-1,2,l),ztau(1,1),.true.,-1,psiWP(1,-1,1))    !WZgA
     &      -s1c(psiZ(1,-1,2),wa(1,1),.true.,-1,psiWPg(1,-1,1,l))      !WgAZ
     &      -s1c(psiA(1,-1,2),ztau(1,1),.true.,-1,psiWPg(1,-1,1,l))    !WgZA
     &      -s1c(psiZ(1,-1,2),wa(1,1),.true.,-1,braketgWP(1,-1,1,l))   !gWAZ
     &      -s1c(psiA(1,-1,2),ztau(1,1),.true.,-1,braketgWP(1,-1,1,l)) !gWZA
     D     )*clr(ifl(2),2,-1)*clr(ifl(2),1,-1)

         m2(1,l) =
     &     (-s1c(braketgWP(1,-1,2,l),atau(1,1),.true.,-1,psiA(1,-1,1)) !AzWg
     &      -s1c(braketgWP(1,-1,2,l),wa(1,1),.true.,-1,psiAZ(1,-1,1))  !zAWg
     &      -s1c(psiWPg(1,-1,2,l),atau(1,1),.true.,-1,psiA(1,-1,1))    !AzgW
     &      -s1c(psiWPg(1,-1,2,l),wa(1,1),.true.,-1,psiAZ(1,-1,1))     !zAgW
     &      -s1c(psiWP(1,-1,2),atau(1,1),.true.,-1,psiAg(1,-1,1,l))    !AgzW
     &      -s1c(psiWP(1,-1,2),wa(1,1),.true.,-1,psiAZg(1,-1,1,l))     !zAgW
     &      -s1c(psiWP(1,-1,2),atau(1,1),.true.,-1,braketgA(1,-1,1,l)) !gAzW
     &      -s1c(psiWP(1,-1,2),wa(1,1),.true.,-1,braketgAZ(1,-1,1,l))  !gzAW
     A     )*clr(ifl(1),1,-1)*clr(ifl(1),1,-1)+
     &     (-s1c(braketgAZ(1,-1,2,l),wp(1,1),.true.,-1,psiA(1,-1,1))   !AWzg
     &      -s1c(psiAZg(1,-1,2,l),wp(1,1),.true.,-1,psiA(1,-1,1))      !AWgz
     &      -s1c(psiAZ(1,-1,2),wp(1,1),.true.,-1,psiAg(1,-1,1,l))      !AgWz
     &      -s1c(psiAZ(1,-1,2),wp(1,1),.true.,-1,braketgA(1,-1,1,l))   !gAWz
     B     )*clr(ifl(2),1,-1)*clr(ifl(1),1,-1)+
     &     (-s1c(braketgA(1,-1,2,l),wp(1,1),.true.,-1,psiAZ(1,-1,1))   !zWAg
     &      -s1c(psiAg(1,-1,2,l),wp(1,1),.true.,-1,psiAZ(1,-1,1))      !zWgA
     &      -s1c(psiA(1,-1,2),wp(1,1),.true.,-1,psiAZg(1,-1,1,l))      !zgWA
     &      -s1c(psiA(1,-1,2),wp(1,1),.true.,-1,braketgAZ(1,-1,1,l))   !gzWA
     C     )*clr(ifl(1),1,-1)*clr(ifl(2),1,-1)+
     &     (-s1c(braketgAZ(1,-1,2,l),wa(1,1),.true.,-1,psiWP(1,-1,1))  !WAzg
     &      -s1c(braketgA(1,-1,2,l),atau(1,1),.true.,-1,psiWP(1,-1,1)) !WzAg
     &      -s1c(psiAZg(1,-1,2,l),wa(1,1),.true.,-1,psiWP(1,-1,1))     !WAgz
     &      -s1c(psiAg(1,-1,2,l),atau(1,1),.true.,-1,psiWP(1,-1,1))    !WzgA
     &      -s1c(psiAZ(1,-1,2),wa(1,1),.true.,-1,psiWPg(1,-1,1,l))     !WgAz
     &      -s1c(psiA(1,-1,2),atau(1,1),.true.,-1,psiWPg(1,-1,1,l))    !WgzA
     &      -s1c(psiAZ(1,-1,2),wa(1,1),.true.,-1,braketgWP(1,-1,1,l))  !gWAz
     &      -s1c(psiA(1,-1,2),atau(1,1),.true.,-1,braketgWP(1,-1,1,l)) !gWzA
     D     )*clr(ifl(2),1,-1)*clr(ifl(2),1,-1)


         mat(1,-1,l,1) = (m1(1,l) + m2(1,l))*clr(3,3,-1) ! ud -> WZA

      enddo

c Calculation of processes qq -> V V g -> W Z A g -------------

c qq -> W A/Z (-> W Z/A) g

      do l=1,2
c W A
         call curr6(-1,braketgA(1,-1,2,l),pga(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,1,l) )         ! Wlt A g
         call curr6(-1,psiAg(1,-1,2,l),pag(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )         ! Wlt g A
         call curr6(-1,psiA(1,-1,2),fqa(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,3,l) )  ! g Wlt A
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m1(1,l) = -dotcc(wpwz(0,1),jwp(0,-1,4,l))
         m1(2,l) = zm2i(3)*dotrc(pwpz(0),wpwz(0,1))
     &              *dotrc(pwpz(0),jwp(0,-1,4,l))

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiWP(1,-1,1),fqp(0,1),jwp(0,-1,1,l) )         ! W Alt g
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiWPg(1,-1,1,l),pWPg(0,1),jwp(0,-1,2,l) )     ! W g Alt
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgWP(1,-1,1,l),pgWP(0,1),jwp(0,-1,3,l) )  ! g W Alt
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m1(1,l) = m1(1,l) -dotcc(aza(0,1),jwp(0,-1,4,l))

c W Z
ccc same current as before
c         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
c     &              psiWP(1,-1,1),fqp(0,1),jwp(0,-1,1,l) )         ! W Zlt g
c         call curr6(-1,psi(1,-1,2),p(0,2),
c     &              psiWPg(1,-1,1,l),pWPg(0,1),jwp(0,-1,2,l) )     ! W g Zlt
c         call curr6(-1,psi(1,-1,2),p(0,2),
c     &              braketgWP(1,-1,1,l),pgWP(0,1),jwp(0,-1,3,l) )  ! g W Zlt
c         do mu = 0,5
c           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
c         enddo
         m2(1,l) = -dotcc(zaz(0,1),jwp(0,-1,4,l))

         call curr6(-1,braketgZ(1,-1,2,l),pgZ(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,1,l) )         ! Wlt Z g
         call curr6(-1,psiZg(1,-1,2,l),pZg(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )         ! Wlt g Z
         call curr6(-1,psiZ(1,-1,2),fqz(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,3,l) )  ! g Wlt Z
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m2(1,l) = m2(1,l) -dotcc(wpwa(0,1),jwp(0,-1,4,l))
         m2(2,l) = zm2i(3)*dotrc(pwpa(0),wpwa(0,1))
     &              *dotrc(pwpa(0),jwp(0,-1,4,l))

c W z
         call curr6(-1,braketgAZ(1,-1,2,l),pgZ(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,1,l) )         ! Wlt z g
         call curr6(-1,psiAZg(1,-1,2,l),pZg(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )         ! Wlt g z
         call curr6(-1,psiAZ(1,-1,2),fqz(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,3,l) )  ! g Wlt z
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m1(1,l) = m1(1,l) -dotcc(wpwa(0,1),jwp(0,-1,4,l))
         m1(2,l) = m1(2,l) + zm2i(3)*dotrc(pwpa(0),wpwa(0,1))
     &              *dotrc(pwpa(0),jwp(0,-1,4,l))

c A W
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiA(1,-1,1),fqa(0,1),jwp(0,-1,1,l) )       ! A Wlt g
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiAg(1,-1,1,l),pAg(0,1),jwp(0,-1,2,l) )    ! A g Wlt 
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgA(1,-1,1,l),pgA(0,1),jwp(0,-1,3,l) ) ! g A Wlt
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m3(1,l) = -dotcc(wpwz(0,1),jwp(0,-1,4,l))
         m3(2,l) = zm2i(3)*dotrc(pwpz(0),wpwz(0,1))
     &              *dotrc(pwpz(0),jwp(0,-1,4,l))

         call curr6(-1,braketgWP(1,-1,2,l),pgWP(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,1,l) )         ! Alt W g
         call curr6(-1,psiWPg(1,-1,2,l),pWPg(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )         ! Alt g W
         call curr6(-1,psiWP(1,-1,2),fqp(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,3,l) )  ! g Alt W
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m3(1,l) = m3(1,l) -dotcc(aza(0,1),jwp(0,-1,4,l))

c Z W
ccc same current as before
c         call curr6(-1,braketgWP(1,-1,2,l),pgWP(0,2),
c     &              psi(1,-1,1),p(0,1),jwp(0,-1,1,l) )         ! Zlt W g
c         call curr6(-1,psiWPg(1,-1,2,l),pWPg(0,2),
c     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l) )         ! Zlt g W
c         call curr6(-1,psiWP(1,-1,2),fqp(0,2),
c     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,3,l) )  ! g Zlt W
c         do mu = 0,5
c           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
c         enddo
         m4(1,l) = -dotcc(zaz(0,1),jwp(0,-1,4,l))

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiZ(1,-1,1),fqz(0,1),jwp(0,-1,1,l) )       ! Z Wlt g
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiZg(1,-1,1,l),pZg(0,1),jwp(0,-1,2,l) )    ! Z g Wlt 
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgZ(1,-1,1,l),pgZ(0,1),jwp(0,-1,3,l) ) ! g Z Wlt
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m4(1,l) = m4(1,l) -dotcc(wpwa(0,1),jwp(0,-1,4,l))
         m4(2,l) = zm2i(3)*dotrc(pwpa(0),wpwa(0,1))
     &              *dotrc(pwpa(0),jwp(0,-1,4,l))

c z W
         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psiAZ(1,-1,1),fqz(0,1),jwp(0,-1,1,l) )       ! z Wlt g
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              psiAZg(1,-1,1,l),pZg(0,1),jwp(0,-1,2,l) )    ! z g Wlt 
         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketgAZ(1,-1,1,l),pgZ(0,1),jwp(0,-1,3,l) ) ! g z Wlt
         do mu = 0,5
           jwp(mu,-1,4,l)=jwp(mu,-1,1,l)+jwp(mu,-1,2,l)+jwp(mu,-1,3,l)
         enddo
         m3(1,l) = m3(1,l) -dotcc(wpwa(0,1),jwp(0,-1,4,l))
         m3(2,l) = m3(2,l) + zm2i(3)*dotrc(pwpa(0),wpwa(0,1))
     &              *dotrc(pwpa(0),jwp(0,-1,4,l))

ccc

         mat(1,-1,l,2) = ( (m1(1,l) + m1(2,l)) * clr(ifl(2),1,-1)
     &                    +(m2(1,l) + m2(2,l)) * clr(ifl(2),2,-1)
     &                    +(m3(1,l) + m3(2,l)) * clr(ifl(1),1,-1)
     &                    +(m4(1,l) + m4(2,l)) * clr(ifl(1),2,-1)
     &                  )*clr(3,3,-1)

      enddo          ! gluon polarization


c Calculation of processes qq -> A/Z g -> W W A g ---------------------

      do l =1,2

         call curr6(-1,psi(1,-1,2),p(0,2),
     &              braketg(1,-1,1,l),pq(0,1),jwp(0,-1,1,l))

         call curr6(-1,braketg(1,-1,2,l),pq(0,2),
     &              psi(1,-1,1),p(0,1),jwp(0,-1,2,l))

         do mu = 0,5
            jwp(mu,-1,4,l) = jwp(mu,-1,1,l) + jwp(mu,-1,2,l)
         enddo

c contract with atowwa/ztowwa tensor

         mw = -dotcc(zwwa(0,1),jwp(0,-1,4,l))
         mat(1,-1,l,3)  = mw * clr(3,3,-1)

      enddo              ! gluon polarization

c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph

      if (ldebugm) then
         
c         do isig = -1,1,2
            do i = 1,3
                  
                  m1kb(1,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(1,isig,i) = 0d0 ! kartesian basis l = 2
                  
                  do lh = -1,1,2
                     madgr(i,isig,lh,1) = 0d0
                  enddo
                  
            enddo
c         enddo
         
c         do isig = -1,1,2
            do lh = -1,1,2
                  do i =1,3
                     
                     madgr(i,isig,lh,1) = tampr(i,isig,lh,1)
                     
                  enddo
            enddo
c         enddo
         
c         do isig = -1,1,2
            do i=1,3
                  
                  m1kb(1,isig,i) = (madgr(i,isig,-1,1) - 
     &                 madgr(i,isig,1,1))/sqrt(2d0)       
                  m2kb(1,isig,i) = -im*(madgr(i,isig,-1,1) + 
     &                 madgr(i,isig,1,1))/sqrt(2d0)
                  
            enddo
c         enddo
         
c print results from checks:

         do i=1,3
c               do isig = -1,-1,2
                  
                  print*,"mat(",1,",",isig,",1,",i,") = ",mat(1,isig,1,i)
                  temp = dreal(mat(1,isig,1,i))**2+dimag(mat(1,isig,1,i))**2
                  print*,"m1kb(",1,",",isig,",",i,") = ",m1kb(1,isig,i)
                  temp1 = dreal(m1kb(1,isig,i))**2+dimag(m1kb(1,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
                  
                  print*,"mat(",1,",",isig,",2,",i,") = ",mat(1,isig,2,i)
                  temp = dreal(mat(1,isig,2,i))**2+dimag(mat(1,isig,2,i))**2
                  print*,"m2kb(",1,",",isig,",",i,") = ",m2kb(1,isig,i)
                  temp1 = dreal(m2kb(1,isig,i))**2+dimag(m2kb(1,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
c               enddo
         enddo
         

      endif ! end of debugging block

c------------------------------------------------------

c sum the graphs, square them and map them onto udwzaj 

         res(1) = 0
c         do isig = -1,1,2
            do l = 1,2
               mm(1,isig,l) = 0
               do i = 1,3
                  mm(1,isig,l) = mm(1,isig,l) +  mat(1,isig,l,i)
               enddo

               res(1)= res(1) + dreal(mm(1,isig,l))**2
     &                        + dimag(mm(1,isig,l))**2

            enddo
c         enddo
         res(1) = res(1)*4d0*fpials(2)  ! C_2*3 is the color factor

      udwzaj = res(1)
      
      return
      end


