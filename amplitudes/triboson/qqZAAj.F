c****************subroutine QQZZAJ *************************************
c
c	Giuseppe Bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c	Last modified: November 2010
C
C  QQZAA calculates the matrix elements**2 for 
C        q1 q1bar    ---->   Z A A g,   
C         Z ---> f1-bar f1
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
c      return uuzaa = |M_real|^2   etc.
c
c       fpials is attached only in the end of the code
c
c       l is the gluon polarization in the kartesian basis (l=1,2)
c       k is the process ID (1:uuzaaj,2:ddzaaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqZAAj(bos,pbar,sign,qbar,gsign,uuzaaj,ddzaaj)

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
      double precision  uuzaaj, ddzaaj, res(2)
      double precision  p(0:3,2+n_v)

      integer  sign(2+n_v), gsign, mu, i, k, isig
      integer  ifl(2), lh, l, bos

      double precision fqze(0:4,2), fqa1(0:4,2), fqa2(0:4,2)
      double precision pq(0:4,2)
      double precision pgZe(0:4,2), pZeg(0:4,2)
      double precision pgA1(0:4,2),pA1g(0:4,2)
      double precision pgA2(0:4,2),pA2g(0:4,2)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psize(2,-1:1,2), psiae(2,-1:1,2)
      double complex psia1(2,-1:1,2), psia2(2,-1:1,2)
      double complex braketg(2,-1:1,2,2)
      double complex braketgZe(2,-1:1,2,2), braketgAe(2,-1:1,2,2)
      double complex braketgA1(2,-1:1,2,2), braketgA2(2,-1:1,2,2)
      double complex psiZeg(2,-1:1,2,2), psiAeg(2,-1:1,2,2)
      double complex psiA1g(2,-1:1,2,2), psiA2g(2,-1:1,2,2)
      double complex ja1(0:5,-1:1,4,2),
     1               ja2(0:5,-1:1,4,2)
      double complex mat(2,-1:1,2,5), madgr(5,-1:1,-1:1,2)
      double complex matvbf(2,-1:1,2,3)

      double complex m1kb(2,-1:1,5), m2kb(2,-1:1,5), temp, temp1
      double complex mm(2,-1:1,2), zero
      double complex ma, mz, m1, m2
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

      logical ldebugm, linit
      data linit /.true./
      save linit
      parameter(ldebugm = .false.)

c---------------------------- BEGIN CODE -----------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:
      fpials(2) = fpi*als(1,1)
      fpials(3) = fpi*als(2,1)
c      fpials(2) = fpi*alfas

c define flavors of external quarks for the 2 NC subprocesses
      if (linit) then
         linit = .false.
         k = 1                  ! uuzaaj: u = 3
         ifl(k) = 3
         k = 2                  ! ddzaaj: d = 4
         ifl(k) = 4
      endif

      do k = 1,2
         do isig = -1,1,2
            do i = 1,5
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
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa1(0,1),
     &              atau1(1,1), psia1(1,isig,2),fqa1(0,2))
         call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qa2(0,1),
     &              atau2(1,1), psia2(1,isig,2),fqa2(0,2))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,1),
     &              ze(1,1), psize(1,isig,1),fqze(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa1(0,1),
     &              atau1(1,1), psia1(1,isig,1),fqa1(0,1))
         call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qa2(0,1),
     &              atau2(1,1), psia2(1,isig,1),fqa2(0,1))
         if (bos.eq.211) then ! only for bos=211 (e-e+aa)
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,qze(0,1),
     &           ae(1,1),psiae(1,isig,2),fqze(0,2))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,qze(0,1),
     &           ae(1,1), psiae(1,isig,1),fqze(0,1))
         endif
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
            call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &           qze(0,1),ze(1,1),braketgZe(1,isig,1,l),pgZe(0,1))
            call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &           qze(0,1),ze(1,1),braketgZe(1,isig,2,l),pgZe(0,2))
            call ket2r(psiZe(1,isig,1),.false.,fqze(0,1),isig,q,eps(0,l),
     &           psiZeg(1,isig,1,l),pZeg(0,1))
            call bra2r(psiZe(1,isig,2),.false.,fqze(0,2),isig,q,eps(0,l),
     &           psiZeg(1,isig,2,l),pZeg(0,2))
            if (bos.eq.211) then ! only for bos=211 (e-e+aa)
               call ket2c(braketg(1,isig,1,l),.false.,pq(0,1),isig,
     &              qze(0,1),ae(1,1),braketgAe(1,isig,1,l),pgZe(0,1))
               call bra2c(braketg(1,isig,2,l),.false.,pq(0,2),isig,
     &              qze(0,1),ae(1,1),braketgAe(1,isig,2,l),pgZe(0,2))
               call ket2r(psiAe(1,isig,1),.false.,fqze(0,1),isig,q,eps(0,l),
     &              psiAeg(1,isig,1,l),pZeg(0,1))
               call bra2r(psiAe(1,isig,2),.false.,fqze(0,2),isig,q,eps(0,l),
     &              psiAeg(1,isig,2,l),pZeg(0,2))
            endif
         enddo                   !fermion helicity
      enddo                      !gluon polarization

c ********************************************************
      if (bos.eq.211) then
c ********************************************************

c qqbar -> Z A A g
      do l=1,2
         do isig = -1,1,2
           m1=
     1      -s1c(psiaeg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     1      -s1c(braketgae(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     1      -s1c(psiae(1,isig,2),atau1(1,1),.true.,isig,psia2g(1,isig,1,l))
     1      -s1c(psiae(1,isig,2),atau1(1,1),.true.,isig,braketga2(1,isig,1,l))
     2      -s1c(psiaeg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     2      -s1c(braketgae(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     2      -s1c(psiae(1,isig,2),atau2(1,1),.true.,isig,psia1g(1,isig,1,l))
     2      -s1c(psiae(1,isig,2),atau2(1,1),.true.,isig,braketga1(1,isig,1,l))
     3      -s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psiae(1,isig,1))
     3      -s1c(braketga2(1,isig,2,l),atau1(1,1),.true.,isig,psiae(1,isig,1))
     3      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,psiaeg(1,isig,1,l))
     3      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketgae(1,isig,1,l))
     4      -s1c(psia2g(1,isig,2,l),ae(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(braketga2(1,isig,2,l),ae(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(psia2(1,isig,2),ae(1,1),.true.,isig,psia1g(1,isig,1,l))
     4      -s1c(psia2(1,isig,2),ae(1,1),.true.,isig,braketga1(1,isig,1,l))
     5      -s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(braketga1(1,isig,2,l),atau2(1,1),.true.,isig,psiae(1,isig,1))
     5      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,psiaeg(1,isig,1,l))
     5      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketgae(1,isig,1,l))
     6      -s1c(psia1g(1,isig,2,l),ae(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(braketga1(1,isig,2,l),ae(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(psia1(1,isig,2),ae(1,1),.true.,isig,psia2g(1,isig,1,l))
     6      -s1c(psia1(1,isig,2),ae(1,1),.true.,isig,braketga2(1,isig,1,l))

           m2=
     1      -s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(braketga2(1,isig,2,l),atau1(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,psizeg(1,isig,1,l))
     1      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketgze(1,isig,1,l))
     2      -s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psize(1,isig,1))
     2      -s1c(braketga1(1,isig,2,l),atau2(1,1),.true.,isig,psize(1,isig,1))
     2      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,psizeg(1,isig,1,l))
     2      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketgze(1,isig,1,l))
     3      -s1c(psizeg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     3      -s1c(braketgze(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     3      -s1c(psize(1,isig,2),atau1(1,1),.true.,isig,psia2g(1,isig,1,l))
     3      -s1c(psize(1,isig,2),atau1(1,1),.true.,isig,braketga2(1,isig,1,l))
     4      -s1c(psizeg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(braketgze(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(psize(1,isig,2),atau2(1,1),.true.,isig,psia1g(1,isig,1,l))
     4      -s1c(psize(1,isig,2),atau2(1,1),.true.,isig,braketga1(1,isig,1,l))
     5      -s1c(psia2g(1,isig,2,l),ze(1,1),.true.,isig,psia1(1,isig,1))
     5      -s1c(braketga2(1,isig,2,l),ze(1,1),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),ze(1,1),.true.,isig,psia1g(1,isig,1,l))
     5      -s1c(psia2(1,isig,2),ze(1,1),.true.,isig,braketga1(1,isig,1,l))
     6      -s1c(psia1g(1,isig,2,l),ze(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(braketga1(1,isig,2,l),ze(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(psia1(1,isig,2),ze(1,1),.true.,isig,psia2g(1,isig,1,l))
     6      -s1c(psia1(1,isig,2),ze(1,1),.true.,isig,braketga2(1,isig,1,l))
           
           mat(1,isig,l,1) =  m1 * clr(3,1,isig)**3               !uu->AAA
     1                      + m2 * clr(3,2,isig)*clr(3,1,isig)**2 !uu->AAZ

           mat(2,isig,l,1) =  m1 * clr(4,1,isig)**3               !dd->AAA
     1                      + m2 * clr(4,2,isig)*clr(4,1,isig)**2 !dd->AAZ
        enddo
      enddo

c qqbar -> A A/Z (-> Z A) g
      do l = 1,2
         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketga1(1,-1,1,l),pgA1(0,1), ja1(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psia1g(1,-1,1,l),pA1g(0,1), ja1(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psia1(1,-1,1),fqa1(0,1), ja1(0,-1,3,l) )

         call curr6(1,psi(1,-1,2),p(0,2),
     &           braketga2(1,-1,1,l),pgA2(0,1), ja2(0,-1,1,l) )
         call curr6(1,psi(1,-1,2),p(0,2),
     &           psia2g(1,-1,1,l),pA2g(0,1), ja2(0,-1,2,l) )
         call curr6(1,braketg(1,-1,2,l),pq(0,2),
     &           psia2(1,-1,1),fqa2(0,1), ja2(0,-1,3,l) )

         do isig = -1,1,2
            do mu = 0,5
               ja1(mu,isig,4,l) = ja1(mu,isig,1,l) + ja1(mu,isig,2,l)
     &                          + ja1(mu,isig,3,l)
               ja2(mu,isig,4,l) = ja2(mu,isig,1,l) + ja2(mu,isig,2,l)
     &                          + ja2(mu,isig,3,l)
            enddo

c contract with atoza/ztoza tensor
            m1 = -dotcc(aaze2(0,1),ja1(0,isig,4,l))
     1           -dotcc(aaze1(0,1),ja2(0,isig,4,l))
            m2 = -dotcc(zzae2(0,1),ja1(0,isig,4,l))
     1           -dotcc(zzae1(0,1),ja2(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,2) = (m1*clr(ifl(k),1,isig)
     1                         +  m2*clr(ifl(k),2,isig))
     2                         *  clr(ifl(k),1,isig)
            enddo
         enddo
      enddo

c qqbar -> A/Z(-> Z A) A g
      do l = 1,2
         call curr6(1,psia1(1,-1,2),fqa1(0,2),
     &           braketg(1,-1,1,l),pq(0,1), ja1(0,-1,1,l) )
         call curr6(1,psia1g(1,-1,2,l),pA1g(0,2),
     &           psi(1,-1,1),p(0,1), ja1(0,-1,2,l) )
         call curr6(1,braketga1(1,-1,2,l),pgA1(0,2),
     &           psi(1,-1,1),p(0,1), ja1(0,-1,3,l) )

         call curr6(1,psia2(1,-1,2),fqa2(0,2),
     &           braketg(1,-1,1,l),pq(0,1), ja2(0,-1,1,l) )
         call curr6(1,psia2g(1,-1,2,l),pA2g(0,2),
     &           psi(1,-1,1),p(0,1), ja2(0,-1,2,l) )
         call curr6(1,braketga2(1,-1,2,l),pgA2(0,2),
     &           psi(1,-1,1),p(0,1), ja2(0,-1,3,l) )

         do isig = -1,1,2
            do mu = 0,5
               ja1(mu,isig,4,l) = ja1(mu,isig,1,l) + ja1(mu,isig,2,l)
     &                          + ja1(mu,isig,3,l)
               ja2(mu,isig,4,l) = ja2(mu,isig,1,l) + ja2(mu,isig,2,l)
     &                          + ja2(mu,isig,3,l)
            enddo

c contract with atoza/ztoza tensor
            m1 = -dotcc(aaze2(0,1),ja1(0,isig,4,l))
     1           -dotcc(aaze1(0,1),ja2(0,isig,4,l))
            m2 = -dotcc(zzae2(0,1),ja1(0,isig,4,l))
     1           -dotcc(zzae1(0,1),ja2(0,isig,4,l))
            do k = 1,2
               mat(k,isig,l,3) = (m1*clr(ifl(k),1,isig)
     1                         +  m2*clr(ifl(k),2,isig))
     2                         * clr(ifl(k),1,isig)
            enddo
         enddo
      enddo

c qqbar -> A/Z g -> Z A A g ---------------------
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
            ma = -dotcc(azaa(0,1),jqqa(0,isig,l))
            mz = -dotcc(zzaa(0,1),jqqa(0,isig,l))
            do k = 1,2
               mat(k,isig,l,4)  = ma*clr(ifl(k),1,isig)
               mat(k,isig,l,5)  = mz*clr(ifl(k),2,isig)
            enddo
         enddo
      enddo

c ********************************************************
      elseif (bos.eq.511) then
c ********************************************************

c qqbar -> Z A A g
      do l=1,2
        do isig = -1,1,2
           m1=
     1      -s1c(psia2g(1,isig,2,l),atau1(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(braketga2(1,isig,2,l),atau1(1,1),.true.,isig,psize(1,isig,1))
     1      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,psizeg(1,isig,1,l))
     1      -s1c(psia2(1,isig,2),atau1(1,1),.true.,isig,braketgze(1,isig,1,l))
     2      -s1c(psia1g(1,isig,2,l),atau2(1,1),.true.,isig,psize(1,isig,1))
     2      -s1c(braketga1(1,isig,2,l),atau2(1,1),.true.,isig,psize(1,isig,1))
     2      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,psizeg(1,isig,1,l))
     2      -s1c(psia1(1,isig,2),atau2(1,1),.true.,isig,braketgze(1,isig,1,l))
     3      -s1c(psizeg(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     3      -s1c(braketgze(1,isig,2,l),atau1(1,1),.true.,isig,psia2(1,isig,1))
     3      -s1c(psize(1,isig,2),atau1(1,1),.true.,isig,psia2g(1,isig,1,l))
     3      -s1c(psize(1,isig,2),atau1(1,1),.true.,isig,braketga2(1,isig,1,l))
     4      -s1c(psizeg(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(braketgze(1,isig,2,l),atau2(1,1),.true.,isig,psia1(1,isig,1))
     4      -s1c(psize(1,isig,2),atau2(1,1),.true.,isig,psia1g(1,isig,1,l))
     4      -s1c(psize(1,isig,2),atau2(1,1),.true.,isig,braketga1(1,isig,1,l))
     5      -s1c(psia2g(1,isig,2,l),ze(1,1),.true.,isig,psia1(1,isig,1))
     5      -s1c(braketga2(1,isig,2,l),ze(1,1),.true.,isig,psia1(1,isig,1))
     5      -s1c(psia2(1,isig,2),ze(1,1),.true.,isig,psia1g(1,isig,1,l))
     5      -s1c(psia2(1,isig,2),ze(1,1),.true.,isig,braketga1(1,isig,1,l))
     6      -s1c(psia1g(1,isig,2,l),ze(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(braketga1(1,isig,2,l),ze(1,1),.true.,isig,psia2(1,isig,1))
     6      -s1c(psia1(1,isig,2),ze(1,1),.true.,isig,psia2g(1,isig,1,l))
     6      -s1c(psia1(1,isig,2),ze(1,1),.true.,isig,braketga2(1,isig,1,l))
           
           mat(1,isig,l,1) =  m1 * clr(3,2,isig)*clr(3,1,isig)**2 !uu->AAZ
           mat(2,isig,l,1) =  m1 * clr(4,2,isig)*clr(4,1,isig)**2 !dd->AAZ

        enddo
      enddo

c qqbar -> A/Z g -> Z A A g ---------------------
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
            ma = -dotcc(azaa(0,1),jqqa(0,isig,l))
            mz = -dotcc(zzaa(0,1),jqqa(0,isig,l))
            do k = 1,2
               mat(k,isig,l,4)  = ma*clr(ifl(k),1,isig)
               mat(k,isig,l,5)  = mz*clr(ifl(k),2,isig)
            enddo
         enddo
      enddo

      endif
c-----------------------------------------------------
c------------------------------------------------------
c comparison with madgraph
      if (ldebugm) then
         do isig = -1,1,2
            do k = 1,2
               do i =1,3
                  m1kb(k,isig,i) = 0d0 ! kartesian basis l = 1
                  m2kb(k,isig,i) = 0d0 ! kartesian basis l = 2
                  do lh = -1,1,2
                     madgr(i,isig,lh,k) = tampr(i,isig,lh,k)
                  enddo
               enddo
               do lh = 1,2
                  matvbf(k,isig,lh,1) = mat(k,isig,lh,1)
                  matvbf(k,isig,lh,2) = mat(k,isig,lh,2)+mat(k,isig,lh,3)
                  matvbf(k,isig,lh,3) = mat(k,isig,lh,4)+mat(k,isig,lh,5)
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

         do k=1,2
            do isig = -1,-1,2
               do i=1,3
                  
                  print*,"mat(",k,",",isig,",1,",i,") = ",matvbf(k,isig,1,i)
                  temp = dreal(matvbf(k,isig,1,i))**2+dimag(matvbf(k,isig,1,i))**2
                  print*,"m1kb(",k,",",isig,",",i,") = ",m1kb(k,isig,i)
                  temp1 = dreal(m1kb(k,isig,i))**2+dimag(m1kb(k,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
                  
                  print*,"mat(",k,",",isig,",2,",i,") = ",matvbf(k,isig,2,i)
                  temp = dreal(matvbf(k,isig,2,i))**2+dimag(matvbf(k,isig,2,i))**2
                  print*,"m2kb(",k,",",isig,",",i,") = ",m2kb(k,isig,i)
                  temp1 = dreal(m2kb(k,isig,i))**2+dimag(m2kb(k,isig,i))**2
                  print*,"temp/temp1 = ",temp/temp1
               enddo
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
               do i = 1,5
                  mm(k,isig,l) = mm(k,isig,l) +  mat(k,isig,l,i)
               enddo
               res(k)= res(k) + dreal(mm(k,isig,l))**2
     &                        + dimag(mm(k,isig,l))**2
            enddo
         enddo
         res(k) = res(k)*4d0*fpials(2)  ! C_2*3 is the color factor
      enddo
      uuzaaj = res(1)
      ddzaaj = res(2)
      return
      end
