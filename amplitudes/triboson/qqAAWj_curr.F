c****************subroutine QQAAWJ *************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Last modified: September 2010
C
C  QQVAAj_curr calculates the current  
C     (q1 q2    ---->   W+- A A)_mu ,   
C     W+- ---> f1-bar f2
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c       k is the process ID (1:uuwwaj,2:ddwwaj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqAAWj_curr(momin, final, N, bos, fsign)

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

      double complex final(0:3)
      double precision  pbar(0:3,n_p+n_v)
      double precision  momin(0:3,n_p+n_v)

      double precision  p(0:3,n_p+n_v),fq0(0:4,2),fqa1(0:4,2),fqa2(0:4,2)
      double precision  fq0a1(0:4,2),fq0a2(0:4,2),fqa1a2(0:4,2)
      double precision  fq0a1a2(0:4,2)

      integer  sign(n_p+n_v), mu, i, k, N
      integer fsign(n_p+n_v)
      integer  bos


      double precision pwa1(0:4), pwa2(0:4), pwa1a2(0:4)
      double precision eps(0:3,2) ! g in kartesian basis
      double complex psi(2,-1:1,2)


      double complex psiw0(2,-1:1, 2),psia1(2,-1:1, 2),
     1               psia2(2,-1:1, 2)
      double complex psiw0a1(2,-1:1, 2),psiw0a2(2,-1:1, 2),
     1               psia1a2(2,-1:1, 2)
      double complex psiw0a1a2(2,-1:1, 2)


      double complex jcurr(0:5,-1:1)
   
      real*8 dotrr,theta
      double complex im
      double complex dotcc, s1c, dotrc
      double complex zero 
      external dotcc, s1c, dotrc,dotrr,theta
      integer u(2),d(2)
      parameter (im = (0d0,1d0))
      parameter (zero = (0d0,0d0) )

c variable for q_mu q_nu terms in the propagators of the massive vector bosons

      double complex zm2i(2:3)
      save zm2i

c  helicity selection

      INTEGER h
      COMMON /hcount / h

      logical ldebug, ldebugm, linit, christophmods

      data linit /.true./
      save linit
      parameter(ldebugm = .false., ldebug=.false., christophmods=.false.)

c -------------- BEGIN ENTRY --------------------------------

      do i=1,n_p+n_v
        do mu=0,3
          pbar(mu,i) = momin(mu,i)
        enddo
        sign(i) = fsign(i)
      enddo
      goto 1
c ----
      ENTRY UGDAAWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,3)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,2)
      enddo
      sign(1) =  1
      sign(2) =  1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) = -1
      goto 1

c ----
      ENTRY DGUAAWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,3)
        pbar(mu,2) = momin(mu,1)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,2)
      enddo
      sign(1) = -1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) = -1
      goto 1

c ----
      ENTRY UDGAAWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,2)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,3)
      enddo
      sign(1) =  1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) =  1
      goto 1

1     continue

c---------------------------- BEGIN CODE -----------------------

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
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
      endif

c define the internal momenta
 
      do mu = 0,3
         do i = 1,n_p+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwa1(mu) = qw0(mu,N) + qa1(mu,N)
         pwa2(mu) = qw0(mu,N) + qa2(mu,N)
         pwa1a2(mu) = qw0(mu,N) + qa1(mu,N) + qa2(mu,N)
      enddo
      pwa1(4) = pwa1(0)**2-pwa1(1)**2-pwa1(2)**2-pwa1(3)**2
      pwa2(4) = pwa2(0)**2-pwa2(1)**2-pwa2(2)**2-pwa2(3)**2
      pwa1a2(4) = pwa1a2(0)**2-pwa1a2(1)**2-pwa1a2(2)**2-pwa1a2(3)**2

c----------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

c get other bra and ket vectors

c first boson to <bra| and |ket>
c W
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qw0(0,N),w0(1,N),
     &              psiw0(1,-1,2),fq0(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qw0(0,N),w0(1,N),
     &              psiw0(1,-1,1),fq0(0,1))
c A1
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,2),fqa1(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa1(0,N),
     &     atau1(1,N),psia1(1,-1,1),fqa1(0,1))
c A2
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,2),fqa2(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,qa2(0,N),
     &     atau2(1,N),psia2(1,-1,1),fqa2(0,1))
      
c second boson to <bra| and |ket>
c W A1
cc WA1 via current
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,pwa1(0),wpwa1(0,N),
     &              psiw0a1(1,-1,2),fq0a1(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,pwa1(0),wpwa1(0,N),
     &              psiw0a1(1,-1,1),fq0a1(0,1))
cc W A1
      call bra2cadd(psiw0(1,-1,2),.false.,fq0(0,2),-1,qa1(0,N),
     &     atau1(1,N),psiw0a1(1,-1,2),fq0a1(0,2),clr(u(k),1,-1))
      call ket2cadd(psiw0(1,-1,1),.false.,fq0(0,1),-1,qa1(0,N),
     &     atau1(1,N),psiw0a1(1,-1,1),fq0a1(0,1),clr(d(k),1,-1))
cc A1 W
      call bra2cadd(psia1(1,-1,2),.false.,fqa1(0,2),-1,qw0(0,N),
     &     w0(1,N),psiw0a1(1,-1,2),fq0a1(0,2),clr(d(k),1,-1))
      call ket2cadd(psia1(1,-1,1),.false.,fqa1(0,1),-1,qw0(0,N),
     &     w0(1,N),psiw0a1(1,-1,1),fq0a1(0,1),clr(u(k),1,-1))
c W A2
cc WA2 via current
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,pwa2(0),wpwa2(0,N),
     &              psiw0a2(1,-1,2),fq0a2(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,pwa2(0),wpwa2(0,N),
     &              psiw0a2(1,-1,1),fq0a2(0,1))
cc W A2
      call bra2cadd(psiw0(1,-1,2),.false.,fq0(0,2),-1,qa2(0,N),
     &     atau2(1,N),psiw0a2(1,-1,2),fq0a2(0,2),clr(u(k),1,-1))
      call ket2cadd(psiw0(1,-1,1),.false.,fq0(0,1),-1,qa2(0,N),
     &     atau2(1,N),psiw0a2(1,-1,1),fq0a2(0,1),clr(d(k),1,-1))
cc A2 W
      call bra2cadd(psia2(1,-1,2),.false.,fqa2(0,2),-1,qw0(0,N),
     &     w0(1,N),psiw0a2(1,-1,2),fq0a2(0,2),clr(d(k),1,-1))
      call ket2cadd(psia2(1,-1,1),.false.,fqa2(0,1),-1,qw0(0,N),
     &     w0(1,N),psiw0a2(1,-1,1),fq0a2(0,1),clr(u(k),1,-1))
c A1 A2
cc A1 A2
      call bra2c(psia1(1,-1,2),.false.,fqa1(0,2),-1,qa2(0,N),
     &     atau2(1,N),psia1a2(1,-1,2),fqa1a2(0,2))
      call ket2c(psia1(1,-1,1),.false.,fqa1(0,1),-1,qa2(0,N),
     &     atau2(1,N),psia1a2(1,-1,1),fqa1a2(0,1))
cc A2 A1
      call bra2cadd(psia2(1,-1,2),.false.,fqa2(0,2),-1,qa1(0,N),
     &     atau1(1,N),psia1a2(1,-1,2),fqa1a2(0,2),1d0)
      call ket2cadd(psia2(1,-1,1),.false.,fqa2(0,1),-1,qa1(0,N),
     &     atau1(1,N),psia1a2(1,-1,1),fqa1a2(0,1),1d0)

c third boson to <bra| and |ket>
cc WA1A2 via current
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,pwa1a2(0),wpwaa(0,N),
     &              psiw0a1a2(1,-1,2),fq0a1a2(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,pwa1a2(0),wpwaa(0,N),
     &              psiw0a1a2(1,-1,1),fq0a1a2(0,1))
cc WA1 A2
      call bra2cadd(psiw0a1(1,-1,2),.false.,fq0a1(0,2),-1,qa2(0,N),
     &     atau2(1,N),psiw0a1a2(1,-1,2),fq0a1a2(0,2),clr(u(k),1,-1))
      call ket2cadd(psiw0a1(1,-1,1),.false.,fq0a1(0,1),-1,qa2(0,N),
     &     atau2(1,N),psiw0a1a2(1,-1,1),fq0a1a2(0,1),clr(d(k),1,-1))
cc A2 WA1
      call bra2cadd(psia2(1,-1,2),.false.,fqa2(0,2),-1,pwa1(0),
     &     wpwa1(0,N),psiw0a1a2(1,-1,2),fq0a1a2(0,2),clr(d(k),1,-1))
      call ket2cadd(psia2(1,-1,1),.false.,fqa2(0,1),-1,pwa1(0),
     &     wpwa1(0,N),psiw0a1a2(1,-1,1),fq0a1a2(0,1),clr(u(k),1,-1))
cc WA2 A1
      call bra2cadd(psiw0a2(1,-1,2),.false.,fq0a2(0,2),-1,qa1(0,N),
     &     atau1(1,N),psiw0a1a2(1,-1,2),fq0a1a2(0,2),clr(u(k),1,-1))
      call ket2cadd(psiw0a2(1,-1,1),.false.,fq0a2(0,1),-1,qa1(0,N),
     &     atau1(1,N),psiw0a1a2(1,-1,1),fq0a1a2(0,1),clr(d(k),1,-1))
cc A1 WA2
      call bra2cadd(psia1(1,-1,2),.false.,fqa1(0,2),-1,pwa2(0),
     &     wpwa2(0,N),psiw0a1a2(1,-1,2),fq0a1a2(0,2),clr(d(k),1,-1))
      call ket2cadd(psia1(1,-1,1),.false.,fqa1(0,1),-1,pwa2(0),
     &     wpwa2(0,N),psiw0a1a2(1,-1,1),fq0a1a2(0,1),clr(u(k),1,-1))
cc A1A2 W
      call bra2cadd(psia1a2(1,-1,2),.false.,fqa1a2(0,2),-1,qw0(0,N),
     &     w0(1,N),psiw0a1a2(1,-1,2),fq0a1a2(0,2),clr(d(k),1,-1)**2)
      call ket2cadd(psia1a2(1,-1,1),.false.,fqa1a2(0,1),-1,qw0(0,N),
     &     w0(1,N),psiw0a1a2(1,-1,1),fq0a1a2(0,1),clr(u(k),1,-1)**2)

c----------------------------------------------------------------------

c Build current
c (number of vector bosons attached)
c <0|3>
         call curr6(-1,psi(1,-1,2),p(0,2),
     &          psiw0a1a2(1,-1,1),fq0a1a2(0,1),jcurr(0,-1) )
c <1|2>
cc <W|A1A2>
         call curr6addm(-1,psiw0(1,-1,2),fq0(0,2),
     &          psia1a2(1,-1,1),fqa1a2(0,1),jcurr(0,-1),
     &            clr(u(k),1,-1)**2,clr(u(k),1,1)**2)
cc <A1|WA2>
         call curr6addm(-1,psia1(1,-1,2),fqa1(0,2),
     &          psiw0a2(1,-1,1),fq0a2(0,1),jcurr(0,-1),
     &            clr(d(k),1,-1),clr(d(k),1,1))
cc <A2|WA1>
         call curr6addm(-1,psia2(1,-1,2),fqa2(0,2),
     &          psiw0a1(1,-1,1),fq0a1(0,1),jcurr(0,-1),
     &            clr(d(k),1,-1),clr(d(k),1,1))
c <2|1>
cc <A1A2|W>
         call curr6addm(-1,psia1a2(1,-1,2),fqa1a2(0,2),
     &          psiw0(1,-1,1),fq0(0,1),jcurr(0,-1),
     &            clr(d(k),1,-1)**2,clr(d(k),1,1)**2)
cc <WA2|A1>
         call curr6addm(-1,psiw0a2(1,-1,2),fq0a2(0,2),
     &          psia1(1,-1,1),fqa1(0,1),jcurr(0,-1),
     &            clr(u(k),1,-1),clr(u(k),1,1))
cc <WA1|A2>
         call curr6addm(-1,psiw0a1(1,-1,2),fq0a1(0,2),
     &          psia2(1,-1,1),fqa2(0,1),jcurr(0,-1),
     &            clr(u(k),1,-1),clr(u(k),1,1))
c <0|3>
         call curr6add(-1,psiw0a1a2(1,-1,2),fq0a1a2(0,2),
     &          psi(1,-1,1),p(0,1),jcurr(0,-1) )
      
c----------------------------------------------------------------------

c store current in output variable
c add remaining coupling factors g_s and g_Wff

      do i=0,3
        final(i) = -jcurr(i,-1) * gg(1) * clr(1,3,-1)
      enddo

      return
      end

c======================================================================

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculate squared matrix element with arguments suitable for dipoles 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine UGDAAWj_msq(momin, final, N, bos, ps_number)

      implicit none

#include "VBFNLO/utilities/global.inc"

      integer ps_number

      double precision momin(0:3,n_p+n_v), final
      double precision pbar(0:3,n_p+n_v), qbar(0:3)
      double precision avgfac
      integer N, bos, sign(7), gsign
      integer mu, nlo

      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,3)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,2)
        qbar(mu) = pbar(mu,7)
      enddo
      sign(1) =  1
      sign(2) =  1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) = -1
      gsign = sign(7)
      avgfac = 1d0/4d0/3d0/8d0
      goto 1

c ----
      ENTRY DGUAAWj_msq(momin, final, N, bos, ps_number)
      do mu=0,3
        pbar(mu,1) = momin(mu,3)
        pbar(mu,2) = momin(mu,1)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,2)
        qbar(mu) = pbar(mu,7)
      enddo
      sign(1) = -1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) = -1
      gsign = sign(7)
      avgfac = 1d0/4d0/3d0/8d0
      goto 1

c ----
      ENTRY UDGAAWj_msq(momin, final, N, bos, ps_number)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,2)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
        pbar(mu,5) = momin(mu,6)
        pbar(mu,6) = momin(mu,7)
        pbar(mu,7) = momin(mu,3)
        qbar(mu) = pbar(mu,7)
      enddo
      sign(1) =  1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
      sign(5) =  1
      sign(6) =  1
      sign(7) =  1
      gsign = sign(7)
      avgfac = 1d0/4d0/3d0/3d0
      goto 1

1     continue
      nlo = 0
      call qqAAWj(pbar,sign, qbar, gsign, nlo, N, final, -bos,ps_number)
      final = final * avgfac

      return
      end

