c****************subroutine QQWWJ *************************************
c
c       Michael Rauch, <rauch@particle.uni-karlsruhe.de>
c       Last modified: December 2012
C
C  QQGWW_curr calculates the current  
C     (q1 q2    ---->   W+ W- g)_mu ,   
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c       fpials is attached only in the end of the code
c
c       k is the process ID (1:uuwwj,2:ddwwj)
c       isig is the helicity of partons 1 and 2
c
c*********************************************************************

      subroutine qqgWW_curr(final,N,p1,p2,p3,fsign1,fsign2,gsign)

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

      double complex final(0:3,-1:1,2)
      double precision p1(0:3), p2(0:3), p3(0:3)
      integer fsign1, fsign2, gsign

      double precision  pbar(0:3,n_p+n_v), qbar(0:3)
      double precision  p(0:3,n_p+n_v)

      integer  sign(n_p+n_v), mu, i, isig, N
      integer fsign(n_p+n_v)


      double complex psi(2,-1:1,2),psiv(2,-1:1,4,4),psivv(2,-1:1,4,4,4)
      double complex gvv12(1:6,-1:1,3:4)
      double precision pmomv(0:4,4,4),pmomvv(0:4,4,4,4)
      double precision res(2)
   
      real*8 dotrr,theta
      double precision contract_grjcj
      double complex im
      double complex dotcc, s1c, dotrc
      double complex zero 
      external dotcc, s1c, dotrc,dotrr,theta, contract_grjcj

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


      do mu=0,3
        pbar(mu,1) = p1(mu)
        pbar(mu,2) = p2(mu)
        qbar(mu)   = p3(mu)
      enddo
      sign(1) = fsign1
      sign(2) = fsign2

c -------------- BEGIN CODE --------------------------------

c initialize & precompute stuff needed below:

c fix strong coupling gs**2 for the two quarks:

      fpials(2) = fpi*als(1,N)
      fpials(3) = fpi*als(2,N)
c      fpials(2) = fpi*alfas

c---------------------------- BEGIN CODE -----------------------

c define flavors of external quarks for the 2 NC subprocesses

      if (linit) then
         linit = .false.
         zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))
      endif

c define the internal momenta
 
      do mu = 0,3
         do i = 1,2
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c----------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

!
!  First construct the fermion line which has everything attached
!
!  Add particles to psi
!
!  in psiv() etc: a=1, z=2, w+=3, w-=4
!  order in arguments is as in math. expression <bra|V...   ...V|ket>
!
        do isig=-1,1,2
! first V to <bra| or |ket>
          call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &      qww(0,N),aww(0,N),psiv(1,isig,2,1),pmomv(0,2,1))
          call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &      qww(0,N),zww(0,N),psiv(1,isig,2,2),pmomv(0,2,2))
          call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &      qww(0,N),aww(0,N),psiv(1,isig,1,1),pmomv(0,1,1))
          call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &      qww(0,N),zww(0,N),psiv(1,isig,1,2),pmomv(0,1,2))
          if (isig.eq.-1) then ! W only to left
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qp(0,N),wp(1,N),psiv(1,isig,2,3),pmomv(0,2,3))
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qm(0,N),wm(1,N),psiv(1,isig,2,4),pmomv(0,2,4))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qp(0,N),wp(1,N),psiv(1,isig,1,3),pmomv(0,1,3))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qm(0,N),wm(1,N),psiv(1,isig,1,4),pmomv(0,1,4))

! second V to <bra|V or V|ket> ! only W's
            call bra2c(psiv(1,isig,2,4),.false.,pmomv(0,2,4),isig,
     &        qp(0,N),wp(1,N),psivv(1,isig,2,4,3),pmomvv(0,2,4,3))
            call bra2c(psiv(1,isig,2,3),.false.,pmomv(0,2,3),isig,
     &        qm(0,N),wm(1,N),psivv(1,isig,2,3,4),pmomvv(0,2,3,4))
            call ket2c(psiv(1,isig,1,4),.false.,pmomv(0,1,4),isig,
     &        qp(0,N),wp(1,N),psivv(1,isig,1,3,4),pmomvv(0,1,3,4))
            call ket2c(psiv(1,isig,1,3),.false.,pmomv(0,1,3),isig,
     &        qm(0,N),wm(1,N),psivv(1,isig,1,4,3),pmomvv(0,1,4,3))
          endif

        enddo

c----------------------------------------------------------------------

! initialize g*12 
        do isig=-1,1,2
        do i=3,4
        do mu=1,4
          gvv12(mu,isig,i) = 0d0
        enddo
          gvv12(5,isig,i) = dcmplx(
     &      p(0,1)-p(0,2)-qp(0,N)-qm(0,N),
     &      p(3,1)-p(3,2)-qp(3,N)-qm(3,N))
          gvv12(6,isig,i) = dcmplx(
     &      p(1,1)-p(1,2)-qp(1,N)-qm(1,N),
     &      p(2,1)-p(2,2)-qp(2,N)-qm(2,N))
        enddo
        enddo

! add currents 
c both W
cc u quark
ccc W
        call curr6addm(-1,psivv(1,-1,2,4,3),pmomvv(0,2,4,3),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,4),pmomv(0,2,4),psiv(1,-1,1,3),pmomv(0,1,3),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),p(0,2),psivv(1,-1,1,4,3),pmomvv(0,1,4,3),
     &    gvv12(1,-1,3),clr(3,3,-1)**2,clr(3,3,+1)**2)
ccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
        call curr6addm(+1,psi(1,-1,2),p(0,2),psiv(1,-1,1,1),pmomv(0,1,1),
     &    gvv12(1,-1,3),-clr(3,1,-1),-clr(3,1,+1))
ccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
        call curr6addm(+1,psi(1,-1,2),p(0,2),psiv(1,-1,1,2),pmomv(0,1,2),
     &    gvv12(1,-1,3),-clr(3,2,-1),-clr(3,2,+1))
cc d quark
ccc W
        call curr6addm(-1,psivv(1,-1,2,3,4),pmomvv(0,2,3,4),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psiv(1,-1,2,3),pmomv(0,2,3),psiv(1,-1,1,4),pmomv(0,1,4),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
        call curr6addm(-1,psi(1,-1,2),p(0,2),psivv(1,-1,1,3,4),pmomvv(0,1,3,4),
     &    gvv12(1,-1,4),clr(3,3,-1)**2,clr(3,3,+1)**2)
ccc A
        call curr6addm(+1,psiv(1,-1,2,1),pmomv(0,2,1),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
        call curr6addm(+1,psi(1,-1,2),p(0,2),psiv(1,-1,1,1),pmomv(0,1,1),
     &    gvv12(1,-1,4),-clr(4,1,-1),-clr(4,1,+1))
ccc Z
        call curr6addm(+1,psiv(1,-1,2,2),pmomv(0,2,2),psi(1,-1,1),p(0,1),
     &    gvv12(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
        call curr6addm(+1,psi(1,-1,2),p(0,2),psiv(1,-1,1,2),pmomv(0,1,2),
     &    gvv12(1,-1,4),-clr(4,2,-1),-clr(4,2,+1))
      
c----------------------------------------------------------------------

c store current in output variable
c (0,0,i) contains integrated cross section

      do i=1,2
        final(0,0,i) = 0d0
        do isig=-1,1,2
          do mu=0,3
            final(mu,isig,i) = gvv12(mu+1,isig,i+2)*sqrt(fpials(2))
          enddo
          final(0,0,i) = final(0,0,i) - contract_grjcj(final(0,isig,i),final(0,isig,i))
        enddo
        final(0,0,i) = final(0,0,i)*4d0 ! colour factor
      enddo

      if (ldebug) then
        do mu=0,3
          pbar(mu,1) = p1(mu)
          pbar(mu,2) = p2(mu)
          qbar(mu)   = p3(mu)
        enddo
        fsign(1) = fsign1
        fsign(2) = fsign2

c compare with standard LO ME
        call qqVVj(pbar,fsign,qbar,gsign, 0, N, res(1), res(2), 1)

        do i=1,2
          write(*,*) "ME comparison: ", N, i, dble(final(0,0,i)), res(i), dble(final(0,0,i))/res(i)
        enddo
        read(*,*)
      endif

      return
      end

c======================================================================

