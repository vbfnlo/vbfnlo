c****************subroutine QQZZJ *************************************
c
c       Michael Rauch, <michael.rauch@kit.edu>
c       Last modified: October 2014
C
C  QQGZAZA_curr calculates the current  
C     (q1 q2    ---->   Z/A Z/A g)_mu ,   
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

      subroutine qqgZAZA_curr(final,N,p1,p2,p3,fsign1,fsign2,gsign,bos)

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

      integer  sign(n_p+n_v), mu, i, k, isig, N
      integer fsign(n_p+n_v)
      integer  l, bos

      double complex psi(2,-1:1,2),psiv(2,-1:1,4,2,4),psivv(2,-1:1,4,3,4)
      double complex gvv12(1:6,-1:1,3:4)
      double precision pmomv(0:4,4,2,4),pmomvv(0:4,4,3,4)
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

c variables for ZZ,ZA,AA 
      integer typ, typl
      integer cpp(4), cpm(4)
      double complex cepsp(6,4), cepsm(6,4)
      logical lhaszvv


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

c define the internal momenta
 
      do mu = 0,3
         do i = 1,2
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
      enddo

c----------------------------------------------------------------------

c map polarisation vectors and indices onto the arrays, so we can use do-loops
c makes debugging way easier
c typ: 1:AA, 2: ZA, 3: AZ, 4: ZZ (so can just stop at right index for AA and ZA)
      if (linit) then
        zm2i(3) = 1d0/dcmplx(xm2(3),-xmg(3))

        cpp(1) = 1
        cpp(2) = 2
        cpp(3) = 1
        cpp(4) = 2
        cpm(1) = 1
        cpm(2) = 1
        cpm(3) = 2
        cpm(4) = 2
        if (bos.eq.22) then
          typ = 4
          lhaszvv = .true.
        else if (bos.eq.21) then
          typ = 2
          lhaszvv = .true.
        else if (bos.eq.11) then
          typ = 1
          lhaszvv = .false.
        else
          write(*,*) "Error in ggZZg: wrong bos ",bos
          stop
        endif
c        if (narrowwidth.eq.0) then
          typl = 1
c        else
c          typl = typ
c        endif
        linit = .false.
      endif
      do mu=1,6
        cepsp(mu,1) = atau1(mu,N)
        cepsp(mu,2) = wp(mu,N)
        cepsp(mu,3) = atau1(mu,N)
        cepsp(mu,4) = wp(mu,N)
        cepsm(mu,1) = atau2(mu,N)
        cepsm(mu,2) = atau2(mu,N)
        cepsm(mu,3) = wm(mu,N)
        cepsm(mu,4) = wm(mu,N)
      enddo
      if (bos.eq.11) then
        do mu=0,3
          qp(mu,N) = qa1(mu,N)
          qm(mu,N) = qa2(mu,N)
        enddo
      endif

c----------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

!
!  First construct the fermion line which has everything attached
!
!  Add particles to psi
!
!
        do isig=-1,1,2
! first V to <bra| or |ket>
          do k=typ,typl,-1
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qp(0,N),cepsp(1,k),psiv(1,isig,2,1,k),pmomv(0,2,1,k))
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qm(0,N),cepsm(1,k),psiv(1,isig,2,2,k),pmomv(0,2,2,k))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qp(0,N),cepsp(1,k),psiv(1,isig,1,1,k),pmomv(0,1,1,k))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qm(0,N),cepsm(1,k),psiv(1,isig,1,2,k),pmomv(0,1,2,k))
! second V to <bra|V or V|ket> 
            call bra2c(psiv(1,isig,2,2,k),.false.,pmomv(0,2,2,k),isig,
     &        qp(0,N),cepsp(1,k),psivv(1,isig,2,2,k),pmomvv(0,2,2,k))
            call bra2c(psiv(1,isig,2,1,k),.false.,pmomv(0,2,1,k),isig,
     &        qm(0,N),cepsm(1,k),psivv(1,isig,2,1,k),pmomvv(0,2,1,k))
            call ket2c(psiv(1,isig,1,2,k),.false.,pmomv(0,1,2,k),isig,
     &        qp(0,N),cepsp(1,k),psivv(1,isig,1,1,k),pmomvv(0,1,1,k))
            call ket2c(psiv(1,isig,1,1,k),.false.,pmomv(0,1,1,k),isig,
     &        qm(0,N),cepsm(1,k),psivv(1,isig,1,2,k),pmomvv(0,1,2,k))
          enddo
! joint VV tensor to <bra| or |ket> - last index in psivv/pmomvv denotes boson index
          if (lhaszvv) then
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qaa(0,N),wz(1,N),psivv(1,isig,2,3,2),pmomvv(0,2,3,2))
            call bra2c(psi(1,isig,2),.true.,p(0,2),isig,
     &        qaa(0,N),wa(1,N),psivv(1,isig,2,3,1),pmomvv(0,2,3,1))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qaa(0,N),wz(1,N),psivv(1,isig,1,3,2),pmomvv(0,1,3,2))
            call ket2c(psi(1,isig,1),.true.,p(0,1),isig,
     &        qaa(0,N),wa(1,N),psivv(1,isig,1,3,1),pmomvv(0,1,3,1))
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
c both bosons
        do i=3,4 ! fermion flavour
          do k=typ,typl,-1
            do l=1,2 ! 2 permutations of attaching V_e V_mu
cc VV|g|
              call curr6addm(+1,psivv(1,-1,2,l,k),pmomvv(0,2,l,k),
     &          psi(1,-1,1),p(0,1),gvv12(1,-1,i),
     &          clr(i,cpp(k),-1)*clr(i,cpm(k),-1),
     &          clr(i,cpp(k),+1)*clr(i,cpm(k),+1))
cc |g|VV
              call curr6addm(+1,psi(1,-1,2),p(0,2),
     &          psivv(1,-1,1,l,k),pmomvv(0,1,l,k),gvv12(1,-1,i),
     &          clr(i,cpp(k),-1)*clr(i,cpm(k),-1),
     &          clr(i,cpp(k),+1)*clr(i,cpm(k),+1))
            enddo
cc V_e|g|V_mu
            call curr6addm(+1,psiv(1,-1,2,1,k),pmomv(0,2,1,k),
     &        psiv(1,-1,1,2,k),pmomv(0,1,2,k),gvv12(1,-1,i),
     &        clr(i,cpp(k),-1)*clr(i,cpm(k),-1),
     &        clr(i,cpp(k),+1)*clr(i,cpm(k),+1))
cc V_mu|g|V_e
            call curr6addm(+1,psiv(1,-1,2,2,k),pmomv(0,2,2,k),
     &        psiv(1,-1,1,1,k),pmomv(0,1,1,k),gvv12(1,-1,i),
     &        clr(i,cpp(k),-1)*clr(i,cpm(k),-1),
     &        clr(i,cpp(k),+1)*clr(i,cpm(k),+1))
          enddo
          if (lhaszvv) then
            do k=1,2 ! typ
cc (VV)|g|
              call curr6addm(+1,psivv(1,-1,2,3,k),pmomvv(0,2,3,k),
     &          psi(1,-1,1),p(0,1),gvv12(1,-1,i),
     &          -clr(i,k,-1),
     &          -clr(i,k,+1))
cc |g|(VV)
              call curr6addm(+1,psi(1,-1,2),p(0,2),
     &          psivv(1,-1,1,3,k),pmomvv(0,1,3,k),gvv12(1,-1,i),
     &          -clr(i,k,-1),
     &          -clr(i,k,+1))
            enddo
          endif
        enddo
      
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
        if (bos.eq.22) then
          call qqZZj(pbar,fsign,qbar,gsign, 0, N, res(1), res(2), 1)
        else if (bos.eq.21) then
          call qqZAj(pbar,fsign,qbar,gsign, 0, N, res(1), res(2), 1)
        else if (bos.eq.11) then
          call qqAAj(pbar,fsign,qbar,gsign, 0, N, res(1), res(2), 1)
        endif

        do i=1,2
          write(*,*) "ME comparison: ", N, i, dble(final(0,0,i)), res(i), dble(final(0,0,i))/res(i)
        enddo
        read(*,*)
      endif

      return
      end

c======================================================================

