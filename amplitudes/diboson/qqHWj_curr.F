c****************subroutine QQHWJ *************************************
c
c      adapted from Michael Rauch QQAAWJ_curr
c       Last modified: September 2010
C
C  QQVHj_curr calculates the current  
C     (q1 q2    ---->   W+- H)_mu ,   
C     W+- ---> f1-bar f2
C
C  and crossing related processes. Pauli interference terms for
c  identical fermions are neglected.
c
c*********************************************************************

      subroutine qqHWj_curr(momin, final, N, bos, fsign)

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"

c     electroweak couplings are taken from KOPPLN

      double precision  clr, xm2, xmg, b
      COMMON /BKOPOU/   CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6)

c variables for the main part of the program      

      ! array length 5 means n_p+n_v
      ! TODO: remove or include leptons consistently
      double complex final(0:3)
      double precision pbar(0:3,6)
      double precision momin(0:3,5)
      double precision p(0:3,5)
      double precision fqwh(0:4,2)

      integer sign(5), mu, i, N
      integer fsign(5)
      integer bos

      double precision pwh(0:4)
      double complex psi(2,-1:1,2)
      double complex psiwh(2,-1:1, 2)
      double complex jcurr(0:5,-1:1)
   
      logical ldebug, ldebugm
      parameter(ldebugm = .false., ldebug=.false.)

c -------------- BEGIN ENTRY --------------------------------

      do i=1,5 !n_p+n_v
        do mu=0,3
          pbar(mu,i) = momin(mu,i)
        enddo
        sign(i) = fsign(i)
      enddo
      goto 1
c ----
      ENTRY UGDHWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,3)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6) ! Higgs, undefined in Wjj case
        pbar(mu,5) = momin(mu,2)
      enddo
      sign(1) =  1
      sign(2) =  1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) = -1
      goto 1

c ----
      ENTRY DGUHWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,3)
        pbar(mu,2) = momin(mu,1)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6)
        pbar(mu,5) = momin(mu,2)
      enddo
      sign(1) = -1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) = -1
      goto 1

c ----
      ENTRY UDGHWj_curr(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,2)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6)
        pbar(mu,5) = momin(mu,3)
      enddo
      sign(1) =  1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) =  1
      goto 1

1     continue

c---------------------------- BEGIN CODE -----------------------

c define the internal momenta

      if(ldebug) then
         print*, 'in qqHWj_curr.F'
         print*, "pbar ", pbar
      endif
      do mu = 0,3
         do i = 1,5 !n_p+n_v
            p(mu,i) = pbar(mu,i)*sign(i)
         enddo
         pwh(mu) = qw0(mu,N) + qh(mu,N)
!         if (bos.eq.36 .or. bos.eq.46) then
!         elseif (bos.eq.3 .or. bos.eq.4) then
!            pwh(mu) = qw0(mu,N) 
!         endif
      enddo
      pwh(4) = pwh(0)**2-pwh(1)**2-pwh(2)**2-pwh(3)**2

c----------------------------------------------------------------------

c get the external quark spinors (including factor sqrt(2E) )

      call psi0m(2,pbar(0,1),sign(1),psi)

!c get other bra and ket vectors


      if(ldebug) then
         print*, "qw0", qw0
         print*, "qh", qh
         print*, "pwh", pwh
!         print*, "wtowh", wtowh(0:5,N)
      endif
cc WH via current
      call bra2c(psi(1,-1,2),.true.,p(0,2),-1,pwh(0),wtowh(0,N),
     &              psiwh(1,-1,2),fqwh(0,2))
      call ket2c(psi(1,-1,1),.true.,p(0,1),-1,pwh(0),wtowh(0,N),
     &              psiwh(1,-1,1),fqwh(0,1))


c----------------------------------------------------------------------

c Build current
c (number of vector bosons attached)
c <0|3>
         call curr6(-1,psi(1,-1,2),p(0,2),
     &          psiwh(1,-1,1),fqwh(0,1),jcurr(0,-1) )
c <0|3>
         call curr6add(-1,psiwh(1,-1,2),fqwh(0,2),
     &          psi(1,-1,1),p(0,1),jcurr(0,-1) )
      
c----------------------------------------------------------------------

c store current in output variable
c add remaining coupling factors g_s and g_Wff

       do i=0,3
         final(i) = -jcurr(i,-1) * gg(1) * clr(1,3,-1)
      enddo

      if(ldebug) then
         print*, "final ", final
         print*, "psiwh", psiwh
         print*, "fqwh", fqwh
      endif

      return
      end

c======================================================================

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculate squared matrix element with arguments suitable for dipoles 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine UGDHWj_msq(momin, final, N, bos)

      implicit none

#include "VBFNLO/utilities/global.inc"

      ! 5 or 6 = n_p + n_v in the following lines
      double precision momin(0:3,5), final !, final2
      double precision pbar(0:3,6), qbar(0:4)
      double precision avgfac
      integer N, bos, sign(6), gsign
      integer mu, nlo

      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,3)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6)
        pbar(mu,5) = momin(mu,2)
        qbar(mu) = pbar(mu,5)
      enddo
      sign(1) =  1
      sign(2) =  1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) = -1
      gsign = sign(5)
      avgfac = 1d0/4d0/3d0/8d0
      goto 1

c ----
      ENTRY DGUHWj_msq(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,3)
        pbar(mu,2) = momin(mu,1)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6)
        pbar(mu,5) = momin(mu,2)
        qbar(mu) = pbar(mu,5)
      enddo
      sign(1) = -1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) = -1
      gsign = sign(5)
      avgfac = 1d0/4d0/3d0/8d0
      goto 1

c ----
      ENTRY UDGHWj_msq(momin, final, N, bos)
      do mu=0,3
        pbar(mu,1) = momin(mu,1)
        pbar(mu,2) = momin(mu,2)
        pbar(mu,3) = momin(mu,4)
        pbar(mu,4) = momin(mu,5)
!        pbar(mu,5) = momin(mu,6)
        pbar(mu,5) = momin(mu,3)
        qbar(mu) = pbar(mu,5)
      enddo
      sign(1) =  1
      sign(2) = -1
      sign(3) =  1
      sign(4) = -1
!      sign(5) =  1
      sign(5) =  1
      gsign = sign(5)
      avgfac = 1d0/4d0/3d0/3d0 ! 1/4 helicity of quarks 
      goto 1

1     continue
      nlo = 0

      ! 6th entry is dummy argument
      ! qqHWj only uses pbar(:,1:3) and qbar
      sign(6)=1
      pbar(0:3,6) = (/0d0,0d0,0d0,0d0/)

      ! ignore g mass entry
      qbar(4) = 0d0

      call qqHWj(pbar, sign, qbar, gsign, nlo, N, final, -bos)
c      call qqWHj(pbar, sign, qbar, gsign, final2, bos)
c      if (final/final2 -1 .gt. 1d-5) then
c         print*, "in qqHWj_curr: ", final, final2
c      endif
      final = final * avgfac

      return
      end

