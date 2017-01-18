c*********************  subroutine m2s_qqWA *********************************        
c
c       Adapted from m2s_VV by Vera Hankele
c
c       Johannes Bellm, <johannes@particle.physik.uni-karlsruhe.de>
c	Initial version: Juli 2011
c	Last modified:
c
c***************************************************************************
      real*8 function m2s_qqWA(
     &      bos,  !in:  XX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-, X=6 H
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )
      implicit none
c
c declare input/output variables
c
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"

      integer bos, nlo
      real*8 xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c  helicity selection
c      integer jsig, jsig1, jsig3, jsig5
c      common /chelsum/ jsig,jsig1,jsig3,jsig5


c  helicity selection
      INTEGER h, it 
      COMMON /hcount / h


      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

c declare external functions

      real*8 dotrr
      integer FL_VV
      external dotrr
      external FL_VV

c declare local variables

      real*8 q_sf
      integer init/0/, I, mu
      save init

c declare variables for summation process

      INTEGER ires

      integer FSIGN(2+n_v), physToDiag(2), nmin, nmax
      save nmax
c
c store contributions from subprocess ip in res(ip) where
      real*8 res(maxnumsubproc)
c
c declare parton distribution variables
c
      real*8 pdf(-6:6,2)
c
c define program switches
c
      logical ldebug
      parameter (ldebug=.false.)		!output debug information

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

      if ( init .eq. 0 ) then
         if(bos.eq.31) then
            write(6,*) " "
            write(6,*) "W+A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if(bos.eq.41) then
            write(6,*) " "
            write(6,*) "W-A amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else
            write(6,*) " wrong call of m2s_qqWA, bos = ",bos
            write(6,*) "-----------------------------------------------"
         endif
         call DIBOSwtowa_reset

         call printnfl(.true.)
!          h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif


c...Les Houches interface

c lepton helicities for Z -> ta- ta+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity combination (h=1:2)
        h = mod(h,2) + 1
	it = (-1)**(h+1)

        if ((lha.or.hepmc) .and. (.not. doNLO)) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=it
        endif

c for WA->2l+A precalculate W->WA leptonic tensors

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3)       !ve
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4)       !e+          
      CALL VXXXXX(v(0,3,1),ZERO ,it,1,an(1,1))  !W(1,5)       !A


      do mu = 0,3
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+
         qa(mu,1) = v(mu,3,1)
      enddo
      qp(4,1) = dotrr(qp(0,1),qp(0,1))


c leptonic tensors
       if (with_anom) then ! anomalous gauge boson couplings    noanomal
                if (bos.eq.31) then
                      call DIBOSwptowpa_an(v,1,wtowa(0,1))
                else if (bos.eq.41) then
                      call DIBOSwmtowma_an(v,1,wtowa(0,1))
                endif
       else  ! SM
               if (bos.eq.31) then
                     call DIBOSwptowpa(v,1,wtowa(0,1))
               else if (bos.eq.41) then
                     call DIBOSwmtowma(v,1,wtowa(0,1))
               endif
       endif


c later you need the an: wp: w+->e+ve

      CALL JIOXXX(wep(1,1),wve(1,1),GWF ,WMASS,WWIDTH,wp(1,1)) !!W+


c if desired, print output information

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f10.3, 2x) )
 20      format( " v(", i1, ") = ", 4(f10.3, 2x) )
c
         write(6,*) " "
         do i = 1, n_v
            write(6,20) i, v(0,i,1), v(1,i,1), v(2,i,1), v(3,i,1)
         end do
         print*," parton momenta "
         do i = 1, n_p-1
            write(6,10) i, p(0,i), p(1,i), p(2,i), p(3,i)
         end do
c
         write(6,*) " "
         write(6,*) "xi(1), xi(2) =", xi
         write(6,*) "mu_f =",sqrt(mufsq(1,1))
      end if

c
c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
c
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

        
      DO IRES = 1,NMAX
         RES(IRES)=0d0
      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VV(fsign,-1)            !reset counter for subprocesses to 0

C*******************  q1 q2 ---> WA   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1 
      fsign(2) = -1 
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1


      call qq_WA(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax) 


C*******************  q1<->q2 ---> WA   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WA(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax)

c*****************  end of process evaluation  *********************


      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO


c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
      endif

      m2s = m2s*2d0     ! factor 2 for random helicity summation

      m2s_qqWA = m2s
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WA
c
c*****************************************************************

      subroutine qq_WA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      real*8 p(0:3,max_p-1), v(0:3,n_v), xi(nx)
      real*8 pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
c
c wbf_zh calls the amplitude square routines 
c             qqWA    for qq --> VA --> 2 leptons + A
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      double precision N                              ! color factors
      parameter(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      real*8 udwa, udwa2
      real*8 pbar(0:3,2+n_v), dotrr, polcol
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV


      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      udwa = 0d0 
      udwa2 = 0d0


      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo


         do mu = 0,3             ! kinematics for WA --> 2 lepton + A decay
            pbar(mu,3) = v(mu,1) ! ve
            pbar(mu,4) = v(mu,2) ! e+
            pbar(mu,5) = v(mu,3) ! A
         enddo

#ifdef WITH_MADGRAPH
         if(ldebug) then

            if (bos.eq.31) then
               call sudwpa(pbar,fsign,udwa2)
            else if (bos.eq.41) then
               call sduwma(pbar,fsign,udwa2)
            endif
         endif
#endif

         call QQWA(pbar,fsign,nlo,1,udwa, bos)

#ifdef WITH_MADGRAPH
         if (ldebug) then
            if(abs(udwa/udwa2 -1).ge.1.0d-6) then
               print*,"udwa =   ",udwa
               print*,"udwa2 =  ",udwa2
            endif
         endif
#endif


      NCmatrixelt(0,0)=udwa


      do if1=2,(nfl/2)*2,2
         if (bos.eq.31) then
            iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
            iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
         else if (bos.eq.41) then
            iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
            iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
         endif
         k=FL_VV(iflav,1)

         res(k)=pdf(sign1*iflav(1),1)
     &         *pdf(sign2*iflav(2),2)
     &         *NCmatrixelt(0,0)*polcol

      enddo 

      nmax = FL_VV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WA
c
c*****************************************************************


      subroutine DIBOSwtowa_reset
      implicit none
#include "VBFNLO/utilities/global.inc"
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,3
         do j = 1,3
            wtowa(mu,j) = zero            
         enddo
      enddo
      print*," W to WA  tensors initialized to ",zero 
      return
      end

      
