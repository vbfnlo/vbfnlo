*****************************************************************************
c
c    begin function m2s_qqAA
c
c*****************************************************************************
c
c       michael rauch, <rauch@particle.physik.uni-karlsruhe.de>
c       Adapted from m2s_qqAAA
c       Last modified Feb 2011 (rauch)
c
c     this function computes |ME|^2 for qq -> A A 
c     at LO and virtual corrections via call of 'qqAA.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqAA(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )

      use globalvars, only: ldoblha, lglobalprint

      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
#include "tensor.inc"

c declare input/output variables
      INTEGER bos, nlo
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions
      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VV

c  helicity selection
      INTEGER h
      COMMON /hcount / h

      DOUBLE PRECISION zero
      PARAMETER (zero=0d0)

      DOUBLE PRECISION q_sf
      INTEGER init/0/, I, mu
      SAVE init

c declare variables for summation process
      INTEGER ires

      INTEGER FSIGN(2+n_v), physToDiag(2), nmin, nmax
      SAVE nmax

c store contributions from subprocess ip in res(ip) where
      DOUBLE PRECISION res(maxnumsubproc)

c declare parton distribution variables
      DOUBLE PRECISION pdf(-6:6,2)

c define program switches
      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c ------------- BEGIN CODE -------------------

      if ( (init .eq. 0).and.(.not.ldoblha) ) then
         
         if (lglobalprint) then
            write(6,*) " "
            write(6,*) "AA amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif

         call printnfl(.false.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else if ((init.eq.0).and.ldoblha) then
         init = 1
         h = 1
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif
      
      do mu = 0,3
         qa1(mu,1) = v(mu,1,1)! A1
         qa2(mu,1) = v(mu,2,1)! A2
         qaa(mu,1) = v(mu,1,1) + v(mu,2,1)! A1+A2
      enddo

      qa1(4,1) = dotrr(qa1(0,1),qa1(0,1))
      qa2(4,1) = dotrr(qa2(0,1),qa2(0,1))
      qaa(4,1) = dotrr(qaa(0,1),qaa(0,1))

c photon polarisation vectors
c select helicity: h = random number for photon helicity (h=1:2)

      if (ldoblha) then
         h = blha_ranhelcomb
      end if
     
      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         helicity(1)= ie
         helicity(2)= iu
      end if
      
c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)= ie
           helicity(2)= iu
        endif

      CALL VXXXXX(qa1(0,1),ZERO,ie,1,atau1(1,1)) !A1
      CALL VXXXXX(qa2(0,1),ZERO,iu,1,atau2(1,1)) !A2

      if (ldoblha) then
         physToDiag(1)=blha_physdiagmap(1,blha_cursubproc,blha_curproc)
         physToDiag(2)=blha_physdiagmap(2,blha_cursubproc,blha_curproc)
         fsign(1) = blha_fsign(1, blha_cursubproc, blha_curproc)
         fsign(2) = blha_fsign(2, blha_cursubproc, blha_curproc)
         fsign(3) = 1
         fsign(4) = 1

         call qq_AA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &        res,nmin,nmax)
         return
      endif

      
      if ((nlo .eq. -8) .or. (nlo.eq.1)) 
     &  call aa_htoaa(v(0,1,1),1,w4l(h,1),0) !H -> A A
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
         do i = 1, max_p-1
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Now sum over the subprocesses contributing to A A production

      nmax = FL_VV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar1 ---> A A   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = 1

      call qq_AA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)

C*******************  q1bar q1 ---> A A   **********************

      if (nlo .ne. -8) then
cc gg is symmetric -> would be double-counting
      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_AA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,
     &     res,nmin,nmax)
      endif

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
      m2s = m2s * 4d0    ! factor 4 for random helicity summation
      m2s_qqAA = m2s
      RETURN
      END
c*****************************************************************************
c
c    end function m2s_qqAA
c
c*****************************************************************************
c******************************************************************
c
c   begin subroutine qq_AA
c
c*****************************************************************

      subroutine qq_AA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                 pdf,res,nmin,nmax)

      use globalvars, only: ldoblha
     
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,n_p), v(0:3,n_v), xi(nx)
      double precision pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax

c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       leptons and photon momenta
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c  OUTPUT: res(k)
c
c  In and Output:
c          nmin, nmax      range of process numbers for this call
c
      double precision N                              ! color factors
      parameter(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      double precision uuaa, ddaa, uuaa2, ddaa2
      double precision pbar(0:3,2+n_v), dotrr, polcol, polcolg
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VV
      external dotrr, FL_VV

      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      uuaa = 0d0
      uuaa2 = 0d0
      ddaa = 0d0
      ddaa2 = 0d0

      nmin = nmax+1
      do i = 1,2
         diagToPhys(physToDiag(i)) = i
      enddo

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))
      polcolg = 1d0/(4d0*(N**2-1)**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo

      do mu = 0,3             ! kinematics for AA
         pbar(mu,3) = v(mu,1) ! A1
         pbar(mu,4) = v(mu,2) ! A2
      enddo

      if ((nlo .eq. -8).or.((nlo.eq.1).and.(physToDiag(1).eq.1))) then
c gluon-fusion loop contributions
        call GGAA(pbar,fsign,nlo,1,uuaa)
        iflav(1)=0
        iflav(2)=0
        k=FL_VV(iflav,1)
  
        res(k)=pdf(sign1*iflav(1),1)
     &        *pdf(sign2*iflav(2),2)
     &        *uuaa*polcolg

      endif
      if (nlo .ne. -8) then
#ifdef WITH_MADGRAPH
        if (ldebug) then
           call suu_aa(pbar,fsign,uuaa2)
           call sdd_aa(pbar,fsign,ddaa2)
        endif
#endif

        call QQAA(pbar,fsign,nlo,1,uuaa,ddaa)

#ifdef WITH_MADGRAPH
        if (ldebug) then
           print*,"uuaa =   ",uuaa
           print*,"uuaa2 =  ",uuaa2
           print*,"ddaa =   ",ddaa
           print*,"ddaa2 =  ",ddaa2
           read(*,*)
        endif
#endif

        if (ldoblha) then
           return
        endif

        NCmatrixelt(0,0)=uuaa
        NCmatrixelt(1,1)=ddaa

        do if1=1,nfl
              iflav(1)=if1*fsign(physToDiag(1))
              iflav(2)=if1*fsign(physToDiag(2))
              k=FL_VV(iflav,1)
  
              res(k)=pdf(sign1*iflav(1),1)
     &              *pdf(sign2*iflav(2),2)
     &              *NCmatrixelt(mod(if1,2),mod(if1,2))*polcol
  
        enddo 

      endif

      nmax = FL_VV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_AA
c
c*****************************************************************

