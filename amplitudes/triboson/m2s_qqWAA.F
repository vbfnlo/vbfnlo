*****************************************************************************
c
c    begin function m2s_qqWAA
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.uni-karlsruhe.de>
c       Adapted from m2s_qqWWA
c       Last modified December 2009 (giuseppe)
c
c       this function computes |ME|^2 for qq-> ve(ve~) e+(-) A A ("W+-AA")
c	at LO and virtual corrections via call of 'qqWAA.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqWAA(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      ps_number,
     &      NW
     &                        )
      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"

c declare input/output variables

      INTEGER bos, nlo, NW, ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions

      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

c declare local variables

c  helicity selection
 
      INTEGER h
      COMMON /hcount / h

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

      DOUBLE PRECISION q_sf, qq
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

      if ( init .eq. 0 ) then

         if (bos.eq.311) then
            write(6,*) " "
            write(6,*) "AAW+ amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         elseif(bos.eq.411) then
            write(6,*) " "
            write(6,*) "AAW- amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         endif
         call vtowaa_reset

         call printnfl(.true.)
         h = 1
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif

      do mu = 0,3
         qw0(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+-
         qa1(mu,1) = v(mu,3,1)             !A1
         qa2(mu,1) = v(mu,4,1)             !A2
         qaa(mu,1) = qa1(mu,1)+qa2(mu,1)
         qaaw(mu,1) = qaa(mu,1)+qw0(mu,1)
      enddo

      qw0(4,1) = dotrr(qw0(0,1),qw0(0,1))
      qa1(4,1) = dotrr(qa1(0,1),qa1(0,1))
      qa2(4,1) = dotrr(qa2(0,1),qa2(0,1))
      qaa(4,1) = dotrr(qaa(0,1),qaa(0,1))
      qaaw(4,1) = dotrr(qaaw(0,1),qaaw(0,1))

c for WAA->2lAA precalculate leptonic tensors

c lepton spinors and W+, A polarisation vectors
 
c select helicity: h = random number for photon helicity (h=1:2)

      h = mod(h,4) + 1
      ie = (-1)**(h+1)
      iu = (-1)**((h-1)/2)

c Les Houches interface
      if((lha.or.hepmc).and.(.not.doNLO)) then
         helicity(1) = -1
         helicity(2) =  1
         helicity(3) =  iu
         helicity(4) =  it
      endif

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1))           !ve 
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1))           !e+     
      CALL VXXXXX(qa1(0,1),ZERO ,ie,1,atau1(1,1))          !A1
      CALL VXXXXX(qa2(0,1),ZERO ,iu,1,atau2(1,1))          !A2
      CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, w0(1,1))!W+- 

#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
            qq = qaa(0,1)*qw0(0,1) - qaa(1,1)*qw0(1,1)-
     &           qaa(2,1)*qw0(2,1) - qaa(3,1)*qw0(3,1)
            xw=dotrc(qaa(0,1),w0(1,1))/qq
            do mu = 0,3
               w0p(mu+1,1) = w0(mu+1,1) - xw*qw0(mu,1)
               atau1p(mu+1,1) = atau1(mu+1,1)
               atau2p(mu+1,1) = atau2(mu+1,1)
            enddo
            do mu = 5,6
               w0p(mu,1) = w0(mu,1)
               atau1p(mu,1) = atau1(mu,1)
               atau2p(mu,1) = atau2(mu,1)
            enddo
	     xatau1=0
	     xatau2=0
         endif
#endif


      if (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant WAA mass
         call anomal_formfactor(qaaw(0,1),qw0(0,1),qa1(0,1),qa2(0,1))


c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         if (bos.eq.311) then
            call wptowanew_anomal(v(0,1,1),1,h,wpwa1(0,1),1,NW)!W+ -> W+ A -> ve e+ a
            call wptowanew_anomal(v(0,1,1),1,h,wpwa2(0,1),2,NW)!W+ -> W+ A -> ve e+ a
            call wptowaa_anomal(v(0,1,1),1,h,wpwaa(0,1),0,NW)  !W+ -> ve e+ a a
         elseif (bos.eq.411) then
            call wmtowanew_anomal(v(0,1,1),1,h,wmwa1(0,1),1,NW)!W- -> W- A -> mu- vm~ a
            call wmtowanew_anomal(v(0,1,1),1,h,wmwa2(0,1),2,NW)!W- -> W- A -> mu- vm~ a
            call wmtowaa_anomal(v(0,1,1),1,h,wmwaa(0,1),0,NW)  !W- -> ve e+ a a
         endif


      else !SM

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         if (bos.eq.311) then
            call wptowanew(v(0,1,1),1,h,wpwa1(0,1),1,NW)!W+ -> W+ A -> ve e+ a
            call wptowanew(v(0,1,1),1,h,wpwa2(0,1),2,NW)!W+ -> W+ A -> ve e+ a
            call wptowaa(v(0,1,1),1,h,wpwaa(0,1),0,NW)  !W+ -> ve e+ a a
         elseif (bos.eq.411) then
            call wmtowanew(v(0,1,1),1,h,wmwa1(0,1),1,NW)!W- -> W- A -> mu- vm~ a
            call wmtowanew(v(0,1,1),1,h,wmwa2(0,1),2,NW)!W- -> W- A -> mu- vm~ a
            call wmtowaa(v(0,1,1),1,h,wmwaa(0,1),0,NW)  !W- -> ve e+ a a
         endif

      endif !anom or SM

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

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VVV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W W A   **************************

c   physToDiag(ext.momentum label) = Feynman diagram label

      physToDiag(1)=1     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=2

      fsign(1) = 1
      fsign(2) = -1
      fsign(3) = 1
      fsign(4) = -1
      fsign(5) = 1
      fsign(6) = 1

      call qq_WAA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)

C*******************  q1bar q2 ---> W W A   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WAA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc).and..not.doNLO) then
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

      m2s_qqWAA = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWAA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WAA
c
c*****************************************************************

      subroutine qq_WAA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax,ps_number

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
      DOUBLE PRECISION N                              ! color factors
      PARAMETER(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      DOUBLE PRECISION udwaa, duwaa, udwaa2, duwaa2
      DOUBLE PRECISION pbar(0:3,2+n_v), dotrr, polcol
      DOUBLE PRECISION NCmatrixelt(0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(4), diagToPhys(2), FL_VVV
      EXTERNAL dotrr, FL_VVV


      LOGICAL ldebug
      PARAMETER (ldebug=.false.)
      PARAMETER (zero = (0d0,0d0) )
      INTEGER i,if1,k,mu

      udwaa = 0d0
      duwaa = 0d0
      udwaa2 = 0d0
      duwaa2 = 0d0

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

      do mu = 0,3             ! kinematics for WAA
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! A1
         pbar(mu,6) = v(mu,4) ! A2
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         if (bos.eq.311) then
            call sud_waa(pbar,fsign,udwaa2)
         else if (bos.eq.411) then
            call sdu_waa(pbar,fsign,duwaa2)
         endif
      endif
#endif

      if (bos.eq.311) then
         call QQAAW(pbar,fsign,nlo,1,udwaa,bos,ps_number)
      elseif (bos.eq.411) then
         call QQAAW(pbar,fsign,nlo,1,duwaa,bos,ps_number)
      endif

      if (ldebug) then
         if (bos.eq.311) then
            print*,"udwaa =   ",udwaa
#ifdef WITH_MADGRAPH
            print*,"udwaa2 =   ",udwaa2
#endif
         elseif (bos.eq.411) then
            print*,"duwaa =   ",duwaa
#ifdef WITH_MADGRAPH
            print*,"duwaa2 =   ",duwaa2
#endif
         endif
         read(*,*)
      endif

      NCmatrixelt(0)=udwaa
      NCmatrixelt(1)=duwaa

      do if1=1,(nfl/2)*2,2
         if (bos.eq.311) then
            iflav(1)=(if1-physToDiag(1)+2)*fsign(physToDiag(1))
            iflav(2)=(if1+physToDiag(1)-1)*fsign(physToDiag(2))
            k=FL_VVV(iflav,1)
            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *NCmatrixelt(0)*polcol
         elseif (bos.eq.411) then
            iflav(1)=(if1+physToDiag(1)-1)*fsign(physToDiag(1))
            iflav(2)=(if1-physToDiag(1)+2)*fsign(physToDiag(2))
            k=FL_VVV(iflav,1)
            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *NCmatrixelt(1)*polcol
         endif


      enddo 

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WAA
c
c*****************************************************************

      subroutine vtowaa_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            wpwa1(mu,j) = zero
            wpwa2(mu,j) = zero
            wpwaa(mu,j) = zero
            wmwa1(mu,j) = zero
            wmwa2(mu,j) = zero
            wmwaa(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
