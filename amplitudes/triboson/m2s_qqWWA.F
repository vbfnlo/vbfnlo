*****************************************************************************
c
c    begin function m2s_qqWWA
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c       Adapted from m2s_qqWWZ by V. Hankele
c       Last modified May 2008 (giuseppe)
c
c       this function computes |ME|^2 for qq-> ve e+ mu-vm~ A ("WWA")
c	at LO and virtual corrections via call of 'qqWWA.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqWWA(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn    !in:  random number array
     &                        )
      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)
      complex*16 mg2_tens(0:5), mg2_tens2(0:5)
      double precision sumtest

c declare input/output variables

      INTEGER bos, nlo
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

         write(6,*) " "
         write(6,*) "W+W-A amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowwa_reset

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
         qp(mu,1) = v(mu,1,1)+v(mu,2,1)   !W+
         qm(mu,1) = v(mu,3,1)+v(mu,4,1)   !W-
         qz(mu,1) = v(mu,5,1)             !A
         qww(mu,1) = qp(mu,1) + qm(mu,1)
         qwwz(mu,1) = qww(mu,1) + qz(mu,1)
      enddo

      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qm(4,1) = dotrr(qm(0,1),qm(0,1))
      qz(4,1) = dotrr(qz(0,1),qz(0,1))
      qww(4,1) = dotrr(qww(0,1),qww(0,1))
      qwwz(4,1) = dotrr(qwwz(0,1),qwwz(0,1))

c for WWA->4l precalculate A->WWA,Z->WWA etc..

c lepton spinors and W+, W-, A polarisation vectors
 
c select helicity: h = random number for photon helicity (h=1:2)

      h = mod(h,2) + 1
      it = (-1)**(h+1)

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))          !ve 
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))          !e+        
      CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu(1,1)) !W(1,5))          !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm(1,1)) !W(1,6))          !vm~
      CALL VXXXXX(qz(0,1),ZERO ,it,1,atau(1,1)) !W(1,7))          !A
      CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, wp(1,1))    !W(1,8) W+
      CALL JIOXXX(wvm(1,1),wmu(1,1),GWF,WMASS,WWIDTH, wm(1,1))    !W(1,9) W-

#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c wpp=wp-xp*qp, wmp=wm-xm*qm, ztaup=ztau-xz*qz, ataup=atau-xa*qz
c better choice
            qq = qp(0,1)*qww(0,1) - qp(1,1)*qww(1,1)-
     &           qp(2,1)*qww(2,1) - qp(3,1)*qww(3,1)
            xp = dotrc(qww(0,1),wp(1,1))/qq
            qq = qm(0,1)*qww(0,1) - qm(1,1)*qww(1,1)-
     &           qm(2,1)*qww(2,1) - qm(3,1)*qww(3,1)
            xm = dotrc(qww(0,1),wm(1,1))/qq
            qq = qz(0,1)*qww(0,1) - qz(1,1)*qww(1,1)-
     &           qz(2,1)*qww(2,1) - qz(3,1)*qww(3,1)
            xa=0d0
            do mu = 0,3
               wpp(mu+1,1) = wp(mu+1,1) - xp*qp(mu,1)
               wmp(mu+1,1) = wm(mu+1,1) - xm*qm(mu,1)
               ataup(mu+1,1) = atau(mu+1,1)
            enddo
            do mu = 5,6
               wpp(mu,1) = wp(mu,1)
               wmp(mu,1) = wm(mu,1)
               ataup(mu,1) = atau(mu,1)
            enddo

         endif
#endif


c leptonic tensors:
      SELECT CASE(procid)
      CASE(WPhadWMA, WPWMhadA)

         call anomal_formfactor(qwwz(0,1),qp(0,1),qm(0,1),qz(0,1))
         call WWA_atoww_had(v(0,1,1),1,aww(0,1),0)   !A -> W+ W- -> ve e+ mu- vm~
         call WWA_ztoww_had(v(0,1,1),1,zww(0,1),0)   !Z -> W+ W- -> ve e+ mu- vm~
         call wptowa_had(v(0,1,1),1,h,wpwa(0,1),0)     !W+ -> W+ A -> ve e+ a
         call wmtowa_had(v(0,1,1),1,h,wmwa(0,1),0)     !W- -> W- A -> mu- vm~ a
         call az_towwa_had(v(0,1,1),1,h,awwa(0,1),zwwa(0,1),0) !A/Z -> WWA -> ve e+ mu- vm~ a


      CASE DEFAULT

      if (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant WWA mass
         call anomal_formfactor(qwwz(0,1),qp(0,1),qm(0,1),qz(0,1))

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         call TRIBOSatoww_anomal(v(0,1,1),1,aww(0,1),0)   !A -> W+ W- -> ve e+ mu- vm~
         call TRIBOSztoww_anomal(v(0,1,1),1,zww(0,1),0)   !Z -> W+ W- -> ve e+ mu- vm~
         call wptowa_anomal(v(0,1,1),1,h,wpwa(0,1),0)     !W+ -> W+ A -> ve e+ a
         call wmtowa_anomal(v(0,1,1),1,h,wmwa(0,1),0)     !W- -> W- A -> mu- vm~ a
         call az_towwa_anomal(v(0,1,1),1,h,awwa(0,1),zwwa(0,1),0) !A/Z -> WWA -> ve e+ mu- vm~ a

  
      else  ! SM

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         call TRIBOSatoww(v(0,1,1),1,aww(0,1),0)   !A -> W+ W- -> ve e+ mu- vm~
         call TRIBOSztoww(v(0,1,1),1,zww(0,1),0)   !Z -> W+ W- -> ve e+ mu- vm~
         call wptowa(v(0,1,1),1,h,wpwa(0,1),0)     !W+ -> W+ A -> ve e+ a
         call wmtowa(v(0,1,1),1,h,wmwa(0,1),0)     !W- -> W- A -> mu- vm~ a
!    old tensor routines
!       call atowwa(v(0,1,1),1,h,awwa(0,1),0)     !A -> WWA -> ve e+ mu- vm~ a
!       call ztowwa(v(0,1,1),1,h,zwwa(0,1),0)     !Z -> WWA -> ve e+ mu- vm~ a
         call az_towwa(v(0,1,1),1,h,awwa(0,1),zwwa(0,1),0) !A/Z -> WWA -> ve e+ mu- vm~ a

      endif ! anom or SM

      END SELECT


#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         print*, ""
         call wptowa_mg2_had(v(0,1,1),1,h,mg2_tens,0)     !W+ -> W+ A -> ve e+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wpwa(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor w+ to w+a:", sumtest
         call wmtowa_mg2_had(v(0,1,1),1,h,mg2_tens,0)     !W+ -> W+ A -> ve e+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wmwa(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor w- to w-a:", sumtest
         call atowwa_mg2_had(v(0,1,1),1,h,mg2_tens,0)     !A -> WWA -> ve e+ mu- vm~ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( awwa(i,1) / mg2_tens(i) - 1)
           ! print*, awwa(i,1) , mg2_tens(i)
         enddo
         print*, "tensor a to wwa:", sumtest
         call ztowwa_mg2_had(v(0,1,1),1,h,mg2_tens2,0)     !Z -> WWA -> ve e+ mu- vm~ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zwwa(i,1) / mg2_tens2(i) - 1)
           ! print*, zwwa(i,1) , mg2_tens(i)
         enddo
         print*, "tensor z to wwa:", sumtest
         read(*,*)
      endif
#endif


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
      fsign(6) = -1
      fsign(7) = 1

      call qq_WWA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)

C*******************  q1bar q2 ---> W W A   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WWA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)

c*****************  end of process evaluation  *********************

      SELECT CASE(process)
      CASE(WPhadWMA)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPWMhadA)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      END SELECT

      if (init.eq.1) then
         init = init+1
      endif

      m2s = 0
      DO IRES = 1,NMAX
         m2s = m2s + RES(IRES)
      ENDDO

c...Les Houches interface - the most propable subprocess   
      if ((lha.or.hepmc).and..not.doNLO) then
         helicity(1) = -1
         helicity(2) =  1
         helicity(3) = -1
         helicity(4) =  1
         helicity(5) =  it

         i=0
         weight=0.d0
         rnumb=RandomNumber()
         do while((i.le.nmax).and.(weight.le.rnumb*m2s))
            i=i+1
            weight=weight+res(i)
            iprocess=i
         enddo
         if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
            rnumb=RandomNumber()
            SELECT CASE(process)
            CASE(WPhadWMA)
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            CASE(WPWMhadA)
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            END SELECT
         endif
      endif

      m2s = m2s * 2d0    ! factor 2 for random helicity summation

      m2s_qqWWA = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWWA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WWA
c
c*****************************************************************

      subroutine qq_WWA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax

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
      DOUBLE PRECISION uuwwa, ddwwa, uuwwa2, ddwwa2
      DOUBLE PRECISION pbar(0:3,2+n_v), dotrr, polcol
      DOUBLE PRECISION NCmatrixelt(0:1,0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(4), diagToPhys(2), FL_VVV
      EXTERNAL dotrr, FL_VVV


      LOGICAL ldebug
      PARAMETER (ldebug=.false.)
      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)
      PARAMETER (zero = (0d0,0d0) )
      INTEGER i,if1,k,mu

      uuwwa = 0d0
      uuwwa2 = 0d0
      ddwwa = 0d0
      ddwwa2 = 0d0

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

      do mu = 0,3             ! kinematics for WWA
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! vm~
         pbar(mu,7) = v(mu,5) ! A
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         call suu_wwa(pbar,fsign,uuwwa2)
         call sdd_wwa(pbar,fsign,ddwwa2)
      endif
#endif

      call QQVVA(pbar,fsign,nlo,1,uuwwa, ddwwa)

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         call wptowa_mg2_had(v(0,1),1,1,wpwa(0,1),0)     !W+ -> W+ A -> ve e+ a
         call wmtowa_mg2_had(v(0,1),1,1,wmwa(0,1),0)     !W- -> W- A -> mu- vm~ a
         call atowwa_mg2_had(v(0,1),1,1,awwa(0,1),0)     !A -> WWA -> ve e+ mu- vm~ a
         call ztowwa_mg2_had(v(0,1),1,1,zwwa(0,1),0)     !Z -> WWA -> ve e+ mu- vm~ a

         call QQVVA(pbar,fsign,nlo,1,uuwwa2, ddwwa2)

         call wptowa_had(v(0,1),1,1,wpwa(0,1),0)     !W+ -> W+ A -> ve e+ a
         call wmtowa_had(v(0,1),1,1,wmwa(0,1),0)     !W- -> W- A -> mu- vm~ a
         call az_towwa_had(v(0,1),1,1,awwa(0,1),zwwa(0,1),0) !A/Z -> WWA -> ve e+ mu- vm~ a

         print*,"1 - uuwwa/uuwwa2 =   ",abs(1-uuwwa/uuwwa2)
         print*,"1 - ddwwa/ddwwa2 =   ",abs(1-ddwwa/ddwwa2)
         read(*,*)
      endif
      if (ldebug) then
         print*,"uuwwa =   ",uuwwa
         print*,"uuwwa2 =  ",uuwwa2
         print*,"ddwwa =   ",ddwwa
         print*,"ddwwa2 =  ",ddwwa2
      endif
#endif

      NCmatrixelt(0,0)=uuwwa
      NCmatrixelt(1,1)=ddwwa

      do if1=1,(nfl/2)*2
            iflav(1)=if1*fsign(physToDiag(1))
            iflav(2)=if1*fsign(physToDiag(2))
            k=FL_VVV(iflav,1)

            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *NCmatrixelt(mod(if1,2),mod(if1,2))*polcol

      enddo 

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WWA
c
c*****************************************************************

      subroutine vtowwa_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            aww(mu,j)  = zero
            zww(mu,j)  = zero
            awwa(mu,j) = zero
            zwwa(mu,j) = zero
            wpwa(mu,j) = zero
            wmwa(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
