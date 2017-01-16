*****************************************************************************
c
c    begin function m2s_qqWZA
c
c*****************************************************************************
c
c       michael rauch, <rauch@particle.physik.uni-karlsruhe.de>
c       Adapted from m2s_qqWWA by G. Bozzi
c       Adapted from m2s_qqWWZ by V. Hankele
c       Last modified Jul 2009 (rauch)
c
c       this function computes |ME|^2 for qq-> ve e+ mu- mu+ A ("W+ZA")
c                                     and qq-> ve~ e- mu- mu+ A ("W-ZA")
c	at LO and virtual corrections via call of 'qqWZA.f'
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqWZA(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      ps_number
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

c declare input/output variables

      INTEGER bos, nlo,ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s
      DOUBLE PRECISION q34(0:4,3) !check if I really need these

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

c variables for hadronic decays
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      integer N_gen_W
      real*8 fac_W
      external mjj2
      
      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)
      complex*16 mg2_tens(0:5)
      double precision sumtest


c define program switches

      LOGICAL ldebug
      PARAMETER (ldebug=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c ------------- BEGIN CODE -------------------

      if ( init .eq. 0 ) then
         if(bos.eq.321) then
            write(6,*) " "
            write(6,*) "W+ZA amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if(bos.eq.421) then
            write(6,*) " "
            write(6,*) "W-ZA amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else
            write(6,*) " wrong call of m2s_qqWZA, bos = ",bos
            write(6,*) "-----------------------------------------------"
         endif
         call vtowza_reset

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
         qz(mu,1) = v(mu,3,1)+v(mu,4,1)   !Z 
         qa(mu,1) = v(mu,5,1)             !A
         qwz(mu,1) = qp(mu,1) + qz(mu,1)
         qwza(mu,1) = qwz(mu,1) + qa(mu,1)
      enddo

      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qz(4,1) = dotrr(qz(0,1),qz(0,1))
      qa(4,1) = dotrr(qa(0,1),qa(0,1))
      qwz(4,1) = dotrr(qwz(0,1),qwz(0,1))
      qwza(4,1) = dotrr(qwza(0,1),qwza(0,1))

c for WZA->4l precalculate W->WZA etc..

c lepton spinors and W+, Z, A polarisation vectors
 
c select helicity: h = random number for lepton & photon helicity (h=1:4)

      h = mod(h,4) + 1
      iu = (-1)**(h+1)
      it = (-1)**((h-1)/2)

      if ((lha.or.hepmc) .and. (.not. doNLO)) then
         helicity(1)=-1
         helicity(2)= 1
         helicity(3)=-iu
         helicity(4)= iu
         helicity(5)= it
      end if

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))          !ve 
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))          !e+        
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,1,wmu(1,1)) !W(1,5))          !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,iu,-1,wmp(1,1)) !W(1,6))          !mu+
      CALL VXXXXX(qa(0,1),ZERO ,it,1,wa(1,1)) !W(1,7))          !A
      CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, wp(1,1))    !W(1,8) W+

      SELECT CASE(procid)
      CASE(WPZhadA, WMZhadA)
        ! for general output to all combinations: up-type first
        if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
          ! set couplings for Z hadronic decay into uubar
          call setZtouu
        endif
        CALL JIOXXX(wmp(1,1),wmu(1,1),GZ_ZF,ZMASS,ZWIDTH, ztau(1,1))  !Z(1,9) Z
        CALL JIOXXX(wmp(1,1),wmu(1,1),GZ_AF,ZERO ,ZERO  , atau(1,1))  !Z(1,9) A(Z)

      CASE DEFAULT
        CALL JIOXXX(wmp(1,1),wmu(1,1),GZL,ZMASS,ZWIDTH, ztau(1,1))  !Z(1,9) Z
        CALL JIOXXX(wmp(1,1),wmu(1,1),GAL,ZERO ,ZERO  , atau(1,1))  !Z(1,9) A(Z)
      END SELECT

#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c wpp=wp-xp*qp, wmp=wm-xm*qm, ztaup=ztau-xz*qz, ataup=atau-xa*qz
c choice: WZ
            qq = qp(0,1)*qwz(0,1) - qp(1,1)*qwz(1,1)-
     &           qp(2,1)*qwz(2,1) - qp(3,1)*qwz(3,1)
            xp = dotrc(qwz(0,1),wp(1,1))/qq
            qq = qz(0,1)*qwz(0,1) - qz(1,1)*qwz(1,1)-
     &           qz(2,1)*qwz(2,1) - qz(3,1)*qwz(3,1)
            xz = dotrc(qwz(0,1),ztau(1,1))/qq
            xaz = dotrc(qwz(0,1),atau(1,1))/qq
            xa=0d0
            do mu = 0,3
               wpp(mu+1,1) = wp(mu+1,1) - xp*qp(mu,1)
               ztaup(mu+1,1) = ztau(mu+1,1) - xz*qz(mu,1)
               ataup(mu+1,1) = atau(mu+1,1) - xaz*qz(mu,1)
               wap(mu+1,1) = wa(mu+1,1)
            enddo
            do mu = 5,6
               wpp(mu,1) = wp(mu,1)
               ztaup(mu,1)  = ztau(mu,1)
               ataup(mu,1)  = atau(mu,1)
               wap(mu,1) = wa(mu,1)
            enddo

         endif
#endif


      SELECT CASE(procid)
      CASE(WPhadZA, WPZhadA, WMhadZA, WMZhadA)

         call anomal_formfactor(qwza(0,1),qp(0,1),qz(0,1),qa(0,1))
         if (bos.eq.321) then
            call TRIBOSwptowz_had(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa_had(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza_had(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz_had(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa_had(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza_had(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza_had(v(0,1,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza_had(v(0,1,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a


      CASE DEFAULT

      if (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant WZA mass
         call anomal_formfactor(qwza(0,1),qp(0,1),qz(0,1),qa(0,1))

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         if (bos.eq.321) then
            call TRIBOSwptowz_anomal(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa_anomal(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza_anomal(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz_anomal(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa_anomal(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza_anomal(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza_anomal(v(0,1,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza_anomal(v(0,1,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a


      else !SM

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

         if (bos.eq.321) then
            call TRIBOSwptowz(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza(v(0,1,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza(v(0,1,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a

      endif !anom or SM

      END SELECT


#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then

         print*, ""

         if (bos.eq.321) then

            call TRIBOSwptowz_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !W -> W+ Z -> ve e+ mu- mu+
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wpwz(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w+ to wz:", sumtest

            call TRIBOSwptowa_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !W -> W+ A -> ve e+ a
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wpwa(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w+ to wa:", sumtest

            call wptowza_mg2_had(v(0,1,1),1,h,mg2_tens,0)    !W -> WZA -> ve e+ mu- mu+ a
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zwwa(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w+ to wza:", sumtest

         else if (bos.eq.421) then

            call TRIBOSwmtowz_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !W -> W- Z -> ve~ e- mu- mu+
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wpwz(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w- to wz:", sumtest

            call TRIBOSwmtowa_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !W -> W- A -> ve~ e- a
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wpwa(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w- to wa:", sumtest

            call wmtowza_mg2_had(v(0,1,1),1,h,mg2_tens,0)    !W -> WZA -> ve~ e- mu- mu+ a
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zwwa(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w- to wza:", sumtest

         endif

         call TRIBOSatoza_mg2_had(v(0,1,1),1,h,mg2_tens,0)   !A -> Z A -> mu- mu+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( aza(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to za:", sumtest

         call TRIBOSztoza_mg2_had(v(0,1,1),1,h,mg2_tens,0)   !Z -> Z A -> mu- mu+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zaz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor z to za:", sumtest

         read(*,*)
      endif
#endif


c if desired, print output information

      if ( ldebug ) then
 10      format( " p(", i1, ") = ", 4(f32.16, 2x) )
 20      format( " v(", i1, ") = ", 4(f32.16, 2x) )
c 10.3
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

C*******************  q1 qbar2 ---> W Z A   **************************

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

      call qq_WZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)

C*******************  q1bar q2 ---> W Z A   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WPZhadA, WMZhadA)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        CALL JIOXXX(wmp(1,1),wmu(1,1),GZ_ZF,ZMASS,ZWIDTH, ztau(1,1))  !Z(1,9) Z
        CALL JIOXXX(wmp(1,1),wmu(1,1),GZ_AF,ZERO ,ZERO  , atau(1,1))  !Z(1,9) A(Z)

#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c wpp=wp-xp*qp, wmp=wm-xm*qm, ztaup=ztau-xz*qz, ataup=atau-xa*qz
c choice: WZ
            qq = qp(0,1)*qwz(0,1) - qp(1,1)*qwz(1,1)-
     &           qp(2,1)*qwz(2,1) - qp(3,1)*qwz(3,1)
            xp = dotrc(qwz(0,1),wp(1,1))/qq
            qq = qz(0,1)*qwz(0,1) - qz(1,1)*qwz(1,1)-
     &           qz(2,1)*qwz(2,1) - qz(3,1)*qwz(3,1)
            xz = dotrc(qwz(0,1),ztau(1,1))/qq
            xaz = dotrc(qwz(0,1),atau(1,1))/qq
            xa=0d0
            do mu = 0,3
               wpp(mu+1,1) = wp(mu+1,1) - xp*qp(mu,1)
               ztaup(mu+1,1) = ztau(mu+1,1) - xz*qz(mu,1)
               ataup(mu+1,1) = atau(mu+1,1) - xaz*qz(mu,1)
               wap(mu+1,1) = wa(mu+1,1)
            enddo
            do mu = 5,6
               wpp(mu,1) = wp(mu,1)
               ztaup(mu,1)  = ztau(mu,1)
               ataup(mu,1)  = atau(mu,1)
               wap(mu,1) = wa(mu,1)
            enddo

         endif
#endif

         if (bos.eq.321) then
            call TRIBOSwptowz_had(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa_had(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza_had(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz_had(v(0,1,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa_had(v(0,1,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza_had(v(0,1,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza_had(v(0,1,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza_had(v(0,1,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a


C*******************  q1 qbar2 ---> W Z A   **************************

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

      call qq_WZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)

C*******************  q1bar q2 ---> W Z A   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax,ps_number)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif


      SELECT CASE(procid)
      CASE(WPhadZA, WMhadZA)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPZhadA, WMZhadA)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
        DO IRES = 1,NMAX
          if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
             if (ires.le.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_up      ! factor fac_Z_up for up-type with all combinations
             elseif (ires.gt.(NMAX/2)) then
                res(IRES) = res(IRES) * fac_Z_down    ! factor fac_Z_down for down-type with all combinations
             endif
             if (mod(nmax,2).ne.0) then
               print*, "Something went wrong with the hadronic decays!"
               stop
             endif
          elseif (mod(abs(finalquarks(1)),2).eq.0) then
             res(IRES) = res(IRES) * fac_Z_up         ! up-type final states
          else
             res(IRES) = res(IRES) * fac_Z_down       ! down-type final states
          endif
        ENDDO
      END SELECT


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
         SELECT CASE(process)
         CASE(WPhadZA)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(WMhadZA)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(WPZhadA, WMZhadA)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               if (i.le.(nmax/2)) then    ! up-type
                 finalquarks_psp(1) =  2 + 2* INT(rnumb*N_gen_up)
                 finalquarks_psp(2) = -2 - 2* INT(rnumb*N_gen_up)
               else                       ! down-type
                 finalquarks_psp(1) =  1 + 2* INT(rnumb*N_gen_down)
                 finalquarks_psp(2) = -1 - 2* INT(rnumb*N_gen_down)
               endif
            endif
         END SELECT
      endif

      m2s = m2s * 4d0    ! factor 4 for random helicity summation

      m2s_qqWZA = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWZA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_WZA
c
c*****************************************************************

      subroutine qq_WZA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax,ps_number)
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
      INTEGER ps_number

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
      DOUBLE PRECISION udwza, udwza2 
      DOUBLE PRECISION pbar(0:3,2+n_v), dotrr, polcol
      DOUBLE PRECISION NCmatrixelt(0:1,0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(4), diagToPhys(2), FL_VVV
      EXTERNAL dotrr, FL_VVV


      LOGICAL ldebug
      PARAMETER (ldebug=.false.)
      PARAMETER (zero = (0d0,0d0) )
      INTEGER i,if1,k,mu
      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

c  helicity selection, for debugging
      INTEGER h
      COMMON /hcount / h

      udwza = 0d0
      udwza2 = 0d0

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

      do mu = 0,3             ! kinematics for WZA
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! mu+
         pbar(mu,7) = v(mu,5) ! A
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         if (bos.eq.321) then
            call sud_wza(pbar,fsign,udwza2)
         else if (bos.eq.421) then
            call sdu_wza(pbar,fsign,udwza2)
         endif
      endif
#endif

      call qqWZA(pbar,fsign,nlo,1,udwza,bos,ps_number)

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then

         if (bos.eq.321) then
            call TRIBOSwptowz_mg2_had(v(0,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa_mg2_had(v(0,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza_mg2_had(v(0,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz_mg2_had(v(0,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa_mg2_had(v(0,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza_mg2_had(v(0,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza_mg2_had(v(0,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza_mg2_had(v(0,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a

         call qqWZA(pbar,fsign,nlo,1,udwza2,bos,ps_number)

         if (bos.eq.321) then
            call TRIBOSwptowz_had(v(0,1),1,h,wpwz(0,1),0)  !W -> W+ Z -> ve e+ mu- mu+
            call TRIBOSwptowa_had(v(0,1),1,h,wpwa(0,1),0)  !W -> W+ A -> ve e+ a
            call wptowza_had(v(0,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve e+ mu- mu+ a
         else if (bos.eq.421) then
            call TRIBOSwmtowz_had(v(0,1),1,h,wpwz(0,1),0)  !W -> W- Z -> ve~ e- mu- mu+
            call TRIBOSwmtowa_had(v(0,1),1,h,wpwa(0,1),0)  !W -> W- A -> ve~ e- a
            call wmtowza_had(v(0,1),1,h,zwwa(0,1),0)    !W -> WZA -> ve~ e- mu- mu+ a
         endif
         call TRIBOSatoza_had(v(0,1),1,h,aza(0,1),0)   !A -> Z A -> mu- mu+ a
         call TRIBOSztoza_had(v(0,1),1,h,zaz(0,1),0)   !Z -> Z A -> mu- mu+ a

         print*,"1 - udwza/udwza2 =   ",abs(1-udwza/udwza2)
         read(*,*)
      endif
#endif

      if (ldebug) then
         print*,"udwza =   ",udwza
#ifdef WITH_MADGRAPH
         print*,"udwza2 =  ",udwza2
#endif
         read(*,*)
      endif

      NCmatrixelt(0,1)=udwza

      do if1=2,(nfl/2)*2,2
         if (bos.eq.321) then
            iflav(1)=( if1 - physToDiag(1) + 1 )*fsign(physToDiag(1))
            iflav(2)=( if1 - physToDiag(2) + 1 )*fsign(physToDiag(2))
         else if (bos.eq.421) then
            iflav(1)=( if1 + physToDiag(1) - 2 )*fsign(physToDiag(1))
            iflav(2)=( if1 + physToDiag(2) - 2 )*fsign(physToDiag(2))
         endif
         k=FL_VVV(iflav,1)

         res(k)=pdf(sign1*iflav(1),1)
     &         *pdf(sign2*iflav(2),2)
     &         *NCmatrixelt(mod(if1,2),mod(if1+1,2))*polcol

      enddo 

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WZA
c
c*****************************************************************

      subroutine vtowza_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            wpwz(mu,j)  = zero
            wpwa(mu,j)  = zero
            aza(mu,j)   = zero
            zaz(mu,j)   = zero
            zwwa(mu,j)  = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
