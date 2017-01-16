*****************************************************************************
c
c    begin function m2s_qqWWZ
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c       
c
c       this function computes |ME|^2 for qq-> ve e+ mu-vm~ ta- ta+ ("WWZ")
c	at LO and virtual corrections via call of "qqWWZ.f"
c
c	For the decay Z -> ta- ta+ random helicity summation is employed
c
c*****************************************************************************

      DOUBLE PRECISION FUNCTION m2s_qqWWZ(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      NW    !in:  = 1:  NW Approximation, = 0: Full calculation
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
      INTEGER bos, nlo, NW
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions
      DOUBLE PRECISION dotrr
      DOUBLE COMPLEX dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

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

c define program switches and debugging variables
      LOGICAL ldebug, ldebug_ltens
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebug_ltens=.false.)
      complex*16 mg2_tens(0:5)
      double precision sumtest

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c variables for hadronic decays
      integer N_gen_W
      real*8 fac_W, mjj2
      external mjj2
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down

c ------------- BEGIN CODE -------------------

      if ( init .eq. 0 ) then

         write(6,*) " "
         write(6,*) "W+W-Z amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowwz_reset

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
         qz(mu,1) = v(mu,5,1)+v(mu,6,1)   !Z
         qww(mu,1) = qp(mu,1) + qm(mu,1)
         qwwz(mu,1) = qww(mu,1) + qz(mu,1)
      enddo

      qp(4,1) = dotrr(qp(0,1),qp(0,1))
      qm(4,1) = dotrr(qm(0,1),qm(0,1))
      qz(4,1) = dotrr(qz(0,1),qz(0,1))
      qww(4,1) = dotrr(qww(0,1),qww(0,1))
      qwwz(4,1) = dotrr(qwwz(0,1),qwwz(0,1))

c for WWZ->6l precalculate A->WWZ,Z->WWZ etc..

c lepton spinors and W+, W-, Z polarisation vectors
c lepton helicities for Z -> ta- ta+ not fixed
c -> sum over all possible helicities in |M|**2
 
c select helicity: h ... random number for lepton helicity combination (h=1:2)
        h = mod(h,2) + 1
	it = (-1)**(h+1)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-1
           helicity(2)= 1
           helicity(3)=-1
           helicity(4)= 1
           helicity(5)=-it
           helicity(6)= it
        endif

      CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve(1,1)) !W(1,3))          !ve 
      CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep(1,1)) !W(1,4))          !e+        
      CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu(1,1)) !W(1,5))          !mu-      
      CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm(1,1)) !W(1,6))          !vm~
      CALL OXXXXX(v(0,5,1),LMASS ,-it,1,ztm(1,1)) !W(1,7))        !ta-      
      CALL IXXXXX(v(0,6,1),LMASS ,it,-1,ztp(1,1)) !W(1,8))        !ta+

      CALL JIOXXX(wep(1,1),wve(1,1),GWF,WMASS,WWIDTH, wp(1,1))  !W(1,9)  W+
      CALL JIOXXX(wvm(1,1),wmu(1,1),GWF,WMASS,WWIDTH, wm(1,1))  !W(1,10) W-

      SELECT CASE(procid)
      CASE(WWZhad)
            ! for general output to all combinations: up-type first
            if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
              ! set couplings for Z hadronic decay into uubar
              call setZtouu
            endif
            CALL JIOXXX(ztp(1,1),ztm(1,1),GZ_AF,ZERO ,ZERO  ,atau(1,1)) !W(1,11) A
            CALL JIOXXX(ztp(1,1),ztm(1,1),GZ_ZF,ZMASS,ZWIDTH,ztau(1,1)) !W(1,14) Z

      CASE DEFAULT
            CALL JIOXXX(ztp(1,1),ztm(1,1),GAL,ZERO ,ZERO  ,atau(1,1)) !W(1,11) A
            CALL JIOXXX(ztp(1,1),ztm(1,1),GZL,ZMASS,ZWIDTH,ztau(1,1)) !W(1,14) Z

      END SELECT

c************reduction of pentagon contribution *************
#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c wpp=wp-xp*qp, wmp=wm-xm*qm, ztaup=ztau-xz*qz, ataup=atau-xa*qz
c in eps~ * q_WW = 0
            qq = qp(0,1)*qww(0,1)-qp(1,1)*qww(1,1)-qp(2,1)*qww(2,1)
     &          -qp(3,1)*qww(3,1)
            xp = dotrc(qww(0,1),wp(1,1))/qq
            qq = qm(0,1)*qww(0,1)-qm(1,1)*qww(1,1)-qm(2,1)*qww(2,1)
     &          -qm(3,1)*qww(3,1)
            xm = dotrc(qww(0,1),wm(1,1))/qq
            qq = qz(0,1)*qww(0,1)-qz(1,1)*qww(1,1)-qz(2,1)*qww(2,1)
     &          -qz(3,1)*qww(3,1)
            xz = dotrc(qww(0,1),ztau(1,1))/qq
            xa = dotrc(qww(0,1),atau(1,1))/qq
            do mu = 0,3
               wpp(mu+1,1) = wp(mu+1,1) - xp*qp(mu,1)
               wmp(mu+1,1) = wm(mu+1,1) - xm*qm(mu,1)
               ztaup(mu+1,1) = ztau(mu+1,1) - xz*qz(mu,1)
               ataup(mu+1,1) = atau(mu+1,1) - xa*qz(mu,1)
            enddo
            do mu = 5,6
               wpp(mu,1) = wp(mu,1)
               wmp(mu,1) = wm(mu,1)
               ztaup(mu,1) = ztau(mu,1)
               ataup(mu,1) = atau(mu,1)
            enddo
         endif
#endif
c*****************************************************************

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

      SELECT CASE(procid)
      CASE(WPhadWMZ, WPWMhadZ, WWZhad)
            call anomal_formfactor(qwwz(0,1),qp(0,1),qm(0,1),qz(0,1))
            call TRIBOSatoww_had(v(0,1,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
            call TRIBOSztoww_had(v(0,1,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
            call wptowz_had(v(0,1,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
            call wmtowz_had(v(0,1,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
            call atowwz_had(v(0,1,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
            call ztowwz_had(v(0,1,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
      CASE DEFAULT
            if (with_anom) then ! anomalous gauge boson couplings
c             using global form factor for all tensors of one phase space point
c             this ensures proper cancellations for anomalous contributions
c             energy scale is invariant WWZ mass
              call anomal_formfactor(qwwz(0,1),qp(0,1),qm(0,1),qz(0,1))
              call TRIBOSatoww_anomal(v(0,1,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
              call TRIBOSztoww_anomal(v(0,1,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
              call wptowz_anomal(v(0,1,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
              call wmtowz_anomal(v(0,1,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
              call atowwz_anomal(v(0,1,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
              call ztowwz_anomal(v(0,1,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
            else if (with_kk) then ! Kaluza Klein Scenario
#ifdef WITH_KK
              call TRIBOSatoww(v(0,1,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
              call TRIBOSztoww(v(0,1,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
              call wptowz(v(0,1,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
              call wmtowz(v(0,1,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
              call atowwz_kk(v(0,1,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
              call ztowwz_kk(v(0,1,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
#endif
            else ! SM scenario
               call TRIBOSatoww(v(0,1,1),1,aww(0,1),NW) !A -> W+ W- -> ve e+ mu- vm~
               call TRIBOSztoww(v(0,1,1),1,zww(0,1),NW) !Z -> W+ W- -> ve e+ mu- vm~
               call wptowz(v(0,1,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
               call wmtowz(v(0,1,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
               call atowwz(v(0,1,1),1,h,awwz(0,1),NW) !A->WWZ -> ve e+ mu- vm~ ta- ta+
               call ztowwz(v(0,1,1),1,h,zwwz(0,1),NW) !Z->WWZ-> ve e+ mu- vm~ ta- ta+
            endif
      END SELECT

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         print*, ""
         call ATOWWZ_mg_had(v(0,1,1),1,h,mg2_tens,NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( awwz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to wwz:", sumtest
         call ZTOWWZ_mg_had(v(0,1,1),1,h,mg2_tens,NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zwwz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor z to wwz:", sumtest
         call TRIBOSatoww_mg_had(v(0,1,1),1,mg2_tens,NW)   !A -> W+ W- -> ve e+ mu- vm~
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( aww(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to ww:", sumtest
         call TRIBOSztoww_mg_had(v(0,1,1),1,mg2_tens,NW)   !Z -> W+ W- -> ve e+ mu- vm~
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zww(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor z to ww:", sumtest
         call wptowz_mg_had(v(0,1,1),1,h,mg2_tens,NW) !W+ -> W+ Z -> ve e+ ta- ta+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wpwz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor w+ to wz:", sumtest
         call wmtowz_mg_had(v(0,1,1),1,h,mg2_tens,NW) !W- -> W- Z -> mu- vm~ ta- ta+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( wmwz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor w- to wz:", sumtest
         read(*,*)
      endif
#endif

c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VVV(fsign,-1)    !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> W W Z   **************************

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
      fsign(8) = -1

      call qq_WWZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW)

C*******************  q1bar q2 ---> W W Z   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WWZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(WWZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        CALL JIOXXX(ztp(1,1),ztm(1,1),GZ_AF,ZERO ,ZERO  ,atau(1,1)) !W(1,11) A
        CALL JIOXXX(ztp(1,1),ztm(1,1),GZ_ZF,ZMASS,ZWIDTH,ztau(1,1)) !W(1,14) Z
        call TRIBOSatoww_had(v(0,1,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
        call TRIBOSztoww_had(v(0,1,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
        call wptowz_had(v(0,1,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
        call wmtowz_had(v(0,1,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
        call atowwz_had(v(0,1,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
        call ztowwz_had(v(0,1,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+

c************reduction of pentagon contribution *************
#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c wpp=wp-xp*qp, wmp=wm-xm*qm, ztaup=ztau-xz*qz, ataup=atau-xa*qz
c in eps~ * q_WW = 0
            qq = qp(0,1)*qww(0,1)-qp(1,1)*qww(1,1)-qp(2,1)*qww(2,1)
     &          -qp(3,1)*qww(3,1)
            xp = dotrc(qww(0,1),wp(1,1))/qq
            qq = qm(0,1)*qww(0,1)-qm(1,1)*qww(1,1)-qm(2,1)*qww(2,1)
     &          -qm(3,1)*qww(3,1)
            xm = dotrc(qww(0,1),wm(1,1))/qq
            qq = qz(0,1)*qww(0,1)-qz(1,1)*qww(1,1)-qz(2,1)*qww(2,1)
     &          -qz(3,1)*qww(3,1)
            xz = dotrc(qww(0,1),ztau(1,1))/qq
            xa = dotrc(qww(0,1),atau(1,1))/qq
            do mu = 0,3
               wpp(mu+1,1) = wp(mu+1,1) - xp*qp(mu,1)
               wmp(mu+1,1) = wm(mu+1,1) - xm*qm(mu,1)
               ztaup(mu+1,1) = ztau(mu+1,1) - xz*qz(mu,1)
               ataup(mu+1,1) = atau(mu+1,1) - xa*qz(mu,1)
            enddo
            do mu = 5,6
               wpp(mu,1) = wp(mu,1)
               wmp(mu,1) = wm(mu,1)
               ztaup(mu,1) = ztau(mu,1)
               ataup(mu,1) = atau(mu,1)
            enddo
         endif
#endif
c*****************************************************************


C*******************  q1 qbar2 ---> W W Z   **************************

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
      fsign(8) = -1

      call qq_WWZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW)

C*******************  q1bar q2 ---> W W Z   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_WWZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(WPhadWMZ)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPWMhadZ)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,3,1),v(0,4,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WWZhad)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,5,1),v(0,6,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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
         SELECT CASE(process)
         CASE(WPhadWMZ)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(WPWMhadZ)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(WWZhad)
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

      m2s = m2s * 2d0     ! factor 2 for random helicity summation

      m2s_qqWWZ = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqWWZ
c
c*****************************************************************************

c*****************************************************************************
c
c   begin subroutine qq_WWZ
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c       
c
c
c       qq_WWZ calls the amplitude square routines 
c             qqWWZ    for qq --> WWZ --> 6 leptons
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c          NW              NW = 1: NW approximation, NW = 0: full calc.
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
c*****************************************************************************

      SUBROUTINE qq_WWZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax, NW)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
      INTEGER NW

      DOUBLE PRECISION N                              ! color factors
      PARAMETER(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      DOUBLE PRECISION uuwwz, ddwwz, uuwwz2, ddwwz2
      DOUBLE PRECISION pbar(0:3,2+n_v), dotrr, polcol
      DOUBLE PRECISION NCmatrixelt(0:1,0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(2), diagToPhys(2), FL_VVV
      EXTERNAL dotrr, FL_VVV


      LOGICAL ldebug, ldebug_ltens
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebug_ltens=.false.)
      PARAMETER (zero = (0d0,0d0) )
c  helicity selection, for debugging
      INTEGER h
      COMMON /hcount / h
      INTEGER i,if1,k,mu

      uuwwz = 0d0
      uuwwz2 = 0d0
      ddwwz = 0d0
      ddwwz2 = 0d0

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

      do mu = 0,3             ! kinematics for WWZ --> 6 lepton decay
         pbar(mu,3) = v(mu,1) ! ve
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! vm~
         pbar(mu,7) = v(mu,5) ! tau-
         pbar(mu,8) = v(mu,6) ! tau+
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         call suu_wwz(pbar,fsign,uuwwz2,nw)
         call sdd_wwz(pbar,fsign,ddwwz2,nw)
      endif
#endif

#ifdef WITH_KK
      if (.not. mgtest) then
         call QQWWZ(pbar,fsign,nlo,1,uuwwz,ddwwz,NW)
      else
         call suu_wwz_mhlm(pbar,fsign,uuwwz)
         call sdd_wwz_mhlm(pbar,fsign,ddwwz)
      endif
#else
      call QQWWZ(pbar,fsign,nlo,1,uuwwz,ddwwz,NW)
#endif

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         call TRIBOSatoww_mg_had(v(0,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
         call TRIBOSztoww_mg_had(v(0,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
         call wptowz_mg_had(v(0,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
         call wmtowz_mg_had(v(0,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
         call atowwz_mg_had(v(0,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
         call ztowwz_mg_had(v(0,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
         call QQWWZ(pbar,fsign,nlo,1,uuwwz2,ddwwz2,NW)
         call TRIBOSatoww_had(v(0,1),1,aww(0,1),NW)   !A -> W+ W- -> ve e+ mu- vm~
         call TRIBOSztoww_had(v(0,1),1,zww(0,1),NW)   !Z -> W+ W- -> ve e+ mu- vm~
         call wptowz_had(v(0,1),1,h,wpwz(0,1),NW) !W+ -> W+ Z -> ve e+ ta- ta+
         call wmtowz_had(v(0,1),1,h,wmwz(0,1),NW) !W- -> W- Z -> mu- vm~ ta- ta+
         call atowwz_had(v(0,1),1,h,awwz(0,1),NW) !A -> WWZ -> ve e+ mu- vm~ ta- ta+
         call ztowwz_had(v(0,1),1,h,zwwz(0,1),NW) !Z -> WWZ -> ve e+ mu- vm~ ta- ta+
         print*,"1 - uuwwz/uuwwz2 =   ",abs(1-uuwwz/uuwwz2)
         print*,"1 - ddwwz/ddwwz2 =   ",abs(1-ddwwz/ddwwz2)
         read(*,*)
      endif
      if (ldebug) then
         print*,"uuwwz =   ",uuwwz
         print*,"uuwwz2 =  ",uuwwz2
         print*,"ddwwz =   ",ddwwz
         print*,"ddwwz2 =  ",ddwwz2
         read(*,*)
      endif
#endif

      NCmatrixelt(0,0)=uuwwz
      NCmatrixelt(1,1)=ddwwz

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
c   end subroutine qq_WWZ
c
c*****************************************************************

c******************************************************************
c
c   begin function FL_VVV
c
c*****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c       
c  this function fills hepup color and flavor assignments for
c  VVV processes with two incoming quarks only
c
c  id1,2 are flavor identifiers for incoming quarks
c  colstruc = is irrelevant here
c
c  color structure of these interactions is trivial
c  and colstruc was included only to make this program 
c  more consistent with the others.
c
c  Note that colstruc is also a flag:
c  if colstruc=-1, then we are resetting.
c
c  The color structure is always the same:
c  color always flows from id1->id2
c  no gluons are involved.
c
c********************************************************************
 
      INTEGER function FL_VVV(iflav,colstruc)

      IMPLICIT NONE

      INTEGER iflav(2),colstruc    ! input for color and flavor assignment

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
      
      INTEGER lkup
      COMMON /localkup/ lkup(numParticles,7,maxNumSubProc)
 
      INTEGER listposition
      INTEGER numresets
      SAVE listposition
      SAVE numresets

      DATA numresets /-1/
      DATA listposition /0/

      if (colstruc.eq.-1) then  ! we are resetting
c       write(*,*) "we are resetting"
        listposition=0
        numresets=numresets+1
      else if (colstruc.eq.-2) then
c       don"t do anything, just return the number of subprocesses.
      else
c increment the counter regardless of whether or not it's the 1st time through.
        listposition=listposition+1

c       we fill info for this subprocess,if required
        if(numresets.eq.0 .and. (lha.or.hepmc)) then
          call fillColoredPartons_VVV(iflav,
     1                               listposition+numdecay)
        endif       ! numresets.eq.0
      endif         
      FL_VVV=listposition
      end
c******************************************************************************
c
c   end function FL_VVV
c
c******************************************************************************

c*****************************************************************************
c
c    begin  subroutine fillColoredPartons_VVV
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c
c  assigns values to the variables in the common block localHEPUP
c  in particular, this subroutine assigns values to those variables that 
c  will be stored in the lookup tables generated by writeHEPUPtable.  
c  As the name suggests, this routine only stores the information for the 
c  colored partons.  Particles without color will be dealt with in the 
c  subroutine fillColorless.
c
c*****************************************************************************

      SUBROUTINE fillColoredPartons_VVV(iflav,listposition)

      IMPLICIT NONE

c  include files
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/process.inc"

c  input parameters
      INTEGER iflav(2),id1,id2,listposition

c  local variables

      
      id1=iflav(1)
      id2=iflav(2)
   
      select case(process)
      case(WWZ,ZZWP,ZZWM,ZZZ,WWWM,WWWP,
     &     WPhadWMZ, WPWMhadZ, WWZhad, 
     &     ZZhadWP, ZZhadWM, ZZWPhad, ZZWMhad,
     &     WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM,
     &     ZZZhad) 
         lnup(listposition)=numParticles+3
      case (WWA, ZZA, WPZA, WMZA,
     &      WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
         lnup(listposition) = numParticles+1
      case (WPAA, WMAA, ZAA, ZnAA)
         lnup(listposition) = numParticles-1
      case (AAA)
         lnup(listposition) = numParticles-3
      end select

      listup(1,listposition)=-1   !incoming quarks
      listup(2,listposition)=-1
 
      lidup(1,listposition)=id1     ! the 2 quarks
      lidup(2,listposition)=id2

      if (id1.gt.0) then
        licolup(1,1,listposition)=501
        licolup(2,1,listposition)=0
        licolup(1,2,listposition)=0
        licolup(2,2,listposition)=501
      else
        licolup(1,1,listposition)=0
        licolup(2,1,listposition)=501
        licolup(1,2,listposition)=501
        licolup(2,2,listposition)=0
      endif

      end  ! fillColoredPartons_VVV
c*****************************************************************************
c
c    end subroutine fillColoredPartons_VVV
c
c*****************************************************************************

c*****************************************************************************
c
c    begin subroutine vtowwz_reset
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2006 March 23
c	Last modified:  2009 January 28
c       
c  this subroutines sets all leptonic tensors needed in the calculation 
c  of pp -> WWZ to zero.
c  
c*****************************************************************************

      SUBROUTINE vtowwz_reset

      IMPLICIT NONE

#include "tensor.inc"

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))
      INTEGER j,mu

      do mu = 0,5
         do j = 1,3
            aww(mu,j)  = zero
            zww(mu,j)  = zero
            awwz(mu,j) = zero
            zwwz(mu,j) = zero
            wpwz(mu,j) = zero
            wmwz(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
c*****************************************************************************
c
c    end subroutine vtowwz_reset
c
c*****************************************************************************
