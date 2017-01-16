*****************************************************************************
c
c    begin function m2s_qqZZZ
c
c*****************************************************************************

c	Vera Hankele, <vera@particle.physik.uni-karlsruhe.de>
c	Initial version:  September 2007
c	Last modified: June 2009
c
c       this function computes |ME|^2 for qq -> ZZZ
c	at LO and virtual corrections via call of "qqZZZ.F" and "qqZZZ1.F"
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqZZZ(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,    !in:  random number array
     $      NW,     !in:  = 1:  NW Approximation, = 0: Full calculation
     &      ps_number                  )
      IMPLICIT NONE

c declare global variables

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/process.inc"

c for hadronic decay
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2
      complex*16 mg2_tens(0:5), mg2_tens2(0:5), mg2_tens3(0:5)
      double precision sumtest

c declare input/output variables

      INTEGER bos, nlo, NW,ps_number
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions

      DOUBLE PRECISION dotrr
      INTEGER FL_VVV
      EXTERNAL dotrr
      EXTERNAL FL_VVV

c helicity selection

      INTEGER h
      COMMON /hcount/ h

c declare local variables

      COMPLEX*16 zero
      PARAMETER (zero=(0d0,0d0))

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
      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

c variables for leshouches interface
      DOUBLE PRECISION weight,rnumb,RandomNumber

c ------------- BEGIN CODE -------------------

      if ( init .eq. 0 ) then

         write(6,*) " "
         write(6,*) "ZZZ amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtozzz_reset

         call printnfl(.false.)
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

c lepton spinors and Z, Z, Z polarisation vectors
c lepton helicities for Z -> l- l+ not fixed
c -> sum over all possible helicities in |M|**2

c select helicity: h ... random number for lepton helicity
c                        combination (h=1:8)

        h = mod(h,8) + 1
        ie = sign(1,4-h)
        iu = (-1)**(h+1)
        it = (-1)**((h-1)/2)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5)=-it
           helicity(6)= it
        endif

      CALL OXXXXX(v(0,1,1),ZERO ,-ie,+1,zem(1,1)) !W(1,3))       !e-
      CALL IXXXXX(v(0,2,1),ZERO ,+ie,-1,zep(1,1)) !W(1,4))       !e+
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,+1,zmm(1,1)) !W(1,5))       !mu-
      CALL IXXXXX(v(0,4,1),ZERO ,+iu,-1,zmp(1,1)) !W(1,6))       !mu+
      CALL OXXXXX(v(0,5,1),ZERO ,-it,+1,ztm(1,1)) !W(1,7))       !ta-
      CALL IXXXXX(v(0,6,1),ZERO ,+it,-1,ztp(1,1)) !W(1,8))       !ta+

      SELECT CASE(procid)
      CASE(ZZZhad)
         ! for general output to all combinations: up-type first
         if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
           ! set couplings for Z hadronic decay into uubar
           call setZtouu
         endif
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO , ZERO ,ae(1,1))   !W(1,9)  A
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH,ze(1,1))   !W(1,9)  Z

      CASE DEFAULT
        CALL JIOXXX(zep(1,1),zem(1,1),GAL,ZERO , ZERO ,ae(1,1))   !W(1,9)  A
        CALL JIOXXX(zep(1,1),zem(1,1),GZL,ZMASS,ZWIDTH,ze(1,1))   !W(1,9)  Z
      END SELECT

      CALL JIOXXX(zmp(1,1),zmm(1,1),GAL,ZERO ,ZERO  ,amu(1,1))  !W(1,10) A
      CALL JIOXXX(zmp(1,1),zmm(1,1),GZL,ZMASS,ZWIDTH,zmu(1,1))  !W(1,10) Z
      CALL JIOXXX(ztp(1,1),ztm(1,1),GAL,ZERO ,ZERO  ,atau(1,1)) !W(1,11) A
      CALL JIOXXX(ztp(1,1),ztm(1,1),GZL,ZMASS,ZWIDTH,ztau(1,1)) !W(1,14) Z

      do mu = 0,3
         qze(mu,1) = v(mu,1,1)+v(mu,2,1)   !Z
         qzmu(mu,1) = v(mu,3,1)+v(mu,4,1)   !Z
         qztau(mu,1) = v(mu,5,1)+v(mu,6,1)   !Z
         qzzz(mu,1) = qze(mu,1) + qzmu(mu,1) + qztau(mu,1)
      enddo

      qze(4,1) = dotrr(qze(0,1),qze(0,1))
      qzmu(4,1) = dotrr(qzmu(0,1),qzmu(0,1))
      qztau(4,1) = dotrr(qztau(0,1),qztau(0,1))
      qzzz(4,1) = dotrr(qzzz(0,1),qzzz(0,1))


      SELECT CASE(procid)
      CASE(ZZZhad)
            call anomal_formfactor(qzzz(0,1),qze(0,1),qzmu(0,1),qztau(0,1))
c leptonic tensors
c all diagrams
            if (ps_number.eq.1) then
               call atozzz_had(v(0,1,1),1,h,azzz(0,1),NW)   !A -> ZZZ -> 6 leptons
               call ztozzz_had(v(0,1,1),1,h,zzzz(0,1),NW)   !Z -> ZZZ -> 6 leptons
               call zzz_atozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
               call zzz_atozz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,azezt(0,1))
               call atozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,azmzt(0,1))  ! 2nd and 3rd Z: no had decay
               call zzz_ztozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
               call zzz_ztozz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,zzezt(0,1))
               call ztozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,zzmzt(0,1))  ! 2nd and 3rd Z: no had decay
            endif
c only higgs contributions
            call ztozzzHiggs_had(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations

      CASE DEFAULT
      if (with_anom) then ! anomalous gauge boson couplings
c       using global form factor for all tensors of one phase space point
c       this ensures proper cancellations for anomalous contributions
c       energy scale is invariant ZZZ mass
         call anomal_formfactor(qzzz(0,1),qze(0,1),qzmu(0,1),qztau(0,1))

c leptonic tensors
         if (NW.eq.1) then
c NW approximation (only Higgs and neutral quartic-couplings contribution exists in the leptonic tensors)
            call ztozzz_anomal(v(0,1,1),1,h,zzzz(0,1),NW) !Z -> ZZZ -> 6 leptons
            call atozzz_anomal(v(0,1,1),1,h,azzz(0,1),NW) !A -> ZZZ -> 6 leptons
c only higgs contributions
            call ztozzzHiggs_anomal(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         else
c all diagrams
            if (ps_number.eq.1) then
               call atozzz_anomal(v(0,1,1),1,h,azzz(0,1),NW)   !A -> ZZZ -> 6 leptons
               call ztozzz_anomal(v(0,1,1),1,h,zzzz(0,1),NW)   !Z -> ZZZ -> 6 leptons
               call atozz_anomal(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
               call atozz_anomal(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,azezt(0,1))
               call atozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,azmzt(0,1))
               call ztozz_anomal(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
               call ztozz_anomal(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,zzezt(0,1))
               call ztozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,zzmzt(0,1))
            endif
c only higgs contributions
            call ztozzzHiggs_anomal(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         endif

      else ! SM

c leptonic tensors
         if (NW.eq.1) then
c NW approximation (only Higgs contribution exists in the leptonic tensors)
            call ztozzz(v(0,1,1),1,h,zzzz(0,1),NW) !Z -> ZZZ -> 6 leptons
c only higgs contributions
            call ztozzzHiggs(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         else
c all diagrams
            if (ps_number.eq.1) then
               call atozzz(v(0,1,1),1,h,azzz(0,1))   !A -> ZZZ -> 6 leptons
               call ztozzz(v(0,1,1),1,h,zzzz(0,1),NW) !Z -> ZZZ -> 6 leptons
               call atozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
               call atozz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,azezt(0,1))
               call atozz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,azmzt(0,1))
               call ztozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
               call ztozz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,zzezt(0,1))
               call ztozz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,zzmzt(0,1))
            endif
c only higgs contributions
            call ztozzzHiggs(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         endif

      endif !anom couplings or SM
      END SELECT ! procid


#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         print*, ""
         print*, "after ltens-definition, mgtest vs original:"
         if (ps_number.eq.1) then
            call atozzz_mg2_had(v(0,1,1),1,h,mg2_tens,NW)   !A -> ZZZ -> 6 leptons
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( azzz(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor a to zzz:", sumtest
            call ztozzz_mg2_had(v(0,1,1),1,h,mg2_tens,NW)   !Z -> ZZZ -> 6 leptons
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzzz(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor z to zzz:", sumtest
            call zzz_atozz_mg2_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,mg2_tens)
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( azezm(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor a to zz 1:", sumtest
            call zzz_atozz_mg2_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,mg2_tens)
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( azezt(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor a to zz 2:", sumtest
            call zzz_ztozz_mg2_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,mg2_tens)
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzezm(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor z to zz 1:", sumtest
            call zzz_ztozz_mg2_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,mg2_tens)
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzezt(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor z to zz 2:", sumtest
         endif
c only higgs contributions
         call ztozzzHiggs_had(v(0,1,1),1,h,mg2_tens,mg2_tens2,mg2_tens3,NW) !Z -> H Z, H -> Z Z for all permutations
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzzz_higgs1(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor z to zzz via higgs 1:", sumtest
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzzz_higgs2(i,1) / mg2_tens2(i) - 1)
            enddo
            print*, "tensor z to zzz via higgs 2:", sumtest
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( zzzz_higgs3(i,1) / mg2_tens3(i) - 1)
            enddo
            print*, "tensor z to zzz via higgs 3:", sumtest
         read(*,*)
      endif
#endif


c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1))        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2))        ! f_a(x2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Now sum over the subprocesses contributing to Z production

      nmax = FL_VVV(fsign,-1)    !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z Z   **************************

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

      call qq_ZZZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW,ps_number)

C*******************  q1bar q2 ---> W W Z   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW,ps_number)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZZhad)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO , ZERO ,ae(1,1))   !W(1,9)  A
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH,ze(1,1))   !W(1,9)  Z

            if (ps_number.eq.1) then
               call atozzz_had(v(0,1,1),1,h,azzz(0,1),NW)   !A -> ZZZ -> 6 leptons
               call ztozzz_had(v(0,1,1),1,h,zzzz(0,1),NW)   !Z -> ZZZ -> 6 leptons
               call zzz_atozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
               call zzz_atozz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,azezt(0,1))
               call atozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,azmzt(0,1))  ! 2nd and 3rd Z: no had decay
               call zzz_ztozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
               call zzz_ztozz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),-ie,-it,zzezt(0,1))
               call ztozz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),-iu,-it,zzmzt(0,1))  ! 2nd and 3rd Z: no had decay
            endif
c only higgs contributions
            call ztozzzHiggs_had(v(0,1,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations


C*******************  q1 qbar2 ---> Z Z Z   **************************

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

      call qq_ZZZ(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW,ps_number)

C*******************  q1bar q2 ---> W W Z   **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW,ps_number)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(ZZZhad)
        call hadDecayFactor_Z(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)
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
         CASE(ZZZhad)
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

      m2s = m2s * 8.0d0 !factor for random helicity summation

      m2s_qqZZZ = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqZZZ
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZZZ
c
c*****************************************************************

      subroutine qq_ZZZ(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax, NW,ps_number)

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      INTEGER physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax, NW
      Integer ps_number
c
c wbf_zh calls the amplitude square routines 
c             qqVVV    for qq --> VVV --> 6 leptons
c
c  INPUT:  p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c          NW              NW = 1: NW approximation, NW = 0: full calc.
c
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
      DOUBLE PRECISION N                              ! color factors
      PARAMETER(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "tensor.inc"

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

c  helicity selection, for debugging
      INTEGER h
      COMMON /hcount / h

c local variables
      DOUBLE PRECISION uuzzz, ddzzz, uuzzz2, ddzzz2
      DOUBLE PRECISION uuzzz_mg, ddzzz_mg, uuzzz2_mg, ddzzz2_mg
      DOUBLE PRECISION pbar(0:3,2+n_v), dotrr, polcol
      DOUBLE PRECISION NCmatrixelt(0:1,0:1)
      DOUBLE COMPLEX zero
      INTEGER iflav(4), diagToPhys(2), FL_VVV
      EXTERNAL dotrr, FL_VVV


      LOGICAL ldebug
      PARAMETER (ldebug=.false.)
      PARAMETER (zero = (0d0,0d0) )
      INTEGER i,if1,k,mu

      uuzzz = 0d0
      uuzzz2 = 0d0
      ddzzz = 0d0
      ddzzz2 = 0d0

      uuzzz_mg = 0d0
      uuzzz2_mg = 0d0
      ddzzz_mg = 0d0
      ddzzz2_mg = 0d0

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
         pbar(mu,3) = v(mu,1) ! e-
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! mu+
         pbar(mu,7) = v(mu,5) ! ta-
         pbar(mu,8) = v(mu,6) ! ta+
      enddo

c call of the subroutine to calculate matrix element squared
c last entry = 1: only qqbar -> ZZZ; = 0 all diagrams

      if (ps_number.eq.1) then ! subtract Higgs contr.
         call qqZZZ(pbar,fsign,nlo,1,uuzzz,ddzzz,NW)
         call qqZZZ1(pbar,fsign,nlo,1,uuzzz2,ddzzz2,NW)
         uuzzz = uuzzz - uuzzz2
         ddzzz = ddzzz - ddzzz2
      else if (ps_number.eq.2) then ! only Higgs contr.
         call qqZZZ1(pbar,fsign,nlo,1,uuzzz,ddzzz,NW)
      else ! complete |M|^2
         call qqZZZ(pbar,fsign,nlo,1,uuzzz,ddzzz,NW)
      endif

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
        ie = sign(1,4-h)
        iu = (-1)**(h+1)
        it = (-1)**((h-1)/2)
         if (ps_number.eq.1) then
            call atozzz_mg2_had(v(0,1),1,h,azzz(0,1),NW)   !A -> ZZZ -> 6 leptons
            call ztozzz_mg2_had(v(0,1),1,h,zzzz(0,1),NW)   !Z -> ZZZ -> 6 leptons
            call zzz_atozz_mg2_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
            call zzz_atozz_mg2_had(v(0,1),v(0,2),v(0,5),v(0,6),-ie,-it,azezt(0,1))
            call zzz_ztozz_mg2_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
            call zzz_ztozz_mg2_had(v(0,1),v(0,2),v(0,5),v(0,6),-ie,-it,zzezt(0,1))
         endif
         call ztozzzHiggs_mg2_had(v(0,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         if (ps_number.eq.1) then ! subtract Higgs contr.
            call qqZZZ(pbar,fsign,nlo,1,uuzzz_mg,ddzzz_mg,NW)
            call qqZZZ1(pbar,fsign,nlo,1,uuzzz2_mg,ddzzz2_mg,NW)
            uuzzz_mg = uuzzz_mg - uuzzz2_mg
            ddzzz_mg = ddzzz_mg - ddzzz2_mg
         else if (ps_number.eq.2) then ! only Higgs contr.
            call qqZZZ1(pbar,fsign,nlo,1,uuzzz_mg,ddzzz_mg,NW)
         else ! complete |M|^2
            call qqZZZ(pbar,fsign,nlo,1,uuzzz_mg,ddzzz_mg,NW)
         endif
         if (ps_number.eq.1) then
            call atozzz_had(v(0,1),1,h,azzz(0,1),NW)   !A -> ZZZ -> 6 leptons
            call ztozzz_had(v(0,1),1,h,zzzz(0,1),NW)   !Z -> ZZZ -> 6 leptons
            call zzz_atozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
            call zzz_atozz_had(v(0,1),v(0,2),v(0,5),v(0,6),-ie,-it,azezt(0,1))
            call zzz_ztozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
            call zzz_ztozz_had(v(0,1),v(0,2),v(0,5),v(0,6),-ie,-it,zzezt(0,1))
         endif
         call ztozzzHiggs_had(v(0,1),1,h,zzzz_higgs1(0,1),zzzz_higgs2(0,1), 
     1              zzzz_higgs3(0,1),NW) !Z -> H Z, H -> Z Z for all permutations
         print*,"uuzzz    =  ",uuzzz
         print*,"uuzzz_mg =  ",uuzzz_mg, "  ratio:", abs(uuzzz/uuzzz_mg - 1)
         print*,"ddzzz    =  ",ddzzz
         print*,"ddzzz_mg =  ",ddzzz_mg, "  ratio:", abs(ddzzz/ddzzz_mg - 1)
         print*,"uuzzz2    = ",uuzzz2
         print*,"uuzzz2_mg = ",uuzzz2_mg, "  ratio:", abs(uuzzz2/uuzzz2_mg - 1)
         print*,"ddzzz2    = ",ddzzz2
         print*,"ddzzz2_mg = ",ddzzz2_mg, "  ratio:", abs(ddzzz2/ddzzz2_mg - 1)
         read(*,*)
      endif
      if (ldebug) then
         call suu_zzz(pbar,fsign,uuzzz2,NW)
         call sdd_zzz(pbar,fsign,ddzzz2,NW)
         print*,"uuzzz =   ",uuzzz
         print*,"uuzzz2 =  ",uuzzz2
         print*,"uuzzz2/uuzzz - 1 = ",uuzzz2/uuzzz - 1
         print*,"ddzzz =   ",ddzzz
         print*,"ddzzz2 =  ",ddzzz2
         print*,"ddzzz/ddzzz2 - 1 = ",ddzzz/ddzzz2 - 1
      endif
#endif

      NCmatrixelt(0,0)=uuzzz
      NCmatrixelt(1,1)=ddzzz

      do if1=1,nfl
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
c   end subroutine qq_ZZZ
c
c*****************************************************************

      subroutine vtozzz_reset
      implicit none
#include "tensor.inc"
      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            azzz(mu,j)  = zero
            zzzz(mu,j)  = zero
            azezm(mu,j) = zero
            azezt(mu,j) = zero
            azmzt(mu,j) = zero
            zzezm(mu,j) = zero
            zzezt(mu,j) = zero
            zzmzt(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
