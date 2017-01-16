*****************************************************************************
c
c    begin function m2s_qqZZA
c
c*****************************************************************************
c
c       giuseppe bozzi, <giuseppe@particle.physik.uni-karlsruhe.de>
c       Adapted from m2s_qqZZZ by V. Hankele
c       Last modified December 2008 (giuseppe)
c
c       this function computes |ME|^2 for qq-> e- e+ mu- mu+ A ("ZZA")
c	at LO and virtual corrections via call of 'qqZZA.f'
c
c
c       no reduction of pentagon contribution as in the other tri-boson
c       codes is performed because it drops out in the end due to gauge
c       invariance.
c
c ------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION m2s_qqZZA(
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
      integer N_gen_up, N_gen_down
      real*8 fac_Z_up, fac_Z_down, mjj2
      external mjj2

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)
      complex*16 mg2_tens(0:5)
      double precision sumtest

c declare input/output variables
      INTEGER bos, nlo
      DOUBLE PRECISION xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s

c declare external functions
      DOUBLE PRECISION dotrr
      double complex dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

c  helicity selection

      INTEGER h
      COMMON /hcount / h

      DOUBLE COMPLEX zero
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
         write(6,*) "ZZA amplitude square information:"
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtozza_reset

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
      
      do mu = 0,3
         qze(mu,1) = v(mu,1,1)+v(mu,2,1)  ! A,Z(e)
         qzmu(mu,1) = v(mu,3,1)+v(mu,4,1) ! A,Z(m)
         qz(mu,1) = v(mu,5,1)             ! A
         qzz(mu,1) = qze(mu,1) + qzmu(mu,1)             !ZZ 4-momentum
         qzza(mu,1) = qze(mu,1) + qzmu(mu,1) + qz(mu,1) !ZZA 4-momentum
      enddo
      
      qze(4,1) = dotrr(qze(0,1),qze(0,1))
      qzmu(4,1) = dotrr(qzmu(0,1),qzmu(0,1))
      qz(4,1) = dotrr(qz(0,1),qz(0,1))
      qzz(4,1) = dotrr(qzz(0,1),qzz(0,1))
      qzza(4,1) = dotrr(qzza(0,1),qzza(0,1))
      qwwz(4,1) = qzza(4,1)
      
c     for ZZA->4l precalculate A->ZZA,Z->ZZA etc..
      
c     lepton spinors and Z,A polarisation vectors
      
c     select helicity: h = random number for lepton+photon helicity (h=1:8)
      
      h = mod(h,8) + 1

      ie = sign(1,4-h)
      if(bos .eq. 221) then
         iu = (-1)**(h+1) !muons
      else if(bos .eq. 251) then
         iu = 1 ! only left-handed neutrinos
      end if
      it = (-1)**((h-1)/2)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5) =  it
        endif
      
      CALL OXXXXX(v(0,1,1),ZERO ,-ie,+1,zem(1,1)) !e- 
      CALL IXXXXX(v(0,2,1),ZERO ,+ie,-1,zep(1,1)) !e+        
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,+1,zmm(1,1)) !mu- or v      
      CALL IXXXXX(v(0,4,1),ZERO ,+iu,-1,zmp(1,1)) !mu+ or v~
      CALL VXXXXX(qz(0,1),ZERO ,it,1,atau(1,1)) !A

      SELECT CASE(procid)
      CASE(ZZhadA)
        ! for general output to all combinations: up-type first
        if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
          ! set couplings for Z hadronic decay into uubar
          call setZtouu
        endif
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO,ZERO, ae(1,1))    !A*e
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH, ze(1,1)) !Z*e

      CASE DEFAULT
        CALL JIOXXX(zep(1,1),zem(1,1),GAL,ZERO,ZERO, ae(1,1))    !A*e
        CALL JIOXXX(zep(1,1),zem(1,1),GZL,ZMASS,ZWIDTH, ze(1,1)) !Z*e
      END SELECT

      if(bos .eq. 221) then
         CALL JIOXXX(zmp(1,1),zmm(1,1),GAL,ZERO,ZERO, amu(1,1))   !A*m
         CALL JIOXXX(zmp(1,1),zmm(1,1),GZL,ZMASS,ZWIDTH, zmu(1,1))!Z*m
      else if(bos .eq. 251) then
         CALL JIOXXX(zmp(1,1),zmm(1,1),GZN,ZMASS,ZWIDTH, zmu(1,1))!Z*v
         amu(1:4,1) = 0
         amu(5:6,1) = zmu(5:6,1)
      endif
c************ reduction of pentagon contribution *************
#ifdef WITH_NLO
      if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c     reduction to "true" pentagons (see paper)
         qq = qze(0,1)*qzz(0,1) - qze(1,1)*qzz(1,1)-
     &        qze(2,1)*qzz(2,1) - qze(3,1)*qzz(3,1)
         xze = dotrc(qzz(0,1),ze(1,1))/qq
         xae = dotrc(qzz(0,1),ae(1,1))/qq
         qq = qzmu(0,1)*qzz(0,1) - qzmu(1,1)*qzz(1,1)-
     &        qzmu(2,1)*qzz(2,1) - qzmu(3,1)*qzz(3,1)
         xzmu = dotrc(qzz(0,1),zmu(1,1))/qq
         xamu = dotrc(qzz(0,1),amu(1,1))/qq
         qq = qz(0,1)*qzz(0,1) - qz(1,1)*qzz(1,1)-
     &        qz(2,1)*qzz(2,1) - qz(3,1)*qzz(3,1)
         xaon = dotrc(qzz(0,1),atau(1,1))/qq

cc----------------------------------------------------------
         xze = (0.0d0,0.0d0)
         xae = (0.0d0,0.0d0)
         xzmu = (0.0d0,0.0d0)
         xamu = (0.0d0,0.0d0)
         xaon = (0.0d0,0.0d0)
c------------------------------------------------------------

         do mu = 0,3
            zepp(mu+1,1) = ze(mu+1,1) - xze*qze(mu,1)
            zmup(mu+1,1) = zmu(mu+1,1) - xzmu*qzmu(mu,1)
            aep(mu+1,1) = ae(mu+1,1) - xae*qze(mu,1)
            amup(mu+1,1) = amu(mu+1,1) - xamu*qzmu(mu,1)
            ataup(mu+1,1) = atau(mu+1,1) - xaon*qz(mu,1)
         enddo
         do mu = 5,6
            zepp(mu,1) = ze(mu,1)
            zmup(mu,1) = zmu(mu,1)
            aep(mu,1) = ae(mu,1)
            amup(mu,1) = amu(mu,1)
            ataup(mu,1) = atau(mu,1)
         enddo         

      endif
#endif
c*****************************************************************


c leptonic tensors:
      SELECT CASE(procid)
      CASE(ZZhadA)
         call anomal_formfactor(qzza(0,1),qze(0,1),qzmu(0,1),qz(0,1))
         call zza_atozz_had(v(0,1,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+ 
         call zza_ztozz_had(v(0,1,1),1,zzztens(0,1),0)      !Z -> Z Z -> e- e+ mu- mu+
         call zza_atozae_had(v(0,1,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
         call zza_ztozae_had(v(0,1,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
         call zza_atozam_anomal(v(0,1,1),1,h,aazm(0,1),0)  !A -> Z A -> mu- mu+ a
         call zza_ztozam_anomal(v(0,1,1),1,h,zzam(0,1),0)  !Z -> Z A -> mu- mu+ a
         call zza_atozza_had(v(0,1,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
         call zza_ztozza_had(v(0,1,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a

      CASE DEFAULT

      if (with_anom) then ! anomalous gauge boson couplings
c        using global form factor for all tensors of one phase space point
c        this ensures proper cancellations for anomalous contributions
c        energy scale is invariant ZZA mass
         call anomal_formfactor(qzza(0,1),qze(0,1),qzmu(0,1),qz(0,1))

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)
         call zza_atozz_anomal(v(0,1,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+
         call zza_ztozz_anomal(v(0,1,1),1,zzztens(0,1),0)      !Z -> Z Z -> e- e+ mu- mu+
         call zza_atozae_anomal(v(0,1,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
         call zza_ztozae_anomal(v(0,1,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
         call zza_atozam_anomal(v(0,1,1),1,h,aazm(0,1),0)  !A -> Z A -> mu- mu+ a
         call zza_ztozam_anomal(v(0,1,1),1,h,zzam(0,1),0)  !Z -> Z A -> mu- mu+ a
         call zza_atozza_anomal(v(0,1,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
         call zza_ztozza_anomal(v(0,1,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a


      else  ! SM

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)
         if(bos .eq. 221) then !mu- mu+
            call zza_atozz(v(0,1,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+
            call zza_ztozz(v(0,1,1),1,zzztens(0,1),0)  !Z -> Z Z -> e- e+ mu- mu+
            call zza_atozae(v(0,1,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
            call zza_ztozae(v(0,1,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
            call zza_atozam(v(0,1,1),1,h,aazm(0,1),0)  !A -> Z A -> mu- mu+ a
            call zza_ztozam(v(0,1,1),1,h,zzam(0,1),0)  !Z -> Z A -> mu- mu+ a
            call zza_atozza(v(0,1,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
            call zza_ztozza(v(0,1,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a
         else if(bos .eq. 251) then !v v~
            call zzna_atozz(v(0,1,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ v v~
            call zzna_ztozz(v(0,1,1),1,zzztens(0,1),0)  !Z -> Z Z -> e- e+ v v~
            call zza_atozae(v(0,1,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
            call zza_ztozae(v(0,1,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
            aazm(0:3,1) = 0
            zzam(0:3,1) = 0
            call zzna_atozza(v(0,1,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ v v~ a
            call zzna_ztozza(v(0,1,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ v v~ a
         endif

      endif ! anom or SM

      end select


#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         print*, ""
         call zza_atozz_mg2_had(v(0,1,1),1,mg2_tens,0)      !A -> Z Z -> e- e+ mu- mu+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( azz(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to zz:", sumtest
         call zza_ztozz_mg2_had(v(0,1,1),1,mg2_tens,0)      !Z -> Z Z -> e- e+ mu- mu+
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zzztens(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to zz:", sumtest
         call zza_atozae_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !A -> Z A -> e- e+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( aaze(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to za:", sumtest
         call zza_ztozae_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !Z -> Z A -> e- e+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zzae(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor z to za:", sumtest
         call zza_atozza_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !A -> ZZA -> e- e+ mu- mu+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( azza(i,1) / mg2_tens(i) - 1)
           print*, azza(i,1) , mg2_tens(i)
         enddo
         print*, "tensor a to zza:", sumtest
         call zza_ztozza_mg2_had(v(0,1,1),1,h,mg2_tens,0)  !Z -> ZZA -> e- e+ mu- mu+ a
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zzza(i,1) / mg2_tens(i) - 1)
           print*, zzza(i,1) , mg2_tens(i)
         enddo
         print*, "tensor z to zza:", sumtest
         read(*,*)
      endif
#endif


c     distributions in the incoming (anti)protons.
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) ) ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) ) ! f_a(x2)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
C     Now sum over the subprocesses contributing to Z production
      
      nmax = FL_VVV(fsign,-1) !reset counter for subprocesses to 0

C*******************  q1 qbar2 ---> Z Z A   **************************

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

      call qq_ZZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)

C*******************  q1bar q2 ---> Z Z A   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZhadA)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd

        CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO,ZERO, ae(1,1))    !A*e
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH, ze(1,1)) !Z*e

c************ reduction of pentagon contribution *************
#ifdef WITH_NLO
      if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c     reduction to "true" pentagons (see paper)
         qq = qze(0,1)*qzz(0,1) - qze(1,1)*qzz(1,1)-
     &        qze(2,1)*qzz(2,1) - qze(3,1)*qzz(3,1)
         xze = dotrc(qzz(0,1),ze(1,1))/qq
         xae = dotrc(qzz(0,1),ae(1,1))/qq
         qq = qzmu(0,1)*qzz(0,1) - qzmu(1,1)*qzz(1,1)-
     &        qzmu(2,1)*qzz(2,1) - qzmu(3,1)*qzz(3,1)
         xzmu = dotrc(qzz(0,1),zmu(1,1))/qq
         xamu = dotrc(qzz(0,1),amu(1,1))/qq
         qq = qz(0,1)*qzz(0,1) - qz(1,1)*qzz(1,1)-
     &        qz(2,1)*qzz(2,1) - qz(3,1)*qzz(3,1)
         xaon = dotrc(qzz(0,1),atau(1,1))/qq

cc----------------------------------------------------------
         xze = (0.0d0,0.0d0)
         xae = (0.0d0,0.0d0)
         xzmu = (0.0d0,0.0d0)
         xamu = (0.0d0,0.0d0)
         xaon = (0.0d0,0.0d0)
c------------------------------------------------------------

         do mu = 0,3
            zepp(mu+1,1) = ze(mu+1,1) - xze*qze(mu,1)
            zmup(mu+1,1) = zmu(mu+1,1) - xzmu*qzmu(mu,1)
            aep(mu+1,1) = ae(mu+1,1) - xae*qze(mu,1)
            amup(mu+1,1) = amu(mu+1,1) - xamu*qzmu(mu,1)
            ataup(mu+1,1) = atau(mu+1,1) - xaon*qz(mu,1)
         enddo
         do mu = 5,6
            zepp(mu,1) = ze(mu,1)
            zmup(mu,1) = zmu(mu,1)
            aep(mu,1) = ae(mu,1)
            amup(mu,1) = amu(mu,1)
            ataup(mu,1) = atau(mu,1)
         enddo         

      endif
#endif
c*****************************************************************

         call zza_atozz_had(v(0,1,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+
         call zza_ztozz_had(v(0,1,1),1,zzztens(0,1),0)      !Z -> Z Z -> e- e+ mu- mu+
         call zza_atozae_had(v(0,1,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
         call zza_ztozae_had(v(0,1,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
         call zza_atozam_anomal(v(0,1,1),1,h,aazm(0,1),0)  !A -> Z A -> mu- mu+ a
         call zza_ztozam_anomal(v(0,1,1),1,h,zzam(0,1),0)  !Z -> Z A -> mu- mu+ a
         call zza_atozza_had(v(0,1,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
         call zza_ztozza_had(v(0,1,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a


C*******************  q1 qbar2 ---> Z Z A   **************************

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

      call qq_ZZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)

C*******************  q1bar q2 ---> Z Z A   **********************

      physToDiag(1)=2     !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZA(xi,p,v,physToDiag,fsign,bos,nlo,pdf,res,nmin,nmax)


c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select

c*****************  end of process evaluation  *********************

      if (init.eq.1) then
         init = init+1
      endif

      SELECT CASE(procid)
      CASE(ZZhadA)
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
         CASE(ZZhadA)
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
      
      if (bos .eq. 221) then
        m2s = m2s * 8d0           ! factor 2 for random helicity summation
      else if (bos .eq. 251) then
        m2s = m2s * 4d0           ! factor 2 for random helicity summation - except Z->nu nu~
      endif
      
      m2s_qqZZA = m2s
      
      RETURN
      END
      
c*****************************************************************************
c
c    end function m2s_qqZZA
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZZA
c
c*****************************************************************

      subroutine qq_ZZA(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax)
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
      double precision uuzza, ddzza, uuzza2, ddzza2
      double precision pbar(0:3,3+n_v), dotrr, polcol
      double precision NCmatrixelt(0:1,0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VVV
      external dotrr, FL_VVV

      logical ldebug
      parameter (ldebug=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu

      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)

c  helicity selection, for debugging
      INTEGER h
      COMMON /hcount / h

      uuzza = 0d0
      uuzza2 = 0d0
      ddzza = 0d0
      ddzza2 = 0d0

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

      do mu = 0,3             ! kinematics for ZZA
         pbar(mu,3) = v(mu,1) ! e-
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! mu+
         pbar(mu,7) = v(mu,5) ! A
         pbar(mu,8) = 0.0d0   !additional momentum not used for VVA processes
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
        if (bos .eq. 221) then
         call suu_zza(pbar,fsign,uuzza2)
         call sdd_zza(pbar,fsign,ddzza2)
        else if (bos .eq. 251) then
         call suu_zzna(pbar,fsign,uuzza2)
         call sdd_zzna(pbar,fsign,ddzza2)
        endif
      endif
#endif

      call QQZZA(pbar,fsign,nlo,1,uuzza, ddzza)

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         call zza_atozz_mg2_had(v(0,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+
         call zza_ztozz_mg2_had(v(0,1),1,zzztens(0,1),0)      !Z -> Z Z -> e- e+ mu- mu+
         call zza_atozae_mg2_had(v(0,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
         call zza_ztozae_mg2_had(v(0,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
         call zza_atozza_mg2_had(v(0,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
         call zza_ztozza_mg2_had(v(0,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a
         call QQZZA(pbar,fsign,nlo,1,uuzza2, ddzza2)
         call zza_atozz_had(v(0,1),1,azz(0,1),0)      !A -> Z Z -> e- e+ mu- mu+
         call zza_ztozz_had(v(0,1),1,zzztens(0,1),0)      !Z -> Z Z -> e- e+ mu- mu+
         call zza_atozae_had(v(0,1),1,h,aaze(0,1),0)  !A -> Z A -> e- e+ a
         call zza_ztozae_had(v(0,1),1,h,zzae(0,1),0)  !Z -> Z A -> e- e+ a
         call zza_atozza_had(v(0,1),1,h,azza(0,1),0)  !A -> ZZA -> e- e+ mu- mu+ a
         call zza_ztozza_had(v(0,1),1,h,zzza(0,1),0)  !Z -> ZZA -> e- e+ mu- mu+ a
         print*,"1 - uuzza/uuzza2 =   ",abs(1-uuzza/uuzza2)
         print*,"1 - ddzza/ddzza2 =   ",abs(1-ddzza/ddzza2)
         read(*,*)
      endif
      if (ldebug) then
         print*,"uuzza =   ",uuzza
         print*,"uuzza2 =  ",uuzza2
         print*,"ddzza =   ",ddzza
         print*,"ddzza2 =  ",ddzza2
         read(*,*)
      endif
#endif

      NCmatrixelt(0,0)=uuzza
      NCmatrixelt(1,1)=ddzza

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
c   end subroutine qq_ZZA
c
c*****************************************************************

      subroutine vtozza_reset
      implicit none
#include "tensor.inc"
      double complex zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            azz(mu,j)  = zero
            zzztens(mu,j)  = zero
            aaze(mu,j)  = zero
            zzae(mu,j)  = zero
            aazm(mu,j)  = zero
            zzam(mu,j)  = zero
            azza(mu,j) = zero
            zzza(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
