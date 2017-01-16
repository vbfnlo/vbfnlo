*****************************************************************************
c
c    begin function m2s_qqZZW
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008 January 12
c	Last modified:  2009 January 28
c
c       this function computes |ME|^2 for qq -> ZZW+-
c	at LO and virtual corrections via call of "qqZZW.f"
c
c*****************************************************************************

      DOUBLE PRECISION FUNCTION m2s_qqZZW(
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
      double complex dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

c helicity selection
      INTEGER h
      COMMON /hcount/ h

c declare local variables
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

         if(bos.eq.223) then
            write(6,*) " "
            write(6,*) "ZZW+ amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else if (bos.eq.224) then
            write(6,*) " "
            write(6,*) "ZZW- amplitude square information:"
            write(6,*) "-----------------------------------------------"
            write(6,*) " "
         else
            write(6,*) " wrong call of m2s_qqZZW, bos = ",bos
            write(6,*) "-----------------------------------------------"
         endif
         call vtozzw_reset

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

c lepton spinors and Z, Z, W polarisation vectors
c lepton helicities for Z -> l- l+ not fixed
c -> sum over all possible helicities in |M|**2

c select helicity: h ... random number for lepton helicity
c                        combination (h=1:4)

        h = mod(h,4) + 1
        ie = sign(1,2-h)
        iu = (-1)**(h+1)

c...Les Houches interface
        if ((lha.or.hepmc).and..not.doNLO) then
           helicity(1)=-ie
           helicity(2)= ie
           helicity(3)=-iu
           helicity(4)= iu
           helicity(5)=-1
           helicity(6)= 1
        endif

      do mu = 0,3
         qze(mu,1) = v(mu,1,1)+v(mu,2,1)   !Z
         qzmu(mu,1) = v(mu,3,1)+v(mu,4,1)  !Z
         qwbos(mu,1) = v(mu,5,1)+v(mu,6,1)   !W
         qzz(mu,1) = qze(mu,1)+qzmu(mu,1)   !ZZ
         qzzw(mu,1) = qze(mu,1) + qzmu(mu,1) + qwbos(mu,1)
      enddo
      qze(4,1) = dotrr(qze(0,1),qze(0,1))
      qzmu(4,1) = dotrr(qzmu(0,1),qzmu(0,1))
      qwbos(4,1) = dotrr(qwbos(0,1),qwbos(0,1))
      qzz(4,1) = dotrr(qzz(0,1),qzz(0,1))
      qzzw(4,1) = dotrr(qzzw(0,1),qzzw(0,1))

      CALL OXXXXX(v(0,1,1),ZERO ,-ie,+1,zem(1,1)) !e-
      CALL IXXXXX(v(0,2,1),ZERO ,+ie,-1,zep(1,1)) !e+
      CALL OXXXXX(v(0,3,1),ZERO ,-iu,+1,zmm(1,1)) !mu-
      CALL IXXXXX(v(0,4,1),ZERO ,+iu,-1,zmp(1,1)) !mu+
      CALL OXXXXX(v(0,5,1),ZERO , -1,+1,ztm(1,1)) !vt !tau-
      CALL IXXXXX(v(0,6,1),ZERO ,  1,-1,ztp(1,1)) !ta+ !!vta~

      SELECT CASE(procid)
      CASE(ZZhadWP,ZZhadWM)
         ! for general output to all combinations: up-type first
         if ((finalquarks(1).eq.93 .or. finalquarks(1).eq.94)) then
           ! set couplings for Z hadronic decay into uubar
           call setZtouu
         endif
         CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO , ZERO ,ae(1,1))   !A
         CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH,ze(1,1))   !Z

      CASE DEFAULT
         CALL JIOXXX(zep(1,1),zem(1,1),GAL,ZERO , ZERO ,ae(1,1))   !A
         CALL JIOXXX(zep(1,1),zem(1,1),GZL,ZMASS,ZWIDTH,ze(1,1))   !Z

      END SELECT
      CALL JIOXXX(zmp(1,1),zmm(1,1),GAL,ZERO ,ZERO  ,amu(1,1))  !A
      CALL JIOXXX(zmp(1,1),zmm(1,1),GZL,ZMASS,ZWIDTH,zmu(1,1))  !Z

      CALL JIOXXX(ztp(1,1),ztm(1,1),GWF,WMASS,WWIDTH,wbos(1,1))!W+ !W-


c************reduction of pentagon contribution *************
#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c ze=zepp-xze*qze, zmu=zmup-xzmu*qmu, wbos=wbosp-xw*qwbos
c ae=aep-xae*qze, amu=amup-xamu*qmu
c first choice
            qq = qze(0,1)*qzz(0,1)-qze(1,1)*qzz(1,1)- qze(2,1)*qzz(2,1)
     &          -qze(3,1)*qzz(3,1)
            xze = dotrc(qzz(0,1),ze(1,1))/qq
            xae = dotrc(qzz(0,1),ae(1,1))/qq
            qq = qzmu(0,1)*qzz(0,1)-qzmu(1,1)*qzz(1,1)-qzmu(2,1)*qzz(2,1)
     &          -qzmu(3,1)*qzz(3,1)
            xzmu = dotrc(qzz(0,1),zmu(1,1))/qq
            xamu = dotrc(qzz(0,1),amu(1,1))/qq
            qq = qwbos(0,1)*qzz(0,1)-qwbos(1,1)*qzz(1,1)-qwbos(2,1)*qzz(2,1)
     &          -qwbos(3,1)*qzz(3,1)
            xw = dotrc(qzz(0,1),wbos(1,1))/qq
c--------------------------------------------------------
            do mu = 0,3
               zepp(mu+1,1) = ze(mu+1,1) - xze*qze(mu,1)
               zmup(mu+1,1) = zmu(mu+1,1) - xzmu*qzmu(mu,1)
               aep(mu+1,1) = ae(mu+1,1) - xae*qze(mu,1)
               amup(mu+1,1) = amu(mu+1,1) - xamu*qzmu(mu,1)
               wbosp(mu+1,1) = wbos(mu+1,1) - xw*qwbos(mu,1)
            enddo
            do mu = 5,6
               zepp(mu,1) = ze(mu,1)
               aep(mu,1) = ae(mu,1)
               zmup(mu,1) = zmu(mu,1)
               amup(mu,1) = amu(mu,1)
               wbosp(mu,1) = wbos(mu,1)
            enddo

         endif
#endif
c*****************************************************************

c leptonic tensors, last entry defines whether the non resonant diagrams
c contribute ( = 0: all diagrams, = 1: no non resonant diagrams)

      SELECT CASE(procid)
      CASE(ZZhadWP, ZZhadWM)
              call anomal_formfactor(qzzw(0,1),qze(0,1),qzmu(0,1),qwbos(0,1))
              if (bos.eq.223) then
                 call wptozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              endif
              call zzw_atozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
              call zzw_ztozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
      CASE(ZZWPhad, ZZWMhad)
              call anomal_formfactor(qzzw(0,1),qze(0,1),qzmu(0,1),qwbos(0,1))
              if (bos.eq.223) then
                 call wptozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_had(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_had(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              endif
              call zzw_atozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
              call zzw_ztozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
      CASE DEFAULT
        if (with_anom) then ! anomalous gauge boson couplings
c             using global form factor for all tensors of one phase space point
c             this ensures proper cancellations for anomalous contributions
c             energy scale is invariant ZZW mass
              call anomal_formfactor(qzzw(0,1),qze(0,1),qzmu(0,1),qwbos(0,1))
              if (bos.eq.223) then
                 call wptozzw_anomal(v(0,1,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_anomal(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw_anomal(v(0,1,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_anomal(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              endif
              if (NW.eq.0) then
                 call atozz_anomal(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
                 call ztozz_anomal(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
              endif
        else if (with_kk) then ! kk scenario
#ifdef WITH_KK
              if (bos.eq.223) then
                 call wptozzw_kk(v(0,1,1),1,h,wzzw(0,1),NW)
                 call zzw_wptowz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw_kk(v(0,1,1),1,h,wzzw(0,1),NW)
                 call zzw_wmtowz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              endif
              if (NW.eq.0) then
                 call atozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
                 call ztozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
              endif
#endif
        else  ! SM
              if (bos.eq.223) then
                 call wptozzw(v(0,1,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw(v(0,1,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)
              endif
              if (NW.eq.0) then
                 call atozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
                 call ztozz(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))
              endif
        endif
      END SELECT

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         if (bos.eq.223) then
            call wptozzw_mg_had(v(0,1,1),1,h,mg2_tens,NW) ! W+ -> ZZW
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wzzw(i,1) / mg2_tens(i) - 1)
              print*, wzzw(i,1) , mg2_tens(i), abs(wzzw(i,1) / mg2_tens(i))
            enddo
            print*, ""
            print*, "tensor w+ to zzw:", sumtest
            call zzw_wptowz_mg_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,mg2_tens,NW)
            sumtest = 0
            do i = 0,5
              sumtest = sumtest + abs( wzw1(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w+ to zw:", sumtest
         else if (bos.eq.224) then
            call wmtozzw_mg_had(v(0,1,1),1,h,mg2_tens,NW) ! W- -> ZZW
            sumtest = 0
            do i = 0,3
              sumtest = sumtest + abs( wzzw(i,1) / mg2_tens(i) - 1)
              print*, wzzw(i,1) , mg2_tens(i), abs(wzzw(i,1) / mg2_tens(i))
            enddo
            print*, ""
            print*, "tensor w- to zzw:", sumtest
            call zzw_wmtowz_mg_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,mg2_tens,NW)
            sumtest = 0
            do i = 0,3
              sumtest = sumtest + abs( wzw1(i,1) / mg2_tens(i) - 1)
            enddo
            print*, "tensor w- to zw:", sumtest
         endif
         call zzw_atozz_mg_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,mg2_tens)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( azezm(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor a to zz:", sumtest
         call zzw_ztozz_mg_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,mg2_tens)
         sumtest = 0
         do i = 0,5
           sumtest = sumtest + abs( zzezm(i,1) / mg2_tens(i) - 1)
         enddo
         print*, "tensor z to zz:", sumtest
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

C*******************  q1 qbar2 ---> Z Z W   **************************

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

      call qq_ZZW(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW)

C*******************  q1bar q2 ---> Z Z W  **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZW(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW)


c**************  end of process evaluation part one ******************


! now run again for down-type decay products if 93 93 / 94 94 is 
! requested and the z decays hadronically.

      SELECT CASE(procid)
      CASE(ZZhadWP, ZZhadWM)
       if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
        ! set couplings for Z hadronic decay into ddbar
        call setZtodd
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_AF,ZERO , ZERO ,ae(1,1))   !A
        CALL JIOXXX(zep(1,1),zem(1,1),GZ_ZF,ZMASS,ZWIDTH,ze(1,1))   !Z
        if (bos.eq.223) then
           call wptozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
           call zzw_wptowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
           call zzw_wptowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
        else if (bos.eq.224) then
           call wmtozzw_had(v(0,1,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
           call zzw_wmtowz_had(v(0,1,1),v(0,2,1),v(0,5,1),v(0,6,1),1,-ie,wzw1(0,1),NW)
           call zzw_wmtowz_anomal(v(0,3,1),v(0,4,1),v(0,5,1),v(0,6,1),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
        endif
        call zzw_atozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,azezm(0,1))
        call zzw_ztozz_had(v(0,1,1),v(0,2,1),v(0,3,1),v(0,4,1),-ie,-iu,zzezm(0,1))

c************reduction of pentagon contribution *************
#ifdef WITH_NLO
         if (abs(nlo).eq.5 .or. nlo.eq.-4) then
c ze=zepp-xze*qze, zmu=zmup-xzmu*qmu, wbos=wbosp-xw*qwbos
c ae=aep-xae*qze, amu=amup-xamu*qmu
c first choice
            qq = qze(0,1)*qzz(0,1)-qze(1,1)*qzz(1,1)- qze(2,1)*qzz(2,1)
     &          -qze(3,1)*qzz(3,1)
            xze = dotrc(qzz(0,1),ze(1,1))/qq
            xae = dotrc(qzz(0,1),ae(1,1))/qq
            qq = qzmu(0,1)*qzz(0,1)-qzmu(1,1)*qzz(1,1)-qzmu(2,1)*qzz(2,1)
     &          -qzmu(3,1)*qzz(3,1)
            xzmu = dotrc(qzz(0,1),zmu(1,1))/qq
            xamu = dotrc(qzz(0,1),amu(1,1))/qq
            qq = qwbos(0,1)*qzz(0,1)-qwbos(1,1)*qzz(1,1)-qwbos(2,1)*qzz(2,1)
     &          -qwbos(3,1)*qzz(3,1)
            xw = dotrc(qzz(0,1),wbos(1,1))/qq
c--------------------------------------------------------
            do mu = 0,3
               zepp(mu+1,1) = ze(mu+1,1) - xze*qze(mu,1)
               zmup(mu+1,1) = zmu(mu+1,1) - xzmu*qzmu(mu,1)
               aep(mu+1,1) = ae(mu+1,1) - xae*qze(mu,1)
               amup(mu+1,1) = amu(mu+1,1) - xamu*qzmu(mu,1)
               wbosp(mu+1,1) = wbos(mu+1,1) - xw*qwbos(mu,1)
            enddo
            do mu = 5,6
               zepp(mu,1) = ze(mu,1)
               aep(mu,1) = ae(mu,1)
               zmup(mu,1) = zmu(mu,1)
               amup(mu,1) = amu(mu,1)
               wbosp(mu,1) = wbos(mu,1)
            enddo

         endif
#endif
c*****************************************************************

C*******************  q1 qbar2 ---> Z Z W   **************************

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

      call qq_ZZW(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax, NW)

C*******************  q1bar q2 ---> Z Z W  **********************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) = -1

      call qq_ZZW(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax, NW)

c***  end of extra process evaluation for hadronic decay  ***

       endif
      end select


c*****************  end of process evaluation  *********************


      SELECT CASE(procid)
      CASE(ZZWPhad, ZZWMhad)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,5,1),v(0,6,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(ZZhadWM, ZZhadWP)
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
         CASE(ZZWPhad)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            endif
         CASE(ZZWMhad)
            if (finalquarks(1).eq.93 .or. finalquarks(1).eq.94) then
               rnumb=RandomNumber()
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            endif
         CASE(ZZhadWP, ZZhadWM)
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

      m2s = m2s * 4.0d0 !factor for random helicity summation

      m2s_qqZZW = m2s
      
      RETURN
      END

c*****************************************************************************
c
c    end function m2s_qqZZW
c
c*****************************************************************************

c******************************************************************
c
c   begin subroutine qq_ZZW
c
c*****************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008 January 15
c	Last modified:  2009 January 28
c       
c
c       qq_ZZW calls the amplitude square routines 
c             qqZZW    for qq --> ZZW --> 6 leptons
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
c*****************************************************************************

      SUBROUTINE qq_ZZW(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax, NW)
      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      DOUBLE PRECISION p(0:3,n_p), v(0:3,n_v), xi(nx)
      DOUBLE PRECISION pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax
      integer NW

      DOUBLE PRECISION N                              ! color factors
      parameter(N=3d0)
c
c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/process.inc"
#include "tensor.inc"

c local variables
      double precision udzzw, duzzw, udzzw2, duzzw2
      double precision pbar(0:3,2+n_v), dotrr, polcol
      double precision NCmatrixelt(0:1)
      double complex zero
      integer iflav(4), diagToPhys(2), FL_VVV
      external dotrr, FL_VVV

      logical ldebug
      parameter (ldebug=.false.)
      LOGICAL ldebug_ltens
      PARAMETER (ldebug_ltens=.false.)
c  helicity selection, for debugging
      INTEGER h
      COMMON /hcount / h
      parameter (zero = (0d0,0d0) )

      integer i,if1,k,mu

      udzzw = 0.0d0
      udzzw2 = 0.0d0
      duzzw = 0.0d0
      duzzw2 = 0.0d0

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

      do mu = 0,3             ! kinematics for ZZW --> 6 lepton decay
         pbar(mu,3) = v(mu,1) ! e-
         pbar(mu,4) = v(mu,2) ! e+
         pbar(mu,5) = v(mu,3) ! mu-
         pbar(mu,6) = v(mu,4) ! mu+
         pbar(mu,7) = v(mu,5) ! vt - !ta-
         pbar(mu,8) = v(mu,6) ! ta+  !vt~
      enddo

#ifdef WITH_MADGRAPH
      if (ldebug) then
         if (bos.eq.223) then
            call sud_zzw(pbar,fsign,udzzw2)
         else if (bos.eq.224) then
            call sdu_zzw(pbar,fsign,duzzw2)
         endif
      endif
#endif

c call of the subroutine to calculate matrix element squared
      if (bos.eq.223) then
#ifdef WITH_KK
         if (.not.mgtest) then
            call qqZZW(pbar,fsign,nlo,1,udzzw,bos, NW)
         else
            call sud_zzw_mhlm(pbar,fsign,udzzw)
         endif
#else
         call qqZZW(pbar,fsign,nlo,1,udzzw,bos, NW)
#endif
      else if (bos.eq.224) then
#ifdef WITH_KK
         if (.not. mgtest) then
            call qqZZW(pbar,fsign,nlo,1,duzzw,bos, NW)
         else
            call sdu_zzw_mhlm(pbar,fsign,duzzw)
         endif
#else
         call qqZZW(pbar,fsign,nlo,1,duzzw,bos, NW)
#endif
      endif

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
         ie = sign(1,2-h)
         iu = (-1)**(h+1)
         SELECT CASE(procid)
         CASE(ZZhadWP, ZZhadWM)
              if (bos.eq.223) then
                 call wptozzw_mg_had(v(0,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_mg_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_anomal(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              else if (bos.eq.224) then
                 call wmtozzw_mg_had(v(0,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_mg_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_anomal(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              endif
              call zzw_atozz_mg_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
              call zzw_ztozz_mg_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
         CASE(ZZWPhad, ZZWMhad)
              if (bos.eq.223) then
                 call wptozzw_mg_had(v(0,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_mg_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_mg_had(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw_mg_had(v(0,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_mg_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_mg_had(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)
              endif
              call zzw_atozz_mg_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
              call zzw_ztozz_mg_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
         end select
         if (bos.eq.223) then
            call qqZZW(pbar,fsign,nlo,1,udzzw2,bos, NW)
            print*,"1 - udzzw2/udzzw =   ",abs(1-udzzw2/udzzw)
         else if (bos.eq.224) then
            call qqZZW(pbar,fsign,nlo,1,duzzw2,bos, NW)
            print*,"1 - duzzw2/duzzw =   ",abs(1-duzzw2/duzzw)
         endif
         SELECT CASE(procid)
         CASE(ZZhadWP, ZZhadWM)
              if (bos.eq.223) then
                 call wptozzw_had(v(0,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_anomal(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_anomal(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)       ! no hadronic decay of the second Z
              endif
              call zzw_atozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
              call zzw_ztozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
         CASE(ZZWPhad, ZZWMhad)
              if (bos.eq.223) then
                 call wptozzw_had(v(0,1),1,h,wzzw(0,1),NW) ! W+ -> ZZW
                 call zzw_wptowz_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wptowz_had(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)
              else if (bos.eq.224) then
                 call wmtozzw_had(v(0,1),1,h,wzzw(0,1),NW) ! W- -> ZZW
                 call zzw_wmtowz_had(v(0,1),v(0,2),v(0,5),v(0,6),1,-ie,wzw1(0,1),NW)
                 call zzw_wmtowz_had(v(0,3),v(0,4),v(0,5),v(0,6),1,-iu,wzw2(0,1),NW)
              endif
              call zzw_atozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,azezm(0,1))
              call zzw_ztozz_had(v(0,1),v(0,2),v(0,3),v(0,4),-ie,-iu,zzezm(0,1))
         end select
         read(*,*)
      endif
      if (ldebug) then
         print*,"udzzw=   ",udzzw
         print*,"udzzw2 =  ",udzzw2
         print*,"udzzw/udzzw2 = ",udzzw/udzzw2
         print*,"udzzw2/udzzw - 1 = ",udzzw2/udzzw - 1
      endif
#endif

      NCmatrixelt(0)=udzzw
      NCmatrixelt(1)=duzzw

      do if1 = 1,(nfl/2)*2,2
         if (bos.eq.223) then
            iflav(1) = (if1-physToDiag(1)+2)*fsign(physToDiag(1))
            iflav(2) = (if1+physToDiag(1)-1)*fsign(physToDiag(2))
         else if (bos.eq.224) then
            iflav(1) = (if1+physToDiag(1)-1)*fsign(physToDiag(1))
            iflav(2) = (if1-physToDiag(1)+2)*fsign(physToDiag(2))
         endif
 
         k=FL_VVV(iflav,1)

         res(k)=pdf(sign1*iflav(1),1)
     &         *pdf(sign2*iflav(2),2)
     &         *NCmatrixelt(bos-223)*polcol
      enddo

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_ZZW
c
c*****************************************************************

c*****************************************************************************
c
c    begin subroutine vtozzw_reset
c
c*****************************************************************************
c
c       Vera Hankele, <vera@particle.uni-karlsruhe.de>
c	Initial version:  2008 January 20
c	Last modified:  2009 January 28
c       
c  this subroutines sets all leptonic tensors needed in the calculation 
c  of pp -> ZZW to zero.
c  
c*****************************************************************************

      SUBROUTINE vtozzw_reset

      IMPLICIT NONE

#include "tensor.inc"

      complex*16 zero
      parameter (zero=(0d0,0d0))
      integer j,mu

      do mu = 0,5
         do j = 1,3
            wzzw(mu,j) = zero
            wzw1(mu,j) = zero
            wzw2(mu,j) = zero
            azezm(mu,j) = zero
            zzezm(mu,j) = zero
         enddo
      enddo
      print*," precalculated tensors initialized to ",zero 
      return
      end
c*****************************************************************************
c
c    end subroutine vtozzw_reset
c
c*****************************************************************************
