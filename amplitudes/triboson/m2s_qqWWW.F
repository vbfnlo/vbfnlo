*****************************************************************************
c
c    begin function m2s_qqWWW
c
c*****************************************************************************
c
c       Stefan Prestel, <sprestel@particle.uni-karlsruhe.de>
c	Initial version:  
c	Last modified by Vera Hankele:  2009 January 28
c
c this function computes |ME|^2 for qq-> ve e+ mu-vm~ vt/ta- ta+/vt~ ("WWW")
c at LO and virtual corrections via call of 'qqWWWm' or 'qqWWWp' depending on
c overall charge
c
c*****************************************************************************

      DOUBLE PRECISION FUNCTION m2s_qqWWW(
     &      bos,  !in:  XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
     &      nlo,  !in:  Born only for nlo=0, +interference for 1
     &      xi,   !in:  Feynman x parameters of incoming fermions
     &      p,    !in:  momenta associated with external fermions
     &      v,    !in:  momenta associated with gauge boson dec. fermions
     &      rn,   !in:  random number array
     &      NW    !in:  NW=1: Narrow width approximation,
                  !     NW=0: full calculation
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

c declare input/output variables
      INTEGER bos, nlo
      double precision xi(nx), p(0:3,max_p), v(0:3,n_v,3), rn(1), m2s
      INTEGER NW

c declare external functions
      double precision dotrr
      double complex dotrc
      INTEGER FL_VVV
      EXTERNAL dotrr, dotrc
      EXTERNAL FL_VVV

c declare local variables
      double complex zero
      PARAMETER (zero=(0d0,0d0))

      double precision q_sf, qq
      INTEGER init/0/, I, mu
      SAVE init

c declare variables for summation process
      INTEGER ires

      INTEGER FSIGN(2+n_v), physToDiag(2), nmin, nmax
      SAVE nmax

c variables for leshouches interface
      double precision weight,rnumb,RandomNumber

c store contributions from subprocess ip in res(ip) where
      double precision res(maxnumsubproc)


c declare parton distribution variables
      double precision pdf(-6:6,2)

c define program switches
      LOGICAL ldebug, ldebug_ltens
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebug_ltens=.false.)
      complex*16 mg2_tens(0:5)
      double precision sumtest
      double precision swmomenta(0:3,6)

c --------------------------- BEGIN CODE --------------------------------------

      if ( init .eq. 0 ) then

         write(6,*) " "
         if( bos.eq.343 ) then
            write(6,*) "W+W-W+ -Amplitude square information:"
         else if(bos.eq.434) then
            write(6,*) "W-W+W- -Amplitude square information:"
         else 
            write(6,*) "Process not implemented, check boson identifier"
         endif
         write(6,*) "-----------------------------------------------"
         write(6,*) " "
         call vtowww_reset

         call printnfl(.true.)
         if (NW.eq.0) then
            print*,"Calculation of amplitude with all off-shell contributions"
         else
            print*,"Calculation of amplitude in narrow width approximation"
         endif
         init = 1
         do i = 1,maxnumsubproc
            res(i) = 0
         enddo
         print*
      else
         do i = 1,nmax
            res(i) = 0
         enddo
      endif

c Les Houches interface
      if((lha.or.hepmc).and.(.not.doNLO)) then
         helicity(1) = -1
         helicity(2) =  1
         helicity(3) = -1
         helicity(4) =  1
         helicity(5) = -1
         helicity(6) =  1
      endif

c switch momenta for W-W+W- so that electron spinor is always created out of v(0,1,j)
      if(bos.eq.434) then
         do mu=0,3
         swmomenta(mu,1) = v(mu,3,1) ! ve
         swmomenta(mu,2) = v(mu,4,1) ! e+
         swmomenta(mu,3) = v(mu,1,1) ! mu-
         swmomenta(mu,4) = v(mu,2,1) ! vm~
         enddo
         do i=1,4
            do mu=0,3
                v(mu,i,1) = swmomenta(mu,i)
            enddo
         enddo
      endif

c lepton spinors and W+, W-, W+/- polarisation vectors

      if(bos.eq.343) then  ! lepton-wavefunctions for u d~ --> w+ w- w+
         CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve1(1,1)) !W(1,3)         !ve
         CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep1(1,1)) !W(1,4)         !e+
         CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu1(1,1)) !W(1,5)         !mu-
         CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm1(1,1)) !W(1,6)         !vm~
         CALL OXXXXX(v(0,5,1),ZERO ,-1,1,wvt(1,1) ) !W(1,7)         !vtau
         CALL IXXXXX(v(0,6,1),ZERO ,1,-1,wtp(1,1)) !W(1,8)          !ta+
         CALL JIOXXX(wep1(1,1),wve1(1,1),GWF,WMASS,WWIDTH, wp1(1,1)) !W+ #1
         CALL JIOXXX(wvm1(1,1),wmu1(1,1),GWF,WMASS,WWIDTH, w1m(1,1)) !W- #1
         CALL JIOXXX(wtp(1,1) ,wvt(1,1) ,GWF,WMASS,WWIDTH, wp2(1,1)) !W+ #2
         do mu = 0,3
            qp1(mu,1)     = v(mu,1,1)+v(mu,2,1)       !W+ #1
            q1m(mu,1)     = v(mu,3,1)+v(mu,4,1)       !W-
            qp2(mu,1)     = v(mu,5,1)+v(mu,6,1)       !W+ #2
            qwp1wm(mu,1)  = qp1(mu,1)+ q1m(mu,1)      !momentum of W+#1 and W-
            qwp2wm(mu,1)  = qp2(mu,1)+ q1m(mu,1)      !momentum of W+#2 and W-
            qwpwmwp(mu,1) = qwp1wm(mu,1) + qp2(mu,1)  !momentum of W+W-W+
         enddo
         qp1(4,1)     = dotrr(qp1(0,1),qp1(0,1))
         q1m(4,1)     = dotrr(q1m(0,1),q1m(0,1))
         qp2(4,1)     = dotrr(qp2(0,1),qp2(0,1))
         qwp1wm(4,1)  = dotrr(qwp1wm(0,1),qwp1wm(0,1))
         qwp2wm(4,1)  = dotrr(qwp2wm(0,1),qwp2wm(0,1))
         qwpwmwp(4,1) = dotrr(qwpwmwp(0,1),qwpwmwp(0,1))

      else if (bos.eq.434) then ! lepton-wavefunctions for u~d --> w+ w- w-
         CALL OXXXXX(v(0,1,1),ZERO ,-1,1,wve2(1,1)) !ve
         CALL IXXXXX(v(0,2,1),ZERO ,1,-1,wep2(1,1)) !e+
         CALL OXXXXX(v(0,3,1),ZERO ,-1,1,wmu2(1,1)) !mu-
         CALL IXXXXX(v(0,4,1),ZERO ,1,-1,wvm2(1,1)) !vm~
         CALL OXXXXX(v(0,5,1),ZERO,-1,1,wtm(1,1) )  ! tau-
         CALL IXXXXX(v(0,6,1),ZERO,1,-1,wvtbar(1,1))! vtau~
         CALL JIOXXX(wep2(1,1),wve2(1,1),GWF,WMASS,WWIDTH, w1p(1,1)) !W+ #1
         CALL JIOXXX(wvm2(1,1),wmu2(1,1),GWF,WMASS,WWIDTH, wm1(1,1)) !W- #1
         CALL JIOXXX(wvtbar(1,1),wtm(1,1),GWF,WMASS,WWIDTH,wm2(1,1)) !W- #2
         do mu = 0,3
            q1p(mu,1)     = v(mu,1,1)+v(mu,2,1)       !W+ #2
            qm1(mu,1)     = v(mu,3,1)+v(mu,4,1)       !W- #1
            qm2(mu,1)     = v(mu,5,1)+v(mu,6,1)       !W- #2
            qwpwm1(mu,1)  = q1p(mu,1)+ qm1(mu,1)      !momentum of W+ and W-#1
            qwpwm2(mu,1)  = q1p(mu,1)+ qm2(mu,1)      !momentum of W+ and W-#2
            qwpwmwm(mu,1) = qwpwm1(mu,1) + qm2(mu,1)  !momentum of W+W-W-
         enddo
         q1p(4,1)     = dotrr(q1p(0,1),q1p(0,1))
         qm1(4,1)     = dotrr(qm1(0,1),qm1(0,1))
         qm2(4,1)     = dotrr(qm2(0,1),qm2(0,1))
         qwpwm1(4,1)  = dotrr(qwpwm1(0,1),qwpwm1(0,1))
         qwpwm2(4,1)  = dotrr(qwpwm2(0,1),qwpwm2(0,1))
         qwpwmwm(4,1) = dotrr(qwpwmwm(0,1),qwpwmwm(0,1))

      endif

#ifdef WITH_NLO
      if (abs(nlo).eq.5 .or. nlo.eq.-4) then
         if(bos.eq.343) then
            qq  = dotrr(qp1(0,1),qwpwmwp(0,1))
            xp1 = dotrc(qwpwmwp(0,1),wp1(1,1))/qq
            qq  = dotrr(qp2(0,1),qwpwmwp(0,1))
            xp2 = dotrc(qwpwmwp(0,1),wp2(1,1))/qq
            qq  = dotrr(q1m(0,1),qwpwmwp(0,1))
            x1m = dotrc(qwpwmwp(0,1),w1m(1,1))/qq
            do mu = 0,3
               wpp1(mu+1,1) = wp1(mu+1,1) - xp1*qp1(mu,1)
               wpp2(mu+1,1) = wp2(mu+1,1) - xp2*qp2(mu,1)
               w1mp(mu+1,1) = w1m(mu+1,1) - x1m*q1m(mu,1)
            enddo
            do mu = 5,6
               wpp1(mu,1) = wp1(mu,1)
               wpp2(mu,1) = wp2(mu,1)
               w1mp(mu,1) = w1m(mu,1)
            enddo
         else
            qq  = dotrr(qm1(0,1),qwpwmwm(0,1))
            xm1 = dotrc(qwpwmwm(0,1),wm1(1,1))/qq
            qq  = dotrr(qm2(0,1),qwpwmwm(0,1))
            xmm2= dotrc(qwpwmwm(0,1),wm2(1,1))/qq 
            qq  = dotrr(q1p(0,1),qwpwmwm(0,1))
            x1p = dotrc(qwpwmwm(0,1),w1p(1,1))/qq
            do mu = 0,3
               wmm1(mu+1,1) = wm1(mu+1,1) - xm1*qm1(mu,1)
               wmm2(mu+1,1) = wm2(mu+1,1) - xmm2*qm2(mu,1)
               w1pm(mu+1,1) = w1p(mu+1,1) - x1p*q1p(mu,1)
            enddo
            do mu = 5,6
               wmm1(mu,1) = wm1(mu,1)
               wmm2(mu,1) = wm2(mu,1)
               w1pm(mu,1) = w1p(mu,1)
            enddo
         endif
      endif   ! end nlo = +5,-5,-4 block
#endif

c calculate leptonic tensors

      SELECT CASE(procid)
      CASE(WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM)
        if(bos.eq.343) then
           call anomal_formfactor(qwpwmwp(0,1),qp1(0,1),q1m(0,1),qp2(0,1))
           call WWW_atoww1_had(v(0,1,1),1,343,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww2_had(v(0,1,1),1,343,WWW_aww(0,1,2),NW)           !A -> mu- vm~ vt   tau+
           call WWW_ztoww1_had(v(0,1,1),1,343,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww2_had(v(0,1,1),1,343,WWW_zww(0,1,2),NW)           !Z -> mu- vm~ vt   tau+
           call WWW_wptowww_had(v(0,1,1),1,1,wptowww(0,1),NW)              !W+ -> ve e+ mu- vm~ vt ta+
        else
           call anomal_formfactor(qwpwmwm(0,1),q1p(0,1),qm1(0,1),qm2(0,1))
           call WWW_atoww1_had(v(0,1,1),1,434,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww3_had(v(0,1,1),1,434,WWW_aww(0,1,3),NW)           !A -> ve  e+  tau- vt~
           call WWW_ztoww1_had(v(0,1,1),1,434,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww3_had(v(0,1,1),1,434,WWW_zww(0,1,3),NW)           !Z -> ve  e+  tau- vt~
           call WWW_wmtowww_had(v(0,1,1),1,1,wmtowww(0,1),NW)              !W- -> ve e+ mu- vm~ ta- vt~
        endif

      CASE DEFAULT
      if (with_anom) then ! anomalous gauge boson couplings
c       using global form factor for all tensors of one phase space point
c       this ensures proper cancellations for anomalous contributions
c       energy scale is invariant WWW mass
        if(bos.eq.343) then
           call anomal_formfactor(qwpwmwp(0,1),qp1(0,1),q1m(0,1),qp2(0,1))
           call WWW_atoww1_anomal(v(0,1,1),1,343,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww2_anomal(v(0,1,1),1,343,WWW_aww(0,1,2),NW)           !A -> mu- vm~ vt   tau+
           call WWW_ztoww1_anomal(v(0,1,1),1,343,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww2_anomal(v(0,1,1),1,343,WWW_zww(0,1,2),NW)           !Z -> mu- vm~ vt   tau+
           call WWW_wptowww_anomal(v(0,1,1),1,1,wptowww(0,1),NW)              !W+ -> ve e+ mu- vm~ vt ta+
        else
           call anomal_formfactor(qwpwmwm(0,1),q1p(0,1),qm1(0,1),qm2(0,1))
           call WWW_atoww1_anomal(v(0,1,1),1,434,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww3_anomal(v(0,1,1),1,434,WWW_aww(0,1,3),NW)           !A -> ve  e+  tau- vt~
           call WWW_ztoww1_anomal(v(0,1,1),1,434,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww3_anomal(v(0,1,1),1,434,WWW_zww(0,1,3),NW)           !Z -> ve  e+  tau- vt~
           call WWW_wmtowww_anomal(v(0,1,1),1,1,wmtowww(0,1),NW)              !W- -> ve e+ mu- vm~ ta- vt~
        endif
      else if (with_kk) then ! kk scenario
#ifdef WITH_KK
        if(bos.eq.343) then
           call WWW_atoww1(v(0,1,1),1,343,WWW_aww(0,1,1),NW)           !A -> W+ W- -> ve  e+  mu-  vm~
           call WWW_atoww2(v(0,1,1),1,343,WWW_aww(0,1,2),NW)           !A -> W+ W- -> mu- vm~ vt   tau+
           call WWW_ztoww1(v(0,1,1),1,343,WWW_zww(0,1,1),NW)           !Z -> W+ W- -> ve  e+  mu-  vm~
           call WWW_ztoww2(v(0,1,1),1,343,WWW_zww(0,1,2),NW)           !Z -> W+ W- -> mu- vm~ vt   tau+
           call WWW_wptowww_kk(v(0,1,1),1,1,wptowww(0,1),NW)           !W+ -> WWW -> ve e+ mu- vm~ vt ta+
        else ! bos = 434
           call WWW_atoww1(v(0,1,1),1,434,WWW_aww(0,1,1),NW)           !A -> W+ W- -> ve  e+  mu-  vm~
           call WWW_atoww3(v(0,1,1),1,434,WWW_aww(0,1,3),NW)           !A -> W+ W- -> ve  e+  tau- vt~
           call WWW_ztoww1(v(0,1,1),1,434,WWW_zww(0,1,1),NW)           !Z -> W+ W- -> ve  e+  mu-  vm~
           call WWW_ztoww3(v(0,1,1),1,434,WWW_zww(0,1,3),NW)           !Z -> W+ W- -> ve  e+  tau- vt~
           call WWW_wmtowww_kk(v(0,1,1),1,1,wmtowww(0,1),NW)           !W- -> WWW -> ve e+ mu- vm~ ta- vt~
        endif
#endif
      else ! SM
        if(bos.eq.343) then
           call WWW_atoww1(v(0,1,1),1,343,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww2(v(0,1,1),1,343,WWW_aww(0,1,2),NW)           !A -> mu- vm~ vt   tau+
           call WWW_ztoww1(v(0,1,1),1,343,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww2(v(0,1,1),1,343,WWW_zww(0,1,2),NW)           !Z -> mu- vm~ vt   tau+
           call WWW_wptowww(v(0,1,1),1,1,wptowww(0,1),NW)              !W+ -> ve e+ mu- vm~ vt ta+
        else
           call WWW_atoww1(v(0,1,1),1,434,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
           call WWW_atoww3(v(0,1,1),1,434,WWW_aww(0,1,3),NW)           !A -> ve  e+  tau- vt~
           call WWW_ztoww1(v(0,1,1),1,434,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
           call WWW_ztoww3(v(0,1,1),1,434,WWW_zww(0,1,3),NW)           !Z -> ve  e+  tau- vt~
           call WWW_wmtowww(v(0,1,1),1,1,wmtowww(0,1),NW)              !W- -> ve e+ mu- vm~ ta- vt~
        endif

      endif
      END SELECT


#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
        if(bos.eq.343) then
           print*, ""
           call WWW_atoww1_mg2_had(v(0,1,1),1,343,mg2_tens,NW)           !A -> ve  e+  mu-  vm~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_aww(i,1,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
!             print*,  WWW_aww(i,1,1) , mg2_tens(i), abs( (WWW_aww(i,1,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor a to ww 1:", sumtest
           call WWW_atoww2_mg2_had(v(0,1,1),1,343,mg2_tens,NW)           !A -> mu- vm~ vt   tau+
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_aww(i,1,2)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor a to ww 2:", sumtest
           call WWW_ztoww1_mg2_had(v(0,1,1),1,343,mg2_tens,NW)           !Z -> ve  e+  mu-  vm~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_zww(i,1,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor z to ww 1:", sumtest
           call WWW_ztoww2_mg2_had(v(0,1,1),1,343,mg2_tens,NW)           !Z -> mu- vm~ vt   tau+
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_zww(i,1,2)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor z to ww 2:", sumtest
           call WWW_wptowww_mg2_had(v(0,1,1),1,1,mg2_tens,NW)              !W+ -> ve e+ mu- vm~ vt ta+
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (wptowww(i,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
!             print*, wptowww(i,1),mg2_tens(i)
           enddo
           print*, "tensor w+ to www+:", sumtest
        else
           print*, ""
           call WWW_atoww1_mg2_had(v(0,1,1),1,434,mg2_tens,NW)           !A -> ve  e+  mu-  vm~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_aww(i,1,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor a to ww 1:", sumtest
           call WWW_atoww3_mg2_had(v(0,1,1),1,434,mg2_tens,NW)           !A -> ve  e+  tau- vt~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_aww(i,1,3)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor a to ww 3:", sumtest
           call WWW_ztoww1_mg2_had(v(0,1,1),1,434,mg2_tens,NW)           !Z -> ve  e+  mu-  vm~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_zww(i,1,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor z to ww 1:", sumtest
           call WWW_ztoww3_mg2_had(v(0,1,1),1,434,mg2_tens,NW)           !Z -> ve  e+  tau- vt~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (WWW_zww(i,1,3)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor z to ww 3:", sumtest
           call WWW_wmtowww_mg2_had(v(0,1,1),1,1,mg2_tens,NW)              !W- -> ve e+ mu- vm~ ta- vt~
           sumtest = 0
           do i = 0,5
             sumtest = sumtest + abs( (wmtowww(i,1)+1d-200) / (mg2_tens(i)+1d-200) - 1)
           enddo
           print*, "tensor w- to www-:", sumtest
        endif
        read(*,*)
      endif
#endif


c call PDF subroutine in order to determine parton
c distributions in the incoming (anti)protons.
      q_sf = sqrt(mufsq(1,1))

      call pdfproton( xi(1), q_sf, pdf(-6,1) )        ! f_a(x1)
      call pdfproton( xi(2), q_sf, pdf(-6,2) )        ! f_a(x2)

C******************************************************************************

      nmax = FL_VVV(fsign,-1)            !reset counter for subprocesses to 0

c******************************************************************************

C*******************  q1 q2bar ---> W W W   ********************************
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

      call qq_WWW(xi,p,v,physToDiag,fsign,bos,nlo,
     1            pdf,res,nmin,nmax,NW)

C*******************  q2bar q1 ---> W W W   ********************************

      physToDiag(1)=2    !physToDiag(1/2) are labels of incoming quarks
      physToDiag(2)=1

      fsign(1) = 1
      fsign(2) =-1

      call qq_WWW(xi,p,v,physToDiag,fsign,bos, nlo,
     1            pdf,res,nmin,nmax,NW)

c*****************  end of process evaluation  ********************************

!     reverse momenta switch for W-W+W- (in order for a correct event output etc...)
      if(bos.eq.434) then
         do mu=0,3
         swmomenta(mu,1) = v(mu,3,1) ! ve
         swmomenta(mu,2) = v(mu,4,1) ! e+
         swmomenta(mu,3) = v(mu,1,1) ! mu-
         swmomenta(mu,4) = v(mu,2,1) ! vm~
         enddo
         do i=1,4
            do mu=0,3
                v(mu,i,1) = swmomenta(mu,i)
            enddo
         enddo
      endif

      SELECT CASE(process)
      CASE(WPhadWMWP,WMhadWPWM)
        call hadDecayFactor_W(finalquarks(1), abs(mjj2(v(0,1,1),v(0,2,1))), N_gen_W, fac_W)
        DO IRES = 1,NMAX
          res(IRES) = res(IRES) * fac_W            ! factor fac_W for hadronic decay
        ENDDO
      CASE(WPWMhadWP,WMWPhadWM)
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
            CASE(WPhadWMWP, WMWPhadWM)
               finalquarks_psp(1) =  2 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -1 - 2* INT(rnumb*2)
            CASE(WMhadWPWM, WPWMhadWP)
               finalquarks_psp(1) =  1 + 2* INT(rnumb*2)
               finalquarks_psp(2) = -2 - 2* INT(rnumb*2)
            END SELECT
         endif
      endif

      m2s_qqWWW = m2s

      RETURN
      END

c******************************************************************************
c
c    end function m2s_qqWWW
c
c******************************************************************************

c******************************************************************************
c
c   begin subroutine qq_WWW
c
c******************************************************************************
c
c       Stefan Prestel, <sprestel@particle.uni-karlsruhe.de>
c	Initial version:  
c	Last modified by Vera Hankele:  2009 January 28
c       
c       qq_WWW calls the amplitude square routines 
c             qqWWW   for qq --> WWW --> 6 leptons
c
c  INPUT:  xi(1:2)         momentum fraction xi of external partons
c          p(0:3,4)        external physical parton momenta
c          v(0:3,nv)       Z decay momenta, for Higgs production only
c                          the sum q(mu) = v(mu,1)+...+v(mu,nv) is needed
c          physToDiag(4)   physToDiag(ext.mom. label) = Feynman diagram label
c          fsign           sign factors for amplitude calls; see qqZqq
c          bos             XXX Boson identifier,X=1 photon, X=2 Z, X=3,4 W+-
c          nlo             nlo = 0: calculate LO only
c                          nlo = 1: calculate full NLO (Born + interference)
c          pdf             complete pdf's as computed in pdfproton
c          NW              NW=1: Narrow width approximation,
c                          NW=0: full calculation
c  OUTPUT:
c          uuvv(ID)   output in format as in qqZqq(j), but corresponding
c          etc.            to res(*) in m2s_qqZqq
c     or   res(k)   
c
c  In and Output
c    nmin, nmax            range of process numbers for this call
c
c*****************************************************************************

      SUBROUTINE qq_WWW(xi,p,v,physToDiag,fsign,bos, nlo,
     1                  pdf,res,nmin,nmax,NW)
c

      IMPLICIT NONE

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"

      double precision p(0:3,n_p), v(0:3,n_v), xi(nx)
      double precision pdf(-6:6,2), res(maxnumsubproc)
      integer physToDiag(2), fsign(2+n_v), bos, nlo, nmin, nmax,NW

      double precision N                              ! color factors
      parameter(N=3d0)

c alfas, scales etc
#include "VBFNLO/utilities/scales.inc"
#include "tensor.inc"
#include "VBFNLO/utilities/coupl.inc"

c local variables
      double precision udbarWWW, udbarWWW2                 !for process pp --> W+W-W+ --> Leptons
      double precision dubarWWW, dubarWWW2                 !for process pp --> W+W-W- --> Leptons

      double precision pbar(0:3,2+n_v), dotrr, polcol
      double precision Matrixelt(0:1,0:1)
      double complex zero
      integer iflav(2), FL_VVV
      external dotrr, FL_VVV


      LOGICAL ldebug, ldebug_ltens
      PARAMETER (ldebug=.false.)
      PARAMETER (ldebug_ltens=.false.)
      parameter (zero = (0d0,0d0) )
      integer i,if1,k,mu
      integer stopp
      data stopp /1/
      save stopp
      integer counter                   ! to count number of M^2 differing to madgraph by epsilon
      data counter/0/
      save counter

      integer debugpspoints              ! number of phasespace points tested
      parameter (debugpspoints = 2)
      double precision epsilon         ! parameter in which madgraph and handwritten M^2 agree
      parameter ( epsilon = 1.0d-8)

c initialize output

      udbarWWW = 0.0d0
      dubarWWW = 0.0d0

      udbarWWW2 = 0.0d0
      dubarWWW2 = 0.0d0

      nmin = nmax+1

      polcol = 1d0/(4d0*N**2*xi(1)*xi(2))

      do mu = 0,3
         do i = 1,2
            pbar(mu,physToDiag(i))=p(mu,i)
         enddo
      enddo

      do mu = 0,3   ! kinematics for w+w-w+ | w-w+w- --> 6 lepton decay
         pbar(mu,3) = v(mu,1) ! ve    | ve
         pbar(mu,4) = v(mu,2) ! e+    | e+
         pbar(mu,5) = v(mu,3) ! mu-   | mu-
         pbar(mu,6) = v(mu,4) ! vm~   | vm~
         pbar(mu,7) = v(mu,5) ! vt    | ta-
         pbar(mu,8) = v(mu,6) ! tau+  | vt~
      enddo
c choose process via boson identifier bos

      if(bos.eq.343) then
#ifdef WITH_MADGRAPH  /* Madgraph code for u d~ -> ve e+ mu- vm~ vt ta+ */
         if (ldebug) then
             call sudbar_WWW(pbar,fsign,udbarWWW2,NW)
         endif
#endif

#ifdef WITH_KK
         if (.not. mgtest) then
           call QQWWWp(pbar,fsign,nlo,1,udbarWWW,NW)
         else
           call sudbar_WWW_mhlm(pbar,fsign,udbarWWW)     ! Madgraph code for u d~ -> ve e+ mu- vm~ vt ta+
         endif
#else
      call QQWWWp(pbar,fsign,nlo,1,udbarWWW,NW)
#endif

         Matrixelt(0,1)=udbarWWW

         do if1=1,(nfl/2)*2,2
            iflav(1)=( if1 - PhysToDiag(1) + 2 )*fsign(PhysToDiag(1))
            iflav(2)=( if1 - PhysToDiag(2) + 2 )*fsign(PhysToDiag(2))
            k=FL_VVV(iflav,1)

c sign1, sign2 read from global.inc: identifier of hadrons
c        in colliding beams 1 and 2 (proton or antiproton)

            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *Matrixelt(mod(if1+1,2),mod(if1,2))*polcol
         enddo

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
            call WWW_atoww1_mg2_had(v(0,1),1,343,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
            call WWW_atoww2_mg2_had(v(0,1),1,343,WWW_aww(0,1,2),NW)           !A -> mu- vm~ vt   tau+
            call WWW_ztoww1_mg2_had(v(0,1),1,343,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
            call WWW_ztoww2_mg2_had(v(0,1),1,343,WWW_zww(0,1,2),NW)           !Z -> mu- vm~ vt   tau+
            call WWW_wptowww_mg2_had(v(0,1),1,1,wptowww(0,1),NW)              !W+ -> ve e+ mu- vm~ vt ta+
            call QQWWWp(pbar,fsign,nlo,1,udbarWWW2,NW)
            call WWW_atoww1_had(v(0,1),1,343,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
            call WWW_atoww2_had(v(0,1),1,343,WWW_aww(0,1,2),NW)           !A -> mu- vm~ vt   tau+
            call WWW_ztoww1_had(v(0,1),1,343,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
            call WWW_ztoww2_had(v(0,1),1,343,WWW_zww(0,1,2),NW)           !Z -> mu- vm~ vt   tau+
            call WWW_wptowww_had(v(0,1),1,1,wptowww(0,1),NW)              !W+ -> ve e+ mu- vm~ vt ta+
            print*,"udbarWWW =  ",udbarWWW
            print*,"udbarWWW2 = ",udbarWWW2
            print*,"udbarWWW/udbarWWW2 - 1 = ",abs(udbarWWW/udbarWWW2 - 1)
            read(*,*)
         endif
         if (ldebug) then
            if (abs(udbarWWW/udbarWWW2 - 1) .gt. epsilon) then
               print*," "
               print*,"udbarWWW =  ",udbarWWW
               print*,"udbarWWW2 = ",udbarWWW2
               print*,"udbarWWW/udbarWWW2 - 1 = ",udbarWWW/udbarWWW2 - 1
               print*," "
               counter = counter +1
            endif
            if (stopp .eq. debugpspoints) then
               print*,counter, " points in ",debugpspoints,"
     1                differ by number greater than ",epsilon
               print*,"Leaving VBFNLO"
               stop
            endif
            stopp = stopp +1
         endif
#endif

      else if(bos.eq.434) then   !start of ubard-comp. 

#ifdef WITH_MADGRAPH  /* !Madgraph du~ -> ve e+ mu- vm~ ta- vt~ */
         if(ldebug) then
            call sdubar_WWW(pbar,fsign,dubarWWW2,NW)
         endif
#endif

#ifdef WITH_KK
      if (.not. mgtest) then
         call QQWWWm(pbar,fsign,nlo,1,dubarWWW,NW)
      else
         call sdubar_WWW_mhlm(pbar,fsign,dubarWWW2,NW)
      endif
#else
         call QQWWWm(pbar,fsign,nlo,1,dubarWWW,NW)
#endif

         Matrixelt(1,0)=dubarWWW

         do if1=1,(nfl/2)*2,2
            iflav(1)=( if1 + PhysToDiag(1) - 1 )*fsign(PhysToDiag(1))
            iflav(2)=( if1 + PhysToDiag(2) - 1 )*fsign(PhysToDiag(2))
            k=FL_VVV(iflav,1)

            res(k)=pdf(sign1*iflav(1),1)
     &            *pdf(sign2*iflav(2),2)
     &            *Matrixelt(mod(if1,2),mod(if1+1,2))*polcol
         enddo

#ifdef WITH_MADGRAPH
      if (ldebug_ltens) then
            call WWW_atoww1_mg2_had(v(0,1),1,434,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
            call WWW_atoww3_mg2_had(v(0,1),1,434,WWW_aww(0,1,3),NW)           !A -> ve  e+  tau- vt~
            call WWW_ztoww1_mg2_had(v(0,1),1,434,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
            call WWW_ztoww3_mg2_had(v(0,1),1,434,WWW_zww(0,1,3),NW)           !Z -> ve  e+  tau- vt~
            call WWW_wmtowww_mg2_had(v(0,1),1,1,wmtowww(0,1),NW)              !W- -> ve e+ mu- vm~ ta- vt~
            call QQWWWm(pbar,fsign,nlo,1,dubarWWW2,NW)
            call WWW_atoww1_had(v(0,1),1,434,WWW_aww(0,1,1),NW)           !A -> ve  e+  mu-  vm~
            call WWW_atoww3_had(v(0,1),1,434,WWW_aww(0,1,3),NW)           !A -> ve  e+  tau- vt~
            call WWW_ztoww1_had(v(0,1),1,434,WWW_zww(0,1,1),NW)           !Z -> ve  e+  mu-  vm~
            call WWW_ztoww3_had(v(0,1),1,434,WWW_zww(0,1,3),NW)           !Z -> ve  e+  tau- vt~
            call WWW_wmtowww_had(v(0,1),1,1,wmtowww(0,1),NW)              !W- -> ve e+ mu- vm~ ta- vt~
            print*,"dubarWWW =  ",dubarWWW
            print*,"dubarWWW2 = ",dubarWWW2
            print*,"dubarWWW/dubarWWW2 - 1 = ",abs(dubarWWW/dubarWWW2 - 1)
            read(*,*)
         endif
         if (ldebug) then
            if (abs(dubarWWW/dubarWWW2 - 1) .gt. epsilon) then
               print*," "
               print*,"dubarWWW =  ",dubarWWW
               print*,"dubarWWW2 = ",dubarWWW2
               print*,"dubarWWW/dubarWWW2 - 1 = ",dubarWWW/dubarWWW2 - 1
               print*," "
               counter = counter +1
            endif
            if (stopp .eq. debugpspoints) then
               print*,counter, " points in ",debugpspoints,
     1              " differ by number greater than ",epsilon
               print*,"Leaving VBFNLO"
               stop
            endif
            stopp = stopp +1
         endif
#endif

      else
         write(6,*) "Process not implemented, check boson identifier"
      endif !end of dubar-computation block

      nmax = FL_VVV(iflav,-2)

      end
c******************************************************************
c
c   end subroutine qq_WWW
c
c*****************************************************************

c*****************************************************************************
c
c    begin subroutine vtowww_reset
c
c*****************************************************************************
c
c       Stefan Prestel, <sprestel@particle.uni-karlsruhe.de>
c	Initial version:  
c	Last modified by Vera Hankele:  2009 January 28
c       
c  this subroutines sets all leptonic tensors needed in the calculation 
c  of pp -> WWW to zero.
c  
c*****************************************************************************

      SUBROUTINE vtowww_reset

      IMPLICIT NONE

#include "tensor.inc"
      double complex zero
      parameter (zero=(0d0,0d0))
      integer j,jj,mu

      do mu = 0,5
         do j = 1,3
            do jj =1,3
              WWW_aww(mu,j,jj)  = zero
              WWW_zww(mu,j,jj)  = zero
            enddo
            wptowww(mu,j) = zero
            wmtowww(mu,j) = zero
         enddo
      enddo
      print*,' precalculated tensors initialized to ',zero 
      return
      end
c*****************************************************************************
c
c    end subroutine vtowww_reset
c
c*****************************************************************************
