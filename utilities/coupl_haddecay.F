! This subroutine sets all couplings relevant for hadronic boson decays

      SUBROUTINE coupl_haddecay

      implicit none

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"

      ! set all couplings as for leptonic decays
      GWP_ZUN(1) = GZN(1)
      GWP_ZUN(2) = GZN(2)
      GWP_ZDL(1) = GZL(1)
      GWP_ZDL(2) = GZL(2)
      GWP_AUN(1) = dcmplx(0.)
      GWP_AUN(2) = dcmplx(0.)
      GWP_ADL(1) = GAL(1)
      GWP_ADL(2) = GAL(2)
      GWM_ZUN(1) = GZN(1)
      GWM_ZUN(2) = GZN(2)
      GWM_ZDL(1) = GZL(1)
      GWM_ZDL(2) = GZL(2)
      GWM_AUN(1) = dcmplx(0.)
      GWM_AUN(2) = dcmplx(0.)
      GWM_ADL(1) = GAL(1)
      GWM_ADL(2) = GAL(2)
      GZ_ZF(1) = GZL(1)
      GZ_ZF(2) = GZL(2)
      GZ_ZUN(1) = GZN(1)
      GZ_ZUN(2) = GZN(2)
      GZ_AF(1) = GAL(1)
      GZ_AF(2) = GAL(2)
      Ztouu = .false.
      Ztodd = .false.

      SELECT CASE(procid)
      CASE(WPhadWMlep,GFWPhadWMlep,WPhadZ,
     &     WPhadWMjj, WPhadZjj, WPhadWMZ, ZZWPhad, WPhadWMWP, WMWPhadWM, WPhadZJ,
     &     Hjj_WPhadWM, WPhadWPjj,
     &     WMhadWMjj,
     &     WPhadWMA, WPhadZA,
     &     WPhadWMj,GFWPhadWMj)  ! WMWM is calculated using WPWP
        ! set couplings for W+ hadronic decay
        GWP_ZUN(1) = GZU(1)
        GWP_ZUN(2) = GZU(2)
        GWP_ZDL(1) = GZD(1)
        GWP_ZDL(2) = GZD(2)
        GWP_AUN(1) = GAU(1)
        GWP_AUN(2) = GAU(2)
        GWP_ADL(1) = GAD(1)
        GWP_ADL(2) = GAD(2)
      CASE(WPlepWMhad,GFWPlepWMhad,WMhadZ,
     &     WPWMhadjj, WMhadZjj, WPWMhadZ, ZZWMhad, WPWMhadWP, WMhadWPWM, WMhadZJ,
     &     Hjj_WPWMhad,
     &     WPWMhadA, WMhadZA,
     &     WPWMhadj,GFWPWMhadj)
        ! set couplings for W- hadronic decay
        GWM_ZUN(1) = GZU(1)
        GWM_ZUN(2) = GZU(2)
        GWM_ZDL(1) = GZD(1)
        GWM_ZDL(2) = GZD(2)
        GWM_AUN(1) = GAU(1)
        GWM_AUN(2) = GAU(2)
        GWM_ADL(1) = GAD(1)
        GWM_ADL(2) = GAD(2)
      CASE(ZZhad,WPZhad,WMZhad,
     &     WPZhadjj, WMZhadjj, ZZhadjj, WWZhad, ZZhadWP, ZZhadWM, ZZZhad, WPZhadJ, WMZhadJ,
     &     Hjj_ZZhad,
     &     ZZhadA, WPZhadA, WMZhadA)
        if (finalquarks(1).eq.93 .or. mod(abs(finalquarks(1)),2).eq.0) then         ! the mod-part also includes ID=94!
           ! set couplings for Z hadronic decay into uubar
           ! any modifications here also have to go into the relevant m2s-routines!
           GZ_ZF(1) = GZU(1)
           GZ_ZF(2) = GZU(2)
           GZ_AF(1) = GAU(1)
           GZ_AF(2) = GAU(2)
           Ztouu = .true.
        else
           ! set couplings for Z hadronic decay into ddbar
           ! any modifications here also have to go into the relevant m2s-routines!
           GZ_ZF(1) = GZD(1)
           GZ_ZF(2) = GZD(2)
           GZ_AF(1) = GAD(1)
           GZ_AF(2) = GAD(2)
           Ztodd = .true.
        endif
        GZ_ZUN(1) = GZU(1)
        GZ_ZUN(2) = GZU(2)
      END SELECT

      RETURN
      END


      SUBROUTINE setZtouu

      implicit none

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

      GZ_ZF(1) = GZU(1)
      GZ_ZF(2) = GZU(2)
      GZ_AF(1) = GAU(1)
      GZ_AF(2) = GAU(2)
      Ztouu = .true.
      Ztodd = .false.

      RETURN
      END


      SUBROUTINE setZtodd

      implicit none

#include "VBFNLO/utilities/coupl.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

      GZ_ZF(1) = GZD(1)
      GZ_ZF(2) = GZD(2)
      GZ_AF(1) = GAD(1)
      GZ_AF(2) = GAD(2)
      Ztouu = .false.
      Ztodd = .true.

      RETURN
      END


      subroutine hadDecayFactor_Z(finalquarks_in, QSQ_had, N_gen_up, N_gen_down, fac_Z_up, fac_Z_down)

      ! This subroutine checks whether the Q^2 of the hadronically decaying Z / gamma
      ! is large enough to produce the desired final state.
      ! Furthermore, the color factor for hadronic Z decays is applied here.
      ! Alpha_S corrections will be applied here as well in the future.

      ! input:
      ! * finalquarks_in: particle ID of the first decay product
      ! * QSQ_had: Q^2 of the hadronically decaying boson
      !   Q^2 of vector boson is the same in all tilde cinematics, so use this globally 
      !   for one phase space point
      ! output: 
      ! * N_gen_up:   contributing number of up-type quarks to Z decay
      ! * N_gen_down: contributing number of down-type quarks to Z decay
      ! * N_gen_W:    contributing number of generations to W decay
      ! * fac_*:      N_gen * color factor (* alpha_s corrections)

      IMPLICIT NONE

      integer finalquarks, finalquarks_in, N_gen_up, N_gen_down
      double precision QSQ_had, sqrt_s_had
      double precision fac_Z_up, fac_Z_down
      double precision ud_threshold, s_threshold, c_threshold, b_threshold

      integer pdflib, valueLHA
      COMMON/PDFparameters/pdflib, valueLHA

      double precision alphas5, alphas5_hardwired
      external alphas5, alphas5_hardwired

      parameter (
     &  b_threshold  = 8.77  ,
     &  c_threshold  = 2.56  ,
     &  s_threshold  = 1.13 ,
     &  ud_threshold = 0.373  )

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

      finalquarks = abs(finalquarks_in)
      sqrt_s_had = sqrt(QSQ_had)


      ! Is Z decay possible?
      N_gen_up   = 0
      N_gen_down = 0
      select case (finalquarks)
      case (1)
         if (sqrt_s_had.ge.ud_threshold) N_gen_down = 1
      case (2)
         if (sqrt_s_had.ge.ud_threshold) N_gen_up   = 1
      case (3)
         if (sqrt_s_had.ge.s_threshold)  N_gen_down = 1
      case (4)
         if (sqrt_s_had.ge.c_threshold)  N_gen_up   = 1
      case (5)
         if (sqrt_s_had.ge.b_threshold)  N_gen_down = 1
      case (93)
         if (sqrt_s_had.ge.c_threshold) then
            N_gen_up   = 2
            N_gen_down = 2
         elseif (sqrt_s_had.ge.s_threshold) then
            N_gen_up   = 1
            N_gen_down = 2
         elseif (sqrt_s_had.ge.ud_threshold) then
            N_gen_up   = 1
            N_gen_down = 1
         endif
      case (94)
         if (sqrt_s_had.ge.b_threshold) then
            N_gen_up   = 2
            N_gen_down = 3
         elseif (sqrt_s_had.ge.c_threshold) then
            N_gen_up   = 2
            N_gen_down = 2
         elseif (sqrt_s_had.ge.s_threshold) then
            N_gen_up   = 1
            N_gen_down = 2
         elseif (sqrt_s_had.ge.ud_threshold) then
            N_gen_up   = 1
            N_gen_down = 1
         endif
      case default
         print*, "Error: Wrong final state declaration."
         STOP
      end select


      ! include factor 3 from color into fac_*
      fac_Z_up   = N_gen_up   * 3.d0
      fac_Z_down = N_gen_down * 3.d0


      ! include NLO decay corrections here
      if (nlo_semilep.gt.0) then
        if (pdflib.eq.1 .and. .not.doNLO) then   ! lhapdf used and LO production
          fac_Z_up   = fac_Z_up   * (1.d0 + alphas5_hardwired(max(QSQ_had,1d0),1)/pi)
          fac_Z_down = fac_Z_down * (1.d0 + alphas5_hardwired(max(QSQ_had,1d0),1)/pi)
        else
          fac_Z_up   = fac_Z_up   * (1.d0 + alphas5(max(QSQ_had,1d0),1)/pi)
          fac_Z_down = fac_Z_down * (1.d0 + alphas5(max(QSQ_had,1d0),1)/pi)
        endif
      endif

      return
      end


      subroutine hadDecayFactor_W(finalquarks_in, QSQ_had, N_gen_W, fac_W)

      ! This subroutine applies the the color factor for hadronic W decays.
      ! Alpha_S corrections will be applied here as well in the future.

      ! input:
      ! * finalquarks_in: particle ID of the first decay product
      ! * QSQ_had: Q^2 of the hadronically decaying boson
      !   Q^2 of vector boson is the same in all tilde cinematics, so use this globally 
      !   for one phase space point
      ! output: 
      ! * N_gen_W:    contributing number of generations to W decay
      ! * fac_*:      N_gen * color factor (* alpha_s corrections)

      IMPLICIT NONE

      integer finalquarks, finalquarks_in, N_gen_W
      double precision QSQ_had, sqrt_s_had
      double precision fac_W

      integer pdflib, valueLHA
      COMMON/PDFparameters/pdflib, valueLHA
      
      double precision alphas5, alphas5_hardwired
      external alphas5, alphas5_hardwired

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"

      finalquarks = abs(finalquarks_in)
      sqrt_s_had = sqrt(QSQ_had)


      ! How many generations in W decay?
      select case (finalquarks)
      case (93,94)
         N_gen_W = 2    ! Assume W decay always possible
      case default
         N_gen_W = 1    ! Assume W decay always possible
      end select


      ! include factor 3 from color into fac_W
      fac_W = N_gen_W * 3.d0


      ! include NLO decay corrections here
      if (nlo_semilep.gt.0) then
        if (pdflib.eq.1 .and. .not.doNLO) then   ! lhapdf used and LO production
          fac_W = fac_W * (1.d0 + alphas5_hardwired(max(QSQ_had,1d0),1)/pi)
        else
          fac_W = fac_W * (1.d0 + alphas5(max(QSQ_had,1d0),1)/pi)
        endif
      endif

      return
      end


