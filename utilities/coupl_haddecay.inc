c============================================================================
c
c  Define common block containing all coupling constants for hadronic decays.
c
c============================================================================

      double complex GWP_ZUN(2), GWP_ZDL(2), GWP_AUN(2), GWP_ADL(2)
      double complex GWM_ZUN(2), GWM_ZDL(2), GWM_AUN(2), GWM_ADL(2)
      double complex GZ_ZUN(2),  GZ_ZF(2),  GZ_AF(2)
      logical Ztodd, Ztouu

      common /COUPL_HAD/  GWP_ZUN, GWP_ZDL, GWP_AUN, GWP_ADL,
     &                    GWM_ZUN, GWM_ZDL, GWM_AUN, GWM_ADL,
     &                    GZ_ZUN,  GZ_ZF,  GZ_AF

      common /COUPL_HAD_Z/ Ztodd, Ztouu

c============================================================================
c
c  Define common block containing for approx. NLO decay of vector bosons.
c
c============================================================================

      integer nlo_semilep 

      common /NLO_SEMILEP_DECAY/ nlo_semilep

