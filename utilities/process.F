c*************************************************************************  
c     This file contains subroutines for the process initialization.
c
c     It also contains subroutines for assigning the right
c     process identifier, the number of final state leptons and default
c     values for the final state leptons
c     In additon it contains routines to order the leptons in the final state
c     and to assign the right symmetryfactor  
c*************************************************************************  
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE InitProcess
c     SUBROUTINE PrintProcInfo
c     SUBROUTINE InitRealEmissions
c     SUBROUTINE InitQEDEmissions
c     SUBROUTINE proc_assignment(procID,process,N)
c     SUBROUTINE lepton_assignment(N,finallep)
c     SUBROUTINE symmetry_factor
c     SUBROUTINE lepton_gen
c
c*************************************************************************  


c*************************************************************************
      SUBROUTINE InitProcess
c*************************************************************************
c     This routine is intended to do the initialization of parameters that
c     are process dependent. However process dependent coupling parameters
c     should be set in InitCouplings, or shouldn"t they?
c*************************************************************************
          use globalvars, only: ldoblha
      implicit none
#include "global.inc"
#include "process.inc"
      integer lflavr(5:6), i

      photoncut_vbf = .false.       ! print out qsqAmin cut from global.inc

c first set phase space dimension & NLO parameters
      SELECT CASE(process)

#ifdef WITH_VBF 
***********************
* VBF Higgs processes *
***********************
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar)
         n_p = 4
         if (ldoblha .and. (process.eq.Hjj)) then
           n_v = 1
         else
           n_v = 2
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         if (ewcor_switch) NLO_Loops = 2
         N_points(1) = N_points(0)        ! #points for virtual corr.
         N_points(2) = N_points(1)-3      ! #points for ew corns 
         call qqhqqgi(0)
         if (with_anomHiggs) photoncut_vbf = .true.

      CASE(AAjj)        ! VBF pp -> spin-2 (or Higgs) jj -> aa jj
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         photoncut_vbf = .true.
      CASE(Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)        ! VBF pp -> spin-2 jj -> VVjj
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr. 
      CASE(Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu,
     &     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad) ! VBF pp -> Hjj -> VVjj
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         if (ewcor_switch) NLO_Loops = 2
         N_points(1) = N_points(0)        ! #points for virtual corr. 
         N_points(2) = N_points(1)-3      ! #points for ew corns
         call qqhqqgi(0)
         if (with_anomHiggs) photoncut_vbf = .true.
      CASE(HjjjLO, HjjjLO_AA, HjjjLO_mu, HjjjLO_tau, HjjjLO_bbar)
         n_p = 5
         if (ldoblha .and. (process.eq.HjjjLO)) then
           n_v = 1
         else
           n_v = 2
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
         call qqhqqgi(0)
         if (with_anomHiggs) photoncut_vbf = .true.
      CASE(AAjjjLO)         ! VBF pp -> spin-2 (or Higgs) jjj -> aa jjj
         n_p = 5
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
         photoncut_vbf = .true.
      CASE(Sp2jjjLO_WW, Sp2jjjLO_ZZ_ll, Sp2jjjLO_ZZ_lnu)    ! VBF pp -> spin-2 jjj -> VV jjj
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
      CASE(HjjjLO_WW, HjjjLO_ZZ_ll, HjjjLO_ZZ_lnu,
     &     HjjjLO_WPhadWM, HjjjLO_WPWMhad, HjjjLO_ZZhad)
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0
         call qqhqqgi(0)
         if (with_anomHiggs) photoncut_vbf = .true.
      CASE(HHjj)
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
      CASE(HHjjjLO)
         n_p = 5
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false. ! process at lo
         NLO_Loops = 0
      CASE(HHjj_bbtau,HHjj_bbAA)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
      CASE(HHjjjLO_bbtau,HHjjjLO_bbAA)
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false. ! process at lo
         NLO_Loops = 0

********************
* VBF WW processes *
********************
      CASE(WPWMjj, WPhadWMjj, WPWMhadjj, WPWPjj, WPhadWPjj)
         n_p = 4  ! actual number of partons
         n_v = 4  ! actual number of other particles
         PS_dimension = 3*(n_p+n_v-nx)-3  ! phasespace dimension
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 3  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr. 
         N_points(2) = N_points(0) - 4    ! -||- for boxes
         N_points(3) = N_points(0) - 7    ! -||- for pentagons
         call qqhqqi(3)
         if (process.ne.WPWPjj .and. process.ne.WPhadWPjj) photoncut_vbf = .true.
      CASE(WMWMjj, WMhadWMjj)
         n_p = 4  ! actual number of partons
         n_v = 4  ! actual number of other particles
         PS_dimension = 3*(n_p+n_v-nx)-3  ! phasespace dimension
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 3  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr. 
         N_points(2) = N_points(0) - 4    ! -||- for boxes
         N_points(3) = N_points(0) - 7    ! -||- for pentagons
         call qqhqqi(3)
         sign1 = -sign1    ! get W-W-jj from W+W+jj by switching p <--> pbar
         sign2 = -sign2
      CASE(WPWMjjjLO, WPhadWMjjjLO, WPWMhadjjjLO, WPWPjjjLO, WPhadWPjjjLO)
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         call qqhqqi(3)
         if (process.ne.WPWPjjjLO .and. process.ne.WPhadWPjjjLO) photoncut_vbf = .true.
      CASE(WMWMjjjLO, WMhadWMjjjLO)
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         call qqhqqi(3)
         sign1 = -sign1    ! get W-W-jj from W+W+jj by switching p <--> pbar
         sign2 = -sign2
********************
* VBF ZZ processes *
********************
      CASE(ZZjj_ll, ZZjj_lnu, ZZhadjj)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 7
         photoncut_vbf = .true.
      CASE(ZZjjjLO_ll, ZZjjjLO_lnu, ZZhadjjjLO)  
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         photoncut_vbf = .true.
********************
* VBF W+Z processes *
********************
      CASE(WPZjj,WPhadZjj,WPZhadjj)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 7
         photoncut_vbf = .true.
      CASE(WPZjjjLO,WPhadZjjjLO,WPZhadjjjLO)  
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         photoncut_vbf = .true.
********************
* VBF W-Z processes *
********************
      CASE(WMZjj,WMhadZjj,WMZhadjj)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 7
         photoncut_vbf = .true.
      CASE(WMZjjjLO,WMhadZjjjLO,WMZhadjjjLO)  
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
         photoncut_vbf = .true.
********************
* VBF W+A processes *
********************
      CASE(WPAjj)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 7
      CASE(WPAjjjLO)  
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
********************
* VBF W-A processes *
********************
      CASE(WMAjj)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 7
      CASE(WMAjjjLO)  
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
********************
* VBF ZA processes *
********************
      CASE(ZAjj)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 5
         photoncut_vbf = .true.
      CASE(ZAjjjLO)  
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         photoncut_vbf = .true.
      CASE(ZAjj_n)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) - 4
         N_points(3) = N_points(0) - 5
         photoncut_vbf = .true.
      CASE(ZAjjjLO_n)  
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0
         photoncut_vbf = .true.
*******************
* VBF Z processes *
*******************
      CASE(Zjj_l)
         n_p = 4
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 1 
         N_points(1) = N_points(0) 
         lflavr(5) = 2          
         lflavr(6) = 2          
         call qqbqqi(2,lflavr)
         photoncut_vbf = .true.
      CASE(Zjj_nu)
         n_p = 4
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 1 
         N_points(1) = N_points(0) 
         lflavr(5) = 1          
         lflavr(6) = 1          
         call qqbqqi(2,lflavr)
         photoncut_vbf = .true.
      CASE(ZjjjLO_l)
         n_p = 5
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         lflavr(5) = 2          
         lflavr(6) = 2          
         call qqbqqi(2,lflavr)
         photoncut_vbf = .true.
      CASE(ZjjjLO_nu)
         n_p = 5
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         lflavr(5) = 1          
         lflavr(6) = 1          
         call qqbqqi(2,lflavr)
         photoncut_vbf = .true.
***********************
* VBF W+/W- processes *
***********************
      CASE(WPjj) ! VBF pp->W+jj->l+nu
         n_p = 4
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.true.
         NLO_Loops = 1 
         N_points(1) = N_points(0) 
         lflavr(5) = 2          
         lflavr(6) = 1          
         call qqbqqi(3,lflavr)
         photoncut_vbf = .true.
      CASE(WMjj) ! VBF pp->W-jj->l-nu
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 1
         N_points(1) = N_points(0)
         lflavr(5) = 1
         lflavr(6) = 2
         call qqbqqi(4,lflavr)
         photoncut_vbf = .true.
      CASE(WPjjjLO)
         n_p = 5
         n_v = 2  
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.
         NLO_Loops = 0 
         lflavr(5) = 2          
         lflavr(6) = 1          
         call qqbqqi(3,lflavr)
         photoncut_vbf = .true.
      CASE(WMjjjLO)
         n_p = 5
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
         lflavr(5) = 1
         lflavr(6) = 2
         call qqbqqi(4,lflavr)
         photoncut_vbf = .true.
*******************
* VBF A processes *
*******************
      CASE(Ajj)
         n_p = 4
         n_v = 1
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         lflavr(5) = 2  ! lflavr is not used at all in qqbqqi for ibos=1 
         lflavr(6) = 2
         call qqbqqi(1,lflavr)
         photoncut_vbf = .true.
      CASE(AjjjLO)
         n_p = 5
         n_v = 1
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
         lflavr(5) = 2   ! lflavr is not used at all in qqbqqi for ibos=1
         lflavr(6) = 2 
         call qqbqqi(1,lflavr) 
         photoncut_vbf = .true.
********************************
* VBF Higgs + photon processes *
********************************    
      CASE(HAjj, HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar)
         n_p = 4
         if (ldoblha .and. (process.eq.HAjj)) then
           n_v = 2
         else
           n_v = 3
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         call Get_HAjj_couplings()
      CASE(HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu)
         n_p = 4
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         call Get_HAjj_couplings()    
      CASE(HAjjjLO, HAjjjLO_AA, HAjjjLO_mu, HAjjjLO_tau, HAjjjLO_bbar)
         n_p = 5
         if (ldoblha .and. (process.eq.HAjjjLO)) then
           n_v = 2
         else
           n_v = 3
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at nlo
         NLO_Loops = 0  ! number of different virtual contributions
         call Get_HAjj_couplings()
      CASE(HAjjjLO_WW, HAjjjLO_ZZ_ll, HAjjjLO_ZZ_lnu)
         n_p = 5
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at nlo
         NLO_Loops = 0  ! number of different virtual contributions
         call Get_HAjj_couplings()
#endif


#ifdef WITH_HJJJ
***********************
* VBF Higgs + 3 jets *
***********************
      CASE(Hjjj,Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar)
         n_p = 5
         if (ldoblha .and. (process.eq.Hjjj)) then
           n_v = 1
         else
           n_v = 2
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         call init_qqhqqg_virt
      CASE(Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu) ! VBF pp -> Hjjj -> VVjjj
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 1  ! number of different virtual contributions
         N_points(1) = N_points(0)        ! #points for virtual corr.
         call init_qqhqqg_virt
      CASE(HjjjjLO,HjjjjLO_AA, HjjjjLO_mu, HjjjjLO_tau, HjjjjLO_bbar)
         n_p = 6
         if (ldoblha .and. (process.eq.HjjjjLO)) then
           n_v = 1
         else
           n_v = 2
         endif
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
      CASE(HjjjjLO_WW, HjjjjLO_ZZ_ll, HjjjjLO_ZZ_lnu) ! VBF pp -> Hjjjj ->VVjjjj
         n_p = 6
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3 
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
#endif


#ifdef WITH_GGF
**********************
* GF Higgs processes *
**********************
      CASE(GFHjj)
         n_p = 4
         n_v = 1
         PS_dimension = 3*(n_p-nx)
         hasNLO=.false.
         NLO_Loops = 0
         call init_gluonfusion    
      CASE(GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar) ! GF pp -> Hjj -> xxjj
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at lo
         NLO_Loops = 0  
         call init_gluonfusion
      CASE(GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
         call init_gluonfusion      
#endif


#ifdef WITH_QCDV
*******************
* QCD V processes *
*******************
      CASE(QCDZjj_l, QCDZjj_nu)
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) -3
         N_points(3) = N_points(0)
      CASE(QCDZjjjLO_l,QCDZjjjLO_nu, QCDWPjjjLO,QCDWMjjjLO)
         n_p = 5
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
      CASE(QCDWPjj,QCDWMjj)
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0) -3
         N_points(3) = N_points(0)
#endif


#ifdef WITH_QCDVV
********************
* QCD VV processes *
********************      
      CASE(QCDAAjj)
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) 
         N_points(2) = N_points(0) -5
         N_points(3) = N_points(0) -1  
         print *,"Sorry - this process is not working!"
         stop
      CASE(QCDWPZjj, QCDWMZjj)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) 
         N_points(2) = N_points(0) -5
         N_points(3) = N_points(0)
      CASE(QCDWPAjj, QCDWMAjj)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) 
         N_points(2) = N_points(0) -6
         N_points(3) = N_points(0) -1
      CASE(QCDZAjj_l)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) 
         N_points(2) = N_points(0) -6
         N_points(3) = N_points(0) -1
      CASE(QCDZAjj_n)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) 
         N_points(2) = N_points(0) -6
         N_points(3) = N_points(0) -1
      CASE(QCDZZjj_ll, QCDZZjj_lnu)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0) +0
         N_points(2) = N_points(0) -5
         N_points(3) = N_points(0) -1
         N_points(4) = N_points(0) -10
! !          N_points(5) = N_points(0) -8
      CASE(QCDWPWPjj,QCDWMWMjj)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
	 N_points(1) = N_points(0) 
	 N_points(2) = N_points(0) -4
         N_points(3) = N_points(0) -3
      CASE(QCDWPZjjjLO, QCDWMZjjjLO,QCDWPWPjjjLO,QCDWMWMjjjLO,QCDZZjjjLO_ll,QCDZZjjjLO_lnu)
         n_p = 5
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.false.
         NLO_Loops = 0
      CASE(QCDWPAjjjLO, QCDWMAjjjLO)
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.false.
         NLO_Loops = 0
      CASE(QCDZAjjjLO_l,QCDZAjjjLO_n)
         n_p = 5
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.false.
         NLO_Loops = 0
      CASE(QCDAAjjjLO)
         n_p = 5
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-3+1
         hasNLO=.false.
         NLO_Loops = 0
#endif


#ifdef WITH_DIBOSON
*********************
* W+W- processes *
*********************
      CASE(WW,WPhadWMlep,WPlepWMhad)       ! pp -> W+W-
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         hasFLoops=.true.
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(WWjLO,WPhadWMlepjLO,WPlepWMhadjLO)       ! pp -> W+W- j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
      CASE(GFWW,GFWPhadWMlep,GFWPlepWMhad)       ! pp -> W+W-
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0

****************
* WZ processes *
****************
      CASE(WPZ,WPhadZ,WPZhad)       ! pp -> W+Z
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 1
         N_points(1) = N_points(0)
      CASE(WPZjLOb,WPhadZjLO,WPZhadjLO)       ! pp -> W+Z j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
      CASE(WMZ,WMhadZ,WMZhad)       ! pp -> W-Z
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 1
         N_points(1) = N_points(0)
      CASE(WMZjLOb,WMhadZjLO,WMZhadjLO)       ! pp -> W-Z j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0

****************
* ZZ processes *
****************
      CASE(ZZ,ZZhad)       ! pp -> ZZ
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         hasFLoops=.true.
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(ZZjLO,ZZhadjLO)       ! pp -> ZZ j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
      CASE(GFZZ,GFZZhad)       ! pp -> ZZ 
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0

****************
* WA processes *
****************
      CASE(WPA)       ! pp -> W+A
         n_p = 2
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         NLO_Loops = 1
         N_points(1) = N_points(0)
      CASE(WPAjLOb)       ! pp -> W+A j
         n_p = 3
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0
      CASE(WMA)       ! pp -> W-A
         n_p = 2
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         NLO_Loops = 1
         N_points(1) = N_points(0)
      CASE(WMAjLOb)       ! pp -> W-A j
         n_p = 3
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0


****************
* WH processes *
****************
      !N_points(0): LO number of integration points
      !N_points(1): virtual corrections
      !N_points(-1): real emissions

      ! pp -> W+H

      CASE(
     7      WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar,
     7      WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar 
     &    )
         n_p = 2
         if (ldoblha .and. (process == WPH .or. process == WMH)) then
            n_v = 3
         else
            n_v = 4 !4 for higgsdecay and vbfnlo-internal (fake decay products)
         endif
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.true.
         NLO_LOOPS = 1
         N_points(1) = N_points(0)
      CASE(      
     7      WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7      WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu
     &    )
         n_p = 2
         n_v = 6 
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.true.
         NLO_LOOPS = 1
         N_points(1) = N_points(0)

      ! pp -> W+H j
      CASE(
     7      WPHjLOb, WPHjLOb_AA, WPHjLOb_mu, WPHjLOb_tau, WPHjLOb_bbar,
     7      WMHjLOb, WMHjLOb_AA, WMHjLOb_mu, WMHjLOb_tau, WMHjLOb_bbar
     &    )
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.
         NLO_LOOPS = 0
      CASE(      
     7      WPHjLOb_WW, WPHjLOb_ZZ_ll,  WPHjLOb_ZZ_lnu,
     7      WMHjLOb_WW, WMHjLOb_ZZ_ll,  WMHjLOb_ZZ_lnu
     &    )
         n_p = 3
         n_v = 6 
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.
         NLO_LOOPS = 0

****************
* W processes  *
****************
      CASE(WP_only,WM_only)       ! pp -> W+-
         n_p = 2
         n_v = 2
         PS_dimension = 4 ! 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         NLO_LOOPS = 1
         N_points(1) = N_points(0)
         !N_points(0): LO number of integration points
         !N_points(1): virtual corrections
         !N_points(-1): real emissions
      CASE(WPjLOb,WMjLOb)   ! pp -> W+- j
         n_p = 3
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.
         NLO_LOOPS = 0


****************
* ZA processes *
****************
      CASE(ZA)       ! pp -> ZA
         n_p = 2
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         hasFLoops=.true.
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(ZAjLO)       ! pp -> ZA j
         n_p = 3
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0
      CASE(GFZA)       ! pp -> ZA 
         n_p = 2
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0

****************
* AA processes *
****************
      CASE(AA)       ! pp -> AA
         n_p = 2
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         hasFLoops=.true.
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(AAjLO)       ! pp -> AA j
         n_p = 3
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0
      CASE(GFAA)       ! pp -> AA 
         n_p = 2
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0

#endif


#ifdef WITH_DIBOSONJET
      CASE(WWJ,WPHADWMJ,WPWMHADJ)
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         hasFLoops=.true.
         NLO_Loops = 4
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
         N_points(3) = N_points(0)-6
         N_points(4) = N_points(0)-6
      CASE(WWJJLO,WPHADWMJJLO,WPWMHADJJLO)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
         NLO_Loops = 0
      CASE(GFWWj,GFWPHADWMJ,GFWPWMHADJ)
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
         NLO_Loops = 0
      CASE(WMAJ,WPAJ)
         n_p = 3
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(WMAJLO,WPAJLO)
         n_p = 4
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
      CASE(WMZJ,WMhadZJ,WMZhadJ,WPZJ,WPhadZJ,WPZhadJ)
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
      CASE(WMZJLO,WMhadZJjLO,WMZhadJjLO,WPZJLO,WPhadZJjLO,WPZhadJjLO)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
      CASE(ZZJ)
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.  ! process at nlo
         hasFLoops=.true.
         NLO_Loops = 4
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
         N_points(3) = N_points(0)-6
         N_points(4) = N_points(0)-12
      CASE(ZZJJLO)
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
         NLO_Loops = 0
      CASE(GFZZj)
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.  ! process at LO
         NLO_Loops = 0

*********************
* W+- H j processes *
*********************
      CASE(
     7      WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar,
     7      WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar 
     &    )
         n_p = 3
         if (ldoblha .and. (process == WPHJ .or. process == WMHJ)) then
            n_v = 3
         else
            n_v = 4 !4 for higgsdecay and vbfnlo-internal (fake decay products)
         endif
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.true.  ! process at nlo
         !NLO_Loops = 1
         if(dotoploops) then
            NLO_Loops = 3 ! vertex + box + toploops
         else
            NLO_Loops = 2 
         endif
         N_points(1) = N_points(0)
         N_points(2) = N_points(1)-2 ! real emission points
         ! TODO: set toploop points         
         N_points(3) = N_points(1) 
      CASE(
     7      WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7      WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &    )
         n_p = 3
         n_v = 6 !h->VV->4l + w->2l
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.true.  ! process at nlo
         !NLO_Loops = 1
         if(dotoploops) then
            NLO_Loops = 3 ! vertex + box + toploops
         else
            NLO_Loops = 2 
         endif
         N_points(1) = N_points(0)
         N_points(2) = N_points(1)-2 ! real emission points
         N_points(3) = N_points(1) 


      CASE(
     7      WPHJJLO, WPHJJLO_AA, WPHJJLO_mu, WPHJJLO_tau, WPHJJLO_bbar,
     7      WMHJJLO, WMHJJLO_AA, WMHJJLO_mu, WMHJJLO_tau, WMHJJLO_bbar
     &    )
         n_p = 4
         if (ldoblha .and. (process == WPHJJLO .or. process == WMHJJLO)) then
            n_v = 3
         else
            n_v = 4 !4 for higgsdecay and vbfnlo-internal (fake decay products)
         endif
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.
      CASE(      
     7      WPHJJLO_WW, WPHJJLO_ZZ_ll,  WPHJJLO_ZZ_lnu,
     7      WMHJJLO_WW, WMHJJLO_ZZ_ll,  WMHJJLO_ZZ_lnu
     &    )
         n_p = 4
         n_v = 6 
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.


*********************
* W+- j processes *
*********************
      CASE(WMJ,WPJ)
         n_p = 3
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.true.  ! process at nlo
         NLO_Loops = 2
         N_points(1) = N_points(0)
         N_points(2) = N_points(1)-2 ! real emission points
      CASE(WMJJLO,WPJJLO)
         n_p = 4
         n_v = 2
         PS_dimension = 3*(n_p+n_v-nx-1)
         hasNLO=.false.  ! process at LO
#endif


#ifdef WITH_TRIBOSON
*********************
* W+W-Z processes *
*********************
      CASE(WWZ,WPhadWMZ,WPWMhadZ,WWZhad)       ! pp -> W+W-Z
         n_p = 2
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WWZjLO,WPhadWMZjLO,WPWMhadZjLO,WWZhadjLO)     ! pp -> W+W-Z j
         n_p = 3
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* ZZZ processes *
*********************
      CASE(ZZZ,ZZZhad)       ! pp -> ZZZ
         n_p = 2
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(ZZZjLO,ZZZhadjLO)     ! pp -> ZZZ j
         n_p = 3
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0
*********************
* ZZW processes *
*********************
      CASE(ZZWP, ZZhadWP, ZZWPhad, ZZWM, ZZhadWM, ZZWMhad)       ! pp -> ZZW
         n_p = 2
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(ZZWPjLO, ZZhadWPjLO, ZZWPhadjLO, ZZWMjLO, ZZhadWMjLO, ZZWMhadjLO)     ! pp -> ZZW j
         n_p = 3
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
********************
* WWW processes    *
********************
      CASE(WWWP, WPhadWMWP, WPWMhadWP, WWWM, WMhadWPWM, WMWPhadWM)          ! pp -> WWW
         n_p = 2
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-4
         N_points(3) = N_points(0)-4
      CASE(WWWPjLO, WPhadWMWPjLO, WPWMhadWPjLO, WWWMjLO, WMhadWPWMjLO, WMWPhadWMjLO)    ! pp -> WWWj
         n_p = 3
         n_v = 6
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W+W-A processes *
*********************
      CASE(WWA, WPhadWMA, WPWMhadA)       ! pp -> W+W-A
         n_p = 2
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WWAjLO, WPhadWMAjLO, WPWMhadAjLO)     ! pp -> W+W-A j
         n_p = 3
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* ZZA processes *
*********************
      CASE(ZZA, ZZhadA, ZZnA)       ! pp -> ZZA
         n_p = 2
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(ZZAjLO, ZZhadAjLO, ZZnAjLO)     ! pp -> ZZA j
         n_p = 3
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W+ZA processes *
*********************
      CASE(WPZA, WPhadZA, WPZhadA)       ! pp -> W+ZA
         n_p = 2
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WPZAjLO, WPhadZAjLO, WPZhadAjLO)     ! pp -> W+ZA j
         n_p = 3
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W-ZA processes *
*********************
      CASE(WMZA, WMhadZA, WMZhadA)       ! pp -> W-ZA
         n_p = 2
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WMZAjLO, WMhadZAjLO, WMZhadAjLO)     ! pp -> W-ZA j
         n_p = 3
         n_v = 5
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W+AA processes *
*********************
      CASE(WPAA)       ! pp -> W+AA
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WPAAjLO)     ! pp -> W+AA j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W-AA processes *
*********************
      CASE(WMAA)       ! pp -> W-AA
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(WMAAjLO)     ! pp -> W-AA j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* ZAA processes *
*********************
      CASE(ZAA,ZnAA)       ! pp -> ZAA
         n_p = 2
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(ZAAjLO,ZnAAjLO)     ! pp -> ZAA j
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* AAA processes *
*********************
      CASE(AAA)       ! pp -> AAA
         n_p = 2
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.true.
         NLO_Loops = 3
         N_points(1) = N_points(0)
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
      CASE(AAAjLO)     ! pp -> AAA j
         n_p = 3
         n_v = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
         hasNLO=.false.
         NLO_Loops = 0
#endif

#ifdef WITH_TRIBOSONJET
*********************
* W+AAj processes *
*********************
      CASE(WPAAj)       ! pp -> W+AAj
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 5
         N_points(1) = N_points(0) +1
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
         N_points(4) = N_points(0)-7
         N_points(5) = N_points(0)-5
      CASE(WPAAjjLO)     ! pp -> W+AA jj
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0
*********************
* W-AAj processes *
*********************
      CASE(WMAAj)       ! pp -> W-AAj
         n_p = 3
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.true.
         NLO_Loops = 5
         N_points(1) = N_points(0) +1
         N_points(2) = N_points(0)-2
         N_points(3) = N_points(0)-4
         N_points(4) = N_points(0)-7
         N_points(5) = N_points(0)-5
      CASE(WMAAjjLO)     ! pp -> W-AA jj
         n_p = 4
         n_v = 4
         PS_dimension = 3*(n_p+n_v-nx)-3
         hasNLO=.false.
         NLO_Loops = 0


#endif

      CASE DEFAULT
         print*,"In initprocess no initialization of phase space dimension."
         print*,"Maybe the desired process was not included in configure-option enable_processes!"
         stop 2
      END SELECT

c set n_kin, the number of different kinematics to n_kin = 1
c for LO and virtual parts of the calculation
      n_kin = 1
      n_qed = 1


c now set particle ids
      SELECT CASE(procID)
************************************
* VBF/QCD Higgs, WW & ZZ processes *
************************************
      CASE(Hjj) ! VBF pp -> Hjj 
         particle_IDs(1)= 100  ! no real particle
         particle_IDs(2)= -100 ! no real particle
      CASE(Hjj_AA, AAjj,QCDAAjj) ! VBF pp -> Hjj (for AAjj: spin-2 jj)-> gammagammajj
         particle_IDs(1)= 22 ! photon
         particle_IDs(2)= 22 ! photon
      CASE(Hjj_mu) ! VBF pp -> Hjj -> mu+mu-jj
         particle_IDs(1)= -13  ! mu+
         particle_IDs(2)= 13 ! mu-
      CASE(Hjj_tau) ! VBF pp -> Hjj -> tau+tau-jj 
         particle_IDs(1)= -15 ! tau+
         particle_IDs(2)= 15  ! tau-
      CASE(Hjj_bbar) ! VBF pp -> Hjj -> bbjj
         particle_IDs(1)= 5  ! b
         particle_IDs(2)= -5 ! bbar
      CASE(Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu, WPWMjj, WPWPjj, WMWMjj, 
     1       ZZjj_ll, ZZjj_lnu, Zjj_l, Zjj_nu, WPjj, WMjj, WPZjj, WMZjj, 
     2       Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo 
      CASE(WPAjj, WMAjj,QCDWPAjj,QCDWMAjj,QCDZAjj_l,QCDZAjj_n, ZAjj, ZAjj_n)
	 do i = 1, N_lepton
	    particle_IDs(i) = finallep(i)
         enddo 
	 particle_IDs(3) = 22
      ! vbf processes where 1st boson decays hadronically
      CASE(WPhadWMjj, ZZhadjj, WPhadWPjj, WMhadWMjj,
     &     WPhadZjj, WMhadZjj, Hjj_WPhadWM, Hjj_ZZhad)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! vbf processes where 2nd boson decays hadronically
      CASE(WPWMhadjj, Hjj_WPWMhad,
     &     WPZhadjj, WMZhadjj)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      CASE(HHjj)
         particle_IDs(1) = 25 ! higgs
         particle_IDs(2) = 25 ! higgs
      CASE(HHjj_bbtau) ! VBF pp -> HHjj -> bbtau+tau-jj
         particle_IDs(1) = 5   ! b
         particle_IDs(2) = -5  ! bbar
         particle_IDs(3) = -15 ! tau+
         particle_IDs(4) = 15  ! tau-
      CASE(HHjj_bbAA) ! VBF pp -> HHjj -> bbAAjj
         particle_IDs(1) = 5   ! b
         particle_IDs(2) = -5  ! bbar
         particle_IDs(3) = 22  ! photon
         particle_IDs(4) = 22  ! photon
      CASE(Ajj)
         particle_IDs(1)= 22      
      CASE(HAjj)                ! VBF pp -> HAjj 
         particle_IDs(1)= 100  ! no real particle
         particle_IDs(2)= -100 ! no real particle
         particle_IDs(3)= 22
      CASE(HAjj_AA) ! VBF pp -> HAjj -> A A A jj
         particle_IDs(1)= 22 ! photon
         particle_IDs(2)= 22 ! photon
         particle_IDs(3)= 22
      CASE(HAjj_mu) ! VBF pp -> HAjj -> A mu+mu-jj
         particle_IDs(1)= -13  ! mu+
         particle_IDs(2)= 13 ! mu-
         particle_IDs(3)= 22
      CASE(HAjj_tau) ! VBF pp -> HAjj -> A tau+tau-jj 
         particle_IDs(1)= -15 ! tau+
         particle_IDs(2)= 15  ! tau-
         particle_IDs(3)= 22
      CASE(HAjj_bbar) ! VBF pp -> HAjj -> A bbjj
         particle_IDs(1)= 5  ! b
         particle_IDs(2)= -5 ! bbar
         particle_IDs(3)= 22
      CASE(HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
         particle_IDs(5)= 22
**********************
* VBF Higgs + 3 jets *
**********************
      CASE(Hjjj) ! VBF pp -> Hjjj
         particle_IDs(1)= 100  ! no real particle
         particle_IDs(2)= -100 ! no real particle
      CASE(Hjjj_AA)              ! VBF pp -> Hjjj -> gammagammajj
         particle_IDs(1)= 22 ! photon
         particle_IDs(2)= 22 ! photon
      CASE(Hjjj_mu) ! VBF pp -> Hjjj -> mu+mu-jjj
         particle_IDs(1)= -13  ! mu+
         particle_IDs(2)= 13 ! mu-
      CASE(Hjjj_tau) ! VBF pp -> Hjjj -> tau+tau-jjj 
         particle_IDs(1)= -15 ! tau+
         particle_IDs(2)= 15  ! tau-
      CASE(Hjjj_bbar) ! VBF pp -> Hjjj -> bbjjj
         particle_IDs(1)= 5  ! b
         particle_IDs(2)= -5 ! bbar
      CASE(Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
******************************************************
* GF Higgs processes and other QCD induced processes *
*****************************************************
      CASE(GFHjj)
         particle_IDs(1) = 25
      CASE(GFHjj_AA)                ! GF pp -> Hjj -> gammagammajj
         particle_IDs(1)= 22 ! photon
         particle_IDs(2)= 22 ! photon
      CASE(GFHjj_mu) ! GF pp -> Hjj -> mu+mu-jj
         particle_IDs(1)= -13  ! mu+
         particle_IDs(2)= 13 ! mu-
      CASE(GFHjj_tau) ! GF pp -> Hjj -> tau+tau-jj 
         particle_IDs(1)= -15 ! tau+
         particle_IDs(2)= 15  ! tau-
      CASE(GFHjj_bbar) ! GF pp -> Hjj -> bbjj
         particle_IDs(1)= 5  ! b
         particle_IDs(2)= -5 ! bbar
      CASE(GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu, 
     -        QCDZjj_l, QCDZjj_nu,
     -        QCDWPjj,QCDWMjj,QCDWPZjj,QCDWMZjj,
     -        QCDWPWPjj,QCDWMWMjj,QCDZZjj_ll,QCDZZjj_lnu)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
*********************
* VV and VVV processes *
*********************
      CASE(WW, GFWW, WWZ, ZZZ, ZZWP, ZZWM, WWWP, WWWM,
     -     WPZ, WMZ, ZZ, GFZZ)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
      ! diboson processes where 1st boson decays hadronically
      CASE(WPhadWMlep, GFWPhadWMlep,
     &     ZZhad, GFZZhad,
     &     WPhadZ, WMhadZ)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! diboson processes where 2nd boson decays hadronically
      CASE(WPlepWMhad, GFWPlepWMhad,
     &     WPZhad, WMZhad)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      ! triboson processes where 1st boson decays hadronically
      CASE(WPhadWMZ, ZZhadWP, ZZhadWM, WPhadWMWP, WMhadWPWM, ZZZhad)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
         Particle_IDs(5) = finallep(3)
         Particle_IDs(6) = finallep(4)
      ! triboson processes where 2nd boson decays hadronically
      CASE(WPWMhadZ, WPWMhadWP, WMWPhadWM)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
         Particle_IDs(5) = finallep(3)
         Particle_IDs(6) = finallep(4)
      ! triboson processes where 3rd boson decays hadronically
      CASE(WWZhad, ZZWPhad, ZZWMhad)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finallep(3)
         Particle_IDs(4) = finallep(4)
         Particle_IDs(5) = finalquarks(1)
         Particle_IDs(6) = finalquarks(2)
      ! triboson processes with photon where 1st boson decays hadronically
      CASE(WPhadWMA, ZZhadA, WPhadZA, WMhadZA)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
         Particle_IDs(5) = 22
      ! triboson processes with photon where 2nd boson decays hadronically
      CASE(WPWMhadA, WPZhadA, WMZhadA)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
         Particle_IDs(5) = 22
      CASE(
     1     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu,
     1     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu,
     1     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     1     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu
     2     )
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
         ! IDs for Higgs and the decay products handled later
      CASE(WWA, ZZA, ZZnA, WPZA, WMZA,WPA,WMA, ZA, GFZA)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
         Particle_IDs(N_lepton+1) = 22
      CASE(WPAA,WMAA,ZAA,ZnAA,WPAAj,WMAAj)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
         Particle_IDs(N_lepton+1) = 22
         Particle_IDs(N_lepton+2) = 22
      CASE(AA,GFAA)
         Particle_IDs(1) = 22
         Particle_IDs(2) = 22
      CASE(AAA)
         Particle_IDs(1) = 22
         Particle_IDs(2) = 22
         Particle_IDs(3) = 22
*********************
* VVj processes     *
*********************
      CASE(WWJ, GFWWj, WMZJ, WPZJ, ZZJ, GFZZj)
         do i=1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
      CASE(WMAJ, WPAJ) 
         do i=1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo
         Particle_IDs(3) = 22
      ! diboson + jet processes where 1st boson decays hadronically
      CASE(WPHADWMJ,GFWPHADWMJ,WMhadZJ, WPhadZJ)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! diboson + jet processes where 2nd boson decays hadronically
      CASE(WPWMHADJ,GFWPWMHADJ,WMZhadJ, WPZhadJ)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)

*********************
* V + nj processes  *
*********************
      CASE(WP_only, WM_only, WPJ, WMJ)
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
         enddo

************************
*    Default target    *
************************
      CASE DEFAULT
         print*,"No particle ID assignment for this process "
         stop
      END SELECT
      

****************************************************************
* Set particle IDs for Higgs decay products in WH(j) processes *
****************************************************************

      SELECT CASE(procID)
      CASE(WPH, WPHJ, WMH, WMHJ) ! pp -> H 
         particle_IDs(N_lepton+1)= 100  ! no real particle
         particle_IDs(N_lepton+2)= -100 ! no real particle
      CASE(WPH_AA, WPHJ_AA, WMH_AA, WMHJ_AA) ! pp -> H -> gammagamma
         particle_IDs(N_lepton+1)= 22 ! photon
         particle_IDs(N_lepton+2)= 22 ! photon
      CASE(WPH_mu, WPHJ_mu, WMH_mu, WMHJ_mu) ! pp -> H -> mu+mu-
         particle_IDs(N_lepton+1)= -13  ! mu+
         particle_IDs(N_lepton+2)= 13 ! mu-
      CASE(WPH_tau, WPHJ_tau, WMH_tau, WMHJ_tau) ! pp -> H -> tau+tau- 
         particle_IDs(N_lepton+1)= -15 ! tau+
         particle_IDs(N_lepton+2)= 15  ! tau-
      CASE(WPH_bbar, WPHJ_bbar, WMH_bbar, WMHJ_bbar) ! pp -> H -> bb
         particle_IDs(N_lepton+1)= 5  ! b
         particle_IDs(N_lepton+2)= -5 ! bbar
      END SELECT
      
      RETURN
      END


c*************************************************************************
      SUBROUTINE PrintProcInfo(printrefonly)
c*************************************************************************
c     Prints out the information in "procinfo.dat" depending on which
c     process is selected.
c*************************************************************************
      implicit none
#include "global.inc"
#include "process.inc"

      logical printrefonly

      character*250 filename
      character*8 istr 
      character*250 line
      integer iunit
      logical lglobalprint, lreference

      character*250 GetInputPath
      external GetInputPath

      iunit=10
      lglobalprint = .false.

c Open file
      filename = "procinfo.dat"
      call GetInputFile(iunit,filename,*77)

      if (LOplusjet) then
        write (istr,"(a3,I4.4)") "###",procID
      else
        write (istr,"(a2,I4.4)") "##",procID
      endif

      
      if (.not. printrefonly) then
        print *,""
        print *,"      information on the selected process "
        print *,"-----------------------------------------------"
      endif
      
c read lines and print infos
 50   CONTINUE
      read(iunit,"(a100)",end=99,err=99) line    
      if (index(line,"##").eq.1) lglobalprint = .false.
      if (INDEX(line,istr).eq.1) then 
         lglobalprint = .true.
         lreference = .false.
         goto 50
      endif
      if (lglobalprint .and. .not.printrefonly) print *, trim(line)
c get "References:" lines separately if printrefonly
      if (INDEX(line,"Reference").eq.1) then
        lreference = .true.
        if (printrefonly) line = line((index(line,":")+2):len_trim(line))
        if (lglobalprint .and. printrefonly) then
          write(*,*)'-----------------------------------------------------------------'
          write(*,*)'Process-specific references (please cite as well if appropriate):'
        endif
      endif
      if (lglobalprint .and. printrefonly .and. lreference) then
        print *, trim(line)
      endif

      GOTO 50
 99   CONTINUE
      
      close(iunit)
      RETURN
 77   CONTINUE
      print *
      print *," Error in reading ",trim(filename)
      print *
      RETURN
      END


#ifdef WITH_NLO
c*************************************************************************
      SUBROUTINE InitRealEmissions
c*************************************************************************
c     This reinitializes certain parameters that change for the calculation
c     of the real emissions (e.g. the phase-space dimension)
c*************************************************************************
      implicit none

#include "global.inc"
#include "process.inc"

   
      SELECT CASE(process)
#ifdef WITH_VBF
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW,
     1     Hjj_ZZ_ll, Hjj_ZZ_lnu,
     2     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad)
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         call qqhqqgi(0)     
      CASE(AAjj, Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)   
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(Ajj)
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3 
      CASE(HAjj, HAjj_AA, HAjj_mu, 
     2     HAjj_tau, HAjj_bbar, HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu)
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         call Get_HAjj_couplings()
      CASE(WPWMjj,WPhadWMjj,WPWMhadjj,WPWPjj,WMWMjj,
     &     WPhadWPjj, WMhadWMjj)  ! note: W-W-jj is obtained via W+W+jj.  
         n_p = 5                  !    p and pbar have already been switched
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         call qqhqqgi(3)
      CASE(ZZjj_ll, ZZjj_lnu, ZZhadjj)   
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(WPZjj,WPhadZjj,WPZhadjj,WMZjj,WMhadZjj,WMZhadjj)   
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(HHjj, HHjj_bbtau, HHjj_bbAA)
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(WPAjj,WMAjj,ZAjj,ZAjj_n)   
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
#endif

#ifdef WITH_HJJJ
      CASE(Hjjj, Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_bbar, Hjjj_WW,
     1     Hjjj_ZZ_ll, Hjjj_ZZ_lnu)
         n_p = 6
         n_kin = 25
         PS_dimension = 3*(n_p+n_v-nx)-3
#endif

      CASE(Zjj_l, Zjj_nu, WPjj, WMjj)
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3   
      CASE(WW, WWZ, ZZWP, ZZWM, WWWP, WWWM, WWA, ZZA, ZZnA, WPZA, WMZA, 
     1     WPAA, WMAA, ZAA, ZnAA, WPZ, WMZ, ZZ,
     2     WPhadWMlep,WPlepWMhad,ZZhad,
     3     WPhadZ,WPZhad,WMhadZ,WMZhad,
     4     WPhadWMZ,WPWMhadZ,WWZhad,ZZhadWP,ZZWPhad,ZZhadWM,ZZWMhad,
     5     WPhadWMWP,WPWMhadWP,WMhadWPWM,WMWPhadWM,
     6     WPhadWMA, WPWMhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA, ZZhadA)
         n_p = 3
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(WPA, WMA, ZA, AA, 
     &     WP_only, WM_only,
     1     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu,
     1     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu
     2    )
         n_p = 3
         n_kin = 3 ! needed randomness for real emission
         PS_dimension = 3*(n_p+n_v-nx)-2
      CASE(WWJ,WMAJ,WPAJ,WMZJ,WPZJ,ZZJ,
     &     WPHADWMJ,WPWMHADJ,
     &     WMhadZJ,WMZhadJ,WPhadZJ,WPZhadJ,
     &     WPJ, WMJ,
     1     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     1     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu
     &   )
         n_p = 4
         n_kin = 11
         PS_dimension = 3*(n_p+n_v-nx)-3
      CASE(ZZZ, AAA, ZZZhad)
         n_p = 3
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-2
      CASE(WPAAj, WMAAj)
         n_p = 4
         n_kin = 11
         PS_dimension = 3*(n_p+n_v-nx)-3
#ifdef WITH_QCDVV
      CASE(QCDWPZjj,QCDWMZjj,QCDWPWPjj,QCDWMWMjj,QCDZZjj_ll,QCDZZjj_lnu
     &,QCDWPAjj,QCDWMAjj,QCDZAjj_l,QCDZAjj_n,QCDAAjj)
         n_p = 5
         n_kin =16
         PS_dimension = 3*(n_p+n_v-nx)-3+1
#endif

#ifdef WITH_QCDV
      CASE(QCDWPjj,QCDWMjj,QCDZjj_l,QCDZjj_nu)
         n_p = 5
         n_kin =16
         PS_dimension = 3*(n_p+n_v-nx)-3+1
#endif
      CASE DEFAULT
         print*,"No definition of the PS dimension "//
     &          "for the real emission"
         stop
      END SELECT

      RETURN
      END


c*************************************************************************
      SUBROUTINE InitQEDEmissions
c*************************************************************************
c     This reinitializes certain parameters that change for the calculation
c     of the qed emissions (e.g. the phase-space dimension)
c*************************************************************************
      implicit none

#include "global.inc"
#include "process.inc"

   
      SELECT CASE(process)
#ifdef WITH_VBF
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW,
     1     Hjj_ZZ_ll, Hjj_ZZ_lnu,
     2     Hjj_ZZhad, Hjj_WPhadWM, Hjj_WPWMhad)
         n_qed = 13
         n_p = 5
         n_kin = 3
         PS_dimension = 3*(n_p+n_v-nx)-3
         call qqhqqgi(0)     
#endif
      CASE DEFAULT
         print*,"No definition of the PS dimension "//
     &          "for the real emission"
         stop
      END SELECT

      RETURN
      END


c*************************************************************************
      SUBROUTINE InitFLoop
c*************************************************************************
c     This reinitializes certain parameters that change for the calculation
c     of the fermion loop 
c*************************************************************************
      implicit none

#include "global.inc"
#include "process.inc"

      SELECT CASE(process)
      CASE(WW,ZZ,ZA,AA,WPhadWMlep,WPlepWMhad,ZZhad)
c nothing at the moment
      CASE DEFAULT
c nothing at the moment
      END SELECT

      RETURN
      END
c ifdef with_nlo
#endif


c*************************************************************************  
      subroutine proc_assignment
c*************************************************************************  
c     INPUT: procID       number to identify process
c
c     OUTPUT:   common blocks in process.inc are filled with
c               N_lepton: number of expected leptons in the final state
c               process:  process identifier
c               finallep_default:  default numbers for final state leptons
c     
c*************************************************************************  
          use globalvars, only: ldoblha
          use readinput, only: read_logical
      implicit none

      logical kkswitch, anomswitch
     
#include "process.inc"
#include "global.inc"
#include "cuts.inc"


      if (ldoblha) then
        kkswitch = .false.
        anomswitch = .false.
      else
        call read_logical("KK_MOD",kkswitch,.false.)
        call read_logical("ANOM_CPL",anomswitch,.false.) 
      endif

      Loops_sub_LO=1
      Loops_sub_NLO=1
      PS_Loops = 1

      ! for hadronic decay: default is all combinations of four light quarks
      N_decayquarks = 0
      finalquarks_default(1)=  93
      finalquarks_default(2)=  93

c ********************************

      SELECT CASE(procID)
      CASE (QCDAAjj) 
         print*," "
         print*,"This process is not yet implemented (procID = ",procID,") !"
         print*," "
         stop
      END SELECT

c ********************************

      SELECT CASE(procID)
      CASE (WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ)
         print*," "
         print*,"Hadronic decay implementation not available for this process ID (procID = ",procID,")"
         print*," "
         stop
      END SELECT

c ********************************
c set common options for classes of processes 
c group different decay channels
c      
      SELECT CASE(procID)
      CASE(
     1     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu,
     1     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     2     WP_only
     &    )
         N_lepton = 2                     ! for WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu this is overwritten later
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(
     1     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu,
     1     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu,
     2     WM_only      
     &    )
         N_lepton = 2                     ! for WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu this is overwritten later
         finallep_default(1)=  11
         finallep_default(2)= -12
      END SELECT


c ********************************
c set common options for classes of processes 
c group different number of subprocesses
c      
      SELECT CASE(procID)
      CASE(
     1     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     1     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu
     &    )
         if (LOplusjet) then
            Loops_sub_LO=1 ! 2 to seperate 2 quarkline and 1 quarkline diagrams
         endif
         Loops_sub_NLO=1
      END SELECT

c   set variable process, number of leptons in the final state and
c   default values for the final state leptons
      SELECT CASE(procID)
      CASE(Hjj)
         if (LOplusjet) process = HjjjLO
         if (.not.LOplusjet) process = Hjj
         N_lepton = 0
      CASE(Hjj_AA)
         if (LOplusjet) process = HjjjLO_AA
         if (.not.LOplusjet) process = Hjj_AA
         N_lepton = 0
      CASE(AAjj)
         if (LOplusjet) process = AAjjjLO
         if (.not.LOplusjet) process = AAjj
         N_lepton = 0
      CASE(Sp2jj_WW)
         if (LOplusjet) process = Sp2jjjLO_WW
         if (.not.LOplusjet) process = Sp2jj_WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)= -12
         finallep_default(4)=  11
      CASE(Sp2jj_ZZ_ll)
         if (LOplusjet) process = Sp2jjjLO_ZZ_ll
         if (.not.LOplusjet) process = Sp2jj_ZZ_ll
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)= -11
         finallep_default(4)=  11
      CASE(Sp2jj_ZZ_lnu)
         if (LOplusjet) process = Sp2jjjLO_ZZ_lnu
         if (.not.LOplusjet) process = Sp2jj_ZZ_lnu
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(Hjj_mu)
         if (LOplusjet) process = HjjjLO_mu
         if (.not.LOplusjet) process = Hjj_mu
         N_lepton = 0
      CASE(Hjj_tau)
         if (LOplusjet) process = HjjjLO_tau
         if (.not.LOplusjet) process = Hjj_tau
         N_lepton = 0
      CASE(Hjj_bbar)
         if (LOplusjet) process = HjjjLO_bbar
         if (.not.LOplusjet) process = Hjj_bbar
         N_lepton = 0
      CASE(Hjj_WW)
         if (LOplusjet) process = HjjjLO_WW
         if (.not.LOplusjet) process = Hjj_WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(Hjj_WPhadWM)
         if (LOplusjet) process = HjjjLO_WPhadWM
         if (.not.LOplusjet) process = Hjj_WPhadWM
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(Hjj_WPWMhad)
         if (LOplusjet) process = HjjjLO_WPWMhad
         if (.not.LOplusjet) process = Hjj_WPWMhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(Hjj_ZZ_ll)
         if (LOplusjet) process = HjjjLO_ZZ_ll
         if (.not.LOplusjet) process = Hjj_ZZ_ll
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(Hjj_ZZhad)
         if (LOplusjet) process = HjjjLO_ZZhad
         if (.not.LOplusjet) process = Hjj_ZZhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
      CASE(Hjj_ZZ_lnu)
         if (LOplusjet) process = HjjjLO_ZZ_lnu
         if (.not.LOplusjet) process = Hjj_ZZ_lnu
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(HHjj)
         if (LOplusjet) process = HHjjjLO
         if (.not.LOplusjet) process = HHjj
         N_lepton = 0
      CASE(HHjj_bbtau)
         if (LOplusjet) process = HHjjjLO_bbtau
         if (.not.LOplusjet) process = HHjj_bbtau
         N_lepton = 0
      CASE(HHjj_bbAA)
         if (LOplusjet) process = HHjjjLO_bbAA
         if (.not.LOplusjet) process = HHjj_bbAA
         N_lepton = 0
      CASE(HAjj)
         if (LOplusjet) process = HAjjjLO
         if (.not.LOplusjet) process = HAjj
         N_lepton = 0
      CASE(HAjj_AA)
         if (LOplusjet) process = HAjjjLO_AA
         if (.not.LOplusjet) process = HAjj_AA
         N_lepton = 0
      CASE(HAjj_mu)
         if (LOplusjet) process = HAjjjLO_mu
         if (.not.LOplusjet) process = HAjj_mu
         N_lepton = 0
      CASE(HAjj_tau)
         if (LOplusjet) process = HAjjjLO_tau
         if (.not.LOplusjet) process = HAjj_tau
         N_lepton = 0
      CASE(HAjj_bbar)
         if (LOplusjet) process = HAjjjLO_bbar
         if (.not.LOplusjet) process = HAjj_bbar
         N_lepton = 0
      CASE(HAjj_WW)
         if (LOplusjet) process = HAjjjLO_WW
         if (.not.LOplusjet) process = HAjj_WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(HAjj_ZZ_ll)
         if (LOplusjet) process = HAjjjLO_ZZ_ll
         if (.not.LOplusjet) process = HAjj_ZZ_ll
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(HAjj_ZZ_lnu)
         if (LOplusjet) process = HAjjjLO_ZZ_lnu
         if (.not.LOplusjet) process = HAjj_ZZ_lnu
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(Ajj)
         if (LOplusjet) process = AjjjLO
         if (.not.LOplusjet) process = Ajj
         N_lepton = 0
      CASE(Zjj_l)
         if (LOplusjet) process = ZjjjLO_l
         if (.not.LOplusjet) process = Zjj_l
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)=  11
      CASE(Zjj_nu)
         if (LOplusjet) process = ZjjjLO_nu
         if (.not.LOplusjet) process = Zjj_nu
         N_lepton = 2
         finallep_default(1)= -12
         finallep_default(2)=  12
      CASE(WPjj)
         if (LOplusjet) process = WPjjjLO
         if (.not.LOplusjet) process = WPjj
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)=  12
      CASE(WMjj)
         if (LOplusjet) process = WMjjjLO
         if (.not.LOplusjet) process = WMjj
         N_lepton = 2
         finallep_default(1)= -12
         finallep_default(2)=  11
      CASE(WPWMjj)
         if (LOplusjet) process = WPWMjjjLO
         if (.not.LOplusjet) process = WPWMjj
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)= -12
         finallep_default(4)=  11
      CASE(WPhadWMjj)
         if (LOplusjet) process = WPhadWMjjjLO
         if (.not.LOplusjet) process = WPhadWMjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -12
         finallep_default(2)=  11
      CASE(WPWMhadjj)
         if (LOplusjet) process = WPWMhadjjjLO
         if (.not.LOplusjet) process = WPWMhadjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(WPWPjj)
         if (LOplusjet) process = WPWPjjjLO
         if (.not.LOplusjet) process = WPWPjj
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)= 12
         finallep_default(3)= -11
         finallep_default(4)= 12
      CASE(WPhadWPjj)
         if (LOplusjet) process = WPhadWPjjjLO
         if (.not.LOplusjet) process = WPhadWPjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -11
         finallep_default(2)= 12
      CASE(WMWMjj)
         if (LOplusjet) process = WMWMjjjLO
         if (.not.LOplusjet) process = WMWMjj
         N_lepton = 4
         finallep_default(1)= 11
         finallep_default(2)= -12
         finallep_default(3)= 11
         finallep_default(4)= -12
      CASE(WMhadWMjj)
         if (LOplusjet) process = WMhadWMjjjLO
         if (.not.LOplusjet) process = WMhadWMjj
         N_decayquarks = 2
         N_lepton = 2
         finallep_default(1)= 11
         finallep_default(2)= -12
      CASE(WPZjj)
         if (LOplusjet) process = WPZjjjLO
         if (.not.LOplusjet) process = WPZjj
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  12
         finallep_default(3)= -11
         finallep_default(4)=  11
      CASE(WPhadZjj)
         if (LOplusjet) process = WPhadZjjjLO
         if (.not.LOplusjet) process = WPhadZjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -11
         finallep_default(2)=  11
      CASE(WPZhadjj)
         if (LOplusjet) process = WPZhadjjjLO
         if (.not.LOplusjet) process = WPZhadjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -11
         finallep_default(2)=  12
      CASE(WMZjj)
         if (LOplusjet) process = WMZjjjLO
         if (.not.LOplusjet) process = WMZjj
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -12
         finallep_default(3)= -11
         finallep_default(4)=  11
      CASE(WMhadZjj)
         if (LOplusjet) process = WMhadZjjjLO
         if (.not.LOplusjet) process = WMhadZjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -11
         finallep_default(2)=  11
      CASE(WMZhadjj)
         if (LOplusjet) process = WMZhadjjjLO
         if (.not.LOplusjet) process = WMZhadjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(WPAjj)
         if (LOplusjet) process = WPAjjjLO
         if (.not.LOplusjet) process = WPAjj
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)=  12
         PS_Loops = 2
      CASE(WMAjj)
         if (LOplusjet) process = WMAjjjLO
         if (.not.LOplusjet) process = WMAjj
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
         PS_Loops = 2
      CASE(ZAjj)
         if (LOplusjet) process = ZAjjjLO
         if (.not.LOplusjet) process = ZAjj
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)=  11
         PS_Loops = 2
      CASE(ZAjj_n)
         if (LOplusjet) process = ZAjjjLO_n
         if (.not.LOplusjet) process = ZAjj_n
         N_lepton = 2
         finallep_default(1)=  12
         finallep_default(2)= -12
         PS_Loops = 1
      CASE(ZZjj_ll)
         if (LOplusjet) process = ZZjjjLO_ll
         if (.not.LOplusjet) process = ZZjj_ll
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)= -11
         finallep_default(4)=  11
      CASE(ZZhadjj)
         if (LOplusjet) process = ZZhadjjjLO
         if (.not.LOplusjet) process = ZZhadjj
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)= -11
         finallep_default(2)=  11
      CASE(ZZjj_lnu)
         if (LOplusjet) process = ZZjjjLO_lnu
         if (.not.LOplusjet) process = ZZjj_lnu
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(Hjjj)
         if (LOplusjet) process = HjjjjLO
         if (.not.LOplusjet) process = Hjjj
         N_lepton = 0
      CASE(Hjjj_AA)
         if (LOplusjet) process = HjjjjLO_AA
         if (.not.LOplusjet) process = Hjjj_AA
         N_lepton = 0
      CASE(Hjjj_mu)
         if (LOplusjet) process = HjjjjLO_mu
         if (.not.LOplusjet) process = Hjjj_mu
         N_lepton = 0
      CASE(Hjjj_tau)
         if (LOplusjet) process = HjjjjLO_tau
         if (.not.LOplusjet) process = Hjjj_tau
         N_lepton = 0
      CASE(Hjjj_bbar)
         if (LOplusjet) process = HjjjjLO_bbar
         if (.not.LOplusjet) process = Hjjj_bbar
         N_lepton = 0
      CASE(Hjjj_WW)
         if (LOplusjet) process = HjjjjLO_WW
         if (.not.LOplusjet) process = Hjjj_WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(Hjjj_ZZ_ll)
         if (LOplusjet) process = HjjjjLO_ZZ_ll
         if (.not.LOplusjet) process = Hjjj_ZZ_ll
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(Hjjj_ZZ_lnu)
         if (LOplusjet) process = HjjjjLO_ZZ_lnu
         if (.not.LOplusjet) process = Hjjj_ZZ_lnu
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(GFHjj)
         process = GFHjj
         N_lepton = 0
      CASE(GFHjj_AA)
         process = GFHjj_AA
         N_lepton = 0
      CASE(GFHjj_mu)
         process = GFHjj_mu
         N_lepton = 0
      CASE(GFHjj_tau)
         process = GFHjj_tau
         N_lepton = 0
      CASE(GFHjj_bbar)
         process = GFHjj_bbar
         N_lepton = 0
      CASE(GFHjj_WW)
         process = GFHjj_WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(GFHjj_ZZ_ll)
         process = GFHjj_ZZ_ll
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(GFHjj_ZZ_lnu)
         process = GFHjj_ZZ_lnu
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  12
         finallep_default(4)= -12
      CASE(QCDZjj_l)
         if (.not.LOplusjet) then
            process = QCDZjj_l
            Loops_sub_LO = 7
         else
            process = QCDZjjjLO_l
         endif
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)= 11
      CASE(QCDZjj_nu)
         if (.not.LOplusjet) then
            process = QCDZjj_nu
            Loops_sub_LO = 7
         else
            process = QCDZjjjLO_nu
         endif
         N_lepton = 2
         finallep_default(2)= -12
         finallep_default(1)= 12
      CASE(QCDWPjj)
         if (.not.LOplusjet) then
	    process = QCDWPjj
	    Loops_sub_LO = 8
	    Loops_sub_NLO= 1
	 else
	    process = QCDWPjjjLO
	 endif
         N_lepton = 2
         finallep_default(2)= 12
         finallep_default(1)= -11
         PS_Loops = 1
      CASE(QCDWMjj)
         if (.not.LOplusjet) then
	    process = QCDWMjj
	    Loops_sub_LO = 8
	    Loops_sub_NLO= 1
	 else
	    process = QCDWMjjjLO
	 endif
         N_lepton = 2
         finallep_default(2)= -12
         finallep_default(1)= 11
         PS_Loops = 1
      CASE(QCDWPZjj)
         if (.not.LOplusjet) then
	    process = QCDWPZjj
	    Loops_sub_LO = 8
	    Loops_sub_NLO= 1
	 else
	    process = QCDWPZjjjLO
	 endif
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)= 12
         finallep_default(3)= -11
         finallep_default(4)= 11
         PS_Loops = 1
      CASE(QCDWMZjj)
         if (.not.LOplusjet) then
	    process = QCDWMZjj
	    Loops_sub_LO = 8
	 else
	    process = QCDWMZjjjLO
	 endif
         N_lepton = 4
         finallep_default(1)= 11
         finallep_default(2)= -12
         finallep_default(3)= 11
         finallep_default(4)= -11
         PS_Loops = 1
      CASE(QCDWPAjj)
         if (.not.LOplusjet) then
            process = QCDWPAjj
            Loops_sub_LO = 8
            Loops_sub_NLO= 1
         else
            process = QCDWPAjjjLO
         endif
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)= 12
         PS_Loops = 2
      CASE(QCDWMAjj)
         if (.not.LOplusjet) then
            process = QCDWMAjj
            Loops_sub_LO = 8
         else
            process = QCDWMAjjjLO
         endif
         N_lepton = 2
         finallep_default(1)= 11
         finallep_default(2)= -12
         PS_Loops = 2
      CASE(QCDZAjj_l)
         if (.not.LOplusjet) then
            process = QCDZAjj_l
            Loops_sub_LO = 7
            Loops_sub_NLO= 1
         else
            process = QCDZAjjjLO_l
         endif
         N_lepton = 2
         finallep_default(1)= -11
         finallep_default(2)= 11
         PS_Loops = 2
      CASE(QCDZAjj_n)
         if (.not.LOplusjet) then
            process = QCDZAjj_n
            Loops_sub_LO = 7
            Loops_sub_NLO= 1
         else
            process = QCDZAjjjLO_n
         endif
         N_lepton = 2
         finallep_default(1)= 12
         finallep_default(2)= -12
         PS_Loops = 1
      CASE(QCDAAjj)
         if (.not.LOplusjet) then
            process = QCDAAjj
            Loops_sub_LO = 7
         else
            process = QCDAAjjjLO
         endif
         N_lepton = 0
         PS_Loops = 1
      CASE(QCDZZjj_ll)
         if (.not.LOplusjet) then
            process = QCDZZjj_ll
            Loops_sub_LO = 7
         else
            process = QCDZZjjjLO_ll
         endif
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)= -11
         finallep_default(4)=  11
         PS_Loops = 1
      CASE(QCDZZjj_lnu)
         if (.not.LOplusjet) then
            process = QCDZZjj_lnu
            Loops_sub_LO = 7
         else
            process = QCDZZjjjLO_lnu
         endif
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)=  11
         finallep_default(3)=  12
         finallep_default(4)= -12
         PS_Loops = 1
      CASE(QCDWPWPjj)
         if (.not.LOplusjet) then
	    process = QCDWPWPjj
	    Loops_sub_LO = 3
	    Loops_sub_NLO = 1
	 else
	    process = QCDWPWPjjjLO
	 endif
         N_lepton = 4
         finallep_default(1)= -11
         finallep_default(2)= 12
         finallep_default(3)= -11
         finallep_default(4)= 12
         PS_Loops = 1
      CASE(QCDWMWMjj)
         if (.not.LOplusjet) then
            process = QCDWMWMjj
            Loops_sub_LO = 3
            Loops_sub_NLO = 1
         else
            process = QCDWMWMjjjLO
         endif
         N_lepton = 4
         finallep_default(1)= 11
         finallep_default(2)= -12
         finallep_default(3)= 11
         finallep_default(4)= -12
         PS_Loops = 1
      CASE(WW)
         if (LOplusjet) process = WWjLO
         if (.not.LOplusjet) process = WW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(WPhadWMlep)
         if (LOplusjet) process = WPhadWMlepjLO
         if (.not.LOplusjet) process = WPhadWMlep
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(WPlepWMhad)
         if (LOplusjet) process = WPlepWMhadjLO
         if (.not.LOplusjet) process = WPlepWMhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(GFWW)
         process = GFWW
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(GFWPhadWMlep)
         process = GFWPhadWMlep
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(GFWPlepWMhad)
         process = GFWPlepWMhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(WPZ)
         if (LOplusjet) process = WPZjLOb
         if (.not.LOplusjet) process = WPZ
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
         N_lepton = 4
      CASE(WPhadZ)
         if (LOplusjet) process = WPhadZjLO
         if (.not.LOplusjet) process = WPhadZ
         finallep_default(1)=  11
         finallep_default(2)= -11
         N_decayquarks = 2
         N_lepton = 2
      CASE(WPZhad)
         if (LOplusjet) process = WPZhadjLO
         if (.not.LOplusjet) process = WPZhad
         finallep_default(1)=  12
         finallep_default(2)= -11
         N_lepton = 2
         N_decayquarks = 2
      CASE(WMZ)
         if (LOplusjet) process = WMZjLOb
         if (.not.LOplusjet) process = WMZ
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -12
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(WMhadZ)
         if (LOplusjet) process = WMhadZjLO
         if (.not.LOplusjet) process = WMhadZ
         N_decayquarks = 2
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
      CASE(WMZhad)
         if (LOplusjet) process = WMZhadjLO
         if (.not.LOplusjet) process = WMZhad
         N_decayquarks = 2
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(ZZ)
         if (LOplusjet) process = ZZjLO
         if (.not.LOplusjet) process = ZZ
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(ZZhad)
         if (LOplusjet) process = ZZhadjLO
         if (.not.LOplusjet) process = ZZhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
      CASE(GFZZ)
         process = GFZZ
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(GFZZhad)
         process = GFZZhad
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
      CASE(WPA)
         if (LOplusjet) process = WPAjLOb
         if (.not.LOplusjet) process = WPA
         N_lepton = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
         PS_Loops = 2
      CASE(WMA)
         if (LOplusjet) process = WMAjLOb
         if (.not.LOplusjet) process = WMA
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
         PS_Loops = 2

         !Higgsstrahlung
      CASE(WPH)
         if (LOplusjet) process = WPHjLOb
         if (.not.LOplusjet) process = WPH
      CASE(WPH_AA)
         if (LOplusjet) process = WPHjLOb_AA
         if (.not.LOplusjet) process = WPH_AA
      CASE(WPH_mu)
         if (LOplusjet) process = WPHjLOb_mu
         if (.not.LOplusjet) process = WPH_mu
      CASE(WPH_tau)
         if (LOplusjet) process = WPHjLOb_tau
         if (.not.LOplusjet) process = WPH_tau
      CASE(WPH_bbar)
         if (LOplusjet) process = WPHjLOb_bbar
         if (.not.LOplusjet) process = WPH_bbar
      CASE(WPH_WW)
         if (LOplusjet) process = WPHjLOb_WW
         if (.not.LOplusjet) process = WPH_WW
         N_lepton = 6
         finallep_default(3)=  12
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -12
      CASE(WPH_ZZ_ll)
         if (LOplusjet) process = WPHjLOb_ZZ_ll
         if (.not.LOplusjet) process = WPH_ZZ_ll
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -11
      CASE(WPH_ZZ_lnu)
         if (LOplusjet) process = WPHjLOb_ZZ_lnu
         if (.not.LOplusjet) process = WPH_ZZ_lnu
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  12
         finallep_default(6)= -12
      CASE(WMH)
         if (LOplusjet) process = WMHjLOb
         if (.not.LOplusjet) process = WMH
      CASE(WMH_AA)
         if (LOplusjet) process = WMHjLOb_AA
         if (.not.LOplusjet) process = WMH_AA
      CASE(WMH_mu)
         if (LOplusjet) process = WMHjLOb_mu
         if (.not.LOplusjet) process = WMH_mu
      CASE(WMH_tau)
         if (LOplusjet) process = WMHjLOb_tau
         if (.not.LOplusjet) process = WMH_tau
      CASE(WMH_bbar)
         if (LOplusjet) process = WMHjLOb_bbar
         if (.not.LOplusjet) process = WMH_bbar
      CASE(WMH_WW)
         if (LOplusjet) process = WMHjLOb_WW
         if (.not.LOplusjet) process = WMH_WW
         N_lepton = 6
         finallep_default(3)=  12
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -12
      CASE(WMH_ZZ_ll)
         if (LOplusjet) process = WMHjLOb_ZZ_ll
         if (.not.LOplusjet) process = WMH_ZZ_ll
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -11
      CASE(WMH_ZZ_lnu)
         if (LOplusjet) process = WMHjLOb_ZZ_lnu
         if (.not.LOplusjet) process = WMH_ZZ_lnu
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  12
         finallep_default(6)= -12
      CASE(WP_only)
         if (LOplusjet) process = WPjLOb
         if (.not.LOplusjet) process = WP_only
      CASE(WM_only)
         if (LOplusjet) process = WMjLOb
         if (.not.LOplusjet) process = WM_only

      CASE(ZA)
         if (LOplusjet) process = ZAjLO
         if (.not.LOplusjet) process = ZA
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
         PS_Loops = 2
      CASE(GFZA)
         process = GFZA
         N_lepton = 2
         finallep_default(1)=  11
         finallep_default(2)= -11
      CASE(AA)
         if (LOplusjet) process = AAjLO
         if (.not.LOplusjet) process = AA
         N_lepton = 0
      CASE(GFAA)
         process = GFAA
         N_lepton = 0
      CASE(WWJ)
         if (LOplusjet) process = WWJJLO
         if (.not.LOplusjet) process = WWJ
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(WPHADWMJ)
         if (LOplusjet) process = WPHADWMJJLO
         if (.not.LOplusjet) process = WPHADWMJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(WPWMHADJ)
         if (LOplusjet) process = WPWMHADJJLO
         if (.not.LOplusjet) process = WPWMHADJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(GFWWJ)
         process = GFWWJ
         N_lepton = 4
         finallep_default(1)=  12
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -12
      CASE(GFWPHADWMJ)
         process = GFWPHADWMJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  11
         finallep_default(2)= -12
      CASE(GFWPWMHADJ)
         process = GFWPWMHADJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1)=  12
         finallep_default(2)= -11
      CASE(WMAJ)
         if (LOplusjet) process = WMAJLO
         if (.not.LOplusjet) process = WMAJ
         N_lepton = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         PS_Loops = 2
      CASE(WPAJ)
         if (LOplusjet) process = WPAJLO
         if (.not.LOplusjet) process = WPAJ
         N_lepton = 2
         finallep_default(1) = -11
         finallep_default(2) =  12
         PS_Loops = 2
      CASE(WMZJ)
         if (LOplusjet) process = WMZJLO
         if (.not.LOplusjet) process = WMZJ
         N_lepton = 4
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  11
         finallep_default(4) = -11
      CASE(WMhadZJ)
         if (LOplusjet) process = WMhadZJLO
         if (.not.LOplusjet) process = WMhadZJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
      CASE(WMZhadJ)
         if (LOplusjet) process = WMZhadJLO
         if (.not.LOplusjet) process = WMZhadJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
      CASE(WPZJ)
         if (LOplusjet) process = WPZJLO
         if (.not.LOplusjet) process = WPZJ
         N_lepton = 4
         finallep_default(1) = -11
         finallep_default(2) =  12
         finallep_default(3) = -11
         finallep_default(4) =  11
      CASE(ZZJ)
         if (LOplusjet) process = ZZJJLO
         if (.not.LOplusjet) process = ZZJ
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(GFZZJ)
         process = GFZZj
         N_lepton = 4
         finallep_default(1)=  11
         finallep_default(2)= -11
         finallep_default(3)=  11
         finallep_default(4)= -11
      CASE(WPhadZJ)
         if (LOplusjet) process = WPhadZJjLO
         if (.not.LOplusjet) process = WPhadZJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) = -11
         finallep_default(2) =  11
      CASE(WPZhadJ)
         if (LOplusjet) process = WPZhadJLO
         if (.not.LOplusjet) process = WPZhadJ
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) = -11
         finallep_default(2) =  12

      CASE(WMJ)
         if (LOplusjet) then
            process = WMJJLO
            Loops_sub_LO=1! 2 to seperate 2 quarkline and 1 quarkline diagrams
         endif
         if (.not.LOplusjet) process = WMJ
         N_lepton = 2
         finallep_default(1) = 11
         finallep_default(2) = -12
         Loops_sub_NLO=1! 2 to seperate 2 quarkline and 1 quarkline diagrams

      CASE(WPJ)
         if (LOplusjet) then
            process = WPJJLO
            Loops_sub_LO=1! 2 to seperate 2 quarkline and 1 quarkline diagrams
         endif
         if (.not.LOplusjet) process = WPJ
         N_lepton = 2
         finallep_default(1) =  12 
         finallep_default(2) =  -11
         Loops_sub_NLO=1! 2 to seperate 2 quarkline and 1 quarkline diagrams

      CASE(WPHJ)
         if (LOplusjet) process = WPHJJLO
         if (.not.LOplusjet) process = WPHJ
      CASE(WPHJ_AA)
         if (LOplusjet) process = WPHJJLO_AA
         if (.not.LOplusjet) process = WPHJ_AA
      CASE(WPHJ_mu)
         if (LOplusjet) process = WPHJJLO_mu
         if (.not.LOplusjet) process = WPHJ_mu
      CASE(WPHJ_tau)
         if (LOplusjet) process = WPHJJLO_tau
         if (.not.LOplusjet) process = WPHJ_tau
      CASE(WPHJ_bbar)
         if (LOplusjet) process = WPHJJLO_bbar
         if (.not.LOplusjet) process = WPHJ_bbar
      CASE(WPHJ_WW)
         if (LOplusjet) process = WPHJJLO_WW
         if (.not.LOplusjet) process = WPHJ_WW
         N_lepton = 6
         finallep_default(3)=  12
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -12
      CASE(WPHJ_ZZ_ll)
         if (LOplusjet) process = WPHJJLO_ZZ_ll
         if (.not.LOplusjet) process = WPHJ_ZZ_ll
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -11
      CASE(WPHJ_ZZ_lnu)
         if (LOplusjet) process = WPHJJLO_ZZ_lnu
         if (.not.LOplusjet) process = WPHJ_ZZ_lnu
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  12
         finallep_default(6)= -12
      CASE(WMHJ)
         if (LOplusjet) process = WMHJJLO
         if (.not.LOplusjet) process = WMHJ
      CASE(WMHJ_AA)
         if (LOplusjet) process = WMHJJLO_AA
         if (.not.LOplusjet) process = WMHJ_AA
      CASE(WMHJ_mu)
         if (LOplusjet) process = WMHJJLO_mu
         if (.not.LOplusjet) process = WMHJ_mu
      CASE(WMHJ_tau)
         if (LOplusjet) process = WMHJJLO_tau
         if (.not.LOplusjet) process = WMHJ_tau
      CASE(WMHJ_bbar)
         if (LOplusjet) process = WMHJJLO_bbar
         if (.not.LOplusjet) process = WMHJ_bbar
      CASE(WMHJ_WW)
         if (LOplusjet) process = WMHJJLO_WW
         if (.not.LOplusjet) process = WMHJ_WW
         N_lepton = 6
         finallep_default(3)=  12
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -12
      CASE(WMHJ_ZZ_ll)
         if (LOplusjet) process = WMHJJLO_ZZ_ll
         if (.not.LOplusjet) process = WMHJ_ZZ_ll
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  11
         finallep_default(6)= -11
      CASE(WMHJ_ZZ_lnu)
         if (LOplusjet) process = WMHJJLO_ZZ_lnu
         if (.not.LOplusjet) process = WMHJ_ZZ_lnu
         N_lepton = 6
         finallep_default(3)=  11
         finallep_default(4)= -11
         finallep_default(5)=  12
         finallep_default(6)= -12

      CASE(WWZ)
         if (LOplusjet) process = WWZjLO
         if (.not.LOplusjet) process = WWZ
         N_lepton = 6
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
         finallep_default(5) =  11
         finallep_default(6) = -11
      CASE(WPhadWMZ)
         if (LOplusjet) process = WPhadWMZjLO
         if (.not.LOplusjet) process = WPhadWMZ
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  11
         finallep_default(4) = -11
      CASE(WPWMhadZ)
         if (LOplusjet) process = WPWMhadZjLO
         if (.not.LOplusjet) process = WPWMhadZ
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
      CASE(WWZhad)
         if (LOplusjet) process = WWZhadjLO
         if (.not.LOplusjet) process = WWZhad
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
      CASE(ZZZ)
         if (LOplusjet) process = ZZZjLO
         if (.not.LOplusjet) process = ZZZ
         N_lepton = 6
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         finallep_default(5) =  11
         finallep_default(6) = -11
         PS_Loops = 2
      CASE(ZZZhad)
         if (LOplusjet) process = ZZZhadjLO
         if (.not.LOplusjet) process = ZZZhad
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         PS_Loops = 2
      CASE(ZZWP)
         if (LOplusjet) process = ZZWPjLO
         if (.not.LOplusjet) process = ZZWP
         N_lepton = 6
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         finallep_default(5) =  12
         finallep_default(6) = -11
      CASE(ZZhadWP)
         if (LOplusjet) process = ZZhadWPjLO
         if (.not.LOplusjet) process = ZZhadWP
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  12
         finallep_default(4) = -11
      CASE(ZZWPhad)
         if (LOplusjet) process = ZZWPhadjLO
         if (.not.LOplusjet) process = ZZWPhad
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
      CASE(ZZWM)
         if (LOplusjet) process = ZZWMjLO
         if (.not.LOplusjet) process = ZZWM
         N_lepton = 6
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         finallep_default(5) =  11
         finallep_default(6) = -12
      CASE(ZZhadWM)
         if (LOplusjet) process = ZZhadWMjLO
         if (.not.LOplusjet) process = ZZhadWM
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
      CASE(ZZWMhad)
         if (LOplusjet) process = ZZWMhadjLO
         if (.not.LOplusjet) process = ZZWMhad
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
      CASE(WWWP)
         if (LOplusjet) process = WWWPjLO
         if (.not.LOplusjet) process = WWWP
         N_lepton = 6
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
         finallep_default(5) =  12
         finallep_default(6) = -11
      CASE(WPhadWMWP)
         if (LOplusjet) process = WPhadWMWPjLO
         if (.not.LOplusjet) process = WPhadWMWP
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  12
         finallep_default(4) = -11
      CASE(WPWMhadWP)
         if (LOplusjet) process = WPWMhadWPjLO
         if (.not.LOplusjet) process = WPWMhadWP
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  12
         finallep_default(4) = -11
      CASE(WWWM)
         if (LOplusjet) process = WWWMjLO
         if (.not.LOplusjet) process = WWWM
         N_lepton = 6
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  12
         finallep_default(4) = -11
         finallep_default(5) =  11
         finallep_default(6) = -12
      CASE(WMhadWPWM)
         if (LOplusjet) process = WMhadWPWMjLO
         if (.not.LOplusjet) process = WMhadWPWM
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
      CASE(WMWPhadWM)
         if (LOplusjet) process = WMWPhadWMjLO
         if (.not.LOplusjet) process = WMWPhadWM
         N_lepton = 4
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  11
         finallep_default(4) = -12
      CASE(WWA)
         if (LOplusjet) process = WWAjLO
         if (.not.LOplusjet) process = WWA
         N_lepton = 4
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -12
         PS_Loops = 3
      CASE(WPhadWMA)
         if (LOplusjet) process = WPhadWMAjLO
         if (.not.LOplusjet) process = WPhadWMA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         PS_Loops = 3
      CASE(WPWMhadA)
         if (LOplusjet) process = WPWMhadAjLO
         if (.not.LOplusjet) process = WPWMhadA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         PS_Loops = 3
      CASE(ZZA)
         if (LOplusjet) process = ZZAjLO
         if (.not.LOplusjet) process = ZZA
         N_lepton = 4
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         PS_Loops = 3
      CASE(ZZhadA)
         if (LOplusjet) process = ZZhadAjLO
         if (.not.LOplusjet) process = ZZhadA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         PS_Loops = 3
      CASE(ZZnA)
         if (LOplusjet) process = ZZnAjLO
         if (.not.LOplusjet) process = ZZnA
         N_lepton = 4
         finallep_default(1) =  11
         finallep_default(2) = -11
         finallep_default(3) =  12
         finallep_default(4) = -12
         PS_Loops = 3
      CASE(WPZA)
         if (LOplusjet) process = WPZAjLO
         if (.not.LOplusjet) process = WPZA
         N_lepton = 4
         finallep_default(1) =  12
         finallep_default(2) = -11
         finallep_default(3) =  11
         finallep_default(4) = -11
         PS_Loops = 3
      CASE(WPhadZA)
         if (LOplusjet) process = WPhadZAjLO
         if (.not.LOplusjet) process = WPhadZA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         PS_Loops = 3
      CASE(WPZhadA)
         if (LOplusjet) process = WPZhadAjLO
         if (.not.LOplusjet) process = WPZhadA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         PS_Loops = 3
      CASE(WMZA)
         if (LOplusjet) process = WMZAjLO
         if (.not.LOplusjet) process = WMZA
         N_lepton = 4
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  11
         finallep_default(4) = -11
         PS_Loops = 3
      CASE(WMhadZA)
         if (LOplusjet) process = WMhadZAjLO
         if (.not.LOplusjet) process = WMZhadA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         PS_Loops = 3
      CASE(WMZhadA)
         if (LOplusjet) process = WMZhadAjLO
         if (.not.LOplusjet) process = WMZhadA
         N_lepton = 2
         N_decayquarks = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         finallep_default(3) =  11
         finallep_default(4) = -11
         PS_Loops = 3
      CASE(WPAA)
         if (LOplusjet) process = WPAAjLO
         if (.not.LOplusjet) process = WPAA
         N_lepton = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         IF(anomswitch) THEN
         PS_Loops = 5
         Else
         PS_Loops = 4
         ENDIF
      CASE(WMAA)
         if (LOplusjet) process = WMAAjLO
         if (.not.LOplusjet) process = WMAA
         N_lepton = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         IF(anomswitch) THEN
         PS_Loops = 5
         Else
         PS_Loops = 4
         ENDIF
      CASE(ZAA)
         if (LOplusjet) process = ZAAjLO
         if (.not.LOplusjet) process = ZAA
         N_lepton = 2
         finallep_default(1) =  11
         finallep_default(2) = -11
         IF(anomswitch) THEN
         PS_Loops = 5
         Else
         PS_Loops = 4
         ENDIF
      CASE(ZnAA)
         if (LOplusjet) process = ZnAAjLO
         if (.not.LOplusjet) process = ZnAA
         N_lepton = 2
         finallep_default(1) =  12
         finallep_default(2) = -12
         IF(anomswitch) THEN
         PS_Loops = 4
         Else
         PS_Loops = 1
         ENDIF
      CASE(AAA)
         if (LOplusjet) process = AAAjLO
         if (.not.LOplusjet) process = AAA
         N_lepton = 0
         IF(anomswitch) THEN
         PS_Loops = 4
         Else
         PS_Loops = 1
         ENDIF
      CASE(WPAAj)
         if (LOplusjet) then
           process = WPAAjjLO
           Loops_sub_LO=2
         endif
         if (.not.LOplusjet) then
           process = WPAAj
           Loops_sub_LO=6
         endif
         N_lepton = 2
         finallep_default(1) =  12
         finallep_default(2) = -11
         IF(anomswitch) THEN
         PS_Loops = 5
         Else
         PS_Loops = 4
         ENDIF
         Loops_sub_NLO=2
      CASE(WMAAj)
         if (LOplusjet) then
           process = WMAAjjLO
           Loops_sub_LO=2
         endif
         if (.not.LOplusjet) then
           process = WMAAj
           Loops_sub_LO=6
         endif
         N_lepton = 2
         finallep_default(1) =  11
         finallep_default(2) = -12
         IF(anomswitch) THEN
         PS_Loops = 5
         Else
         PS_Loops = 4
         ENDIF
         Loops_sub_NLO=2
      CASE DEFAULT
         print*," "
         print*,"Please select a valid process ID (procID = ",procID,")"
         print*," "
         stop
      END SELECT


c minimal number of jets in the final state for the different processes
      SELECT CASE(procID)
      CASE(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_WW, Hjj_ZZ_ll,
     1     Hjj_ZZ_lnu, HHjj, HHjj_bbtau, HHjj_bbAA, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj, WPWMjj, WPWPjj, 
     2     WMWMjj, ZZjj_ll, ZZjj_lnu,
     3     WPZjj, WMZjj, HAjj, HAjj_AA, HAjj_mu, HAjj_tau,
     4     HAjj_WW, HAjj_ZZ_ll, HAjj_ZZ_lnu, AAjj,QCDAAjj,
     5     QCDZZjj_ll,QCDZZjj_lnu,QCDWPZjj,QCDWMZjj,QCDWPWPjj,QCDWMWMjj,
     6     QCDZjj_l,QCDZjj_nu,QCDWPjj, QCDWMjj, Sp2jj_WW, Sp2jj_ZZ_ll,
     7     Sp2jj_ZZ_lnu,WPAjj, WMAjj,QCDWPAjj,QCDWMAjj,QCDZAjj_l,QCDZAjj_n,ZAjj,ZAjj_n)
         if (LOplusjet) njet_min = 3
         if (.not.LOplusjet) njet_min = 2
      CASE(WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj,
     &     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, WPhadWPjj, WMhadWMjj)
         if (LOplusjet) njet_min = 5
         if (.not.LOplusjet) njet_min = 4
      CASE(Hjj_bbar,HAjj_bbar)
         if (LOplusjet) njet_min = 5
         if (.not.LOplusjet) njet_min = 4
      CASE(Hjjj,Hjjj_AA, Hjjj_mu, Hjjj_tau, Hjjj_WW, 
     1        Hjjj_ZZ_ll,Hjjj_ZZ_lnu)
         if (LOplusjet) njet_min = 4
         if (.not.LOplusjet) njet_min = 3
      CASE(Hjjj_bbar)
         if (LOplusjet) njet_min = 6
         if (.not.LOplusjet) njet_min = 5
      CASE(GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, 
     -        GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 2
      CASE(GFHjj_bbar)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 4
      CASE(WW, WPZ, WMZ, ZZ, WPA, WMA, ZA, AA,
     &     WWZ, ZZZ, ZZWP, ZZWM, WWWP, WWWM, WWA, ZZA, ZZnA,
     &     WPZA, WMZA, WPAA, WMAA, ZAA, ZnAA, AAA,
     &     WP_only, WM_only,
     &     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu,
     &     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu
     &     )    
         if (LOplusjet) njet_min = 1
         if (.not.LOplusjet) njet_min = 0
      CASE(WPH_bbar, WMH_bbar) !force detected jets in H-> b bbar final state
         if (LOplusjet) njet_min = 3
         if (.not.LOplusjet) njet_min = 2
      CASE(WPhadWMlep,WPlepWMhad,ZZhad,
     &     WPhadZ,WPZhad,WMhadZ,WMZhad,
     &     WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, 
     &     WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
         if (LOplusjet) njet_min = 3
         if (.not.LOplusjet) njet_min = 2
      CASE(GFWW, GFZZ, GFZA, GFAA)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 0
      CASE(GFWWj, GFZZj)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 1
      CASE(GFWPhadWMlep,GFWPlepWMhad,GFZZhad)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 2
      CASE(GFWPHADWMJ, GFWPWMHADj)
         if (LOplusjet) then
            print*,"No LO process with additional Jet implemented"
            print*,"Process without additional jet is calculated"
            LOplusjet = .false.
         endif
         njet_min = 3
      CASE(WWJ, WPZJ, WMZJ, WMAJ, WPAJ, ZZJ,
     &     WPJ, WMJ,
     3     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     3     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu,
     &     WPAAj, WMAAj)
         if (LOplusjet) njet_min = 2
         if (.not.LOplusjet) njet_min = 1
      CASE( WPHJ_bbar, WMHJ_bbar) ! Higgsdecay b quarks are also included in jet algorithm
         if (LOplusjet) njet_min = 4
         if (.not.LOplusjet) njet_min = 3
      CASE(WPHADWMJ, WPWMHADJ, WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ)
         if (LOplusjet) njet_min = 4
         if (.not.LOplusjet) njet_min = 3
      END SELECT

      return
      end


c*************************************************************************  
      subroutine lepton_assignment
c*************************************************************************
c     input and output is transferred via the common blocks in process.inc
c  
c     INPUT:    N_lepton:   number of leptons in the final state
c               finallep:   final state lepton numbers requested by the user
c               finallep_default: default final state lepton numbers
c
c     OUTPUT:   finallep:   ordered and checked lepton numbers         
c     
c*************************************************************************  
      implicit none

#include "global.inc"
#include "process.inc"
#include "mssm.inc"

c local variables
      integer nl, ordered(7), i, j, k, lep(7)
      integer i2

      nl = N_lepton
      do i = 1, N_lepton
         if (finallep(i) .ge. 98) then            
            alllep_flag=.true.
            if (finallep(i) .eq. 98) lepgens = 2
            if (finallep(i) .eq. 99) lepgens = 3
         end if
         lep(i) = finallep(i)
         ordered(i) = 0
      enddo

c leptons number 1, 3, ... 
      do k = 1, N_lepton,2
         i = 1
         do while (i < nl+1)
            if ( ((-1*lep(i)*finallep_default(k)).le.0) .and.
     1         ((abs(lep(i)).eq.abs(finallep_default(k)))   .or.
     2          (abs(lep(i)).eq.abs(finallep_default(k))+2) .or.
     3          (abs(lep(i)).eq.abs(finallep_default(k))+4)) ) then
c now check if we can find its partner
               i2 = 1
               do while (i2 < nl+1)
                  if ( (i.ne.i2) .and.
     &                 ((-1*lep(i2)*finallep_default(k+1)).le.0) .and.
     1                 (abs(abs(lep(i2))-abs(lep(i))).le.1)       .and.
     2                (mod(max(abs(lep(i2)),abs(lep(i))),2).eq.0 .or.
     3                 abs(lep(i2)).eq.abs(lep(i)))   .and.
     4               ((abs(lep(i2)).eq.abs(finallep_default(k+1)))   .or.
     5                (abs(lep(i2)).eq.abs(finallep_default(k+1))+2) .or.
     6                (abs(lep(i2)).eq.abs(finallep_default(k+1))+4)) ) then
c got it ...
                     ordered(k) = lep(i)
                     ordered(k+1) = lep(i2)
c remove first lepton from list
                     do j = 1, i-1
                        lep(j) = lep(j)
                     enddo
                     do j = i, nl
                        lep(j) = lep(j+1)
                        if (j .eq. i2) i2 = i2-1
                     enddo
                     nl = nl - 1
                     i = nl+1
c remove second lepton from list
                     do j = 1, i2-1
                        lep(j) = lep(j)
                     enddo
                     do j = i2, nl
                        lep(j) = lep(j+1)
                     enddo
                     nl = nl - 1
                     i2 = nl+1
                  else
                     i2 = i2+1
                  endif
               enddo
               i = i+1

            else
               i = i+1
            endif
         enddo
      enddo

      j = 1
      nl = N_lepton
      do while(j<nl+1)
         if (ordered(j).eq.0) then
            do i = 1, N_lepton
               finallep(i) = finallep_default(i)
            enddo
            if (.not. alllep_flag) then
               print*," "
               print*,"  Numbers for the final state leptons are not"
               print*,"  valid default values are used: ",finallep_default
               print*," "
            end if
            nl = 0
         else
            do i = 1, N_lepton
               finallep(i) = ordered(i)
            enddo
            j = j+1
          endif
       enddo
       ! the following code is unused
      ! do i = 1, N_lepton
      ! ! note: not valid for processes with hadronic decays:
      ! ! here index of leptons and index of v do not correspond...
      ! ! would have to add a shift (at the moment for WPlepWMhad +2)
      ! ! but atm this is only used in WAAj...
      !    select case(finallep(i))
      !    case(11,-11)
      !       particle_mass(i) = ME
      !    case(13,-13)
      !       particle_mass(i) = MM
      !    case(15,-15)
      !       particle_mass(i) = ML
      !    case default
      !       particle_mass(i) = 0d0
      !    end select
      ! enddo

c$$$      do i = 1, N_lepton
c$$$         write(*,*)'i, finallep(i) =', i, finallep(i)
c$$$      end do
c$$$      stop

      return 
      end


c*************************************************************************  
      subroutine quark_assignment
c*************************************************************************
c     input and output is transferred via the common blocks in process.inc
c  
c     INPUT:    N_decayquarks:       number of quarks in the final state
c               finalquarks:         final state quark numbers requested by the user
c               finalquarks_default: default final state quark numbers
c
c     OUTPUT:   finalquarks:         ordered and checked quark numbers         
c     
c*************************************************************************  
      implicit none

#include "global.inc"
#include "process.inc"

c local variables
      integer ordered(2), i, quarks(2), temp

      ! no process with semileptonic decays
      if (N_decayquarks.ne.2) return

      do i = 1, 2
         quarks(i) = finalquarks(i)
         ordered(i) = 0
      enddo

      ! all combinations is okay for every process
      if (quarks(1).eq.93 .or. quarks(1).eq.94) then
        quarks(2)=quarks(1)
        finalquarks(2)=finalquarks(1)
        return
      endif

      select case (procID)
      ! W+ decay
      case (WPhadWMlep,GFWPhadWMlep,WPhadZ,
     &      WPhadWMjj, WPhadZjj, WPhadWMZ, ZZWPhad, WPhadWMWP, WMWPhadWM, WPhadZJ,
     &      WPhadWMa, WPhadZA,
     &      Hjj_WPhadWM, WPhadWPjj, WPhadWMj, GFWPhadWMj)
         if ( (quarks(1)+quarks(2).eq.1) .and. (abs(quarks(1))+abs(quarks(2)).le.7) 
     &        .and. (mod(quarks(1)**2+quarks(2)**2,5).eq.0) ) then
            ordered(1)=max(quarks(1),quarks(2))
            ordered(2)=min(quarks(1),quarks(2))
         endif
      ! W- decay
      case (WPlepWMhad,GFWPlepWMhad,WMhadZ,
     &      WPWMhadjj, WMhadZjj, WPWMhadZ, ZZWMhad, WPWMhadWP, WMhadWPWM, WMhadZJ,
     &      WPWMhadA, WMhadZA,
     &      Hjj_WPWMhad, WMhadWMjj, WPWMhadj, GFWPWMhadj)
         if ( (quarks(1)+quarks(2).eq.-1) .and. (abs(quarks(1))+abs(quarks(2)).le.7) 
     &        .and. (mod(quarks(1)**2+quarks(2)**2,5).eq.0) ) then
            ordered(1)=max(quarks(1),quarks(2))
            ordered(2)=min(quarks(1),quarks(2))
         endif
      ! Z decay
      case (ZZhad,GFZZhad,WPZhad,WMZhad,
     &      WPZhadjj, WMZhadjj, ZZhadjj, WWZhad, ZZhadWP, ZZhadWM, ZZZhad, WPZhadJ, WMZhadJ,
     &      ZZhadA, WPZhadA, WMZhadA,
     &      Hjj_ZZhad)
         if ( (quarks(1)+quarks(2).eq.0) .and. (abs(quarks(1))+abs(quarks(2)).le.10) 
     &        .and. (abs(quarks(1))+abs(quarks(2)).ne.0)  ) then
            ordered(1)=max(quarks(1),quarks(2))
            ordered(2)=min(quarks(1),quarks(2))
         endif
      case default
         print*, "Hadronic Decay requested for a process where it's not implemented!"
         stop
      end select

      if (ordered(1).eq.0 .or. ordered(2).eq.0) then
         finalquarks(1) = finalquarks_default(1)
         finalquarks(2) = finalquarks_default(2)
         print*," "
         print*,"  Numbers for the final state quarks from boson decays are not valid,"
         print*,"  default values are used: ",finalquarks_default
      else
         finalquarks(1) = ordered(1)
         finalquarks(2) = ordered(2)
      endif

      ! the expected order for the LHE output is quark antiquark (generalized for leptonic
      ! and hadronic decay: particle antiparticle), but for several processes
      ! it is switched. therefore we will switch them as well.
      ! for IDs 93, 94 this switch is also in the lesHouches.F.in and has to be adjusted when
      ! changes occur !!!
      select case (procID)
      case (WPWMhadjj, WPhadZjj, WPZhadjj, WMZhadjj, ZZhadjj, WPhadZJ, WPZhadJ, WPhadWPjj)
        temp = finalquarks(1)
        finalquarks(1) = finalquarks(2)
        finalquarks(2) = temp
      end select

      return 
      end


c*************************************************************************  
      subroutine symmetry_factor
c*************************************************************************
c     input and output is transferred via the common blocks in process.inc
c  
c     INPUT:    N_lepton:   number of leptons in the final state
c               finallep:   final state leptons
c
c     OUTPUT:   symmetryfactor:  symmetry factor for identical particles
c                                in the final state, modified by appropriate
c                                factor if LEPTONS = 99 (any leptons) is 
c                                chosen in vbfnlo.dat    
c     
c*************************************************************************  
      implicit none

#include "process.inc"
#include "global.inc"

      symmetryfactor = 1.0d0

      SELECT CASE(procID) 
      CASE(Zjj_l,Zjj_nu,WPjj,WMjj,WPA,WMA,ZA,GFZA,WMAJ,WPAJ,WMAJJ,WPAJJ,
     &        QCDWPAJJ,QCDWMAjj,QCDZAjj_l,QCDZAjj_n,ZAjj,ZAjj_n,! only all leps
     &        QCDZjj_l,QCDZjj_nu,QCDWPjj,QCDWMjj)
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 3d0   
            if (lepgens .eq. 2) symmetryfactor = 2d0   
         end if
      CASE(Hjj_WW,Hjj_ZZ_lnu,Hjjj_WW,Hjjj_ZZ_lnu,HAjj_WW,HAjj_ZZ_lnu,
     &        WPWMjj,ZZjj_lnu,WPZjj,WMZjj,WW,WPZ,WMZ,WWA,WPZA,WMZA,GFWW,
     &        GFHjj_WW,GFHjj_ZZ_lnu,WMZJ,WPZJ,       ! only all leps
     &        WPH_WW, WPHjLOb_WW,
     &        WMH_WW, WMHjLOb_WW,
     &        WPHJ_WW, WPHJJLO_WW,
     &        WMHJ_WW, WMHJJLO_WW,
     &        WPH_ZZ_lnu, WPHjLOb_ZZ_lnu,    
     &        WMH_ZZ_lnu, WMHjLOb_ZZ_lnu,    
     &        WPHJ_ZZ_lnu, WPHJJLO_ZZ_lnu,    
     &        WMHJ_ZZ_lnu, WMHJJLO_ZZ_lnu,   
     &        QCDWPZjj,QCDWMZjj,QCDZZjj_lnu)
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 9d0  
            if (lepgens .eq. 2) symmetryfactor = 4d0  
         end if
      CASE(WPhadWMjj,WPWMhadjj,WPhadZjj,WPZhadjj,WMhadZjj,WMZhadjj,
     &     WPhadWMlep,WPlepWMhad,WPhadZ,WPZhad,WMhadZ,WMZhad,GFWPhadWMlep,GFWPlepWMhad,
     &     Hjj_WPhadWM, Hjj_WPWMhad,
     &     WMhadZJ,WMZhadJ,WPhadZJ,WPZhadJ,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA) ! only all leps
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 3d0  
            if (lepgens .eq. 2) symmetryfactor = 2d0  
         end if
      CASE(WWZ) ! only all leps
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 27d0  
            if (lepgens .eq. 2) symmetryfactor = 8d0  
         end if
      CASE(WPhadWMZ, WPWMhadZ, WWZhad) ! only all leps
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 9d0  
            if (lepgens .eq. 2) symmetryfactor = 4d0  
         end if
      CASE(Hjj_ZZ_ll,Hjjj_ZZ_ll,HAjj_ZZ_ll,GFHjj_ZZ_ll,
     &     ZZjj_ll,WPWPjj,WMWMjj,ZZ,GFZZ,ZZA, Sp2jj_ZZ_ll, ZZJ, GFZZJ,
     &     QCDWPWPjj,QCDWMWMjj,QCDZZjj_ll)
         if( (finallep(1).eq.finallep(3)) .and.    !ZZ(W) -> 4 e + ..and 
     1        (finallep(2).eq.finallep(4)) ) then  !ZZ(W) -> 4 mu + ..cases
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 4.5d0  
            if (lepgens .eq. 2) symmetryfactor = 2d0  
         end if
      CASE(
     &     WPH_ZZ_ll, WPHjLOb_ZZ_ll,    
     &     WMH_ZZ_ll, WMHjLOb_ZZ_ll,    
     &     WPHJ_ZZ_ll, WPHJJLO_ZZ_ll,    
     &     WMHJ_ZZ_ll, WMHJJLO_ZZ_ll    
     &    )
         if( (finallep(3).eq.finallep(5)) .and.    !ZZ(W) -> 4 e + ..and 
     1        (finallep(4).eq.finallep(6)) ) then  !ZZ(W) -> 4 mu + ..cases
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 4.5d0  
            if (lepgens .eq. 2) symmetryfactor = 2d0  
         end if
      CASE(ZZhadjj,WPhadWPjj,WMhadWMjj,ZZhad,GFZZhad, Hjj_ZZhad)
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 3d0  
            if (lepgens .eq. 2) symmetryfactor = 2d0  
         end if
      CASE(ZZWP,ZZWM)
         if( (finallep(1).eq.finallep(3)) .and.    !ZZ(W) -> 4 e + ..and 
     1        (finallep(2).eq.finallep(4)) ) then  !ZZ(W) -> 4 mu + ..cases
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 13.5d0 
            if (lepgens .eq. 2) symmetryfactor = 4d0 
         end if
      CASE(WWWP,WWWM)
         if( (finallep(1).eq.finallep(5)) .and.
     1        (finallep(2).eq.finallep(6)) ) then
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 13.5d0 
            if (lepgens .eq. 2) symmetryfactor = 4d0 
         end if 
      CASE(WPWMhadWP,WMWPhadWM,ZZWPhad,ZZWMhad,ZZZhad)
         if( (finallep(1).eq.finallep(3)) .and.
     1        (finallep(2).eq.finallep(4)) ) then
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 4.5d0 
            if (lepgens .eq. 2) symmetryfactor = 2d0 
         end if 
      CASE(WPhadWMWP,WMhadWPWM,ZZhadWP,ZZhadWM)
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 9d0 
            if (lepgens .eq. 2) symmetryfactor = 4d0 
         end if 
      CASE(ZZZ)
         if( (finallep(1).eq.finallep(3)) .and.
     1       (finallep(3).eq.finallep(5)) .and.
     2       (finallep(2).eq.finallep(4)) .and.
     3       (finallep(4).eq.finallep(6)) ) then
            symmetryfactor = 1.0d0/6.0d0
         elseif( ( (finallep(1).eq.finallep(3)) .and.
     1             (finallep(2).eq.finallep(4)) ) .or.
     2           ( (finallep(1).eq.finallep(5)) .and.
     3             (finallep(2).eq.finallep(6)) ) .or.
     4           ( (finallep(3).eq.finallep(5)) .and.
     5             (finallep(4).eq.finallep(6)) ) ) then
            symmetryfactor = 0.5d0
         endif
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 4.5d0
            if (lepgens .eq. 2) symmetryfactor = 4d0/3d0
         end if
      CASE(ZZnA) 
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 9d0  
            if (lepgens .eq. 2) symmetryfactor = 4d0  
         end if
      CASE(AA,GFAA,AAjj,QCDAAjj)
         symmetryfactor = 0.5d0
      CASE(WPAA,WMAA,ZAA,ZnAA,WPAAj,WMAAj)
         symmetryfactor = 0.5d0
         if (alllep_flag) then
            if (lepgens .eq. 3) symmetryfactor = 1.5d0
            if (lepgens .eq. 2) symmetryfactor = 1d0
         end if
      CASE(AAA)
         symmetryfactor = 1.0d0/6.0d0
      CASE(HHjj)
         symmetryfactor = 0.5d0
      END SELECT

c$$$      write(*,*)'symmetryfactor =', symmetryfactor
c$$$      write(*,*)'finallep(1) =', finallep(1)
c$$$      write(*,*)'finallep(2) =', finallep(2)
c$$$      write(*,*)'finallep(3) =', finallep(3)
c$$$      write(*,*)'finallep(4) =', finallep(4)
c$$$      stop

      return 
      end


********************************************************************************
********************************************************************************

c Subroutine that randomly generates the lepton numbers if 'ALL LEPTONS'
c (LEPTONS = 99 in vbfnlo.dat) has been chosen

      subroutine lepton_gen

      implicit none

#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/global.inc"

      double precision rn, RandomNumber
      integer i
      external RandomNumber

      ! default leptons are always first generation!
      do i = 1, N_lepton, 2
         rn = RandomNumber()
         finallep(i)   = sign(abs(finallep_default(i)  ) + 2*int(rn*lepgens), finallep_default(i)  )
         finallep(i+1) = sign(abs(finallep_default(i+1)) + 2*int(rn*lepgens), finallep_default(i+1))
      enddo


** Setting final particle ID for LHE output

      SELECT CASE (procID)
      ! vbf processes where 1st boson decays hadronically
      CASE(WPhadWMjj, ZZhadjj,
     &     WPhadZjj, WMhadZjj,
     &     WPhadWPjj,WMhadWMjj,
     &     Hjj_WPhadWM, Hjj_ZZhad)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! vbf processes where 2nd boson decays hadronically
      CASE(WPWMhadjj, Hjj_WPWMhad,
     &     WPZhadjj, WMZhadjj)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      ! diboson processes where 1st boson decays hadronically
      CASE(WPhadWMlep, GFWPhadWMlep,
     &     ZZhad, GFZZhad,
     &     WPhadZ, WMhadZ)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! diboson processes where 2nd boson decays hadronically
      CASE(WPlepWMhad, GFWPlepWMhad,
     &     WPZhad, WMZhad)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      ! triboson processes where 1st boson decays hadronically
      CASE(WPhadWMZ, ZZhadWP, ZZhadWM, WPhadWMWP, WMhadWPWM, ZZZhad)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
         Particle_IDs(5) = finallep(3)
         Particle_IDs(6) = finallep(4)
      ! triboson processes where 2nd boson decays hadronically
      CASE(WPWMhadZ, WPWMhadWP, WMWPhadWM)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
         Particle_IDs(5) = finallep(3)
         Particle_IDs(6) = finallep(4)
      ! triboson processes where 3rd boson decays hadronically
      CASE(WWZhad, ZZWPhad, ZZWMhad)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finallep(3)
         Particle_IDs(4) = finallep(4)
         Particle_IDs(5) = finalquarks(1)
         Particle_IDs(6) = finalquarks(2)
      ! triboson processes with photon where 1st boson decays hadronically
      CASE(WPhadWMA, ZZhadA, WPhadZA, WMhadZA)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! triboson processes with photon where 2nd boson decays hadronically
      CASE(WPWMhadA, WPZhadA, WMZhadA)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      ! diboson + jet processes where 1st boson decays hadronically
      CASE(WPhadWMJ, GFWPhadWMJ, WMhadZJ, WPhadZJ)
         Particle_IDs(1) = finalquarks(1)
         Particle_IDs(2) = finalquarks(2)
         Particle_IDs(3) = finallep(1)
         Particle_IDs(4) = finallep(2)
      ! diboson + jet processes where 2nd boson decays hadronically
      CASE(WPWMhadJ, GFWPWMhadJ, WMZhadJ, WPZhadJ)
         Particle_IDs(1) = finallep(1)
         Particle_IDs(2) = finallep(2)
         Particle_IDs(3) = finalquarks(1)
         Particle_IDs(4) = finalquarks(2)
      CASE DEFAULT
         do i = 1, N_lepton
            Particle_IDs(i) = finallep(i)
            ! write(*,*)'i, particle_ids(i) =', i, particle_ids(i)
         enddo 
      end select

      end


********************************************************************************
********************************************************************************

*** Final print out of basic process information

      subroutine printFinalProc

      implicit none

#include "process.inc"
#include "global.inc"
#include "cuts.inc"
#include "coupl_haddecay.inc"

      character*50 lep(n_lepton+N_decayquarks), speclep, lepout, pin1, pin2, pin
      character*100 pout
      character*50 pouts1, pouts2

      integer i


c FINAL STATE LEPTONS      
      do i = 1, n_lepton + N_decayquarks
         if (alllep_flag) then  ! multiple leptons
            SELECT CASE (particle_ids(i))
            CASE (11,13,15)
               if (lepgens .eq. 3) lep(i) = "l-"
               if (lepgens .eq. 2) lep(i) = "e/mu-"
            CASE (-11,-13,-15)
               if (lepgens .eq. 3) lep(i) = "l+"
               if (lepgens .eq. 2) lep(i) = "e/mu+"
            CASE (12,14,16)
               if (lepgens .eq. 3) lep(i) = "vl"
               if (lepgens .eq. 2) lep(i) = "ve/mu"
            CASE (-12,-14,-16)
               if (lepgens .eq. 3) lep(i) = "vl~"
               if (lepgens .eq. 2) lep(i) = "ve/mu~"
            END SELECT
         else   ! specific lepton picked by user
          SELECT CASE (particle_ids(i))
            CASE (11)
               lep(i) = "e-"
            CASE (-11)
               lep(i) = "e+"
            CASE (12)
               lep(i) = "ve"
            CASE (-12)
               lep(i) = "ve~"
            CASE (13)
               lep(i) = "mu-"
            CASE (-13)
               lep(i) = "mu+"
            CASE (14)
               lep(i) = "vmu"
            CASE (-14)
               lep(i) = "vmu~"
            CASE (15)
               lep(i) = "mtau-"
            CASE (-15)
               lep(i) = "mtau+"
            CASE (16)
               lep(i) = "vtau"
            CASE (-16)
               lep(i) = "vtau~"
            END SELECT
         end if

         SELECT CASE (particle_ids(i))
         CASE (1)
            lep(i) = "d"
         CASE (-1)
            lep(i) = "d~"
         CASE (2)
            lep(i) = "u"
         CASE (-2)
            lep(i) = "u~"
         CASE (3)
            lep(i) = "s"
         CASE (-3)
            lep(i) = "s~"
         CASE (4)
            lep(i) = "c"
         CASE (-4)
            lep(i) = "c~"
         CASE (5)
            lep(i) = "b"
         CASE (-5)
            lep(i) = "b~"
         CASE (93)
            if (mod(i,2).eq.1) then
               lep(i) = "q"
            else
               lep(i) = "q~"
            endif
         CASE (94)
            if (mod(i,2).eq.1) then
               lep(i) = "q/b"
            else
               lep(i) = "q~/b~"
            endif
         CASE (100)
            lep(i) = ""
         END SELECT
      end do

      lepout = ' '
      do i = 1, n_lepton + N_decayquarks
         speclep = lep(i)
         lepout = TRIM(lepout)//" "//
     &        TRIM(speclep)
      end do
c      write(*,*)'lepout =', lepout


c INCOMING PARTICLES
      SELECT CASE(procID)
      CASE(WMWMjj, WMhadWMjj)
        if (sign1 .eq. 1) then
           pin1 = "pbar"
        else
           pin1 = "p"
        end if
        if (sign2 .eq. 1)then
           pin2 = "pbar"
        else
           pin2 = "p"
        end if
      CASE DEFAULT
        if (sign1 .eq. 1) then
           pin1 = "p"
        else
           pin1 = "pbar"
        end if
        if (sign2 .eq. 1)then
           pin2 = "p"
        else
           pin2 = "pbar"
        end if
      END SELECT
      pin = TRIM(pin1)//" "//TRIM(pin2)
      
c prefill pouts{1,2,3} to reduce amount of code below
c for WPH:      
c pouts1 hold H decay information
c pouts2 is switch W+ W-      
      SELECT CASE(procID)
      CASE(WPH, WPHJ, WMH, WMHJ)
         pouts1 = " H"
      CASE(WPH_AA, WPHJ_AA, WMH_AA, WMHJ_AA)
         pouts1 = " A A"
      CASE(WPH_mu, WPHJ_mu, WMH_mu, WMHJ_mu)
         pouts1 = " mu+ mu-"
      CASE(WPH_tau, WPHJ_tau, WMH_tau, WMHJ_tau)
         pouts1 = " tau+ tau-"
      CASE(WPH_bbar, WPHJ_bbar, WMH_bbar, WMHJ_bbar)
         pouts1 = " b bbar"
      CASE(WPH_WW, WPHJ_WW, WMH_WW, WMHJ_WW)
         pouts1 = " W+ W-"
      CASE(WPH_ZZ_ll, WPHJ_ZZ_ll, WMH_ZZ_ll, WMHJ_ZZ_ll, WPH_ZZ_lnu, WPHJ_ZZ_lnu, WMH_ZZ_lnu, WMHJ_ZZ_lnu)
         pouts1 = " Z Z"
      END SELECT

      SELECT CASE(procID)
      CASE(
     3     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     3     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     &     WP_only, WPj
     &    )
         pouts2 = "W+"
      CASE(
     3     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu,
     3     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu,
     &     WM_only, WMj
     &    )
         pouts2 = "W-"
      END SELECT
      

c THE PROCESS ITSELF
      SELECT CASE(procID)
c HIGGS via VBF and gluon fusion
      CASE(Hjj,GFHjj)
         pout = 'H jj'
         if (LOplusjet) pout = "Hjjj"
      CASE(Hjj_AA,GFHjj_AA)
         pout = "H jj --> AA jj"
         if (LOplusjet) pout = "H jjj --> AA jjj"
      CASE(Hjj_mu,GFHjj_mu)
         pout = "H jj --> mu+ mu- jj" 
         if (LOplusjet) pout = "H jjj --> mu+ mu- jjj"          
      CASE(Hjj_tau,GFHjj_tau)
         pout = "H jj --> tau+ tau- jj" 
         if (LOplusjet) pout = "H jjj --> tau+ tau- jjj"
      CASE(Hjj_bbar,GFHjj_bbar)
         pout = "H jj --> b bbar jj"
         if (LOplusjet) pout = "H jjj --> b bbar jjj"
      CASE(Hjj_WW,GFHjj_WW, Hjj_WPhadWM, Hjj_WPWMhad)
         pout = "H jj --> W+ W- jj --> "//trim(lepout)//
     &        " jj"
         if (LOplusjet) pout = "H jjj --> W+ W- jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(Hjj_ZZ_ll,Hjj_ZZ_lnu,GFHjj_ZZ_ll,GFHjj_ZZ_lnu, Hjj_ZZhad)
         pout = "H jj --> Z Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = "H jjj --> Z Z jjj --> "//
     &        trim(lepout)//" jjj"

      CASE(Sp2jj_WW)
       if (with_spin2) then
         pout = "spin2 jj --> W+ W- jj --> "//trim(lepout)//
     &        " jj"
         if (LOplusjet) pout = "spin2 jjj --> W+ W- jjj --> "//
     &        trim(lepout)//" jjj"
        else
         pout = "H jj --> W+ W- jj --> "//trim(lepout)//
     &        " jj"
         if (LOplusjet) pout = "H jjj --> W+ W- jjj --> "//
     &        trim(lepout)//" jjj"
        endif
      CASE(Sp2jj_ZZ_ll,Sp2jj_ZZ_lnu)
       if (with_spin2) then
         pout = "spin2 jj --> Z Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = "spin2 jjj --> Z Z jjj --> "//
     &        trim(lepout)//" jjj"
       else
         pout = "H jj --> Z Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = "H jjj --> Z Z jjj --> "//
     &        trim(lepout)//" jjj"
       endif

c HIGGS PAIR VIA VBF
      CASE(HHjj)
         pout = 'H H jj'
         if (LOplusjet) pout = " H H jjj"
      CASE(HHjj_bbtau)
         pout = "H H jj --> b bbar tau+ tau- jj"
         if (LOplusjet) pout = "H H jjj --> b bbar tau+ tau- jjj"
      CASE(HHjj_bbAA)
         pout = "H H jj --> b bbar A A jj" 
         if (LOplusjet) pout = "H H jjj --> b bbar A A jjj"

c HIGGS + JET VIA VBF
      CASE(Hjjj)
         pout = 'H jjj'
         if (LOplusjet) pout = "H jjjj"
      CASE(Hjjj_AA)
         pout = "H jjj --> A A jjj"
         if (LOplusjet) pout = "H jjjj --> AA jjjj"
      CASE(Hjjj_mu)
         pout = "H jjj --> mu+ mu- jjj" 
         if (LOplusjet) pout = "H jjjj --> mu+ mu- jjjj"          
      CASE(Hjjj_tau)
         pout = "H jjj --> tau+ tau- jjj" 
         if (LOplusjet) pout = "H jjjj --> tau+ tau- jjjj"
      CASE(Hjjj_bbar)
         pout = "H jjj --> b bbar jjj"
         if (LOplusjet) pout = "H jjjj --> b bbar jjjj"
      CASE(Hjjj_WW)
         pout = "H jjj --> W+ W- jjj --> "//trim(lepout)//
     &        " jjj"
         if (LOplusjet) pout = "H jjjj --> W+ W- jjjj --> "//
     &        trim(lepout)//" jjjj"
      CASE(Hjjj_ZZ_ll,Hjjj_ZZ_lnu)
         pout = "H jjj --> Z Z jjj --> "//trim(lepout)//
     &        " jjj"
         if (LOplusjet) pout = "H jjjj --> Z Z jjjj --> "//
     &        trim(lepout)//" jjjj"

c SINGLE VECTOR BOSON via VBF
      CASE(Zjj_l, Zjj_nu)
         pout = " Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " Z jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WPjj)
         pout = " W+ jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W+ jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WMjj)
         pout = " W- jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W- jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(Ajj)
         pout = " A jj"
         if (LOplusjet) pout = " A jjj"

c DIBOSON PRODUCTION via VBF
      CASE(WPWMjj, WPhadWMjj, WPWMhadjj)
         pout = " W+ W- jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W+ W- jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(ZZjj_ll, ZZjj_lnu, ZZhadjj)
         pout = " Z Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " Z Z jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WPZjj, WPhadZjj, WPZhadjj)
         pout = " W+ Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W+ Z jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WMZjj, WMhadZjj, WMZhadjj)
         pout = " W- Z jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W- Z jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WPAjj)
         pout = " W+A jj --> "//trim(lepout)//" A jj"
         if (LOplusjet) pout = " W+A jjj --> "//
     &        trim(lepout)//" A jjj"
      CASE(WMAjj)
         pout = " W-A jj --> "//trim(lepout)//" A jj"
         if (LOplusjet) pout = " W-A jjj --> "//
     &        trim(lepout)//" A jjj"
      CASE(ZAjj,ZAjj_n)
         pout = " Z A jj --> "//trim(lepout)//" A jj"
         if (LOplusjet) pout = " Z A jjj --> "//
     &        trim(lepout)//" A jjj"
      CASE(AAjj)
         pout = " spin2 jj --> A A jj"
         if (LOplusjet) pout = " spin2 jj --> A A jjj"
      CASE(WPWPjj, WPhadWPjj)
         pout = " W+ W+ jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W+ W+ jjj --> "//
     &        trim(lepout)//" jjj"
      CASE(WMWMjj, WMhadWMjj)
         pout = " W- W- jj --> "//trim(lepout)//" jj"
         if (LOplusjet) pout = " W- W- jjj --> "//
     &        trim(lepout)//" jjj"

c DIBOSON PRODUCTION
      CASE(WW, GFWW, WPhadWMlep, WPlepWMhad, GFWPhadWMlep, GFWPlepWMhad)
         pout = " W+ W- --> "//trim(lepout)
         if (LOplusjet) pout = " W+ W- j --> "//
     &        trim(lepout)//" j"
      CASE(WPZ, WPhadZ, WPZhad)
         pout = " W+ Z --> "//trim(lepout)
         if (LOplusjet) pout = " W+ Z j --> "//
     &        trim(lepout)//" j"
      CASE(WMZ, WMhadZ, WMZhad)
         pout = " W- Z --> "//trim(lepout)
         if (LOplusjet) pout = " W- Z j --> "//
     &        trim(lepout)//" j"
      CASE(ZZ, GFZZ, ZZhad, GFZZhad)
         pout = " Z Z --> "//trim(lepout)
         if (LOplusjet) pout = " Z Z j --> "//
     &        trim(lepout)//" j"
      CASE(WPA)
         pout = " W+ A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " W+ A j --> "//
     &        trim(lepout)//" A j"
      CASE(WMA)
         pout = " W- A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " W- A j --> "//
     &        trim(lepout)//" j"

      ! Higgsstrahlung
      CASE(
     &   WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar,
     &   WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar
     &    )
         pout = " "//TRIM(pouts2)//" H --> "//TRIM(lepout)//TRIM(pouts1)
         if (LOplusjet) pout = " "//TRIM(pouts2)//"H j --> "//
     &        TRIM(lepout)//TRIM(pouts1)//" j"
      CASE(
     3   WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     3   WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu
     &    )
         pout = " "//TRIM(pouts2)//" H --> "
     &    //TRIM(pouts2)//" "//TRIM(pouts1)//" --> "//TRIM(lepout)
         if (LOplusjet) pout = " "//TRIM(pouts2)//"H j --> "
     &    //TRIM(pouts2)//" "//TRIM(pouts1)//" j --> "//
     &        TRIM(lepout)//" j"

      CASE(ZA, GFZA)
         pout = " Z A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " Z A j --> "//
     &        trim(lepout)//" j"
      CASE(AA, GFAA)
         pout = " A A "
         if (LOplusjet) pout = " AA j"

c single W production
      CASE(WP_only, WM_only)
         pout = " "//TRIM(pouts2)//" --> "//TRIM(lepout)//""
         if (LOplusjet) pout = " "//TRIM(pouts2)//" j --> "//
     &        TRIM(lepout)//" j"

      CASE(WPj, WMj)
         pout = " "//TRIM(pouts2)//" j --> "//TRIM(lepout)//" j"
         if (LOplusjet) pout = " "//TRIM(pouts2)//" jj --> "//
     &        TRIM(lepout)//" jj"

c TRIBOSON PRODUCTION
      CASE(WWZ, WPhadWMZ, WPWMhadZ, WWZhad)
         pout = " W+ W- Z --> "//trim(lepout)
         if (LOplusjet) pout = " W+ W- Z j --> "//
     &        trim(lepout)//" j"
      CASE(ZZWP, ZZhadWP, ZZWPhad)
         pout = " Z Z W+ --> "//trim(lepout)
         if (LOplusjet) pout = " Z Z W+ j --> "//
     &        trim(lepout)//" j"
      CASE(ZZWM, ZZhadWM, ZZWMhad)
         pout = " Z Z W- --> "//trim(lepout)
         if (LOplusjet) pout = " Z Z W- j --> "//
     &        trim(lepout)//" j"
      CASE(WWWP, WPhadWMWP, WPWMhadWP)
         pout = " W+ W- W+ --> "//trim(lepout)
         if (LOplusjet) pout = " W+ W- W+ j --> "//
     &        trim(lepout)//" j"
      CASE(WWWM, WMhadWPWM, WMWPhadWM)
         pout = " W+ W- W- --> "//trim(lepout)
         if (LOplusjet) pout = " W+ W- W- j --> "//
     &        trim(lepout)//" j"
      CASE(ZZZ, ZZZhad)
         pout = " ZZZ --> "//trim(lepout)
         if (LOplusjet) pout = " ZZZ j --> "//
     &        trim(lepout)//" j"
      CASE(WWA, WPhadWMA, WPWMhadA)
         pout = " W+ W- A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " W+ W- A j --> "//
     &        trim(lepout)//" A j"
      CASE(ZZA, ZZhadA, ZZnA)
         pout = " ZZA --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " ZZA j --> "//
     &        trim(lepout)//" A j"
      CASE(WPZA, WPhadZA, WPZhadA)
         pout = " W+ Z A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " W+ Z A j --> "//
     &        trim(lepout)//" A j"
      CASE(WMZA, WMhadZA, WMZhadA)
         pout = " W- Z A --> "//trim(lepout)//" A"
         if (LOplusjet) pout = " W- Z A j --> "//
     &        trim(lepout)//" A j"
      CASE(WPAA)
         pout = " W+ A A --> "//trim(lepout)//" A A"
         if (LOplusjet) pout = " W+ A A j --> "//
     &        trim(lepout)//" A A j"
      CASE(WMAA)
         pout = " W- A A --> "//trim(lepout)//" A A"
         if (LOplusjet) pout = " W- A A j --> "//
     &        trim(lepout)//" A A j"
      CASE(ZAA,ZnAA)
         pout = " Z A A --> "//trim(lepout)//" A A"
         if (LOplusjet) pout = " Z A A j --> "//
     &        trim(lepout)//" A A j"
      CASE(AAA)
         pout = " A A A "
         if (LOplusjet) pout = " A A A j "

c TRIBOSON PLUS JET
      CASE(WPAAj)
         pout = " W+ A A j --> "//trim(lepout)//" A A j"
         if (LOplusjet) pout = " W+ A A jj --> "//
     &        trim(lepout)//" A A jj"
      CASE(WMAAj)
         pout = " W- A A j --> "//trim(lepout)//" A A j"
         if (LOplusjet) pout = " W- A A jj --> "//
     &        trim(lepout)//" A A jj"

c HIGGS PLUS PHOTON via VBF
      CASE(HAjj)
         pout = 'H A jj'
         if (LOplusjet) pout = "H A jjj"
      CASE(HAjj_AA)
         pout = "H A jj --> A A A jj"
         if (LOplusjet) pout = "H A jjj --> A A A jjj"
      CASE(HAjj_mu)
         pout = "H A jj --> mu+ mu- A jj" 
         if (LOplusjet) pout = "H A jjj --> mu+ mu- A jjj"          
      CASE(HAjj_tau)
         pout = "H A jj --> tau+ tau- A jj" 
         if (LOplusjet) pout = "H A jjj --> tau+ tau- A jjj"
      CASE(HAjj_bbar)
         pout = "H A jj --> b bbar A jj"
         if (LOplusjet) pout = "H A jjj --> b bbar A jjj"
      CASE(HAjj_WW)
         pout = "H A jj --> W+ W- A jj --> "//
     &        trim(lepout)//" A jj"
         if (LOplusjet) pout = "H A jjj --> W+ W- A jjj --> "//
     &        trim(lepout)//" A jjj"
      CASE(HAjj_ZZ_ll,HAjj_ZZ_lnu)
         pout = "H A jj --> Z Z A jj --> "//trim(lepout)//
     &        " A jj"
         if (LOplusjet) pout = "H A jjj --> Z Z A jjj --> "//
     &        trim(lepout)//" A jjj"

c DIBOSON PLUS JET
      CASE(WPZj)
         pout = " W+ Z j --> "//trim(lepout)//" j"
         if (LOplusjet) pout = " W+ Z jj --> "//
     &        trim(lepout)//" jj"
      CASE(WMZj)
         pout = " W- Z j --> "//trim(lepout)//" j"
         if (LOplusjet) pout = " W- Z jj --> "//
     &        trim(lepout)//" jj"
      CASE(WPAj)
         pout = " W+ A j --> "//trim(lepout)//" A j"
         if (LOplusjet) pout = " W+ A jj --> "//
     &        trim(lepout)//" A jj"
      CASE(WMAj)
         pout = " W- A j --> "//trim(lepout)//" A j"
         if (LOplusjet) pout = " W- A jj --> "//
     &        trim(lepout)//" A jj"
      CASE(WWj,WPHADWMJ,WPWMHADJ,GFWWj,GFWPHADWMJ,GFWPWMHADJ)
         pout = " W+ W- j --> "//trim(lepout)//" j"
         if (LOplusjet) pout = " W+ W- jj --> "//
     &        trim(lepout)//" jj"
      CASE(ZZj,GFZZj)
         pout = " Z Z j --> "//trim(lepout)//" j"
         if (LOplusjet) pout = " Z Z jj --> "//
     &        trim(lepout)//" jj"

      CASE(
     &   WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar,
     &   WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar
     &    )
         pout = " "//TRIM(pouts2)//" H j --> "//TRIM(lepout)//TRIM(pouts1)//" j"
         if (LOplusjet) pout = " "//TRIM(pouts2)//"H jj --> "//
     &        TRIM(lepout)//TRIM(pouts1)//" jj"
      CASE(
     3   WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     3   WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &    )
         pout = " "//TRIM(pouts2)//" H j --> "
     &    //TRIM(pouts2)//" "//TRIM(pouts1)//" j --> "//TRIM(lepout)//" jj"
         if (LOplusjet) pout = " "//TRIM(pouts2)//"H jj --> "
     &    //TRIM(pouts2)//" "//TRIM(pouts1)//" jj --> "//
     &        TRIM(lepout)//" jj"

c QCD induced Vjj production
      CASE(QCDWPjj)
         pout = " W+ jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W+ jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"
      CASE(QCDWMjj)
         pout = " W- jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W- jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"

      CASE(QCDZjj_l, QCDZjj_nu)
         pout = " Z jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " Z jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"


c QCD induced VVjj production
      CASE(QCDAAjj)
         pout = " A A jj (QCD induced)"
         if (LOplusjet) pout = " A A jjj (QCD induced)"
      CASE(QCDZZjj_ll, QCDZZjj_lnu)
         pout = " Z Z jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " Z Z jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"
      CASE(QCDWPZjj)
         pout = " W+ Z jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W+ Z jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"
      CASE(QCDWMZjj)
         pout = " W- Z jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W- Z jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"
      CASE(QCDWPAjj)
         pout = " W+ A jj --> "//trim(lepout)//" A jj (QCD induced)"
         if (LOplusjet) pout = " W+ A jjj --> "//
     &        trim(lepout)//" A jjj (QCD induced)"
      CASE(QCDWMAjj)
         pout = " W- A jj --> "//trim(lepout)//" A jj (QCD induced)"
         if (LOplusjet) pout = " W- A jjj --> "//
     &        trim(lepout)//" A jjj (QCD induced)"
      CASE(QCDZAjj_l)
         pout = " Z A jj --> "//trim(lepout)//" A jj (QCD induced)"
         if (LOplusjet) pout = " Z A jjj --> "//
     &        trim(lepout)//" A jjj (QCD induced)"
      CASE(QCDZAjj_n)
         pout = " Z A jj --> "//trim(lepout)//" A jj (QCD induced)"
         if (LOplusjet) pout = " Z A jjj --> "//
     &        trim(lepout)//" A jjj (QCD induced)"
      CASE(QCDWPWPjj)
         pout = " W+ W+ jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W+ W+ jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"
      CASE(QCDWMWMjj)
         pout = " W- W- jj --> "//trim(lepout)//" jj (QCD induced)"
         if (LOplusjet) pout = " W- W- jjj --> "//
     &        trim(lepout)//" jjj (QCD induced)"


      END SELECT


      write(*,*)'  '
      if ((procID .ge. 4000) .and. (procID .le. 4999)) then
         write(*,*)' GLUON FUSION (1-loop):'
      end if

      write(*,*)' PROCESS: ', procID, "  :   ", trim(pin), 
     &     " --> ", trim(pout)

      if (with_anom .or. with_anomHiggs) 
     &     write(*,*)'   with anomalous couplings'
      if (semileptonic) then
        select case (nlo_semilep)
        case (0)
          write(*,*)'   with leading order hadronic decay of one vector boson'
        case (1)
          write(*,*)'   with approximative next-to leading order hadronic decay of one vector boson'
        end select
      endif
      if (with_kk) write(*,*)'   running in Kaluza-Klein mode'
      if (with_spin2) write(*,*)'   with spin-two particles'
      if (model .eq. 2) write(*,*)'   running in the MSSM'

      write(*,*)'  '


      end
