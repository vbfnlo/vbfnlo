c*************************************************************************
      SUBROUTINE InitGlobalParameters
c*************************************************************************
c     Initialize the global parameters which are defined in "global.inc"
c*************************************************************************

      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/coupl_haddecay.inc"
     
       integer i, lhapdfused    

      double precision alphas5_hardwired
      external alphas5_hardwired

c default values for the grid names for the moment here
       character*50 gridname1_default(1:max_ps), 
     &              gridname2_default(1:max_ps),
     &              gridname3_default(1:max_ps),
     &              gridname4_default(1:max_ps)


c Assign values to the global parameters here.
c Read in values from 'vbfnlo.dat'.

      doBorn=.true.
      doVirtuals=.false.
      doEmissions=.false.
      doFLoops=.false.
      
      call loadfile("vbfnlo.dat",.true.)
      print *," "
      print *,"       Information on global parameters"
      print *,"-----------------------------------------------"
      call read_int("PROCESS",procID,100)


      if ( (procID.ge.4100) .and. (procID.le.4999) ) then
         print*,"  SELECTED PROCESS NUMBER NOT APPROPRIATE" //
     &        "FOR PROGRAM VBFNLO. "
         print*,"  FOR GLUON FUSION PROCESSES PLEASE RUN THE" //
     &        "PROGRAM GGFLO! "
         STOP
      endif

      ! switch whether W, Z, H should be forced on-shell
      !   0: full width dependence
      !   1: force on-shell
      narrowwidth = 0

c     Validity information for the NW implementation can be found
c     in the developer wiki. 
c     It's not officially supported and has several limitations!
!      SELECT CASE(procID)
!      CASE(WWZ, ZZZ, ZZWP, ZZWM, WWWP, WWWM)                 ! implemented in SM and anom
!         call read_int("NARROW_WIDTH",narrowwidth,0)
!      CASE(WPAA, WMAA, WPAAj, WMAAj, 
!     &     WP_only, WM_only, WPj, WMJ,           ! implemented in SM
!     7     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
!     7     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu,
!     7     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
!     7     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
!     &    )
!         call read_int("NARROW_WIDTH",narrowwidth,0)
!      CASE(ZAA)                                              ! implemented only in PS generator for SM
!         call read_int("NARROW_WIDTH",narrowwidth,0)
!      CASE(ZA,GFZA, WMA, WPA, ZZ, GFZZ, WMZ, WPZ, WW, GFWW)  ! implemented only in PS generator
!         call read_int("NARROW_WIDTH",narrowwidth,0)
!      END SELECT

      call read_logical("LOPROCESS_PLUS_JET",LOplusjet,.false.)
c      call proc_assignment(procID)
      call proc_assignment
      alllep_flag = .false.
      if (N_lepton.gt.0) then
         call read_IntList("LEPTONS",N_lepton,finallep,finallep_default)
         call lepton_assignment
      endif
      if (N_decayquarks.gt.0) then
         call read_IntList("DECAY_QUARKS",N_decayquarks,finalquarks,finalquarks_default)
         call quark_assignment
      endif
      call symmetry_factor

      call read_int("LO_ITERATIONS",N_iterations1,4)
      call read_int("NLO_ITERATIONS",N_iterations2,4)
      call read_int("LO_POINTS",N_points(0),20)
      call read_int("NLO_POINTS",N_points(-1),20)
      call read_logical("NLO_SWITCH",doNLO,.false.)
      call read_logical("EWCOR_SWITCH",ewcor_switch,.false.)
      call read_int("FERMIONLOOP",floops,3)

      if (N_decayquarks.gt.0) then
        call read_int("NLO_SEMILEP_DECAY",nlo_semilep,0)
      else
        nlo_semilep = 0
      endif

      call read_logical("LHA_SWITCH",lha,.false.)   
      call read_logical("HEPMC_SWITCH",hepmc,.false.)   
#ifndef WITH_HEPMC
      if (hepmc) then
        write(*,*)'Sorry, HepMC output not included.'
        write(*,*)'Reconfigure with --with-hepmc '
        write(*,*)'  or switch off HEPMC_SWITCH'
        hepmc = .false.
        stop
      endif
#endif
      if (lha.or.hepmc) then
         if (ewcor_switch .or. doNLO) then
            write(*,*)'Sorry, we cannot output LHA or HepMC event files'
            write(*,*)'when we are not working at leading order!'
            lha = .false.
            hepmc = .false.
            stop
         end if
         if (lha) call read_string("LHA_FILE",lhaFILE,"event.lhe")
         if (hepmc) call read_string("HEPMC_FILE",hepmcFILE,
     -        "event.hepmc")
         call read_logical("UNWEIGHTING_SWITCH",unweighting,.true.)
         desiredEventCount = 0
         if (unweighting) call read_int("DESIRED_EVENT_COUNT",desiredEventCount,0)
         partialUnweight = .false.
         if (desiredEventCount.gt.0) call read_logical("PARTIAL_UNWEIGHTING",partialUnweight,.false.)
         call read_logical("TAUMASS",taumass_switch,.false.)
      else
         unweighting = .false.
         desiredEventCount = 0
         partialUnweight = .false.
         taumass_switch = .false.
      end if

       do i = 1,PS_Loops
          gridname1_default(i) = "grid2"
          gridname2_default(i) = "grid3"
          gridname3_default(i) = "grid4"
          gridname4_default(i) = "grid5"
       enddo

      call read_stringList("LO_GRID",PS_Loops,gridname1,gridname1_default)
      call read_stringList("NLO_GRID",PS_Loops,gridname2,gridname2_default)
      call read_stringList("PHTN_GRID",PS_Loops,gridname3,gridname3_default)
      call read_stringList("FLOOP_GRID",PS_Loops,gridname4,gridname4_default)

      call read_Real("ECM",ecm,13000d0)
      call read_int("BEAM1",sign1,1)
      call read_int("BEAM2",sign2,1)

      call read_int("ID_MUF",id_muf,0)
      call read_int("ID_MUR",id_mur,0)
      call read_Real("MUF_USER",muf_user,100d0)
      call read_Real("MUR_USER",mur_user,100d0)
      call read_Real("XIF",xif,1d0)
      call read_Real("XIR",xir,1d0)
      
c read in some hidden options 
      call read_logical_hidden("BORN_OFF", born_off,.false.)
      call read_logical_hidden("VIRTUALS_OFF", virtuals_off,.false.) 
      call read_logical_hidden("EMISSIONS_OFF", emissions_off,.false.)
      call read_logical_hidden("REDUCE_SUB_STAT", reduceSubStat,.true.)

      if (LOplusjet) then
         doNLO = .false.
         ewcor_switch = .false.
      endif

#ifdef WITH_NLO
      if (doNLO.and.(.not.LOplusjet)) then
         doVirtuals=.true.
         doEmissions=.true.
         doFLoops=.true.
      endif
      
      doBorn=doBorn .and. .not. born_off.or..not.doNLO
      doVirtuals=doVirtuals .and. .not. virtuals_off
      doEmissions=doEmissions .and. .not. emissions_off
#else
      print *," "
      print *," NLO features are disabled. "
      print *," Processes are available only at LO level."
      if (doNLO.and.(.not.LOplusjet)) then
         print *," "
         print *," Error: selected NLO process is not available."
         stop
      endif
#endif

      call closefile

      if (nlo_semilep.gt.0 .and. .not.doNLO) then
        call loadfile("vbfnlo.dat",.false.)
        call read_int("PDF_SWITCH",lhapdfused,0)
        call closefile
        if (lhapdfused.eq.1) then
          ! with lhapdf simultaneous use of LO and NLO pdf is not possible.
          ! NLO effects in semileptonic decays should be with NLO pdfs
          ! => use internal alphas-routine
          print *," "
          print *," Using internal alpha_s routine for NLO effects in semileptonic decays as"
          print *," LHAPDF does not provide NLO and LO alpha_s simultaneously:"
          print *," NLO alpha_s(mZ = 91.187 GeV) =", alphas5_hardwired(91.187d0**2,1)
          print *," "
        endif
      endif

      RETURN
      END


c*************************************************************************
      SUBROUTINE InitGlobalParameters_GGFLO
c*************************************************************************
c     Initialize the global parameters which are defined in "global.inc"
c*************************************************************************

      implicit none
#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/lha.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/process.inc"
     
c Assign values to the global parameters here.
c Read in values from 'vbfnlo.dat'.

      doBorn=.true.
      doVirtuals=.false.
      doEmissions=.false.
      doFLoops=.true.
      infloops = .true.

      call loadfile("vbfnlo.dat",.true.)
      print *," "
      print *,"       Information on global parameters"
      print *,"-----------------------------------------------"
      call read_int("PROCESS",procID,100)
      if ( (procID.lt.4100) .or. (procID.gt.4999) ) then
c     if you change this condition also change the correlated if statement in InitGlobalParameters for vbfnlo!
         print*,"  SELECTED PROCESS NUMBER NOT APPROPRIATE"//
     &        " FOR PROGRAM GGFLO. "
         print*,"  THIS PROGRAM ONLY WORKS WITH GLUON FUSION"//
     &        " PROCESSES. "
         print*,"  PLEASE RUN THE PROGRAM VBFNLO! "
         STOP
      endif
      narrowwidth = 0
c      call read_logical("LOPROCESS_PLUS_JET",LOplusjet,.false.)
      call proc_assignment
      if (N_lepton.gt.0) then
         call read_IntList("LEPTONS",N_lepton,finallep,finallep_default)
         call lepton_assignment
      endif
      if (N_decayquarks.gt.0) then
         call read_IntList("DECAY_QUARKS",N_decayquarks,finalquarks,finalquarks_default)
         call quark_assignment
      endif
      call symmetry_factor

      call read_int("LO_ITERATIONS",N_iterations1,4)
      call read_int("NLO_ITERATIONS",N_iterations2,4)
      call read_int("LO_POINTS",N_points(0),15)
      call read_int("NLO_POINTS",N_points(-1),15)
      call read_logical("NLO_SWITCH",doNLO,.false.)
      call read_int("FERMIONLOOP",floops,3)
      if ( (floops .eq. 0) .and. 
     &    ( (procID.ge.4300) .and. (procID.le.4379) ) ) then
        write(*,*)'Diboson production via gluon-fusion requested'
        write(*,*)'but fermionic loops switched off.'
        write(*,*)'Result will be zero.'
        stop
      endif
      call read_string("LO_GRID",gridname1,"grid2")
      call read_string("NLO_GRID",gridname2,"grid3")
      call read_logical("LHA_SWITCH",lha,.false.)
      call read_logical("HEPMC_SWITCH",hepmc,.false.)   
#ifndef WITH_HEPMC
      if (hepmc) then
        write(*,*)'Sorry, HepMC output not included.'
        write(*,*)'Reconfigure with --with-hepmc '
        write(*,*)'or switch off HEPMC_SWITCH'
        hepmc = .false.
        stop
      endif
#endif
      if (lha) call read_string("LHA_FILE",lhaFILE,"event.lhe")
      if (hepmc) call read_string("HEPMC_FILE",hepmcFILE,"event.hepmc")
      if (lha.or.hepmc) then
        call read_logical("UNWEIGHTING_SWITCH",unweighting,.true.)
        call read_logical("TAUMASS",taumass_switch,.false.)
        desiredEventCount = 0
        if (unweighting) call read_int("DESIRED_EVENT_COUNT",desiredEventCount,0)
        partialUnweight = .false.
        if (desiredEventCount.gt.0) call read_logical("PARTIAL_UNWEIGHTING",partialUnweight,.false.)
      else
        unweighting = .false.
        desiredEventCount = 0
        partialUnweight = .false.
        taumass_switch = .false.
      endif
     
      call read_Real("ECM",ecm,14000d0)
      call read_int("BEAM1",sign1,1)
      call read_int("BEAM2",sign2,1)

      call read_int("ID_MUF",id_muf,0)
      call read_int("ID_MUR",id_mur,0)
      call read_Real("MUF_USER",muf_user,100d0)
      call read_Real("MUR_USER",mur_user,100d0)
      call read_Real("XIF",xif,1d0)
      call read_Real("XIR",xir,1d0)

      if (LOplusjet) then
         doNLO = .false.
      endif

#ifdef WITH_NLO
      if (doNLO.and.(.not.LOplusjet)) then
         doVirtuals=.true.
         doEmissions=.true.
      endif
#else
      print *," "
      print *," NLO features are disabled. "
      print *," Processes are available only at LO level."
      if (doNLO.and.(.not.LOplusjet)) then
         print *," "
         print *," Error: selected NLO process is not available."
         stop
      endif
#endif

      call closefile

      RETURN
      END


***********************************
