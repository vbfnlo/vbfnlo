 
****************************************************************************
*   This file contains functions and subroutines for the cuts.
c***************************************************************************
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     FUNCTION Cuts(jets, njets, leptons, nl,invisible, ninv, photons, nphotons, L)
c     FUNCTION Basic_Cuts(partons, npartons, jets, njets, leptons, nl, invisible, ninv, photons, nphotons)
c     FUNCTION VBF_Cuts(jets, njets, leptons, nl, photons, nphotons)
c     FUNCTION JetVeto_Cuts(jets, njets, leptons, nl, photons, nphotons)
c     SUBROUTINE DefineLeptons(v, Lmax, leptons, n_leptons)
c     SUBROUTINE DefinePhotons(v, Lmax, photons, n_photons)
c     SUBROUTINE DefineInvisible(v, Lmax, invisible, n_invisible)
c
c***************************************************************************
*   Last modified: 02.02.2011
****************************************************************************
c*************************************************************************
      SUBROUTINE InitCuts
          use globalvars, only: lglobalprint
          use readinput
c*************************************************************************
c     Reads in the values for the cuts.
c*************************************************************************
      implicit none
#include "global.inc"
#include "process.inc"
#include "cuts.inc"
#include "lha.inc"

      integer i, n_lcharged, n_photon

c     logical variable that checks if process is a VBF process
c     (common block tagjetdef in cuts.inc)
      SELECT CASE(procID)
      CASE(Hjj, Hjj_AA, AAjj, Hjj_mu, Hjj_tau, Hjj_bbar, Hjj_WW,
     &     Hjj_ZZ_ll, Hjj_ZZ_lnu, Ajj, Zjj_l, Zjj_nu, WPjj, WMjj,
     &     WPWMjj, WPWPjj, WMWMjj, ZZjj_ll, ZZjj_lnu, WPZjj, WMZjj, 
     &     WPAjj, WMAjj, ZAjj, ZAjj_n, HAjj, 
     &     HAjj_AA, HAjj_mu, HAjj_tau, HAjj_bbar, HAjj_WW, 
     &     HAjj_ZZ_ll, HAjj_ZZ_lnu,
     &     HHjj, HHjj_bbtau, HHjj_bbAA,
     &     GFHjj, GFHjj_AA, GFHjj_mu, GFHjj_tau, GFHjj_bbar,
     &     GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu, 
     &     QCDWPWPjj,QCDWMWMjj,QCDWPZjj,QCDWMZjj,QCDZZjj_ll,QCDZZjj_lnu,QCDAAjj,
     &     QCDWPjj, QCDWPAjj, QCDWMAjj,QCDWMjj, QCDZjj_l,QCDZAjj_l,QCDZAjj_n,
     &     QCDZjj_nu,
     &     Hjjj, Hjjj_AA, Hjjj_mu,
     &     Hjjj_tau, Hjjj_bbar, Hjjj_WW, Hjjj_ZZ_ll, Hjjj_ZZ_lnu,
     &     WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj,
     &     ZZhadjj, WPhadWPjj, WMhadWMjj,
     &     Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad, 
     &     Sp2jj_WW, Sp2jj_ZZ_ll, Sp2jj_ZZ_lnu)
         vbfprocess = .true.
      case default
         vbfprocess = .false.
      END SELECT

c     logical variable that checks if process has semileptonic decays
c     (common block semilepcuts in cuts.inc)
      SELECT CASE(procID)
      CASE(WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj,
     &    ZZhadjj, WPhadWPjj, WMhadWMjj, WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP,
     &    ZZWPhad, ZZhadWM, ZZWMhad, WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM,
     &    ZZZhad, WPhadWMlep, WPlepWMhad, GFWPhadWMlep, GFWPlepWMhad, ZZhad,
     &    GFZZhad, WPhadZ, WPZhad, WMhadZ, WMZhad,
     &    Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad,
     &    WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     &    WPHADWMJ,WPWMHADJ,WMhadZJ,WMZhadJ,WPhadZJ,WPZhadJ)
         semileptonic = .true.
      case default
         semileptonic = .false.
      END SELECT

c     logical variable that checks if process had a hadronically decaying Z
c     (common block techcuts in process.inc)
      SELECT CASE(procID)
      CASE(WPZhadjj, WMZhadjj, ZZhadjj, WWZhad, ZZhadWP, ZZhadWM, 
     &     ZZZhad, ZZhad, GFZZhad, WPZhad, WMZhad, Hjj_ZZhad,
     &     ZZhadA, WPZhadA, WMZhadA,WMZhadJ,WPZhadJ)
         zdecay_had = .true.           ! true for all processes with hadronically decaying Z 
      case default
         zdecay_had = .false.
      end select

      call loadfile("cuts.dat",.false.)

c jet cuts

      do i=1,max_jets
         ptmin_def_default(i)=20.d0
      enddo
      call read_realList("PT_JET_MIN",max_jets,ptmin_def,ptmin_def_default,.true.)
      call read_real("Y_P_MAX",ymax_part,5d0)
      call read_Real("PGENKTJET",jetgenktp,1d0) ! 1:kT, 0:C/A, -1:anti-kT 
      call read_real("Y_JET_MAX",ymax_def,4.5d0)
      call read_real("RJJ_MIN",rjjmin,0.8d0)

c leptons cuts
      n_lcharged = 0
      do i=1,max_v
         ptlmin_def_default(i)=10.d0
*         ptlmax_def_default(i)=ecm
          if ( (abs(particle_IDs(i)).eq.11) ! check if it is a lepton
     &         .or.(abs(particle_IDs(i)).eq.13)
     &         .or.(abs(particle_IDS(i)).eq.15) ! taus!!
     &         ) then
            n_lcharged = n_lcharged + 1
          endif
      enddo
      call read_real("RLL_MIN",rll_min,0.4d0)
      call read_real("RLL_MAX",rll_max,50.0d0)
      call read_real("RJL_MIN",rjl_min,0.6d0)
      call read_real("Y_L_MAX", ymax_l,2.5d0)
      call read_realList("PT_L_MIN",max_v,ptlmin,ptlmin_def_default,.true.)
*      call read_realList("PT_L_MAX",max_v,ptlmax,ptlmax_def_default,.true.)
      call read_real("MLL_MIN",mllmin,15d0)
      call read_real("MLL_MAX",mllmax,14000d0)
      call read_logical("MLL_OSONLY",mllosonly,.true.)
*      call read_logical("NEUTRINO_CUT",lneutrinoCut,.false.)
*      call read_real("PT_IN_MIN", ptinmin, 10d0)
      call read_real("PTMISS_MIN", ptmiss_min, 0d0)
*      call read_real("MZA_MIN", mza_min, 0d0)
!      ! Cut on neutrino pT, default 0 to be compatible with old
!      ! configuration files (which don't contain PT_IN_MIN)
!      call read_real("PT_IN_MIN", ptinmin, 0d0)

c photon cuts
      n_photon = 0
      do i=1,max_v
         ptlmin_def_default(i)=20.d0
*         ptlmax_def_default(i)=ecm
          if ( (abs(particle_IDs(i)).eq.22) ! check if it is a photon
     &         ) then
            n_photon = n_photon + 1
          endif
      enddo
      call read_real("RGG_MIN",rgg_min,0.6d0)
      call read_real("RGG_MAX",rgg_max,50.0d0)
      call read_real("RJG_MIN",rjg_min,0.6d0)
      call read_real("RLG_MIN",rlg_min,0.6d0)
      call read_real("Y_G_MAX",ymax_g,1.5d0)
      call read_realList("PT_G_MIN",max_v,ptgmin,ptlmin_def_default,.true.)
*      call read_realList("PT_G_MAX",max_v,ptgmax,ptlmax_def_default,.true.)
      call read_real("PHISOLCUT",isol_phopart,0.7d0)
      call read_real("EFISOLCUT",eff_phopart,1d0)
      call read_real("MLG_MIN",mlAmin,0d0)
      call read_real("MLG_MAX",mlAmax,14000d0)
*      call read_int("PHOTONHEL",photon_hel,0)
*      call read_real("MT_MIN",MTmin,0d0)
*      call read_real("MT_MAX",MTmax,14000d0)

c VBF cuts - only applied to tagging jets
      call read_real("ETAJJ_MIN",etajjmin, 0d0)
      call read_logical("YSIGN",ysign,.false.)
      call read_logical("LRAPIDGAP",lrapidgap,.false.)
      call read_real("DELY_JL",dely_jlep, 0d0)
      call read_logical("GRAPIDGAP",grapidgap,.false.)
      call read_real("DELY_JG",dely_jphot, 0d0)
      call read_real("MDIJ_MIN",mdijmin, 600d0)
      call read_real("MDIJ_MAX",mdijmax,14000d0)
      ! call read_real("VBFQSQAMIN",qsqAmin,4d0)   ! set in global.inc

c Jet veto cuts
      call read_logical("JVETO",lveto,.false.)
      call read_real("DELY_JVETO",dely_jveto,0d0)
      call read_real("YMAX_VETO",ymax_veto,4.5d0)
      call read_real("PTMIN_VETO",ptmin_veto,10d0)

*      print *,"min m(ZA) (only for(ZAjj)   (M_ZA_MIN): ", mza_min
      
c semileptonic decay cuts
      if (semileptonic) then
        call read_int("SINGLE_DECAYJET",semilep_fatjet,0)
        if (semilep_fatjet.lt.0 .or. semilep_fatjet.gt.2) then
           print*, "ERROR:"
           print*, "For SINGLE_DECAYJET only values 0, 1 or 2 are valid!"
           STOP
        endif
        if (LOplusjet .and. semilep_fatjet.eq.2) then
           print*, "CLASH: SINGLE_DECAYJET = 2 in cuts.dat and LOPROCESS_PLUS_JET = true"
           print*, "in vbfnlo.dat is not allowed to be set at at the same time!"
           STOP
        endif
        call read_int("DEF_TAGJET",def_tagjet,1)
        if (def_tagjet.le.0 .or. def_tagjet.gt.4) then
           print*, "ERROR:"
           print*, "For DEF_TAGJET only values 1 - 4 are valid!"
           STOP
        endif
        if (def_tagjet.eq.2) ysign = .true.
        if(ptmin_tag_1.lt.ptmin_tag_2) then
           print*,"Warning: PTMIN_TAG_1 is lower than PTMIN_TAG_2!"
        endif
        if (DEF_TAGJET.eq.2) then
          call read_real("ETA_CENTRAL",eta_central,2d0)
        else
          eta_central=-1d0
        endif
        call read_real("PTMIN_TAG_1",ptmin_tag_1,20d0)
        call read_real("PTMIN_TAG_2",ptmin_tag_2,20d0)
        call read_logical("HARD_CENTRAL",hard_central,.false.)
        call read_real("PTMIN_CENTRAL",ptmin_central,20d0)
        call read_real("QSQAMIN_ZDEC",qsqAmin_had,100d0)
        call read_int("RECONST_HAD_V",semilep_reconst_v_had,0)
        if (semilep_reconst_v_had.lt.0 .or. semilep_reconst_v_had.gt.3) then
           print*, "ERROR:"
           print*, "For RECONST_HAD_V only values 0, 1, 2 or 3 are valid!"
           STOP
        endif
        call read_real("V_MASS_RANGE",semilep_v_mass_range, 20d0)
        call read_logical("VBFCUTS_ALWAYS",semilepAlwaysVBFCuts,.false.)
      else
        semilepAlwaysVBFCuts = .false.
        semilep_fatjet = 0
        def_tagjet = 1
      endif


      ! minimal number of jets: set in process.F to njet_min
      ! min_jets is used as minimum number of jets within cut routine and defineJets
      ! can be different from njet_min for semileptonic decays
      ! for hadronic z-decays: additional photon cut in phasespace.F
      min_jets=njet_min 
      if (semilep_fatjet.eq.2) min_jets = min_jets - 1

      call closefile

      do i=1,max_jets-1
         if(ptmin_def(i).lt.ptmin_def(i+1)) then
            print*,"Warning: PT_JET_MIN is not ordered; higher values"
            print*,"should come first"
         endif
      enddo
      do i=1,n_lcharged-1
         if(ptlmin(i).lt.ptlmin(i+1)) then
            print*,"Warning: PT_L_MIN is not ordered; higher values"
            print*,"should come first"
         endif
      enddo
      do i=1,n_photon-1
         if(ptgmin(i).lt.ptgmin(i+1)) then
            print*,"Warning: PT_G_MIN is not ordered; higher values"
            print*,"should come first"
         endif
      enddo

c print information on cuts
      if (lglobalprint) then
      print *,""
      print *,"         information on cut parameters"
      print *,"-----------------------------------------------"
      print *,""
      ! print max(1, njet_min) jets for LO, njet_min+1 for NLO
      if (donlo) then
        print "(A,5(G9.4e1))"," minimum jet pt        (PT_JET_MIN):    ",(ptmin_def(i), i=1,njet_min+1)
      else
        print "(A,5(G9.4e1))"," minimum jet pt        (PT_JET_MIN):    ",(ptmin_def(i), i=1,max(1,njet_min))
      endif
      print *,"maximum jet rapidity                 (Y_JET_MAX): ",ymax_def
      print *,"maximum parton pseudorapidity          (Y_P_MAX): ",ymax_part
      print *,"minimum jet-jet R-separation           (RJJ_MIN): ",rjjmin
      if (semilep_fatjet.eq.1) then
        ! min_jets-1 partially allowed
        print *,"minimum number of jets                (NJET_MIN): ",min_jets, " (",min_jets-1,")"
      else
        print *,"minimum number of jets                (NJET_MIN): ",min_jets
      endif
      print *,"single jet from hadronic decay (SINGLE_DECAYJET): ",semilep_fatjet
      print *,"exponent of generalised kT algorithm (PGENKTJET): ",jetgenktp
      print *,"minimum lepton-lepton R-separation     (RLL_MIN): ",rll_min
      print *,"maximum lepton-lepton R-separation     (RLL_MAX): ",rll_max
      print *,"minimum jet-lepton R-separation        (RJL_MIN): ",rjl_min
      print *,"minimum jet-photon R-separation        (RJG_MIN): ",rjg_min
      print *,"minimum lepton-photon R-separation     (RLG_MIN): ",rlg_min   
      print *,"minimum photon-photon R-separation     (RGG_MIN): ",rgg_min
      print *,"maximum photon-photon R-separation     (RGG_MAX): ",rgg_max
      print *,"maximum lepton rapidity                (Y_L_MAX): ",ymax_l
      print "(A,5(G9.4e1))"," minimum lepton pt       (PT_L_MIN):    ",(ptlmin(i), i=1,n_lcharged)
*      print "(A,5(G9.4e1))"," maximum lepton pt       (PT_L_MAX):    ",(ptlmax(i), i=1,n_lcharged)
*      print *,"minimum invisible pt                  (PT_IN_MIN): ",ptinmin
      print *,"minimum invariant dilepton mass        (MLL_MIN): ",mllmin
      print *,"maximum invariant dilepton mass        (MLL_MAX): ",mllmax
      print *,"m_ll cuts only on opposite-sign l   (MLL_OSONLY): ",mllosonly
*      print *,"treat neutrino as lepton:         (NEUTRINO_CUT): ",lneutrinoCut
      print *,"maximum photon rapidity                (Y_G_MAX): ",ymax_g
      print "(A,5(G9.4e1))"," minimum photon pt       (PT_G_MIN):    ",(ptgmin(i), i=1,n_photon)
*      print "(A,5(G9.4e1))"," maximum photon pt       (PT_G_MAX):    ",(ptgmax(i), i=1,n_photon)
      print *,"photon isolation cut                 (PHISOLCUT): ",isol_phopart
      print *,"efficiency of photon isolation cut   (EFISOLCUT): ",eff_phopart
*      print *,"photon helicity                      (PHOTONHEL): ",photon_hel
      print *,"minimum invariant lepton-photon mass   (MLG_MIN): ",mlAmin
      print *,"maximum invariant lepton-photon mass   (MLG_MAX): ",mlAmax
      print *,"minimal missing transverse momentum (PTMISS_MIN): ",ptmiss_min
*      print *,"minimum Wgamma transverse mass          (MT_MIN): ",MTmin
*      print *,"maximum Wgamma transverse mass          (MT_MAX): ",MTmax
      print *,"minimum jet rapidity separation      (ETAJJ_MIN): ",etajjmin
      print *,"tagging jets in opposite hemispheres     (YSIGN):   ",ysign
      print *,"leptons fall inside rapidity gap     (LRAPIDGAP):   ",lrapidgap
      print *,"min leptons y-dist from tagging jets   (DELY_JL): ",dely_jlep
      print *,"photons fall inside rapidity gap     (GRAPIDGAP):   ",grapidgap
      print *,"min photons y-dist from tagging jets   (DELY_JG): ",dely_jphot
      print *,"dijet min mass cut on tag jets        (MDIJ_MIN): ",mdijmin
      print *,"dijet max mass cut on tag jets        (MDIJ_MAX): ",mdijmax
      print *,"veto criteria for jets                   (JVETO):   ",lveto
      print *,"minimum veto-tag y-dist             (DELY_JVETO): ",dely_jveto
      print *,"maximum |y| for veto jet             (YMAX_VETO): ",ymax_veto
      print *,"minimum pT for veto jet             (PTMIN_VETO): ",ptmin_veto
      if (semileptonic) then
         print *,"definition of tagging jets          (DEF_TAGJET): ",def_tagjet
         if (DEF_TAGJET.eq.2) print *,'max |y| for "central region"       (ETA_CENTRAL): ',eta_central
         print *,'minimal pt for harder tagging jet  (PTMIN_TAG_1): ',ptmin_tag_1
         print *,'minimal pt for softer tagging jet  (PTMIN_TAG_2): ',ptmin_tag_2
         if (hard_central) print *,'min pt for jet in central region (PTMIN_CENTRAL): ',ptmin_central
         print *,'mass reconst. cut for V -> had   (RECONST_HAD_V): ',semilep_reconst_v_had
         print *,'mass range for RECONST_HAD_V      (V_MASS_RANGE): ',semilep_v_mass_range
      endif
      print *,""


      if ((qsqAmin.gt.0d0 .and. photoncut_vbf) .or. 
     &    (qsqAmin_had.gt.0d0 .and. zdecay_had .and. procID.ne.Hjj_ZZhad)) then
         print *,""
         print *,"         information on technical cuts"
         print *,"-----------------------------------------------"
         print*," "
      endif
      if (qsqAmin.gt.0d0 .and. photoncut_vbf) then
         print 171," Photon exchange with virtuality Q^2 < ",qsqAmin,
     1   " GeV^2 is eliminated. "
         print*," Recover by calculating  q gamma ---> q g V (V)   etc."
         print*," "
      endif
      if (qsqAmin_had.gt.0d0 .and. zdecay_had .and. procID.ne.Hjj_ZZhad) then
         print*," "
         print 171," Hadronically decaying Z bosons / photons with Q^2 < ",qsqAmin_had,
     1   " GeV^2 are neglected. "
      endif
 171  format(a,f10.1,a)

      endif 

      END

c*************************************************************************
      LOGICAL FUNCTION Cuts(jets, njets, leptons, nl, invisible, ninv, 
     -     photons, nphotons, L)
c*************************************************************************
c     This is the general cut routine. It is intended to be the interface
c     to any basic, specialized or user-defined cut routines. 
c*************************************************************************
c     INPUT 
c     jets      : the array with the jet information as it is returned by the 
c                 jetdefinition routine.
c     njets     : number of visible jets
c     leptons   : the array with the lepton information as it is returned by the 
c                 defineleptons routine.
c     nl        : number of visible leptons
c     invisible : the array with the invisible particle information as it is 
c                  returned by the  defineinvisible routine.
c     ninv      : number of invisible particles
c     photons   : the array with the photon information as it is returned by the 
c                 definephotons routine.
c     nphotons  : number of visible photons
c     L         : id for momentum set in jets, leptons, invisible, photons
c
c     OUTPUT
c     Cuts      : .true. or .false. depending on whether the event passes the cuts.
c*************************************************************************

      IMPLICIT NONE
#include "global.inc"
#include "process.inc"
#include "cuts.inc"

c     input variables
      real*8 partons(0:7,max_jets,max_kin)
      real*8 jets(0:7,max_jets,max_kin) 
      real*8 leptons(0:8,max_v,max_kin), invisible(0:8,max_v,max_kin)
      real*8 photons(0:7,max_v,max_kin)   
      integer npartons(max_kin), ninv(max_kin)
      integer L, njets(max_kin), nl(max_kin), nphotons(max_kin)

      integer jetc

c     external functions
      logical Basic_Cuts, VBF_Cuts, JetVeto_Cuts
      external Basic_Cuts, VBF_Cuts, JetVeto_Cuts

      common/partons/partons,npartons

      ! initialize jetc, tag1, tag2, so it's also set for non-vbf processes
      jetc = 0
      tag1 = 0
      tag2 = 0

c define tagging jets for VBF processes
      if (vbfprocess .or. semilepAlwaysVBFCuts) call DefineTaggingJets(jets(0,1,L), njets(L))


c call of function Basic_Cuts for all processes
      Cuts = Basic_Cuts( partons(0,1,L), npartons(L), jets(0,1,L), 
     1                   njets(L), leptons(0,1,L), nl(L), 
     2                   invisible(0,1,L), ninv(L),
     3                   photons(0,1,L), nphotons(L), decayPartonsCombined(L) )

c call of vbf_cuts only for vbf processes
      if ((vbfprocess .or. semilepAlwaysVBFCuts) .and. Cuts) then
         Cuts = VBF_Cuts( jets(0,1,L), njets(L), leptons(0,1,L),
     1                      nl(L), photons(0,1,L), nphotons(L), jetc)
      endif


c call of jet veto for all processes
c has to be called after VBF_Cuts as this routine returns the ID of the
c additional central jet for the VBF processes
      if (Cuts) then
         Cuts = JetVeto_Cuts( jets(0,1,L), njets(L), leptons(0,1,L),
     1                          nl(L), photons(0,1,L), nphotons(L), jetc)
      endif

      RETURN
      END


c*************************************************************************
      SUBROUTINE DefineTaggingJets(jets, njets)
c*************************************************************************
      implicit none
#include "global.inc"
#include "cuts.inc"
#include "koppln.inc"
#include "process.inc"

c     input variables
      real*8 jets(0:7,max_jets) 
      integer njets

      integer i, j
      integer tempjet

      real*8 mjj2, vectorbosonmass, m_inv, mjj
      integer exclude1, exclude2
      external mjj2

      ! define tagging jets. Usually it's the two hardest jets, 
      ! except for semileptonic decay processes
      if (semileptonic) then
         select case (def_tagjet)
         case (1)   ! largest pt
            tag1 = 1
            tag2 = 2
         case (2)   ! largest pt outside central region, opposite hemisphere
            tag1 = 0
            tag2 = 0
            do i = 1,njets
               if (jets(6,i).gt.eta_central) then
                  tag1 = i
                  exit
               endif
            enddo
            do i = 1,njets
               if (jets(6,i).lt.(-eta_central)) then
                  tag2 = i
                  exit
               endif
            enddo
            if (tag1*tag2.eq.0) return
            if (jets(5,tag2).gt.jets(5,tag1)) then
               tempjet = tag1
               tag1 = tag2
               tag2 = tempjet
            endif
         case (3)   ! largest eta separation, pt_tag1 > pt_tag2
            tag1 = 1
            tag2 = 1
            do i = 2,njets
               if (jets(6,i).gt.jets(6,tag1)) tag1 = i
               if (jets(6,i).lt.jets(6,tag2)) tag2 = i
            enddo
            if (tag1.eq.1 .and. tag2.eq.1) then
               tag1 = 1
               tag2 = 2
            endif
            if (jets(5,tag2).gt.jets(5,tag1)) then
               tempjet = tag1
               tag1 = tag2
               tag2 = tempjet
            endif
         case (4)   ! largest pT, but exclude jet compatible with mV
            if (zdecay_had) then
              vectorbosonmass = xmz
            else
              vectorbosonmass = xmw
            endif
            m_inv = 1d20
            exclude1 = 0
            exclude2 = 0
            ! dijet case
            do i = 1,njets-1
               do j = i+1,njets
                  mjj = sqrt(abs(mjj2(jets(0,i),jets(0,j))))
                  if (abs(mjj-vectorbosonmass).lt.abs(m_inv-vectorbosonmass)) then
                     m_inv = mjj
                     exclude1 = i
                     exclude2 = j
                  endif
               enddo
            enddo
            ! monojet case
            if (semilep_fatjet.gt.0) then
               do i = 1,njets
                  mjj = sqrt(abs(mjj2(jets(0,i),jets(0,i)))/4d0)
                  if (abs(mjj-vectorbosonmass).lt.abs(m_inv-vectorbosonmass)) then
                     m_inv = mjj
                     exclude1 = i
                     exclude2 = 0
                  endif
               enddo
            endif
            tag1 = 0
            tag2 = 0
            do i = 1,njets
               if (i.ne.exclude1 .and. i.ne.exclude2) then
                  tag1 = i
                  exit
               endif
            enddo
            do i = (tag1+1),njets
               if (i.ne.exclude1 .and. i.ne.exclude2) then
                  tag2 = i
                  exit
               endif
            enddo
         end select
      else
         tag1 = 1
         tag2 = 2
      endif

      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION BASIC_CUTS(partons, npartons, jets, njets,
     1                  leptons, nl, invisible, ninv, photons, nphotons,
     2                  fatDecayJet)
c*************************************************************************
c     This routine offers basic cuts for all processes.
c*************************************************************************
      IMPLICIT NONE
#include "global.inc"
#include "cuts.inc"
#include "koppln.inc"
#include "process.inc"

c     input variables
      real*8 partons(0:7,max_jets)
      real*8 jets(0:7,max_jets) 
      real*8 leptons(0:8,max_v)
      real*8 invisible(0:8,max_v)      
      real*8 photons(0:7,max_v)
      integer npartons, njets, nl, nphotons, ninv

      integer i, j,  k, kk, kmax
      real*8 mjj2,mll,mlA,mjj

      real*8 mom(0:3), vectorbosonmass
      logical fatDecayJet

      logical pass_cuts, had_rec, ltemp

      double precision RSEPS
      external RSEPS
c     photon isolation stuff

      double precision delta(npartons)
      double precision dlist(npartons,nphotons)
      double precision ptlist(npartons,nphotons)
      double precision flist(npartons,nphotons)
      double precision spt(npartons,nphotons)
c pt-ordered leptons and photons
      double precision ptlepton(max_v)
      double precision ptphoton(max_v)
      real*8 mass2
      external mass2
  
      logical ldebug
      parameter (ldebug=.false.)
  
      pass_cuts = .true.


      if(ldebug) then
         print*, 'in cuts.F:'
         print*, 'j1', jets(:,1)
         print*, 'j2', jets(:,2)
         print*, 'l1', leptons(:,1)
         print*, 'l2', leptons(:,2)
         print*, 'l3', leptons(:,3)
         print*, 'l4', leptons(:,4)
      endif


c   JET CUTS
c=================

c...require minimum min_jets jets in final state
      if (semilep_fatjet.eq.1 .and. fatDecayJet) then
        pass_cuts = pass_cuts .and. ( njets.ge.(min_jets-1) )
      else
        pass_cuts = pass_cuts .and. ( njets.ge.min_jets )
      endif

c...jet-jet r-separation
      do i=1,njets-1
         do j = i+1, njets
            pass_cuts=pass_cuts.and.
     &      ( RSEPS(jets(6,i),jets(7,i),jets(6,j),jets(7,j)).gt.rjjmin)
         enddo
      enddo


c   LEPTON CUTS
c=================
c...maximum lepton rapidity and minimum lepton pt
      do i=1, n_v
        ptlepton(i) = 0d0
      enddo
      do i=1, nl
        do j=1,i
          if (ptlepton(j).lt.leptons(5,i)) then
            do k=i-1,j,-1
              ptlepton(k+1)=ptlepton(k)
            enddo
            ptlepton(j) = leptons(5,i)
            exit
          endif
        enddo
      enddo
      do i=1, nl
         pass_cuts=pass_cuts.and. 
     &   ( (abs(leptons(6,i)).lt.ymax_l).and.(ptlepton(i).gt.ptlmin(i)))
*     &   .and.(ptlepton(i).le.ptlmax(i)) )  
      enddo
!c apply lepton cuts also to neutrinos for debugging
!      do i=1, n_v
!        ptlepton(i) = 0d0
!      enddo
!      do i=1, ninv
!        do j=1,i
!          if (ptlepton(j).lt.invisible(5,i)) then
!            do k=i-1,j,-1
!              ptlepton(k+1)=ptlepton(k)
!            enddo
!            ptlepton(j) = invisible(5,i)
!            exit
!          endif
!        enddo
!      enddo
!      do i=1, ninv
!         pass_cuts=pass_cuts.and. 
!     &   ( (abs(invisible(6,i)).lt.ymax_l).and.(ptlepton(i).gt.ptlmin(i)))
!*     &   .and.(ptlepton(i).le.ptlmax(i)) )  
!      enddo
!c end

c...cut on minimum pt of invisible particles
c$$$      do i=1, ninv
c$$$         pass_cuts=pass_cuts.and. 
c$$$     &   ( (invisible(5,i).gt.ptinmin) )
c$$$     
c$$$      enddo

c...lepton-lepton r-separation min and max 
      do i=1, nl-1
         do j = i+1, nl
            pass_cuts=pass_cuts.and.
     -           ( RSEPS(leptons(6,i),leptons(7,i),leptons(6,j),
     -           leptons(7,j)).gt.rll_min .and. RSEPS(leptons(6,i),
     -           leptons(7,i),leptons(6,j),leptons(7,j)).lt.rll_max)
         enddo
      enddo

c... cut on invariant dilepton mass
      do i = 1,nl-1
         do j = i+1,nl
            if ( (leptons(8,i)*leptons(8,j).lt.0) .or. !opposite sign leptons only
     &           ( .not. mllosonly) ) then
               mll = mjj2(leptons(0,i),leptons(0,j))
               pass_cuts=pass_cuts.and.(mll.gt.mllmin**2).and.
     -              (mll.lt.mllmax**2)
            endif
         enddo
      enddo


c   PHOTON CUTS
c=================
c...maximum photon rapidity and minimum photon pt
      do i=1, n_v
        ptphoton(i) = 0d0
      enddo
      do i=1, nphotons
        do j=1,i
          if (ptphoton(j).lt.photons(5,i)) then
            do k=i-1,j,-1
              ptphoton(k+1)=ptphoton(k)
            enddo
            ptphoton(j) = photons(5,i)
            exit
          endif
        enddo
      enddo
      do i=1, nphotons
         pass_cuts=pass_cuts.and. 
     &   ( (abs(photons(6,i)).lt.ymax_g).and.(ptphoton(i).gt.ptgmin(i)))
*     &   .and.(ptphoton(i).le.ptgmax(i)) )  
      enddo

c...photon-photon r-separation min and max
      do i=1, nphotons-1
         do j = i+1, nphotons
            pass_cuts=pass_cuts.and.
     -           ( RSEPS(photons(6,i),photons(7,i),photons(6,j),
     -           photons(7,j)).gt.rgg_min .and. RSEPS(photons(6,i),
     -           photons(7,i),photons(6,j),photons(7,j)).lt.rgg_max)
         enddo
      enddo

c...transverse mass from hep-ph/9305314

c     transverse mass only useful for Wgamma production!!!
c$$$      if(nphotons.eq.1)then
c$$$               MT2 = (sqrt(mjj2(leptons(0,1),photons(0,1))+ 
c$$$     &                (photons(1,1)+leptons(1,1))**2 +(photons(2,1)+leptons(2,1))**2)+
c$$$     &                 invisible(5,1))**2     
c$$$     &                -(photons(1,1)+leptons(1,1)+invisible(1,1))**2
c$$$     &               -(photons(2,1)+leptons(2,1)+invisible(2,1))**2
c$$$!               print*,MT2, mla, MTmin, MTmax 
c$$$               pass_cuts=pass_cuts.and.(MT2.gt.MTmin**2).and.(MT2.le.MTmax**2)
c$$$      endif

!     c...NEW photon isolation cut according to Frixione,PLB 429(1998),369

      ltemp=.true.
!       print*
!       print*
!       print*
!       print*,'cuts.F'
!       print*, "nphotons",nphotons, "npartons",npartons
       do i=1, nphotons
         k=0
         do j=1, npartons
            delta(j)=RSEPS(photons(6,i),photons(7,i),partons(6,j),
     -           partons(7,j))
!             print*,"j",j, npartons, delta(j)
!             print*,"partons(6,j)",partons(6,j)
!             print*,"partons(7,j)",partons(7,j)
!             print*, "parton", partons(:,j)
            if (delta(j).lt.isol_phopart) then
               k=k+1
               dlist(k,i)=delta(j)
               ptlist(k,i)=partons(5,j)
               flist(k,i)=(1.-dcos(dlist(k,i)))/(1.-dcos(isol_phopart))
            endif
         enddo
         kmax=k
         do k=1,kmax
            spt(k,i)=0d0
            do kk=1,kmax
               if (dlist(kk,i).le.dlist(k,i)) then
                  spt(k,i)=spt(k,i)+ptlist(kk,i)
               endif
            enddo
!             ltemp=ltemp
            pass_cuts=pass_cuts
     &           .and.(spt(k,i).le.eff_phopart*photons(5,i)*flist(k,i))

         enddo
       enddo
!          print*, "pass_cuts in cuts", ltemp,eff_phopart,isol_phopart
!          print*


c   ADDITIONAL CUTS
c=================================
c...jet-lepton r-separation
      do i=1, nl
         do j=1, njets
            pass_cuts=pass_cuts.and.
     &           ( RSEPS(leptons(6,i),leptons(7,i),jets(6,j),
     -           jets(7,j)).gt.rjl_min )
         enddo
      enddo

c...jet-photon r-separation
      do i=1, nphotons
         do j=1, njets
            pass_cuts=pass_cuts.and.
     &           ( RSEPS(photons(6,i),photons(7,i),jets(6,j),
     -           jets(7,j)).gt.rjg_min )
         enddo
      enddo

c...lepton-photon r-separation
      do i=1, nphotons
         do j=1, nl
            pass_cuts=pass_cuts.and.
     &           ( RSEPS(photons(6,i),photons(7,i),leptons(6,j),
     -           leptons(7,j)).gt.rlg_min )
         enddo
      enddo

c...minimal invariant lepton+photon mass for all combinations of charged leptons and photons
      do i = 1,nl
         do j = 1,nphotons
               mlA = mjj2(leptons(0,i),photons(0,j))
               pass_cuts=pass_cuts.and.(mlA.gt.mlAmin**2).and.
     &             (mlA.le.mlAmax**2)
         enddo
      enddo

c...minimal missing transverse momentum
      do j = 1,2  ! only transverse elements relevant
         mom(j) = 0d0
         do i = 1,nl
            mom(j) = mom(j) + leptons(j,i)
         enddo
         do i = 1,nphotons
            mom(j) = mom(j) + photons(j,i)
         enddo
         do i = 1,njets
            mom(j) = mom(j) + jets(j,i)
         enddo
      enddo
c Alternative using transverse momenta of invisible particles (neutrinos, stable BSM particles like LSP)
c      do j = 1,2  ! only transverse elements relevant
c         mom(j) = 0d0
c         do i = 1,ninv
c            mom(j) = mom(j) + invisible(j,i)
c         enddo
c      enddo
      pass_cuts = pass_cuts .and. (sqrt(mom(1)**2+mom(2)**2).ge.ptmiss_min)


c...reconstruct mass of hadronically decaying jet

      if (semilep_reconst_v_had.gt.0) then

        had_rec=.false.
        if (zdecay_had) then
          vectorbosonmass = xmz
        else
          vectorbosonmass = xmw
        endif

c... dijet case
        if (semilep_reconst_v_had.eq.1 .or. semilep_reconst_v_had.eq.3) then
           do i = 1,njets-1
              do j = i+1,njets
                 if(i.ne.tag1 .and. i.ne.tag2 .and. j.ne.tag1 .and. j.ne.tag2) then !do not consider tagging jets
                    mjj = sqrt(abs(mjj2(jets(0,i),jets(0,j))))
                    if (abs(mjj-vectorbosonmass).lt.semilep_v_mass_range) had_rec=.true.
                 endif
              enddo
           enddo
        endif

c... monojet case
        if (semilep_reconst_v_had.eq.2 .or. semilep_reconst_v_had.eq.3) then
           do i = 1,njets
              if(i.ne.tag1 .and. i.ne.tag2) then !do not consider tagging jets
                 mjj = sqrt(abs(mjj2(jets(0,i),jets(0,i)))/4d0)
                 if (abs(mjj-vectorbosonmass).lt.semilep_v_mass_range) had_rec=.true.
              endif
           enddo
        endif

        pass_cuts=pass_cuts.and.had_rec

      endif

cfc      select case(procID)
cfc      case(ZAjj,ZAjjjLO)
cfc! m_za_min
cfc      temp(:3) = leptons(:3,1)+leptons(:3,2)+photons(:3,1)
cfc      temp(4) = mass2(temp(0))
cfc      temp(4)=sqrt(temp(4))
cfc      pass_cuts = pass_cuts .and. temp(4).ge.mza_min
cfc      case(ZAjj_n,ZAjjjLO_n)
cfc! m_za_min
cfc      temp(:3) = invisible(:3,1)+invisible(:3,2)+photons(:3,1)
cfc      temp(4) = mass2(temp(0))
cfc      temp(4)=sqrt(temp(4))
cfc      pass_cuts = pass_cuts .and. temp(4).ge.mza_min
cfcc      if (.not.temp(4).ge.mza_min) print*, "pass_cuts2", temp(4)
cfc      end select


      BASIC_CUTS = pass_cuts
      
      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION VBF_Cuts(jets, njets, leptons, nl, photons, 
     -     nphotons, jetc)
c*************************************************************************
c     This routine offers vbf specific cuts, that is cuts for processes 
c     with 2 or more jets.
c     jetc returns the number of the additional hard central jet if requested, else 0
c*************************************************************************
      IMPLICIT NONE
#include "global.inc"
#include "cuts.inc"

c     input variables
      real*8 jets(0:7,max_jets) 
      integer njets, nl, nphotons
      real*8 leptons(0:8,max_v)
      real*8 photons(0:7,max_v)

      integer i
      real*8 ylepmax,ylepmin,m2,mjj2
      real*8 yphotmax, yphotmin

      integer jetc

      logical pass_cuts

      double precision RSEPS
      external RSEPS

      logical Basic_Cuts
      external Basic_cuts

      
      pass_cuts = .true.


c check if two tagging jets have been found
      if (tag1*tag2.eq.0) then
         VBF_cuts = .false.
         return
      endif

c...jet-jet rapidity separation
      if (njets.ge.2) then
        pass_cuts=pass_cuts.and.
     &        (abs(jets(6,tag1)-jets(6,tag2)).gt.etajjmin)
      endif

c...jets #1 and #2 must have opposite sign rapidity
      if (ysign) then
         if (njets.ge.2) then
            pass_cuts=pass_cuts.and. ( (jets(6,tag1)*jets(6,tag2)).le.0d0 )
         endif
      endif

c...leptons fall inside rapidity gap
      if(lrapidgap) then
         if (njets.ge.2) then
            ylepmax = max(jets(6,tag1),jets(6,tag2)) - dely_jlep
            ylepmin = min(jets(6,tag1),jets(6,tag2)) + dely_jlep
            do i=1,nl
               pass_cuts=pass_cuts.and.(leptons(6,i).lt.ylepmax.and.
     &              leptons(6,i).gt.ylepmin)      
            enddo
         endif
      endif

c...photons fall inside rapidity gap
      if(grapidgap) then
         if (njets.ge.2) then
            yphotmax = max(jets(6,tag1),jets(6,tag2)) - dely_jphot
            yphotmin = min(jets(6,tag1),jets(6,tag2)) + dely_jphot
            do i=1,nphotons
               pass_cuts=pass_cuts.and.(photons(6,i).lt.yphotmax.and.
     &              photons(6,i).gt.yphotmin)      
            enddo
         endif
      endif

c...  tagging jets invariant dijet mass cut
      if (njets.ge.2) then
         m2=mjj2(jets(0,tag1),jets(0,tag2))
         pass_cuts=pass_cuts.and.(m2.ge.mdijmin**2.and.m2.le.mdijmax**2)
      endif

c...extra vbf cuts for semileptonic processes
      if (semileptonic) then

c...require hard central jet
         jetc=0
         if (hard_central) then
            if (DEF_TAGJET.eq.2) then
               do i = 1,njets
                  if ( abs(jets(6,i)).le.eta_central .and. jets(5,i).ge.ptmin_central
     &                .and. i.ne.tag1 .and. i.ne.tag2 ) then
                     jetc = i
                     exit
                  endif
               enddo
            else
               do i = 1,njets
                  if ( jets(6,i).le.max(jets(6,tag1),jets(6,tag2))
     &                .and. jets(6,i).ge.min(jets(6,tag1),jets(6,tag2))
     &                .and. jets(5,i).ge.ptmin_central
     &                .and. i.ne.tag1 .and. i.ne.tag2 ) then
                     jetc = i
                     exit
                  endif
               enddo
            endif
            pass_cuts=pass_cuts.and.(jetc.gt.0)
         endif

c...minimal pt for tagging jets
         pass_cuts=pass_cuts.and.(jets(5,tag1).gt.ptmin_tag_1)
         pass_cuts=pass_cuts.and.(jets(5,tag2).gt.ptmin_tag_2)

      endif

      VBF_cuts = pass_cuts
      
      RETURN
      END


c*************************************************************************
      LOGICAL FUNCTION JetVeto_Cuts(jets, njets, leptons, nl, photons, 
     -     nphotons, jetc)
c*************************************************************************
c     This routine offers jet veto cuts, that is cuts on additional jets
c     apart from the standard number of the LO process.
c*************************************************************************
      IMPLICIT NONE
#include "global.inc"
#include "process.inc"
#include "cuts.inc"

c     input variables
      integer jetc
      real*8 jets(0:7,max_jets) 
      integer njets, nl, nphotons
      real*8 leptons(0:8,max_v)
      real*8 photons(0:7,max_v)

      integer i, nveto, maxnjet
      real*8 yvetomax,yvetomin


      logical pass_cuts
      
      pass_cuts = .true.

C... veto jet cut      
      nveto=0
      if(lveto) then
         ! for vbf processes veto additional central jets
         if (vbfprocess) then

            yvetomax = max(jets(6,tag1),jets(6,tag2)) - dely_jveto
            yvetomin = min(jets(6,tag1),jets(6,tag2)) + dely_jveto

            if(njets.gt.min_jets)then
               do i=1,njets
                  if(jets(6,i).le.yvetomax.and.
     &               jets(6,i).ge.yvetomin.and.
     &               abs(jets(6,i)).le.ymax_veto.and.
     &               jets(5,i).gt.ptmin_veto.and.
     &               i.ne.tag1 .and. i.ne.tag2 .and. i.ne.jetc) 
     &                  nveto=nveto+1            
               enddo
               if (hard_central) nveto = nveto + 1
               nveto = nveto - min_jets + 2
            endif

         ! for all other processes just veto the additional jet
         else

            maxnjet = min_jets

            if(njets.gt.min_jets)then
               do i=maxnjet+1,njets
                  if(abs(jets(6,i)).le.ymax_veto.and.
     &               jets(5,i).gt.ptmin_veto) 
     &                  nveto=nveto+1            
               enddo
            endif

         endif

         pass_cuts=pass_cuts.and.(nveto.le.0)

      endif
      
      JetVeto_cuts = pass_cuts
      
      RETURN
      END


c*************************************************************************
      SUBROUTINE DefineLeptons(v, Lmax, leptons, n_leptons)
c*************************************************************************
c     Similar to "DefineJets" in jetdefinition.f this routine extracts
c     the leptons out of "v" into "leptons". Only e, mu, tau particles are
c     taken. In addition pt, y and phi information is calculated and stored.
c*************************************************************************
      implicit none
#include "global.inc"
#include "cuts.inc"
      real*8 v(0:3,max_v,max_kin)
      real*8 leptons(0:8,max_v,max_kin)
      integer n_leptons(max_kin), Lmax
      real*8 pt,y,phi
      integer i, l, mu

      do l=1, Lmax

         do i=1, max_v
            do mu=0,8
               leptons(mu,i,l) = 0d0
            enddo
         enddo

         n_leptons(l)=0
c     get leptons out of decay_particles     
         do i=1, n_v
            if ( (abs(particle_IDs(i)).eq.11) ! check if it is a lepton
     &           .or.(abs(particle_IDs(i)).eq.13)
     &           .or.(abs(particle_IDS(i)).eq.15) ! taus!!
cfc     &       .or.lneutrinoCut .and.(
cfc     &               (abs(particle_IDs(i)).eq.12) ! check if it is a lepton
cfc     &           .or.(abs(particle_IDs(i)).eq.14)
cfc     &           .or.(abs(particle_IDS(i)).eq.16))
     &           ) then
               n_leptons(l)=n_leptons(l)+1
               do mu = 0,3      ! transfer to "leptons"
                  leptons(mu,n_leptons(l),l) = v(mu,i,l)
               enddo
               leptons(8,n_leptons(l),l) = particle_IDs(i)
            endif 
         enddo
c     calculate rapidity, pt, phi-angle for the leptons
         do i=1, n_leptons(l)
            call lego( leptons(0,i,l), pt, y, phi )
            leptons(4,i,l) = 0.0 ! mass  
            leptons(5,i,l) = pt  ! transverse momentum  
            leptons(6,i,l) = y   ! rapidity 
            leptons(7,i,l) = phi ! azimuthal angle
         end do
         
      enddo

      RETURN
      END


c*************************************************************************
      SUBROUTINE DefinePhotons(v, Lmax, photons, n_photons)
c*************************************************************************
c     Similar to "DefineJets" in jetdefinition.f this routine extracts
c     the photons out of "v" into "photons". In addition pt, y and phi 
c     information is calculated and stored.
c*************************************************************************
      implicit none
#include "global.inc"
      real*8 v(0:3,max_v,max_kin)
      real*8 photons(0:7,max_v,max_kin)
      integer n_photons(max_kin), Lmax
      real*8 pt,y,phi
      integer i, l, mu

      do l=1, Lmax

         do i=1, max_v
            do mu=0,7
               photons(mu,i,l) = 0d0
            enddo
         enddo

         n_photons(l)=0
c     get photons out of decay_particles     
         do i=1, n_v
            if ( (abs(particle_IDs(i)).eq.22) ! check if it is a photon
     &           ) then
               n_photons(l)=n_photons(l)+1
               do mu = 0,3      ! transfer to "photons"
                  photons(mu,n_photons(l),l) = v(mu,i,l)
               enddo
            endif 
         enddo
c     calculate rapidity, pt, phi-angle for the leptons
         do i=1, n_photons(l)
            call lego(photons(0,i,l), pt, y, phi )
            photons(4,i,l) = 0.0 ! mass  
            photons(5,i,l) = pt ! transverse momentum  
            photons(6,i,l) = y  ! rapidity 
            photons(7,i,l) = phi ! azimuthal angle
         end do
         
      enddo

      RETURN
      END
     
c*************************************************************************
      SUBROUTINE DefineInvisible(v, Lmax, invisible, n_invisible)
c*************************************************************************
c     Similar to "DefineJets" in jetdefinition.f this routine extracts
c     the invisible particles out of "v" into "invisible". In addition pt,
c     y and phi information is calculated and stored.
c*************************************************************************
      implicit none
#include "global.inc"
#include "cuts.inc"
      real*8 v(0:3,max_v,max_kin)
      real*8 invisible(0:8,max_v,max_kin)
      integer n_invisible(max_kin), Lmax
      real*8 pt,y,phi
      integer i, l, mu

      do l=1, Lmax
         do i=1, max_v
            do mu=0,8
               invisible(mu,i,l) = 0d0
            enddo
         enddo

         n_invisible(l)=0
c     get invisible particles out of decay_particles     
         do i=1, n_v
            if ( (abs(particle_IDs(i)).eq.12) ! check if it is a neutrino
     &           .or.(abs(particle_IDs(i)).eq.14)
     &           .or.(abs(particle_IDS(i)).eq.16)
cfc   if ( (.not.lneutrinoCut) .and.((abs(particle_IDs(i)).eq.12) ! check if it is a neutrino
cfc     &           .or.(abs(particle_IDs(i)).eq.14)
cfc     &           .or.(abs(particle_IDS(i)).eq.16))
     &           ) then
               n_invisible(l)=n_invisible(l)+1
               do mu = 0,3      ! transfer to "invisible"
                  invisible(mu,n_invisible(l),l) = v(mu,i,l)
               enddo
               invisible(8,n_invisible(l),l) = particle_IDs(i)
            endif 
         enddo
	 

c     calculate rapidity, pt, phi-angle for the invisible particles
         do i=1, n_invisible(l)
            call lego( invisible(0,i,l), pt, y, phi )
            invisible(4,i,l) = 0.0 ! mass  
            invisible(5,i,l) = pt  ! transverse momentum  
            invisible(6,i,l) = y   ! rapidity 
            invisible(7,i,l) = phi ! azimuthal angle
         end do
      enddo    
      

      RETURN
      END
 
