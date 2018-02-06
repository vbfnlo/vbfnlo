!***************************************************************************
!  This is an include-file of the VBFNLO main program.
!  It contains declarations of global variables.
!***************************************************************************
!  Author(s) : G. Klaemke
!***************************************************************************
!  Last modified: 15.09.2006
!***************************************************************************

! the process identification number
      integer procID

! switch to enable (default) or disable the calculation of toploops
! so far only used in WHj@NLO
      logical dotoploops
      parameter (dotoploops = .true.)

! enables built-in ckm values for processes supporting CKM
! so far only W, Wj, WH, WHj
      logical luseckm
      parameter (luseckm = .false.)

! the number of iterations of the integration
      integer N_iterations1, N_iterations2
! the number of points/shots in the last iteration = 2**N_points
! an array is used because there are several loops (LO, virtual contrib., real emissions)
! -1       : real emissions loop
!  0       : LO loop
!  1,2,... : the various virtual contributions
      integer N_points(-1:5)

! the center of mass energy of the collider
      double precision ecm

! center of mass energy^2 in electroweak processes used for form factors and kmatrix unitarisation
      double precision ss

! phasespace dimensionality - this is process dependend!!

      integer*4 PS_dimension

! length of the monaco random number array. (Must always be larger than PS_dimension.)

      integer max_PS_dim
      parameter (max_PS_dim = 30)

! the maximum number of 4-momenta stored in 'partons'
! here, the partons are the incoming quarks/gluons + the 2 tagging jets
! + the real emission (for NLO only)

      integer max_p
      parameter (max_p=6)

! the maximum number of additional particles like a higgs boson, or any decay products.
! additional/decay particles means: everything which is not a parton involved in
! the basic QCD 2->2(+1)+X process.

      integer max_v
      parameter (max_v=10)

! number of internal x parameters
      integer nx
      parameter( nx=2 )

! the maximum number of jets
      integer max_jets
      parameter (max_jets=max_p+max_v-nx)

! max number of momentum configurations for real emission processes
      integer max_kin
      parameter (max_kin=25)

! max number of different phase spaces
      integer max_ps
      parameter (max_ps = 5)

! max number of NLO contributions
      integer max_nlo
      parameter (max_nlo = 11)

! maximal number of resonances that can be put into one random number
! (see MultiResonances in phasespace/ps_tools.F)
      integer max_res, max_multires
      parameter (max_res=20, max_multires=10)

! actual number of entries in partons
      integer n_p

! actual number of entries in decay_particles
      integer n_v

! actual number of kinematics
      integer n_kin, n_qed
      
      
      integer photon_hel
      common/ZAjjhel/ photon_hel


!     The particle identification code according to the Monte Carlo particle numbering scheme.
!     reminder :  e- =  11, nue    =  12, mu- =  13, numu    =  14, tau- =  15, nutau  =  16,
!                 e+ = -11, nuebar = -12, mu+ = -13, numubar = -14, tau+ = -15, nutaub = -16,
!             photon = 22,      Z0 = 23,   W+ = 24,       h0 = 25,    A0 = 36,      H+ = 37
!     in addition, we identify a general parton (quark/gluon) with id=93
      integer particle_IDs(max_v)
! file names for input and output of monaco grid information
      character*50 gridname1(1:max_ps), gridname2(1:max_ps), gridname3(1:max_ps), gridname4(1:max_ps)

!     proton (+1) / antiproton (-1) beams
      integer sign1, sign2

!     is NLO implemented for the process?
      logical hasNLO, hasFLoops

!     switches for born/virtual contribution/real emission
      logical doBorn, doVirtuals, doEmissions, doFLoops, inFLoops


!* Flags, switches, parameters for running in SUSY and with ew
!* corrections
! switches:
!  loopSqr_switch = whether |M_{loop}|^{2} is calculated
!  ewcor_switch = whether ew corns are included
!  FeynH_switch =  whether FeynHiggs is used
!  slha_switch = whether a SLHA file is used
      logical loopSqr_switch, ewcor_switch, FeynH_switch, slha_switch
      logical loopSqr_switchIN
! flags:
!  model = SM (1) and MSSM (2)
!  Higgs type  = SM (0), h0 (1), HH (2), A0 (3)
!  ewapprox = which approximations the user has selected for the ew corns
!             (0 = no approximations)
!  sector = which particles are included in the ew corns
!  ewscheme = method by which LO em coupling is calculated
      integer model, HiggsType, ewapprox, sector, ewscheme, floops, MHloops


!     number of loops for different virtual contributions
      integer NLO_Loops

!     number of loops for different phase space generators
      integer PS_Loops, Loops_sub_LO,Loops_sub_NLO,sub_number
!      commom /subprocnu/sub_number

! For VBF Higgs, flag determining manner in which the Higgs is distributed
!   0 = on-shell; 1 = Breit-Wigner; 2 = Passarino
      integer higgsscheme

! put everything in global common blocks: split into integer, logical etc.
      common /cglobalD/ ecm
      common /cglobalI/ ewapprox, model, HiggsType, ewscheme, floops
      common /cglobalI/ sector, MHloops, procID, higgsscheme
      common /cglobalI/ N_iterations1, N_iterations2, N_points, n_p, n_v
      common /cglobalI/ sign1, sign2, NLO_Loops, PS_dimension, particle_IDs
      common /cglobalI/ PS_Loops, sub_number,Loops_sub_LO,Loops_sub_NLO
      common /cglobalL/ FeynH_switch, slha_switch, ewcor_switch
      common /cglobalL/ loopSqr_switch, loopSqr_switchIN
      common /cglobalL/ hasNLO, hasFLoops
      common /cglobalL/ doBorn, doVirtuals, doEmissions, doFLoops
      common /cglobalL/ inFLoops
      common /cglobalC/ gridname1, gridname2, gridname3, gridname4


      common /cglobal2/ n_kin, n_qed

      common /kmatrix/ ss

! distribute finite terms for Born*virtual+real emission = c*|Born|**2 + ....
! with c = -9/2 - pi^2/3 = creal + cvirt
      double precision creal, cvirt, pi, pi2o3p9h
      parameter (pi = 3.141592653589793d0, pi2o3p9h = pi**2/3d0+9d0/2d0)
      parameter (cvirt = pi**2/3d0-7d0, creal = -pi2o3p9h-cvirt) ! for carlo

! and for pp -> H jjj.
      double precision crealH3j, cvirtH3j
      parameter (cvirtH3j = pi**2/3d0-8d0, crealH3j = 2d0*pi**2/3d0-13d0/2d0)

! and for pp -> VV (V)
      double precision crealVV, cvirtVV, pi2o3m4
      parameter (pi2o3m4=2d0*(pi**2/3d0-4d0))
      parameter (cvirtVV = 4/3d0*pi**2-8d0,crealVV = pi2o3m4-cvirtVV)

! fix number of external quark flavors to sum over; for CC processes both
! up- and down-type quarks of the highest generation need to be allowed
! (i.e. the largest even number less than or equal to nfl is used)
! For vbf processes that involve both CC and NC channels, by switching vbfNFLb
!  to .true. and nfl=5 external bottom quarks will be included in the NCs but
!  not in the CCs
! For other processes nfl=5 allows external bottom quarks if the process has
!  no CC components
! Set nfl=4 to have no external bottom quarks anywhere
      integer nfl, nflVBF
      logical vbfNFLb
      common /vbfNFL/ nfl, nflVBF, vbfnflb

! switch to enable external b bbar contributions
      logical lincludebb
      parameter (lincludebb=.true.)
! switch to enable external b bbar contributions only if they are final
!   and recombined into one jet
!   this is ignored if lincludebb is true
      logical lincludebbfinal
      parameter (lincludebbfinal=.true.)
! gets set by cuts routine to signal that final state particles
! are combined (there is only one jet)
      logical lincludebbfinaldo


! flag to signal whether two quarks in the final state are recombined into a jet
      logical lqqrecombined

      common /cglobalL2/ lqqrecombined, lincludebbfinaldo

! qsqAmin is the minimum virtuality for t-channel photon exchange.
! Subamplitudes will be set to zero if photon virtuality falls below
! qsqAmin
      double precision qsqAmin
      parameter (qsqAmin=4d0)
!      common /photondivVBF/qsqAmin


! qsqAmin_had is the minimum virtuality for a hadronically decaying photon
! in the semileptonic processes.
! Amplitude will be set to zero if Q^2_gamma < qsqAmin_had
      double precision qsqAmin_had
      common /photondiv/qsqAmin_had
!      parameter (qsqAmin_had=100d0)


! flag for warnings in matrixelement routines
      logical lwarn
      parameter (lwarn = .false.)

      integer igauge
      parameter (igauge=2)
! igauge allows to select different ways to insure EW gauge invariance
! in the face of finite W and Z width effects for EW W and Z production:
!
!   igauge = 1       overall factor scheme
!            2       complex mass scheme, m**2-->m**2-i*m*Gamma
!                    in all propagators
!            3       for ibos=2 only: naive Z Breit Wigner
!           -1       width = 0, for testing purposes only

      logical with_anom
      common /anom_switch/ with_anom
      logical with_anomHiggs
      common /anomHiggs/ with_anomHiggs
      logical with_kk
      common /kk_switch/ with_kk
      logical with_spin2
      common /spin2_switch/ with_spin2


      integer max_sub
      parameter (max_sub=20)

! instabilities for pentagons
      integer pentcountm, pentcount2m,pentcountps(max_ps),pentcount2ps(max_ps)
      integer pentcountsub(max_sub),pentcount2sub(max_sub),pentcountm_Q,pentcountsub_Q(max_sub)
      common /instpent/ pentcountm, pentcount2m, pentcountps, pentcount2ps,pentcountsub,pentcount2sub
      common /instpent/ pentcountm_Q,pentcountsub_Q
      save /instpent/


! instabilities for hexagons
      integer hexcountm, hexcount2m,hexcountps(max_ps),hexcount2ps(max_ps)
      integer hexcountsub(max_sub),hexcount2sub(max_sub),hexcountm_Q,hexcountsub_Q(max_sub)
      common /insthex/ hexcountm, hexcount2m, hexcountps, hexcount2ps,hexcountsub,hexcount2sub
      common /insthex/ hexcountm_Q, hexcountsub_Q
      save /insthex/


! instabilities for boxes
      integer boxcountm, boxcount2m,boxcountps(max_ps),boxcount2ps(max_ps)
      integer boxcountsub(max_sub),boxcount2sub(max_sub),boxcountm_Q,boxcountsub_q(max_sub)
      common /instbox/ boxcountm, boxcount2m, boxcountps, boxcount2ps,boxcountsub,boxcount2sub
      common /instbox/ boxcountm_Q,boxcountsub_q
      save /instbox/


! instabilities for boxes
      integer box2countm, box2count2m,box2countm_Q,box2countps(max_ps),box2count2ps(max_ps)
      integer box2countsub(max_sub),box2count2sub(max_sub),box2countsub_q(max_sub)
      common /instbox2/ box2countm, box2count2m, box2countps, box2count2ps,box2countsub,box2count2sub
      common /instbox2/ box2countm_q,box2countsub_q
      save /instbox2/

! instabilities per subprocess
      real*8 instab_error
      common /cinstab_error/ instab_error

! data files
      logical lreaddata, lwritedata
      common /rwdata/ lreaddata, lwritedata

! store contributions split up in subproc, phase space and gnlo
! gnlo = -1: LO
! gnlo = NLO_Loops+1: RE
! sum up contributions in index 0
      integer final_Ncall
      real*8 final_xsec, final_sdev2,final_time
      common /final_sub/ final_xsec(-1:max_sub,-1:max_ps,-1:max_nlo),final_sdev2(-1:max_sub,-1:max_ps,-1:max_nlo)
      common /final_sub/ final_time(-1:max_sub,-1:max_ps,-1:max_nlo)
      common /final_sub_Int / final_Ncall(-1:max_sub,-1:max_ps,-1:max_nlo)

! default values from GramDets (GGFLO)
      real*8 GRAMCTO0,GRAMDTO0

      common/GRAMTO/GRAMCTO0,GRAMDTO0

! output file for cross section data
      character*50 XSECfile
      common/xsecfile/ XSECfile

! hidden debug switches
      logical born_off, virtuals_off, emissions_off,reduceSubStat
      common /nlo_switch / born_off, virtuals_off, emissions_off,reduceSubStat

! flag to output timing information
      logical showtiming
      common/showtimingcom/ showtiming
