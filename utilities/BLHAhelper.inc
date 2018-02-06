!*************************************************************************
!
!     This file contains variable declarations for the
!     Binoth Les Houches Accords
!     (arXiv:1001.1307 and arXiv:1308.3462)
!
!*************************************************************************
!
!     Author: Michael Rauch <michael.rauch@kit.edu>, Sep 2013
!     Last changed: Mar 2015 by Michael Rauch
!
!*************************************************************************

! variables which store set entries - initialize in BLHA_initialize

       double precision blha_wwidth, blha_zwidth, blha_hwidth
       double precision blha_sin2w, blha_vev, blha_gf
       double precision blha_xmw, blha_xmz, blha_xmh
       double precision blha_alpha, blha_ewfac, blha_alphas
       double precision blha_CA, blha_CF
       double precision blha_gammaQuark, blha_gammaGluon
       double precision blha_tgammaQuark, blha_tgammaGluon
       double precision blha_KQuark, blha_KGluon
       integer blha_Nc
       integer blha_ewrenormscheme
       logical blha_recomp, blha_DR, blha_ranhelsum
       logical blha_couplingsoff, blha_anomcoupl

       common/BLHA_parameters/ blha_wwidth, blha_zwidth, blha_hwidth
       common/BLHA_parameters/ blha_sin2w, blha_vev, blha_gf
       common/BLHA_parameters/ blha_xmw, blha_xmz, blha_xmh
       common/BLHA_parameters/ blha_alpha, blha_ewfac, blha_alphas
       common/BLHA_parameters/ blha_CA, blha_CF
       common/BLHA_parameters/ blha_gammaQuark, blha_gammaGluon
       common/BLHA_parameters/ blha_tgammaQuark, blha_tgammaGluon
       common/BLHA_parameters/ blha_KQuark, blha_KGluon
       common/BLHA_parameters/ blha_Nc
       common/BLHA_parameters/ blha_ewrenormscheme
       common/BLHA_parameters/ blha_recomp, blha_DR, blha_ranhelsum
       common/BLHA_parameters/ blha_couplingsoff, blha_anomcoupl


       integer max_blhaproc, max_subproc
       parameter (max_blhaproc = 10000) ! TODO, use allocate to reduce size
       parameter (max_subproc = 6)

       integer blha_numproc                                           ! number of processes
       integer blha_numsubproc(max_blhaproc)                          ! number of subprocesses in process
       integer blha_numptcl(max_blhaproc)                             ! number of particles in process
       integer blha_numparton(max_blhaproc)                           ! number of partons in process
       integer blha_numelweak(max_blhaproc)                           ! number of elweak in process
       integer blha_amptype(max_blhaproc)                             ! type of amplitude (tree,loop,cctree,sctree)
       integer blha_procsubproc(max_subproc,max_blhaproc)             ! process ID for subprocess
       logical blha_lojsubproc(max_subproc,max_blhaproc)              ! subprocess is LO+j one
       integer blha_idsubproc(max_subproc,max_blhaproc)               ! internal sub-id for subprocess
       integer blha_multsubproc(max_subproc,max_blhaproc)             ! internal sub-id for subprocess
       integer blha_particlemap(max_p+max_v,max_subproc,max_blhaproc) ! momenta mapping MatchBox to phys
       integer blha_invmap(max_p+max_v,max_subproc,max_blhaproc)      ! momenta mapping phys to MatchBox
       integer blha_physdiagmap(max_p+max_v,max_subproc,max_blhaproc) ! momenta mapping phys to diag
       integer blha_diagphysmap(max_p+max_v,max_subproc,max_blhaproc) ! momenta mapping diag to phys
       integer blha_fsign(max_p+max_v,max_subproc,max_blhaproc)       ! fsign mapping
       integer blha_alphasorder(max_blhaproc)                         ! order in alpha_s
       integer blha_alphaorder(max_blhaproc)                          ! order in alpha_ew

       common/BLHA_procmap/ blha_numproc, blha_numsubproc
       common/BLHA_procmap/ blha_numptcl, blha_numparton, blha_numelweak
       common/BLHA_procmap/ blha_amptype
       common/BLHA_procmap/ blha_procsubproc, blha_idsubproc
       common/BLHA_procmap/ blha_multsubproc
       common/BLHA_procmap/ blha_particlemap, blha_invmap
       common/BLHA_procmap/ blha_physdiagmap, blha_diagphysmap
       common/BLHA_procmap/ blha_fsign
       common/BLHA_procmap/ blha_alphasorder, blha_alphaorder
       common/BLHA_procmap/ blha_lojsubproc


       double precision blha_amp(2*(max_p+max_v)**2), blha_scale ! TODO: reduce array size
       double precision blha_helrand
       integer blha_ranhelcomb
       integer blha_curproc, blha_cursubproc, blha_lastprocId
       logical blha_lastLOplusJet

       common/BLHA_thisproc/ blha_amp, blha_scale
       common/BLHA_thisproc/ blha_helrand, blha_ranhelcomb
       common/BLHA_thisproc/ blha_curproc, blha_cursubproc, blha_lastprocID
       common/BLHA_thisproc/ blha_lastLOplusJet

       integer blha_pssubproc(max_blhaproc)                           ! last used subproc for phasespace
       integer blha_pssubstep(max_blhaproc)                           ! steps with different
       integer blha_numbosons(max_blhaproc)                           ! number of fermion pairs
       integer blha_bosons(3,max_subproc,max_blhaproc)                ! boson type of fermion pairs

       common/BLHA_phasespace/ blha_pssubproc, blha_pssubstep, blha_numbosons, blha_bosons

       integer blha_lastPhaseSpace

       common/BLHA_thisPhaseSpace/ blha_lastPhaseSpace

       double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT
       double precision BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU
       double precision BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
       COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT
       COMMON /BRANCH/ BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU
       COMMON /BRANCH/ BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

