!*************************************************************************  
!
!     This file contains all functions defined in the
!     Binoth Les Houches Accords 
!     (arXiv:1001.1307 and arXiv:1308.3462)
!
!*************************************************************************  
!
!     Author: Michael Rauch <michael.rauch@kit.edu>
!     Initial version: Sep 2013
!     Last changed: Jun 2015 by Michael Rauch
!
!     TODO: split into useful subparts
!
!*************************************************************************  
!   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
!*************************************************************************
!     SUBROUTINE OLP_SetParameter_VBFNLO(para, re, im, ierr)
!     SUBROUTINE OLP_GetParameter_VBFNLO(para, re, im, ierr)
!     SUBROUTINE OLP_EvalSubProcess_VBFNLO(bproc, pp, mu, alphas, rval)
!     SUBROUTINE OLP_EvalSubProcess2_VBFNLO(bproc, pp, mu, rval, acc)
!     SUBROUTINE OLP_PhaseSpacePoint(proc, r, pp, weight)
!     SUBROUTINE OLP_Polvec(p, q, eps)
!     SUBROUTINE OLP_Info(version, message)
!     SUBROUTINE VBFNLO_BLHA2Amp(pp)
!     LOGICAL FUNCTION VBFNLO_ClosestOnshell(v,proc,subproc)
!     SUBROUTINE BLHA_initialize()
!     SUBROUTINE BLHA_start()
!     SUBROUTINE BLHA_dorecomp()
!     SUBROUTINE VBFNLO_SetupProcess(nparticles,pdgprocess,orderAlphas,orderAlpha,amptype,procok)
!     SUBROUTINE BLHA_setDR(DR)
!     SUBROUTINE BLHA_setNf(Nf)
!     SUBROUTINE BLHA_setNc(Nc)
!     SUBROUTINE BLHA_setCouplingsOff()
!     SUBROUTINE BLHA_setEWRenormScheme(scheme)
!     SUBROUTINE BLHA_cctree(i,j,res)
!     SUBROUTINE BLHA_sctree(i,j,res)
!     SUBROUTINE BLHA_error(errmsg,filename,lineno)
!     SUBROUTINE BLHA_amptypeerror(amptype,filename,lineno)
!*************************************************************************  

      recursive SUBROUTINE BLHA_initialize()
          use globalvars, only: lglobalprint, ldoscales, ldoblha, seed
          use readinput, only: inputpath, pdfpath
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"
      integer pdflib
      COMMON/PDFparameters/pdflib

! definitions for second Higgs
      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU, &
                       BH2GAM, BH2GAMZ, XMH2, XGH2,  &
                       sinba, cosba, ch2ww, ch2zz
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU, &
                       BH2GAM, BH2GAMZ, XMH2, XGH2,  &
                       sinba, cosba, ch2ww, ch2zz

      logical blha_runinit
      data blha_runinit /.true./

      double precision alphas5_hardwired
      external alphas5_hardwired
 
      if (blha_runinit) then
        blha_runinit = .false.
        blha_numproc = 0
        blha_lastprocID = 0
        blha_Nc = 3
        nfl = 5
        blha_pssubproc = 1

        ldoscales = .false.
        lglobalprint = .false.
        ldoblha = .true.

        with_kk=.false.
        with_anomHiggs = .false.
        with_anom = .false.
        with_spin2 = .false.
        model = 0
        ewcor_switch = .false.
        higgsmix = 0
        higgsscheme = 0
        pdflib = 0
        inputpath=""
        pdfpath=""
        sub_number = 1

        blha_wwidth = -1
        blha_zwidth = -1
        blha_hwidth = -1
        blha_xmw    = -1
        blha_xmz    = -1
        blha_xmh    = -1
        blha_sin2w  = -1
        blha_vev    = -1
        blha_gf     = -1
        blha_alpha  = -1
        blha_ewfac  = -1
        blha_alphas = -1

        seed = 7782
        call InitRandomNumbers()

! default values of parameters
        xmt = 172.4d0
        xmb = 4.855d0
        xmc = 1.65d0
        xmtau = 1.77684d0
        gf = 1.16637d-5
        xmw = 80.398d0
        xmz = 91.1876d0
        sin2w = 1-xmw**2/xmz**2
        xmh = 126.d0
        alfa = -1
        alfas = alphas5_hardwired(xmz**2,0)
! default masses for anomalous Higgs couplings
        Mf(1,1) = 0d0
        Mf(1,2) = 0d0
        Mf(1,3) = 0d0
        Mf(2,1) = 0.00051099891D0
        Mf(2,2) = 0.105658367D0
        Mf(2,3) = xmtau
        Mf(3,1) = 0.0024D0
        Mf(3,2) = xmc
        Mf(3,3) = xmt
        Mf(4,1) = 0.00475D0
        Mf(4,2) = 0.104D0
        Mf(4,3) = xmb
! second Higgs parameters
        sinba = 1d0
        cosba = 0d0

        blha_ranhelsum = .false.
        blha_helrand = -1
        blha_couplingsoff = .false.
        blha_ewrenormscheme = 0
        blha_anomcoupl = .false.

        blha_recomp = .true.
        call BLHA_dorecomp()

! just a default -- change via OLP_SetParameter_VBFNLO
        procID = Hjj

        call proc_assignment()
        call InitProcess()

        call InitPhaseSpace()
        blha_lastPhaseSpace = procId

        call InitGaugeTest()
      endif

      end

!*************************************************************************  

      SUBROUTINE BLHA_start()
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      blha_numproc = 0
      call BLHA_setDR(0)
      call BLHA_setNf(5)
      call BLHA_setNc(3)

      return
      end

!*************************************************************************  

      SUBROUTINE OLP_SetParameter_VBFNLO(para, re, im, ierr)
          use readinput, only: usedefaults
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      character(*) para
      double precision re, im
      integer ierr

      call BLHA_initialize()

!cc mass
      if ((trim(para).eq."mass(1)") .or. &
          (trim(para).eq."mass(2)") .or. &
          (trim(para).eq."mass(3)") .or. &
          (trim(para).eq."mass(11)") .or. &
          (trim(para).eq."mass(12)") .or. &
          (trim(para).eq."mass(13)") .or. &
          (trim(para).eq."mass(14)") .or. &
          (trim(para).eq."mass(15)") .or. &
          (trim(para).eq."mass(21)") .or. &
          (trim(para).eq."mass(22)")) then
! must be massless
        if (re.ne.0d0) then
          ierr = 0
          write(*,*) "VBFNLO error in OLP_SetParameter: Parameter", &
            trim(para)
          write(*,*) "  Particle must be massless."
        else
          ierr = 1
        endif
        return
      elseif (trim(para).eq."mass(4)"  .or. &
              trim(para).eq."charm_mass") then ! charm
        ierr = 1
        if (re .eq. xmc) return
        xmc = re
      elseif (trim(para).eq."mass(5)" .or. &
              trim(para).eq."bottom_mass") then ! bottom
        ierr = 1
        if (re .eq. xmb) return
        xmb = re
      elseif (trim(para).eq."mass(6)" .or. &
              trim(para).eq."top_mass") then ! top
        ierr = 1
        if (re .eq. xmt) return
        xmt = re
      elseif (trim(para).eq."mass(16)" .or. &
              trim(para).eq."tau_mass") then ! tau
        ierr = 1
        if (re .eq. xmtau) return
        xmtau = re
      elseif (trim(para).eq."mass(23)" .or. &
              trim(para).eq."Z_mass") then ! Z
        ierr = 1
        if (re .eq. blha_xmz) return
        blha_xmz = re
      elseif (trim(para).eq."mass(24)" .or. &
              trim(para).eq."W_mass") then ! W
        ierr = 1
        if (re .eq. blha_xmw) return
        blha_xmw = re
      elseif (trim(para).eq."mass(25)" .or. &
              trim(para).eq."H_mass") then ! H
        ierr = 1
        if (re .eq. blha_xmh) return
        blha_xmh = re
      elseif (para(1:min(len_trim(para),5)).eq."mass(") then
! unknown particle
        ierr = 2
        return
!cc width
      elseif ((trim(para).eq."width(1)") .or. &
              (trim(para).eq."width(2)") .or. &
              (trim(para).eq."width(3)") .or. &
              (trim(para).eq."width(11)") .or. &
              (trim(para).eq."width(12)") .or. &
              (trim(para).eq."width(13)") .or. &
              (trim(para).eq."width(14)") .or. &
              (trim(para).eq."width(15)") .or. &
              (trim(para).eq."width(21)") .or. &
              (trim(para).eq."width(22)")) then
! must have zero width
        if (re.ne.0d0) then
          ierr = 0
          write(*,*) "VBFNLO error in OLP_SetParameter: Parameter",trim(para)
          write(*,*) "  Particle is massless and must have zero width."
        else
          ierr = 1
        endif
        return
      elseif (trim(para).eq."width(4)" .or. &
              trim(para).eq."charm_width") then ! charm
        ierr = 2
        return
      elseif (trim(para).eq."width(5)" .or. &
              trim(para).eq."bottom_width") then ! bottom
        ierr = 2
        return
      elseif (trim(para).eq."width(6)" .or. &
              trim(para).eq."top_width") then ! top
        ierr = 2
        return
      elseif (trim(para).eq."width(16)" .or. &
              trim(para).eq."tau_width") then ! tau
        ierr = 2
        return
      elseif (trim(para).eq."width(23)" .or. &
              trim(para).eq."Z_width") then ! Z
        ierr = 1
        if (re .eq. blha_zwidth) return
        blha_zwidth = re
      elseif (trim(para).eq."width(24)" .or. &
              trim(para).eq."W_width") then ! W
        ierr = 1
        if (re .eq. blha_wwidth) return
        blha_wwidth = re
      elseif (trim(para).eq."width(25)" .or. &
              trim(para).eq."H_width") then ! H
        ierr = 1
        if (re .eq. blha_hwidth) return
        blha_hwidth = re
      elseif ((para(1:min(len_trim(para),6)).eq."width(")) then
! unknown particle
        ierr = 2
        return
!cc CKM matrix 
      elseif ((trim(para).eq."VV12") .or.  &
              (trim(para).eq."VV13") .or.  &
              (trim(para).eq."VV21") .or.  &
              (trim(para).eq."VV23") .or.  &
              (trim(para).eq."VV31") .or.  &
              (trim(para).eq."VV32")) then  
        if ((re.ne.0d0) .or. (im.ne.0d0)) then
          ierr = 0
          write(*,*) "VBFNLO error in OLP_SetParameter: Parameter", &
            trim(para)
          write(*,*) "  Non-diagonal CKM matrix not supported."
        else
          ierr = 1
        endif
        return
      elseif ((trim(para).eq."VV11") .or.  &
              (trim(para).eq."VV22") .or.  &
              (trim(para).eq."VV33")) then  
        if ((re.ne.1d0) .or. (im.ne.0d0)) then
          ierr = 0
          write(*,*) "VBFNLO error in OLP_SetParameter: Parameter",trim(para)
          write(*,*) "  Non-diagonal CKM matrix not supported."
        else
          ierr = 1
        endif
        return
!cc EW parameters
      elseif (trim(para).eq."sw2" .or. &
              trim(para).eq."sin_th_2") then
        ierr = 1
        if (re .eq. blha_sin2w) return
        blha_sin2w = re
      elseif (trim(para).eq."vev") then
        ierr = 1
        if (re .eq. blha_vev) return
        blha_vev = re
      elseif (trim(para).eq."Gf") then
        ierr = 1
        if (re .eq. blha_gf) return
        blha_gf = re
      elseif (trim(para).eq."alpha") then
        ierr = 1
        if (re .eq. blha_alpha) return
        blha_alpha = re
      elseif (trim(para).eq."ewfactor") then
! must be a positive number
        if (re.le.0d0) then
          ierr = -2
          return
        endif
        ierr = 1
        if (re .eq. blha_ewfac) return
        blha_ewfac = re
      elseif (trim(para).eq."alphas") then
        ierr = 1
        if (re .eq. blha_alphas) return
        blha_alphas = re
!cc Process ID
      elseif (trim(para).eq."process") then
        process = re
        procID = re
        ierr = 1
        call InitProcess()
        return
!cc Phase space dimension (# of random numbers required)
!cc -- not a settable parameter
      elseif (trim(para).eq."PSdimension") then
        ierr = -2
        return
!cc centre-of-mass energy
      elseif (trim(para).eq."sqrtS") then
        ierr = 1
        if (re .eq. ecm) return
        ecm = re
        blha_lastPhaseSpace = 0
!cc random helicity summation
      elseif (trim(para).eq."ranhelsum") then
        if (re .eq. 0) then
          blha_ranhelsum = .false.
        else
          blha_ranhelsum = .true.
        endif
        ierr = 1
        return
      elseif (trim(para).eq."HelicityRN") then
        blha_helrand = re
        ierr = 1
        return
!cc anomalous couplings
      elseif (trim(para).eq."anomcoupl") then
        if (re .eq. 0) then
          blha_anomcoupl = .false.
          blha_lastprocID = 0
          ierr = 1
        else
          blha_anomcoupl = .true.
          blha_lastprocID = 0
          ierr = 1
          ! check for existence if file
          call InitAnomCouplings(blha_anomcoupl)
          if (usedefaults) ierr = 0
        endif
        return
!cc number of colours
      elseif (trim(para).eq."Nc") then
        call BLHA_setNc(int(re))
        ierr = 1
        return
!cc number of light flavours
      elseif (trim(para).eq."Nf") then
        call BLHA_setNf(int(re))
        ierr = 1
        return
!cc Unknown parameter
      else
        ierr = -2
        return
      endif

!cc parameter changed -- needs recomputation
      blha_recomp = .true.

      end

!*************************************************************************  

      SUBROUTINE OLP_GetParameter_VBFNLO(para, re, im, ierr)
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      DOUBLE PRECISION CLR,XM2,XMG,Bk,Vk,Ak
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),Bk(6,6,6), &
           Vk(4,5),Ak(4,5)

      character(*) para
      double precision re, im
      integer ierr

! use default values
      re = 0d0
      im = 0d0

! make sure parameters are up-to-date
      call BLHA_dorecomp()

!cc mass
      if ((trim(para).eq."mass(1)") .or. &
          (trim(para).eq."mass(2)") .or. &
          (trim(para).eq."mass(3)") .or. &
          (trim(para).eq."mass(11)") .or. &
          (trim(para).eq."mass(12)") .or. &
          (trim(para).eq."mass(13)") .or. &
          (trim(para).eq."mass(14)") .or. &
          (trim(para).eq."mass(15)") .or. &
          (trim(para).eq."mass(21)") .or. &
          (trim(para).eq."mass(22)")) then
! must be massless
        re = 0d0
        ierr = 1
      elseif (trim(para).eq."mass(4)"  .or. &
              trim(para).eq."charm_mass") then ! charm
        re = xmc 
        ierr = 1
      elseif (trim(para).eq."mass(5)" .or. &
              trim(para).eq."bottom_mass") then ! bottom
        re = xmb
        ierr = 1
      elseif (trim(para).eq."mass(6)" .or. &
              trim(para).eq."top_mass") then ! top
        re = xmt
        ierr = 1
      elseif (trim(para).eq."mass(16)" .or. &
              trim(para).eq."tau_mass") then ! tau
        re = xmtau
        ierr = 1
      elseif (trim(para).eq."mass(23)" .or. &
              trim(para).eq."Z_mass") then ! Z
        re = xmz 
        ierr = 1
      elseif (trim(para).eq."mass(24)" .or. &
              trim(para).eq."W_mass") then ! W
        re = xmw
        ierr = 1
      elseif (trim(para).eq."mass(25)" .or. &
              trim(para).eq."H_mass") then ! H
        re = xmh
        ierr = 1
      elseif (para(1:min(len_trim(para),5)).eq."mass(") then
! unknown particle
        ierr = 2
!cc width
      elseif ((trim(para).eq."width(1)") .or. &
              (trim(para).eq."width(2)") .or. &
              (trim(para).eq."width(3)") .or. &
              (trim(para).eq."width(11)") .or. &
              (trim(para).eq."width(12)") .or. &
              (trim(para).eq."width(13)") .or. &
              (trim(para).eq."width(14)") .or. &
              (trim(para).eq."width(15)") .or. &
              (trim(para).eq."width(21)") .or. &
              (trim(para).eq."width(22)")) then
! must have zero width
        re = 0d0
        ierr = 1
      elseif (trim(para).eq."width(4)" .or. &
              trim(para).eq."charm_width") then ! charm
        ierr = 2
      elseif (trim(para).eq."width(5)" .or. &
              trim(para).eq."bottom_width") then ! bottom
        ierr = 2
      elseif (trim(para).eq."width(6)" .or. &
              trim(para).eq."top_width") then ! top
        ierr = 2
      elseif (trim(para).eq."width(16)" .or. &
              trim(para).eq."tau_width") then ! tau
        ierr = 2
      elseif (trim(para).eq."width(23)" .or. &
              trim(para).eq."Z_width") then ! Z
        re = xmg(2)/sqrt(xm2(2))
        ierr = 1
      elseif (trim(para).eq."width(24)" .or. &
              trim(para).eq."W_width") then ! W
        re = xmg(3)/sqrt(xm2(3))
        ierr = 1
      elseif (trim(para).eq."width(25)" .or. &
              trim(para).eq."H_width") then ! H
        re = xmg(6)/sqrt(xm2(6))
        ierr = 1
      elseif ((para(1:min(len_trim(para),6)).eq."width(")) then
! unknown particle
        ierr = 2
!cc CKM matrix 
      elseif ((trim(para).eq."VV12") .or.  &
              (trim(para).eq."VV13") .or.  &
              (trim(para).eq."VV21") .or.  &
              (trim(para).eq."VV23") .or.  &
              (trim(para).eq."VV31") .or.  &
              (trim(para).eq."VV32")) then  
        re = 0d0
        ierr = 1
      elseif ((trim(para).eq."VV11") .or.  &
              (trim(para).eq."VV22") .or.  &
              (trim(para).eq."VV33")) then  
        re = 1d0
        ierr = 1
!cc EW parameters
      elseif (trim(para).eq."sw2" .or. &
              trim(para).eq."sin_th_2") then
        re = sin2w
        ierr = 1
      elseif (trim(para).eq."vev") then
        re = 1d0/sqrt(sqrt(2d0)*gf) 
        ierr = 1
      elseif (trim(para).eq."Gf") then
        re = gf
        ierr = 1
      elseif (trim(para).eq."alpha") then
        re = alfa 
        ierr = 1
      elseif (trim(para).eq."ewfactor") then
        re = blha_ewfac 
        ierr = 1
      elseif (trim(para).eq."alphas") then
        re = alfas 
        ierr = 1
!cc Process ID
      elseif (trim(para).eq."process") then
        re = process
        ierr = 1
!cc Phase space dimension (# of random numbers required)
      elseif (trim(para).eq."PSdimension") then
        re = PS_dimension
        ierr = 1
!cc centre-of-mass energy
      elseif (trim(para).eq."sqrtS") then
        re = ecm 
        ierr = 1
!cc random helicity summation
      elseif (trim(para).eq."ranhelsum") then
        if (blha_ranhelsum) then
          re = 1
        else
          re = 0
        endif
        ierr = 1
      elseif (trim(para).eq."HelicityRN") then
        re = abs(blha_helrand) 
        ierr = 1
!cc anomalous couplings
      elseif (trim(para).eq."anomcoupl") then
        if (blha_anomcoupl) then
          re = 1
        else
          re = 0
        endif
        ierr = 1
!cc number of colours
      elseif (trim(para).eq."Nc") then
        re = blha_Nc
        ierr = 1
!cc number of light flavours
      elseif (trim(para).eq."Nf") then
        re = nfl
        ierr = 1
!cc Unknown parameter
      else
        ierr = 2
      endif

      end

!*************************************************************************  

      SUBROUTINE BLHA_dorecomp()
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/scales.inc"
#include "VBFNLO/utilities/mssm.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision e, g2, s, c, z, w, q, g

! nothing to do
      if (.not.blha_recomp) return

! is already initialized?
      call BLHA_initialize()

! now -- here comes the tricky part
! let's see what has been set as parameters and what we can make out of it

! copy all input parameters
      if (blha_xmw .gt. 0d0) xmw = blha_xmw
      if (blha_xmz .gt. 0d0) xmz = blha_xmz
      if (blha_xmh .gt. 0d0) xmh = blha_xmh
      if (blha_sin2w .gt. 0d0) sin2w = blha_sin2w
      if (blha_gf .gt. 0d0) GF = blha_gf
      if (blha_alphas .gt. 0d0) alfas = blha_alphas

! MW - MZ - sw2
      if ((blha_xmw .gt. 0d0) .and. (blha_xmz .gt. 0d0) .and.  &
          (blha_sin2w .ge. 0d0)) then
! everything ok
      elseif ((blha_xmw .gt. 0d0) .and. (blha_xmz .gt. 0d0)) then
        sin2w =  1 - xmw**2/xmz**2
      elseif ((blha_xmw .gt. 0d0) .and. (blha_sin2w .gt. 0d0)) then
        xmz =  xmw/sqrt(1-sin2w)
      elseif ((blha_xmz .gt. 0d0) .and. (blha_sin2w .gt. 0d0)) then
        xmw =  xmz*sqrt(1-sin2w)
      elseif (blha_xmw .gt. 0d0) then !+default xmz
        sin2w =  1 - xmw**2/xmz**2
      elseif (blha_xmz .gt. 0d0) then !+default sin2w
        xmw =  xmz*sqrt(1-sin2w)
      elseif (blha_sin2w .gt. 0d0) then !+default xmz
        xmw =  xmz*sqrt(1-sin2w)
!      else ! happy with defaults
      endif
      
! GF - alpha
      if ( (blha_alpha .gt. 0d0) .and. &
           ( (blha_gf .le. 0d0) .or. &
             (blha_ewrenormscheme .eq. 1) .or. &
             (blha_ewrenormscheme .eq. 2) ) &
         ) then
        gf = PI*blha_alpha/sqrt(2d0)/xmw**2/sin2w
! GF - vev
      else if ((blha_vev .gt. 0d0) .and. (blha_gf .le. 0d0)) then
        gf =  1/(sqrt(2d0)*blha_vev**2)
      endif

      als(1,1) = alfas
      als(2,1) = alfas
      als(3,1) = alfas

      G = SQRT(4*PI*ALFAS)
      S = SQRT(SIN2W)
      C = SQRT(1.d0 - SIN2W)
      G2 = SQRT(8.d0*GF/SQRT(2.d0))*XMZ*C
      if (blha_ewfac.gt.0d0) G2 = G2*sqrt(blha_ewfac)
      Z = G2/4.d0/C
      W = G2/SQRT(8.d0)
      Q = G2*C
      if ( (blha_alpha .gt. 0d0) .and. .not. (blha_ewrenormscheme .eq. 3) ) then
        E = sqrt(4.d0*PI*blha_alpha)
        if (blha_ewfac.gt.0d0) E = E*sqrt(blha_ewfac)
      else
        E = G2*S
      endif
      ALFA = E**2/(4.d0*PI)

! W,Z,H decays
! BRs are not actually used, but set to avoid recalculation in koppln
      call clearwidths()
      if (blha_wwidth .gt. 0d0) then
        xgw = blha_wwidth
        BWNE = 10.86d-2
        BWUD = 33.71d-2
        BWTB = 0d0
      endif
      if (blha_zwidth .gt. 0d0) then
        xgz = blha_zwidth
        BZNN =  6.67d-2
        BZEE =  3.37d-2
        BZUU = 11.60d-2
        BZDD = 15.60d-2
        BZTT =  0d0
      endif
      if (blha_hwidth .gt. 0d0) then
        xgh = blha_hwidth
        BHWW   = 2.19d-1
        BHZZ   = 2.72d-2
        BHGG   = 8.54d-2
        BHTT   = 0d0
        BHBB   = 5.72d-1
        BHCC   = 2.89d-2
        BHTAU  = 6.27d-2
        BHMU   = 2.18d-4
        BHGAM  = 2.28d-3
        BHGAMZ = 1.56d-3
      endif

! masses
      if (xmtau .gt. 0d0) Mf(2,3) = xmtau
      if (xmc .gt. 0d0) Mf(3,2) = xmc
      if (xmt .gt. 0d0) Mf(3,3) = xmt
      if (xmb .gt. 0d0) Mf(4,3) = xmb

      call setparams(g2)
      call InitAnomCouplings(blha_anomcoupl)
      call koppln(2,e,g2,s,c,z,w,q,g)
      call ctrans(xmb)
      call coupl_haddecay()

      blha_lastprocID = 0

      blha_recomp = .false.

      end

!*************************************************************************  

      SUBROUTINE VBFNLO_SetupProcess(nparticles,pdgprocess,orderAlphas, &
         orderAlpha, amptype, procok)
!*************************************************************************
!     Map particles to processes
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer nparticles
      integer pdgprocess(nparticles)
      integer orderAlphas, orderAlpha, amptype, procok
! local variables
      integer i,j
      integer nparton, pdgparton(nparticles), mapparton(nparticles)
      integer nelweak, pdgelweak(nparticles), mapelweak(nparticles)
! external
      logical isparton
      external isparton

      if (blha_numproc.ge.max_blhaproc) then
        procok = -2
        return
      endif
      blha_numproc = blha_numproc+1
      blha_numsubproc(blha_numproc) = 0

      if (.not.isparton(pdgprocess(1)) .or.  &
          .not.isparton(pdgprocess(2))) then
        procok = 0
        return
      endif

      nparton = 2
      pdgparton(1) = pdgprocess(1)
      mapparton(1) = 1
      pdgparton(2) = pdgprocess(2)
      mapparton(2) = 2
      nelweak = 0

      do i=3,nparticles
        if (isparton(pdgprocess(i))) then
          nparton = nparton+1
          pdgparton(nparton) = pdgprocess(i)
          mapparton(nparton) = i
        else
          nelweak = nelweak+1
          pdgelweak(nelweak) = pdgprocess(i)
          mapelweak(nelweak) = i
        endif
      enddo

! Figure out which process we are in
      if ((nparton .eq. 4) .and. (orderAlphas .eq. 0)) then
! VBF
        call MomMapping_VBF4(pdgparton,mapparton)
      else if ((nparton .eq. 5) .and. (orderAlphas .eq. 1)) then
! VBF+j
        call MomMapping_VBF5(pdgparton,mapparton)
      else if ((nparton .eq. 6) .and. (orderAlphas .eq. 2)) then
! VBF+2j
        call MomMapping_VBF6(pdgparton,mapparton)
      else if ((nparton .eq. 2) .and. (orderAlphas .eq. 0)) then
! QCD+0j
        call MomMapping_QCD(pdgparton, mapparton, nparton)
      else if ((nparton .eq. 3) .and. (orderAlphas .eq. 1)) then 
! QCD+1j
        call MomMapping_QCD(pdgparton, mapparton, nparton)
      else if ((nparton .eq. 4) .and. (orderAlphas .eq. 2)) then
! QCD+2j
        call MomMapping_QCD4(pdgparton,mapparton)
      else if ((nparton .eq. 5) .and. (orderAlphas .eq. 3)) then
! QCD+3j
        call MomMapping_QCD5(pdgparton,mapparton)
      else
! not yet supported
        procok = 0
        return
      endif
      call MomMapping_EW(nparton,nelweak,pdgelweak,mapelweak, &
                         orderAlphas,orderAlpha)
! create diagphysmap
      do i=1,nparton
        do j=1,blha_numsubproc(blha_numproc)
          blha_diagphysmap(blha_physdiagmap(i,j,blha_numproc),j,blha_numproc) = i
        enddo
      enddo
! create invmap
      do i=1,nparton+nelweak
        do j=1,blha_numsubproc(blha_numproc)
          blha_invmap(blha_particlemap(i,j,blha_numproc),j,blha_numproc) = i
        enddo
      enddo

      blha_numptcl(blha_numproc)   = nparticles
      blha_numparton(blha_numproc) = nparton
      blha_numelweak(blha_numproc) = nelweak
      blha_amptype(blha_numproc)   = abs(amptype)

      if (blha_numsubproc(blha_numproc) .lt. 0) then
        procok = 0
      else
        procok = blha_numproc
      endif
      return
      end

!*************************************************************************  

      SUBROUTINE OLP_EvalSubProcess_VBFNLO(bproc, pp, mu, alphas, rval)
!*************************************************************************
!     Evaluate amplitude
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer bproc
      double precision pp(0:4,blha_numptcl(bproc))
      double precision mu
      double precision alphas
      double precision rval(*)

! local variables
      double precision acc

      blha_alphas = alphas
      blha_recomp = .true.

      call OLP_EvalSubProcess2_VBFNLO(bproc, pp, mu, rval, acc)

      return
      end

!*************************************************************************  

      SUBROUTINE OLP_EvalSubProcess2_VBFNLO(bproc, pp, mu, rval, acc)
!*************************************************************************
!     Evaluate amplitude
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer bproc
      double precision pp(0:4,blha_numptcl(bproc))
      double precision mu
      double precision rval(*)
      double precision acc

! local variables
      integer i,j,nrval,numhelsum
      logical linfnan

! external
      double precision RandomNumber
      external RandomNumber

      blha_curproc = bproc
      blha_scale = mu

      if (blha_amptype(blha_curproc).eq.2) then !cctree
        nrval = blha_numptcl(blha_curproc)* &
                (blha_numptcl(blha_curproc)-1)/2
      else if (blha_amptype(blha_curproc).eq.3) then !sctree
        nrval = 2*blha_numptcl(blha_curproc)**2
      else
        nrval = 4
      endif

      do i=1,nrval
        rval(i) = 0d0
      enddo
      acc     = 0d0
 
! make sure parameters are up-to-date
      call BLHA_dorecomp()

      if (blha_ranhelsum) then
! helicity random number not set from MC -> generate
        if (blha_helrand .lt. 0) then
          blha_helrand = RandomNumber()
        endif
      endif

      do i=1,blha_numsubproc(blha_curproc)
        blha_cursubproc=i
        procId=blha_procsubproc(blha_cursubproc,blha_curproc)
        LOplusJet=blha_lojsubproc(blha_cursubproc,blha_curproc)
        if ((procId .ne. blha_lastprocId) .or.  &
            (LOplusJet .neqv. blha_lastLOplusJet)) then
          call proc_assignment()
          call InitProcess()
          blha_lastprocId = procId
          blha_lastLOplusJet = LOplusJet
        endif

        do j=1,nrval
          blha_amp(j) = 0d0
        enddo

        call VBFNLO_BLHA2Amp(pp)

        if (blha_ranhelsum) then
          numhelsum = blha_multsubproc(blha_cursubproc,blha_curproc)
        else 
          numhelsum = 1
        endif
        do j=1,nrval
          rval(j) = rval(j) + blha_amp(j)*numhelsum
        enddo

      enddo

      if (blha_ranhelsum) then
! invalidate as used
        blha_helrand = -blha_helrand
      endif

! remove couplings if requested
      if (blha_couplingsoff) then
        do j=1,nrval
          rval(j) = rval(j) &
                    /(blha_alphas)**blha_alphasorder(blha_curproc) &
                    /(blha_alpha)**blha_alphaorder(blha_curproc)
        enddo
        if (blha_amptype(blha_curproc).eq.1) then !loop
          do j=1,3 ! epsilon^{-2,-1,0} coefficients have one alphas more
            rval(j) = rval(j)/(blha_alphas)
          enddo
        endif
      endif

! check for Inf and NaN
      linfnan = .false.
      do j=1,nrval
        if (rval(j) .ne. rval(j)) linfnan = .true.    ! NaN
        if (rval(j)+1 .eq. rval(j)) linfnan = .true.  ! Inf
      enddo
      if (linfnan) rval(1:nrval) = 0

      return
      end     

!*************************************************************************  

      SUBROUTINE OLP_PhaseSpacePoint(proc, rpsnum, r, pp, weight)
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer proc
      double precision rpsnum
      double precision r(*)
      double precision pp(0:3,*)
      double precision weight

! local variables
      integer subproc
      double precision p(0:3,max_p,max_kin)
      double precision x(2)
      double precision v(0:3,max_v,max_kin)
      integer i,mu,ps_number

! external
      logical Choose_PS, VBFNLO_ClosestOnshell
      external Choose_PS, VBFNLO_ClosestOnshell


! do something only if subprocesses exist
      if (blha_numsubproc(proc) .lt. 1) then
        weight = 0
        return
      endif

! cycle through subprocesses
      blha_pssubproc(proc) = &
        mod( blha_pssubproc(proc)+blha_pssubstep(proc)-1,  &
             blha_numsubproc(proc) ) + 1
      subproc = blha_pssubproc(proc)

! re-initialize process if necessary
      procId=blha_procsubproc(subproc,proc)
      LOplusJet=blha_lojsubproc(subproc,proc)
      if ((procId .ne. blha_lastprocId) .or.  &
          (LOplusJet .neqv. blha_lastLOplusJet)) then
        call proc_assignment()
        call InitProcess()
        blha_lastprocId = procId
        blha_lastLOplusJet = LOplusJet
      endif
! re-initialize phase space if necessary
      if ((blha_procsubproc(subproc,proc) .ne.  &
             blha_lastPhaseSpace)) then
! make sure masses, etc. are up-to-date
        call blha_dorecomp()
        call InitPhaseSpace()
        blha_lastPhaseSpace = blha_procsubproc(subproc,proc)
      endif

! get phase-space point
      weight = 1
      ps_number = 1+int(rpsnum*PS_Loops)
      call phasespace(r, p, x, v, ps_number, weight)
      weight = weight*PS_Loops
      if (.not. Choose_PS(ps_number, v, 1)) weight = 0d0

! check phase-space selection -- skip if only one ew subproc
      if ( (weight .gt. 0d0) .and.                                 &
           (blha_numsubproc(proc) .gt. blha_pssubstep(proc)) .and. &
           .not. VBFNLO_ClosestOnshell(v,proc,subproc))            &
        weight = 0d0
      weight = weight * blha_numsubproc(proc)/blha_pssubstep(proc)

! remove flux factor and conversion into fb
      if (weight .gt. 0d0) then
        weight = weight*2d0*x(1)*x(2)*ecm**2
        weight = weight/3.89379304d11
      endif

! map partons back
      do i=1,blha_numparton(proc)
        do mu=0,3
          pp(mu,blha_particlemap(i,subproc,proc)) &
            = p(mu,i,1)
        enddo
      enddo
! map elweak back
      do i=1,blha_numelweak(proc)
        do mu=0,3
          pp(mu,blha_particlemap(i+blha_numparton(proc), &
               subproc,proc)) &
            = v(mu,i,1) 
        enddo
      enddo

      return
      end

!*************************************************************************  

      SUBROUTINE OLP_Polvec(p, q, eps)
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision p(0:3)
      double precision pbar(0:4)
      double precision q(0:3)
      double precision eps(0:7)

! local variables
      integer i
      double complex ceps(0:3)

      pbar(0:3) = p(0:3)
      pbar(4) = 0d0

! epsilon(p,+) 
      call helvec(pbar,-1,+1,ceps)

      do i=0,3
        eps(2*i)   = dreal(ceps(i))
        eps(2*i+1) = dimag(ceps(i))
      enddo

      return
      end

!*************************************************************************  

      SUBROUTINE OLP_Info(version, message)
        use VBFNLOVersion, only: &
          vbfnloversionstring, setVersion, vbfnloreference
!*************************************************************************
!     Pass parameters
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      character*15 version
      character*255 message

      call setVersion
      version = vbfnloversionstring
      message =  &
        trim(vbfnloreference(1))//achar(10)// &
        trim(vbfnloreference(2))//achar(10)// &
        'and process-specific references'//achar(10)// &
        trim(vbfnloreference(3))

      return
      end

!*************************************************************************  

      SUBROUTINE VBFNLO_BLHA2Amp(pp)
!*************************************************************************
!     Maps BLHA subprocess to internal m2s routine
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision pp(0:4,blha_numptcl(blha_curproc))
! local variables
      double precision p(0:3,max_p,max_kin)
      double precision v(0:3,max_v,max_kin)
      integer i, mu, k, numhelsum, ranhelcomb
! amplitude arguments
      double precision x(nx) 
      double precision rnd(max_PS_dim)
      integer nlo, ps_number
! external
      double precision dotrr, mjj2
      external dotrr, mjj2


      if (blha_amptype(blha_curproc) .eq. 1) then
        nlo = 1
        doVirtuals = .true.
      else
        nlo = 0
        doVirtuals = .false.
      endif
      ps_number = 1
      x(1) = 0.5d0 ! dummy
      x(2) = 0.5d0 ! dummy

! map partons to p
      do i=1,blha_numparton(blha_curproc)
        do mu=0,3
          p(mu,i,1) = pp(mu,blha_particlemap( &
                            i,blha_cursubproc,blha_curproc))
        enddo
      enddo
! clear others to zero -- in case they are used somewhere in the code
      do i=blha_numparton(blha_curproc)+1,max_p
        do mu=0,3
          p(mu,i,1) = 0d0
        enddo
      enddo
      do k=2,max_kin
        do i=1,max_p
          do mu=0,3
            p(mu,i,k) = 0d0
          enddo
        enddo
      enddo
! map elweak to v
      do i=1,blha_numelweak(blha_curproc)
        do mu=0,3
          v(mu,i,1) = pp(mu,blha_particlemap( &
                            i+blha_numparton(blha_curproc), &
                            blha_cursubproc,blha_curproc))
        enddo
      enddo
! clear others to zero -- in case they are used somewhere in the code
      do i=blha_numelweak(blha_curproc)+1,max_v
        do mu=0,3
          v(mu,i,1) = 0d0
        enddo
      enddo
      do k=2,max_kin
        do i=1,max_v
          do mu=0,3
            v(mu,i,k) = 0d0
          enddo
        enddo
      enddo

      if (blha_ranhelsum) then
        numhelsum = 1
        ranhelcomb =  &
          int(blha_multsubproc(blha_cursubproc,blha_curproc) &
              *blha_helrand)
      else 
        numhelsum = blha_multsubproc(blha_cursubproc,blha_curproc)
        ranhelcomb = 0
      endif
      do i=1,numhelsum
        blha_ranhelcomb = ranhelcomb+i
        call Amplitude(rnd,p,x,v,nlo,ps_number)
      enddo

      return
      end

!*************************************************************************  

      LOGICAL FUNCTION VBFNLO_ClosestOnshell(v,proc,subproc)
!*************************************************************************
!     check for subprocess closest to on-shell
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      double precision v(0:3,max_v,max_kin)
      integer proc,subproc
! local variables
      integer bos,sub,smallest
      double precision minv(0:max_subproc)
      character*255 errmsg
! external
      double precision mjj2
      external mjj2

! get product of inverse Breit-Wigners for each subproc
      do sub=1,blha_numsubproc(proc),blha_pssubstep(proc)
        minv(sub) = 1
        do bos=1,blha_numbosons(proc)
          if (abs(blha_bosons(bos,sub,proc)) .eq. 23) then
            minv(sub) = minv(sub) *                                  &
              ( (mjj2(v(0,blha_invmap(                               &
                            blha_particlemap(                        &
                              blha_numparton(proc)+2*bos-1,          &
                              subproc,proc),                         &
                            sub,proc)-blha_numparton(proc)           &
                         ,1),                                        &
                      v(0,blha_invmap(                               &
                            blha_particlemap(                        &
                              blha_numparton(proc)+2*bos,            &
                              subproc,proc),                         &
                            sub,proc)-blha_numparton(proc)           &
                         ,1))                                        &
                 -xmz**2)**2 + (xmz*xgz)**2 )
          else if (abs(blha_bosons(bos,sub,proc)) .eq. 24) then
            minv(sub) = minv(sub) *                                  &
              ( (mjj2(v(0,blha_invmap(                               &
                            blha_particlemap(                        &
                              blha_numparton(proc)+2*bos-1,          &
                              subproc,proc),                         &
                            sub,proc)-blha_numparton(proc)           &
                         ,1),                                        &
                      v(0,blha_invmap(                               &
                            blha_particlemap(                        &
                              blha_numparton(proc)+2*bos,            &
                              subproc,proc),                         &
                            sub,proc)-blha_numparton(proc)           &
                         ,1))                                        &
                 -xmw**2)**2 + (xmw*xgw)**2 )
          else
            write(errmsg,'(A,I2)')  &
             "Unknown boson in phase-space: ", blha_bosons(bos,sub,proc)
            call BLHA_error(errmsg,__FILE__,__LINE__)
          endif
        enddo
      enddo

! look for smallest
      minv(0) = minv(1)
      smallest = 1
      do sub=1+blha_pssubstep(proc),blha_numsubproc(proc),blha_pssubstep(proc)
        if ( minv(sub) .lt. minv(0) ) then
          smallest = sub
          minv(0) = minv(sub)
        endif
      enddo

      if (smallest .eq. subproc) then
        VBFNLO_ClosestOnshell = .true.
      else
        VBFNLO_ClosestOnshell = .false.
      endif

      return
      end

!*************************************************************************

      SUBROUTINE BLHA_setDR(DR)
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer DR

      if (DR.ne.0) then
        blha_DR=.true.
      else
        blha_DR=.false.
      endif
      call BLHA_setNc(blha_Nc)

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_setNf(Nf)
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer Nf

      nfl = Nf
      call BLHA_setNc(blha_Nc)

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_setNc(Nc)
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer Nc

      blha_Nc=Nc
! N==-1 means N_c->Infinity
      if (blha_Nc .eq. -1) then
        blha_CA = 3
        blha_CF = 3/2d0
      else if (blha_Nc .eq. 0) then
        blha_CA = 0
        blha_CF = 0
        write(*,*) "VBFNLO Warning: BLHA_setNc called with Nc==0"
      else
        blha_CA = blha_Nc
        blha_CF = (blha_Nc**2-1)/(2d0*blha_Nc)
      endif
      blha_gammaQuark = 3/2d0*blha_CF
      blha_gammaGluon = 11/6d0*blha_CA-1/3d0*nfl
      blha_KQuark = (7/2d0-pi**2/6d0)*blha_CF
      blha_KGluon = (67/18d0-pi**2/6d0)*blha_CA-5/9d0*nfl
      if (blha_DR) then
        blha_tgammaQuark = -blha_CF/2d0  
        blha_tgammaGluon = -blha_CA/6d0
      else
        blha_tgammaQuark = 0
        blha_tgammaGluon = 0
      endif

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_setCouplingsOff()
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      blha_couplingsoff = .true.

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_setEWRenormScheme(scheme)
!*************************************************************************
!     Initialize variables
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer scheme

      blha_ewrenormscheme = scheme

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_cctree(i,j,res)
!*************************************************************************
!     calculates the index of the colour-correlation matrix
!     for colour correlations of particle i and j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer i,j
      double precision res
! local variables
      integer iindex,jindex,ccindex

      iindex = blha_particlemap( &
                 blha_diagphysmap(i,blha_cursubproc,blha_curproc), &
                 blha_cursubproc,blha_curproc)-1
      jindex = blha_particlemap( &
                 blha_diagphysmap(j,blha_cursubproc,blha_curproc), &
                 blha_cursubproc,blha_curproc)-1

      if (iindex<jindex) then
        ccindex = 1+iindex+jindex*(jindex-1)/2
      else
        ccindex = 1+jindex+iindex*(iindex-1)/2
      endif
      blha_amp(ccindex) = blha_amp(ccindex) + res

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_sctree(i,j,res)
!*************************************************************************
!     calculates the index of the colour-correlation matrix
!     for colour correlations of index i and j
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      integer i,j
      double complex res
! local variables
      integer iindex,jindex,scindex

      iindex = blha_particlemap( &
                 blha_diagphysmap(i,blha_cursubproc,blha_curproc), &
                 blha_cursubproc,blha_curproc)-1
      jindex = blha_particlemap( &
                 blha_diagphysmap(j,blha_cursubproc,blha_curproc), &
                 blha_cursubproc,blha_curproc)-1

      scindex = 1+2*iindex+2*blha_numptcl(blha_curproc)*jindex
      blha_amp(scindex)   = blha_amp(scindex)   + dreal(res)
      blha_amp(scindex+1) = blha_amp(scindex+1) + dimag(res)

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_error(errmsg,filename,lineno)
!*************************************************************************
!     abort with an error message
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      character(*) errmsg, filename
      integer lineno

      write(*,*) "Unrecoverable error in VBFNLO BLHA interface:"
      write(*,*) "file",filename(1:len(filename)),",line",lineno,":"
      write(*,*) errmsg(1:len(errmsg))
      stop

      return
      end

!*************************************************************************  

      SUBROUTINE BLHA_amptypeerror(amptype,filename,lineno)
!*************************************************************************
!     abort due to wrong amplitude type in function
!*************************************************************************
      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/BLHAhelper.inc"

      character(*) filename
      integer amptype, lineno
! local variables
      character*255 errmsg

      write(errmsg,'(A,I2)')  &
       "Function called with wrong amplitude type number ", amptype
      call BLHA_error(errmsg,filename,lineno)

      return
      end

!*************************************************************************  

