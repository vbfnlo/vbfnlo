c*************************************************************************
      SUBROUTINE InitCouplings
c*************************************************************************
c     This is the interface to the koppln-routine for the initialization
c     of coupling constants. 
c*************************************************************************
       use globalvars, only: lglobalprint
       use readinput
       use VBFNLO_HIGGSTO, only: VBFNLO_HIGGSTO_gridHt

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/process.inc"
#include "VBFNLO/utilities/ckm.inc"

** Switches for anomalous couplings, kaluza-klein
      logical anom,kkswitch, spin2switch

** Used for determination of slha input file
      character*250 GetInputPath
      external GetInputPath
      character*250 path


** Used for determining the strong coupling constant:
      double precision alphas5
      external alphas5


** Function which determines whether a parameter has already been set with a 
** reasonable value
      logical replace
      external replace

** 1/Alfa_QED: this is the default read-in value, but is immediately converted
** to alfa
      double precision invalfa

** Coupling parameters
      double precision e, g2, s, c, z, w, q, g

** Branching ratios and widths
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
c definitions for second Higgs
      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      integer gsba, gcba
      logical with_secondhiggs, use_secondhiggs
      common /COUPL_TWOHIGGS/ gsba, gcba,
     &                        with_secondhiggs, use_secondhiggs
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
c temporary variables
      double precision tBHWW,tBHZZ,tBHGG,tBHTT,tBHBB,tBHCC,tBHTAU,tBHMU,
     &                 tBHGAM, tBHGAMZ, tXMH, tXGH, tsinba

      double precision srt
      external srt

      call clearwidths
      
      call loadfile("vbfnlo.dat",.true.)

      if (lglobalprint) then
      print *," "
      print *,"              Physics parameters"
      print *,"-----------------------------------------------"
      endif

*sophy: read model (SM/MSSM) and HiggsType (SM/h0/H0/A0)
      call read_Int("HTYPE",HiggsType,0)
      call read_Int("MODEL",model,1)
      call read_logical("ANOM_CPL",anom,.false.) 
      call read_logical("KK_MOD",kkswitch,.false.)
      call read_logical("SPIN2",spin2switch,.false.)      


c set ckm values, this could be an input at some point

      if(luseckm) then
!         call read_realList("CKM",9,vckm,vckm,.false.)
       vckm = reshape( (/0.97428,  0.22535,  0.0,
     &                   0.22535,  0.97428,  0.0,
     &                   0.0,      0.0,      1.0  /),(/3,3/) )

        vckm = vckm**2
        ! later only squared matrix elements are needed for folding
        ! with pdfs
      else
         vckm = reshape( (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 /),(/3,3/) )
      endif

      call closefile

      if (HiggsType .ne. 0) then
         call loadfile("susy.dat",.true.)
*sophy: determine whether FeynHiggs should be used
         call read_logical("FEYNH_SWITCH",FeynH_switch,.false.)
*sophy: determine whether we're working from a slha file
         call read_logical("SLHA_SWITCH",slha_switch,.true.)

         if ((.not. FeynH_switch) .and. (.not. slha_switch)) then
            write(*,*)'Please set either FeynH_switch to .true.'
            write(*,*)'or SLHA_SWITCH to .true. when working with'
            write(*,*)'MSSM-type Higgs bosons (HTYPE = 1, 2, 3)'
            stop
         end if

         if (slha_switch) then
*sophy: setting SLHA file path name
            call read_string("SLHA_FILE",slhaFILE,"sps1a.slha")
            path = GetInputPath()
            slhaFILE = trim(path)//"/"//slhaFILE
*         write(*,*)'Reading from SLHA file =', slhaFILE
*sophy:  reading SLHA file. Argument is debugging
            call readSLHA(1)
         end if
         call closefile
      else
         FeynH_switch = .false.
         slha_switch = .false.
      end if


      call loadfile("vbfnlo.dat",.true.)

*sophy: reading higgs mass if it hasn't been set by slha
      if (HiggsType .eq. 0) then
         call read_Real("HMASS",xmh,126d0)
         if (model .eq. 2) then
            write(*,*)'When working in the MSSM, please choose'
            write(*,*)'HTYPE = 1, 2 or 3'
            stop
         end if
      end if

* reading higgs width
      if (replace(xgh,3) .or. (HiggsType .eq. 0)) then
         if (higgsscheme .ne. 2) then   ! higgsscheme=2=CPS
           call read_Real("HWIDTH",xgh,-999d0)
c get width from HXSWG Table
           if (xgh .eq. -998d0) call VBFNLO_HIGGSTO_gridHt(XMH,XGH)
         endif
      end if

* reading W width
      if (replace(xgw,3)) then
        call read_Real("WWIDTH",xgw,-999d0,.true.)
      end if
* reading Z width
      if (replace(xgz,3)) then
        call read_Real("ZWIDTH",xgz,-999d0,.true.)
      end if

c read in additional parameters for two-Higgs model
      if (model .eq. 3) then
        with_secondhiggs = .true.
        use_secondhiggs = .false.
        call read_Real("H2MASS",xmh2,126d0)
        if (higgsscheme .ne. 2) then   ! higgsscheme=2=CPS
          call read_Real("H2WIDTH",xgh2,-999d0)
c get width from HXSWG Table
          if (xgh2 .eq. -998d0) call VBFNLO_HIGGSTO_gridHt(XMH2,XGH2)
        endif
        call read_Real("SIN2BA",sinba,1d0)
        call read_Real("COS2BA",cosba,-999d0)
        if (replace(cosba,3)) then
          cosba = srt(1-sinba)
        else
          cosba = srt(cosba)
        endif
        sinba = srt(sinba)
      else
        with_secondhiggs = .false.
        sinba = 1d0
        cosba = 0d0
      endif

*sophy: reading input quark and tau masses if they haven't been set by slha file
      if (replace(xmt,0)) call read_real("TOPMASS",xmt,172.4d0)
      if (replace(xmb,0)) call read_real("BOTTOMMASS",xmb,4.855d0)
      if (replace(xmc,0)) call read_real("CHARMMASS",xmc,1.65d0)
      if (replace(xmtau,0)) call read_real("TAU_MASS",xmtau,1.77684d0)

*sophy: determing electroweak scheme, reading in ew values if they're not 
* set by slha file.  All the lovely details of EWSCHEME are on the wiki.
      call read_Int("EWSCHEME",ewscheme,3)

      Select Case (EWSCHEME)
      Case(1)
         if (replace(gf,0)) call read_Real("FERMI_CONST",gf,
     -        1.16637d-5)
         if (replace(alfa,0)) then
            call read_Real("INVALFA",invalfa,-1D0)
            alfa = 1D0/invalfa
         end if
         if (alfa .eq. -1D0) call read_Real("ALFA",alfa,
     -        7.7552841117d-3)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
         xmw = -1.d0  ! xmw and sin2w are calculated in this scheme
         sin2w = -1.d0
      CASE(2)
         if (replace(gf,0)) call read_Real("FERMI_CONST",
     -        gf,1.16637d-5)
         if (replace(sin2w,0)) call read_Real("SIN2W",sin2w,
     -        0.222646d0)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
         xmw = -1.d0  ! xmw and alfa are calculated in this scheme
         alfa = -1.d0
      CASE(3)
         if (replace(gf,0)) call read_Real("FERMI_CONST",gf,
     -        1.16637d-5)
         if (replace(xmw,0)) call read_Real("WMASS",xmw,80.398d0)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
         alfa = -1.d0 ! alfa and sin2w are calculated in this scheme
         sin2w = -1.d0
      CASE(4)  ! Nothing is calculated in this scheme: all are inputs
         if (replace(gf,0)) call read_Real("FERMI_CONST",gf,
     -        1.16637d-5)
         if (replace(alfa,0)) then
            call read_Real("INVALFA",invalfa,-1D0) 
            alfa = 1D0/invalfa
         end if
         if (alfa .eq. -1D0) call read_Real("ALFA",alfa,
     -        7.2973525376d-3)
         call read_Real("SIN2W",sin2w,0.222646d0)
         if (replace(xmw,0)) call read_Real("WMASS",xmw,80.398d0)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
      CASE(5)  ! Input value of alfa_qed should be alfa(MZ)
         if (replace(gf,0)) call read_Real("FERMI_CONST",gf,
     -        1.16637d-5)
         if (replace(alfa,0)) then
            call read_Real("INVALFA",invalfa,-1D0) 
            alfa = 1D0/invalfa
         end if
         if (alfa .eq. -1D0) call read_Real("ALFA",alfa,
     -        7.7552841117d-3)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
         if (replace(xmw,0)) call read_Real("WMASS",xmw,80.398d0)
      CASE(6)  ! Input value of alfa_qed should be alfa(0)
         if (replace(gf,0)) call read_Real("FERMI_CONST",gf,
     -        1.16637d-5)
         if (replace(alfa,0)) then
            call read_Real("INVALFA",invalfa,-1D0) 
            alfa = 1D0/invalfa
         end if
         if (alfa .eq. -1D0) call read_Real("ALFA",alfa,
     -        7.2973525376d-3)
         if (replace(xmz,0)) call read_Real("ZMASS",xmz,91.1876d0)
         if (replace(xmw,0)) call read_Real("WMASS",xmw,80.398d0)
      Case Default
         print*,"Invalid choice for EWSCHEME"
      End Select


** Setting 'basic' alfas - alfas at MZ at NLO, according to whichever PDF
** we're using    
      alfas = alphas5(xmz**2D0,1)


*sophy: setting 'sector': which corrections are included
      if (ewcor_switch) then

*sophy: read in value of electroweak corrections and any approximations wanted
         call read_Int("EW_APPROX",ewapprox,5)

         if (ewapprox .eq. 1) then
            sector = 1
         else if (ewapprox .eq. 2) then
            sector = 2
         else 
            if (model .eq. 1) then  ! SM
               sector = 3
            else  ! MSSM
               sector = ewapprox
               if (ewapprox .eq. 0) then
                  sector = 6
               end if
            end if
         end if
      else
         sector = 0
      end if
      

*sophy: killing run if ew cor have been selected with non VBF Higgs production
      if (ewcor_switch .and. (procID .gt. 109 .and. procID.ne.1010)) then
      Print*, "********************************"//
     &        "*********************************" 
      Print*, "Clash: Electroweak corrections  "//
     &        "are only included for VBF Higgs "// 
     &        "production! "
      Print*, "********************************"//
     &        "*********************************" 
      STOP
      endif


      if ((procID .eq. 191) .and. (.not. spin2switch)) then
         Print*, "***************************************"//
     &        "***************************************" 
         Print*, "Sorry: process 191 is only available "//
     &           "as a spin-2 process!" 
         Print*, "If you want ppbar -> Hjj -> AAjj please "//
     &           "use process = 101!"  
         Print*, "***************************************"//
     &        "***************************************" 
         STOP
      endif

      if (anom .and. kkswitch
     &     .or. anom .and. spin2switch .or. 
     &    kkswitch .and. spin2switch) then
      Print*, "***************************************"//
     &        "***************************************" 
      Print*, "Clash: Two of ANOM_CPL, KK_MOD and Spin2 "//
     &        "in vbfnlo.dat are set to 'true'!" 
      Print*, "***************************************"//
     &        "***************************************" 
      STOP
      endif

** sophy: kill commands if electroweak corrections are chosen as well as
**        anomalous couplings or KK
      if (anom .and. ewcor_switch) then
         print *, "********************************"//
     &            "*********************************" 
         print*, "Clash: Both ANOM_CPL and EWCOR_SWITCH "//
     &           "in vbfnlo.dat are set to 'true'!" 
         print*, "********************************"//
     &           "*********************************" 
         stop
      endif
      if (kkswitch .and. ewcor_switch) then
         print *, "********************************"//
     &            "*********************************" 
         print*, "Clash: Both KK_MOD and EWCOR_SWITCH "//
     &           "in vbfnlo.dat are set to 'true'!" 
         print*, "********************************"//
     &           "*********************************" 
         stop
      endif
** sophy: kill command if anomalous couplings and MSSM both
**        chosen
      if (anom .and. (model .eq. 2)) then
         print *, "********************************"//
     &            "*********************************" 
         print*, "Clash: We cannot run with anomalous"//
     &           "couplings in the MSSM"
         print*, "********************************"//
     &           "*********************************" 
         stop
      endif
** sophy: kill command if kaluza-klein and MSSM both
**        chosen
      if (kkswitch .and. (model .eq. 2)) then
         print *, "********************************"//
     &            "*********************************" 
         print*, "Clash: KK_MOD and the MSSM have"//
     &           "both been chosen!"
         print*, "********************************"//
     &           "*********************************" 
         stop
      endif

      if (with_secondhiggs .and. 
     &  (procID .lt. 200 .or. procID.ge.300)) then
      Print*, "********************************"//
     &        "*********************************" 
      Print*, "Two-Higgs model only implemented "//
     &        "for VBF -> VV (V in W,Z) "// 
     &        "production! "
      Print*, "********************************"//
     &        "*********************************" 
      STOP
      endif


** Writing out SM / MSSM info 
      if (lglobalprint) then
      if (model .eq. 1) then
         write(*,*)'  '
         write(*,*)'Running in the SM'
      else if (model .eq. 2) then
         write(*,*)'  '
         write(*,*)'Running in the MSSM'
      else if (model .eq. 3) then
         write(*,*)'  '
         write(*,*)'Running in the Two-Higgs model'
      else
         write(*,*)'  '
         write(*,*)'Please choose either MODEL = 1 (SM), 2 (MSSM) '//
     &             'or 3 (Two-Higgs)'
         stop
      end if
      endif

      call closefile

      with_kk=.false.
      with_spin2 = .false.


** Setting the electroweak parameters
      call setEWpara(e,g2,s,c,z,w,q,g)

      call InitAnomCouplings(anom)

      if (model.eq.3) then
c save Higgs1 parameters to temporary and overwrite with Higgs2 parameters
        tBHWW   = BHWW   
        tBHZZ   = BHZZ   
        tBHGG   = BHGG   
        tBHTT   = BHTT   
        tBHBB   = BHBB   
        tBHCC   = BHCC   
        tBHTAU  = BHTAU  
        tBHMU   = BHMU   
        tBHGAM  = BHGAM  
        tBHGAMZ = BHGAMZ 
        tXMH    = XMH    
        tXGH    = XGH    
        tsinba  = sinba  
        BHWW   = BH2WW   
        BHZZ   = BH2ZZ   
        BHGG   = BH2GG   
        BHTT   = BH2TT   
        BHBB   = BH2BB   
        BHCC   = BH2CC   
        BHTAU  = BH2TAU  
        BHMU   = BH2MU   
        BHGAM  = BH2GAM  
        BHGAMZ = BH2GAMZ 
        XMH    = XMH2    
        XGH    = XGH2    
        sinba  = cosba  
c call KOPPLN
        CALL KOPPLN(0,e,g2,s,c,z,w,q,g)
c move result into Higgs2 parameters and copy back
        BH2WW   = BHWW  
        BH2ZZ   = BHZZ  
        BH2GG   = BHGG  
        BH2TT   = BHTT  
        BH2BB   = BHBB  
        BH2CC   = BHCC  
        BH2TAU  = BHTAU 
        BH2MU   = BHMU  
        BH2GAM  = BHGAM 
        BH2GAMZ = BHGAMZ
        XMH2    = XMH    
        XGH2    = XGH    
        ch2ww   = ahvv(1,3,4)
        ch2zz   = ahvv(1,2,2)

        BHWW   = tBHWW   
        BHZZ   = tBHZZ   
        BHGG   = tBHGG   
        BHTT   = tBHTT   
        BHBB   = tBHBB   
        BHCC   = tBHCC   
        BHTAU  = tBHTAU  
        BHMU   = tBHMU   
        BHGAM  = tBHGAM  
        BHGAMZ = tBHGAMZ 
        XMH    = tXMH    
        XGH    = tXGH    
        sinba  = tsinba  
      endif
      CALL KOPPLN(2,e,g2,s,c,z,w,q,g)
      call ctrans(xmb)

      if (kkswitch) then
#ifdef WITH_KK
         Select Case (procID)
         Case(WPWMjj, WPZjj, WMZjj, ZZjj_ll, ZZjj_lnu,
     -           WWZ, ZZWP, ZZWM, WWWP, WWWM)
            with_kk= kkswitch
            call kk_coupl_inp
	   CASE DEFAULT
         print *," Error: KK mode not available for this process!"
	   stop		
         END Select
#else
         print *," Error: KK mode not available."
         print *,' reconfigure with "--enable-kk" and then recomplile.'
         stop 2
#endif
      end if

 
      if (spin2switch) then
         Select Case (procID)
         Case(ZZjj_lnu, ZZjj_ll, WPWMjj, WPZjj, WMZjj, AAjj, Sp2jj_WW,
     &        Sp2jj_ZZ_ll, AAjjjLO, Sp2jjjLO_WW, Sp2jjjLO_ZZ_ll, 
     &        WPWMjjjLO, WPZjjjLO, WMZjjjLO, ZZjjjLO_ll, ZZjjjLO_lnu,
     &        Sp2jj_ZZ_lnu, Sp2jjjLO_ZZ_lnu)
            with_spin2= spin2switch 
            call read_spin2coupl
            call spin2width
         CASE DEFAULT
            print *," Error: spin-2 couplings not available for this process!"
            stop		
         END Select
      else if (procID .eq. AAjj) then
         print *,"Sorry! Process 240 (pp->AAjj) only runs in spin-2 models"
      end if

!     call coupling settings for hadronic boson decays
      CALL coupl_haddecay

      RETURN
      END 


c*************************************************************************
      SUBROUTINE InitAnomCouplings(anom)
c*************************************************************************
c     This routine reads in anomalous couplings if requested 
c*************************************************************************

      implicit none

#include "VBFNLO/utilities/global.inc"
#include "VBFNLO/utilities/koppln.inc"
#include "VBFNLO/utilities/process.inc"

      logical anom

** Dummy flag to enable the output of a warning if anomalous couplings are 
** requested for a process that does not support them
      logical anom_called


      anom_called = .false.

      with_anomHiggs = .false.
      with_anom = .false.

#ifdef WITH_VBF
      if (anom) then
         Select Case (procID)
         Case(Hjj, Hjj_AA, Hjj_mu, Hjj_tau, Hjj_bbar,  
     -           Hjj_WW, Hjj_ZZ_ll, Hjj_ZZ_lnu,
     -           Hjj_WPhadWM, Hjj_WPWMhad, Hjj_ZZhad)
            with_anomHiggs=anom
            call read_anomHVV
            anom_called = .true.
         END SELECT
      endif
#endif


#ifdef WITH_GGF
      if (anom) then
         Select Case (procID)
         Case(GFHjj_WW, GFHjj_ZZ_ll, GFHjj_ZZ_lnu)
            with_anomHiggs=anom
            call read_anomHVV
            anom_called = .true.
         END SELECT
      endif
#endif


#ifdef WITH_VBF
      if (anom) then
         Select Case (procID)
         Case(WPWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,WPWPjj,WMWMjj,WPAjj,WMAjj,ZAjj,ZAjj_n,
     &        WPhadWMjj, WPWMhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, ZZhadjj, WPhadWPjj, WMhadWMjj,
     &        Ajj, Zjj_l, Zjj_nu, WPjj, WMjj,WPAjjjLO,WMAjjjLO,ZAjjjLO)
            call read_anomVcouplings
            anom_called = .true.
            with_anom = anom
            with_anomHiggs = anom
         END SELECT
      endif
#endif


#ifdef WITH_TRIBOSON
      if (anom) then
         Select Case (procID)
         Case(WWZ,ZZZ,ZZWP,ZZWM,WWWP,WWWM,WWA,ZZA,ZZnA,WPZA,WMZA,WPAA,WMAA,ZAA,ZnAA,AAA,
     &        WPAAj,WMAAj,
     &        WPhadWMZ, WPWMhadZ, WWZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, 
     &        WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM, ZZZhad,
     &        WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA)
            call read_anomVcouplings
            anom_called = .true.
            with_anomHiggs = anom
            with_anom = anom
         END SELECT
      endif
#endif

#ifdef WITH_DIBOSON
      if (anom) then
c anomalous vector couplings
         Select Case (procID)
         Case(WW,WWjLO,WPZ,WMZ,WPA,WMA,WPZjLOb,
     &        WPhadWMlep,WPlepWMhad,
     &        WPhadZ,WPZhad,WMhadZ,WMZhad,
     &       WMZjLOb,WPAjLOb,WMAjLOb,
     &       ZA,ZZ,ZAjLO,ZZjLO)
            with_anom = anom
            call read_anomVcouplings
            anom_called = .true.
         END SELECT
c anomalous Higgs couplings
         Select Case (procID)
         Case(WW,ZZ,ZA,AA,GFWW,GFZZ,GFZA,GFAA,
     &        WPhadWMlep,WPlepWMhad,GFWPhadWMlep,GFWPlepWMhad,
     &        ZZhad,GFZZhad,
     &        WWj,ZZj,GFWWj,GFZZj)
            with_anomHiggs=anom
            call read_anomHVV
            anom_called = .true.
         END SELECT
      endif
#endif

#ifdef WITH_DIBOSONJET
      if (anom) then
         Select Case (procID)
         Case(WMAJ,WPAJ,WMZJ,WPZJ,
     &        WPhadZJ, WPZhadJ, WMhadZJ, WMZhadJ)
            with_anom=anom
            call read_anomVcouplings
            anom_called = .true.
         END SELECT
      endif
#endif


      if (anom) then
         Select Case (procID)
         Case(
     1     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll, WPH_ZZ_lnu,
     1     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll, WMH_ZZ_lnu,
     1     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll, WPHJ_ZZ_lnu,
     1     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll, WMHJ_ZZ_lnu
     1       )
            call read_anomVcouplings
            anom_called = .true.
            with_anomHiggs = anom
            with_anom = anom
         END SELECT
      endif


      if (anom .and. (.not. anom_called)) then
         write(*,*)'  '
         write(*,*)'WARNING!'
         write(*,*)'Anomalous couplings have been selected'
         write(*,*)'with a process that does not accept them.'
         write(*,*)'Standard couplings will be used.'
         write(*,*)'   '
      end if

      if (.not.with_anom) call set_anomVcouplings_zero     ! set the anomalous vector boson couplings to zero

      return
      end

c*************************************************************************


C========= subroutine KOPPLN  ==============================================
C
C  Adapted from RKWW by R.Kleiss
C
C  H --> gamma gamma decay added by D. Rainwater,
C        from HHG, Mar 97
c  #included D5 operators by TFigy, June 04
C  H --> b bbar decay modified by D. Rainwater,
C        from HHG, Mar 98
C
c  H --> V*V* -> all added by T.Figy, Oct 2003
c
C  Calculate the fermion-fermion-boson and three-boson couplings
C  in the standard model, without colour factors.
C  also calculate mass**2 and mass*width of the bosons.
C  V   = FFB vector coupling
C  A   = FFB axial-vector coupling
C  CLR = FFB left- and righthanded couplings
C  B   = BBB coupling
C  The fermion indices are:  1: neutrino of electron (muon,tau)
C                            2: electron (muon,tau)
C                            3: up quark (charm,top)
C                            4: down quark (strange,bottom)
C  The boson indices are     1: photon
C                            2: Z0 boson
C                            3: W+ boson
C                            4: W- boson
C                            5: gluon
C                            6: Higgs
C  The helicity indices are  -1: lefthanded
C                            +1: righthanded
C  No Kobayashi-Maskawa mixing is implemented here.
C  The W+- are identified as outgoing from the BBB vertex.
C  XM2 = boson mass**2
C  XMG = boson mass*width
C  Also the branching ratios of W and Z decay are calculated in the SM
C
C  The fundamental parameters:
C  ALFA = The Q.E.D. fine structure constant
C  SIN2W = sin**2 of the weak mixing angle
C  ALFAS = the QCD fine structure constant
C  XMZ = mass of the Z0 boson in GeV
C  XMW = mass of the W+- bosons in GeV
C  XMH = mass of the Higgs boson in GeV
C  XMT = mass of the top quark in GeV
C
C  printout is steered by INFO:
C   <=0 :   no printout
C   <=1 :   printout of couplings only
C    =2 :   printout of input parameters and branching ratios
C   else:   printout of couplings and W/Z branching ratios
C
C===========================================================================
C

*** IMPORTANT NOTE!
*** In order to avoid hard-wiring values of the electroweak parameters when
*** using anomalous couplings, the setting of these parameters has been moved to
*** a new subroutine, 'setEWpara', which can be found in utilities/parameters.F


      SUBROUTINE KOPPLN(INFO,e,g2,s,c,z,w,q,g)
          use globalvars, only: lglobalprint
          use VBFNLO_HIGGSTO, only: VBFNLO_HIGGSTO_pole

      implicit none

      double complex  f_tau_t, f_tau_b, f_tau_W, F_t, F_b, F_W, dum
      double precision  beta0

      double precision  zero
      parameter ( zero = 0.0d0 )
      integer  info, i, j, k
      logical ldebug
c
c INPUT
c
#include "koppln.inc"
#include "global.inc"
#include "mssm.inc"
c
c OUTPUT
c
      double precision XM2s(6),XMGs(6),
     & BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     & BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,BHGAM
      COMMON /BKOPOUshort/ XM2s,XMGs
 

** clrCT is the counterterm for qqV
      double complex clrCT(3:4,2:4,-1:1)
      common /bkopouCT/ clrCT
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)

c     ahvv(i,v1,v2) i = 1,2,3 are coefficients of eq 1 of 
c     Physics Letters B 591,297
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL
      real*8 treefacW, treefacZ, loopfac 
      common/lhcoup/ treefacW, treefacZ, loopfac
      double precision g5hvv(2,1:4,1:4),lambda5
      common/hcoupl/g5hvv,lambda5
c     form factor for hvv couplings
      integer ffac
      double precision m2ff, mff 
      logical lff
      common/ formfacmass/ m2ff,mff,ffac,lff
C
C local variables
C

** coupling parameters, widths, ..
      double precision e,g,s,c,g2,z,w,q, 
     & gwh, gzh, decw, gwne, gwud, gwtb, decz, gznn, gzee, gzuu, gzdd,
     & gztt, xgz,xgw,xlw,xlz,
     & ghww,ghzz,ghtt,betat,qcdt,
     & ghcc,ghbb,ghtau,ghgg,ghgam,tau_t,eta_p,eta_m,tau_b,tau_w,xgh
      double precision xmMU,betaMU,ghmu,gamt,alphas5,alphas
      double precision A_12_t, A_12_b
      external alphas5, alphas
      double precision BHGAMZ,GHGAMZ
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      double precision xmax,xmin
      double complex II_1,II_2,A_sum,A_w,A_t,A_b
      double precision cos2w,tan2w,cotw,lamb_b,lamb_t,lamb_w

c parameters needed for calculation of H -> B BBAR
      double precision alphas_H, alphas_B, BMASS_POLE, BMASS_B
      double precision alphas_C, CMASS_POLE, CMASS_C
      double precision DELTA_qq, DELTA_H, C1, C2
      integer Nf

c function to calculate phase space reduction for heavy fermionic decays
      double precision heavy
      external heavy


** Function which determines whether a parameter has already been set with a 
** reasonable value
      logical replace
      external replace

      double precision hww_width,gaus
      external hww_width, gaus

      double precision relerr, abserr

c
c mass dependence for H--->W*W,Z*Z (From Gunion et al., Higgs Hunters Guide)
c
c      FH(X) = - (1.d0 - X**2)*(23.5d0*X**2 - 6.5d0 + 1.d0/X**2) 
c     &        - 3.d0*(1.d0 - 6.d0*X**2 + 4.d0*X**4)*LOG(X)
c     &        + 3.d0*(1.d0 - 8.d0*X**2 + 20.d0*X**4)/
c     &                                   SQRT(4.d0*X**2 - 1.d0)
c     &        * ACOS( (3.d0*X**2 - 1.d0)/(2.d0*X**3) )
c      PI = 4.d0*ATAN(1.d0)

      ldebug = .false.   ! alfa.eq.0.122431


c     initialize ahvv to zero
      A_sum=(0.0d0,0.0d0)
      dum =(0.0d0,0.0d0)
      if(.not.with_anomHiggs) then
         do i = 1,4
            do j = 1,4
               ahvv(1,i,j) = (0.0d0,0.0d0)
               ahvv(2,i,j) = (0.0d0,0.0d0)
               ahvv(3,i,j) = (0.0d0,0.0d0)
               g5hvv(2,i,j) =0.0d0
               g5hvv(1,i,j) =0.0d0
            enddo
         enddo
         treefacW = 1.0d0
         treefacZ = 1.0d0
         loopfac = 1.0d0
         lambda5 = 1D0
         lff=.false.
      endif


c     initialize ahvv to zero
      if(with_anomHiggs) call fillanomhcoupl(1,A_sum,dum,loopfac)

c     fill ahvv coefficients h-w-w and h-z-z
      if(with_anomHiggs) call fillanomhcoupl(2,A_sum,dum,loopfac)


c THE FFB VECTOR COUPLING
c
* photon-fermion coupling
      V(1,1) = zero
      V(2,1) = -E
      V(3,1) = 2.d0*E/3.d0
      V(4,1) = -E/3.d0
* Z-fermion coupling
      V(1,2) =  Z
      V(2,2) = -Z*(1.d0 - 4.d0*SIN2W)
      V(3,2) =  Z*(1.d0 - 8.d0/3.d0*SIN2W)
      V(4,2) = -Z*(1.d0 - 4.d0/3.d0*SIN2W)

* W-fermion coupling
      DO 1 K = 1,4
         DO 1 J = 3,4
            V(K,J) = W
 1    continue
* gluon-fermion coupling
      V(1,5) = zero
      V(2,5) = zero
      V(3,5) = G
      V(4,5) = G

c
c THE FFB AXIAL-VECTOR COUPLING
c
* photon/gluon - fermion coupling
      DO 2 K = 1,4
         DO 2 J = 1,5,4
            A(K,J) = zero
 2    continue
* Z-fermion coupling
      A(1,2) = -Z
      A(2,2) =  Z
      A(3,2) = -Z
      A(4,2) =  Z
* W-fermion coupling
      do K = 1,4
         do J = 3,4
            A(K,J) = -W
         enddo
      enddo

* switching sin-theta-w conventions for the mssm
      if (model .eq. 2) then
         do k = 1, 4
            do j = 2, 4
               V(k,j) = -V(k,j)
               A(k,j) = -A(k,j)
            end do
         end do
      end if

c
c THE FFB LEFT- AND RIGHTHANDED COUPLINGS
c
      DO K = 1,4
         DO J = 1,5
            CLR(K,J,0) = zero
            CLR(K,J,-1)= V(K,J)-A(K,J)
            CLR(K,J,1) = V(K,J)+A(K,J)
         end do
      end do
      if (ldebug) then
         print*," electron L coupling "
         print*, (clr(2,i,-1),i=1,5)
         print*," electron R coupling "
         print*, (clr(2,i,+1),i=1,5)
      endif

* The ffb left and right handed counterterm couplings:  These are needed if
* we're calculating ew corrections.  If we're including corrections from non-
* (s)fermions, an effective coupling (which includes loop as well as counterterm
* diagrams) is used instead.  It's called from amplitudes/vvjj/qqhqq.F so that
* the momentum dependence can be taken into account
      if (ewcor_switch) then
         call qqV_ct(clrCT,model,sector)
      end if

c
c THE BBB COUPLINGS
c
      DO I = 1,6
         DO J = 1,6
            DO K = 1,6
               B(I,J,K) = zero
            end do
         end do
      end do
* W-W-photon coupling
      B(3,4,1) = -E
      B(4,1,3) = -E
      B(1,3,4) = -E
      B(1,4,3) =  E
      B(4,3,1) =  E
      B(3,1,4) =  E
* W-W-Z coupling
      B(3,4,2) = -Q
      B(4,2,3) = -Q
      B(2,3,4) = -Q
      B(2,4,3) =  Q
      B(4,3,2) =  Q
      B(3,2,4) =  Q
* altering sin-theta-w conventions for mssm
      if (model .eq. 2) then
         B(3,4,2) = Q
         B(4,2,3) = Q
         B(2,3,4) = Q
         B(2,4,3) = -Q
         B(4,3,2) = -Q
         B(3,2,4) = -Q
      end if         
* gluon coupling
      B(5,5,5) =  G
c
c Higgs couplings to W and Z
c
      GWH = G2
      GZH = GWH/C**2

      
** The triboson processes with 3 massive gauge bosons contain Higgs 
** contributions.  There should be no problem when working in the MSSM if the
** Higgs is SM-like, but out of the decoupling regime additional corrections
** would need to be made.  When working in the MSSM we thus forcibly reset
** the higgs type to 1 (h0) and use an effective Higgs mixing angle
      if (model .eq. 2) then
         if ((procID .ge. 400) .and. (procID .le. 450)) then

            HiggsType = 1

            XMH = MHiggs_loop(1)
            MHiggs(1) = XMH
            MH = XMH
            MH2 = MH**2
            
            propLevel = 0
            SBA = SBAeff           

            write(*,*)'  '
            write(*,*)'NOTE: This triboson process includes SM-like'
            write(*,*)'Higgs contributions.  This is only a good'
            write(*,*)'approximation in the decoupling region of '
            write(*,*)'the MSSM.'
            write(*,*)'Consequently, we have set HTYPE = 1 (h0)'
            write(*,*)'and use an effective Higgs mixing angle to'
            write(*,*)'calculate the Higgs couplings.'
            write(*,*)'Sin(beta-alpha) =', SBA
            write(*,*)'  '

         end if
      else if (model .eq. 3) then
        GWH = GWH*sinba
        GZH = GZH*sinba
      end if


* Setting susy Higgs couplings :
      if (model .eq. 2) then
         if ((propLevel .eq. 1) .or. (propLevel .eq.2)) then
            GWH = GWH*(Hgs(HiggsType,1)*SBA + Hgs(HiggsType,2)*CBA)
            GZH = GZH*(Hgs(HiggsType,1)*SBA + Hgs(HiggsType,2)*CBA)
         else
            if (HiggsType .eq. 1) then
               GWH = GWH*SBA
               GZH = GZH*SBA
            else if (HiggsType .eq. 2) then
               GWH = GWH*CBA
               GZH = GZH*CBA
            else if (HiggsType .eq. 3) then
               GWH = 0D0
               GZH = 0D0
            end if
         end if
      end if


* Z-Z-Higgs coupling
      B(2,2,6) = GZH
      B(2,6,2) = GZH
      B(6,2,2) = GZH
* W-W-Higgs coupling
      B(3,4,6) = GWH
      B(4,3,6) = GWH
      B(3,6,4) = GWH
      B(4,6,3) = GWH
      B(6,3,4) = GWH
      B(6,4,3) = GWH

*      write(*,*)'  '
*      write(*,*)'Higgs-Z-Z coupling =', GZH
*      write(*,*)'Higgs-W-W coupling =', GWH

      
* a1 are sm hzz and hww couplings
      ahvv(1,2,2) = treefacZ*B(6,2,2)*xmw 
      do i=3,4
         do j=3,4
            ahvv(1,i,j) = treefacW*B(6,i,j)*xmw 
         enddo
      enddo
c
c Calculate the total decay widths and branching ratios of the heavy bosons if 
c we're not using a slha file, or if they're not included in the slha file
c
c W first
c
      DECW = XMW/12.d0/PI

      if (replace(BWNE,0)) then
         GWNE = DECW*( V(1,3)**2 + A(1,3)**2 )
      else
         GWNE = BWNE*XGW
      end if

      if (replace(BWUD,0)) then
         GWUD = DECW*( V(3,3)**2 + A(3,3)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GWUD = BWUD*XGW
      end if

      if (replace(BWTB,0)) then
         GWTB = zero
         IF(XMT.GT.XMW) GOTO 31
         GWTB = GWUD*HEAVY( V(3,3)/A(3,3), XMT/XMW, zero )
      else
         GWTB = BWTB*XGW
      end if

 31   if (replace(XGW,3)) then
         XGW  = 3.d0*GWNE + 2.d0*GWUD + GWTB
      end if

      BWNE = GWNE/XGW
      BWUD = GWUD/XGW
      BWTB = GWTB/XGW
c
c Z next
c
      DECZ = XMZ/12.d0/PI

      if (replace(BZNN,0)) then
         GZNN = DECZ*( V(1,2)**2 + A(1,2)**2 )
      else
         GZNN = BZNN*XGZ
      end if

      if (replace(BZEE,0)) then
         GZEE = DECZ*( V(2,2)**2 + A(2,2)**2 )
      else
         GZEE = BZEE*XGZ
      end if

      if (replace(BZUU,0)) then
         GZUU = DECZ*( V(3,2)**2 + A(3,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GZUU = BZUU*XGZ
      end if

      if (replace(BZDD,0)) then
         GZDD = DECZ*( V(4,2)**2 + A(4,2)**2 )*3.d0*(1.d0 + ALFAS/PI)
      else
         GZDD = BZDD*XGZ
      end if

      if (replace(BZTT,0)) then
         GZTT = zero
         IF(XMT.GT.XMZ/2.d0) GOTO 32
         GZTT = GZUU*HEAVY( V(3,2)/A(3,2) , XMT/XMZ , XMT/XMZ )
      else
         GZTT = BZTT*XGZ
      end if

 32   if (replace(XGZ,3)) then
         XGZ  = 3.d0*GZNN + 3.d0*GZEE + 2.d0*GZUU + 3.d0*GZDD + GZTT
      end if
      BZNN = GZNN/XGZ
      BZEE = GZEE/XGZ
      BZUU = GZUU/XGZ
      BZDD = GZDD/XGZ
      BZTT = GZTT/XGZ

c
c     mass and width info 
c
      XM2(1) = zero
      XM2(2) = XMZ**2
      XM2(3) = XMW**2
      XM2(4) = XM2(3)
      XM2(5) = zero
      xm2(6) = xmh**2
c
      XMG(1) = zero
      XMG(2) = XMZ*XGZ
      XMG(3) = XMW*XGW
      XMG(4) = XMG(3)
      XMG(5) = zero
c     cannot compute xmg(6) at this point
c
c dummy xm2 and xmg for common block BKOPOUshort 
c
      do i = 1,6
         xm2s(i) = xm2(i)
      end do
      do i = 1,5
         xmgs(i) = xmg(i)
      end do     


c******************************************************************
c Next the Higgs decays to WW, ZZ, bb, cc, tau tau and gamma gamma.
c These are either taken from the slha file, from FH, or are 
c calculated below
c******************************************************************

** can now set the CPS Higgs width
      if (higgsscheme .eq. 2) then   ! higgsscheme=2=CPS
        if (replace(xgh,3) .or. (HiggsType .eq. 0)) then
           call VBFNLO_HIGGSTO_pole(xmh,xgh)
        endif
        if (model .eq. 3) then
          call VBFNLO_HIGGSTO_pole(xmh2,xgh2)
        endif
      endif

** Kill the code if we have a very low higgs mass and are not working in the
** complex MSSM
      IF ((XMH.LT.60.d0) .and. (higgsmix .ne. 3)) THEN
         XMH = 60.d0
         print*," Higgs mass below 60 GeV is not valid: m_H = ",xmh
      ENDIF

      XLZ = (XMZ/XMH)**2
      XLW = (XMW/XMH)**2


c     Higgs decays to WW and ZZ
c     subroutines below use  BKOPOUshort common block
c **********  First h -> Z*Z* **********************
      if (replace(bhzz,3) .or. replace(xgh,3)) then

         call hzz_width_flavor_sum(ghzz)
         if (replace(xgh,3) .and. (.not. replace(bhzz,3))) then
            xgh = bhzz*ghzz
         end if
      else
         ghzz = bhzz*xgh
      end if


c*********** Next h -> W*W* *************************
      if (replace(bhww,3) .or. replace(xgh,3)) then
         call calZ(0.0d0,xmin,4)      
         call calZ(xmh**2,xmax,4)
         if (ldebug) print*," xmin,xmax = ",xmin,xmax
        
         relerr = 1d-6
         abserr = 1d-10
         ghww = gaus(hww_width,xmin,xmax,relerr,abserr)   !partial width for h-> V*V* -> all

         if(ghww.lt.0.0d0) ghww = 0.0d0
         if (replace(xgh,3) .and. (.not. replace(bhww,3))) then
            xgh = bhww*ghww
         end if
      else 
         ghww = bhww*xgh
      end if


** Calculating running bottom mass
      if (xmb .gt. 0d0) then
        BMASS_POLE = xmb
      else ! sensible default
        BMASS_POLE = 4.855d0
      endif

c the coefficient fuction c is evaluated within the five flavour 
c approximation for the whole mass range. For m_Higgs > m_top,
c the deviation to the six flavour scheme is less the 1%
      Nf=5 ! TODO: combine with global nf for consistency

      alphas_B=alphas5(BMASS_POLE**2,1)
      alphas_H=alphas5(xmh**2,1)
         
c calculate the running MSbar bottom mass mb(m_b) = BMASS_B
      BMASS_B=BMASS_POLE*(1.0d0-4.d0/3.d0*alphas_B/PI 
     &        +(1.0414d0*Nf-14.3323d0)*alphas_B**2/PI**2)

c calculate the running MSbar bottom mass mb(m_H) = BMASS_H
      alphas_B = alphas5(BMASS_B**2,1)

      C1=(23.0d0/6.0d0*alphas_H/PI)**(12.d0/23.d0)*(1.d0+1.17549d0
     &     *alphas_H/PI+1.50071d0*alphas_H**2/PI**2+0.172478d0
     &     *alphas_H**3/PI**3)
      C2=(23d0/6.d0*alphas_B/PI)**(12.d0/23d0)*(1.d0+1.17549d0
     &     *alphas_B/PI+1.50071d0*alphas_B**2/PI**2+0.172478d0
     &     *alphas_B**3/PI**3)

      BMASS_H=BMASS_B*C1/C2


c*********** Next h -> b bbar ************************
      if (replace(bhbb,3) .or. replace(xgh,3)) then
c Now calculate the H -> bbar width
         DELTA_qq = 5.67d0*alphas_H/PI + (35.94-1.36*Nf)*
     &        alphas_H**2/PI**2 + (164.14-25.77*Nf+0.26*Nf**2)*
     &        alphas_H**3/PI**3

         DELTA_H = alphas_H**2/PI**2 *(1.57-2.0d0/3.0d0*
     -        LOG(xmh**2/xmt**2) + 1.0d0/9.0d0*
     -        (LOG(BMASS_H**2/xmh**2))**2)

         GHBB = 3.0d0*GF/(4.0d0*PI*SQRT(2.0d0)) * XMH * BMASS_H**2
     -        *( 1 + DELTA_qq + DELTA_H )
         if (replace(xgh,3) .and. (.not. replace(bhbb,3))) then
            xgh = bhbb*ghbb
         end if
      else 
         ghbb = bhbb*xgh
      end if  ! end of Higgs to b bar width calculation


c*********** Next h -> c cbar ************************
      if (replace(bhcc,3) .or. replace(xgh,3)) then
         if (xmc .gt. 0d0) then
           CMASS_POLE = xmc
         else ! sensible default
           CMASS_POLE = 1.65d0
         endif
c determine number of flavors with mass smaller than the Higgs mass
c and the running alphas according to that number
         if (xmh.gt.xmt) then
            Nf = 6
         else
            Nf = 5
         endif
         alphas_H = alphas5(xmh**2,1)
         alphas_C = alphas5(CMASS_POLE**2,1)

c calculate the running MSbar charm mass mc(m_c) = CMASS_C
         CMASS_C = CMASS_POLE*(1.0d0-4.d0/3.d0*alphas_C/PI 
     1        + (1.0414*Nf-14.3323)*alphas_C**2/PI**2)

c calculate the running MSbar charm mass mc(m_H) = CMASS_H
         C1=(23.0d0/6.0d0*alphas_H/PI)**(12.d0/23.d0)*(1.d0+1.175d0*
     -        alphas_H/PI+1.501d0*alphas_H**2/PI**2+0.1725d0*
     &        alphas_H**3/PI**3)
         C2=(23d0/6.d0*alphas_C/PI)**(12.d0/23d0)*(1.d0+1.175d0*
     &        alphas_C/PI+1.501d0*alphas_C**2/PI**2+0.1725d0*
     &        alphas_C**3/PI**3)

         CMASS_H = CMASS_C*C1/C2

c Now calculate the H -> ccbar width
         DELTA_qq = 5.67d0*alphas_H/PI + (35.94-1.36*Nf)*
     &        alphas_H**2/PI**2 + (164.14-25.77*Nf+0.26*Nf**2)*
     &        alphas_H**3/PI**3

         DELTA_H = alphas_H**2/PI**2 *(1.57-2.0d0/3.0d0*
     &        LOG(xmh**2/xmt**2)+1.0d0/9.0d0*
     &        (LOG(CMASS_H**2/xmh**2))**2)

         GHCC = 3.0d0*GF/(4.0d0*PI*SQRT(2.0d0)) * XMH * CMASS_H**2
     &        *( 1 + DELTA_qq + DELTA_H )
         if (replace(xgh,3) .and. (.not. replace(bhcc,3))) then
            xgh = bhcc*ghcc
         end if
      else 
         ghcc = bhcc*xgh
      end if  ! end of higgs to c cbar width calculation


c*********** Next h -> t tbar ************************
      if (replace(bhtt,3) .or. replace(xgh,3)) then
         IF ( XMH.GT.2.d0*XMT ) THEN
            BETAT = SQRT(1.d0 - 4.d0*(XMT/XMH)**2)
            beta0 = 7.d0        !11.d0 - 12.d0/3.d0
            QCDT = ( LOG(2.d0*XMT/0.2d0) / 
     &           LOG(XMH/0.2d0) )**(8.d0/beta0)*( 1.d0 + 3.d0*ALFAS/PI )
            GHTT = 3.d0*GF*XMT**2/(4.d0*PI*SQRT(2.d0)) * BETAT**3 * 
     &           XMH * QCDT
         ELSE
            GHTT = zero
         ENDIF
         if (replace(xgh,3) .and. (.not. replace(bhtt,3))) then
            xgh = bhtt*ghtt
         end if
      else 
         ghtt = bhtt*xgh
      end if  ! end of higgs to t tbar 


** Higgs to tau-plus tau-minus decay
      if (replace(bhtau,0) .or. replace(xgh,3)) then
*         XMtau = 1.77684d0      !  and  H ----> tau+ tau-
*         XMtau = ML
         BETAt = SQRT(1.d0 - 4.d0*(XMtau/XMH)**2)
         GHtau = GF*XMtau**2/(4.d0*PI*SQRT(2.d0)) * BETAt**3 * XMH 
         if (replace(xgh,3) .and. (.not. replace(bhtau,3))) then
            xgh = bhtau*ghtau
         end if
      else 
         ghtau = bhtau*xgh
      end if


** Higgs to muons decay      
      if (replace(bhmu,3) .or. replace(xgh,3)) then
*         XMmu = 0.106d0         !  and  H ----> mu+ mu-
         xmmu = MM
         BETAmu = SQRT(1.d0 - 4.d0*(XMmu/XMH)**2)
         GHmu = GF*XMmu**2/(4.d0*PI*SQRT(2.d0)) * BETAmu**3 * XMH 
         if (replace(xgh,3) .and. (.not. replace(bhmu,3))) then
            xgh = bhmu*ghmu
         end if
      else if (xgh .gt. 0d0) then
         ghmu = bhmu*xgh
      end if


** Higgs to gg 
** Leading-order is exact result -- see e.g. Djouadi, hep-ph/0503172
      tau_t = xmh**2/(4d0*xmt**2)
      if (tau_t .le. 1d0) then
        A_12_t = asin(sqrt(tau_t))**2
      else
        A_12_t = -( 
     &    log( (1+sqrt(1-1/tau_t))/(1-sqrt(1-1/tau_t)) ) - dcmplx(0d0,pi)
     &    )**2/4d0
      endif
      A_12_t = 2*(tau_t+(tau_t-1)*A_12_t)/tau_t**2

      tau_b = xmh**2/(4d0*BMASS_POLE**2)
      if (tau_b .le. 1d0) then
        A_12_b = asin(tau_b)**2
      else
        A_12_b = -( 
     &    log( (1+sqrt(1-1/tau_b))/(1-sqrt(1-1/tau_b)) ) - dcmplx(0d0,pi)
     &    )**2/4d0
      endif
      A_12_b = 2*(tau_b+(tau_b-1)*A_12_b)/tau_b**2

      GHgg = GF*alphas_H**2*XMH**3/(36d0*sqrt(2d0)*pi**3)
      GHgg = GHgg*abs(3d0/4d0*(A_12_t+A_12_b))**2

      if (tau_t .le. 1d0) then
c higher-order corrections -- see eg. Schreck, Steinhauser, arXiv:0708.0916
      GHgg = GHgg * (1d0 + 
     &  alphas_H/pi*(17.9167d0+9.574d0*tau_t+5.571d0*tau_t**2
     &              +3.533d0*tau_t**3+2.395d0*tau_t**4+1.708d0*tau_t**4)
     & +(alphas_H/pi)**2*(156.808d0+109.365d0*tau_t+74.434d0*tau_t**2)
     & +(alphas_H/pi)**3*(467.684d0)
     & )
      endif

c for  H --> gamma gamma  I use p. 25 of Higgs Hunter"s Guide,
c using t,b and W loops only
c Physics Letters B 318,1 eq 12 with cpodd d5 operator
c      GHGAM = (alfa*G2/xmw)**2 * (xmh/pi)**3 / 1024.d0
      if (replace(bhgam,3) .or. replace(xgh,3)) then
         ghgam = alfa*g2/(4.0d0*pi*xmw)
         tau_t = 4.d0*(xmt/xmh)**2
         if ( tau_t.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_t)
            eta_m = 1.d0 - sqrt(1.d0 - tau_t)
            f_tau_t = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_t = ( asin(sqrt(1.d0/tau_t)) )**2
         end if
         F_t = -2.d0*tau_t* (1.d0 + (1.d0 - tau_t)*f_tau_t)

         tau_b = 4.d0*(BMASS_POLE/xmh)**2
         if ( tau_b.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_b)
            eta_m = 1.d0 - sqrt(1.d0 - tau_b)
            f_tau_b = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_b = ( asin(sqrt(1.d0/tau_b)) )**2
         end if
         F_b = -2.d0*tau_b* (1.d0 + (1.d0 - tau_b)*f_tau_b)
         
         tau_W = 4.d0*(xmw/xmh)**2
         if ( tau_W.lt.1.d0 ) then
            eta_p = 1.d0 + sqrt(1.d0 - tau_W)
            eta_m = 1.d0 - sqrt(1.d0 - tau_W)
            f_tau_W = -0.25d0 * ( log(eta_p/eta_m) + pi*(0d0,-1.d0) )**2
         else
            f_tau_W = ( asin(sqrt(1.d0/tau_W)) )**2
         end if
         F_W = 2.d0 + 3.d0*tau_W*(1.d0 + (2.d0 - tau_W)*f_tau_W)
      
         dum = 3.d0*(4.d0*F_t + F_b)/9.d0 + F_W*sinba
         
         dum = ghgam*dum
         ahvv(2,1,1) = dum

         if(with_anomHiggs) call fillanomhcoupl(2,A_sum,dum,loopfac)

         ghgam = xmh**3/(64.0d0*pi) * 
     &        ( dble(ahvv(2,1,1))**2+dimag(ahvv(2,1,1))**2+
     &        dble(ahvv(3,1,1))**2+dimag(ahvv(3,1,1))**2 )
c      GHGAM = GHGAM * ( dble(dum)**2 + dimag(dum)**2 ) 
         if (replace(xgh,3) .and. (.not. replace(bhgam,3))) then
            xgh = bhgam*ghgam
         end if
      else if (xgh .gt. 0d0) then
         ghgam = bhgam*xgh
      end if


c for H -> gamma Z (on-shell)
c cp even and odd d5 operators
c     Using Higgs Hunters guide
c     t,b ,W loops are include in SM calculation
c     compute At,Ab,Aw
      if (replace(bhgamz,3) .or. replace(xgh,3)) then
         cos2w = 1.0d0 - sin2w
         tan2w = sin2w/cos2w
         cotw = sqrt(1.0d0/tan2w)
c     
c     b quark loop
         tau_b = 4.d0*(BMASS_POLE/xmh)**2
         lamb_b = 4.0d0*(BMASS_POLE/xmz)**2
         call compints(tau_b,lamb_b,II_1,II_2)
c   
         A_b = -3.0d0 *2.0d0*V(4,1)*V(4,2)*(2.0d0*c/(g2*E))*(II_1 -II_2)
         A_b = A_b/(s*c) 
c     top quark loop
         tau_t = 4.d0*(xmt/xmh)**2
         lamb_t = 4.0d0*(xmt/xmz)**2
         call compints(tau_t,lamb_t,II_1,II_2)
c   
         A_t = -3.0d0 *2.0d0*V(3,1)*V(3,2)*(2.0d0*c/(g2*E))*(II_1 -II_2)
         A_t = A_t/(s*c) 
c     W loop
         tau_w = 4.0d0*(xmw/xmh)**2
         lamb_w = 4.0d0*(xmw/xmz)**2
         call compints(tau_w,lamb_w,II_1,II_2)
         A_w = 4.0d0*(3.0d0-tan2w)*II_2
         A_w = A_w + ((1.0d0+2.0d0/tau_w)*tan2w
     $        -(5.0d0+2.0d0/tau_w))*II_1
         A_w = -cotw*A_w*sinba
         A_sum = A_w + A_t + A_b
         A_sum = A_sum * alfa*g2/(4.0d0*pi*xmw)
         ahvv(2,1,2) = -A_sum 
         ahvv(2,2,1) = -A_sum

         if(with_anomHiggs) call fillanomhcoupl(2,A_sum,dum,loopfac)
c     
         if(xmh.gt.xmz) then
            call hzgamma_width(ghgamz)
         else
            ghgamz = 0.0d0
         endif

         if (replace(xgh,3) .and. (.not. replace(bhgamz,3))) then
            xgh = bhgamz*ghgamz
         end if
      else 
         ghgamz = bhgamz*xgh
      end if  ! end of Higgs decay to gamma Z

c H total decay width and branching ratios
      if (replace(xgh,3)) then
         XGH = GHWW + GHZZ + GHBB + GHTT + GHCC + GHTAU + GHMU 
     $        + GHGG + GHGAM + GHGAMZ
      end if
c
      BHWW  = GHWW/XGH
      BHZZ  = GHZZ/XGH
      BHGG  = GHGG/XGH
      BHTT  = GHTT/XGH
      BHBB  = GHBB/XGH
      BHTAU = GHTAU/XGH
      BHMU  = GHMU/XGH
      BHCC  = GHCC/XGH
      BHGAM = GHGAM/XGH
      BHGAMZ = GHGAMZ/XGH      

** Checking that the relevant branching ratio is less than one
      call checkBR


** top width: this is calculated if we're not using a slha file or if it's not 
** present in the slha file
      if (replace(gamt,0)) then      
         gamt = gf*xmt**3/(8d0*pi*sqrt(2d0))*(1-(xmw/xmt)**2)**2 *
     -        (1+2d0*(xmw/xmt)**2)*
     -        (1 - 2d0*alphas5(xmt**2,1)/(3d0*pi)*(2*pi**2/3d0 - 2.5d0))
      end if


c width info for higgs
c
      XMG(6) = XMH*XGH 
c
c dummy xmg for commbon block BKOPOUshort for higgs
c
      xmgs(6) = xmg(6)      
c

*** Initialising scales
      call InitScales

      if (lglobalprint) then

      print *," "
      print *,"  Information on couplings, masses and widths"
      print *,"==============================================="
      print *," "
      print *,"  SM coupling parameters "
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F14.12)") "G_Fermi",": ",GF
      write(*,"(T4,A,T36,A,F10.8)") "SIN**2 of the weak mixing angle",
     &     ": ",SIN2W
      write(*,"(T4,A,T36,A,G13.8)") "Q.E.D. fine structure constant",
     &     ": 1/",1.D0/ALFA
      write(*,"(T4,A,T36,A,F10.8)") "SU(2) coupling constant g",
     &     ": ",e/s
      write(*,"(T4,A,T36,A,F10.8)") "U(1) coupling constant g'",
     &     ": ",e/c

      write(*,"(T4,A,T36,A,F10.8)") "Q.C.D. LO alfas(MZ)",
     &     ": ",ALFAS_LO
      write(*,"(T4,A,T36,A,F10.8)") "Q.C.D. NLO alfas(MZ)",
     &     ": ",ALFAS 
c output, which is not so relevant to the user
      if (.false.) then
         write(*,*) " "
         write(*,12) C,S,Z,W,Q
   12 FORMAT(/" C,S,Z,W,Q =",5F10.6)
      write(6,*) " "
      write(6,13) ((V(J,K),K=1,5),J=1,4)
      write(6,14) ((A(J,K),K=1,5),J=1,4)
      write(6,15) ((CLR(J,K,-1),K=1,5),J=1,4)
      write(6,16) ((CLR(J,K,1),K=1,5),J=1,4)
      write(6,17) (((B(I,J,K),K=1,6),J=1,6),I=1,6)
   13 FORMAT(" FFB VECTOR COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   14 FORMAT(" FFB AXIAL-VECTOR COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   15 FORMAT(" FFB LEFTHANDED COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   16 FORMAT(" FFB RIGHTHANDED COUPLING ( ROW=F, COL= B )"/,4(5F10.6/))
   17 FORMAT(" BBB NON-ABELIAN COUPLINGS"/,6( 6(6F10.6/)/) )
      write(6,*) " "
      write(6,19) (XM2(J),XMG(J),J=1,6)
   19 FORMAT(/," MASS**2 AND MASS*WIDTH VALUES"/,6(2F10.2/))
      endif

      print *," "
      print *,"     particle masses"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F9.4)") "mass of the Z0 boson",": ",XMZ
      write(*,"(T4,A,T36,A,F9.4)") "mass of the W+/- boson",": ",XMW
      write(*,"(T4,A,T36,A,F9.4)") "mass of the Higgs boson",": ",XMH
      write(*,"(T4,A,T36,A,F9.4)") "mass of the top quark",": ",XMT
      write(*,"(T4,A,T36,A,F9.4)") "mass of the bottom quark at MH",
     &     ": ",BMASS_H
      write(*,"(T3,A)") "(All other fermions have mass equal ZERO)"
      IF (INFO.LE.1) RETURN
      print *," "
      print *,"     particle widths"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Z0",": ",
     &     XGZ," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of W+ or W-",": ",
     &     XGW," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of Higgs",": ",
     &     XGH," GeV"
      if (model.eq.3)
     &  write(*,"(T4,A,T36,A,F10.6,A)") "total widths of second Higgs",": ",
     &       XGH2," GeV"
      write(*,"(T4,A,T36,A,F10.6,A)") "total widths of top quark",": ",
     &     gamt," GeV"

      IF (with_anomHiggs) THEN
         print*," "         
         PRINT*," partial Higgs decay width "
         print*,"---------------------------"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> WW",": ",ghww," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> ZZ",": ",ghzz," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> bottom anti-bottom",
     &        ": ", ghbb," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> top anti-top",
     &        ": ",ghtt," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> charm anti-charm",
     &        ": ",  ghcc," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> tau anti-tau",": ",
     &        ghtau," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> mu anti-mu",": ",
     &        ghmu," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> gluon gluon",": ",
     &        ghgg," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> gamma gamma",": ",
     &        ghgam," GeV"
         write(*,"(T4,A,T36,A,D12.6,A)") "H -> Z gamma",": ",
     &        ghgamz," GeV"
      ENDIF

      print *," "
      print *,"    branching ratios"
      print*,"---------------------------"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> lepton anti-lepton",": ",
     &     BWNE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "W -> massless quark anti-quark",": ", BWUD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "W -> top bottom",": ",
     &     BWTB*100.d0," %"
      print *," "
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> neutrino anti-neutrino",": ",BZNN*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> charged lepton anti-lepton",": ",BZEE*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> up-type quark anti-quark",": ",BZUU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") 
     &     "Z -> down-type quark anti-quark",": ",BZDD*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "Z -> top anti-top",": ",
     &     BZTT*100.d0," %"
      print *," "
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> WW",": ",bhww*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> ZZ",": ",bhzz*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> bottom anti-bottom",": ",
     &     BHBB*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> top anti-top",": ",
     &     BHTT*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> charm anti-charm",": ",
     &     BHcc*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> tau anti-tau",": ",
     &     BHTAU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> mu anti-mu",": ",
     &     BHMU*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> gluon gluon",": ",
     &     BHGG*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> gamma gamma",": ",
     &     BHGAM*100.d0," %"
      write(*,"(T4,A,T36,A,F10.7,A)") "H -> Z gamma",": ",
     &     BHGAMZ*100.d0," %"
      print *," "
      if (model.eq.3) then
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> WW",": ",bh2ww*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> ZZ",": ",bh2zz*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> bottom anti-bottom",": ",
     &       BH2BB*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> top anti-top",": ",
     &       BH2TT*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> charm anti-charm",": ",
     &       BH2cc*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> tau anti-tau",": ",
     &       BH2TAU*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> mu anti-mu",": ",
     &       BH2MU*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> gluon gluon",": ",
     &       BH2GG*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> gamma gamma",": ",
     &       BH2GAM*100.d0," %"
        write(*,"(T4,A,T36,A,F10.7,A)") "H2 -> Z gamma",": ",
     &       BH2GAMZ*100.d0," %"
        print *," "
      endif

      endif ! lprint

      if (ldebug) stop

      RETURN

      END
      

********************************************************************************
********************************************************************************
c                         
c     The subroutine below calculate the partial width for H-> V*V* semi-
c     analytically.
      
      function hww_width_amp(xx,yy)
      implicit none

      double precision lamb1,lamb2,xx,yy,hww_width_amp,
     $     fun5de,fun5do,funsm5de,funsm
c     
#include "global.inc"
#include "mssm.inc"

      double precision g5hvv(2,1:4,1:4),lambda5
      common/hcoupl/g5hvv,lambda5
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision xmv, xmh
      double precision lambda,constants
      double precision Q2(2)
      integer id,bos,bos1,bos2
      double complex ahvv(3,4,4), ahvvL(3,4,4)
      common/tensorhvv/ ahvv, ahvvL

      ! this is no longer used for z bosons
      bos1 = 3
      bos2 = 4

      xmv = sqrt(xm2s(bos1))
      xmh = sqrt(xm2s(6))
c      print*,"mh=",xmh
      call calQ2(xx,Q2(1),4)
c      Q2(1) = xx 
      lamb1 = Q2(1)/xmh**2
c      
      call calQ2(yy,Q2(2),4)
c      Q2(2) = yy
      lamb2 = Q2(2)/xmh**2
c
      call lambda1(1.0d0,lamb1,lamb2,lambda) !call 1 time for all
      funsm = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1       10.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     
      constants = dble(ahvv(1,bos1,bos2))**2 + 
     1     dimag(ahvv(1,bos1,bos2))**2
      constants = constants * xmh**3/((PI**3) * 64.0d0 * xmv**4)

      funsm = funsm * constants

  
c     cp odd contribution
      fun5do = lambda * (1.0d0 + lamb1**2 + lamb2**2 - 
     1     2.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
      constants = g5hvv(2,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5do = fun5do * constants
      
c     cp even contribution 
      fun5de = lambda * (1.0d0 + lamb1**2 + lamb2**2 + 
     1     4.0d0*lamb1*lamb2 -  2.0d0 * lamb1 -  2.0d0 * lamb2)
     2     *lamb1*lamb2
      constants = g5hvv(1,bos1,bos2)**2 * 
     $     xmh**7/( 8.0d0 * xmv**4 * (PI**3) * (lambda5**2))
      fun5de = fun5de*constants

c     interefence of SM and CP even 
      funsm5de = lambda * (1.0d0 - lamb1 - lamb2) 
     2     *lamb1*lamb2
c
c id = 1 for photon, 2 for Z , 3 for W , 4 for W
c
      constants = -g5hvv(1,bos1,bos2)*dble(ahvv(1,bos1,bos2))* 
     $     3.0d0*xmh**5/( 8.0d0 * xmv**4 * (PI**3) * lambda5)
      funsm5de = funsm5de * constants


      hww_width_amp = funsm + funsm5de + fun5do + fun5de !sum 
      end


      subroutine calQ2(x,Q2,id)
      implicit none
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision x,Q2
      integer id
c      for W bosons
      Q2 = xm2s(id) + xmgs(id) * dtan(x)
      end


      subroutine calZ(QQ2,Z,id)
      implicit none
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs
      double precision QQ2,Z
      integer id

c      print*,"mv=",sqrt(xm2s(id))
      Z = datan((QQ2 - xm2s(id))/xmgs(id))
      end


********************************************************************************
********************************************************************************

      subroutine lambda1(x,y,z,lambda)
      implicit none
      double precision x,y,z,lambda2, lambda

      lambda2 = x**2 + y**2 + z**2 - 2.0d0*x*y - 2.0d0*y*z - 2.0d0*x*z
      if (lambda2.gt.0d0) then
         lambda = sqrt(lambda2)
      else
         lambda = 0d0
      endif
      end


      function hww_width_1d(y)
      implicit none
      double precision hww_width_1d,y,hww_width_amp,x
      common/x_ghww/ x

      hww_width_1d = hww_width_amp(x,y)
      end


      function hww_width(xx)
      implicit none
      double precision hww_width
      double precision xx,y1,y2,x,hww_width_1d
      common /x_ghww/ x
      double precision gaus
      external hww_width_1d
      external gaus
      double precision Q2,q2max
      double precision  xmh
      double precision XM2s(6),XMGs(6)
      COMMON /BKOPOUshort/ XM2s,XMGs

      double precision relerr, abserr

      xmh = sqrt(xm2s(6))
      x = xx
      call calQ2(x,Q2,4)
      q2max = (xmh - sqrt(Q2))**2
      call calZ(q2max,y2,4) 
      call calZ(0.0d0,y1,4)

      relerr = 1d-6
      abserr = 1d-10

      hww_width = gaus(hww_width_1d,y1,y2,relerr,abserr)
      ! this routine is no longer used for z bosons
      ! if(id.eq.2) then          !symmetry factor of h->ZZ
      !    h =ss/2.0d0
      ! else
      !    h = ss
      ! endif

      end


! this is a numerical integration function \int_a^b f
! it's based on cernlibs gaus(f,a,b,eeps) and has some overlap with
! the quadpack routines
!
! eeps > 0: absolute precision
! eeps < 0: relative precision
      recursive double precision function gaus(f,a,b,relerr,abserr)
      implicit none

      double precision, intent(in) :: a, b
      double precision, intent(in) :: relerr, abserr
      double precision :: AA, BB, C1, C2, DELTA
      double precision S8, S16, U, Y
      integer i
      integer count

      double precision f
      external f

      double precision, parameter :: numacc = 1.0e-19
      double precision, parameter, dimension(12) :: W = (
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/),
     & X = (
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/)

      count = 0

      delta=numacc*abs(A-B)
      gaus=0.d0
      AA=A
      Y=B-A
      do while (abs(y)>= delta)
          s16=0
          s8=1e20
          do while (abs(S16-S8) > max(relerr*abs(S16), abserr))
              BB=AA+Y
              count = count + 1
              if (count .ge. 150000) then
                  write(*,*)'Sorry!  We are having problems calculating a '
                  write(*,*)'decay width.  Try altering the anomalous couplings'
                  write(*,*)'or putting the width in by-hand using a SLHA file'
                  stop
              end if

              C1=0.5d0*(AA+BB)
              C2=C1-AA
              S8=0.d0
              DO I=1,4
                  S8=S8+W(I)*(F(C1+x(i)*c2)+F(C1-x(i)*c2))
              end do
              S16=0.d0
              DO I=5,12
                  U=X(I)*C2
                  S16=S16+W(I)*(F(C1+U)+F(C1-U))
              end do
              S8=S8*C2
              S16=S16*C2

              Y=0.5d0*Y
              if (abs(y) <= delta) return ! this should not happen for normal integrands
          enddo
          GAUS=GAUS+S16
          AA=BB
          Y=B-AA
      enddo

      END


********************************************************************************
********************************************************************************
C This subroutine computes the parametric integrals in the 
c Higgs Hunters Guide on pg 29 (eq. 2.24)

      subroutine compints(a,b,I1,I2)
      implicit none
      double precision a,b
      double complex I1,I2
      double complex fci,g
      external fci
      external g

      I1 = a*b/(2.0d0*(a-b)) + a**2 * b**2/(2.0*(a-b)**2) * (fci(a)-fci(b))
     $     + a**2 * b/((a-b)**2) *(g(a) - g(b))
      I2 = -a*b/(2.0d0*(a -b))* (fci(a) - fci(b))
      return

      end


********************************************************************************
********************************************************************************
c
** Function g from eq. 2.24 of the Higgs Hunters' Guide (p29), used in the
** calculation of H -> Z photon

      double complex function g(aa)
      implicit none
      double precision aa,nup,num
      double precision pi
      parameter(pi =3.141592653589793d0) 

      if(aa.ge.1.0d0) then
         g = sqrt(aa -1.0d0)*ASIN(1.0d0/sqrt(aa))
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         g = 0.5d0* sqrt(1.0d0 - aa)*(log(nup/num)*(1.d0,0.d0)
     $        + pi*(0.d0,-1.d0))
      endif
      return

      end


********************************************************************************
********************************************************************************
c
** Function from eq. 2.19 of Higgs Hunters' Guide (p26), used in calculation of
** H -> Z photon

      double complex function fci(aa)
      implicit none
      double precision aa,nup,num
      double precision pi
      parameter(pi =3.141592653589793d0) 
      
      if(aa.ge.1.0d0) then
         fci = (ASIN(sqrt(1.0d0/aa)))**2
      
      else
         nup = 1.0d0 + sqrt(1.0d0 - aa)
         num = 1.0d0 - sqrt(1.0d0 - aa)
         fci = -0.25d0*(log(nup/num)
     $        + pi*(0.d0,-1.d0))**2
      endif
      return

      end


********************************************************************************
********************************************************************************
c	all couplings, widths, and masses are computed from 
c	values stored in common/ bkopout/
c
c	edited by Barbara Jager on March 10, 2005

	subroutine ctrans(xmb)
          use globalvars, only: ldoblha
c     some coupling parameters are modified here. (This should be corrected.)
c***********************************************************************
c 	This subroutine sets up the couplings, masses and widths 
c	needed by HELAS as they are given in koppln.f
c***********************************************************************

	implicit none      
 
c common blocks which need to be filled for HELAS (masses, widths, couplings):	
#include "coupl.inc"
#include "scales.inc"
#include "global.inc"
#include "BLHAhelper.inc"
			
      double precision g2,ee2,sw,cw,vv,sc2
      logical loutput
      parameter (loutput = .false.)
	
c     for linit=0, alfas is filled according to the input value
c     for linit=1, alfas is filled according to scales.inc input, als(1,1).
      integer linit
      data linit /0/
      save linit
c common block with input values (only needed for top mass xmt):

      double precision alphas, xmt, alfa, xmz, xmw,s2w, xmh, gf, alfas_lo
      common /bkopin/ alphas, xmt, alfa, xmz, xmw, s2w, xmh, gf, alfas_lo
      double precision xmb
	
c common block which has been filled in koppln.f and contains 
c FFB, BBB, and Higgs couplings (but not in form needed by HELAS):
      double complex clrCT
      common /bkopouCT/ clrCT(3:4,2:4,-1:1)
      DOUBLE PRECISION CLR,XM2,XMG,B,V,A
      COMMON /BKOPOU/ CLR(4,5,-1:1),XM2(6),XMG(6),B(6,6,6),V(4,5),A(4,5)
      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XMH2, XGH2, 
     &                 sinba, cosba, ch2ww, ch2zz
      

      if ((linit.eq.0).or.ldoblha) then
c fill basic parameters:
	G      =  B(5,5,5)
	GG(1)  = -G
	GG(2)  = -G
	alfas  =  G**2/(pi*4.d0)  ! strong coupling at fixed ren. scale MZ
        linit = linit+1
      else if (linit.eq.1) then
	G      =  sqrt(als(1,1)*4*pi)
	GG(1)  = -G
	GG(2)  = -G
	alfas  =  G**2/(pi*4.d0)  ! strong coupling at selected ren. scale

        return
      endif

	ee     = -CLR(2,1,-1)
	ee2    =  ee**2
	alpha  =  ee2/(pi*4.d0)
	
	g2     = sqrt(2.d0)*CLR(2,3,-1)
	cw     = B(2,4,3)/g2
	sin2w  = 1.d0 - cw**2
	sw     = sqrt(sin2w)
	
	
c fill masses:

	hmass = sqrt(xm2(6))
	wmass = sqrt(xm2(3))
	zmass = sqrt(xm2(2))
	amass = 0.d0
	tmass = xmt
c	tmass = 175.d0
        bmass = xmb
c	bmass = 0.d0
c	bmass = 4.5d0	! value given in qqhqq.f
c	bmass = 4.7d0	! value from switchmom.f here 	 
	cmass = 0.d0
	lmass = 0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!        bmass = 0.0
	
c  fill widths:			

	hwidth = xmg(6)/hmass 
	wwidth = xmg(3)/wmass
	zwidth = xmg(2)/zmass
!        print*, "wwidth",wwidth
!        print*, "zwidth",zwidth
	awidth = 0.d0	
	lwidth =  0d0
	twidth = 1.6d0	! no top emerges in WW-prod. -> value here "arbitrary"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!	wwidth = 2.06
!        wwidth = 2.10457
!        twidth = 1.6511d0
        
	
c  auxiliary variables:	
	
	vv  = 2.d0*wmass*sw/ee
	sc2 = sin2w*(1.d0-sin2w)
	
c  fill gauge couplings: 
c
c	fermion-fermion-vector couplings: 

	gal(1) = -dcmplx(clr(2,1,-1))	! photon
	gal(2) = -dcmplx(clr(2,1, 1))
	gau(1) = -dcmplx(clr(3,1,-1))	
	gau(2) = -dcmplx(clr(3,1, 1))	
	gad(1) = -dcmplx(clr(4,1,-1))	
	gad(2) = -dcmplx(clr(4,1, 1))	
	
	gwf(1) = -dcmplx(clr(1,3,-1))	! W+- boson (same for all fermions)
	gwf(2) = -dcmplx(clr(1,3, 1))	! zero
	
	gzn(1) = -dcmplx(clr(1,2,-1))	! Z boson
	gzn(2) = -dcmplx(clr(1,2, 1))	! zero	
	gzl(1) = -dcmplx(clr(2,2,-1))
	gzl(2) = -dcmplx(clr(2,2, 1))	
	gzu(1) = -dcmplx(clr(3,2,-1))
	gzu(2) = -dcmplx(clr(3,2, 1))	
	gzd(1) = -dcmplx(clr(4,2,-1))
	gzd(2) = -dcmplx(clr(4,2, 1))
	
c  	vector boson couplings:

	gw   = ee/sw
	gwwa = -dcmplx(b(3,4,1))
	gwwz = -dcmplx(b(3,4,2))

!         print*, gwwz, gwwa, gw, wmass
!         STOP
	
c	gauge - higgs boson couplings:

	gwwh  = dcmplx(b(3,4,6)*wmass)	
	gzzh  = dcmplx(b(2,2,6)*wmass)	
	ghhh  = dcmplx(-hmass**2/vv*3d0)
	gwwhh = dcmplx( ee2/sin2w*0.5d0 )
        gzzhh = dcmplx( ee2/sc2*0.5d0)
        ghhhh = ghhh/vv
        if (with_secondhiggs) then
	h2mass   = xmh2
	h2width  = xgh2
	gwwh2  = ch2ww
	gzzh2  = ch2zz
	gsba   = sinba
	gcba   = cosba
        endif
	

c-----------------------------------------------

c for reference:
c	info on couplings returned by ctrans:

      if (loutput) then
	write(6,*) "couplings returned by ctrans:"
	write(6,*) "alpha,ee, sin2w, alfas,g,gg:"
	write(6,*)  alpha,ee, sin2w, alfas,g,gg
	write(6,*) 
	write(6,*) "hmass, wmass, zmass, amass:"
	write(6,*)  hmass, wmass, zmass, amass
	write(6,*) 
	write(6,*) "tmass, bmass, lmass, cmass:"
	write(6,*)  tmass, bmass, lmass, cmass
	write(6,*) 
	write(6,*) "hwidth, wwidth, zwidth:" 
	write(6,*)  hwidth, wwidth, zwidth
	write(6,*) 
	write(6,*) "twidth, lwidth, awidth:"
	write(6,*)  twidth, lwidth, awidth
	write(6,*) 
	write(6,*) "gal   , gad   , gau   , gwf:"
	write(6,*)  gal   , gad   , gau   , gwf
	write(6,*) 
	write(6,*) "gzn   , gzl   , gzd   , gzu:"
	write(6,*)  gzn   , gzl   , gzd   , gzu
	write(6,*) 
	write(6,*) "gw, gwwa, gwwz:"
	write(6,*)  gw, gwwa, gwwz
	write(6,*) 
	write(6,*) "gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh:"
	write(6,*)  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh
	write(6,*) 
      endif
	
c--------------------------------------------------------------------------	

	return
	end


********************************************************************************
********************************************************************************
********************************************************************************

** sophy@particle.uni-karlsruhe.de
** This subroutine is going to calculate the counterterm coupling
** for the qqV vertex as long as 'sector' is less than 3 (i.e. only fermions
** and sfermions are included in the electroweak loops)

      subroutine qqV_ct(clrCT,model,sector)

      implicit none

**  The counterterm couplings qqV      
      double complex clrCT(3:4,2:4,-1:1)
      double complex fcplCTfm(4,6)
      double complex clrCTfm(3:4,2:4,-1:1)

      common /FMclr/ clrCTfm, fcplCTfm

      integer model,sector

#include "mssm.inc"   


** FERMION counterterm couplings!
** First the coupling udW: right handed
      clrCTfm(3,3,1) = 0D0

** coupling udW: left handed
      clrCTfm(3,3,-1) = (dZet1fm - (dSW1fm/SW) + (dZW1fm/2.D0)) 
      clrCTfm(3,3,-1) = clrCTfm(3,3,-1)*EL/(SQRT(2.D0)*SW)

      if (model .eq. 2) then         
         clrCTfm(3,3,-1) = (dZet1fm + (dSW1fm/SW) + (dZW1fm/2.D0)) 
         clrCTfm(3,3,-1) = -clrCTfm(3,3,-1)*EL/(SQRT(2.D0)*SW)
      end if

      clrCTfm(4,4,-1) = clrCTfm(3,3,-1)
      clrCTfm(4,4,1) = clrCTfm(3,3,1)
      

** Now the coupling uuZ: left handed
      clrCTfm(3,2,-1) = EL*((gL3*dZZZ1fm/2.D0)+dgL3fm-(dZAZ1fm/3.D0))

** coupling uuZ: right handed
      clrCTfm(3,2,1) = EL*((gR3*dZZZ1fm/2.D0)+dgR3fm-(dZAZ1fm/3.D0))


** Now the coupling ddZ: left handed
      clrCTfm(4,2,-1) = EL*((gL4*dZZZ1fm/2.D0)+dgL4fm+(dZAZ1fm/6.D0))

** coupling ddZ: right handed
      clrCTfm(4,2,1) = EL*((gR4*dZZZ1fm/2.D0)+dgR4fm+(dZAZ1fm/6.D0))

      
      write(*,*) '   '
      write(*,*) 'Counterterm qqV couplings FOR FERMIONS ...'
      write(*,*) 'udWR=', clrCTfm(3,3,1)
      write(*,*) 'udWL=', clrCTfm(3,3,-1)
      
      write(*,*) 'uuZR=', clrCTfm(3,2,1)
      write(*,*) 'uuZL=', clrCTfm(3,2,-1)
      write(*,*) 'ddZR=', clrCTfm(4,2,1)
      write(*,*) 'ddZL=', clrCTfm(4,2,-1)


      clrCT(3,3,1) = clrCTfm(3,3,1)
      clrCT(3,3,-1) = clrCTfm(3,3,-1) 
      clrCT(4,4,-1) = clrCT(3,3,-1)
      clrCT(4,4,1) = clrCT(3,3,1)
      clrCT(3,2,-1) = clrCTfm(3,2,-1)
      clrCT(3,2,1) = clrCTfm(3,2,1)
      clrCT(4,2,-1) = clrCTfm(4,2,-1)
      clrCT(4,2,1) = clrCTfm(4,2,1)


      end


********************************************************************************
c
c PHASE SPACE REDUCTION FACTOR FOR HEAVY FERMIONIC DECAYS
c
      double precision function heavy(X,Y1,Y2)

      implicit none

      double precision X, Y1, Y2
      
      HEAVY = ( 1.d0 - .5d0*(Y1**2 + Y2**2) - 
     &           0.5d0*(Y1**2 - Y2**2)**2 +
     &           3.d0*Y1*Y2*((X**2 - 1.d0)/(X**2 + 1.d0)) )*
     &           SQRT((1.d0 - Y1**2 - Y2**2)**2 - 4.d0*Y1**2*Y2**2)

      end 


********************************************************************************

