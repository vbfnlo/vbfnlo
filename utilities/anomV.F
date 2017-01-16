********************************************************************************
********************************************************************************
***                                                                          ***
*** anomV.F                                                                  ***
*** 1 June 2012                                                              ***
*** sophy@particle.uni-karlsruhe                                             ***
*** Last modified by Oliver Schlimpert: June 2012                            ***
***                                                                          ***
*** Legacy NOTE: this is a consolidation of the 4 previous anomalous boson   ***
*** couplings read-in routines                                               ***
***                                                                          ***
********************************************************************************
********************************************************************************

      subroutine read_anomVcouplings
        use globalvars, only: lglobalprint
        use readinput
c reads in initial set of anomalous coupling constants 
c from file anomV.dat. The coupling constants
c are stored in common blocks in anomV_cplg.inc      

      implicit none

** Dummy parameter used to check consistency of Delta kappa_Z and _photon
      double precision consistency, massscale2

** Dummy parameters for filling HVV couplings
      double precision treefacW, treefacZ, loopfac, dum
      common /lhcoup/ treefacW, treefacZ, loopfac

#include "an_couplings.inc"   
#include "mssm.inc"
#include "global.inc"
#include "process.inc"


      if (lglobalprint) then
        print *," "
        print *,"  Information on anomalous coupling parameters  "
        print *,"------------------------------------------------"
        call loadfile("anomV.dat",.true.)
      else
        call loadfile("anomV.dat",.false.)
      endif

** Read switches controlling whether we use a formfactor:
      call read_logical("FORMFAC",formfact,.false.)
      call read_logical("KMATRIX",kmatrix,.false.)
      call read_logical("FORMFAC_IND",formfac_ind,.false.)
** universal formfactors:
      if (formfact) then
         call read_real("FFMASSSCALE",ffmassscale2,2000d0)
         if (ffmassscale2 .lt. 1d-3) then
            write(*,*)'The input value of FFMASSSCALE is too small'
            write(*,*)"We will use FFMASSSCALE = 2000 GeV instead"
            ffmassscale2 = 2000d0
         end if            
         ffmassscale2 = ffmassscale2**2
         call read_int("FFEXP",ffexponent,2)
         if (ffexponent .lt. 0) then
            write(*,*)"You cannot use a negative exponent for the"
            write(*,*)"formfactor.  We will instead use"
            write(*,*)"FFEXP = 2"
            ffexponent = 2
         end if
         if (.not. formfac_ind) then
            massscale2FB = ffmassscale2
            massscale2FW = ffmassscale2
            massscale2FWWW = ffmassscale2
            ffexpFWWW = ffexponent
            ffexpFW = ffexponent
            ffexpFB = ffexponent

            massscale2L = ffmassscale2
            massscale2KZ = ffmassscale2
            massscale2KA = ffmassscale2
            massscale2G = ffmassscale2
            ffexpL = ffexponent
            ffexpKZ = ffexponent
            ffexpKA = ffexponent
            ffexpG = ffexponent
         end if
      end if


** Now we move on to the dim-6 operators that parametrise the WWZ and WWgamma
** couplings:
      call read_int("TRIANOM",trianom,1)

      if (trianom .eq. 1) then

         call read_real("FWWW",fwww_0,0d0)
         call read_real("FW",fw_0,0d0)
         call read_real("FB",fb_0,0d0)

         if (formfac_ind .and. formfact) then
            call read_real("MASS_SCALE_FWWW",massscale2FWWW,2000D0)
            massscale2FWWW = massscale2FWWW**2
            if (massscale2FWWW .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FWWW is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FWWW = ffmassscale2
               write(*,*)"MASS_SCALE_FWWW =", sqrt(massscale2FWWW)
            end if            
            call read_int("FFEXP_FWWW",ffexpFWWW,2)
            if (ffexpFWWW .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFWWW = ffexponent 
               write(*,*)"FFEXP_FWWW =", ffexpFWWW
            end if

            call read_real("MASS_SCALE_FW",massscale2FW,2000D0)
            massscale2FW = massscale2FW**2
            if (massscale2FW .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FW is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FW = ffmassscale2
               write(*,*)"MASS_SCALE_FW =", sqrt(massscale2FW)
            end if            
            call read_int("FFEXP_FW",ffexpFW,2) 
            if (ffexpFW .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFW = ffexponent 
               write(*,*)"FFEXP_FW =", ffexpFW
            end if

            call read_real("MASS_SCALE_FB",massscale2FB,2000D0)
            massscale2FB = massscale2FB**2
            if (massscale2FB .lt. 1d-3) then
               write(*,*)'The input value MASS_SCALE_FB is too small'
               write(*,*)"We will use the universal formfactor scale."
               massscale2FB = ffmassscale2
               write(*,*)"MASS_SCALE_FB =", sqrt(massscale2FB)
            end if            
            call read_int("FFEXP_FB",ffexpFB,2)
            if (ffexpFB .lt. 0) then
               write(*,*)"You cannot use a negative exponent for the"
               write(*,*)"formfactor.  We will instead use value from"
               write(*,*)"the universal formfactor"
               ffexpFB = ffexponent 
               write(*,*)"FFEXP_FB =", ffexpFB
            end if
         end if

** setting trianom=2 parametrisation
         aDkappa0_0 = 0.5d0*MW2*(fb_0 + fw_0)
         zDkappa0_0 = 0.5d0*MZ2*(CW2*fw_0 - SW2*fb_0)
         lambda0_0 = 3d0*MW2*EL*EL*fwww_0/(2d0*SW2)
         zDg0_0 = 0.5d0*MZ2*fw_0
** formfactors: note formfactor for kappa0(z/gamma) is not a simple conversion 
**              and is therefore done in the formfactor routine
         massscale2L = massscale2FWWW
         ffexpL = ffexpFWWW
         massscale2G = massscale2FW
         ffexpG = ffexpFW


      else if (trianom .eq. 2) then

         call read_real("LAMBDA0",lambda0_0,0D0)
         call read_real("ZDELTAG1",zDg0_0,0D0)
         call read_real("ZDELTAKAPPA0",zDkappa0_0,0D0)
         call read_real("ADELTAKAPPA0",aDkappa0_0,0D0)

** If zdeltag, zdeltakappa0 or adeltakappa0 is zero, we set it to be consistent
         if ((zDkappa0_0.eq.0D0) .and. 
     -        ((zDg0_0.ne.0D0) .or. (aDkappa0_0.ne.0D0)))then
            zDkappa0_0 = zDg0_0-SW2/CW2*aDkappa0_0
            write(*,*)'Note!' 
            write(*,*)'ZDELTAKAPPA0 = ZDELTAG0-SW2/CW2*ADELTAKAPPA0' 
            write(*,*)'Thus setting ZDELTAKAPPA0 = ', zDkappa0_0
         else if ((zDg0_0.eq.0D0) .and. 
     -        ((zDkappa0_0.ne.0D0) .or. (aDkappa0_0.ne.0D0))) then
            zDg0_0 = zDkappa0_0 + SW2/CW2*aDkappa0_0
            write(*,*)'Note!' 
            write(*,*)'ZDELTAG0 = ZDELTAKAPPA0+SW2/CW2*ADELTAKAPPA0'
            write(*,*)'Thus setting ZDELTAG0 = ', zDg0_0
         else if ((aDkappa0_0.eq.0D0) .and. 
     -        ((zDkappa0_0.ne.0D0) .or. (zDg0_0.ne.0D0))) then
            aDkappa0_0 = CW2/SW2*(zDg0_0-zDkappa0_0)
            write(*,*)'Note!' 
            write(*,*)'ADELTAKAPPA0 = CW2/SW2*(ZDELTAG0-ZDELTAKAPPA0)'
            write(*,*)'Thus setting ADELTAKAPPA0 = ', aDkappa0_0
         end if

** Consistency check on kappa0, lambda0, g
         consistency = (zDkappa0_0/(SW2*MZ2)) + (aDkappa0_0/MW2)
         consistency = consistency - (zDg0_0/(SW2*MZ2))
         consistency = sqrt(consistency**2)
         if (consistency .gt. 1D-7) then
            aDkappa0_0 = CW2/SW2*(zDg0_0-zDkappa0_0)
            write(*,*)'WARNING! The values for ZDELTAKAPPA0,'
            write(*,*)'ADELTAKAPPA0 AND ZDELTAG0 are not consistent.'
            write(*,*)'We will use ZDELTAKAPPA0 AND ZDELTAG0 as input.'
            write(*,*)'This gives ADELTAKAPPA0 =', aDkappa0_0
            write(*,*)'If you want to use ADELTAKAPPA0 as input'
            write(*,*)'set either ZDELTAKAPPA0 or ZDELTAG0 to 0.'
         end if

** Setting individual formfactors
         if (formfact .and. formfac_ind) then
            call read_real("MASS_SCALE_LAMBDA",massscale2L,2000D0)
            call read_int("FFEXP_LAMBDA",ffexpL,2)
            massscale2L = massscale2L**2
            call read_real("MASS_SCALE_G",massscale2G,2000D0)
            call read_int("FFEXP_G",ffexpG,2)
            massscale2G = massscale2G**2
            call read_real("MASS_SCALE_ZKAPPA",massscale2KZ,2000D0)
            call read_int("FFEXP_ZKAPPA",ffexpKZ,2)
            massscale2KZ = massscale2KZ**2
            call read_real("MASS_SCALE_AKAPPA",massscale2KA,2000D0)
            call read_int("FFEXP_AKAPPA",ffexpKA,2)
            massscale2KA = massscale2KA**2

* Check for formfactor mass scales set to zero
            if (massscale2L .eq. 0d0) then
               write(*,*)"You have input 0 GeV mass scales for lambda "
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", ffmassscale2
               massscale2L = ffmassscale2
            end if
               
            if ((massscale2KA .eq. 0d0) .and. (massscale2KZ .eq. 0d0)
     &           .and. (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for all "
               write(*,*)"individual formfactors."
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_ZKAPPA and MASS_SCALE_G"
               massscale2KZ = ffmassscale2
               massscale2G = ffmassscale2
            else if ((massscale2KA .eq. 0d0) .and. 
     &              (massscale2KZ .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_ZKAPPA and MASS_SCALE_AKAPPA"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_ZKAPPA"
               massscale2KZ = ffmassscale2
            else if ((massscale2KA .eq. 0d0) .and.
     &           (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_AKAPPA and MASS_SCALE_G"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_G"
               massscale2G = ffmassscale2
            else if ((massscale2KZ .eq. 0d0) .and.
     &           (massscale2G .eq. 0d0)) then
               write(*,*)"You have input 0 GeV mass scales for"
               write(*,*)"both MASS_SCALE_ZKAPPA and MASS_SCALE_G"
               write(*,*)"We will use the universal formfactor"
               write(*,*)"mass scale =", sqrt(ffmassscale2)
               write(*,*)"for MASS_SCALE_G"
               massscale2G = ffmassscale2
            end if

** Checking for negative exponents in formfactor
            if (ffexpKA .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKA = 2"
               ffexpKA = 2
            end if
            if (ffexpKZ .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpKZ = 2
            end if
            if (ffexpG .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpG = 2
            end if
            if (ffexpL .lt. 0) then
               write(*,*)"You cannot use negative exponents in the"
               write(*,*)"formfactors.  We will instead use"
               write(*,*)"FFEXPKZ = 2"
               ffexpL = 2
            end if

** Consistency checks and settings:
            if ((massscale2KA .eq. 0d0) .and. 
     &           (aDkappa0_0 .ne. 0d0)) then 
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/(SW2*aDkappa0_0)
               massscale2 = 1000d0/(massscale2**
     &              dble(-1d0/ffexpKA) - 1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KA = massscale2
                  write(*,*)'MASS_SCALE_AKAPPA is set to', 
     &                 sqrt(massscale2KA)
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_G and MASS_SCALE_ZKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            else if ((massscale2KZ .eq. 0d0) .and.
     &              (zDkappa0_0 .ne. 0d0)) then
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - (SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA)))/zDkappa0_0
               massscale2 = 1000d0/(massscale2**dble(-1d0/ffexpKZ) - 
     &              1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KZ = massscale2
                  write(*,*)'MASS_SCALE_ZKAPPA is set to', 
     &                 sqrt(abs(massscale2KZ))
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_G and MASS_SCALE_AKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            else if ((massscale2G .eq. 0d0) .and. 
     &              (zDg0_0 .ne. 0d0)) then
               massscale2 = ((SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA)) + zDkappa0_0/((1d0 + 
     &              1000d0/massscale2KZ)**dble(ffexpKZ)))/zDg0_0
               massscale2 = 1000d0/(massscale2**dble(-1d0/ffexpG) - 
     &              1d0)  
               if (massscale2G .gt. 1d0) then
                  massscale2G = massscale2
                  write(*,*)'MASS_SCALE_G is set to', 
     &                 sqrt(abs(massscale2G))
                  write(*,*)'for consistency'
               else
                  write(*,*)"Sorry! We cannot make the input scales"
                  write(*,*)"MASS_SCALE_ZKAPPA and MASS_SCALE_AKAPPA"
                  write(*,*)"consistent.  Please change a value or"
                  write(*,*)"use a universal formfactor."
                  stop
               end if
            end if

            consistency = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/SW2
            consistency = consistency - aDkappa0_0/((1d0 + 1000d0/
     &              massscale2KA)**dble(ffexpKA))
            if (abs(consistency) .gt. 1d-7) then
               massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &              dble(ffexpG)) - zDkappa0_0/((1d0 + 1000d0/
     &              massscale2KZ)**dble(ffexpKZ)))*CW2/(SW2*aDkappa0_0)
               massscale2 = 1000d0/(massscale2**
     &              dble(-1d0/ffexpKA) - 1d0)
               if (massscale2 .gt. 1d0) then
                  massscale2KA = massscale2
                  write(*,*)"The mass scales MASS_SCALE_G,MASS_SCALE_KZ"
                  write(*,*)"and MASS_SCALE_KA are inconsistent."
                  write(*,*)"We will use MASS_SCALE_G AND MASS_SCALE_KZ"
                  write(*,*)"as inputs, which gives"
                  write(*,*)'MASS_SCALE_AKAPPA =', sqrt(massscale2KA)
                  write(*,*)"If you want to use MASS_SCALE_AKAPPA as an"
                  write(*,*)"input, set either MASS_SCALE_G or "
                  write(*,*)"MASS_SCALE_ZKAPPA to zero."
               else
                  massscale2 = (zDg0_0/((1d0 + 1000d0/massscale2G)**
     &                 dble(ffexpG)) - (SW2*aDkappa0_0/CW2)/((1d0 + 
     &                 1000d0/massscale2KA)**dble(ffexpKA)))/zDkappa0_0
                  massscale2 = 1000d0/(massscale2**
     &                 dble(-1d0/ffexpKZ) - 1d0)
                  if (massscale2 .gt. 1d0) then
                     massscale2KZ = massscale2
                     write(*,*)"The mass scales MASS_SCALE_G,"
                     write(*,*)"MASS_SCALE_KZ and MASS_SCALE_KA are" 
                     write(*,*)" inconsistent.  We will use "
                     write(*,*) "MASS_SCALE_G AND MASS_SCALE_KA as"
                     write(*,*)"inputs, which gives"
                     write(*,*)'MASS_SCALE_ZKAPPA =', sqrt(massscale2KZ)
                     write(*,*)"If you want to use MASS_SCALE_AKAPPA "
                     write(*,*)"as an input, set either MASS_SCALE_G"
                     write(*,*)"or MASS_SCALE_ZKAPPA to zero."
                  else
                     massscale2 = ((SW2*aDkappa0_0/CW2)/((1d0 + 1000d0/
     &                    massscale2KA)**dble(ffexpKA)) + zDkappa0_0/ 
     &                    ((1d0 + 1000d0/massscale2KZ)**dble(ffexpKZ)))/
     &                    zDg0_0
                     massscale2 = 1000d0/(massscale2**
     &                    dble(-1d0/ffexpG) - 1d0)  
                     if (massscale2 .gt. 1d0) then
                        massscale2G = massscale2
                        write(*,*)"The mass scales MASS_SCALE_G,"
                        write(*,*)"MASS_SCALE_KZ and MASS_SCALE_KA are" 
                        write(*,*)"inconsistent.  We will use "
                        write(*,*) "MASS_SCALE_KZ AND MASS_SCALE_KA as"
                        write(*,*)"inputs, which gives"
                        write(*,*)'MASS_SCALE_G =', sqrt(massscale2G)
                        write(*,*)"If you want to use MASS_SCALE_G "
                        write(*,*)"as an input, set either "
                        write(*,*)"MASS_SCALE_AKAPPA or "
                        write(*,*)"MASS_SCALE_ZKAPPA to zero."
                     else
                        write(*,*)"Sorry! We cannot make the input"
                        write(*,*)"formfactor mass scales MASS_SCALE_G,"
                        write(*,*)"MASS_SCALE_AKAPPA and"
                        write(*,*)"MASS_SCALE_ZKAPPA consistent. Please"
                        write(*,*)"change one or more values or use a"
                        write(*,*)"universal formfactor."
                        stop
                     end if
                  end if
               end if
                  
            end if               

         end if


** Setting trianom=1 parameterisation
         fwww_0 = 2D0*SW2*lambda0_0/(3D0*MW2*EL*EL)
         fw_0 = 2D0*zDg0_0/MZ2
         fb_0 = (CW2*zDg0_0 - zDkappa0_0)*2D0/(SW2*MZ2)
* Formfactors: note formfactor for FB is not a simple conversion and is thus
*              set in the formfactor subroutine
         massscale2FWWW = massscale2L
         ffexpFWWW = ffexpL
         massscale2FW = massscale2G
         ffexpFW = ffexpG

      else 
         write(*,*)'Invalid entry for TRIANOM:'
         write(*,*)'gauge boson anomalous coupling parametrisation'
         write(*,*)'Please set TRIANOM = 1 or 2'
         stop
      end if


** Some inputs are only used for specific processes
      SELECT CASE(procID)
      CASE(WPWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,WPWPjj,WMWMjj,WPAjj,WMAjj,ZAjj,ZAjj_n,
     &     WWZ,ZZZ,ZZWP,ZZWM,WWWP,WWWM,WWA,ZZA,ZZnA,WPZA,WMZA,WPAA,WMAA,ZAA,ZnAA,AAA,
     &     WPAAj,WMAAj,WW,GFWW,WPZ,WMZ,WPA,WMA,
     &     WPhadWMjj, WPWMhadjj, ZZhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, WPhadWPjj, WMhadWMjj,
     &     WPhadWMZ, WPWMhadZ, WWZhad, ZZZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, WPhadWMWP, WPWMhadWP, 
     &     WMhadWPWM, WMWPhadWM, WPhadWMlep, WPlepWMhad, GFWPhadWMlep, GFWPlepWMhad, WPhadZ, WPZhad, WMhadZ, WMZhad,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     7     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu,
     7     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &     )

** dim-8 operators for triboson and VVjj processes: 
         call read_real("FS0",fs0_0,0.0d0)
         call read_real("FS1",fs1_0,0.0d0)
         call read_real("FS2",fs2_0,0.0d0)
         call read_real("FT0",ft0_0,0.0d0)
         call read_real("FT1",ft1_0,0.0d0)
         call read_real("FT2",ft2_0,0.0d0)
         call read_real("FT5",ft5_0,0.0d0)
         call read_real("FT6",ft6_0,0.0d0)
         call read_real("FT7",ft7_0,0.0d0)
         call read_real("FT8",ft8_0,0.0d0)
         call read_real("FT9",ft9_0,0.0d0)
         call read_real("FM0",fm0_0,0.0d0)
         call read_real("FM1",fm1_0,0.0d0)
         call read_real("FM2",fm2_0,0.0d0)
         call read_real("FM3",fm3_0,0.0d0)
         call read_real("FM4",fm4_0,0.0d0)
         call read_real("FM5",fm5_0,0.0d0)
         call read_real("FM6",fm6_0,0.0d0)
         call read_real("FM7",fm7_0,0.0d0)
** dim-6 operators which are only used for diboson, triboson and VVjj production
c         call read_real("FBW",fbw_0,0.0d0) ! fbw introduces B-W3 mixing
c         call read_real("FDW",fdw_0,0.0d0) ! if the contributions to the gauge-boson propagators are neglected, fdw violates gauge invariance!
         fbw_0 = 0.0d0
         fdw_0 = 0.0d0
         call read_real("FWW",fww_0,0.0d0)
         call read_real("FBB",fbb_0,0.0d0)
         call read_real("FBWT",fbwtilde_0,0.0d0)
         call read_real("FDWT", fdwtilde_0,0.0d0)   
         call read_real("FWWT",fwwtilde_0,0.0d0)
         call read_real("FBBT",fbbtilde_0,0.0d0)
         call read_real("FWT",fwtilde_0,0.0d0)
         call read_real("FBT",fbtilde_0,0.0d0)
         call read_real("FWWWT", fwwwtilde_0,0.0d0)
c         call read_real("FPHI1", fphi1_0,0.0d0) ! fphi1 contributes to the Z boson mass
         fphi1_0 = 0.0d0
      END SELECT

      
      call closefile
*      print *," "


** converting to HVV coupling notation, for Higgs width calculation
      treefacZ = 1d0
      treefacW = 1d0
      loopfac = 1d0
      dum = 0d0

      Select Case(ProcID)
      CASE(WPWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,WPWPjj,WMWMjj,WPAjj,WMAjj,ZAjj,ZAjj_n,
     &     WWZ,ZZZ,ZZWP,ZZWM,WWWP,WWWM,WWA,ZZA,ZZnA,WPZA,WMZA,WPAA,WMAA,ZAA,ZnAA,AAA,
     &     WPAAj,WMAAj,GFWW,WW,WPZ,WMZ,WPA,WMA,
     &     WPhadWMjj, WPWMhadjj, ZZhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, WPhadWPjj, WMhadWMjj,
     &     WPhadWMZ, WPWMhadZ, WWZhad, ZZZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, WPhadWMWP, WPWMhadWP, 
     &     WMhadWPWM, WMWPhadWM, GFWPhadWMlep, GFWPlepWMhad, WPhadWMlep, WPlepWMhad, WPhadZ, WPZhad, WMhadZ, WMZhad,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     7     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu,
     7     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &     )

         call anomH_convert(4,4,(trianom .eq. 2),(trianom .eq. 1),
     &        dum, dum, dum, dum, zDg0_0, aDkappa0_0, dum,
     &        fww_0, fwwtilde_0, fbb_0, fbbtilde_0, fw_0, fb_0, fbtilde_0, 
     &        fbw_0, fbwtilde_0, fwtilde_0, fphi1_0)
      
      CASE DEFAULT

         call anomH_convert(4,4,(trianom .eq. 2),(trianom .eq. 1),
     &        dum, dum, dum, dum, zDg0_0, aDkappa0_0, dum,
     &        dum, dum, dum, dum, fw_0, fb_0, dum,
     &        dum, dum, dum, dum)

      END SELECT


** Setting formfactors
      fwww=fwww_0
      fw=fw_0
      fb=fb_0

      lambda0 = lambda0_0
      zDg0 = zDg0_0
      zDkappa0 = zDkappa0_0
      aDkappa0 = aDkappa0_0

      Select Case(ProcID)
      CASE(WPWMjj,ZZjj_ll,ZZjj_lnu,WPZjj,WMZjj,WPWPjj,WMWMjj,WPAjj,WMAjj,ZAjj,ZAjj_n,
     &     WWZ,ZZZ,ZZWP,ZZWM,WWWP,WWWM,WWA,ZZA,ZZnA,WPZA,WMZA,WPAA,WMAA,ZAA,ZnAA,AAA,
     &     WPAAj,WMAAj,WW,WPZ,WMZ,WPA,WMA,
     &     WPhadWMjj, WPWMhadjj, ZZhadjj, WPhadZjj, WPZhadjj, WMhadZjj, WMZhadjj, WPhadWPjj, WMhadWMjj,
     &     WPhadWMZ, WPWMhadZ, WWZhad, ZZZhad, ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, WPhadWMWP, WPWMhadWP, 
     &     WMhadWPWM, WMWPhadWM, WPhadWMlep, WPlepWMhad, WPhadZ, WPZhad, WMhadZ, WMZhad,
     &     WPhadWMA, WPWMhadA, ZZhadA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     7     WPH, WPH_AA, WPH_mu, WPH_tau, WPH_bbar, WPH_WW, WPH_ZZ_ll,  WPH_ZZ_lnu,
     7     WMH, WMH_AA, WMH_mu, WMH_tau, WMH_bbar, WMH_WW, WMH_ZZ_ll,  WMH_ZZ_lnu,
     7     WPHJ, WPHJ_AA, WPHJ_mu, WPHJ_tau, WPHJ_bbar, WPHJ_WW, WPHJ_ZZ_ll,  WPHJ_ZZ_lnu,
     7     WMHJ, WMHJ_AA, WMHJ_mu, WMHJ_tau, WMHJ_bbar, WMHJ_WW, WMHJ_ZZ_ll,  WMHJ_ZZ_lnu
     &     )
! don't need to list GFWW here as this process only uses the anom HVV couplings which were filled with anomH_convert

         fs0=fs0_0
         fs1=fs1_0
         fs2=fs2_0
         ft0=ft0_0
         ft1=ft1_0
         ft2=ft2_0
         ft5=ft5_0
         ft6=ft6_0
         ft7=ft7_0
         ft8=ft8_0
         ft9=ft9_0
         fm0=fm0_0
         fm1=fm1_0
         fm2=fm2_0
         fm3=fm3_0
         fm4=fm4_0
         fm5=fm5_0
         fm6=fm6_0
         fm7=fm7_0

         fbw = fbw_0
         fdw = fdw_0
         fww = fww_0
         fbb = fbb_0
         fbwtilde = fbwtilde_0
         fdwtilde = fdwtilde_0
         fwwtilde = fwwtilde_0
         fbbtilde = fbbtilde_0
         fwtilde = fwtilde_0
         fbtilde = fbtilde_0
         fwwwtilde = fwwwtilde_0
         fphi1=fphi1_0
      END SELECT


      end


      subroutine set_anomVcouplings_zero
c     this subroutine sets all relevant couplings to SM values,
c     so that the anomalous routines can also be used for SM calculations
c     when anom=false

#include "an_couplings.inc"   

      ! form factor for diboson, triboson, VBF, diboson+jet
      formfact=.false.

      ! couplings for triboson
      fs0=0
      fs1=0
      fs2=0
      ft0=0
      ft1=0
      ft2=0
      ft5=0
      ft6=0
      ft7=0
      ft8=0
      ft9=0
      fm0=0
      fm1=0
      fm2=0
      fm3=0
      fm4=0
      fm5=0
      fm6=0
      fm7=0

      ! couplings for diboson, triboson, VBF
      fwww=0
      fw=0
      fb=0

      ! couplings for VBF
      fbw = 0
      fdw = 0
      fww = 0
      fbb = 0
      fbwtilde = 0
      fdwtilde = 0
      fwwtilde = 0
      fbbtilde = 0
      fwtilde = 0
      fbtilde = 0
      fwwwtilde = 0

      ! couplings for diboson+jet
      lambda0 = 0
      zDg0 = 0
      zDkappa0 = 0
      aDkappa0 = 0


      return
      end

********************************************************************************
********************************************************************************

