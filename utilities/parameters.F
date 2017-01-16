********************************************************************************
********************************************************************************
**                                                                           ***
******************************************************************           ***
*** sophy@particle.uni-karlsruhe.de                            ***           ***
***                                                            ***           ***
*** This set of subroutines sets all of the parameters used    ***           ***
*** in the calculation of electroweak corrections in the code: ***           ***
*** i.e. it calls FeynHiggs, and calculates Higgs masses and   ***           ***
*** renormalisation constants.                                 ***           ***
***                                                            ***           ***
*** VERY IMPORTANT NOTE!!!                                     ***           ***
*** The implementation of the user input EWSCHEME was changed  ***           ***
*** for version 2.5.0 of VBFNLO.  In order to switch back to   ***           ***
*** the old implementation, comment out the NEW subroutine     ***           ***
*** 'setEWpara' and uncomment the OLD subroutine 'setEWpara',  ***           ***
*** which are the first and second subroutines in this file    ***           ***
*** respectively.                                              ***           ***
*** For more details, see the wiki.                            ***           ***
******************************************************************           ***
**                                                                           ***
********************************************************************************
********************************************************************************

*** NEW v 2.5.0 subroutine that sets electroweak parameters according to 
*** EWSCHEME, and e and g2, which control the couplings.

      subroutine setEWpara(e,g2,s,c,z,w,q,g)

      implicit none

      double precision e, g2, s, c, z, w, q, g, xmwsq


#include "koppln.inc"
#include "global.inc"
#include "mssm.inc"


** ewscheme = 1: MW2 and SW2 are set according to the inputs alfa, gfermi and
** MZ.  NOTE ALSO:  using this scheme can lead to odd values of the W mass.  
** A warning is output if this is the case.
      if (ewscheme .eq. 1) then
         XMWSQ = XMZ**2/2.d0 + SQRT(XMZ**4/4.d0 - PI*ALFA*XMZ**2/
     -        SQRT(2.d0)/GF) 
         SIN2W = 1.d0 - XMWSQ/XMZ**2
         XMW = XMZ*SQRT(1.d0 - SIN2W)

** ewscheme = 2: alfa_qed and W mass are set according to the inputs gfermi,
** MZ and SW2
      else if (ewscheme .eq. 2) then
         XMW = XMZ*SQRT(1.d0 - SIN2W)

** ewschemes 3, 5 and 6: SW2 is set according to the input MW and MZ
      else if ((ewscheme .eq. 3) .or. (ewscheme .ge. 5)) then
        SIN2W = 1.d0 - (XMW/XMZ)**2
      end if


      G = SQRT(4*PI*ALFAS)
      S = SQRT(SIN2W)
      C = SQRT(1.d0 - SIN2W)
   

** G2 = E/S IS CALCULATED FROM GFERMI AND THE Z MASS
      G2 = SQRT(8.d0*GF/SQRT(2.d0))*XMZ*C
 
** For ewscheme = 1,5,6, all couplings are set according to the input alfa.
      if ((ewscheme .eq. 1) .or. (ewscheme .ge. 5)) then
         E = SQRT(4.d0*PI*ALFA)
         G2 = E/S

** For ewscheme 4, photon couplings are set according to the input alfa, and 
** all other couplings are set according to the input gfermi.  NOTE: This can
** lead to problems with gauge invariance if alfa and gfermi are not consistent!
      else if (ewscheme .eq. 4) then
         E = SQRT(4.d0*PI*ALFA)

** For ewscheme 2 and 3, all couplings are set according to the input gfermi
      else if ((ewscheme .eq. 2) .or. (ewscheme .eq. 3)) then
         E = G2*S
         ALFA = E**2/(4.d0*PI)
      end if


      Z = G2/4.d0/C
      W = G2/SQRT(8.d0)
      Q = G2*C


* call to subroutine setting parameters 
*      if ((HiggsType .ne. 0) .or. ewcor_switch) then
      call setparams(G2)
*      end if


** WARNING MESSAGES:
      if (ewscheme .eq. 4) then
         write(*,*)'  '
         write(*,*)'WARNING! Note that the inputs:'
         write(*,*)'  ALFA, FERMI_CONST, WMASS, ZMASS and SIN2W '
         write(*,*)'are not independent quantities.  If they are not'
         write(*,*)'consistent, problems with gauge invariance may'
         write(*,*)'arise.'
         write(*,*)'In this scheme, the photon couplings are set'
         write(*,*)'according to the input ALFA, and all other'
         write(*,*)'couplings are set according to FERMI_CONST.'
         write(*,*)'  '
      end if

      if (XMW .lt. 80D0) then
         write(*,*)'  '
         write(*,*)'WARNING!! MW =', XMW
         write(*,*)'This is very low!'
         write(*,*)'  '
      end if
    

      end


********************************************************************************
********************************************************************************

c$$$*** OLD (pre v 2.5.0) subroutine that sets electroweak parameters according 
c$$$*** to EWSCHEME, and e and g2, which control the couplings.
c$$$
c$$$      subroutine setEWpara(e,g2,s,c,z,w,q,g)
c$$$
c$$$      implicit none
c$$$
c$$$      double precision e, g2, s, c, z, w, q, g, xmwsq
c$$$
c$$$** Flag that will warn if there are problems with consistency - i.e. fermion
c$$$** and higgs couplings set according to different inputs
c$$$      integer iflag
c$$$
c$$$#include "koppln.inc"
c$$$#include "global.inc"
c$$$#include "mssm.inc"
c$$$
c$$$
c$$$      iflag = 0
c$$$
c$$$
c$$$* for 1/alpha_QED < 136 the Weak angle is calculated from the input
c$$$* values of alpha(assumed at scale m_W**2), m_Z and G_Fermi
c$$$
c$$$** setting MW and sin_theta_w according to the input variable ewscheme
c$$$      if (ewscheme .eq. 1) then
c$$$         if (alfa .gt. 1d0/136d0) then
c$$$            XMWSQ = XMZ**2/2.d0 + SQRT(XMZ**4/4.d0 - 
c$$$     -           PI*ALFA*XMZ**2/SQRT(2.d0)/GF) 
c$$$            SIN2W = 1.d0 - XMWSQ/XMZ**2
c$$$            XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$         else
c$$$** NOTE: if alfa_qed is input as less than or equal to alfa(0),sin^2 theta_w
c$$$**       is fixed at 0.2312.
c$$$            SIN2W = 0.2312d0
c$$$            XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$         end if
c$$$      else if (ewscheme .eq. 2) then
c$$$         XMW = XMZ*SQRT(1.d0 - SIN2W)
c$$$      else if ((ewscheme .eq. 3) .or. (ewscheme .ge. 5)) then
c$$$        SIN2W = 1.d0 - (XMW/XMZ)**2
c$$$      end if
c$$$
c$$$
c$$$      G = SQRT(4*PI*ALFAS)
c$$$      S = SQRT(SIN2W)
c$$$      C = SQRT(1.d0 - SIN2W)
c$$$   
c$$$c
c$$$c G2 = E/S IS CALCULATED FROM GFERMI AND THE Z MASS
c$$$c
c$$$      G2 = SQRT(8.d0*GF/SQRT(2.d0))*XMZ*C
c$$$      IF ( ALFA.GT.1.d0/136.d0 .or. ewscheme .ge. 5) THEN
c$$$** ewscheme 1 or 4 if alfa > alfa(0), and ewscheme 5 and 6
c$$$        E = SQRT(4.d0*PI*ALFA)
c$$$        G2 = E/S
c$$$      ELSEIF ( ALFA.GT.1.d0/138.d0 ) THEN
c$$$** ewscheme 1 or 4 if alfa = alfa(0): e calculated from (input) alfa, G2 from 
c$$$** g_fermi. WARNING! gauge cancellations may no longer work!
c$$$        E = SQRT(4.d0*PI*ALFA)
c$$$        iflag = 1
c$$$      ELSE
c$$$** ewscheme 2 or 3 or if alfa < alfa(0): e, G2 and alfa all calculated from 
c$$$** g_fermi. 
c$$$        E = G2*S
c$$$        ALFA = E**2/(4.d0*PI)
c$$$        if (ewscheme .eq. 1 .or. ewscheme .eq. 4) then
c$$$           iflag = 2
c$$$        end if
c$$$      ENDIF
c$$$
c$$$
c$$$      Z = G2/4.d0/C
c$$$      W = G2/SQRT(8.d0)
c$$$      Q = G2*C
c$$$
c$$$
c$$$* call to subroutine setting parameters 
c$$$*      if ((HiggsType .ne. 0) .or. ewcor_switch) then
c$$$      call setparams(G2)
c$$$*      end if
c$$$
c$$$
c$$$      if (iflag .eq. 1) then
c$$$         write(*,*)"WARNING! The Higgs coupling is set from G_FERMI,"
c$$$         write(*,*)"and the electron charge is set from the input ALFA."
c$$$         write(*,*)"Gauge cancellations may no longer work"
c$$$      else if (iflag .eq. 2) then
c$$$         write(*,*)"NOTE that all couplings are set from G_FERMI, not"
c$$$         write(*,*)"from the input ALFA"
c$$$      end if
c$$$
c$$$
c$$$      end


********************************************************************************
********************************************************************************

***  Subroutine that sets all (non-ewscheme) parameters 

      subroutine setparams(G2)
          use globalvars, only: lglobalprint, ldoblha
          use readinput

      implicit none


*** The following three parameters are used by FeynHiggs.  They are set in the
*** subroutine 'bench_set':

** This is a FeynHiggs flag: t1_CplxApprox determines how the two-loop 
** corrections are treated in the presence of complex parameters
      integer t1_CplxApprox

** The soft-SUSY breaking parameters in the sfermion sector
      double precision MSusy

** Scale information for FeynHiggs
      double precision scalefactor


** These do-loop counters are used to set squark and other SUSY parameters 
      integer sfe1, sfe2, gen, type
      integer neu1, neu2, cha1, cha2

** Flag identifying which of the neutral/charged currents is being considered
      integer cur

** Variables used for choosing the charge renormalisation scheme
      integer charge
      double precision Alfa0, AlfaMZ, AlfaGF
      double precision delA0, delAmz, delAgf
      double precision G2

** Flags for the PDF set which is used
      integer pdf

** Flag for whether to include bottom Yukawa coupling, and value of Delta MB
** from FeynHiggs
      logical delMB_switch
      double complex DeltaMB


** By (annoying) convention, the sign of the su(2) covariant derivative is 
** different in the SM and MSSM. This factor takes care of this fact, and 
** ensures everything's consistent.
      double precision su2sign


** Dummy parameters used in calculation of the tree level higgs sector
      double precision avgMH2, deltaMH2
      double precision Mh02tree, MHH2tree, MA02tree, MHp2tree


** Dummy parameters: these are used to calculate the SUSY parameters (sfermion,
** chargino, neutralino sectors) if they haven't been set by either a SLHA file
** or FeynHiggs
      logical first
      double precision MSfcalc(2,4,3), MChaCalc(2), MNeuCalc(4)
      double complex USfcalc(2,2,4,3), UChaCalc(2,2), VChaCalc(2,2)
      double complex ZNeuCalc(4,4)


** Function which determines whether a parameter has already been set with a 
** reasonable value (for real and complex parameters)
      logical replace
      external replace

** Debugging flag
      integer debug
      
      
** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "koppln.inc"
#include "process.inc"

      ! set debug
      debug = 0


      if (model .eq. 2) then  ! We're working in the MSSM
** Reading input from file susy.dat unless we've already read from a SLHA
** file.
         call loadfile("susy.dat",.true.)

         write(*,*)" "
         write(*,*)"              Global SUSY parameters"
         write(*,*)"-----------------------------------------------"
     
         call read_Int("PROPLEVEL", propLevel, 1)
         call read_logical("DELMB_SWITCH",delMB_switch,.false.) 
         call read_Int("MH_LOOPS", MHloops, 0)
         call read_logical("LOOPSQR_SWITCH",loopsqr_switchIN,.true.) 
         call closefile
      end if

** Read the PDF set specified (default is hard-wired cteq) and delta alfa 
** (unless it's been set by slha)
      if ((sector .ge. 3) .or. (ewscheme .eq. 6)) then
         call loadfile("vbfnlo.dat",.true.)
         call read_Int("PDF_SWITCH",pdf,0)
         if (replace(delAlfa,0)) then
            call read_Real("DEL_ALFA",delAlfa,0.059047686d0)
         end if
         call closefile
      end if
         
** Messages that record the input details
      if (model .eq. 2) then   ! We're working in the MSSM

** Which Higgs is being produced
         write(*,*)'  '
         if (HiggsType .eq. 1) then
            write(*,*)'Producing the light CP-even higgs'
         else if (HiggsType .eq. 2) then
            write(*,*)'Producing the heavy CP-even higgs'
         else if (HiggsType .eq. 3) then
            write(*,*)'Producing the CP-odd higgs'
         end if

         write(*,*)'  '
         if (propLevel .ne. 0) then
#ifndef WITH_LT
            if (.not. FeynH_switch) then
               write(*,*)'Sorry!  We cannot use the Higgs propagator'
               write(*,*)'Z matrix if neither FeynHiggs nor LoopTools'
               write(*,*)'is being used'
               write(*,*)'Switching to propLevel = 0'
               propLevel = 0
            end if
#endif
         end if
            
         write(*,*)'  '
         if (propLevel .eq. 0) then
            write(*,*)'No Higgs propagator corrections are included'
            write(*,*)'Instead, we use an effective Higgs mixing angle'
            write(*,*)'NOTE! This is not necessarily a good'
            write(*,*)'approximation in many areas of the MSSM'
         else if (propLevel .eq. 1) then
            write(*,*)'Higgs propagator corrections included' 
            write(*,*)'at leading order'
         end if
         write(*,*)'  '
         if (ewcor_switch) then
            if (propLevel .eq. 2) then
               write(*,*)'Higgs propagator corrections included at' 
               write(*,*)'leading order and loop level'
            else if (propLevel .eq. 3) then
               write(*,*)'Higgs propagator corrections included as a' 
               write(*,*)'loop correction'
            end if

** Which electroweak corrections are included
            if (ewapprox .eq. 0) then
               write(*,*)'SORRY, the complete electroweak corrections'
               write(*,*)'in the MSSM are not yet included.  Set'
               write(*,*)'EWAPPROX to 5 for the most-complete MSSM'
               write(*,*)'corrections presently available.'
               stop
            else if (ewapprox .eq. 1) then
               write(*,*)'Electroweak corrections from the third ' 
               write(*,*)'generation (s)quarks are included'
            else if (ewapprox .eq. 2) then
               write(*,*)'Electroweak corrections from (s)fermions' 
               write(*,*)'are included'
            else if (ewapprox .eq. 3) then
               write(*,*)'All SM-type corrections are included,'
               write(*,*)'together with sfermion corrections'
            else if (ewapprox .eq. 4) then
               write(*,*)'Complete MSSM corrections to the hVV'
               write(*,*)'vertex are included, together with'
               write(*,*)'SM-type and sfermion corrections'
               write(*,*)'elsewhere'
            else if (ewapprox .eq. 5) then
               write(*,*)'Complete MSSM corrections to the hVV'
               write(*,*)'vertex, VV self energy and qqV vertex'
               write(*,*)'are included, together with SM-type'
               write(*,*)'boxes and pentagons'
            end if

** Checking mrst2004qed is chosen iff full electroweak corrections are
** wanted
            if (sector .ge. 3) then
               if (pdf .ne. 2) then
                  write(*,*)' '
                  write(*,*)'Note that photon-induced processes can' 
                  write(*,*)'only be included if the pdfset MRST2004QED' 
                  write(*,*)'is chosen (PDF_SWITCH=2 in vbfnlo.dat).'
               end if

** Output for internal higgs mass tree/corrected
               write(*,*)' '
               if (MHloops .eq. 0) then
                  write(*,*)'Higgs masses within electroweak loops are'
                  write(*,*)'set to their tree-level values'
               else if (MHloops .eq. 1) then
                  write(*,*)'Higgs masses within electroweak loops are'
                  write(*,*)'set to their corrected values'
               else
                  write(*,*)'Please input MH_LOOPS = 0 (tree level '
                  write(*,*)'Higgs masses in loops) or 1 (corrected '
                  write(*,*)'Higgs masses in loops) in susy.dat'
                  stop
            end if
            end if

** Output for loop^2 contributions
            write(*,*)' '
            if (loopsqr_switchIN) then
               write(*,*)'Squared electroweak corrections from fermions'
               write(*,*)'will be included if'
               write(*,*)'     |M_{(s)fermion loop}| > 0.1 |M_{born}|'
            else
               write(*,*)'Squared electroweak corrections from fermions'
               write(*,*)'are neglected'
            end if

         end if

      else   !  working in the SM

** Standard Model info messages
** Higgs mass message
         if (lglobalprint) write(*,*)'  '
         if (FeynH_switch .or. slha_switch) then
            if (HiggsType .eq. 1) then
               write(*,*)'MH = mass of the light cp-even higgs'
            else if (HiggsType .eq. 2) then
               write(*,*)'MH = mass of the heavy cp-even higgs'
            else if (HiggsType .eq. 3) then
               write(*,*)'MH = mass of the cp-odd higgs'
            end if
         else
            if (lglobalprint) write(*,*)'Higgs mass =', XMH
         end if

** Which electroweak corrections are included
         if (lglobalprint) write(*,*)'  '
         if (sector .eq. 1) then
            write(*,*)'Electroweak corrections from the third' 
            write(*,*)'generation quarks are included'
         else if (sector .eq. 2) then
            write(*,*)'Electroweak corrections from fermions' 
            write(*,*)'are included'
         else if (sector .ge. 3) then
            write(*,*)'All electroweak corrections are included'
** Checking mrst2004qed is chosen iff full electroweak corrections are
** wanted
            if (pdf .ne. 2) then
               write(*,*)'Note that photon-induced processes can only' 
               write(*,*)'be included if the pdfset MRST2004QED is' 
               write(*,*)'chosen (PDF_SWITCH=2 in vbfnlo.dat).'
            end if
         end if
      end if


** Electroweak sector - PDG (ntbk 7, p87)
*      MZ = 91.1876D0
*      MW = 80.398D0
** Matching 0710.4749 
*      MZ = 91.1535D0
*      MW = 80.397D0

** Setting gauge boson parameters to match input (read either from SLHA file - 
** in slha_read.F, or vbfnlo.dat - in koppln.F)
      MZ = DBLE(xmz)
      MW = DBLE(xmw)

      MZ2 = MZ**2D0
      MW2 = MW**2D0

      SW2 = Sin2W
      SW = sqrt(SW2)
      CW2 = 1 - SW2
      CW = sqrt(CW2)

*      CW2 = MW2/MZ2
*      SW2 = 1 - CW2
*      CW = sqrt(CW2)
*      SW = sqrt(SW2)


      if (ewcor_switch) then
         if (((MW2/MZ2)-CW2) .gt. 0.001D0) then
            write(*,*)'  '
            write(*,*)'STOPPING CODE!!!'
            write(*,*)'The choice of ewscheme means that'
            write(*,*)'MW2/MZ2 is not equal to CW2.'
            write(*,*)'The electroweak corrections will'
            write(*,*)'not be finite.'
            stop
         end if
      end if

   
** Light fermion masses - set according to the PDG if they haven't been set in 
** the slha file.  c,, t and b masses will have been set by this point, either 
** by a slha file or via the .dat input files
      if (replace(Mf(2,1),0)) then
         ME = 0.00051099891D0
      else
         ME = Mf(2,1)
      end if
      if (replace(Mf(3,1),0)) then
         MsU = 0.0024D0
      else
         MsU = Mf(3,1)
      end if
      if (replace(Mf(4,1),0)) then
         MD = 0.00475D0
      else
         MD = Mf(4,1)
      end if
      if (replace(Mf(2,2),0)) then
         MM = 0.105658367D0
      else
         MM = Mf(2,2)
      end if
      if (replace(Mf(4,2),0)) then
         MS = 0.104D0
      else
         MS = Mf(4,2)
      end if

** b,t,c and tau masses - input (either slha or vbfnlo.dat)
      if (ldoblha) then
      MC = Mf(3,2)
      MB = Mf(4,2)
      MT = Mf(3,3)
      ML = Mf(2,3)
      else
      MC = DBLE(xmc)
      MB = DBLE(xmb)
      MT = DBLE(xmt)
      ML = DBLE(xmtau)
      endif

** Quark masses - PDG
*      MC = 1.27D0 
*      MT = 172.6D0   
*      MB = 4.2D0  
** Matching 0804.2676
*      MT = 170.9D0
*      MB = 4.7D0
** Matching 0710.4749
*      MsU = 0.0066D0
*      MC = 1.2D0
*      MT = 174.3D0      
*      MD = 0.0066D0
*      MS = 0.150D0
*      MB = 4.3D0
*      ME = 0.000510998927D0
*      MM = 0.105658369D0
*      ML = 1.77699D0

      MU2 = MsU**2D0
      MC2 = MC**2D0
      MT2 = MT**2D0
      MD2 = MD**2D0
      MS2 = MS**2D0
      MB2 = MB**2D0


** Choosing charge renormalisation scheme based on value of AlfaQED and ewscheme
*      AlfaQED = (G2**2D0)*SW2/(4D0*Pi)
      AlfaQED = Alfa
      AlfaQED2 = AlfaQED**2D0
      gfermi = GF

      if (ewcor_switch) then

         if (ewscheme .eq. 3) then
** define LO alfa via gfermi: 
            charge = 2
         else if (ewscheme .eq. 5) then
** define LO alfa as alfa(MZ)
            charge = 1
         else if (ewscheme .eq. 6) then
** define LO alfa as alfa(0)
            charge = 0
         else 
** in other schemes, the the code chooses whether alfa(0), alfa(MZ) or alfa(GF)
** is more appropriate
            AlfaGF = 1.41421356237309504880168872421D0
            AlfaGF = (AlfaGF/Pi*gfermi*MW2*SW2)
         
            Alfa0 = 1/137.03599981D0
            
            AlfaMZ = Alfa0/(1D0-delAlfa)

            delA0 = (AlfaQED -Alfa0)**2D0
            delAmz = (AlfaQED - AlfaMZ)**2D0
            delAgf = (AlfaQED - AlfaGF)**2D0

            if ((delA0 .lt. delAmz) .and. (delA0 .lt. delAgf)) then
               charge = 0
            else if ((delAmz .lt. delA0) .and. (delAmz .lt. delAgf)) then
               charge = 1
            else if ((delAgf .lt. delA0) .and. (delAgf .lt. delAmz)) then
               charge = 2
            end if

            write(*,*)'  '
            write(*,*)'Note that the code has chosen a charge'
            write(*,*)'renormalisation scheme that best fits the input'
            write(*,*)'electroweak parameters.  For better results,'
            write(*,*)'choose ewscheme = 3, 5 or 6.'
            
         end if

         write(*,*)'  '
         if (charge .eq. 0) then
            write(*,*)'Charge renormalisation: alfa(0)'
         else if (charge .eq. 1) then
            write(*,*)'Charge renormalisation: alfa(MZ)'
         else if (charge .eq. 2) then
            write(*,*)'Charge renormalisation: alfa(GF)'
            if ((model .eq. 2)) then
               if ((ewapprox .eq. 0) .or. (ewapprox .ge. 4)) then
                  write(*,*)'NOTE THAT, when working in the full MSSM,'
                  write(*,*)'and parametrising LO by gfermi, an'
                  write(*,*)'approximation is used such that the finite'
                  write(*,*)'Delta R does not include contributions'
                  write(*,*)'charginos and neutralinos'
               end if
            end if
         end if
      end if

     
      EL = sqrt(4D0*Pi*AlfaQED)
*      write(*,*)'Alfa =', AlfaQED


** These are the values that the charges and L/R axial vector couplings can
** take for up-type and down-type fermions.  The variable 'su2sign' takes care
** of the always-fun sin-theta-w conventions
      if (model .eq. 2) then
         su2sign = -1D0
      else
         su2sign = 1D0
      end if

      Qu = 2/3D0
      Qd = -1/3D0
      
      gL3 = su2sign*((1/2D0)-(2*SW2/3D0))/(SW*CW)
      gR3 = -su2sign*(2*SW/(3D0*CW))
      gL4 = su2sign*(-(1/2D0)+(SW2/3D0))/(SW*CW)
      gR4 = su2sign*SW/(3D0*CW)


** Setting renormalisation scales
      scalefactor = 1D0
      massQ = 1D-5
#ifdef WITH_LT
      if (ewcor_switch) call renSCALES(scalefactor,1)
#else
      if (ewcor_switch) then
         write(*,*)'Sorry, we cannot calculate the electroweak'
         write(*,*)'corrections if LoopTools has not been linked.'
         stop
      end if
#endif


** Now setting the SUSY parameters according to the chosen benchmark scenario 
** parameter: 'bench' (set in susy.dat) or the slha file
      if (FeynH_switch .or. slha_switch) then
         call bench_set(higgsmix, TB, MH, MA0, MHp,
     &        M3SL, M3SE, M3SQ, M3SU, M3SD,
     &        M2SL, M2SE, M2SQ, M2SU, M2SD,
     &        M1SL, M1SE, M1SQ, M1SU, M1SD,
     &        MUE,
     &        Af(2,3,3), Af(3,3,3), Af(4,3,3),
     &        Af(2,2,2), Af(3,2,2), Af(4,2,2),
     &        Af(2,1,1), Af(3,1,1), Af(4,1,1),
     &        M_1, M_2, M_3, MSusy, slha_switch, bench, debug)
         if (higgsmix .eq. 2) t1_CplxApprox = 0
         if (higgsmix .eq. 3) t1_CplxApprox = 1
      end if

      if (model .eq. 2) then  ! working in the MSSM
         if ((procID .gt. 4120) .or. (procID .lt. 3999)) then 
* i.e. not gluon fusion
            if (higgsmix .ne. 3 .and. HiggsType .eq. 3) then
               write(*,*)'The CP-odd Higgs boson is not produced at'
               write(*,*)'tree level unless we are working in the'
               write(*,*)'complex MSSM.'
*               if (.not. ewcor_switch) then
                  stop
*               else
*                  doBorn = .false.
*               end if
            else if (higgsmix .eq. 3 .and. HiggsType .eq. 3) then
               if (propLevel .eq. 0 .or. propLevel .eq. 3) then
                  write(*,*)'The CP-odd Higgs boson is not produced at'
                  write(*,*)'tree level unless we include propagator'
                  write(*,*)'corrections at LO and are working in the'
                  write(*,*)'complex MSSM.'
*                  if (.not. ewcor_switch) then
                    stop
*                  else 
*                     doBorn = .false.
*                  end if
               end if
            end if
         end if  ! end of gluon fusion
      end if  ! end of mssm


** FeynHiggs parameter setting:
      if (FeynH_switch .and. (HiggsType .ne. 0)) then

** FeynHiggs parameter setting:
#ifdef WITH_FH
         call FHpara(t1_CplxApprox, MSusy, scalefactor, DeltaMB, debug)
#else
         write(*,*)'Sorry, you cannot set FeynH_switch to TRUE'
         write(*,*)'if you have not linked to FeynHiggs.'
         write(*,*)'Use --with-FEYNHIGGS=/path/to/feynhiggs'
         stop
#endif

         if (delMB_switch) then
            write(*,*)'  '
            write(*,*)'Corrections to the bottom Yukawa coupling'
            write(*,*)'are included, with Delta MB =', DeltaMB
         end if

         MH = DBLE(MHiggs(HiggsType))
         MH2 = MH**2D0

** Setting Higgs mass as used elsewhere in the code
         XMH = DBLE(MHiggs(HiggsType))


** Setting 'safe' corrected Higgs masses
         do type = 1, 4
            MHiggs_loop(type) = MHiggs(type)
         end do


** Correcting the bottom Yukawa coupling, according to the flag delMB_switch 
** (set in susy.dat)
         if (delMB_switch) then
            MB = MB/abs(1D0 + DeltaMB)
            MB2 = MB**2
         end if

      else  !  if we're not using FeynHiggs
         MH = XMH
         MH2 = MH**2D0
      end if

      if (lglobalprint) then
         write(*,*)'  '
         if (higgsscheme .eq. 0) then
            write(*,*) 'Higgs is produced on-shell'
         else if (higgsscheme .eq. 1) then
            write(*,*) 'Higgs is produced with a Breit-Wigner distribution'
         else if (higgsscheme .eq. 2) then
            write(*,*) 'Higgs is produced with a distribution a la Passarino'
         end if
      end if

    
** Deubgging messages
      if (debug .ne. 0) then
         write(*,*)'    '
         write(*,*)'SM parameters ...'
         write(*,"(A,F7.3)")'MS = ', MS
         write(*,"(A,F7.3)")'MC = ', MC
         write(*,"(A,F7.3)")'MB = ', MB
         write(*,"(A,F7.3)")'MT = ', MT
         write(*,"(A,F7.3)")'MW = ', MW
         write(*,"(A,F7.3)")'MZ = ', MZ
      
         if (model .eq. 2) then
            write(*,*)'      '
            write(*,*)'SUSY parameters ...'
            write(*,"(A,F6.3)") "Tan[beta] = ", TB
            write(*,"(A,F8.3)") "MA0 = ", MA0
            write(*,"(A,F8.3)") "MHp =", MHp
            write(*,"(A,'('F8.3,',',F8.3')')") "MUE = ", MUE
            write(*,"(A,'('F8.3,',',F8.3')')") "At = ", Af(3,3,3)
            write(*,"(A,F8.3)") "MSusy = ", MSusy
            write(*,"(A,'('F8.3,',',F8.3')')") "M_1 = ", M_1
            write(*,"(A,'('F8.3,',',F8.3')')") "M_2 = ", M_2
            write(*,"(A,'('F8.3,',',F8.3')')") "M_3 = ", M_3
         end if
      end if


** Quark masses:  Here I reset all first and second generation quark masses to
** a single parameter massQ, which controls the collinear divergences.  massQ
** is set in the file renormalisation.F.
      MsU = massQ
      MC = massQ
      MD = massQ
      MS = massQ
 
      MU2 = MsU**2D0
      MC2 = MC**2D0

      MD2 = MD**2D0
      MS2 = MS**2D0

      ME2 = ME**2D0
      MM2 = MM**2D0
      ML2 = ML**2D0


** Setting the fermion masses
      Mf(1,1) = 0D0
      Mf(1,2) = 0D0
      Mf(1,3) = 0D0

      Mf(2,1) = ME
      Mf(2,2) = MM
      Mf(2,3) = ML

      Mf(3,1) = MsU
      Mf(3,2) = MC
      Mf(3,3) = MT

      Mf(4,1) = MD
      Mf(4,2) = MS
      Mf(4,3) = MB

      do type = 1, 4
         do gen = 1, 3
            Mf2(type, gen) = Mf(type, gen)**2D0
         end do
      end do


** Setting the SUSY parameters:
      if (model .eq. 2) then

** Dealing with the higgs sector in SUSY
         Mh0 = MHiggs(1)
         MHH = MHiggs(2)
         MA0 = MHiggs(3)
         MHp = MHiggs(4)

         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2

         TB2 = TB**2D0
         CB2 = 1/(1 + TB2)
         SB2 = 1 - CB2
         CB = sqrt(CB2)
         SB = CB*TB
         S2B = 2*TB*CB2
         C2B = CB2 - SB2   
         
         SA2eff = SAeff**2D0
         CA2eff = 1 - SA2eff
         CAeff = sqrt(CA2eff)  
         S2Aeff = 2*SAeff*CAeff
         C2Aeff = CA2eff - SA2eff
         
         SABeff = SAeff*CB + CAeff*SB
         CABeff = CAeff*CB - SAeff*SB
         SBAeff = SB*CAeff - CB*SAeff
         CBAeff = CB*CAeff + SB*SAeff


** Calculating tree level Higgs sector
         if (higgsmix .eq. 2) then !  real MSSM
            MHiggs_tree(3) = MHiggs_loop(3) ! MA is the input mass
            MA0 = MHiggs_tree(3)
            MA02 = MA0**2
            MA02tree = MA02
            
            MHp2tree = MA02 + MW2
            
            avgMH2 = 0.5D0*((MA02) + MZ2)
            deltaMH2 = MA0*MZ*C2B         
         else                   ! complex MSSM
            MHiggs_tree(4) = MHiggs_loop(4) ! MHp is the input mass
            MHp = MHiggs_tree(4)
            MHp2 = MHp**2
            MHp2tree = MHp2

            MA02tree = MHp2 - MW2

            avgMH2 = 0.5D0*(MHp2 - MW2 + MZ2)
            deltaMH2 = sqrt(MHp2-MW2)*MZ*C2B
         end if

         deltaMH2 = sqrt((avgMH2 - deltaMH2)*(avgMH2 + deltaMH2))
         Mh02tree = avgMH2 - deltaMH2
         MHH2tree = avgMH2 + deltaMH2

         if (higgsmix .eq. 3) then  ! complex MSSM
            MA02tree = (MHp**2) - MW2
            MHp2tree = MHp**2
         else  !  real MSSM
            MHp2tree = (MA0**2) + MW2
            MA02tree = MA0**2
         end if

         TAtree = (-((MA0**2D0) + MZ2)*SB*CB)/
     -        (MZ2*CB2+(MA0**2D0)*SB2-Mh02tree)
         CA2tree = 1/(1 + TAtree**2D0)
         SA2tree = 1 - CA2tree
         SAtree = -sqrt(SA2tree)
         CAtree = SAtree/TAtree 
         S2Atree = 2*SAtree*CAtree
         C2Atree = CA2tree - SA2tree

         SBAtree = SB*CAtree - CB*SAtree
         CBAtree = CB*CAtree + SB*SAtree 
         SABtree = SAtree*CB + CAtree*SB
         CABtree = CAtree*CB - SAtree*SB

** Setting 'safe' tree level Higgs masses
         MHiggs_tree(1) = sqrt(Mh02tree)
         MHiggs_tree(2) = sqrt(MHH2tree)
         MHiggs_tree(3) = sqrt(MA02tree)
         MHiggs_tree(4) = sqrt(MHp2tree)
         
  
** Setting the Higgs mixing angle
         if (propLevel .eq. 0) then
* If no propagator-type corrections are used, we use the effective Higgs mixing
* angle
            SA = SAeff
            SA2 = SA2eff
            CA2 = CA2eff
            CA = CAeff  
            S2A = S2Aeff
            C2A = C2Aeff
         
            SAB = SABeff
            CAB = CABeff
            SBA = SBAeff
            CBA = CBAeff
         else  
* We use the tree level higgs mixing angle alfa when propagator-type corrections
* are used
            SA = SAtree
            SA2 = SA2tree
            CA2 = CA2tree
            CA = CAtree  
            S2A = S2Atree
            C2A = C2Atree
         
            SAB = SABtree
            CAB = CABtree
            SBA = SBAtree
            CBA = CBAtree
         end if
      

** Remember, this is only if we're working in the MSSM

** Now, sfermion parameters.  We need an internal (tree-level) calculation of 
** these if:
**   -  present in SLHA, no FH and calculating Z factors.  In this case, we
**      leave the MSf/USf as set by SLHA, but use tree level (internal) values
**      to calculate the Z factors.
**   -  not present in SLHA, no FH and calculating ew corrections or Z factors
**   -  not present in SLHA, no FH and calculating gluon fusion

         if (.not. FeynH_switch) then
            call SfermionSector(MSfcalc,USfcalc,debug)
         end if

         first = .false.
         do sfe1 = 1,2
            do type = 1,4
               do gen = 1,3
                  if (replace(MSf(sfe1,type,gen),1)) then 
                     MSf(sfe1,type,gen) = MSfcalc(sfe1,type,gen)
                     first = .true.
                  end if
                  MSf2(sfe1,type,gen)= MSf(sfe1,type,gen)**2D0
                  
                  do sfe2 = 1,2
                     if (replace(DBLE(USf(sfe1,sfe2,type,gen)),
     -                    1)) then
                        USf(sfe1,sfe2,type,gen) = 
     -                       USfcalc(sfe1,sfe2,type,gen)
                        first = .true.
                     end if
                     USfC(sfe1,sfe2,type,gen) =
     -                    DCONJG(USf(sfe1,sfe2,type,gen))
                  end do                  
               end do
            end do
         end do

** writing out sfermion info
         write(*,*)'  '
         write(*,*)'Sfermion sector ..'
         if (first) write(*,*)'(Internal calculation)'
         do gen = 1, 3
            write(*,"(A,I1,A,F9.2)")'MSf(1,',gen,') = ', 
     -           MSf(1,1,gen)
         end do
         do type = 2, 4
            do gen = 1, 3
               write(*,"(A,I1,A,I1,A,F9.5,A,F9.5)")'MSf(',type,',',
     -              gen,') = ', MSf(1,type,gen),',', MSf(2,type,gen)
               
               write(*,"(A,I1,A,I1,A,F9.5,A,F9.5,A)")'USf(1,1,',
     -              type,',',gen,') = (', DBLE(USf(1,1,type,gen)),',', 
     -              DIMAG(USf(1,1,type,gen)),')'
               write(*,"(A,I1,A,I1,A,F9.5,A,F9.5,A)")'USf(1,2,',
     -              type,',',gen,') = (', DBLE(USf(1,2,type,gen)),',', 
     -              DIMAG(USf(1,2,type,gen)),')'
               write(*,"(A,I1,A,I1,A,F9.5,A,F9.5,A)")'USf(2,1,',
     -              type,',',gen,') = (', DBLE(USf(2,1,type,gen)),',', 
     -              DIMAG(USf(2,1,type,gen)),')'
               write(*,"(A,I1,A,I1,A,F9.5,A,F9.5,A)")'USf(2,2,',
     -              type,',',gen,') = (', DBLE(USf(2,2,type,gen)),',', 
     -              DIMAG(USf(2,2,type,gen)),')'
            end do
         end do


         MUEC = DCONJG(MUE)
            
         do type = 1, 4
            do gen = 1, 3
               do sfe1 = 1, 3
                  AfC(type,gen,gen) = DCONJG(Af(type,gen,gen)) 
               end do
            end do
         end do


** The chargino and neutralino parameters are only needed if we want to 
** calculate the Z propagator matrices, or if we want to calculate ew
** corrections and are including non (s)fermions in the loops
         if ((ewcor_switch .and. (sector .gt. 2)) .or. 
     -        ((propLevel .ne. 0) .and. (.not. FeynH_switch))) then

** We need to perform an internal calculation of the chargino and neutralino
** parameters if:
**      - we're using a SLHA file (with and without FeynHiggs)

            if (slha_switch) then
               call CharginoSector(MChaCalc,UChaCalc,VChaCalc,debug)

               do cha1 = 1, 2
                  MCha(cha1) = MChaCalc(cha1)
                  
                  do cha2 = 1, 2
                     UCha(cha1,cha2) = UChaCalc(cha1,cha2)
                     VCha(cha1,cha2) = VChaCalc(cha1,cha2)
                  end do
               end do
 
** Writing out chargino information
               write(*,*)'  '
               write(*,*)'Chargino sector ..'
               write(*,"(A,F9.3,A,F9.3)")'MCha', MCha(1),',', MCha(2)
               write(*,"(A,F8.5,A,F8.5,A)")'UCha(1,1) = (',
     -              DBLE(UCha(1,1)),',',DIMAG(UCha(1,1)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'UCha(1,2) = (',
     -              DBLE(UCha(1,2)),',',DIMAG(UCha(1,2)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'UCha(2,1) = (',
     -              DBLE(UCha(2,1)),',',DIMAG(UCha(2,1)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'UCha(2,2) = (',
     -              DBLE(UCha(2,2)),',',DIMAG(UCha(2,2)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'VCha(1,1) = (',
     -              DBLE(VCha(1,1)),',',DIMAG(VCha(1,1)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'VCha(1,2) = (',
     -              DBLE(VCha(1,2)),',',DIMAG(VCha(1,2)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'VCha(2,1) = (',
     -              DBLE(VCha(2,1)),',',DIMAG(VCha(2,1)),')'
               write(*,"(A,F8.5,A,F8.5,A)")'VCha(2,2) = (',
     -              DBLE(VCha(2,2)),',',DIMAG(VCha(2,2)),')'

** And now the neutralinos, which work in the same way as the charginos  
               call NeutralinoSector(MNeuCalc,ZNeuCalc,debug)

               do neu1 = 1, 4
                  MNeu(neu1) = MNeuCalc(neu1)
               
                  do neu2 = 1, 4
                     ZNeu(neu1,neu2) = ZNeuCalc(neu1,neu2)
                  end do
               end do

** Writing out neutralino information
               write(*,*)'  '
               write(*,*)'Neutralino sector ..'
               if (.not. first) write(*,*)'(Internal calculation)'
               do neu1 = 1, 4
                  write(*,"(A,I1,A,F9.3)")'MNeu(', neu1,') = ', 
     -                 MNeu(neu1)
                  do neu2 = 1, 4
                     write(*,"(A,I1,A,I1,A,F8.5,A,F8.5,A)")'ZNeu(',neu1,
     -                    ',',neu2,') = (', DBLE(ZNeu(neu1,neu2)),',',
     -                    DIMAG(ZNeu(neu1,neu2)),')'
                  end do
               end do


            end if ! end of if slha_switch

         end if  ! end of need chargino/neutralino parameters 


         do cha1 = 1, 2
            MCha2(cha1) = MCha(cha1)**2
            
            do cha2 = 1, 2
               UChaC(cha1,cha2) = DCONJG(UCha(cha1,cha2))
               VChaC(cha1,cha2) = DCONJG(VCha(cha1,cha2))
            end do
         end do
         
         do neu1 = 1, 4
            MNeu2(neu1) = MNeu(neu1)**2   
            
            do neu2 = 1, 4
               ZNeuC(neu1,neu2) = DCONJG(ZNeu(neu1,neu2))
            end do
         end do

      end if   ! end of SUSY only


** Setting to tree level Higgs masses for the calculation of the renormalisation
** constants
      if ((model .eq. 2) .and. (MHloops .eq. 0)) then
         do type = 1, 4
            MHiggs(type) = MHiggs_tree(type)
            MHiggs2(type) = MHiggs(type)**2
         end do
         
         Mh0 = MHiggs(1)
         MHH = MHiggs(2)
         MA0 = MHiggs(3)
         MHp = MHiggs(4)

         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2

         MH = MHiggs(HiggsType)
         MH2 = MHiggs2(HiggsType)
      end if
** End of setting tree level Higgs masses  


** Calculating the renormalization constants. 
#ifdef WITH_LT      
      if (ewcor_switch) call ren_con(charge)
#endif

** Next: outputting the SUSY Higgs info and - if we're not using FeynHiggs - 
** calculating the Z propagator matrix if requested
      if (model .eq. 2) then
         write(*,*)'  '

         if (FeynH_switch) then
            write(*,*)'FeynHiggs calculation of Higgs parameters:'
         else if (slha_switch) then  
** If we're not using FH, we need to caclulate the Z propagator matrix, unless
** we're using an effective (higgs-mixing) alpha approximation
            if (propLevel .ne. 0) then
               if (higgsmix .eq. 3) then
                  write(*,*)'Sorry, we can not calculate the Z matrix'
                  write(*,*)'internally in the complex MSSM yet.'
                  write(*,*)'You can either use an effective Higgs'
                  write(*,*)'mixing angle alpha (propLevel = 0),'
                  write(*,*)'or link to FeynHiggs.'
                  stop
               end if
#ifdef WITH_LT
               call slha_zmatrix(0) ! argument > 0 gives debugging info
#else
               write(*,*)'Sorry, we cannot calculate the Higgs-'
               write(*,*)'propagator-type corrections internally'
               write(*,*)'without LoopTools.  Either link to LoopTools'
               write(*,*)'or FeynHiggs, or use propLevel = 0.'
               stop
#endif
            end if  !  End of propLevel = 0
            write(*,*)'SLHA Higgs parameters:'
         end if  

** writing out Higgs info
         do cur = 1, 3
            write(*,"(A,I1,A,F8.3)")'MHiggs(',cur,') = ', 
     -           MHiggs_loop(cur)
            if (propLevel .ne. 0) then
               write(*,"(A,I1,A,'('F10.7,','F10.7')')")
     -              'propagator corn: Hgs(',cur,',h) = ', Hgs(cur,1)
               write(*,"(A,I1,A,'('F10.7,','F10.7')')")
     -              'propagator corn: Hgs(',cur,',H) = ', Hgs(cur,2)
               write(*,"(A,I1,A,'('F10.7,','F10.7')')")
     -              'propagator corn: Hgs(',cur,',A) = ', Hgs(cur,3)
            end if
         end do
         write(*,"(A,F8.5)")'Effective sin(beta-alpha) = ', SBAeff
      end if  !  end of MSSM only


** Resetting to corrected Higgs masses
      if ((model .eq. 2) .and. (MHloops .eq. 0)) then
         do type = 1, 4
            MHiggs(type) = MHiggs_loop(type)
            MHiggs2(type) = MHiggs(type)**2
         end do
         
         Mh0 = MHiggs_loop(1)
         MHH = MHiggs_loop(2)
         MA0 = MHiggs_loop(3)
         MHp = MHiggs_loop(4)
         
         Mh02 = Mh0**2
         MHH2 = MHH**2
         MA02 = MA0**2
         MHp2 = MHp**2

         MH = MHiggs(HiggsType)
         MH2 = MHiggs2(HiggsType)
      end if
** End of resetting Higgs masses

      
** These are the values for the "renormalisation constant" of the L/R couplings
** can take for up-type and down-type fermions
      if (ewcor_switch) then
         dgR3 = gR3*(dZet1 + su2sign/(CW2*SW)*dSW1)
         dgR4 = gR4*(dZet1 + su2sign/(CW2*SW)*dSW1)
         dgL3 = su2sign/(2*SW*CW) * (dZet1 + (SW2 - CW2)/
     -        (su2sign*CW2*SW)*dSW1) + dgR3
         dgL4 = -su2sign/(2*SW*CW) *(dZet1 + (SW2 - CW2)/
     -        (su2sign*CW2*SW)*dSW1) + dgR4
** and for fermions:
         dgR3fm = gR3*(dZet1fm + su2sign/(CW2*SW)*dSW1fm)
         dgR4fm = gR4*(dZet1fm + su2sign/(CW2*SW)*dSW1fm)
         dgL3fm = su2sign/(2*SW*CW) * (dZet1fm + (SW2 - CW2)/
     -        (su2sign*CW2*SW)*dSW1fm) + dgR3fm
         dgL4fm = -su2sign/(2*SW*CW) *(dZet1fm + (SW2 - CW2)/
     -        (su2sign*CW2*SW)*dSW1fm) + dgR4fm
      else
         dgR3 = 0D0
         dgR4 = 0D0
         dgL3 = 0D0
         dgL4 = 0D0
      end if


** Setting the mass for each current when all externals are quarks (rather than
** antiquarks).  See notebook 7, page 1 for notation.  Labelling: [particle, 
** current].
      do cur = 1, 6
         if (cur .eq. 1) then
            massCur(1,cur) = MsU
            massCur(2,cur) = MD
            massCur(3,cur) = MC
            massCur(4,cur) = MS
         else if (cur .eq. 2) then
            massCur(1,cur) = MsU
            massCur(2,cur) = MD
            massCur(3,cur) = MS
            massCur(4,cur) = MC
         else if (cur .eq. 3) then
            massCur(1,cur) = MD
            massCur(2,cur) = MsU
            massCur(3,cur) = MC
            massCur(4,cur) = MS
         else if (cur .eq. 4) then
            massCur(1,cur) = MD
            massCur(2,cur) = MsU
            massCur(3,cur) = MS
            massCur(4,cur) = MC         
         else if (cur .eq. 5) then
            massCur(1,cur) = MsU
            massCur(2,cur) = MD
            massCur(3,cur) = MS
            massCur(4,cur) = MC
         else if (cur.eq. 6) then
            massCur(1,cur) = MD
            massCur(2,cur) = MsU
            massCur(3,cur) = MC
            massCur(4,cur) = MS
         end if

         massCur(5,cur) = MH
         
      end do


** Setting the charges of each current [labelling: (current, particle)] for all
** quarks (qk = 1).  
      do cur = 1, 6
         if (cur .eq. 1) then
            chargeCur(cur,1) = Qu
            chargeCur(cur,2) = Qd
            chargeCur(cur,3) = Qu
            chargeCur(cur,4) = Qd
         else if (cur .eq. 2 .or. cur .eq. 5) then
            chargeCur(cur,1) = Qu
            chargeCur(cur,2) = Qd
            chargeCur(cur,3) = Qd
            chargeCur(cur,4) = Qu
         else if (cur .eq. 3 .or. cur .eq. 6) then
            chargeCur(cur,1) = Qd
            chargeCur(cur,2) = Qu
            chargeCur(cur,3) = Qu
            chargeCur(cur,4) = Qd
         else if (cur .eq. 4) then
            chargeCur(cur,1) = Qd
            chargeCur(cur,2) = Qu
            chargeCur(cur,3) = Qd
            chargeCur(cur,4) = Qu
         end if
         chargeCur(cur,5) = 0
      end do


** Setting the generic values of the susy couplings (se ntbk, 29-7-10)
      if (model .eq. 2) then
** Setting values for neutralino couplings: see ntbk 29-7-10
         do neu1 = 1, 4
            do gen = 1, 2
               SNeu1u(neu1,gen) = Mf(3,gen)*ZNeu(neu1,4)/SB
               SNeu1d(neu1,gen) = Mf(4,gen)*ZNeu(neu1,3)/CB
            end do
         end do

         SNeu2u = 1D0
         SNeu2d = -1D0

         SNeu3u = 2D0
         SNeu3d = -1D0
** End of setting values for neutralino couplings


** Setting values for chargino couplings: see ntbk 29-7-10
         do cha1 = 1, 2
            do gen = 1, 2
               SCha1u(cha1,gen) = Mf(3,gen)*VCha(cha1,2)/SB
               SCha1d(cha1,gen) = Mf(4,gen)*UCha(cha1,2)/CB
            end do
            SCha2u(cha1) = UCha(cha1,1)
            SCha2d(cha1) = VCha(cha1,1)
         end do
** End of setting values for chargino couplings


** Setting values for chargino-neutralino-W couplings
         do cha1 = 1, 2
            do neu1 = 1, 4
               cnw(1,cha1,neu1) = VChaC(cha1,1)*ZNeu(neu1,2) -
     -              (VChaC(cha1,2)*ZNeu(neu1,4)/sqrt(2D0))
               cnw(2,cha1,neu1) = UCha(cha1,1)*ZNeuC(neu1,2) +
     -              (UCha(cha1,2)*ZNeuC(neu1,3)/sqrt(2D0))
               cnw(3,cha1,neu1) = VCha(cha1,1)*ZNeuC(neu1,2) -
     -              (VCha(cha1,2)*ZNeuC(neu1,4)/sqrt(2D0))
               cnw(4,cha1,neu1) = UChaC(cha1,1)*ZNeu(neu1,2) +
     -              (UChaC(cha1,2)*ZNeu(neu1,3)/sqrt(2D0))
            end do
         end do
** End of setting values for chargino-neutralino-W couplings


** Setting values for chargino-chargino-Z couplings
         do cha1 = 1, 2
            do cha2 = 1, 2
               zcc(cha1,cha2,1) = -((UChaC(cha1,2)*UCha(cha2,2)/2D0) +
     -              UChaC(cha1,1)*UCha(cha2,1))
               zcc(cha1,cha2,2) = -((VChaC(cha1,2)*VCha(cha2,2)/2D0) +
     -              VChaC(cha1,1)*VCha(cha2,1))
            end do
         end do
** End of setting values for chargino-chargino-Z couplings


** Setting values for higgs-squark-squark couplings
         SAqku = 1D0/TB
         SAqkd = TB

         do gen = 1, 2 
            Shqku(gen) = -(Af(3,gen,gen)*CA + MUEC*SA)/SB
            Shqkd(gen) = (Af(4,gen,gen)*SA + MUEC*CA)/CB
         end do
         Shqk1u = -CA/SB
         Shqk1d = SA/CB

         do gen = 1, 2
            SHHu(gen) = (CA*MUEC - SA*Af(3,gen,gen))/SB
            SHHd(gen) = (SA*MUEC - CA*Af(4,gen,gen))/CB
         end do
         SHH1u = -SA/SB
         SHH1d = -CA/CB
** End of setting values for higgs-squark-squark couplings         
         

      end if  ! end of setting values of generic SUSY couplings


**** END of parameters section  ****


      end


*******************************************************************************
*******************************************************************************

*** This subroutine sets the parameters according to a benchmark scenario, 
*** which is set by the parameter 'bench' (set in susy.dat).  Note that if
*** slha_switch is set to .true., this overrides the choice of 'bench'.

*** Benchmarks 1-4 from: hep-ph/0202167
*** CPX scenario from: hep-ex/0602042
*** SPS parameters from:  http://www.ippp.dur.ac.uk/~georg/sps/

      subroutine bench_set(higgsmix, TB, MH, MA0, MHp,
     &     M3SL, M3SE, M3SQ, M3SU, M3SD,
     &     M2SL, M2SE, M2SQ, M2SU, M2SD,
     &     M1SL, M1SE, M1SQ, M1SU, M1SD,
     &     MUE,
     &     Atau, At, Ab, Amu, Ac, As, Ae, Au, Ad,
     &     M_1, M_2, M_3, MSusy, slha_switch, bench, debug)
          use readinput

      implicit none

** Higgs parameters
      double precision TB, MA0, MHp, MH

** Soft-susy breaking parameters in the sfermion sector
      double precision MSusy
      double precision M3SL, M2SL, M1SL
      double precision M3SE, M2SE, M1SE
      double precision M3SQ, M2SQ, M1SQ
      double precision M3SU, M2SU, M1SU
      double precision M3SD, M2SD, M1SD

** Trilinear coupling parameters      
      double complex Atau, At, Ab
      double complex Amu, Ac, As
      double complex Ae, Au, Ad
        
** Higgs mixing parameter and gaugino mass parameters
      double complex MUE, M_1, M_2, M_3

** Debuggin flag and flag controlling which benchmark we're working in
      integer debug, bench

** FeynHiggs flag, controls whether real or complex parameters are considered
      integer higgsmix

** Switch determining whether we're using a SLHA file or not
      logical slha_switch

** Function which determines whether a parameter has already been set with a 
** reasonable value
      logical replace
      external replace


** Reading input from file susy.dat: default values are SPS1a
      call loadfile("susy.dat",.true.)

      write(*,*)" "
      write(*,*)"              SUSY physics parameters"
      write(*,*)"-----------------------------------------------"

      if (.not. slha_switch) then
         call read_Int("BENCH", bench, 1)

         if ((bench .lt. 0) .or. (bench .gt. 19) .or. 
     -        ((bench .gt. 5) .and. (bench .lt. 10))) then
            write(*,*)'REPLACE choice of input BENCH!'
            write(*,*)'Please input BENCH = 0..5 or 10..19'
            stop
         end if
      end if

      if (replace(TB,0)) call read_Real("TANB", TB, 10D0)
      if (replace(MA0,0)) call read_Real("MASSA0", MA0, 400D0)
      if (replace(MHp,0)) call read_Real("MASSHP", MHp, -1D0)

      if (replace(DBLE(higgsmix),0)) then
         if (MA0 .lt. 0D0) then ! complex MSSM
            higgsmix = 3
         else if (MHp .lt. 0D0) then ! real MSSM
            higgsmix = 2
         end if
      end if
      if ((MA0 .lt. 0D0) .and. (MHp .lt. 0D0)) then
         write(*,*)'Please enter a value for either the cp-odd Higgs'
         write(*,*)'mass or the charged Higgs mass.'
         stop
      else
         if (replace(DBLE(higgsmix),0)) higgsmix = 2
      end if

** user defined scenario: reading input - default values are SPS1a
      if (slha_switch .or. (bench .eq. 0)) then
         if (replace(M3SQ,0)) call read_Real("M3SQ", M3SQ,495.9D0)
         if (replace(M3SU,0)) call read_Real("M3SU", M3SU,424.8D0)
         if (replace(M3SD,0)) call read_Real("M3SD", M3SD,516.9D0)
         if (replace(M3SL,0)) call read_Real("M3SL", M3SL,195.8D0)
         if (replace(M3SE,0)) call read_Real("M3SE", M3SE,133.6D0)
         if (replace(M2SQ,0)) call read_Real("M2SQ", M2SQ,539.9D0)
         if (replace(M2SU,0)) call read_Real("M2SU", M2SU,521.7D0)
         if (replace(M2SD,0)) call read_Real("M2SD", M2SD,519.5D0)
         if (replace(M2SL,0)) call read_Real("M2SL", M2SL,196.6D0)
         if (replace(M2SE,0)) call read_Real("M2SE", M2SE,136.2D0)

         if (replace(DBLE(mue),0)) call read_cplx("MUE", MUE, (352.4D0,0D0))
         if (replace(DBLE(M_2),0)) call read_cplx("M_2", M_2, (192.7D0,0D0))
         if (replace(DBLE(M_1),0)) call read_cplx("M_1", M_1, (99.1D0,0D0))
         if (replace(DBLE(M_3),0)) call read_cplx("M_3", M_3, (595.2D0,0D0))

         if ((slha_switch) .and. .not. (replace(DBLE(At),0) .or. 
     -        replace(DBLE(Ab),0) .or. replace(DBLE(Atau),0))) then
            if (replace(DBLE(At),0)) then
               if (replace(DBLE(Ab),0)) then ! neither At nor Ab are set
                  At = Atau
                  Ab = Atau
               else if (replace(DBLE(Atau),0)) then ! neither At nor Atau set
                  At = Ab
                  Atau = Ab
               else ! only At is not set
                  At = Ab
               end if
            else if (replace(DBLE(Ab),0)) then
               if (replace(DBLE(Atau),0)) then ! neither Ab nor Atau are set
                  Ab = At
                  Atau = At
               else if (replace(DBLE(At),0)) then ! neither Ab nor At are set
                  At = Atau
                  Ab = Atau
               else ! only Ab is not set
                  Ab = At
               end if
            else if (replace(DBLE(Atau),0)) then
               if (replace(DBLE(At),0)) then ! neither Atau nor At are set
                  At = Ab
                  Atau = Ab
               else if (replace(DBLE(Ab),0)) then ! neither Atau nor Ab are set
                  Ab = At
                  Atau = At
               else ! only Atau is not set
                  Atau = Ab
               end if
            end if                  
         else
            if (replace(DBLE(At),0)) call read_cplx("AT", At, (-510.0D0,0D0))
            if (replace(DBLE(Ab),0)) call read_cplx("AB", Ab, (-772.7D0,0D0))
            if (replace(DBLE(Atau),0)) call read_cplx("ATAU", Atau,
     -           (-254.2D0,0D0))
         end if
   
         if (replace(DBLE(Amu),0)) Amu = Atau
         if (replace(DBLE(Ac),0)) Ac = At
         if (replace(DBLE(As),0)) As = Ab
         if (replace(DBLE(Ae),0)) Ae = Atau
         if (replace(DBLE(Au),0)) Au = At
         if (replace(DBLE(Ad),0)) Ad = Ab

         MSusy = M3SQ
         if (replace(M1SL,0)) M1SL = M2SL
         if (replace(M1SE,0)) M1SE = M2SE
         if (replace(M1SQ,0)) M1SQ = M2SQ
         if (replace(M1SU,0)) M1SU = M2SU
         if (replace(M1SD,0)) M1SD = M2SD 

         MSusy = M3SQ

      end if
    

      call closefile
      write(*,*)'  '
      

** completely user input scenario!
      if ((bench .eq. 0) .or. slha_switch) return

        
      if (bench .eq. 1) then
         
         write(*,*) "Mh max scenario"
         
         MSusy = 1000D0
         
         MUE = 200D0
         
         Atau = (2D0*MSusy)+(MUE/TB)
         
         M_2 = 200D0
         M_3 = 0.8D0*MSusy       
         
         higgsmix = 2
         
         
      else if (bench .eq. 2) then
         
         write(*,*) "no-mixing scenario"
         
         MSusy = 2000D0
         
         MUE = 200D0
         
         Atau = MUE/TB
         
         M_2 = 200D0
         M_3 = 0.8D0*MSusy       
         
         higgsmix = 2
         
         
      else if (bench .eq. 3) then
         
         write(*,*) "gluophobic scenario"
         
         MSusy = 350D0
         
         MUE = 300D0
         
         Atau = -750D0+(MUE/TB)
         
         M_2 = 300D0
         M_3 = 500D0        
         
         higgsmix = 2
         
         
      else if (bench .eq. 4) then
         
         write(*,*) "small alpha-eff scenario"
         
         MSusy = 800D0
         
         MUE = 2.5D0*MSusy
         
         Atau = -1100D0+(MUE/TB)
         
         M_2 = 500D0
         M_3 = 500D0       
         
         higgsmix = 2
         
         
      else if (bench .eq. 5) then
         
         write(*,*) "CPX scenario"
         
         MSusy = 500D0
         
         MUE = 2000D0
         
         Atau = (0D0,900D0)
         
         M_2 = 200D0
         M_3 = 1000D0*Exp((0D0,1.5708D0))       
         
         higgsmix = 3

      end if
      

      M_1 = (0D0,0D0)

      M3SL = MSusy
      M2SL = M3SL
      M1SL = M2SL
      M3SE = MSusy
      M2SE = M3SE
      M1SE = M2SE
      M3SQ = MSusy
      M2SQ = M3SQ
      M1SQ = M2SQ
      M3SU = MSusy
      M2SU = M3SU
      M1SU = M2SU
      M3SD = MSusy
      M2SD = M3SD
      M1SD = M2SD
      
      At = Atau
      Ab = Atau
      Amu = Atau
      Ac = At
      As = Ab
      Ae = Atau
      Au = At
      Ad = Ab


      if (bench .eq. 10) then

         write(*,*)'benchmark SPS1a ...' 

         M3SQ = 495.9D0
         M3SU = 424.8D0
         M3SD = 516.9D0
         M3SL = 195.8D0
         M3SE = 133.6D0
         M2SQ = 539.9D0
         M2SU = 521.7D0
         M2SD = 519.5D0
         M2SL = 196.6D0
         M2SE = 136.2D0

         At = -510.0D0
         Ab = -772.7D0
         Atau = -254.2D0

         MA0 = 393.6D0
         MHp = -1
         TB = 10D0

         MUE = 352.4D0
         M_2 = 192.7D0
         M_1 = 99.1D0
         M_3 = 595.2D0

      else if(bench .eq. 11) then

         write(*,*)'benchmark SPS1b ...'  

         M3SQ = 762.5D0
         M3SU = 670.7D0
         M3SD = 780.3D0
         M3SL = 323.8D0
         M3SE = 218.6D0
         M2SQ = 836.2D0
         M2SU = 807.5D0
         M2SD = 803.9D0
         M2SL = 334.0D0
         M2SE = 248.3D0

         At = -729.3D0
         Ab = -987.4D0
         Atau = -195.8D0

         MA0 =  525.5D0
         MHp = -1
         TB = 30D0

         MUE = 495.6D0
         M_2 = 310.9D0
         M_1 = 162.8D0
         M_3 = 916.1D0
         
      else if (bench .eq. 12) then

         write(*,*)'benchmark SPS2 ...'  

         M3SQ = 1295.3D0
         M3SU = 998.5D0
         M3SD = 1519.9D0
         M3SL = 1449.6D0
         M3SE = 1438.9D0
         M2SQ = 1533.6D0
         M2SU = 1530.5D0
         M2SD = 1530.3D0
         M2SL = 1455.6D0
         M2SE = 1451.0D0

         At = -563.7D0
         Ab = -797.2D0
         Atau = -187.8D0

         MA0 =  1443.0D0
         MHp = -1
         TB = 10D0

         MUE = 124.8D0
         M_2 = 234.1D0
         M_1 = 120.4D0
         M_3 = 784.4D0
         
      else if (bench .eq. 13) then

         write(*,*)'benchmark SPS3 ...'  

         M3SQ = 760.7D0
         M3SU = 661.2D0
         M3SD = 785.6D0
         M3SL = 282.4D0
         M3SE = 170D0
         M2SQ = 818.3D0
         M2SU = 792.6D0
         M2SD = 788.9D0
         M2SL = 283.3D0
         M2SE = 173.0D0

         At = -733.5D0
         Ab = -1042.2D0
         Atau = -246.1D0

         MA0 =  572.4D0
         MHp = -1
         TB = 10D0

         MUE = 508.6D0
         M_2 = 311.4D0
         M_1 = 162.8D0
         M_3 = 914.3D0
         
      else if (bench .eq. 14) then

         write(*,*)'benchmark SPS4 ...'  

         M3SQ = 640.1D0
         M3SU = 556.8D0
         M3SD = 673.4D0
         M3SL = 394.7D0
         M3SE = 289.5D0
         M2SQ = 732.2D0
         M2SU = 716.0D0
         M2SD = 713.9D0
         M2SL = 445.9D0
         M2SE = 414.2D0

         At = -552.2D0
         Ab = -729.5D0
         Atau = -102.3D0

         MA0 =  404.4D0
         MHp = -1
         TB = 50D0

         MUE = 377D0
         M_2 = 233.2D0
         M_1 = 120.8D0
         M_3 = 721D0
         

      else if (bench .eq. 15) then

         write(*,*)'benchmark SPS5 ...'  

         M3SQ = 535.2D0
         M3SU = 360.5D0
         M3SD = 620.5D0
         M3SL = 250.1D0
         M3SE = 180.9D0
         M2SQ = 643.9D0
         M2SU = 625.4D0
         M2SD = 622.9D0
         M2SL = 252.2D0
         M2SE = 186.8D0

         At = -905.6D0
         Ab = -1671.4D0
         Atau = -1179.3D0

         MA0 =  693.9D0
         MHp = -1
         TB = 5D0

         MUE = 639.8D0
         M_2 = 234.6D0
         M_1 = 121.4D0
         M_3 = 710.3D0
         
      else if (bench .eq. 16) then

         write(*,*)'benchmark SPS6 ...'  

         M3SQ = 591.2D0
         M3SU = 517.0D0
         M3SD = 619.0D0
         M3SL = 259.7D0
         M3SE = 230.5D0
         M2SQ = 641.3D0
         M2SU = 629.3D0
         M2SD = 621.8D0
         M2SL = 260.7D0
         M2SE = 232.8D0

         At = -570.0D0
         Ab = -811.3D0
         Atau = -213.4D0

         MA0 =  463.0D0
         MHp = -1
         TB = 10D0

         MUE = 393.9D0
         M_2 = 232.1D0
         M_1 = 195.9D0
         M_3 = 708.5D0

      else if (bench .eq. 17) then

         write(*,*)'benchmark SPS7 ...'  

         M3SQ = 836.3D0
         M3SU = 780.1D0
         M3SD = 826.9D0
         M3SL = 256.8D0
         M3SE = 117.6D0
         M2SQ = 861.3D0
         M2SU = 831.3D0
         M2SD = 828.6D0
         M2SL = 257.2D0
         M2SE = 119.7D0

         At = -319.4D0
         Ab = -350.5D0
         Atau = -39.0D0

         MA0 =  377.9D0
         MHp = -1
         TB = 15D0

         MUE = 300D0
         M_2 = 326.8D0
         M_1 = 168.6D0
         M_3 = 926.0D0

      else if (bench .eq. 18) then

         write(*,*)'benchmark SPS8 ...'  

         M3SQ = 1042.7D0
         M3SU = 952.7D0
         M3SD = 1025.5D0
         M3SL = 352.8D0
         M3SE = 167.2D0
         M2SQ = 1081.6D0
         M2SU = 1033.8D0
         M2SD = 1029.0D0
         M2SL = 353.5D0
         M2SE = 170.4D0

         At = -296.7D0
         Ab = -330.3D0
         Atau = -36.7D0

         MA0 =  514.5D0
         MHp = -1
         TB = 15D0

         MUE = 398.3D0
         M_2 = 271.8D0
         M_1 = 140.0D0
         M_3 = 820.5D0

      else if (bench .eq. 19) then

         write(*,*)'benchmark SPS9 ...'  

         M3SQ = 1111.6D0
         M3SU = 1003.2D0
         M3SD = 1231.7D0
         M3SL = 307.4D0
         M3SE = 281.2D0
         M2SQ = 1219.2D0
         M2SU = 1227.9D0
         M2SD = 1237.6D0
         M2SL = 316.2D0
         M2SE = 300.0D0

         At = -350.3D0
         Ab = 216.4D0
         Atau = 1162.4D0

         MA0 =  911.7D0
         MHp = -1
         TB = 10D0

         MUE = 869.9D0
         M_2 = -175.5D0
         M_1 = -550.6D0
         M_3 = 1275.2D0

      end if


      if (bench .ge. 6) then         
* For the SPS points, trilinear coupling parameters of first and second 
* generation are set to the values of the third generation
         Amu = Atau
         Ac = At
         As = Ab
         Ae = Atau
         Au = At
         Ad = Ab

         MSusy = M3SQ
         M1SL = M2SL
         M1SE = M2SE
         M1SQ = M2SQ
         M1SU = M2SU
         M1SD = M2SD 
         
         higgsmix = 2
      end if 

      
      end


********************************************************************************
********************************************************************************

*** A function which determines whether a parameter needs to be set from the 
*** .dat input files, or via an internal calculation.  

      logical function replace(para,calc)
          use globalvars, only: ldoblha

      implicit none

** The parameter which (perhaps) needs to be set
      double precision para

** What type of parameter it is:
**   0: an input parameter (4-3-11 for full list)
**   1: sfermion masses and mixings
      integer calc 

#include "global.inc"
#include "mssm.inc"
#include "BLHAhelper.inc"
#include "VBFNLO/amplitudes/ggf/ggf_process.inc"

      
      replace = .false.


      if (calc .eq. 0) then

         if (slha_switch .or. ldoblha) then
            if (para .eq. -999D0) replace = .true.
         else
            replace = .true.
         end if

      else if (calc .eq. 1) then ! sfermion masses and mixings

         if ((model .eq. 2) .and. (.not. FeynH_switch)) then

            if ((procID .ge. 4100) .and. (procID .le. 4107)) then !gluon fusion
               if (higgs_mix .eq. 2) then
                  replace = .true.
               else if (.not. slha_switch) then
                  replace = .true.
               else if (para .eq. -999D0) then
                  replace = .true.
               end if
            else if (propLevel .ne. 0) then
               replace = .true.
            else if (ewcor_switch) then
               if (.not. slha_switch) then
                  replace = .true.
               else if (para .eq. -999D0) then
                  replace = .true.
               end if
            end if
                  
         end if

      else if (calc .eq. 2) then ! chargino and neutralino parameters
         
         if (FeynH_switch) then
            if (para .eq. -999D0) replace = .true.
         else
            replace = .true.
         end if

      else if (calc .eq. 3) then ! Higgs widths

c$$$         if (FeynH_switch) then
c$$$            if (para .eq. -999D0) replace = .true.
c$$$         else
c$$$            if (slha_switch) then
c$$$               if (para .eq. -999D0) replace = .true.
c$$$            else
c$$$               replace = .true.
c$$$            end if
c$$$         end if
         if (para .eq. -999d0) replace = .true.

      else

         write(*,*)'SOPHY! You have a problem with REPLACE.'
         stop

      end if         


      end


********************************************************************************
********************************************************************************

*** Subroutine to initialise branching ratios and total widths to -999

      subroutine clearwidths

      implicit none


      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

      double precision BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XGH2, sinba, cosba
      COMMON /BRANCH2/ BH2WW,BH2ZZ,BH2GG,BH2TT,BH2BB,BH2CC,BH2TAU,BH2MU,
     &                 BH2GAM, BH2GAMZ, XGH2, sinba, cosba


** widths:
      xgw = -999d0
      xgz = -999d0
      xgh = -999d0
      gamt = -999d0


** w branching ratios
      bwne = -999d0
      bwud = -999d0
      bwtb = -999d0


** z branching ratios
      bznn = -999d0
      bzee = -999d0
      bzuu = -999d0
      bzdd = -999d0
      bztt = -999d0


** higgs branching ratios
      bhww = -999d0
      bhzz = -999d0
      bhgg = -999d0
      bhtt = -999d0
      bhbb = -999d0
      bhcc = -999d0
      bhtau = -999d0
      bhmu = -999d0
      bhgam = -999d0
      bhgamz = -999d0

** Higgs2 branching ratios
      BH2WW   = -999d0
      BH2ZZ   = -999d0
      BH2GG   = -999d0
      BH2TT   = -999d0
      BH2BB   = -999d0
      BH2CC   = -999d0
      BH2TAU  = -999d0
      BH2MU   = -999d0
      BH2GAM  = -999d0
      BH2GAMZ = -999d0
      XGH2    = -999d0

      end


********************************************************************************
********************************************************************************

*** This subroutine checks that the relevant branching ratio for a particular
*** process is not greater than one.  Only resonant Higgs bosons are considered.

      subroutine checkBR

      implicit none

#include "global.inc"
#include "process.inc"


      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT


      SELECT CASE(procID)
      CASE(Hjj_AA, WPH_AA, WPHJ_AA, WMH_AA, WMHJ_AA, Hjjj_AA, HAjj_AA, GFHjj_AA, AA, GFAA)
         if (BHGAM .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> 2 photons'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHGAM
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(WPAA, WMAA, ZAA, ZnAA)
         if (BHGAM .gt. 1D0 .and. with_anom) then      ! relevant due to anomalous AZH, AAH couplings
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> 2 photons'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHGAM
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(ZA, GFZA)
         if (BHGAMZ .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> photon + Z'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHGAMZ
            stop
         end if
      CASE(WPZA, WMZA, WPhadZA, WPZhadA, WMhadZA, WMZhadA,
     &     ZAA, ZnAA)
         if (BHGAMZ .gt. 1D0 .and. with_anom) then      ! relevant due to anomalous AZH, AAH couplings
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> photon + Z'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHGAMZ
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(Hjj_mu, WPH_mu, WPHJ_mu, WMH_mu, WMHJ_mu, Hjjj_mu, HAjj_mu, GFHjj_mu)
         if (BHMU .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> 2 muons'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHMU
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(Hjj_tau, WPH_tau, WPHJ_tau, WMH_tau, WMHJ_tau, Hjjj_tau, HAjj_tau, GFHjj_tau)
         if (BHTAU .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> tau+ tau-'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHTAU
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(Hjj_bbar, WPH_bbar, WPHJ_bbar, WMH_bbar, WMHJ_bbar, Hjjj_bbar, HAjj_bbar, GFHjj_bbar)
         if (BHBB .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> b bbar'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHBB
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(Hjj_WW, WPH_WW, WPHJ_WW, WMH_WW, WMHJ_WW, Hjjj_WW, HAjj_WW, GFHjj_WW, 
     -        WPWMjj, 
     -        WWZ, WWWP, WWWM, WW, GFWW,
     -        WPhadWMlep,WPlepWMhad,
     -        GFWPhadWMlep,GFWPlepWMhad,
     -        WPhadWMjj, WPWMhadjj, WPhadWMZ, WPWMhadZ, WWZhad,
     -        WPhadWMWP, WPWMhadWP, WMhadWPWM, WMWPhadWM,
     -        Hjj_WPhadWM, Hjj_WPWMhad)
         if (BHWW .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> W+ W-'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHWW
            stop
         end if
      CASE(WWA, WPhadWMA, WPWMhadA)
         if (BHWW .gt. 1D0 .and. with_anom) then      ! relevant due to anomalous AZH, AAH couplings
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> W+ W-'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHWW
            stop
         end if
      END SELECT

      SELECT CASE(procID)
      CASE(Hjj_ZZ_ll, Hjj_ZZ_lnu, Hjjj_ZZ_ll, Hjjj_ZZ_lnu, 
     &    WPH_ZZ_ll, WPHJ_ZZ_ll, WMH_ZZ_ll, WMHJ_ZZ_ll,
     &    WPH_ZZ_lnu, WPHJ_ZZ_lnu, WMH_ZZ_lnu, WMHJ_ZZ_lnu,
     -        HAjj_ZZ_ll, HAjj_ZZ_lnu, GFHjj_ZZ_ll, GFHjj_ZZ_lnu,  
     -        ZZjj_ll, ZZjj_lnu, 
     -        ZZWP, ZZWM, ZZZ, ZZ, GFZZ, ZZhad, GFZZhad,
     -        ZZhadWP, ZZWPhad, ZZhadWM, ZZWMhad, ZZhadjj, ZZZhad, Hjj_ZZhad)
         if (BHZZ .gt. 1D0) then
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> Z Z'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHZZ
            stop
         end if
      CASE(ZZA, ZZhadA, ZZnA)
         if (BHZZ .gt. 1D0 .and. with_anom) then      ! relevant due to anomalous AZH, AAH couplings
            write(*,*)'The branching ratio:'
            write(*,*)'    Higgs --> Z Z'
            write(*,*)'is greater than one.'
            write(*,*)'branching ratio =', BHZZ
            stop
         end if
      END SELECT


      end



