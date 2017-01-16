********************************************************************************
********************************************************************************
*** 14th January 2011                                                        ***
*** slha_read.F                                                              ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contains subroutines to read all parameters in from a SLHA     ***
*** file and calculate the SUSY parameters (e.g. sfermion masses) if         ***
*** FeynHiggs is not being used.  It's mostly stolen from Thomas Hahn's      ***
*** SLHALib-2.2 (hep–ph/0605049).                                            ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine reads in all parameters and decays from a SLHA file.  If a 
*** parameter is not present it sets its value to -999, which later signifies
*** to the code that this parameter still needs to be set.

      subroutine readSLHA(debug)

      implicit none

#include "SLHA.h"
#include "mssm.inc"
#include "koppln.inc"
#include "global.inc"

      integer error
      double complex slhadata(nslhadata)

** do loop parameters
      integer type, gen1, gen2, sfe1, sfe2, cha1, cha2, neu1, neu2

** inverse of em coupling Alfa at MZ, Alfa(0)
      double precision invAlfaMZ, Alfa0

** effective Higgs mixing angle 
      double precision Aeff

** sfermion mixing matrices in 'extended' format, and dummy parameters used
** in conversion extended -> non-extended format
      double complex UASf(6,6,4)
      integer asfe1, asfe2


** Total widths and branching ratios
      double precision BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT
      COMMON /BRANCH/ BWNE,BWUD,BWTB,BZNN,BZEE,BZUU,BZDD,BZTT,
     &                BHWW,BHZZ,BHGG,BHTT,BHBB,BHCC,BHTAU,BHMU,
     &                BHGAM, BHGAMZ, XGW, XGZ, XGH, GAMT

** PDG index of higgs
      integer higgs_pdg


** Function which determines whether a parameter has already been set with a 
** reasonable value
      logical replace
      external replace

** Function that grabs the total width or branching ratio from the SLHA input
      double precision VSLHADecay
      external VSLHADecay

** Debugging flag
      integer debug

      
      call VSLHAClear(slhadata)


** Reading in SLHA data
      call VSLHAReadFile(error,slhadata,slhaFILE,0)

      write(*,*)'  '


** Checking whether we're working with CP violation
      if (DBLE(ModSel_CPV) .eq. 2D0) then
         write(*,*)'Effects of CP violation are included!'
         higgsmix = 3
      else
         write(*,*)'No CP violation effects are included!'
         higgsmix = 2
      end if


** ew gauge boson masses, gfermi:
      GF = SMInputs_GF
      xmz = SMInputs_MZ
      xmw = Mass_MW

** setting the em coupling
      invAlfaMZ = SMInputs_invAlfaMZ
      if (invAlfaMZ .ne. -999D0) then
         Alfa = 1/invAlfaMZ

         if (ewscheme .eq. 6) then
** note that, if we're working with ewscheme=6 we want to parameterise LO by
** alfa(0), but the slha file contains only alfa(MZ).  We therefore calculate
** Delta Alfa from the PDG value of Alfa(0) and the input value of Alfa(MZ)
            Alfa0 = 1/137.03599981D0
            DelAlfa = 1D0 - (Alfa0/Alfa)
            Alfa = Alfa0
         else
            DelAlfa = -999D0
         end if
      else
         Alfa = -999D0
      end if


** fermion masses:
      do type = 1, 4
         do gen1 = 1, 3
            Mf(type,gen1) = SMInputs_Mf(type,gen1)
         end do
      end do

      xmc = Mf(3,2)
      xmb = Mf(4,3)
      xmt = Mf(3,3)
      xmtau = Mf(2,3)


** trilinear coupling parameters:
      do type = 2, 4
         do gen1 = 1, 3
            do gen2 = 1, 3
               Af(type,gen1,gen2) = Af_Af(gen1,gen2,type)
            end do
         end do
      end do


** sfermion masses
      do sfe1 = 1, 2
         do type = 1, 4
            do gen1 = 1, 3
               MSf(sfe1,type,gen1) = Mass_MSf(sfe1,type,gen1)
            end do
         end do
      end do


** sfermion mixings:  (first in extended format, from slha, then convert)
      do sfe1 = 1, 6
         do sfe2 = 1, 6
            do type = 1, 4
               UASf(sfe1,sfe2,type) = ASfMix_UASf(sfe1,sfe2,type)
            end do
         end do
      end do

      do sfe1 = 1, 2
         do sfe2 = 1, 2
            do type = 1, 4
               do gen1 = 1,3
                  asfe1 = gen1 + 3*(sfe1-1)
                  asfe2 = gen1 + 3*(sfe2-1)
                  USf(sfe1,sfe2,type,gen1) = UASf(asfe1,asfe2,type)
               end do
            end do
         end do
      end do  

      do sfe1 = 1, 2
         do sfe2 = 1, 2
            do type = 1, 4
               if (replace(DBLE(USf(sfe1,sfe2,type,3)),0)) then
                  USf(sfe1,sfe2,type,3) = SfMix_USf(sfe1,sfe2,type)
               end if
            end do
         end do
      end do


** Higgs sector:
      MHiggs_loop(1) = Mass_Mh0
      MHiggs_loop(2) = Mass_MHH
      MHiggs_loop(3) = Mass_MA0
      if (replace(MHiggs_loop(3),0)) MHiggs_loop(3) = ExtPar_MA0
      if (replace(MHiggs_loop(3),0)) MHiggs_loop(3) = 
     -     sqrt(ExtPar_MA02)
      MHiggs_loop(4) = Mass_MHp
      if (replace(MHiggs_loop(4),0)) MHiggs_loop(4) = ExtPar_MHp

      XMH = MHiggs_loop(HiggsType)
      Mh0 = MHiggs_loop(1)
      MHH = MHiggs_loop(2)
      MA0 = MHiggs_loop(3)
      MHp = MHiggs_loop(4)
      
      TB = HMix_TB
      if (replace(TB,0)) TB = ExtPar_TB
      if (replace(TB,0)) TB = MinPar_TB
      
      MUE = ExtPar_MUE
      if (replace(DBLE(MUE),0)) MUE = HMix_MUE

      Aeff = Alpha_Alpha
      SAeff = SIN(Aeff)
      

** Chargino masses and mixings:
      do cha1 = 1, 2
         MCha(cha1) = Mass_MCha(cha1)
         do cha2 = 1, 2
            UCha(cha1,cha2) = UMix_UCha(cha1,cha2)
            VCha(cha1,cha2) = VMix_VCha(cha1,cha2)
         end do
      end do


** Neutralino masses and mixings:
      do neu1 = 1, 4
         MNeu(neu1) = Mass_MNeu(neu1)
         do neu2 = 1, 4
            ZNeu(neu1,neu2) = NMix_ZNeu(neu1,neu2)
         end do
      end do    


** Soft SUSY-breaking parameters in the sfermion sector
      M1SL = ExtPar_MSL(1)
      if (replace(M1SL,0)) M1SL = MSoft_MSL(1)
      M2SL = ExtPar_MSL(2)
      if (replace(M2SL,0)) M2SL = MSoft_MSL(2)
      M3SL = ExtPar_MSL(3)
      if (replace(M3SL,0)) M3SL = MSoft_MSL(3)

      M1SE = ExtPar_MSE(1)
      if (replace(M1SE,0)) M1SE = MSoft_MSE(1)
      M2SE = ExtPar_MSE(2)
      if (replace(M2SE,0)) M2SE = MSoft_MSE(2)
      M3SE = ExtPar_MSE(3)
      if (replace(M3SE,0)) M3SE = MSoft_MSE(3)

      M1SQ = ExtPar_MSQ(1)
      if (replace(M1SQ,0)) M1SQ = MSoft_MSQ(1)
      M2SQ = ExtPar_MSQ(2)
      if (replace(M2SQ,0)) M2SQ = MSoft_MSQ(2)
      M3SQ = ExtPar_MSQ(3)
      if (replace(M3SQ,0)) M3SQ = MSoft_MSQ(3)

      M1SU = ExtPar_MSU(1)
      if (replace(M1SU,0)) M1SU = MSoft_MSU(1)
      M2SU = ExtPar_MSU(2)
      if (replace(M2SU,0)) M2SU = MSoft_MSU(2)
      M3SU = ExtPar_MSU(3)
      if (replace(M3SU,0)) M3SU = MSoft_MSU(3)

      M1SD = ExtPar_MSD(1)
      if (replace(M1SD,0)) M1SD = MSoft_MSD(1)
      M2SD = ExtPar_MSD(2)
      if (replace(M2SD,0)) M2SD = MSoft_MSD(2)
      M3SD = ExtPar_MSD(3)
      if (replace(M3SD,0)) M3SD = MSoft_MSD(3)


** Gaugino mass parameters
      M_1 = ExtPar_M1
      if (replace(DBLE(M_1),0)) M_1 = MSoft_M1
      M_2 = ExtPar_M2
      if (replace(DBLE(M_2),0)) M_2 = MSoft_M2
      M_3 = ExtPar_M3
      if (replace(DBLE(M_3),0)) M_3 = MSoft_M3


** Debugging info:
      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'  '
         write(*,*)'Parameters as set by SLHA file'
         write(*,*)'NOTE that if the parameter is -999D0, we could not'
         write(*,*)'find it in the SLHA file, and will be overwritten'

         write(*,*)'  '
         write(*,*)'EW sector ..'
         write(*,*)'MW, MZ =', xmw, xmz
         write(*,*)'gfermi =', GF
         write(*,*)'Alfa(MZ) =', Alfa

         write(*,*)'  '
         write(*,*)'Fermion masses ..'
         do type = 1, 4
            do gen1 = 1, 3
               write(*,*)'Mf(',type,gen1,') =', Mf(type,gen1)
            end do
         end do

         write(*,*)'  '
         write(*,*)'Trilinear coupling parameters..'
         do type = 2, 4
            do gen1 = 1, 3
               do gen2 = 1, 3
                  write(*,*)'Af(', type,gen1,gen2,') = ',
     -                 Af(type,gen1,gen2)
               end do
            end do
         end do

         write(*,*)'  '
         write(*,*)'Sfermion masses ..'
         do gen1 = 1, 3
            write(*,*)'MSf(1,',gen1,') =', MSf(1,1,gen1)
         end do
         do type = 2, 4
            do gen1 = 1, 3
               write(*,*)'MSf(',type,gen1,') =', MSf(1,type,gen1), 
     -              MSf(2,type,gen1)
            end do
         end do

         write(*,*)'  '
         write(*,*)'Sfermion mixing ..'
         do type = 1, 4
            do gen1 = 1, 3
               do sfe1 = 1, 2
                  write(*,*)'USf(1/2,',sfe1,type,gen1,') =', 
     -                 USf(1,sfe1,type,gen1), USf(2,sfe1,type,gen1)
               end do
            end do
         end do

         write(*,*)'  '
         write(*,*)'Higgs sector ..'
         do type = 1, 4
            write(*,*)'MHiggs(', type,') = ', MHiggs_loop(type)
         end do
         write(*,*)'Tan beta =', TB
         write(*,*)'MUE =', MUE
         write(*,*)'SAeff =', SAeff
      
         write(*,*)'  '
         write(*,*)'Chargino masses and mixings ..'
         do cha1 = 1, 2
            write(*,*)'MCha(',cha1,') =', MCha(cha1)
            write(*,*)'UCha(',cha1,'1/2) =', UCha(cha1,1), UCha(cha1,2)
            write(*,*)'VCha(',cha1,'1/2) =', VCha(cha1,1), VCha(cha1,2)
         end do

         write(*,*)'  '
         write(*,*)'Neutralino masses and mixings ..'
         do neu1 = 1, 4
            write(*,*)'MNeu(', neu1,') =', MNeu(neu1)
            do neu2 = 1, 4
               write(*,*)'ZNeu(',neu1,neu2,') =', ZNeu(neu1,neu2)
            end do
         end do    

         write(*,*)'  '
         write(*,*)'Soft SUSY-breaking parameters ..'
         write(*,*)'M1SL =', M1SL
         write(*,*)'M2SL =', M2SL
         write(*,*)'M3SL =', M3SL
         write(*,*)'M1SE =', M1SE
         write(*,*)'M2SE =', M2SE
         write(*,*)'M3SE =', M3SE
         write(*,*)'M1SQ =', M1SQ
         write(*,*)'M2SQ =', M2SQ
         write(*,*)'M3SQ =', M3SQ
         write(*,*)'M1SU =', M1SU
         write(*,*)'M2SU =', M2SU
         write(*,*)'M3SU =', M3SU
         write(*,*)'M1SD =', M1SD
         write(*,*)'M2SD =', M2SD
         write(*,*)'M3SD =', M3SD

         write(*,*)'  '
         write(*,*)'Gaugino mass parameters ..'
         write(*,*)'M_1 =', M_1
         write(*,*)'M_2 =', M_2
         write(*,*)'M_3 =', M_3
      end if


** Now for the decay widths and branching ratios
** NOTE: VSLHADecay(slhadata, parentID, no. children, child ID1, child ID2)

** Z widths and decays
      XGZ = VSLHADecay(slhadata,23,0,0,0,0,0)
      BZNN = VSLHADecay(slhadata,23,2,12,-12,0,0) + ! decay to neutrinos
     -     VSLHADecay(slhadata,23,2,14,-14,0,0) + 
     -     VSLHADecay(slhadata,23,2,16,-16,0,0)
      if (BZNN .lt. 0D0) BZNN = -999D0
      BZEE = VSLHADecay(slhadata,23,2,11,-11,0,0) + ! decay to leptons
     -     VSLHADecay(slhadata,23,2,13,-13,0,0) + 
     -     VSLHADecay(slhadata,23,2,15,-15,0,0)
      if (BZEE .lt. 0D0) BZEE = -999D0
      BZUU = VSLHADecay(slhadata,23,2,2,-2,0,0) + ! decay to u-type quarks
     -     VSLHADecay(slhadata,23,2,4,-4,0,0) + 
     -     VSLHADecay(slhadata,23,2,6,-6,0,0)
      if (BZUU .lt. 0D0) BZUU = -999D0
      BZDD = VSLHADecay(slhadata,23,2,1,-1,0,0) + ! decay to d-type quarks
     -     VSLHADecay(slhadata,23,2,3,-3,0,0) + 
     -     VSLHADecay(slhadata,23,2,5,-5,0,0)
      if (BZDD .lt. 0D0) BZDD = -999D0
      BZTT = VSLHADecay(slhadata,23,2,6,-6,0,0)    ! decay to tops
      if (BZTT .lt. 0D0) BZTT = -999D0


** W widths and decays
      XGW = VSLHADecay(slhadata,24,0,0,0,0,0)
      BWNE = VSLHADecay(slhadata,24,2,11,-12,0,0) +  ! decay to lepton + 
                                                           !   neutrino
     -     VSLHADecay(slhadata,24,2,13,-14,0,0)+ 
     -     VSLHADecay(slhadata,24,2,15,-16,0,0)
      if (BWNE .lt. 0D0) BWNE = -999D0
      BWUD = VSLHADecay(slhadata,24,2,1,-2,0,0) ! decay to up and down
      BWTB = VSLHADecay(slhadata,24,2,5,-6,0,0) ! decay to top and bottom


** Top width
      GAMT = VSLHADecay(slhadata,6,0,0,0,0,0)


** Higgs decays:
      if (HiggsType .eq. 1) then
         higgs_pdg = 25
      else if (HiggsType .eq. 2) then
         higgs_pdg = 35
      else if (HiggsType .eq. 3) then
         higgs_pdg = 36
      end if
      XGH = VSLHADecay(slhadata,higgs_pdg,0,0,0,0,0)
      BHZZ = VSLHADecay(slhadata,higgs_pdg,2,23,23,0,0) ! decay to ZZ
      BHWW = VSLHADecay(slhadata,higgs_pdg,2,24,-24,0,0) ! decay to WW
      BHBB = VSLHADecay(slhadata,higgs_pdg,2,5,-5,0,0) ! decay to bottoms
      BHTT = VSLHADecay(slhadata,higgs_pdg,2,6,-6,0,0) ! decay to tops
      BHTAU = VSLHADecay(slhadata,higgs_pdg,2,15,-15,0,0) ! decay to taus
      BHCC = VSLHADecay(slhadata,higgs_pdg,2,4,-4,0,0) ! decay to charms
      BHMU = VSLHADecay(slhadata,higgs_pdg,2,13,-13,0,0) ! decay to muons
      BHGG = VSLHADecay(slhadata,higgs_pdg,2,21,21,0,0) ! decay to gluons
      BHGAM = VSLHADecay(slhadata,higgs_pdg,2,22,22,0,0) ! decay to ptns
      BHGAMZ = VSLHADecay(slhadata,higgs_pdg,2,22,23,0,0) ! decay to ptn+Z


      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'  '
         write(*,*)'Widths and branching ratios as set by SLHA file'
         write(*,*)'NOTE that if the width/ratio is -999D0, we could'
         write(*,*)'not find it in the SLHA file, and it will be'
         write(*,*)'overwritten later.'

         write(*,*)'  '
         write(*,*)'Z decays ..'
         write(*,*)'Z boson width =', XGZ
         write(*,*)'Z -> neutrinos BR =', BZNN
         write(*,*)'Z -> leptons BR =', BZEE
         write(*,*)'Z -> up-type quarks BR =', BZUU
         write(*,*)'Z -> down-type quarks BR =', BZDD
         write(*,*)'Z -> top quarks BR =', BZTT

         write(*,*)'  '
         write(*,*)'W decays ..'
         write(*,*)'W boson width =', XGW
         write(*,*)'W -> lepton+neutrino BR =', BWNE
         write(*,*)'Z -> up + down BR =', BWUD
         write(*,*)'Z -> top + bottom BR =', BWTB

         write(*,*)'  '
         write(*,*)'top width =', GAMT

         write(*,*)'  '
         write(*,*)'Higgs decays ..'
         write(*,*)'Higgs boson width =', XGH
         write(*,*)'H -> ZZ BR =', BHZZ
         write(*,*)'H -> WW BR =', BHWW
         write(*,*)'H -> top quarks BR =', BHTT
         write(*,*)'H -> bottom quarks BR =', BHBB
         write(*,*)'H -> charm quarks BR =', BHCC
         write(*,*)'H -> taus BR =', BHTAU
         write(*,*)'H -> muons BR =', BHMU
         write(*,*)'H -> gluons BR =', BHGG
         write(*,*)'H -> photons BR =', BHGAM
         write(*,*)'H -> photon + Z BR =', BHGAMZ    
         write(*,*)'  '     
      end if


      end


********************************************************************************
********************************************************************************

********************************************************************************
********************************************************************************
***                                                                          ***
*** This final selection of subroutines calculates the sfermion, chargino    ***
*** and neutralino sectors at tree level.  They are used if FeynHiggs is not ***
*** being used and the SLHA file from which we're determining the SUSY       ***
*** parameters does not contain them.  They are adapted from various parts   ***
*** of FormCalc (Thomas Hahn).                                               ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** Sfermion sector: Firstly, we consider the sfermion masses and mixings

      subroutine SfermionSector(MSfcalc,USfcalc,debug)

      implicit none

#include "mssm.inc"

** The calculated sfermion masses and mixings: this is the OUTPUT
      double precision MSfcalc(2,4,3)
      double complex USfcalc(2,2,4,3)

** Do-loop parameters: fermion type and generation
      integer type, gen

** Dummy parameters 
      double precision i3part, qpart
      double precision DSf(2,4) 
      double complex Xf(4,3,3)
      double precision M12sq, delta, t, c, m1, m2
      double complex s
      double precision msq(2)
      double precision CTB

** Components of the sfermion mass matrix
      double precision M11, M22
      double complex M12

** Electron charge
      double precision Qe
      parameter (Qe = -1D0)
** SUSY breaking masses, and their squares
      double precision MSS(2,2:4,3), MSS2(2,2:4,3)

** Debugging flag
      integer debug


** Setting SUSY breaking masses
      MSS(1,2,1) = M1SL
      MSS(1,2,2) = M2SL
      MSS(1,2,3) = M3SL

      MSS(2,2,1) = M1SE
      MSS(2,2,2) = M2SE
      MSS(2,2,3) = M3SE

      MSS(1,3,1) = M1SQ
      MSS(1,3,2) = M2SQ
      MSS(1,3,3) = M3SQ

      MSS(2,3,1) = M1SU
      MSS(2,3,2) = M2SU
      MSS(2,3,3) = M3SU

      MSS(2,4,1) = M1SD
      MSS(2,4,2) = M2SD
      MSS(2,4,3) = M3SD


** Some dummy parameters      
      i3part = 0.5D0*MZ2*C2B
      qpart = SW2*MZ2*C2B

      DSf(1,1) = i3part
      DSf(2,1) = 0
      DSf(2,2) = Qe*qpart
      DSf(1,2) = -i3part - DSf(2,2)
      DSf(2,3) = Qu*qpart
      DSf(1,3) = +i3part - DSf(2,3)
      DSf(2,4) = Qd*qpart
      DSf(1,4) = -i3part - DSf(2,4)


** Setting Xf:
      do type = 2, 4
         do gen = 1, 3
            if (type .eq. 3) then
               CTB = 1D0/TB
            else
               CTB = TB
            end if
            Xf(type,gen,gen) = Af(type,gen,gen) - DCONJG(MUE)*CTB
         end do
      end do


      do gen = 1, 3
         MSS(1,4,gen) = MSS(1,3,gen)
         
         MSS2(1,2,gen) = MSS(1,2,gen)**2
         MSS2(2,2,gen) = MSS(2,2,gen)**2
         MSS2(1,3,gen) = MSS(1,3,gen)**2
         MSS2(2,3,gen) = MSS(2,3,gen)**2
         MSS2(1,4,gen) = MSS(1,4,gen)**2
         MSS2(2,4,gen) = MSS(2,4,gen)**2

** Firstly: sneutrinos
         if((MSS2(1,2,gen) + DSf(1,1)) .gt. 1D300) then
            msq(1) = 0D0
            msq(2) = MSS2(1,2,gen) + DSf(1,1)
            USfcalc(1,1,1,gen) = 0
            USfcalc(2,2,1,gen) = 0
            USfcalc(1,2,1,gen) = 1
            USfcalc(2,1,1,gen) = 1
         else
            msq(1) = MSS2(1,2,gen) + DSf(1,1)
            msq(2) = 0D0
            USfcalc(1,1,1,gen) = 1
            USfcalc(2,2,1,gen) = 1
            USfcalc(1,2,1,gen) = 0
            USfcalc(2,1,1,gen) = 0
         endif

         if( msq(1) .lt. 0 ) then
            write(*,*) "Ooops! Negative mass squares in sfermion sector"
            stop
         endif

         MSfcalc(1,1,gen) = sqrt(msq(1))
         MSfcalc(2,1,gen) = sqrt(msq(2))

         
         do type = 2, 4

** Components of sfermion mass matrix:
            M11 = MSS2(1,type,gen) + DSf(1,type) + Mf2(type,gen)
            M22 = MSS2(2,type,gen) + DSf(2,type) + Mf2(type,gen)
            M12 = Mf(type,gen)*DCONJG(Xf(type,gen,gen))
            
** Calculating the sfermion masses and elements of the mixing matrix
            M12sq = DBLE(M12)**2 + DIMAG(M12)**2
            t = 0.5D0*(M11 - M22)
            t = 1/(t + sign(sqrt(t**2 + M12sq), t))
      
            delta = t*M12sq
            m1 = M11 + delta
            m2 = M22 - delta

            c = 1/sqrt(delta*t + 1)
            s = c*t*M12

** Setting sfermion masses and mixings such that m1 < m2
            if( m1 .gt. m2 ) then
               msq(1) = m2
               msq(2) = m1
               USfcalc(1,1,type,gen) = -DCONJG(s)
               USfcalc(2,2,type,gen) = s
               USfcalc(1,2,type,gen) = c
               USfcalc(2,1,type,gen) = c
            else
               msq(1) = m1
               msq(2) = m2
               USfcalc(1,1,type,gen) = c
               USfcalc(2,2,type,gen) = c
               USfcalc(1,2,type,gen) = s
               USfcalc(2,1,type,gen) = -DCONJG(s)
            endif

            MSfcalc(1,type,gen) = sqrt(msq(1))
            MSfcalc(2,type,gen) = sqrt(msq(2))
         enddo

      end do


      if (debug .ne. 0) then
         do type = 1, 4
            do gen = 1, 3
               write(*,*)' '
               write(*,*)'type, gen =', type, gen
               write(*,*)'MSfcalc =', MSfcalc(1,type,gen), 
     -              MSfcalc(2,type,gen)
               write(*,*)'USfcalc 11 =', USFcalc(1,1,type,gen)
               write(*,*)'USfcalc 12 =', USFcalc(1,2,type,gen)
               write(*,*)'USfcalc 21 =', USFcalc(2,1,type,gen)
               write(*,*)'USfcalc 22 =', USFcalc(2,2,type,gen)
            end do
         end do
      end if


      end


********************************************************************************
********************************************************************************

*** Now we focus on the chargino sector, calculating the masses and mixings

      subroutine CharginoSector(MChaCalc,UChaCalc,VChaCalc,debug)

      implicit none

** The calculated values of the chargino mass and mixing matrices
      double precision MChaCalc(2)
      double complex UChaCalc(2,2), VChaCalc(2,2)

** Components of chargino mass matrix
      double complex chaM(2,2), chaMinit(2,2)


** Dummy parameters used to diagonalise the matrices
      integer p, q, j, pi(3)
      double precision red, off, thresh
      double precision t, dv, dw, xv, xw, invc
      double complex x, y, ssv, ssw, tv, tw, evp, evq, f
      double complex A(3,3), V(2,3), W(3,3)   
      integer sweep

      double precision sq
      double complex c

** A function which checks whether the matrix neuM has changed: true if it has,
** false if it's stayed the same
      logical matCheck
      external matCheck

** Debugging flag
      integer debug

#include "mssm.inc"


** The function sq, used in the matrix
      sq(c) = DBLE(c*DCONJG(c)) !  a function


** Setting the components of the chargino mass matrix
      chaM(1,1) = M_2
      chaM(2,1) = sqrt(2D0)*MW*CB
      chaM(1,2) = sqrt(2D0)*MW*SB
      chaM(2,2) = MUE


** The following calculates the mass eigenvalues and the mixing matrices:
      do p = 1, 2
         do q = 1, 2
	    A(q,p) = 0
	    V(q,p) = 0
	    W(q,p) = 0
         enddo
         V(p,p) = 1
         W(p,p) = 1
      enddo

      do p = 1, 2
         do q = 1, 2
            A(q,p) = chaM(q,p)
         enddo
      enddo
      
      red = .01D0/2**4

      
      do p = 1, 2
         do q = 1, 2
            chaMinit(p,q) = chaM(p,q)
         end do
      end do


      do sweep = 1, 50

** This if loop runs on the first run through the 'sweep' do loop or if the
** matrix chaM has changed
** NOTE: see ntbk 19-1-11 for structure write-up
         if ((sweep .eq. 1) .or. (matCheck(chaM,chaMinit,2))) then
            off = 0
            do q = 1, 2
               do p = 1, q - 1
                  off = off + sq(A(p,q)) + sq(A(q,p))
               enddo
            enddo
         end if


         do p = 1, 2
            do q = 1, 2
               chaMinit(p,q) = chaM(p,q)
            end do
         end do
         

         if (off .gt. (2D0**(-102))) then

            thresh = 0
            if( sweep .lt. 4 ) thresh = off*red

            off = sq(A(1,2)) + sq(A(2,1))
            if( sweep .gt. 4 .and. off .lt.
     &           (2D0**(-102))*max(sq(A(1,1)), sq(A(2,2))) ) then
               A(1,2) = 0
               A(2,1) = 0
            else
               if( off .gt. thresh ) then
                  xv = DBLE((A(1,1) - A(2,2))*DCONJG(A(1,1) + A(2,2)))
                  xw = DBLE((A(1,2) - A(2,1))*DCONJG(A(1,2) + A(2,1)))
                  dv = .5D0*(xv + xw)
                  dw = .5D0*(xv - xw)
               
                  tv = DCONJG(A(1,1))*A(2,1) + A(2,2)*DCONJG(A(1,2))
                  xv = sqrt(dv**2 + sq(tv))
                  tw = DCONJG(A(1,1))*A(1,2) + A(2,2)*DCONJG(A(2,1))
                  xw = sqrt(dw**2 + sq(tw))
                  
                  t = sign(1D0, min(abs(dv + xv), abs(dw + xw)) -
     &                 min(abs(dv - xv), abs(dw - xw)))
                  
                  tv = tv/(dv + t*xv)
                  invc = sqrt(1 + sq(tv))
                  ssv = tv/invc
                  tv = tv/(invc + 1)
                  
                  tw = tw/(dw + t*xw)
                  invc = sqrt(1 + sq(tw))
                  ssw = tw/invc
                  tw = tw/(invc + 1)
                  
                  x = A(1,1)
                  y = A(2,1)
                  evp = invc*(x + DCONJG(ssv)*(y - tv*x))
                  x = A(1,2)
                  y = A(2,2)
                  evq = invc*(y - ssv*(x + DCONJG(tv)*y))
                  
                  do j = 1, 2
                     x = A(j,1)
                     y = A(j,2)
                     A(j,1) = x + DCONJG(ssw)*(y - tw*x)
                     A(j,2) = y - ssw*(x + DCONJG(tw)*y)
                     x = A(1,j)
                     y = A(2,j)
                     A(1,j) = x + DCONJG(ssv)*(y - tv*x)
                     A(2,j) = y - ssv*(x + DCONJG(tv)*y)
                  enddo
                  
                  A(1,1) = evp
                  A(2,1) = 0
                  A(1,2) = 0
                  A(2,2) = evq
                  
                  do j = 1, 2
                     x = V(j,1)
                     y = V(j,2)
                     V(j,1) = x + ssv*(y - DCONJG(tv)*x)
                     V(j,2) = y - DCONJG(ssv)*(x + tv*y)
                  enddo
                  
                  do j = 1, 2
                     x = W(j,1)
                     y = W(j,2)
                     W(j,1) = x + ssw*(y - DCONJG(tw)*x)
                     W(j,2) = y - DCONJG(ssw)*(x + tw*y)
                  enddo
               endif
            endif

            if (sweep .eq. 50) then
               write(*,*) "Bad convergence in chargino sector"
            end if
         end if  
      enddo


* make the diagonal elements nonnegative
      do p = 1, 2
         MChaCalc(p) = abs(A(p,p))
         if( (MChaCalc(p) .gt. (2D0**(-52))) .and. 
     -        (MChaCalc(p) .ne. DBLE(A(p,p))) ) then
	    f = A(p,p)/MChaCalc(p)
	    do q = 1, 2
               W(q,p) = W(q,p)*f
	    enddo
         endif
      enddo

* sort the singular values
      do p = 1, 2
         pi(p) = p
      enddo
      do p = 1, 2
         j = p
         t = MChaCalc(p)
         do q = p + 1, 2
            if( (t - MChaCalc(q)) .gt. 0 ) then
               j = q
               t = MChaCalc(q)
            endif
         enddo

         MChaCalc(j) = MChaCalc(p)
         MChaCalc(p) = t

         q = pi(j)
         pi(j) = pi(p)

         do j = 1, 2
	    UChaCalc(p,j) = V(j,q)
	    VChaCalc(p,j) = W(j,q)
         enddo
      enddo


      end


********************************************************************************
********************************************************************************

*** This subroutine looks at the neutralino sector, calculating the masses
*** and mixings

      subroutine NeutralinoSector(MNeuCalc, ZNeuCalc, debug)

      implicit none

** The calculated values of the neutralino mass and mixing matrix
      double precision MNeuCalc(4)
      double complex ZNeuCalc(4,4)

** The components of the neutralino mass matrix
      double complex neuM(4,4), neuMinit(4,4)

** A function which checks whether the matrix neuM has changed: true if it has,
** false if it's stayed the same
      logical matCheck
      external matCheck


** These parameters are used when finding the mass eigenvalues and components of
** the mixing matrix      
      integer p, q, j
      double precision red, off, thresh
      double precision sqp, sqq, t, invc
      double complex f, x, y
      double complex ev(2,4)
      
      integer sweep

      double precision sq
      double complex c

** Debugging flag
      integer debug


#include "mssm.inc"


      sq(c) = DBLE(c*DCONJG(c))  ! a function


* GUT relation (Higgs Hunters Guide p. 287)  !  SOPHY CHECK THIS ISN'T ELSEWHERE
      if (ABS(M_1) .le. 0D0) M_1 = 5/3D0*SW2/CW2*M_2


** Setting the components of the neutralino mass matrix
      neuM(1,1) = M_1
      neuM(2,2) = M_2
      neuM(3,1) = -MZ*SW*CB
      neuM(1,3) = neuM(3,1)
      neuM(4,1) = MZ*SW*SB
      neuM(1,4) = neuM(4,1)
      neuM(3,2) = MZ*CW*CB
      neuM(2,3) = neuM(3,2)
      neuM(4,2) = -MZ*CW*SB
      neuM(2,4) = neuM(4,2)
      neuM(4,3) = -MUE
      neuM(3,4) = neuM(4,3)
      neuM(2,1) = 0
      neuM(1,2) = 0
      neuM(3,3) = 0
      neuM(4,4) = 0


      do p = 1, 4
         do q = 1, 4
            neuMinit(p,q) = neuM(p,q)
         end do
      end do


** The following calculates the mass eigenvalues and the mixing matrices:
      do p = 1, 4
         ev(1,p) = 0
         ev(2,p) = neuM(p,p)
      enddo

      do p = 1, 4
         do q = 1, 4
	    ZNeuCalc(q,p) = 0
         enddo
         ZNeuCalc(p,p) = 1
      enddo

      red = .04D0/4**4


** Setting initial neuM
      do p = 1, 4
         do q = 1, 4
            neuMinit(p,q) = neuM(p,q)
         end do
      end do


      do sweep = 1, 50

** This first if loop runs on the first iteration of the sweep loop OR if the
** matrix neuM has been changed in the previous iteration of the sweep loop
** NOTE: see ntbk 19-1-11 for structure write-up
         if ((sweep .eq. 1) .or. matCheck(neuMinit,neuM,4)) then
            off = 0
            do q = 2, 4
               do p = 1, q - 1
                  off = off + sq(neuM(p,q))
               enddo
            enddo
         end if


         do p = 1, 4
            do q = 1, 4
               neuMinit(p,q) = neuM(p,q)
            end do
         end do      


         if(off .gt. (2D0**(-103))) then
         
            thresh = 0
            if( sweep .lt. 4 ) thresh = off*red
         
            do q = 2, 4
               do p = 1, q - 1

                  off = sq(neuM(p,q))
                  sqp = sq(ev(2,p))
                  sqq = sq(ev(2,q))

                  if( (sweep .gt. 4) .and.
     &                 (off .lt. (2D0**(-103))*max(sqp, sqq))) then
                     neuM(p,q) = 0
                  else
                     if( off .gt. thresh ) then

                        t = 0.5D0*abs(sqp - sqq)

                        if( t .eq. 0 ) then
                           f = 1
                           if( sqp .ne. 0 ) f = sqrt(ev(2,q)/ev(2,p))
                        else
                           f = sign(1D0, sqp - sqq)*
     &                          (ev(2,q)*DCONJG(neuM(p,q)) +
     &                          DCONJG(ev(2,p))*neuM(p,q))
                        endif

                        t = t + sqrt(t**2 + DBLE(f)**2 + DIMAG(f)**2)
                        f = f/t
                        
                        ev(1,p) = ev(1,p) + neuM(p,q)*DCONJG(f)
                        ev(2,p) = neuM(p,p) + ev(1,p)
                        ev(1,q) = ev(1,q) - neuM(p,q)*f
                        ev(2,q) = neuM(q,q) + ev(1,q)
                        
                        t = DBLE(f)**2 + DIMAG(f)**2
                        invc = sqrt(t + 1)
                        f = f/invc
                        t = t/(invc*(invc + 1))
                     
                        do j = 1, p - 1
                           x = neuM(j,p)
                           y = neuM(j,q)
                           neuM(j,p) = x + (DCONJG(f)*y - t*x)
                           neuM(j,q) = y - (f*x + t*y)
                        enddo
                     
                        do j = p + 1, q - 1
                           x = neuM(p,j)
                           y = neuM(j,q)
                           neuM(p,j) = x + (DCONJG(f)*y - t*x)
                           neuM(j,q) = y - (f*x + t*y)
                        enddo
                        
                        do j = q + 1, 4
                           x = neuM(p,j)
                           y = neuM(q,j)
                           neuM(p,j) = x + (DCONJG(f)*y - t*x)
                           neuM(q,j) = y - (f*x + t*y)
                        enddo

                        neuM(p,q) = 0
                        
                        do j = 1, 4
                           x = ZNeuCalc(p,j)
                           y = ZNeuCalc(q,j)
                           ZNeuCalc(p,j) = x + (f*y - t*x)
                           ZNeuCalc(q,j) = y - (DCONJG(f)*x + t*y)
                        enddo
                     endif
                  endif
               enddo
            enddo

            do p = 1, 4
               ev(1,p) = 0
               neuM(p,p) = ev(2,p)
            enddo
            
            if (sweep .eq. 50) then
               write(*,*) "Bad convergence in the neutralino sector"
            end if

         end if  

      enddo


* make the diagonal elements nonnegative
      do p = 1, 4
         MNeuCalc(p) = abs(neuM(p,p))

         if((MNeuCalc(p) .gt. (2D0**(-52))) .and. 
     -        (MNeuCalc(p) .ne. DBLE(neuM(p,p))) ) then
	    f = sqrt(neuM(p,p)/MNeuCalc(p))
	    do q = 1, 4
               ZNeuCalc(p,q) = ZNeuCalc(p,q)*f
	    enddo
         endif

      enddo
              
* sort the eigenvalues
      do p = 1, 3
         j = p
         t = MNeuCalc(p)
         do q = p + 1, 4
	    if( (t - MNeuCalc(q)) .gt. 0 ) then
               j = q
               t = MNeuCalc(q)
	    endif
         enddo
         
         if( j .ne. p ) then
	    MNeuCalc(j) = MNeuCalc(p)
	    MNeuCalc(p) = t
	    do q = 1, 4
               x = ZNeuCalc(p,q)
               ZNeuCalc(p,q) = ZNeuCalc(j,q)
               ZNeuCalc(j,q) = x
	    enddo
         endif
      enddo


      end


********************************************************************************
********************************************************************************

*** If any of the elements of the matrices mat1 and mat2 are different, this
*** function returns a 'true' value, otherwise it's 'false'

      logical function matCheck(mat1,mat2,n)

      implicit none

** Size of matrices and do loop parameters
      integer n, i, j

** The two matrices we're checking and the absolute difference between elements
      double complex mat1(n,n), mat2(n,n)
      double precision matDif


      matCheck = .false.


      do i = 1, n
         do j = 1, n
            matDif = ABS(mat1(i,j)-mat2(i,j))
            if (matDif .gt. 0.001D0) then
               matCheck = .true.
            end if
         end do
      end do


      end


********************************************************************************

********************************************************************************
********************************************************************************
***                                                                          ***
*** This set of subroutines does the technical parts of reading in from a    ***
*** SLHA file.  They are taken, with a few small modifications, from         ***
*** SLHALib-2.2, by Thomas Hahn.                                             ***
*** The SLHA format is described in hep-ph/0311123 and 0801.0045 [hep-ph].   ***
*** An example SLHA file (sps1a.slha) is provided with the input files.      ***
***                                                                          ***
********************************************************************************
********************************************************************************

** The main subroutine reads the SLHA file, and files all the information in 
** the array sd

      subroutine VSLHAReadFile(error, slhain, filename, abort)

      implicit none


#include "PDG.h"
#include "SLHAIntern.h"

      integer error, abort      
      double complex slhain(nslhadata)
      character*(*) filename

      double complex slhadata(nslhadata)
      double precision sd(2*nslhadata)
      equivalence (slhadata,sd)

      integer VSLHATokens
      external VSLHATokens

      integer unit, block, im, n, n1, n2, n3, pdgno
      integer decay, parent, parent_id, i, j, p, id(10)
      double precision blockQ, val
      character*80 text
      
      integer lineno
      character*250 line
      common /slhaline/ lineno, line

      integer m
      parameter (m = 2**30 - 1)
      
#define Range(i,min,max) \
iand(max-i,m).le.max-min

#define Vector(x,n1,max) \
if( Range(n1,1,max) ) x(n1)

#define Matrix(x,n1,n2,max) \
if( Range(n1,1,max).and.Range(n2,1,max) ) x(n1,n2)

#define RMatrix(x,n1,n2,max1,max2) \
if( Range(n1,1,max1).and.Range(n2,1,max2) ) x(n1,n2)

#define Triple(x,n1,n2,n3,max) \
if( Range(n1,1,max).and.Range(n2,1,max).and.Range(n3,1,max) ) x(n1,n2,n3)

#define Scan \
read(line, *, err=997, end=997)

#include "SLHAReadBlocks.h"

      slhadata = slhain

      error = 0
      block = 0
      pdgno = 0
      decay = 1


*      if( filename(1:lnblnk(filename)) .eq. "-" ) then
      if( trim(filename) .eq. "-" ) then
         unit = 5
         if( lineno .ne. 0 ) goto 501
      else
         open(10, file=filename, status="old", err=996)
         unit = 10
         lineno = 0
      endif

 100  format(A)
 101  format(A, I5)
 102  format(A, A, A, I5)
 103  format(A, A)
      
 500  continue
      lineno = lineno + 1
      read(unit, 100, end=999) line
 501  if( VSLHATokens(line) .eq. 0 ) goto 500
      
      if( line(1:6) .eq. "BLOCK " ) then
         i = 7
         im = 1
         if( line(7:8) .eq. "IM" ) then
	    im = 2
	    i = 9
         endif
         n = index(line(i:), " ")
         
         do block = 1, nblocks
	    if( line(i:i+n-1) .eq. blockname(block)(1:n) ) then
               blockQ = 0
               if( line(i+n:i+n+1) .eq. "Q=" )
     &              read(line(i+n+2:), *) blockQ
               goto 500
	    endif
         enddo
         print 102, "Unknown block ", line(7:n+3+2*i),
     &        " in input line", lineno
         block = -1
         goto 500
      endif
      
      if( line(1:6) .eq. "DECAY " ) then
         if( decay + 2 .ge. LengthDecay ) then
	    print 101, "Out of memory in input line", lineno
	    block = -1
         else
	    block = tecys
	    read(line(7:), *, err=997, end=997) parent_id, val
	    parent = decay
	    decay = decay + 2
	    Decay(parent) = Decay_Entry(parent_id, 2)
	    Decay(parent + 1) = val
         endif
         goto 500
      endif
      
      goto (
     &     modse, sminp, minpa, extpa, qxtpa, nmrun, mass, tmass,
     &     nmix, umix, vmix, staum, stopm, sbotm,
     &     alfa, talfa, hmix, gauge, msoft,
     &     ae, au, ad, ye, yu, yd,
     &     lllei, llqdi, luddi, llle, llqd, ludd,
     &     tllei, tlqdi, tuddi, tlle, tlqd, tudd,
     &     rki, rk, rdi, rd, rvevi, rvev, rmlhi, rmlh,
     &     rnmix, rumix, rvmix, rhmix, ramix, rlmix,
     &     vckmi, vckm, umnsi, umns,
     &     msl2i, mse2i, msq2i, msu2i, msd2i,
     &     msl2, mse2, msq2, msu2, msd2, 
     &     tei, tui, tdi, te, tu, td,
     &     snmix, slmix, usmix, tsmix, ssmix, samix,
     &     hcmix, nnmix, nhmix, namix,
*     &     prcob, tecys
     &     prcob, spinf, tcinf, tecys
     &     ) block
      
      if( block .eq. 0 ) then
* depending on the abort flag, one can stop reading here and
* assume the file is in a different format...
         if( abort .ne. 0 ) goto 998
* ... or issue a warning and continue
         print 101, "Superfluous text in input line", lineno
      endif

      goto 500

 996  print 103, "Cannot open ", filename
      error = 1
      return

 997  print 101, "Syntax error in input line", lineno
      goto 500

modse	Scan n, val
      if( n .eq. 1  ) ModSel_Model = val
      if( n .eq. 3  ) ModSel_Content = val
      if( n .eq. 4  ) ModSel_RPV = val
      if( n .eq. 5  ) ModSel_CPV = val
      if( n .eq. 6  ) ModSel_FV = val
      if( n .eq. 11 ) ModSel_GridPts = val
      if( n .eq. 12 ) ModSel_Qmax = val
      if( n .eq. 21 ) then
         pdgno = pdgno + 1
         if( pdgno .le. 5 ) ModSel_PDG(pdgno) = val
      endif
      goto 500

sminp	Scan n, val
      if( n .eq. 1  ) SMInputs_invAlfaMZ = val
      if( n .eq. 2  ) SMInputs_GF = val
      if( n .eq. 3  ) SMInputs_AlfasMZ = val
      if( n .eq. 4  ) SMInputs_MZ = val
      if( n .eq. 5  ) SMInputs_Mb = val
      if( n .eq. 6  ) SMInputs_Mt = val
      if( n .eq. 7  ) SMInputs_Mtau = val
      if( n .eq. 8  ) SMInputs_Mnu3 = val
      if( n .eq. 11 ) SMInputs_Me = val
      if( n .eq. 12 ) SMInputs_Mnu1 = val
      if( n .eq. 13 ) SMInputs_Mmu = val
      if( n .eq. 14 ) SMInputs_Mnu2 = val
      if( n .eq. 21 ) SMInputs_Md = val
      if( n .eq. 22 ) SMInputs_Mu = val
      if( n .eq. 23 ) SMInputs_Ms = val
      if( n .eq. 24 ) SMInputs_Mc = val
      goto 500

minpa	Scan n, val
      if( n .eq. 1 ) MinPar_M0 = val
      if( n .eq. 2 ) MinPar_M12 = val
      if( n .eq. 3 ) MinPar_TB = val
      if( n .eq. 4 ) MinPar_signMUE = val
      if( n .eq. 5 ) MinPar_A = val
      if( n .eq. 6 ) MinPar_cgrav = val
      goto 500
      
extpa	Scan n, val
      if( n .eq. 0  ) ExtPar_Q = val
      if( n .eq. 1  ) ExtPar_M1 = val
      if( n .eq. 2  ) ExtPar_M2 = val
      if( n .eq. 3  ) ExtPar_M3 = val
      if( n .eq. 11 ) ExtPar_At = val
      if( n .eq. 12 ) ExtPar_Ab = val
      if( n .eq. 13 ) ExtPar_Atau = val
      if( n .eq. 21 ) ExtPar_MHd2 = val
      if( n .eq. 22 ) ExtPar_MHu2 = val
      if( n .eq. 23 ) ExtPar_MUE = val
      if( n .eq. 24 ) ExtPar_MA02 = val
      if( n .eq. 25 ) ExtPar_TB = val
      if( n .eq. 26 ) ExtPar_MA0 = val
      if( n .eq. 27 ) ExtPar_MHp = val
      if( Range(n, 31, 33) ) ExtPar_MSL(n - 30) = val
      if( Range(n, 34, 36) ) ExtPar_MSE(n - 33) = val
      if( Range(n, 41, 43) ) ExtPar_MSQ(n - 40) = val
      if( Range(n, 44, 46) ) ExtPar_MSU(n - 43) = val
      if( Range(n, 47, 49) ) ExtPar_MSD(n - 46) = val
      if( Range(n, 51, 53) ) ExtPar_N5(n - 50) = val
      if( n .eq. 61 ) ExtPar_lambda = val
      if( n .eq. 62 ) ExtPar_kappa = val
      if( n .eq. 63 ) ExtPar_Alambda = val
      if( n .eq. 64 ) ExtPar_Akappa = val
      if( n .eq. 65 ) ExtPar_lambdaS = val
      if( n .eq. 66 ) ExtPar_xiF = val
      if( n .eq. 67 ) ExtPar_xiS = val
      if( n .eq. 68 ) ExtPar_MUEprime = val
      if( n .eq. 69 ) ExtPar_mS2prime = val
      if( n .eq. 70 ) ExtPar_mS2 = val
      goto 500

qxtpa	Scan n, val
      if( n .eq. 1  ) QExtPar_QM1 = val
      if( n .eq. 2  ) QExtPar_QM2 = val
      if( n .eq. 3  ) QExtPar_QM3 = val
      if( n .eq. 11 ) QExtPar_QAt = val
      if( n .eq. 12 ) QExtPar_QAb = val
      if( n .eq. 13 ) QExtPar_QAtau = val
      if( n .eq. 21 ) QExtPar_QMHd2 = val
      if( n .eq. 22 ) QExtPar_QMHu2 = val
      if( n .eq. 23 ) QExtPar_QMUE = val
      if( n .eq. 24 ) QExtPar_QMA02 = val
      if( n .eq. 25 ) QExtPar_QTB = val
      if( n .eq. 31 ) QExtPar_QMSL = val
      if( n .eq. 34 ) QExtPar_QMSE = val
      if( n .eq. 41 ) QExtPar_QMSQ = val
      if( n .eq. 44 ) QExtPar_QMSU = val
      if( n .eq. 47 ) QExtPar_QMSD = val
      goto 500

nmrun	Scan n, val
      NMSSMRun_Q = blockQ
      if( n .eq. 1 ) NMSSMRun_lambda = val
      if( n .eq. 2 ) NMSSMRun_kappa = val
      if( n .eq. 3 ) NMSSMRun_Alambda = val
      if( n .eq. 4 ) NMSSMRun_Akappa = val
      if( n .eq. 5 ) NMSSMRun_lambdaS = val
      if( n .eq. 6 ) NMSSMRun_xiF = val
      if( n .eq. 7 ) NMSSMRun_xiS = val
      if( n .eq. 8 ) NMSSMRun_MUEprime = val
      if( n .eq. 9 ) NMSSMRun_mS2prime = val
      if( n .eq. 10 ) NMSSMRun_mS2 = val
      goto 500

mass	Scan n, val
      if( n .eq. PDG_nu_e        ) Mass_Mf(1,1) = val
      if( n .eq. PDG_electron    ) Mass_Mf(2,1) = val
      if( n .eq. PDG_up          ) Mass_Mf(3,1) = val
      if( n .eq. PDG_down        ) Mass_Mf(4,1) = val
      if( n .eq. PDG_nu_mu       ) Mass_Mf(1,2) = val
      if( n .eq. PDG_muon        ) Mass_Mf(2,2) = val
      if( n .eq. PDG_charm       ) Mass_Mf(3,2) = val
      if( n .eq. PDG_strange     ) Mass_Mf(4,2) = val
      if( n .eq. PDG_nu_tau      ) Mass_Mf(1,3) = val
      if( n .eq. PDG_tau         ) Mass_Mf(2,3) = val
      if( n .eq. PDG_top         ) Mass_Mf(3,3) = val
      if( n .eq. PDG_bottom      ) Mass_Mf(4,3) = val
      if( n .eq. PDG_snu_e1      ) Mass_MSf(1,1,1) = val
      if( n .eq. PDG_snu_e2      ) Mass_MSf(2,1,1) = val
      if( n .eq. PDG_selectron1  ) Mass_MSf(1,2,1) = val
      if( n .eq. PDG_selectron2  ) Mass_MSf(2,2,1) = val
      if( n .eq. PDG_sup1        ) Mass_MSf(1,3,1) = val
      if( n .eq. PDG_sup2        ) Mass_MSf(2,3,1) = val
      if( n .eq. PDG_sdown1      ) Mass_MSf(1,4,1) = val
      if( n .eq. PDG_sdown2      ) Mass_MSf(2,4,1) = val
      if( n .eq. PDG_snu_mu1     ) Mass_MSf(1,1,2) = val
      if( n .eq. PDG_snu_mu2     ) Mass_MSf(2,1,2) = val
      if( n .eq. PDG_smuon1      ) Mass_MSf(1,2,2) = val
      if( n .eq. PDG_smuon2      ) Mass_MSf(2,2,2) = val
      if( n .eq. PDG_scharm1     ) Mass_MSf(1,3,2) = val
      if( n .eq. PDG_scharm2     ) Mass_MSf(2,3,2) = val
      if( n .eq. PDG_sstrange1   ) Mass_MSf(1,4,2) = val
      if( n .eq. PDG_sstrange2   ) Mass_MSf(2,4,2) = val
      if( n .eq. PDG_snu_tau1    ) Mass_MSf(1,1,3) = val
      if( n .eq. PDG_snu_tau2    ) Mass_MSf(2,1,3) = val
      if( n .eq. PDG_stau1       ) Mass_MSf(1,2,3) = val
      if( n .eq. PDG_stau2       ) Mass_MSf(2,2,3) = val
      if( n .eq. PDG_stop1       ) Mass_MSf(1,3,3) = val
      if( n .eq. PDG_stop2       ) Mass_MSf(2,3,3) = val
      if( n .eq. PDG_sbottom1    ) Mass_MSf(1,4,3) = val
      if( n .eq. PDG_sbottom2    ) Mass_MSf(2,4,3) = val
      if( n .eq. PDG_Z           ) Mass_MZ = val
      if( n .eq. PDG_W           ) Mass_MW = val
      if( n .eq. PDG_h0          ) Mass_Mh0 = val
      if( n .eq. PDG_HH          ) Mass_MHH = val
      if( n .eq. PDG_A0          ) Mass_MA0 = val
      if( n .eq. PDG_Hp          ) Mass_MHp = val
      if( n .eq. PDG_neutralino1 ) Mass_MNeu(1) = val
      if( n .eq. PDG_neutralino2 ) Mass_MNeu(2) = val
      if( n .eq. PDG_neutralino3 ) Mass_MNeu(3) = val
      if( n .eq. PDG_neutralino4 ) Mass_MNeu(4) = val
      if( n .eq. PDG_neutralino5 ) Mass_MNeu(5) = val
      if( n .eq. PDG_chargino1   ) Mass_MCha(1) = val
      if( n .eq. PDG_chargino2   ) Mass_MCha(2) = val
      if( n .eq. PDG_gluino      ) Mass_MGl = val
      if( n .eq. PDG_gravitino   ) Mass_MGrav = val
      goto 500
      
tmass	Scan n, val
      if( n .eq. PDG_h0 ) DMass_DeltaMh0 = val
      if( n .eq. PDG_HH ) DMass_DeltaMHH = val
      if( n .eq. PDG_A0 ) DMass_DeltaMA0 = val
      if( n .eq. PDG_Hp ) DMass_DeltaMHp = val
      goto 500

nmix	Scan n1, n2, val
      Matrix(NMix_ZNeu, n1,n2, 4) = val
      goto 500

umix	Scan n1, n2, val
      Matrix(UMix_UCha, n1,n2, 2) = val
      goto 500

vmix	Scan n1, n2, val
      Matrix(VMix_VCha, n1,n2, 2) = val
      goto 500

staum	Scan n1, n2, val
      Matrix(StauMix_USf, n1,n2, 2) = val
      goto 500

stopm	Scan n1, n2, val
      Matrix(StopMix_USf, n1,n2, 2) = val
      goto 500

sbotm	Scan n1, n2, val
      Matrix(SbotMix_USf, n1,n2, 2) = val
      goto 500

alfa	Scan Alpha_Alpha
      goto 500
        
talfa	Scan DAlpha_DeltaAlpha
      goto 500

hmix	Scan n, val
      HMix_Q = blockQ
      if( n .eq. 1 ) HMix_MUE = val
      if( n .eq. 2 ) HMix_TB = val
      if( n .eq. 3 ) HMix_VEV = val
      if( n .eq. 4 ) HMix_MA02 = val
      goto 500

gauge	Scan n, val
      Gauge_Q = blockQ
      if( n .eq. 1 ) Gauge_g1 = val
      if( n .eq. 2 ) Gauge_g2 = val
      if( n .eq. 3 ) Gauge_g3 = val
      goto 500

msoft	Scan n, val
      MSoft_Q = blockQ
      if( n .eq. 1  ) MSoft_M1 = val
      if( n .eq. 2  ) MSoft_M2 = val
      if( n .eq. 3  ) MSoft_M3 = val
      if( n .eq. 21 ) MSoft_MHd2 = val
      if( n .eq. 22 ) MSoft_MHu2 = val
      if( Range(n, 31, 33) ) MSoft_MSL(n - 30) = val
      if( Range(n, 34, 36) ) MSoft_MSE(n - 33) = val
      if( Range(n, 41, 43) ) MSoft_MSQ(n - 40) = val
      if( Range(n, 44, 46) ) MSoft_MSU(n - 43) = val
      if( Range(n, 47, 49) ) MSoft_MSD(n - 46) = val
      goto 500

ae	Scan n1, n2, val
      Ae_Q = blockQ
      Matrix(Ae_Af, n1,n2, 3) = val
      goto 500

au	Scan n1, n2, val
      Au_Q = blockQ
      Matrix(Au_Af, n1,n2, 3) = val
      goto 500

ad	Scan n1, n2, val
      Ad_Q = blockQ
      Matrix(Ad_Af, n1,n2, 3) = val
      goto 500

ye	Scan n1, n2, val
      Ye_Q = blockQ
      Matrix(Ye_Yf, n1,n2, 3) = val
      goto 500

yu	Scan n1, n2, val
      Yu_Q = blockQ
      Matrix(Yu_Yf, n1,n2, 3) = val
      goto 500
      
yd	Scan n1, n2, val
      Yd_Q = blockQ
      Matrix(Yd_Yf, n1,n2, 3) = val
      goto 500
      
lllei	Scan n1, n2, n3, val
      Triple(RVLamLLEIn_lamLLE, n1,n2,n3, 3) = val
      goto 500
      
llqdi	Scan n1, n2, n3, val
      Triple(RVLamLQDIn_lamLQD, n1,n2,n3, 3) = val
      goto 500

luddi	Scan n1, n2, n3, val
      Triple(RVLamUDDIn_lamUDD, n1,n2,n3, 3) = val
      goto 500

llle	Scan n1, n2, n3, val
      RVLamLLE_Q = blockQ
      Triple(RVLamLLE_lamLLE, n1,n2,n3, 3) = val
      goto 500

llqd	Scan n1, n2, n3, val
      RVLamLQD_Q = blockQ
      Triple(RVLamLQD_lamLQD, n1,n2,n3, 3) = val
      goto 500

ludd	Scan n1, n2, n3, val
      RVLamUDD_Q = blockQ
      Triple(RVLamUDD_lamUDD, n1,n2,n3, 3) = val
      goto 500

tllei	Scan n1, n2, n3, val
      Triple(RVTLLEIn_TLLE, n1,n2,n3, 3) = val
      goto 500

tlqdi	Scan n1, n2, n3, val
      Triple(RVTLQDIn_TLQD, n1,n2,n3, 3) = val
      goto 500

tuddi	Scan n1, n2, n3, val
      Triple(RVTUDDIn_TUDD, n1,n2,n3, 3) = val
      goto 500
      
tlle	Scan n1, n2, n3, val
      RVTLLE_Q = blockQ
      Triple(RVTLLE_TLLE, n1,n2,n3, 3) = val
      goto 500

tlqd	Scan n1, n2, n3, val
      RVTLQD_Q = blockQ
      Triple(RVTLQD_TLQD, n1,n2,n3, 3) = val
      goto 500
      
tudd	Scan n1, n2, n3, val
      RVTUDD_Q = blockQ
      Triple(RVTUDD_TUDD, n1,n2,n3, 3) = val
      goto 500

rki	Scan n, val
      Vector(RVKappaIn_kappa, n, 3) = val
      goto 500

rk	Scan n, val
      RVKappa_Q = blockQ
      Vector(RVKappa_kappa, n, 3) = val
      goto 500

rdi	Scan n, val
      Vector(RVDIn_D, n, 3) = val
      goto 500

rd	Scan n, val
      RVD_Q = blockQ
      Vector(RVD_D, n, 3) = val
      goto 500
      
rvevi	Scan n, val
      Vector(RVSnVEVIn_VEV, n, 3) = val
      goto 500

rvev	Scan n, val
      RVSnVEV_Q = blockQ
      Vector(RVSnVEV_VEV, n, 3) = val
      goto 500

rmlhi	Scan n, val
      Vector(RVM2LH1In_M2LH1, n, 3) = val
      goto 500

rmlh	Scan n, val
      RVM2LH1_Q = blockQ
      Vector(RVM2LH1_M2LH1, n, 3) = val
      goto 500

rnmix	Scan n1, n2, val
      Matrix(RVNMix_ZNeu, n1,n2, 7) = val
      goto 500

rumix	Scan n1, n2, val
      Matrix(RVUMix_UCha, n1,n2, 5) = val
      goto 500

rvmix	Scan n1, n2, val
      Matrix(RVVMix_VCha, n1,n2, 5) = val
      goto 500

rhmix	Scan n1, n2, val
      Matrix(RVHMix_UH, n1,n2, 5) = val
      goto 500

ramix	Scan n1, n2, val
      RMatrix(RVAMix_UA, n1,n2, 4,5) = val
      goto 500

rlmix	Scan n1, n2, val
      RMatrix(RVLMix_CLep, n1,n2, 7,8) = val
      goto 500

vckmi	Scan n, val
      if( n .eq. 1 ) VCKMIn_lambda = val
      if( n .eq. 2 ) VCKMIn_A = val
      if( n .eq. 3 ) VCKMIn_rhobar = val
      if( n .eq. 4 ) VCKMIn_etabar = val
      goto 500

vckm	Scan n1, n2, val
      VCKM_Q = blockQ
      Matrix(VCKM_VCKM, n1,n2, 3) = val
      goto 500
      
umnsi	Scan n, val
      if( n .eq. 1 ) UPMNSIn_theta12 = val
      if( n .eq. 2 ) UPMNSIn_theta23 = val
      if( n .eq. 3 ) UPMNSIn_theta13 = val
      if( n .eq. 4 ) UPMNSIn_delta13 = val
      if( n .eq. 5 ) UPMNSIn_alpha1 = val
      if( n .eq. 6 ) UPMNSIn_alpha2 = val
      goto 500

umns	Scan n1, n2, val
      UPMNS_Q = blockQ
      Matrix(UPMNS_UPMNS, n1,n2, 3) = val
      goto 500
      
msl2i	Scan n1, n2, val
      Matrix(MSL2In_MSL2, n1,n2, 3) = val
      goto 500

mse2i	Scan n1, n2, val
      Matrix(MSE2In_MSE2, n1,n2, 3) = val
      goto 500

msq2i	Scan n1, n2, val
      Matrix(MSQ2In_MSQ2, n1,n2, 3) = val
      goto 500

msu2i	Scan n1, n2, val
      Matrix(MSU2In_MSU2, n1,n2, 3) = val
      goto 500

msd2i	Scan n1, n2, val
      Matrix(MSD2In_MSD2, n1,n2, 3) = val
      goto 500

msl2	Scan n1, n2, val
      MSL2_Q = blockQ
      Matrix(MSL2_MSL2, n1,n2, 3) = val
      goto 500

mse2	Scan n1, n2, val
      MSE2_Q = blockQ
      Matrix(MSE2_MSE2, n1,n2, 3) = val
      goto 500
      
msq2	Scan n1, n2, val
      MSQ2_Q = blockQ
      Matrix(MSQ2_MSQ2, n1,n2, 3) = val
      goto 500

msu2	Scan n1, n2, val
      MSU2_Q = blockQ
      Matrix(MSU2_MSU2, n1,n2, 3) = val
      goto 500

msd2	Scan n1, n2, val
      MSD2_Q = blockQ
      Matrix(MSD2_MSD2, n1,n2, 3) = val
      goto 500

tei	Scan n1, n2, val
      Matrix(TeIn_Tf, n1,n2, 3) = val
      goto 500

tui	Scan n1, n2, val
      Matrix(TuIn_Tf, n1,n2, 3) = val
      goto 500

tdi	Scan n1, n2, val
      Matrix(TdIn_Tf, n1,n2, 3) = val
      goto 500

te	Scan n1, n2, val
      Te_Q = blockQ
      Matrix(Te_Tf, n1,n2, 3) = val
      goto 500

tu	Scan n1, n2, val
      Tu_Q = blockQ
      Matrix(Tu_Tf, n1,n2, 3) = val
      goto 500

td	Scan n1, n2, val
      Td_Q = blockQ
      Matrix(Td_Tf, n1,n2, 3) = val
      goto 500

snmix	Scan n1, n2, val
      Matrix(SnuMix_UASf, n1,n2, 3) = val
      goto 500

slmix	Scan n1, n2, val
      Matrix(SelMix_UASf, n1,n2, 6) = val
      goto 500

usmix	Scan n1, n2, val
      Matrix(USqMix_UASf, n1,n2, 6) = val
      goto 500

tsmix	Scan n1, n2, val
      Matrix(DSqMix_UASf, n1,n2, 6) = val
      goto 500

ssmix	Scan n1, n2, val
      Matrix(SnSMix_US, n1,n2, 3) = val
      goto 500

samix	Scan n1, n2, val
      Matrix(SnAMix_UA, n1,n2, 3) = val
      goto 500

hcmix	Scan n1, n2, val
      Matrix(CVHMix_UH, n1,n2, 4) = val
      goto 500

nnmix	Scan n1, n2, val
      Matrix(NMNMix_ZNeu, n1,n2, 5) = val
      goto 500

nhmix	Scan n1, n2, val
      Matrix(NMHMix_UH, n1,n2, 3) = val
      goto 500

namix	Scan n1, n2, val
      Matrix(NMAMix_UA, n1,n2, 3) = val
      goto 500

prcob	Scan n, val
      if( n .eq. 1  ) PrecObs_DeltaRho = val
      if( n .eq. 2  ) PrecObs_MWMSSM = val
      if( n .eq. 3  ) PrecObs_MWSM = val
      if( n .eq. 4  ) PrecObs_SW2effMSSM = val
      if( n .eq. 5  ) PrecObs_SW2effSM = val
      if( n .eq. 11 ) PrecObs_gminus2mu = val
      if( n .eq. 21 ) PrecObs_EDMeTh = val
      if( n .eq. 22 ) PrecObs_EDMn = val
      if( n .eq. 23 ) PrecObs_EDMHg = val
      if( n .eq. 31 ) PrecObs_bsgammaMSSM = val
      if( n .eq. 32 ) PrecObs_bsgammaSM = val
      if( n .eq. 33 ) PrecObs_DeltaMsMSSM = val
      if( n .eq. 34 ) PrecObs_DeltaMsSM = val
      goto 500

spinf	Scan n, text
      goto 500

tcinf	Scan n, text
      goto 500

tecys	Scan val, n, (id(i), i = 1, n)
      if( val .eq. 0 ) goto 500
      if( decay + n + 1 .ge. LengthDecay ) then
         print 101, "Out of memory in input line", lineno
         block = -1
         goto 500
      endif

      do i = 1, n
         p = i
         do j = i + 1, n
	    if( id(j) .lt. id(p) ) p = j
         enddo
         Decay(decay) = Decay_Entry(id(p), n + 1)
         decay = decay + 1
         id(p) = id(i)
      enddo
      Decay(decay) = val
      decay = decay + 1
      Decay(parent) = Decay_Entry(parent_id, decay - parent)
      goto 500

 998  error = 2
      slhain = slhadata
      if( unit .eq. 5 ) return

 999  close(unit)

      slhain = slhadata

      end


********************************************************************************
********************************************************************************

** This function returns the selected branching ratio
**     parent_id = PDG id of decaying particle
**     nchildren = number of particles parent decays into
**     childX_id = PDG id of decay product

      double precision function VSLHADecay(slhain, parent_id,
     &     nchildren, child1_id, child2_id, child3_id, child4_id)

      implicit none

#include "SLHAIntern.h"

      double complex slhain(nslhadata)
      integer parent_id
      integer nchildren, child1_id, child2_id, child3_id, child4_id

      integer decay, id(4), i, j, p, nextparent, nextchild
      
      double complex slhadata(nslhadata)
      double precision sd(2*nslhadata)
      equivalence (slhadata,sd)

      slhadata = slhain

** The BR is set to -999 if it's not present in the SLHA file     
      VSLHADecay = invalid


*      VSLHADecay = 0
      
      decay = 1


      do while( Decay_Id(decay) .ne. parent_id )
         if( Decay(decay) .eq. invalid ) return
         decay = decay + Decay_Next(decay)
      enddo

      if( nchildren .eq. 0 ) goto 1
      if( nchildren .ge. 1 ) id(1) = child1_id
      if( nchildren .ge. 2 ) id(2) = child2_id
      if( nchildren .ge. 3 ) id(3) = child3_id
      if( nchildren .ge. 4 ) id(4) = child4_id
      
      do i = 1, nchildren
         p = i
         do j = i + 1, nchildren
	    if( id(j) .lt. id(p) ) p = j
         enddo
         if( p .ne. i ) then
	    j = id(i)
	    id(i) = id(p)
	    id(p) = j
         endif
      enddo

      nextparent = decay + Decay_Next(decay)
      decay = decay + 2

      do while( decay .lt. nextparent )
         nextchild = decay + Decay_Next(decay)
         if( decay + nchildren + 1 .eq. nextchild ) then
	    if( Decay_Id(decay) .eq. id(1) ) then
               if( decay + 2 .eq. nextchild ) goto 1
               decay = decay + 1
               if( Decay_Id(decay) .eq. id(2) ) then
                  if( decay + 2 .eq. nextchild ) goto 1
                  decay = decay + 1
                  if( Decay_Id(decay) .eq. id(3) ) then
                     if( decay + 2 .eq. nextchild ) goto 1
                     decay = decay + 1
                     if( Decay_Id(decay) .eq. id(4) ) then
                        if( decay + 2 .eq. nextchild ) goto 1
                     endif
                  endif
               endif
	    endif
         endif
         decay = nextchild
      enddo

      slhain = slhadata

      return

 1    VSLHADecay = Decay(decay + 1)

      slhain = slhadata

      end


********************************************************************************
********************************************************************************

      integer function VSLHATokens(s)

      implicit none

      character*(*) s
      integer si, di, c

      VSLHATokens = 0
      di = 0
      do si = 1, len(s)
         c = ichar(s(si:si))
         if( c .eq. ichar('#') ) goto 1
         if( c .eq. ichar(' ') .or. c .eq. 9 ) then
	    if( di .ne. 0 .and.
     &           s(di:di) .ne. ' ' .and. s(di:di) .ne. '=' ) then
               VSLHATokens = VSLHATokens + 1
               di = di + 1
               s(di:di) = ' '
	    endif
         else
	    if( c .ne. ichar('=') .or.
     &           di .eq. 0 .or. s(di:di) .ne. ' ' )
     &           di = di + 1
	    if( c .gt. 96 ) c = iand(c, 223)
	    s(di:di) = char(c)
         endif
      enddo
      
 1    do si = di + 1, len(s)
         s(si:si) = ' '
      enddo


      end


********************************************************************************
********************************************************************************

** This subroutine clears the values of the slhadata array

      subroutine VSLHAClear(slhadata)

      implicit none

#include "SLHADefs.h"

      double complex slhadata(*)

      integer i


      do i = 1, OffsetDecays
         slhadata(i) = invalid
      enddo

      do i = OffsetDecays + 1, nslhadata
         slhadata(i) = DCMPLX(invalid, invalid)
      enddo


      end


********************************************************************************
