********************************************************************************
********************************************************************************
*** renormalisation.F                                                        ***
*** sophy@particle.uni-karlsruhe.de, 26-11-2010                              ***
***                                                                          ***
*** This file is all to do with renormalisation.  It sets the                ***
*** renormalisation scales with LoopTools, and calculates the                ***
*** renormalisation constants for whichever model we're working on.  It can  ***
*** also calculate the Higgs propagator Z matrix.  This file is only         ***
*** compiled if LoopTools has been linked.                                   ***
********************************************************************************
********************************************************************************

*** Subroutine to set the renormalisation scales using LoopTools

      subroutine renSCALES(scalefactor,debug)

      implicit none

      double precision scalefactor

** Debugging flag
      integer debug

#include "mssm.inc"
#include "global.inc"
#include "looptools.h"


** Setting the renormalisation scales.  Please note that mudim = mu^2 and 
** setlambda sets the photon mass SQUARED.
      call setdelta(0D0)
      call setmudim(MT2)
      call setlambda(1D0)
      massQ = 1d-5

      Divergence = getdelta()
      massG = getlambda()
      scalefactor = sqrt(getmudim()/MT2)

** Outputting renormalisation scales for reference purposes
      if (ewcor_switch) then
         write(*,*)'  '
         write(*,*)'Renormalisation scales ...'
         write(*,*)'delta =', getdelta()
         write(*,*)'mudim =', getmudim()
         write(*,*)'lambda =', getlambda()
         write(*,*)'massQ =', massQ
      end if


      end


********************************************************************************
********************************************************************************

*** This subroutine calls the appropriate (i.e. sm/mssm) subroutine to 
*** calculate the renormalisation constants

      subroutine ren_con(charge)

      implicit none

** Flag for charge renormalisation
      integer charge

#include "global.inc"
#include "mssm.inc"


      if (model .eq. 1) then             ! SM
         call ren_conSM(charge,1)
      else if (model .eq. 2) then        ! MSSM
** Calculating renconFULL for use in the hVV vertex when ewapprox = 4 is used
         if ((ewapprox .eq. 4) .or. 
     -        ((propLevel .ne. 0) .and. (.not. feynH_switch))) then
            call ren_conSUSY(charge,1,.true.)
         end if
         call ren_conSUSY(charge,1,.false.)
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the renormalisation constants in the standard 
*** model.  See notebook 5, page 28 for definitions, but be wary of 
*** sin-theta-w conventions.

      subroutine ren_conSM(charge,debug)

      implicit none


** em coupling
      double precision Alfa, Alfa2

** Particle do loop parameters
      integer g3ini, Gen3

** Flag for charge renormalisation
      integer charge

** Dummy variables used to calculate dZet1 (otherwise known as the bane of my
** life).  The charge renormalisation is dependent on the parameter 'charge',
** which is set in parameters.F according to the vbfnlo.dat flag 'ewscheme'
      double complex dZAAlt1, dZeR
      double complex dZAAlt1fm, dZeRfm

** Dummy parameters
      double complex dup2, dup3, dup4, dup5, dup6

** Debugging flag
      integer debug

** This is needed for looptools. 
#include "looptools.h"

** Declaration of parameters
#include "mssm.inc"
#include "global.inc"


** Clearing the values of the renormalization constants
      dMWsq1 = (0D0,0D0)
      dMZsq1 = (0D0,0D0)
      dSW1 = (0D0,0D0)
      dZet1 = (0D0,0D0)
      dZW1 = (0D0,0D0)
      dZZZ1 = (0D0,0D0)
      dZH1 = (0D0,0D0)
      dZh0 = (0D0,0D0)
      dZHH = (0D0,0D0)
      dZh0HH = (0D0,0D0)
      dMAZsq1 = (0D0,0D0)
      dTB1 = (0D0,0D0)
      dZAA1 = (0D0,0D0)
      dZZA1 = (0D0,0D0)
      dZAZ1 = (0D0,0D0)
      dZfLu = (0D0,0D0)
      dZfLd = (0D0,0D0)
      dZfLc = (0D0,0D0)
      dZfLs = (0D0,0D0)
      dZfRu = (0D0,0D0)
      dZfRd = (0D0,0D0)
      dZfRc = (0D0,0D0)
      dZfRs = (0D0,0D0)


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Dummy parameters: clearing their values
      dZAAlt1 = (0D0,0D0)
      dZeR = (0D0,0D0)


** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in koppln.F and is 1 for generation 3 quarks, 
** 2 for all fermions and 3 for all particles
      if (sector .eq. 1) then
         g3ini = 3
      else 
         g3ini = 1
      end if


** FERMIONIC CONTRIBUTIONS
**    For the quarks 
      do Gen3 = g3ini, 3

         dMWsq1 = dMWsq1 - (3*Alfa)/(4.D0*Pi*SW2)*
     -        DBLE(A0(Mf2(4,Gen3))) + (3*Alfa)/(2.D0*Pi*SW2)*
     -        DBLE(B0i(bb00,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa*MW2)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(4.D0*Pi*SW2)*
     -        (DBLE(B0i(bb0,MW2,Mf2(3,Gen3),Mf2(4,Gen3)))*Mf2(3,Gen3))

         dMZsq1 = dMZsq1 - (Alfa*(9 + SW2*(-24 + 32*SW2)))/
     -        (24.D0*CW2*Pi*SW2)*DBLE(A0(Mf2(3,Gen3))) - 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(A0(Mf2(4,Gen3))) - 
     -        (Alfa*(-9 + (24 - 32*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -        (Alfa*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -        (Alfa*MZ2*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -        (Alfa*MZ2*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))*Mf2(3,Gen3))-
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3))


         if (charge .ne. 2) then ! i.e. LO via alfa(0) or alfa(MZ)

            if (Gen3 .eq. 3) then
               dZAA1 = dZAA1 + (4*Alfa)/(3.D0*Pi)*
     -              DBLE(B0i(bb1,0.D0,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -              (8*Alfa)/(3.D0*Pi)*DBLE(B0i(dbb00,0.D0,Mf2(3,Gen3),
     -              Mf2(3,Gen3))) 
            end if

         else  ! LO via GF

            dZeR = dZeR - Alfa/(8.D0*MW2*Pi*SW2)*(-3*A0(Mf2(4,Gen3)) + 
     -           6*B0i(bb00,0.D0,Mf2(3,Gen3),Mf2(4,Gen3)) - 
     -           3*B0i(bb0,0.D0,Mf2(3,Gen3),Mf2(4,Gen3))*Mf2(3,Gen3))
         end if


         dZAZ1 = dZAZ1 - 2/MZ2*(-((Alfa*(-3 + 8*SW2))/(6.D0*CW*Pi)*
     -        DBLE(A0(Mf2(3,Gen3)))/SW) - (Alfa*(-3 + 4*SW2))/
     -        (12.D0*CW*Pi)*DBLE(A0(Mf2(4,Gen3)))/SW - 
     -        (Alfa*(3 - 8*SW2))/(3.D0*CW*Pi)*
     -        DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -        (Alfa*(3 - 4*SW2))/(6.D0*CW*Pi)*
     -        DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))/SW - 
     -        (Alfa*MZ2*(-3 + 8*SW2))/(6.D0*CW*Pi)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -        (Alfa*MZ2*(-3 + 4*SW2))/(12.D0*CW*Pi)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))/SW)

         dZZA1 = dZZA1 + 2/MZ2*(-((Alfa*(-3 + 8*SW2))/
     -        (6.D0*CW*Pi)*DBLE(A0(Mf2(3,Gen3)))/SW) - 
     -        (Alfa*(-3 + 4*SW2))/(12.D0*CW*Pi)*
     -        DBLE(A0(Mf2(4,Gen3)))/SW - 
     -        (Alfa*(3 - 8*SW2))/(3.D0*CW*Pi)*
     -        DBLE(B0i(bb00,0.D0,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -        (Alfa*(3 - 4*SW2))/(6.D0*CW*Pi)*
     -        DBLE(B0i(bb00,0.D0,Mf2(4,Gen3),Mf2(4,Gen3)))/SW)

         dZH1 = dZH1 + (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -        (DBLE(B0i(bb1,MH2,Mf2(3,Gen3),Mf2(3,Gen3)))*Mf2(3,Gen3))+
     -        (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -        (MH2*DBLE(B0i(dbb1,MH2,Mf2(3,Gen3),Mf2(3,Gen3)))*
     -        Mf2(3,Gen3)) + (3*Alfa)/(2.D0*MW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MH2,Mf2(3,Gen3),Mf2(3,Gen3)))*
     -        Mf2(3,Gen3)**2) + (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -        (DBLE(B0i(bb1,MH2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3))+
     -        (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -        (MH2*DBLE(B0i(dbb1,MH2,Mf2(4,Gen3),Mf2(4,Gen3)))*
     -        Mf2(4,Gen3)) + (3*Alfa)/(2.D0*MW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MH2,Mf2(4,Gen3),Mf2(4,Gen3)))*
     -        Mf2(4,Gen3)**2)

         dZW1 = dZW1 + (3*Alfa)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(2.D0*Pi*SW2)*
     -        DBLE(B0i(dbb00,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa*MW2)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(dbb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa)/(4.D0*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MW2,Mf2(3,Gen3),Mf2(4,Gen3)))*Mf2(3,Gen3))

         dZZZ1 = dZZZ1 + (Alfa*(9 + SW2*(-24 + 32*SW2)))/
     -        (24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (Alfa*(-9 + (24 - 32*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (Alfa*MZ2*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*MZ2*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))*Mf2(3,Gen3))
     -        + (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3))

      end do  


** charge renormalisation: parametrising by alfa(0) or alfa(MZ)
      if (charge .ne. 2) then

         do Gen3 = 1, 3

** light up-type quarks:
            if (Gen3 .ne. 3) then
               dZAAlt1 = dZAAlt1 - 1/MZ2*(-((4*Alfa)/(3.D0*Pi)*
     -              DBLE(A0(Mf2(3,Gen3)))) + (8*Alfa)/(3.D0*Pi)*
     -              DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -              (4*Alfa*MZ2)/(3.D0*Pi)*
     -              DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))))
            end if

** down-type quarks                
            dZAAlt1 = dZAAlt1 - 1/MZ2*(-Alfa/(3.D0*Pi)*
     -           DBLE(A0(Mf2(4,Gen3))) + (2*Alfa)/(3.D0*Pi)*
     -           DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -           (Alfa*MZ2)/(3.D0*Pi)*
     -           DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))))

** leptons
            dZAAlt1 = dZAAlt1 - 1/MZ2*(-(Alfa/Pi*
     -           DBLE(A0(Mf2(2,Gen3)))) + (2*Alfa)/Pi*
     -           DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     -           (Alfa*MZ2)/Pi*DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),
     -           Mf2(2,Gen3))))

** If we're working with only third generation quarks, the remaining
** part of Delta Alfa must be calculated.  See ntbk, 30-7-10
            if (sector .eq. 1) then

               if (Gen3 .ne. 3) then
** light up-type quarks
                  dZAA1 = dZAA1 - (4*Alfa)/(3.D0*Pi)*DBLE(B0i(bb1,0.D0,
     -                 Mf2(3,Gen3),Mf2(3,Gen3))) + (8*Alfa)/(3.D0*Pi)*
     -                 DBLE(B0i(dbb00,0.D0,Mf2(3,Gen3),Mf2(3,Gen3))) 
                  
** down-type quarks
                  dZAA1 = dZAA1 - Alfa/(3.D0*Pi)*DBLE(B0i(bb1,0.D0,
     -                 Mf2(4,Gen3),Mf2(4,Gen3))) + (2*Alfa)/(3.D0*Pi)*
     -                 DBLE(B0i(dbb00,0.D0,Mf2(4,Gen3),Mf2(4,Gen3)))
               end if

** leptons
               dZAA1 = dZAA1 - Alfa/Pi*DBLE(B0i(bb1,0.D0,Mf2(2,Gen3),
     -              Mf2(2,Gen3))) + (2*Alfa)/Pi*DBLE(B0i(dbb00,0.D0,
     -              Mf2(2,Gen3),Mf2(2,Gen3)))

            end if

         end do

      end if


**    For the leptons
      if (sector .ne. 1) then 

         dMZsq1 = dMZsq1 + Alfa/(8.D0*CW2*Pi*SW2)*
     -        (6*DBLE(B0i(bb00,MZ2,0.D0,0.D0)) - 
     -        3*MZ2*DBLE(B0i(bb1,MZ2,0.D0,0.D0)))
         
         dZZZ1 = dZZZ1 + Alfa/(8.D0*CW2*Pi*SW2)*
     -        (3*DBLE(B0i(bb1,MZ2,0.D0,0.D0)) - 
     -        6*DBLE(B0i(dbb00,MZ2,0.D0,0.D0)) + 
     -        3*MZ2*DBLE(B0i(dbb1,MZ2,0.D0,0.D0)))

         do Gen3 = 1, 3

            dMWsq1 = dMWsq1 - Alfa/(4.D0*Pi*SW2)*DBLE(A0(Mf2(2,Gen3))) + 
     -           Alfa/(2.D0*Pi*SW2)*
     -           DBLE(B0i(bb00,MW2,0.D0,Mf2(2,Gen3))) - (Alfa*MW2)/
     -           (4.D0*Pi*SW2)*DBLE(B0i(bb1,MW2,0.D0,Mf2(2,Gen3)))

            dMZsq1 = dMZsq1 - (Alfa*(1 + SW2*(-4 + 8*SW2)))/
     -           (8.D0*CW2*Pi*SW2)*DBLE(A0(Mf2(2,Gen3))) - 
     -           (Alfa*(-1 + (4 - 8*SW2)*SW2))/(4.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     -           (Alfa*MZ2*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))*Mf2(2,Gen3))


            if (charge .eq. 2) then  !  LO via GF
               dZeR = dZeR - Alfa/(8.D0*MW2*Pi*SW2)*(-A0(Mf2(2,Gen3)) + 
     -              2*B0i(bb00,0.D0,0.D0,Mf2(2,Gen3)))
            end if


            dZAZ1 = dZAZ1 - 2/MZ2*(-((Alfa*(-1 + 4*SW2))/(4.D0*CW*Pi)*
     -           DBLE(A0(Mf2(2,Gen3)))/SW) - 
     -           (Alfa*(1 - 4*SW2))/(2.D0*CW*Pi)*
     -           DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))/SW - 
     -           (Alfa*MZ2*(-1 + 4*SW2))/(4.D0*CW*Pi)*
     -           DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))/SW)

            dZZA1 = dZZA1 + 2/MZ2*(-((Alfa*(-1 + 4*SW2))/
     -           (4.D0*CW*Pi)*DBLE(A0(Mf2(2,Gen3)))/SW) - 
     -           (Alfa*(1 - 4*SW2))/(2.D0*CW*Pi)*
     -           DBLE(B0i(bb00,0.D0,Mf2(2,Gen3),Mf2(2,Gen3)))/SW)


            dZH1 = dZH1 + Alfa/(4.D0*MW2*Pi*SW2)*
     -           (DBLE(B0i(bb1,MH2,Mf2(2,Gen3),Mf2(2,Gen3)))*
     -           Mf2(2,Gen3)) + Alfa/(4.D0*MW2*Pi*SW2)*
     -           (MH2*DBLE(B0i(dbb1,MH2,Mf2(2,Gen3),Mf2(2,Gen3)))*
     -           Mf2(2,Gen3)) + Alfa/(2.D0*MW2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MH2,Mf2(2,Gen3),Mf2(2,Gen3)))*
     -           Mf2(2,Gen3)**2)

            dZW1 = dZW1 + Alfa/(4.D0*Pi*SW2)*
     -           DBLE(B0i(bb1,MW2,0.D0,Mf2(2,Gen3)))-Alfa/(2.D0*Pi*SW2)* 
     -           DBLE(B0i(dbb00,MW2,0.D0,Mf2(2,Gen3))) + (Alfa*MW2)/
     -           (4.D0*Pi*SW2)*DBLE(B0i(dbb1,MW2,0.D0,Mf2(2,Gen3)))

            dZZZ1 = dZZZ1 + (Alfa*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2*Pi*
     -           SW2)*DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           (Alfa*(-1 + (4 - 8*SW2)*SW2))/(4.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(dbb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           (Alfa*MZ2*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(dbb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))*Mf2(2,Gen3))
         end do
      end if


** fermion values
      dMWsq1fm = dMWsq1
      dMZsq1fm = dMZsq1
      dZW1fm = dZW1
      dZZZ1fm = dZZZ1
      dZH1fm = dZH1
      dMAZsq1fm = dMAZsq1
      dZAA1fm = dZAA1
      dZZA1fm = dZZA1
      dZAZ1fm = dZAZ1
      dZAAlt1fm = dZAAlt1
      dZeRfm = dZeR
      dZh0fm = dZh0
      dZHHfm = dZHH
      dZh0HHfm = dZh0HH
      dTB1fm = dTB1


** Bosonic contributions
      if (sector .eq. 3) then        
         dMWsq1 = dMWsq1 - (Alfa*MW2)/(6.D0*Pi*SW2) - 
     -        (Alfa*MW2*(1 - CW2 - SW2))/(2.D0*Pi*SW2) + 
     -        Alfa/(16.D0*Pi*SW2)*DBLE(A0(MH2)) - 
     -        (Alfa*(-7 + 4*SW2))/(8.D0*Pi*SW2)*DBLE(A0(MW2)) + 
     -        (Alfa*(1/4.D0 + CW2))/(4.D0*Pi*SW2)*DBLE(A0(MZ2)) - 
     -        (Alfa*MW2)/Pi*DBLE(B0i(bb0,MW2,0.D0,MW2)) + 
     -        (Alfa*MW2)/(4.D0*Pi*SW2)*DBLE(B0i(bb0,MW2,MH2,MW2)) + 
     -        ((-5*Alfa*CW2*MW2)/(4.D0*Pi*SW2) - 
     -        (Alfa*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(bb0,MW2,MW2,MZ2)) - 
     -        (2*Alfa)/Pi*DBLE(B0i(bb00,MW2,0.D0,MW2)) - 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb00,MW2,MH2,MW2)) - 
     -        (Alfa*(1/4.D0 + 2*CW2))/(Pi*SW2)*
     -        DBLE(B0i(bb00,MW2,MW2,MZ2)) - 
     -        (Alfa*MW2)/(2.D0*Pi)*DBLE(B0i(bb1,MW2,0.D0,MW2)) - 
     -        (Alfa*CW2*MW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MW2,MW2,MZ2)) 

         dMZsq1 = dMZsq1 - (Alfa*CW2*MZ2)/(6.D0*Pi*SW2) - 
     -        (Alfa*(-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/(8.D0*CW2*Pi*SW2)*
     -        DBLE(A0(MW2)) + Alfa/(16.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(MH2)) + DBLE(A0(MZ2))) + 
     -        (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*
     -        DBLE(B0i(bb0,MZ2,MH2,MZ2)) + 
     -        ((-5*Alfa*CW2*MZ2)/(4.D0*Pi*SW2) - 
     -        (Alfa*MW2*(CW2**2 - SW2**2))/(2.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(bb0,MZ2,MW2,MW2)) - 
     -        Alfa/(4.D0*CW2*Pi*SW2)*DBLE(B0i(bb00,MZ2,MH2,MZ2)) - 
     -        (Alfa*(CW2*(9*CW2 - 2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,MW2,MW2)) - 
     -        (Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MZ2,MW2,MW2))


         if (charge .ne. 2) then  ! LO via alfa(0), alfa(MZ)
            dZAA1 = dZAA1 + Alfa*(1/(6.D0*Pi) + 
     -           5/(4.D0*Pi)*DBLE(B0i(bb0,0.D0,MW2,MW2)) + 
     -           1/(2.D0*Pi)*DBLE(B0i(bb1,0.D0,MW2,MW2)) + 
     -           3/Pi*DBLE(B0i(dbb00,0.D0,MW2,MW2)))
         else  !  LO via GF
            dZeR = dZeR - 1/(2.D0*MW2)*(-(Alfa*MW2*(1 - CW2 - SW2))/
     -           (2.D0*Pi*SW2) + Alfa/(16.D0*Pi*SW2)*A0(MH2) - 
     -           (Alfa*(-7 + 4*SW2))/(8.D0*Pi*SW2)*A0(MW2) + 
     -           (Alfa*(1/4.D0 + CW2))/(4.D0*Pi*SW2)*A0(MZ2) + 
     -           (Alfa*MW2)/(4.D0*Pi)*B0i(bb0,0.D0,0.D0,MW2) + 
     -           (Alfa*MW2)/(4.D0*Pi*SW2)*B0i(bb0,0.D0,MH2,MW2) - 
     -           (Alfa*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2)*
     -           B0i(bb0,0.D0,MW2,MZ2) - 
     -           (2*Alfa)/Pi*B0i(bb00,0.D0,0.D0,MW2) - 
     -           Alfa/(4.D0*Pi*SW2)*B0i(bb00,0.D0,MH2,MW2) - (Alfa*
     -           (1/4.D0 + 2*CW2))/(Pi*SW2)*B0i(bb00,0.D0,MW2,MZ2))
         end if


         dZAZ1 = dZAZ1 - (2/MZ2*((Alfa*CW*MZ2)/(6.D0*Pi)/SW + 
     -        DBLE(B0i(bb0,MZ2,MW2,MW2))*((Alfa*MW2)/(2.D0*
     -        CW*Pi)/SW + (5*Alfa*CW*MZ2)/(4.D0*Pi)/SW) - 
     -        (Alfa*(5*CW2 - SW2))/(4.D0*CW*Pi)*DBLE(A0(MW2))/SW - 
     -        (Alfa*(-5*CW2 + SW2))/(2.D0*CW*Pi)*
     -        DBLE(B0i(bb00,MZ2,MW2,MW2))/SW + 
     -        (Alfa*CW*MZ2)/(2.D0*Pi)*DBLE(B0i(bb1,MZ2,MW2,MW2))/SW))

         dZZA1 = dZZA1 + 2/MZ2*(-((Alfa*(5*CW2 - SW2))/(4.D0*CW*Pi)*
     -        DBLE(A0(MW2))/SW) + 
     -        (Alfa*MW2)/(2.D0*CW*Pi)*DBLE(B0i(bb0,0.D0,MW2,MW2))/SW - 
     -        (Alfa*(-5*CW2 + SW2))/(2.D0*CW*Pi)*
     -        DBLE(B0i(bb00,0.D0,MW2,MW2))/SW)

         dZH1 = dZH1 - (DBLE(B0i(dbb0,MH2,MW2,MW2))*
     -        (-(Alfa/(8.D0*Pi*SW2)*MH2) + 
     -        Alfa/(16.D0*MW2*Pi*SW2)*(MH2**2 + 12*MW2**2))) + 
     -        Alfa/(8.D0*Pi*SW2)*DBLE(B0i(bb0,MH2,MW2,MW2)) + 
     -        Alfa/(16.D0*CW2*Pi*SW2)*DBLE(B0i(bb0,MH2,MZ2,MZ2)) - 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb1,MH2,MW2,MW2)) - 
     -        Alfa/(8.D0*CW2*Pi*SW2)*DBLE(B0i(bb1,MH2,MZ2,MZ2)) - 
     -        (9*Alfa)/(32.D0*MW2*Pi*SW2)*
     -        (MH2**2*DBLE(B0i(dbb0,MH2,MH2,MH2))) - 
     -        Alfa/(4.D0*Pi*SW2)*(MH2*DBLE(B0i(dbb1,MH2,MW2,MW2))) - 
     -        Alfa/(8.D0*CW2*Pi*SW2)*(MH2*DBLE(B0i(dbb1,MH2,MZ2,MZ2))) - 
     -        DBLE(B0i(dbb0,MH2,MZ2,MZ2))*
     -        (-(Alfa/(16.D0*CW2*Pi*SW2)*MH2) + 
     -        Alfa/(32.D0*CW2**2*MW2*Pi*SW2)*
     -        (16*MW2**2 + CW2*(-4*MW2*MZ2 + CW2*MH2**2)))

         dZW1 = dZW1 + Alfa/(6.D0*Pi*SW2) + 
     -        (5*Alfa)/(4.D0*Pi)*DBLE(B0i(bb0,MW2,0.D0,MW2)) + 
     -        (5*Alfa*CW2)/(4.D0*Pi*SW2)*DBLE(B0i(bb0,MW2,MW2,MZ2)) + 
     -        Alfa/(2.D0*Pi)*DBLE(B0i(bb1,MW2,0.D0,MW2)) + 
     -        (Alfa*CW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MW2,MW2,MZ2)) + 
     -        (Alfa*MW2)/Pi*DBLE(B0i(dbb0,MW2,0.D0,MW2)) - 
     -        (Alfa*MW2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb0,MW2,MH2,MW2)) + 
     -        ((5*Alfa*CW2*MW2)/(4.D0*Pi*SW2) + 
     -        (Alfa*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(dbb0,MW2,MW2,MZ2)) + 
     -        (2*Alfa)/Pi*DBLE(B0i(dbb00,MW2,0.D0,MW2)) + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(dbb00,MW2,MH2,MW2)) + 
     -        (Alfa*(1/4.D0 + 2*CW2))/(Pi*SW2)*
     -        DBLE(B0i(dbb00,MW2,MW2,MZ2)) + 
     -        (Alfa*MW2)/(2.D0*Pi)*DBLE(B0i(dbb1,MW2,0.D0,MW2)) + 
     -        (Alfa*CW2*MW2)/(2.D0*Pi*SW2)*DBLE(B0i(dbb1,MW2,MW2,MZ2))

         dZZZ1 = dZZZ1 + (Alfa*CW2)/(6.D0*Pi*SW2) + 
     -        (5*Alfa*CW2)/(4.D0*Pi*SW2)*DBLE(B0i(bb0,MZ2,MW2,MW2)) + 
     -        (Alfa*CW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MZ2,MW2,MW2)) - 
     -        (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*
     -        DBLE(B0i(dbb0,MZ2,MH2,MZ2)) + 
     -        ((5*Alfa*CW2*MZ2)/(4.D0*Pi*SW2) + 
     -        (Alfa*MW2*(CW2**2 - SW2**2))/(2.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(dbb0,MZ2,MW2,MW2)) + 
     -        Alfa/(4.D0*CW2*Pi*SW2)*DBLE(B0i(dbb00,MZ2,MH2,MZ2)) + 
     -        (Alfa*(CW2*(9*CW2 - 2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,MW2,MW2)) + 
     -        (Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*DBLE(B0i(dbb1,MZ2,MW2,MW2))


** fermion field renormalisation
         dZfLu = Alfa/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MU2,MD2,MW2)) + 
     -        (Alfa*MU2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MU2,MD2,MW2)) + 
     -        Alfa/(2.D0*Pi)*(DBLE(B0i(bb1,MU2,MU2,MZ2))*gL3**2) + 
     -        (2*Alfa*MU2)/Pi*(DBLE(B0i(dbb0,MU2,MU2,MZ2))*gL3*
     -        gR3) + (Alfa*MU2)/(2.D0*Pi)*
     -        (DBLE(B0i(dbb1,MU2,MU2,MZ2))*(gL3**2 + gR3**2))+
     -        Alfa/(8.D0*Pi*SW2)*
     -        (1 + SW2*(2*gL3**2 + 2*Qu**2)) - 
     -        Alfa/(2.D0*Pi)*((DBLE(B0i(bb0,MU2,0.D0,MU2)) + 
     -        DBLE(B0i(bb1,MU2,0.D0,MU2)))*Qu**2) + 
     -        (Alfa*MU2)/Pi*(DBLE(B0i(dbb0,MU2,0.D0,MU2))*Qu**2) - 
     -        (Alfa*MU2)/Pi*(DBLE(B0i(dbb1,MU2,0.D0,MU2))*Qu**2)

         dZfLd = Alfa/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MD2,MU2,MW2)) + 
     -        (Alfa*MD2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MD2,MU2,MW2)) + 
     -        Alfa/(2.D0*Pi)*(DBLE(B0i(bb1,MD2,MD2,MZ2))*gL4**2) + 
     -        (2*Alfa*MD2)/Pi*(DBLE(B0i(dbb0,MD2,MD2,MZ2))*gL4*
     -        gR4) + (Alfa*MD2)/(2.D0*Pi)*
     -        (DBLE(B0i(dbb1,MD2,MD2,MZ2))*(gL4**2 + gR4**2))+
     -        Alfa/(8.D0*Pi*SW2)*
     -        (1 + SW2*(2*gL4**2 + 2*Qd**2)) - 
     -        Alfa/(2.D0*Pi)*((DBLE(B0i(bb0,MD2,0.D0,MD2)) + 
     -        DBLE(B0i(bb1,MD2,0.D0,MD2)))*Qd**2) + 
     -        (Alfa*MD2)/Pi*(DBLE(B0i(dbb0,MD2,0.D0,MD2))*Qd**2) - 
     -        (Alfa*MD2)/Pi*(DBLE(B0i(dbb1,MD2,0.D0,MD2))*Qd**2)

         dZfLc = Alfa/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MC2,MS2,MW2)) + 
     -        (Alfa*MC2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MC2,MS2,MW2)) + 
     -        Alfa/(2.D0*Pi)*(DBLE(B0i(bb1,MC2,MC2,MZ2))*gL3**2) + 
     -        (2*Alfa*MC2)/Pi*(DBLE(B0i(dbb0,MC2,MC2,MZ2))*gL3*
     -        gR3) + (Alfa*MC2)/(2.D0*Pi)*
     -        (DBLE(B0i(dbb1,MC2,MC2,MZ2))*(gL3**2 + gR3**2))+
     -        Alfa/(8.D0*Pi*SW2)*
     -        (1 + SW2*(2*gL3**2 + 2*Qu**2)) - 
     -        Alfa/(2.D0*Pi)*((DBLE(B0i(bb0,MC2,0.D0,MC2)) + 
     -        DBLE(B0i(bb1,MC2,0.D0,MC2)))*Qu**2) + 
     -        (Alfa*MC2)/Pi*(DBLE(B0i(dbb0,MC2,0.D0,MC2))*Qu**2) - 
     -        (Alfa*MC2)/Pi*(DBLE(B0i(dbb1,MC2,0.D0,MC2))*Qu**2)

         dZfLs = Alfa/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MS2,MC2,MW2)) + 
     -        (Alfa*MS2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MS2,MC2,MW2)) + 
     -        Alfa/(2.D0*Pi)*(DBLE(B0i(bb1,MS2,MS2,MZ2))*gL4**2) + 
     -        (2*Alfa*MS2)/Pi*(DBLE(B0i(dbb0,MS2,MS2,MZ2))*gL4*
     -        gR4) + (Alfa*MS2)/(2.D0*Pi)*
     -        (DBLE(B0i(dbb1,MS2,MS2,MZ2))*(gL4**2 + gR4**2))+
     -        Alfa/(8.D0*Pi*SW2)*
     -        (1 + SW2*(2*gL4**2 + 2*Qd**2)) - 
     -        Alfa/(2.D0*Pi)*((DBLE(B0i(bb0,MS2,0.D0,MS2)) + 
     -        DBLE(B0i(bb1,MS2,0.D0,MS2)))*Qd**2) + 
     -        (Alfa*MS2)/Pi*(DBLE(B0i(dbb0,MS2,0.D0,MS2))*Qd**2) - 
     -        (Alfa*MS2)/Pi*(DBLE(B0i(dbb1,MS2,0.D0,MS2))*Qd**2)


** right-hand fermions field renormalisation
         dup5 = DBLE(B0i(bb0,MU2,0.D0,MU2)) +DBLE(B0i(bb1,MU2,0.D0,MU2))

         dZfRu = Alfa/(9.D0*CW2*Pi) - 
     -        (2*Alfa)/(9.D0*Pi)*dup5 + 
     -        (2*Alfa*SW2)/(9.D0*CW2*Pi)*DBLE(B0i(bb1,MU2,MU2,MZ2)) + 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb0,MU2,0.D0,MU2)) - 
     -        (2*Alfa*MU2*(3 - 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MU2,MU2,MZ2)) - 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb1,MU2,0.D0,MU2)) + 
     -        (Alfa*MU2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MU2,MD2,MW2)) + 
     -        (Alfa*MU2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MU2,MU2,MZ2))

         dup5 = DBLE(B0i(bb0,MD2,0.D0,MD2)) +DBLE(B0i(bb1,MD2,0.D0,MD2))

         dZfRd = Alfa/(36.D0*CW2*Pi) - Alfa/(18.D0*Pi)*dup5 + 
     -        (Alfa*SW2)/(18.D0*CW2*Pi)*DBLE(B0i(bb1,MD2,MD2,MZ2)) + 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb0,MD2,0.D0,MD2)) - 
     -        (Alfa*MD2*(3 - 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MD2,MD2,MZ2)) - 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb1,MD2,0.D0,MD2)) + 
     -        (Alfa*MD2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MD2,MD2,MZ2)) + 
     -        (Alfa*MD2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MD2,MU2,MW2))

         dup6 = DBLE(B0i(bb0,MC2,0.D0,MC2)) +DBLE(B0i(bb1,MC2,0.D0,MC2))

         dZfRc = Alfa/(9.D0*CW2*Pi) - (2*Alfa)/(9.D0*Pi)*dup6 + 
     -        (2*Alfa*SW2)/(9.D0*CW2*Pi)*DBLE(B0i(bb1,MC2,MC2,MZ2)) + 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb0,MC2,0.D0,MC2)) - 
     -        (2*Alfa*MC2*(3 - 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MC2,MC2,MZ2)) - 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb1,MC2,0.D0,MC2)) + 
     -        (Alfa*MC2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MC2,MC2,MZ2)) + 
     -        (Alfa*MC2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MC2,MS2,MW2))

         dup6 = DBLE(B0i(bb0,MS2,0.D0,MS2)) +DBLE(B0i(bb1,MS2,0.D0,MS2))

         dZfRs = Alfa/(36.D0*CW2*Pi) - Alfa/(18.D0*Pi)*dup6 + 
     -        (Alfa*SW2)/(18.D0*CW2*Pi)*DBLE(B0i(bb1,MS2,MS2,MZ2)) + 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb0,MS2,0.D0,MS2)) - 
     -        (Alfa*MS2*(3 - 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MS2,MS2,MZ2)) - 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb1,MS2,0.D0,MS2)) + 
     -        (Alfa*MS2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MS2,MC2,MW2)) + 
     -        (Alfa*MS2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MS2,MS2,MZ2))

      end if  ! End of bosonic contributions


      dSW1 = CW2/(2.D0*MW2*MZ2)*(-(MZ2*dMWsq1) + MW2*dMZsq1)/SW


** Charge renormalisation: see notebook 8, page 67.  
      if (charge .ne. 2) then
         dZet1 = 1/(2.D0*CW)*(-(dZZA1*SW) + CW*(-dZAA1 - dZAAlt1))
      else if (charge .eq. 2) then
         dZet1 = dZeR + 1/(2.D0*MW2)*((dMWsq1*SW + 2*MW2*dSW1))/SW
         if (sector .eq. 3) then
         dZet1 = dZet1 - 1/(2.D0*CW*MW2)*(MW2*dZZA1)/SW - Alfa*(6.D0 + 
     -           (7.D0 - 4*SW2)*LOG(CW2)/(2.D0*SW2))/(8.D0*Pi*SW2)
         end if
      end if

      if (charge .eq. 0) then
         dZet1 = dZet1 + (delAlfa/2.D0)
*         dZet1 = dZet1 + (0.02758D0 + 0.031498D0)/2.D0
      end if


** fermion values of dSW1 and dZet1
      dSW1fm = CW2/(2.D0*MW2*MZ2)*(-(MZ2*dMWsq1fm) + MW2*dMZsq1fm)/SW
* Charge renormalisation: see notebook 8, page 67.  
      if (charge .ne. 2) then
         dZet1fm = 1/(2.D0*CW)*(-(dZZA1fm*SW) + CW*(-dZAA1fm - dZAAlt1fm))
      else if (charge .eq. 2) then
         dZet1fm = dZeRfm + 1/(2.D0*MW2)*((dMWsq1fm*SW + 2*MW2*dSW1fm))/SW
      end if
      if (charge .eq. 0) then
         dZet1fm = dZet1fm + (delAlfa/2.D0)
      end if


      if (debug .ne. 0) then
         write(*,*)'  '
         write(*,*)'Renormalisation constants ...'
         write(*,*)'dMWsq1 =', dMWsq1
         write(*,*)'dMZsq1 =', dMZsq1
         write(*,*)'dSW1 =', dSW1
         write(*,*)'dZet1 =', dZet1
         write(*,*)'dZW1 =', dZW1
         write(*,*)'dZZZ1 =', dZZZ1
         write(*,*)'dZH1 =', dZH1
         write(*,*)'dZAA1 =', dZAA1
         write(*,*)'dZAAlt1 =', dZAAlt1
         write(*,*)'dZZA1 =', dZZA1 
         write(*,*)'dZAZ1 =', dZAZ1
         write(*,*)'dZfLu =', dZfLu
         write(*,*)'dZfLd =', dZfLd
         write(*,*)'dZfLc =', dZfLc
         write(*,*)'dZfLs =', dZfLs
         write(*,*)'dZfRu =', dZfRu
         write(*,*)'dZfRd =', dZfRd
         write(*,*)'dZfRc =', dZfRc
         write(*,*)'dZfRs =', dZfRs
      
 
** Debugging output: fermion values 
         write(*,*)'         '
         write(*,*)'FERMION VALUES!:'
         write(*,*)'dMWsq1 =', dMWsq1fm
         write(*,*)'dMZsq1 =', dMZsq1fm
         write(*,*)'dSW1 =', dSW1fm
         write(*,*)'dZet1 =', dZet1fm
         write(*,*)'dZW1 =', dZW1fm
         write(*,*)'dZZZ1 =', dZZZ1fm
         write(*,*)'dZH1 =', dZH1fm
         write(*,*)'dZAA1 =', dZAA1fm
         write(*,*)'dZAAlt1 =', dZAAlt1fm
         write(*,*)'dZZA1 =', dZZA1fm
         write(*,*)'dZAZ1 =', dZAZ1fm
      end if


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the renormalisation constants in SUSY - for 
*** definitions, see notebook 5, page 28

      subroutine ren_conSUSY(charge,debug,full)

      implicit none


** Do loop parameters 
      integer Gen3, Sfe3, Sfe4, g3ini
      integer Cha3, Cha4, Neu3, Neu4

** Dummy parameters
      double complex dup4, dup7, dup8, dup21, dup22
      double complex chi9, chi10, dup5
      double complex dup13, dup14, dup15, dup16
      double complex dup18, dup19
      double complex eta24, eta25, eta26, eta27, eta28, eta29, eta30
      double complex eta31, eta32, eta33, eta34, eta35, eta36, eta37
      double complex eta38, eta39, eta40, eta41, eta42, eta43, eta44
      double complex eta45, eta46
      double complex dup11, dup12, dup40

      double complex dup23, dup24, dup25, dup26, dup27, dup28
      double complex dup29, dup30, dup31, dup32, dup33, dup35
      double complex dup36, dup37, dup38, dup39, dup41
      double complex dup42, dup43, dup44, dup45, dup46
      double complex dup47, dup48, dup49, dup50, dup51
      double complex dup52, dup53, dup54, dup55, dup56
      double complex dup57, dup58, dup59, dup60

      double precision IndexDelta(4,4)

      double complex dZAAlt1, dZeR
      double complex dZAAlt1fm, dZeRfm
      double complex dMWsq1sm, dSW1sm, dZZA1sm, dZeRsm
      double complex dZZA1chi, dZAA1chi


      double precision Alfa, Alfa2

      integer charge

** Flag determining whether contributions from neutralinos/charginos are 
** calculated and flag determining whether we need to add the suffix 'full'
** to the renormalisation constant names
      logical chiRUN, full

** Debugging flag
      integer debug


** Model parameter declarations
#include "mssm.inc"   
#include "global.inc"   

** This is needed for looptools. 
#include "looptools.h"


** Clearing the values of the renormalization constants
      dMWsq1 = (0D0,0D0)
      dMZsq1 = (0D0,0D0)
      dSW1 = (0D0,0D0)
      dZet1 = (0D0,0D0)
      dZW1 = (0D0,0D0)
      dZZZ1 = (0D0,0D0)
      dZH1 = (0D0,0D0)
      dZh0 = (0D0,0D0)
      dZHH = (0D0,0D0)
      dZh0HH = (0D0,0D0)
      dMAZsq1 = (0D0,0D0)
      dTB1 = (0D0,0D0)
      dZAA1 = (0D0,0D0)
      dZZA1 = (0D0,0D0)
      dZAZ1 = (0D0,0D0)
      dZfLu = (0D0,0D0)
      dZfLd = (0D0,0D0)
      dZfLc = (0D0,0D0)
      dZfLs = (0D0,0D0)
      dZfRu = (0D0,0D0)
      dZfRd = (0D0,0D0)
      dZfRc = (0D0,0D0)
      dZfRs = (0D0,0D0)

** Neutralinos are included in the renormalisation constant calculations if
** we are running the full MSSM (ewapprox = 0), if we are running the full MSSM
** at hVV, VV and qqV with only SM-type boxes (ewapprox = 5) or we are 
** calculating 'rencon_FULL', used in the calculation of the Higgs vertex while
** other contributions are SM-type or sfermions (ewapprox = 4)
      chiRUN  = full .or. (sector .ge. 5)

** Sets the correct starting value for the loops over generation for the quarks.
** The parameter 'sector' is set in input.dat (line 9) and is 1 for generation 3
** quarks, 2 for all fermions and 3 for all particles.  
      if ((sector .eq. 1) .and. (.not. full)) then
         g3ini = 3
      else
         g3ini = 1
      end if

      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      do Neu3 = 1, 4
         do Neu4 = 1, 4
            if (Neu3 .eq. Neu4) then
               IndexDelta(Neu3,Neu4) = 1D0
            else
               IndexDelta(Neu3,Neu4) = 0D0
            end if
         end do
      end do


      dZAAlt1 = (0D0,0D0)
      dZeR = (0D0,0D0)
      
      dMWsq1sm = (0D0,0D0)
      dSW1sm = (0D0,0D0)
      dZZA1sm = (0D0,0D0)
      dZeRsm = (0D0,0D0)
      dZZA1chi = (0D0,0D0)
      dZAA1chi = (0D0,0D0)


** We start with contributions from the (s)quark sector

** quarks
      do Gen3 = g3ini, 3

** W boson mass renormalisation: on-shell
         dMWsq1 = dMWsq1 - (3*Alfa)/(4.D0*Pi*SW2)*
     -        DBLE(A0(Mf2(4,Gen3))) + (3*Alfa)/(2.D0*Pi*SW2)*
     -        DBLE(B0i(bb00,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa*MW2)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(4.D0*Pi*SW2)*
     -        (DBLE(B0i(bb0,MW2,Mf2(3,Gen3),Mf2(4,Gen3)))*Mf2(3,Gen3)) 

** Z boson mass renormalisation: on-shell
         dMZsq1 = dMZsq1 - (Alfa*(9 + SW2*(-24 + 32*SW2)))/
     -        (24.D0*CW2*Pi*SW2)*DBLE(A0(Mf2(3,Gen3))) - 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(A0(Mf2(4,Gen3))) - 
     -        (Alfa*(-9 + (24 - 32*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -        (Alfa*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -        (Alfa*MZ2*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -        (Alfa*MZ2*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))*Mf2(3,Gen3))-
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3))
            
** Z-photon mixing
            dZZA1 = dZZA1 + 2/MZ2*(-((Alfa*(3 - 8*SW2))/(6.D0*CW*Pi)*
     -           DBLE(A0(Mf2(3,Gen3)))/SW) - (Alfa*(3 - 4*SW2))/
     -           (12.D0*CW*Pi)*DBLE(A0(Mf2(4,Gen3)))/SW - 
     -           (Alfa*(-3 + 8*SW2))/(3.D0*CW*Pi)*
     -           DBLE(B0i(bb00,0.D0,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -           (Alfa*(-3 + 4*SW2))/(6.D0*CW*Pi)*
     -           DBLE(B0i(bb00,0.D0,Mf2(4,Gen3),Mf2(4,Gen3)))/SW)

** Charge renormalisation  
         if (charge .ne. 2) then ! parametrising lo by alpha(MZ) or alpha(0) 

** Top quark :
            if (Gen3 .eq. 3) then
               dZAA1 = dZAA1 + (4*Alfa)/(3.D0*Pi)*
     -              DBLE(B0i(bb1,0.D0,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -              (8*Alfa)/(3.D0*Pi)*DBLE(B0i(dbb00,0.D0,
     -              Mf2(3,Gen3),Mf2(3,Gen3))) 
            end if

         else   ! charge renormalisation: parametrising lo by GF
         
            dZeR = dZeR - Alfa/(8.D0*MW2*Pi*SW2)*
     -           (-3*A0(Mf2(4,Gen3)) + 
     -           6*B0i(bb00,0.D0,Mf2(3,Gen3),Mf2(4,Gen3)) - 
     -           3*B0i(bb0,0.D0,Mf2(3,Gen3),Mf2(4,Gen3))*Mf2(3,Gen3))

         end if   ! end of charge renormalisation


** Higgs field renormalisation: DR bar
         dZh0 = dZh0 - Alfa/(8.D0*MW2*Pi*SW2)*(DBLE(Divergence)*
     -        (3*CA2*CB2*Mf2(3,Gen3) + 3*SA2*SB2*Mf2(4,Gen3)))/
     -        (CB2*SB2)         

         dZHH = dZHH - Alfa/(8.D0*MW2*Pi*SW2)*(DBLE(Divergence)*
     -        (3*CB2*SA2*Mf2(3,Gen3) + 3*CA2*SB2*Mf2(4,Gen3)))/
     -        (CB2*SB2)    

** W boson field renormalisation: on-shell
         dZW1 = dZW1 + (3*Alfa)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) - 
     -        (3*Alfa)/(2.D0*Pi*SW2)*
     -        DBLE(B0i(dbb00,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa*MW2)/(4.D0*Pi*SW2)*
     -        DBLE(B0i(dbb1,MW2,Mf2(3,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa)/(4.D0*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MW2,Mf2(3,Gen3),Mf2(4,Gen3)))*Mf2(3,Gen3))


** Z boson field renormalisation: on-shell
         dZZZ1 = dZZZ1 + (Alfa*(9 + SW2*(-24 + 32*SW2)))/
     -        (24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (Alfa*(-9 + (24 - 32*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*(-9 + (12 - 8*SW2)*SW2))/(12.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (Alfa*MZ2*(9 + SW2*(-24 + 32*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) + 
     -        (Alfa*MZ2*(9 + SW2*(-12 + 8*SW2)))/(24.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -        (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))*Mf2(3,Gen3))
     -        + (3*Alfa)/(8.D0*CW2*Pi*SW2)*
     -        (DBLE(B0i(dbb0,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))*Mf2(4,Gen3)) 

** photon-Z field renormalisation
         dZAZ1 = dZAZ1 - 2/MZ2*(-((Alfa*(3 - 8*SW2))/(6.D0*CW*Pi)*
     -        DBLE(A0(Mf2(3,Gen3)))/SW) - (Alfa*(3 - 4*SW2))/
     -        (12.D0*CW*Pi)*DBLE(A0(Mf2(4,Gen3)))/SW - 
     -        (Alfa*(-3 + 8*SW2))/(3.D0*CW*Pi)*
     -        DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -        (Alfa*(-3 + 4*SW2))/(6.D0*CW*Pi)*
     -        DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))/SW - 
     -        (Alfa*MZ2*(3 - 8*SW2))/(6.D0*CW*Pi)*
     -        DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3)))/SW - 
     -        (Alfa*MZ2*(3 - 4*SW2))/(12.D0*CW*Pi)*
     -        DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3)))/SW)

      end do


** Contributions from squarks 
      do Sfe3 = 1,2
         do Gen3 = g3ini, 3            
            dMWsq1 = dMWsq1 + (3*Alfa)/(8.D0*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,3,Gen3)))*USf(Sfe3,1,3,Gen3)*
     -           USfC(Sfe3,1,3,Gen3)) + (3*Alfa)/(8.D0*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,4,Gen3)))*USf(Sfe3,1,4,Gen3)*
     -           USfC(Sfe3,1,4,Gen3))

            dMZsq1 = dMZsq1 + Alfa/(24.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,3,Gen3)))*((3 - 4*SW2)**2*
     -           (USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 
     -           16*SW2**2*(USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3))))+
     -           Alfa/(24.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,4,Gen3)))*((3 - 2*SW2)**2*
     -           (USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)) + 
     -           4*SW2**2*(USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))) 
               
                  
            dup21 = SW2*(4*USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3)+ 
     -           4*USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)) - 
     -           3*USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3)                 
            dup22 = SW2*(2*USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)+ 
     -           2*USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)) - 
     -           3*USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)
            dZZA1 = dZZA1 + 2/MZ2*(-(Alfa/(6.D0*CW*Pi)*
     -           (dup21*DBLE(A0(MSf2(Sfe3,3,Gen3))))/SW) - 
     -           Alfa/(12.D0*CW*Pi)*
     -           (dup22*DBLE(A0(MSf2(Sfe3,4,Gen3))))/SW + 
     -           Alfa/(3.D0*CW*Pi)*
     -           (dup21*DBLE(B0i(bb00,0.D0,MSf2(Sfe3,3,Gen3),
     -           MSf2(Sfe3,3,Gen3))))/SW + Alfa/(6.D0*CW*Pi)*
     -           (dup22*DBLE(B0i(bb00,0.D0,MSf2(Sfe3,4,Gen3),
     -           MSf2(Sfe3,4,Gen3))))/SW)

** Charge renormalisation 
            if (charge .ne. 2) then ! LO via alpha(0), alpha(MZ)
               dZAA1 = dZAA1 + (4*Alfa)/(3.D0*Pi)*
     -              DBLE(B0i(dbb00,0.D0,MSf2(Sfe3,3,Gen3),
     -              MSf2(Sfe3,3,Gen3))) + Alfa/(3.D0*Pi)*
     -              DBLE(B0i(dbb00,0.D0,MSf2(Sfe3,4,Gen3),
     -              MSf2(Sfe3,4,Gen3)))
            else  ! LO via GF
               dZeR = dZeR - Alfa/(16.D0*MW2*Pi*SW2)*
     -              (3*A0(MSf2(Sfe3,3,Gen3))*USf(Sfe3,1,3,Gen3)*
     -              USfC(Sfe3,1,3,Gen3) + 3*A0(MSf2(Sfe3,4,Gen3))*
     -              USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)) 
            end if

            dup7 = SW2*(4*USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -           4*USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)) - 
     -           3*USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3)            
            dup8 = SW2*(2*USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -           2*USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)) - 
     -           3*USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3)

            dZAZ1 = dZAZ1 - 2/MZ2*
     -           (-(Alfa/(6.D0*CW*Pi)*
     -           (dup7*DBLE(A0(MSf2(Sfe3,3,Gen3))))/SW) - 
     -           Alfa/(12.D0*CW*Pi)*
     -           (dup8*DBLE(A0(MSf2(Sfe3,4,Gen3))))/SW + 
     -           Alfa/(3.D0*CW*Pi)*
     -           (dup7*DBLE(B0i(bb00,MZ2,MSf2(Sfe3,3,Gen3),
     -           MSf2(Sfe3,3,Gen3))))/SW + 
     -           Alfa/(6.D0*CW*Pi)*
     -           (dup8*DBLE(B0i(bb00,MZ2,MSf2(Sfe3,4,Gen3),
     -           MSf2(Sfe3,4,Gen3))))/SW)
         enddo
      enddo


** More contributions from squarks
      do Sfe4 = 1,2
         do Sfe3 = 1,2
            do Gen3 = g3ini, 3               
               dMWsq1 = dMWsq1 - (3*Alfa)/(2.D0*Pi*SW2)*
     -              (DBLE(B0i(bb00,MW2,MSf2(Sfe3,3,Gen3),
     -              MSf2(Sfe4,4,Gen3)))*
     -              USf(Sfe3,1,3,Gen3)*USf(Sfe4,1,4,Gen3)*
     -              USfC(Sfe3,1,3,Gen3)*USfC(Sfe4,1,4,Gen3))  

               dMZsq1 = dMZsq1 - 
     -              Alfa/(12.D0*CW2*Pi*SW2)*(DBLE(B0i(bb00,MZ2,
     -              MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -              (SW2*(4*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -              4*USf(Sfe4,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)) - 
     -              3*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3))*
     -              (SW2*(4*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3) + 
     -              4*USf(Sfe3,2,3,Gen3)*USfC(Sfe4,2,3,Gen3)) - 
     -              3*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3))) - 
     -              Alfa/(12.D0*CW2*Pi*SW2)*(DBLE(B0i(bb00,MZ2,
     -              MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -              (SW2*(2*USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -              2*USf(Sfe4,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)) - 
     -              3*USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3))*
     -              (SW2*(2*USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3) + 
     -              2*USf(Sfe3,2,4,Gen3)*USfC(Sfe4,2,4,Gen3)) - 
     -              3*USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3)))

** Charge renormalisation
               if (charge .eq. 2) then ! LO via GF
                  dZeR = 1/(4.D0*MW2*Pi*SW2)*(4*MW2*Pi*SW2*dZeR + 
     -                 3*Alfa*(B0i(bb00,0.D0,MSf2(Sfe3,3,Gen3),
     -                 MSf2(Sfe4,4,Gen3))*USf(Sfe3,1,3,Gen3)*
     -                 USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,3,Gen3)*
     -                 USfC(Sfe4,1,4,Gen3))) 
               end if
          
               dZW1 = dZW1 + (3*Alfa)/(2.D0*Pi*SW2)*
     -              (DBLE(B0i(dbb00,MW2,MSf2(Sfe3,3,Gen3),
     -              MSf2(Sfe4,4,Gen3)))*USf(Sfe3,1,3,Gen3)*
     -              USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,3,Gen3)*
     -              USfC(Sfe4,1,4,Gen3))  
          
               dZZZ1 = dZZZ1 + Alfa/(12.D0*CW2*Pi*SW2)*
     -              (DBLE(B0i(dbb00,MZ2,MSf2(Sfe3,3,Gen3),
     -              MSf2(Sfe4,3,Gen3)))*
     -              (SW2*(4*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -              4*USf(Sfe4,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)) - 
     -              3*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,1,3,Gen3))*
     -              (SW2*(4*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3) + 
     -              4*USf(Sfe3,2,3,Gen3)*USfC(Sfe4,2,3,Gen3)) - 
     -              3*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,1,3,Gen3))) + 
     -              Alfa/(12.D0*CW2*Pi*SW2)*
     -              (DBLE(B0i(dbb00,MZ2,MSf2(Sfe3,4,Gen3),
     -              MSf2(Sfe4,4,Gen3)))*
     -              (SW2*(2*USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -              2*USf(Sfe4,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)) - 
     -              3*USf(Sfe4,1,4,Gen3)*USfC(Sfe3,1,4,Gen3))*
     -              (SW2*(2*USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3) + 
     -              2*USf(Sfe3,2,4,Gen3)*USfC(Sfe4,2,4,Gen3)) - 
     -              3*USf(Sfe3,1,4,Gen3)*USfC(Sfe4,1,4,Gen3)))  
            enddo
         enddo
      enddo


** Charge renormalisation:  parametrising lo by alpha(MZ) or alpha(0)
      if (charge .ne. 2) then    

         do Gen3 = 1, 3

** Down-type quarks:
            dZAAlt1 = dZAAlt1 - 
     -           1/MZ2*(-Alfa/(3.D0*Pi)*DBLE(A0(Mf2(4,Gen3))) + 
     -           (2*Alfa)/(3.D0*Pi)*
     -           DBLE(B0i(bb00,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))) - 
     -           (Alfa*MZ2)/(3.D0*Pi)*
     -           DBLE(B0i(bb1,MZ2,Mf2(4,Gen3),Mf2(4,Gen3))))
            
** light up-type quarks
            if (Gen3 .ne. 3) then
               dZAAlt1 = dZAAlt1 - 
     -              1/MZ2*(-((4*Alfa)/(3.D0*Pi)*DBLE(A0(Mf2(3,Gen3)))) + 
     -              (8*Alfa)/(3.D0*Pi)*
     -              DBLE(B0i(bb00,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))) - 
     -              (4*Alfa*MZ2)/(3.D0*Pi)*
     -              DBLE(B0i(bb1,MZ2,Mf2(3,Gen3),Mf2(3,Gen3))))
            end if

** leptons
            dZAAlt1 = dZAAlt1 - 
     -           1/MZ2*(-(Alfa/Pi*DBLE(A0(Mf2(2,Gen3)))) + 
     -           (2*Alfa)/Pi*DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),
     -           Mf2(2,Gen3))) - 
     -           (Alfa*MZ2)/Pi*DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),
     -           Mf2(2,Gen3))))


** If we're only considering loop contributions from third generation
** (s)quarks, we need to calculate the rest of Delta Alfa.  See ntbk, 
** 30-7-10.
            if ((sector .eq. 1) .and. (.not. full)) then

** ! (nb: if gen=3, dZAA1 contribution has already been calculated)
               if (Gen3 .ne. 3) then 
** light up-type quarks
                  dZAA1 = dZAA1 - 
     -                 (4*Alfa)/(3.D0*Pi)*DBLE(B0i(bb1,0.D0,Mf2(3,Gen3),
     -                 Mf2(3,Gen3))) + 
     -                 (8*Alfa)/(3.D0*Pi)*DBLE(B0i(dbb00,0.D0,Mf2(3,Gen3),
     -                 Mf2(3,Gen3))) 
                  
** down-type quarks
                  dZAA1 = dZAA1 - Alfa/(3.D0*Pi)*
     -                 DBLE(B0i(bb1,0.D0,Mf2(4,Gen3),Mf2(4,Gen3))) + 
     -                 (2*Alfa)/(3.D0*Pi)*DBLE(B0i(dbb00,0.D0,Mf2(4,Gen3),
     -                 Mf2(4,Gen3)))
               end if

** leptons
               dZAA1 = dZAA1 - Alfa/Pi*DBLE(B0i(bb1,0.D0,Mf2(2,Gen3),
     -              Mf2(2,Gen3))) + (2*Alfa)/Pi*DBLE(B0i(dbb00,0.D0,
     -              Mf2(2,Gen3),Mf2(2,Gen3)))
               
            end if  ! end of if sector = 1
            
         end do

      end if
** end of charge renormalisation
      

** (S)LEPTONIC CONTRIBUTIONS
      if ((sector .ne. 1) .or. full) then

** Firstly, the leptons
        dMZsq1 = dMZsq1 + Alfa/(8.D0*CW2*Pi*SW2)*
     -        (6*DBLE(B0i(bb00,MZ2,0.D0,0.D0)) - 
     -        3*MZ2*DBLE(B0i(bb1,MZ2,0.D0,0.D0)))

         dZZZ1 = dZZZ1 + Alfa/(8.D0*CW2*Pi*SW2)*
     -        (3*DBLE(B0i(bb1,MZ2,0.D0,0.D0)) - 
     -        6*DBLE(B0i(dbb00,MZ2,0.D0,0.D0)) + 
     -        3*MZ2*DBLE(B0i(dbb1,MZ2,0.D0,0.D0)))

         do Gen3 = 1,3
            dMWsq1 = dMWsq1 - Alfa/(4.D0*Pi*SW2)*DBLE(A0(Mf2(2,Gen3))) + 
     -           Alfa/(2.D0*Pi*SW2)*DBLE(B0i(bb00,MW2, 
     -           0.D0,Mf2(2,Gen3))) - (Alfa*MW2)/(4.D0*Pi*SW2)*
     -           DBLE(B0i(bb1,MW2,0.D0,Mf2(2,Gen3)))

           dMZsq1 = dMZsq1 - (Alfa*(1 + SW2*(-4 + 8*SW2)))/
     -          (8.D0*CW2*Pi*SW2)*DBLE(A0(Mf2(2,Gen3))) - 
     -          (Alfa*(-1 + (4 - 8*SW2)*SW2))/(4.D0*CW2*Pi*SW2)*
     -          DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     -          (Alfa*MZ2*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2*Pi*SW2)*
     -          DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) - 
     -          Alfa/(8.D0*CW2*Pi*SW2)*
     -          (DBLE(B0i(bb0,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))*Mf2(2,Gen3))
               
           dZZA1 = dZZA1 + 2/MZ2*(-((Alfa*(1 - 4*SW2))/
     -          (4.D0*CW*Pi)*DBLE(A0(Mf2(2,Gen3)))/SW) - 
     -          (Alfa*(-1 + 4*SW2))/(2.D0*CW*Pi)*
     -          DBLE(B0i(bb00,0.D0,Mf2(2,Gen3),Mf2(2,Gen3)))/SW)
           
** Charge renormalisation 
            if (charge .eq. 2) then  ! LO via GF
               dZeR = dZeR - Alfa/(16.D0*MW2*Pi*SW2)*
     -              (-2*A0(Mf2(2,Gen3)) + A0(MSf2(1,1,Gen3)) + 
     -              4*B0i(bb00,0.D0,0.D0,Mf2(2,Gen3)))               
            end if

            dZh0 = 1/(8.D0*MW2*Pi*SW2)*(8*MW2*Pi*SW2*(CB2*dZh0) - 
     -           Alfa*(SA2*DBLE(Divergence)*Mf2(2,Gen3)))/CB2

            dZHH = 1/(8.D0*MW2*Pi*SW2)*(8*MW2*Pi*SW2*(CB2*dZHH) - 
     -           Alfa*(CA2*DBLE(Divergence)*Mf2(2,Gen3)))/CB2 

            dZW1 = dZW1 + Alfa/(4.D0*Pi*SW2)*
     -           DBLE(B0i(bb1,MW2,0.D0,Mf2(2,Gen3))) - Alfa/(2.D0*Pi*
     -           SW2)*DBLE(B0i(dbb00,MW2,0.D0,Mf2(2,Gen3))) +(Alfa*MW2)/
     -           (4.D0*Pi*SW2)*DBLE(B0i(dbb1,MW2,0.D0,Mf2(2,Gen3)))

            dZZZ1 = dZZZ1 + (Alfa*(1 + SW2*(-4 + 8*SW2)))/
     -           (8.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           (Alfa*(-1 + (4 - 8*SW2)*SW2))/(4.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(dbb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           (Alfa*MZ2*(1 + SW2*(-4 + 8*SW2)))/(8.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(dbb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))*Mf2(2,Gen3))

            dZAZ1 = dZAZ1 - 2/MZ2*
     -           (-((Alfa*(1 - 4*SW2))/(4.D0*CW*Pi)*
     -           DBLE(A0(Mf2(2,Gen3)))/SW) - 
     -           (Alfa*(-1 + 4*SW2))/(2.D0*CW*Pi)*
     -           DBLE(B0i(bb00,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))/SW - 
     -           (Alfa*MZ2*(1 - 4*SW2))/(4.D0*CW*Pi)*
     -           DBLE(B0i(bb1,MZ2,Mf2(2,Gen3),Mf2(2,Gen3)))/SW)
         end do


** now contributions from sleptons
         do Gen3 = 1,3
            dMWsq1 = dMWsq1 + 
     -           Alfa/(8.D0*Pi*SW2)*DBLE(A0(MSf2(1,1,Gen3))) 

           dMZsq1 = dMZsq1 + 
     -          Alfa/(8.D0*CW2*Pi*SW2)*DBLE(A0(MSf2(1,1,Gen3))) - 
     -          Alfa/(4.D0*CW2*Pi*SW2)*
     -          DBLE(B0i(bb00,MZ2,MSf2(1,1,Gen3),MSf2(1,1,Gen3))) 

            dZZZ1 = dZZZ1 + Alfa/(4.D0*CW2*Pi*SW2)*
     -           DBLE(B0i(dbb00,MZ2,MSf2(1,1,Gen3),MSf2(1,1,Gen3))) 
         enddo
         
         do Sfe3 = 1,2
            do Gen3 = 1,3               
               dMWsq1 = dMWsq1 + 
     -              Alfa/(8.D0*Pi*SW2)*(DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -              USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) -  Alfa/
     -              (2.D0*Pi*SW2)*(DBLE(B0i(bb00,MW2,MSf2(1,1,Gen3),
     -              MSf2(Sfe3,2,Gen3)))*USf(Sfe3,1,2,Gen3)*
     -              USfC(Sfe3,1,2,Gen3))

              dMZsq1 = dMZsq1 + Alfa/(8.D0*CW2*Pi*SW2)*
     -             (DBLE(A0(MSf2(Sfe3,2,Gen3)))*((1 - 2*SW2)**2*
     -             (USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -             4*SW2**2*(USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3))))
                  
              dup11 = SW2*(2*USf(Sfe3,1,2,Gen3)*
     -             USfC(Sfe3,1,2,Gen3) + 
     -             2*USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)) - 
     -             USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)
                  
              dZZA1 = dZZA1 + 2/MZ2*(-(Alfa/(4.D0*CW*Pi)*
     -             (dup11*DBLE(A0(MSf2(Sfe3,2,Gen3))))/SW) + 
     -             Alfa/(2.D0*CW*Pi)*
     -             (dup11*DBLE(B0i(bb00,0.D0,MSf2(Sfe3,2,Gen3),
     -             MSf2(Sfe3,2,Gen3))))/SW)

** Charge renormalisation
              if (charge .ne. 2) then ! LO via alfa(0), alfa(MZ)
                 dZAA1 = dZAA1 + Alfa/Pi*DBLE(B0i(dbb00,0.D0,
     -                MSf2(Sfe3,2,Gen3),MSf2(Sfe3,2,Gen3)))
              else  !  LO via GF
                  dZeR = dZeR - Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (A0(MSf2(Sfe3,2,Gen3))*USf(Sfe3,1,2,Gen3)*
     -                 USfC(Sfe3,1,2,Gen3) - 4*
     -                 B0i(bb00,0.D0,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3))*
     -                 USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3))
               end if

               dZW1 = dZW1 + Alfa/(2.D0*Pi*SW2)*(DBLE(B0i(dbb00,
     -              MW2,MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3)))*
     -              USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3))
   
               dup4 = SW2*(2*USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -              2*USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)) - 
     -              USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)               

               dZAZ1 = dZAZ1 - 2/MZ2*
     -              (-(Alfa/(4.D0*CW*Pi)*
     -              (dup4*DBLE(A0(MSf2(Sfe3,2,Gen3))))/SW) + 
     -              Alfa/(2.D0*CW*Pi)*
     -              (dup4*DBLE(B0i(bb00,MZ2,MSf2(Sfe3,2,Gen3),
     -              MSf2(Sfe3,2,Gen3))))/SW)
            enddo
         enddo


** More contributions from sleptons
        do Sfe4 = 1,2
           do Sfe3 = 1,2
              do Gen3 = 1,3                 
                 dMZsq1 = dMZsq1- Alfa/(4.D0*CW2*Pi*SW2)*(DBLE(B0i(bb00,
     -                MZ2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                (SW2*(2*USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -                2*USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)) - 
     -                USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3))*
     -                (SW2*(2*USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     -                2*USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3)) - 
     -                USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3)))

                  dZZZ1 = dZZZ1 + Alfa/(4.D0*CW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb00,MZ2,MSf2(Sfe3,2,Gen3),
     -                 MSf2(Sfe4,2,Gen3)))*
     -                 (SW2*(2*USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -                 2*USf(Sfe4,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)) - 
     -                 USf(Sfe4,1,2,Gen3)*USfC(Sfe3,1,2,Gen3))*
     -                 (SW2*(2*USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     -                 2*USf(Sfe3,2,2,Gen3)*USfC(Sfe4,2,2,Gen3)) - 
     -                 USf(Sfe3,1,2,Gen3)*USfC(Sfe4,1,2,Gen3)))
               enddo
            enddo
         enddo

** END OF (s)leptonic contributions
      end if


** fermion values
      dMWsq1fm = dMWsq1
      dMZsq1fm = dMZsq1
      dZW1fm = dZW1
      dZZZ1fm = dZZZ1
      dZH1fm = dZH1
      dMAZsq1fm = dMAZsq1
      dZAA1fm = dZAA1
      dZZA1fm = dZZA1
      dZAZ1fm = dZAZ1
      dZAAlt1fm = dZAAlt1
      dZeRfm = dZeR
      dZh0fm = dZh0
      dZHHfm = dZHH
      dZh0HHfm = dZh0HH
      dTB1fm = dTB1


** CONTRIBUTIONS FROM SM-TYPE BOSONS AND THEIR SUPERPARTNERS
      if ((sector .ge. 3) .or. full) then

         dMWsq1 = dMWsq1 - (Alfa*MW2)/(6.D0*Pi*SW2) - (Alfa*MW2*(1 - 
     -        CW2 - SW2))/(2.D0*Pi*SW2) - 
     -        (Alfa*(-7 + 4*SW2))/(8.D0*Pi*SW2)*DBLE(A0(MW2)) + 
     -        (Alfa*(1/4.D0 + CW2))/(4.D0*Pi*SW2)*DBLE(A0(MZ2)) - 
     -        (Alfa*MW2)/Pi*DBLE(B0i(bb0,MW2,0.D0,MW2)) + 
     -        ((-5*Alfa*CW2*MW2)/
     -        (4.D0*Pi*SW2) - (Alfa*MW2*(2*CW2**2 - SW2**2))/
     -        (4.D0*CW2*Pi*SW2))*DBLE(B0i(bb0,MW2,MW2,MZ2)) - 
     -        (2*Alfa)/Pi*DBLE(B0i(bb00,MW2,0.D0,MW2)) - 
     -        (Alfa*(1/4.D0 + 2*CW2))/
     -        (Pi*SW2)*DBLE(B0i(bb00,MW2,MW2,MZ2)) - 
     -        (Alfa*MW2)/(2.D0*Pi)*DBLE(B0i(bb1,MW2,0.D0,MW2)) - 
     -        (Alfa*CW2*MW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MW2,MW2,MZ2))
         
         dMWsq1 = dMWsq1 + Alfa/(16.D0*Pi*SW2)*
     -        (DBLE(A0(MA02)) + DBLE(A0(Mh02)) + DBLE(A0(MHH2))) + 
     -        Alfa/(8.D0*Pi*SW2)*DBLE(A0(MHp2)) + (Alfa*MW2)/
     -        (4.D0*Pi*SW2)*(SBA**2*DBLE(B0i(bb0,MW2,Mh02,MW2))) + 
     -        (Alfa*MW2)/(4.D0*Pi*SW2)*(CBA**2*
     -        DBLE(B0i(bb0,MW2,MHH2,MW2))) - 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb00,MW2,MA02,MHp2)) - 
     -        Alfa/(4.D0*Pi*SW2)*(SBA**2*(DBLE(B0i(bb00,MW2,Mh02,MW2)) + 
     -        DBLE(B0i(bb00,MW2,MHH2,MHp2)))) - Alfa/(4.D0*Pi*SW2)*
     -        (CBA**2*(DBLE(B0i(bb00,MW2,Mh02,MHp2)) + 
     -        DBLE(B0i(bb00,MW2,MHH2,MW2)))) 


         dMZsq1 = dMZsq1 - (Alfa*CW2*MZ2)/(6.D0*Pi*SW2) - (Alfa*
     -        (-SW2**2 + CW2*(-9*CW2 + 2*SW2)))/(8.D0*CW2*Pi*SW2)*
     -        DBLE(A0(MW2)) + Alfa/(16.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(MZ2))) + ((-5*Alfa*CW2*MZ2)/(4.D0*Pi*SW2) - 
     -        (Alfa*MW2*(CW2**2 - SW2**2))/(2.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(bb0,MZ2,MW2,MW2)) - (Alfa*(CW2*(9*CW2 - 
     -        2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,MW2,MW2)) - 
     -        (Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MZ2,MW2,MW2))

         dMZsq1 = dMZsq1 - (Alfa*
     -        (-SW2**2 + CW2*(-CW2 + 2*SW2)))/(8.D0*CW2*Pi*SW2)*
     -        DBLE(A0(MHp2)) + Alfa/(16.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(MA02)) + DBLE(A0(Mh02)) + DBLE(A0(MHH2))) + 
     -        (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*
     -        (SBA**2*DBLE(B0i(bb0,MZ2,Mh02,MZ2))) + 
     -        (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*(CBA**2*
     -        DBLE(B0i(bb0,MZ2,MHH2,MZ2))) - Alfa/(4.D0*CW2*Pi*SW2)*
     -        (SBA**2*(DBLE(B0i(bb00,MZ2,MA02,MHH2)) + 
     -        DBLE(B0i(bb00,MZ2,Mh02,MZ2)))) - Alfa/(4.D0*CW2*Pi*
     -        SW2)*(CBA**2*(DBLE(B0i(bb00,MZ2,MA02,Mh02)) + 
     -        DBLE(B0i(bb00,MZ2,MHH2,MZ2)))) - 
     -        (Alfa*(CW2*(CW2 - 2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb00,MZ2,MHp2,MHp2)) 

         if (charge .ne. 2) then ! LO via alfa(0) or alfa(MZ)
            dZAA1 = dZAA1 + Alfa*(1/(6.D0*Pi) + 
     -           5/(4.D0*Pi)*DBLE(B0i(bb0,0.D0,MW2,MW2)) + 
     -           1/(2.D0*Pi)*DBLE(B0i(bb1,0.D0,MW2,MW2)) + 
     -           1/Pi*DBLE(B0i(dbb00,0.D0,MHp2,MHp2)) + 
     -           3/Pi*DBLE(B0i(dbb00,0.D0,MW2,MW2)))
         else ! LO via GF
            dZeR = dZeR - Alfa/(32.D0*MW2*Pi*SW2)*(A0(MA02) + A0(Mh02) + 
     -           A0(MHH2)) - Alfa/(16.D0*MW2*Pi*SW2)*A0(MHp2) - 
     -           (Alfa*(7 - 4*CW2 - 4*SW2))/(16.D0*MW2*Pi*SW2)*A0(MW2) - 
     -           (Alfa*(1 + 12*CW2))/(32.D0*MW2*Pi*SW2)*A0(MZ2) - 
     -           Alfa/(8.D0*Pi)*B0i(bb0,0.D0,0.D0,MW2) - 
     -           Alfa/(8.D0*Pi*SW2)*(SBA**2*B0i(bb0,0.D0,Mh02,MW2)) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA**2*B0i(bb0,0.D0,MHH2,MW2)) - 
     -           (Alfa*(-2*CW2**2*MZ2 + MW2*SW2**2))/
     -           (8.D0*CW2*MW2*Pi*SW2)*B0i(bb0,0.D0,MW2,MZ2) + 
     -           Alfa/(MW2*Pi)*B0i(bb00,0.D0,0.D0,MW2) + 
     -           Alfa/(8.D0*MW2*Pi*SW2)*B0i(bb00,0.D0,MA02,MHp2) + 
     -           Alfa/(8.D0*MW2*Pi*SW2)*(SBA**2*
     -           (B0i(bb00,0.D0,Mh02,MW2) + B0i(bb00,0.D0,MHH2,MHp2))) +
     -           Alfa/(8.D0*MW2*Pi*SW2)*(CBA**2*
     -           (B0i(bb00,0.D0,Mh02,MHp2) + B0i(bb00,0.D0,MHH2,MW2))) +
     -           (Alfa*(1 + 8*CW2))/(8.D0*MW2*Pi*SW2)*
     -           B0i(bb00,0.D0,MW2,MZ2)
         end if
         
         dZZA1 = dZZA1 -((Alfa*(-CW2 + SW2))/
     -        (2.D0*CW*MZ2*Pi)*DBLE(A0(MHp2))/SW) - (Alfa*
     -        (-5*CW2 + SW2))/(2.D0*CW*MZ2*Pi)*DBLE(A0(MW2))/SW - 
     -        (Alfa*MW2)/(CW*MZ2*Pi)*DBLE(B0i(bb0,0.D0,MW2,MW2))/SW - 
     -        (Alfa*(CW2 - SW2))/(CW*MZ2*Pi)*
     -        DBLE(B0i(bb00,0.D0,MHp2,MHp2))/SW - 
     -        (Alfa*(5*CW2 - SW2))/(CW*MZ2*Pi)*
     -        DBLE(B0i(bb00,0.D0,MW2,MW2))/SW


         eta34 = Alfa/(8.D0*CW2*Pi*SW2)*(SBA**2*DBLE(Divergence) + 
     -        CBA**2*(DBLE(Divergence) + 2*CW2*DBLE(Divergence)) + 
     -        2*CW2*(SBA**2*DBLE(Divergence)))
         
         dZh0 = dZh0 + eta34         
         dZHH = dZHH + eta34


         dZW1 = dZW1 + Alfa/(6.D0*Pi*SW2) + 
     -        (5*Alfa)/(4.D0*Pi)*DBLE(B0i(bb0,MW2,0.D0,MW2)) + 
     -        (5*Alfa*CW2)/(4.D0*Pi*SW2)*DBLE(B0i(bb0,MW2,MW2,MZ2)) + 
     -        Alfa/(2.D0*Pi)*DBLE(B0i(bb1,MW2,0.D0,MW2)) + 
     -        (Alfa*CW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MW2,MW2,MZ2)) + 
     -        (Alfa*MW2)/Pi*DBLE(B0i(dbb0,MW2,0.D0,MW2)) + 
     -        ((5*Alfa*CW2*MW2)/(4.D0*Pi*SW2) + 
     -        (Alfa*MW2*(2*CW2**2 - SW2**2))/(4.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(dbb0,MW2,MW2,MZ2)) + 
     -        (2*Alfa)/Pi*DBLE(B0i(dbb00,MW2,0.D0,MW2)) + 
     -        (Alfa*(1/4.D0 + 2*CW2))/(Pi*SW2)*
     -        DBLE(B0i(dbb00,MW2,MW2,MZ2)) + 
     -        (Alfa*MW2)/(2.D0*Pi)*DBLE(B0i(dbb1,MW2,0.D0,MW2)) + 
     -        (Alfa*CW2*MW2)/(2.D0*Pi*SW2)*DBLE(B0i(dbb1,MW2,MW2,MZ2))

         dZW1 = dZW1 - (Alfa*MW2)/(4.D0*Pi*SW2)*
     -        (SBA**2*DBLE(B0i(dbb0,MW2,Mh02,MW2))) - 
     -        (Alfa*MW2)/(4.D0*Pi*SW2)*
     -        (CBA**2*DBLE(B0i(dbb0,MW2,MHH2,MW2))) + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(dbb00,MW2,MA02,MHp2)) + 
     -        Alfa/(4.D0*Pi*SW2)*(SBA**2*
     -        (DBLE(B0i(dbb00,MW2,Mh02,MW2)) + 
     -        DBLE(B0i(dbb00,MW2,MHH2,MHp2)))) + 
     -        Alfa/(4.D0*Pi*SW2)*(CBA**2*
     -        (DBLE(B0i(dbb00,MW2,Mh02,MHp2)) + 
     -        DBLE(B0i(dbb00,MW2,MHH2,MW2)))) 


         dZZZ1 = dZZZ1 + (Alfa*CW2)/(6.D0*Pi*SW2) + 
     -        (5*Alfa*CW2)/(4.D0*Pi*SW2)*DBLE(B0i(bb0,MZ2,MW2,MW2)) + 
     -        (Alfa*CW2)/(2.D0*Pi*SW2)*DBLE(B0i(bb1,MZ2,MW2,MW2)) + 
     -        ((5*Alfa*CW2*MZ2)/(4.D0*Pi*SW2) + 
     -        (Alfa*MW2*(CW2**2 - SW2**2))/(2.D0*CW2*Pi*SW2))*
     -        DBLE(B0i(dbb0,MZ2,MW2,MW2)) + 
     -        (Alfa*(CW2*(9*CW2 - 2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,MW2,MW2)) + 
     -        (Alfa*CW2*MZ2)/(2.D0*Pi*SW2)*DBLE(B0i(dbb1,MZ2,MW2,MW2))

         dZZZ1 = dZZZ1 - (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*
     -        (SBA**2*DBLE(B0i(dbb0,MZ2,Mh02,MZ2))) - 
     -        (Alfa*MW2)/(4.D0*CW2**2*Pi*SW2)*
     -        (CBA**2*DBLE(B0i(dbb0,MZ2,MHH2,MZ2))) + 
     -        Alfa/(4.D0*CW2*Pi*SW2)*
     -        (SBA**2*(DBLE(B0i(dbb00,MZ2,MA02,MHH2)) + 
     -        DBLE(B0i(dbb00,MZ2,Mh02,MZ2)))) + 
     -        Alfa/(4.D0*CW2*Pi*SW2)*
     -        (CBA**2*(DBLE(B0i(dbb00,MZ2,MA02,Mh02)) + 
     -        DBLE(B0i(dbb00,MZ2,MHH2,MZ2)))) + 
     -        (Alfa*(CW2*(CW2 - 2*SW2) + SW2**2))/(4.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb00,MZ2,MHp2,MHp2)) 


         dZAZ1 = dZAZ1 + (Alfa*CW)/(3.D0*Pi)/SW + 
     -        (Alfa*(-CW2 + SW2))/(2.D0*CW*MZ2*Pi)*DBLE(A0(MHp2))/SW - 
     -        (Alfa*(5*CW2 - SW2))/(2.D0*CW*MZ2*Pi)*DBLE(A0(MW2))/SW - 
     -        (Alfa*(-CW2 + SW2))/(CW*MZ2*Pi)*
     -        DBLE(B0i(bb00,MZ2,MHp2,MHp2))/SW - 
     -        (Alfa*(-5*CW2 + SW2))/(CW*MZ2*Pi)*
     -        DBLE(B0i(bb00,MZ2,MW2,MW2))/SW + 
     -        (Alfa*CW)/Pi*DBLE(B0i(bb1,MZ2,MW2,MW2))/SW - 
     -        2/MZ2*(DBLE(B0i(bb0,MZ2,MW2,MW2))*
     -        (-((Alfa*MW2)/(2.D0*CW*Pi)/SW) - 
     -        (5*Alfa*CW*MZ2)/(4.D0*Pi)/SW))


** fermion field renormalisation : on shell
         dup25 = DBLE(B0i(bb0,MU2,0.D0,MU2))+DBLE(B0i(bb1,MU2,0.D0,MU2))

         dZfLu = (Alfa*(3 - 4*SW2)**2)/(144.D0*CW2*Pi*SW2) + 
     -        (Alfa*(9 + 8*SW2))/(72.D0*Pi*SW2) - 
     -        (2*Alfa)/(9.D0*Pi)*dup25 + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb1,MU2,MD2,MW2)) + 
     -        (Alfa*(3 - 4*SW2)**2)/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MU2,MU2,MZ2)) + 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb0,MU2,0.D0,MU2)) + 
     -        (2*Alfa*MU2*(-3 + 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MU2,MU2,MZ2)) - 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb1,MU2,0.D0,MU2)) + 
     -        (Alfa*MU2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MU2,MD2,MW2)) + 
     -        (Alfa*MU2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MU2,MU2,MZ2))

         dZfRu = Alfa/(9.D0*CW2*Pi) - (2*Alfa)/(9.D0*Pi)*dup25 + 
     -        (2*Alfa*SW2)/(9.D0*CW2*Pi)*DBLE(B0i(bb1,MU2,MU2,MZ2)) + 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb0,MU2,0.D0,MU2)) + 
     -        (2*Alfa*MU2*(-3 + 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MU2,MU2,MZ2)) - 
     -        (4*Alfa*MU2)/(9.D0*Pi)*DBLE(B0i(dbb1,MU2,0.D0,MU2)) + 
     -        (Alfa*MU2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MU2,MD2,MW2)) + 
     -        (Alfa*MU2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MU2,MU2,MZ2))

         dup25 = DBLE(B0i(bb0,MD2,0.D0,MD2))+DBLE(B0i(bb1,MD2,0.D0,MD2))
         
         dZfLd = (Alfa*(3 - 2*SW2)**2)/(144.D0*CW2*Pi*SW2) + 
     -        (Alfa*(9 + 2*SW2))/(72.D0*Pi*SW2) - Alfa/(18.D0*Pi)*dup25 + 
     -        (Alfa*(3 - 2*SW2)**2)/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MD2,MD2,MZ2)) + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb1,MD2,MU2,MW2)) + 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb0,MD2,0.D0,MD2)) + 
     -        (Alfa*MD2*(-3 + 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MD2,MD2,MZ2)) - 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb1,MD2,0.D0,MD2)) + 
     -        (Alfa*MD2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MD2,MD2,MZ2)) + 
     -        (Alfa*MD2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MD2,MU2,MW2))

         dZfRd = Alfa/(36.D0*CW2*Pi) - Alfa/(18.D0*Pi)*dup25 + 
     -        (Alfa*SW2)/(18.D0*CW2*Pi)*DBLE(B0i(bb1,MD2,MD2,MZ2)) + 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb0,MD2,0.D0,MD2)) + 
     -        (Alfa*MD2*(-3 + 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MD2,MD2,MZ2)) - 
     -        (Alfa*MD2)/(9.D0*Pi)*DBLE(B0i(dbb1,MD2,0.D0,MD2)) + 
     -        (Alfa*MD2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MD2,MD2,MZ2)) + 
     -        (Alfa*MD2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MD2,MU2,MW2))

         dup25 = DBLE(B0i(bb0,MC2,0.D0,MC2))+DBLE(B0i(bb1,MC2,0.D0,MC2))

         dZfLc = (Alfa*(3 - 4*SW2)**2)/(144.D0*CW2*Pi*SW2) + 
     -        (Alfa*(9 + 8*SW2))/(72.D0*Pi*SW2) - 
     -        (2*Alfa)/(9.D0*Pi)*dup25 + 
     -        (Alfa*(3 - 4*SW2)**2)/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MC2,MC2,MZ2)) + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb1,MC2,MS2,MW2)) + 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb0,MC2,0.D0,MC2)) + 
     -        (2*Alfa*MC2*(-3 + 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MC2,MC2,MZ2)) - 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb1,MC2,0.D0,MC2)) + 
     -        (Alfa*MC2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MC2,MC2,MZ2)) + 
     -        (Alfa*MC2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MC2,MS2,MW2))

         dZfRc = Alfa/(9.D0*CW2*Pi) - (2*Alfa)/(9.D0*Pi)*dup25 + 
     -        (2*Alfa*SW2)/(9.D0*CW2*Pi)*DBLE(B0i(bb1,MC2,MC2,MZ2)) + 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb0,MC2,0.D0,MC2)) + 
     -        (2*Alfa*MC2*(-3 + 4*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MC2,MC2,MZ2)) - 
     -        (4*Alfa*MC2)/(9.D0*Pi)*DBLE(B0i(dbb1,MC2,0.D0,MC2)) + 
     -        (Alfa*MC2*(9 + SW2*(-24 + 32*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MC2,MC2,MZ2)) + 
     -        (Alfa*MC2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MC2,MS2,MW2))

         dup25 = DBLE(B0i(bb0,MS2,0.D0,MS2))+DBLE(B0i(bb1,MS2,0.D0,MS2))

         dZfLs = (Alfa*(3 - 2*SW2)**2)/(144.D0*CW2*Pi*SW2) + 
     -        (Alfa*(9 + 2*SW2))/(72.D0*Pi*SW2) - Alfa/(18.D0*Pi)*dup25 + 
     -        Alfa/(4.D0*Pi*SW2)*DBLE(B0i(bb1,MS2,MC2,MW2)) + 
     -        (Alfa*(3 - 2*SW2)**2)/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(bb1,MS2,MS2,MZ2)) + 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb0,MS2,0.D0,MS2)) + 
     -        (Alfa*MS2*(-3 + 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MS2,MS2,MZ2)) - 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb1,MS2,0.D0,MS2)) + 
     -        (Alfa*MS2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MS2,MC2,MW2)) + 
     -        (Alfa*MS2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MS2,MS2,MZ2))

         dZfRs = Alfa/(36.D0*CW2*Pi) - Alfa/(18.D0*Pi)*dup25 + 
     -        (Alfa*SW2)/(18.D0*CW2*Pi)*DBLE(B0i(bb1,MS2,MS2,MZ2)) + 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb0,MS2,0.D0,MS2)) + 
     -        (Alfa*MS2*(-3 + 2*SW2))/(9.D0*CW2*Pi)*
     -        DBLE(B0i(dbb0,MS2,MS2,MZ2)) - 
     -        (Alfa*MS2)/(9.D0*Pi)*DBLE(B0i(dbb1,MS2,0.D0,MS2)) + 
     -        (Alfa*MS2)/(4.D0*Pi*SW2)*DBLE(B0i(dbb1,MS2,MC2,MW2)) + 
     -        (Alfa*MS2*(9 + SW2*(-12 + 8*SW2)))/(72.D0*CW2*Pi*SW2)*
     -        DBLE(B0i(dbb1,MS2,MS2,MZ2))

         
         dZeRsm = dZeR
         dMWsq1sm = dMWsq1
         dSW1sm = -(CW2/(2.D0*MW2*MZ2)*(-(MZ2*dMWsq1) + MW2*dMZsq1)/SW)
         dZZA1sm = dZZA1


** neutralino and chargino contributions
         if (chiRUN) then

** W boson mass renormalisation
            do Neu3 = 1,4
               do Cha3 = 1,2                  
                  eta24 = (Sqrt(2.D0)*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
     -                 UChaC(Cha3,2)*ZNeu(Neu3,3))/Sqrt(2.D0)
                  eta25 = (Sqrt(2.D0)*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -                 VChaC(Cha3,2)*ZNeu(Neu3,4))/Sqrt(2.D0)
                  eta26 = (Sqrt(2.D0)*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
     -                 UCha(Cha3,2)*ZNeuC(Neu3,3))/Sqrt(2.D0)
                  eta27 = (Sqrt(2.D0)*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
     -                 VCha(Cha3,2)*ZNeuC(Neu3,4))/Sqrt(2.D0)
                  
                  dMWsq1 = dMWsq1 - Alfa/(2.D0*Pi*SW2)*((eta24* 
     -                 eta26 + eta25*eta27)*DBLE(A0(MNeu2(Neu3)))) + 
     -                 Alfa/(Pi*SW2)*((eta24*eta26 + eta25*eta27)*
     -                 DBLE(B0i(bb00,MW2,MCha2(Cha3),MNeu2(Neu3)))) - 
     -                 (Alfa*MW2)/(2.D0*Pi*SW2)*
     -                 ((eta24*eta26 + eta25*eta27)*
     -                 DBLE(B0i(bb1,MW2,MCha2(Cha3),MNeu2(Neu3)))) + 
     -                 Alfa/(2.D0*Pi*SW2)*(DBLE(B0i(bb0,MW2,MCha2(Cha3),
     -                 MNeu2(Neu3)))*(-(eta25*eta27*MCha2(Cha3)) + eta26*
     -                 eta27*MCha(Cha3)*MNeu(Neu3) + eta24*(-(eta26*
     -                 MCha2(Cha3)) + eta25*MCha(Cha3)*MNeu(Neu3))))   
               enddo
            enddo


** Z boson mass renormalisation
            do Cha4 = 1,2
               do Cha3 = 1,2               
                  eta28 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2)) - 
     -                 UCha(Cha4,1)*UChaC(Cha3,1)               
                  eta29 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     -                 UCha(Cha3,1)*UChaC(Cha4,1)               
                  eta30 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     -                 VCha(Cha4,1)*VChaC(Cha3,1)               
                  eta31 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2)) - 
     -                 VCha(Cha3,1)*VChaC(Cha4,1)
               
                  dMZsq1 = dMZsq1 - Alfa/(2.D0*CW2*Pi*SW2)*((eta28*
     -                 eta29 + eta30*eta31)*DBLE(A0(MCha2(Cha4)))) + 
     -                 Alfa/(CW2*Pi*SW2)*((eta28*eta29 + eta30*eta31)*
     -                 DBLE(B0i(bb00,MZ2,MCha2(Cha3),MCha2(Cha4)))) - 
     -                 (Alfa*MZ2)/(2.D0*CW2*Pi*SW2)*
     -                 ((eta28*eta29 + eta30*eta31)*
     -                 DBLE(B0i(bb1,MZ2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 Alfa/(2.D0*CW2*Pi*SW2)*
     -                 (DBLE(B0i(bb0,MZ2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (eta29*eta31*MCha(Cha3)*MCha(Cha4) - 
     -                 eta30*eta31*MCha2(Cha3) + eta28*(eta30*MCha(Cha3)*
     -                 MCha(Cha4) - eta29*MCha2(Cha3))))               
               enddo
            enddo
         
            do Neu4 = 1,4
               do Neu3 = 1,4               
                  eta32 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     -                 ZNeu(Neu4,4)*ZNeuC(Neu3,4)               
                  eta33 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     -                 ZNeu(Neu3,4)*ZNeuC(Neu4,4)
                  
                  dMZsq1 = dMZsq1 - Alfa/(8.D0*CW2*Pi*SW2)*
     -                 (eta32*eta33*DBLE(A0(MNeu2(Neu4)))) + 
     -                 Alfa/(4.D0*CW2*Pi*SW2)*(eta32*eta33*
     -                 DBLE(B0i(bb00,MZ2,MNeu2(Neu3),MNeu2(Neu4)))) -
     -                 (Alfa*MZ2)/(8.D0*CW2*Pi*SW2)*(eta32*eta33*
     -                 DBLE(B0i(bb1,MZ2,MNeu2(Neu3),MNeu2(Neu4))))+
     -                 DBLE(B0i(bb0,MZ2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (Alfa/(16.D0*CW2*Pi*SW2)*
     -                 ((-eta32**2 - eta33**2)*MNeu(Neu3)*MNeu(Neu4)) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(eta32*eta33*MNeu2(Neu3)))
               enddo
            enddo            


** Charge renormalisation
*            if (charge .ne. 2) then  !  LO via alfa(0), alfa(MZ)
               do Cha3 = 1,2
                  dZAA1chi = dZAA1chi + Alfa/Pi*DBLE(B0i(bb1,0.D0,
     -                 MCha2(Cha3),MCha2(Cha3))) - (2*Alfa)/Pi*
     -                 DBLE(B0i(dbb00,0.D0,MCha2(Cha3),MCha2(Cha3)))
               enddo
               dZAA1 = dZAA1 + dZAA1chi
*            else  ! LO via GF
c$$$               do Neu3 = 1,4
c$$$                  do Cha3 = 1,2
c$$$                     dup35 = (Sqrt(2.D0)*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
c$$$     -                    UChaC(Cha3,2)*ZNeu(Neu3,3))/Sqrt(2.D0)
c$$$                     dup36 = (Sqrt(2.D0)*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
c$$$     -                    VChaC(Cha3,2)*ZNeu(Neu3,4))/Sqrt(2.D0)
c$$$                     dup37 = (Sqrt(2.D0)*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
c$$$     -                    UCha(Cha3,2)*ZNeuC(Neu3,3))/Sqrt(2.D0)
c$$$                     dup38 = (Sqrt(2.D0)*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
c$$$     -                    VCha(Cha3,2)*ZNeuC(Neu3,4))/Sqrt(2.D0)
c$$$
c$$$                     dZeR = dZeR - 1/(2.D0*MW2)*
c$$$     -                    (-(Alfa/(2.D0*Pi*SW2)*
c$$$     -                    ((dup35*dup37 + dup36*dup38)*A0(MCha2(Cha3)))) + 
c$$$     -                    Alfa/(Pi*SW2)*((dup35*dup37 + dup36*dup38)*
c$$$     -                    B0i(bb00,0.D0,MCha2(Cha3),MNeu2(Neu3))) + 
c$$$     -                    Alfa/(2.D0*Pi*SW2)*
c$$$     -                    (B0i(bb0,0.D0,MCha2(Cha3),MNeu2(Neu3))*
c$$$     -                    (dup37*dup38*MCha(Cha3)*MNeu(Neu3) - 
c$$$     -                    dup36*dup38*MNeu2(Neu3) + 
c$$$     -                    dup35*(dup36*MCha(Cha3)*MNeu(Neu3) - 
c$$$     -                    dup37*MNeu2(Neu3)))))
c$$$                  enddo
c$$$               enddo
*            end if
            
            do Cha3 = 1,2
               dup40 = 4*SW2 - 2*UCha(Cha3,1)*UChaC(Cha3,1) - 
     -              UCha(Cha3,2)*UChaC(Cha3,2) - 2*VCha(Cha3,1)*
     -              VChaC(Cha3,1) - VCha(Cha3,2)*VChaC(Cha3,2)
               
               dZZA1chi = dZZA1chi + 2/MZ2*
     -              (Alfa/(4.D0*CW*Pi)*(dup40*DBLE(A0(MCha2(Cha3))))/SW - 
     -              Alfa/(2.D0*CW*Pi)*(dup40*DBLE(B0i(bb00,0.D0,
     -              MCha2(Cha3),MCha2(Cha3))))/SW)
            enddo
            dZZA1 = dZZA1 + dZZA1chi


** Higgs field renormalisation
            do Cha4 = 1,2
               do Cha3 = 1,2
               
                  dZh0= dZh0 - Alfa/(8.D0*Pi*SW2)*(DBLE(Divergence)*
     -                 ((SA*UCha(Cha4,2)*VCha(Cha3,1) - CA*UCha(Cha4,1)*
     -                 VCha(Cha3,2))*(SA*UChaC(Cha4,2)*VChaC(Cha3,1) - 
     -                 CA*UChaC(Cha4,1)*VChaC(Cha3,2)) + (SA*
     -                 UCha(Cha3,2)*VCha(Cha4,1) - CA*UCha(Cha3,1)*
     -                 VCha(Cha4,2))*(SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2))))
               
                  dZHH = dZHH - Alfa/(8.D0*Pi*SW2)*(DBLE(Divergence)*
     -                 ((CA*UCha(Cha4,2)*VCha(Cha3,1) + SA*UCha(Cha4,1)*
     -                 VCha(Cha3,2))*(CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 SA*UChaC(Cha4,1)*VChaC(Cha3,2)) + (CA*
     -                 UCha(Cha3,2)*VCha(Cha4,1) + SA*UCha(Cha3,1)*
     -                 VCha(Cha4,2))*(CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -                 SA*UChaC(Cha3,1)*VChaC(Cha4,2))))                  
               enddo
            enddo
            
            do Neu4 = 1,4
               do Neu3 = 1,4               
                  eta35 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)               
                  eta36 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  
                  dZh0 = dZh0 - Alfa/(16.D0*CW2*Pi*SW2)*(DBLE(Divergence)*
     -                 (eta35*SA*ZNeu(Neu3,3) + CA*eta35*ZNeu(Neu3,4) + 
     -                 (-(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1))*
     -                 (SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)))*
     -                 (eta36*SA*ZNeuC(Neu3,3) +CA*eta36*ZNeuC(Neu3,4) + 
     -                 (-(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1))*
     -                 (SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4))))
                  
                  dZHH = dZHH - Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (DBLE(Divergence)*(CA*ZNeu(Neu3,3)*
     -                 (-(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)) + 
     -                 ZNeu(Neu3,4)*
     -                 (CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)) + 
     -                 (-(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1))*
     -                 (CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)))*
     -                 (CA*ZNeuC(Neu3,3)*(-(CW*ZNeuC(Neu4,2)) + 
     -                 SW*ZNeuC(Neu4,1)) + ZNeuC(Neu3,4)*
     -                 (CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)) + 
     -                 (-(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1))*
     -                 (CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4))))
               enddo
            enddo


** W boson field renormalisation
            do Neu3 = 1,4
               do Cha3 = 1,2               
                  eta37 = (Sqrt(2.D0)*UChaC(Cha3,1)*ZNeu(Neu3,2) + 
     -                 UChaC(Cha3,2)*ZNeu(Neu3,3))/Sqrt(2.D0)              
                  eta38 = (Sqrt(2.D0)*VChaC(Cha3,1)*ZNeu(Neu3,2) - 
     -                 VChaC(Cha3,2)*ZNeu(Neu3,4))/Sqrt(2.D0)               
                  eta39 = (Sqrt(2.D0)*UCha(Cha3,1)*ZNeuC(Neu3,2) + 
     -                 UCha(Cha3,2)*ZNeuC(Neu3,3))/Sqrt(2.D0)               
                  eta40 = (Sqrt(2.D0)*VCha(Cha3,1)*ZNeuC(Neu3,2) - 
     -                 VCha(Cha3,2)*ZNeuC(Neu3,4))/Sqrt(2.D0)
                  
                  dZW1 = dZW1 + Alfa/(2.D0*Pi*SW2)*
     -                 ((eta37*eta39 + eta38*eta40)*
     -                 DBLE(B0i(bb1,MW2,MCha2(Cha3),MNeu2(Neu3)))) - 
     -                 Alfa/(Pi*SW2)*((eta37*eta39 + eta38*eta40)*
     -                 DBLE(B0i(dbb00,MW2,MCha2(Cha3),MNeu2(Neu3)))) + 
     -                 (Alfa*MW2)/(2.D0*Pi*SW2)*
     -                 ((eta37*eta39 + eta38*eta40)*
     -                 DBLE(B0i(dbb1,MW2,MCha2(Cha3),MNeu2(Neu3)))) - 
     -                 Alfa/(2.D0*Pi*SW2)*(DBLE(B0i(dbb0,MW2,MCha2(Cha3),
     -                 MNeu2(Neu3)))*(-(eta38*eta40*MCha2(Cha3)) + 
     -                 eta39*eta40*MCha(Cha3)*MNeu(Neu3) + 
     -                 eta37*(-(eta39*MCha2(Cha3)) + 
     -                 eta38*MCha(Cha3)*MNeu(Neu3))))               
               enddo
            enddo


** Z boson field renormalisation
            do Cha4 = 1,2
               do Cha3 = 1,2               
                  eta41 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(UCha(Cha4,2)*UChaC(Cha3,2)) - 
     -                 UCha(Cha4,1)*UChaC(Cha3,1)               
                  eta42 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(UCha(Cha3,2)*UChaC(Cha4,2)) - 
     -                 UCha(Cha3,1)*UChaC(Cha4,1)               
                  eta43 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(VCha(Cha4,2)*VChaC(Cha3,2)) - 
     -                 VCha(Cha4,1)*VChaC(Cha3,1)               
                  eta44 = SW2*IndexDelta(Cha3,Cha4) - 
     -                 1/2.D0*(VCha(Cha3,2)*VChaC(Cha4,2)) - 
     -                 VCha(Cha3,1)*VChaC(Cha4,1)
                  
                  dZZZ1 = dZZZ1 + Alfa/(2.D0*CW2*Pi*SW2)*
     -                 ((eta41*eta42 + eta43*eta44)*
     -                 DBLE(B0i(bb1,MZ2,MCha2(Cha3),MCha2(Cha4)))) - 
     -                 Alfa/(CW2*Pi*SW2)*((eta41*eta42 + eta43*eta44)*
     -                 DBLE(B0i(dbb00,MZ2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 (Alfa*MZ2)/(2.D0*CW2*Pi*SW2)*
     -                 ((eta41*eta42 + eta43*eta44)*
     -                 DBLE(B0i(dbb1,MZ2,MCha2(Cha3),MCha2(Cha4)))) - 
     -                 Alfa/(2.D0*CW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MZ2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (eta42*eta44*MCha(Cha3)*MCha(Cha4) - eta43*eta44*
     -                 MCha2(Cha3) + eta41*(eta43*MCha(Cha3)*MCha(Cha4) - 
     -                 eta42*MCha2(Cha3))))                  
               enddo
            enddo
         
            do Neu4 = 1,4
               do Neu3 = 1,4
                  
                  eta45 = ZNeu(Neu4,3)*ZNeuC(Neu3,3) - 
     -                 ZNeu(Neu4,4)*ZNeuC(Neu3,4)               
                  eta46 = ZNeu(Neu3,3)*ZNeuC(Neu4,3) - 
     -                 ZNeu(Neu3,4)*ZNeuC(Neu4,4)
                  
                  dZZZ1 = dZZZ1 + Alfa/(8.D0*CW2*Pi*SW2)*(eta45*eta46*
     -                 DBLE(B0i(bb1,MZ2,MNeu2(Neu3),MNeu2(Neu4))))-
     -                 Alfa/(4.D0*CW2*Pi*SW2)*(eta45*eta46*
     -                 DBLE(B0i(dbb00,MZ2,MNeu2(Neu3),MNeu2(Neu4)))) +
     -                 (Alfa*MZ2)/(8.D0*CW2*Pi*SW2)*(eta45*eta46*
     -                 DBLE(B0i(dbb1,MZ2,MNeu2(Neu3),MNeu2(Neu4)))) -
     -                 DBLE(B0i(dbb0,MZ2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (Alfa/(16.D0*CW2*Pi*SW2)*
     -                 ((-eta45**2 - eta46**2)*MNeu(Neu3)*MNeu(Neu4)) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(eta45*eta46*MNeu2(Neu3))) 
               enddo
            enddo


** photon-Z field renormalisation
            do Cha3 = 1,2
               dup15 = 4*SW2 - 2*UCha(Cha3,1)*UChaC(Cha3,1) - 
     -              UCha(Cha3,2)*UChaC(Cha3,2) - 2*VCha(Cha3,1)*
     -              VChaC(Cha3,1) - VCha(Cha3,2)*VChaC(Cha3,2)

               dZAZ1 = dZAZ1 - 2/MZ2*(Alfa/(4.D0*CW*Pi)*(dup15*
     -              DBLE(A0(MCha2(Cha3))))/SW- Alfa/(2.D0*CW*Pi)*(dup15*
     -              DBLE(B0i(bb00,MZ2,MCha2(Cha3),MCha2(Cha3))))/SW + 
     -              (Alfa*MZ2)/(4.D0*CW*Pi)*(dup15*
     -              DBLE(B0i(bb1,MZ2,MCha2(Cha3),MCha2(Cha3))))/SW)
            enddo


** fermion field renormalisation
            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfLu = dZfLu + Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MU2,MCha2(Cha3),MSf2(Sfe3,4,1)))*
     -                 (2*MW*(CB*UChaC(Cha3,1)*USf(Sfe3,1,4,1)) - 
     -                 MD*(Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*
     -                 (2*MW*(CB*UCha(Cha3,1)*USfC(Sfe3,1,4,1)) - 
     -                 MD*(Sqrt(2.D0)*UCha(Cha3,2)*USfC(Sfe3,2,4,1))))/
     -                 CB2 - (Alfa*MU2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MU2,MCha2(Cha3),MSf2(Sfe3,4,1)))*
     -                 (-(MU2*(CB2*USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)*
     -                 VCha(Cha3,2)*VChaC(Cha3,2))) + SB2*(-(2*MW2*(CB2*
     -                 UChaC(Cha3,1)*USf(Sfe3,1,4,1))) + MD*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*
     -                 UCha(Cha3,1)*USfC(Sfe3,1,4,1) + SB2*(MD*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,1)*USf(Sfe3,1,4,1)) - MD2*
     -                 (UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*UCha(Cha3,2)*
     -                 USfC(Sfe3,2,4,1)))/(CB2*SB2) + (Alfa*MU2)/(8.D0*
     -                 MW2*Pi*SW2)*(DBLE(B0i(dbb0,MU2,MCha2(Cha3),
     -                 MSf2(Sfe3,4,1)))*MCha(Cha3)*(MW*(CB*Sqrt(2.D0)*
     -                 UCha(Cha3,1)*USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)*
     -                 VCha(Cha3,2)) - MD*(UCha(Cha3,2)*USf(Sfe3,1,4,1)*
     -                 USfC(Sfe3,2,4,1)*VCha(Cha3,2)) + (MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,1)*USf(Sfe3,1,4,1)) - MD*
     -                 (UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*
     -                 USfC(Sfe3,1,4,1)*VChaC(Cha3,2)))/(CB*SB)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup26 = 3*CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)
                  dup27 = 4*MW*(SB*SW*USf(Sfe3,2,3,1)*ZNeu(Neu3,1)) - 
     -                 3*CW*MsU*(USf(Sfe3,1,3,1)*ZNeu(Neu3,4))
                  dup28 = 3*CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)
                  dup29 = 4*MW*(SB*SW*USfC(Sfe3,2,3,1)*ZNeuC(Neu3,1)) - 
     -                 3*CW*MsU*(USfC(Sfe3,1,3,1)*ZNeuC(Neu3,4))

                  dZfLu = dZfLu + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MU2,MNeu2(Neu3),MSf2(Sfe3,3,1)))*
     -                 (MW*(dup26*SB*USfC(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USfC(Sfe3,2,3,1)*ZNeu(Neu3,4)))*
     -                 (MW*(dup28*SB*USf(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USf(Sfe3,2,3,1)*ZNeuC(Neu3,4))))/SB2 - 
     -                 (Alfa*MU2)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MU2,MNeu2(Neu3),MSf2(Sfe3,3,1)))*
     -                 (-(dup27*dup29)-(MW*(dup26*SB*USfC(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USfC(Sfe3,2,3,1)*ZNeu(Neu3,4)))*
     -                 (MW*(dup28*SB*USf(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USf(Sfe3,2,3,1)*ZNeuC(Neu3,4)))))/SB2 + 
     -                 (Alfa*MsU)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MU2,MNeu2(Neu3),MSf2(Sfe3,3,1)))*
     -                 (dup27*(MW*(dup26*SB*USfC(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USfC(Sfe3,2,3,1)*ZNeu(Neu3,4))) + 
     -                 dup29*(MW*(dup28*SB*USf(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USf(Sfe3,2,3,1)*ZNeuC(Neu3,4))))*
     -                 MNeu(Neu3))/SB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfRu = dZfRu - 
     -                 (Alfa*MU2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MU2,MCha2(Cha3),MSf2(Sfe3,4,1)))*
     -                 (-(MU2*(CB2*USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)*
     -                 VCha(Cha3,2)*VChaC(Cha3,2))) + SB2*(-(2*MW2*(CB2*
     -                 UChaC(Cha3,1)*USf(Sfe3,1,4,1))) + MD*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*
     -                 UCha(Cha3,1)*USfC(Sfe3,1,4,1) + SB2*(MD*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,1)*USf(Sfe3,1,4,1)) - 
     -                 MD2*(UChaC(Cha3,2)*USf(Sfe3,2,4,1)))*
     -                 UCha(Cha3,2)*USfC(Sfe3,2,4,1)))/(CB2*SB2) + 
     -                 (Alfa*MU2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MU2,MCha2(Cha3),MSf2(Sfe3,4,1)))*
     -                 USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)*VCha(Cha3,2)*
     -                 VChaC(Cha3,2))/SB2 + 
     -                 (Alfa*MU2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MU2,MCha2(Cha3),MSf2(Sfe3,4,1)))*
     -                 MCha(Cha3)*(MW*
     -                 (CB*Sqrt(2.D0)*UCha(Cha3,1)*USf(Sfe3,1,4,1)*
     -                 USfC(Sfe3,1,4,1)*VCha(Cha3,2)) -MD*(UCha(Cha3,2)*
     -                 USf(Sfe3,1,4,1)*USfC(Sfe3,2,4,1)*VCha(Cha3,2)) + 
     -                 (MW*(CB*Sqrt(2.D0)*UChaC(Cha3,1)*
     -                 USf(Sfe3,1,4,1)) - MD*(UChaC(Cha3,2)*
     -                 USf(Sfe3,2,4,1)))*USfC(Sfe3,1,4,1)*
     -                 VChaC(Cha3,2)))/(CB*SB)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup30 = 3*CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)
                  dup31 = 4*MW*(SB*SW*USf(Sfe3,2,3,1)*ZNeu(Neu3,1)) - 
     -                 3*CW*MsU*(USf(Sfe3,1,3,1)*ZNeu(Neu3,4))
                  dup32 = 3*CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)
                  dup33 = 4*MW*(SB*SW*USfC(Sfe3,2,3,1)*ZNeuC(Neu3,1)) - 
     -                 3*CW*MsU*(USfC(Sfe3,1,3,1)*ZNeuC(Neu3,4))

                  dZfRu = dZfRu + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (dup31*dup33*DBLE(B0i(bb1,MU2,MNeu2(Neu3),
     -                 MSf2(Sfe3,3,1))))/SB2 - 
     -                 (Alfa*MU2)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MU2,MNeu2(Neu3),MSf2(Sfe3,3,1)))*
     -                 (-(dup31*dup33) - 
     -                 (MW*(dup30*SB*USfC(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USfC(Sfe3,2,3,1)*ZNeu(Neu3,4)))*
     -                 (MW*(dup32*SB*USf(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USf(Sfe3,2,3,1)*ZNeuC(Neu3,4)))))/SB2 + 
     -                 (Alfa*MsU)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MU2,MNeu2(Neu3),MSf2(Sfe3,3,1)))*
     -                 (dup31*(MW*(dup30*SB*USfC(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USfC(Sfe3,2,3,1)*ZNeu(Neu3,4))) + 
     -                 dup33*(MW*(dup32*SB*USf(Sfe3,1,3,1)) + 
     -                 3*CW*MsU*(USf(Sfe3,2,3,1)*ZNeuC(Neu3,4))))*
     -                 MNeu(Neu3))/SB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfLd = dZfLd + Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 (2*MW*(SB*USfC(Sfe3,1,3,1)*VCha(Cha3,1)) - 
     -                 MsU*(Sqrt(2.D0)*USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*
     -                 (2*MW*(SB*USf(Sfe3,1,3,1)*VChaC(Cha3,1)) - 
     -                 MsU*(Sqrt(2.D0)*USf(Sfe3,2,3,1)*VChaC(Cha3,2))))/
     -                 SB2 + (Alfa*MD2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 MCha(Cha3)*((MW*(SB*Sqrt(2.D0)*USfC(Sfe3,1,3,1)*
     -                 VCha(Cha3,1)) - MsU*(USfC(Sfe3,2,3,1)*
     -                 VCha(Cha3,2)))*UCha(Cha3,2)*USf(Sfe3,1,3,1) + 
     -                 (MW*(SB*Sqrt(2.D0)*USf(Sfe3,1,3,1)*VChaC(Cha3,1)) - 
     -                 MsU*(USf(Sfe3,2,3,1)*VChaC(Cha3,2)))*
     -                 UChaC(Cha3,2)*USfC(Sfe3,1,3,1)))/(CB*SB) - 
     -                 (Alfa*MD2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 (-(MD2*(SB2*UCha(Cha3,2)*UChaC(Cha3,2)*
     -                 USf(Sfe3,1,3,1)*USfC(Sfe3,1,3,1))) + CB2*((-(2*
     -                 MW2*(SB2*USfC(Sfe3,1,3,1)*VCha(Cha3,1))) + 
     -                 MsU*MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*
     -                 USf(Sfe3,1,3,1)*VChaC(Cha3,1) + 
     -                 (MsU*MW*(SB*Sqrt(2.D0)*USfC(Sfe3,1,3,1)*
     -                 VCha(Cha3,1)) - 
     -                 MU2*(USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*
     -                 USf(Sfe3,2,3,1)*VChaC(Cha3,2))))/(CB2*SB2)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup26 = 3*CW*ZNeu(Neu3,2) - SW*ZNeu(Neu3,1)
                  dup27 = 3*CW*ZNeuC(Neu3,2) - SW*ZNeuC(Neu3,1)

                  dZfLd = dZfLd + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MD2,MNeu2(Neu3),MSf2(Sfe3,4,1)))*
     -                 (MW*(CB*dup26*USfC(Sfe3,1,4,1)) - 
     -                 3*CW*MD*(USfC(Sfe3,2,4,1)*ZNeu(Neu3,3)))*
     -                 (MW*(CB*dup27*USf(Sfe3,1,4,1)) - 
     -                 3*CW*MD*(USf(Sfe3,2,4,1)*ZNeuC(Neu3,3))))/CB2 + 
     -                 (Alfa*MD)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MD2,MNeu2(Neu3),MSf2(Sfe3,4,1)))*
     -                 ((2*MW*(CB*SW*USf(Sfe3,2,4,1)*ZNeu(Neu3,1)) + 
     -                 3*CW*MD*(USf(Sfe3,1,4,1)*ZNeu(Neu3,3)))*
     -                 (MW*(CB*dup26*USfC(Sfe3,1,4,1)) - 
     -                 3*CW*MD*(USfC(Sfe3,2,4,1)*ZNeu(Neu3,3))) + 
     -                 (MW*(CB*dup27*USf(Sfe3,1,4,1)) - 
     -                 3*CW*MD*(USf(Sfe3,2,4,1)*ZNeuC(Neu3,3)))*
     -                 (2*MW*(CB*SW*USfC(Sfe3,2,4,1)*ZNeuC(Neu3,1)) + 
     -                 3*CW*MD*(USfC(Sfe3,1,4,1)*ZNeuC(Neu3,3))))*
     -                 MNeu(Neu3))/CB2 + 
     -                 (Alfa*MD)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MD2,MNeu2(Neu3),MSf2(Sfe3,4,1)))*
     -                 ((MD*((MW2*(CB2*(SW2*ZNeuC(Neu3,1) - 
     -                 3*CW*(SW*ZNeuC(Neu3,2)))*ZNeu(Neu3,1))-3*MW2*
     -                 (CB2*(CW*(SW*ZNeuC(Neu3,1)) - 
     -                 3*CW2*ZNeuC(Neu3,2))*ZNeu(Neu3,2)) + 
     -                 9*CW2*MD2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,1,4,1)) + 9*MD2*MW*
     -                 (CB*(CW*(SW*ZNeuC(Neu3,1)) - CW2*ZNeuC(Neu3,2))*
     -                 USfC(Sfe3,2,4,1)*ZNeu(Neu3,3)))*
     -                 USf(Sfe3,1,4,1) + 
     -                 (MD*((4*MW2*SW2*(CB2*ZNeu(Neu3,1)*ZNeuC(Neu3,1)) + 
     -                 9*CW2*MD2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,2,4,1)) + 9*MD2*MW*
     -                 (CB*(CW*(SW*ZNeu(Neu3,1)) - CW2*ZNeu(Neu3,2))*
     -                 USfC(Sfe3,1,4,1)*ZNeuC(Neu3,3)))*
     -                 USf(Sfe3,2,4,1)))/CB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfRd = dZfRd + (Alfa*MD2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 UCha(Cha3,2)*UChaC(Cha3,2)*USf(Sfe3,1,3,1)*
     -                 USfC(Sfe3,1,3,1))/CB2 + 
     -                 (Alfa*MD2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 MCha(Cha3)*((MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,1,3,1)*VCha(Cha3,1)) - 
     -                 MsU*(USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*UCha(Cha3,2)*
     -                 USf(Sfe3,1,3,1) + 
     -                 (MW*(SB*Sqrt(2.D0)*USf(Sfe3,1,3,1)*VChaC(Cha3,1)) - 
     -                 MsU*(USf(Sfe3,2,3,1)*VChaC(Cha3,2)))*
     -                 UChaC(Cha3,2)*USfC(Sfe3,1,3,1)))/(CB*SB) - 
     -                 (Alfa*MD2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MD2,MCha2(Cha3),MSf2(Sfe3,3,1)))*
     -                 (-(MD2*(SB2*UCha(Cha3,2)*UChaC(Cha3,2)*
     -                 USf(Sfe3,1,3,1)*USfC(Sfe3,1,3,1))) + CB2*((-(2*
     -                 MW2*(SB2*USfC(Sfe3,1,3,1)*VCha(Cha3,1))) + 
     -                 MsU*MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*
     -                 USf(Sfe3,1,3,1)*VChaC(Cha3,1) + 
     -                 (MsU*MW*(SB*Sqrt(2.D0)*USfC(Sfe3,1,3,1)*
     -                 VCha(Cha3,1)) - 
     -                 MU2*(USfC(Sfe3,2,3,1)*VCha(Cha3,2)))*
     -                 USf(Sfe3,2,3,1)*VChaC(Cha3,2))))/(CB2*SB2)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup28 = 2*MW*(CB*SW*USf(Sfe3,2,4,1)*ZNeu(Neu3,1)) + 
     -                 3*CW*MD*(USf(Sfe3,1,4,1)*ZNeu(Neu3,3))
                  dup29 = 2*MW*(CB*SW*USfC(Sfe3,2,4,1)*ZNeuC(Neu3,1)) + 
     -                 3*CW*MD*(USfC(Sfe3,1,4,1)*ZNeuC(Neu3,3))

                  dZfRd = dZfRd + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (dup28*dup29*DBLE(B0i(bb1,MD2,MNeu2(Neu3),
     -                 MSf2(Sfe3,4,1))))/CB2 + 
     -                 (Alfa*MD)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MD2,MNeu2(Neu3),MSf2(Sfe3,4,1)))*
     -                 (dup28*(MW*(CB*USfC(Sfe3,1,4,1)*
     -                 (3*CW*ZNeu(Neu3,2) - SW*ZNeu(Neu3,1))) - 
     -                 3*CW*MD*(USfC(Sfe3,2,4,1)*ZNeu(Neu3,3))) + 
     -                 dup29*(MW*(CB*USf(Sfe3,1,4,1)*
     -                 (3*CW*ZNeuC(Neu3,2) - SW*ZNeuC(Neu3,1))) - 
     -                 3*CW*MD*(USf(Sfe3,2,4,1)*ZNeuC(Neu3,3))))*
     -                 MNeu(Neu3))/CB2 + 
     -                 (Alfa*MD)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MD2,MNeu2(Neu3),MSf2(Sfe3,4,1)))*
     -                 ((MD*((MW2*(CB2*(SW2*ZNeuC(Neu3,1) - 
     -                 3*CW*(SW*ZNeuC(Neu3,2)))*ZNeu(Neu3,1))-3*MW2*
     -                 (CB2*(CW*(SW*ZNeuC(Neu3,1)) - 
     -                 3*CW2*ZNeuC(Neu3,2))*ZNeu(Neu3,2)) + 
     -                 9*CW2*MD2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,1,4,1)) + 9*MD2*MW*
     -                 (CB*(CW*(SW*ZNeuC(Neu3,1)) - CW2*ZNeuC(Neu3,2))*
     -                 USfC(Sfe3,2,4,1)*ZNeu(Neu3,3)))*
     -                 USf(Sfe3,1,4,1) + 
     -                 (MD*((4*MW2*SW2*(CB2*ZNeu(Neu3,1)*ZNeuC(Neu3,1)) + 
     -                 9*CW2*MD2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,2,4,1)) + 9*MD2*MW*
     -                 (CB*(CW*(SW*ZNeu(Neu3,1)) - CW2*ZNeu(Neu3,2))*
     -                 USfC(Sfe3,1,4,1)*ZNeuC(Neu3,3)))*
     -                 USf(Sfe3,2,4,1)))/CB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfLc = dZfLc + Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 (2*MW*(CB*UChaC(Cha3,1)*USf(Sfe3,1,4,2)) - 
     -                 MS*(Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 (2*MW*(CB*UCha(Cha3,1)*USfC(Sfe3,1,4,2)) - 
     -                 MS*(Sqrt(2.D0)*UCha(Cha3,2)*USfC(Sfe3,2,4,2))))/
     -                 CB2 - (Alfa*MC2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 (-(MC2*(CB2*USf(Sfe3,1,4,2)*USfC(Sfe3,1,4,2)*
     -                 VCha(Cha3,2)*VChaC(Cha3,2))) + SB2*(-(2*MW2*(CB2*
     -                 UChaC(Cha3,1)*USf(Sfe3,1,4,2))) + MS*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 UCha(Cha3,1)*USfC(Sfe3,1,4,2) + 
     -                 SB2*(MS*MW*(CB*Sqrt(2.D0)*UChaC(Cha3,1)*
     -                 USf(Sfe3,1,4,2)) - 
     -                 MS2*(UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 UCha(Cha3,2)*USfC(Sfe3,2,4,2)))/(CB2*SB2) + 
     -                 (Alfa*MC2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 MCha(Cha3)*(MW*
     -                 (CB*Sqrt(2.D0)*UCha(Cha3,1)*USf(Sfe3,1,4,2)*
     -                 USfC(Sfe3,1,4,2)*VCha(Cha3,2)) -MS*(UCha(Cha3,2)*
     -                 USf(Sfe3,1,4,2)*USfC(Sfe3,2,4,2)*VCha(Cha3,2)) + 
     -                 (MW*(CB*Sqrt(2.D0)*UChaC(Cha3,1)*USf(Sfe3,1,4,2)) - 
     -                 MS*(UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 USfC(Sfe3,1,4,2)*VChaC(Cha3,2)))/(CB*SB)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup26 = 3*CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)
                  dup27 = 4*MW*(SB*SW*USf(Sfe3,2,3,2)*ZNeu(Neu3,1)) - 
     -                 3*CW*MC*(USf(Sfe3,1,3,2)*ZNeu(Neu3,4))
                  dup28 = 3*CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)
                  dup29 = 4*MW*(SB*SW*USfC(Sfe3,2,3,2)*ZNeuC(Neu3,1)) - 
     -                 3*CW*MC*(USfC(Sfe3,1,3,2)*ZNeuC(Neu3,4))

                  dZfLc = dZfLc + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MC2,MNeu2(Neu3),MSf2(Sfe3,3,2)))*
     -                 (MW*(dup26*SB*USfC(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USfC(Sfe3,2,3,2)*ZNeu(Neu3,4)))*
     -                 (MW*(dup28*SB*USf(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USf(Sfe3,2,3,2)*ZNeuC(Neu3,4))))/SB2 - 
     -                 (Alfa*MC2)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MC2,MNeu2(Neu3),MSf2(Sfe3,3,2)))*
     -                 (-(dup27*dup29) - 
     -                 (MW*(dup26*SB*USfC(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USfC(Sfe3,2,3,2)*ZNeu(Neu3,4)))*
     -                 (MW*(dup28*SB*USf(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USf(Sfe3,2,3,2)*ZNeuC(Neu3,4)))))/SB2 + 
     -                 (Alfa*MC)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MC2,MNeu2(Neu3),MSf2(Sfe3,3,2)))*
     -                 (dup27*(MW*(dup26*SB*USfC(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USfC(Sfe3,2,3,2)*ZNeu(Neu3,4))) + 
     -                 dup29*(MW*(dup28*SB*USf(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USf(Sfe3,2,3,2)*ZNeuC(Neu3,4))))*
     -                 MNeu(Neu3))/SB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfRc = dZfRc - (Alfa*MC2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 (-(MC2*(CB2*USf(Sfe3,1,4,2)*USfC(Sfe3,1,4,2)*
     -                 VCha(Cha3,2)*VChaC(Cha3,2))) + SB2*(-(2*MW2*(CB2*
     -                 UChaC(Cha3,1)*USf(Sfe3,1,4,2))) + MS*MW*(CB*
     -                 Sqrt(2.D0)*UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 UCha(Cha3,1)*USfC(Sfe3,1,4,2) + 
     -                 SB2*(MS*MW*(CB*Sqrt(2.D0)*UChaC(Cha3,1)*
     -                 USf(Sfe3,1,4,2)) - 
     -                 MS2*(UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 UCha(Cha3,2)*USfC(Sfe3,2,4,2)))/(CB2*SB2) + 
     -                 (Alfa*MC2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 USf(Sfe3,1,4,2)*USfC(Sfe3,1,4,2)*VCha(Cha3,2)*
     -                 VChaC(Cha3,2))/SB2 + 
     -                 (Alfa*MC2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MC2,MCha2(Cha3),MSf2(Sfe3,4,2)))*
     -                 MCha(Cha3)*(MW*
     -                 (CB*Sqrt(2.D0)*UCha(Cha3,1)*USf(Sfe3,1,4,2)*
     -                 USfC(Sfe3,1,4,2)*VCha(Cha3,2)) -MS*(UCha(Cha3,2)*
     -                 USf(Sfe3,1,4,2)*USfC(Sfe3,2,4,2)*VCha(Cha3,2)) + 
     -                 (MW*(CB*Sqrt(2.D0)*UChaC(Cha3,1)*USf(Sfe3,1,4,2)) - 
     -                 MS*(UChaC(Cha3,2)*USf(Sfe3,2,4,2)))*
     -                 USfC(Sfe3,1,4,2)*VChaC(Cha3,2)))/(CB*SB)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup30 = 3*CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)
                  dup31 = 4*MW*(SB*SW*USf(Sfe3,2,3,2)*ZNeu(Neu3,1)) - 
     -                 3*CW*MC*(USf(Sfe3,1,3,2)*ZNeu(Neu3,4))
                  dup32 = 3*CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)
                  dup33 = 4*MW*(SB*SW*USfC(Sfe3,2,3,2)*ZNeuC(Neu3,1)) - 
     -                 3*CW*MC*(USfC(Sfe3,1,3,2)*ZNeuC(Neu3,4))

                  dZfRc = dZfRc + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (dup31*dup33*DBLE(B0i(bb1,MC2,MNeu2(Neu3),
     -                 MSf2(Sfe3,3,2))))/SB2 - 
     -                 (Alfa*MC2)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MC2,MNeu2(Neu3),MSf2(Sfe3,3,2)))*
     -                 (-(dup31*dup33) - 
     -                 (MW*(dup30*SB*USfC(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USfC(Sfe3,2,3,2)*ZNeu(Neu3,4)))*
     -                 (MW*(dup32*SB*USf(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USf(Sfe3,2,3,2)*ZNeuC(Neu3,4)))))/SB2 + 
     -                 (Alfa*MC)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MC2,MNeu2(Neu3),MSf2(Sfe3,3,2)))*
     -                 (dup31*(MW*(dup30*SB*USfC(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USfC(Sfe3,2,3,2)*ZNeu(Neu3,4))) + 
     -                 dup33*(MW*(dup32*SB*USf(Sfe3,1,3,2)) + 
     -                 3*CW*MC*(USf(Sfe3,2,3,2)*ZNeuC(Neu3,4))))*
     -                 MNeu(Neu3))/SB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfLs = dZfLs + Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 (2*MW*(SB*USfC(Sfe3,1,3,2)*VCha(Cha3,1)) - 
     -                 MC*(Sqrt(2.D0)*USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*
     -                 (2*MW*(SB*USf(Sfe3,1,3,2)*VChaC(Cha3,1)) - MC*
     -                 (Sqrt(2.D0)*USf(Sfe3,2,3,2)*VChaC(Cha3,2))))/SB2 + 
     -                 (Alfa*MS2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 MCha(Cha3)*((MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,1,3,2)*VCha(Cha3,1)) - 
     -                 MC*(USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*UCha(Cha3,2)*
     -                 USf(Sfe3,1,3,2) + (MW*(SB*Sqrt(2.D0)*
     -                 USf(Sfe3,1,3,2)*VChaC(Cha3,1)) - 
     -                 MC*(USf(Sfe3,2,3,2)*VChaC(Cha3,2)))*
     -                 UChaC(Cha3,2)*USfC(Sfe3,1,3,2)))/(CB*SB) - 
     -                 (Alfa*MS2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 (-(MS2*(SB2*UCha(Cha3,2)*UChaC(Cha3,2)*
     -                 USf(Sfe3,1,3,2)*USfC(Sfe3,1,3,2)))+CB2*((-(2*MW2*
     -                 (SB2*USfC(Sfe3,1,3,2)*VCha(Cha3,1))) + MC*MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*
     -                 USf(Sfe3,1,3,2)*VChaC(Cha3,1) + 
     -                 (MC*MW*(SB*Sqrt(2.D0)*USfC(Sfe3,1,3,2)*
     -                 VCha(Cha3,1)) - 
     -                 MC2*(USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*
     -                 USf(Sfe3,2,3,2)*VChaC(Cha3,2))))/(CB2*SB2)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup26 = 3*CW*ZNeu(Neu3,2) - SW*ZNeu(Neu3,1)
                  dup27 = 3*CW*ZNeuC(Neu3,2) - SW*ZNeuC(Neu3,1)

                  dZfLs = dZfLs + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MS2,MNeu2(Neu3),MSf2(Sfe3,4,2)))*
     -                 (MW*(CB*dup26*USfC(Sfe3,1,4,2)) - 
     -                 3*CW*MS*(USfC(Sfe3,2,4,2)*ZNeu(Neu3,3)))*
     -                 (MW*(CB*dup27*USf(Sfe3,1,4,2)) - 
     -                 3*CW*MS*(USf(Sfe3,2,4,2)*ZNeuC(Neu3,3))))/CB2 + 
     -                 (Alfa*MS)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MS2,MNeu2(Neu3),MSf2(Sfe3,4,2)))*
     -                 ((2*MW*(CB*SW*USf(Sfe3,2,4,2)*ZNeu(Neu3,1)) + 
     -                 3*CW*MS*(USf(Sfe3,1,4,2)*ZNeu(Neu3,3)))*
     -                 (MW*(CB*dup26*USfC(Sfe3,1,4,2)) - 
     -                 3*CW*MS*(USfC(Sfe3,2,4,2)*ZNeu(Neu3,3))) + 
     -                 (MW*(CB*dup27*USf(Sfe3,1,4,2)) - 
     -                 3*CW*MS*(USf(Sfe3,2,4,2)*ZNeuC(Neu3,3)))*
     -                 (2*MW*(CB*SW*USfC(Sfe3,2,4,2)*ZNeuC(Neu3,1)) + 
     -                 3*CW*MS*(USfC(Sfe3,1,4,2)*ZNeuC(Neu3,3))))*
     -                 MNeu(Neu3))/CB2 + 
     -                 (Alfa*MS)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MS2,MNeu2(Neu3),MSf2(Sfe3,4,2)))*
     -                 ((MS*((MW2*(CB2*(SW2*ZNeuC(Neu3,1) - 
     -                 3*CW*(SW*ZNeuC(Neu3,2)))*ZNeu(Neu3,1))-3*MW2*
     -                 (CB2*(CW*(SW*ZNeuC(Neu3,1)) - 
     -                 3*CW2*ZNeuC(Neu3,2))*ZNeu(Neu3,2)) + 
     -                 9*CW2*MS2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,1,4,2)) + 9*MS2*MW*
     -                 (CB*(CW*(SW*ZNeuC(Neu3,1)) - CW2*ZNeuC(Neu3,2))*
     -                 USfC(Sfe3,2,4,2)*ZNeu(Neu3,3)))*
     -                 USf(Sfe3,1,4,2) + 
     -                 (MS*((4*MW2*SW2*(CB2*ZNeu(Neu3,1)*ZNeuC(Neu3,1)) + 
     -                 9*CW2*MS2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,2,4,2)) + 9*MS2*MW*
     -                 (CB*(CW*(SW*ZNeu(Neu3,1)) - CW2*ZNeu(Neu3,2))*
     -                 USfC(Sfe3,1,4,2)*ZNeuC(Neu3,3)))*
     -                 USf(Sfe3,2,4,2)))/CB2
               enddo
            enddo

            do Sfe3 = 1,2
               do Cha3 = 1,2
                  dZfRs = dZfRs + (Alfa*MS2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 UCha(Cha3,2)*UChaC(Cha3,2)*USf(Sfe3,1,3,2)*
     -                 USfC(Sfe3,1,3,2))/CB2 + 
     -                 (Alfa*MS2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 MCha(Cha3)*((MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,1,3,2)*VCha(Cha3,1)) - 
     -                 MC*(USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*UCha(Cha3,2)*
     -                 USf(Sfe3,1,3,2) + 
     -                 (MW*(SB*Sqrt(2.D0)*USf(Sfe3,1,3,2)*VChaC(Cha3,1)) - 
     -                 MC*(USf(Sfe3,2,3,2)*VChaC(Cha3,2)))*
     -                 UChaC(Cha3,2)*USfC(Sfe3,1,3,2)))/(CB*SB) - 
     -                 (Alfa*MS2)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MS2,MCha2(Cha3),MSf2(Sfe3,3,2)))*
     -                 (-(MS2*(SB2*UCha(Cha3,2)*UChaC(Cha3,2)*
     -                 USf(Sfe3,1,3,2)*USfC(Sfe3,1,3,2))) + 
     -                 CB2*((-(2*MW2*(SB2*USfC(Sfe3,1,3,2)*
     -                 VCha(Cha3,1))) + MC*MW*
     -                 (SB*Sqrt(2.D0)*USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*
     -                 USf(Sfe3,1,3,2)*VChaC(Cha3,1) + 
     -                 (MC*MW*(SB*Sqrt(2.D0)*USfC(Sfe3,1,3,2)*
     -                 VCha(Cha3,1)) - 
     -                 MC2*(USfC(Sfe3,2,3,2)*VCha(Cha3,2)))*
     -                 USf(Sfe3,2,3,2)*VChaC(Cha3,2))))/(CB2*SB2)
               enddo
            enddo

            do Sfe3 = 1,2
               do Neu3 = 1,4
                  dup28 = 2*MW*(CB*SW*USf(Sfe3,2,4,2)*ZNeu(Neu3,1)) + 
     -                 3*CW*MS*(USf(Sfe3,1,4,2)*ZNeu(Neu3,3))
                  dup29 = 2*MW*(CB*SW*USfC(Sfe3,2,4,2)*ZNeuC(Neu3,1)) + 
     -                 3*CW*MS*(USfC(Sfe3,1,4,2)*ZNeuC(Neu3,3))

                  dZfRs = dZfRs + Alfa/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (dup28*dup29*DBLE(B0i(bb1,MS2,MNeu2(Neu3),
     -                 MSf2(Sfe3,4,2))))/CB2 + 
     -                 (Alfa*MS)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MS2,MNeu2(Neu3),MSf2(Sfe3,4,2)))*
     -                 (dup28*(MW*(CB*USfC(Sfe3,1,4,2)*
     -                 (3*CW*ZNeu(Neu3,2) - SW*ZNeu(Neu3,1))) - 
     -                 3*CW*MS*(USfC(Sfe3,2,4,2)*ZNeu(Neu3,3))) + 
     -                 dup29*(MW*(CB*USf(Sfe3,1,4,2)*
     -                 (3*CW*ZNeuC(Neu3,2) - SW*ZNeuC(Neu3,1))) - 
     -                 3*CW*MS*(USf(Sfe3,2,4,2)*ZNeuC(Neu3,3))))*
     -                 MNeu(Neu3))/CB2 + 
     -                 (Alfa*MS)/(72.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb1,MS2,MNeu2(Neu3),MSf2(Sfe3,4,2)))*
     -                 ((MS*((MW2*(CB2*(SW2*ZNeuC(Neu3,1) - 
     -                 3*CW*(SW*ZNeuC(Neu3,2)))*ZNeu(Neu3,1))-3*MW2*
     -                 (CB2*(CW*(SW*ZNeuC(Neu3,1)) - 
     -                 3*CW2*ZNeuC(Neu3,2))*ZNeu(Neu3,2)) + 
     -                 9*CW2*MS2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,1,4,2)) + 9*MS2*MW*
     -                 (CB*(CW*(SW*ZNeuC(Neu3,1)) - CW2*ZNeuC(Neu3,2))*
     -                 USfC(Sfe3,2,4,2)*ZNeu(Neu3,3)))*
     -                 USf(Sfe3,1,4,2) + (MS*((4*MW2*SW2*(CB2*
     -                 ZNeu(Neu3,1)*ZNeuC(Neu3,1)) + 
     -                 9*CW2*MS2*(ZNeu(Neu3,3)*ZNeuC(Neu3,3)))*
     -                 USfC(Sfe3,2,4,2)) + 9*MS2*MW*
     -                 (CB*(CW*(SW*ZNeu(Neu3,1)) - CW2*ZNeu(Neu3,2))*
     -                 USfC(Sfe3,1,4,2)*ZNeuC(Neu3,3)))*
     -                 USf(Sfe3,2,4,2)))/CB2
               enddo
            enddo

** End of contributions from chargino/neutralinos
         end if


** End of contributions from SM-type bosons and their superpartners
      end if


** Electroweak mixing angle renormalisation
      dSW1 = -(CW2/(2.D0*MW2*MZ2)*(-(MZ2*dMWsq1) + MW2*dMZsq1)/SW)


*** CHARGE renormalisation
      if (charge .ne. 2) then
         
** parametrising lo by alpha(MZ)
         dZet1 = 1/(2.D0*CW)*(dZZA1*SW + CW*(-dZAA1 - dZAAlt1))

** parametrising lo by alpha(0)         
         if (charge .eq. 0) then
            dZet1 = dZet1 + (0.5D0*delAlfa)
         end if
                  
      else if (charge .eq. 2) then
** Parametrising lo by GF        

         dZet1 = dZeR + 1/(2.D0*MW2)*(dMWsq1*SW - 2*MW2*dSW1)/SW
         if ((sector .ge. 3) .or. full) then
            dZet1 = dZet1 + 1/(2.D0*CW)*dZZA1/SW - Alfa*(6.D0 + 
     -           (7.D0 - 4*SW2)*LOG(CW2)/(2.D0*SW2))/(8.D0*Pi*SW2)

            if (chiRUN) then
               dZet1 = dZeRsm + 1/(2.D0*MW2)*(dMWsq1sm*SW - 
     -              2*MW2*dSW1sm)/SW
               dZet1 = dZet1 + 1/(2.D0*CW)*dZZA1sm/SW - Alfa*(6.D0 + 
     -           (7.D0 - 4*SW2)*LOG(CW2)/(2.D0*SW2))/(8.D0*Pi*SW2)
               dZet1 = dZet1 - (dZAA1chi/2D0) + (SW*dZZA1chi/(2D0*CW))
            end if
         end if

** NOTE: At present, if parametrising LO by GF and working in the full MSSM
**       we use an approximate charge renormalisation that is given by:
**         dZet1(0) - 0.5 Delta R^{SM-type + sfms}
**       see ntbk, 21-08-2010

      end if
** End of charge renormalisation


** Higgs field renormalisation
      dZh0HH = (CA*(dZh0 - dZHH)*SA)/C2A

** Beta renormalisation
      dTB1 = -(1/2.D0*(-dZh0 + dZHH)/C2A)

** Z -> CP-odd A0 boson mass mixing term
      dMAZsq1 = -(2*(0D0,1D0)*MZ*(CB*dTB1*SB))


***** FERMION VALUES *****
*** Electroweak mixing angle renormalisation
      dSW1fm = -(CW2/(2.D0*MW2*MZ2)*(-(MZ2*dMWsq1fm) + MW2*dMZsq1fm)/SW)

*** CHARGE renormalisation
      if (charge .ne. 2) then
         
** parametrising lo by alpha(MZ)
         dZet1fm = 1/(2.D0*CW)*(dZZA1fm*SW + CW*(-dZAA1fm - dZAAlt1fm))

** parametrising lo by alpha(0)         
         if (charge .eq. 0) then
            dZet1fm = dZet1fm + (0.5D0*delAlfa)
         end if
                  
      else if (charge .eq. 2) then
** Parametrising lo by GF        

         dZet1fm = dZeRfm + 1/(2.D0*MW2)*(dMWsq1fm*SW - 2*MW2*dSW1fm)/SW

      end if
** End of charge renormalisation


** Higgs field renormalisation
      dZh0HHfm = (CA*(dZh0fm - dZHHfm)*SA)/C2A

** Beta renormalisation
      dTB1fm = -(1/2.D0*(-dZh0fm + dZHHfm)/C2A)

** Z -> CP-odd A0 boson mass mixing term
      dMAZsq1fm = -(2*(0D0,1D0)*MZ*(CB*dTB1fm*SB))
***** END FERMION VALUES *****


** Renaming (with suffix full) all necessary renormalisation constants
      if (full) then
         dMWsq1full = dMWsq1
         dMZsq1full = dMZsq1
         dSW1full = dSW1
         dZet1full = dZet1
         dZW1full = dZW1
         dZZZ1full = dZZZ1
         dZh0full = dZh0
         dZHHfull = dZHH
         dZh0HHfull = dZh0HH
         dTB1full = dTB1
         dZZA1full = dZZA1
      end if
      
 
** Debugging output
      if (debug .ne. 0) then
         write(*,*)'         '
         if (full) write(*,*)'Ren Cons in the FULL MSSM:'
         write(*,*)'Renormalisation constants at Mh0 =', MHiggs_loop(1)
         write(*,*)'dMWsq1 =', dMWsq1
         write(*,*)'dMZsq1 =', dMZsq1
         write(*,*)'dSW1 =', dSW1
         write(*,*)'dZW1 =', dZW1
         write(*,*)'dZZZ1 =', dZZZ1
         write(*,*)'dZh0 =', dZh0
         write(*,*)'dZHH =', dZHH
         write(*,*)'dZh0HH =', dZh0HH
         write(*,*)'dTB1 =', dTB1
         write(*,*)'dMAZsq1 =', dMAZsq1
         write(*,*)'dZAA1 =', dZAA1
         write(*,*)'dZAZ1 =', dZAZ1
         write(*,*)'dZZA1 =', dZZA1
         write(*,*)'dZet1 =', dZet1
         write(*,*)'dZfLu =', dZfLu
         write(*,*)'dZfLd =', dZfLd
         write(*,*)'dZfLc =', dZfLc
         write(*,*)'dZfLs =', dZfLs
         write(*,*)'dZfRu =', dZfRu
         write(*,*)'dZfRd =', dZfRd
         write(*,*)'dZfRc =', dZfRc
         write(*,*)'dZfRs =', dZfRs
      end if
      
 
** Debugging output: fermion values 
      if (debug .ne. 0) then
         write(*,*)'         '
         write(*,*)'FERMION VALUES!:'
         write(*,*)'Renormalisation constants at Mh0 =', MHiggs_loop(1)
         write(*,*)'dMWsq1 =', dMWsq1fm
         write(*,*)'dMZsq1 =', dMZsq1fm
         write(*,*)'dSW1 =', dSW1fm
         write(*,*)'dZW1 =', dZW1fm
         write(*,*)'dZZZ1 =', dZZZ1fm
         write(*,*)'dZh0 =', dZh0fm
         write(*,*)'dZHH =', dZHHfm
         write(*,*)'dZh0HH =', dZh0HHfm
         write(*,*)'dTB1 =', dTB1fm
         write(*,*)'dMAZsq1 =', dMAZsq1fm
         write(*,*)'dZAA1 =', dZAA1fm
         write(*,*)'dZAZ1 =', dZAZ1fm
         write(*,*)'dZZA1 =', dZZA1fm
         write(*,*)'dZet1 =', dZet1fm
      end if


      end


********************************************************************************
********************************************************************************


********************************************************************************
********************************************************************************
***                                                                          ***
*** The next set of subroutines contains the code needed in order to         ***
*** calculate the Z propagator matrix (at one-loop level) if FeynHiggs is    ***
*** not being used.  As such, it contains (renormalised) self energies and   ***
*** some renormalisation constants.  The code was generated using FeynArts   ***
*** and FormCalc, using reference hep-ph/0611326.                            ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine will calculate the Higgs propagator Z matrices if we're not
*** using FeynHiggs

      subroutine slha_zmatrix(debug)

      implicit none

** higgs mass and field renormalisation constants
      double complex dMHiggs(3,3), dZHiggs(3,3)


** renormalised self energies and derivative: (storage arrays and functions)
** labelling: rse(i,j,k) = RSE_{ij}(M_{k}^{2})
      double complex rse(3,3,3), drse(3,3,3)
      double complex rseFN, drseFN
      external rseFN
      external drseFN


** do loop parameters Z_ij
      integer i, j, k, noti, nj, nk


** dummy parameters, used to calculate individual terms in the Z factors
      double complex termB, termC, termD, termE
      double complex termBprime, termCprime, termDprime, termEprime
      double complex num, denom

** components of Z factors
      double complex Znorm(3), Zprop(3,3)

      integer debug
      

#include "mssm.inc"

      
** Clearing the values of all arrays and dummy parameters
      termB = (0D0,0D0)
      termC = (0D0,0D0)
      termD = (0D0,0D0)
      termE = (0D0,0D0)
      termBprime = (0D0,0D0)
      termCprime = (0D0,0D0)
      termDprime = (0D0,0D0)
      termEprime = (0D0,0D0)

      do i = 1, 3
         Znorm(i) = (0D0,0D0)
         do j = 1, 3
            Hgs(i,j) = (0D0,0D0)
            Zprop(i,j) = (0D0,0D0)
            do k = 1, 3
               rse(i,j,k) = (0D0,0D0)
               rse(i,j,k) = (0D0,0D0)
            end do
         end do
      end do


** Calculating the renormalisation constants needed for the renormalised self
** energies
      call Zfac_rencon(dMHiggs,dZHiggs,debug)
      

** Self energies required for both the MSSM with real and with complex 
** parameters
      rse(1,1,1) = rseFN(dMHiggs,dZHiggs,1,1,1,debug)
      rse(1,1,2) = rseFN(dMHiggs,dZHiggs,1,1,2,debug)
      rse(1,2,1) = rseFN(dMHiggs,dZHiggs,1,2,1,debug)
      rse(1,2,2) = rseFN(dMHiggs,dZHiggs,1,2,2,debug)
      rse(2,1,1) = rse(1,2,1)
      rse(2,1,2) = rse(1,2,2)
      rse(2,2,1) = rseFN(dMHiggs,dZHiggs,2,2,1,debug)
      rse(2,2,2) = rseFN(dMHiggs,dZHiggs,2,2,2,debug)

      drse(1,1,1) = drseFN(dZHiggs,1,1,1,debug)
      drse(1,1,2) = drseFN(dZHiggs,1,1,2,debug)
      drse(1,2,1) = drseFN(dZHiggs,1,2,1,debug)
      drse(1,2,2) = drseFN(dZHiggs,1,2,2,debug)  
      drse(2,1,1) = drse(1,2,1)
      drse(2,1,2) = drse(1,2,2)
      drse(2,2,1) = drseFN(dZHiggs,2,2,1,debug)
      drse(2,2,2) = drseFN(dZHiggs,2,2,2,debug)


** AND NOW, the calculation you've all been waiting for ..

      if (higgsmix .eq. 2) then !!  MSSM with REAL parameters
    
         do i = 1, 2

            if (i .eq. 1) noti = 2
            if (i .eq. 2) noti = 1

            Znorm(i) = MHiggs_loop(i)**2 - MHiggs_tree(noti)**2 + 
     -           rse(noti,noti,i)
            Znorm(i) = Znorm(i)*2D0*rse(i,noti,i)*drse(i,noti,i)
            Znorm(i) = Znorm(i) - (rse(i,noti,i)**2)*(1D0+
     -           drse(noti,noti,i))
            Znorm(i) = Znorm(i)/((MHiggs_loop(i)**2 - 
     -           MHiggs_tree(noti)**2 + rse(noti,noti,i))**2)
            
            Znorm(i) = 1D0 + DBLE(drse(i,i,i)) - DBLE(Znorm(i))
            Znorm(i) = 1D0/Znorm(i)
            
            do j = 1, 2
               if (j .eq. i) then
                  Zprop(i,j) = (1D0,0D0)
               else
                  Zprop(i,j) = MHiggs_loop(i)**2 - MHiggs_tree(j)**2 +
     -                 rse(j,j,i)
                  Zprop(i,j) = -rse(i,j,i)/Zprop(i,j)
               end if
               
               Hgs(i,j) = sqrt(Znorm(i))*Zprop(i,j)
               
            end do
            
         end do
         
         Hgs(1,3) = (0D0,0D0)
         Hgs(2,3) = (0D0,0D0)
         Hgs(3,3) = (1D0,0D0)
         Hgs(3,2) = (0D0,0D0)
         Hgs(3,1) = (0D0,0D0)
         
      else                      !!  MSSM with COMPLEX parameters

         write(*,*)'Sorry, we can not calculate the Z factors in'
         write(*,*)'the MSSM with complex parameters yet.  Please'
         write(*,*)'use FeynHiggs or set propLevel = 0.'
         stop

** First, setting the self energies that are not needed for the MSSM with real
** parameters
         rse(1,1,3) = rseFN(dMHiggs,dZHiggs,1,1,3,debug)
         rse(1,2,3) = rseFN(dMHiggs,dZHiggs,1,2,3,debug)
         rse(2,1,3) = rse(1,2,3)
         rse(2,2,3) = rseFN(dMHiggs,dZHiggs,2,2,3,debug)

         rse(3,1,1) = rseFN(dMHiggs,dZHiggs,3,1,1,debug)
         rse(3,1,2) = rseFN(dMHiggs,dZHiggs,3,1,2,debug)
         rse(3,1,3) = rseFN(dMHiggs,dZHiggs,3,1,3,debug)
         rse(1,3,1) = rse(3,1,1)
         rse(1,3,2) = rse(3,1,2)
         rse(1,3,3) = rse(3,1,3)
         rse(3,2,1) = rseFN(dMHiggs,dZHiggs,3,2,1,debug)
         rse(3,2,2) = rseFN(dMHiggs,dZHiggs,3,2,2,debug)
         rse(3,2,3) = rseFN(dMHiggs,dZHiggs,3,2,3,debug)
         rse(2,3,1) = rse(3,2,1)
         rse(2,3,2) = rse(3,2,2)
         rse(2,3,3) = rse(3,2,3)
         rse(3,3,1) = rseFN(dMHiggs,dZHiggs,3,3,1,debug)
         rse(3,3,2) = rseFN(dMHiggs,dZHiggs,3,3,2,debug)
         rse(3,3,3) = rseFN(dMHiggs,dZHiggs,3,3,3,debug) 

         drse(1,1,3) = drseFN(dZHiggs,1,1,3,debug)
         drse(1,2,3) = drseFN(dZHiggs,1,2,3,debug)
         drse(2,1,3) = drse(1,2,3)
         drse(2,2,3) = drseFN(dZHiggs,2,2,3,debug)

         drse(3,1,1) = drseFN(dZHiggs,3,1,1,debug)
         drse(3,1,2) = drseFN(dZHiggs,3,1,2,debug)
         drse(3,1,3) = drseFN(dZHiggs,3,1,3,debug)
         drse(1,3,1) = drse(3,1,1)
         drse(1,3,2) = drse(3,1,2)
         drse(1,3,3) = drse(3,1,3)
         drse(3,2,1) = drseFN(dZHiggs,3,2,1,debug)
         drse(3,2,2) = drseFN(dZHiggs,3,2,2,debug)
         drse(3,2,3) = drseFN(dZHiggs,3,2,3,debug)
         drse(2,3,1) = drse(3,2,1)
         drse(2,3,2) = drse(3,2,2)
         drse(2,3,3) = drse(3,2,3)
         drse(3,3,1) = drseFN(dZHiggs,3,3,1,debug)
         drse(3,3,2) = drseFN(dZHiggs,3,3,2,debug)
         drse(3,3,3) = drseFN(dZHiggs,3,3,3,debug)     

         do i = 1, 3
            
            if (i .eq. 1) then
               nj = 2
               nk = 3
            else if (i .eq. 2) then
               nj = 3
               nk = 1
            else if (i .eq. 3) then
               nj = 1
               nk = 2
            end if

            termB = 2D0*(MHiggs_loop(i)**2 + rse(i,nj,i))*
     -           (MHiggs_loop(i)**2 + rse(nj,nk,i))*
     -           (MHiggs_loop(i)**2 + rse(nk,i,i))
            termBprime = (MHiggs_loop(i)**2 + rse(i,nj,i))*
     -           (1D0 + drse(nk,i,i))
            termBprime = termBprime + ((1D0+drse(nj,nk,i))*
     -           (MHiggs_loop(i)**2 + rse(nk,i,i)))
            termBprime = termBprime*2D0*(MHiggs_loop(i)**2 +
     -           rse(i,nj,i))
            termBprime = termBprime + 2D0*(1D0 + drse(i,nj,i))*
     -           (MHiggs_loop(i)**2 + rse(nj,nk,i))*
     -           (MHiggs_loop(i)**2 + rse(nk,i,i))
            
            termC = ((MHiggs_loop(i)**2 + rse(nk,i,i))**2)*
     -           (MHiggs_loop(i)**2 - MHiggs_tree(nj)**2 +
     -           rse(nj,nj,i))
            termCprime = 2D0*(MHiggs_loop(i)**2 + rse(nk,i,i))*
     -           (1D0 + drse(nk,i,i))*(MHiggs_loop(i)**2 -
     -           MHiggs_tree(nj)**2 + rse(nj,nj,i))
            termCprime = termCprime + ((MHiggs_loop(i)**2 +
     -           rse(nk,i,i))**2)*(1D0 + drse(nj,nj,i))
            
            termD = ((MHiggs_loop(i)**2 + rse(i,nj,i))**2)*
     -           (MHiggs_loop(i)**2 - MHiggs_tree(nk)**2 +
     -           rse(nk,nk,i))
            termDprime = 2D0*(MHiggs_loop(i)**2 + rse(i,nj,i))*
     -           (1D0 + drse(i,nj,i))*(MHiggs_loop(i)**2 -
     -           MHiggs_tree(nk)**2 + rse(nk,nk,i))
            termDprime = termDprime + ((MHiggs_loop(i)**2 +
     -           rse(i,nj,i))**2)*(1D0 + drse(nk,nk,i))
            
            termE = ((MHiggs_loop(i)**2 - MHiggs_tree(nj)**2 +
     -           rse(nj,nj,i))*(MHiggs_loop(i)**2 - 
     -           MHiggs_tree(nk)**2 + rse(nk,nk,i))) - 
     -           (MHiggs_loop(i)**2 + rse(nj,nk,i))**2
            termEprime = (MHiggs_loop(i)**2 - MHiggs_tree(nj)**2 +
     -           rse(nj,nj,i))*(1D0 + drse(nk,nk,i))
            termEprime = (MHiggs_loop(i)**2 - MHiggs_tree(nk)**2 +
     -           rse(nk,nk,i))*(1D0 + drse(nj,nj,i)) + termEprime
            termEprime = -2D0*(MHiggs_loop(i)**2 + rse(nj,nk,i))*
     -           (1D0 + drse(nj,nk,i)) + termEprime
            
            Znorm(i) = termE*(termBprime - termCprime - termDprime) -
     -           termEprime*(termB - termC - termD)
            Znorm(i) = drse(i,i,i) - Znorm(i)*(0D0,1D0)/(termE**2)
            Znorm(i) = 1D0/(1D0 + DBLE(Znorm(i)))

            do j = 1, 3

               if (j .eq. i) then
                  Zprop(i,j) = (1D0,0D0)
               else

** setting k, the index that's neither i or j
                  if (i .eq. 1) then
                     if (j .eq. 2) k = 3
                     if (j .eq. 3) k = 2
                  else if (i .eq. 2) then
                     if (j .eq. 1) k = 3
                     if (j .eq. 3) k = 1
                  else if (i .eq. 3) then
                     if (j .eq. 1) k = 2
                     if (j .eq. 2) k = 1
                  end if

                  num = MHiggs_loop(i)**2 - MHiggs_tree(k)**2 + 
     -                 rse(k,k,i)
                  num = num*rse(i,j,i) - rse(j,k,i)*rse(k,i,i)

                  denom = MHiggs_loop(i)**2 - MHiggs_tree(j)**2 +
     -                 rse(j,j,i)
                  denom = (MHiggs_loop(i)**2 - MHiggs_tree(k)**2 +
     -                 rse(k,k,i))*denom
                  denom = rse(j,k,i)**2 - denom

                  Zprop(i,j) = num/denom

               end if
                  
               Hgs(i,j) = sqrt(Znorm(i))*Zprop(i,j)

            end do
         end do

      end if


** Debugging info: Output the renormalised self energies
      if (debug .ne. 0) then
         write(*,*)'  '
         do i = 1, 3
            do j = 1, 3
               write(*,*)'Self energy ', i, j
               do k = 1, 3
                  if (k .eq. 1) write(*,*)'At Mh0'
                  if (k .eq. 2) write(*,*)'At MHH'
                  if (k .eq. 3) write(*,*)'At MA0'
                  write(*,*)'rse =', rse(i,j,k)
                  write(*,*)'drse =', drse(i,j,k)
               end do
            end do
         end do
      end if


      do i = 1, 3
         do j = 1, 3
            ZHiggs(i,j) = Hgs(i,j)
         end do
      end do


      end


********************************************************************************
********************************************************************************

*** This subroutine sets or calculates the Higgs renormalisation constants
*** needed in order to calculate the renormalised self energies

      subroutine Zfac_rencon(dMHiggs,dZHiggs,debug)

      implicit none

** higgs mass, field and tadpole renormalisation constants
      double complex dMHp
      double complex dMHiggs(3,3), dZHiggs(3,3), dTadH(3)

** do loop variables
      integer i, j

** self energy A0 -> A0 function
      double complex seA0
      external seA0

** do loop parameters (run over e.g. chargino index)
      integer Gen2, Gen3, Cha2, Cha3, Sfe2, Sfe3, Sfe4, Neu2, Neu3

** dummy parameters
      double complex dup107
      double complex dup2, dup3, dup4, dup5, dup6, dup7, dup8, dup9

** em coupling constant
      double precision Alfa, Alfa2

** dummy parameters
      double complex dup58, dup59, dup60, dup61, dup62, dup63, dup64
      double complex dup65, dup66, dup67, dup68, dup69, dup70
      double complex dup71, dup72, dup73, dup74, dup75, dup76, dup77
      double complex dup78, dup79, dup80, dup81, dup82, dup83
      double precision ZHiggsD(3,3), UHiggsD(3,3)
      double complex nup1, nup2, nup3, nup4, nup5
      double complex nup6, nup7, nup8, nup9, nup10
      double complex nup11, nup12, nup13, nup14, nup15
      double complex nup16, nup17, nup18, nup19, nup20
      double complex nup21, nup22, nup23, nup24, nup25
      double complex nup26, nup27, nup28, nup29, nup30
      double complex nup31, nup32, nup33, nup34, nup35
      double complex nup36, nup37, nup38, nup39, nup40
      double complex nup41, nup42, nup43, nup44, nup45

      integer debug


#include "mssm.inc"
#include "global.inc"
#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      do cha3 = 1, 3
         do cha2 = 1, 3
            if (cha3 .eq. cha2) then
               ZHiggsD(cha3,cha2) = (1D0,0D0)
               UHiggsD(cha3,cha2) = (1D0,0D0)
            else
               ZHiggsD(cha3,cha2) = (0D0,0D0)
               UHiggsD(cha3,cha2) = (0D0,0D0)
            end if
         end do
      end do


** For safety purposes, clearing the values of the renormalisation constants
      do i = 1, 3
         dTadH(i) = (0D0,0D0)
         do j = 1, 3
            dZHiggs(i,j) = (0D0,0D0)
            dMHiggs(i,j) = (0D0,0D0)
         end do
      end do


** Firstly, the tadpoles  
** h0 tadpole:
      do Gen2 = 1,3
         dTadH(1) = dTadH(1) + 1/Pi**2*
     -        (EL*(MZ/(32.D0*CW)*(SAB*DBLE(A0(MSf2(1,1,Gen2)))) + 
     -        1/MW*(3/8.D0*(CA*DBLE(A0(Mf2(3,Gen2)))*Mf2(3,Gen2))/
     -        SB + (SA*(-(1/8.D0*
     -        (DBLE(A0(Mf2(2,Gen2)))*Mf2(2,Gen2))) - 
     -        3/8.D0*(DBLE(A0(Mf2(4,Gen2)))*Mf2(4,Gen2))))/CB)))/SW
      enddo

      do Sfe2 = 1,2
         do Gen2 = 1,3
            dTadH(1) = dTadH(1) + 1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,2,Gen2)))*
     -           (USf(Sfe2,1,2,Gen2)*
     -           (CW*((CA*MUE + SA*AfC(2,Gen2,Gen2))*Mf(2,Gen2)*
     -           USfC(Sfe2,2,2,Gen2)) + 
     -           (MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -           2*CW*(SA*Mf2(2,Gen2)))*USfC(Sfe2,1,2,Gen2)) + 
     -           USf(Sfe2,2,2,Gen2)*
     -           (CW*((CA*MUEC + SA*Af(2,Gen2,Gen2))*Mf(2,Gen2)*
     -           USfC(Sfe2,1,2,Gen2)) + 
     -           2*(-(MW*MZ*SW2*(CB*SAB)) + CW*(SA*Mf2(2,Gen2)))*
     -           USfC(Sfe2,2,2,Gen2))))/(CB*SW) - 
     -           1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,3,Gen2)))*
     -           (USf(Sfe2,1,3,Gen2)*
     -           (3*CW*((MUE*SA + CA*AfC(3,Gen2,Gen2))*Mf(3,Gen2)*
     -           USfC(Sfe2,2,3,Gen2)) + 
     -           (MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -           6*CW*(CA*Mf2(3,Gen2)))*USfC(Sfe2,1,3,Gen2)) + 
     -           USf(Sfe2,2,3,Gen2)*
     -           (3*CW*((MUEC*SA + CA*Af(3,Gen2,Gen2))*Mf(3,Gen2)*
     -           USfC(Sfe2,1,3,Gen2)) + 
     -           2*(-(2*MW*MZ*SW2*(SAB*SB)) + 
     -           3*CW*(CA*Mf2(3,Gen2)))*USfC(Sfe2,2,3,Gen2))))/
     -           (SB*SW) + 1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,4,Gen2)))*
     -           (USf(Sfe2,1,4,Gen2)*
     -           (3*CW*((CA*MUE + SA*AfC(4,Gen2,Gen2))*Mf(4,Gen2)*
     -           USfC(Sfe2,2,4,Gen2)) + 
     -           (MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -           6*CW*(SA*Mf2(4,Gen2)))*USfC(Sfe2,1,4,Gen2)) + 
     -           USf(Sfe2,2,4,Gen2)*
     -           (3*CW*((CA*MUEC + SA*Af(4,Gen2,Gen2))*Mf(4,Gen2)*
     -           USfC(Sfe2,1,4,Gen2)) - 
     -           2*(MW*MZ*SW2*(CB*SAB) - 3*CW*(SA*Mf2(4,Gen2)))*
     -           USfC(Sfe2,2,4,Gen2))))/(CB*SW)
         enddo
      enddo

      dTadH(1) = dTadH(1) + (MW*(2*CW2*MW2 + MZ2))/(16.D0*CW2*Pi**2)*
     -     (EL*SBA)/SW - MW/(64.D0*CW2*Pi**2)*
     -     (C2B*EL*SAB*DBLE(A0(MA02)))/SW - 
     -     (3*MW)/(64.D0*CW2*Pi**2)*(C2A*EL*SAB*DBLE(A0(Mh02)))/SW + 
     -     MW/(64.D0*CW2*Pi**2)*
     -     (EL*(2*CAB*S2A + C2A*SAB)*DBLE(A0(MHH2)))/SW - 
     -     (3*MW)/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MZ2))*(SBA - 1/6.D0*(C2B*SAB)))/SW - 
     -     MW/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MHp2))*(C2B*SAB + 2*CW2*SBA))/SW - 
     -     MW/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MW2))*(-(C2B*SAB) + 6*CW2*SBA))/SW
      
      do Cha2 = 1,2
         dTadH(1) = dTadH(1) - 1/(8.D0*Pi**2)*
     -        (EL*DBLE(A0(MCha2(Cha2)))*MCha(Cha2)*
     -        (SA*(UCha(Cha2,2)*VCha(Cha2,1) + 
     -        UChaC(Cha2,2)*VChaC(Cha2,1)) + 
     -        CA*(-(UCha(Cha2,1)*VCha(Cha2,2)) - 
     -        UChaC(Cha2,1)*VChaC(Cha2,2))))/(SW*Sqrt(2.D0))
      enddo

      do Neu2 = 1,4
         dup107 = SA*ZNeu(Neu2,3) + CA*ZNeu(Neu2,4)
         dTadH(1) = dTadH(1) + 1/(16.D0*CW*Pi**2)*
     -        (EL*DBLE(A0(MNeu2(Neu2)))*MNeu(Neu2)*
     -        (-(CW*(dup107*ZNeu(Neu2,2))) + 
     -        dup107*SW*ZNeu(Neu2,1) + 
     -        (-(CW*ZNeuC(Neu2,2)) + SW*ZNeuC(Neu2,1))*
     -        (SA*ZNeuC(Neu2,3) + CA*ZNeuC(Neu2,4))))/SW
      enddo


** HH tadpole:
      do Gen2 = 1, 3
         dTadH(2) = dTadH(2) + 1/Pi**2*
     -        (EL*(-(MZ/(32.D0*CW)*(CAB*DBLE(A0(MSf2(1,1,Gen2))))) + 
     -        1/MW*(3/8.D0*(SA*DBLE(A0(Mf2(3,Gen2)))*Mf2(3,Gen2))/
     -        SB + (CA*(1/8.D0*(DBLE(A0(Mf2(2,Gen2)))*Mf2(2,Gen2)) + 
     -        3/8.D0*(DBLE(A0(Mf2(4,Gen2)))*Mf2(4,Gen2))))/CB)))/SW
      enddo

      do Sfe2 = 1,2
         do Gen2 = 1,3
            dTadH(2) = dTadH(2) + 1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,2,Gen2)))*
     -           (USf(Sfe2,1,2,Gen2)*
     -           (CW*((MUE*SA - CA*AfC(2,Gen2,Gen2))*Mf(2,Gen2)*
     -           USfC(Sfe2,2,2,Gen2)) + 
     -           (MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -           2*CW*(CA*Mf2(2,Gen2)))*USfC(Sfe2,1,2,Gen2)) + 
     -           USf(Sfe2,2,2,Gen2)*
     -           (CW*((MUEC*SA - CA*Af(2,Gen2,Gen2))*Mf(2,Gen2)*
     -           USfC(Sfe2,1,2,Gen2)) + 
     -           2*(MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen2)))*
     -           USfC(Sfe2,2,2,Gen2))))/(CB*SW) - 
     -           1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,3,Gen2)))*
     -           (USf(Sfe2,1,3,Gen2)*
     -           (3*CW*((-(CA*MUE) + SA*AfC(3,Gen2,Gen2))*
     -           Mf(3,Gen2)*USfC(Sfe2,2,3,Gen2)) + 
     -           (MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -           6*CW*(SA*Mf2(3,Gen2)))*USfC(Sfe2,1,3,Gen2)) + 
     -           USf(Sfe2,2,3,Gen2)*
     -           (3*CW*((-(CA*MUEC) + SA*Af(3,Gen2,Gen2))*
     -           Mf(3,Gen2)*USfC(Sfe2,1,3,Gen2)) + 
     -           2*(2*MW*MZ*SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen2)))*
     -           USfC(Sfe2,2,3,Gen2))))/(SB*SW) - 
     -           1/(32.D0*CW*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,4,Gen2)))*
     -           (USf(Sfe2,1,4,Gen2)*
     -           (3*CW*((-(MUE*SA) + CA*AfC(4,Gen2,Gen2))*
     -           Mf(4,Gen2)*USfC(Sfe2,2,4,Gen2)) + 
     -           (MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -           6*CW*(CA*Mf2(4,Gen2)))*USfC(Sfe2,1,4,Gen2)) + 
     -           USf(Sfe2,2,4,Gen2)*
     -           (3*CW*((-(MUEC*SA) + CA*Af(4,Gen2,Gen2))*
     -           Mf(4,Gen2)*USfC(Sfe2,1,4,Gen2)) - 
     -           2*(MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen2)))*
     -           USfC(Sfe2,2,4,Gen2))))/(CB*SW)
         enddo
      enddo

      dTadH(2) = dTadH(2) + (MW*(2*CW2*MW2 + MZ2))/(16.D0*CW2*Pi**2)*
     -     (CBA*EL)/SW + MW/(64.D0*CW2*Pi**2)*
     -     (C2B*CAB*EL*DBLE(A0(MA02)))/SW + 
     -     MW/(64.D0*CW2*Pi**2)*
     -     (EL*(C2A*CAB - 2*S2A*SAB)*DBLE(A0(Mh02)))/SW - 
     -     (3*MW)/(64.D0*CW2*Pi**2)*(C2A*CAB*EL*DBLE(A0(MHH2)))/SW - 
     -     (3*MW)/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MZ2))*(CBA + 1/6.D0*(C2B*CAB)))/SW - 
     -     MW/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MHp2))*(-(C2B*CAB) + 2*CW2*CBA))/SW - 
     -     MW/(32.D0*CW2*Pi**2)*
     -     (EL*DBLE(A0(MW2))*(C2B*CAB + 6*CW2*CBA))/SW

      do Cha2 = 1,2
         dTadH(2) = dTadH(2) + 1/(8.D0*Pi**2)*
     -        (EL*DBLE(A0(MCha2(Cha2)))*MCha(Cha2)*
     -        (CA*(UCha(Cha2,2)*VCha(Cha2,1) + 
     -        UChaC(Cha2,2)*VChaC(Cha2,1)) + 
     -        SA*(UCha(Cha2,1)*VCha(Cha2,2) + 
     -        UChaC(Cha2,1)*VChaC(Cha2,2))))/(SW*Sqrt(2.D0))
      enddo

      do Neu2 = 1,4
         dTadH(2) = dTadH(2) - 1/(16.D0*CW*Pi**2)*
     -        (EL*DBLE(A0(MNeu2(Neu2)))*MNeu(Neu2)*
     -        ((-(CW*(CA*ZNeu(Neu2,3))) + CW*(SA*ZNeu(Neu2,4)))*
     -        ZNeu(Neu2,2) + SW*ZNeu(Neu2,1)*
     -        (CA*ZNeu(Neu2,3) - SA*ZNeu(Neu2,4)) + 
     -        (-(CW*ZNeuC(Neu2,2)) + SW*ZNeuC(Neu2,1))*
     -        (CA*ZNeuC(Neu2,3) - SA*ZNeuC(Neu2,4))))/SW
      enddo


* A0 tadpole:
      do Gen2 = 1,3 ! Sophy! Why no leptons?
         dTadH(3) = dTadH(3) + (3*cI)/(8.D0*MW*Pi**2)*
     -        (EL*((DBLE(A0(Mf2(3,Gen2)))*Mf2(3,Gen2))/TB + 
     -        TB*DBLE(A0(Mf2(4,Gen2)))*Mf2(4,Gen2)))/SW
      enddo

      do Sfe2 = 1,2
         do Gen2 = 1,3
            dTadH(3) = dTadH(3) + cI/(32.D0*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,2,Gen2)))*Mf(2,Gen2)*
     -           (-((MUEC + TB*Af(2,Gen2,Gen2))*USf(Sfe2,2,2,Gen2)*
     -           USfC(Sfe2,1,2,Gen2)) + 
     -           (MUE + TB*AfC(2,Gen2,Gen2))*USf(Sfe2,1,2,Gen2)*
     -           USfC(Sfe2,2,2,Gen2)))/SW + 
     -           (3*cI)/(32.D0*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,3,Gen2)))*Mf(3,Gen2)*
     -           (-(MUEC*TB*USf(Sfe2,2,3,Gen2)*USfC(Sfe2,1,3,Gen2)) - 
     -           Af(3,Gen2,Gen2)*USf(Sfe2,2,3,Gen2)*
     -           USfC(Sfe2,1,3,Gen2) + 
     -           MUE*TB*USf(Sfe2,1,3,Gen2)*USfC(Sfe2,2,3,Gen2) + 
     -           AfC(3,Gen2,Gen2)*USf(Sfe2,1,3,Gen2)*
     -           USfC(Sfe2,2,3,Gen2)))/(SW*TB) + 
     -           (3*cI)/(32.D0*MW*Pi**2)*
     -           (EL*DBLE(A0(MSf2(Sfe2,4,Gen2)))*Mf(4,Gen2)*
     -           (-((MUEC + TB*Af(4,Gen2,Gen2))*USf(Sfe2,2,4,Gen2)*
     -           USfC(Sfe2,1,4,Gen2)) + 
     -           (MUE + TB*AfC(4,Gen2,Gen2))*USf(Sfe2,1,4,Gen2)*
     -           USfC(Sfe2,2,4,Gen2)))/SW
         enddo
      enddo

      do Cha2 = 1,2
         dTadH(3) = dTadH(3) + cI/(8.D0*Pi**2)*
     -        (EL*DBLE(A0(MCha2(Cha2)))*MCha(Cha2)*
     -        (SB*(UCha(Cha2,2)*VCha(Cha2,1) - 
     -        UChaC(Cha2,2)*VChaC(Cha2,1)) + 
     -        CB*(UCha(Cha2,1)*VCha(Cha2,2) - 
     -        UChaC(Cha2,1)*VChaC(Cha2,2))))/(SW*Sqrt(2.D0))
      enddo

      do Neu2 = 1,4
         dTadH(3) = dTadH(3) - cI/(16.D0*CW*Pi**2)*
     -        (EL*DBLE(A0(MNeu2(Neu2)))*MNeu(Neu2)*
     -        ((-(CW*(SB*ZNeu(Neu2,3))) + CW*(CB*ZNeu(Neu2,4)))*
     -        ZNeu(Neu2,2) + SW*ZNeu(Neu2,1)*
     -        (SB*ZNeu(Neu2,3) - CB*ZNeu(Neu2,4)) - 
     -        (-(CW*ZNeuC(Neu2,2)) + SW*ZNeuC(Neu2,1))*
     -        (SB*ZNeuC(Neu2,3) - CB*ZNeuC(Neu2,4))))/SW
      enddo


** Now, the field renormalisation constants
      dZHiggs(1,1) = dZh0full
      dZHiggs(1,2) = dZh0HHfull
      dZHiggs(2,1) = dZh0HHfull
      dZHiggs(2,2) = dZHHfull

      if (higgsmix .eq. 3) then ! additional value needed for complex MSSM
         dZHiggs(3,3) = SB2*dZHHfull + CB2*dZh0full
      end if


** Lastly, the mass renormalisation constants
      if (higgsmix .eq. 2) then   ! MSSM with real parameters
         
         dMHiggs(3,3) = DBLE(seA0(3))

      else  ! MSSM with complex parameters

         dMHp = -(Alfa*MW2)/(4.D0*Pi*SW2) - 
     -        (Alfa*MZ2*(CW2 - SW2)**2)/(8.D0*CW2*Pi*SW2) + 
     -        Alfa/(32.D0*CW2*Pi*SW2)*(C2B**2*DBLE(A0(MA02))) + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*(C2B**2*DBLE(A0(MHp2))) - 
     -        Alfa/Pi*(MHp2*DBLE(B0i(bb0,MHp2,0.D0,MHp2))) + 
     -        Alfa/(16.D0*Pi*SW2)*
     -        ((-MA02 - MHp2 + MW2)*DBLE(B0i(bb0,MHp2,MA02,MW2))) - 
     -        (Alfa*(CW2 - SW2)**2)/(8.D0*CW2*Pi*SW2)*
     -        (MHp2*DBLE(B0i(bb0,MHp2,MHp2,MZ2))) - 
     -        Alfa/(2.D0*Pi)*(MHp2*DBLE(B0i(bb1,MHp2,0.D0,MHp2))) + 
     -        Alfa/(8.D0*Pi*SW2)*(MHp2*DBLE(B0i(bb1,MHp2,MA02,MW2))) + 
     -        Alfa/(8.D0*Pi*SW2)*(CBA**2*MHp2*
     -        DBLE(B0i(bb1,MHp2,Mh02,MW2))) + 
     -        Alfa/(8.D0*Pi*SW2)*(MHp2*SBA**2*
     -        DBLE(B0i(bb1,MHp2,MHH2,MW2))) + 
     -        (Alfa*(CW2 - SW2)**2)/(8.D0*CW2*Pi*SW2)*
     -        (MHp2*DBLE(B0i(bb1,MHp2,MHp2,MZ2))) + 
     -        DBLE(B0i(bb0,MHp2,Mh02,MW2))*
     -        (-(Alfa/(16.D0*Pi*SW2)*(CBA**2*Mh02)) - 
     -        Alfa/(16.D0*Pi*SW2)*(CBA**2*MHp2) + 
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -        (-(S2B*SAB) + CW2*CBA)**2) + 
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MHp2,MHH2,MHp2))*(C2B*CAB - 2*CW2*CBA)**2)+
     -        Alfa/(32.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(Mh02))*(SW2*(C2A*C2B) + CW2*(1 - S2A*S2B))) + 
     -        Alfa/(32.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(MHH2))*(-(SW2*(C2A*C2B)) + CW2*(1 + S2A*S2B)))+
     -        DBLE(B0i(bb0,MHp2,MHH2,MW2))*
     -        (-(Alfa/(16.D0*Pi*SW2)*(MHH2*SBA**2)) - 
     -        Alfa/(16.D0*Pi*SW2)*(MHp2*SBA**2) + 
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*(CAB*S2B - CW2*SBA)**2)+
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MHp2,Mh02,MHp2))*(C2B*SAB + 2*CW2*SBA)**2)+
     -        DBLE(A0(MW2))*(Alfa/(16.D0*CW2*Pi*SW2)*
     -        (-C2B**2 + S2B**2) + (7*Alfa)/(16.D0*Pi*SW2)*(1 - 1/7.D0*
     -        CBA**2 - 1/7.D0*SBA**2))+
     -        (3*Alfa)/(16.D0*CW2*Pi*SW2)*
     -        (DBLE(A0(MZ2))*((CW2 - SW2)**2 + 
     -        1/6.D0*(-(SW2*C2B**2) + CW2*(1 + S2B**2))))

         do Gen3 = 1,3
            dup2 = Mf2(3,Gen3) + TB2**2*Mf2(4,Gen3)

            dMHp = dMHp - (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -           (dup2*DBLE(A0(Mf2(4,Gen3))))/TB2 - 
     -           (3*Alfa)/(4.D0*MW2*Pi*SW2)*(dup2*MHp2*
     -           DBLE(B0i(bb1,MHp2,Mf2(3,Gen3),Mf2(4,Gen3))))/
     -           TB2 - Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(1,1,Gen3)))*
     -           ((MW2 - 2*CW2*MW2)*C2B - 2*CW2*(TB2*Mf2(2,Gen3)))) + 
     -           Alfa/(4.D0*MW2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHp2,Mf2(3,Gen3),Mf2(4,Gen3)))*
     -           (-3*dup2*Mf2(3,Gen3) - 6*TB2*Mf2(3,Gen3)*Mf2(4,Gen3)))
     -           /TB2
         enddo

         do Neu3 = 1,4
            do Cha3 = 1,2
               dup3 = -(1/CW*(UCha(Cha3,2)*
     -              (CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)))/Sqrt(2.D0)) + 
     -              UCha(Cha3,1)*ZNeu(Neu3,3)
               dup4 = 1/CW*(VCha(Cha3,2)*
     -              (CW*ZNeu(Neu3,2) + SW*ZNeu(Neu3,1)))/Sqrt(2.D0) + 
     -              VCha(Cha3,1)*ZNeu(Neu3,4)
               dup5 = -(1/CW*(UChaC(Cha3,2)*
     -              (CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)))/Sqrt(2.D0)) + 
     -              UChaC(Cha3,1)*ZNeuC(Neu3,3)
               dup6 = 1/CW*(VChaC(Cha3,2)*
     -              (CW*ZNeuC(Neu3,2) + SW*ZNeuC(Neu3,1)))/Sqrt(2.D0) + 
     -              VChaC(Cha3,1)*ZNeuC(Neu3,4)

               dMHp = dMHp + Alfa/(2.D0*Pi*SW2)*
     -              ((-(CB2*dup4*dup6) - dup3*dup5*SB2)*
     -              DBLE(A0(MCha2(Cha3)))) + 
     -              Alfa/(2.D0*Pi*SW2)*((CB2*dup4*dup6*MHp2 + 
     -              dup3*dup5*MHp2*SB2)*
     -              DBLE(B0i(bb1,MHp2,MCha2(Cha3),MNeu2(Neu3)))) + 
     -              Alfa/(2.D0*Pi*SW2)*(DBLE(B0i(bb0,MHp2,MCha2(Cha3),
     -              MNeu2(Neu3)))*
     -              (CB*(-(dup3*dup4*SB*MCha(Cha3)*MNeu(Neu3)) - 
     -              dup5*dup6*SB*MCha(Cha3)*MNeu(Neu3)) + 
     -              CB2*dup4*dup6*(MHp2 - MNeu2(Neu3)) + 
     -              dup3*dup5*(MHp2*SB2 - SB2*MNeu2(Neu3))))
            enddo
         enddo

         do Sfe3 = 1,2
            do Gen3 = 1,3
               dMHp = dMHp + Alfa/(8.D0*MW2*Pi*SW2)*(DBLE(B0i(bb0,MHp2,
     -              MSf2(1,1,Gen3),MSf2(Sfe3,2,Gen3)))*
     -              ((MW2*S2B - TB*Mf2(2,Gen3))*USf(Sfe3,1,2,Gen3) - 
     -              (MUEC + TB*Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -              USf(Sfe3,2,2,Gen3))*
     -              ((MW2*S2B - TB*Mf2(2,Gen3))*USfC(Sfe3,1,2,Gen3) - 
     -              (MUE + TB*AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -              USfC(Sfe3,2,2,Gen3))) - 
     -              Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -              (DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -              (MW2*(C2B*USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3)) - 
     -              2*(MW2*SW2*C2B + CW2*(TB2*Mf2(2,Gen3)))*
     -              USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3))) + 
     -              Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -              (DBLE(A0(MSf2(Sfe3,3,Gen3)))*
     -              (TB2*((MW2 + 2*CW2*MW2)*C2B + 
     -              6*CW2*(TB2*Mf2(4,Gen3)))*USf(Sfe3,1,3,Gen3)*
     -              USfC(Sfe3,1,3,Gen3) + 
     -              2*(-(2*MW2*SW2*(C2B*TB2)) + 3*CW2*Mf2(3,Gen3))*
     -              USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)))/TB2 - 
     -              Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -              (DBLE(A0(MSf2(Sfe3,4,Gen3)))*
     -              ((C2B*(-(MW2*TB2) + 4*CW2*MW2*TB2) - 
     -              6*CW2*Mf2(3,Gen3))*USf(Sfe3,1,4,Gen3)*
     -              USfC(Sfe3,1,4,Gen3) - 
     -              2*TB2*(MW2*SW2*C2B + 3*CW2*(TB2*Mf2(4,Gen3)))*
     -              USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))/TB2
            enddo
         enddo

         do Sfe4 = 1,2
            do Sfe3 = 1,2
               do Gen3 = 1,3
                  dup7 = -(MW2*(S2B*TB)) + Mf2(3,Gen3) + TB2*Mf2(4,Gen3)

                  dMHp = dMHp + (3*Alfa)/(8.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(bb0,MHp2,MSf2(Sfe3,3,Gen3),
     -                 MSf2(Sfe4,4,Gen3)))*((-(dup7*USf(Sfe4,1,4,Gen3))- 
     -                 TB*(MUEC + TB*Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                 USf(Sfe4,2,4,Gen3))*USfC(Sfe3,1,3,Gen3) - 
     -                 Mf(3,Gen3)*(MUE*TB*USf(Sfe4,1,4,Gen3) + 
     -                 AfC(3,Gen3,Gen3)*USf(Sfe4,1,4,Gen3) + 
     -                 Mf(4,Gen3)*USf(Sfe4,2,4,Gen3) + 
     -                 TB2*Mf(4,Gen3)*USf(Sfe4,2,4,Gen3))*
     -                 USfC(Sfe3,2,3,Gen3))*
     -                 ((-(dup7*USf(Sfe3,1,3,Gen3)) - 
     -                 MUEC*TB*Mf(3,Gen3)*USf(Sfe3,2,3,Gen3) - 
     -                 Af(3,Gen3,Gen3)*Mf(3,Gen3)*USf(Sfe3,2,3,Gen3))*
     -                 USfC(Sfe4,1,4,Gen3) - 
     -                 Mf(4,Gen3)*(TB*(MUE + TB*AfC(4,Gen3,Gen3))*
     -                 USf(Sfe3,1,3,Gen3) + 
     -                 (1 + TB2)*Mf(3,Gen3)*USf(Sfe3,2,3,Gen3))*
     -                 USfC(Sfe4,2,4,Gen3)))/TB2
               enddo
            enddo
         enddo

         dMHiggs(3,3) = dMHp - dMWsq1

         dMHiggs(1,3) = -EL*dTadH(3)*SBA/(2D0*MZ*SW*CW)
         dMHiggs(3,1) = dMHiggs(1,3)
         dMHiggs(2,3) = -EL*dTadH(3)*CBA/(2D0*MZ*SW*CW)
         dMHiggs(3,2) = dMHiggs(2,3)

      end if


      dMHiggs(1,1) = dMHiggs(3,3)*CBA*CBA + dMZsq1full*SAB*SAB 
      dMHiggs(1,1) = dMHiggs(1,1) + EL*(dTadH(2)*CBA*SBA*SBA - 
     -     dTadH(1)*SBA*(1D0 + CBA*CBA))/
     -     (2D0*MZ*SW*CW) 
      dMHiggs(1,1) = dMHiggs(1,1) + dTB1full*S2B*
     -     (-SBA*CBA*MHiggs_tree(3)**2 + SAB*CAB*MZ2)

      dMHiggs(1,2) = -SBA*CBA*dMHiggs(3,3) - SAB*CAB*dMZsq1full
      dMHiggs(1,2) = dMHiggs(1,2) - EL*(dTadH(2)*(SBA**3) + 
     -     dTadH(1)*(CBA**3))/(2D0*MZ*SW*CW)
      dMHiggs(1,2) = dMHiggs(1,2) - dTB1full*SB*CB*
     -     ((CBA*CBA - SBA*SBA)*MHiggs_tree(3)**2 + 
     -     (CAB*CAB - SAB*SAB)*MZ2)

      dMHiggs(2,1) = dMHiggs(1,2)

      dMHiggs(2,2) = dMHiggs(3,3)*SBA*SBA + dMZsq1full*CAB*CAB 
      dMHiggs(2,2) = dMHiggs(2,2) + EL*(dTadH(1)*CBA*CBA*SBA - 
     -     dTadH(2)*CBA*(1D0 + SBA*SBA))/
     -     (2D0*MZ*SW*CW) 
      dMHiggs(2,2) = dMHiggs(2,2) + dTB1full*S2B*
     -     (SBA*CBA*MHiggs_tree(3)**2 - SAB*CAB*MZ2)


      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'Renormalisation constants for Z matrix calculation'
         write(*,*)'dMh0 =', dMHiggs(1,1)
         write(*,*)'dMh0HH =', dMHiggs(1,2)
         write(*,*)'dMHH =', dMHiggs(2,2)
         write(*,*)'dMA0 =', dMHiggs(3,3)
         write(*,*)'dMHp =', dMHp
         write(*,*)'dZh0 =', dZHiggs(1,1)
         write(*,*)'dZh0HH =', dZHiggs(1,2)
         write(*,*)'dZHH =', dZHiggs(2,2)
         write(*,*)'dTh0 =', dTadH(1)
         write(*,*)'dTHH =', dTadH(2)
         write(*,*)'dMZsq1full =', dMZsq1full
      end if


      end


********************************************************************************
********************************************************************************

*** This function calculates the renormalised self energy

      double complex function rseFN(dMHiggs,dZHiggs,i,j,k,debug)  

      implicit none

** higgs mass and field renormalisation constants
      double complex dMHiggs(3,3), dZHiggs(3,3)

** rse labelling: rseFN(i,j,k) = rse_{ij} at (M_{k}^{2})
      integer i, j, k, dummyi, dummyj

** self energy of the Higgs combination ijk
      double complex seFN

** self energy A0 -> A0 function
      double complex seA0
      external seA0

** em coupling
      double precision Alfa, Alfa2

** do loop variables
      integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

** dummy parameters
      double complex dup124, dup125, dup126, dup127, dup128, dup129
      double complex dup130, dup131, dup132, dup133, dup134, dup135
      double complex dup136, dup137, dup138, dup139, dup140, dup141
      double complex dup142, dup143, dup144, dup145, dup146, dup147
      double complex dup148, dup149, dup150, dup151, dup152, dup153
      double complex dup154, dup155, dup156, dup157, dup158, dup159
      double complex dup160, dup161, dup162, dup163, dup164, dup165
      double complex dup166, dup167, dup168, dup169, dup170, dup171
      double complex dup172, dup173, dup174, dup175, dup176, dup177
      double complex dup178, dup179, dup180, dup181, dup182, dup183
      double complex dup184, dup185, dup186, dup187, dup188, dup189
      double complex dup190, dup191, dup192, dup193, dup194, dup195
      double complex dup196, dup197, dup198, dup199, dup200, dup201
      double complex dup202, dup203, dup204, dup205, dup206, dup207
      double complex dup208, dup209, dup210, dup211, dup212, dup213
      double complex dup214

      integer debug

#include "mssm.inc"
#include "global.inc"
#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Clearing the self energy value      
      seFN = (0D0,0D0)


** Right, we're going to set the dummy indices, such that dummyi is less than
** dummyj
      if (j .lt. i) then
         dummyi = j
         dummyj = i
      else
         dummyi = i
         dummyj = j
      end if


** Firstly, we calculate the required self energy
      if (dummyi .eq. 1) then

         if (dummyj .eq. 1) then  !  This is h0-h0 self energy

            do Gen3 = 1, 3
               seFN = seFN + Alfa/(Pi*SW2)*(1/CW2*(-(1/16.D0*(C2A*
     -              DBLE(A0(MSf2(1,1,Gen3))))) + MZ2/16.D0*(SAB**2*
     -              DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(1,1,Gen3),
     -              MSf2(1,1,Gen3))))) + 1/MW2*((CA2*(-(3/4.D0*
     -              ((DBLE(A0(Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/SB2 + 
     -              (SA2*(-(1/4.D0*((DBLE(A0(Mf2(2,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3))) - 
     -              1/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) - 3/4.D0*
     -              ((DBLE(A0(Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB2))
            enddo

            do Sfe3 = 1,2
               do Gen3 = 1,3
                  seFN = seFN - Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -                 ((MW2*(-1 + 2*SW2)*(C2A*CB2) - 
     -                 2*CW2*(SA2*Mf2(2,Gen3)))*USf(Sfe3,1,2,Gen3)*
     -                 USfC(Sfe3,1,2,Gen3) - 
     -                 2*(MW2*SW2*(C2A*CB2) + CW2*(SA2*Mf2(2,Gen3)))*
     -                 USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))/CB2 + 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,3,Gen3)))*
     -                 ((MW2*(-3 + 4*SW2)*(C2A*SB2) + 
     -                 6*CW2*(CA2*Mf2(3,Gen3)))*USf(Sfe3,1,3,Gen3)*
     -                 USfC(Sfe3,1,3,Gen3) + 2*(-(2*MW2*SW2*(C2A*SB2)) + 
     -                 3*CW2*(CA2*Mf2(3,Gen3)))*USf(Sfe3,2,3,Gen3)*
     -                 USfC(Sfe3,2,3,Gen3)))/SB2 - 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,4,Gen3)))*
     -                 ((MW2*(-3 + 2*SW2)*(C2A*CB2) - 
     -                 6*CW2*(SA2*Mf2(4,Gen3)))*USf(Sfe3,1,4,Gen3)*
     -                 USfC(Sfe3,1,4,Gen3) - 
     -                 2*(MW2*SW2*(C2A*CB2) + 3*CW2*(SA2*Mf2(4,Gen3)))*
     -                 USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))/CB2
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dup152 = -(MW*MZ*SW2*(CB*SAB)) + 
     -                    CW*(SA*Mf2(2,Gen3))
                     dup153 = -(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3))
                     dup154 = MW*MZ*SW2*(CB*SAB) - 3*CW*(SA*Mf2(4,Gen3))
                     dup155 = MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 2*CW*
     -                    (SA*Mf2(2,Gen3))
                     dup156 = MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 6*CW*
     -                    (CA*Mf2(3,Gen3))
                     dup157 = MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 6*CW*
     -                    (SA*Mf2(4,Gen3))

                     seFN = seFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + dup155*
     -                    USfC(Sfe3,1,2,Gen3)) + USf(Sfe4,2,2,Gen3)*
     -                    (CW*((CA*MUEC + SA*Af(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -                    2*dup152*USfC(Sfe3,2,2,Gen3)))*
     -                    (USf(Sfe3,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,2,2,Gen3)) + 
     -                    dup155*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + 
     -                    2*dup152*USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    dup156*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*
     -                    (3*CW*((MUEC*SA + CA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 
     -                    2*dup153*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*
     -                    (3*CW*((MUE*SA + CA*AfC(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,2,3,Gen3)) + 
     -                    dup156*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*
     -                    (3*CW*((MUEC*SA + CA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,1,3,Gen3)) + 
     -                    2*dup153*USfC(Sfe4,2,3,Gen3))))/SB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*
     -                    (3*CW*((CA*MUE + SA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,2,4,Gen3)) + 
     -                    dup157*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 
     -                    2*dup154*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*
     -                    (3*CW*((CA*MUE + SA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,2,4,Gen3)) + 
     -                    dup157*USfC(Sfe4,1,4,Gen3)) + 
     -                    USf(Sfe3,2,4,Gen3)*
     -                    (3*CW*((CA*MUEC + SA*Af(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,1,4,Gen3)) - 
     -                    2*dup154*USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

            dup141 = 1 - 1/4.D0*CBA**2 - 1/4.D0*SBA**2

            seFN = seFN + Alfa/(32.D0*CW2*Pi*SW2)*
     -           (C2A*C2B*DBLE(A0(MA02))) + 
     -           (3*Alfa)/(32.D0*CW2*Pi*SW2)*(C2A**2*DBLE(A0(Mh02))) - 
     -           Alfa/(32.D0*CW2*Pi*SW2)*((1 - 3*S2A**2)*
     -           DBLE(A0(MHH2))) + (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2B**2*SAB**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MA02))) + 
     -           (9*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2A**2*SAB**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,Mh02))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((C2A*CAB - 2*S2A*SAB)**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,MHH2))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           ((2*CAB*S2A + C2A*SAB)**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHH2,MHH2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA**2*MHp2) + (Alfa*MW2)/(8.D0*
     -           CW2**2*Pi*SW2)*(-(S2B*SAB) + CW2*CBA)**2) + Alfa/
     -           (16.D0*CW2*Pi*SW2)*(DBLE(A0(MHp2))*(SW2*(C2A*C2B) + 
     -           CW2*(1 - S2A*S2B))) + Alfa/(16.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MW2))*(-(SW2*(C2A*C2B)) + 
     -           CW2*(1 + 8*dup141 + S2A*S2B))) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2))*(-(Alfa/
     -           (16.D0*CW2*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) + 
     -           Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(CW2*(CBA**2*MA02)) + MW2*(S2B**2*SAB**2))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (C2B*SAB + 2*CW2*SBA)**2) + Alfa/(4.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MZ2))*(1 - 1/8.D0*(C2A*C2B) - 1/4.D0*CBA**2 - 
     -           1/4.D0*SBA**2)) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (C2B**2*SAB**2 + 12*CW2**2*SBA**2)) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -           SBA**2)) + Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (MW2*(C2B**2*SAB**2) + 14*MW2*SBA**2 - 
     -           2*CW2*MZ2*SBA**2)) - Alfa/(8.D0*CW2**2*Pi*SW2)*
     -           (2*MW2*SBA**2 + CW2*(MZ2 + CW2*MW2*(2 + 4*SBA**2)))
            
            seFN = seFN + Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2)))

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup142 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup143 = SA*UCha(Cha3,2)*VCha(Cha4,1) - 
     -                 CA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup144 = SA*UChaC(Cha4,2)*VChaC(Cha3,1) - 
     -                 CA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup145 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  seFN = seFN - Alfa/(4.D0*Pi*SW2)*((dup142*dup144 + 
     -                 dup143*dup145)*DBLE(A0(MCha2(Cha4)))) - 
     -                 Alfa/(4.D0*Pi*SW2)*((dup142*dup144 + dup143*
     -                 dup145)*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 Alfa/(4.D0*Pi*SW2)*(DBLE(B0i(bb0,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (-(dup144*dup145*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup143*dup145*MCha2(Cha3) + 
     -                 dup142*(-(dup143*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup144*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup146 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup147 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup148 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup149 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup150 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup151 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)

                  seFN = seFN - Alfa/(8.D0*CW2*Pi*SW2)*
     -                 (DBLE(A0(MNeu2(Neu4)))*
     -                 (dup146*dup148 + dup147*SA*ZNeu(Neu3,3) + 
     -                 CA*dup147*ZNeu(Neu3,4))*
     -                 (dup149*dup151 + dup150*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup150*ZNeuC(Neu3,4))) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*
     -                 (dup146*dup148 + dup147*SA*ZNeu(Neu3,3) + 
     -                 CA*dup147*ZNeu(Neu3,4))*
     -                 (dup149*dup151 + dup150*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup150*ZNeuC(Neu3,4))) + 
     -                 DBLE(B0i(bb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*(-(Alfa/(8.D0*CW2*Pi*SW2)*
     -                 (MNeu2(Neu3)*
     -                 (dup146*dup148 + dup147*SA*ZNeu(Neu3,3) + 
     -                 CA*dup147*ZNeu(Neu3,4))*
     -                 (dup149*dup151 + dup150*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup150*ZNeuC(Neu3,4)))) - 
     -                 Alfa/(16.D0*CW2*Pi*SW2)*(MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup146*dup148 + dup147*SA*ZNeu(Neu3,3) + 
     -                 CA*dup147*ZNeu(Neu3,4))**2 + 
     -                 (dup149*dup151 + dup150*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup150*ZNeuC(Neu3,4))**2)))
               enddo
            enddo

         else if (dummyj .eq. 2) then  !  This is the h0-HH self energy

            do Gen3 = 1,3
               seFN = seFN + Alfa/(Pi*SW2)*(-(1/(16.D0*CW2)*
     -              (S2A*DBLE(A0(MSf2(1,1,Gen3))) + MZ2*(CAB*SAB*
     -              DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(1,1,Gen3),
     -              MSf2(1,1,Gen3)))))) + 1/MW2*(CA*SA*((-(3/4.D0*
     -              ((DBLE(A0(Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2))/SB2 + 
     -              (1/4.D0*((DBLE(A0(Mf2(2,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3)) + 
     -              1/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) + 3/4.D0*
     -              ((DBLE(A0(Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) + 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2))/CB2)))
            enddo

            do Sfe3 = 1,2
               do Gen3 = 1,3
                  seFN = seFN - Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (S2A*DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -                 ((MW2*(-1 + 2*SW2)*CB2 + CW2*Mf2(2,Gen3))*
     -                 USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -                 (-(2*MW2*SW2*CB2) + CW2*Mf2(2,Gen3))*
     -                 USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))/CB2 + 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (S2A*DBLE(A0(MSf2(Sfe3,3,Gen3)))*
     -                 ((MW2*(-3 + 4*SW2)*SB2 + 3*CW2*Mf2(3,Gen3))*
     -                 USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -                 (-(4*MW2*SW2*SB2) + 3*CW2*Mf2(3,Gen3))*
     -                 USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)))/SB2 - 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (S2A*DBLE(A0(MSf2(Sfe3,4,Gen3)))*
     -                 ((MW2*(-3 + 2*SW2)*CB2 + 3*CW2*Mf2(4,Gen3))*
     -                 USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -                 (-(2*MW2*SW2*CB2) + 3*CW2*Mf2(4,Gen3))*
     -                 USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))/CB2
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     seFN = seFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                    2*CW*(SA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 2*(-(MW*MZ*SW2*
     -                    (CB*SAB)) + CW*(SA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*(USf(Sfe3,1,2,Gen3)*
     -                    (CW*((MUE*SA - CA*AfC(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe4,2,2,Gen3)) + 
     -                    (MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3)))*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*(CW*((MUEC*SA - CA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + 
     -                    2*(MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3)))*
     -                    USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                    6*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,1,3,Gen3)) + 
     -                    2*(-(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*
     -                    (3*CW*((-(CA*MUE) + SA*AfC(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,2,3,Gen3)) + 
     -                    (MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3)))*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*
     -                    (3*CW*((-(CA*MUEC) + SA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,1,3,Gen3)) + 2*(2*MW*MZ*
     -                    SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen3)))*
     -                    USfC(Sfe4,2,3,Gen3))))/SB2

                     seFN = seFN - Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*(3*CW*((CA*MUE + SA*
     -                    AfC(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                    6*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*SW2*(CB*SAB) - 
     -                    3*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*
     -                    (3*CW*((-(MUE*SA) + CA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3)))*USfC(Sfe4,1,4,Gen3)) + 
     -                    USf(Sfe3,2,4,Gen3)*
     -                    (3*CW*((-(MUEC*SA) + CA*Af(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,1,4,Gen3)) - 2*(MW*MZ*
     -                    SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen3)))*
     -                    USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

            seFN = seFN -((Alfa*(1 + 2*CW2**2)*MW2)/(4.D0*CW2**2*Pi*
     -           SW2)*(CBA*SBA)) + Alfa/(32.D0*CW2*Pi*SW2)*
     -           (C2B*S2A*DBLE(A0(MA02))) + 
     -           (3*Alfa)/(32.D0*CW2*Pi*SW2)*(C2A*S2A*DBLE(A0(Mh02))) - 
     -           (3*Alfa)/(32.D0*CW2*Pi*SW2)*(C2A*S2A*DBLE(A0(MHH2))) - 
     -           Alfa/(32.D0*CW2*Pi*SW2)*(C2B*S2A*DBLE(A0(MZ2))) - 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2B**2*CAB*SAB*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MA02))) - 
     -           (3*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2A*SAB*(C2A*CAB - 2*S2A*SAB)*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,Mh02))) - 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((-(C2A**2*CAB*SAB) + 4*CAB*S2A**2*SAB + 
     -           2*C2A*S2A*(-CAB**2 + SAB**2))*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,MHH2))) - 
     -           (3*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2A*CAB*(2*CAB*S2A + C2A*SAB)*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHH2,MHH2))) - 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) - 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2))) + 
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MHp2))*(SW2*(C2B*S2A) + CW2*(C2A*S2B))) - 
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           (DBLE(A0(MW2))*(SW2*(C2B*S2A) + CW2*(C2A*S2B))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (-(C2B*CAB) + 2*CW2*CBA)*(C2B*SAB + 2*CW2*SBA)) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(CBA*
     -           MHiggs_loop(k)**2*SBA)) + Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(C2B**2*CAB*SAB)) + 
     -           CBA*(14*MW2*SBA - 2*CW2*MZ2*SBA))) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (-(C2B**2*CAB*SAB) + 12*CW2**2*(CBA*SBA))) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2))*
     -           (Alfa/(16.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA) + 
     -           Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(CAB*S2B**2*SAB)) + CW2*(CBA*MA02*SBA))) + 
     -           Alfa/(8.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (CAB*(CW2*MW2*(CBA*S2B) - MW2*(S2B**2*SAB)) + 
     -           CW2*MW2*(S2B*SAB*SBA) + CW2**2*(CBA*(MHiggs_loop(k)**2*
     -           SBA + MHp2*SBA - MW2*SBA))))
            
            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup158 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup159 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 SA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup160 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 SA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup161 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  seFN = seFN + Alfa/(4.D0*Pi*SW2)*((dup158*dup160 + 
     -                 dup159*dup161)*DBLE(A0(MCha2(Cha4)))) + 
     -                 Alfa/(4.D0*Pi*SW2)*((dup158*dup160 + dup159*
     -                 dup161)*MHiggs_loop(k)**2*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -                 (DBLE(B0i(bb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(dup160*dup161*MCha(Cha3)*
     -                 MCha(Cha4) + dup159*dup161*MCha2(Cha3) + 
     -                 dup158*(dup159*MCha(Cha3)*MCha(Cha4) + 
     -                 dup160*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup162 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup163 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup164 = CA*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup165 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup166 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup167 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup168 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup169 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup170 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup171 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)
                  dup172 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)
                  dup173 = -((dup162*dup167 + CA*dup163*ZNeu(Neu3,3) + 
     -                 dup165*ZNeu(Neu3,4))*
     -                 (dup168*dup171 + dup169*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup169*ZNeuC(Neu3,4))) + 
     -                 (-(dup162*dup166) + dup165*ZNeu(Neu3,3) + 
     -                 dup164*ZNeu(Neu3,4))*
     -                 (dup168*dup172 + CA*dup169*ZNeuC(Neu3,3) + 
     -                 dup170*ZNeuC(Neu3,4))

                  seFN = seFN - Alfa/(16.D0*CW2*Pi*SW2)*(dup173*
     -                 DBLE(A0(MNeu2(Neu4)))) - 
     -                 Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (dup173*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(bb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*(-(Alfa/(16.D0*CW2*Pi*SW2)*(dup173*
     -                 MNeu2(Neu3))) - Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (MNeu(Neu3)*MNeu(Neu4)*
     -                 ((-(dup162*dup166) + dup165*ZNeu(Neu3,3) + 
     -                 dup164*ZNeu(Neu3,4))*
     -                 (dup162*dup167 + CA*dup163*ZNeu(Neu3,3) + 
     -                 dup165*ZNeu(Neu3,4)) - 
     -                 (dup168*dup171 + dup169*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup169*ZNeuC(Neu3,4))*
     -                 (dup168*dup172 + CA*dup169*ZNeuC(Neu3,3) + 
     -                 dup170*ZNeuC(Neu3,4)))))
               enddo
            enddo

         else if (dummyj .eq. 3) then  !  This is the h0-A0 self energy

            do Gen3 = 1,3
               seFN = seFN + (Alfa*cI)/(MW2*Pi*SW2)*
     -              ((CA*(-(3/4.D0*((DBLE(A0(Mf2(3,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/(SB*TB) + 
     -              (SA*TB*(3/4.D0*((DBLE(A0(Mf2(4,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) + 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB)
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     seFN = seFN + (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    Mf(2,Gen3)*(USf(Sfe4,1,2,Gen3)*(CW*((CA*MUE + 
     -                    SA*AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                    2*CW*(SA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 2*(-(MW*MZ*SW2*
     -                    (CB*SAB)) + CW*(SA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*
     -                    (-((MUEC + TB*Af(2,Gen3,Gen3))*
     -                    USf(Sfe3,2,2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     -                    (MUE + TB*AfC(2,Gen3,Gen3))*
     -                    USf(Sfe3,1,2,Gen3)*USfC(Sfe4,2,2,Gen3)))/CB - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    Mf(3,Gen3)*(USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*
     -                    SA + CA*AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                    6*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,1,3,Gen3)) + 
     -                    2*(-(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,2,3,Gen3)))*
     -                    (-(MUEC*TB*USf(Sfe3,2,3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3)) - Af(3,Gen3,Gen3)*
     -                    USf(Sfe3,2,3,Gen3)*USfC(Sfe4,1,3,Gen3) + MUE*
     -                    TB*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -                    AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)))/(SB*TB) + 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    Mf(4,Gen3)*(USf(Sfe4,1,4,Gen3)*(3*CW*((CA*
     -                    MUE + SA*AfC(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                    6*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*SW2*(CB*SAB) - 
     -                    3*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,2,4,Gen3)))*
     -                    (-((MUEC + TB*Af(4,Gen3,Gen3))*
     -                    USf(Sfe3,2,4,Gen3)*USfC(Sfe4,1,4,Gen3)) + 
     -                    (MUE + TB*AfC(4,Gen3,Gen3))*
     -                    USf(Sfe3,1,4,Gen3)*USfC(Sfe4,2,4,Gen3)))/CB
                  enddo
               enddo
            enddo
            
            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup124 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup125 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 CB*UCha(Cha3,1)*VCha(Cha4,2)
                  dup126 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 CB*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup127 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  seFN = seFN + (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup124*dup126) + dup125*dup127)*
     -                 DBLE(A0(MCha2(Cha4)))) + (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup124*dup126) + dup125*dup127)*
     -                 MHiggs_loop(k)**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -                 MCha2(Cha3),MCha2(Cha4)))) + (Alfa*cI)/(4.D0*Pi*
     -                 SW2)*(DBLE(B0i(bb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(-(dup126*dup127*MCha(Cha3)*
     -                 MCha(Cha4)) + dup125*dup127*MCha2(Cha3) + 
     -                 dup124*(dup125*MCha(Cha3)*MCha(Cha4) - 
     -                 dup126*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup128 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup129 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup130 = CA*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup131 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup132 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup133 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup134 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
                  dup135 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup136 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup137 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
                  dup138 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)
                  dup139 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)
                  dup140 = (dup128*dup134 + dup129*SB*ZNeu(Neu3,3) + 
     -                 dup131*ZNeu(Neu3,4))*
     -                 (dup135*dup138 + dup136*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup136*ZNeuC(Neu3,4)) + 
     -                 (-(dup128*dup133) + dup132*ZNeu(Neu3,3) + 
     -                 dup130*ZNeu(Neu3,4))*
     -                 (dup135*dup139 + dup136*SB*ZNeuC(Neu3,3) + 
     -                 dup137*ZNeuC(Neu3,4))

                  seFN = seFN + (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup140*DBLE(A0(MNeu2(Neu4)))) + 
     -                 (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup140*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(bb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*((Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup140*MNeu2(Neu3)) - (Alfa*cI)/(16.D0*CW2*Pi*
     -                 SW2)*(MNeu(Neu3)*MNeu(Neu4)*
     -                 ((-(dup128*dup133) + dup132*ZNeu(Neu3,3) + 
     -                 dup130*ZNeu(Neu3,4))*
     -                 (dup128*dup134 + dup129*SB*ZNeu(Neu3,3) + 
     -                 dup131*ZNeu(Neu3,4)) + 
     -                 (dup135*dup138 + dup136*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup136*ZNeuC(Neu3,4))*
     -                 (dup135*dup139 + dup136*SB*ZNeuC(Neu3,3) + 
     -                 dup137*ZNeuC(Neu3,4)))))
               enddo
            enddo

         end if

      else if (dummyi .eq. 2) then

         if (dummyj .eq. 2) then  ! This is the HH-HH self energy

            do Gen3 = 1,3
               seFN = seFN + Alfa/(Pi*SW2)*(1/(16.D0*CW2)*
     -              (C2A*DBLE(A0(MSf2(1,1,Gen3))) + MZ2*(CAB**2*
     -              DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(1,1,Gen3),
     -              MSf2(1,1,Gen3))))) + 1/MW2*((SA2*(-(3/4.D0*
     -              ((DBLE(A0(Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/SB2 + (CA2*
     -              (-(1/4.D0*((DBLE(A0(Mf2(2,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3))) - 
     -              1/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) - 3/4.D0*
     -              ((DBLE(A0(Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB2))
            enddo

            do Sfe3 = 1,2
               do Gen3 = 1,3
                  seFN = seFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -                 ((MW2*(-1 + 2*SW2)*(C2A*CB2) + 
     -                 2*CW2*(CA2*Mf2(2,Gen3)))*USf(Sfe3,1,2,Gen3)*
     -                 USfC(Sfe3,1,2,Gen3) + 
     -                 2*(-(MW2*SW2*(C2A*CB2)) + CW2*(CA2*Mf2(2,Gen3)))*
     -                 USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3)))/CB2 + 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,3,Gen3)))*
     -                 ((MW2*(3 - 4*SW2)*(C2A*SB2) + 
     -                 6*CW2*(SA2*Mf2(3,Gen3)))*USf(Sfe3,1,3,Gen3)*
     -                 USfC(Sfe3,1,3,Gen3) + 2*(2*MW2*SW2*(C2A*SB2) + 3*
     -                CW2*(SA2*Mf2(3,Gen3)))*
     -                 USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)))/SB2 + 
     -                 Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                 (DBLE(A0(MSf2(Sfe3,4,Gen3)))*
     -                 ((MW2*(-3 + 2*SW2)*(C2A*CB2) + 
     -                 6*CW2*(CA2*Mf2(4,Gen3)))*USf(Sfe3,1,4,Gen3)*
     -                 USfC(Sfe3,1,4,Gen3) - 
     -                 2*(MW2*SW2*(C2A*CB2) - 3*CW2*(CA2*Mf2(4,Gen3)))*
     -                 USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))/CB2
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dup203 = MUEC*SA - CA*Af(2,Gen3,Gen3)
                     dup204 = -(CA*MUEC) + SA*Af(3,Gen3,Gen3)
                     dup205 = -(MUEC*SA) + CA*Af(4,Gen3,Gen3)
                     dup206 = MUE*SA - CA*AfC(2,Gen3,Gen3)
                     dup207 = -(CA*MUE) + SA*AfC(3,Gen3,Gen3)
                     dup208 = -(MUE*SA) + CA*AfC(4,Gen3,Gen3)
                     dup209 = MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3))
                     dup210 = 2*MW*MZ*SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen3))
                     dup211 = MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen3))
                     dup212 = MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3))
                     dup213 = MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3))
                     dup214 = MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3))

                     seFN = seFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*
     -                    (CW*(dup206*Mf(2,Gen3)*USfC(Sfe3,2,2,Gen3)) + 
     -                    dup212*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*
     -                    (CW*(dup203*Mf(2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -                    2*dup209*USfC(Sfe3,2,2,Gen3)))*
     -                    (USf(Sfe3,1,2,Gen3)*
     -                    (CW*(dup206*Mf(2,Gen3)*USfC(Sfe4,2,2,Gen3)) + 
     -                    dup212*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*
     -                    (CW*(dup203*Mf(2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     -                    2*dup209*USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*(dup207*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + dup213*
     -                    USfC(Sfe3,1,3,Gen3)) + USf(Sfe4,2,3,Gen3)*(3*
     -                    CW*(dup204*Mf(3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 
     -                    2*dup210*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*(3*CW*(dup207*Mf(3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)) + 
     -                    dup213*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*(3*CW*(dup204*Mf(3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3)) + 
     -                    2*dup210*USfC(Sfe4,2,3,Gen3))))/SB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*(3*CW*(dup208*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    dup214*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*(dup205*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 
     -                    2*dup211*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*(3*CW*(dup208*Mf(4,Gen3)*
     -                    USfC(Sfe4,2,4,Gen3)) + 
     -                    dup214*USfC(Sfe4,1,4,Gen3)) + 
     -                    USf(Sfe3,2,4,Gen3)*
     -                    (3*CW*(dup205*Mf(4,Gen3)*USfC(Sfe4,1,4,Gen3)) - 
     -                    2*dup211*USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

            dup141 = 1 - 1/4.D0*CBA**2 - 1/4.D0*SBA**2

            seFN = seFN - (Alfa/(32.D0*CW2*Pi*SW2)*
     -           (C2A*C2B*DBLE(A0(MA02)))) - Alfa/(32.D0*CW2*Pi*SW2)*
     -           ((1 - 3*S2A**2)*DBLE(A0(Mh02))) + 
     -           (3*Alfa)/(32.D0*CW2*Pi*SW2)*(C2A**2*DBLE(A0(MHH2))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2B**2*CAB**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MA02))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           ((C2A*CAB - 2*S2A*SAB)**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,Mh02))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((2*CAB*S2A + C2A*SAB)**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,MHH2))) + 
     -           (9*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2A**2*CAB**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHH2,MHH2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (C2B*CAB - 2*CW2*CBA)**2) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (C2B**2*CAB**2 + 12*CW2**2*CBA**2)) - 
     -           Alfa/(8.D0*CW2**2*Pi*SW2)*(CW2*(2*CW2*MW2 + MZ2) + 
     -           (2*MW2 + 4*CW2**2*MW2)*CBA**2)+
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))*(-(Alfa/
     -           (16.D0*CW2*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) + 
     -           Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (MW2*(C2B**2*CAB**2) + (14*MW2 - 2*CW2*MZ2)*CBA**2))+
     -           Alfa/(16.D0*CW2*Pi*SW2)*(DBLE(A0(MW2))*(SW2*(C2A*C2B) + 
     -           CW2*(1 + 8*dup141 - S2A*S2B))) + 
     -           Alfa/(16.D0*CW2*Pi*SW2)*(DBLE(A0(MHp2))*(-(SW2*
     -           (C2A*C2B)) + CW2*(1 + S2A*S2B)))+
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) - 
     -           Alfa/(8.D0*Pi*SW2)*(MHp2*SBA**2) + (Alfa*MW2)/(8.D0*
     -           CW2**2*Pi*SW2)*(CAB*S2B - CW2*SBA)**2)+
     -           Alfa/(4.D0*CW2*Pi*SW2)*(DBLE(A0(MZ2))*(1 + 1/8.D0*
     -           (C2A*C2B) - 1/4.D0*CBA**2 - 1/4.D0*SBA**2)) + 
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) - 
     -           Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(CAB**2*S2B**2)) + CW2*(MA02*SBA**2)))
            
            seFN = seFN + Alfa/(4.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2)))

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup191 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     -                 SA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup192 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 SA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup193 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 SA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup194 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -                 SA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  seFN = seFN - Alfa/(4.D0*Pi*SW2)*((dup191*dup193 + 
     -                 dup192*dup194)*DBLE(A0(MCha2(Cha4)))) - 
     -                 Alfa/(4.D0*Pi*SW2)*((dup191*dup193 + dup192*
     -                 dup194)*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 Alfa/(4.D0*Pi*SW2)*(DBLE(B0i(bb0,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (-(dup193*dup194*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup192*dup194*MCha2(Cha3) + 
     -                 dup191*(-(dup192*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup193*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup195 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup196 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup197 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup198 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup199 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup200 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup201 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup202 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

                  seFN = seFN - Alfa/(8.D0*CW2*Pi*SW2)*
     -                 (DBLE(A0(MNeu2(Neu4)))*
     -                 (dup195*dup198 + CA*dup196*ZNeu(Neu3,3) + 
     -                 dup197*ZNeu(Neu3,4))*
     -                 (dup199*dup202 + CA*dup200*ZNeuC(Neu3,3) + 
     -                 dup201*ZNeuC(Neu3,4))) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*
     -                 (dup195*dup198 + CA*dup196*ZNeu(Neu3,3) + 
     -                 dup197*ZNeu(Neu3,4))*
     -                 (dup199*dup202 + CA*dup200*ZNeuC(Neu3,3) + 
     -                 dup201*ZNeuC(Neu3,4))) + DBLE(B0i(bb0,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (-(Alfa/(8.D0*CW2*Pi*SW2)*(MNeu2(Neu3)*
     -                 (dup195*dup198 + CA*dup196*ZNeu(Neu3,3) + 
     -                 dup197*ZNeu(Neu3,4))*
     -                 (dup199*dup202 + CA*dup200*ZNeuC(Neu3,3) + 
     -                 dup201*ZNeuC(Neu3,4)))) - 
     -                 Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup195*dup198 + CA*dup196*ZNeu(Neu3,3) + 
     -                 dup197*ZNeu(Neu3,4))**2 + 
     -                 (dup199*dup202 + CA*dup200*ZNeuC(Neu3,3) + 
     -                 dup201*ZNeuC(Neu3,4))**2)))
               enddo
            enddo

         else if (dummyj .eq. 3) then  !  This is the HH-A0 self energy

            do Gen3 = 1,3
               seFN = seFN + (Alfa*cI)/(MW2*Pi*SW2)*
     -              ((SA*(-(3/4.D0*((DBLE(A0(Mf2(3,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/(SB*TB) + 
     -              (CA*TB*(-(3/4.D0*((DBLE(A0(Mf2(4,Gen3))) + 
     -              MHiggs_loop(k)**2*
     -              DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB)
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     seFN = seFN + (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    Mf(2,Gen3)*(USf(Sfe4,1,2,Gen3)*(CW*((MUE*SA - 
     -                    CA*AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((MUEC*SA - CA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 
     -                    2*(MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*(-((MUEC + TB*
     -                    Af(2,Gen3,Gen3))*USf(Sfe3,2,2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + (MUE + TB*
     -                    AfC(2,Gen3,Gen3))*USf(Sfe3,1,2,Gen3)*
     -                    USfC(Sfe4,2,2,Gen3)))/CB - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*(DBLE(B0i(bb0,
     -                    MHiggs_loop(k)**2,MSf2(Sfe3,3,Gen3),
     -                    MSf2(Sfe4,3,Gen3)))*Mf(3,Gen3)*
     -                    (USf(Sfe4,1,3,Gen3)*
     -                    (3*CW*((-(CA*MUE) + SA*AfC(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*
     -                    (3*CW*((-(CA*MUEC) + SA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 2*(2*MW*MZ*
     -                    SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen3)))*
     -                    USfC(Sfe3,2,3,Gen3)))*(-(MUEC*TB*
     -                    USf(Sfe3,2,3,Gen3)*USfC(Sfe4,1,3,Gen3)) - 
     -                    Af(3,Gen3,Gen3)*USf(Sfe3,2,3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3) + MUE*TB*
     -                    USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -                    AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)))/(SB*TB) - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    Mf(4,Gen3)*(USf(Sfe4,1,4,Gen3)*
     -                    (3*CW*((-(MUE*SA) + CA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*
     -                    (3*CW*((-(MUEC*SA) + CA*Af(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*
     -                    SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen3)))*
     -                    USfC(Sfe3,2,4,Gen3)))*(-((MUEC + TB*
     -                    Af(4,Gen3,Gen3))*USf(Sfe3,2,4,Gen3)*
     -                    USfC(Sfe4,1,4,Gen3)) + (MUE + TB*
     -                    AfC(4,Gen3,Gen3))*USf(Sfe3,1,4,Gen3)*
     -                    USfC(Sfe4,2,4,Gen3)))/CB
                  enddo
               enddo
            enddo

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup174 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     -                 SA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup175 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 CB*UCha(Cha3,1)*VCha(Cha4,2)
                  dup176 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 CB*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup177 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -                 SA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  seFN = seFN - (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup174*dup176) + dup175*dup177)*
     -                 DBLE(A0(MCha2(Cha4)))) - 
     -                 (Alfa*cI)/(4.D0*Pi*SW2)*((-(dup174*dup176) + 
     -                 dup175*dup177)*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 (Alfa*cI)/(4.D0*Pi*SW2)*(DBLE(B0i(bb0,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (dup176*dup177*MCha(Cha3)*MCha(Cha4) - 
     -                 dup175*dup177*MCha2(Cha3) + 
     -                 dup174*(-(dup175*MCha(Cha3)*MCha(Cha4)) + 
     -                 dup176*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup178 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup179 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup180 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup181 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup182 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
                  dup183 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup184 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup185 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup186 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
                  dup187 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup188 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)
                  dup189 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)
                  dup190 = (dup178*dup183 + CA*dup179*ZNeu(Neu3,3) + 
     -                 dup181*ZNeu(Neu3,4))*
     -                 (dup184*dup188 + dup185*SB*ZNeuC(Neu3,3) + 
     -                 dup186*ZNeuC(Neu3,4)) - 
     -                 (dup178*dup182 + dup179*SB*ZNeu(Neu3,3) + 
     -                 dup180*ZNeu(Neu3,4))*
     -                 (dup184*dup189 + CA*dup185*ZNeuC(Neu3,3) + 
     -                 dup187*ZNeuC(Neu3,4))

                  seFN = seFN + (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup190*DBLE(A0(MNeu2(Neu4)))) + 
     -                 (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup190*MHiggs_loop(k)**2*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(bb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*((Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup190*MNeu2(Neu3)) - (Alfa*cI)/(16.D0*CW2*Pi*
     -                 SW2)*(MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup178*dup182 + dup179*SB*ZNeu(Neu3,3) + 
     -                 dup180*ZNeu(Neu3,4))*
     -                 (dup178*dup183 + CA*dup179*ZNeu(Neu3,3) + 
     -                 dup181*ZNeu(Neu3,4)) - 
     -                 (dup184*dup188 + dup185*SB*ZNeuC(Neu3,3) + 
     -                 dup186*ZNeuC(Neu3,4))*
     -                 (dup184*dup189 + CA*dup185*ZNeuC(Neu3,3) + 
     -                 dup187*ZNeuC(Neu3,4)))))
               enddo
            enddo

         end if
  
      else  !  This is the A0-A0 self energy

         if ((higgsmix .eq. 2) .and. (k .eq. 3)) then
            seFN = dMHiggs(3,3)
         else
            seFN = seA0(k)
         end if

      end if


** And finally, the renormalised self energy
      if (i .eq. j) then  ! hh, HH or AA
         rseFN = dZHiggs(i,i)*(MHiggs_loop(k)**2 - MHiggs_tree(i)**2)
         rseFN = rseFN + seFN - dMHiggs(i,i)
      else 
         rseFN = 0.5D0*(MHiggs_tree(i)**2 + MHiggs_tree(j)**2)
         rseFN = dZHiggs(i,j)*(MHiggs_loop(k)**2 - rseFN)
         rseFN = rseFN + seFN - dMHiggs(i,j)
      end if

      
      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'rse: i,j,k =', i,j,k
         write(*,*)'rse =', rseFN
      end if


      end


********************************************************************************
********************************************************************************

*** This function calculates the derivative of the renormalised self energy

      double complex function drseFN(dZHiggs,i,j,k,debug)

      implicit none

** higgs mass and field renormalisation constants
      double complex dZHiggs(3,3)

** rse labelling: drseFN(i,j,k) = drse_{ij} at (M_{k}^{2})
      integer i, j, k, dummyi, dummyj

** derivative of self energy function
      double complex dseFN

** em coupling
      double precision Alfa, Alfa2

** do loop parameters
      integer Cha2, Cha3, Cha4, Gen2, Gen3, Neu2, Neu3, Neu4, Sfe2
      integer Sfe3, Sfe4

** dummy parameters
      double complex dup10, dup100, dup101, dup102, dup103, dup104
      double complex dup105, dup106, dup11, dup12
      double complex dup13, dup14, dup15
      double complex dup16, dup17, dup18, dup19, dup2, dup20, dup21
      double complex dup22, dup23, dup24, dup25, dup26, dup27, dup28
      double complex dup29, dup3, dup30, dup31, dup32, dup33, dup34
      double complex dup35, dup36, dup37, dup38, dup39, dup4, dup40
      double complex dup41, dup42, dup43, dup44, dup45, dup46, dup47
      double complex dup48, dup49, dup5, dup50, dup51, dup52, dup53
      double complex dup54, dup55, dup56, dup57, dup58, dup59, dup6
      double complex dup60, dup61, dup62, dup63, dup64, dup65, dup66
      double complex dup67, dup68, dup69, dup7, dup70, dup71, dup72
      double complex dup73, dup74, dup75, dup76, dup77, dup78, dup79
      double complex dup8, dup80, dup81, dup82, dup83, dup84, dup85
      double complex dup86, dup87, dup88, dup89, dup9, dup90, dup91
      double complex dup92, dup93, dup94, dup95, dup96, dup97, dup98
      double complex dup99

      integer debug

#include "mssm.inc"
#include "looptools.h"
#include "global.inc"


** Clearing the value of dseFN
      dseFN = (0D0,0D0)


      Alfa = AlfaQED
      Alfa2 = AlfaQED2
        

** Right, we're going to set the dummy indices, such that dummyi is less than
** dummyj
      if (j .lt. i) then
         dummyi = j
         dummyj = i
      else
         dummyi = i
         dummyj = j
      end if

** Firstly, we calculate the required self energy
      if (dummyi .eq. 1) then

         if (dummyj .eq. 1) then  !  This is h0-h0 self energy

            dseFN = dseFN - (Alfa/(16.D0*CW2*Pi*SW2)*
     -           (CBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2)))) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))) - 
     -           Alfa/(8.D0*Pi*SW2)*(SBA**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))) - 
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           (SBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (CBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (SBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2B**2*SAB**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MA02))) + 
     -           (9*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2A**2*SAB**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,Mh02))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((C2A*CAB - 2*S2A*SAB)**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,MHH2))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           ((2*CAB*S2A + C2A*SAB)**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHH2,MHH2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA**2*MHp2) + (Alfa*MW2)/(8.D0*
     -           CW2**2*Pi*SW2)*(-(S2B*SAB) + CW2*CBA)**2
     -           ) + DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(CBA**2*
     -           MHiggs_loop(k)**2)) + Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(CW2*(CBA**2*MA02)) + MW2*(S2B**2*SAB**2))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (C2B*SAB + 2*CW2*SBA)**2) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (C2B**2*SAB**2 + 12*CW2**2*SBA**2)) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MZ2,MZ2))*(-(Alfa/
     -           (16.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) + 
     -           Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (MW2*(C2B**2*SAB**2) + 14*MW2*SBA**2 - 
     -           2*CW2*MZ2*SBA**2))

            dseFN = dseFN + Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -           SBA**2*DBLE(B0i(dbb1,MHiggs_loop(k)**2,MZ2,MZ2)))

            do Gen3 = 1,3
               dseFN = dseFN + Alfa/(Pi*SW2)*(MZ2/(16.D0*CW2)*
     -              (SAB**2*DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -              MSf2(1,1,Gen3),MSf2(1,1,Gen3)))) + 
     -              1/MW2*((CA2*(-(3/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/SB2 + 
     -              (SA2*(-(1/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3))) - 
     -              1/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) - 
     -              3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB2))
            enddo

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup34 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup35 = SA*UCha(Cha3,2)*VCha(Cha4,1) - 
     -                 CA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup36 = SA*UChaC(Cha4,2)*VChaC(Cha3,1) - 
     -                 CA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup37 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  dseFN = dseFN - Alfa/(4.D0*Pi*SW2)*((dup34*dup36 + 
     -                 dup35*dup37)*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -                 MCha2(Cha3),MCha2(Cha4)))) - Alfa/(4.D0*Pi*SW2)*
     -                 ((dup34*dup36 + dup35*dup37)*MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(-(dup36*dup37*MCha(Cha3)*
     -                 MCha(Cha4)) - dup35*dup37*MCha2(Cha3) + 
     -                 dup34*(-(dup35*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup36*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup38 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup39 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup40 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup41 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup42 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup43 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)

                  dseFN = dseFN - Alfa/(8.D0*CW2*Pi*SW2)*(DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (dup38*dup40 + dup39*SA*ZNeu(Neu3,3) + 
     -                 CA*dup39*ZNeu(Neu3,4))*
     -                 (dup41*dup43 + dup42*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup42*ZNeuC(Neu3,4))) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*(dup38*dup40 + dup39*SA*
     -                 ZNeu(Neu3,3) + CA*dup39*ZNeu(Neu3,4))*
     -                 (dup41*dup43 + dup42*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup42*ZNeuC(Neu3,4))) + DBLE(B0i(dbb0,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (-(Alfa/(8.D0*CW2*Pi*SW2)*(MNeu2(Neu3)*
     -                 (dup38*dup40 + dup39*SA*ZNeu(Neu3,3) + 
     -                 CA*dup39*ZNeu(Neu3,4))*
     -                 (dup41*dup43 + dup42*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup42*ZNeuC(Neu3,4)))) - Alfa/(16.D0*CW2*Pi*
     -                 SW2)*(MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup38*dup40 + dup39*SA*ZNeu(Neu3,3) + 
     -                 CA*dup39*ZNeu(Neu3,4))**2 + 
     -                 (dup41*dup43 + dup42*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup42*ZNeuC(Neu3,4))**2)))
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dup44 = -(MW*MZ*SW2*(CB*SAB)) + CW*(SA*Mf2(2,Gen3))
                     dup45 = -(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3))
                     dup46 = MW*MZ*SW2*(CB*SAB) - 3*CW*(SA*Mf2(4,Gen3))
                     dup47 = MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                    2*CW*(SA*Mf2(2,Gen3))
                     dup48 = MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                    6*CW*(CA*Mf2(3,Gen3))
                     dup49 = MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                    6*CW*(SA*Mf2(4,Gen3))

                     dseFN = dseFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + 
     -                    dup47*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 
     -                    2*dup44*USfC(Sfe3,2,2,Gen3)))*
     -                    (USf(Sfe3,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,2,2,Gen3)) + 
     -                    dup47*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + 
     -                    2*dup44*USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    dup48*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,1,3,Gen3)) + 
     -                    2*dup45*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)) + 
     -                    dup48*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3)) + 
     -                    2*dup45*USfC(Sfe4,2,3,Gen3))))/SB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*(3*CW*((CA*MUE + SA*
     -                    AfC(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    dup49*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 
     -                    2*dup46*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*(3*CW*((CA*MUE + SA*
     -                    AfC(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe4,2,4,Gen3)) + 
     -                    dup49*USfC(Sfe4,1,4,Gen3)) + 
     -                    USf(Sfe3,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe4,1,4,Gen3)) - 
     -                    2*dup46*USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

         else if (dummyj .eq. 2) then  !  This is the h0-HH self energy

            dseFN = dseFN + Alfa/(16.D0*CW2*Pi*SW2)*
     -           (CBA*SBA*DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(8.D0*Pi*SW2)*(CBA*SBA*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA*SBA*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))) - 
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           (CBA*SBA*DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))) - 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (CBA*SBA*DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) - 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*SBA*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (CBA*SBA*DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2))) - 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2B**2*CAB*SAB*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MA02))) - 
     -           (3*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2A*SAB*(C2A*CAB - 2*S2A*SAB)*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,Mh02))) - 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((-(C2A**2*CAB*SAB) + 4*CAB*S2A**2*SAB + 
     -           2*C2A*S2A*(-CAB**2 + SAB**2))*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,MHH2))) - 
     -           (3*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           (C2A*CAB*(2*CAB*S2A + C2A*SAB)*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHH2,MHH2))) - 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MA02,MZ2))) - 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (-(C2B*CAB) + 2*CW2*CBA)*(C2B*SAB + 2*CW2*SBA)) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MZ2,MZ2))*(-(Alfa/
     -           (16.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA)) + 
     -           Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(C2B**2*CAB*SAB)) + 
     -           CBA*(14*MW2*SBA - 2*CW2*MZ2*SBA))) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (-(C2B**2*CAB*SAB) + 12*CW2**2*(CBA*SBA))) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MZ2))*
     -           (Alfa/(16.D0*CW2*Pi*SW2)*(CBA*MHiggs_loop(k)**2*SBA) + 
     -           Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(CAB*S2B**2*SAB)) + CW2*(CBA*MA02*SBA))) + 
     -           Alfa/(8.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (CAB*(CW2*MW2*(CBA*S2B) - MW2*(S2B**2*SAB)) + 
     -           CW2*MW2*(S2B*SAB*SBA) + CW2**2*(CBA*(MHiggs_loop(k)**2*
     -           SBA + MHp2*SBA - MW2*SBA))))
            
            dseFN = dseFN + Alfa/(8.D0*CW2*Pi*SW2)*(CBA*
     -           MHiggs_loop(k)**2*SBA*DBLE(B0i(dbb1,MHiggs_loop(k)**2,
     -           MZ2,MZ2)))

            do Gen3 = 1,3
               dseFN = dseFN + Alfa/(Pi*SW2)*(-(MZ2/(16.D0*CW2)*
     -              (CAB*SAB*DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -              MSf2(1,1,Gen3),MSf2(1,1,Gen3))))) + 
     -              1/MW2*(CA*SA*((-(3/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2))/SB2 + 
     -              (1/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3)) + 
     -              1/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) + 
     -              3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) + 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2))/CB2)))
            enddo


            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup50 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup51 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 SA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup52 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 SA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup53 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  dseFN = dseFN + Alfa/(4.D0*Pi*SW2)*((dup50*dup52 + 
     -                 dup51*dup53)*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -                 MCha2(Cha3),MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -                 ((dup50*dup52 + dup51*dup53)*MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(dup52*dup53*MCha(Cha3)*
     -                 MCha(Cha4) + dup51*dup53*MCha2(Cha3) + 
     -                 dup50*(dup51*MCha(Cha3)*MCha(Cha4) + 
     -                 dup52*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup54 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup55 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup56 = CA*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup57 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup58 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup59 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup60 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup61 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup62 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup63 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)
                  dup64 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)
                  dup65 = -((dup54*dup59 + CA*dup55*ZNeu(Neu3,3) + 
     -                 dup57*ZNeu(Neu3,4))*
     -                 (dup60*dup63 + dup61*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup61*ZNeuC(Neu3,4))) + 
     -                 (-(dup54*dup58) + dup57*ZNeu(Neu3,3) + 
     -                 dup56*ZNeu(Neu3,4))*
     -                 (dup60*dup64 + CA*dup61*ZNeuC(Neu3,3) + 
     -                 dup62*ZNeuC(Neu3,4))

                  dseFN = dseFN - Alfa/(16.D0*CW2*Pi*SW2)*(dup65*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4))))-Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (dup65*MHiggs_loop(k)**2*DBLE(B0i(dbb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(dbb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*(-(Alfa/(16.D0*CW2*Pi*SW2)*(dup65*
     -                 MNeu2(Neu3))) - Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (MNeu(Neu3)*MNeu(Neu4)*
     -                 ((-(dup54*dup58) + dup57*ZNeu(Neu3,3) + 
     -                 dup56*ZNeu(Neu3,4))*
     -                 (dup54*dup59 + CA*dup55*ZNeu(Neu3,3) + 
     -                 dup57*ZNeu(Neu3,4)) - 
     -                 (dup60*dup63 + dup61*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup61*ZNeuC(Neu3,4))*
     -                 (dup60*dup64 + CA*dup61*ZNeuC(Neu3,3) + 
     -                 dup62*ZNeuC(Neu3,4)))))
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dseFN = dseFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*(CW*((CA*MUE + SA*
     -                    AfC(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                    2*CW*(SA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((CA*MUEC + SA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 2*(-(MW*MZ*SW2*
     -                    (CB*SAB)) + CW*(SA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*(USf(Sfe3,1,2,Gen3)*
     -                    (CW*((MUE*SA - CA*AfC(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe4,2,2,Gen3)) + 
     -                    (MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3)))*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*(CW*((MUEC*SA - CA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + 
     -                    2*(MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3)))*
     -                    USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                    6*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,1,3,Gen3)) + 
     -                    2*(-(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*
     -                    (3*CW*((-(CA*MUE) + SA*AfC(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,2,3,Gen3)) + 
     -                    (MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3)))*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*
     -                    (3*CW*((-(CA*MUEC) + SA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe4,1,3,Gen3)) + 2*(2*MW*MZ*
     -                    SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen3)))*
     -                    USfC(Sfe4,2,3,Gen3))))/SB2
                     
                     dseFN = dseFN - Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*(3*CW*((CA*MUE + SA*
     -                    AfC(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                    6*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*SW2*(CB*SAB) - 
     -                    3*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*
     -                    (3*CW*((-(MUE*SA) + CA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3)))*USfC(Sfe4,1,4,Gen3)) + 
     -                    USf(Sfe3,2,4,Gen3)*
     -                    (3*CW*((-(MUEC*SA) + CA*Af(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe4,1,4,Gen3)) - 
     -                    2*(MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*
     -                    Mf2(4,Gen3)))*USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

         else if (dummyj .eq. 3) then  !  This is the h0-A0 self energy

            do Gen3 = 1,3
               dseFN = dseFN + (Alfa*cI)/(MW2*Pi*SW2)*((CA*(-(3/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/(SB*TB) + 
     -              (SA*TB*(3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -              Mf2(4,Gen3),Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) + 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB)
            enddo

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup17 = SA*UCha(Cha4,2)*VCha(Cha3,1) - 
     -                 CA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup18 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 CB*UCha(Cha3,1)*VCha(Cha4,2)
                  dup19 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 CB*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup20 = SA*UChaC(Cha3,2)*VChaC(Cha4,1) - 
     -                 CA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  dseFN = dseFN + (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup17*dup19) + dup18*dup20)*DBLE(B0i(bb1,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 (Alfa*cI)/(4.D0*Pi*SW2)*((-(dup17*dup19) + dup18*
     -                 dup20)*MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) + (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(-(dup19*dup20*MCha(Cha3)*
     -                 MCha(Cha4)) + dup18*dup20*MCha2(Cha3) + 
     -                 dup17*(dup18*MCha(Cha3)*MCha(Cha4) - 
     -                 dup19*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup21 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup22 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup23 = CA*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup24 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup25 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup26 = SA*ZNeu(Neu4,3) + CA*ZNeu(Neu4,4)
                  dup27 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
                  dup28 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup29 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup30 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
                  dup31 = SA*ZNeuC(Neu4,3) + CA*ZNeuC(Neu4,4)
                  dup32 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)
                  dup33 = (dup21*dup27 + dup22*SB*ZNeu(Neu3,3) + 
     -                 dup24*ZNeu(Neu3,4))*
     -                 (dup28*dup31 + dup29*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup29*ZNeuC(Neu3,4)) + 
     -                 (-(dup21*dup26) + dup25*ZNeu(Neu3,3) + 
     -                 dup23*ZNeu(Neu3,4))*
     -                 (dup28*dup32 + dup29*SB*ZNeuC(Neu3,3) + 
     -                 dup30*ZNeuC(Neu3,4))

                  dseFN = dseFN + (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup33*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -                 MNeu2(Neu3),MNeu2(Neu4))))+
     -                 (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup33*MHiggs_loop(k)**2*DBLE(B0i(dbb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(dbb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*((Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup33*MNeu2(Neu3)) - 
     -                 (Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (MNeu(Neu3)*MNeu(Neu4)*
     -                 ((-(dup21*dup26) + dup25*ZNeu(Neu3,3) + 
     -                 dup23*ZNeu(Neu3,4))*
     -                 (dup21*dup27 + dup22*SB*ZNeu(Neu3,3) + 
     -                 dup24*ZNeu(Neu3,4)) + 
     -                 (dup28*dup31 + dup29*SA*ZNeuC(Neu3,3) + 
     -                 CA*dup29*ZNeuC(Neu3,4))*
     -                 (dup28*dup32 + dup29*SB*ZNeuC(Neu3,3) + 
     -                 dup30*ZNeuC(Neu3,4)))))
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dseFN = dseFN + (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    Mf(2,Gen3)*(USf(Sfe4,1,2,Gen3)*
     -                    (CW*((CA*MUE + SA*AfC(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(-1 + 2*SW2)*(CB*SAB) + 
     -                    2*CW*(SA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*
     -                    (CW*((CA*MUEC + SA*Af(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -                    2*(-(MW*MZ*SW2*(CB*SAB)) + CW*(SA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*(-((MUEC + TB*
     -                    Af(2,Gen3,Gen3))*USf(Sfe3,2,2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + (MUE + TB*
     -                    AfC(2,Gen3,Gen3))*USf(Sfe3,1,2,Gen3)*
     -                    USfC(Sfe4,2,2,Gen3)))/CB - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),
     -                    MSf2(Sfe4,3,Gen3)))*Mf(3,Gen3)*
     -                    (USf(Sfe4,1,3,Gen3)*(3*CW*((MUE*SA + CA*
     -                    AfC(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(-3 + 4*SW2)*(SAB*SB) + 
     -                    6*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*(3*CW*((MUEC*SA + CA*
     -                    Af(3,Gen3,Gen3))*Mf(3,Gen3)*
     -                    USfC(Sfe3,1,3,Gen3)) + 
     -                    2*(-(2*MW*MZ*SW2*(SAB*SB)) + 
     -                    3*CW*(CA*Mf2(3,Gen3)))*USfC(Sfe3,2,3,Gen3)))*
     -                    (-(MUEC*TB*USf(Sfe3,2,3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3)) - 
     -                    Af(3,Gen3,Gen3)*USf(Sfe3,2,3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3) + 
     -                    MUE*TB*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -                    AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)))/(SB*TB) + 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    Mf(4,Gen3)*(USf(Sfe4,1,4,Gen3)*
     -                    (3*CW*((CA*MUE + SA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CB*SAB) + 
     -                    6*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*(3*CW*((CA*MUEC + SA*
     -                    Af(4,Gen3,Gen3))*Mf(4,Gen3)*
     -                    USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*SW2*(CB*SAB) - 
     -                    3*CW*(SA*Mf2(4,Gen3)))*USfC(Sfe3,2,4,Gen3)))*
     -                    (-((MUEC + TB*Af(4,Gen3,Gen3))*
     -                    USf(Sfe3,2,4,Gen3)*USfC(Sfe4,1,4,Gen3)) + 
     -                    (MUE + TB*AfC(4,Gen3,Gen3))*
     -                    USf(Sfe3,1,4,Gen3)*USfC(Sfe4,2,4,Gen3)))/CB
                  enddo
               enddo
            enddo

         end if  ! end of h0-A0 self energy

      else if (dummyi .eq. 2) then

         if (dummyj .eq. 2) then  ! This is the HH-HH self energy

            do Gen3 = 1,3
               dseFN = dseFN + Alfa/(Pi*SW2)*(MZ2/(16.D0*CW2)*
     -              (CAB**2*DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -              MSf2(1,1,Gen3),MSf2(1,1,Gen3)))) + 1/MW2*((SA2*
     -              (-(3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -              Mf2(3,Gen3),Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/SB2 + 
     -              (CA2*(-(1/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3))))*Mf2(2,Gen3))) - 
     -              1/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -              Mf2(2,Gen3)))*Mf2(2,Gen3)**2) - 
     -              3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3)) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB2))
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dup95 = MUEC*SA - CA*Af(2,Gen3,Gen3)
                     dup96 = -(CA*MUEC) + SA*Af(3,Gen3,Gen3)
                     dup97 = -(MUEC*SA) + CA*Af(4,Gen3,Gen3)
                     dup98 = MUE*SA - CA*AfC(2,Gen3,Gen3)
                     dup99 = -(CA*MUE) + SA*AfC(3,Gen3,Gen3)
                     dup100 = -(MUE*SA) + CA*AfC(4,Gen3,Gen3)
                     dup101 = MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3))
                     dup102 = 2*MW*MZ*SW2*(CAB*SB) + 
     -                    3*CW*(SA*Mf2(3,Gen3))
                     dup103 = MW*MZ*SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen3))
                     dup104 = MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3))
                     dup105 = MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3))
                     dup106 = MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3))
                     
                     dseFN = dseFN + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    (USf(Sfe4,1,2,Gen3)*
     -                    (CW*(dup98*Mf(2,Gen3)*USfC(Sfe3,2,2,Gen3)) + 
     -                    dup104*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*
     -                    (CW*(dup95*Mf(2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -                    2*dup101*USfC(Sfe3,2,2,Gen3)))*
     -                    (USf(Sfe3,1,2,Gen3)*
     -                    (CW*(dup98*Mf(2,Gen3)*USfC(Sfe4,2,2,Gen3)) + 
     -                    dup104*USfC(Sfe4,1,2,Gen3)) + 
     -                    USf(Sfe3,2,2,Gen3)*
     -                    (CW*(dup95*Mf(2,Gen3)*USfC(Sfe4,1,2,Gen3)) + 
     -                    2*dup101*USfC(Sfe4,2,2,Gen3))))/CB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    (USf(Sfe4,1,3,Gen3)*
     -                    (3*CW*(dup99*Mf(3,Gen3)*USfC(Sfe3,2,3,Gen3)) + 
     -                    dup105*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*
     -                    (3*CW*(dup96*Mf(3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 
     -                    2*dup102*USfC(Sfe3,2,3,Gen3)))*
     -                    (USf(Sfe3,1,3,Gen3)*
     -                    (3*CW*(dup99*Mf(3,Gen3)*USfC(Sfe4,2,3,Gen3)) + 
     -                    dup105*USfC(Sfe4,1,3,Gen3)) + 
     -                    USf(Sfe3,2,3,Gen3)*
     -                    (3*CW*(dup96*Mf(3,Gen3)*USfC(Sfe4,1,3,Gen3)) + 
     -                    2*dup102*USfC(Sfe4,2,3,Gen3))))/SB2 + 
     -                    Alfa/(48.D0*CW2*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    (USf(Sfe4,1,4,Gen3)*(3*CW*(dup100*Mf(4,Gen3)*
     -                    USfC(Sfe3,2,4,Gen3)) + 
     -                    dup106*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*
     -                    (3*CW*(dup97*Mf(4,Gen3)*USfC(Sfe3,1,4,Gen3)) - 
     -                    2*dup103*USfC(Sfe3,2,4,Gen3)))*
     -                    (USf(Sfe3,1,4,Gen3)*(3*CW*(dup100*Mf(4,Gen3)*
     -                    USfC(Sfe4,2,4,Gen3)) + dup106*
     -                    USfC(Sfe4,1,4,Gen3)) + USf(Sfe3,2,4,Gen3)*
     -                    (3*CW*(dup97*Mf(4,Gen3)*USfC(Sfe4,1,4,Gen3)) - 
     -                    2*dup103*USfC(Sfe4,2,4,Gen3))))/CB2
                  enddo
               enddo
            enddo

            dseFN = dseFN - (Alfa/(16.D0*CW2*Pi*SW2)*
     -           (SBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MZ2)))) - 
     -           Alfa/(8.D0*Pi*SW2)*(SBA**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2))) - 
     -           Alfa/(8.D0*Pi*SW2)*(CBA**2*
     -           DBLE(B0i(bb0,MHiggs_loop(k)**2,MW2,MW2))) - 
     -           Alfa/(16.D0*CW2*Pi*SW2)*
     -           (CBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MZ2,MZ2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (SBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(SBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA**2*
     -           DBLE(B0i(bb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*
     -           (CBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,MZ2,MZ2))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2B**2*CAB**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MA02))) + 
     -           (Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*
     -           ((C2A*CAB - 2*S2A*SAB)**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,Mh02))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           ((2*CAB*S2A + C2A*SAB)**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,MHH2))) + 
     -           (9*Alfa*MW2)/(32.D0*CW2**2*Pi*SW2)*(C2A**2*CAB**2*
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHH2,MHH2))) + 
     -           Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MA02,MZ2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -           Alfa/(4.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,MW2,MW2))) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MHp2))*
     -           (C2B*CAB - 2*CW2*CBA)**2) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MW2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2)) + 
     -           (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -           (C2B**2*CAB**2 + 12*CW2**2*CBA**2)) + 
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MZ2,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(CBA**2*
     -           MHiggs_loop(k)**2)) + Alfa/(32.D0*CW2**2*Pi*SW2)*
     -           (MW2*(C2B**2*CAB**2) + (14*MW2 - 2*CW2*MZ2)*CBA**2))+
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MW2))*
     -           (-(Alfa/(8.D0*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2)) - 
     -           Alfa/(8.D0*Pi*SW2)*(MHp2*SBA**2) + (Alfa*MW2)/(8.D0*
     -           CW2**2*Pi*SW2)*(CAB*S2B - CW2*SBA)**2)+
     -           DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MZ2))*
     -           (-(Alfa/(16.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -           SBA**2)) - Alfa/(16.D0*CW2**2*Pi*SW2)*
     -           (-(MW2*(CAB**2*S2B**2)) + CW2*(MA02*SBA**2)))

            dseFN = dseFN + Alfa/(8.D0*CW2*Pi*SW2)*(CBA**2*
     -           MHiggs_loop(k)**2*DBLE(B0i(dbb1,MHiggs_loop(k)**2,
     -           MZ2,MZ2)))

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup83 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     -                 SA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup84 = CA*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 SA*UCha(Cha3,1)*VCha(Cha4,2)
                  dup85 = CA*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 SA*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup86 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -                 SA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  dseFN = dseFN - Alfa/(4.D0*Pi*SW2)*((dup83*dup85 + 
     -                 dup84*dup86)*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -                 MCha2(Cha3),MCha2(Cha4)))) - Alfa/(4.D0*Pi*SW2)*
     -                 ((dup83*dup85 + dup84*dup86)*MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))*(-(dup85*dup86*MCha(Cha3)*
     -                 MCha(Cha4)) - dup84*dup86*MCha2(Cha3) + 
     -                 dup83*(-(dup84*MCha(Cha3)*MCha(Cha4)) - 
     -                 dup85*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup87 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup88 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup89 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup90 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup91 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup92 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup93 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup94 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

                  dseFN = dseFN - Alfa/(8.D0*CW2*Pi*SW2)*
     -                 (DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*(dup87*dup90 + CA*dup88*
     -                 ZNeu(Neu3,3) + dup89*ZNeu(Neu3,4))*
     -                 (dup91*dup94 + CA*dup92*ZNeuC(Neu3,3) + 
     -                 dup93*ZNeuC(Neu3,4))) - 
     -                 Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*
     -                 DBLE(B0i(dbb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*
     -                 (dup87*dup90 + CA*dup88*ZNeu(Neu3,3) + 
     -                 dup89*ZNeu(Neu3,4))*
     -                 (dup91*dup94 + CA*dup92*ZNeuC(Neu3,3) + 
     -                 dup93*ZNeuC(Neu3,4))) + DBLE(B0i(dbb0,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))*
     -                 (-(Alfa/(8.D0*CW2*Pi*SW2)*(MNeu2(Neu3)*
     -                 (dup87*dup90 + CA*dup88*ZNeu(Neu3,3) + 
     -                 dup89*ZNeu(Neu3,4))*
     -                 (dup91*dup94 + CA*dup92*ZNeuC(Neu3,3) + 
     -                 dup93*ZNeuC(Neu3,4)))) - Alfa/(16.D0*CW2*Pi*SW2)*
     -                 (MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup87*dup90 + CA*dup88*ZNeu(Neu3,3) + 
     -                 dup89*ZNeu(Neu3,4))**2 + 
     -                 (dup91*dup94 + CA*dup92*ZNeuC(Neu3,3) + 
     -                 dup93*ZNeuC(Neu3,4))**2)))
               enddo
            enddo

         else if (dummyj .eq. 3) then  !  This is the HH-A0 self energy

            do Gen3 = 1,3
               dseFN = dseFN + (Alfa*cI)/(MW2*Pi*SW2)*((SA*(-(3/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3))))*Mf2(3,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -              Mf2(3,Gen3)))*Mf2(3,Gen3)**2)))/(SB*TB) + 
     -              (CA*TB*(-(3/4.D0*
     -              ((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -              DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3))))*Mf2(4,Gen3))) - 
     -              3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -              Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))/CB)
            enddo

            do Cha4 = 1,2
               do Cha3 = 1,2
                  dup66 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     -                 SA*UCha(Cha4,1)*VCha(Cha3,2)
                  dup67 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -                 CB*UCha(Cha3,1)*VCha(Cha4,2)
                  dup68 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -                 CB*UChaC(Cha4,1)*VChaC(Cha3,2)
                  dup69 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -                 SA*UChaC(Cha3,1)*VChaC(Cha4,2)

                  dseFN = dseFN - (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup66*dup68) + dup67*dup69)*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MCha2(Cha3),
     -                 MCha2(Cha4)))) - (Alfa*cI)/(4.D0*Pi*SW2)*
     -                 ((-(dup66*dup68) + dup67*dup69)*
     -                 MHiggs_loop(k)**2*DBLE(B0i(dbb1,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))) + 
     -                 (Alfa*cI)/(4.D0*Pi*SW2)*(DBLE(B0i(dbb0,
     -                 MHiggs_loop(k)**2,MCha2(Cha3),MCha2(Cha4)))*
     -                 (dup68*dup69*MCha(Cha3)*MCha(Cha4) - 
     -                 dup67*dup69*MCha2(Cha3) + 
     -                 dup66*(-(dup67*MCha(Cha3)*MCha(Cha4)) + 
     -                 dup68*MCha2(Cha3))))
               enddo
            enddo

            do Neu4 = 1,4
               do Neu3 = 1,4
                  dup70 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
                  dup71 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
                  dup72 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
                  dup73 = CW*(SA*ZNeu(Neu4,2)) - SA*SW*ZNeu(Neu4,1)
                  dup74 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
                  dup75 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)
                  dup76 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
                  dup77 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
                  dup78 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
                  dup79 = CW*(SA*ZNeuC(Neu4,2)) - SA*SW*ZNeuC(Neu4,1)
                  dup80 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)
                  dup81 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)
                  dup82 = (dup70*dup75 + CA*dup71*ZNeu(Neu3,3) + 
     -                 dup73*ZNeu(Neu3,4))*
     -                 (dup76*dup80 + dup77*SB*ZNeuC(Neu3,3) + 
     -                 dup78*ZNeuC(Neu3,4)) - 
     -                 (dup70*dup74 + dup71*SB*ZNeu(Neu3,3) + 
     -                 dup72*ZNeu(Neu3,4))*
     -                 (dup76*dup81 + CA*dup77*ZNeuC(Neu3,3) + 
     -                 dup79*ZNeuC(Neu3,4))

                  dseFN = dseFN + (Alfa*cI)/(16.D0*CW2*Pi*SW2)*(dup82*
     -                 DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4))))+(Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup82*MHiggs_loop(k)**2*DBLE(B0i(dbb1,
     -                 MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))) + 
     -                 DBLE(B0i(dbb0,MHiggs_loop(k)**2,MNeu2(Neu3),
     -                 MNeu2(Neu4)))*((Alfa*cI)/(16.D0*CW2*Pi*SW2)*
     -                 (dup82*MNeu2(Neu3)) - (Alfa*cI)/(16.D0*CW2*Pi*
     -                 SW2)*(MNeu(Neu3)*MNeu(Neu4)*
     -                 ((dup70*dup74 + dup71*SB*ZNeu(Neu3,3) + 
     -                 dup72*ZNeu(Neu3,4))*
     -                 (dup70*dup75 + CA*dup71*ZNeu(Neu3,3) + 
     -                 dup73*ZNeu(Neu3,4)) - 
     -                 (dup76*dup80 + dup77*SB*ZNeuC(Neu3,3) + 
     -                 dup78*ZNeuC(Neu3,4))*
     -                 (dup76*dup81 + CA*dup77*ZNeuC(Neu3,3) + 
     -                 dup79*ZNeuC(Neu3,4)))))
               enddo
            enddo

            do Sfe4 = 1,2
               do Sfe3 = 1,2
                  do Gen3 = 1,3
                     dseFN = dseFN + (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3)))*
     -                    Mf(2,Gen3)*(USf(Sfe4,1,2,Gen3)*
     -                    (CW*((MUE*SA - CA*AfC(2,Gen3,Gen3))*
     -                    Mf(2,Gen3)*USfC(Sfe3,2,2,Gen3)) + 
     -                    (MW*MZ*(1 - 2*SW2)*(CAB*CB) - 
     -                    2*CW*(CA*Mf2(2,Gen3)))*USfC(Sfe3,1,2,Gen3)) + 
     -                    USf(Sfe4,2,2,Gen3)*(CW*((MUEC*SA - CA*
     -                    Af(2,Gen3,Gen3))*Mf(2,Gen3)*
     -                    USfC(Sfe3,1,2,Gen3)) + 
     -                    2*(MW*MZ*SW2*(CAB*CB) - CW*(CA*Mf2(2,Gen3)))*
     -                    USfC(Sfe3,2,2,Gen3)))*(-((MUEC + TB*
     -                    Af(2,Gen3,Gen3))*USf(Sfe3,2,2,Gen3)*
     -                    USfC(Sfe4,1,2,Gen3)) + (MUE + TB*
     -                    AfC(2,Gen3,Gen3))*USf(Sfe3,1,2,Gen3)*
     -                    USfC(Sfe4,2,2,Gen3)))/CB - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                    Mf(3,Gen3)*(USf(Sfe4,1,3,Gen3)*
     -                    (3*CW*((-(CA*MUE) + SA*AfC(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe3,2,3,Gen3)) + 
     -                    (MW*MZ*(3 - 4*SW2)*(CAB*SB) + 
     -                    6*CW*(SA*Mf2(3,Gen3)))*USfC(Sfe3,1,3,Gen3)) + 
     -                    USf(Sfe4,2,3,Gen3)*
     -                    (3*CW*((-(CA*MUEC) + SA*Af(3,Gen3,Gen3))*
     -                    Mf(3,Gen3)*USfC(Sfe3,1,3,Gen3)) + 2*(2*MW*MZ*
     -                    SW2*(CAB*SB) + 3*CW*(SA*Mf2(3,Gen3)))*
     -                    USfC(Sfe3,2,3,Gen3)))*(-(MUEC*TB*
     -                    USf(Sfe3,2,3,Gen3)*USfC(Sfe4,1,3,Gen3)) - 
     -                    Af(3,Gen3,Gen3)*USf(Sfe3,2,3,Gen3)*
     -                    USfC(Sfe4,1,3,Gen3) + MUE*TB*
     -                    USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -                    AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -                    USfC(Sfe4,2,3,Gen3)))/(SB*TB) - 
     -                    (Alfa*cI)/(16.D0*CW*MW2*Pi*SW2)*
     -                    (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                    MSf2(Sfe3,4,Gen3),MSf2(Sfe4,4,Gen3)))*
     -                    Mf(4,Gen3)*(USf(Sfe4,1,4,Gen3)*
     -                    (3*CW*((-(MUE*SA) + CA*AfC(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,2,4,Gen3)) + 
     -                    (MW*MZ*(-3 + 2*SW2)*(CAB*CB) + 
     -                    6*CW*(CA*Mf2(4,Gen3)))*USfC(Sfe3,1,4,Gen3)) + 
     -                    USf(Sfe4,2,4,Gen3)*
     -                    (3*CW*((-(MUEC*SA) + CA*Af(4,Gen3,Gen3))*
     -                    Mf(4,Gen3)*USfC(Sfe3,1,4,Gen3)) - 2*(MW*MZ*
     -                    SW2*(CAB*CB) - 3*CW*(CA*Mf2(4,Gen3)))*
     -                    USfC(Sfe3,2,4,Gen3)))*(-((MUEC + TB*
     -                    Af(4,Gen3,Gen3))*USf(Sfe3,2,4,Gen3)*
     -                    USfC(Sfe4,1,4,Gen3)) + (MUE + TB*
     -                    AfC(4,Gen3,Gen3))*USf(Sfe3,1,4,Gen3)*
     -                    USfC(Sfe4,2,4,Gen3)))/CB
                  enddo
               enddo
            enddo

         end if  !  This is the end of the HH-A0 self energy
  
      else  !  This is the A0-A0 self energy

         dup2 = 1/CW2*(-(CW2*(CBA**2*Mh02)) + MW2*(S2B**2*SAB**2))
         dup3 = 1/CW2*(-(MW2*(CAB**2*S2B**2)) + CW2*(MHH2*SBA**2))
         dup4 = Alfa/(8.D0*Pi*SW2)*(-MHiggs_loop(k)**2 - MHp2 + MW2)

         dseFN = dseFN + dup4*DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHp2,MW2)) - 
     -        Alfa/(16.D0*CW2*Pi*SW2)*
     -        (CBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,MZ2))) - 
     -        Alfa/(16.D0*CW2*Pi*SW2)*
     -        (SBA**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MHH2,MZ2))) - 
     -        Alfa/(8.D0*Pi*SW2)*
     -        DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2)) + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*
     -        (CBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,Mh02,MZ2))) + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*
     -        (SBA**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,MHH2,MZ2))) + 
     -        Alfa/(4.D0*Pi*SW2)*
     -        DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2)) + 
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -        (C2B**2*SAB**2*
     -        DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,Mh02))) + 
     -        (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -        (C2B**2*CAB**2*
     -        DBLE(B0i(dbb0,MHiggs_loop(k)**2,MA02,MHH2))) + 
     -        Alfa/(16.D0*CW2*Pi*SW2)*
     -        ((dup2 - CBA**2*MHiggs_loop(k)**2)*
     -        DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mh02,MZ2))) + 
     -        Alfa/(16.D0*CW2*Pi*SW2)*
     -        ((-dup3 - MHiggs_loop(k)**2*SBA**2)*
     -        DBLE(B0i(dbb0,MHiggs_loop(k)**2,MHH2,MZ2))) + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*
     -        (CBA**2*MHiggs_loop(k)**2*
     -        DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mh02,MZ2))) + 
     -        Alfa/(8.D0*CW2*Pi*SW2)*
     -        (MHiggs_loop(k)**2*SBA**2*
     -        DBLE(B0i(dbb1,MHiggs_loop(k)**2,MHH2,MZ2))) + 
     -        Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*
     -        DBLE(B0i(dbb1,MHiggs_loop(k)**2,MHp2,MW2)))

         do Gen3 = 1,3
            dseFN = dseFN + Alfa/(MW2*Pi*SW2)*
     -           ((3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -           Mf2(3,Gen3))) + MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -           Mf2(3,Gen3))))*Mf2(3,Gen3)) + 
     -           3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(3,Gen3),
     -           Mf2(3,Gen3)))*Mf2(3,Gen3)**2))/TB2 + 
     -           TB2*(-(1/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -           Mf2(2,Gen3),Mf2(2,Gen3))) + MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -           Mf2(2,Gen3))))*Mf2(2,Gen3))) + 
     -           3/4.D0*((DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -           Mf2(4,Gen3))) + MHiggs_loop(k)**2*
     -           DBLE(B0i(dbb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -           Mf2(4,Gen3))))*Mf2(4,Gen3)) + 
     -           3/2.D0*(DBLE(B0i(dbb0,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -           Mf2(4,Gen3)))*Mf2(4,Gen3)**2)))
         enddo


         do Cha4 = 1,2
            do Cha3 = 1,2
               dup5 = SB*UCha(Cha4,2)*VCha(Cha3,1) + 
     -              CB*UCha(Cha4,1)*VCha(Cha3,2)
               dup6 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -              CB*UCha(Cha3,1)*VCha(Cha4,2)
               dup7 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -              CB*UChaC(Cha4,1)*VChaC(Cha3,2)
               dup8 = SB*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -              CB*UChaC(Cha3,1)*VChaC(Cha4,2)

               dseFN = dseFN - Alfa/(4.D0*Pi*SW2)*((dup5*dup7 + 
     -              dup6*dup8)*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -              MCha2(Cha3),MCha2(Cha4)))) - 
     -              Alfa/(4.D0*Pi*SW2)*((dup5*dup7 + dup6*dup8)*
     -              MHiggs_loop(k)**2*DBLE(B0i(dbb1,MHiggs_loop(k)**2,
     -              MCha2(Cha3),MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -              (DBLE(B0i(dbb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -              MCha2(Cha4)))*(dup7*dup8*MCha(Cha3)*MCha(Cha4) - 
     -              dup6*dup8*MCha2(Cha3) + dup5*(dup6*MCha(Cha3)*
     -              MCha(Cha4) - dup7*MCha2(Cha3))))
            enddo
         enddo

         do Neu4 = 1,4
            do Neu3 = 1,4
               dup9 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
               dup10 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
               dup11 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
               dup12 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
               dup13 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
               dup14 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
               dup15 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
               dup16 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)

               dseFN = dseFN - Alfa/(8.D0*CW2*Pi*SW2)*
     -              (DBLE(B0i(bb1,MHiggs_loop(k)**2,MNeu2(Neu3),
     -              MNeu2(Neu4)))*(dup12*dup9 + dup10*SB*ZNeu(Neu3,3) + 
     -              dup11*ZNeu(Neu3,4))*
     -              (dup13*dup16 + dup14*SB*ZNeuC(Neu3,3) + 
     -              dup15*ZNeuC(Neu3,4))) - Alfa/(8.D0*CW2*Pi*SW2)*
     -              (MHiggs_loop(k)**2*DBLE(B0i(dbb1,MHiggs_loop(k)**2,
     -              MNeu2(Neu3),MNeu2(Neu4)))*
     -              (dup12*dup9 + dup10*SB*ZNeu(Neu3,3) + 
     -              dup11*ZNeu(Neu3,4))*
     -              (dup13*dup16 + dup14*SB*ZNeuC(Neu3,3) + 
     -              dup15*ZNeuC(Neu3,4))) + DBLE(B0i(dbb0,
     -              MHiggs_loop(k)**2,MNeu2(Neu3),MNeu2(Neu4)))*
     -              (-(Alfa/(8.D0*CW2*Pi*SW2)*(MNeu2(Neu3)*
     -              (dup12*dup9 + dup10*SB*ZNeu(Neu3,3) + 
     -              dup11*ZNeu(Neu3,4))*
     -              (dup13*dup16 + dup14*SB*ZNeuC(Neu3,3) + 
     -              dup15*ZNeuC(Neu3,4)))) + Alfa/(16.D0*CW2*Pi*SW2)*
     -              (MNeu(Neu3)*MNeu(Neu4)*
     -              ((dup12*dup9 + dup10*SB*ZNeu(Neu3,3) + 
     -              dup11*ZNeu(Neu3,4))**2 + 
     -              (dup13*dup16 + dup14*SB*ZNeuC(Neu3,3) + 
     -              dup15*ZNeuC(Neu3,4))**2)))
            enddo
         enddo

         do Sfe4 = 1,2
            do Sfe3 = 1,2
               do Gen3 = 1,3
                  dseFN = dseFN - Alfa/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                 MSf2(Sfe3,2,Gen3),
     -                 MSf2(Sfe4,2,Gen3)))*Mf2(2,Gen3)*
     -                 (-((MUEC + TB*Af(2,Gen3,Gen3))*
     -                 USf(Sfe4,2,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     -                 (MUE + TB*AfC(2,Gen3,Gen3))*USf(Sfe4,1,2,Gen3)*
     -                 USfC(Sfe3,2,2,Gen3))*(-((MUEC + TB*
     -                 Af(2,Gen3,Gen3))*USf(Sfe3,2,2,Gen3)*
     -                 USfC(Sfe4,1,2,Gen3)) + 
     -                 (MUE + TB*AfC(2,Gen3,Gen3))*USf(Sfe3,1,2,Gen3)*
     -                 USfC(Sfe4,2,2,Gen3))) - 
     -                 (3*Alfa)/(16.D0*MW2*Pi*SW2)*
     -                 (DBLE(B0i(dbb0,MHiggs_loop(k)**2,
     -                 MSf2(Sfe3,3,Gen3),MSf2(Sfe4,3,Gen3)))*
     -                 Mf2(3,Gen3)*(-(MUEC*TB*USf(Sfe4,2,3,Gen3)*
     -                 USfC(Sfe3,1,3,Gen3)) - 
     -                 Af(3,Gen3,Gen3)*USf(Sfe4,2,3,Gen3)*
     -                 USfC(Sfe3,1,3,Gen3) + 
     -                 MUE*TB*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,2,3,Gen3) + 
     -                 AfC(3,Gen3,Gen3)*USf(Sfe4,1,3,Gen3)*
     -                 USfC(Sfe3,2,3,Gen3))*(-(MUEC*TB*
     -                 USf(Sfe3,2,3,Gen3)*USfC(Sfe4,1,3,Gen3)) - 
     -                 Af(3,Gen3,Gen3)*USf(Sfe3,2,3,Gen3)*
     -                 USfC(Sfe4,1,3,Gen3) + 
     -                 MUE*TB*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -                 AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -                 USfC(Sfe4,2,3,Gen3)))/TB2 - 
     -                 (3*Alfa)/(16.D0*MW2*Pi*SW2)*(DBLE(B0i(dbb0,
     -                 MHiggs_loop(k)**2,MSf2(Sfe3,4,Gen3),
     -                 MSf2(Sfe4,4,Gen3)))*Mf2(4,Gen3)*(-((MUEC + TB*
     -                 Af(4,Gen3,Gen3))*USf(Sfe4,2,4,Gen3)*
     -                 USfC(Sfe3,1,4,Gen3)) + 
     -                 (MUE + TB*AfC(4,Gen3,Gen3))*USf(Sfe4,1,4,Gen3)*
     -                 USfC(Sfe3,2,4,Gen3))*(-((MUEC + TB*
     -                 Af(4,Gen3,Gen3))*USf(Sfe3,2,4,Gen3)*
     -                 USfC(Sfe4,1,4,Gen3)) + 
     -                 (MUE + TB*AfC(4,Gen3,Gen3))*USf(Sfe3,1,4,Gen3)*
     -                 USfC(Sfe4,2,4,Gen3)))
               enddo
            enddo
         enddo

      end if  ! end of A0 -> A0


      drseFN = dseFN + dZHiggs(i,j)


      if (debug .ne. 0) then
         write(*,*)' '
         write(*,*)'drse: i,j,k =', i,j,k
         write(*,*)'drse =', drseFN
      end if


      end


********************************************************************************
********************************************************************************

*** This is a function used to calculate the self energy A0 -> A0.  The argument
*** is the index of the Higgs mass at which the self energy is calculated.  It 
*** is a separate function so that it can be used both in the renormalised self
*** energy calculations and in the renormalisation constant determination.

      double complex function seA0(k)

      implicit none

** Self energy
*      double complex seA0

** Index of Higgs Mass at which the self energy is evaluated
      integer k

** em coupling
      double precision Alfa, Alfa2

** Do loop variables for, e.g. chargino index
      integer Sfe3, Sfe4, Neu3, Neu4, Gen3, Cha3, Cha4

** Dummy parameters
      double complex dup2, dup3, dup4
      double complex dup112, dup113, dup114, dup115
      double complex dup116, dup117, dup118, dup119, dup120
      double complex dup121, dup122, dup123

#include "mssm.inc"
#include "looptools.h"
#include "global.inc"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

** Clearing the value of seA0
      seA0 = (0D0,0D0)

      do Gen3 = 1,3
         seA0 = seA0 - (3*Alfa)/(4.D0*MW2*Pi*SW2)*
     -        (SB2*DBLE(A0(Mf2(3,Gen3)))*Mf2(3,Gen3))/(SB2*TB2)-
     -        (3*Alfa)/(4.D0*MW2*Pi*SW2)*(SB2*MHiggs_loop(k)**2*
     -        DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(3,Gen3),Mf2(3,Gen3)))*
     -        Mf2(3,Gen3))/(SB2*TB2) + (3*Alfa*cI)/(4.D0*MW2*Pi*SW2)*
     -        (cI*SB2*DBLE(A0(Mf2(4,Gen3)))*Mf2(4,Gen3))/CB2 + 
     -        (3*Alfa*cI)/(4.D0*MW2*Pi*SW2)*(cI*SB2*MHiggs_loop(k)**2*
     -        DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(4,Gen3),
     -        Mf2(4,Gen3)))*Mf2(4,Gen3))/CB2 + Alfa/(4.D0*MW2*Pi*SW2)*
     -        (DBLE(B0i(bb0,MHiggs_loop(k)**2,Mf2(4,Gen3),Mf2(4,Gen3)))*
     -        (3*cI*(cI*SB2*Mf2(4,Gen3)**2)+3*SB2*Mf2(4,Gen3)**2))/CB2 

         seA0 = seA0 + Alfa/(Pi*SW2)*(-(1/(16.D0*CW2)*
     -        (C2B*DBLE(A0(MSf2(1,1,Gen3))))) + 
     -        1/MW2*(TB2*(-(1/4.D0*
     -        ((DBLE(A0(Mf2(2,Gen3))) + MHiggs_loop(k)**2*
     -        DBLE(B0i(bb1,MHiggs_loop(k)**2,Mf2(2,Gen3),
     -        Mf2(2,Gen3))))*Mf2(2,Gen3))))))
      enddo

      do Sfe3 = 1,2
         do Gen3 = 1,3
            seA0 = seA0 + Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,2,Gen3)))*
     -           ((MW2*(1 - 2*SW2)*C2B + 2*CW2*(TB2*Mf2(2,Gen3)))*
     -           USf(Sfe3,1,2,Gen3)*USfC(Sfe3,1,2,Gen3) + 
     -           2*(MW2*SW2*C2B + CW2*(TB2*Mf2(2,Gen3)))*
     -           USf(Sfe3,2,2,Gen3)*USfC(Sfe3,2,2,Gen3))) + 
     -           Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,3,Gen3)))*
     -           ((MW2*(C2B*(-3*TB2 + 4*SW2*TB2)) + 6*CW2*Mf2(3,Gen3))*
     -           USf(Sfe3,1,3,Gen3)*USfC(Sfe3,1,3,Gen3) + 
     -           2*(-(2*MW2*SW2*(C2B*TB2)) + 3*CW2*Mf2(3,Gen3))*
     -           USf(Sfe3,2,3,Gen3)*USfC(Sfe3,2,3,Gen3)))/TB2 + 
     -           Alfa/(16.D0*CW2*MW2*Pi*SW2)*
     -           (DBLE(A0(MSf2(Sfe3,4,Gen3)))*
     -           ((MW2*(3 - 2*SW2)*C2B + 6*CW2*(TB2*Mf2(4,Gen3)))*
     -           USf(Sfe3,1,4,Gen3)*USfC(Sfe3,1,4,Gen3) + 
     -           2*(MW2*SW2*C2B + 3*CW2*(TB2*Mf2(4,Gen3)))*
     -           USf(Sfe3,2,4,Gen3)*USfC(Sfe3,2,4,Gen3)))
         enddo
      enddo

      do Sfe4 = 1,2
         do Sfe3 = 1,2
            do Gen3 = 1,3
               seA0 = seA0 - Alfa/(16.D0*MW2*Pi*SW2)*
     -              (DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(Sfe3,2,Gen3),
     -              MSf2(Sfe4,2,Gen3)))*Mf2(2,Gen3)*
     -              (-((MUEC + TB*Af(2,Gen3,Gen3))*USf(Sfe4,2,2,Gen3)*
     -              USfC(Sfe3,1,2,Gen3)) + 
     -              (MUE + TB*AfC(2,Gen3,Gen3))*USf(Sfe4,1,2,Gen3)*
     -              USfC(Sfe3,2,2,Gen3))*
     -              (-((MUEC + TB*Af(2,Gen3,Gen3))*USf(Sfe3,2,2,Gen3)*
     -              USfC(Sfe4,1,2,Gen3)) + 
     -              (MUE + TB*AfC(2,Gen3,Gen3))*USf(Sfe3,1,2,Gen3)*
     -              USfC(Sfe4,2,2,Gen3))) - 
     -              (3*Alfa)/(16.D0*MW2*Pi*SW2)*
     -              (DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(Sfe3,3,Gen3),
     -              MSf2(Sfe4,3,Gen3)))*Mf2(3,Gen3)*
     -              (-(MUEC*TB*USf(Sfe4,2,3,Gen3)*USfC(Sfe3,1,3,Gen3)) - 
     -              Af(3,Gen3,Gen3)*USf(Sfe4,2,3,Gen3)*
     -              USfC(Sfe3,1,3,Gen3) + 
     -              MUE*TB*USf(Sfe4,1,3,Gen3)*USfC(Sfe3,2,3,Gen3) + 
     -              AfC(3,Gen3,Gen3)*USf(Sfe4,1,3,Gen3)*
     -              USfC(Sfe3,2,3,Gen3))*
     -              (-(MUEC*TB*USf(Sfe3,2,3,Gen3)*USfC(Sfe4,1,3,Gen3)) - 
     -              Af(3,Gen3,Gen3)*USf(Sfe3,2,3,Gen3)*
     -              USfC(Sfe4,1,3,Gen3) + 
     -              MUE*TB*USf(Sfe3,1,3,Gen3)*USfC(Sfe4,2,3,Gen3) + 
     -              AfC(3,Gen3,Gen3)*USf(Sfe3,1,3,Gen3)*
     -              USfC(Sfe4,2,3,Gen3)))/TB2 - 
     -              (3*Alfa)/(16.D0*MW2*Pi*SW2)*
     -              (DBLE(B0i(bb0,MHiggs_loop(k)**2,MSf2(Sfe3,4,Gen3),
     -              MSf2(Sfe4,4,Gen3)))*Mf2(4,Gen3)*
     -              (-((MUEC + TB*Af(4,Gen3,Gen3))*USf(Sfe4,2,4,Gen3)*
     -              USfC(Sfe3,1,4,Gen3)) + 
     -              (MUE + TB*AfC(4,Gen3,Gen3))*USf(Sfe4,1,4,Gen3)*
     -              USfC(Sfe3,2,4,Gen3))*
     -              (-((MUEC + TB*Af(4,Gen3,Gen3))*USf(Sfe3,2,4,Gen3)*
     -              USfC(Sfe4,1,4,Gen3)) + 
     -              (MUE + TB*AfC(4,Gen3,Gen3))*USf(Sfe3,1,4,Gen3)*
     -              USfC(Sfe4,2,4,Gen3)))
            enddo
         enddo
      enddo      

      dup2 = 1/CW2*(-(CW2*(CBA**2*Mh02)) + MW2*(S2B**2*SAB**2))
      dup3 = 1/CW2*(-(MW2*(CAB**2*S2B**2)) + CW2*(MHH2*SBA**2))
      dup4 = Alfa/(8.D0*Pi*SW2)*(-MHiggs_loop(k)**2 - MHp2 + MW2)

      seA0 = seA0 - (Alfa*(2*CW2*MW2 + MZ2))/(8.D0*CW2*Pi*SW2) + 
     -     dup4*DBLE(B0i(bb0,MHiggs_loop(k)**2,MHp2,MW2)) + 
     -     (3*Alfa)/(32.D0*CW2*Pi*SW2)*(C2B**2*DBLE(A0(MA02))) + 
     -     Alfa/(32.D0*CW2*Pi*SW2)*(C2A*C2B*DBLE(A0(Mh02))) - 
     -     Alfa/(32.D0*CW2*Pi*SW2)*(C2A*C2B*DBLE(A0(MHH2))) + 
     -     Alfa/(16.D0*CW2*Pi*SW2)*(C2B**2*DBLE(A0(MHp2))) + 
     -     (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -     (C2B**2*SAB**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,Mh02))) + 
     -     (Alfa*MW2)/(16.D0*CW2**2*Pi*SW2)*
     -     (C2B**2*CAB**2*DBLE(B0i(bb0,MHiggs_loop(k)**2,MA02,MHH2))) + 
     -     Alfa/(16.D0*CW2*Pi*SW2)*
     -     ((dup2 - CBA**2*MHiggs_loop(k)**2)*
     -     DBLE(B0i(bb0,MHiggs_loop(k)**2,Mh02,MZ2))) + 
     -     Alfa/(16.D0*CW2*Pi*SW2)*
     -     ((-dup3 - MHiggs_loop(k)**2*SBA**2)*
     -     DBLE(B0i(bb0,MHiggs_loop(k)**2,MHH2,MZ2))) + 
     -     Alfa/(8.D0*CW2*Pi*SW2)*(CBA**2*MHiggs_loop(k)**2*
     -     DBLE(B0i(bb1,MHiggs_loop(k)**2,Mh02,MZ2))) + 
     -     Alfa/(8.D0*CW2*Pi*SW2)*(MHiggs_loop(k)**2*SBA**2*
     -     DBLE(B0i(bb1,MHiggs_loop(k)**2,MHH2,MZ2))) + 
     -     Alfa/(4.D0*Pi*SW2)*(MHiggs_loop(k)**2*
     -     DBLE(B0i(bb1,MHiggs_loop(k)**2,MHp2,MW2))) + 
     -     DBLE(A0(MW2))*((3*Alfa)/(8.D0*Pi*SW2) + 
     -     Alfa/(16.D0*CW2*Pi*SW2)*
     -     (-(SW2*C2B**2) + CW2*(1 + S2B**2))) + 
     -     DBLE(A0(MZ2))*(-(Alfa/(32.D0*CW2*Pi*SW2)*(1 - 3*S2B**2)) + 
     -     Alfa/(4.D0*CW2*Pi*SW2)*(1 - 1/4.D0*CBA**2 - 1/4.D0*SBA**2))

      do Cha4 = 1,2
         do Cha3 = 1,2
            dup112 = SB*UCha(Cha4,2)*VCha(Cha3,1) + 
     -           CB*UCha(Cha4,1)*VCha(Cha3,2)
            dup113 = SB*UCha(Cha3,2)*VCha(Cha4,1) + 
     -           CB*UCha(Cha3,1)*VCha(Cha4,2)
            dup114 = SB*UChaC(Cha4,2)*VChaC(Cha3,1) + 
     -           CB*UChaC(Cha4,1)*VChaC(Cha3,2)
            dup115 = SB*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     -           CB*UChaC(Cha3,1)*VChaC(Cha4,2)

            seA0 = seA0 - Alfa/(4.D0*Pi*SW2)*((dup112*dup114 + dup113*
     -           dup115)*DBLE(A0(MCha2(Cha4)))) - 
     -           Alfa/(4.D0*Pi*SW2)*((dup112*dup114 + dup113*dup115)*
     -           MHiggs_loop(k)**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -           MCha2(Cha3),MCha2(Cha4)))) + Alfa/(4.D0*Pi*SW2)*
     -           (DBLE(B0i(bb0,MHiggs_loop(k)**2,MCha2(Cha3),
     -           MCha2(Cha4)))*(dup114*dup115*MCha(Cha3)*MCha(Cha4) - 
     -           dup113*dup115*MCha2(Cha3) + 
     -           dup112*(dup113*MCha(Cha3)*MCha(Cha4) - 
     -           dup114*MCha2(Cha3))))
         enddo
      enddo

      do Neu4 = 1,4
         do Neu3 = 1,4
            dup116 = -(CW*ZNeu(Neu3,2)) + SW*ZNeu(Neu3,1)
            dup117 = -(CW*ZNeu(Neu4,2)) + SW*ZNeu(Neu4,1)
            dup118 = CB*(CW*ZNeu(Neu4,2) - SW*ZNeu(Neu4,1))
            dup119 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)
            dup120 = -(CW*ZNeuC(Neu3,2)) + SW*ZNeuC(Neu3,1)
            dup121 = -(CW*ZNeuC(Neu4,2)) + SW*ZNeuC(Neu4,1)
            dup122 = CB*(CW*ZNeuC(Neu4,2) - SW*ZNeuC(Neu4,1))
            dup123 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)

            seA0 = seA0 - Alfa/(8.D0*CW2*Pi*SW2)*(DBLE(A0(MNeu2(Neu4)))*
     -           (dup116*dup119 + dup117*SB*ZNeu(Neu3,3) + 
     -           dup118*ZNeu(Neu3,4))*
     -           (dup120*dup123 + dup121*SB*ZNeuC(Neu3,3) + 
     -           dup122*ZNeuC(Neu3,4))) - Alfa/(8.D0*CW2*Pi*SW2)*
     -           (MHiggs_loop(k)**2*DBLE(B0i(bb1,MHiggs_loop(k)**2,
     -           MNeu2(Neu3),MNeu2(Neu4)))*
     -           (dup116*dup119 + dup117*SB*ZNeu(Neu3,3) + 
     -           dup118*ZNeu(Neu3,4))*
     -           (dup120*dup123 + dup121*SB*ZNeuC(Neu3,3) + 
     -           dup122*ZNeuC(Neu3,4))) +DBLE(B0i(bb0,MHiggs_loop(k)**2,
     -           MNeu2(Neu3),MNeu2(Neu4)))*(-(Alfa/(8.D0*CW2*Pi*SW2)*
     -           (MNeu2(Neu3)*(dup116*dup119 + dup117*SB*ZNeu(Neu3,3) + 
     -           dup118*ZNeu(Neu3,4))*
     -           (dup120*dup123 + dup121*SB*ZNeuC(Neu3,3) + 
     -           dup122*ZNeuC(Neu3,4)))) + Alfa/(16.D0*CW2*Pi*SW2)*
     -           (MNeu(Neu3)*MNeu(Neu4)*
     -           ((dup116*dup119 + dup117*SB*ZNeu(Neu3,3) + 
     -           dup118*ZNeu(Neu3,4))**2 + 
     -           (dup120*dup123 + dup121*SB*ZNeuC(Neu3,3) + 
     -           dup122*ZNeuC(Neu3,4))**2)))
         enddo
      enddo


      end


********************************************************************************
