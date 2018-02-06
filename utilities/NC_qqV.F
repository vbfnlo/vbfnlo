*******************************************************************************
*******************************************************************************

**************************************************************************
** NC_qqV.F                                                             **
** sophy@particle.uni-karlsruhe.de                                      **
**                                                                      **
** This set of subroutines calculates the matrix squared elements from  **
** corrections to the qqV vertex for the neutral currents.              **
** Note that, by default, these subroutines are not used.  Instead we   **
** calculate an effective coupling for the qqV vertex, which            **
** incorporates both loop and counterterm diagrams.                     **
**************************************************************************

*******************************************************************************
*******************************************************************************

*** This subroutine calculates the counterterm part of the qqV loop component

      subroutine boxesQct(NC,mass2,gL,gR,q,dgL,dgR,S,T,T24,T14,S34,U,
     &                      Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared, left/right couplings and charge of external particles
      double precision mass2(5), gL(4), gR(4), q(4), dGL(4), dGR(4)

** Fermion field renormalisation constants for specific current NC.  
      double complex dZfL(4), dZfR(4)
      integer i


** Dummy parameters
      double complex dup19, dup20
      double complex dup17, dup18


** Current identifier
      integer NC  


** QED coupling
      double precision Alfa, Alfa2

** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Setting the quark field renormalisation constants for current NC
      dZfL(1) = dZfLu
      dZfL(2) = dZfLd
      dZfL(3) = dZfLc
      dZfL(4) = dZfLs

      dZfR(1) = dZfRu
      dZfR(2) = dZfRd
      dZfR(3) = dZfRc
      dZfR(4) = dZfRs


      if (NC .ge. 3) then
         dZfL(1) = dZfLd
         dZfL(2) = dZfLu
         dZfR(1) = dZfRd
         dZfR(2) = dZfRu
      end if

      if (NC .eq. 2 .or. NC .eq. 4) then
         dZfL(3) = dZfLs
         dZfL(4) = dZfLc
         dZfR(3) = dZfRs
         dZfR(4) = dZfRc
      end if
         

      dup19 = 2*dgL(3) + dZfL(3)*gL(3) + dZZZ1*gL(3) - dZAZ1*q(3)
      dup20 = 2*dgR(3) + dZfR(3)*gR(3) + dZZZ1*gR(3) - dZAZ1*q(3)

      Cloop(1) = Cloop(1) + 1/CW2*(Alfa*(Abb4*
     -     (-(4*MW*Pi*(dup19*EL*gL(1))) - 
     -     4*MW*Pi*(EL*gL(1)*gL(3)*DCONJG(dZfL(3))))) + 
     -     Alfa*(Abb2*(-(4*MW*Pi*(dup19*EL*gR(1))) - 
     -     4*MW*Pi*(EL*gL(3)*gR(1)*DCONJG(dZfL(3))))) + 
     -     Alfa*(Abb3*(-(4*MW*Pi*(dup20*EL*gL(1))) - 
     -     4*MW*Pi*(EL*gL(1)*gR(3)*DCONJG(dZfR(3))))) + 
     -     Alfa*(Abb1*(-(4*MW*Pi*(dup20*EL*gR(1))) - 
     -     4*MW*Pi*(EL*gR(1)*gR(3)*DCONJG(dZfR(3))))))/
     -     (SW*(-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U))

      dup17 = 2*dgL(1) + dZfL(1)*gL(1) + dZZZ1*gL(1) - dZAZ1*q(1)

      dup18 = 2*dgR(1) + dZfR(1)*gR(1) + dZZZ1*gR(1) - dZAZ1*q(1)

      Cloop(1) = Cloop(1) + 1/CW2*(Alfa*(Abb4*
     -     (-(4*MW*Pi*(dup17*EL*gL(3))) - 
     -     4*MW*Pi*(EL*gL(1)*gL(3)*DCONJG(dZfL(1))))) + 
     -     Alfa*(Abb3*(-(4*MW*Pi*(dup17*EL*gR(3))) - 
     -     4*MW*Pi*(EL*gL(1)*gR(3)*DCONJG(dZfL(1))))) + 
     -     Alfa*(Abb2*(-(4*MW*Pi*(dup18*EL*gL(3))) - 
     -     4*MW*Pi*(EL*gL(3)*gR(1)*DCONJG(dZfR(1))))) + 
     -     Alfa*(Abb1*(-(4*MW*Pi*(dup18*EL*gR(3))) - 
     -     4*MW*Pi*(EL*gR(1)*gR(3)*DCONJG(dZfR(1))))))/
     -     (SW*(-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U))      


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set P of loop corrections to the upper qqV
*** vertex for NC currents 13

      subroutine boxesQuP_13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup1, tmp1


** QED coupling
      double precision Alfa, Alfa2

** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup1 = 1/SW2*(gL(2) + gL(1)*(2*SW2*gL(1)**2 + 2*SW2*q(1)**2))

      tmp1 = (2*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*(gR(1)**2 + 
     -     q(1)**2))/SW + (2*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*
     -     (gR(1)**2 + q(1)**2))/SW + Cval(cc0,iintQ4)*(-(4*Alfa2*MW)/
     -     (CW2*SW2)*(AbbSumQ7*EL*gL(2)*gL(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (AbbSumQ5*EL*gL(2)*gR(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb4*EL*gL(2)*gL(3)*T)/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb3*EL*gL(2)*gR(3)*T)/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(3)*(gR(2)*mass2(2) - MW2*gL(2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb3*EL*gR(3)*(gR(2)*mass2(2) - MW2*gL(2)))/SW) + 
     -     1/(CW2*SW2)*(Bval(bb0,iintQ2)*
     -     (-(2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Bval(bb0,iintQ1)*
     -     (Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*gR(1)**2) - 4*MW*q(1)**2)) + 
     -     Alfa2*(Abb4*EL*gL(1)*
     -     (-(4*MW*(gL(1)**2*gL(3))) - 4*MW*(gL(3)*q(1)**2))) + 
     -     Alfa2*(Abb3*EL*gL(1)*
     -     (-(4*MW*(gL(1)**2*gR(3))) - 4*MW*(gR(3)*q(1)**2))) + 
     -     Alfa2*(Abb1*EL*gR(1)*
     -     (-(4*MW*(gR(1)**2*gR(3))) - 4*MW*(gR(3)*q(1)**2)))))/SW+
     -     1/CW2*(Cval(cc0,iintQ3)*
     -     (-(8*Alfa2*MW*(AbbSumQ7*EL*gL(1)*gL(3)*mass2(1)*q(1)**2)) - 
     -     4*Alfa2*MW*(AbbSumQ3*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(AbbSumQ5*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) - 
     -     8*Alfa2*MW*(AbbSumQ1*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     Alfa2*(Abb4*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(1)*q(1)**2) + 
     -     4*MW*(gL(1)*gL(3)*q(1)**2*T))) + 
     -     Alfa2*(Abb2*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)**2) + 
     -     4*MW*(gL(3)*gR(1)*q(1)**2*T))) + 
     -     Alfa2*(Abb3*EL*
     -     (4*MW*(gR(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     4*MW*(gL(1)*gR(3)*q(1)**2*T))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     4*MW*(gR(1)*gR(3)*q(1)**2*T)))))/SW

      tmp1 = tmp1 + Cval(cc0,iintQ5)*
     -     (-(8*Alfa2*MW)/CW2*(AbbSumQ7*EL*gL(1)**3*gL(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumQ3*EL*gL(3)*gR(1)**3*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumQ5*EL*gL(1)**3*gR(3)*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(AbbSumQ1*EL*gR(1)**3*gR(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*T)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*T)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*T)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*T)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)*(gR(1)*mass2(1) - MZ2*gL(1)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)*(gR(1)*mass2(1) - MZ2*gL(1)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)**2*(-(gL(1)*mass2(1)) + MZ2*gR(1)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)*(-(gL(1)*mass2(1)) + MZ2*gR(1)))/SW)+
     -     1/(CW2*SW2)*(Cval(cc00,iintQ4)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ5)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ3)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2)))/SW

      tmp1 = tmp1 + 1/(CW2*SW2)*(Cval(cc1,iintQ4)*
     -     (-6*Alfa2*MW*(AbbSumQ17*EL*gL(2)*gL(3)*mass2(1)) - 
     -     2*Alfa2*MW*(AbbSumQ15*EL*gL(2)*gR(3)*mass2(1)) + 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*T) + 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc1,iintQ5)*
     -     (-12*Alfa2*MW*(AbbSumQ17*EL*gL(1)**3*gL(3)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)**3*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)**3*gR(3)*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)**3*gR(3)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)*T) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3*T) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)*T) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc1,iintQ3)*
     -     (-12*Alfa2*MW*(AbbSumQ17*EL*gL(1)*gL(3)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) - 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2*T)))/SW

      tmp1 = tmp1 + 1/CW2*(Cval(cc2,iintQ3)*
     -     (-12*Alfa2*MW*(AbbSumQ17*EL*gL(1)*gL(3)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) - 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2*T) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2*T)))/SW

      tmp1 = tmp1 + (Alfa2*MW)/CW2*(Abb4*dup1*EL*gL(3))/SW + 
     -     (Alfa2*MW)/CW2*(Abb3*dup1*EL*gR(3))/SW + 
     -     1/(CW2*SW2)*(Cval(cc2,iintQ4)*
     -     (-6*Alfa2*MW*(AbbSumQ17*EL*gL(2)*gL(3)*mass2(1)) - 
     -     2*Alfa2*MW*(AbbSumQ15*EL*gL(2)*gR(3)*mass2(1)) + 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*T) + 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc2,iintQ5)*
     -     (-12*Alfa2*MW*(AbbSumQ17*EL*gL(1)**3*gL(3)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)**3*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)**3*gR(3)*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)**3*gR(3)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)*T) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3*T) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)*T) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3)*T)))/SW

      Cloop(1) = tmp1/((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U)) + Cloop(1)
        
        
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set N of loop corrections to the upper qqV
*** vertex for NC currents 13 

      subroutine boxesQuN_13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup2, dup3


** QED coupling
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup2 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)

      dup3 = 4*mass2(1) + 3*mass2(3) + 2*MH2


      Cloop(1) = Cloop(1) +((Alfa2*CW*MW)/(CW2*SW2**2)*(Abb4*EL*gL(3)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*(Abb3*EL*gR(3)) + 
     -     1/(CW2*SW2**2)*(Bval(bb0,iintQ6)*
     -     (-(2*Alfa2*CW*MW*(Abb4*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 1/(CW2*SW2**2)*
     -     (Cval(cc00,iintQ7)*(-(4*Alfa2*CW*MW*(Abb4*EL*gL(3))) - 
     -     4*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 1/(CW2*SW2**2)*
     -     (Cval(cc2,iintQ7)*
     -     (5*Alfa2*CW*MW*(AbbSumQ26*EL*gL(3)*mass2(1)) + 
     -     3*Alfa2*CW*MW*(AbbSumQ25*EL*gR(3)*mass2(1)) - 
     -     2*Alfa2*CW*MW*(Abb4*EL*gL(3)*T) - 
     -     2*Alfa2*CW*MW*(Abb3*EL*gR(3)*T))) + 
     -     1/(CW2*SW2**2)*(Cval(cc1,iintQ7)*(Alfa2*CW*(Abb4*
     -     (8*MW*(dup2*EL*gL(3)) - 4*MW*(dup3*EL*gL(3)))) + 
     -     Alfa2*CW*(Abb3*
     -     (8*MW*(dup2*EL*gR(3)) - 4*MW*(dup3*EL*gR(3)))))) + 
     -     1/(CW2*SW2**2)*(Cval(cc0,iintQ7)*
     -     (5*Alfa2*CW*MW*(AbbSumQ26*EL*gL(3)*mass2(1)) + 
     -     3*Alfa2*CW*MW*(AbbSumQ25*EL*gR(3)*mass2(1)) + 
     -     Alfa2*CW*(Abb4*EL*gL(3)*(-(2*MW*mass2(2)) - 2*MW*T)) + 
     -     Alfa2*CW*(Abb3*EL*gR(3)*(-(2*MW*mass2(2)) - 2*MW*T)))))/
     -     ((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - S - 
     -     T24 - U))
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set P of loop corrections to the upper qqV
*** vertex for NC currents 24 

      subroutine boxesQuP_24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup1, tmp1


** QED coupling
      double precision Alfa, Alfa2

** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup1 = 1/SW2*(gL(2) + gL(1)*(2*SW2*gL(1)**2 + 2*SW2*q(1)**2))

      tmp1 = -((2*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*(gR(1)**2 + 
     -     q(1)**2))/SW) - (2*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*
     -     (gR(1)**2 + q(1)**2))/SW + 
     -     Cval(cc0,iintQ4)*((4*Alfa2*MW)/(CW2*SW2)*
     -     (AbbSumQ7*EL*gL(2)*gL(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (AbbSumQ5*EL*gL(2)*gR(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb4*EL*gL(2)*gL(3)*T)/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb3*EL*gL(2)*gR(3)*T)/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(3)*(gR(2)*mass2(2) - MW2*gL(2)))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb3*EL*gR(3)*(gR(2)*mass2(2) - MW2*gL(2)))/SW) + 
     -     1/(CW2*SW2)*(Bval(bb0,iintQ2)*
     -     (2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Bval(bb0,iintQ1)*
     -     (Alfa2*(Abb2*EL*gL(3)*gR(1)*(4*MW*gR(1)**2 + 4*MW*q(1)**2)) + 
     -     Alfa2*(Abb4*EL*gL(1)*
     -     (4*MW*(gL(1)**2*gL(3)) + 4*MW*(gL(3)*q(1)**2))) + 
     -     Alfa2*(Abb3*EL*gL(1)*
     -     (4*MW*(gL(1)**2*gR(3)) + 4*MW*(gR(3)*q(1)**2))) + 
     -     Alfa2*(Abb1*EL*gR(1)*
     -     (4*MW*(gR(1)**2*gR(3)) + 4*MW*(gR(3)*q(1)**2)))))/SW + 
     -     1/CW2*(Cval(cc0,iintQ3)*
     -     (8*Alfa2*MW*(AbbSumQ7*EL*gL(1)*gL(3)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ3*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ5*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     8*Alfa2*MW*(AbbSumQ1*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     Alfa2*(Abb4*EL*
     -     (-(4*MW*(gL(3)*gR(1)*mass2(1)*q(1)**2)) - 
     -     4*MW*(gL(1)*gL(3)*q(1)**2*T))) + 
     -     Alfa2*(Abb2*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)**2)) - 
     -     4*MW*(gL(3)*gR(1)*q(1)**2*T))) + 
     -     Alfa2*(Abb3*EL*
     -     (-(4*MW*(gR(1)*gR(3)*mass2(1)*q(1)**2)) - 
     -     4*MW*(gL(1)*gR(3)*q(1)**2*T))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)**2)) - 
     -     4*MW*(gR(1)*gR(3)*q(1)**2*T)))))/SW

      tmp1 = tmp1 + Cval(cc0,iintQ5)*
     -     ((8*Alfa2*MW)/CW2*(AbbSumQ7*EL*gL(1)**3*gL(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumQ3*EL*gL(3)*gR(1)**3*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumQ5*EL*gL(1)**3*gR(3)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumQ1*EL*gR(1)**3*gR(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*T)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*T)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*T)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*T)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)*(gR(1)*mass2(1) - MZ2*gL(1)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)*(gR(1)*mass2(1) - MZ2*gL(1)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)**2*(-(gL(1)*mass2(1)) + MZ2*gR(1)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)*(-(gL(1)*mass2(1)) + MZ2*gR(1)))/SW)+
     -     1/(CW2*SW2)*(Cval(cc00,iintQ4)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(2)*gL(3))) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ5)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ3)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2)))/SW

      tmp1 = tmp1 + 1/(CW2*SW2)*(Cval(cc1,iintQ4)*
     -     (6*Alfa2*MW*(AbbSumQ17*EL*gL(2)*gL(3)*mass2(1)) + 
     -     2*Alfa2*MW*(AbbSumQ15*EL*gL(2)*gR(3)*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*T) - 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc1,iintQ5)*
     -     (12*Alfa2*MW*(AbbSumQ17*EL*gL(1)**3*gL(3)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)**3*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)**3*gR(3)*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)**3*gR(3)*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc1,iintQ3)*
     -     (12*Alfa2*MW*(AbbSumQ17*EL*gL(1)*gL(3)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2*T)))/SW

      tmp1 = tmp1 + 1/CW2*(Cval(cc2,iintQ3)*
     -     (12*Alfa2*MW*(AbbSumQ17*EL*gL(1)*gL(3)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)*mass2(1)*q(1)**2) + 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)*gR(3)*mass2(1)*q(1)**2) + 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)*gR(3)*mass2(1)*q(1)**2) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)**2*T)))/SW

      tmp1 = tmp1 - (Alfa2*MW)/CW2*(Abb4*dup1*EL*gL(3))/SW - 
     -     (Alfa2*MW)/CW2*(Abb3*dup1*EL*gR(3))/SW + 
     -     1/(CW2*SW2)*(Cval(cc2,iintQ4)*
     -     (6*Alfa2*MW*(AbbSumQ17*EL*gL(2)*gL(3)*mass2(1)) + 
     -     2*Alfa2*MW*(AbbSumQ15*EL*gL(2)*gR(3)*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*T) - 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*T)))/SW + 
     -     1/CW2*(Cval(cc2,iintQ5)*
     -     (12*Alfa2*MW*(AbbSumQ17*EL*gL(1)**3*gL(3)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumQ12*EL*gL(3)*gR(1)**3*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumQ15*EL*gL(1)**3*gR(3)*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumQ9*EL*gR(1)**3*gR(3)*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3)*T)))/SW

      Cloop(1) = tmp1/((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U)) + Cloop(1)
        
        
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set N of loop corrections to the upper qqV
*** vertex for NC currents 24

      subroutine boxesQuN_24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup2, dup3


** QED coupling
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup2 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup3 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + (-((Alfa2*CW*MW)/(CW2*SW2**2)*(Abb4*EL*
     -     gL(3))) - (Alfa2*CW*MW)/(CW2*SW2**2)*(Abb3*EL*gR(3)) + 
     -     1/(CW2*SW2**2)*(Bval(bb0,iintQ6)*
     -     (2*Alfa2*CW*MW*(Abb4*EL*gL(3)) + 
     -     2*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 
     -     1/(CW2*SW2**2)*(Cval(cc00,iintQ7)*
     -     (4*Alfa2*CW*MW*(Abb4*EL*gL(3)) + 
     -     4*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 
     -     1/(CW2*SW2**2)*(Cval(cc2,iintQ7)*
     -     (-5*Alfa2*CW*MW*(AbbSumQ26*EL*gL(3)*mass2(1)) - 
     -     3*Alfa2*CW*MW*(AbbSumQ25*EL*gR(3)*mass2(1)) + 
     -     2*Alfa2*CW*MW*(Abb4*EL*gL(3)*T) + 
     -     2*Alfa2*CW*MW*(Abb3*EL*gR(3)*T))) + 1/(CW2*SW2**2)*
     -     (Cval(cc1,iintQ7)*(Alfa2*CW*(Abb4*
     -     (-(8*MW*(dup2*EL*gL(3))) + 4*MW*(dup3*EL*gL(3))))+
     -     Alfa2*CW*(Abb3*(-(8*MW*(dup2*EL*gR(3))) + 
     -     4*MW*(dup3*EL*gR(3)))))) + 
     -     1/(CW2*SW2**2)*(Cval(cc0,iintQ7)*
     -     (-5*Alfa2*CW*MW*(AbbSumQ26*EL*gL(3)*mass2(1)) - 
     -     3*Alfa2*CW*MW*(AbbSumQ25*EL*gR(3)*mass2(1)) + 
     -     Alfa2*CW*(Abb4*EL*gL(3)*(2*MW*mass2(2) + 2*MW*T)) + 
     -     Alfa2*CW*(Abb3*EL*gR(3)*(2*MW*mass2(2) + 2*MW*T)))))
     -     /((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - S - 
     -     T24 - U))

      
      end


*******************************************************************************
*******************************************************************************
 
*** This subroutine calculates the set P of loop corrections to the lower qqV
*** vertex for NC currents 12

      subroutine boxesQlP_12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)
      
      implicit none
      
      
** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup10, dup11, dup12, dup13, dup14, dup4, dup5
      double complex dup6, dup7, dup8, dup9, tmp2

** QED coupling
      double precision Alfa, Alfa2

** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup4 = -(gL(3)*gR(1)) + gL(1)*(gL(3) - gR(3)) + gR(1)*gR(3)
      dup5 = gR(3)*mass2(3) + gL(3)*(-2*mass2(1) - MH2)
      dup6 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + MH2)
      dup7 = gR(3)**3*mass2(3) + gL(3)**3*(-2*mass2(1) - MH2)
      dup8 = -(gL(3)*mass2(3)) + gR(3)*(2*mass2(1) + MH2)
      dup9 = gL(3)*mass2(3) + gR(3)*(2*mass2(1) + MH2)
      dup10 = -(gL(3)**3*mass2(3)) + gR(3)**3*(2*mass2(1) + MH2)
      dup11 = gR(4)*mass2(4) + gL(4)*(2*mass2(1) + MH2 - MW2)
      dup12 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + MH2 - MZ2)
      dup13 = -(gL(3)*mass2(3)) + gR(3)*(-2*mass2(1) - MH2 + MZ2)
      dup14 = 1/SW2*(gL(4) + gL(3)*(2*SW2*gL(3)**2 + 2*SW2*q(3)**2))

      tmp2 = (2*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*(gR(3)**2 + 
     -     q(3)**2))/SW + (2*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*
     -     (gR(3)**2 + q(3)**2))/SW + Cval(cc1,iintQ11)*((2*Alfa2*MW)/
     -     (CW2*SW2)*(Abb3*EL*gL(1)*gL(4)*mass2(3))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*gR(1)*mass2(3))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(S + T24 + U))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(S + T24 + U))/SW) + 
     -     1/(CW2*SW2)*(Bval(bb0,iintQ9)*
     -     (-(2*Alfa2*MW*(Abb4*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/(CW2*SW2)*(Cval(cc00,iintQ11)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(1)*gL(4)) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ12)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/CW2*(Cval(cc00,iintQ10)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(3)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(3)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(3)**2) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(3)**2)))/SW + 
     -     1/CW2*(Bval(bb0,iintQ8)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*gL(3)**2) - 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*gR(3)**2) - 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*gR(3)**2) - 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb2*EL*gL(3)*
     -     (-(4*MW*(gL(3)**2*gR(1))) - 4*MW*(gR(1)*q(3)**2)))))/SW

      tmp2 = tmp2 + Cval(cc2,iintQ11)*
     -     (-(2*Alfa2*MW)/(CW2*SW2)*(Abb3*EL*gL(1)*gL(4)*mass2(3))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*gR(1)*mass2(3))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(S + T24 + U))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(S + T24 + U))/SW) 

      tmp2 = tmp2 + 1/CW2*(Cval(cc2,iintQ10)*(Alfa2*(Abb4*
     -     (-(4*MW*(dup5*EL*gL(1)*q(3)**2)) + 
     -     MW*(EL*gL(1)*gL(3)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb2*(-(4*MW*(dup5*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb3*(4*MW*(dup8*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb1*(4*MW*(dup8*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U))))))/SW

      tmp2 = tmp2 + 1/CW2*(Cval(cc1,iintQ10)*(Alfa2*(Abb4*(-(4*MW*(dup5*
     -     EL*gL(1)*q(3)**2)) + MW*(EL*gL(1)*gL(3)*q(3)**2*(-4*S - 
     -     4*T24 - 4*U)))) + Alfa2*(Abb2*
     -     (-(4*MW*(dup5*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb3*(4*MW*(dup8*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb1*(4*MW*(dup8*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U))))))/
     -     SW + 1/CW2*(Cval(cc0,iintQ10)*
     -     (Alfa2*(AbbSumQ40*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(3)**2) - 
     -     4*MW*(gR(3)*mass2(3)*q(3)**2))) + 
     -     Alfa2*(AbbSumQ35*EL*
     -     (-(4*MW*(gL(3)*gR(1)*mass2(3)*q(3)**2)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(3)**2))) + Alfa2*(Abb4*
     -     (4*MW*(dup6*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gL(3)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb2*(4*MW*(dup6*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb3*(4*MW*(dup9*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb1*(4*MW*(dup9*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(-4*S - 4*T24 - 4*U))))))/SW

      tmp2 = tmp2 + Cval(cc2,iintQ12)*
     -     ((4*Alfa2*MW)/CW2*(Abb3*dup10*EL*gL(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*dup7*EL*gL(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup10*EL*gR(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup7*EL*gR(1))/SW - (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)**3*(S + T24 + U))/SW - (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**3*gR(1)*(S + T24 + U))/SW - (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)**3*(S + T24 + U))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*(S + T24 + U))/SW)
     -     + 1/CW2*(Cval(cc1,iintQ12)*(Alfa2*(Abb4*
     -     (-(4*MW*(dup7*EL*gL(1))) + 
     -     MW*(EL*gL(1)*gL(3)**3*(-4*S - 4*T24 - 4*U)))) + Alfa2*(Abb2*
     -     (-(4*MW*(dup7*EL*gR(1))) + 
     -     MW*(EL*gL(3)**3*gR(1)*(-4*S - 4*T24 - 4*U)))) + Alfa2*(Abb3*
     -     (4*MW*(dup10*EL*gL(1)) + 
     -     MW*(EL*gL(1)*gR(3)**3*(-4*S - 4*T24 - 4*U)))) + Alfa2*(Abb1*
     -     (4*MW*(dup10*EL*gR(1)) + 
     -     MW*(EL*gR(1)*gR(3)**3*(-4*S - 4*T24 - 4*U))))))/SW

      tmp2 = tmp2 + (Alfa2*MW)/CW2*(Abb4*dup14*EL*gL(1))/SW + 
     -     (Alfa2*MW)/CW2*(Abb2*dup14*EL*gR(1))/SW + 
     -     Cval(cc0,iintQ12)*((4*Alfa2*MW)/CW2*
     -     (Abb4*dup12*EL*gL(1)*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup12*EL*gL(3)**2*gR(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup13*EL*gL(1)*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup13*EL*gR(1)*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumQ40*EL*gL(1)*(-gL(3)**3 + gR(3)**3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumQ35*EL*gR(1)*(-gL(3)**3 + gR(3)**3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)**3*(S + T24 + U))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**3*gR(1)*(S + T24 + U))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)**3*(S + T24 + U))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*(S + T24 + U))/SW)
     -     + 1/(CW2*SW2)*(Cval(cc0,iintQ11)*
     -     (2*Alfa2*MW*(AbbSumQ40*EL*gL(1)*gL(4)*mass2(3)) - 
     -     2*Alfa2*MW*(AbbSumQ35*EL*gL(4)*gR(1)*mass2(3)) + 
     -     Alfa2*(Abb4*(2*MW*(dup11*EL*gL(1)) + 
     -     MW*(EL*gL(1)*gL(4)*(-2*S - 2*T24 - 2*U)))) + 
     -     Alfa2*(Abb2*(2*MW*(dup11*EL*gR(1)) + 
     -     MW*(EL*gL(4)*gR(1)*(-2*S - 2*T24 - 2*U))))))/SW

      Cloop(1) = tmp2/((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U)) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set N of loop corrections to the lower qqV
*** vertex for NC currents 12

      subroutine boxesQlN_12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup15, dup16


** QED coupling
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup15 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup16 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = ((Alfa2*CW*MW)/(CW2*SW2**2)*(Abb4*EL*gL(1)) + (Alfa2*
     -     CW*MW)/(CW2*SW2**2)*(Abb2*EL*gR(1)) + Cval(cc2,iintQ14)*
     -     (-(Alfa2*CW*MW)/(CW2*SW2**2)*(AbbSumQ61*EL*gL(1)*mass2(3)) + 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*(AbbSumQ60*EL*gR(1)*mass2(3)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*(Abb4*EL*gL(1)*(2*mass2(1) + 
     -     MH2)) - (2*Alfa2*CW*MW)/(CW2*SW2**2)*(Abb2*EL*gR(1)*(2*
     -     mass2(1) + MH2)) + (2*Alfa2*CW*MW)/(CW2*SW2**2)* (Abb4*EL*
     -     gL(1)*(S + T24 + U)) + (2*Alfa2*CW*MW)/(CW2*SW2**2)*(Abb2*EL*
     -     gR(1)*(S + T24 + U))) + Cval(cc0,iintQ14)*(-(Alfa2*CW*MW)/
     -     (CW2*SW2**2)*(AbbSumQ61*EL*gL(1)*mass2(3)) + (3*Alfa2*CW*MW)/
     -     (CW2*SW2**2)*(AbbSumQ60*EL*gR(1)*mass2(3)) - (2*Alfa2*CW*MW)/
     -     (CW2*SW2**2)*(Abb4*EL*gL(1)*(2*mass2(1) + mass2(4) + MH2)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(2*mass2(1) + mass2(4) + MH2)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb4*EL*gL(1)*(S + T24 + U)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(S + T24 + U))) + 
     -     1/(CW2*SW2**2)*
     -     (Bval(bb0,iintQ13)*
     -     (-(2*Alfa2*CW*MW*(Abb4*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(Abb2*EL*gR(1)))) + 
     -     1/(CW2*SW2**2)*
     -     (Cval(cc00,iintQ14)*
     -     (-(4*Alfa2*CW*MW*(Abb4*EL*gL(1))) - 
     -     4*Alfa2*CW*MW*(Abb2*EL*gR(1)))) + 
     -     1/(CW2*SW2**2)*
     -     (Cval(cc1,iintQ14)*
     -     (Alfa2*CW*(Abb4*
     -     (8*MW*(dup15*EL*gL(1)) - 4*MW*(dup16*EL*gL(1)))) + 
     -     Alfa2*CW*(Abb2*
     -     (8*MW*(dup15*EL*gR(1)) - 4*MW*(dup16*EL*gR(1)))))))/
     -     ((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - S - 
     -     T24 - U))
      
      
      end


*******************************************************************************
*******************************************************************************
      
*** This subroutine calculates the set P of loop corrections to the lower qqV
*** vertex for NC currents 34

      subroutine boxesQlP_34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)
      
      implicit none
      
      
** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup10, dup11, dup12, dup13, dup14, dup4, dup5
      double complex dup6, dup7, dup8, dup9, tmp2

** QED coupling
      double precision Alfa, Alfa2

** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup4 = -(gL(3)*gR(1)) + gL(1)*(gL(3) - gR(3)) + gR(1)*gR(3)
      dup5 = gR(3)*mass2(3) + gL(3)*(-2*mass2(1) - MH2)
      dup6 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + MH2)
      dup7 = gR(3)**3*mass2(3) + gL(3)**3*(-2*mass2(1) - MH2)
      dup8 = -(gL(3)*mass2(3)) + gR(3)*(2*mass2(1) + MH2)
      dup9 = gL(3)*mass2(3) + gR(3)*(2*mass2(1) + MH2)
      dup10 = -(gL(3)**3*mass2(3)) + gR(3)**3*(2*mass2(1) + MH2)
      dup11 = gR(4)*mass2(4) + gL(4)*(2*mass2(1) + MH2 - MW2)
      dup12 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + MH2 - MZ2)
      dup13 = -(gL(3)*mass2(3)) + gR(3)*(-2*mass2(1) - MH2 + MZ2)
      dup14 = 1/SW2*(gL(4) + gL(3)*(2*SW2*gL(3)**2 + 2*SW2*q(3)**2))

      tmp2 = -((2*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*(gR(3)**2 + 
     -     q(3)**2))/SW) - (2*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*
     -     (gR(3)**2 + q(3)**2))/SW + Cval(cc1,iintQ11)*((2*Alfa2*MW)/
     -     (CW2*SW2)*(Abb3*EL*gL(1)*gL(4)*mass2(3))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*gR(1)*mass2(3))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(S + T24 + U))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(S + T24 + U))/SW) + 
     -     1/(CW2*SW2)*(Bval(bb0,iintQ9)*
     -     (2*Alfa2*MW*(Abb4*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/(CW2*SW2)*(Cval(cc00,iintQ11)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(1)*gL(4))) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Cval(cc00,iintQ12)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/CW2*(Cval(cc00,iintQ10)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(3)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(3)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(3)**2) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(3)**2)))/SW + 
     -     1/CW2*(Bval(bb0,iintQ8)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*(4*MW*gL(3)**2 + 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*gR(3)**2 + 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*gR(3)**2 + 4*MW*q(3)**2)) + 
     -     Alfa2*(Abb2*EL*gL(3)*
     -     (4*MW*(gL(3)**2*gR(1)) + 4*MW*(gR(1)*q(3)**2)))))/SW

      tmp2 = tmp2 + Cval(cc2,iintQ11)*
     -     ((2*Alfa2*MW)/(CW2*SW2)*(Abb3*EL*gL(1)*gL(4)*mass2(3))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*gR(1)*mass2(3))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb4*EL*gL(1)*gL(4)*(S + T24 + U))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb2*EL*gL(4)*gR(1)*(S + T24 + U))/SW) 

      tmp2 = tmp2 + 1/CW2*(Cval(cc2,iintQ10)*(Alfa2*(Abb4*
     -     (4*MW*(dup5*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gL(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb2*
     -     (4*MW*(dup5*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb3*
     -     (-(4*MW*(dup8*EL*gL(1)*q(3)**2)) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb1*
     -     (-(4*MW*(dup8*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U))))))/
     -     SW

      tmp2 = tmp2 + 1/CW2*(Cval(cc1,iintQ10)*(Alfa2*(Abb4*
     -     (4*MW*(dup5*EL*gL(1)*q(3)**2) + 
     -     MW*(EL*gL(1)*gL(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb2*
     -     (4*MW*(dup5*EL*gR(1)*q(3)**2) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb3*
     -     (-(4*MW*(dup8*EL*gL(1)*q(3)**2)) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb1*
     -     (-(4*MW*(dup8*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U))))))/
     -     SW + 1/CW2*(Cval(cc0,iintQ10)*
     -     (Alfa2*(AbbSumQ40*EL*gL(1)*
     -     (-(4*MW*(gL(3)*mass2(3)*q(3)**2)) + 
     -     4*MW*(gR(3)*mass2(3)*q(3)**2))) + 
     -     Alfa2*(AbbSumQ35*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(3)**2) - 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(3)**2))) + 
     -     Alfa2*(Abb4*
     -     (-(4*MW*(dup6*EL*gL(1)*q(3)**2)) + 
     -     MW*(EL*gL(1)*gL(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb2*
     -     (-(4*MW*(dup6*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gL(3)*gR(1)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb3*
     -     (-(4*MW*(dup9*EL*gL(1)*q(3)**2)) + 
     -     MW*(EL*gL(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb1*
     -     (-(4*MW*(dup9*EL*gR(1)*q(3)**2)) + 
     -     MW*(EL*gR(1)*gR(3)*q(3)**2*(4*S + 4*T24 + 4*U))))))/
     -     SW

      tmp2 = tmp2 + Cval(cc2,iintQ12)*
     -     (-((4*Alfa2*MW)/CW2*(Abb3*dup10*EL*gL(1))/SW) + 
     -     (4*Alfa2*MW)/CW2*(Abb4*dup7*EL*gL(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup10*EL*gR(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup7*EL*gR(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)**3*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**3*gR(1)*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)**3*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*(S + T24 + U))/SW)
     -     + 1/CW2*(Cval(cc1,iintQ12)*
     -     (Alfa2*(Abb4*
     -     (4*MW*(dup7*EL*gL(1)) + 
     -     MW*(EL*gL(1)*gL(3)**3*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb2*
     -     (4*MW*(dup7*EL*gR(1)) + 
     -     MW*(EL*gL(3)**3*gR(1)*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb3*
     -     (-(4*MW*(dup10*EL*gL(1))) + 
     -     MW*(EL*gL(1)*gR(3)**3*(4*S + 4*T24 + 4*U)))) + 
     -     Alfa2*(Abb1*
     -     (-(4*MW*(dup10*EL*gR(1))) + 
     -     MW*(EL*gR(1)*gR(3)**3*(4*S + 4*T24 + 4*U))))))/SW

      tmp2 = tmp2 - (Alfa2*MW)/CW2*(Abb4*dup14*EL*gL(1))/SW - 
     -     (Alfa2*MW)/CW2*(Abb2*dup14*EL*gR(1))/SW + 
     -     Cval(cc0,iintQ12)*(-((4*Alfa2*MW)/CW2*
     -     (Abb4*dup12*EL*gL(1)*gL(3)**2)/SW) - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup12*EL*gL(3)**2*gR(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup13*EL*gL(1)*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup13*EL*gR(1)*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumQ40*EL*gL(1)*(-gL(3)**3 + gR(3)**3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumQ35*EL*gR(1)*(-gL(3)**3 + gR(3)**3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)**3*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**3*gR(1)*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)**3*(S + T24 + U))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*(S + T24 + U))/SW)
     -     + 1/(CW2*SW2)*(Cval(cc0,iintQ11)*
     -     (-2*Alfa2*MW*(AbbSumQ40*EL*gL(1)*gL(4)*mass2(3)) + 
     -     2*Alfa2*MW*(AbbSumQ35*EL*gL(4)*gR(1)*mass2(3)) + 
     -     Alfa2*(Abb4*
     -     (-(2*MW*(dup11*EL*gL(1))) + 
     -     MW*(EL*gL(1)*gL(4)*(2*S + 2*T24 + 2*U)))) + 
     -     Alfa2*(Abb2*
     -     (-(2*MW*(dup11*EL*gR(1))) + 
     -     MW*(EL*gL(4)*gR(1)*(2*S + 2*T24 + 2*U))))))/SW

      Cloop(1) = tmp2/((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - 
     -     S - T24 - U)) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the set N of loop corrections to the lower qqV
*** vertex for NC currents 34

      subroutine boxesQlN_34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variables
      double complex dup15, dup16


** QED coupling
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup15 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup16 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + (-((Alfa2*CW*MW)/(CW2*SW2**2)*(Abb4*EL*
     -     gL(1))) - (Alfa2*CW*MW)/(CW2*SW2**2)*(Abb2*EL*gR(1)) + 
     -     Cval(cc2,iintQ14)*((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (AbbSumQ61*EL*gL(1)*mass2(3)) - 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (AbbSumQ60*EL*gR(1)*mass2(3)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb4*EL*gL(1)*(2*mass2(1) + MH2)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(2*mass2(1) + MH2)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb4*EL*gL(1)*(S + T24 + U)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(S + T24 + U))) + 
     -     Cval(cc0,iintQ14)*((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (AbbSumQ61*EL*gL(1)*mass2(3)) - 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (AbbSumQ60*EL*gR(1)*mass2(3)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb4*EL*gL(1)*(2*mass2(1) + mass2(4) + MH2)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(2*mass2(1) + mass2(4) + MH2)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb4*EL*gL(1)*(S + T24 + U)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb2*EL*gR(1)*(S + T24 + U))) + 
     -     1/(CW2*SW2**2)*
     -     (Bval(bb0,iintQ13)*
     -     (2*Alfa2*CW*MW*(Abb4*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(Abb2*EL*gR(1)))) + 
     -     1/(CW2*SW2**2)*
     -     (Cval(cc00,iintQ14)*
     -     (4*Alfa2*CW*MW*(Abb4*EL*gL(1)) + 
     -     4*Alfa2*CW*MW*(Abb2*EL*gR(1)))) + 
     -     1/(CW2*SW2**2)*(Cval(cc1,iintQ14)*(Alfa2*CW*(Abb4*
     -     (-(8*MW*(dup15*EL*gL(1))) + 4*MW*(dup16*EL*gL(1))))
     -     + Alfa2*CW*
     -     (Abb2*(-(8*MW*(dup15*EL*gR(1))) + 
     -     4*MW*(dup16*EL*gR(1)))))))/
     -     ((-MZ2 + T)*(2*mass2(1) + 3*mass2(3) + MH2 - MZ2 - S - 
     -     T24 - U))
      
      
      end


********************************************************************************
