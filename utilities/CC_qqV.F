*******************************************************************************
*******************************************************************************

**************************************************************************
** CC_qqV.F                                                             **
** sophy@particle.uni-karlsruhe.de                                      **
**                                                                      **
** This set of subroutines calculates the matrix squared elements from  **
** corrections to the qqV vertex for the charged currents.              **
** Note that, by default, these subroutines are not used.  Instead we   **
** calculate an effective coupling for the qqV vertex, which            **
** incorporates both loop and counterterm diagrams.                     **
**************************************************************************

*******************************************************************************
*******************************************************************************

*** This subroutine calculates the counterterm part of the qqV loop component

      subroutine qqVCT(qk,CC,mass2,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared, left/right couplings and charge of external particles
      double precision mass2(5)

** Fermion field renormalisation constants for specific current CC.  Labelling:
** Particle number i
      double complex dZfL(4), dZfR(4)
      integer i

** Current and quark identifier
      integer CC, qk  


      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** Setting the quark field renormalisation constants for current CC
      dZfL(1) = dZfLu
      dZfL(2) = dZfLd
      dZfL(3) = dZfLc
      dZfL(4) = dZfLs

      if (CC .eq. 1) then
         dZfL(3) = dZfLs
         dZfL(4) = dZfLc
      else if (CC .eq. 2) then
         dZfL(1) = dZfLd
         dZfL(2) = dZfLu
      end if  


** upper line diagram
* up 51
      Cloop(1) = (-((2*Alfa*MW*Pi)/SW2*(Abb1*EL*dZfL(1))/SW) - 
     -     (2*Alfa*MW*Pi)/SW2**2*
     -     (Abb1*EL*(-2*dSW1*SW + 2*SW2*dZet1 + SW2*dZW1))/SW - 
     -     (4*Alfa*MW*Pi)/SW2*
* down 51
     -     (Abb1*EL*(1/2.D0*DCONJG(dZfL(2))))/SW)/
     -     ((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U))
      

** lower line diagram
      Cloop(1) = Cloop(1) + 
     -     (-((2*Alfa*MW*Pi)/SW2*
* charm 51
     -     (Abb1*EL*DCONJG(dZfL(4)))/SW) - 
* strange 51
     -     (2*Alfa*MW*Pi)/SW2*(Abb1*EL*dZfL(3))/SW - 
     -     (2*Alfa*MW*Pi)/SW2**2*
     -     (Abb1*EL*(-2*dSW1*SW + 
     -     SW2*(2*dZet1 + dZW1)))/SW)/
     -     ((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U))


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the upper quark line
*** which have a relative minus sign between currents 5 and 6, for quark 
*** patterns 1 and 3

      subroutine qqVUP_n13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

      
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** set A
      Cloop(1) = (Cval(cc1,iintQu6)*((Alfa2*CW*MW)/SW2**2*(Abb1*EL*
     -     gL(1)*(-mass2(1) + 
     -     mass2(2)))) + Cval(cc1,iintQu7)*(-(Alfa2*MW)/SW2*(Abb1*EL*
     -     (-mass2(1) + mass2(2))*q(1))/SW) + Cval(cc0,iintQu7)*
     -     ((2*Alfa2*
     -     MW)/SW2*(Abb1*EL*(-mass2(1) + mass2(2))*q(1))/SW - (2*Alfa2*
     -     MW)/SW2*(Abb1*EL*q(1)*T)/SW) -(2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*
     -     Bval(bb0,iintQu4))/SW -(2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*
     -     Bval(bb0,iintQu5)) - (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*
     -     Cval(cc00,iintQu6)) - (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*
     -     Cval(cc00,iintQu7))/SW + Cval(cc0,iintQu6)*
     -     (-(3*Alfa2*CW*MW)/SW2**2*(AbbQu4*EL*(-gL(1) + gR(1))*mass(1)*
     -     mass(2)) - (3*Alfa2*CW*MW)/SW2**2*(Abb1*EL*(gR(1)*mass2(1) - 
     -     gL(1)*mass2(2))) - (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*
     -     gL(1)*T)) + Cval(cc2,iintQu6)*((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQu4*EL*gL(1)*mass(1)*mass(2)) - (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*T) + (3*Alfa2*CW*MW)/SW2**2*(Abb1*EL*
     -     gL(1)*(mass2(2) + 2/3.D0*mass2(1)))) + Cval(cc2,iintQu7)*
     -     ((3*Alfa2*
     -     MW)/SW2*(AbbQu4*EL*mass(1)*mass(2)*q(1))/SW - (2*Alfa2*MW)/
     -     SW2*(Abb1*EL*q(1)*T)/SW + (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*
     -     (mass2(2) + 3/2.D0*mass2(1)))/SW) +(Alfa2*
     -     MW)/SW2**2*(Abb1*EL*(SW2*q(1) + CW*(gL(1)*SW)))/SW)/((-MW2 + 
     -     T)*(mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - MW2 - 
     -     S - T24 - U))

** set B
      Cloop(1) = Cloop(1) + 
     -     (Cval(cc1,iintQu9)*((Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(-mass2(1) + mass2(2)))) + 
     -     Cval(cc0,iintQu8)*
     -     ((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*T)/SW) + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Bval(bb0,iintQu4))/SW + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*Bval(bb0,iintQu5)) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Cval(cc00,iintQu8))/SW + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*Cval(cc00,iintQu9)) + 
     -     Cval(cc1,iintQu8)*
     -     (-(Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(2))/SW) + 
     -     Cval(cc0,iintQu9)*
     -     ((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQu4*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(-(gL(2)*mass2(1)) + gR(2)*mass2(2))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*T) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (0*EL*mass(2)*(gR(2) - 1/3.D0*gL(2)))) + 
     -     Cval(cc2,iintQu8)*
     -     (-(3*Alfa2*MW)/SW2*(AbbQu4*EL*mass(1)*mass(2)*q(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*T)/SW - 
     -     (3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*(mass2(2) + 2/3.D0*mass2(1)))/SW) + 
     -     Cval(cc2,iintQu9)*
     -     (-(3*Alfa2*CW*MW)/SW2**2*(AbbQu4*EL*gL(2)*mass(1)*mass(2)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*T) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(mass2(2) + 3/2.D0*mass2(1)))) - 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(2) + CW*(gL(2)*SW)))/SW)/
     -     ((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U))
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the upper quark line
*** which do not have a relative minus sign between currents 5 and 6, for quark 
*** patterns 1 and 3

      subroutine qqVUP_p13(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex tmp1

      
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      tmp1 = (0D0,0D0)


      tmp1 = (Alfa2*MW)/SW2*(Abb1*EL*(gL(1)*gL(2) + q(1)*q(2)))/SW + 
     -     Cval(cc1,iintQu3)*(-(2*Alfa2*MW)/SW2*(AbbQu4*EL*gL(1)*gL(2)
     -     *mass(1)*mass(2))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2)))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*T)/SW) + 
     -     Cval(cc1,iintQu2)*(-(2*Alfa2*MW)/SW2*(AbbQu4*EL*mass(1)
     -     *mass(2)*q(1)*q(2))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(2))*q(1)*q(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*T)/SW) - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(1)*gL(2) + q(1)*q(2))*
     -     Bval(bb0,iintQu1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*Cval(cc00,iintQu2))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*Cval(cc00,iintQu3))/SW + 
     -     Cval(cc2,iintQu2)*(-(2*Alfa2*MW)/SW2*(AbbQu4*EL*mass(1)*
     -     mass(2)*q(1)*q(2))/SW + (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*
     -     q(2)*T)/SW - (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*(mass2(2) + 
     -     1/2.D0*mass2(1)))/SW) + 
     -     Alfa2/SW2*(Abb1*EL*Cval(cc0,iintQu2)*
     -     (-(2*MW*(mass2(1)*q(1)*q(2))) - 2*MW*(mass2(2)*q(1)*q(2)) + 
     -     2*MW*(q(1)*q(2)*T)))/SW + 
     -     1/SW2*(Cval(cc0,iintQu3)*
     -     (Alfa2*(AbbQu4*EL*
     -     (-(2*MW*(gL(1)*gL(2)*mass(1)*mass(2))) + 
     -     2*MW*(gR(1)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(2)*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(2) + 
     -     MW*(-2*MZ2 + 2*T)))))/SW

      tmp1 = tmp1 + Cval(cc2,iintQu3)*(-(2*Alfa2*MW)/SW2*(AbbQu4*EL*
     -     gL(1)*gL(2)*mass(1)*mass(2))/SW + (2*Alfa2*MW)/SW2*(Abb1*EL*
     -     gL(1)*gL(2)*T)/SW - (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*
     -     (mass2(2) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the upper quark line
*** which have a relative minus sign between currents 5 and 6, for quark 
*** patterns 2 and 4

      subroutine qqVUP_n24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

      
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** set A
      Cloop(1) = (Cval(cc1,iintQu7)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(-mass2(1) + mass2(2)))) + 
     -     Cval(cc0,iintQu6)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*T)/SW) - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Bval(bb0,iintQu4))/SW - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*Bval(bb0,iintQu5)) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Cval(cc00,iintQu6))/SW - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*Cval(cc00,iintQu7)) + 
     -     Cval(cc1,iintQu6)*((Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(2))/SW) + 
     -     Cval(cc0,iintQu7)*(-(3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQu3*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(-(gL(2)*mass2(1)) + gR(2)*mass2(2))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*T) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (0*EL*mass(2)*(gR(2) - 1/3.D0*gL(2)))) + 
     -     Cval(cc2,iintQu6)*
     -     ((3*Alfa2*MW)/SW2*(AbbQu3*EL*mass(1)*mass(2)*q(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*T)/SW + 
     -     (3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*(mass2(2) + 2/3.D0*mass2(1)))/SW) + 
     -     Cval(cc2,iintQu7)*
     -     ((3*Alfa2*CW*MW)/SW2**2*(AbbQu3*EL*gL(2)*mass(1)*mass(2)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*T) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(mass2(2) + 3/2.D0*mass2(1)))) + 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(2) + CW*(gL(2)*SW)))/SW)/
     -     ((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U))

** set B
      Cloop(1) = Cloop(1) + (Cval(cc1,iintQu8)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(-mass2(1) + mass2(2)))) + 
     -     Cval(cc1,iintQu9)*((Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1))/SW) + 
     -     Cval(cc0,iintQu9)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*T)/SW) + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Bval(bb0,iintQu4))/SW + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*Bval(bb0,iintQu5)) + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*Cval(cc00,iintQu8)) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Cval(cc00,iintQu9))/SW + 
     -     Cval(cc0,iintQu8)*((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQu3*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(1)*mass2(1) - gL(1)*mass2(2))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*T)) + 
     -     Cval(cc2,iintQu8)*
     -     (-(3*Alfa2*CW*MW)/SW2**2*(AbbQu3*EL*gL(1)*mass(1)*mass(2)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*T) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(mass2(2) + 2/3.D0*mass2(1)))) + 
     -     Cval(cc2,iintQu9)*
     -     (-(3*Alfa2*MW)/SW2*(AbbQu3*EL*mass(1)*mass(2)*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*T)/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*(mass2(2) + 3/2.D0*mass2(1)))/SW) - 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(1) + CW*(gL(1)*SW)))/SW)/
     -    ((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U))
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the upper quark line
*** which do not have a relative minus sign between currents 5 and 6, for quark 
*** patterns 2 and 4

      subroutine qqVUP_p24(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,
     -                   Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex tmp1

      
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      tmp1 = (Alfa2*MW)/SW2*(Abb1*EL*(gL(1)*gL(2) + q(1)*q(2)))/SW + 
     -     Cval(cc1,iintQu3)*(-(2*Alfa2*MW)/SW2*(AbbQu3*EL*gL(1)*gL(2)*
     -     mass(1)*mass(2))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2)))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*T)/SW) + 
     -     Cval(cc1,iintQu2)*(-(2*Alfa2*MW)/SW2*(AbbQu3*EL*mass(1)*
     -     mass(2)*q(1)*q(2))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(2))*q(1)*q(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*T)/SW) - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(1)*gL(2) + q(1)*q(2))*
     -     Bval(bb0,iintQu1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*Cval(cc00,iintQu2))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*Cval(cc00,iintQu3))/
     -     SW + Cval(cc2,iintQu2)*(-(2*Alfa2*MW)/SW2*(AbbQu3*EL*mass(1)*
     -     mass(2)*q(1)*q(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(2)*T)/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(2)*(mass2(2) + 1/2.D0*mass2(1)))/SW) + 
     -     Alfa2/SW2*(Abb1*EL*Cval(cc0,iintQu2)*
     -     (-(2*MW*(mass2(1)*q(1)*q(2))) - 2*MW*(mass2(2)*q(1)*q(2)) + 
     -     2*MW*(q(1)*q(2)*T)))/SW + 1/SW2*(Cval(cc0,iintQu3)*
     -     (Alfa2*(AbbQu3*EL*
     -     (-(2*MW*(gL(1)*gL(2)*mass(1)*mass(2))) + 
     -     2*MW*(gR(1)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(2)*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(2) + 
     -     MW*(-2*MZ2 + 2*T)))))/SW

      tmp1 = tmp1 + Cval(cc2,iintQu3)*
     -     (-(2*Alfa2*MW)/SW2*(AbbQu3*EL*gL(1)*gL(2)*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(2)*T)/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(2)*(mass2(2) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the lower quark line
*** which have a relative minus sign between currents 5 and 6, for quark 
*** patterns 1 and 2

      subroutine qqVDOWN_n12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex tmp2, tmp3

      
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


** set A
      tmp2 = Cval(cc1,iintQl6)*
     -     ((Alfa2*MW)/SW2*(Abb1*EL*(-mass2(3) + mass2(4))*q(4))/SW)-
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Bval(bb0,iintQl4))/SW - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Bval(bb0,iintQl5)) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Cval(cc00,iintQl6))/SW - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Cval(cc00,iintQl7)) + 
     -     Cval(cc0,iintQl6)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + 
     -     2*mass2(4) + MH2)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*(S + T24 + U))/SW) + 
     -     Cval(cc1,iintQl7)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(-mass2(3) + mass2(4)))) + 
     -     Cval(cc0,iintQl7)*(-(3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQl3*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*(S + T24 + U)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(4)*mass2(4) + 
     -     gL(4)*(2/3.D0*mass2(1) + 2/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4) + 2/3.D0*MH2)))) + 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(4) + CW*(gL(4)*SW)))/SW

      tmp2 = tmp2 + Cval(cc2,iintQl7)*((3*Alfa2*CW*MW)/SW2**2*(AbbQl3*
     -     EL*gL(4)*mass(3)*mass(4)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*(S + T24 + U)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + MH2 + 
     -     1/2.D0*mass2(3)))) + Cval(cc2,iintQl6)*((3*Alfa2*MW)/SW2*
     -     (AbbQl3*EL*mass(3)*mass(4)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*(S + T24 + U))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(mass2(1) + mass2(2) + mass2(3) + MH2 - 
     -     1/2.D0*mass2(4)))/SW) 

      Cloop(1) = tmp2/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 

** set B
      tmp3 = Cval(cc1,iintQl9)*
     -    ((Alfa2*MW)/SW2*(Abb1*EL*(-mass2(3) + mass2(4))*q(3))/SW)+
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Bval(bb0,iintQl4))/SW + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Bval(bb0,iintQl5)) + 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Cval(cc00,iintQl8)) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Cval(cc00,iintQl9))/SW + 
     -     Cval(cc0,iintQl9)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + 3*mass2(3) + MH2)*q(3))/
     -     SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(3)*(S + T24 + U))/SW) + Cval(cc1,iintQl8)*
     -     (-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(-mass2(3) + mass2(4)))) + 
     -     Cval(cc0,iintQl8)*((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQl3*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*(S + T24 + U)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(2/3.D0*mass2(1) + 2/3.D0*mass2(2) + 
     -     4/3.D0*mass2(3) - 1/3.D0*mass2(4) + 2/3.D0*MH2)))) - 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(3) + CW*(gL(3)*SW)))/SW

      tmp3 = tmp3 + Cval(cc2,iintQl9)*(-(3*Alfa2*MW)/SW2*(AbbQl3*EL*
     -     mass(3)*mass(4)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*(S + T24 + U))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(3)*(mass2(1) + mass2(2) + MH2 + 
     -     1/2.D0*mass2(3)))/SW) + Cval(cc2,iintQl8)*(-(3*Alfa2*CW*MW)/
     -     SW2**2*(AbbQl3*EL*gL(3)*mass(3)*mass(4)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*(S + T24 + U)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) + mass2(3) + MH2 - 
     -     1/2.D0*mass2(4)))) 

      Cloop(1) = tmp3/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the lower quark line
*** which do not have a relative minus sign between currents 5 and 6, for quark 
*** patterns 1 and 2

      subroutine qqVDOWN_p12(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex dup1, tmp1

      
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup1 = mass2(1) + mass2(2) + mass2(3) - mass2(4) + MH2

      tmp1 = (Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gL(4) + q(3)*q(4)))/SW + 
     -     Cval(cc0,iintQl3)*((2*Alfa2*MW)/SW2*
     -     (AbbQl3*EL*(-(gL(3)*gL(4)) + gR(3)*gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*gL(4)*
     -     (-mass2(1) - mass2(2) - mass2(3) - MH2 + MZ2))/SW-
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*gL(4)*(S + T24 + U))/SW)
     -     - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gL(4) + q(3)*q(4))*Bval(bb0,iintQl1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*q(4)*Cval(cc00,iintQl2))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*gL(4)*Cval(cc00,iintQl3))/SW + 
     -     Alfa2/SW2*(Abb1*EL*Cval(cc0,iintQl2)*
     -     (2*MW*(mass2(1)*q(3)*q(4)) + 2*MW*(mass2(2)*q(3)*q(4)) + 
     -     2*MW*(mass2(3)*q(3)*q(4)) + 2*MW*(MH2*q(3)*q(4)) + 
     -     MW*(q(3)*q(4)*(-2*S - 2*T24 - 2*U))))/SW + 
     -     1/SW2*(Cval(cc1,iintQl3)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*gL(3)*gL(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(3)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(2) + 2*MW*MH2 + 
     -     MW*(-2*S - 2*T24 - 2*U)))))/SW + 
     -     1/SW2*(Cval(cc1,iintQl2)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*mass(3)*mass(4)*q(3)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (2*MW*(mass2(1)*q(3)*q(4)) + 2*MW*(mass2(2)*q(3)*q(4)) + 
     -     2*MW*(MH2*q(3)*q(4)) + 
     -     MW*(q(3)*q(4)*(-2*S - 2*T24 - 2*U))))))/SW

      tmp1 = tmp1 + 
     -     1/SW2*(Cval(cc2,iintQl3)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*gL(3)*gL(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*gL(3)*gL(4)) + 
     -     MW*(EL*gL(3)*gL(4)*(-2*S - 2*T24 - 2*U))))))/SW + 
     -     1/SW2*(Cval(cc2,iintQl2)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*mass(3)*mass(4)*q(3)*q(4)) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(3)*q(4)) + 
     -     MW*(EL*q(3)*q(4)*(-2*S - 2*T24 - 2*U))))))/SW

      Cloop(1) = tmp1/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set loops on the lower quark line
*** which have a relative minus sign between currents 5 and 6, for quark 
*** patterns 3 and 4

      subroutine qqVDOWN_n34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex tmp2, tmp3

      
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      tmp2 = Cval(cc1,iintQl7)*((Alfa2*MW)/SW2*((-3*0 + 0)*
     -     EL*mass(3)*q(3))/SW - 
     -     (Alfa2*MW)/SW2*(Abb1*EL*(-mass2(3) + mass2(4))*q(3))/SW)-
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Bval(bb0,iintQl4))/SW - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Bval(bb0,iintQl5)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Cval(cc00,iintQl6)) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Cval(cc00,iintQl7))/SW + 
     -     Cval(cc0,iintQl7)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + 3*mass2(3) + MH2)*q(3))/
     -     SW + (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*(S + T24 + U))/SW) + 
     -     Cval(cc1,iintQl6)*((Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(-mass2(3) + mass2(4)))) + 
     -     Cval(cc0,iintQl6)*(-(3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQl3*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*(S + T24 + U)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(2/3.D0*mass2(1) + 2/3.D0*mass2(2) + 
     -     4/3.D0*mass2(3) - 1/3.D0*mass2(4) + 2/3.D0*MH2)))) + 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(3) + CW*(gL(3)*SW)))/SW

      tmp2 = tmp2 + Cval(cc2,iintQl7)*((3*Alfa2*MW)/SW2*(AbbQl3*EL*
     -     mass(3)*mass(4)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*(S + T24 + U))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*(mass2(1) + mass2(2) + MH2 + 
     -     1/2.D0*mass2(3)))/SW) + 
     -     Cval(cc2,iintQl6)*((3*Alfa2*CW*MW)/SW2**2*(AbbQl3*EL*gL(3)*
     -     mass(3)*mass(4)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*(S + T24 + U)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) + mass2(3) + MH2 - 
     -     1/2.D0*mass2(4)))) 

      Cloop(1) = tmp2/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 

** set B
      tmp3 = Cval(cc1,iintQl8)*
     -     (-(Alfa2*MW)/SW2*(Abb1*EL*(-mass2(3) + mass2(4))*q(4))/SW)+
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Bval(bb0,iintQl4))/SW + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Bval(bb0,iintQl5)) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Cval(cc00,iintQl8))/SW + 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Cval(cc00,iintQl9)) + 
     -     Cval(cc0,iintQl8)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + 
     -     2*mass2(4) + MH2)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*(S + T24 + U))/SW) + 
     -     Cval(cc1,iintQl9)*((Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(-mass2(3) + mass2(4)))) + 
     -     Cval(cc0,iintQl9)*((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbQl3*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*(S + T24 + U)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*(Abb1*EL*(gR(4)*mass2(4) + 
     -     gL(4)*(2/3.D0*mass2(1) + 2/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4) + 2/3.D0*MH2)))) - 
     -     (Alfa2*MW)/SW2**2*(Abb1*EL*(SW2*q(4) + CW*(gL(4)*SW)))/SW

      tmp3 = tmp3 + Cval(cc2,iintQl9)*(-(3*Alfa2*CW*MW)/SW2**2*(AbbQl3*
     -     EL*gL(4)*mass(3)*mass(4)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*(S + T24 + U)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + MH2 + 
     -     1/2.D0*mass2(3)))) + Cval(cc2,iintQl8)*(-(3*Alfa2*MW)/SW2*
     -     (AbbQl3*EL*mass(3)*mass(4)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*(S + T24 + U))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(mass2(1) + mass2(2) + mass2(3) + MH2 - 
     -     1/2.D0*mass2(4)))/SW) 

      Cloop(1) = tmp3/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the diagram set of loops on the lower quark line
*** which do not have a relative minus sign between currents 5 and 6, for quark 
*** patterns 3 and 4

      subroutine qqVDOWN_p34(mass,mass2,gL,gR,q,S,T,T24,T14,S34,U,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass (and squared), left/right coupling and charge of external particles
      double precision mass(5), mass2(5), gL(4), gR(4), q(4) 

** Dummy variable
      double complex dup1, tmp1

      
      double precision Alfa, Alfa2


** Parameter, abbreviation and integral declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h" 


      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup1 = mass2(1) + mass2(2) + mass2(3) - mass2(4) + MH2

      tmp1 = (Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gL(4) + q(3)*q(4)))/SW + 
     -     Cval(cc0,iintQl3)*((2*Alfa2*MW)/SW2*
     -     (AbbQl3*EL*(-(gL(3)*gL(4)) + gR(3)*gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*gL(4)*
     -     (-mass2(1) - mass2(2) - mass2(3) - MH2 + MZ2))/SW-
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*gL(4)*(S + T24 + U))/SW)
     -     - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gL(4) + q(3)*q(4))*Bval(bb0,iintQl1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*q(4)*Cval(cc00,iintQl2))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*gL(4)*Cval(cc00,iintQl3))/SW + 
     -     Alfa2/SW2*(Abb1*EL*Cval(cc0,iintQl2)*
     -     (2*MW*(mass2(1)*q(3)*q(4)) + 2*MW*(mass2(2)*q(3)*q(4)) + 
     -     2*MW*(mass2(3)*q(3)*q(4)) + 2*MW*(MH2*q(3)*q(4)) + 
     -     MW*(q(3)*q(4)*(-2*S - 2*T24 - 2*U))))/SW + 
     -     1/SW2*(Cval(cc1,iintQl3)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*gL(3)*gL(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(3)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(2) + 2*MW*MH2 + 
     -     MW*(-2*S - 2*T24 - 2*U)))))/SW + 
     -     1/SW2*(Cval(cc1,iintQl2)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*mass(3)*mass(4)*q(3)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (2*MW*(mass2(1)*q(3)*q(4)) + 2*MW*(mass2(2)*q(3)*q(4)) + 
     -     2*MW*(MH2*q(3)*q(4)) + 
     -     MW*(q(3)*q(4)*(-2*S - 2*T24 - 2*U))))))/SW

      tmp1 = tmp1 + 
     -     1/SW2*(Cval(cc2,iintQl3)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*gL(3)*gL(4)*mass(3)*mass(4)) + Alfa2*
     -     (Abb1*(2*MW*(dup1*EL*gL(3)*gL(4)) + 
     -     MW*(EL*gL(3)*gL(4)*(-2*S - 2*T24 - 2*U))))))/SW + 
     -     1/SW2*(Cval(cc2,iintQl2)*
     -     (-2*Alfa2*MW*(AbbQl3*EL*mass(3)*mass(4)*q(3)*q(4)) + 
     -     Alfa2*(Abb1*(2*MW*(dup1*EL*q(3)*q(4)) + 
     -     MW*(EL*q(3)*q(4)*(-2*S - 2*T24 - 2*U))))))/SW

      Cloop(1) = tmp1/((-MW2 + T)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     mass2(4) + MH2 - MW2 - S - T24 - U)) 
      
      
      end


********************************************************************************
