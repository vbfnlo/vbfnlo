********************************************************************************
********************************************************************************
***                                                                          ***
*** CC_boxes.F                                                               ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contains all of the matrix elements for the electroweak        ***
*** box diagrams for the charged currents.                                   ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine calculates the boxesNu part of the loop component for
*** qk = 1 and qk = 3 (qq and qqbar), with photons.

      subroutine boxesNu_p13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5)  

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, tmp1

      integer qk

      double precision Alfa, Alfa2


** Paremter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = MH2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)


      tmp1 = Dval(dd0,iintU7)*
     -     (-(2*Alfa2*MW)/SW2*(Abb1*EL*(mass2(1) + MH2)*q(2))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(2)*(Abb34 + Abb1*(-2*S34 - 2*T)))/SW) 
     -     + (Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(-q(1) + q(2))*Cval(cc0,iintU1))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Cval(cc0,iintU2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Cval(cc0,iintU3))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Dval(dd00,iintU6))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Dval(dd00,iintU7))/SW +
     -     Dval(dd0,iintU6)*((2*Alfa2*MW)/SW2*(Abb1*EL*(mass2(2) + MH2)*
     -     q(1))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(1)*(Abb10 + Abb1*(2*T + 2*T14)))/SW) +
     -     1/(2.D0*SW2)*(Dval(dd0,iintU5)*
     -     (2*Alfa2*MW*(Abb10*EL*q(1))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd0,iintU4)*
     -     (2*Alfa2*MW*(Abb34*EL*q(2))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd1,iintU5)*
     -     (2*Alfa2*MW*(Abb10*EL*q(1)) + 
     -     Alfa2*(Abb1*(2*MW*(dup1*EL*q(1)) - 2*MW*(EL*q(1)*S34))))
     -     )/SW + 1/(2.D0*SW2)*
     -     (Dval(dd1,iintU4)*
     -     (2*Alfa2*MW*(Abb34*EL*q(2)) + 
     -     Alfa2*(Abb1*
     -     (-(2*MW*(dup1*EL*q(2))) + 2*MW*(EL*q(2)*T14)))))/SW

      tmp1 = tmp1 + Dval(dd1,iintU6)*
     -     ((Alfa2*MW)/SW2*
     -     (Abb1*EL*(3*mass2(1) - 2*mass2(2) + MH2)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumU22 + Abb1*(-S34 + 2*T14)))/SW) + 
     -     1/SW2*(Dval(dd12,iintU6)*
     -     (2*Alfa2*MW*(AbbSumU8*EL*q(1))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU8*EL*q(2)))))/SW + 
     -     1/SW2*(Dval(dd1,iintU7)*
     -     (Alfa2*MW*(Abb34*EL*q(2)) + 
     -     Alfa2*(Abb1*EL*(-(MW*(mass2(2)*q(2))) + MW*(q(2)*T14))))
     -     )/SW

      tmp1 = tmp1 + 1/(2.D0*SW2)*
     -     (Dval(dd2,iintU5)*
     -     (2*Alfa2*MW*(Abb10*EL*q(1)) + 
     -     Alfa2*MW*(Abb1*EL*mass2(1)*q(1)) - 
     -     Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(1))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (2*Alfa2*MW*(Abb34*EL*q(2))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd2,iintU4)*
     -     (2*Alfa2*MW*(Abb34*EL*q(2)) + 
     -     Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(2)) - 
     -     Alfa2*MW*(Abb1*EL*mass2(2)*q(2))))/SW + 
     -     Dval(dd2,iintU7)*(-(Alfa2*MW)/SW2*
     -     (Abb1*EL*(-2*mass2(1) + 3*mass2(2) + MH2)*q(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*(AbbSumU12 + Abb1*(-S34 + 1/2.D0*T14)))/SW) + 
     -     1/SW2*(Dval(dd2,iintU6)*
     -     (Alfa2*MW*(Abb10*EL*q(1)) + 
     -     Alfa2*(Abb1*EL*(MW*(mass2(1)*q(1)) - MW*(q(1)*S34)))))/
     -     SW

      tmp1 = tmp1 + Dval(dd3,iintU6)*((3*Alfa2*MW)/SW2*(Abb14*EL*
     -     mass(1)*mass(2)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(1)*(Abb10 + Abb1*(2*T + 2*T14)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*(mass2(2) + MH2 + 5/2.D0*mass2(1)))/SW) + 
     -     Dval(dd3,iintU7)*(-(3*Alfa2*MW)/SW2*(Abb14*EL*mass(1)*
     -     mass(2)*q(2))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(2)*(Abb34 + Abb1*(-2*S34 - 2*T)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*(mass2(1) + MH2 + 5/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Dval(dd23,iintU6)*
     -     (2*Alfa2*MW*(Abb10*EL*q(1))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd3,iintU5)*
     -     (2*Alfa2*MW*(Abb10*EL*q(1)) + 
     -     Alfa2*(Abb1*EL*
     -     (MW*(mass2(1)*q(1)) + MW*(mass2(2)*q(1)) - 
     -     2*MW*(q(1)*S34)))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd3,iintU4)*
     -     (2*Alfa2*MW*(Abb34*EL*q(2)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(MW*(mass2(1)*q(2))) - MW*(mass2(2)*q(2)) + 
     -     2*MW*(q(2)*T14)))))/SW

      Cloop(1) = tmp1/
     -     (mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - 
     -     MW2 - S - T24 - U) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************


*** This subroutine calculates the boxesNu part of the loop component for
*** qk = 2 and qk = 4 (qbarq and qbarqbar), with photons.

      subroutine boxesNu_p24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5)  

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, tmp1

      integer qk

      double precision Alfa, Alfa2


** Paremter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = MH2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)


      tmp1 = Dval(dd0,iintU6)*
     -     ((2*Alfa2*MW)/SW2*(Abb1*EL*(mass2(2) + MH2)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(1)*(Abb35 + Abb1*(2*T + 2*T14)))/SW)
     -     + (Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(-q(1) + q(2))*Cval(cc0,iintU1))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Cval(cc0,iintU2))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Cval(cc0,iintU3))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*Dval(dd00,iintU6))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*Dval(dd00,iintU7))/SW + 
     -     Dval(dd0,iintU7)*(-(2*Alfa2*MW)/SW2*(Abb1*EL*(mass2(1) + 
     -     MH2)*q(2))/SW - (Alfa2*MW)/SW2*
     -     (EL*q(2)*(Abb2 + Abb1*(-2*S34 - 2*T)))/SW) + 
     -     1/(2.D0*SW2)*(Dval(dd0,iintU5)*
     -     (2*Alfa2*MW*(Abb35*EL*q(1))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd0,iintU4)*
     -     (2*Alfa2*MW*(Abb2*EL*q(2))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd1,iintU5)*
     -     (2*Alfa2*MW*(Abb35*EL*q(1)) + 
     -     Alfa2*(Abb1*(2*MW*(dup1*EL*q(1)) - 2*MW*(EL*q(1)*S34))))
     -     )/SW + 1/(2.D0*SW2)*
     -     (Dval(dd1,iintU4)*
     -     (2*Alfa2*MW*(Abb2*EL*q(2)) + 
     -     Alfa2*(Abb1*
     -     (-(2*MW*(dup1*EL*q(2))) + 2*MW*(EL*q(2)*T14)))))/SW

      tmp1 = tmp1 + Dval(dd1,iintU6)*
     -     ((Alfa2*MW)/SW2*
     -     (Abb1*EL*(3*mass2(1) - 2*mass2(2) + MH2)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumU22 + Abb1*(-S34 + 2*T14)))/SW) + 
     -     1/SW2*(Dval(dd12,iintU6)*
     -     (2*Alfa2*MW*(AbbSumU7*EL*q(1))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU7*EL*q(2)))))/SW + 
     -     1/SW2*(Dval(dd1,iintU7)*
     -     (Alfa2*MW*(Abb2*EL*q(2)) + 
     -     Alfa2*(Abb1*EL*(-(MW*(mass2(2)*q(2))) + MW*(q(2)*T14))))
     -     )/SW

      tmp1 = tmp1 + 1/(2.D0*SW2)*(Dval(dd2,iintU5)*
     -     (2*Alfa2*MW*(Abb35*EL*q(1)) + 
     -     Alfa2*MW*(Abb1*EL*mass2(1)*q(1)) - 
     -     Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(1))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd2,iintU4)*
     -     (2*Alfa2*MW*(Abb2*EL*q(2)) + 
     -     Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(2)) - 
     -     Alfa2*MW*(Abb1*EL*mass2(2)*q(2))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (2*Alfa2*MW*(Abb2*EL*q(2))))/SW + 
     -     Dval(dd2,iintU7)*(-(Alfa2*MW)/SW2*
     -     (Abb1*EL*(-2*mass2(1) + 3*mass2(2) + MH2)*q(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*(AbbSumU11 + Abb1*(-S34 + 1/2.D0*T14)))/SW) + 
     -     1/SW2*(Dval(dd2,iintU6)*
     -     (Alfa2*MW*(Abb35*EL*q(1)) + 
     -     Alfa2*(Abb1*EL*(MW*(mass2(1)*q(1)) - MW*(q(1)*S34)))))/
     -     SW

      tmp1 = tmp1 + Dval(dd3,iintU6)*((3*Alfa2*MW)/SW2*(Abb21*EL*
     -     mass(1)*mass(2)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(1)*(Abb35 + Abb1*(2*T + 2*T14)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*(mass2(2) + MH2 + 5/2.D0*mass2(1)))/SW) + 
     -     Dval(dd3,iintU7)*(-(3*Alfa2*MW)/SW2*(Abb21*EL*mass(1)*
     -     mass(2)*q(2))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(2)*(Abb2 + Abb1*(-2*S34 - 2*T)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*(mass2(1) + MH2 + 5/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Dval(dd23,iintU6)*
     -     (2*Alfa2*MW*(Abb35*EL*q(1))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd3,iintU5)*
     -     (2*Alfa2*MW*(Abb35*EL*q(1)) + 
     -     Alfa2*(Abb1*EL*
     -     (MW*(mass2(1)*q(1)) + MW*(mass2(2)*q(1)) - 
     -     2*MW*(q(1)*S34)))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd3,iintU4)*
     -     (2*Alfa2*MW*(Abb2*EL*q(2)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(MW*(mass2(1)*q(2))) - MW*(mass2(2)*q(2)) + 
     -     2*MW*(q(2)*T14)))))/SW

      Cloop(1) = tmp1/
     -     (mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - 
     -     MW2 - S - T24 - U) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNu part of the loop component for
*** qk = 1 and qk = 3 (qq and qqbar), with no photons.

      subroutine boxesNu_nop13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5)  

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup2, dup3, tmp2

      integer qk

      double precision Alfa, Alfa2


** Paremter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = 3*mass2(1) - 2*mass2(2) + MH2
      dup3 = MH2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)

      tmp2 = (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*Cval(cc0,iintU10)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(-gL(1) + gL(2))*Cval(cc0,iintU11)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*Cval(cc0,iintU12)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*Cval(cc0,iintU8)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*Cval(cc0,iintU9)) + 
     -     Dval(dd0,iintU14)*(-(3*Alfa2*CW*MW)/SW2**2*
     -     (Abb14*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(1)*mass2(1) + 
     -     gL(1)*(-(4/3.D0*mass2(2)) - 1/3.D0*MH2))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb34 + Abb1*(-T - 1/2.D0*S34)))) + 
     -     Dval(dd0,iintU13)*((3*Alfa2*CW*MW)/SW2**2*
     -     (Abb14*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(2)*mass2(2) + 
     -     gL(2)*(-(2/3.D0*mass2(1)) - 2/3.D0*mass2(2) - 2/3.D0*MH2))
     -     ) - (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(Abb34 + Abb1*(-S34 - T + 1/2.D0*T14)))) + 
     -     1/(CW*SW2)*(Dval(dd0,iintU15)*
     -     (Alfa2*(Abb1*(MW*(dup3*EL*gL(1)) - MW*(EL*gL(1)*S34)))))+
     -     1/(CW*SW2)*(Dval(dd0,iintU16)*
     -     (Alfa2*(Abb1*(-(MW*(dup3*EL*gL(2))) + MW*(EL*gL(2)*T14))))
     -     )

      tmp2 = tmp2 + (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*Dval(dd00,iintU13)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*Dval(dd00,iintU14)) + 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*Dval(dd00,iintU17)) - 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*Dval(dd00,iintU18)) + 
     -     Dval(dd0,iintU17)*((3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb14*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(2)*(Abb10 + Abb1*(2*T + 2*T14))) + 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*(gR(2)*mass2(2) + 
     -     gL(2)*(-(5/3.D0*mass2(1)) - 2/3.D0*MH2)))) + 
     -     Dval(dd0,iintU18)*(-(3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb14*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(Abb34 + Abb1*(-2*S34 - 2*T))) - 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*(gR(1)*mass2(1) + 
     -     gL(1)*(-(5/3.D0*mass2(2)) - 2/3.D0*MH2)))) + 
     -     Dval(dd1,iintU13)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(-2*mass2(1) + 2*mass2(2) + MH2)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(Abb34 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14))))

      tmp2 = tmp2 + Dval(dd1,iintU15)*
     -     ((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb14*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(1)*mass2(1) + gL(1)*(mass2(2) + 2*MH2))) + 
     -     (Alfa2*MW)/(CW*SW2)*(EL*gL(1)*(Abb10 - Abb1*S34))) + 
     -     Dval(dd1,iintU16)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb14*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(2)*mass2(2) + gL(2)*(mass2(1) + 2*MH2))) + 
     -     (Alfa2*MW)/(CW*SW2)*(EL*gL(2)*(Abb34 + Abb1*T14))) + 
     -     Dval(dd1,iintU14)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(2*mass2(1) - 2*mass2(2) + MH2)) + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb34 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintU17)*
     -     (Alfa2*CW*MW*(AbbSumU22*EL*gL(2)) + 
     -     Alfa2*CW*(Abb1*
     -     (-(MW*(dup2*EL*gL(2))) + MW*(EL*gL(2)*(-S34 + 2*T14)))
     -     )))

      tmp2 = tmp2 + 
     -     1/SW2**2*(Dval(dd11,iintU14)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU8*EL*gL(1))))) + 
     -     1/SW2**2*(Dval(dd11,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU8*EL*gL(2)))) + 
     -     1/SW2**2*(Dval(dd12,iintU14)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU8*EL*gL(1))))) + 
     -     1/SW2**2*(Dval(dd12,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU8*EL*gL(2)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintU18)*
     -     (-(Alfa2*CW*MW*(Abb34*EL*gL(1))) + Alfa2*CW*(Abb1*EL*
     -     gL(1)*(MW*mass2(2) - MW*T14))))

      tmp2 = tmp2 + 
     -     1/SW2**2*(Dval(dd13,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU27*EL*gL(2)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd13,iintU18)*
     -     (-(2*Alfa2*CW*MW*(Abb34*EL*gL(1))) )) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintU18)*
     -     (2*Alfa2*CW*MW*(AbbSumU8*EL*gL(1)) )) + 
     -     1/SW2**2*(Dval(dd13,iintU14)*
     -     (2*Alfa2*CW*MW*(Abb34*EL*gL(1)) )) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintU17)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU8*EL*gL(2))))) + 
     -     1/SW2**2*(Dval(dd2,iintU13)*
     -     (-(Alfa2*CW*MW*(Abb34*EL*gL(2))) + 
     -     Alfa2*CW*(Abb1*EL*gL(2)*(MW*mass2(2) - MW*T14))))

      tmp2 = tmp2 + Dval(dd2,iintU16)*
     -     (-(Alfa2*MW)/(2.D0*CW*SW2)*(Abb14*EL*gR(2)*mass(1)*mass(2)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(2)*T14) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(2)*(MH2 + 1/2.D0*mass2(1)))
     -     ) + Dval(dd2,iintU15)*
     -     ((Alfa2*MW)/(2.D0*CW*SW2)*(Abb14*EL*gR(1)*mass(1)*mass(2)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(1)*S34) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(1)*(MH2 + 1/2.D0*mass2(2)))
     -     ) + Dval(dd2,iintU18)*
     -     ((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*(-2*mass2(1) + 3*mass2(2) + MH2)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(AbbSumU12 + Abb1*(-S34 + 1/2.D0*T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd2,iintU17)*
     -     (-(Alfa2*CW*MW*(Abb10*EL*gL(2))) + Alfa2*CW*
     -     (Abb1*EL*gL(2)*(-(MW*mass2(1)) + MW*S34)))) + 
     -     1/SW2**2*(Dval(dd2,iintU14)*
     -     (Alfa2*CW*MW*(AbbSumU22*EL*gL(1)) + 
     -     Alfa2*CW*(Abb1*
     -     (-(MW*(dup2*EL*gL(1))) + MW*(EL*gL(1)*(-S34 + 2*T14)))
     -     )))

      tmp2 = tmp2 + Dval(dd3,iintU13)*
     -     (-(3*Alfa2*CW*MW)/SW2**2*(Abb14*EL*gL(2)*mass(1)*mass(2)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(mass2(1) + 2*mass2(2) + MH2)) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(4*Abb34 + Abb1*(-2*S34 - 2*T + T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintU17)*
     -     (-(2*Alfa2*CW*MW*(Abb10*EL*gL(2))))) + 
     -     1/SW2**2*(Dval(dd23,iintU13)*
     -     (-(2*Alfa2*CW*MW*(Abb34*EL*gL(2)))))

      tmp2 = tmp2 + 
     -     Dval(dd3,iintU17)*(-(3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb14*EL*gL(2)*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(2)*(Abb10 + Abb1*(2*T + 2*T14))) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*(mass2(2) + MH2 + 5/2.D0*mass2(1)))) + 
     -     Dval(dd3,iintU18)*((3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb14*EL*gL(1)*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(Abb34 + Abb1*(-2*S34 - 2*T))) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*(mass2(1) + MH2 + 5/2.D0*mass2(2)))) + 
     -     1/SW2**2*(Dval(dd33,iintU13)*
     -     (-(2*Alfa2*CW*MW*(Abb34*EL*gL(2))))) + 
     -     1/(2.D0*CW*SW2)*(Dval(dd3,iintU15)*
     -     (2*Alfa2*MW*(Abb1*EL*gL(1)*MH2))) + 
     -     1/(2.D0*CW*SW2)*(Dval(dd3,iintU16)*
     -     (-2*Alfa2*MW*(Abb1*EL*gL(2)*MH2))) + 
     -     Dval(dd3,iintU14)*((3*Alfa2*CW*MW)/SW2**2*(Abb14*EL*gL(1)*
     -     mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(2*mass2(1) + 4*mass2(2) + MH2)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb34 + Abb1*(-T - 1/2.D0*S34))))

      Cloop(1) = tmp2/
     -     (mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - 
     -     MW2 - S - T24 - U) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNu part of the loop component for
*** qk = 2 and qk = 4 (qbarq and qbarqbar), with no photons.

      subroutine boxesNu_nop24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5)  

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup2, dup3, tmp2

      integer qk

      double precision Alfa, Alfa2


** Paremter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = 3*mass2(1) - 2*mass2(2) + MH2
      dup3 = MH2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)

      tmp2 = (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*Cval(cc0,iintU10)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(-gL(1) + gL(2))*Cval(cc0,iintU11)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*Cval(cc0,iintU12)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(2)*Cval(cc0,iintU8)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*Cval(cc0,iintU9)) + 
     -     Dval(dd0,iintU14)*(-(3*Alfa2*CW*MW)/SW2**2*
     -     (Abb21*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(1)*mass2(1) + 
     -     gL(1)*(-(4/3.D0*mass2(2)) - 1/3.D0*MH2))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb2 + Abb1*(-T - 1/2.D0*S34)))) + 
     -     Dval(dd0,iintU13)*((3*Alfa2*CW*MW)/SW2**2*
     -     (Abb21*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*(gR(2)*mass2(2) + 
     -     gL(2)*(-(2/3.D0*mass2(1)) - 2/3.D0*mass2(2) - 2/3.D0*MH2))
     -     ) - (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(Abb2 + Abb1*(-S34 - T + 1/2.D0*T14)))) + 
     -     1/(CW*SW2)*(Dval(dd0,iintU15)*
     -     (Alfa2*(Abb1*(MW*(dup3*EL*gL(1)) - MW*(EL*gL(1)*S34)))))+
     -     1/(CW*SW2)*(Dval(dd0,iintU16)*
     -     (Alfa2*(Abb1*(-(MW*(dup3*EL*gL(2))) + MW*(EL*gL(2)*T14))))
     -     )

      tmp2 = tmp2 + (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*Dval(dd00,iintU13)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(1)*Dval(dd00,iintU14)) + 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*Dval(dd00,iintU17)) - 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*Dval(dd00,iintU18)) + 
     -     Dval(dd0,iintU17)*((3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb21*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(2)*(Abb35 + Abb1*(2*T + 2*T14))) + 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*(gR(2)*mass2(2) + 
     -     gL(2)*(-(5/3.D0*mass2(1)) - 2/3.D0*MH2)))) + 
     -     Dval(dd0,iintU18)*(-(3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb21*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(Abb2 + Abb1*(-2*S34 - 2*T))) - 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*(gR(1)*mass2(1) + 
     -     gL(1)*(-(5/3.D0*mass2(2)) - 2/3.D0*MH2)))) + 
     -     Dval(dd1,iintU13)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(-2*mass2(1) + 2*mass2(2) + MH2)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(Abb2 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14))))

      tmp2 = tmp2 + Dval(dd1,iintU15)*
     -     ((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb21*EL*(-gL(1) + gR(1))*mass(1)*mass(2)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(1)*mass2(1) + gL(1)*(mass2(2) + 2*MH2))) + 
     -     (Alfa2*MW)/(CW*SW2)*(EL*gL(1)*(Abb35 - Abb1*S34))) + 
     -     Dval(dd1,iintU16)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb21*EL*(-gL(2) + gR(2))*mass(1)*mass(2)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(2)*mass2(2) + gL(2)*(mass2(1) + 2*MH2))) + 
     -     (Alfa2*MW)/(CW*SW2)*(EL*gL(2)*(Abb2 + Abb1*T14))) + 
     -     Dval(dd1,iintU14)*(-(Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(2*mass2(1) - 2*mass2(2) + MH2)) + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb2 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintU17)*
     -     (Alfa2*CW*MW*(AbbSumU22*EL*gL(2)) + 
     -     Alfa2*CW*(Abb1*
     -     (-(MW*(dup2*EL*gL(2))) + MW*(EL*gL(2)*(-S34 + 2*T14)))
     -     )))

      tmp2 = tmp2 + 1/SW2**2*(Dval(dd11,iintU14)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU7*EL*gL(1))))) + 
     -     1/SW2**2*(Dval(dd11,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU7*EL*gL(2)))) + 
     -     1/SW2**2*(Dval(dd12,iintU14)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU7*EL*gL(1))))) + 
     -     1/SW2**2*(Dval(dd12,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU7*EL*gL(2)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintU18)*
     -     (-(Alfa2*CW*MW*(Abb2*EL*gL(1))) + 
     -     Alfa2*CW*(Abb1*EL*gL(1)*(MW*mass2(2) - MW*T14))))

      tmp2 = tmp2 + 
     -     1/SW2**2*(Dval(dd13,iintU13)*
     -     (2*Alfa2*CW*MW*(AbbSumU26*EL*gL(2)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd13,iintU18)*
     -     (-(2*Alfa2*CW*MW*(Abb2*EL*gL(1))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintU18)*
     -     (2*Alfa2*CW*MW*(AbbSumU7*EL*gL(1)))) + 
     -     1/SW2**2*(Dval(dd13,iintU14)*
     -     (2*Alfa2*CW*MW*(Abb2*EL*gL(1)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintU17)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU7*EL*gL(2))))) + 
     -     1/SW2**2*(Dval(dd2,iintU13)*
     -     (-(Alfa2*CW*MW*(Abb2*EL*gL(2))) + 
     -     Alfa2*CW*(Abb1*EL*gL(2)*(MW*mass2(2) - MW*T14))))

      tmp2 = tmp2 + Dval(dd2,iintU16)*
     -     (-(Alfa2*MW)/(2.D0*CW*SW2)*(Abb21*EL*gR(2)*mass(1)*mass(2)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(2)*T14) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(2)*(MH2 + 1/2.D0*mass2(1)))
     -     ) + Dval(dd2,iintU15)*
     -     ((Alfa2*MW)/(2.D0*CW*SW2)*(Abb21*EL*gR(1)*mass(1)*mass(2)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(1)*S34) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(1)*(MH2 + 1/2.D0*mass2(2)))
     -     ) + Dval(dd2,iintU18)*
     -     ((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*(-2*mass2(1) + 3*mass2(2) + MH2)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(AbbSumU11 + Abb1*(-S34 + 1/2.D0*T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd2,iintU17)*
     -     (-(Alfa2*CW*MW*(Abb35*EL*gL(2))) + 
     -     Alfa2*CW*(Abb1*EL*gL(2)*(-(MW*mass2(1)) + MW*S34)))) + 
     -     1/SW2**2*(Dval(dd2,iintU14)*
     -     (Alfa2*CW*MW*(AbbSumU22*EL*gL(1)) + 
     -     Alfa2*CW*(Abb1*
     -     (-(MW*(dup2*EL*gL(1))) + MW*(EL*gL(1)*(-S34 + 2*T14)))
     -     )))

      tmp2 = tmp2 + Dval(dd3,iintU13)*
     -     (-(3*Alfa2*CW*MW)/SW2**2*(Abb21*EL*gL(2)*mass(1)*mass(2)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(2)*(mass2(1) + 2*mass2(2) + MH2)) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(2)*(4*Abb2 + Abb1*(-2*S34 - 2*T + T14)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintU17)*
     -     (-(2*Alfa2*CW*MW*(Abb35*EL*gL(2))))) + 
     -     1/SW2**2*(Dval(dd23,iintU13)*
     -     (-(2*Alfa2*CW*MW*(Abb2*EL*gL(2)))))

      tmp2 = tmp2 + Dval(dd3,iintU17)*(-(3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb21*EL*gL(2)*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(2)*(Abb35 + Abb1*(2*T + 2*T14))) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(2)*(mass2(2) + MH2 + 5/2.D0*mass2(1)))) + 
     -     Dval(dd3,iintU18)*((3*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb21*EL*gL(1)*mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(1)*(Abb2 + Abb1*(-2*S34 - 2*T))) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(1)*(mass2(1) + MH2 + 5/2.D0*mass2(2)))) + 
     -     1/SW2**2*(Dval(dd33,iintU13)*
     -     (-(2*Alfa2*CW*MW*(Abb2*EL*gL(2))) )) + 
     -     1/(2.D0*CW*SW2)*(Dval(dd3,iintU15)*
     -     (2*Alfa2*MW*(Abb1*EL*gL(1)*MH2))) + 
     -     1/(2.D0*CW*SW2)*(Dval(dd3,iintU16)*
     -     (-2*Alfa2*MW*(Abb1*EL*gL(2)*MH2))) + 
     -     Dval(dd3,iintU14)*((3*Alfa2*CW*MW)/SW2**2*(Abb21*EL*gL(1)*
     -     mass(1)*mass(2)) + 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(1)*(2*mass2(1) + 4*mass2(2) + MH2)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(Abb2 + Abb1*(-T - 1/2.D0*S34))))

      Cloop(1) = tmp2/
     -     (mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - 
     -     MW2 - S - T24 - U) + Cloop(1)

      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNl part of the loop component for
*** qk = 1 and qk = 2 (qq and qbarq), with photons.

      subroutine boxesNl_p12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, tmp2

      integer qk

      double precision Alfa, Alfa2


** Paremeter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(3) + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)


      tmp2 = Dval(dd0,iintL10)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-3*mass2(1) - 3*mass2(2) - 4*mass2(3) - 
     -     mass2(4) - 3*MH2 + MW2)*q(4))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL3 + 
     -     Abb1*(-4*S - 2*S34 - 4*T - 2*T14 - 2*T24 - 4*U)))
     -     /SW) + Dval(dd0,iintL9)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup1*EL*q(3))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-2*S - 2*U)))/SW) + 
     -     Dval(dd0,iintL7)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(-mass2(1) - mass2(2) - 2*mass2(3) - 
     -     mass2(4) - MH2 + MW2)*q(4))/SW - 
     -     (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(4)*(2*AbbSumL3 + Abb1*(S + T24 + U)))/SW) - 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb1*EL*q(4)*Cval(cc0,iintL3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Cval(cc0,iintL4))/SW + 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb1*EL*q(3)*Cval(cc0,iintL5))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Cval(cc0,iintL6))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Dval(dd00,iintL10))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Dval(dd00,iintL9))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd0,iintL8)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3)))))/SW

      tmp2 = tmp2 + Dval(dd1,iintL10)*((Alfa2*MW)/SW2*
     -     (EL*q(4)*(3*AbbSumL12 + Abb1*(-4*S - 2*T - T24 - 4*U)))/
     -     SW + (3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(3*mass2(3) + MH2 + 4/3.D0*mass2(1) + 
     -     4/3.D0*mass2(2)))/SW) + 
     -     Dval(dd1,iintL7)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(4)*(2*AbbSumL3 + Abb1*(S + 3*T24 + U)))/SW + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*q(4)*(mass2(4) + MH2 + 1/3.D0*mass2(1) + 
     -     1/3.D0*mass2(2) + 2/3.D0*mass2(3)))/SW) + 
     -     Dval(dd1,iintL8)*(-(Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(S34 + T14 + T24)))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (Abb1*EL*q(3)*(mass2(1) + mass2(2) + MH2 - 
     -     1/2.D0*mass2(3) + 3/2.D0*mass2(4)))/SW) + 
     -     Dval(dd1,iintL9)*(-(3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(3)*(MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4)))/SW + 
     -     (3*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL21 + Abb1*(T24 + 1/3.D0*S34 + 1/3.D0*T14)))/
     -     SW)

      tmp2 = tmp2 + Dval(dd2,iintL10)*
     -     ((Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + 2*mass2(2) + 4*mass2(3) + 
     -     mass2(4))*q(4))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL3 + Abb1*(2*S + T24 + 2*U)))/SW) + 
     -     1/SW2*(Dval(dd12,iintL9)*
     -     (-(2*Alfa2*MW*(AbbSumL13*EL*q(3)))))/SW + 

     -     1/SW2*(Dval(dd12,iintL10)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*q(4))))/SW 

      tmp2 = tmp2 + Dval(dd2,iintL9)*
     -     (-((2*Alfa2*MW)/SW2*(Abb1*dup1*EL*q(3))/SW) - 
     -     (Alfa2*MW)/SW2*(EL*q(3)*(AbbSumL16 + Abb1*(-S - T - U)))/SW) + 
     -     Dval(dd2,iintL7)*((Alfa2*MW)/(2.D0*SW2)*(Abb8*EL*mass(3)*
     -     mass(4)*q(4))/SW + (Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + 
     -     3*mass2(4) + MH2)*q(4))/SW - (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(4)*(2*AbbSumL3 + Abb1*(S + T24 + U)))/SW) + 
     -     1/SW2*(Dval(dd23,iintL9)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3)))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd2,iintL8)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3))) - 
     -     3*Alfa2*MW*(Abb1*EL*mass2(3)*q(3)) - 
     -     Alfa2*MW*(Abb8*EL*mass(3)*mass(4)*q(3))))/SW + 
     -     1/SW2*(Dval(dd23,iintL10)*
     -     (-(2*Alfa2*MW*(AbbSumL3*EL*q(4)))))/SW

      tmp2 = tmp2 + Dval(dd3,iintL9)*
     -     ((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3))*q(3))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-2*S - 2*U)))/SW) + 
     -     Dval(dd3,iintL10)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + MH2)*q(4))/SW
     -     + (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL3 + Abb1*(-2*S - 2*T - 2*U)))/SW) + 
     -     Dval(dd3,iintL8)*(-((2*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(3))/SW) - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-S - T - U)))/SW) + 
     -     Dval(dd3,iintL7)*((Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(mass2(3) + mass2(4) + 1/2.D0*mass2(1) + 
     -     1/2.D0*mass2(2)))/SW - 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL3 + Abb1*(T24 + 1/2.D0*S + 1/2.D0*U)))/SW)

      Cloop(1) = tmp2/(-MW2 + T) + Cloop(1)

            
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNl part of the loop component for
*** qk = 3 and qk = 4 (qqbar and qbarqbar), with photons.

      subroutine boxesNl_p34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, tmp2

      integer qk

      double precision Alfa, Alfa2


** Paremeter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(3) + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)

      tmp2 = Dval(dd0,iintL7)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(7*mass2(1) + 7*mass2(2) + 6*mass2(3) + 
     -     3*mass2(4) + 7*MH2 + MW2)*q(4))/SW - 
     -     (Alfa2*MW)/(2.D0*SW2)*(EL*q(4)*(2*AbbSumL1 + 
     -     Abb1*(-3*S - 4*S34 - 4*T - 4*T14 - 3*T24 - 3*U)))
     -     /SW) + Dval(dd0,iintL9)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup1*EL*q(3))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-2*S - 2*U)))/SW) + 
     -     Dval(dd0,iintL10)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) - mass2(2) - 2*mass2(3) - MH2 + 
     -     MW2)*q(4))/SW + (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL1 + Abb1*(-2*S - 2*T - 2*U)))/SW) - 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb1*EL*q(4)*Cval(cc0,iintL3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Cval(cc0,iintL4))/SW + 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb1*EL*q(3)*Cval(cc0,iintL5))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Cval(cc0,iintL6))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(4)*Dval(dd00,iintL10))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(3)*Dval(dd00,iintL9))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd0,iintL8)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3)))))/SW

      tmp2 = tmp2 + Dval(dd1,iintL10)*((Alfa2*MW)/SW2*
     -     (EL*q(4)*(3*AbbSumL12 + Abb1*(-4*S - 2*T - T24 - 4*U)))/
     -     SW + (3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(3*mass2(3) + MH2 + 4/3.D0*mass2(1) + 
     -     4/3.D0*mass2(2)))/SW) + 
     -     Dval(dd1,iintL7)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(4)*(2*AbbSumL1 + Abb1*(S + 3*T24 + U)))/SW + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*q(4)*(mass2(4) + MH2 + 1/3.D0*mass2(1) + 
     -     1/3.D0*mass2(2) + 2/3.D0*mass2(3)))/SW) + 
     -     Dval(dd1,iintL8)*(-(Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(S34 + T14 + T24)))/SW + 
     -     (Alfa2*MW)/SW2*(Abb1*EL*q(3)*(mass2(1) + mass2(2) + MH2 - 
     -     1/2.D0*mass2(3) + 3/2.D0*mass2(4)))/SW) + 
     -     Dval(dd1,iintL9)*(-(3*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(3)*(MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4)))/SW + (3*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL21 + Abb1*(T24 + 1/3.D0*S34 + 1/3.D0*T14)))/
     -     SW)

      tmp2 = tmp2 + Dval(dd2,iintL10)*((Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + 2*mass2(2) + 4*mass2(3) + 
     -     mass2(4))*q(4))/SW - (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL1 + Abb1*(2*S + T24 + 2*U)))/SW) + 
     -     1/SW2*(Dval(dd12,iintL9)*
     -     (-(2*Alfa2*MW*(AbbSumL13*EL*q(3)))))/SW + 
     -     1/SW2*(Dval(dd12,iintL10)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*q(4))))/SW 

      tmp2 = tmp2 + Dval(dd2,iintL9)*
     -     (-((2*Alfa2*MW)/SW2*(Abb1*dup1*EL*q(3))/SW) - (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-S - T - U)))/SW) + 
     -     Dval(dd2,iintL7)*((Alfa2*MW)/(2.D0*SW2)*(Abb27*EL*mass(3)
     -     *mass(4)*q(4))/SW + (Alfa2*MW)/(2.D0*SW2)*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + 
     -     3*mass2(4) + MH2)*q(4))/SW - (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(4)*(2*AbbSumL1 + Abb1*(S + T24 + U)))/SW) + 
     -     1/SW2*(Dval(dd23,iintL9)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3)))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd2,iintL8)*
     -     (-(2*Alfa2*MW*(AbbSumL16*EL*q(3))) - 
     -     3*Alfa2*MW*(Abb1*EL*mass2(3)*q(3)) - 
     -     Alfa2*MW*(Abb27*EL*mass(3)*mass(4)*q(3))))/SW + 
     -     1/SW2*(Dval(dd23,iintL10)*
     -     (-(2*Alfa2*MW*(AbbSumL1*EL*q(4)))))/SW

      tmp2 = tmp2 + Dval(dd3,iintL9)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3))*q(3))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-2*S - 2*U)))/SW) + 
     -     Dval(dd3,iintL8)*(-((2*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(3))/SW) - (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumL16 + Abb1*(-S - T - U)))/SW) - 
     -     (2*Alfa2*MW)/SW2*(0*EL*mass(3)*q(3)*Dval(dd33,iintL9))/SW + 
     -     Dval(dd3,iintL10)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) + MH2)*q(4))/SW
     -     + (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL1 + Abb1*(-2*S - 2*T - 2*U)))/SW) + 
     -     Dval(dd3,iintL7)*((Alfa2*MW)/SW2*
     -     (Abb1*EL*q(4)*(mass2(3) + mass2(4) + 1/2.D0*mass2(1) + 
     -     1/2.D0*mass2(2)))/SW - (Alfa2*MW)/SW2*
     -     (EL*q(4)*(AbbSumL1 + Abb1*(T24 + 1/2.D0*S + 1/2.D0*U)))/SW)

      Cloop(1) = tmp2/(-MW2 + T) + Cloop(1)

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNl part of the loop component for
*** qk = 1 and qk = 2 (qq and qbarq), with no photons.

      subroutine boxesNl_nop12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup2, dup3, dup4, dup5, dup6, dup7, dup8, tmp3

      integer qk

      double precision Alfa, Alfa2


** Paremeter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup3 = 4*S + 2*T + T24 + 4*U
      dup4 = 3*mass2(3) + MH2 + 4/3.D0*mass2(1) + 4/3.D0*mass2(2)
      dup5 = MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4)
      dup6 = 2*AbbSumL3 + Abb1*(-S34 - 2*T - T14 + T24)
      dup7 = AbbSumL21 + Abb1*(T24 + 1/3.D0*S34 + 1/3.D0*T14)
      dup8 = AbbSumL16 + Abb1*(-S - U - 1/2.D0*T24)

      tmp3 = Dval(dd0,iintL20)*((Alfa2*MW)/(CW*SW2)*
     -     (Abb8*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (Abb1*EL*(gR(4)*mass2(4) + gL(4)*(MH2 + MW2/2.D0))) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*T24)) + 
     -     Dval(dd0,iintL21)*(-(2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + mass2(3) + 
     -     mass2(4))) - (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-2*S - 2*U)))) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(4)*Cval(cc0,iintL11)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Cval(cc0,iintL12)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Cval(cc0,iintL13)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*Cval(cc0,iintL14)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*Cval(cc0,iintL15)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(3)*Cval(cc0,iintL16)) + 
     -     Dval(dd0,iintL17)*(-((Alfa2*CW*MW)/SW2**2*(dup6*EL*gL(3))) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) - mass2(3) + MH2))) + 
     -     Dval(dd0,iintL19)*(-(Alfa2*MW)/(CW*SW2)*
     -     (Abb8*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(3)*(S34 + T14 + T24)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(-mass2(1) - mass2(2) - MH2 - 1/2.D0*mass2(3) - 
     -     3/2.D0*mass2(4))))) + 
     -     Dval(dd0,iintL18)*(-((2*Alfa2*CW*MW)/SW2**2*
     -     (dup8*EL*gL(4))) - (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + MW2 + 
     -     5/2.D0*mass2(3) + 1/2.D0*MH2)))

      tmp3 = tmp3 + Dval(dd1,iintL18)*((Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(2*mass2(1) + 2*mass2(2) - mass2(3) + 
     -     3*mass2(4) + MH2)) - (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(4)*(2*AbbSumL16 + Abb1*(S34 + T14 + T24)))) + 
     -     Dval(dd0,iintL22)*((2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*(-3*mass2(1) - 3*mass2(2) - 4*mass2(3) - 
     -     mass2(4) - 3*MH2 + MZ2)) - (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL3 + 
     -     Abb1*(-4*S - 2*S34 - 4*T - 2*T14 - 2*T24 - 4*U))))
     -     + (4*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*Dval(dd00,iintL17)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Dval(dd00,iintL18)) - 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*Dval(dd00,iintL21)) + 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*Dval(dd00,iintL22)) + 
     -     Dval(dd1,iintL17)*((5*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(MH2 + 2/5.D0*mass2(1) + 2/5.D0*mass2(2) + 
     -     1/5.D0*mass2(3) + 3/5.D0*mass2(4))) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumL3 + Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14))))

      tmp3 = tmp3 + Dval(dd1,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(-2*mass2(1) - 2*mass2(2) - 3*mass2(4) - 
     -     2*MH2))) - (Alfa2*MW)/(CW*SW2)*
     -     (EL*gL(3)*(AbbSumL16 + Abb1*(S34 + T14 + T24)))) + 
     -     Dval(dd1,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(4)*mass2(4) + 
     -     gL(4)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     2*mass2(4) + 3*MH2))) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (EL*gL(4)*(2*AbbSumL3 + Abb1*(S + 3*T24 + U)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintL21)*
     -     (3*Alfa2*CW*MW*(Abb1*dup5*EL*gL(4)) - 
     -     3*Alfa2*CW*MW*(dup7*EL*gL(4))))

      tmp3 = tmp3 + 1/SW2**2*(Dval(dd11,iintL17)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3)))) + 
     -     1/SW2**2*(Dval(dd11,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))) + 
     -     1/SW2**2*(Dval(dd12,iintL17)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3)))) + 
     -     1/SW2**2*(Dval(dd12,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintL22)*
     -     (-(3*Alfa2*CW*MW*(AbbSumL12*EL*gL(3))) + 
     -     Alfa2*CW*(Abb1*
     -     (MW*(dup3*EL*gL(3)) - 3*MW*(dup4*EL*gL(3))))))

      tmp3 = tmp3 + 1/SW2**2*(Dval(dd13,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL16*EL*gL(4))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintL22)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3))) )) + 
     -     1/SW2**2*(Dval(dd13,iintL17)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL3*EL*gL(3))) ))+ 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintL21)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))

      tmp3 = tmp3 + Dval(dd2,iintL21)*
     -     ((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-S - T - U))) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(3) + 1/2.D0*mass2(1) + 
     -     1/2.D0*mass2(2)))) + 
     -     Dval(dd2,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-2*gL(4) + gR(4))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*T24) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*(MH2 + 1/2.D0*mass2(3)))
     -     ) + Dval(dd2,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-2*gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(3)*(S34 + T14 + T24)) + 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) + mass2(3) + MH2 + 
     -     3/2.D0*mass2(4)))) + 
     -     1/SW2**2*(Dval(dd2,iintL18)*
     -     (-(3*Alfa2*CW*MW*(AbbSumL12*EL*gL(4))) + 
     -     Alfa2*CW*(Abb1*
     -     (MW*(dup3*EL*gL(4)) - 3*MW*(dup4*EL*gL(4)))))) + 
     -     1/SW2**2*(Dval(dd2,iintL17)*
     -     (3*Alfa2*CW*MW*(Abb1*dup5*EL*gL(3)) - 
     -     3*Alfa2*CW*MW*(dup7*EL*gL(3))))

      tmp3 = tmp3 + Dval(dd2,iintL22)*(-(Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*(2*mass2(1) + 2*mass2(2) + 4*mass2(3) + 
     -     mass2(4))) + (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL3 + Abb1*(2*S + T24 + 2*U)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintL21)*
     -     (2*Alfa2*CW*MW*(AbbSumL16*EL*gL(4)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintL22)*
     -     (2*Alfa2*CW*MW*(AbbSumL3*EL*gL(3)))) 

      tmp3 = tmp3 + Dval(dd3,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(3)*mass2(3) + gL(3)*(-mass2(4) + 2*MH2)))) + 
     -     Dval(dd3,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb8*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(4)*mass2(4) + 3*gL(4)*MH2)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(4)*T24)) + 
     -     Dval(dd3,iintL21)*(-(2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + mass2(3))) - 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-2*S - 2*U)))) + 
     -     Dval(dd3,iintL18)*(-((2*Alfa2*CW*MW)/SW2**2*
     -     (dup8*EL*gL(4))) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(2*mass2(1) + 2*mass2(2) + 7*mass2(3) + 
     -     MH2))) + 
     -     1/SW2**2*(Dval(dd3,iintL17)*
     -     (-(Alfa2*CW*MW*(Abb1*dup2*EL*gL(3))) - 
     -     Alfa2*CW*MW*(dup6*EL*gL(3))))

      tmp3 = tmp3 + Dval(dd3,iintL22)*
     -     (-((2*Alfa2*CW*MW)/(CW2*SW2**2)*(Abb1*dup2*EL*gL(3))) - 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL3 + Abb1*(-2*S - 2*T - 2*U)))) 

      Cloop(1) = tmp3/(-MW2 + T) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesNl part of the loop component for
*** qk = 3 and qk = 4 (qqbar and qbarqbar), with no photons.

      subroutine boxesNl_nop34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup2, dup3, dup4, dup5, dup6, dup7, dup8, tmp3

      integer qk

      double precision Alfa, Alfa2


** Paremeter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup3 = 4*S + 2*T + T24 + 4*U
      dup4 = 3*mass2(3) + MH2 + 4/3.D0*mass2(1) + 4/3.D0*mass2(2)
      dup5 = MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 2/3.D0*mass2(4)
      dup6 = 2*AbbSumL1 + Abb1*(-S34 - 2*T - T14 + T24)
      dup7 = AbbSumL21 + Abb1*(T24 + 1/3.D0*S34 + 1/3.D0*T14)
      dup8 = AbbSumL16 + Abb1*(-S - U - 1/2.D0*T24)

      tmp3 = Dval(dd0,iintL17)*
     -     (-((Alfa2*CW*MW)/SW2**2*(dup6*EL*gL(3))) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) - mass2(3) + MH2)))+
     -     Dval(dd0,iintL20)*((Alfa2*MW)/(CW*SW2)*
     -     (Abb27*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (Abb1*EL*(gR(4)*mass2(4) + gL(4)*(MH2 + MW2/2.D0))) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*T24)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(4)*Cval(cc0,iintL11)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Cval(cc0,iintL12)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Cval(cc0,iintL13)) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*Cval(cc0,iintL14)) - 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*Cval(cc0,iintL15)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(3)*Cval(cc0,iintL16)) + 
     -     Dval(dd0,iintL19)*(-(Alfa2*MW)/(CW*SW2)*
     -     (Abb27*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (Abb1*EL*gL(3)*(2*S + S34 + 2*T + T14 + T24 + 2*U)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(3*mass2(1) + 3*mass2(2) + 3*MH2 + 
     -     7/2.D0*mass2(3) + 1/2.D0*mass2(4))))) + 
     -     Dval(dd0,iintL18)*(-((2*Alfa2*CW*MW)/SW2**2*
     -     (dup8*EL*gL(4))) - (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + MW2 + 
     -     5/2.D0*mass2(3) + 1/2.D0*MH2)))

      tmp3 = tmp3 + Dval(dd1,iintL18)*((Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(2*mass2(1) + 2*mass2(2) - mass2(3) + 
     -     3*mass2(4) + MH2)) - (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(4)*(2*AbbSumL16 + Abb1*(S34 + T14 + T24)))) + 
     -     Dval(dd0,iintL21)*(-(2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + mass2(3) + 
     -     mass2(4))) - (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-2*S - 2*U)))) + 
     -     Dval(dd0,iintL22)*((2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*(-mass2(1) - mass2(2) - 2*mass2(3) - 
     -     MH2 + MZ2)) - (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL1 + Abb1*(-2*S - 2*T - 2*U)))) + 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(3)*Dval(dd00,iintL17)) - 
     -     (4*Alfa2*CW*MW)/SW2**2*(Abb1*EL*gL(4)*Dval(dd00,iintL18)) - 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*Dval(dd00,iintL21)) + 
     -     (4*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*Dval(dd00,iintL22)) + 
     -     Dval(dd1,iintL17)*((5*Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(3)*(MH2 + 2/5.D0*mass2(1) + 2/5.D0*mass2(2) + 
     -     1/5.D0*mass2(3) + 3/5.D0*mass2(4))) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumL1 + Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14))))

      tmp3 = tmp3 + Dval(dd1,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*(gR(3)*mass2(3) + 
     -     gL(3)*(-2*mass2(1) - 2*mass2(2) - 3*mass2(4) - 2*MH2))) - 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (EL*gL(3)*(AbbSumL16 + Abb1*(S34 + T14 + T24)))) + 
     -     Dval(dd1,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*(gR(4)*mass2(4) + 
     -     gL(4)*(mass2(1) + mass2(2) + 2*mass2(3) + 
     -     2*mass2(4) + 3*MH2))) - (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (EL*gL(4)*(2*AbbSumL1 + Abb1*(S + 3*T24 + U)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintL21)*
     -     (3*Alfa2*CW*MW*(Abb1*dup5*EL*gL(4)) - 
     -     3*Alfa2*CW*MW*(dup7*EL*gL(4))))

      tmp3 = tmp3 + 1/SW2**2*(Dval(dd11,iintL17)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3)))) + 
     -     1/SW2**2*(Dval(dd11,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))) + 
     -     1/SW2**2*(Dval(dd12,iintL17)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3)))) + 
     -     1/SW2**2*(Dval(dd12,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd1,iintL22)*
     -     (-(3*Alfa2*CW*MW*(AbbSumL12*EL*gL(3))) + Alfa2*CW*(Abb1*
     -     (MW*(dup3*EL*gL(3)) - 3*MW*(dup4*EL*gL(3))))))

      tmp3 = tmp3 + 1/SW2**2*(Dval(dd13,iintL18)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL16*EL*gL(4))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintL22)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL13*EL*gL(3))))) + 
     -     1/SW2**2*(Dval(dd13,iintL17)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL1*EL*gL(3))))) + 
     -     1/(CW2*SW2**2)*(Dval(dd12,iintL21)*
     -     (2*Alfa2*CW*MW*(AbbSumL13*EL*gL(4))))

      tmp3 = tmp3 + Dval(dd2,iintL21)*((Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-S - T - U))) + 
     -     (2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(3) + 1/2.D0*mass2(1) + 
     -     1/2.D0*mass2(2)))) + 
     -     Dval(dd2,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-2*gL(4) + gR(4))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*T24) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(4)*(MH2 + 1/2.D0*mass2(3)))
     -     ) + Dval(dd2,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-2*gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gL(3)*(S34 + T14 + T24)) + 
     -     (Alfa2*MW)/(CW*SW2)*
     -     (Abb1*EL*gL(3)*(mass2(1) + mass2(2) + mass2(3) + MH2 + 
     -     3/2.D0*mass2(4)))) + 1/SW2**2*(Dval(dd2,iintL17)*
     -     (3*Alfa2*CW*MW*(Abb1*dup5*EL*gL(3)) - 
     -     3*Alfa2*CW*MW*(dup7*EL*gL(3)))) + 1/SW2**2*
     -     (Dval(dd2,iintL18)*
     -     (-(3*Alfa2*CW*MW*(AbbSumL12*EL*gL(4))) + 
     -     Alfa2*CW*(Abb1*
     -     (MW*(dup3*EL*gL(4)) - 3*MW*(dup4*EL*gL(4))))))

      tmp3 = tmp3 + Dval(dd2,iintL22)*(-(Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(3)*(2*mass2(1) + 2*mass2(2) + 4*mass2(3) + 
     -     mass2(4))) + (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL1 + Abb1*(2*S + T24 + 2*U)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintL21)*
     -     (2*Alfa2*CW*MW*(AbbSumL16*EL*gL(4)))) + 
     -     1/(CW2*SW2**2)*(Dval(dd23,iintL22)*
     -     (2*Alfa2*CW*MW*(AbbSumL1*EL*gL(3)))) 

      tmp3 = tmp3 + Dval(dd3,iintL19)*(-(Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-gL(3) + gR(3))*mass(3)*mass(4)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(3)*mass2(3) + gL(3)*(-mass2(4) + 2*MH2)))) + 
     -     Dval(dd3,iintL20)*((Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb27*EL*(-gL(4) + gR(4))*mass(3)*mass(4)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*
     -     (Abb1*EL*(gR(4)*mass2(4) + 3*gL(4)*MH2)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gL(4)*T24)) + 
     -     Dval(dd3,iintL21)*(-(2*Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (Abb1*EL*gL(4)*(mass2(1) + mass2(2) + mass2(3))) - 
     -     (3*Alfa2*CW*MW)/(CW2*SW2**2)*(0*EL*gR(4)*mass(4)) - 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(4)*(AbbSumL16 + Abb1*(-2*S - 2*U)))) + 
     -     Dval(dd3,iintL18)*(-((2*Alfa2*CW*MW)/SW2**2*
     -     (dup8*EL*gL(4))) - (Alfa2*CW*MW)/SW2**2*
     -     (Abb1*EL*gL(4)*(2*mass2(1) + 2*mass2(2) + 7*mass2(3) + 
     -     MH2))) + 1/SW2**2*(Dval(dd3,iintL17)*
     -     (-(Alfa2*CW*MW*(Abb1*dup2*EL*gL(3))) - 
     -     Alfa2*CW*MW*(dup6*EL*gL(3))))

      tmp3 = tmp3 + Dval(dd3,iintL22)*
     -     (-((2*Alfa2*CW*MW)/(CW2*SW2**2)*(Abb1*dup2*EL*gL(3))) - 
     -     (Alfa2*CW*MW)/(CW2*SW2**2)*
     -     (EL*gL(3)*(AbbSumL1 + Abb1*(-2*S - 2*T - 2*U)))) 

      Cloop(1) = tmp3/(-MW2 + T) + Cloop(1)
      
      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesPl part of the loop component for
*** qk = 1 and qk = 2 (qq and qbarq)

      subroutine boxesPl_12_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex tmp1

      integer qk

      double precision Alfa, Alfa2


** Declaration of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      tmp1 = Dval(dd0,iintL2)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gR(3)*gR(4)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*gL(4)*
     -     (-mass2(1) - mass2(2) - mass2(3) + MZ2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + Abb1*(-S - U)))/SW) - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*Cval(cc0,iintL1))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*Dval(dd00,iintL2))/SW + 
     -     1/(CW2*SW2)*(Dval(dd12,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd11,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*gL(3)*gL(4))))/SW + 
     -     Dval(dd1,iintL2)*(-(6*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*
     -     (MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 1/3.D0*mass2(4)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + 
     -     Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW)

      tmp1 = tmp1 + Dval(dd2,iintL2)*
     -     ((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*(mass2(1) + mass2(2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + Abb1*(-S - U)))/SW) + 
     -     1/(CW2*SW2)*(Dval(dd23,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd13,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL39*EL*gL(3)*gL(4))))/SW + 
     -     Dval(dd3,iintL2)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(2*AbbSumL16 + 
     -     Abb1*(-S + S34 + T14 + T24 - U)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*(MH2 + 1/2.D0*mass2(4)))/SW)

      tmp1 = tmp1 + 1/(CW2*SW2)*(Dval(dd33,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW

      Cloop(1) = tmp1/(-MW2 + T) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesPl part of the loop component for
*** qk = 3 and qk = 4 (qqbar and qbarqbar)

      subroutine boxesPl_34_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex tmp1

      integer qk

      double precision Alfa, Alfa2


** Declaration of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      tmp1 = Dval(dd0,iintL2)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gR(3)*gR(4)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*gL(4)*
     -     (-mass2(1) - mass2(2) - mass2(3) + MZ2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + Abb1*(-S - U)))/SW) - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*Cval(cc0,iintL1))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*Dval(dd00,iintL2))/SW + 
     -     1/(CW2*SW2)*(Dval(dd12,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd11,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*gL(3)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd13,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL39*EL*gL(3)*gL(4))))/SW + 
     -     Dval(dd1,iintL2)*(-(6*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*
     -     (MH2 + 1/3.D0*mass2(1) + 1/3.D0*mass2(2) + 
     -     1/3.D0*mass2(3) + 1/3.D0*mass2(4)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + 
     -     Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW)

      tmp1 = tmp1 + Dval(dd2,iintL2)*
     -     ((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*(mass2(1) + mass2(2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(AbbSumL16 + Abb1*(-S - U)))/SW) + 
     -     1/(CW2*SW2)*(Dval(dd33,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd23,iintL2)*
     -     (2*Alfa2*MW*(AbbSumL16*EL*gL(3)*gL(4))))/SW + 
     -     Dval(dd3,iintL2)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(3)*gL(4)*(2*AbbSumL16 + 
     -     Abb1*(-S + S34 + T14 + T24 - U)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*gL(4)*(MH2 + 1/2.D0*mass2(4)))/SW)

      Cloop(1) = tmp1/(-MW2 + T) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesPu part of the loop component for
*** qk = 1 and qk = 3 (qq and qqbar)

      subroutine boxesPu_13_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex tmp3

      integer qk

      double precision Alfa, Alfa2


** Declaration of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      tmp3 = Dval(dd0,iintU20)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*(-(gL(1)*gL(2)) + gR(1)*gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2) + MH2 + MZ2))/
     -     SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb10 + Abb1*(T + T14)))/SW) - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*Cval(cc0,iintU19))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*Dval(dd00,iintU20))/SW + 
     -     1/(CW2*SW2)*(Dval(dd11,iintU20)*
     -     (2*Alfa2*MW*(AbbSumU8*EL*gL(1)*gL(2))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd12,iintU20)*
     -     (2*Alfa2*MW*(Abb10*EL*gL(1)*gL(2))))/SW + 
     -     Dval(dd1,iintU20)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(mass2(1) - mass2(2) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb10 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW)

      tmp3 = tmp3 + Dval(dd2,iintU20)*
     -     (-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(1)*gL(2)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(3*mass2(1) + mass2(2) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb10 + Abb1*(T + T14)))/SW) + 
     -     Dval(dd3,iintU20)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(1)*gL(2)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(2*Abb10 + Abb1*(-S34 + T + T14)))/SW) + 
     -     1/(CW2*SW2)*(Dval(dd23,iintU20)*
     -     (2*Alfa2*MW*(Abb10*EL*gL(1)*gL(2))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd13,iintU20)*
     -     (4*Alfa2*MW*(AbbSumU34*EL*gL(1)*gL(2))))/SW

      tmp3 = tmp3 + 1/(CW2*SW2)*
     -     (Dval(dd33,iintU20)*
     -     (2*Alfa2*MW*(Abb10*EL*gL(1)*gL(2))))/SW

      Cloop(1) = tmp3/
     -     (mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + MH2 - 
     -     MW2 - S - T24 - U) + Cloop(1)


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the boxesPu part of the loop component for
*** qk = 2 and qk = 4 (qbarq and qbarqbar)

      subroutine boxesPu_24_CC(qk,mass,mass2,S,T,T24,T14,S34,U,q,gL,
     -     gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass2(5), mass(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex tmp3

      integer qk

      double precision Alfa, Alfa2


** Declaration of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      tmp3 = Dval(dd0,iintU20)*((2*Alfa2*MW)/(CW2*SW2)*(Abb21*EL*
     -     (-(gL(1)*gL(2)) + gR(1)*gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2) + MH2 + MZ2))/
     -     SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb35 + Abb1*(T + T14)))/SW) - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*Cval(cc0,iintU19))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*Dval(dd00,iintU20))/SW + 
     -     1/(CW2*SW2)*(Dval(dd11,iintU20)*
     -     (2*Alfa2*MW*(AbbSumU7*EL*gL(1)*gL(2))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd12,iintU20)*
     -     (2*Alfa2*MW*(Abb35*EL*gL(1)*gL(2)) - 
     -     2*Alfa2*MW*(0*EL*gL(1)*gL(2)*mass(1))))/SW + 
     -     Dval(dd1,iintU20)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(mass2(1) - mass2(2) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb35 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW)

      tmp3 = tmp3 + Dval(dd2,iintU20)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(1)*gL(2)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(3*mass2(1) + mass2(2) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(Abb35 + Abb1*(T + T14)))/SW) + 
     -     Dval(dd3,iintU20)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(1)*gL(2)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(2)*(2*mass2(1) + mass2(2) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(2)*(2*Abb35 + Abb1*(-S34 + T + T14)))/SW) + 
     -     1/(CW2*SW2)*(Dval(dd23,iintU20)*
     -     (2*Alfa2*MW*(Abb35*EL*gL(1)*gL(2))))/SW + 
     -     1/(CW2*SW2)*(Dval(dd13,iintU20)*
     -     (4*Alfa2*MW*(AbbSumU34*EL*gL(1)*gL(2))))/SW

      tmp3 = tmp3 + 1/(CW2*SW2)*(Dval(dd33,iintU20)*
     -     (2*Alfa2*MW*(Abb35*EL*gL(1)*gL(2)) ))/SW

      Cloop(1) = tmp3/(mass2(1) + mass2(2) + 2*mass2(3) + mass2(4) + 
     -     MH2 - MW2 - S - T24 - U) + Cloop(1)


      end


********************************************************************************
