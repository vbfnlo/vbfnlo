********************************************************************************
********************************************************************************
***                                                                          ***
*** CC_pentagons.F                                                           ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contains all of the matrix elements for the electroweak        ***
*** pentagon diagrams for the charged currents.                              ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine calculates the pentagons with photons part of the loop 
*** component for the pattern q + q -> q + H + q (qk = 1)

      subroutine pentP_p_qq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
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
      double complex dup1, dup2

      double precision Alfa, Alfa2


** Declarations of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup2 = mass2(1) - mass2(2) + mass2(3)


      Cloop(1) = Cloop(1) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Dval(dd0,iintP1))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Dval(dd0,iintP2))/SW + 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Eval(ee00,iintP3))/SW + 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Eval(ee00,iintP4))/SW + 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Eval(ee00,iintP5))/SW + 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Eval(ee00,iintP6))/SW + 
     -     Eval(ee1,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(3*mass2(1) - 2*mass2(2) + mass2(3) + MH2)*
     -     q(1)*q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum67 + Abb1*(-S + T - T14)))/SW) + 
     -     Eval(ee1,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum49 + Abb1*(2*S + T14 + 2*T24)))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*
     -     (mass2(1) + mass2(3) + mass2(4) + MH2 + 
     -     1/2.D0*mass2(2)))/SW) + Alfa2/SW2*(Abb1*Eval(ee0,iintP3)*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + MW*(EL*q(1)*q(4)*(-2*S - 2*T - 
     -     2*T14))))/SW + Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP5)*
     -     (-(2*MW*(mass2(2)*q(2)*q(3))) - 2*MW*(mass2(3)*q(2)*q(3)) + 
     -     2*MW*(q(2)*q(3)*U)))/SW + Alfa2/SW2*(Abb1*Eval(ee0,iintP4)*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))/SW + 
     -     1/SW2*(Eval(ee0,iintP6)*
     -     (Alfa2*(Abb1*EL*
     -     (2*MW*(mass2(1)*q(1)*q(3)) + 2*MW*(mass2(3)*q(1)*q(3)) - 
     -     2*MW*(q(1)*q(3)*S)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP4)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 2*mass2(3) + MH2)*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee1,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum4 + Abb1*(2*S + T24)))/SW) + 
     -     Eval(ee1,iintP5)*(-((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(2)*q(3))/SW) - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum19 + Abb1*(-2*S34 - 2*T - 2*U)))/SW) + 
     -     1/SW2*(Eval(ee11,iintP3)*
     -     (4*Alfa2*MW*(Abb23*EL*q(1)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(4))) + 4*MW*(q(1)*q(4)*T)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP4)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum11 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee12,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum29 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee11,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee11,iintP5)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 2*mass2(3) + MH2)*q(2)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum10 + Abb1*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(AbbSum13 + Abb1*S))/SW) + 
     -     Eval(ee13,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(Abb23 + Abb1*T))/SW) + 
     -     Eval(ee12,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum11 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee12,iintP5)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(2)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum29 + Abb1*(-S + T)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP5)*
     -     (4*Alfa2*MW*(AbbSum32*EL*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(1)*q(2)*q(3))) - 
     -     4*MW*(mass2(2)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*T)))))/SW+
     -     1/SW2*(Eval(ee13,iintP4)*
     -     (4*Alfa2*MW*(AbbSum52*EL*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(1)*q(2)*q(4))) - 
     -     4*MW*(mass2(2)*q(2)*q(4)) + 4*MW*(q(2)*q(4)*T)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee14,iintP5)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 3*mass2(3) + MH2)*q(2)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee14,iintP4)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(3) + MH2)*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee2,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(3*mass2(1) + 3*mass2(3) + MH2)*q(1)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum74 + Abb1*(2*S + S34 + T24)))/SW) + 
     -     1/SW2*(Eval(ee14,iintP3)*
     -     (4*Alfa2*MW*(AbbSum72*EL*q(1)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(4))) + 4*MW*(q(1)*q(4)*T)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum14 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee2,iintP4)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(3) + mass2(4) + MH2)*q(2)*
     -     q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum56 + Abb1*(-S + T - T14 - T24)))/SW)+
     -     Eval(ee22,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee2,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum36 + Abb1*(2*T + U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP6)*
     -     (4*Alfa2*MW*(Abb23*EL*q(1)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*T)))))/SW
     -     + 1/SW2*(Eval(ee22,iintP5)*
     -     (4*Alfa2*MW*(Abb23*EL*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(3))) + 4*MW*(q(2)*q(3)*T)))))/SW
     -     + 1/SW2*(Eval(ee22,iintP4)*
     -     (4*Alfa2*MW*(Abb23*EL*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(4))) + 4*MW*(q(2)*q(4)*T)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee23,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(AbbSum78 + Abb1*S))/SW) + 
     -     Eval(ee23,iintP6)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(Abb23 + Abb1*T))/SW) + 
     -     Eval(ee24,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(3) + MH2)*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/SW2*(Eval(ee23,iintP5)*
     -     (4*Alfa2*MW*(Abb23*EL*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(1)*q(2)*q(3))) - 
     -     4*MW*(mass2(2)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*T)))))/SW+
     -     1/SW2*(Eval(ee24,iintP4)*
     -     (4*Alfa2*MW*(AbbSum61*EL*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(4))) + 4*MW*(q(2)*q(4)*T)))))/SW
     -     + 1/SW2*(Eval(ee23,iintP4)*
     -     (4*Alfa2*MW*(Abb23*EL*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(1)*q(2)*q(4))) - 
     -     4*MW*(mass2(2)*q(2)*q(4)) + 4*MW*(q(2)*q(4)*T)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*mass2(1)*q(1)*q(3)*
     -     Eval(ee33,iintP6))/SW + 
     -     Eval(ee3,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*mass(1)*mass(2)*q(1)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(5*mass2(1) + mass2(2) + mass2(3) + MH2)*q(1)*
     -     q(4))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*(S + T + T14))/SW) + 
     -     Eval(ee3,iintP5)*(-((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-2*mass2(2) + mass2(3))*q(2)*q(3))/SW) + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW) + 
     -     1/SW2*(Eval(ee33,iintP3)*
     -     (4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee24,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum17*EL*q(1)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee24,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum42*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(3)) + 4*MW*(mass2(3)*q(2)*q(3)) - 
     -     4*MW*(q(2)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee3,iintP4)*
     -     (-2*Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(2)*q(4)) + 
     -     Alfa2*(Abb1*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) - mass2(4) + 
     -     MH2)*q(2)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*(S34 + T + U))/SW) + 
     -     Eval(ee4,iintP6)*(-(2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2*(Eval(ee34,iintP5)*
     -     (4*Alfa2*MW*(AbbSum45*EL*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee34,iintP3)*
     -     (4*Alfa2*MW*(AbbSum83*EL*q(1)*q(4)) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee34,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum18*EL*q(1)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee34,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSum63*EL*q(2)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(1)*q(2)*q(4))) - 
     -     4*MW*(mass2(2)*q(2)*q(4)) + 4*MW*(q(2)*q(4)*T)))))/SW+
     -     1/SW2*(Eval(ee4,iintP3)*
     -     (Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S34 - 2*T24 - 2*U))))))/SW + 
     -     1/SW2*(Eval(ee4,iintP5)*
     -     (Alfa2*(Abb1*EL*
     -     (2*MW*(mass2(2)*q(2)*q(3)) - 2*MW*(q(2)*q(3)*U)))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*mass2(3)*q(1)*q(3)*
     -     Eval(ee44,iintP6))/SW + 
     -     1/SW2*(Eval(ee44,iintP5)*
     -     (4*Alfa2*MW*(Abb1*EL*mass2(3)*q(2)*q(3))))/SW

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagons with no photons part of the loop 
*** component for the pattern q + q -> q + H + q (qk = 1)

      subroutine pentP_nop_qq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
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
      double complex dup3, dup4, dup5, dup6, dup7, dup8, dup9

      double precision Alfa, Alfa2


** Declarations of parameters etc.
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup3 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup4 = mass2(1) - mass2(2) + mass2(3)
      dup5 = MH2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(2) + 3/2.D0*mass2(3)
      dup6 = 2*mass2(3) + MH2 + 3/2.D0*mass2(1) + 1/2.D0*mass2(2)
      dup7 = 2*mass2(1) + mass2(2) + 2*mass2(3) + MH2
      dup8 = gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) - 2*mass2(2) - mass2(3) - MH2)
      dup9 = AbbSum97 + Abb1*(T24 + 1/2.D0*S + 1/2.D0*T)

      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(1)*gR(3)*mass2(3) + 
     -     gL(3)*(gR(1)*mass2(1) + gL(1)*(mass2(3) + MH2))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum88 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee0,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-3*mass2(1) - 2*mass2(2) - 5*mass2(3) - 
     -     mass2(4) - 3*MH2)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum154 + 
     -     Abb1*(2*S + S34 + T + T14 + 2*T24 + U)))/SW) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Dval(dd0,iintP10))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(3)*Dval(dd0,iintP7))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(4)*Dval(dd0,iintP8))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Dval(dd0,iintP9))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP14)*((2*Alfa2*MW)/SW2*
     -     (AbbSum175*EL*gL(1)*gL(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(-2*mass2(1) + mass2(2)))/SW)+
     -     Eval(ee0,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-4*mass2(1) + mass2(2) - mass2(3) - MH2)))/
     -     SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee0,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(3)*mass2(3)) + gL(3)*(-4*mass2(1) - 2*mass2(2) - 
     -     4*mass2(3) - 4*MH2))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum116 + 
     -     Abb1*(2*S + 2*S34 + 2*T + 2*T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(2*mass2(1) - mass2(2) + mass2(3)))/
     -     SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP11))/SW + 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP12))/SW + 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP13))/SW + 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP14))/SW + 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP15))/SW + 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP16))/SW + 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP17))/SW + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP18)*
     -     (Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP16)*
     -     (Alfa2*(Abb14*EL*
     -     (-(2*MW*(gL(1)*gL(3)*mass(1)*mass(2))) + 
     -     2*MW*(gL(3)*gR(1)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup8*EL*gL(3)) + 
     -     MW*(EL*gL(1)*gL(3)*(2*S34 + 2*T + 2*U))))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP18))/SW + 
     -     Eval(ee1,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-mass2(1) - mass2(2) - 5*mass2(3) - 3*MH2)))
     -     /SW - (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum159 + 
     -     Abb1*(T24 + 1/3.D0*S + 1/3.D0*T)))/SW) + 
     -     Eval(ee1,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(mass2(1) + mass2(2) + 3*mass2(3) + 3*MH2)))/
     -     SW - (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum93 + 
     -     Abb1*(T24 + 1/3.D0*S + 1/3.D0*T)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-4*mass2(1) + 2*mass2(2) - mass2(3) - MH2)))
     -     /SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum67 + Abb1*(-S + T - T14)))/SW) + 
     -     Eval(ee1,iintP14)*((2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(4)*mass2(4) + 
     -     gL(4)*(-2*mass2(1) - mass2(4))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum178 + 
     -     Abb1*(-S + S34 - T - T14 + T24)))/SW) + 
     -     Eval(ee1,iintP13)*((2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(-2*mass2(1) + 2*mass2(3) - 2*mass2(4) + 
     -     2*MH2)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(3*AbbSum123 + 
     -     Abb1*(-2*T14 + 2*T24 - U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(2*mass2(1) + 2*mass2(3) + mass2(4) + 
     -     2*MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum49 + Abb1*(2*S + T14 + 2*T24)))/SW)
     -     + Eval(ee1,iintP17)*
     -     ((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*gR(1)*mass2(1) + 
     -     gL(1)*(-(gR(4)*mass2(4)) + 
     -     gL(4)*(-2*mass2(1) - 3*mass2(2) - mass2(3) - 
     -     2*MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum19 + Abb1*(-2*S34 - 2*T - 2*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP18)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-3*mass2(1) - 4*mass2(3) - MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum4 + Abb1*(2*S + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP15)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW + 1/SW2*(Eval(ee11,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum95*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (4*MW*mass2(3) + 4*MW*MH2 - 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP13)*
     -     (-(4*Alfa2*MW*(AbbSum126*EL*gL(2)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (4*MW*mass2(3) + 4*MW*MH2 - 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum126*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (4*MW*mass2(3) + 4*MW*MH2 - 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP12)*
     -     (-(4*Alfa2*MW*(AbbSum95*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (4*MW*mass2(3) + 4*MW*MH2 - 4*MW*T24))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP16)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee11,iintP17)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee11,iintP18)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/SW2*(Eval(ee12,iintP11)*
     -     (8*Alfa2*MW*(Abb1*dup5*EL*gL(1)*gL(3)) - 
     -     8*Alfa2*MW*(dup9*EL*gL(1)*gL(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup4*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum29 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee12,iintP14)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(1) + mass2(2) + mass2(3)))/SW+ 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum180 + Abb1*(S + T)))/SW) + 
     -     1/SW2*(Eval(ee12,iintP12)*
     -     (8*Alfa2*MW*(Abb1*dup5*EL*gL(2)*gL(4)) - 
     -     8*Alfa2*MW*(dup9*EL*gL(2)*gL(4))))/SW + 
     -     Eval(ee12,iintP13)*
     -     ((8*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(2)*gL(3))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum129 + 
     -     Abb1*(T24 + 1/2.D0*S + 1/2.D0*T)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP16)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup4*EL*gL(1)*gL(3))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum11 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee12,iintP18)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup4*EL*gL(2)*gL(4))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum11 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee12,iintP17)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup4*EL*gL(1)*gL(4))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum29 + Abb1*(-S + T)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP11)*
     -     (-(8*Alfa2*MW*(AbbSum98*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (8*MW*mass2(3) + 8*MW*MH2 - 8*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP15)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(-mass2(1) + mass2(2)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*(EL*gL(2)*gL(3)*(Abb23 + Abb1*T))/SW
     -     ) + 1/SW2*(Eval(ee13,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum181*EL*gL(1)*gL(4)))))/SW + 
     -     Eval(ee13,iintP13)*
     -     ((8*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(2)*gL(3))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum132 + 
     -     Abb1*(T24 + 1/2.D0*S + 1/2.D0*T)))/SW) + 
     -     Eval(ee13,iintP12)*
     -     ((8*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(2)*gL(4))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum161 + 
     -     Abb1*(T24 + 1/2.D0*S + 1/2.D0*T)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee13,iintP16)*
     -     (4*Alfa2*MW*(AbbSum52*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW
     -     + 1/(CW2*SW2)*(Eval(ee13,iintP17)*
     -     (4*Alfa2*MW*(AbbSum32*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW
      
      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP18)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum13 + Abb1*S))/SW) + 
     -     Eval(ee14,iintP12)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(AbbSum162 - Abb1*T24))/SW
     -     ) + Eval(ee14,iintP11)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup3*EL*gL(1)*gL(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum100 + Abb1*(S + T + T24)))/SW) + 
     -     Eval(ee14,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup3*EL*gL(2)*gL(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum134 + Abb1*(S + T + T24)))/SW) + 
     -     1/SW2*(Eval(ee14,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum182*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*MH2 - 4*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP16)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee14,iintP17)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee14,iintP18)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee14,iintP15)*
     -     (4*Alfa2*MW*(AbbSum72*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-4*mass2(1) - 2*mass2(2) - 8*mass2(3) - 
     -     4*MH2))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(3*AbbSum164 + 
     -     Abb1*(3*S + S34 + T + 3*T24 + U)))/SW) + 
     -     Eval(ee2,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(mass2(1) + 4*mass2(3) + 3*MH2))))/SW - 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum102 + Abb1*(T24 + 2/3.D0*S)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(-6*mass2(1) - 3*mass2(2) - 9*mass2(3) - 
     -     6*MH2))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(3*AbbSum136 + 
     -     Abb1*(4*S + 2*S34 + 2*T + 4*T24 + 2*U)))/SW) + 
     -     Eval(ee2,iintP14)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-4*mass2(1) - 2*mass2(2) - mass2(3) - MH2)))
     -     /SW + (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum184 + 
     -     Abb1*(T + 1/3.D0*S + 1/3.D0*T14)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - mass2(2) + mass2(3) + mass2(4) + 
     -     MH2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum56 + Abb1*(-S + T - T14 - T24)))/SW
     -     ) + Eval(ee2,iintP15)*
     -     (-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(-4*mass2(1) - 3*mass2(3) + mass2(4) - 
     -     MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum74 + Abb1*(2*S + S34 + T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP18)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) + mass2(2) - mass2(3))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum14 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee22,iintP11)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee2,iintP17)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-2*mass2(1) - mass2(2) - mass2(3))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum36 + Abb1*(2*T + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee22,iintP13)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee22,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee22,iintP12)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP16)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW + 1/SW2*(Eval(ee22,iintP14)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee23,iintP11)*
     -     ((8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum105 + Abb1*(T24 + 1/2.D0*S)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP17)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW + 1/(CW2*SW2)*
     -     (Eval(ee22,iintP18)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW + 1/SW2*(Eval(ee23,iintP14)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW
     -     + 1/SW2*(Eval(ee23,iintP13)*
     -     (-(8*Alfa2*MW*(AbbSum142*EL*gL(2)*gL(3))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup6*EL*gL(2)*gL(3)) + 
     -     MW*(EL*gL(2)*gL(3)*(-8*S - 4*T - 8*T24))))))/SW + 
     -     1/SW2*(Eval(ee23,iintP12)*
     -     (-(8*Alfa2*MW*(AbbSum167*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup6*EL*gL(2)*gL(4)) + 
     -     MW*(EL*gL(2)*gL(4)*(-8*S - 4*T - 8*T24))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee23,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(2*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum78 + Abb1*S))/SW) + 
     -     Eval(ee23,iintP18)*
     -     (-(Abb1*EL*gL(2)*gL(4)*(-mass2(1) + mass2(2)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*(EL*gL(2)*gL(4)*(Abb23 + Abb1*T))/SW
     -     ) + Eval(ee24,iintP11)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee24,iintP12)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP16)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW
     -     + 1/(CW2*SW2)*(Eval(ee23,iintP17)*
     -     (4*Alfa2*MW*(Abb23*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee24,iintP13)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee24,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum10 + Abb1*(S + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP16)*
     -     (4*Alfa2*MW*(AbbSum61*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW + 1/SW2*(Eval(ee24,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum187*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*T)))
     -     )/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(3*mass2(3) + 3*MH2)))/SW - 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum107 + Abb1*(T24 + 1/3.D0*S)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP17)*
     -     (-(4*Alfa2*MW*(AbbSum42*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP18)*
     -     (-(4*Alfa2*MW*(AbbSum17*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*gR(2)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-5*mass2(1) - 2*mass2(2) - 7*mass2(3) - 
     -     3*MH2)))/SW - 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum171 + Abb1*(S + T24 + 2/3.D0*T)))/SW)
     -     + Eval(ee3,iintP13)*
     -     (-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(12*mass2(1) + 10*mass2(2) + 12*mass2(3) + 
     -     2*mass2(4) + 10*MH2)))/SW - 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum146 + 
     -     Abb1*(2*S + 2*T + 2*T24 + U + 4/3.D0*S34 + 
     -     2/3.D0*T14)))/SW) + 
     -     1/SW2*(Eval(ee3,iintP14)*
     -     (2*Alfa2*MW*(AbbSum175*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(Abb14*EL*gL(4)*gR(1)*mass(1)*mass(2)) + 
     -     2*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP18)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*gR(2)*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(3*mass2(1) + mass2(3)))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(2)*gL(4)*S)/SW) + 
     -     Eval(ee3,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*
     -     (6*mass2(1) + mass2(2) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee3,iintP17)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*gR(1)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(-mass2(2) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     Eval(ee3,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*
     -     (mass2(1) + 2*mass2(2) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee34,iintP11)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum111 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee34,iintP12)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum173 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee33,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum149 + Abb1*(S + T + T24)))/SW) + 
     -     Eval(ee33,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum172 + Abb1*(S + T + T24)))/SW) + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*mass2(1)*Eval(ee33,iintP18))/SW + 
     -     1/(CW2*SW2)*(Eval(ee33,iintP15)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(1))))/SW + 
     -     1/SW2*(Eval(ee33,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum109*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (4*MW*mass2(3) + 4*MW*MH2 - 4*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee34,iintP13)*
     -     ((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*
     -     (2*mass2(1) + mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum151 + Abb1*(-S - T - T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP17)*
     -     (4*Alfa2*MW*(AbbSum45*EL*gL(1)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP15)*
     -     (4*Alfa2*MW*(AbbSum83*EL*gL(2)*gL(3)) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(1))))/SW + 
     -     1/SW2*(Eval(ee34,iintP14)*
     -     (4*Alfa2*MW*(AbbSum189*EL*gL(1)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP18)*
     -     (-(4*Alfa2*MW*(AbbSum18*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP16)*
     -     (-(4*Alfa2*MW*(AbbSum63*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(2) + 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + 5*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum154 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee4,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + gL(1)*(mass2(3) + MH2)))/SW
     -     - (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum88 + Abb1*(S + T24)))/SW) + 
     -     Eval(ee4,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb13*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-4*mass2(1) - 4*mass2(2) - 4*mass2(3) - 
     -     4*MH2)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum116 + 
     -     Abb1*(2*S + 2*S34 + 2*T + 2*T24 + 3*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(2)*gR(3)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(2) - mass2(3) + 
     -     mass2(4) - MH2)))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S34 + T24 + U))/SW) + 
     -     Eval(ee4,iintP14)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(4)*mass2(4) + 
     -     gL(4)*(-2*mass2(1) + mass2(2) - mass2(4))))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(2*AbbSum191 + 
     -     Abb1*(-S + S34 - T - T14 + T24 + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP16)*
     -     (2*Alfa2*MW*(Abb13*EL*gL(1)*gR(3)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb14*EL*
     -     (-(2*MW*(gL(1)*gL(3)*mass(1)*mass(2))) + 
     -     2*MW*(gL(3)*gR(1)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup8*EL*gL(3)) + 
     -     MW*(EL*gL(1)*gL(3)*(2*S34 + 2*T + 2*U))))))/SW

      Cloop(1) = Cloop(1) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(4)*mass2(3)*
     -     Eval(ee44,iintP12))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*mass2(3)*Eval(ee44,iintP18))/SW + 
     -     Eval(ee4,iintP18)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb13*EL*gL(2)*gR(4)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(2)*gL(4)*S)/SW + 
     -     (6*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(3) + 1/3.D0*mass2(1)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP17)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(3))))/SW + 
     -     1/SW2*(Eval(ee44,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum193*EL*gL(1)*gL(4))) ))/SW + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP17)*
     -     (-2*Alfa2*MW*(Abb13*EL*gL(1)*gR(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (2*MW*mass2(2) + 2*MW*mass2(3) - 2*MW*U))))/SW


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagon part of the loop component for the
*** pattern: qbar + q -> qbar + H + q
       
      subroutine pentP_p_qbarq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
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
      double complex dup1


      double precision Alfa, Alfa2

** Parameter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

     
      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(1) + mass2(2) + mass2(3) + MH2


      Cloop(1) = Cloop(1) + Eval(ee1,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(3) + MH2)*q(1)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum69 + Abb1*(-S - T - T14)))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Dval(dd0,iintP1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Dval(dd0,iintP2))/SW - 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Eval(ee00,iintP3))/SW - 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Eval(ee00,iintP4))/SW - 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Eval(ee00,iintP5))/SW - 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Eval(ee00,iintP6))/SW + 
     -     Eval(ee1,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(3*mass2(2) + 3*mass2(3) + MH2)*q(2)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum48 + Abb1*(T14 + T24 + 2*U)))/SW) + 
     -     Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP6)*
     -     (2*MW*(mass2(1)*q(1)*q(3)) + 2*MW*(mass2(3)*q(1)*q(3)) - 
     -     2*MW*(q(1)*q(3)*S)))/SW + 
     -     Alfa2/SW2*(Abb1*Eval(ee0,iintP3)*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S - 2*T - 2*T14))))/SW + 
     -     Alfa2/SW2*(Abb1*Eval(ee0,iintP4)*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))/SW + 
     -     1/SW2*(Eval(ee0,iintP5)*
     -     (Alfa2*(Abb1*EL*(-(2*MW*(mass2(2)*q(2)*q(3))) - 
     -     2*MW*(mass2(3)*q(2)*q(3)) + 2*MW*(q(2)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP5)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 2*mass2(3) + MH2)*q(2)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee11,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 2*mass2(3) + MH2)*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee1,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + 3*mass2(3) + MH2)*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum33 + Abb1*(T24 + 2*U)))/SW) + 
     -     1/SW2*(Eval(ee11,iintP3)*
     -     (4*Alfa2*MW*(Abb29*EL*q(1)*q(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee1,iintP6)*
     -     (2*Alfa2*MW*(AbbSum1*EL*q(1)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(2)*q(4)*(AbbSum13 + Abb1*U))/SW) + 
     -     Eval(ee12,iintP5)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(2)*q(3)*(AbbSum40 + Abb1*U))/SW) + 
     -     Eval(ee12,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(AbbSum40 + Abb1*U))/SW) + 
     -     Eval(ee11,iintP6)*(-(2*Alfa2*MW)/SW2*((-0)*EL*mass(2)*
     -     mass(4)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 2*mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum9 + Abb1*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP5)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(2)*q(3)*(AbbSum42 - Abb1*U))/SW) + 
     -     Eval(ee13,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(2)*q(4)*(AbbSum53 - Abb1*U))/SW) + 
     -     Eval(ee12,iintP6)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(AbbSum13 + Abb1*U))/SW) + 
     -     1/SW2*(Eval(ee13,iintP6)*
     -     (4*Alfa2*MW*(AbbSum16*EL*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee13,iintP3)*
     -     (4*Alfa2*MW*(Abb29*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee14,iintP3)*
     -     (-(4*Alfa2*MW*(AbbSum74*EL*q(1)*q(4))) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(1)*q(4))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + 4*mass2(2) + mass2(3) + MH2)*q(2)*
     -     q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum59 + Abb1*(S34 + U)))/SW) + 
     -     Eval(ee14,iintP6)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 3*mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee14,iintP5)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 3*mass2(3) + MH2)*q(2)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee14,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee2,iintP3)*((2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum76 + Abb1*(S34 + 2*T24 + 2*U)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*
     -     (mass2(2) + mass2(3) + mass2(4) + MH2 + 
     -     1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) - mass2(2) + mass2(3))*q(1)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum20 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee2,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*q(2)*q(3)*(AbbSum43 + Abb1*U))/SW) + 
     -     Eval(ee22,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + 2*mass2(3) + MH2)*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP6)*
     -     (4*Alfa2*MW*(Abb29*EL*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee22,iintP5)*
     -     (4*Alfa2*MW*(Abb29*EL*q(2)*q(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee23,iintP5)*
     -     (4*Alfa2*MW*(Abb29*EL*q(2)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee23,iintP3)*
     -     (4*Alfa2*MW*(AbbSum79*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee22,iintP4)*
     -     (4*Alfa2*MW*(Abb29*EL*q(2)*q(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee23,iintP4)*
     -     (4*Alfa2*MW*(Abb29*EL*q(2)*q(4)) - 
     -     8*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee24,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/SW2*(Eval(ee23,iintP6)*
     -     (4*Alfa2*MW*(Abb29*EL*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee24,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSum63*EL*q(2)*q(4))) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee3,iintP3)*
     -     (2*Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(1)*q(4)) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S - 2*T - 2*T14))))))/SW + 
     -     1/SW2*(Eval(ee24,iintP6)*
     -     (4*Alfa2*MW*(AbbSum27*EL*q(1)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW+
     -     1/SW2*(Eval(ee24,iintP5)*
     -     (4*Alfa2*MW*(AbbSum45*EL*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(3))) - 
     -     4*MW*(mass2(3)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW) + 
     -     Eval(ee3,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW) - (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*mass2(2)*q(2)*q(3)*Eval(ee33,iintP5))/SW + 
     -     Eval(ee3,iintP4)*(-((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*mass(1)*mass(2)*q(2)*q(4))/SW) - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 5*mass2(2) + mass2(3) + MH2)*q(2)*
     -     q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*(S34 + T + U))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*mass(2)*q(2)*q(4)*(0 - 1/2.D0*0))/SW) + 
     -     1/SW2*(Eval(ee34,iintP6)*
     -     (4*Alfa2*MW*(AbbSum30*EL*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee34,iintP3)*
     -     (-(4*Alfa2*MW*(AbbSum83*EL*q(1)*q(4)))))/SW + 
     -     1/SW2*(Eval(ee33,iintP4)*
     -     (-4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee34,iintP4)*
     -     (4*Alfa2*MW*(AbbSum65*EL*q(2)*q(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee34,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum46*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(3))) - 
     -     4*MW*(mass2(3)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*mass2(3)*q(2)*q(3)*
     -     Eval(ee44,iintP5))/SW + 
     -     Eval(ee4,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) - mass2(4) + 
     -     MH2)*q(1)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*(S + T + T14))/SW) + 
     -     Eval(ee4,iintP6)*(-(2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     Eval(ee4,iintP5)*((2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(3)*(mass2(3) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee44,iintP6)*
     -     (-4*Alfa2*MW*(Abb1*EL*mass2(3)*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee4,iintP4)*
     -     (Alfa2*(Abb1*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S + 2*T14 + 2*T24))))))/SW
     

      end


*******************************************************************************
*******************************************************************************
   
*** This subroutine calculates the pentagon part of the loop component for the
*** pattern: qbar + q -> qbar + H + q
       
      subroutine pentP_nop_qbarq_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,
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
      double complex dup2, dup3, dup4, dup5, dup6


      double precision Alfa, Alfa2

** Parameter etc. declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

     
      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = 2*mass2(1) + 2*mass2(2) - mass2(3) + mass2(4)
      dup3 = gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(2) - mass2(3) - MH2)
      dup4 = gR(3)*mass2(3) + gL(3)*(3*mass2(3) + 3*MH2)
      dup5 = 2*AbbSum101 + Abb1*(-S - S34 - T - T14 + T24)
      dup6 = AbbSum89 + Abb1*(S + 2*T24 + 2*U)


      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP10)*((2*Alfa2*MW)/SW2*(dup6*EL*gL(1)*gL(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(-2*mass2(1) - 2*mass2(2) - 4*mass2(3) - 
     -     2*MH2 + 2*MZ2))))/SW) + 
     -     Eval(ee0,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-mass2(1) - 2*mass2(3) - mass2(4) - 2*MH2)))
     -     /SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(Abb34 + Abb1*(S34 + T - T14 + T24 + U)))/
     -     SW) + (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*Dval(dd0,iintP1))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Dval(dd0,iintP7))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(4)*Dval(dd0,iintP8))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Dval(dd0,iintP9))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(2)*gR(3)*mass2(3) + 
     -     gL(3)*(gR(2)*mass2(2) + gL(2)*(mass2(3) + MH2))))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum127 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee0,iintP13)*((2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum178 + 
     -     Abb1*(2*mass2(1) + 2*mass2(2) + 2*mass2(3) + 
     -     mass2(4) + 2*MH2)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*
     -     (mass2(1) + mass2(3) + MH2 + 3/2.D0*mass2(2) + 
     -     1/2.D0*mass2(4)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP14)*
     -     (Alfa2*(Abb21*EL*
     -     (2*MW*(gL(2)*gL(3)*mass(1)*mass(2)) - 
     -     2*MW*(gL(3)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*
     -     (-(2*MW*(dup3*EL*gL(3))) + 
     -     MW*(EL*gL(2)*gL(3)*(-2*S - 2*T - 2*T14))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee0,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - 4*mass2(2) - mass2(3) - MH2)))/SW
     -     + (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW) - 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP10))/SW - 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP11))/SW - 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP12))/SW - 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP13))/SW - 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP14))/SW - 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP15))/SW - 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP16))/SW + 
     -     Eval(ee0,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*
     -     gL(2)*gL(4)*S)/SW + (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*
     -     (MW2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP16)*
     -     (Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(2*MW*mass2(2)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(3*mass2(1) + mass2(2) - 2*mass2(3) - 2*MH2))
     -     )/SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(3*AbbSum159 + 
     -     Abb1*(-S - T - 2*T14 + 2*T24)))/SW) + 
     -     Eval(ee1,iintP10)*(-(2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(-3*mass2(1) - 3*mass2(2) + 2*mass2(3) - 
     -     2*mass2(4))))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(3*AbbSum96 + 
     -     Abb1*(-S - 2*S34 - T - 2*T14 + 2*T24)))/SW) - 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP17))/SW + 
     -     1/SW2*(Eval(ee1,iintP12)*
     -     (6*Alfa2*MW*(AbbSum132*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb8*EL*gL(2)*
     -     (2*MW*(gL(3)*mass(3)*mass(4)) - 2*MW*(gR(3)*mass(3)*mass(4)))
     -     ) + Alfa2*
     -     (Abb1*(-(2*MW*(dup4*EL*gL(2))) + 
     -     MW*(EL*gL(2)*gL(3)*(6*T24 + 2*U))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(3) - MH2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum69 + Abb1*(-S - T - T14)))/SW) + 
     -     Eval(ee1,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(4)*mass2(4) + gL(4)*(-mass2(3) + MH2)))/
     -     SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum180 + Abb1*(T24 - U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(1)*gR(4)*mass2(4) + 
     -     gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-3*mass2(2) - 4*mass2(3) - MH2))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum33 + Abb1*(T24 + 2*U)))/SW) + 
     -     Eval(ee1,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(-(gL(1)*gR(3)*mass2(3)) + 
     -     gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-4*mass2(2) - 3*mass2(3) + mass2(4) - 
     -     MH2))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum48 + Abb1*(T14 + T24 + 2*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gL(4)*(-3*mass2(1) - mass2(3)) - gR(4)*mass2(4))
     -     ))/SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum1 - 2*Abb1*S))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP14)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(2)*gL(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee11,iintP10)*
     -     (4*Alfa2*MW*(AbbSum99*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP12)*
     -     (4*Alfa2*MW*(AbbSum135*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP13)*
     -     (4*Alfa2*MW*(AbbSum135*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee11,iintP11)*
     -     (4*Alfa2*MW*(AbbSum99*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP15)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee11,iintP16)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee11,iintP17)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup2*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(2*AbbSum136 + 
     -     Abb1*(-S - S34 - T - T14 + T24)))/SW) + 
     -     1/SW2*(Eval(ee12,iintP10)*
     -     (4*Alfa2*MW*(Abb1*dup2*EL*gL(1)*gL(3)) + 
     -     4*Alfa2*MW*(dup5*EL*gL(1)*gL(3))))/SW + 
     -     1/SW2*(Eval(ee12,iintP11)*
     -     (4*Alfa2*MW*(Abb1*dup2*EL*gL(2)*gL(4)) + 
     -     4*Alfa2*MW*(dup5*EL*gL(2)*gL(4))))/SW + 
     -     1/SW2*(Eval(ee12,iintP13)*
     -     (4*Alfa2*MW*(AbbSum183*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(3) - 4*MW*U))))/
     -     SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP15)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(2*mass2(2) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum13 + Abb1*U))/SW) + 
     -     Eval(ee12,iintP14)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(2*mass2(2) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum40 + Abb1*U))/SW) + 
     -     Eval(ee12,iintP16)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(2*mass2(2) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum40 + Abb1*U))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP17)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(2) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum13 + Abb1*U))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee13,iintP14)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(2)*gL(3))))/SW + 
     -     Eval(ee13,iintP10)*
     -     (-(8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(MH2 + 3/2.D0*mass2(3)))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum105 + Abb1*(T24 + 1/2.D0*U)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP12)*
     -     (8*Alfa2*MW*(AbbSum138*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee13,iintP11)*
     -     (8*Alfa2*MW*(AbbSum162*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee13,iintP13)*
     -     (4*Alfa2*MW*(AbbSum184*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(3) - 4*MW*U))))/
     -     SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP11)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(AbbSum163 - Abb1*T24))/SW
     -     ) + Eval(ee13,iintP16)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(2*mass2(2) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum42 - Abb1*U))/SW) + 
     -     Eval(ee13,iintP15)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(2*mass2(2) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum53 - Abb1*U))/SW) + 
     -     Eval(ee14,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum106 + Abb1*(T24 + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee13,iintP17)*
     -     (4*Alfa2*MW*(AbbSum16*EL*gL(2)*gL(4))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP12)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum139 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee14,iintP15)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee14,iintP16)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + 3*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee14,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum74*EL*gL(2)*gL(3))) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee14,iintP13)*
     -     (4*Alfa2*MW*(AbbSum185*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(-(4*MW*MH2) + 4*MW*T24))))/
     -     SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP17)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + 3*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee2,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(-3*mass2(1) - 4*mass2(2) - 9*mass2(3) - 
     -     4*MH2))))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(3*AbbSum109 + 2*Abb1*(S + 2*T24 + 2*U)))/
     -     SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(mass2(2) + 4*mass2(3) + 3*MH2))))/SW + 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum141 + Abb1*(T24 + 2/3.D0*U)))/SW) + 
     -     Eval(ee2,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-(gL(2)*gR(4)*mass2(4)) + 
     -     gL(4)*(gR(2)*mass2(2) + 
     -     gL(2)*(5*mass2(3) + mass2(4) + 3*MH2))))/SW + 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum165 + 
     -     Abb1*(T24 - 1/3.D0*T14 + 2/3.D0*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP13)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-2*mass2(1) - 4*mass2(2) - mass2(3) - MH2)))
     -     /SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(3*AbbSum187 + Abb1*(S + T + T14)))/SW) + 
     -     Eval(ee2,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(2*mass2(2) + 2*mass2(3) + mass2(4) + 
     -     2*MH2))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum76 + Abb1*(S34 + 2*T24 + 2*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) + mass2(2) - mass2(3))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum20 + Abb1*(-S + T)))/SW) + 
     -     Eval(ee2,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-3*mass2(2) - mass2(3))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum43 + Abb1*U))/SW) + 
     -     Eval(ee2,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - 5*mass2(2) - mass2(3) - MH2)))/SW
     -     + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum59 + Abb1*(S34 + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee22,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee22,iintP12)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee22,iintP11)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP13)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee23,iintP10)*
     -     (-(8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum115 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee22,iintP14)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP15)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP16)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP17)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(2)*gL(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(2))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP17)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(2)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP15)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee23,iintP13)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(4)) - 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP16)*
     -     (4*Alfa2*MW*(Abb29*EL*gL(1)*gL(4)) - 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP14)*
     -     (4*Alfa2*MW*(AbbSum79*EL*gL(2)*gL(3))))/SW + 
     -     Eval(ee23,iintP12)*
     -     (-(8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(MH2 + 3/2.D0*mass2(3)))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum145 + Abb1*(T24 + 1/2.D0*U)))/SW) + 
     -     Eval(ee23,iintP11)*
     -     (-(8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(MH2 + 3/2.D0*mass2(3)))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum168 + Abb1*(T24 + 1/2.D0*U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee24,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee24,iintP12)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee24,iintP14)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee24,iintP11)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + 3*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum9 + Abb1*(T24 + U)))/SW) + 
     -     1/SW2*(Eval(ee24,iintP13)*
     -     (4*Alfa2*MW*(AbbSum189*EL*gL(1)*gL(4)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(4*mass2(2) + 8*mass2(3) + 4*MH2)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(3*AbbSum119 + Abb1*(S + 4*T24 + 4*U)))/SW
     -     ) + 1/(CW2*SW2)*
     -     (Eval(ee24,iintP15)*
     -     (-(4*Alfa2*MW*(AbbSum63*EL*gL(1)*gL(3))) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP16)*
     -     (4*Alfa2*MW*(AbbSum45*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(3) + 4*MW*U))))/SW
     -     + 1/(CW2*SW2)*(Eval(ee24,iintP17)*
     -     (4*Alfa2*MW*(AbbSum27*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*gR(2)*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-mass2(1) - 2*mass2(2) - 4*mass2(3) - 
     -     mass2(4) - 4*MH2)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(3*AbbSum171 + 
     -     Abb1*(S34 + T - T14 + 3*T24 + U)))/SW) + 
     -     1/SW2*(Eval(ee3,iintP13)*
     -     (2*Alfa2*MW*(AbbSum178*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(Abb21*EL*gL(4)*gR(1)*mass(1)*mass(2)) - 
     -     6*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee3,iintP12)*
     -     (6*Alfa2*MW*(AbbSum147*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb21*EL*
     -     (2*MW*(gL(2)*gL(3)*mass(1)*mass(2)) - 
     -     2*MW*(gL(3)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb8*EL*gL(2)*
     -     (2*MW*(gL(3)*mass(3)*mass(4)) - 2*MW*(gR(3)*mass(3)*mass(4)))
     -     ) + Alfa2*
     -     (Abb1*(-(2*MW*(dup4*EL*gL(2))) + 
     -     MW*(EL*gL(2)*gL(3)*(6*T24 + 2*U))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee3,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*gR(2)*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(3*mass2(1) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(0*EL*gL(2)*gR(4)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(2)*gL(4)*S)/SW) + 
     -     Eval(ee3,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*
     -     (2*mass2(1) + mass2(2) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee3,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*gR(1)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(3*mass2(2) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     Eval(ee3,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*
     -     (mass2(1) + 6*mass2(2) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee33,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + 2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum122 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee34,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum124 + Abb1*(T24 + U)))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(4)*mass2(2)*
     -     Eval(ee33,iintP13))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*mass2(2)*Eval(ee33,iintP16))/SW + 
     -     1/(CW2*SW2)*(Eval(ee33,iintP15)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee33,iintP12)*
     -     (4*Alfa2*MW*(AbbSum149*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW + 
     -     1/SW2*(Eval(ee33,iintP11)*
     -     (4*Alfa2*MW*(AbbSum173*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 4*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee34,iintP11)*
     -     (-((4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(4)*mass(1))/SW) + 
     -     (2*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(4)*mass(3))/SW + 
     -     (2*Alfa2*MW)/SW2*((-0)*EL*gL(2)*gL(4)*mass(1)*mass(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(4)*mass(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(AbbSum175 + Abb1*T24))/SW
     -     ) + Eval(ee34,iintP12)*
     -     (-((4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(1))/SW) + 
     -     (4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(2)*mass(3))/SW - 
     -     (4*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(4))/SW + 
     -     (8*Alfa2*MW)/SW2*(0*EL*gL(2)*gL(3)*mass(1)*mass(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum151 + Abb1*(-T24 - U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP17)*
     -     (4*Alfa2*MW*(AbbSum30*EL*gL(2)*gL(4)) - 
     -     4*Alfa2*MW*(0*EL*gL(2)*gL(4)*mass(1)) + 
     -     2*Alfa2*MW*((-0)*EL*gL(2)*gL(4)*mass(1)*mass(3))))/SW + 
     -     1/SW2*(Eval(ee34,iintP13)*
     -     (-(4*Alfa2*MW*(AbbSum190*EL*gL(1)*gL(4))) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2)) - 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(4)*mass(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum83*EL*gL(2)*gL(3))) + 
     -     4*Alfa2*MW*(0*EL*gL(2)*gL(3)*mass(1)) + 
     -     4*Alfa2*MW*(0*EL*gL(2)*gL(3)*mass(1)*mass(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP15)*
     -     (4*Alfa2*MW*(AbbSum65*EL*gL(1)*gL(3)) - 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(3)*mass(1)) - 
     -     2*Alfa2*MW*(0*EL*gL(1)*gL(3)*mass(2)) - 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2)) + 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(3)*mass(4)) + 
     -     2*Alfa2*MW*((-0)*EL*gL(1)*gL(3)*mass(2)*mass(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP16)*
     -     (-(4*Alfa2*MW*(AbbSum46*EL*gL(1)*gL(4))) + 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(4)*mass(2)) - 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(4)*mass(3)) + 
     -     4*Alfa2*MW*(0*EL*gL(1)*gL(4)*mass(2)*mass(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP10)*((2*Alfa2*MW)/SW2*(dup6*EL*gL(1)*gL(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-2*mass2(1) - 2*mass2(2) - 2*mass2(3) - 
     -     2*MH2)))/SW) + 
     -     Eval(ee4,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*((-Abb3)*EL*gL(2)*gR(3)*mass(1)*mass(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + gL(2)*(mass2(3) + MH2)))/SW
     -     + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum127 + Abb1*(T24 + U)))/SW) + 
     -     Eval(ee4,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb8*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(Abb34 + Abb1*(S + 2*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(2*mass2(3) + MH2 + 1/2.D0*mass2(1)))/
     -     SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP13)*((4*Alfa2*MW)/SW2*
     -     (AbbSum192*EL*gL(1)*gL(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*(gL(4)*mass2(2) + gR(4)*mass2(4)))/SW) + 
     -     Eval(ee4,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(1)*gR(3)*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) - 2*mass2(2) - mass2(3) + 
     -     mass2(4) - MH2)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S + T14 + T24))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP14)*
     -     (-2*Alfa2*MW*(Abb8*EL*gL(2)*gR(3)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb21*EL*
     -     (2*MW*(gL(2)*gL(3)*mass(1)*mass(2)) - 
     -     2*MW*(gL(3)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*
     -     (-(2*MW*(dup3*EL*gL(3))) + 
     -     MW*(EL*gL(2)*gL(3)*(-2*S - 2*T - 2*T14))))))/SW

      Cloop(1) = Cloop(1) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*mass2(3)*Eval(ee44,iintP16))/SW + 
     -     Eval(ee4,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb8*EL*gL(1)*gR(4)*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW - 
     -     (6*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(3) + 1/3.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee44,iintP11)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(3))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP17)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(3))))/SW + 
     -     1/SW2*(Eval(ee44,iintP13)*
     -     (4*Alfa2*MW*(AbbSum194*EL*gL(1)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP17)*
     -     (2*Alfa2*MW*(Abb8*EL*gL(2)*gR(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW

      
      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagons with photons of the loop 
*** component for the pattern: q + qbar -> q + H + qbar

      subroutine pentP_p_qqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
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
      double complex dup1, dup2

      double precision Alfa, Alfa2

** Parameter declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup2 = 2*mass2(1) + mass2(2) + mass2(3) + MH2


      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum76 + Abb1*(S + 2*T + T14)))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Dval(dd0,iintP1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Dval(dd0,iintP2))/SW - 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Eval(ee00,iintP3))/SW - 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Eval(ee00,iintP4))/SW - 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Eval(ee00,iintP5))/SW - 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Eval(ee00,iintP6))/SW + 
     -     Eval(ee1,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum57 + Abb1*(-T14 - 2*U)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*(mass2(3) + 3/2.D0*mass2(2)))/SW) + 
     -     Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP6)*
     -     (2*MW*(mass2(1)*q(1)*q(3)) + 2*MW*(mass2(3)*q(1)*q(3)) - 
     -     2*MW*(q(1)*q(3)*S)))/SW + 
     -     Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP5)*
     -     (-(2*MW*(mass2(2)*q(2)*q(3))) - 2*MW*(mass2(3)*q(2)*q(3)) + 
     -     2*MW*(q(2)*q(3)*U)))/SW + 
     -     Alfa2/SW2*(Abb1*Eval(ee0,iintP4)*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))/SW + 
     -     1/SW2*(Eval(ee0,iintP3)*(Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S - 2*T - 2*T14))))))/SW + 
     -     1/SW2*(Eval(ee1,iintP5)*
     -     (-(2*Alfa2*MW*(AbbSum34*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*(-(4*MW*(mass2(2)*q(2)*q(3))) - 
     -     4*MW*(mass2(3)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum2 + Abb1*(2*S + T24)))/SW) + 
     -     1/SW2*(Eval(ee11,iintP3)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(1)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(2)*q(1)*q(4)) - 4*MW*(q(1)*q(4)*T)))))/SW + 
     -     1/SW2*(Eval(ee11,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(1)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(3))) - 
     -     4*MW*(mass2(4)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW+
     -     1/SW2*(Eval(ee11,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(3))) - 
     -     4*MW*(mass2(4)*q(2)*q(3)) + 4*MW*(q(2)*q(3)*U)))))/SW+
     -     1/SW2*(Eval(ee11,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(2)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(4))) - 
     -     4*MW*(mass2(4)*q(2)*q(4)) + 4*MW*(q(2)*q(4)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(Abb11 + Abb1*T))/SW) + 
     -     Eval(ee12,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum15 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee12,iintP4)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum15 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee12,iintP5)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(2)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum41 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee12,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum41 + Abb1*(S + 2*T + T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*mass(1)*mass(2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum18 + Abb1*(-S - T - T14)))/SW) + 
     -     Eval(ee14,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(1)*q(4))/SW - (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum79 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee14,iintP4)*(-(4*Alfa2*MW)/SW2*(EL*q(2)*q(4)*
     -     (AbbSum11 - Abb1*U))/SW - (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum44*EL*q(2)*q(3))) - 
     -     4*Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(3)) - 4*MW*(q(2)*q(3)*T)))))/SW + 
     -     1/SW2*(Eval(ee13,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSum60*EL*q(2)*q(4))) - 
     -     4*Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(4)) - 4*MW*(q(2)*q(4)*T)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2) + mass2(3))*q(2)*q(3))/SW-
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum47 + Abb1*(T - U)))/SW) + 
     -     Eval(ee2,iintP3)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 2*mass2(2) + mass2(3) + 
     -     2*mass2(4) + MH2)*q(1)*q(4))/SW - (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum80 + Abb1*(-S34 - T24 - 2*U)))/SW) + 
     -     Eval(ee2,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum64 + 
     -     Abb1*(-S - 2*S34 - T14 - T24 - 2*U)))/SW - (6*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*
     -     (mass2(3) + MH2 + 5/3.D0*mass2(2) + 1/3.D0*mass2(4)))/
     -     SW) + 1/SW2*(Eval(ee14,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(1)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*U)))))/SW
     -     + 1/SW2*(Eval(ee14,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(2)*q(3))) + 4*MW*(q(2)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(AbbSum22 + Abb1*S))/SW) + 
     -     Eval(ee23,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*mass(1)*mass(2)*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum84 + Abb1*(-S - T - T14)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP6)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(1)*q(3)))  + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(2)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*T)))))/SW + 
     -     1/SW2*(Eval(ee22,iintP5)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(2)*q(3))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(2)*q(2)*q(3)) - 4*MW*(q(2)*q(3)*T)))))/SW + 
     -     1/SW2*(Eval(ee22,iintP4)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(2)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(2)*q(2)*q(4)) - 4*MW*(q(2)*q(4)*T)))))/SW + 
     -     1/SW2*(Eval(ee23,iintP4)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(2)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(4)) + 4*MW*(mass2(2)*q(2)*q(4)) - 
     -     4*MW*(q(2)*q(4)*T)))))/SW + 
     -     1/SW2*(Eval(ee22,iintP3)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*q(1)*q(4))) + 
     -     Alfa2*(Abb1*EL*
     -     (-(4*MW*(mass2(2)*q(1)*q(4))) - 
     -     4*MW*(mass2(4)*q(1)*q(4)) + 4*MW*(q(1)*q(4)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP6)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-mass2(1) + mass2(2))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(Abb11 + Abb1*T))/SW) + 
     -     Eval(ee24,iintP4)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup2*EL*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum70 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee3,iintP4)*(-((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*mass(1)*mass(2)*q(2)*q(4))/SW) - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + 5*mass2(2) + mass2(3) + MH2)*q(2)*
     -     q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*(S34 + T + U))/SW) + 
     -     Eval(ee24,iintP3)*(-(4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*
     -     (AbbSum11 - Abb1*U))/SW - (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee24,iintP6)*
     -     (4*Alfa2*MW*(AbbSum29*EL*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee24,iintP5)*
     -     (4*Alfa2*MW*(AbbSum52*EL*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee23,iintP5)*
     -     (-(4*Alfa2*MW*(Abb11*EL*q(2)*q(3))) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(3)) + 4*MW*(mass2(2)*q(2)*q(3)) - 
     -     4*MW*(q(2)*q(3)*T)))))/SW + 1/SW2*(Eval(ee3,iintP3)*
     -     (2*Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(1)*q(4)) + 
     -     Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S - 2*T - 2*T14))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW) + 
     -     Eval(ee4,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) - mass2(4) + 
     -     MH2)*q(1)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*(S + T + T14))/SW) + 
     -     Eval(ee34,iintP3)*((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*mass(1)*mass(2)*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum87 + Abb1*(S + T + T14)))/SW) + 
     -     Eval(ee3,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(2) + mass2(3))*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*mass2(1)*q(1)*q(4)*
     -     Eval(ee33,iintP3))/SW + 1/SW2*(Eval(ee33,iintP6)*
     -     (-4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee34,iintP6)*(4*Alfa2*MW*
     -     (AbbSum30*EL*q(1)*q(3))))/SW + 1/SW2*(Eval(ee34,iintP5)*
     -     (-(4*Alfa2*MW*(Abb24*EL*q(2)*q(3)))))/SW + 
     -     1/SW2*(Eval(ee34,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSum73*EL*q(2)*q(4))) - 
     -     4*Alfa2*MW*(Abb14*EL*mass(1)*mass(2)*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(4)) - 4*MW*(q(2)*q(4)*T)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP6)*(-((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(3)*S)/SW) + (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     Eval(ee4,iintP5)*((2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(3)*(mass2(3) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee44,iintP3)*
     -     (-4*Alfa2*MW*(Abb1*EL*mass2(4)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee44,iintP4)*
     -     (-4*Alfa2*MW*(Abb1*EL*mass2(4)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee4,iintP4)*
     -     (Alfa2*(Abb1*
     -     (-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S + 2*T14 + 2*T24))))))/SW
      

      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagons with no photons of the loop 
*** component for the pattern: q + qbar -> q + H + qbar

      subroutine pentP_nop_qqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
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
      double complex dup3, dup4, dup5, dup6, dup7, dup8, dup9

      double precision Alfa, Alfa2

** Parameter declaration
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"


      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup3 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup4 = mass2(2) + mass2(4) + 1/2.D0*mass2(1)
      dup5 = 2*mass2(1) - mass2(2) + mass2(3) + MH2
      dup6 = mass2(3) + MH2 + 2/3.D0*mass2(1) + 1/3.D0*mass2(2)
      dup7 = mass2(3) + MH2 + 2/3.D0*mass2(1) + 2/3.D0*mass2(2)
      dup8 = 2*mass2(1) + mass2(2) + mass2(3) + MH2
      dup9 = AbbSum101 + Abb1*(-S - 2*S34 - T14 - 2*T24)


      Cloop(1) = Cloop(1) + Eval(ee0,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-mass2(1) - 2*mass2(2) - 3*mass2(3) - MH2)))
     -     /SW - (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum164 + Abb1*(S + 2*T14 + T24)))/SW)+
     -     Eval(ee0,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(gR(3)*mass2(3) + gL(3)*(2*mass2(1) + 2*mass2(2) + 
     -     mass2(3) + 3*MH2 + 2*MZ2))))/SW - (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum90 + 
     -     Abb1*(S + 2*S34 + 2*T + 2*T14 + T24)))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(3)*Dval(dd0,iintP1))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Dval(dd0,iintP7))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(4)*Dval(dd0,iintP8))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Dval(dd0,iintP9))/SW

      Cloop(1) = Cloop(1) + Eval(ee0,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(2) - mass2(3) - MH2)))/
     -     SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee0,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(-(gL(2)*gR(3)*mass2(3)) + 
     -     gL(3)*(gR(2)*mass2(2) + 2*gL(2)*mass2(4))))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum129 - Abb1*U))/SW)+
     -     Eval(ee0,iintP13)*(-(4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(AbbSum191 + Abb1*T))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(2) - 1/2.D0*mass2(1))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee0,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - 4*mass2(2) - mass2(3) - MH2)))/SW
     -     + (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW) - 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP10))/SW - 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP11))/SW - 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP12))/SW - 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP13))/SW - 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP14))/SW + 
     -     Eval(ee0,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*S)/SW + (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*
     -     (MW2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP16)*(Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(2*MW*mass2(2)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) - (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP15))/SW - 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP16))/SW - 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP17))/SW + 
     -     Eval(ee1,iintP10)*((4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum94 + 
     -     Abb1*(-S - 4*S34 - 2*T14 - 3*T24)))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(1/2.D0*(gR(3)*mass2(3)) + 
     -     gL(3)*(mass2(1) + mass2(2) + 3/2.D0*mass2(3) + 
     -     5/2.D0*MH2))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP12)*(-(4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum136 + 
     -     Abb1*(-S - 3*S34 - T14 - 3*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(1/2.D0*(gR(3)*mass2(3))) + 
     -     gL(3)*(-2*mass2(3) - 2*MH2 - 3/2.D0*mass2(1) - 
     -     3/2.D0*mass2(2) + 1/2.D0*mass2(4)))))/SW) + 
     -     Eval(ee1,iintP11)*(-(4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(AbbSum168 + 
     -     Abb1*(-S - 2*S34 - 2*T14 - 3*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(1/2.D0*(gR(4)*mass2(4)) + 
     -     gL(4)*(-2*mass2(1) - mass2(2) - 5/2.D0*mass2(3) - 
     -     3/2.D0*MH2))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP14)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(mass2(1) + mass2(2) + mass2(3) + MH2)))/SW-
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum76 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee1,iintP13)*(-(4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(AbbSum195 + 
     -     Abb1*(-S + S34 - T14 + T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(1)*mass2(1) + 
     -     gL(1)*(gL(4)*(1/2.D0*mass2(1) - 3/2.D0*mass2(2) + 
     -     1/2.D0*mass2(4)) - 1/2.D0*(gR(4)*mass2(4)))))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(1)*gR(3)*mass2(3) + 
     -     gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-4*mass2(2) - 2*mass2(3) - mass2(4)))))/
     -     SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum57 + Abb1*(-T14 - 2*U)))/SW) + 
     -     Eval(ee1,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*
     -     (gR(1)*mass2(1) + gL(1)*(-3*mass2(2) - mass2(3))) - 
     -     gL(1)*gR(4)*mass2(4)))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum34 - 2*Abb1*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(-3*mass2(1) - 4*mass2(3) - MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum2 + Abb1*(2*S + T24)))/SW) + 
     -     Eval(ee11,iintP12)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum140 + Abb1*(S34 + T24)))/SW) + 
     -     Eval(ee11,iintP10)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum99 + Abb1*(S34 + T24)))/SW) + 
     -     Eval(ee11,iintP11)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum99 + Abb1*(S34 + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee11,iintP13)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum140 + Abb1*(S34 + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP14)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(2)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*(4*MW*mass2(2) - 4*MW*T))))/
     -     SW + 1/(CW2*SW2)*(Eval(ee11,iintP15)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW
     -     + 1/(CW2*SW2)*(Eval(ee11,iintP16)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP12)*
     -     ((12*Alfa2*MW)/SW2*(Abb1*dup6*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum141 + 
     -     Abb1*(-S - 2*S34 - T14 - 2*T24)))/SW) + 
     -     1/SW2*(Eval(ee12,iintP10)*
     -     (12*Alfa2*MW*(Abb1*dup6*EL*gL(1)*gL(3)) + 
     -     4*Alfa2*MW*(dup9*EL*gL(1)*gL(3))))/SW + 
     -     1/SW2*(Eval(ee12,iintP11)*
     -     (12*Alfa2*MW*(Abb1*dup6*EL*gL(2)*gL(4)) + 
     -     4*Alfa2*MW*(dup9*EL*gL(2)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP17)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP13)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(1)*gL(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum199 + Abb1*(-S - T14)))/SW) + 
     -     Eval(ee12,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(1)*gL(3))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum15 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee12,iintP14)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(2)*gL(3))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum41 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee12,iintP16)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(1)*gL(4))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum41 + Abb1*(S + 2*T + T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP17)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(2)*gL(4))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum15 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee13,iintP12)*
     -     ((12*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum144 + 
     -     Abb1*(-S - 2*S34 - T14 - 2*T24)))/SW) + Eval(ee13,iintP11)*
     -     ((12*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(AbbSum172 + 
     -     Abb1*(-S - 2*S34 - T14 - 2*T24)))/SW) + Eval(ee13,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum106 + Abb1*(S34 + T24)))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(3) + MH2 + 1/2.D0*mass2(2)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP14)*((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(-mass2(1) + mass2(2)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*(EL*gL(2)*gL(3)*(Abb11 + Abb1*T))/SW
     -     ) + Eval(ee13,iintP17)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup3*EL*gL(2)*gL(4))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum18 + Abb1*(-S - T - T14)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP13)*
     -     (4*Alfa2*MW*(Abb16*EL*gL(1)*gL(4)) + 
     -     4*Alfa2*MW*(Abb14*EL*gL(1)*gL(4)*mass(1)*mass(2)) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee13,iintP15)*
     -     (-(4*Alfa2*MW*(AbbSum60*EL*gL(1)*gL(3))) - 
     -     4*Alfa2*MW*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(4*MW*mass2(1) - 4*MW*T))))/
     -     SW + 1/(CW2*SW2)*(Eval(ee13,iintP16)*
     -     (-(4*Alfa2*MW*(AbbSum44*EL*gL(1)*gL(4))) - 
     -     4*Alfa2*MW*(Abb14*EL*gL(1)*gL(4)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(1) - 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP14)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(2)*gL(3))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum79 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee14,iintP10)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(3) - mass2(4) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum108 + Abb1*(-S34 + T - T24)))/SW)+
     -     Eval(ee14,iintP12)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(3) - mass2(4) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum146 + Abb1*(-S34 + T - T24)))/SW)+
     -     Eval(ee14,iintP13)*((4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum202 + Abb1*(-2*S34 - 2*T24 - U)))/
     -     SW + (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(3) + MH2 + 3/2.D0*mass2(2)))/SW)
     -     + 1/SW2*(Eval(ee14,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum173*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP15)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum11 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee14,iintP16)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*U)))
     -     )/SW + 1/(CW2*SW2)*(Eval(ee14,iintP17)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW -
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(4)*mass2(4) + 
     -     gL(4)*(2*mass2(1) + 2*mass2(2) + 4*mass2(4) + 2*MH2))))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(3*AbbSum176 + 
     -     Abb1*(-S - S34 - T - T24 - 3*U)))/SW) + 
     -     Eval(ee2,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(gR(3)*mass2(3) + 
     -     gL(3)*(-3*mass2(1) - 4*mass2(2) - 4*mass2(3) - 
     -      4*mass2(4) - MH2))))/SW - (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(3*AbbSum113 + Abb1*(-2*S - T24 - 4*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - 3*mass2(2) + mass2(3) + MH2)))/SW
     -     - (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(3*AbbSum206 + Abb1*(-S + 2*T - T14)))/SW)
     -     + Eval(ee2,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(mass2(2) - mass2(3) + 4*mass2(4)))))/SW-
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum149 - 2/3.D0*(Abb1*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(2*mass2(2) + mass2(3) + 3*mass2(4) + MH2))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum80 + Abb1*(-S34 - T24 - 2*U)))/SW)+
     -     Eval(ee2,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(4*mass2(2) + 3*mass2(3) + mass2(4) + 3*MH2))
     -     )/SW - (2*Alfa2*MW)/(CW2*SW2)*(EL*gL(1)*gL(3)*(AbbSum64 + 
     -     Abb1*(-S - 2*S34 - T14 - T24 - 2*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP17)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*(gR(2)*mass2(2) + 
     -     gL(2)*(mass2(1) + mass2(3))))/SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum22 + Abb1*S))/SW) + 
     -     Eval(ee2,iintP16)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) + mass2(3))))/SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum47 + Abb1*(T - U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP10)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/(CW2*SW2)*(Eval(ee22,iintP15)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(4*MW*mass2(2) - 4*MW*T))))/
     -     SW + 1/SW2*(Eval(ee22,iintP13)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(2) - 4*MW*T))))/
     -     SW + 1/(CW2*SW2)*(Eval(ee22,iintP16)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(2) - 4*MW*T))))/
     -     SW + 1/SW2*(Eval(ee22,iintP12)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW
     -     + 1/(CW2*SW2)*(Eval(ee22,iintP14)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW
     -     + 1/SW2*(Eval(ee22,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (-(4*MW*mass2(2)) - 4*MW*mass2(4) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(3)*
     -     (mass2(1) + 3*mass2(2) + mass2(3) + 2*mass2(4) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum117 + Abb1*(-S - T - T14 - 2*U)))/
     -     SW) + Eval(ee23,iintP12)*
     -     (-((8*Alfa2*MW)/SW2*(Abb1*dup4*EL*gL(2)*gL(3))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum153 + Abb1*(T + 2*U)))/SW) + 
     -     Eval(ee23,iintP11)*
     -     (-((8*Alfa2*MW)/SW2*(Abb1*dup4*EL*gL(2)*gL(4))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum179 + Abb1*(T + 2*U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP17)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(4*MW*mass2(2) - 4*MW*T))))/
     -     SW + 1/SW2*(Eval(ee23,iintP13)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(2) - 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP17)*((4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(-mass2(1) + mass2(2)))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*(EL*gL(2)*gL(4)*(Abb11 + Abb1*T))/SW
     -     ) + Eval(ee23,iintP14)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup3*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum84 + Abb1*(-S - T - T14)))/SW) + 
     -     Eval(ee24,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum11 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP15)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(3))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (4*MW*mass2(1) + 4*MW*mass2(2) - 4*MW*T))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP16)*
     -     (-(4*Alfa2*MW*(Abb11*EL*gL(1)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(2) - 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP13)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(1)*gL(4))/SW) + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum208 + Abb1*(S + T14)))/SW) + 
     -     Eval(ee24,iintP15)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup8*EL*gL(1)*gL(3))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum70 + Abb1*(S + 2*T + T14)))/SW) + 
     -     Eval(ee24,iintP12)*
     -     (-(4*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum11 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     Eval(ee24,iintP14)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum11 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(4) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee24,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum11*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/(CW2*SW2)*(Eval(ee24,iintP16)*
     -     (4*Alfa2*MW*(AbbSum52*EL*gL(1)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP17)*
     -     (4*Alfa2*MW*(AbbSum29*EL*gL(2)*gL(4))))/SW + 
     -     Eval(ee3,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-3*gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(1/2.D0*(gL(1)*gR(3)*mass2(3)) + 
     -     gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(2) - mass2(4) - 1/2.D0*mass2(3) + 
     -     1/2.D0*MH2))))/SW - (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum120 + 
     -     Abb1*(S34 - U - 1/2.D0*S + 1/2.D0*T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-2*gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(2*AbbSum183 + 
     -     Abb1*(S34 + 3*T - T14 + 3*U)))/SW - (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(1/2.D0*(gR(4)*mass2(4)) + 
     -     gL(4)*(2*mass2(2) + 3/2.D0*mass2(1) - 
     -      1/2.D0*mass2(3) + 3/2.D0*mass2(4) + 1/2.D0*MH2)))
     -     )/SW) + Eval(ee3,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-3*gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gL(3)*(mass2(1) + mass2(2) + 2*mass2(4)) - 
     -     1/2.D0*(gR(3)*mass2(3)))))/SW + (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum156 + Abb1*(T + 3/2.D0*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(3)*
     -     (mass2(1) + 4*mass2(2) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW) + 
     -     Eval(ee3,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-2*gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(AbbSum191 + Abb1*T))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + 1/2.D0*mass2(1)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee3,iintP14)*(Alfa2*(Abb14*EL*
     -     (2*MW*(gL(2)*gL(3)*mass(1)*mass(2)) + 
     -     2*MW*(gL(3)*gR(2)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (2*MW*mass2(2) + 2*MW*mass2(3) + 2*MW*MH2 + 
     -     MW*(-2*S - 2*T - 2*T14)))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee3,iintP16)*
     -     (-2*Alfa2*MW*(Abb14*EL*gL(4)*gR(1)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(2*MW*mass2(2)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee33,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum122 + Abb1*(S34 + T24)))/SW) + 
     -     Eval(ee33,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + mass2(2) + mass2(4)))/SW+
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum159 + Abb1*(T + U)))/SW) + 
     -     Eval(ee33,iintP11)*
     -     ((4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + mass2(2) + mass2(4)))/SW+
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum186 + Abb1*(T + U)))/SW) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*mass2(1)*Eval(ee33,iintP14))/SW + 
     -     1/(CW2*SW2)*(Eval(ee3,iintP17)*
     -     (2*Alfa2*MW*(Abb14*EL*gL(4)*gR(2)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP10)*
     -     (-(4*Alfa2*MW)/SW2*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(3)*
     -     (mass2(1) + mass2(3) - mass2(4) + MH2))/SW + (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum124 + Abb1*(-S34 - T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee33,iintP17)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(1))))/SW + 
     -     Eval(ee34,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum161 + Abb1*(T + U)))/SW - 
     -     (8*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(3)*
     -     (mass2(4) + 1/2.D0*mass2(1) + 1/2.D0*mass2(2)))/SW) + 
     -     1/SW2*(Eval(ee34,iintP13)*
     -     (-(4*Alfa2*MW*(AbbSum210*EL*gL(1)*gL(4))) + 
     -     4*Alfa2*MW*(Abb14*EL*gL(1)*gL(4)*mass(1)*mass(2)) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee34,iintP11)*
     -     (4*Alfa2*MW*(AbbSum188*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(-(4*MW*mass2(2)) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP14)*
     -     ((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup3*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum87 + Abb1*(S + T + T14)))/SW) + 
     -     Eval(ee4,iintP10)*((2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-2*mass2(1) - 2*mass2(2) - mass2(3) - 
     -     4*mass2(4) - MH2)))/SW - (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum90 + Abb1*(-S - T24 - 2*U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP17)*
     -     (4*Alfa2*MW*(AbbSum30*EL*gL(2)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP16)*
     -     (-(4*Alfa2*MW*(Abb24*EL*gL(1)*gL(4)))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP15)*
     -     (-(4*Alfa2*MW*(AbbSum73*EL*gL(1)*gL(3))) - 
     -     4*Alfa2*MW*(Abb14*EL*gL(1)*gL(3)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(4*MW*mass2(1) - 4*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP11)*(-(2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) - mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum164 + Abb1*(-S - T24)))/SW) + 
     -     Eval(ee4,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 4*gL(2)*mass2(4)))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum129 - Abb1*U))/SW)+
     -     Eval(ee4,iintP13)*(-(4*Alfa2*MW)/SW2*
     -     (Abb14*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(2*AbbSum212 + 
     -     Abb1*(S - S34 + T14 - T24 - U)))/SW - (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(4)*gR(1)*mass2(1) + gL(1)*(gL(4)*
     -     (-2*mass2(2) + 1/2.D0*mass2(1) + 1/2.D0*mass2(4))
     -     - 1/2.D0*(gR(4)*mass2(4)))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee4,iintP14)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*gR(3)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(2) - mass2(3) + 
     -     2*mass2(4) - MH2)))/SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee4,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb14*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*gR(3)*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) - 2*mass2(2) - mass2(3) - 
     -     mass2(4) - MH2)))/SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S + T14 + T24))/SW)

      Cloop(1) = Cloop(1) + Eval(ee44,iintP13)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(2*mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum215 + Abb1*(-S34 - T24 - U)))/SW)+ 
     -     Eval(ee4,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*gR(4)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(2)*gL(4)*S)/SW + 
     -     (6*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(3) + 1/3.D0*mass2(1)))/SW) + 
     -     1/SW2*(Eval(ee44,iintP10)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP15)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(4))))/SW + 
     -     1/SW2*(Eval(ee44,iintP12)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP14)*
     -     (-4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP16)*
     -     (-2*Alfa2*MW*(Abb27*EL*gL(1)*gR(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (-(2*MW*mass2(2)) - 2*MW*mass2(3) + 2*MW*U))))/SW
      

      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagon with photons part of the loop 
*** component for the pattern qbar + qbar -> qbar + H +qbar 

      subroutine pentP_p_qbarqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,gL,
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
      double complex dup1

      double precision Alfa, Alfa2


** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup1 = mass2(1) + mass2(2) + mass2(3) + MH2


      Cloop(1) = Cloop(1) + Eval(ee1,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(4*mass2(1) - mass2(2) + mass2(3) + MH2)*q(1)*
     -     q(4))/SW - (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum71 + Abb1*(S + T14)))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Dval(dd0,iintP1))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Dval(dd0,iintP2))/SW + 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(4)*Eval(ee00,iintP3))/SW + 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*Eval(ee00,iintP4))/SW + 
     -     (12*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*Eval(ee00,iintP5))/SW + 
     -     (16*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*Eval(ee00,iintP6))/SW + 
     -     Eval(ee1,iintP4)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(2) + mass2(3) + 
     -     2*mass2(4) + MH2)*q(2)*q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum57 + Abb1*(-2*S - T14 - T24)))/SW) + 
     -     Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP6)*
     -     (2*MW*(mass2(1)*q(1)*q(3)) + 2*MW*(mass2(3)*q(1)*q(3)) - 
     -     2*MW*(q(1)*q(3)*S)))/SW + 
     -     Alfa2/SW2*(Abb1*Eval(ee0,iintP3)*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S - 2*T - 2*T14))))/SW + 
     -     Alfa2/SW2*(Abb1*EL*Eval(ee0,iintP5)*
     -     (-(2*MW*(mass2(2)*q(2)*q(3))) - 2*MW*(mass2(3)*q(2)*q(3)) + 
     -     2*MW*(q(2)*q(3)*U)))/SW + 1/SW2*(Eval(ee0,iintP4)*
     -     (Alfa2*(Abb1*(-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee1,iintP5)*((2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum27+Abb1*(2*S34 + 2*T + T24 + 2*U)))/SW - 
     -     (6*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*
     -     (mass2(3) + MH2 + 2/3.D0*mass2(1) + 2/3.D0*mass2(2)))/
     -     SW) + 1/SW2*(Eval(ee11,iintP3)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(1)*q(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee1,iintP6)*
     -     (2*Alfa2*MW*(AbbSum2*EL*q(1)*q(3)) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee11,iintP5)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(2)*q(3)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(3)) + 4*MW*(mass2(4)*q(2)*q(3)) - 
     -     4*MW*(q(2)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee11,iintP4)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(2)*q(4)) + 
     -     Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(4)) + 4*MW*(mass2(4)*q(2)*q(4)) - 
     -     4*MW*(q(2)*q(4)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum12 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee12,iintP5)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(2)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum38 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee12,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(1)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum38 + Abb1*(-S34 - T - U)))/SW) + 
     -     1/SW2*(Eval(ee11,iintP6)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(1)*q(3)) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(4)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP6)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(3)*(AbbSum12 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee13,iintP5)*(-((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(2)*q(3))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*mass(1)*mass(2)*q(2)*q(3))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum41 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee13,iintP4)*(-((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(2)*q(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*mass(1)*mass(2)*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum61 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee14,iintP3)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(1)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum76 + Abb1*(S34 + T + U)))/SW) + 
     -     1/SW2*(Eval(ee13,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum15*EL*q(1)*q(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee13,iintP3)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(1)*q(4)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(2)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum64 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee14,iintP4)*((4*Alfa2*MW)/SW2*(EL*q(2)*q(4)*(AbbSum9 - 
     -     Abb1*S))/SW + (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(2)*q(4)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     Eval(ee2,iintP3)*((2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSum77 + Abb1*(-2*S - S34)))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*(mass2(3) + 3/2.D0*mass2(1)))/SW) + 
     -     1/SW2*(Eval(ee14,iintP6)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(1)*q(3)) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/SW2*(Eval(ee14,iintP5)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(2)*q(3)) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(2)*q(3)) - 4*MW*(q(2)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*q(1)*q(3)*(AbbSum19 + Abb1*S))/SW) + 
     -     Eval(ee2,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) - mass2(2) + mass2(3))*q(2)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSum45 + Abb1*(-T - U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP6)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(1)*q(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(1)*q(3))))/SW + 
     -     1/SW2*(Eval(ee22,iintP5)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(2)*q(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee23,iintP3)*
     -     (-(4*Alfa2*MW*(AbbSum81*EL*q(1)*q(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(4)) + 
     -     4*Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee22,iintP4)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(2)*q(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee22,iintP3)*
     -     (4*Alfa2*MW*(AbbSum9*EL*q(1)*q(4)) + Alfa2*(Abb1*EL*
     -     (4*MW*(mass2(1)*q(1)*q(4)) + 4*MW*(mass2(4)*q(1)*q(4)) - 
     -     4*MW*(q(1)*q(4)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP3)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*mass(1)*mass(2)*q(1)*q(4))/SW + (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(5*mass2(1) + mass2(2) + mass2(3) + MH2)*q(1)*
     -     q(4))/SW - (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*(S + T + T14))/SW) + 
     -     Eval(ee24,iintP4)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(2) + mass2(3) + MH2)*q(2)*q(4))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum68 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee24,iintP3)*((4*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(AbbSum9 - 
     -     Abb1*S))/SW + (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(4)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2*(Eval(ee23,iintP6)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(1)*q(3)))))/SW + 
     -     1/SW2*(Eval(ee24,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSum24*EL*q(1)*q(3)))))/SW + 
     -     1/SW2*(Eval(ee23,iintP5)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(2)*q(3))) + 
     -     8*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee24,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSum51*EL*q(2)*q(3)))))/SW + 
     -     1/SW2*(Eval(ee23,iintP4)*
     -     (-(4*Alfa2*MW*(Abb22*EL*q(2)*q(4))) + 
     -     8*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(4))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP6)*((2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/SW) + 
     -     Eval(ee34,iintP4)*(-((4*Alfa2*MW)/SW2*
     -     (Abb1*dup1*EL*q(2)*q(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*mass(1)*mass(2)*q(2)*q(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(4)*(AbbSum69 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee3,iintP5)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(-2*mass2(2) + mass2(3))*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(3)*U)/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*mass2(2)*q(2)*q(4)*
     -     Eval(ee33,iintP4))/SW + 1/SW2*(Eval(ee34,iintP6)*
     -     (-(4*Alfa2*MW*(Abb34*EL*q(1)*q(3)))))/SW + 
     -     1/SW2*(Eval(ee33,iintP5)*
     -     (4*Alfa2*MW*(Abb1*EL*mass2(2)*q(2)*q(3))))/SW + 
     -     1/SW2*(Eval(ee34,iintP5)*
     -     (4*Alfa2*MW*(AbbSum52*EL*q(2)*q(3)) ))/SW + 
     -     1/SW2*(Eval(ee34,iintP3)*
     -     (4*Alfa2*MW*(AbbSum87*EL*q(1)*q(4)) + 
     -     4*Alfa2*MW*(Abb1*EL*mass2(1)*q(1)*q(4)) + 
     -     4*Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee3,iintP4)*
     -     (-2*Alfa2*MW*(Abb21*EL*mass(1)*mass(2)*q(2)*q(4)) + 
     -     Alfa2*(Abb1*(-(2*MW*(dup1*EL*q(2)*q(4))) + 
     -     MW*(EL*q(2)*q(4)*(2*S34 + 2*T + 2*U))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP4)*(-(2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(mass2(1) + mass2(2) + mass2(3) - mass2(4) + 
     -     MH2)*q(2)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*q(2)*q(4)*(S34 + T + U))/SW) + 
     -     Eval(ee4,iintP6)*(-(2*Alfa2*MW)/SW2*(Abb1*EL*q(1)*q(3)*S)/
     -     SW + (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2*(Eval(ee44,iintP3)*
     -     (4*Alfa2*MW*(Abb1*EL*mass2(4)*q(1)*q(4))))/SW + 
     -     1/SW2*(Eval(ee44,iintP4)*
     -     (4*Alfa2*MW*(Abb1*EL*mass2(4)*q(2)*q(4))))/SW + 
     -     1/SW2*(Eval(ee4,iintP3)*(Alfa2*(Abb1*
     -     (2*MW*(dup1*EL*q(1)*q(4)) + 
     -     MW*(EL*q(1)*q(4)*(-2*S34 - 2*T24 - 2*U))))))/SW + 
     -     1/SW2*(Eval(ee4,iintP5)*(Alfa2*(Abb1*EL*
     -     (2*MW*(mass2(2)*q(2)*q(3)) - 2*MW*(q(2)*q(3)*U)))))/SW


      end


*******************************************************************************
*******************************************************************************

*** This subroutine calculates the pentagon with no photons part of the loop 
*** component for the pattern qbar + qbar -> qbar + H +qbar 

      subroutine pentP_nop_qbarqbar_CC(mass,mass2,S,T,T24,T14,S34,U,q,
     -     gL,gR,Cloop)

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
      double complex dup2, dup3, dup4, dup5, dup6, dup7, dup8

      double precision Alfa, Alfa2


** Parameter declarations
#include "mssm.inc"
#include "global.inc"
#include "abbr_CC.h"
#include "integrals_CC.h"

#include "looptools.h"

      Alfa = AlfaQED
      Alfa2 = Alfa**2


      dup2 = mass2(1) + mass2(2) + mass2(3) + MH2
      dup3 = mass2(2) + mass2(3) - mass2(4) + MH2
      dup4 = mass2(1) + mass2(2) + mass2(3) - mass2(4) + MH2
      dup5 = mass2(1) - mass2(2) + mass2(3) + mass2(4) + MH2
      dup6 = 2*mass2(1) - mass2(2) + mass2(3) + mass2(4) + MH2
      dup7 = 3*mass2(1) + mass2(2) + mass2(3) + 2*mass2(4) + MH2
      dup8 = 2*AbbSum101 + Abb1*(-S + T14 + T24)


      Cloop(1) = Cloop(1) + Eval(ee0,iintP12)*((2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*(gL(4)*mass2(1) + gR(4)*mass2(4)))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(Abb9 - Abb1*S))/SW) + 
     -     Eval(ee0,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(-(gL(1)*gR(3)*mass2(3)) + 
     -     gL(3)*(gR(1)*mass2(1) + 2*gL(1)*mass2(4))))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum89 - Abb1*S))/SW)-
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Dval(dd0,iintP10))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(3)*Dval(dd0,iintP7))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(1)*gL(4)*Dval(dd0,iintP8))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Dval(dd0,iintP9))/SW

      Cloop(1) = Cloop(1) + Eval(ee0,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-4*mass2(1) + mass2(2) - mass2(3) - MH2)))/
     -     SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee0,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(-2*mass2(1) - mass2(3) - mass2(4) - MH2))
     -     ))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum125 + Abb1*(-S - S34 - T + T14)))/
     -     SW) + Eval(ee0,iintP14)*(-(4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum187 + Abb1*T))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(-2*mass2(2) + 1/2.D0*mass2(1))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee0,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(2*mass2(1) - mass2(2) + mass2(3)))/
     -     SW + (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     Eval(ee0,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) - 2*mass2(2) - mass2(3) - MH2)))/
     -     SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW) + 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP11))/SW + 
     -     (16*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP12))/SW + 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP13))/SW + 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP14))/SW + 
     -     (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*Eval(ee00,iintP15))/SW + 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*Eval(ee00,iintP16))/SW + 
     -     1/(CW2*SW2)*(Eval(ee0,iintP18)*
     -     (Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + (12*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*Eval(ee00,iintP17))/SW + 
     -     (16*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*Eval(ee00,iintP18))/SW + 
     -     Eval(ee1,iintP11)*((4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(3*AbbSum96 + Abb1*(-S + 2*T14 + 2*T24)))/
     -     SW + (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(gL(3)*
     -     (-mass2(1) - mass2(2) - mass2(3) + 
     -      mass2(4) - MH2) - 1/2.D0*(gR(3)*mass2(3)))))/SW
     -     ) + Eval(ee1,iintP12)*(-(4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(3*AbbSum167 + 
     -     Abb1*(-S + 2*T14 + 2*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(4)*(gR(2)*mass2(2) + gL(2)*(mass2(3) + MH2)) - 
     -     1/2.D0*(gL(2)*gR(4)*mass2(4))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP14)*(-(4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum191 + 
     -     Abb1*(S - S34 + 3*T14 + 2*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(4)*gR(1)*mass2(1) + 
     -     gL(1)*(1/2.D0*(gR(4)*mass2(4)) + 
     -     gL(4)*(3/2.D0*mass2(1) - 3/2.D0*mass2(2) + 
     -      1/2.D0*mass2(3) + 1/2.D0*mass2(4) + 3/2.D0*MH2)))
     -     )/SW) + Eval(ee1,iintP13)*(-(4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(1/2.D0*(gR(3)*mass2(3)) + 
     -     gL(3)*(1/2.D0*mass2(1) - 1/2.D0*mass2(2) + 
     -      1/2.D0*mass2(3) - 1/2.D0*mass2(4) + 1/2.D0*MH2)))
     -     )/SW + (6*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum131 + 
     -     Abb1*(T14 - 1/3.D0*S - 1/3.D0*S34 + 2/3.D0*T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(3*mass2(1) - mass2(2) + mass2(3) + MH2)))/SW
     -     - (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum71 + Abb1*(S + T14)))/SW) + 
     -     Eval(ee1,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(2*mass2(1) + mass2(3) + 3*mass2(4) + MH2))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum57 + Abb1*(-2*S - T14 - T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP18)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gL(4)*(-3*mass2(1) - mass2(3)) - gR(4)*mass2(4))
     -     ))/SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum2 - 2*Abb1*S))/SW) + 
     -     Eval(ee1,iintP17)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*(gL(4)*gR(1)*mass2(1) + 
     -     gL(1)*(gR(4)*mass2(4) + 
     -     gL(4)*(-2*mass2(1) - 3*mass2(2) - 4*mass2(3) - 3*MH2))))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(EL*gL(1)*gL(4)*(AbbSum27 + 
     -     Abb1*(2*S34 + 2*T + T24 + 2*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee11,iintP11)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum100 + Abb1*(T14 + T24)))/SW) + 
     -     Eval(ee11,iintP12)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum100 + Abb1*(T14 + T24)))/SW) + 
     -     Eval(ee11,iintP13)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum135 + Abb1*(T14 + T24)))/SW) + 
     -     Eval(ee11,iintP14)*(-(4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum135 + Abb1*(T14 + T24)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP15)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(2)*gL(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP16)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 1/SW2*(Eval(ee12,iintP11)*
     -     (-(4*Alfa2*MW*(Abb1*dup3*EL*gL(1)*gL(3))) + 
     -     4*Alfa2*MW*(dup8*EL*gL(1)*gL(3))))/SW + 
     -     1/SW2*(Eval(ee12,iintP12)*
     -     (-(4*Alfa2*MW*(Abb1*dup3*EL*gL(2)*gL(4))) + 
     -     4*Alfa2*MW*(dup8*EL*gL(2)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP17)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee11,iintP18)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee12,iintP14)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup6*EL*gL(1)*gL(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum195 + Abb1*(-S - T14 - T24)))/SW)+
     -     Eval(ee12,iintP13)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup3*EL*gL(2)*gL(3))/SW) + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(2*AbbSum137 + Abb1*(-S + T14 + T24)))/SW)
     -     + Eval(ee12,iintP15)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum38 + Abb1*(-S34 - T - U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP11)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup4*EL*gL(1)*gL(3))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(2*AbbSum105 + Abb1*(-S + T14 + T24)))/SW)
     -     + Eval(ee12,iintP16)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum12 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee12,iintP18)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum12 + Abb1*(-S34 - T - U)))/SW) + 
     -     Eval(ee12,iintP17)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum38 + Abb1*(-S34 - T - U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP14)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(1)*gL(4))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(1)*gL(4)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum198 + Abb1*(-S - T14 - T24)))/SW)+
     -     Eval(ee13,iintP16)*
     -     (-((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup2*EL*gL(1)*gL(3))/SW) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum61 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee13,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum142 + Abb1*(T14 + T24)))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(3) + MH2 + 1/2.D0*mass2(1)))/SW)
     -     + Eval(ee13,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum171 + Abb1*(T14 + T24)))/SW - 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(3) + MH2 + 1/2.D0*mass2(1)))/SW)
     -     + 1/(CW2*SW2)*(Eval(ee13,iintP15)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(2)*gL(3)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP11)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup4*EL*gL(1)*gL(3))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum107 + Abb1*(-T14 - T24)))/SW) + 
     -     Eval(ee14,iintP13)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup4*EL*gL(2)*gL(3))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum144 + Abb1*(-T14 - T24)))/SW) + 
     -     Eval(ee13,iintP17)*
     -     (-((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup2*EL*gL(1)*gL(4))/SW) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(1)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum41 + Abb1*(S34 + T + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee13,iintP18)*
     -     (-(4*Alfa2*MW*(AbbSum15*EL*gL(2)*gL(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb21*EL*gL(2)*gL(4)*mass(1)*mass(2))))/SW + 
     -     1/SW2*(Eval(ee14,iintP12)*
     -     (4*Alfa2*MW*(AbbSum173*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(4*MW*mass2(1) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP15)*
     -     (-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum76 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee14,iintP16)*((4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum9 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     Eval(ee14,iintP14)*(-(4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum199 - Abb1*(S + 2*T14 + 2*T24)))/SW
     -     - (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW)
     -     + 1/(CW2*SW2)*(Eval(ee14,iintP17)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(1)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(4)*(4*MW*mass2(1) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(-(gR(4)*mass2(4)) + 
     -     gL(4)*(-2*mass2(1) + mass2(3) - mass2(4) + MH2))
     -     ))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(3*AbbSum175 + Abb1*(-2*S + T14 + T24)))/
     -     SW) + Eval(ee2,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*gR(1)*mass2(1) + 
     -     gL(1)*(-(gR(3)*mass2(3)) + 
     -     gL(3)*(mass2(1) - mass2(3) + 4*mass2(4)))))/SW+
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum110 - 2/3.D0*(Abb1*S)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee14,iintP18)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(4*MW*mass2(1) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(2*mass2(1) + 5*mass2(2) + 4*mass2(3) + 5*MH2))))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(3*AbbSum148 + 
     -     Abb1*(2*S34 + 2*T + 4*T14 + 3*T24 + 2*U)))/SW) + 
     -     Eval(ee2,iintP14)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(3*mass2(1) - 5*mass2(2) + mass2(3) + MH2)))/
     -     SW - (6*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(AbbSum204 + 
     -     Abb1*(-(1/3.D0*S) - 1/3.D0*T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*(gL(3)*gR(2)*mass2(2) + 
     -     gL(2)*(gR(3)*mass2(3) + 
     -     gL(3)*(-4*mass2(1) - 2*mass2(3) - mass2(4)))))/
     -     SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum77 + Abb1*(-2*S - S34)))/SW) + 
     -     Eval(ee2,iintP16)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + gL(1)*(mass2(3) + MH2)))/SW
     -     - (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum64 + Abb1*(-S34 - T - U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP18)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(2)*mass2(2) + 
     -     gL(2)*(mass2(1) + mass2(3))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(4)*(AbbSum19 + Abb1*S))/SW) + 
     -     Eval(ee2,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(4)*(gR(1)*mass2(1) + 
     -     gL(1)*(mass2(1) - 2*mass2(2) + mass2(3))))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(4)*(AbbSum45 + Abb1*(-T - U)))/SW) + 
     -     1/SW2*(Eval(ee22,iintP11)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(1)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 1/(CW2*SW2)*(Eval(ee22,iintP16)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee22,iintP14)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee22,iintP13)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP15)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(3)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(3)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW + 
     -     1/SW2*(Eval(ee22,iintP12)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (4*MW*mass2(1) + 4*MW*mass2(4) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(3))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(2*AbbSum152 + 
     -     Abb1*(-2*S - S34 - T - U)))/SW) + Eval(ee23,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb1*dup7*EL*gL(2)*gL(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(2*AbbSum178 + 
     -     Abb1*(-2*S - S34 - T - U)))/SW) + Eval(ee23,iintP11)*
     -     (-(4*Alfa2*MW)/SW2*(Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW + 
     -     (8*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum114 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP17)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee22,iintP18)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(2)*gL(4))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(4)*mass2(2))))/SW

      Cloop(1) = Cloop(1) + 1/(CW2*SW2)*(Eval(ee23,iintP18)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(2)*gL(4)))))/SW + 
     -     Eval(ee24,iintP11)*
     -     ((4*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum9 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP16)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(3))) + 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee23,iintP14)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(4))) + 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP17)*
     -     (-(4*Alfa2*MW*(Abb22*EL*gL(1)*gL(4))) + 
     -     8*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee23,iintP15)*
     -     (-(4*Alfa2*MW*(AbbSum81*EL*gL(2)*gL(3))) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2)) ))/SW + 
     -     1/SW2*(Eval(ee24,iintP12)*
     -     (4*Alfa2*MW*(AbbSum9*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(4*MW*mass2(1) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP14)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup6*EL*gL(1)*gL(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum207 + Abb1*(-S - T14 - T24)))/SW)+
     -     Eval(ee24,iintP16)*(-(4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(mass2(2) + mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum68 + Abb1*(S34 + T + U)))/SW) + 
     -     Eval(ee24,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum9 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(4) + 1/2.D0*mass2(1)))/SW) + 
     -     Eval(ee24,iintP15)*((4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(2)*gL(3)*(AbbSum9 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(mass2(4) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 1/(CW2*SW2)*(Eval(ee24,iintP18)*
     -     (-(4*Alfa2*MW*(AbbSum24*EL*gL(2)*gL(4)))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee24,iintP17)*
     -     (-(4*Alfa2*MW*(AbbSum51*EL*gL(1)*gL(4)))))/SW + 
     -     Eval(ee3,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-3*gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (6*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(3)*(AbbSum117 - Abb1*S))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*(gL(3)*(gR(1)*mass2(1) + 2*gL(1)*mass2(4)) - 
     -     1/2.D0*(gL(1)*gR(3)*mass2(3))))/SW) + 
     -     Eval(ee3,iintP12)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-2*gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(4) + gR(4))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(2)*gL(4)*(3*AbbSum180 - 
     -     Abb1*(3*S + 2*S34 + 2*T + 2*U)))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(2)*mass2(2) + 
     -     gL(2)*(gL(4)*(-2*mass2(2) - mass2(3) - mass2(4) - MH2 - 
     -     5/2.D0*mass2(1)) - 1/2.D0*(gR(4)*mass2(4)))))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP15)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(gL(2) + gR(2))*mass(1)*mass(2))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(2)*gL(3)*
     -     (4*mass2(1) + mass2(2) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S + T + T14))/SW) + 
     -     Eval(ee3,iintP14)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-2*gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(4)*(AbbSum187 + Abb1*T))/SW + 
     -     (12*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*(mass2(2) + 1/6.D0*mass2(1)))/SW) + 
     -     Eval(ee3,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-3*gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(1/2.D0*(gL(2)*gR(3)*mass2(3)) + 
     -     gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(3/2.D0*mass2(3) + 3/2.D0*MH2))))/SW + 
     -     (6*Alfa2*MW)/SW2*(EL*gL(2)*gL(3)*(AbbSum155 + 
     -     Abb1*(T24 + 4/3.D0*T14 + 1/3.D0*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(4)*gR(1)*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(-3*mass2(2) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     1/(CW2*SW2)*(Eval(ee3,iintP18)*
     -     (2*Alfa2*MW*(Abb21*EL*gL(4)*gR(2)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2*(Eval(ee33,iintP11)*
     -     (4*Alfa2*MW*(AbbSum120*EL*gL(1)*gL(3)) - 
     -     4*Alfa2*MW*(Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(4*MW*mass2(4) - 4*MW*S))))/
     -     SW + 1/(CW2*SW2)*(Eval(ee3,iintP16)*(Alfa2*(Abb21*EL*
     -     (-(2*MW*(gL(1)*gL(3)*mass(1)*mass(2))) - 
     -     2*MW*(gL(3)*gR(1)*mass(1)*mass(2)))) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) - 2*MW*MH2 + 
     -     MW*(2*S34 + 2*T + 2*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee33,iintP13)*
     -     ((4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*gL(2)*gL(3)*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum157 + Abb1*(T14 + T24)))/SW) + 
     -     Eval(ee33,iintP12)*
     -     ((4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(4)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(2)*gL(4)*(mass2(3) + MH2))/SW + (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(4)*(AbbSum182 + Abb1*(T14 + T24)))/SW) + 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*mass2(2)*Eval(ee33,iintP16))/SW + 
     -     1/SW2*(Eval(ee33,iintP14)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee33,iintP17)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(1)*gL(4)*mass2(2))))/SW + 
     -     1/SW2*(Eval(ee34,iintP12)*
     -     (-(4*Alfa2*MW*(AbbSum184*EL*gL(2)*gL(4))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(4*MW*mass2(1) - 4*MW*S))))/
     -     SW + 1/SW2*(Eval(ee34,iintP11)*
     -     (-(4*Alfa2*MW*(AbbSum122*EL*gL(1)*gL(3))) - 
     -     4*Alfa2*MW*(Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2)) + 
     -     Alfa2*(Abb1*EL*gL(1)*gL(3)*(8*MW*mass2(4) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP13)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup3*EL*gL(2)*gL(3))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum159 + Abb1*(T14 + T24)))/SW) + 
     -     Eval(ee34,iintP14)*
     -     (-((4*Alfa2*MW)/SW2*(Abb1*dup5*EL*gL(1)*gL(4))/SW) + 
     -     (4*Alfa2*MW)/SW2*(Abb21*EL*gL(1)*gL(4)*mass(1)*mass(2))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum209 + Abb1*(S + T14 + T24)))/SW) + 
     -     Eval(ee34,iintP16)*
     -     (-((4*Alfa2*MW)/(CW2*SW2)*(Abb1*dup2*EL*gL(1)*gL(3))/SW) - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(1)*gL(3)*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/(CW2*SW2)*
     -     (EL*gL(1)*gL(3)*(AbbSum69 + Abb1*(-S34 - T - U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP18)*
     -     (-(4*Alfa2*MW*(Abb34*EL*gL(2)*gL(4)))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP17)*
     -     (4*Alfa2*MW*(AbbSum52*EL*gL(1)*gL(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee34,iintP15)*
     -     (4*Alfa2*MW*(AbbSum87*EL*gL(2)*gL(3)) + 
     -     4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb21*EL*gL(2)*gL(3)*mass(1)*mass(2))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP11)*((2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(1)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(3)*(gR(1)*mass2(1) + 4*gL(1)*mass2(4)))/SW + 
     -     (2*Alfa2*MW)/SW2*(EL*gL(1)*gL(3)*(AbbSum89 - Abb1*S))/SW)+
     -     Eval(ee4,iintP13)*(-(2*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (Abb27*EL*gL(2)*(-gL(3) + gR(3))*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(4*mass2(1) + 4*mass2(2) + 5*mass2(3) + 
     -     5*MH2)))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(AbbSum125 + 
     -     Abb1*(2*S + 2*S34 + 2*T + 4*T14 + 3*T24 + U)))/SW
     -     ) + 1/SW2*(Eval(ee4,iintP12)*
     -     (2*Alfa2*MW*(Abb9*EL*gL(2)*gL(4)) + 
     -     Alfa2*(Abb27*EL*gL(2)*(-(2*MW*(gL(4)*mass(3)*mass(4))) + 
     -     2*MW*(gR(4)*mass(3)*mass(4)))) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*(2*MW*mass2(1) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP15)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(2) + gR(2))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(2)*gR(3)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(3)*(gR(2)*mass2(2) + 
     -     gL(2)*(-2*mass2(1) - mass2(2) - mass2(3) - 
     -     mass2(4) - MH2)))/SW - (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(2)*gL(3)*(S34 + T24 + U))/SW) + 
     -     Eval(ee4,iintP14)*(-(4*Alfa2*MW)/SW2*
     -     (Abb21*EL*gL(4)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (4*Alfa2*MW)/SW2*(Abb1*EL*(gL(4)*gR(1)*mass2(1) + 
     -     gL(1)*(gL(4)*(-3*mass2(2) + 3/2.D0*mass2(1)) + 
     -     1/2.D0*(gR(4)*mass2(4)))))/SW - (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum210 + 
     -     Abb1*(-S + S34 - T14 + U + 1/2.D0*T)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee4,iintP17)*(-(2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(4)*(-mass2(2) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*gR(4)*mass(3)*mass(4))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(1)*gL(4)*U)/SW) + 
     -     Eval(ee4,iintP16)*((2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb21*EL*gL(3)*(-gL(1) + gR(1))*mass(1)*mass(2))/SW - 
     -     (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb27*EL*gL(1)*gR(3)*mass(3)*mass(4))/SW + 
     -     (2*Alfa2*MW)/(CW2*SW2)*(Abb1*EL*gL(3)*(gR(1)*mass2(1) + 
     -     gL(1)*(-mass2(1) - 2*mass2(2) - mass2(3) + 
     -     2*mass2(4) - MH2)))/SW + (2*Alfa2*MW)/(CW2*SW2)*
     -     (Abb1*EL*gL(1)*gL(3)*(S34 + T + U))/SW) + 
     -     1/SW2*(Eval(ee44,iintP11)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee4,iintP18)*
     -     (2*Alfa2*MW*(Abb27*EL*gL(2)*gR(4)*mass(3)*mass(4)) + 
     -     Alfa2*(Abb1*EL*gL(2)*gL(4)*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee44,iintP14)*((4*Alfa2*MW)/SW2*
     -     (Abb1*EL*gL(1)*gL(4)*
     -     (2*mass2(2) + mass2(3) + mass2(4) + MH2))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSum213 + Abb1*(S34 + T + U)))/SW) + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP16)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(1)*gL(3)*mass2(4))))/SW + 
     -     1/SW2*(Eval(ee44,iintP13)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(4))))/SW + 
     -     1/(CW2*SW2)*(Eval(ee44,iintP15)*
     -     (4*Alfa2*MW*(Abb1*EL*gL(2)*gL(3)*mass2(4))))/SW

      end


*******************************************************************************
*******************************************************************************

