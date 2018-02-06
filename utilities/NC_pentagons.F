********************************************************************************
********************************************************************************
***                                                                          ***
*** NC_pentagons.F                                                           ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contains all of the matrix elements for the electroweak        ***
*** pentagon diagrams for the neutral currents.                              ***
***                                                                          ***
********************************************************************************
********************************************************************************

** This subroutine is diagram set P1 for quark + quark -> quark + H + quark 

      subroutine pentP1_qq(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, dup2, dup3, dup4, dup5
      double complex dup6, dup7, dup8, dup9
      double complex tmp1, tmp2, tmp3, tmp4, tmp5, tmp6


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup1 = gL(3)*(gL(3) - 2*gR(3)) + gR(3)**2
      dup2 = AbbSumP239 + Abb1*(2*S + T14 + 2*T24)
      dup3 = AbbSumP247 + Abb2*(-2*S - T14 - 2*T24)
      dup4 = AbbSumP253 + Abb3*(-2*S - T14 - 2*T24)      
      dup5 = AbbSumP255 + Abb4*(2*S + T14 + 2*T24)
        
      tmp1 = 1/CW2*(-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2* 
     -     gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))))/SW

      tmp1 = tmp1 + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(S34 + T + U))/SW

      tmp2 = 1/CW2*(8*Alfa2*MW*(AbbSumP257*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP242*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP252*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP240*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP254*EL*gL(1)**2*gL(3)**2*mass(3)) + 
     -     12*Alfa2*MW*(AbbSumP256*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP246*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP249*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     12*Alfa2*MW*(AbbSumP241*EL*gR(1)**2*gR(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp2 = tmp2 + 1/CW2*
     -     (8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*MH2) + 
     -     Alfa2*(-(4*MW*(dup5*EL*gL(1)**2*gL(3)**2)) + 
     -     4*MW*(dup3*EL*gL(3)**2*gR(1)**2) + 
     -     4*MW*(dup4*EL*gL(1)**2*gR(3)**2) - 
     -     4*MW*(dup2*EL*gR(1)**2*gR(3)**2)))/SW

      tmp3 = -(20*Alfa2*MW)/CW2*
     -     (AbbSumP155*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP140*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP150*EL*gL(1)**2*gR(3)**2*mass2(1))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (20*Alfa2*MW)/CW2*
     -     (AbbSumP137*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP147*EL*gL(1)**2*
     -     (-gL(3)**2 + gR(3)**2)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP138*EL*gR(1)**2*
     -     (-gL(3)**2 + gR(3)**2)*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*MH2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*MH2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*MH2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW
      
      tmp3 = tmp3 - (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP135 + Abb1*(-2*S34 - 2*T - 2*U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP153 + Abb4*(-2*S34 - 2*T - 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP144 + Abb2*(2*S34 + 2*T + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP151 + Abb3*(2*S34 + 2*T + 2*U)))/SW

      tmp4 = 1/CW2*(12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(AbbSumP36*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP20*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     16*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(3)) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp4 = tmp4 + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP7 + Abb1*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP35 + Abb4*(2*S + T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP21 + Abb2*(-S - 1/2.D0*T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP30 + Abb3*(-S - 1/2.D0*T24)))/SW

      Cloop(1) = Cloop(1) + tmp1*Eval(ee0,iintP10) + 
     -     1/CW2*((Dval(dd0,iintP3) + Dval(dd0,iintP4))*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Dval(dd0,iintP1) + Dval(dd0,iintP2))*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee0,iintP11)*(Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP12)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP5)*
     -     (4*Alfa2*MW*(AbbSumP401*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP396*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP399*EL*gL(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP394*EL*gR(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP7)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/CW2*(Eval(ee0,iintP6)*
     -     (12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee0,iintP8)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + tmp2*Eval(ee1,iintP10) + 
     -     1/CW2*((Eval(ee00,iintP12) + Eval(ee00,iintP5))*
     -     (24*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP10) + Eval(ee00,iintP8))*
     -     (32*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) + 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP11) + Eval(ee00,iintP6))*
     -     (24*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*((Eval(ee00,iintP7) + Eval(ee00,iintP9))*
     -     (32*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee0,iintP9)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + tmp3*Eval(ee1,iintP12) + 
     -     Eval(ee1,iintP11)*(-(20*Alfa2*MW)/CW2*
     -     (AbbSumP155*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP140*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP150*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (20*Alfa2*MW)/CW2*
     -     (AbbSumP137*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP147*EL*gL(1)*(-gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP138*EL*gR(1)*(-gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP135 + Abb1*(-2*S34 - 2*T - 2*U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP153 + Abb4*(-2*S34 - 2*T - 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP144 + Abb2*(2*S34 + 2*T + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP151 + Abb3*(2*S34 + 2*T + 2*U)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*mass(1)*q(1)*q(3)*(0 - 3/2.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*mass(3)*q(1)*q(3)*(0 + 1/2.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*mass(3)*q(1)*q(3)*(0 + 1/2.D0*0))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*mass(1)*q(1)*q(3)*(0 - 3/2.D0*0))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*dup1*EL*gL(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*dup1*EL*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP409 + Abb2*(S - S34 + T + T14 - T24)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP413 + Abb3*(S - S34 + T + T14 - T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP404 + Abb1*(-S + S34 - T - T14 + T24)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP416 + Abb4*(-S + S34 - T - T14 + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP6)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP333*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP326*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP330*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP325*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW +
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP324 + Abb1*(-S + T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP328 + Abb2*(-S + T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP331 + Abb3*(-S + T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP332 + Abb4*(-S + T - T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + tmp4*Eval(ee1,iintP8) + 
     -     Eval(ee1,iintP7)*((12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW + (16*Alfa2*MW)/CW2*
     -     (AbbSumP36*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP20*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP26*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP7 + Abb1*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP35 + 
     -     Abb4*(2*S + T24)))/SW + (8*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)*
     -     q(1)*q(3)*(AbbSumP21 + Abb2*(-S - 1/2.D0*T24)))/SW + (8*
     -     Alfa2*MW)/CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP30 + Abb3*
     -     (-S - 1/2.D0*T24)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee1,iintP9)*
     -     (8*Alfa2*MW*(AbbSumP257*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP242*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP252*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP240*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP254*EL*gL(1)*gL(3)*mass(3)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP256*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP246*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP249*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP241*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3)) + 
     -     Alfa2*(-(4*MW*(dup5*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     4*MW*(dup3*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*MW*(dup4*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     4*MW*(dup2*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP32*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee11,iintP12)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP55 + Abb4*(S + T24)))/SW) + 
     -     1/CW2*(Eval(ee11,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP423*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb34*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb42*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP419*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (8*MW*mass2(3) + 8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (8*MW*mass2(3) + 8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (8*MW*mass2(3) + 8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (8*MW*mass2(3) + 8*MW*MH2 - 8*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee11,iintP10) + Eval(ee11,iintP8))*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP55 + Abb4*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee11,iintP7) + Eval(ee11,iintP9))*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP12)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP167 + Abb1*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP169 + Abb2*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP171 + Abb3*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP172 + Abb4*(-S + T)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP426 + Abb1*(S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP428 + Abb2*(S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP430 + Abb3*(S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP431 + Abb4*(S + T)))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP11) + Eval(ee12,iintP6))*
     -     ((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP167 + Abb1*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP169 + Abb2*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP171 + Abb3*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP172 + Abb4*(-S + T)))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP10) + Eval(ee12,iintP8))*
     -     ((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP59 + Abb1*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP61 + Abb2*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP64 + Abb3*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP65 + Abb4*(-S + T)))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP7) + Eval(ee12,iintP9))*
     -     ((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP59 + Abb1*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP61 + Abb2*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP64 + Abb3*(-S + T)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP65 + Abb4*(-S + T)))/SW) + 
     -     1/CW2*(Eval(ee13,iintP10)*
     -     (8*Alfa2*MW*(AbbSumP265*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP262*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(AbbSumP264*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP260*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee13,iintP12)*
     -     (8*Alfa2*MW*(AbbSumP178*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP175*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T))))/SW + 
     -     1/CW2*(Eval(ee13,iintP11)*
     -     (8*Alfa2*MW*(AbbSumP178*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP175*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/
     -     SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee13,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP435*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(AbbSumP434*EL*gR(1)**2*gR(3)**2)))
     -     /SW + 1/CW2*(Eval(ee13,iintP6)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T)))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP7)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP68 + Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP71 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP74 - Abb3*S))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP75 + Abb4*S))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP8)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP68 + Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP71 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP74 - Abb3*S))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP75 + Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee13,iintP9)*
     -     (8*Alfa2*MW*(AbbSumP265*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP262*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP264*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP260*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/
     -     SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP10)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP55 + Abb4*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP11)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP12)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP55 + Abb4*(S + T24)))/SW) + 
     -     1/CW2*(Eval(ee14,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP441*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb34*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb42*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP437*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (8*MW*MH2 - 8*MW*T24)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*(8*MW*MH2 - 8*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee14,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP340*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP339*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP7)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP8)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP55 + Abb4*(S + T24)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee14,iintP9)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP10)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP285*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP278*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP282*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP277*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP276 + Abb1*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP280 + Abb2*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP283 + Abb3*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP284 + Abb4*(S34 + 2*T + U)))/
     -     SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP11)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP189 + Abb1*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP169 + Abb2*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP171 + Abb3*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP193 + Abb4*(2*T + U)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP12)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP8*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP189 + Abb1*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP169 + Abb2*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP171 + Abb3*(2*T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP193 + Abb4*(2*T + U)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP5)*(-(24*Alfa2*MW)/CW2*
     -     (AbbSumP453*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP446*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP450*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP445*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP444 + Abb1*(T + 1/3.D0*S + 1/3.D0*T14)))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP448 + Abb2*(T + 1/3.D0*S + 1/3.D0*T14)))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP451 + Abb3*(T + 1/3.D0*S + 1/3.D0*T14)))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP452 + Abb4*(T + 1/3.D0*S + 1/3.D0*T14)))/SW)


	dup6 = gL(1)*(gL(1) - 2*gR(1)) + gR(1)**2
	dup7 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
	dup8 = 2*mass2(3) + MH2 + 5/2.D0*mass2(1)
	dup9 = 4*mass2(1) + 3*mass2(3) + 2*MH2

        tmp5 = (16*Alfa2*MW)/CW2*
     -     (AbbSumP285*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP278*EL*gL(3)*gR(1)*mass2(1)*
     -     q(1)*q(3))/SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP282*EL*gL(1)*gR(3)*mass2(1)*
     -     q(1)*q(3))/SW + (16*Alfa2*MW)/CW2*
     -     (AbbSumP277*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW + (8*Alfa2*MW)/CW2*(AbbSumP354*EL*
     -     gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP346*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP349*EL*gL(1)*gR(3)*mass2(3)*
     -     q(1)*q(3))/SW + (8*Alfa2*MW)/CW2*
     -     (AbbSumP343*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*
     -     mass2(3)*q(1)*q(3))/SW + (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP342 + Abb1*(2*S + S34 + T24)))/SW +
     -     (8*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP347 + Abb2*(-S - 1/2.D0*S34 - 1/2.D0*T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP351 + Abb3*(-S - 1/2.D0*S34 - 1/2.D0*T24)))/SW

        tmp5 = tmp5 - (4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP353 + Abb4*(2*S + S34 + T24)))/SW

        tmp6 = 1/CW2*(-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2* 
     -     mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))))/SW

        tmp6 = tmp6 + 1/CW2*
     -     (Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*(4*S34 + 4*T + 4*U)) + 
     -     Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*(4*S34 + 4*T + 4*U)))/SW

       Cloop(1) = Cloop(1) + tmp5*Eval(ee2,iintP6) + 
     -     Eval(ee2,iintP7)*((4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP84 + Abb1*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP88 + Abb2*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP92 + Abb3*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP94 + Abb4*(-S + T)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP8)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP23*dup6*EL*gL(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*dup6*EL*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP84 + Abb1*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP88 + Abb2*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP92 + Abb3*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP94 + Abb4*(-S + T)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee2,iintP9)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP285*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP278*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP282*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP277*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP276 + Abb1*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP280 + Abb2*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP283 + Abb3*(S34 + 2*T + U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP284 + Abb4*(S34 + 2*T + U)))/
     -     SW) + 1/CW2*((Eval(ee11,iintP6) + Eval(ee22,iintP11))*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/SW

       Cloop(1) = Cloop(1) + (Eval(ee11,iintP11) + Eval(ee22,iintP6))*
     -     ((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW) + 
     -     1/CW2*((Eval(ee22,iintP12) + Eval(ee22,iintP5))*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T))))/SW

       Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP10) + Eval(ee22,iintP8))*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T))))/SW + 
     -     1/CW2*((Eval(ee22,iintP7) + Eval(ee22,iintP9))*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee23,iintP10)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T))))/SW + 
     -     1/CW2*(Eval(ee23,iintP11)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/
     -     SW

       Cloop(1) = Cloop(1) + Eval(ee23,iintP6)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP357 + Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(Abb32 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(Abb40 - Abb3*S))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP362 + Abb4*S))/SW) + 
     -     1/CW2*((Eval(ee23,iintP12) + Eval(ee23,iintP5))*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee23,iintP8)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW + 
     -     1/CW2*(Eval(ee23,iintP7)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T)))/SW + 
     -     1/CW2*(Eval(ee23,iintP9)*
     -     (8*Alfa2*MW*(Abb64*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/
     -     SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP10)*
     -     (8*Alfa2*MW*(AbbSumP294*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP292*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(AbbSumP293*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP291*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP11)*
     -     (-(8*Alfa2*MW*(AbbSumP205*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*Alfa2*MW*(Abb32*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb40*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP204*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S)))))/SW

       Cloop(1) = Cloop(1) + Eval(ee24,iintP12)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP204 + Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(Abb32 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(Abb40 - Abb3*S))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP205 + Abb4*S))/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP456*EL*gL(1)**2*gL(3)**2)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb63*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP455*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) + 8*MW*T))))/SW

       Cloop(1) = Cloop(1) + Eval(ee24,iintP6)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*
     -     q(3))/SW + (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*(mass2(1) +
     -      mass2(3) + MH2)*q(1)*q(3))/SW + (8*Alfa2*MW)/CW2*(Abb3*EL*
     -     gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW + (8*
     -     Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + 
     -     MH2)*q(1)*q(3))/SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(-S - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(-S - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP46 + Abb1*(S + T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP55 + Abb4*(S + T24)))/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP7)*
     -     (-(8*Alfa2*MW*(AbbSumP106*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     4*Alfa2*MW*(AbbSumP104*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP105*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP103*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S)))))/SW

       Cloop(1) = Cloop(1) + Eval(ee24,iintP8)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP103 + Abb1*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP104 - 2*Abb2*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP105 - 2*Abb3*S))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP106 + Abb4*S))/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP9)*
     -     (8*Alfa2*MW*(AbbSumP294*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP292*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP293*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP291*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP10)*
     -     (-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

       Cloop(1) = Cloop(1) + Eval(ee3,iintP11)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW)

       Cloop(1) = Cloop(1) + Eval(ee3,iintP12)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(-mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*U)/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP5)*
     -     (4*Alfa2*MW*(AbbSumP401*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP396*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP399*EL*gL(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP394*EL*gR(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW

       Cloop(1) = Cloop(1) + Eval(ee3,iintP6)*((28*Alfa2*MW)/CW2*
     -     (AbbSumP381*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP375*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP379*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (28*Alfa2*MW)/CW2*
     -     (AbbSumP373*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(S + T + T14))/SW)

       Cloop(1) = Cloop(1) + Eval(ee3,iintP7)*((4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW)

       Cloop(1) = Cloop(1) + Eval(ee3,iintP8)*(-(4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*S)/SW)

       Cloop(1) = Cloop(1) + Eval(ee3,iintP9)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(S34 + T + U))/SW) +
     -     1/CW2*(Eval(ee33,iintP6)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee33,iintP8)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee33,iintP7)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee34,iintP10)*
     -     (-(8*Alfa2*MW*(AbbSumP304*EL*gL(1)**2*gL(3)**2)) + 
     -     8*Alfa2*MW*(AbbSumP301*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(AbbSumP303*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP299*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(16*MW*mass2(1)) + 8*MW*T))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP12)*
     -     (8*Alfa2*MW*(AbbSumP211*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP210*EL*gR(1)**2*gR(3)**2)))/SW
     -     + 1/CW2*(Eval(ee34,iintP5)*
     -     (8*Alfa2*MW*(AbbSumP461*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP460*EL*gR(1)**2*gR(3)**2)))/
     -     SW + 1/CW2*(Eval(ee34,iintP11)*
     -     (8*Alfa2*MW*(AbbSumP211*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP210*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP385*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP384*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP7)*
     -     (-(8*Alfa2*MW*(AbbSumP112*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     4*Alfa2*MW*(Abb32*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(Abb40*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP111*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S)))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP8)*
     -     (-(8*Alfa2*MW*(AbbSumP112*EL*gL(1)**2*gL(3)**2)) + 
     -     4*Alfa2*MW*(Abb32*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(Abb40*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP111*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (8*MW*mass2(1) + 8*MW*mass2(3) - 8*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (8*MW*mass2(1) + 8*MW*mass2(3) - 8*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) + 8*MW*mass2(3) - 8*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) + 8*MW*mass2(3) - 8*MW*S))))/SW

       Cloop(1) = Cloop(1) + tmp6*Eval(ee4,iintP10) + 
     -     1/CW2*(Eval(ee34,iintP9)*
     -     (-(8*Alfa2*MW*(AbbSumP304*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*Alfa2*MW*(AbbSumP301*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP303*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP299*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(16*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*T)))))/
     -     SW + 1/CW2*(Eval(ee4,iintP11)*(Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U)))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP12)*
     -     (-4*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*U))))/SW

       Cloop(1) = Cloop(1) + Eval(ee4,iintP5)*(-((8*Alfa2*MW)/CW2*
     -     (Abb4*dup8*EL*gL(1)**2*gL(3)**2)/SW) - 
     -     (8*Alfa2*MW)/CW2*(Abb2*dup8*EL*gL(3)**2*gR(1)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*dup8*EL*gL(1)**2*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*dup8*EL*gR(1)**2*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP463 + Abb1*(-S - T - T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP467 + Abb4*(-S - T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP396 + Abb2*(2*S + 2*T + 2*T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP399 + Abb3*(2*S + 2*T + 2*T14)))/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP6)*
     -     (12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U))))))/SW

       Cloop(1) = Cloop(1) + Eval(ee4,iintP7)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*mass(3)*q(1)*q(3)*(0 + 1/2.D0*0))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*mass(3)*q(1)*q(3)*(0 + 1/2.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

       Cloop(1) = Cloop(1) + Eval(ee4,iintP8)*(-(4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*S)/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW)

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP9)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP8*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

       Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee44,iintP12)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP11)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP471*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(AbbSumP469*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*(16*MW*(dup7*EL*gL(1)**2*gL(3)**2) - 
     -     8*MW*(dup9*EL*gL(1)**2*gL(3)**2))) + Alfa2*(Abb2*
     -     (16*MW*(dup7*EL*gL(3)**2*gR(1)**2) - 
     -     8*MW*(dup9*EL*gL(3)**2*gR(1)**2))) + Alfa2*(Abb3*
     -     (16*MW*(dup7*EL*gL(1)**2*gR(3)**2) - 
     -     8*MW*(dup9*EL*gL(1)**2*gR(3)**2))) + Alfa2*(Abb1*
     -     (16*MW*(dup7*EL*gR(1)**2*gR(3)**2) - 
     -     8*MW*(dup9*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee44,iintP8)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP7)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P1 for:
**     antiquark + quark -> antiquark + H + quark 

      subroutine pentP1_qbarq(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, dup2, dup3, dup4, dup5
      double complex dup6, dup7, dup8, dup9, dup10
      double complex dup11, dup12
      double complex tmp1, tmp2, tmp3, tmp4


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup1 = gL(3)*(gL(3) - 2*gR(3)) + gR(3)**2
      dup2 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup3 = MH2 + 5/2.D0*mass2(1) + 3/2.D0*mass2(3)
      dup4 = 4*mass2(1) + 3*mass2(3) + 2*MH2
      dup5 = MZ2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)
      dup6 = AbbSumP243 + Abb2*(U + 1/2.D0*T14 + 1/2.D0*T24)
      dup7 = AbbSumP249 + Abb3*(U + 1/2.D0*T14 + 1/2.D0*T24)

      tmp1 = 1/CW2*(-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*
     -     mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP157*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP146*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP150*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     16*Alfa2*MW*(AbbSumP139*EL*gR(1)**2*gR(3)**2*mass2(3)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp1 = tmp1 - (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP137 + Abb1*(T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP155 + Abb4*(T24 + 2*U)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP147 + Abb2*(U + 1/2.D0*T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP153 + Abb3*(U + 1/2.D0*T24)))/SW

      tmp2 = (8*Alfa2*MW)/CW2*(dup6*EL*gL(3)*gR(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP253*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP238*EL*gL(3)*gR(1)*mass2(1)*q(1)*
     -     q(3))/SW + (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP248*EL*gL(1)*gR(3)*mass2(1)*q(1)*
     -     q(3))/SW - (16*Alfa2*MW)/CW2*
     -     (AbbSumP236*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP15*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW - (8*Alfa2*MW)/CW2*(AbbSumP252*EL*
     -     gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + (4*Alfa2*MW)/CW2*
     -     (AbbSumP242*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP245*EL*gL(1)*gR(3)*mass2(3)*q(1)*
     -     q(3))/SW - (8*Alfa2*MW)/CW2*
     -     (AbbSumP237*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*
     -     mass2(3)*q(1)*q(3))/SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP234 + Abb1*(T14 + T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP250 + Abb4*(T14 + T24 + 2*U)))/SW 

      tmp2 = tmp2 + (8*Alfa2*MW)/CW2*(dup7*EL*gL(1)*gR(3)*q(1)*q(3))/SW

      tmp3 = 1/CW2*(-16*Alfa2*MW*(AbbSumP253*EL*gL(1)**2*gL(3)**2*
     -     mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP238*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP248*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP236*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP252*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP242*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP245*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP237*EL*gR(1)**2*gR(3)**2*mass2(3)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp3 = tmp3 + (8*Alfa2*MW)/CW2*
     -     (dup6*EL*gL(3)**2*gR(1)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(dup7*EL*gL(1)**2*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP234 + Abb1*(T14 + T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP250 + Abb4*(T14 + T24 + 2*U)))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Dval(dd0,iintP3)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*(Dval(dd0,iintP2)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd0,iintP1)*
     -     (Alfa2*(Abb4*EL*gL(1)*
     -     (8*MW*(gL(1)*gL(3)**2) + 8*MW*(gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*
     -     (8*MW*(gL(3)*gR(1)**2) + 8*MW*(gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*
     -     (8*MW*(gL(1)*gR(3)**2) + 8*MW*(gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*
     -     (8*MW*(gR(1)*gR(3)**2) + 8*MW*(gR(3)*q(1)*q(3))))))/SW + 
     -     1/CW2*(Eval(ee0,iintP10)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP11)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP4)*
     -     (4*Alfa2*MW*(AbbSumP397*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP390*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP396*EL*gL(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP387*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(8*MW*(dup3*EL*gL(1)**2*gL(3)**2)) + 
     -     4*MW*(dup4*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(8*MW*(dup3*EL*gL(3)**2*gR(1)**2)) + 
     -     4*MW*(dup4*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(8*MW*(dup3*EL*gL(1)**2*gR(3)**2)) + 
     -     4*MW*(dup4*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(8*MW*(dup3*EL*gR(1)**2*gR(3)**2)) + 
     -     4*MW*(dup4*EL*gR(1)**2*gR(3)**2))) + 
     -     Alfa2*(0*EL*gL(3)*
     -     (8*MW*(gL(3)*gR(1)**2*mass(1)*mass(3)) - 
     -     8*MW*(gR(1)**2*gR(3)*mass(1)*mass(3)))) + 
     -     Alfa2*(0*EL*gL(1)**2*
     -     (-(8*MW*(gL(3)*gR(3)*mass(1)*mass(3))) + 
     -     8*MW*(gR(3)**2*mass(1)*mass(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP6)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/CW2*(Eval(ee0,iintP5)*
     -     (12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP7)*(Alfa2*(Abb4*
     -     (8*MW*(dup5*EL*gL(1)**2*gL(3)**2) - 
     -     4*MW*(EL*gL(1)**2*gL(3)**2*S))) + Alfa2*(Abb2*
     -     (8*MW*(dup5*EL*gL(3)**2*gR(1)**2) - 
     -     4*MW*(EL*gL(3)**2*gR(1)**2*S))) + Alfa2*(Abb3*
     -     (8*MW*(dup5*EL*gL(1)**2*gR(3)**2) - 
     -     4*MW*(EL*gL(1)**2*gR(3)**2*S))) + Alfa2*(Abb1*
     -     (8*MW*(dup5*EL*gR(1)**2*gR(3)**2) - 
     -     4*MW*(EL*gR(1)**2*gR(3)**2*S)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP8)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP9)*
     -     (-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee00,iintP7) + Eval(ee00,iintP9))*
     -     (-(24*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP11) + Eval(ee00,iintP4))*
     -     (-(32*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) - 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP6) + Eval(ee00,iintP8))*
     -     (-(24*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*((Eval(ee00,iintP10) + Eval(ee00,iintP5))*
     -     (-(32*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee1,iintP10)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP15*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW - (16*Alfa2*MW)/CW2*
     -     (AbbSumP157*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP146*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP150*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP139*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP137 +
     -      Abb1*(T24 + 2*U)))/SW+(4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*
     -     q(3)*(AbbSumP155 + Abb4*(T24 + 2*U)))/SW+ (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP147 + Abb2*(U + 1/2.D0*
     -     T24)))/SW + (8*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP153 + Abb3*(U + 1/2.D0*T24)))/SW)

      Cloop(1) = Cloop(1) + tmp1*Eval(ee1,iintP11) + 
     -     Eval(ee1,iintP4)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP32*dup1*EL*gL(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*dup1*EL*gR(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP399 + Abb1*(T24 - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP405 + Abb2*(T24 - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP410 + Abb3*(T24 - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP412 + Abb4*(T24 - U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP5)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP325*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP318*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP322*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP317*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP316 + Abb1*(-S - T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP320 + Abb2*(-S - T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP323 + Abb3*(-S - T - T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP324 + Abb4*(-S - T - T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee1,iintP6)*
     -     (4*Alfa2*MW*(AbbSumP36*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP20*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP30*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP8*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP26*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) - 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP32*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP10*EL*
     -     (-(4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3))) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee1,iintP7)*((12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP26*EL*gL(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP10*EL*gR(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP8 - 2*Abb1*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP20 - 2*Abb2*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP30 - 2*Abb3*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP36 - 2*Abb4*S))/SW)

      Cloop(1) = Cloop(1) + tmp2*Eval(ee1,iintP8) + tmp3*
     -     Eval(ee1,iintP9) + Eval(ee11,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP54 + Abb4*(T24 + U)))/SW) + 
     -     1/CW2*(Eval(ee11,iintP4)*
     -     (8*Alfa2*MW*(AbbSumP420*EL*gL(1)**2*gL(3)**2) - 
     -     4*Alfa2*MW*(AbbSumP417*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP419*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP415*EL*gR(1)**2*gR(3)**2)+ 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(3)) - 8*MW*MH2 + 8*MW*T24))))/SW

      Cloop(1) = Cloop(1) + (Eval(ee11,iintP6) + Eval(ee11,iintP8))*
     -     (-(8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*(mass2(1) + 2*
     -     mass2(3) + MH2)*q(1)*q(3))/SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*
     -     q(3))/SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*
     -     q(3))/SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*
     -     q(3))/SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee11,iintP7) + Eval(ee11,iintP9))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP169 + Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP171 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP173 + Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP174 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP4)*((8*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP422 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP424 - Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP426 - Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP427 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP10) + Eval(ee12,iintP5))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP169 + Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP171 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP173 + Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP174 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP6) + Eval(ee12,iintP8))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP58 + Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP61 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP64 + Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP65 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP7) + Eval(ee12,iintP9))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP58 + Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP61 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP64 + Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP65 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP10)*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP177 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP179 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP181 + Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP182 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP177 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP179 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP181 + Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP182 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP4)*((8*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP429 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP431 - Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP433 - Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP434 - Abb4*U))/SW) + 
     -     1/CW2*(Eval(ee13,iintP5)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup2*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup2*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup2*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup2*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee13,iintP7)*
     -     (8*Alfa2*MW*(AbbSumP71*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP68*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup2*EL*gL(1)**2*gL(3)**2)) + 
     -     8*MW*(dup4*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup2*EL*gL(3)**2*gR(1)**2)) + 
     -     8*MW*(dup4*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup2*EL*gL(1)**2*gR(3)**2)) + 
     -     8*MW*(dup4*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup2*EL*gR(1)**2*gR(3)**2)) + 
     -     8*MW*(dup4*EL*gR(1)**2*gR(3)**2)))))/SW + 
     -     1/CW2*(Eval(ee13,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP71*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP68*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup2*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup2*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup2*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup2*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup4*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP8)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP257 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(Abb72 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(Abb74 + Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP262 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP9)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP257 - Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(Abb72 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(Abb74 + Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP262 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP10)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP54 + Abb4*(T24 + U)))/SW) + 
     -     1/CW2*(Eval(ee14,iintP4)*
     -     (8*Alfa2*MW*(AbbSumP441*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP438*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP440*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP435*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*MH2) + 8*MW*T24)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*MH2) + 8*MW*T24)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*MH2) + 8*MW*T24)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*MH2) + 8*MW*T24))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee14,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP334*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*Alfa2*MW*(AbbSumP332*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP333*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP331*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP6)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 3*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP7)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 3*mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP8)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP272*EL*gL(1)*gL(3)*mass(1)*mass(3)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP9)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP272*EL*gL(1)**2*gL(3)**2*mass(1)*mass(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee2,iintP10)*
     -     (4*Alfa2*MW*(AbbSumP196*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP191*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP194*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP189*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) + 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee2,iintP11)*
     -     (4*Alfa2*MW*(AbbSumP196*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP191*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP194*EL*gL(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP189*EL*gR(1)**2*gR(3)**2) - 
     -     12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP4)*(-(24*Alfa2*MW)/CW2*
     -     (AbbSumP453*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP446*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP450*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP445*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(3*AbbSumP444 + Abb1*(S + T + T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(3*AbbSumP448 + Abb2*(S + T + T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(3*AbbSumP451 + Abb3*(S + T + T14)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(3*AbbSumP452 + Abb4*(S + T + T14)))/
     -     SW)


      dup8 = gL(1)*(gL(1) - 2*gR(1)) + gR(1)**2
      dup9 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup10 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      tmp4 = 1/CW2*(-12*Alfa2*MW*(AbbSumP37*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*MH2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*MH2) + MW*(4*S + 4*T14 + 4*T24))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*MH2) + MW*(4*S + 4*T14 + 4*T24))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*MH2) + MW*(4*S + 4*T14 + 4*T24))))/SW

      tmp4 = tmp4 + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(S + T14 + T24))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (AbbSumP352*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP339*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP347*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP337*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP351*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP342*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP345*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP338*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP336 + Abb1*(S34 + 2*T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP343 + Abb2*(S34 + 2*T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP348 + Abb3*(S34 + 2*T24 + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP350 + Abb4*(S34 + 2*T24 + 2*U)))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*mass(1)*q(1)*q(3)*(0 - 2/3.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*mass(1)*q(1)*q(3)*(0 - 3/2.D0*0))/
     -     SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP6)*((4*Alfa2*MW)/CW2*(Abb4*EL*
     -     gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP82 + Abb1*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP61 + Abb2*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP64 + Abb3*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP89 + Abb4*(-S + T)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP7)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP23*dup8*EL*gL(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*dup8*EL*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP82 + Abb1*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP61 + Abb2*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP64 + Abb3*(-S + T)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP89 + Abb4*(-S + T)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP8)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP253*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP238*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP248*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP236*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW -
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP274 + Abb1*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP276 + Abb2*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP278 + Abb3*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP279 + Abb4*(S34 + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP9)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP253*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP238*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP248*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP236*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP274 + Abb1*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP276 + Abb2*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP278 + Abb3*(S34 + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP279 + Abb4*(S34 + U)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*((Eval(ee22,iintP11) + Eval(ee22,iintP4))*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2)- 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*((Eval(ee11,iintP5) + Eval(ee22,iintP10))*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + (Eval(ee11,iintP10) + Eval(ee22,iintP5))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + 2*mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW) + 
     -     1/CW2*((Eval(ee22,iintP6) + Eval(ee22,iintP8))*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP7) + Eval(ee22,iintP9))*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*((Eval(ee23,iintP11) + Eval(ee23,iintP4))*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) - 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee23,iintP10)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee23,iintP6)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (8*MW*(dup10*EL*gL(1)*gL(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup10*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup10*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup10*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW + 
     -     1/CW2*(Eval(ee23,iintP5)*
     -     (8*Alfa2*MW*(AbbSumP359*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP356*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP358*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP354*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (8*MW*(dup10*EL*gL(1)*gL(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup10*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup10*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup10*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee23,iintP8)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee23,iintP7)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (8*MW*(dup10*EL*gL(1)**2*gL(3)**2) - 
     -     16*MW*(dup9*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup10*EL*gL(3)**2*gR(1)**2) - 
     -     16*MW*(dup9*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup10*EL*gL(1)**2*gR(3)**2) - 
     -     16*MW*(dup9*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup10*EL*gR(1)**2*gR(3)**2) - 
     -     16*MW*(dup9*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee23,iintP9)*
     -     (8*Alfa2*MW*(Abb55*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) - 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee24,iintP10)*
     -     (8*Alfa2*MW*(AbbSumP207*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP205*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP206*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP204*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP204 + Abb1*U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP205 + 2*Abb2*U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP206 + 2*Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP207 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP4)*
     -     (8*Alfa2*MW*(AbbSumP458*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP456*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(AbbSumP457*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP455*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP45 + Abb1*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP50 + Abb2*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP53 + Abb3*(T24 + U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP54 + Abb4*(T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee24,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP102*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb72*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb74*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP101*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP7)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP101 + Abb1*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(Abb72 + Abb2*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(Abb74 + Abb3*U))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP102 + Abb4*U))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP8)*
     -     (-(8*Alfa2*MW*(AbbSumP286*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP285*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP9)*
     -     (-(8*Alfa2*MW*(AbbSumP286*EL*gL(1)**2*gL(3)**2)) + 
     -     8*Alfa2*MW*(Abb59*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb18*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP285*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))
     -     /SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP10)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP11)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*U)/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP4)*
     -     (4*Alfa2*MW*(AbbSumP397*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP390*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP396*EL*gL(1)**2*gR(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP387*EL*gR(1)**2*gR(3)**2) - 
     -     12*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     12*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     12*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP6)*((4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW) + 
     -     Eval(ee3,iintP5)*((12*Alfa2*MW)/CW2*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(S + T + T14))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(S + T + T14))/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP7)*(-(4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*S)/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP8)*(-((28*Alfa2*MW)/CW2*
     -     (AbbSumP295*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW) - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP289*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP293*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (28*Alfa2*MW)/CW2*
     -     (AbbSumP287*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(S34 + T + U))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(S34 + T + U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*mass(1)*q(1)*q(3)*(0 - 1/2.D0*0))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP9)*
     -     (-28*Alfa2*MW*(AbbSumP295*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP23*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP289*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP293*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP15*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     28*Alfa2*MW*(AbbSumP287*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee33,iintP11) + Eval(ee33,iintP4))*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee33,iintP10)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee33,iintP5)*(Alfa2*(Abb4*
     -     (-(8*MW*(dup10*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(8*MW*(dup10*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(8*MW*(dup10*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(8*MW*(dup10*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW + 
     -     1/CW2*(Eval(ee33,iintP6)*(Alfa2*(Abb4*
     -     (-(8*MW*(dup10*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(8*MW*(dup10*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(8*MW*(dup10*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(8*MW*(dup10*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     16*MW*(dup9*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee33,iintP9)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee33,iintP8)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee33,iintP7)*(Alfa2*(Abb4*
     -     (-(8*MW*(dup10*EL*gL(1)**2*gL(3)**2)) + 
     -     16*MW*(dup9*EL*gL(1)**2*gL(3)**2))) + Alfa2*(Abb2*
     -     (-(8*MW*(dup10*EL*gL(3)**2*gR(1)**2)) + 
     -     16*MW*(dup9*EL*gL(3)**2*gR(1)**2))) + Alfa2*(Abb3*
     -     (-(8*MW*(dup10*EL*gL(1)**2*gR(3)**2)) + 
     -     16*MW*(dup9*EL*gL(1)**2*gR(3)**2))) + Alfa2*(Abb1*
     -     (-(8*MW*(dup10*EL*gR(1)**2*gR(3)**2)) + 
     -     16*MW*(dup9*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP10)*
     -     (-(8*Alfa2*MW*(AbbSumP213*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     4*Alfa2*MW*(Abb72*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(Abb74*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP212*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP11)*
     -     (-(8*Alfa2*MW*(AbbSumP213*EL*gL(1)**2*gL(3)**2)) + 
     -     4*Alfa2*MW*(Abb72*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(Abb74*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP212*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(8*MW*mass2(1)) - 8*MW*mass2(3) + 8*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(8*MW*mass2(1)) - 8*MW*mass2(3) + 8*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) - 8*MW*mass2(3) + 8*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(8*MW*mass2(1)) - 8*MW*mass2(3) + 8*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP4)*
     -     (-(8*Alfa2*MW*(AbbSumP466*EL*gL(1)**2*gL(3)**2)) + 
     -     8*Alfa2*MW*(AbbSumP463*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(AbbSumP465*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP461*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/
     -     SW + 1/CW2*(Eval(ee34,iintP5)*
     -     (-(8*Alfa2*MW*(AbbSumP378*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*Alfa2*MW*(AbbSumP375*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP377*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP373*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP7)*
     -     (8*Alfa2*MW*(AbbSumP108*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP107*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (8*MW*(dup10*EL*gL(1)**2*gL(3)**2) - 
     -     16*MW*(dup9*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup10*EL*gL(3)**2*gR(1)**2) - 
     -     16*MW*(dup9*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup10*EL*gL(1)**2*gR(3)**2) - 
     -     16*MW*(dup9*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup10*EL*gR(1)**2*gR(3)**2) - 
     -     16*MW*(dup9*EL*gR(1)**2*gR(3)**2)))))/SW + 
     -     1/CW2*(Eval(ee34,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP108*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP107*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (8*MW*(dup10*EL*gL(1)*gL(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup10*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup10*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup10*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     16*MW*(dup9*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP9)*
     -     (8*Alfa2*MW*(AbbSumP299*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP298*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW
     -     + 1/CW2*(Eval(ee34,iintP8)*
     -     (8*Alfa2*MW*(AbbSumP299*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP298*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP10)*((4*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)*gL(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee4,iintP11)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*U)/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW-
     -     (12*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW-
     -     (12*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW-
     -     (12*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP4)*
     -     (8*Alfa2*MW*(AbbSumP475*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP470*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP473*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP467*EL*gR(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))))/
     -     SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP5)*
     -     (12*Alfa2*MW*(AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP6)*(-(4*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*mass(3)*q(1)*q(3)*(0 - 1/2.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*mass(3)*q(1)*q(3)*(0 - 1/2.D0*0))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP7)*
     -     (4*Alfa2*MW*(AbbSumP32*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + tmp4*Eval(ee4,iintP9) + 
     -     1/CW2*(Eval(ee44,iintP11)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP10)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee4,iintP8)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP37*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP29*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP9*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP23*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP15*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24))))))/SW


      dup11 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup12 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + Eval(ee44,iintP4)*(-((16*Alfa2*MW)/CW2*
     -     (Abb4*dup11*EL*gL(1)**2*gL(3)**2)/SW) + (8*Alfa2*MW)/CW2*
     -     ((AbbSumP483 + Abb4*dup12)*EL*gL(1)**2*gL(3)**2)/SW - 
     -     (16*Alfa2*MW)/CW2*(Abb2*dup11*EL*gL(3)**2*gR(1)**2)/SW - 
     -     (16*Alfa2*MW)/CW2*(Abb3*dup11*EL*gL(1)**2*gR(3)**2)/SW - 
     -     (16*Alfa2*MW)/CW2*(Abb1*dup11*EL*gR(1)**2*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP477 + Abb1*dup12)*EL*gR(1)**2*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP479 + 
     -     Abb2*(8*mass2(1) + 6*mass2(3) + 4*MH2)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP481 + 
     -     Abb3*(8*mass2(1) + 6*mass2(3) + 4*MH2)))/SW) +
     -     1/CW2*(Eval(ee44,iintP6)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee44,iintP7)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW 


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P1 for:
**     quark + antiquark -> quark + H + antiquark 

      subroutine pentP1_qqbar(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, dup2, dup3, dup4, dup5
      double complex dup6, dup7, dup8, dup9, dup10
      double complex dup11, dup12, dup13, dup14, dup15
      double complex dup16
      double complex tmp1, tmp2, tmp3


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup1 = gL(1)*(gL(1) - 2*gR(1)) + gR(1)**2
      dup2 = gL(3)*(gL(3) - 2*gR(3)) + gR(3)**2
      dup3 = MZ2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)
      dup4 = AbbSumP191 + Abb2*(-S - 2*T - T14)
      dup5 = AbbSumP65 + Abb2*(-S - 2*T - T14)
      dup6 = AbbSumP193 + Abb3*(-S - 2*T - T14)
      dup7 = AbbSumP68 + Abb3*(-S - 2*T - T14)

      tmp1 = 1/CW2*(12*Alfa2*MW*(AbbSumP40*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(AbbSumP41*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP20*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     8*Alfa2*MW*(AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP27*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     16*Alfa2*MW*(AbbSumP9*EL*gR(1)**2*gR(3)**2*mass2(3)) + 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp1 = tmp1 + (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP8 + Abb1*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP21 + Abb2*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP32 + Abb3*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP39 + Abb4*(2*S + T24)))/SW

      tmp2 = -(16*Alfa2*MW)/CW2*
     -     (AbbSumP281*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP268*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP277*EL*gL(1)**2*gR(3)**2*mass2(1))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP266*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP282*EL*gL(1)**2*gL(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP271*EL*gL(3)**2*gR(1)**2*mass2(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP273*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP265*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP264 + Abb1*(-T14 - 2*U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP272 + Abb2*(-T14 - 2*U)))/SW

      tmp2 = tmp2 - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP278 + Abb3*(-T14 - 2*U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP280 + Abb4*(-T14 - 2*U)))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Dval(dd0,iintP3)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*(Dval(dd0,iintP2)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd0,iintP1)*
     -     (Alfa2*(Abb4*EL*gL(1)*
     -     (8*MW*(gL(1)*gL(3)**2) + 8*MW*(gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*
     -     (8*MW*(gL(3)*gR(1)**2) + 8*MW*(gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*
     -     (8*MW*(gL(1)*gR(3)**2) + 8*MW*(gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*
     -     (8*MW*(gR(1)*gR(3)**2) + 8*MW*(gR(3)*q(1)*q(3))))))/SW + 
     -     1/CW2*(Eval(ee0,iintP10)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee0,iintP11)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSumP404*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP398*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP403*EL*gL(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(AbbSumP393*EL*gR(1)**2*gR(3)**2) + 
     -     12*Alfa2*MW*(AbbSumP406*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP396*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP401*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumP395*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee0,iintP6)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/CW2*(Eval(ee0,iintP5)*
     -     (12*Alfa2*MW*(AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP7)*(Alfa2*(Abb4*
     -     (8*MW*(dup3*EL*gL(1)**2*gL(3)**2) - 
     -     4*MW*(EL*gL(1)**2*gL(3)**2*S))) + 
     -     Alfa2*(Abb2*
     -     (8*MW*(dup3*EL*gL(3)**2*gR(1)**2) - 
     -     4*MW*(EL*gL(3)**2*gR(1)**2*S))) + 
     -     Alfa2*(Abb3*
     -     (8*MW*(dup3*EL*gL(1)**2*gR(3)**2) - 
     -     4*MW*(EL*gL(1)**2*gR(3)**2*S))) + 
     -     Alfa2*(Abb1*
     -     (8*MW*(dup3*EL*gR(1)**2*gR(3)**2) - 
     -     4*MW*(EL*gR(1)**2*gR(3)**2*S)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP8)*(-(12*Alfa2*MW*
     -     (AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP9)*
     -     (-12*Alfa2*MW*(AbbSumP40*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*((Eval(ee00,iintP7) + 
     -     Eval(ee00,iintP9))*
     -     (-(24*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP11) + Eval(ee00,iintP4))*
     -     (-(32*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) - 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP6) + Eval(ee00,iintP8))*
     -     (-(24*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*((Eval(ee00,iintP10) + Eval(ee00,iintP5))*
     -     (-(32*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee1,iintP10)*
     -     (-(4*Alfa2*MW*(AbbSumP176*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP165*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP172*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP159*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*U) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP33*EL*gL(1)*
     -     (-(4*MW*(gL(3)*mass2(3)*q(1)*q(3))) - 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP168*EL*gL(1)*
     -     (-(4*MW*(gL(3)*mass2(3)*q(1)*q(3))) + 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP13*EL*
     -     (-(4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP160*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee1,iintP11)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumP40*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP17*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP31*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP10*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP168*EL*gL(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP160*EL*gR(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP159 - 2*Abb1*U))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP165 - 2*Abb2*U))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP172 - 2*Abb3*U))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP176 - 2*Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP4)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP25*dup1*EL*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*dup1*EL*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP33*dup2*EL*gL(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*dup2*EL*gR(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP408 + Abb1*(-S + S34 - T14 + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP415 + Abb2*(-S + S34 - T14 + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP421 + Abb3*(-S + S34 - T14 + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP423 + Abb4*(-S + S34 - T14 + T24)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP5)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP346*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP339*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP342*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP337*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW + (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP340 + Abb2*(-S - 2*T - T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP343 + Abb3*(-S - 2*T - T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP335 + 
     -     Abb1*(S + 2*T + T14)))/SW - (4*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*
     -     q(1)*q(3)*(AbbSumP344 + Abb4*(S + 2*T + T14)))/SW)

      Cloop(1) = Cloop(1) + tmp1*Eval(ee1,iintP7) + 
     -     Eval(ee1,iintP6)*((12*Alfa2*MW)/CW2*
     -     (AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP41*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP20*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP27*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP8 + Abb1*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP21 + Abb2*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP32 + Abb3*(2*S + T24)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP39 + Abb4*(2*S + T24)))/SW)

      Cloop(1) = Cloop(1) + tmp2*Eval(ee1,iintP9) + 
     -     Eval(ee1,iintP8)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP281*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP268*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP277*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP266*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP282*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP271*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP273*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP265*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP264 + Abb1*(-T14 - 2*U)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP272 + Abb2*(-T14 - 2*U)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP278 + Abb3*(-T14 - 2*U)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP280 + Abb4*(-T14 - 2*U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee11,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP58 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee11,iintP4)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP426 + Abb1*(S34 + T24)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP430 + Abb2*(S34 + T24)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP432 + Abb3*(S34 + T24)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP433 + Abb4*(S34 + T24)))/SW)+
     -     1/CW2*((Eval(ee11,iintP6) + Eval(ee11,iintP8))*
     -     (-(8*Alfa2*MW*(AbbSumP58*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumP54*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP57*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP50*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee11,iintP7) + Eval(ee11,iintP9))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP58 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP11)*
     -     ((8*Alfa2*MW)/CW2*(dup4*EL*gL(3)**2*gR(1)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(dup6*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP188 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP194 + Abb4*(S + 2*T + T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + Eval(ee12,iintP4)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP434 + Abb1*(-S - T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP437 + Abb2*(-S - T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP439 + Abb3*(-S - T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP440 + Abb4*(-S - T14)))/SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP10) + Eval(ee12,iintP5))*
     -     ((8*Alfa2*MW)/CW2*(dup4*EL*gL(3)*gR(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(dup6*EL*gL(1)*gR(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP188 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP194 + Abb4*(S + 2*T + T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP6) + Eval(ee12,iintP8))*
     -     ((8*Alfa2*MW)/CW2*(dup5*EL*gL(3)*gR(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(dup7*EL*gL(1)*gR(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP62 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP69 + Abb4*(S + 2*T + T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP7) + Eval(ee12,iintP9))*
     -     ((8*Alfa2*MW)/CW2*(dup5*EL*gL(3)**2*gR(1)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(dup7*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP62 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP69 + Abb4*(S + 2*T + T14)))/
     -     SW) + 1/CW2*(Eval(ee13,iintP10)*
     -     (-(8*Alfa2*MW*(AbbSumP202*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP199*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP201*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP196*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T) + 
     -     Alfa2*(AbbSumP98*EL*
     -     (8*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP96*EL*
     -     (-(8*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) + 
     -     8*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP11)*(-(8*Alfa2*MW)/CW2*
     -     (AbbSumP98*EL*gL(3)**2*(-gL(1)**2 + gR(1)**2)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP96*EL*(-gL(1)**2 + gR(1)**2)*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP196 + Abb1*T))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP199 + 2*Abb2*T))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP201 + 2*Abb3*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP202 + Abb4*T))/SW) + 
     -     1/CW2*(Eval(ee13,iintP4)*
     -     (8*Alfa2*MW*(Abb65*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(Abb58*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(Abb62*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb51*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (8*MW*(gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*MW*(gL(3)**2*gR(1)**2*mass2(1)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (8*MW*(gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*MW*(gR(1)**2*gR(3)**2*mass2(1))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP6)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP72 + Abb1*
     -     (-S - T - T14)))/SW + (8*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)*q(1)*
     -     q(3)*(AbbSumP78 + Abb2*(-S - T - T14)))/SW + (8*Alfa2*MW)/
     -     CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP82 + Abb3*(-S - T - 
     -     T14)))/SW + (8*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP83 + Abb4*(-S - T - T14)))/SW) + 1/CW2*
     -     (Eval(ee13,iintP5)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T)))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP7)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP72 + Abb1*(-S - T - T14)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP78 + Abb2*(-S - T - T14)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP82 + Abb3*(-S - T - T14)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP83 + Abb4*(-S - T - T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee13,iintP9)*(-(8*Alfa2*MW)/CW2*
     -     (AbbSumP98*EL*gL(3)**2*(-gL(1)**2 + gR(1)**2)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP96*EL*(-gL(1)**2 + gR(1)**2)*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP285 + Abb1*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(Abb55 + Abb2*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(Abb23 + Abb3*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP289 + Abb4*T))/SW) + 
     -     1/CW2*(Eval(ee13,iintP8)*
     -     (-(8*Alfa2*MW*(AbbSumP289*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP285*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T) + 
     -     Alfa2*(AbbSumP98*EL*
     -     (8*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP96*EL*
     -     (-(8*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) + 
     -     8*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee14,iintP10)*
     -     (-(8*Alfa2*MW*(AbbSumP58*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumP54*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP57*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP50*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(0*EL*gL(1)*gL(3)*mass(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP11)*
     -     (-(8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP58 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP4)*((8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP447 + Abb1*(-2*S34 - 2*T24 - U)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP454 + Abb4*(-2*S34 - 2*T24 - U)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP451 + Abb2*(2*S34 + 2*T24 + U)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP453 + Abb3*(2*S34 + 2*T24 + U)))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP5)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP350 + Abb2*(-2*S - 4*T - 2*T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP351 + Abb3*(-2*S - 4*T - 2*T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP349 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP352 + Abb4*(S + 2*T + T14)))/
     -     SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee14,iintP6)*
     -     (-(8*Alfa2*MW*(AbbSumP58*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumP54*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP57*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP50*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP7)*(-(8*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP58 - Abb4*U))/SW)

      Cloop(1) = Cloop(1) + Eval(ee14,iintP8)*(-(8*Alfa2*MW)/CW2*(EL*
     -     gR(1)*gR(3)*q(1)*q(3)*(AbbSumP50 - Abb1*U))/SW - (8*Alfa2*
     -     MW)/CW2*(EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP57 - Abb3*
     -     U))/SW - (8*Alfa2*MW)/CW2*(EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP58 - Abb4*U))/SW - (16*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*
     -     gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW - (16*Alfa2*
     -     MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*
     -     mass2(1)))/SW - (16*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*
     -     q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW - (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/
     -     SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP10)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP208 + Abb1*(T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP215 + Abb4*(T - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP212 + Abb2*(-T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP214 + Abb3*(-T + U)))/SW) + 
     -     Eval(ee14,iintP9)*(-(8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP58 - Abb4*U))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP11)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP25*dup1*EL*gL(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*dup1*EL*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP208 + Abb1*(T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP215 + Abb4*(T - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP212 + Abb2*(-T + U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP214 + Abb3*(-T + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP4)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (3*AbbSumP456 + Abb1*(-S + 2*T - T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP459 + Abb2*(-S + 2*T - T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP461 + Abb3*(-S + 2*T - T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (3*AbbSumP462 + Abb4*(-S + 2*T - T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP5)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP355 + Abb1*(S + T + T14 - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP358 + Abb2*(S + T + T14 - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP361 + Abb3*(S + T + T14 - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*
     -     (AbbSumP362 + Abb4*(S + T + T14 - U)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee2,iintP6)*
     -     (-(4*Alfa2*MW*(AbbSumP108*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     4*Alfa2*MW*(AbbSumP100*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP105*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP94*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP98*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP96*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S)))))/SW


      dup8 = S - S34 + T14 - T24 - U
      dup9 = gR(1)*mass2(1) + gL(1)*(-mass2(3) - MH2)
      dup10 = gR(1)**2*mass2(1) + gL(1)**2*(-mass2(3) - MH2)
      dup11 = -(gL(1)*mass2(1)) + gR(1)*(mass2(3) + MH2)
      dup12 = -(gL(1)**2*mass2(1)) + gR(1)**2*(mass2(3) + MH2)
      dup13 = AbbSumP65 + Abb2*(-S - 2*T - T14)
      dup14 = AbbSumP68 + Abb3*(-S - 2*T - T14)
      dup15 = AbbSumP299 + Abb1*(-S34 + T - U)
      dup16 = AbbSumP303 + Abb4*(-S34 + T - U)

      tmp3 = -(12*Alfa2*MW)/CW2*
     -     (AbbSumP40*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP17*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP31*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/
     -     SW - (12*Alfa2*MW)/CW2*
     -     (AbbSumP10*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (2*mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (2*mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*(S + T14 + T24))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(S + T14 + T24))/SW

      tmp3 = tmp3 + 1/CW2*
     -     (Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*S + 4*T14 + 4*T24)) + 
     -     Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*(4*S + 4*T14 + 4*T24))
     -     )/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP7)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP98*EL*gL(3)**2*(-gL(1)**2 + gR(1)**2)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP96*EL*(-gL(1)**2 + gR(1)**2)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP94 + Abb1*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP100 - Abb2*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP105 - Abb3*S))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP108 + Abb4*S))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP8)*(-((4*Alfa2*MW)/CW2*
     -     (Abb2*dup11*EL*gL(3)*q(1)*q(3))/SW) + 
     -     (4*Alfa2*MW)/CW2*(Abb4*dup9*EL*gL(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(dup16*EL*gL(1)*gL(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup11*EL*gR(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup9*EL*gR(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(dup15*EL*gR(1)*gR(3)*q(1)*q(3))/SW -
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP65 + Abb2*(S34 - T + U)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP68 + Abb3*(S34 - T + U)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP9)*((4*Alfa2*MW)/CW2*
     -     (Abb4*dup10*EL*gL(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup12*EL*gL(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(dup16*EL*gL(1)**2*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup10*EL*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup12*EL*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(dup15*EL*gR(1)**2*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP65 + Abb2*(S34 - T + U)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP68 + Abb3*(S34 - T + U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP11) + Eval(ee22,iintP4))*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T))))/SW + 
     -     1/CW2*((Eval(ee11,iintP5) + Eval(ee22,iintP10))*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP6) + Eval(ee22,iintP8))*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T)))))/SW + 
     -     1/CW2*((Eval(ee11,iintP10) + Eval(ee22,iintP5))*
     -     (-(8*Alfa2*MW*(AbbSumP58*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumP54*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP57*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP50*EL*gR(1)*gR(3)*q(1)*q(3))+ 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(8*MW*(mass2(1)*q(1)*q(3))) - 
     -     8*MW*(mass2(3)*q(1)*q(3)) + 8*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP7) + Eval(ee22,iintP9))*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (8*MW*mass2(1) - 8*MW*T))))/SW + 
     -     1/CW2*(Eval(ee23,iintP10)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee23,iintP11) + Eval(ee23,iintP4))*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP5)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP365 + Abb1*(-S - T - T14)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP368 + Abb2*(-S - T - T14)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP370 + Abb3*(-S - T - T14)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP371 + Abb4*(-S - T - T14)))/
     -     SW) + 1/CW2*(Eval(ee23,iintP6)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T)))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee23,iintP7)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW + 
     -     1/CW2*(Eval(ee23,iintP9)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (16*MW*mass2(1) - 8*MW*T))))/SW + 
     -     1/CW2*(Eval(ee23,iintP8)*
     -     (-(8*Alfa2*MW*(Abb63*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb49*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (16*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*T)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP11)*
     -     (8*Alfa2*MW*(AbbSumP229*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP227*EL*gL(3)**2*gR(1)**2) + 
     -     4*Alfa2*MW*(AbbSumP228*EL*gL(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP226*EL*gR(1)**2*gR(3)**2)))/SW
     -     + 1/CW2*(Eval(ee24,iintP10)*
     -     (8*Alfa2*MW*(AbbSumP229*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP227*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP228*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP226*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP4)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP465 + Abb2*(-2*S - 2*T14)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP466 + Abb3*(-2*S - 2*T14)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP464 + Abb1*(S + T14)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP467 + Abb4*(S + T14)))/SW)

      Cloop(1) = Cloop(1) + Eval(ee24,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP50 - Abb1*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP54 - Abb2*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP57 - Abb3*U))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP58 - Abb4*U))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)
     -     + 1/CW2*(Eval(ee24,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP120*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP119*EL*gR(1)*gR(3)*q(1)*q(3))))
     -     /SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP8)*((8*Alfa2*MW)/CW2*
     -     (dup13*EL*gL(3)*gR(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(dup14*EL*gL(1)*gR(3)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(3*mass2(1) + mass2(3) + MH2)*q(1)*q(3))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP309 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP310 + Abb4*(S + 2*T + T14)))/
     -     SW) + 1/CW2*(Eval(ee24,iintP7)*
     -     (8*Alfa2*MW*(AbbSumP120*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP119*EL*gR(1)**2*gR(3)**2)))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP9)*((8*Alfa2*MW)/CW2*
     -     (dup13*EL*gL(3)**2*gR(1)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(dup14*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP309 + Abb1*(S + 2*T + T14)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP310 + Abb4*(S + 2*T + T14)))/
     -     SW) + Eval(ee3,iintP10)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP11)*
     -     (-4*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP4)*
     -     (-(4*Alfa2*MW*(AbbSumP404*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP398*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP403*EL*gL(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(AbbSumP393*EL*gR(1)**2*gR(3)**2) + 
     -     12*Alfa2*MW*(AbbSumP406*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP396*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP401*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumP395*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP5)*(Alfa2*(AbbSumP25*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP6)*((4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP7)*
     -     (4*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP8)*
     -     (-20*Alfa2*MW*(AbbSumP317*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP314*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP315*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     20*Alfa2*MW*(AbbSumP312*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP9)*
     -     (-20*Alfa2*MW*(AbbSumP317*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP314*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP315*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     20*Alfa2*MW*(AbbSumP312*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) +1/CW2*(Eval(ee33,iintP7)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW +
     -     1/CW2*(Eval(ee33,iintP5)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee33,iintP6)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee34,iintP11)*
     -     (-(8*Alfa2*MW*(Abb80*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(Abb71*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(Abb73*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb78*EL*gR(1)**2*gR(3)**2)))/SW
     -     + 1/CW2*(Eval(ee34,iintP10)*
     -     (-(8*Alfa2*MW*(Abb80*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(Abb71*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb73*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb78*EL*gR(1)*gR(3)*q(1)*q(3))))/SW
     -     + 1/CW2*(Eval(ee34,iintP4)*
     -     (-(8*Alfa2*MW*(AbbSumP475*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP473*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP474*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP472*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (8*MW*(gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*MW*(gL(3)**2*gR(1)**2*mass2(1)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (8*MW*(gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*MW*(gR(1)**2*gR(3)**2*mass2(1))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP5)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP85*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP76*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP81*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP74*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP78 + Abb2*(-2*S - 2*T - 2*T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP82 + Abb3*(-2*S - 2*T - 2*T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP385 + Abb1*(S + T + T14)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP386 + Abb4*(S + T + T14)))/SW
     -     ) + 1/CW2*(Eval(ee34,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP128*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP127*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP7)*
     -     (8*Alfa2*MW*(AbbSumP128*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP127*EL*gR(1)**2*gR(3)**2)))/SW+
     -     1/CW2*(Eval(ee34,iintP8)*
     -     (-(8*Alfa2*MW*(AbbSumP320*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb55*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb23*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP319*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*T) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*T) + 
     -     Alfa2*(AbbSumP98*EL*
     -     (8*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     8*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP96*EL*
     -     (-(8*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) + 
     -     8*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP9)*(-(8*Alfa2*MW)/CW2*
     -     (AbbSumP98*EL*gL(3)**2*(-gL(1)**2 + gR(1)**2)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP96*EL*(-gL(1)**2 + gR(1)**2)*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP319 + Abb1*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(Abb55 + Abb2*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(Abb23 + Abb3*T))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP320 + Abb4*T))/SW) + 
     -     Eval(ee4,iintP10)*((4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*
     -     q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP11)*
     -     (-4*Alfa2*MW*(AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP4)*((4*Alfa2*MW)/CW2*
     -     ((2*AbbSumP486 + Abb4*dup8)*EL*gL(1)**2*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     ((2*AbbSumP476 + Abb1*dup8)*EL*gR(1)**2*gR(3)**2)/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP406*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP25*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP396*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP401*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP14*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP395*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP481 + Abb2*(-S + S34 - T14 + T24 + U)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP484 + Abb3*(-S + S34 - T14 + T24 + U)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP5)*
     -     (12*Alfa2*MW*(AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP6)*(-(4*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee4,iintP7)*(-(4*Alfa2*MW)/CW2*
     -     (AbbSumP33*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP13*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*S)/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW+
     -     (12*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + 1/3.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + tmp3*Eval(ee4,iintP9) + 
     -     1/CW2*(Eval(ee4,iintP8)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP40*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP17*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP31*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP10*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP25*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP14*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S + 4*T14 + 4*T24))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee44,iintP4)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(2*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(2*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(2*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(2*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*(EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP488 + Abb1*(-S34 - T24 - U)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP495 + Abb4*(-S34 - T24 - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP492 + Abb2*(2*S34 + 2*T24 + 2*U)))/SW - 
     -     (4*Alfa2*MW)/CW2*(EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP494 + Abb3*(2*S34 + 2*T24 + 2*U)))/SW) + 
     -     1/CW2*(Eval(ee44,iintP7)*
     -     (8*Alfa2*MW*(0*EL*gL(1)**2*gL(3)**2*mass(3)) + 
     -     8*Alfa2*MW*(0*EL*gR(1)**2*gR(3)**2*mass(3))))/SW+
     -     1/CW2*(Eval(ee44,iintP5)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee44,iintP9)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP8)*
     -     (-8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P1 for:
**     antiquark + antiquark -> antiquark + H + antiquark 

      subroutine pentP1_qbarqbar(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** Dummy variables
      double complex dup1, dup2, dup3, dup4, dup5
      double complex dup6, dup7, dup8, dup9, dup10
      double complex dup11, dup12, dup13, dup14, dup15
      double complex dup16, dup17, dup18, dup19, dup20
      double complex dup21, dup22, dup23, dup24, dup25
      double complex dup26
      double complex tmp1, tmp2, tmp3, tmp4


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2

*      write(*,*)'new NC qk=4 pent1'


      dup1 = gL(1)*(gL(1) - 2*gR(1)) + gR(1)**2
      dup2 = gL(3)*(gL(3) - 2*gR(3)) + gR(3)**2
      dup3 = S - S34 + 3*T14 + 2*T24
      dup4 = 2*S34 + 2*T + T24 + 2*U
      dup5 = AbbSumP161 + Abb2*(S34 + T + U + 1/2.D0*T24)
      dup6 = AbbSumP170 + Abb3*(S34 + T + U + 1/2.D0*T24)
      dup7 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup8 = 4*mass2(1) + 3*mass2(3) + 2*MH2
      dup9 = gR(1)**2*mass2(1) - 3/2.D0*(gL(1)**2*MH2)
      dup10 = gR(1)**2*mass2(1) + gL(1)**2*(-2*mass2(3) - MH2)
      dup11 = gR(1)**2*mass2(1) + gL(1)**2*(-mass2(3) - MH2)
      dup12 = gR(1)**2*MH2 - 2/3.D0*(gL(1)**2*mass2(1))
      dup13 = -(gL(1)**2*mass2(1)) + gR(1)**2*(mass2(3) + MH2)
      dup14 = -(gL(1)**2*mass2(1)) + gR(1)**2*(2*mass2(3) + MH2)
      dup15 = AbbSumP271 + Abb1*(-S + S34 + T + U)
      dup16 = 2*AbbSumP276 + Abb2*(-S + S34 + T + U)
      dup17 = 2*AbbSumP280 + Abb3*(-S + S34 + T + U)
      dup18 = AbbSumP281 + Abb4*(-S + S34 + T + U)

      tmp1 = 1/CW2*(-20*Alfa2*MW*(AbbSumP177*EL*gL(1)**2*gL(3)**2*
     -     mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP156*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP169*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     20*Alfa2*MW*(AbbSumP152*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP176*EL*gL(1)**2*gL(3)**2*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP160*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     8*Alfa2*MW*(AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) - 
     -     4*Alfa2*MW*(AbbSumP166*EL*gL(1)**2*gR(3)**2*mass2(3)) - 
     -     16*Alfa2*MW*(AbbSumP153*EL*gR(1)**2*gR(3)**2*mass2(3)) - 
     -     12*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) - 
     -     12*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) - 
     -     12*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2))/SW

      tmp1 = tmp1 + 1/CW2*
     -     (4*Alfa2*MW*(AbbSumP174*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(Abb4*dup4*EL*gL(1)**2*gL(3)**2) + 
     -     4*Alfa2*MW*(AbbSumP150*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(8*MW*(dup5*EL*gL(3)**2*gR(1)**2) + 
     -     8*MW*(dup6*EL*gL(1)**2*gR(3)**2)) + 
     -     Alfa2*(Abb1*(4*MW*(dup4*EL*gR(1)**2*gR(3)**2) - 
     -     12*MW*(EL*gR(1)**2*gR(3)**2*MH2))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Dval(dd0,iintP3) + Dval(dd0,iintP4))*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Dval(dd0,iintP1) + Dval(dd0,iintP2))*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP10)*
     -     (-12*Alfa2*MW*(AbbSumP35*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP12)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(1)) - 4*MW*mass2(3) + 4*MW*U))))/SW
     -     + 1/CW2*(Eval(ee0,iintP11)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSumP408*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP400*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP407*EL*gL(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(AbbSumP395*EL*gR(1)**2*gR(3)**2) + 
     -     12*Alfa2*MW*(AbbSumP409*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP398*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP406*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     16*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     12*Alfa2*MW*(AbbSumP396*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee0,iintP7)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) - 
     -     4*MW*(q(1)*q(3)*S)))))/SW + 
     -     1/CW2*(Eval(ee0,iintP6)*
     -     (12*Alfa2*MW*(AbbSumP35*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee0,iintP8)*
     -     (Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee00,iintP12) + Eval(ee00,iintP5))*
     -     (24*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP10) + Eval(ee00,iintP8))*
     -     (32*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2) + 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2)))/SW + 
     -     1/CW2*((Eval(ee00,iintP11) + Eval(ee00,iintP6))*
     -     (24*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     24*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*((Eval(ee00,iintP7) + Eval(ee00,iintP9))*
     -     (32*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     32*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee0,iintP9)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP35*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(3)*q(1)*q(3))) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee1,iintP10)*
     -     (-8*Alfa2*MW*(AbbSumP211*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP202*EL*gL(3)**2*gR(1)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP208*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) - 
     -     8*Alfa2*MW*(AbbSumP201*EL*gR(1)**2*gR(3)**2*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP28*dup2*EL*gL(1)**2*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP9*dup2*EL*gR(1)**2*mass2(3)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*MH2) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*MH2) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*MH2) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*MH2) + 
     -     Alfa2*(4*MW*(dup18*EL*gL(1)**2*gL(3)**2) + 
     -     4*MW*(dup16*EL*gL(3)**2*gR(1)**2) + 
     -     4*MW*(dup17*EL*gL(1)**2*gR(3)**2) + 
     -     4*MW*(dup15*EL*gR(1)**2*gR(3)**2))))/SW

      Cloop(1) = Cloop(1) + tmp1*Eval(ee1,iintP12) + 
     -     1/CW2*(Eval(ee1,iintP11)*
     -     (4*Alfa2*MW*(AbbSumP174*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP150*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     20*Alfa2*MW*(AbbSumP177*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP156*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP169*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     20*Alfa2*MW*(AbbSumP152*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(AbbSumP176*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP160*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP166*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     16*Alfa2*MW*(AbbSumP153*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3)) - 
     -     12*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3)) + 
     -     Alfa2*(8*MW*(dup5*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*MW*(dup6*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP28*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP9*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*
     -     (4*MW*(dup4*EL*gL(1)*gL(3)*q(1)*q(3)) - 
     -     12*MW*(EL*gL(1)*gL(3)*MH2*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (4*MW*(dup4*EL*gR(1)*gR(3)*q(1)*q(3)) - 
     -     12*MW*(EL*gR(1)*gR(3)*MH2*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee1,iintP5)*(-((12*Alfa2*MW)/CW2*
     -     (Abb2*dup12*EL*gL(3)**2)/SW) + 
     -     (8*Alfa2*MW)/CW2*(Abb4*dup9*EL*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     ((AbbSumP432 + Abb4*dup3)*EL*gL(1)**2*gL(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     ((AbbSumP422 - Abb2*dup3)*EL*gL(3)**2*gR(1)**2)/SW - 
     -     (12*Alfa2*MW)/CW2*(Abb1*dup12*EL*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*dup9*EL*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     ((AbbSumP429 - Abb3*dup3)*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     ((AbbSumP414 + Abb1*dup3)*EL*gR(1)**2*gR(3)**2)/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP433*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP421*EL*gL(3)**2*gR(1)**2*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP426*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP415*EL*gR(1)**2*gR(3)**2*mass2(3))/SW)

      Cloop(1) = Cloop(1) + Eval(ee1,iintP6)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP339*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP334*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP337*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP333*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP332 + Abb1*(S + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP193 + Abb2*(S + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP195 + Abb3*(S + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP338 + Abb4*(S + T14)))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee1,iintP7)*
     -     (4*Alfa2*MW*(AbbSumP34*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP16*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP27*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP4*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP35*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP23*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) - 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP28*EL*gL(1)*
     -     (4*MW*(gL(3)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP6*EL*
     -     (-(4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3))) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP9*EL*
     -     (4*MW*(gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee1,iintP8)*((12*Alfa2*MW)/CW2*
     -     (AbbSumP35*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP12*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP26*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP5*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP23*EL*gL(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP6*EL*gR(1)**2*(-gL(3)**2 + gR(3)**2)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP4 - 2*Abb1*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP16 - 2*Abb2*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP27 - 2*Abb3*S))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP34 - 2*Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee1,iintP9)*
     -     (-8*Alfa2*MW*(AbbSumP211*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP202*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP208*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP201*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))- 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*MH2*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*MH2*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*MH2*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*MH2*q(1)*q(3)) + 
     -     Alfa2*(4*MW*(dup18*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     4*MW*(dup16*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     4*MW*(dup17*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     4*MW*(dup15*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP12)*
     -     ((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP52 - Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee11,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP437 + Abb1*(T14 + T24)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP441 + Abb2*(T14 + T24)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP443 + Abb3*(T14 + T24)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP444 + Abb4*(T14 + T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee11,iintP10) + Eval(ee11,iintP8))*
     -     ((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP52 - Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee11,iintP7) + Eval(ee11,iintP9))*
     -     (8*Alfa2*MW*(AbbSumP52*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP48*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP51*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP43*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP12)*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP191 + Abb1*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP193 + Abb2*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP195 + Abb3*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP196 + Abb4*(-S34 - T - U)))/
     -     SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee12,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + 2*mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP447 + Abb1*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP400 + Abb2*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP407 + Abb3*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP450 + Abb4*(-S - T14 - T24)))/SW)

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee12,iintP11) + Eval(ee12,iintP6))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP191 + Abb1*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP193 + Abb2*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP195 + Abb3*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP196 + Abb4*(-S34 - T - U)))/
     -     SW)

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee12,iintP10) + Eval(ee12,iintP8))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP56 + Abb1*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP59 + Abb2*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP62 + Abb3*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP63 + Abb4*(-S34 - T - U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     (Eval(ee12,iintP7) + Eval(ee12,iintP9))*
     -     (-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP56 + Abb1*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP59 + Abb2*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP62 + Abb3*(-S34 - T - U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP63 + Abb4*(-S34 - T - U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP10)*
     -     (-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP286 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP288 + Abb2*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP290 + Abb3*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP291 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP11)*
     -     (-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP200 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP205 + Abb2*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP209 + Abb3*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP210 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP12)*
     -     (-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP200 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP205 + Abb2*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP209 + Abb3*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP210 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP5)*((8*Alfa2*MW)/CW2*
     -     (Abb4*dup10*EL*gL(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*dup14*EL*gL(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*dup10*EL*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*dup14*EL*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP452 + Abb1*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP400 + Abb2*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP407 + Abb3*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP456 + Abb4*(-S - T14 - T24)))/SW) + 
     -     1/CW2*(Eval(ee13,iintP6)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup7*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup8*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup7*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup8*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup7*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup8*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup7*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup8*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP7)*(-((8*Alfa2*MW)/CW2*
     -     ((AbbSumP74 - Abb4*dup8)*EL*gL(1)*gL(3)*q(1)*q(3))/SW) - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP70 - Abb2*dup8)*EL*gL(3)*gR(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP73 - Abb3*dup8)*EL*gL(1)*gR(3)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP66 - Abb1*dup8)*EL*gR(1)*gR(3)*q(1)*q(3))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP75*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP68*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP72*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP67*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP8)*(-((8*Alfa2*MW)/CW2*
     -     ((AbbSumP74 - Abb4*dup8)*EL*gL(1)**2*gL(3)**2)/SW) - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP70 - Abb2*dup8)*EL*gL(3)**2*gR(1)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP73 - Abb3*dup8)*EL*gL(1)**2*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP66 - Abb1*dup8)*EL*gR(1)**2*gR(3)**2)/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP75*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP68*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP72*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP67*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(MH2 + 3/2.D0*mass2(3)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee13,iintP9)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP286 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP288 + Abb2*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP290 + Abb3*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP291 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP10)*
     -     ((8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP52 - Abb4*S))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee14,iintP11)*
     -     (8*Alfa2*MW*(AbbSumP52*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP48*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP51*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP43*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP12)*
     -     ((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP52 - Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP5)*(-(8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP459 - Abb1*(S + 2*T14 + 2*T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP463 - Abb2*(S + 2*T14 + 2*T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP465 - Abb3*(S + 2*T14 + 2*T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP466 - Abb4*(S + 2*T14 + 2*T24)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*
     -     (mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP6)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP193 + Abb2*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP195 + Abb3*(-S34 - T - U)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP345 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP346 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee14,iintP7)*
     -     (8*Alfa2*MW*(AbbSumP52*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP48*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP51*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP43*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) - 8*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP8)*((8*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP52 - Abb4*S))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee14,iintP9)*((8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP52 - Abb4*S))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP10)*((4*Alfa2*MW)/CW2*
     -     (Abb4*dup11*EL*gL(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup13*EL*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup11*EL*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup13*EL*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP302 + Abb1*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP304 + Abb2*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP306 + Abb3*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP307 + Abb4*(-S34 - T - U)))/
     -     SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP11)*(-(4*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP219 + Abb1*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP221 + Abb2*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP224 + Abb3*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP226 + Abb4*(-T - U)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP12)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP19*dup1*EL*gL(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*dup1*EL*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP219 + Abb1*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP221 + Abb2*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP224 + Abb3*(-T - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP226 + Abb4*(-T - U)))/SW)


      dup19 = -2*S + 2*S34 + T - 2*T14 + 2*U
      dup20 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup21 = 2*mass2(1) + 2*mass2(3) + MH2
      dup22 = 4*mass2(1) + 3*mass2(3) + 2*MH2
      dup23 = gR(1)*mass2(1) + gL(1)*(-mass2(3) - MH2)
      dup24 = gR(1)**2*mass2(1) + gL(1)**2*(-2*mass2(3) - MH2)
      dup25 = -(gL(1)*mass2(1)) + gR(1)*(mass2(3) + MH2)
      dup26 = -(gL(1)**2*mass2(1)) + gR(1)**2*(2*mass2(3) + MH2)

      tmp2 = -(12*Alfa2*MW)/CW2*
     -     (AbbSumP35*EL*gL(1)**2*gL(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP12*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP26*EL*gL(1)**2*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/
     -     SW - (12*Alfa2*MW)/CW2*
     -     (AbbSumP5*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-mass2(3) + MH2))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*(S34 + T + U))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(S34 + T + U))/SW

      tmp2 = tmp2 + 1/CW2*
     -     (Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*(4*S34 + 4*T + 4*U)) + 
     -     Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*(4*S34 + 4*T + 4*U)))/
     -     SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP5)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*
     -     (AbbSumP469 + Abb1*(-(1/3.D0*S) - 1/3.D0*T14)))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP471 + Abb2*(-(1/3.D0*S) - 1/3.D0*T14)))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP473 + Abb3*(-(1/3.D0*S) - 1/3.D0*T14)))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*
     -     (AbbSumP474 + Abb4*(-(1/3.D0*S) - 1/3.D0*T14)))/SW)

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP6)*((16*Alfa2*MW)/CW2*
     -     (AbbSumP365*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP352*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP360*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP350*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*(gL(1) + gR(1))*gR(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP364*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP355*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP358*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP351*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP28*EL*gL(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gR(1)*(gL(3) + gR(3))*mass2(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP349 + Abb1*(-2*S - S34)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP356 + Abb2*(-2*S - S34)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP361 + Abb3*(-2*S - S34)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP363 + Abb4*(-2*S - S34)))/SW)

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee2,iintP7)*
     -     (-(4*Alfa2*MW*(AbbSumP96*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP90*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP94*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP85*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP88*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP86*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP8)*((8*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP88*EL*gL(3)**2*(-gL(1)**2 + gR(1)**2)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP86*EL*(-gL(1)**2 + gR(1)**2)*gR(3)**2*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP85 + Abb1*S))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP90 + Abb2*S))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP94 + Abb3*S))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP96 + Abb4*S))/SW)

      Cloop(1) = Cloop(1) + Eval(ee2,iintP9)*((4*Alfa2*MW)/CW2*
     -     (Abb4*dup23*EL*gL(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup25*EL*gL(3)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup23*EL*gR(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup25*EL*gR(3)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(3)*(gL(1) + gR(1))*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumP11*EL*(gL(1) + gR(1))*gR(3)*
     -     mass2(1)*q(1)*q(3))/SW - (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP302 + Abb1*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP304 + Abb2*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP306 + Abb3*(-S34 - T - U)))/
     -     SW - (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP307 + Abb4*(-S34 - T - U)))/
     -     SW) + 1/CW2*((Eval(ee11,iintP6) + Eval(ee22,iintP11))*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP12) + Eval(ee22,iintP5))*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*((Eval(ee11,iintP11) + Eval(ee22,iintP6))*
     -     (8*Alfa2*MW*(AbbSumP52*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP48*EL*gL(3)*gR(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP51*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP43*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (8*MW*(mass2(1)*q(1)*q(3)) + 8*MW*(mass2(3)*q(1)*q(3)) - 
     -     8*MW*(q(1)*q(3)*S)))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee22,iintP10) + Eval(ee22,iintP8))*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee23,iintP10)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*((Eval(ee22,iintP7) + Eval(ee22,iintP9))*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee23,iintP12) + Eval(ee23,iintP5))*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee23,iintP11)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee23,iintP6)*(-((8*Alfa2*MW)/CW2*
     -     ((AbbSumP371 - Abb4*dup22)*EL*gL(1)*gL(3)*q(1)*q(3))/SW) - 
     -     (8*Alfa2*MW)/CW2*
     -     ((Abb70 - Abb2*dup22)*EL*gL(3)*gR(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((Abb72 - Abb3*dup22)*EL*gL(1)*gR(3)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     ((AbbSumP368 - Abb1*dup22)*EL*gR(1)*gR(3)*q(1)*q(3))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP75*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP68*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP72*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (24*Alfa2*MW)/CW2*
     -     (AbbSumP67*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(MH2 + 3/2.D0*mass2(3)))/SW) + 
     -     1/CW2*(Eval(ee23,iintP7)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee23,iintP9)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     16*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee23,iintP8)*
     -     (-(8*Alfa2*MW*(Abb74*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)**2*gR(1)**2) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb64*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gL(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)**2*gR(1)**2)) + 
     -     8*MW*(dup22*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP10)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(1) + mass2(3) + MH2))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP312 + Abb2*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP313 + Abb3*(-S34 - T - U)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP311 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP314 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee24,iintP12)*
     -     (-(8*Alfa2*MW*(AbbSumP237*EL*gL(1)**2*gL(3)**2)) - 
     -     8*Alfa2*MW*(AbbSumP236*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gL(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)**2*gR(1)**2)) + 
     -     8*MW*(dup22*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gR(1)**2*gR(3)**2)))))/SW + 
     -     1/CW2*(Eval(ee24,iintP11)*
     -     (-(8*Alfa2*MW*(AbbSumP237*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumP236*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP5)*((8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(3*mass2(1) + mass2(3) + MH2))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP476 + Abb1*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP400 + Abb2*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP407 + Abb3*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP477 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + Eval(ee24,iintP6)*((8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP43 - Abb1*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP48 - Abb2*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP51 - Abb3*S))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP52 - Abb4*S))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (16*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)
     -     + 1/CW2*(Eval(ee24,iintP7)*
     -     (-(8*Alfa2*MW*(AbbSumP110*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP108*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP109*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(AbbSumP107*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee24,iintP8)*
     -     (-(8*Alfa2*MW*(AbbSumP110*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP108*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP109*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(AbbSumP107*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gL(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)**2*gR(1)**2)) + 
     -     8*MW*(dup22*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP9)*(-(8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(1) + mass2(3) + MH2)*q(1)*q(3))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*(AbbSumP312 + Abb2*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*(AbbSumP313 + Abb3*(-S34 - T - U)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP311 + Abb1*(S34 + T + U)))/SW
     -     + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP314 + Abb4*(S34 + T + U)))/SW
     -     )

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP10)*(Alfa2*(AbbSumP19*
     -     EL*(gL(1)*(-(4*MW*(gL(1)*gL(3)**2*mass2(1))) - 
     -     4*MW*(gL(3)**2*gR(1)*mass2(1))) - 
     -     4*MW*(gL(3)**2*gR(1)**2*mass2(1)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (gL(1)*(-(4*MW*(gL(1)*gR(3)**2*mass2(1))) - 
     -     4*MW*(gR(1)*gR(3)**2*mass2(1))) - 
     -     4*MW*(gR(1)**2*gR(3)**2*mass2(1)))) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (-(4*MW*mass2(3)) - 4*MW*MH2 + 
     -     MW*(4*S34 + 4*T + 4*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP11)*(-(4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(-2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*U)/SW)

      Cloop(1) = Cloop(1) + Eval(ee3,iintP12)*(-(4*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(-3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(-3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(-3*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(-3*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*U)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*U)/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP5)*
     -     (-(4*Alfa2*MW*(AbbSumP408*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(AbbSumP400*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(AbbSumP407*EL*gL(1)**2*gR(3)**2) - 
     -     4*Alfa2*MW*(AbbSumP395*EL*gR(1)**2*gR(3)**2) + 
     -     28*Alfa2*MW*(AbbSumP484*EL*gL(1)**2*gL(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP480*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(AbbSumP483*EL*gL(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     28*Alfa2*MW*(AbbSumP478*EL*gR(1)**2*gR(3)**2*mass2(1)) - 
     -     4*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*T) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*T) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*T) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*T)))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP7)*((4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW) + 
     -     1/CW2*(Eval(ee3,iintP6)*
     -     (20*Alfa2*MW*(AbbSumP177*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP156*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP169*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     20*Alfa2*MW*(AbbSumP152*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) + 4*MW*(mass2(3)*q(1)*q(3)) + 
     -     4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S - 4*T - 4*T14))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee3,iintP8)*
     -     (4*Alfa2*MW*(AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1)) + 
     -     4*Alfa2*MW*(AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee33,iintP10)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW + 
     -     1/CW2*(Eval(ee33,iintP11)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee3,iintP9)*(Alfa2*(AbbSumP19*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(mass2(1)*q(1)*q(3))) - 
     -     4*MW*(mass2(3)*q(1)*q(3)) - 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*((Eval(ee33,iintP12) + Eval(ee33,iintP5))*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(1)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(1))))/SW +
     -     1/CW2*(Eval(ee33,iintP9)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP10)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)**2*gR(1)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(mass2(3) + MH2))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP319 + Abb1*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP320 + Abb4*(-S34 - T - U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP205 + Abb2*(2*S34 + 2*T + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP209 + Abb3*(2*S34 + 2*T + 2*U)))/SW) + 
     -     1/CW2*(Eval(ee34,iintP11)*
     -     (8*Alfa2*MW*(AbbSumP245*EL*gL(1)*gL(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP244*EL*gR(1)*gR(3)*q(1)*q(3))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP5)*((8*Alfa2*MW)/CW2*
     -     (Abb4*dup24*EL*gL(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb2*dup26*EL*gL(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*dup24*EL*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*dup26*EL*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*
     -     (AbbSumP400 + Abb2*(-S - T14 - T24)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*
     -     (AbbSumP407 + Abb3*(-S - T14 - T24)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP487 + Abb1*(S + T14 + T24)))/
     -     SW + (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP488 + Abb4*(S + T14 + T24)))/
     -     SW) + 1/CW2*(Eval(ee34,iintP12)*
     -     (8*Alfa2*MW*(AbbSumP245*EL*gL(1)**2*gL(3)**2) + 
     -     8*Alfa2*MW*(AbbSumP244*EL*gR(1)**2*gR(3)**2)))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP7)*
     -     (-(8*Alfa2*MW*(Abb81*EL*gL(1)*gL(3)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(Abb50*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(Abb55*EL*gL(1)*gR(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb79*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)*gL(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gL(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)*gR(1)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(3)*gR(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gL(1)*gR(3)*q(1)*q(3)))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)*gR(3)*q(1)*q(3))) + 
     -     8*MW*(dup22*EL*gR(1)*gR(3)*q(1)*q(3))))))/SW + 
     -     1/CW2*(Eval(ee34,iintP6)*
     -     (8*Alfa2*MW*(AbbSumP386*EL*gL(1)*gL(3)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb70*EL*gL(3)*gR(1)*q(1)*q(3)) - 
     -     8*Alfa2*MW*(Abb72*EL*gL(1)*gR(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumP385*EL*gR(1)*gR(3)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (8*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     8*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (8*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     8*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3))))))/SW

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee34,iintP8)*
     -     (-(8*Alfa2*MW*(Abb81*EL*gL(1)**2*gL(3)**2)) - 
     -     4*Alfa2*MW*(Abb50*EL*gL(3)**2*gR(1)**2) - 
     -     4*Alfa2*MW*(Abb55*EL*gL(1)**2*gR(3)**2) - 
     -     8*Alfa2*MW*(Abb79*EL*gR(1)**2*gR(3)**2) + 
     -     Alfa2*(Abb4*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gL(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gL(3)**2))) + 
     -     Alfa2*(Abb2*
     -     (-(16*MW*(dup20*EL*gL(3)**2*gR(1)**2)) + 
     -     8*MW*(dup22*EL*gL(3)**2*gR(1)**2))) + 
     -     Alfa2*(Abb3*
     -     (-(16*MW*(dup20*EL*gL(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gL(1)**2*gR(3)**2))) + 
     -     Alfa2*(Abb1*
     -     (-(16*MW*(dup20*EL*gR(1)**2*gR(3)**2)) + 
     -     8*MW*(dup22*EL*gR(1)**2*gR(3)**2)))))/SW

      Cloop(1) = Cloop(1) + tmp2*Eval(ee4,iintP10) + 
     -     Eval(ee34,iintP9)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumP211*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP202*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP208*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP201*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*(mass2(3) + MH2)*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)*q(1)*q(3)*(AbbSumP319 + Abb1*(-S34 - T - U)))/
     -     SW - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)*q(1)*q(3)*(AbbSumP320 + Abb4*(-S34 - T - U)))/
     -     SW + (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)*q(1)*q(3)*
     -     (AbbSumP205 + Abb2*(2*S34 + 2*T + 2*U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)*q(1)*q(3)*
     -     (AbbSumP209 + Abb3*(2*S34 + 2*T + 2*U)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP11)*
     -     (Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(1)*q(1)*q(3)) - 4*MW*(q(1)*q(3)*U)))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP12)*((4*Alfa2*MW)/CW2*
     -     (AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)**2*gL(3)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**2*gR(1)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)**2*gR(3)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)**2*gR(3)**2*(-mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**2*gL(3)**2*U)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**2*gR(1)**2*U)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**2*gR(3)**2*U)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**2*gR(3)**2*U)/SW)

      Cloop(1) = Cloop(1) + Eval(ee4,iintP5)*(-((4*Alfa2*MW)/CW2*
     -     ((AbbSumP496 + Abb2*dup19)*EL*gL(3)**2*gR(1)**2)/SW) - 
     -     (4*Alfa2*MW)/CW2*
     -     ((AbbSumP499 + Abb3*dup19)*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP409*EL*gL(1)**2*gL(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP19*EL*gL(1)*gL(3)**2*gR(1)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP398*EL*gL(3)**2*gR(1)**2*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumP406*EL*gL(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumP11*EL*gL(1)*gR(1)*gR(3)**2*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumP396*EL*gR(1)**2*gR(3)**2*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP491 + 1/2.D0*(Abb1*dup19)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP501 + 1/2.D0*(Abb4*dup19)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP6)*
     -     (12*Alfa2*MW*(AbbSumP35*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (-(4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (-(4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3))) - 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (4*MW*(mass2(3)*q(1)*q(3)) + 4*MW*(MH2*q(1)*q(3)) + 
     -     MW*(q(1)*q(3)*(-4*S34 - 4*T24 - 4*U))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP7)*(-(4*Alfa2*MW)/CW2*(Abb4*
     -     EL*gL(1)*gL(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*S)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*S)/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)*gR(1)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)*q(1)*q(3)*(mass2(3) + 1/2.D0*mass2(1)))/SW)

      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee4,iintP8)*
     -     (4*Alfa2*MW*(AbbSumP28*EL*gL(1)**2*gL(3)*gR(3)*mass2(3)) + 
     -     4*Alfa2*MW*(AbbSumP9*EL*gL(3)*gR(1)**2*gR(3)*mass2(3)) + 
     -     Alfa2*(Abb4*EL*gL(1)**2*gL(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb2*EL*gL(3)**2*gR(1)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb3*EL*gL(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S)) + 
     -     Alfa2*(Abb1*EL*gR(1)**2*gR(3)**2*
     -     (4*MW*mass2(1) + 4*MW*mass2(3) - 4*MW*S))))/SW

      Cloop(1) = Cloop(1) + 
     -     1/CW2*(Eval(ee44,iintP10)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**2*gL(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**2*gR(1)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**2*gR(3)**2*mass2(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**2*gR(3)**2*mass2(3))))/SW + 
     -     1/CW2*(Eval(ee4,iintP9)*
     -     (-(12*Alfa2*MW*
     -     (AbbSumP35*EL*gL(1)*gL(3)*mass2(1)*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumP12*EL*gL(3)*gR(1)*mass2(1)*q(1)*q(3)) - 
     -     4*Alfa2*MW*(AbbSumP26*EL*gL(1)*gR(3)*mass2(1)*q(1)*q(3)) - 
     -     12*Alfa2*MW*(AbbSumP5*EL*gR(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     Alfa2*(AbbSumP19*EL*
     -     (4*MW*(gL(1)*gL(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gL(3)*gR(1)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(AbbSumP11*EL*
     -     (4*MW*(gL(1)*gR(3)*mass2(1)*q(1)*q(3)) + 
     -     4*MW*(gR(1)*gR(3)*mass2(1)*q(1)*q(3)))) + 
     -     Alfa2*(Abb4*EL*gL(1)*gL(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb2*EL*gL(3)*gR(1)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb3*EL*gL(1)*gR(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U)))) + 
     -     Alfa2*(Abb1*EL*gR(1)*gR(3)*
     -     (-(4*MW*(MH2*q(1)*q(3))) + 
     -     MW*(q(1)*q(3)*(4*S34 + 4*T + 4*U))))))/SW

      Cloop(1) = Cloop(1) + Eval(ee44,iintP5)*((8*Alfa2*MW)/CW2*
     -     (Abb4*dup21*EL*gL(1)**2*gL(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*dup21*EL*gL(3)**2*gR(1)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*dup21*EL*gL(1)**2*gR(3)**2)/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*dup21*EL*gR(1)**2*gR(3)**2)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*gR(3)**2*(AbbSumP505 + Abb1*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*gR(1)**2*(AbbSumP205 + Abb2*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gR(3)**2*(AbbSumP209 + Abb3*(S34 + T + U)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*gL(3)**2*(AbbSumP510 + Abb4*(S34 + T + U)))/SW) 


      Cloop(1) = Cloop(1) + 1/CW2*(Eval(ee44,iintP6)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW + 
     -     1/CW2*(Eval(ee44,iintP9)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)*mass2(3)*q(1)*q(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)*mass2(3)*q(1)*q(3))))/SW


      end


********************************************************************************
********************************************************************************
***                                                                          ***
********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 1 and 4 for:
**     quark + quark -> quark + H + quark 

      subroutine pentP2_qq14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      Cloop(1) = Cloop(1) + Eval(ee12,iintP16)*
     -     ((2*Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP172 + Abb4*(-S + T)))/SW
     -     ) + Eval(ee12,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP431 + Abb4*(S + T)))/SW)
     -     + Eval(ee11,iintP16)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP55 + Abb4*(S + T24)))/SW
     -     ) + Eval(ee1,iintP15)*
     -     (-(2*Alfa2*MW)/SW2**2*(Abb4*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(AbbSumP32*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP416 + Abb4*(-S + S34 - T - T14 + T24)))/SW)
     -     + Eval(ee1,iintP16)*
     -     (-(5*Alfa2*MW)/SW2**2*(AbbSumP155*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(AbbSumP147*EL*mass2(3))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP153 + Abb4*(-2*S34 - 2*T - 2*U)))/SW) + 
     -     Eval(ee0,iintP16)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-mass2(1) + 2*mass2(2) + mass2(3)))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) - (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Dval(dd0,iintP13) + Dval(dd0,iintP14)))/SW + 
     -     (6*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP15) + Eval(ee00,iintP16)))/SW + 
     -     Eval(ee0,iintP15)*((Alfa2*MW)/SW2**2*(AbbSumP401*EL)/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(2) - 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee11,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP423*EL)) + Alfa2*(Abb4*EL*
     -     (2*MW*mass2(3) + 2*MW*MH2 - 2*MW*T24))))/SW

      Cloop(1) = Cloop(1) + Eval(ee14,iintP16)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP55 + Abb4*(S + T24)))/SW
     -     ) + Eval(ee2,iintP16)*
     -     (-(3*Alfa2*MW)/SW2**2*(AbbSumP37*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP193 + Abb4*(2*T + U)))/SW)
     -     + 1/SW2**2*(Eval(ee13,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP435*EL))))/SW + Eval(ee2,iintP15)*
     -     (-(6*Alfa2*MW)/SW2**2*(AbbSumP453*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (3*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP452 + Abb4*(T + 1/3.D0*S + 1/3.D0*T14)))/SW) + 
     -     1/SW2**2*(Eval(ee24,iintP16)*
     -     (-(2*Alfa2*MW*(AbbSumP205*EL)) + Alfa2*(Abb4*EL*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2**2*((Eval(ee22,iintP15) + Eval(ee22,iintP16))*
     -     (2*Alfa2*MW*(Abb64*EL) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee24,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP456*EL))+ 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*((Eval(ee23,iintP15) + Eval(ee23,iintP16))*
     -     (2*Alfa2*MW*(Abb64*EL) + 
     -     Alfa2*(Abb4*EL*(-(4*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP16)*
     -     (2*Alfa2*MW*(AbbSumP178*EL) + 
     -     Alfa2*(Abb4*EL*(-(4*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee14,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP441*EL)) + 
     -     Alfa2*(Abb4*EL*(2*MW*MH2 - 2*MW*T24))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP16)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-mass2(1) + mass2(3)))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) + 
     -     Eval(ee4,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP467 + Abb4*(-S - T - T14)))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(3) + MH2 + 5/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee34,iintP16)*
     -     (2*Alfa2*MW*(AbbSumP211*EL)))/SW + 
     -     1/SW2**2*(Eval(ee3,iintP15)*
     -     (Alfa2*MW*(AbbSumP401*EL) + 
     -     Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP15)*
     -     (2*Alfa2*MW*(AbbSumP461*EL)))/SW + 
     -     Eval(ee44,iintP15)*(-(2*Alfa2*MW)/SW2**2*(EL*(AbbSumP471 + 
     -     Abb4*(4*mass2(1) + 3*mass2(3) + 2*MH2)))/SW + 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + MH2 + 3/2.D0*mass2(3)))/SW) + 
     -     1/SW2**2*(Eval(ee44,iintP16)*
     -     (2*Alfa2*MW*(Abb4*EL*mass2(3))))/SW + 
     -     1/SW2**2*(Eval(ee4,iintP16)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*U))))
     -     /SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 2 and 3 for:
**     quark + quark -> quark + H + quark 

      subroutine pentP2_qq23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      Cloop(1) = Cloop(1) + Eval(ee13,iintP13)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP75 + Abb4*S))/SW) + 
     -     (Eval(ee12,iintP13) + Eval(ee12,iintP14))*
     -     ((2*Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP65 + Abb4*(-S + T)))/SW)
     -     + (Eval(ee11,iintP13) + Eval(ee11,iintP14))*
     -     ((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP55 + Abb4*(S + T24)))/SW
     -     ) + Eval(ee1,iintP13)*
     -     ((3*Alfa2*MW)/SW2**2*(AbbSumP37*EL*mass2(1))/SW + 
     -     (4*Alfa2*MW)/SW2**2*(AbbSumP36*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW - 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP35 + Abb4*(2*S + T24)))/SW
     -     ) + Eval(ee1,iintP14)*
     -     ((2*Alfa2*MW)/SW2**2*(AbbSumP257*EL*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(AbbSumP254*EL*mass(3))/SW + 
     -     (3*Alfa2*MW)/SW2**2*(AbbSumP256*EL*mass2(3))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP255 + Abb4*(2*S + T14 + 2*T24)))/SW) + 
     -     (8*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP13) + Eval(ee00,iintP14)))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP13)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*(Eval(ee0,iintP14)*
     -     (-3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP13)*((Alfa2*MW)/SW2**2*(AbbSumP23*EL*mass2(1))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP94 + Abb4*(-S + T)))/SW)+
     -     Eval(ee14,iintP14)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP55 + Abb4*(S + T24)))/SW
     -     ) + Eval(ee14,iintP13)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 3*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP55 + Abb4*(S + T24)))/SW
     -     ) + Eval(ee2,iintP14)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP285*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP284 + Abb4*(S34 + 2*T + U)))/SW) + 
     -     1/SW2**2*(Eval(ee23,iintP13)*
     -     (2*Alfa2*MW*(Abb64*EL) + 2*Alfa2*MW*(Abb4*EL*T)
     -     ))/SW + 1/SW2**2*(Eval(ee24,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP106*EL)) + Alfa2*(Abb4*EL*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2**2*((Eval(ee22,iintP13) + Eval(ee22,iintP14))*
     -     (2*Alfa2*MW*(Abb64*EL) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP14)*(2*Alfa2*MW*(Abb64*EL) + 
     -     Alfa2*(Abb4*EL*(-(4*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP265*EL) + 
     -     Alfa2*(Abb4*EL*(-(4*MW*mass2(1)) + 2*MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP13)*((Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3)))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW) + (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*mass2(1)*Eval(ee33,iintP13))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*mass2(3)*Eval(ee44,iintP13))/SW + 
     -     Eval(ee4,iintP13)*(-(Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW + 
     -     (3*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + 1/3.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee34,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP112*EL)) + Alfa2*(Abb4*EL*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2**2*(Eval(ee24,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP294*EL) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP14)*
     -     (-(2*Alfa2*MW*(AbbSumP304*EL)) + 
     -     Alfa2*(Abb4*EL*(-(4*MW*mass2(1)) + 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee3,iintP14)*
     -     (-3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/
     -     SW + 1/SW2**2*(Eval(ee4,iintP14)*
     -     (-3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/
     -     SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 1 and 4 for:
**     antiquark + quark -> antiquark + H + quark 

      subroutine pentP2_qbarq14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** dummy parameters
      double complex dup13, dup14

** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup13 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup14 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + Eval(ee0,iintP15)*(-(5*Alfa2*MW)/SW2**2*
     -     (AbbSumP485*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-2*mass2(2) + mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(S34 + T + U))/SW) + 
     -     (Eval(ee12,iintP14) + Eval(ee12,iintP15))*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP65 + Abb4*U))/SW) + 
     -     (Eval(ee11,iintP14) + Eval(ee11,iintP15))*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 2*mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP54 + Abb4*(T24 + U)))/SW
     -     ) + Eval(ee1,iintP15)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP253*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(AbbSumP252*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP250 + Abb4*(T14 + T24 + 2*U)))/SW) + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Dval(dd0,iintP12) + Dval(dd0,iintP13)))/SW - 
     -     (6*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP14) + Eval(ee00,iintP15)))/SW + 
     -     Eval(ee0,iintP14)*(-(Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(MW2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)))/SW) + 
     -     1/SW2**2*(Eval(ee1,iintP14)*
     -     (Alfa2*MW*(AbbSumP36*EL) + 
     -     3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + 
     -     Alfa2*MW*(AbbSumP26*EL*mass2(3)) - 
     -     2*Alfa2*MW*(Abb4*EL*S)))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP71*EL) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup13*EL)) + 2*MW*(dup14*EL))))
     -     )/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP14)*((Alfa2*MW)/SW2**2*
     -     (AbbSumP23*EL*mass2(1))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP89 + Abb4*(-S + T)))/SW)+
     -     Eval(ee13,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP262 - Abb4*U))/SW) + 
     -     Eval(ee2,iintP15)*(-(4*Alfa2*MW)/SW2**2*(AbbSumP253*EL*
     -     mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP279 + Abb4*(S34 + U)))/SW)
     -     + Eval(ee14,iintP15)*
     -     ((4*Alfa2*MW)/SW2**2*(AbbSumP272*EL*mass(1)*mass(3))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP54 + Abb4*(T24 + U)))/SW
     -     ) + Eval(ee14,iintP14)*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 3*mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP54 + Abb4*(T24 + U)))/SW
     -     ) + 1/SW2**2*((Eval(ee22,iintP14) + Eval(ee22,iintP15))*
     -     (2*Alfa2*MW*(Abb55*EL) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP15)*
     -     (2*Alfa2*MW*(Abb55*EL) - 
     -     4*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP14)*
     -     (2*Alfa2*MW*(Abb55*EL) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup13*EL)) + 2*MW*(dup14*EL))))
     -     )/SW + 1/SW2**2*
     -     (Eval(ee24,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP102*EL) + 
     -     Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP14)*((Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3)))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW) + 
     -     1/SW2**2*(Eval(ee33,iintP15)*
     -     (-2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP15)*
     -     (2*Alfa2*MW*(AbbSumP299*EL) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee24,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP286*EL)) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee44,iintP14)*
     -     (-2*Alfa2*MW*(Abb4*EL*mass2(3))))/SW + 
     -     1/SW2**2*(Eval(ee33,iintP14)*
     -     (Alfa2*(Abb4*(4*MW*(dup13*EL) - 2*MW*(dup14*EL)))))/
     -     SW + 1/SW2**2*(Eval(ee34,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP108*EL) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup13*EL)) + 2*MW*(dup14*EL))))
     -     )/SW + 1/SW2**2*
     -     (Eval(ee4,iintP14)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*
     -     (Eval(ee4,iintP15)*
     -     (-3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*(-(MW*MH2) + MW*(S + T14 + T24)))))/
     -     SW + 1/SW2**2*(Eval(ee3,iintP15)*
     -     (-7*Alfa2*MW*(AbbSumP295*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/
     -     SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 2 and 3 for:
**     antiquark + quark -> antiquark + H + quark 

      subroutine pentP2_qbarq23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** dummy parameters
      double complex dup13


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      dup13 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + Eval(ee1,iintP12)*((Alfa2*MW)/SW2**2*
     -     (AbbSumP32*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP412 + Abb4*(T24 - U)))/SW)
     -     + Eval(ee12,iintP13)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP174 + Abb4*U))/SW) + 
     -     Eval(ee11,iintP13)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 2*mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP54 + Abb4*(T24 + U)))/SW
     -     ) + Eval(ee1,iintP13)*
     -     (-(3*Alfa2*MW)/SW2**2*(AbbSumP37*EL*mass2(1))/SW - 
     -     (4*Alfa2*MW)/SW2**2*(AbbSumP157*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP155 + Abb4*(T24 + 2*U)))/SW) - 
     -     (8*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP12) + Eval(ee00,iintP13)))/SW + 
     -     Eval(ee0,iintP12)*((Alfa2*MW)/SW2**2*
     -     ((AbbSumP397 + Abb4*dup13)*EL)/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(MH2 + 5/2.D0*mass2(1) + 3/2.D0*mass2(3)))/SW) + 
     -     1/SW2**2*(Eval(ee11,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP420*EL) + 
     -     Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(3)) - 2*MW*MH2 + 2*MW*T24))))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP13)*
     -     (Alfa2*(Abb4*EL*
     -     (-(MW*mass2(1)) - MW*mass2(3) + MW*U))))/SW + 
     -     1/SW2**2*(Eval(ee12,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP427*EL) + 
     -     Alfa2*(Abb4*EL*(2*MW*mass2(3) - 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP12)*(-(6*Alfa2*MW)/SW2**2*
     -     (AbbSumP453*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(3*AbbSumP452 + Abb4*(S + T + T14)))/SW) + 
     -     Eval(ee13,iintP13)*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3)))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP182 - Abb4*U))/SW) + 
     -     Eval(ee14,iintP13)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 3*mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(EL*(AbbSumP54 + Abb4*(T24 + U)))/SW
     -     ) + 1/SW2**2*((Eval(ee22,iintP12) + Eval(ee22,iintP13))*
     -     (2*Alfa2*MW*(Abb55*EL) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*((Eval(ee23,iintP12) + Eval(ee23,iintP13))*
     -     (2*Alfa2*MW*(Abb55*EL) - 
     -     4*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee24,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP458*EL) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee14,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP441*EL) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*MH2) + 2*MW*T24))))/SW + 
     -     1/SW2**2*(Eval(ee2,iintP13)*
     -     (Alfa2*MW*(AbbSumP196*EL) - 
     -     3*Alfa2*MW*(AbbSumP37*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*(-(MW*mass2(3)) + MW*U))))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP434*EL) + 
     -     Alfa2*(Abb4*EL*(2*MW*mass2(3) - 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee3,iintP13)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3)))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*mass2(1)*
     -     (Eval(ee33,iintP12) + Eval(ee33,iintP13)))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*mass2(3)*Eval(ee44,iintP13))/SW + 
     -     Eval(ee4,iintP13)*((Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW - 
     -     (3*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + 1/3.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee3,iintP12)*
     -     (Alfa2*MW*(AbbSumP397*EL) - 
     -     3*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP12)*
     -     (-(2*Alfa2*MW*(AbbSumP466*EL)) - 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee4,iintP12)*
     -     (2*Alfa2*MW*(AbbSumP475*EL) - 
     -     Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     Eval(ee44,iintP12)*
     -     ((2*Alfa2*MW)/SW2**2*((AbbSumP483 + Abb4*dup13)*EL)/SW - 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + MH2 + 3/2.D0*mass2(3)))/SW) + 
     -     1/SW2**2*(Eval(ee34,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP213*EL)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) + 2*MW*U))))/SW
     -     + 1/SW2**2*(Eval(ee24,iintP13)*
     -     (2*Alfa2*MW*(AbbSumP207*EL) + 
     -     Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) + 2*MW*U))))/SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 1 and 4 for:
**     quark + antiquark -> quark + H + antiquark 

      subroutine pentP2_qqbar14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** dummy parameters
      double complex dup13, dup14


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      Cloop(1) = Cloop(1) + Eval(ee13,iintP14)*
     -     ((4*Alfa2*MW)/SW2**2*(AbbSumP85*EL*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP83 + Abb4*(-S - T - T14)))/SW) + 
     -     (Eval(ee12,iintP14) + Eval(ee12,iintP15))*
     -     ((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP69 + Abb4*(S + 2*T + T14)))/SW) + 
     -     Eval(ee1,iintP14)*((3*Alfa2*MW)/SW2**2*(AbbSumP40*EL*
     -     mass2(1))/SW + 
     -     (4*Alfa2*MW)/SW2**2*(AbbSumP41*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW - 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP39 + Abb4*(2*S + T24)))/SW
     -     ) + Eval(ee1,iintP15)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP281*EL*mass2(1))/SW - 
     -     (3*Alfa2*MW)/SW2**2*(AbbSumP282*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP280 + Abb4*(-T14 - 2*U)))/SW) + 
     -     Eval(ee0,iintP15)*(-(5*Alfa2*MW)/SW2**2*(AbbSumP317*EL* 
     -     mass2(1))/SW - (Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-2*mass2(2) + mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(S34 + T + U))/SW) + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Dval(dd0,iintP12) + Dval(dd0,iintP13)))/SW - 
     -     (6*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP14) + Eval(ee00,iintP15)))/SW + 
     -     Eval(ee0,iintP14)*(-(Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(MW2 + 1/2.D0*mass2(1) + 1/2.D0*mass2(3)))/SW) + 
     -     1/SW2**2*((Eval(ee11,iintP14) + Eval(ee11,iintP15))*
     -     (-(2*Alfa2*MW*(AbbSumP58*EL)) + Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee24,iintP15)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP310 + Abb4*(S + 2*T + T14)))/SW) + 
     -     Eval(ee2,iintP15)*((Alfa2*MW)/SW2**2*(Abb2*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP303 + Abb4*(-S34 + T - U)))/SW) + 
     -     Eval(ee14,iintP15)*
     -     (-(2*Alfa2*MW)/SW2**2*(EL*(AbbSumP58 - Abb4*U))/SW - 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee24,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP120*EL)))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP14)*
     -     (-(2*Alfa2*MW*(Abb63*EL)) - 
     -     2*Alfa2*MW*(Abb4*EL*T)))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP289*EL)) + 
     -     2*Alfa2*MW*(AbbSumP98*EL*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*T)))/SW + 
     -     1/SW2**2*(Eval(ee2,iintP14)*
     -     (-(Alfa2*MW*(AbbSumP108*EL)) + 
     -     Alfa2*MW*(AbbSumP98*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*(MW*mass2(3) - MW*S))))/SW + 
     -     1/SW2**2*((Eval(ee22,iintP14) + Eval(ee22,iintP15))*
     -     (Alfa2*(Abb4*EL*(2*MW*mass2(1) - 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP15)*
     -     (Alfa2*(Abb4*EL*(4*MW*mass2(1) - 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee14,iintP14)*
     -     (-(2*Alfa2*MW*(AbbSumP58*EL)) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee4,iintP15)*(-(3*Alfa2*MW)/SW2**2*
     -     (AbbSumP40*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(2*mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(S + T14 + T24))/SW) + 
     -     Eval(ee4,iintP14)*(-(Alfa2*MW)/SW2**2*(Abb4*EL*S)/SW + 
     -     (3*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + 1/3.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee33,iintP14)*
     -     (-2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP14)*
     -     (2*Alfa2*MW*(AbbSumP128*EL)))/SW + 
     -     1/SW2**2*(Eval(ee44,iintP15)*
     -     (-2*Alfa2*MW*(Abb4*EL*mass2(3))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP15)*
     -     (-(2*Alfa2*MW*(AbbSumP320*EL)) + 
     -     2*Alfa2*MW*(AbbSumP98*EL*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*T)))/SW + 
     -     1/SW2**2*(Eval(ee3,iintP14)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*(Eval(ee3,iintP15)*
     -     (-5*Alfa2*MW*(AbbSumP317*EL*mass2(1)) + Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 2 and 3 for:
**     quark + antiquark -> quark + H + antiquark 

      subroutine pentP2_qqbar23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2


      Cloop(1) = Cloop(1) + Eval(ee12,iintP12)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP440 + Abb4*(-S - T14)))/SW) + 
     -     Eval(ee12,iintP13)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP194 + Abb4*(S + 2*T + T14)))/SW) + 
     -     Eval(ee11,iintP12)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP433 + Abb4*(S34 + T24)))/SW) + 
     -     Eval(ee1,iintP12)*((2*Alfa2*MW)/SW2**2*(AbbSumP25*EL*
     -     mass2(1))/SW - (Alfa2*MW)/SW2**2*(AbbSumP33*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP423 + Abb4*(-S + S34 - T14 + T24)))/SW) - 
     -     (8*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP12) + Eval(ee00,iintP13)))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP12)*
     -     (-(Alfa2*MW*(AbbSumP404*EL)) + 
     -     3*Alfa2*MW*(AbbSumP406*EL*mass2(1)) - Alfa2*MW*(Abb4*EL*T)))/
     -     SW + 1/SW2**2*(Eval(ee1,iintP13)*
     -     (-(Alfa2*MW*(AbbSumP176*EL)) - 
     -     3*Alfa2*MW*(AbbSumP40*EL*mass2(1)) - 
     -     Alfa2*MW*(AbbSumP168*EL*mass2(3)) + 
     -     2*Alfa2*MW*(Abb4*EL*U)))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP13)*
     -     (Alfa2*(Abb4*EL*
     -     (-(MW*mass2(1)) - MW*mass2(3) + MW*U))))/SW + 
     -     1/SW2**2*(Eval(ee11,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP58*EL)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(2*MW*mass2(1)) - 2*MW*mass2(3) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee2,iintP12)*((2*Alfa2*MW)/SW2**2*
     -     (AbbSumP85*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(3*AbbSumP462 + Abb4*(-S + 2*T - T14)))/SW) + 
     -     Eval(ee24,iintP12)*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP467 + Abb4*(S + T14)))/SW) + 
     -     Eval(ee2,iintP13)*((2*Alfa2*MW)/SW2**2*
     -     (0*EL*mass(1))/SW + 
     -     (Alfa2*MW)/SW2**2*(AbbSumP25*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP215 + Abb4*(T - U)))/SW)+
     -     Eval(ee14,iintP12)*((2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP454 + Abb4*(-2*S34 - 2*T24 - U)))/SW + 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*(Eval(ee13,iintP12)*
     -     (2*Alfa2*MW*(Abb65*EL) + 
     -     2*Alfa2*MW*(AbbSumP25*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee13,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP202*EL)) + 
     -     2*Alfa2*MW*(AbbSumP98*EL*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*T)))/SW + 
     -     1/SW2**2*((Eval(ee22,iintP12) + Eval(ee22,iintP13))*
     -     (-(2*Alfa2*MW*(Abb63*EL)) + 
     -     Alfa2*(Abb4*EL*(2*MW*mass2(1) - 2*MW*T))))/SW + 
     -     1/SW2**2*((Eval(ee23,iintP12) + Eval(ee23,iintP13))*
     -     (-(2*Alfa2*MW*(Abb63*EL)) + 
     -     Alfa2*(Abb4*EL*(4*MW*mass2(1) - 2*MW*T))))/SW + 
     -     1/SW2**2*(Eval(ee14,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP58*EL)) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(1)) + 2*MW*U))))/SW

      Cloop(1) = Cloop(1) + Eval(ee44,iintP12)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP495 + Abb4*(-S34 - T24 - U)))/SW) + 
     -     Eval(ee4,iintP12)*((3*Alfa2*MW)/SW2**2*(AbbSumP406*EL*
     -     mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(2*AbbSumP486 + Abb4*(S - S34 + T14 - T24 - U)))/SW
     -     ) + 
     -     1/SW2**2*(Eval(ee34,iintP12)*
     -     (-(2*Alfa2*MW*(AbbSumP475*EL)) + 
     -     2*Alfa2*MW*(AbbSumP25*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP13)*
     -     (-(2*Alfa2*MW*(Abb80*EL))))/SW + 
     -     1/SW2**2*(Eval(ee24,iintP13)*
     -     (2*Alfa2*MW*(AbbSumP229*EL)))/SW + 
     -     1/SW2**2*(Eval(ee3,iintP12)*
     -     (-(Alfa2*MW*(AbbSumP404*EL)) + 
     -     3*Alfa2*MW*(AbbSumP406*EL*mass2(1)) - Alfa2*MW*(Abb4*EL*T)))/
     -     SW + 1/SW2**2*(Eval(ee3,iintP13)*
     -     (Alfa2*(Abb4*EL*
     -     (-(MW*mass2(1)) - MW*mass2(3) + MW*U))))/SW + 
     -     1/SW2**2*(Eval(ee4,iintP13)*
     -     (Alfa2*(Abb4*EL*
     -     (-(MW*mass2(1)) - MW*mass2(3) + MW*U))))/SW


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 1 and 4 for:
**     antiquark + antiquark -> antiquark + H + antiquark 

      subroutine pentP2_qbarqbar14(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** dummy parameters
      double complex dup13, dup14


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

*      write(*,*)'new NC qk=4 cur 14 pent2'
      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      Cloop(1) = Cloop(1) + Eval(ee12,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP450 + Abb4*(-S - T14 - T24)))/SW) + 
     -     Eval(ee11,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP444 + Abb4*(T14 + T24)))/SW) + 
     -     Eval(ee1,iintP15)*((2*Alfa2*MW)/SW2**2*(Abb2*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(AbbSumP433*EL*mass2(3))/SW - 
     -     (3*Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP432 + Abb4*(S - S34 + 3*T14 + 2*T24)))/SW)+
     -     Eval(ee0,iintP16)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-mass2(1) + 2*mass2(2) + mass2(3)))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) + 
     -     Eval(ee1,iintP16)*(-(5*Alfa2*MW)/SW2**2*(AbbSumP177*EL*
     -     mass2(1))/SW - 
     -     (4*Alfa2*MW)/SW2**2*(AbbSumP176*EL*mass2(3))/SW - 
     -     (3*Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP174 + Abb4*(2*S34 + 2*T + T24 + 2*U)))/SW)-
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Dval(dd0,iintP13) + Dval(dd0,iintP14)))/SW + 
     -     (6*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP15) + Eval(ee00,iintP16)))/SW + 
     -     1/SW2**2*(Eval(ee11,iintP16)*
     -     (2*Alfa2*MW*(AbbSumP52*EL) + 
     -     Alfa2*(Abb4*EL*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP15)*
     -     (-(Alfa2*MW*(AbbSumP408*EL)) + 
     -     5*Alfa2*MW*(AbbSumP512*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*(-(2*MW*mass2(2)) - MW*T))))/SW

      Cloop(1) = Cloop(1) + Eval(ee13,iintP15)*
     -     ((2*Alfa2*MW)/SW2**2*(Abb2*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(2*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP456 + Abb4*(-S - T14 - T24)))/SW) + 
     -     Eval(ee2,iintP16)*((Alfa2*MW)/SW2**2*(AbbSumP19*EL*mass2(1))/
     -     SW - (Alfa2*MW)/SW2**2*(Abb4*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(EL*(AbbSumP226 + Abb4*(-T - U)))/SW)+
     -     Eval(ee12,iintP16)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP196 + Abb4*(-S34 - T - U)))/SW) + 
     -     Eval(ee13,iintP16)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP211*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP210 + Abb4*(S34 + T + U)))/SW) + 
     -     Eval(ee14,iintP15)*(-(2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP466 - Abb4*(S + 2*T14 + 2*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + MH2 + 3/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*((Eval(ee22,iintP15) + Eval(ee22,iintP16))*
     -     (-(2*Alfa2*MW*(Abb74*EL)) + 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     Eval(ee2,iintP15)*( (2*Alfa2*MW)/SW2**2*(AbbSumP211*EL*
     -     mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW - 
     -     (3*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP474 + Abb4*(-(1/3.D0*S) - 1/3.D0*T14)))/SW) + 
     -     1/SW2**2*(Eval(ee14,iintP16)*
     -     (2*Alfa2*MW*(AbbSumP52*EL) + 
     -     Alfa2*(Abb4*EL*(2*MW*mass2(1) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP15)*
     -     ((2*Alfa2*MW)/SW2**2*(Abb2*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(2*mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP488 + Abb4*(S + T14 + T24)))/SW) + 
     -     Eval(ee4,iintP16)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-mass2(1) + mass2(3)))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) + 
     -     Eval(ee3,iintP16)*(-(Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(-3*mass2(1) + mass2(3)))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*U)/SW) + 
     -     Eval(ee24,iintP15)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(3*mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP477 + Abb4*(S34 + T + U)))/SW) + 
     -     1/SW2**2*((Eval(ee33,iintP15) + Eval(ee33,iintP16))*
     -     (2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*((Eval(ee23,iintP15) + Eval(ee23,iintP16))*
     -     (-(2*Alfa2*MW*(Abb74*EL)) + 
     -     4*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP16)*
     -     (2*Alfa2*MW*(AbbSumP245*EL)))/SW + 
     -     Eval(ee24,iintP16)*(-(2*Alfa2*MW)/SW2**2*(EL*(AbbSumP237 - 
     -     Abb4*(4*mass2(1) + 3*mass2(3) + 2*MH2)))/SW - 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + MH2 + 3/2.D0*mass2(3)))/SW) + 
     -     Eval(ee4,iintP15)*((3*Alfa2*MW)/SW2**2*(AbbSumP409*EL*
     -     mass2(1))/SW - (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP501 + Abb4*(-S + S34 - T14 + U + 1/2.D0*T)))/
     -     SW) + 1/SW2**2*(Eval(ee3,iintP15)*
     -     (-(Alfa2*MW*(AbbSumP408*EL)) + 
     -     7*Alfa2*MW*(AbbSumP484*EL*mass2(1)) - 
     -     Alfa2*MW*(Abb4*EL*T)))/SW

      Cloop(1) = Cloop(1) + Eval(ee44,iintP15)*((2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(2*mass2(1) + 2*mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP510 + Abb4*(S34 + T + U)))/SW) 


      end


********************************************************************************
********************************************************************************

** This subroutine is diagram set P2 for currents 2 and 3 for:
**     antiquark + antiquark -> antiquark + H + antiquark 

      subroutine pentP2_qbarqbar23(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,
     -     Cloop)

      implicit none


** Momentum invariants
      double precision S, T, T14, S34, U, T24

** Loop level component
      double complex Cloop(1)

** Mass squared of external particles
      double precision mass(5), mass2(5) 

** charge and left coupling for particle i
      double precision q(5), gL(4), gR(4)

** dummy parameters
      double complex dup27, dup28


** QED constant
      double precision Alfa, Alfa2

** Parameter, abbreviation and integral declaration
#include "mssm.inc"
#include "global.inc"

#include "abbr_NC.h"
#include "integrals_NC.h"

#include "looptools.h"

*      write(*,*)'new NC qk=4 cur 23 pent2'
      
      Alfa = AlfaQED
      Alfa2 = AlfaQED2

      dup27 = 2*mass2(1) + MH2 + 3/2.D0*mass2(3)
      dup28 = 4*mass2(1) + 3*mass2(3) + 2*MH2

      Cloop(1) = Cloop(1) + (Eval(ee12,iintP13) + Eval(ee12,iintP14))*
     -     (-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP63 + Abb4*(-S34 - T - U)))/SW) + 
     -     Eval(ee1,iintP14)*(-(2*Alfa2*MW)/SW2**2*(AbbSumP211*EL*
     -     mass2(1))/SW + 
     -     (Alfa2*MW)/SW2**2*(AbbSumP28*EL*mass2(3))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*MH2)/SW + 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP281 + Abb4*(-S + S34 + T + U)))/SW) + 
     -     (8*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(Eval(ee00,iintP13) + Eval(ee00,iintP14)))/SW + 
     -     Eval(ee13,iintP13)*
     -     (-((2*Alfa2*MW)/SW2**2*
     -     ((AbbSumP74 - Abb4*dup28)*EL)/SW) - 
     -     (6*Alfa2*MW)/SW2**2*(AbbSumP75*EL*mass2(1))/SW - 
     -     (4*Alfa2*MW)/SW2**2*(Abb4*EL*(MH2 + 3/2.D0*mass2(3)))/SW)
     -     + 1/SW2**2*(Eval(ee1,iintP13)*
     -     (Alfa2*MW*(AbbSumP34*EL) + 
     -     3*Alfa2*MW*(AbbSumP35*EL*mass2(1)) + 
     -     Alfa2*MW*(AbbSumP23*EL*mass2(3)) - 
     -     2*Alfa2*MW*(Abb4*EL*S)))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP13)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*
     -     ((Eval(ee11,iintP13) + Eval(ee11,iintP14))*
     -     (2*Alfa2*MW*(AbbSumP52*EL) + 
     -     Alfa2*(Abb4*EL*
     -     (2*MW*mass2(1) + 2*MW*mass2(3) - 2*MW*S))))/SW + 
     -     1/SW2**2*(Eval(ee0,iintP14)*
     -     (-3*Alfa2*MW*(AbbSumP35*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/SW

      Cloop(1) = Cloop(1) + 
     -     Eval(ee2,iintP14)*((Alfa2*MW)/SW2**2*(Abb2*EL*mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW - 
     -     (Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP307 + Abb4*(-S34 - T - U)))/SW) + 
     -     Eval(ee13,iintP14)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP211*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP291 + Abb4*(S34 + T + U)))/SW) + 
     -     Eval(ee14,iintP14)*
     -     ((2*Alfa2*MW)/SW2**2*(EL*(AbbSumP52 - Abb4*S))/SW + 
     -     (4*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(3) + 1/2.D0*mass2(1)))/SW) + 
     -     1/SW2**2*((Eval(ee22,iintP13) + Eval(ee22,iintP14))*
     -     (-(2*Alfa2*MW*(Abb74*EL)) + 
     -     2*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP14)*
     -     (-(2*Alfa2*MW*(Abb74*EL)) + 
     -     4*Alfa2*MW*(Abb4*EL*mass2(1))))/SW + 
     -     1/SW2**2*(Eval(ee23,iintP13)*
     -     (-(2*Alfa2*MW*(Abb74*EL)) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup27*EL)) + 2*MW*(dup28*EL))))
     -     )/SW + 1/SW2**2*
     -     (Eval(ee24,iintP13)*
     -     (-(2*Alfa2*MW*(AbbSumP110*EL)) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup27*EL)) + 2*MW*(dup28*EL))))
     -     )/SW + 1/SW2**2*
     -     (Eval(ee2,iintP13)*
     -     (-(Alfa2*MW*(AbbSumP96*EL)) + 
     -     Alfa2*MW*(AbbSumP88*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*(MW*mass2(3) - MW*S))))/SW + 
     -     1/SW2**2*(Eval(ee14,iintP13)*
     -     (2*Alfa2*MW*(AbbSumP52*EL) + 
     -     Alfa2*(Abb4*EL*(2*MW*mass2(1) - 2*MW*S))))/SW

      Cloop(1) = Cloop(1) + Eval(ee34,iintP14)*
     -     (-(4*Alfa2*MW)/SW2**2*(AbbSumP211*EL*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2**2*(Abb4*EL*(mass2(3) + MH2))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP320 + Abb4*(-S34 - T - U)))/SW) + 
     -     Eval(ee4,iintP14)*(-(3*Alfa2*MW)/SW2**2*(AbbSumP35*EL*
     -     mass2(1))/SW - 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(-mass2(3) + MH2))/SW + 
     -     (Alfa2*MW)/SW2**2*(Abb4*EL*(S34 + T + U))/SW) + 
     -     Eval(ee24,iintP14)*(-(2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*(mass2(1) + mass2(3) + MH2))/SW + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (EL*(AbbSumP314 + Abb4*(S34 + T + U)))/SW) + 
     -     (2*Alfa2*MW)/SW2**2*
     -     (Abb4*EL*mass2(1)*Eval(ee33,iintP14))/SW - 
     -     (2*Alfa2*MW)/SW2**2*
     -     (0*EL*mass(3)*Eval(ee44,iintP13))/SW + 
     -     1/SW2**2*(Eval(ee44,iintP14)*
     -     (2*Alfa2*MW*(Abb4*EL*mass2(3))))/SW + 
     -     1/SW2**2*(Eval(ee34,iintP13)*
     -     (-(2*Alfa2*MW*(Abb81*EL)) + 
     -     Alfa2*(Abb4*(-(4*MW*(dup27*EL)) + 2*MW*(dup28*EL))))
     -     )/SW + 1/SW2**2*(Eval(ee3,iintP13)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*
     -     (Eval(ee4,iintP13)*
     -     (Alfa2*(Abb4*EL*(MW*mass2(1) + MW*mass2(3) - MW*S))))
     -     /SW + 1/SW2**2*
     -     (Eval(ee3,iintP14)*
     -     (-Alfa2*MW*(AbbSumP19*EL*mass2(1)) + 
     -     Alfa2*(Abb4*EL*
     -     (-(MW*mass2(3)) - MW*MH2 + MW*(S34 + T + U)))))/SW


      end


********************************************************************************
