********************************************************************************
********************************************************************************
***                                                                          ***
*** NC_boxes.F                                                               ***
*** sophy@particle.uni-karlsruhe.de                                          ***
***                                                                          ***
*** This file contains all of the matrix elements for the electroweak        ***
*** box diagrams for the neutral currents.                                   ***
***                                                                          ***
********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams N attached to the lower quark line
*** for currents 1 and 2

      subroutine boxesBlN_12(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7 
      double complex tmp1, tmp2


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

      dup1 = 1/(CW*SW2)*
     -     (SW2*mass2(1) + SW2*mass2(3) + CW2*(mass2(1) + mass2(3)))
      dup2 = 1/(CW*SW2)*(SW2*mass2(4) + CW2*
     -     (8*mass2(1) + 4*mass2(3) + 4*mass2(4)))
      dup3 = 1/(CW*SW2)*(2*SW2*mass2(1) + 
     -     CW2*(-2*mass2(1) - 3*mass2(3) - 3*MH2) + SW2*MH2)
      dup4 = 1/(3.D0*CW*SW2)*
     -     (CW2*(-8*mass2(1) - 4*mass2(3) - 2*MH2) + 3*SW2*MH2)
      dup5 = 1/(CW*SW2)*
     -     (SW2*MH2 + CW2*(4*mass2(1) + 7*mass2(3) + MH2))
      dup6 = 1/(3.D0*CW*SW2)*
     -     (2*SW2*mass2(1) + 3*SW2*MH2 + 
     -     CW2*(16*mass2(1) + 18*mass2(3) + 6*MH2))
      dup7 = 1/(CW*SW2)*
     -     (MW2*SW2 + 2*SW2*MH2 + 
     -     CW2*(8*mass2(1) + 10*mass2(3) + 2*MH2 + 4*MW2))

      tmp1 = (Cval(cc0,iintL1) + Cval(cc0,iintL2))*
     -     (-((Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1))) - 
     -     (Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb2*EL*gR(1)))+
     -     Dval(dd0,iintL3)*
     -     (-((Alfa2*MW)/(2.D0*SW2)*(Abb4*dup7*EL*gL(1))) - 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb2*dup7*EL*gR(1)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb3*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gR(1)*mass2(3)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(1)*T24) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb2*EL*gR(1)*T24) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL28 + Abb2*(-S - U))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - U)))) + 
     -     Dval(dd0,iintL4)*((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL33*EL*gL(1)) - 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb4*dup2*EL*gL(1)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL28*EL*gR(1)) - 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb2*dup2*EL*gR(1)) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL97*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL96*EL*gR(1)*mass2(3)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(1)*(S + U)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb2*EL*gR(1)*(S + U))) + 
     -     1/SW2**2*((Dval(dd00,iintL3) + Dval(dd00,iintL4))*
     -     (-(4*Alfa2*CW*MW*(Abb4*EL*gL(1))) - 
     -     4*Alfa2*CW*MW*(Abb2*EL*gR(1))))

      tmp1 = tmp1 + Dval(dd1,iintL3)*
     -     (-((3*Alfa2*MW)/(2.D0*SW2)*(Abb4*dup4*EL*gL(1))) - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb2*dup4*EL*gR(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb3*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(1)*mass2(3)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(2*AbbSumL28 + Abb2*(S34 + T14))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(2*AbbSumL33 + Abb4*(S34 + T14))) + 
     -     (Alfa2*MW*(-4*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) + 
     -     (Alfa2*MW*(-4*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*T24)) + 
     -     Dval(dd1,iintL4)*(-((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbSumL114*EL*gL(1))) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL33*EL*gL(1)) - 
     -     (Alfa2*MW)/SW2*(Abb4*dup3*EL*gL(1)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumL112*EL*gR(1)) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL28*EL*gR(1)) - 
     -     (Alfa2*MW)/SW2*(Abb2*dup3*EL*gR(1)) - 
     -     (3*Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL102*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL100*EL*gR(1)*mass2(3)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(1)*(S34 + T14)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb2*EL*gR(1)*(S34 + T14)) + 
     -     (Alfa2*MW*(-3*CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) + 
     -     (Alfa2*MW*(-3*CW2 + SW2))/(CW*SW2**2)*(Abb2*EL*gR(1)*T24)
     -     ) + 1/SW2**2*(Dval(dd11,iintL3)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL48*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(AbbSumL45*EL*gR(1)))) + 
     -     1/SW2**2*(Dval(dd12,iintL3)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL48*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(AbbSumL45*EL*gR(1))))

      tmp1 = tmp1 + Dval(dd2,iintL4)*
     -     ((2*Alfa2*MW)/SW2*(Abb4*dup1*EL*gL(1)) + 
     -     (2*Alfa2*MW)/SW2*(Abb2*dup1*EL*gR(1)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (EL*gR(1)*(AbbSumL28 + Abb2*(-S - T - U))) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - T - U)))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (0*EL*(-gL(1) + gR(1))*mass(1)*Dval(dd22,iintL4)) + 
     -     1/SW2**2*(Dval(dd13,iintL3)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL33*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(AbbSumL28*EL*gR(1)))) + 
     -     Dval(dd2,iintL3)*((Alfa2*MW)/(CW*SW2)*(AbbSumL75*EL*gL(1)) - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb4*dup6*EL*gL(1)) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL72*EL*gR(1)) - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb2*dup6*EL*gR(1)) - 
     -     (2*Alfa2*MW)/(CW*SW2)*(AbbSumL74*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL70*EL*gR(1)*mass2(3)) + 
     -     (Alfa2*MW*(2*CW2 + 3*SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) + 
     -     (Alfa2*MW*(2*CW2 + 3*SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*T24) + 
     -     (Alfa2*MW*(8*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*(S + U)) + 
     -     (Alfa2*MW*(8*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*(S + U)) - (3*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL86 - 2/3.D0*(Abb2*T))) - 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL87 - 2/3.D0*(Abb4*T)))) + 1/SW2**2*
     -     (Dval(dd12,iintL4)*
     -     (2*Alfa2*CW*MW*(AbbSumL48*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(AbbSumL45*EL*gR(1))))

      tmp1 = tmp1 + Dval(dd3,iintL3)*
     -     (-((Alfa2*MW)/SW2*(Abb4*dup5*EL*gL(1))) - 
     -     (Alfa2*MW)/SW2*(Abb2*dup5*EL*gR(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL77*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL76*EL*gR(1)*mass2(3)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(1)*T24) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb2*EL*gR(1)*T24) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL28 + Abb2*(-S - U))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - U)))) + 
     -     Dval(dd3,iintL4)*((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL33*EL*gL(1)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL28*EL*gR(1)) + 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL97*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL96*EL*gR(1)*mass2(3)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb4*EL*gL(1)*(2*mass2(1) + mass2(3))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb2*EL*gR(1)*(2*mass2(1) + mass2(3))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(1)*(S + U)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb2*EL*gR(1)*(S + U))) + 
     -     1/SW2**2*(Dval(dd23,iintL4)*
     -     (2*Alfa2*CW*MW*(AbbSumL33*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(AbbSumL28*EL*gR(1)))) 

      tmp2 = -((3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL2*EL*q(1)*
     -     (Cval(cc0,iintL1) + Cval(cc0,iintL2)))/SW) - (4*Alfa2*MW)/
     -     SW2*(AbbSumL2*EL*q(1)*(Dval(dd00,iintL3) + 
     -     Dval(dd00,iintL4)))/SW + (2*Alfa2*MW)/SW2*
     -     (AbbSumL13*EL*q(1)*Dval(dd12,iintL4))/SW - 
     -     (2*Alfa2*MW)/SW2*(0*EL*mass(3)*q(1)*
     -     (Dval(dd11,iintL4) + Dval(dd22,iintL3)))/SW + 
     -     Dval(dd2,iintL3)*((2*Alfa2*MW)/SW2*
     -     (AbbSumL64*EL*mass2(3)*q(1))/SW - 
     -     (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(1)*(2*AbbSumL65 + 6*AbbSumL80 + 
     -     AbbSumL2*(-7*S - 4*T + T24 - 7*U)))/SW - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL2*EL*q(1)*(6*mass2(3) + MH2 + 
     -     14/3.D0*mass2(1)))/SW) + Dval(dd1,iintL4)*
     -     ((3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL93*EL*mass2(3)*q(1))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(1)*(3*AbbSumL107 + AbbSumL7 + 
     -     AbbSumL2*(2*S34 + 2*T14 + 4*T24)))/SW + (4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*q(1)*(mass2(1) + MH2 + 3/4.D0*mass2(3)))/SW)+
     -     Dval(dd1,iintL3)*
     -     (-((Alfa2*MW)/(2.D0*SW2)*(AbbSumL1*EL*mass2(3)*q(1))/SW) - 
     -     (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(1)*(8*AbbSumL7 + AbbSumL2*(4*S34 + 4*T14 + 5*T24)))/
     -     SW + (5*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*(MH2 + 8/5.D0*mass2(1) + 4/5.D0*mass2(3)))/
     -     SW) + Dval(dd0,iintL4)*
     -     (-(Alfa2*MW)/SW2*(AbbSumL90*EL*mass2(3)*q(1))/SW - 
     -     (2*Alfa2*MW)/SW2*(EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL2*EL*q(1)*
     -     (mass2(4) + 8/3.D0*mass2(1) + 4/3.D0*mass2(3)))/SW) + 
     -     Dval(dd0,iintL3)*(-((Alfa2*MW)/SW2*
     -     (AbbSumL1*EL*mass2(3)*q(1))/SW) - (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW - 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL2*EL*q(1)*
     -     (MW2 + 8/3.D0*mass2(1) + 10/3.D0*mass2(3)))/SW) + 
     -     1/SW2*(Dval(dd11,iintL3)*
     -     (-(2*Alfa2*MW*(AbbSumL13*EL*q(1)))))/SW + 
     -     1/SW2*(Dval(dd12,iintL3)*
     -     (-(2*Alfa2*MW*(AbbSumL13*EL*q(1)))))/SW + 
     -     1/SW2*(Dval(dd13,iintL3)*
     -     (-(2*Alfa2*MW*(AbbSumL7*EL*q(1)))))/SW

      tmp2 = tmp2 + Dval(dd3,iintL4)*
     -     (-(Alfa2*MW)/SW2*(AbbSumL90*EL*mass2(3)*q(1))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumL2*EL*(2*mass2(1) + mass2(3))*q(1))/
     -     SW - (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     (2*Alfa2*MW)/SW2*(AbbSumL7*EL*q(1)*Dval(dd23,iintL4))/SW + 
     -     Dval(dd3,iintL3)*((Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL66*EL*mass2(3)*q(1))/SW - (2*Alfa2*MW)/SW2*(EL*q(1)*
     -     (AbbSumL7 + AbbSumL2*(-S - U)))/SW - (7*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*q(1)*(mass2(3) + 4/7.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/(-MZ2 + T) + tmp2/T + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams N attached to the lower quark line
*** for currents 3 and 4

      subroutine boxesBlN_34(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup1, dup2, dup3, dup4, dup5, dup6, dup7
      double complex tmp1, tmp2


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

      dup1 = 1/(CW*SW2)*
     -     (SW2*mass2(1) + SW2*mass2(3) + CW2*(mass2(1) + mass2(3)))
      dup2 = 1/(CW*SW2)*
     -     (SW2*mass2(4) + CW2*
     -     (8*mass2(1) + 4*mass2(3) + 4*mass2(4)))
      dup3 = 1/(CW*SW2)*
     -     (2*SW2*mass2(1) + 
     -     CW2*(-2*mass2(1) - 3*mass2(3) - 3*MH2) + SW2*MH2)
      dup4 = 1/(3.D0*CW*SW2)*
     -     (CW2*(-8*mass2(1) - 4*mass2(3) - 2*MH2) + 3*SW2*MH2)
      dup5 = 1/(CW*SW2)*
     -     (SW2*MH2 + CW2*(4*mass2(1) + 7*mass2(3) + MH2))
      dup6 = 1/(3.D0*CW*SW2)*
     -     (2*SW2*mass2(1) + 3*SW2*MH2 + 
     -     CW2*(16*mass2(1) + 18*mass2(3) + 6*MH2))
      dup7 = 1/(CW*SW2)*
     -     (MW2*SW2 + 2*SW2*MH2 + 
     -     CW2*(8*mass2(1) + 10*mass2(3) + 2*MH2 + 4*MW2))

      tmp1 = (Cval(cc0,iintL1) + Cval(cc0,iintL2))*
     -     ((Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb4*EL*gL(1)) + 
     -     (Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb2*EL*gR(1)))+
     -     Dval(dd0,iintL3)*
     -     ((Alfa2*MW)/(2.D0*SW2)*(Abb4*dup7*EL*gL(1)) + 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb2*dup7*EL*gR(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb3*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb1*EL*gR(1)*mass2(3)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(1)*T24) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb2*EL*gR(1)*T24) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL27 + Abb2*(-S - U))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - U)))) + 
     -     Dval(dd0,iintL4)*(-((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL33*EL*gL(1))) + 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb4*dup2*EL*gL(1)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL27*EL*gR(1)) + 
     -     (Alfa2*MW)/(2.D0*SW2)*(Abb2*dup2*EL*gR(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL97*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL95*EL*gR(1)*mass2(3)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(1)*(S + U)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb2*EL*gR(1)*(S + U))) + 
     -     1/SW2**2*((Dval(dd00,iintL3) + Dval(dd00,iintL4))*
     -     (4*Alfa2*CW*MW*(Abb4*EL*gL(1)) + 
     -     4*Alfa2*CW*MW*(Abb2*EL*gR(1))))

      tmp1 = tmp1 + Dval(dd1,iintL3)*
     -     ((3*Alfa2*MW)/(2.D0*SW2)*(Abb4*dup4*EL*gL(1)) + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb2*dup4*EL*gR(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb3*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(1)*mass2(3)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(2*AbbSumL27 + Abb2*(S34 + T14))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(2*AbbSumL33 + Abb4*(S34 + T14))) - 
     -     (Alfa2*MW*(-4*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) - 
     -     (Alfa2*MW*(-4*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*T24)) + 
     -     Dval(dd1,iintL4)*((3*Alfa2*CW*MW)/SW2**2*
     -     (AbbSumL116*EL*gL(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL33*EL*gL(1)) + 
     -     (Alfa2*MW)/SW2*(Abb4*dup3*EL*gL(1)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumL114*EL*gR(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL27*EL*gR(1)) + 
     -     (Alfa2*MW)/SW2*(Abb2*dup3*EL*gR(1)) + 
     -     (3*Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL101*EL*gL(1)*mass2(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL100*EL*gR(1)*mass2(3)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(1)*(S34 + T14)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb2*EL*gR(1)*(S34 + T14)) - 
     -     (Alfa2*MW*(-3*CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) - 
     -     (Alfa2*MW*(-3*CW2 + SW2))/(CW*SW2**2)*(Abb2*EL*gR(1)*T24)
     -     ) + 1/SW2**2*(Dval(dd11,iintL3)*
     -     (2*Alfa2*CW*MW*(AbbSumL47*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(AbbSumL44*EL*gR(1))))

      tmp1 = tmp1 + 1/SW2**2*(Dval(dd13,iintL3)*
     -     (2*Alfa2*CW*MW*(AbbSumL33*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(AbbSumL27*EL*gR(1)))) + 
     -     Dval(dd2,iintL3)*(-((Alfa2*MW)/(CW*SW2)*
     -     (AbbSumL74*EL*gL(1))) + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb4*dup6*EL*gL(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL71*EL*gR(1)) + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*(Abb2*dup6*EL*gR(1)) + 
     -     (2*Alfa2*MW)/(CW*SW2)*(AbbSumL72*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL70*EL*gR(1)*mass2(3)) - 
     -     (Alfa2*MW*(2*CW2 + 3*SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*T24) - 
     -     (Alfa2*MW*(2*CW2 + 3*SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*T24) - 
     -     (Alfa2*MW*(8*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(1)*(S + U)) - 
     -     (Alfa2*MW*(8*CW2 + SW2))/(2.D0*CW*SW2**2)*
     -     (Abb2*EL*gR(1)*(S + U)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL86 - 2/3.D0*(Abb2*T))) + 
     -     (3*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL87 - 2/3.D0*(Abb4*T)))) + 
     -     1/SW2**2*(Dval(dd12,iintL4)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL47*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(AbbSumL44*EL*gR(1)))) + 
     -     1/SW2**2*(Dval(dd12,iintL3)*
     -     (2*Alfa2*CW*MW*(AbbSumL47*EL*gL(1)) + 
     -     2*Alfa2*CW*MW*(AbbSumL44*EL*gR(1))))

      tmp1 = tmp1 + Dval(dd2,iintL4)*
     -     (-((2*Alfa2*MW)/SW2*(Abb4*dup1*EL*gL(1))) - 
     -     (2*Alfa2*MW)/SW2*(Abb2*dup1*EL*gR(1)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (EL*gR(1)*(AbbSumL27 + Abb2*(-S - T - U))) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - T - U)))) + 
     -     Dval(dd3,iintL3)*((Alfa2*MW)/SW2*(Abb4*dup5*EL*gL(1)) + 
     -     (Alfa2*MW)/SW2*(Abb2*dup5*EL*gR(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL76*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL75*EL*gR(1)*mass2(3)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(1)*T24) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb2*EL*gR(1)*T24) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(1)*(AbbSumL27 + Abb2*(-S - U))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(1)*(AbbSumL33 + Abb4*(-S - U))) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (0*EL*mass(1)*mass(3)*
     -     (-(SW2*gL(1)) + CW2*(2*gL(1) - 2*gR(1)) + SW2*gR(1)))) + 
     -     1/SW2**2*(Dval(dd23,iintL4)*
     -     (-(2*Alfa2*CW*MW*(AbbSumL33*EL*gL(1))) - 
     -     2*Alfa2*CW*MW*(AbbSumL27*EL*gR(1)))) 

      tmp1 = tmp1 + Dval(dd3,iintL4)*
     -     (-((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL33*EL*gL(1))) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (AbbSumL27*EL*gR(1)) - 
     -     (Alfa2*MW)/(CW*SW2)*(AbbSumL97*EL*gL(1)*mass2(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumL95*EL*gR(1)*mass2(3)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb4*EL*gL(1)*(2*mass2(1) + mass2(3))) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (Abb2*EL*gR(1)*(2*mass2(1) + mass2(3))) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(1)*(S + U)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb2*EL*gR(1)*(S + U))) 

      tmp2 = (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*(Cval(cc0,iintL1) + Cval(cc0,iintL2)))/SW + 
     -     (4*Alfa2*MW)/SW2*(AbbSumL2*EL*q(1)*
     -     (Dval(dd00,iintL3) + Dval(dd00,iintL4)))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumL13*EL*q(1)*Dval(dd12,iintL4))/SW + 
     -     Dval(dd2,iintL3)*(-((2*Alfa2*MW)/SW2*
     -     (AbbSumL63*EL*mass2(3)*q(1))/SW) + 
     -     (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(1)*(2*AbbSumL65 + 6*AbbSumL79 + 
     -     AbbSumL2*(-7*S - 4*T + T24 - 7*U)))/SW + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*(6*mass2(3) + MH2 + 14/3.D0*mass2(1)))/SW
     -     ) + Dval(dd1,iintL4)*
     -     (-(3*Alfa2*MW)/(2.D0*SW2)*(AbbSumL92*EL*mass2(3)*q(1))/SW + 
     -     (Alfa2*MW)/SW2*
     -     (EL*q(1)*(3*AbbSumL108 + AbbSumL7 + 
     -     AbbSumL2*(2*S34 + 2*T14 + 4*T24)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*q(1)*(mass2(1) + MH2 + 3/4.D0*mass2(3)))/SW)+
     -     Dval(dd1,iintL3)*
     -     ((Alfa2*MW)/(2.D0*SW2)*(AbbSumL1*EL*mass2(3)*q(1))/SW + 
     -     (Alfa2*MW)/(2.D0*SW2)*
     -     (EL*q(1)*(8*AbbSumL7 + AbbSumL2*(4*S34 + 4*T14 + 5*T24)))/
     -     SW - (5*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*(MH2 + 8/5.D0*mass2(1) + 4/5.D0*mass2(3)))/
     -     SW) + Dval(dd0,iintL4)*
     -     ((Alfa2*MW)/SW2*(AbbSumL90*EL*mass2(3)*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*
     -     (mass2(4) + 8/3.D0*mass2(1) + 4/3.D0*mass2(3)))/SW) + 
     -     Dval(dd0,iintL3)*((Alfa2*MW)/SW2*
     -     (AbbSumL1*EL*mass2(3)*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW + 
     -     (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL2*EL*q(1)*
     -     (MW2 + 8/3.D0*mass2(1) + 10/3.D0*mass2(3)))/SW) + 
     -     1/SW2*(Dval(dd11,iintL3)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*q(1))))/SW + 
     -     1/SW2*(Dval(dd12,iintL3)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*q(1))))/SW + 
     -     1/SW2*(Dval(dd13,iintL3)*
     -     (2*Alfa2*MW*(AbbSumL7*EL*q(1))))/SW

      tmp2 = tmp2 + Dval(dd3,iintL4)*
     -     ((Alfa2*MW)/SW2*(AbbSumL90*EL*mass2(3)*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*(2*mass2(1) + mass2(3))*q(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) - 
     -     (2*Alfa2*MW)/SW2*(AbbSumL7*EL*q(1)*Dval(dd23,iintL4))/SW + 
     -     Dval(dd3,iintL3)*(-((Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumL66*EL*mass2(3)*q(1))/SW) + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW + 
     -     (7*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*q(1)*(mass2(3) + 4/7.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/(-MZ2 + T) + tmp2/T + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams P attached to the lower quark line
*** for currents 1 and 2

      subroutine boxesBlP_12(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup10, dup11, dup12, dup13, dup14, dup15, dup16
      double complex dup17, dup18, dup19, dup20, dup8, dup9, tmp3
      double complex tmp4


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

      dup8 = -(gL(3)*gR(1)) + gL(1)*(gL(3) - gR(3)) + gR(1)*gR(3)
      dup9 = -(gL(3)**3*gR(1)) + gR(1)*gR(3)**3 + gL(1)*(gL(3)**3 - 
     -     gR(3)**3)
      dup10 = MH2 + 2/3.D0*mass2(1) + 2/3.D0*mass2(3)
      dup11 = gR(4)*mass2(4) + gL(4)*(2*mass2(1) + mass2(3) - MW2)
      dup12 = -2*mass2(1) - 2*mass2(3) + MZ2
      dup13 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + mass2(3) - MZ2)
      dup14 = -(gL(3)*mass2(3)) + gR(3)*(-2*mass2(1) - mass2(3) + MZ2)
      dup15 = AbbSumL28 + Abb2*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup16 = AbbSumL33 + Abb4*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup17 = AbbSumL7 + AbbSumL2*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup18 = AbbSumL28 + Abb2*(-S - U - 1/2.D0*T)
      dup19 = AbbSumL33 + Abb4*(-S - U - 1/2.D0*T)
      dup20 = AbbSumL7 + AbbSumL2*(-S - U - 1/2.D0*T)

      tmp3 = Dval(dd0,iintL8)*
     -     ((4*Alfa2*MW)/CW2*(Abb4*dup13*EL*gL(1)*gL(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup13*EL*gL(3)**2*gR(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup14*EL*gL(1)*gR(3)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup14*EL*gR(1)*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*(AbbSumL24 + Abb1*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**3*gR(1)*(AbbSumL28 + Abb2*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*(AbbSumL30 + Abb3*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)**3*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     Dval(dd0,iintL7)*((2*Alfa2*MW)/SW2*
     -     (Abb4*dup11*EL*gL(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb2*dup11*EL*gR(1))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*gL(4)*gR(1)*(AbbSumL28 + Abb2*(-S - U)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/SW2*(Cval(cc0,iintL5)*
     -     (-(2*Alfa2*MW*(Abb4*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/SW2*(Dval(dd00,iintL7)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(1)*gL(4)) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Cval(cc0,iintL6)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3)) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/CW2*(Dval(dd00,iintL8)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW

      tmp3 = tmp3 + Dval(dd11,iintL8)*
     -     ((4*Alfa2*MW)/CW2*(AbbSumL48*EL*gL(1)*gL(3)**3)/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumL45*EL*gL(3)**3*gR(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumL46*EL*gL(1)*gR(3)**3)/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumL43*EL*gR(1)*gR(3)**3)/SW) + 
     -     1/SW2*(Dval(dd11,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL48*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(AbbSumL45*EL*gL(4)*gR(1))))/SW + 
     -     Dval(dd1,iintL8)*(-((12*Alfa2*MW)/CW2*
     -     (Abb4*dup10*EL*gL(1)*gL(3)**3)/SW) + 
     -     (8*Alfa2*MW)/CW2*(dup16*EL*gL(1)*gL(3)**3)/SW - 
     -     (12*Alfa2*MW)/CW2*(Abb2*dup10*EL*gL(3)**3*gR(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(dup15*EL*gL(3)**3*gR(1))/SW - 
     -     (12*Alfa2*MW)/CW2*(Abb3*dup10*EL*gL(1)*gR(3)**3)/SW - 
     -     (12*Alfa2*MW)/CW2*(Abb1*dup10*EL*gR(1)*gR(3)**3)/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)**3*
     -     (AbbSumL24 + Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW+
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)**3*
     -     (AbbSumL30 + Abb3*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW)+
     -     1/SW2*(Dval(dd1,iintL7)*
     -     (-(6*Alfa2*MW*(Abb4*dup10*EL*gL(1)*gL(4))) - 
     -     6*Alfa2*MW*(Abb2*dup10*EL*gL(4)*gR(1)) + 
     -     Alfa2*(4*MW*(dup16*EL*gL(1)*gL(4)) + 
     -     4*MW*(dup15*EL*gL(4)*gR(1)))))/SW

      tmp3 = tmp3 + Dval(dd2,iintL7)*
     -     ((4*Alfa2*MW)/SW2*(Abb4*EL*gL(1)*gL(4)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/SW2*(Abb2*EL*gL(4)*gR(1)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(4)*gR(1)*(AbbSumL28 + Abb2*(-S - U)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/SW2*(Dval(dd12,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(AbbSumL28*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Dval(dd12,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL28*EL*gL(3)**3*gR(1)) + 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL24*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/CW2*(Dval(dd13,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL56*EL*gL(1)*gL(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL54*EL*gL(3)**3*gR(1)) + 
     -     4*Alfa2*MW*(AbbSumL55*EL*gL(1)*gR(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL53*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/SW2*(Dval(dd13,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL56*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(AbbSumL54*EL*gL(4)*gR(1))))/SW

      tmp3 = tmp3 + Dval(dd2,iintL8)*
     -     ((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)**3*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**3*gR(1)*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)**3*mass2(1))/SW + 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*(AbbSumL24 + Abb1*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**3*gR(1)*(AbbSumL28 + Abb2*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*(AbbSumL30 + Abb3*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)**3*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/SW2*(Dval(dd23,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(AbbSumL28*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Dval(dd23,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL28*EL*gL(3)**3*gR(1)) + 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL24*EL*gR(1)*gR(3)**3)))/SW

      tmp3 = tmp3 + Dval(dd3,iintL7)*
     -     ((4*Alfa2*MW)/SW2*(dup19*EL*gL(1)*gL(4))/SW + 
     -     (4*Alfa2*MW)/SW2*(dup18*EL*gL(4)*gR(1))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + mass2(3)))/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + mass2(3)))/SW) + 
     -     Dval(dd3,iintL8)*((8*Alfa2*MW)/CW2*
     -     (dup19*EL*gL(1)*gL(3)**3)/SW + (8*Alfa2*MW)/CW2*(dup18*EL*
     -     gL(3)**3*gR(1))/SW + (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)**3*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**3*gR(1)*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)**3*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)**3*(2*mass2(1) + mass2(3)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)*gR(3)**3*(AbbSumL24 + Abb1*(-S - U -
     -      1/2.D0*T)))/SW + (8*Alfa2*MW)/CW2*(EL*gL(1)*gR(3)**3*
     -     (AbbSumL30 + Abb3*(-S - U - 1/2.D0*T)))/SW) + 1/SW2*
     -     (Dval(dd33,iintL7)*(2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(AbbSumL28*EL*gL(4)*gR(1))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd33,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL28*EL*gL(3)**3*gR(1)) + 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) + 
     -     4*Alfa2*MW*(AbbSumL24*EL*gR(1)*gR(3)**3)))/SW

      tmp4 = Dval(dd11,iintL8)*
     -     ((4*Alfa2*MW)/CW2*(AbbSumL13*EL*gL(3)**2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumL11*EL*gR(3)**2*q(1)*q(3))/SW) + 
     -     Dval(dd0,iintL8)*(-((4*Alfa2*MW)/CW2*
     -     (AbbSumL2*dup12*EL*gL(3)**2*q(1)*q(3))/SW) - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL1*dup12*EL*gR(3)**2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*(AbbSumL4 + AbbSumL1*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*q(1)*q(3)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     Dval(dd0,iintL7)*(-(2*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*(-2*mass2(1) - mass2(3) - mass2(4) + MW2)*
     -     q(1)*q(4))/SW + (2*Alfa2*MW)/SW2*(EL*q(1)*q(4)*(AbbSumL7 + 
     -     AbbSumL2*(-S - U)))/SW) - (2*Alfa2*MW)/SW2*(AbbSumL2*EL*q(1)*
     -     q(4)*Cval(cc0,iintL5))/SW + (4*Alfa2*MW)/SW2*(AbbSumL2*EL*
     -     q(1)*q(4)*Dval(dd00,iintL7))/SW + 1/CW2*(Cval(cc0,iintL6)*
     -     (-(4*Alfa2*MW*(AbbSumL2*EL*gL(3)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumL1*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd00,iintL8)*
     -     (8*Alfa2*MW*(AbbSumL2*EL*gL(3)**2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumL1*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd1,iintL7)*
     -     (-(6*Alfa2*MW*(AbbSumL2*dup10*EL*q(1)*q(4))) + 
     -     4*Alfa2*MW*(dup17*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Dval(dd11,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL13*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Dval(dd12,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4))))/SW + 
     -     Dval(dd1,iintL8)*(-((12*Alfa2*MW)/CW2*
     -     (AbbSumL2*dup10*EL*gL(3)**2*q(1)*q(3))/SW) + 
     -     (8*Alfa2*MW)/CW2*(dup17*EL*gL(3)**2*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*(AbbSumL1*dup10*EL*gR(3)**2*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gR(3)**2*q(1)*q(3)*
     -     (AbbSumL4 + AbbSumL1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW)

      tmp4 = tmp4 + Dval(dd3,iintL7)*
     -     ((4*Alfa2*MW)/SW2*(dup20*EL*q(1)*q(4))/SW + (4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*(2*mass2(1) + mass2(3))*q(1)*q(4))/SW) + 
     -     Dval(dd2,iintL8)*((8*Alfa2*MW)/CW2*
     -     (AbbSumL2*EL*gL(3)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumL1*EL*gR(3)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*(AbbSumL4 + AbbSumL1*(-S - U)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*q(1)*q(3)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     Dval(dd2,iintL7)*((4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*mass2(1)*q(1)*q(4))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     1/CW2*(Dval(dd13,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL17*EL*gL(3)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumL16*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL17*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Dval(dd23,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4))))/SW

      tmp4 = tmp4 + 1/CW2*(Dval(dd33,iintL8)*
     -     (4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd33,iintL7)*
     -     (2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4))))/SW + 
     -     Dval(dd3,iintL8)*((8*Alfa2*MW)/CW2*
     -     (dup20*EL*gL(3)**2*q(1)*q(3))/SW + (8*Alfa2*MW)/CW2*
     -     (AbbSumL2*EL*gL(3)**2*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW+
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumL1*EL*gR(3)**2*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW+
     -     (8*Alfa2*MW)/CW2*(EL*gR(3)**2*q(1)*q(3)*
     -     (AbbSumL4 + AbbSumL1*(-S - U - 1/2.D0*T)))/SW)

      Cloop(1) = tmp3/(-MZ2 + T) + tmp4/T + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams P attached to the lower quark line
*** for currents 3 and 4 

      subroutine boxesBlP_34(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup10, dup11, dup12, dup13, dup14, dup15, dup16
      double complex dup17, dup18, dup8, dup9, tmp3, tmp4


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

      dup8 = MH2 + 2/3.D0*mass2(1) + 2/3.D0*mass2(3)
      dup9 = gR(4)*mass2(4) + gL(4)*(2*mass2(1) + mass2(3) - MW2)
      dup10 = -2*mass2(1) - 2*mass2(3) + MZ2
      dup11 = gR(3)*mass2(3) + gL(3)*(2*mass2(1) + mass2(3) - MZ2)
      dup12 = -(gL(3)*mass2(3)) + gR(3)*(-2*mass2(1) - mass2(3) + MZ2)
      dup13 = AbbSumL27 + Abb2*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup14 = AbbSumL33 + Abb4*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup15 = AbbSumL7 + AbbSumL2*(T24 + 1/2.D0*S34 + 1/2.D0*T14)
      dup16 = AbbSumL27 + Abb2*(-S - U - 1/2.D0*T)
      dup17 = AbbSumL33 + Abb4*(-S - U - 1/2.D0*T)
      dup18 = AbbSumL7 + AbbSumL2*(-S - U - 1/2.D0*T)

      tmp3 = Dval(dd0,iintL8)*
     -     (-((4*Alfa2*MW)/CW2*(Abb4*dup11*EL*gL(1)*gL(3)**2)/SW) - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup11*EL*gL(3)**2*gR(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup12*EL*gL(1)*gR(3)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup12*EL*gR(1)*gR(3)**2)/SW -
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*(AbbSumL23 + Abb1*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**3*gR(1)*(AbbSumL27 + Abb2*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*(AbbSumL30 + Abb3*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)**3*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     Dval(dd0,iintL7)*(-((2*Alfa2*MW)/SW2*
     -     (Abb4*dup9*EL*gL(1))/SW) - 
     -     (2*Alfa2*MW)/SW2*(Abb2*dup9*EL*gR(1))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(4)*gR(1)*(AbbSumL27 + Abb2*(-S - U)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/SW2*(Cval(cc0,iintL5)*
     -     (2*Alfa2*MW*(Abb4*EL*gL(1)*gL(4)) + 
     -     2*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/SW2*(Dval(dd00,iintL7)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(1)*gL(4))) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Cval(cc0,iintL6)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/CW2*(Dval(dd00,iintL8)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)*gL(3)**3)) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)**3*gR(1)) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)*gR(3)**3) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)*gR(3)**3)))/SW

      tmp3 = tmp3 + Dval(dd11,iintL8)*
     -     (-((4*Alfa2*MW)/CW2*(AbbSumL47*EL*gL(1)*gL(3)**3)/SW) - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL44*EL*gL(3)**3*gR(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL46*EL*gL(1)*gR(3)**3)/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL43*EL*gR(1)*gR(3)**3)/SW) + 
     -     1/SW2*(Dval(dd11,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL47*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(AbbSumL44*EL*gL(4)*gR(1))))/SW + 
     -     1/SW2*(Dval(dd12,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(AbbSumL27*EL*gL(4)*gR(1))))/SW + 
     -     Dval(dd1,iintL8)*(-((8*Alfa2*MW)/CW2*
     -     (dup14*EL*gL(1)*gL(3)**3)/SW) + 
     -     (12*Alfa2*MW)/CW2*(Abb4*dup8*EL*gL(1)*gL(3)**3)/SW - 
     -     (8*Alfa2*MW)/CW2*(dup13*EL*gL(3)**3*gR(1))/SW + 
     -     (12*Alfa2*MW)/CW2*(Abb2*dup8*EL*gL(3)**3*gR(1))/SW + 
     -     (12*Alfa2*MW)/CW2*(Abb3*dup8*EL*gL(1)*gR(3)**3)/SW + 
     -     (12*Alfa2*MW)/CW2*(Abb1*dup8*EL*gR(1)*gR(3)**3)/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*
     -     (AbbSumL23 + Abb1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*
     -     (AbbSumL30 + Abb3*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW)+
     -     1/SW2*(Dval(dd1,iintL7)*
     -     (6*Alfa2*MW*(Abb4*dup8*EL*gL(1)*gL(4)) + 
     -     6*Alfa2*MW*(Abb2*dup8*EL*gL(4)*gR(1)) + 
     -     Alfa2*(-(4*MW*(dup14*EL*gL(1)*gL(4))) - 
     -     4*MW*(dup13*EL*gL(4)*gR(1)))))/SW

      tmp3 = tmp3 + Dval(dd13,iintL8)*
     -     (-((4*Alfa2*MW)/CW2*(AbbSumL56*EL*gL(1)*gL(3)**3)/SW) - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL54*EL*gL(3)**3*gR(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL55*EL*gL(1)*gR(3)**3)/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL53*EL*gR(1)*gR(3)**3)/SW) + 
     -     Dval(dd2,iintL7)*(-((4*Alfa2*MW)/SW2*
     -     (Abb4*EL*gL(1)*gL(4)*mass2(1))/SW) - 
     -     (4*Alfa2*MW)/SW2*(Abb2*EL*gL(4)*gR(1)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(4)*gR(1)*(AbbSumL27 + Abb2*(-S - U)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(1)*gL(4)*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/CW2*(Dval(dd12,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3)) - 
     -     4*Alfa2*MW*(AbbSumL27*EL*gL(3)**3*gR(1)) - 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) - 
     -     4*Alfa2*MW*(AbbSumL23*EL*gR(1)*gR(3)**3)))/SW + 
     -     1/SW2*(Dval(dd13,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL56*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(AbbSumL54*EL*gL(4)*gR(1))))/SW

      tmp3 = tmp3 + Dval(dd2,iintL8)*
     -     (-((8*Alfa2*MW)/CW2*(Abb4*EL*gL(1)*gL(3)**3*mass2(1))/SW) - 
     -     (8*Alfa2*MW)/CW2*(Abb2*EL*gL(3)**3*gR(1)*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb3*EL*gL(1)*gR(3)**3*mass2(1))/SW - 
     -     (8*Alfa2*MW)/CW2*(Abb1*EL*gR(1)*gR(3)**3*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*(AbbSumL23 + Abb1*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**3*gR(1)*(AbbSumL27 + Abb2*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*(AbbSumL30 + Abb3*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gL(3)**3*(AbbSumL33 + Abb4*(-S - U)))/SW) + 
     -     1/SW2*(Dval(dd23,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(AbbSumL27*EL*gL(4)*gR(1))))/SW + 
     -     1/CW2*(Dval(dd23,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3)) - 
     -     4*Alfa2*MW*(AbbSumL27*EL*gL(3)**3*gR(1)) - 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) - 
     -     4*Alfa2*MW*(AbbSumL23*EL*gR(1)*gR(3)**3)))/SW

      tmp3 = tmp3 + Dval(dd3,iintL7)*
     -     (-((4*Alfa2*MW)/SW2*(dup17*EL*gL(1)*gL(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*(dup16*EL*gL(4)*gR(1))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb4*EL*gL(1)*gL(4)*(2*mass2(1) + mass2(3)))/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (Abb2*EL*gL(4)*gR(1)*(2*mass2(1) + mass2(3)))/SW) + 
     -     Dval(dd3,iintL8)*(-((8*Alfa2*MW)/CW2*
     -     (dup17*EL*gL(1)*gL(3)**3)/SW) - 
     -     (8*Alfa2*MW)/CW2*(dup16*EL*gL(3)**3*gR(1))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb4*EL*gL(1)*gL(3)**3*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb2*EL*gL(3)**3*gR(1)*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb3*EL*gL(1)*gR(3)**3*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (Abb1*EL*gR(1)*gR(3)**3*(2*mass2(1) + mass2(3)))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(1)*gR(3)**3*(AbbSumL23 + Abb1*(-S - U - 1/2.D0*T)))/SW
     -     - (8*Alfa2*MW)/CW2*
     -     (EL*gL(1)*gR(3)**3*(AbbSumL30 + Abb3*(-S - U - 1/2.D0*T)))/SW
     -     ) + 1/SW2*(Dval(dd33,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(4))) - 
     -     2*Alfa2*MW*(AbbSumL27*EL*gL(4)*gR(1))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd33,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL33*EL*gL(1)*gL(3)**3)) - 
     -     4*Alfa2*MW*(AbbSumL27*EL*gL(3)**3*gR(1)) - 
     -     4*Alfa2*MW*(AbbSumL30*EL*gL(1)*gR(3)**3) - 
     -     4*Alfa2*MW*(AbbSumL23*EL*gR(1)*gR(3)**3)))/SW

      tmp4 = Dval(dd11,iintL8)*
     -     (-((4*Alfa2*MW)/CW2*(AbbSumL13*EL*gL(3)**2*q(1)*q(3))/SW) - 
     -     (4*Alfa2*MW)/CW2*(AbbSumL12*EL*gR(3)**2*q(1)*q(3))/SW) + 
     -     Dval(dd0,iintL8)*((4*Alfa2*MW)/CW2*
     -     (AbbSumL2*dup10*EL*gL(3)**2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumL1*dup10*EL*gR(3)**2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*(AbbSumL4 + AbbSumL1*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*q(1)*q(3)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     Dval(dd0,iintL7)*((2*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*(-2*mass2(1) - mass2(3) - mass2(4) + MW2)*
     -     q(1)*q(4))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     (2*Alfa2*MW)/SW2*(AbbSumL2*EL*q(1)*q(4)*Cval(cc0,iintL5))/SW - 
     -     (4*Alfa2*MW)/SW2*(AbbSumL2*EL*q(1)*q(4)*Dval(dd00,iintL7))/SW + 
     -     1/CW2*(Cval(cc0,iintL6)*
     -     (4*Alfa2*MW*(AbbSumL2*EL*gL(3)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumL1*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd00,iintL8)*
     -     (-(8*Alfa2*MW*(AbbSumL2*EL*gL(3)**2*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumL1*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd1,iintL7)*
     -     (-(4*Alfa2*MW*(dup15*EL*q(1)*q(4))) + 
     -     6*Alfa2*MW*(AbbSumL2*dup8*EL*q(1)*q(4))))/SW + 
     -     1/SW2*(Dval(dd11,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL13*EL*q(1)*q(4)))))/SW + 
     -     1/SW2*(Dval(dd12,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4)))))/SW + 
     -     Dval(dd1,iintL8)*(-((8*Alfa2*MW)/CW2*
     -     (dup15*EL*gL(3)**2*q(1)*q(3))/SW) + 
     -     (12*Alfa2*MW)/CW2*(AbbSumL2*dup8*EL*gL(3)**2*q(1)*q(3))/SW + 
     -     (12*Alfa2*MW)/CW2*(AbbSumL1*dup8*EL*gR(3)**2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*
     -     (AbbSumL4 + AbbSumL1*(T24 + 1/2.D0*S34 + 1/2.D0*T14)))/SW
     -     )

      tmp4 = tmp4 + Dval(dd3,iintL7)*
     -     (-((4*Alfa2*MW)/SW2*(dup18*EL*q(1)*q(4))/SW) - 
     -     (4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*(2*mass2(1) + mass2(3))*q(1)*q(4))/SW) + 
     -     Dval(dd2,iintL8)*(-((8*Alfa2*MW)/CW2*
     -     (AbbSumL2*EL*gL(3)**2*mass2(1)*q(1)*q(3))/SW) - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumL1*EL*gR(3)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*(AbbSumL4 + AbbSumL1*(-S - U)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)**2*q(1)*q(3)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     Dval(dd2,iintL7)*(-((4*Alfa2*MW)/SW2*
     -     (AbbSumL2*EL*mass2(1)*q(1)*q(4))/SW) - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(1)*q(4)*(AbbSumL7 + AbbSumL2*(-S - U)))/SW) + 
     -     1/CW2*(Dval(dd13,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL17*EL*gL(3)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumL16*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL17*EL*q(1)*q(4)))))/SW + 
     -     1/SW2*(Dval(dd23,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4)))))/SW

      tmp4 = tmp4 + 1/CW2*(Dval(dd33,iintL8)*
     -     (-(4*Alfa2*MW*(AbbSumL7*EL*gL(3)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumL4*EL*gR(3)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd33,iintL7)*
     -     (-(2*Alfa2*MW*(AbbSumL7*EL*q(1)*q(4)))))/SW + 
     -     Dval(dd3,iintL8)*(-((8*Alfa2*MW)/CW2*
     -     (dup18*EL*gL(3)**2*q(1)*q(3))/SW) - 
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumL2*EL*gL(3)**2*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (AbbSumL1*EL*gR(3)**2*(2*mass2(1) + mass2(3))*q(1)*q(3))/SW-
     -     (8*Alfa2*MW)/CW2*
     -     (EL*gR(3)**2*q(1)*q(3)*
     -     (AbbSumL4 + AbbSumL1*(-S - U - 1/2.D0*T)))/SW)

      Cloop(1) = tmp3/(-MZ2 + T) + tmp4/T + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams N attached to the upper quark line
*** for currents 1 and 3

      subroutine boxesBuN_13(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup1, dup2, tmp1, tmp2


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

      dup1 = 1/(2.D0*CW*SW2)*(SW2*mass2(1) + SW2*mass2(2) + CW2*(4*
     -     mass2(2) - 2*MH2) + 2*SW2*MH2)
      dup2 = 2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U

      tmp1 = (Cval(cc0,iintU1) + Cval(cc0,iintU2))*(-((Alfa2*MW*(4*CW2 + 
     -     SW2))/(2.D0*CW*SW2**2)*(Abb4*EL*gL(3))) - 
     -     (Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb3*EL*gR(3)))+
     -     Dval(dd0,iintU3)*
     -     (-((Alfa2*MW)/SW2*(Abb4*dup1*EL*gL(3))) - 
     -     (Alfa2*MW)/SW2*(Abb3*dup1*EL*gR(3)) + 
     -     (6*Alfa2*CW*MW)/SW2**2*(AbbSumU62*EL*gL(3)*mass2(1)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU60*EL*gR(3)*mass2(1)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*S34) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*S34) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gR(3)*(Abb62 - Abb3*T)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gL(3)*(Abb82 - Abb4*T))) + 
     -     1/SW2**2*((Dval(dd00,iintU3) + Dval(dd00,iintU4))*
     -     (-(4*Alfa2*CW*MW*(Abb4*EL*gL(3))) - 
     -     4*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 
     -     Dval(dd0,iintU4)*(-((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb82*EL*gL(3))) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb62*EL*gR(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb2*EL*gL(3)*mass2(1)) + 
     -     (7*Alfa2*CW*MW)/SW2**2*(AbbSumU31*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(3)*mass2(1)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU29*EL*gR(3)*mass2(1)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*(S34 + T)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*(S34 + T)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*(Abb4*EL*gL(3)*(SW2*mass2(2) + 
     -     CW2*(4*mass2(2) - 4*MH2))) - (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*mass2(2) + 
     -     CW2*(4*mass2(2) - 4*MH2))))

      tmp1 = tmp1 + Dval(dd1,iintU4)*
     -     ((Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU14*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gR(3)*mass2(1)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(EL*gR(3)*(Abb62 + Abb3*T14)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*(EL*gL(3)*(Abb82 + Abb4*T14))) + 
     -     Dval(dd1,iintU3)*(-(Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(3)*
     -     MH2) - (Alfa2*MW)/(CW*SW2**2)*(Abb3*EL*gR(3)*MH2) + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(3)*(Abb62 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14))) + 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(Abb82 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))) + 
     -     1/SW2**2*(Dval(dd11,iintU3)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(AbbSumU37*EL*gR(3)))) + 
     -     1/SW2**2*(Dval(dd12,iintU4)*
     -     (2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3)) + 2*Alfa2*CW*MW*
     -     (AbbSumU37*EL*gR(3)))) + 1/SW2**2*(Dval(dd13,iintU3)*
     -     (2*Alfa2*CW*MW*(Abb82*EL*gL(3)) + 
     -     2*Alfa2*CW*MW*(Abb62*EL*gR(3)))) + 1/SW2**2*
     -     (Dval(dd12,iintU3)*(-(2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(AbbSumU37*EL*gR(3))))

      tmp1 = tmp1 + 
     -     Dval(dd2,iintU3)*(-((Alfa2*MW)/(CW*SW2)*(Abb45*EL*gL(3))) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb41*EL*gR(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU14*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gR(3)*mass2(1)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb4*EL*gL(3)*S34) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb3*EL*gR(3)*S34) + 
     -     (Alfa2*CW*MW)/SW2**2*(EL*gR(3)*(AbbSumU70 + 2*Abb3*T14)) + 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumU72 + 2*Abb4*T14)) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*MH2 + CW2*(mass2(1) + MH2))) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*MH2 + CW2*(mass2(1) + MH2)))) + 
     -     Dval(dd2,iintU4)*(-((Alfa2*MW)/(CW*SW2)*(Abb82*EL*gL(3))) - 
     -     (Alfa2*MW)/(CW*SW2)*(Abb62*EL*gR(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU14*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gR(3)*mass2(1)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(3)*(AbbSumU41 - Abb3*S34)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumU43 - Abb4*S34)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb4*EL*gL(3)*T14) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb3*EL*gR(3)*T14) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*MH2 + CW2*(mass2(1) + MH2))) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*MH2 + CW2*(mass2(1) + MH2)))) + 
     -     1/SW2**2*(Dval(dd13,iintU4)*
     -     (-(2*Alfa2*CW*MW*(Abb82*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(Abb62*EL*gR(3))))

      tmp1 = tmp1 + Dval(dd3,iintU4)*
     -     (-((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb82*EL*gL(3))) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb62*EL*gR(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb2*EL*gL(3)*mass2(1)) + 
     -     (7*Alfa2*CW*MW)/SW2**2*(AbbSumU31*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(3)*mass2(1)) + 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU29*EL*gR(3)*mass2(1)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*MH2) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*MH2) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*(S34 + T)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*(S34 + T))) + 
     -     Dval(dd3,iintU3)*((6*Alfa2*CW*MW)/SW2**2*(AbbSumU62*EL*gL(3)*
     -     mass2(1)) + (3*Alfa2*CW*MW)/SW2**2*(AbbSumU60*EL*gR(3)*
     -     mass2(1)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb4*EL*gL(3)*S34) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb3*EL*gR(3)*S34) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gR(3)*(Abb62 - Abb3*T)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gL(3)*(Abb82 - Abb4*T)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*mass2(1) - 2*CW2*MH2 + 2*SW2*MH2))-
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*mass2(1) - 2*CW2*MH2 + 2*SW2*MH2)))

      tmp2 = Dval(dd2,iintU4)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     ((-2*AbbSumU1 + AbbSumU5)*EL*mass2(1)*q(3))/SW + 
     -     (Alfa2*MW)/SW2*(EL*q(3)*(2*AbbSumU25 + AbbSumU4 + 
     -     AbbSumU1*(-2*S34 + 2*T14)))/SW) - (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumU1*EL*q(3)*(Cval(cc0,iintU1) + Cval(cc0,iintU2)))/SW - 
     -     (4*Alfa2*MW)/SW2*(AbbSumU1*EL*q(3)*
     -     (Dval(dd00,iintU3) + Dval(dd00,iintU4)))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumU23*EL*q(3)*Dval(dd11,iintU3))/SW + 
     -     Dval(dd2,iintU3)*((Alfa2*MW)/(2.D0*SW2)*((-2*AbbSumU1 + 
     -     AbbSumU5)*EL*mass2(1)*q(3))/SW + (Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU50 + AbbSumU58 + 
     -     AbbSumU1*(-2*S34 + 2*T14)))/SW) + 
     -     Dval(dd0,iintU4)*((Alfa2*MW)/(2.D0*SW2)*
     -     ((14*AbbSumU19 + AbbSumU3)*EL*mass2(1)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 - 3/4.D0*mass2(2)))/SW) + 
     -     Dval(dd0,iintU3)*((6*Alfa2*MW)/SW2*
     -     (AbbSumU53*EL*mass2(1)*q(3))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 + 1/4.D0*mass2(1) - 3/4.D0*mass2(2)))/
     -     SW) + 1/SW2*(Dval(dd13,iintU4)*
     -     (-(2*Alfa2*MW*(AbbSumU4*EL*q(3)))))/SW + 
     -     1/SW2*(Dval(dd12,iintU3)*
     -     (-(2*Alfa2*MW*(AbbSumU23*EL*q(3)))))/SW + 
     -     1/SW2*(Dval(dd12,iintU4)*
     -     (2*Alfa2*MW*(AbbSumU23*EL*q(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintU3)*
     -     (2*Alfa2*MW*(AbbSumU4*EL*q(3))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd1,iintU4)*
     -     (2*Alfa2*MW*(AbbSumU1*EL*mass2(1)*q(3)) - 
     -     Alfa2*MW*(AbbSumU5*EL*mass2(1)*q(3))))/SW + 
     -     Dval(dd1,iintU3)*((4*Alfa2*MW)/SW2*(EL*q(3)*(AbbSumU4 + 
     -     AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW)

      tmp2 = tmp2 + Dval(dd3,iintU4)*((Alfa2*MW)/(2.D0*SW2)*
     -     ((14*AbbSumU19 + AbbSumU3)*EL*mass2(1)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW) + 
     -     Dval(dd3,iintU3)*((6*Alfa2*MW)/SW2*
     -     (AbbSumU53*EL*mass2(1)*q(3))/SW + (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 + 1/4.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/(dup2 - MZ2) + tmp2/dup2 + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams N attached to the upper quark line
*** for currents 2 and 4

      subroutine boxesBuN_24(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup1, dup2, tmp1, tmp2

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

      dup1 = 1/(2.D0*CW*SW2)*(SW2*mass2(1) + SW2*mass2(2) + CW2*(4*
     -     mass2(2) - 2*MH2) + 2*SW2*MH2)
      dup2 = 2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U

      tmp1 = (Cval(cc0,iintU1) + Cval(cc0,iintU2))*
     -     ((Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb4*EL*gL(3)) + 
     -     (Alfa2*MW*(4*CW2 + SW2))/(2.D0*CW*SW2**2)*(Abb3*EL*gR(3)))+
     -     Dval(dd0,iintU3)*
     -     ((Alfa2*MW)/SW2*(Abb4*dup1*EL*gL(3)) + 
     -     (Alfa2*MW)/SW2*(Abb3*dup1*EL*gR(3)) - 
     -     (6*Alfa2*CW*MW)/SW2**2*(AbbSumU64*EL*gL(3)*mass2(1)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU63*EL*gR(3)*mass2(1)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*S34) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*S34) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gR(3)*(Abb62 - Abb3*T)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gL(3)*(Abb82 - Abb4*T))) + 
     -     1/SW2**2*((Dval(dd00,iintU3) + Dval(dd00,iintU4))*
     -     (4*Alfa2*CW*MW*(Abb4*EL*gL(3)) + 
     -     4*Alfa2*CW*MW*(Abb3*EL*gR(3)))) + 
     -     Dval(dd0,iintU4)*((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb82*EL*gL(3)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb62*EL*gR(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb2*EL*gL(3)*mass2(1)) - 
     -     (7*Alfa2*CW*MW)/SW2**2*(AbbSumU29*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(3)*mass2(1)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU27*EL*gR(3)*mass2(1)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*(S34 + T)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*(S34 + T)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*mass2(2) + 
     -     CW2*(4*mass2(2) - 4*MH2))) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*mass2(2) + 
     -     CW2*(4*mass2(2) - 4*MH2))))

      tmp1 = tmp1 + Dval(dd1,iintU4)*
     -     (-(Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU11*EL*gR(3)*mass2(1)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(EL*gR(3)*(Abb62 + Abb3*T14)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(EL*gL(3)*(Abb82 + Abb4*T14))) + 
     -     Dval(dd1,iintU3)*((Alfa2*MW)/(CW*SW2**2)*(Abb4*EL*gL(3)*MH2) + 
     -     (Alfa2*MW)/(CW*SW2**2)*(Abb3*EL*gR(3)*MH2) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(3)*(Abb62 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14))) - 
     -     (4*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(Abb82 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))) + 
     -     1/SW2**2*(Dval(dd11,iintU3)*
     -     (2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3)) + 
     -     2*Alfa2*CW*MW*(AbbSumU37*EL*gR(3)))) + 
     -     1/SW2**2*(Dval(dd13,iintU3)*
     -     (-(2*Alfa2*CW*MW*(Abb82*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(Abb62*EL*gR(3)))) + 1/SW2**2*
     -     (Dval(dd12,iintU3)*
     -     (2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3)) + 
     -     2*Alfa2*CW*MW*(AbbSumU37*EL*gR(3)))) + 1/SW2**2*
     -     (Dval(dd12,iintU4)*
     -     (-(2*Alfa2*CW*MW*(AbbSumU39*EL*gL(3))) - 
     -     2*Alfa2*CW*MW*(AbbSumU37*EL*gR(3)))) 

      tmp1 = tmp1 + 
     -     Dval(dd2,iintU4)*((Alfa2*MW)/(CW*SW2)*(Abb82*EL*gL(3)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb62*EL*gR(3)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU11*EL*gR(3)*mass2(1)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(3)*(AbbSumU41 - Abb3*S34)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumU43 - Abb4*S34)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*T14) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*T14) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*MH2 + CW2*(mass2(1) + MH2))) - 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*MH2 + CW2*(mass2(1) + MH2)))) + 
     -     Dval(dd2,iintU3)*((Alfa2*MW)/(CW*SW2)*(Abb45*EL*gL(3)) + 
     -     (Alfa2*MW)/(CW*SW2)*(Abb41*EL*gR(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU12*EL*gL(3)*mass2(1)) - 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(AbbSumU11*EL*gR(3)*mass2(1)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*S34) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*S34) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (EL*gR(3)*(AbbSumU72 + 2*Abb3*T14)) - 
     -     (Alfa2*CW*MW)/SW2**2*
     -     (EL*gL(3)*(AbbSumU74 + 2*Abb4*T14)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*MH2 + CW2*(mass2(1) + MH2))) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*MH2 + CW2*(mass2(1) + MH2)))) + 
     -     1/SW2**2*(Dval(dd13,iintU4)*
     -     (2*Alfa2*CW*MW*(Abb82*EL*gL(3)) + 
     -     2*Alfa2*CW*MW*(Abb62*EL*gR(3))))

      tmp1 = tmp1 + Dval(dd3,iintU4)*
     -     ((Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb82*EL*gL(3)) + 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*(Abb62*EL*gR(3)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb2*EL*gL(3)*mass2(1)) - 
     -     (7*Alfa2*CW*MW)/SW2**2*(AbbSumU29*EL*gL(3)*mass2(1)) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2)*(Abb1*EL*gR(3)*mass2(1)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU27*EL*gR(3)*mass2(1)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*MH2) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*MH2) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb4*EL*gL(3)*(S34 + T)) + 
     -     (2*Alfa2*CW*MW)/SW2**2*(Abb3*EL*gR(3)*(S34 + T))) + 
     -     Dval(dd3,iintU3)*(-(6*Alfa2*CW*MW)/SW2**2*(AbbSumU64*EL*
     -     gL(3)*mass2(1)) - 
     -     (3*Alfa2*CW*MW)/SW2**2*(AbbSumU63*EL*gR(3)*mass2(1)) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb4*EL*gL(3)*S34) - 
     -     (Alfa2*MW*(-CW2 + SW2))/(CW*SW2**2)*
     -     (Abb3*EL*gR(3)*S34) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gR(3)*(Abb62 - Abb3*T)) - 
     -     (2*Alfa2*CW*MW)/SW2**2*(EL*gL(3)*(Abb82 - Abb4*T)) + 
     -     (Alfa2*MW)/(CW*SW2**2)*
     -     (0*EL*mass(1)*mass(3)*
     -     (-(SW2*gL(3)) + CW2*(2*gL(3) - 2*gR(3)) + SW2*gR(3))) + 
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb4*EL*gL(3)*(SW2*mass2(1) - 2*CW2*MH2 + 2*SW2*MH2))+
     -     (Alfa2*MW)/(2.D0*CW*SW2**2)*
     -     (Abb3*EL*gR(3)*(SW2*mass2(1) - 2*CW2*MH2 + 2*SW2*MH2)))

      tmp2 = Dval(dd2,iintU4)*((Alfa2*MW)/(2.D0*SW2)*
     -     ((-2*AbbSumU1 + AbbSumU5)*EL*mass2(1)*q(3))/SW - 
     -     (Alfa2*MW)/SW2*(EL*q(3)*(2*AbbSumU22 + AbbSumU4 + 
     -     AbbSumU1*(-2*S34 + 2*T14)))/SW) + (3*Alfa2*MW)/(2.D0*SW2)*
     -     (AbbSumU1*EL*q(3)*(Cval(cc0,iintU1) + Cval(cc0,iintU2)))/SW + 
     -     (4*Alfa2*MW)/SW2*(AbbSumU1*EL*q(3)*
     -     (Dval(dd00,iintU3) + Dval(dd00,iintU4)))/SW + 
     -     (2*Alfa2*MW)/SW2*(AbbSumU20*EL*q(3)*Dval(dd11,iintU3))/SW + 
     -     Dval(dd2,iintU3)*(-((Alfa2*MW)/(2.D0*SW2)*
     -     ((-2*AbbSumU1 + AbbSumU5)*EL*mass2(1)*q(3))/SW) - 
     -     (Alfa2*MW)/SW2*(EL*q(3)*(AbbSumU50 + AbbSumU60 + 
     -     AbbSumU1*(-2*S34 + 2*T14)))/SW) + 
     -     Dval(dd0,iintU4)*(-(Alfa2*MW)/(2.D0*SW2)*
     -     ((14*AbbSumU16 + AbbSumU3)*EL*mass2(1)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 - 3/4.D0*mass2(2)))/SW) + 
     -     Dval(dd0,iintU3)*(-((6*Alfa2*MW)/SW2*
     -     (AbbSumU55*EL*mass2(1)*q(3))/SW) - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 + 1/4.D0*mass2(1) - 3/4.D0*mass2(2)))/
     -     SW) + 1/SW2*(Dval(dd13,iintU4)*
     -     (2*Alfa2*MW*(AbbSumU4*EL*q(3))))/SW + 
     -     1/SW2*(Dval(dd12,iintU4)*
     -     (-(2*Alfa2*MW*(AbbSumU20*EL*q(3)))))/SW + 
     -     1/SW2*(Dval(dd12,iintU3)*
     -     (2*Alfa2*MW*(AbbSumU20*EL*q(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintU3)*
     -     (-(2*Alfa2*MW*(AbbSumU4*EL*q(3)))))/SW + 
     -     1/(2.D0*SW2)*(Dval(dd1,iintU4)*
     -     (-2*Alfa2*MW*(AbbSumU1*EL*mass2(1)*q(3)) + 
     -     Alfa2*MW*(AbbSumU5*EL*mass2(1)*q(3))))/SW + 
     -     Dval(dd1,iintU3)*(-(4*Alfa2*MW)/SW2*(EL*q(3)*(AbbSumU4 + 
     -     AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW)

      tmp2 = tmp2 + Dval(dd3,iintU4)*(-((Alfa2*MW)/(2.D0*SW2)*
     -     ((14*AbbSumU16 + AbbSumU3)*EL*mass2(1)*q(3))/SW) - 
     -     (2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW) + 
     -     Dval(dd3,iintU3)*(-((6*Alfa2*MW)/SW2*
     -     (AbbSumU55*EL*mass2(1)*q(3))/SW) - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(3)*(AbbSumU4 + AbbSumU1*(-S34 - T)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(3)*(MH2 + 1/4.D0*mass2(1)))/SW)

      Cloop(1) = tmp1/(dup2 - MZ2) + tmp2/dup2 + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams P attached to the upper quark line
*** for currents 1 and 3

      subroutine boxesBuP_13(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup10, dup3, dup4, dup5, dup6, dup7, dup8, dup9
      double complex tmp3, tmp4


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

      dup3 = -(gL(3)*gR(1)) + gL(1)*(gL(3) - gR(3)) + gR(1)*gR(3)
      dup4 = gR(2)*mass2(2) + gL(2)*(-MH2 - MW2)
      dup5 = gR(1)*mass2(1) + gL(1)*(-MH2 - MZ2)
      dup6 = -(gL(1)*mass2(1)) + gR(1)*(MH2 + MZ2)
      dup7 = 2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U
      dup8 = 2*Abb41 + Abb3*(-S34 + T + T14)
      dup9 = 2*Abb45 + Abb4*(-S34 + T + T14)
      dup10 = 2*AbbSumU50 + AbbSumU1*(-S34 + T + T14)

      tmp3 = Dval(dd0,iintU8)*
     -     ((4*Alfa2*MW)/CW2*(Abb4*dup5*EL*gL(1)**2*gL(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup6*EL*gL(3)*gR(1)**2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup5*EL*gL(1)**2*gR(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup6*EL*gR(1)**2*gR(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU104*EL*gL(1)**3*gL(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU98*EL*gL(3)*gR(1)**3*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU103*EL*gL(1)**3*gR(3)*mass2(1))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU96*EL*gR(1)**3*gR(3)*mass2(1))/SW + (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(Abb24 + Abb1*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(Abb33 + Abb2*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gR(3)*(Abb41 + Abb3*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     Dval(dd0,iintU7)*((2*Alfa2*MW)/SW2*(Abb4*dup4*EL*gL(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb3*dup4*EL*gR(3))/SW - 
     -     (6*Alfa2*MW)/SW2*(AbbSumU104*EL*gL(2)*gL(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumU103*EL*gL(2)*gR(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(T + T14)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     1/SW2*(Cval(cc0,iintU5)*
     -     (-(2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/SW2*(Dval(dd00,iintU7)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Cval(cc0,iintU6)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3))) - 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) - 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd00,iintU8)*
     -     (8*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)) + 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) + 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) + 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW + 
     -     Dval(dd1,iintU8)*(- (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*
     -     MH2)/SW - (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)**3*gR(3)*
     -     (Abb24 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)**3*
     -     (Abb33 + Abb2*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**3*gR(3)*
     -     (Abb41 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**3*gL(3)*
     -     (Abb45 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW) + 
     -     Dval(dd1,iintU7)*(-(2*Alfa2*MW)/SW2*(Abb4*EL*gL(2)*gL(3)* 
     -     MH2)/SW - (2*Alfa2*MW)/SW2*(Abb3*EL*gL(2)*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW + (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW) + 1/SW2*(Dval(dd11,iintU7)*
     -     (2*Alfa2*MW*(AbbSumU39*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(AbbSumU37*EL*gL(2)*gR(3))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd11,iintU8)*
     -     (4*Alfa2*MW*(AbbSumU39*EL*gL(1)**3*gL(3)) + 
     -     4*Alfa2*MW*(AbbSumU115*EL*gL(3)*gR(1)**3) + 
     -     4*Alfa2*MW*(AbbSumU37*EL*gL(1)**3*gR(3)) + 
     -     4*Alfa2*MW*(AbbSumU114*EL*gR(1)**3*gR(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (4*Alfa2*MW*(AbbSumU121*EL*gL(2)*gL(3)) + 
     -     4*Alfa2*MW*(AbbSumU120*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintU8)*
     -     (8*Alfa2*MW*(AbbSumU121*EL*gL(1)**3*gL(3)) + 
     -     8*Alfa2*MW*(AbbSumU119*EL*gL(3)*gR(1)**3) + 
     -     8*Alfa2*MW*(AbbSumU120*EL*gL(1)**3*gR(3)) + 
     -     8*Alfa2*MW*(AbbSumU117*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + 1/SW2*(Dval(dd2,iintU7)*
     -     (-6*Alfa2*MW*(AbbSumU104*EL*gL(2)*gL(3)*mass2(1)) - 
     -     2*Alfa2*MW*(AbbSumU103*EL*gL(2)*gR(3)*mass2(1)) - 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*MH2) - 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*MH2) + 
     -     Alfa2*(2*MW*(dup9*EL*gL(2)*gL(3)) + 
     -     2*MW*(dup8*EL*gL(2)*gR(3)))))/SW + 
     -     1/CW2*(Dval(dd13,iintU8)*
     -     (4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3)) + 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) + 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) + 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd2,iintU8)*
     -     ((4*Alfa2*MW)/CW2*(dup9*EL*gL(1)**3*gL(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(dup8*EL*gL(1)**3*gR(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU104*EL*gL(1)**3*gL(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU98*EL*gL(3)*gR(1)**3*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU103*EL*gL(1)**3*gR(3)*mass2(1))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU96*EL*gR(1)**3*gR(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(2*Abb24 + Abb1*(-S34 + T + T14)))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(2*Abb33 + Abb2*(-S34 + T + T14)))/SW)+
     -     1/SW2*(Dval(dd22,iintU7)*
     -     (2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd3,iintU7)*
     -     (-(8*Alfa2*MW)/SW2*(AbbSumU132*EL*gL(2)*gL(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumU131*EL*gL(2)*gR(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb4*EL*gL(2)*gL(3)*MH2)/SW - 
     -     (2*Alfa2*MW)/SW2*(Abb3*EL*gL(2)*gR(3)*MH2)/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(T + T14)))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     1/SW2*(Dval(dd23,iintU7)*
     -     (2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Dval(dd22,iintU8)*
     -     (4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3)) + 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) + 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) + 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintU8)*
     -     (4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3)) + 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) + 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) + 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd3,iintU8)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumU132*EL*gL(1)**3*gL(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU128*EL*gL(3)*gR(1)**3*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU131*EL*gL(1)**3*gR(3)*mass2(1))/SW - 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumU126*EL*gR(1)**3*gR(3)*mass2(1))/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(Abb24 + Abb1*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(Abb33 + Abb2*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gR(3)*(Abb41 + Abb3*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) 

      tmp4 = Dval(dd0,iintU7)*
     -     (-(6*Alfa2*MW)/SW2*(AbbSumU79*EL*mass2(1)*q(2)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*(-mass2(2) + MH2 + MW2)*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW) + 
     -     Dval(dd0,iintU8)*(-(12*Alfa2*MW)/CW2*
     -     (AbbSumU79*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU75*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU1*EL*gL(1)**2*(-mass2(1) + MH2 + MZ2)*q(1)*q(3))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU3*EL*gR(1)**2*(-mass2(1) + MH2 + MZ2)*q(1)*q(3))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*q(1)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*(AbbSumU77 + AbbSumU3*(T + T14)))/SW)-
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(2)*q(3)*Cval(cc0,iintU5))/SW + 
     -     (4*Alfa2*MW)/SW2*(AbbSumU1*EL*q(2)*q(3)*Dval(dd00,iintU7))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU23*EL*q(2)*q(3)*Dval(dd11,iintU7))/SW + 
     -     1/CW2*(Cval(cc0,iintU6)*
     -     (-(4*Alfa2*MW*(AbbSumU1*EL*gL(1)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumU3*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd00,iintU8)*
     -     (8*Alfa2*MW*(AbbSumU1*EL*gL(1)**2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumU3*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd11,iintU8)*
     -     (4*Alfa2*MW*(AbbSumU23*EL*gL(1)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumU84*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (4*Alfa2*MW*(AbbSumU87*EL*q(2)*q(3))))/SW + 
     -     Dval(dd1,iintU7)*(-(2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(2)* 
     -     q(3))/SW + (4*Alfa2*MW)/SW2*(EL*q(2)*q(3)*(AbbSumU50 + 
     -     AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW) + 
     -     Dval(dd1,iintU8)*(-(4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2* 
     -     MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW + 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**2*q(1)*q(3)*
     -     (AbbSumU50 + AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW+
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)**2*q(1)*q(3)*
     -     (AbbSumU77 + AbbSumU3*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW)

      tmp4 = tmp4 + Dval(dd2,iintU8)*
     -     ((4*Alfa2*MW)/CW2*(dup10*EL*gL(1)**2*q(1)*q(3))/SW - 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU79*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU75*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*
     -     (2*AbbSumU77 + AbbSumU3*(-S34 + T + T14)))/SW) + 
     -     1/CW2*(Dval(dd13,iintU8)*
     -     (4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumU77*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintU8)*
     -     (8*Alfa2*MW*(AbbSumU87*EL*gL(1)**2*q(1)*q(3)) + 
     -     8*Alfa2*MW*(AbbSumU86*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd22,iintU8)*
     -     (4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumU77*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd23,iintU7)*
     -     (2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3))))/SW + 
     -     1/SW2*(Dval(dd22,iintU7)*
     -     (2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3))))/SW + 
     -     1/SW2*(Dval(dd2,iintU7)*
     -     (2*Alfa2*MW*(dup10*EL*q(2)*q(3)) - 
     -     6*Alfa2*MW*(AbbSumU79*EL*mass2(1)*q(2)*q(3)) - 
     -     2*Alfa2*MW*(AbbSumU1*EL*MH2*q(2)*q(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintU8)*
     -     (4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumU77*EL*gR(1)**2*q(1)*q(3))))/SW

      tmp4 = tmp4 + Dval(dd3,iintU7)*
     -     (-(8*Alfa2*MW)/SW2*(AbbSumU94*EL*mass2(1)*q(2)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW) + 
     -     Dval(dd3,iintU8)*(-(16*Alfa2*MW)/CW2*
     -     (AbbSumU94*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU92*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*q(1)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*(AbbSumU77 + AbbSumU3*(T + T14)))/SW)+
     -     (2*Alfa2*MW)/SW2*
     -     (0*EL*mass(1)*q(2)*q(3)*Dval(dd33,iintU7))/SW 

      Cloop(1) = tmp3/(dup7 - MZ2) + tmp4/dup7 + Cloop(1)


      end


********************************************************************************
********************************************************************************

*** This subroutine calculates box diagrams P attached to the upper quark line
*** for currents 2 and 4

      subroutine boxesBuP_24(mass,mass2,S,T,T24,T14,S34,U,q,gL,gR,Cloop)

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
      double complex dup10, dup3, dup4, dup5, dup6, dup7, dup8, dup9
      double complex tmp3, tmp4


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

      dup3 = -(gL(3)*gR(1)) + gR(1)*gR(3) + gL(1)*(-gL(3) + gR(3))
      dup4 = gR(2)*mass2(2) + gL(2)*(-MH2 - MW2)
      dup5 = gR(1)*mass2(1) + gL(1)*(-MH2 - MZ2)
      dup6 = -(gL(1)*mass2(1)) + gR(1)*(MH2 + MZ2)
      dup7 = 2*mass2(1) + 3*mass2(3) + MH2 - S - T24 - U
      dup8 = 2*Abb41 + Abb3*(-S34 + T + T14)
      dup9 = 2*Abb45 + Abb4*(-S34 + T + T14)
      dup10 = 2*AbbSumU50 + AbbSumU1*(-S34 + T + T14)

      tmp3 = Dval(dd0,iintU8)*
     -     (-((4*Alfa2*MW)/CW2*(Abb4*dup5*EL*gL(1)**2*gL(3))/SW) + 
     -     (4*Alfa2*MW)/CW2*(Abb2*dup6*EL*gL(3)*gR(1)**2)/SW - 
     -     (4*Alfa2*MW)/CW2*(Abb3*dup5*EL*gL(1)**2*gR(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*dup6*EL*gR(1)**2*gR(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU107*EL*gL(1)**3*gL(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU100*EL*gL(3)*gR(1)**3*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU105*EL*gL(1)**3*gR(3)*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU99*EL*gR(1)**3*gR(3)*mass2(1))/SW -
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(Abb24 + Abb1*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(Abb33 + Abb2*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gR(3)*(Abb41 + Abb3*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     Dval(dd0,iintU7)*(-((2*Alfa2*MW)/SW2*
     -     (Abb4*dup4*EL*gL(3))/SW) - 
     -     (2*Alfa2*MW)/SW2*(Abb3*dup4*EL*gR(3))/SW + 
     -     (6*Alfa2*MW)/SW2*(AbbSumU107*EL*gL(2)*gL(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(AbbSumU105*EL*gL(2)*gR(3)*mass2(1))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(T + T14)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     1/SW2*(Cval(cc0,iintU5)*
     -     (2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)) + 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/SW2*(Dval(dd00,iintU7)*
     -     (-(4*Alfa2*MW*(Abb4*EL*gL(2)*gL(3))) - 
     -     4*Alfa2*MW*(Abb3*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Cval(cc0,iintU6)*
     -     (4*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3)) + 
     -     4*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) + 
     -     4*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) + 
     -     4*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd00,iintU8)*
     -     (-(8*Alfa2*MW*(Abb4*EL*gL(1)**3*gL(3))) - 
     -     8*Alfa2*MW*(Abb2*EL*gL(3)*gR(1)**3) - 
     -     8*Alfa2*MW*(Abb3*EL*gL(1)**3*gR(3)) - 
     -     8*Alfa2*MW*(Abb1*EL*gR(1)**3*gR(3))))/SW + 
     -     Dval(dd1,iintU8)*((4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*
     -     MH2)/SW + (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)**3*gR(3)*
     -     (Abb24 + Abb1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gL(3)*gR(1)**3*
     -     (Abb33 + Abb2*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**3*gR(3)*
     -     (Abb41 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**3*gL(3)*
     -     (Abb45 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW) + 
     -     Dval(dd1,iintU7)*((2*Alfa2*MW)/SW2*(Abb4*EL*gL(2)*gL(3)*MH2)/
     -     SW + (2*Alfa2*MW)/SW2*(Abb3*EL*gL(2)*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW - (4*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(-(1/2.D0*S34) + 1/2.D0*T14)))/
     -     SW) + 1/SW2*(Dval(dd11,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU39*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(AbbSumU37*EL*gL(2)*gR(3))))/SW

      tmp3 = tmp3 + 1/CW2*(Dval(dd11,iintU8)*
     -     (-(4*Alfa2*MW*(AbbSumU39*EL*gL(1)**3*gL(3))) - 
     -     4*Alfa2*MW*(AbbSumU118*EL*gL(3)*gR(1)**3) - 
     -     4*Alfa2*MW*(AbbSumU37*EL*gL(1)**3*gR(3)) - 
     -     4*Alfa2*MW*(AbbSumU117*EL*gR(1)**3*gR(3))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (-(2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (-(4*Alfa2*MW*(AbbSumU124*EL*gL(2)*gL(3))) - 
     -     4*Alfa2*MW*(AbbSumU123*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintU8)*
     -     (-(8*Alfa2*MW*(AbbSumU124*EL*gL(1)**3*gL(3))) - 
     -     8*Alfa2*MW*(AbbSumU122*EL*gL(3)*gR(1)**3) - 
     -     8*Alfa2*MW*(AbbSumU123*EL*gL(1)**3*gR(3)) - 
     -     8*Alfa2*MW*(AbbSumU120*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + 1/SW2*(Dval(dd2,iintU7)*
     -     (6*Alfa2*MW*(AbbSumU107*EL*gL(2)*gL(3)*mass2(1)) + 
     -     2*Alfa2*MW*(AbbSumU105*EL*gL(2)*gR(3)*mass2(1)) + 
     -     2*Alfa2*MW*(Abb4*EL*gL(2)*gL(3)*MH2) + 
     -     2*Alfa2*MW*(Abb3*EL*gL(2)*gR(3)*MH2) + 
     -     Alfa2*(-(2*MW*(dup9*EL*gL(2)*gL(3))) - 
     -     2*MW*(dup8*EL*gL(2)*gR(3)))))/SW + 
     -     1/CW2*(Dval(dd13,iintU8)*
     -     (-(4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3))) - 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) - 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) - 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd2,iintU8)*
     -     (-((4*Alfa2*MW)/CW2*(dup9*EL*gL(1)**3*gL(3))/SW) - 
     -     (4*Alfa2*MW)/CW2*(dup8*EL*gL(1)**3*gR(3))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU107*EL*gL(1)**3*gL(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU100*EL*gL(3)*gR(1)**3*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU105*EL*gL(1)**3*gR(3)*mass2(1))/SW + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU99*EL*gR(1)**3*gR(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(2*Abb24 + Abb1*(-S34 + T + T14)))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(2*Abb33 + Abb2*(-S34 + T + T14)))/SW)+
     -     1/SW2*(Dval(dd22,iintU7)*
     -     (-(2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd3,iintU7)*
     -     ((8*Alfa2*MW)/SW2*(AbbSumU135*EL*gL(2)*gL(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(AbbSumU133*EL*gL(2)*gR(3)*mass2(1))/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb4*EL*gL(2)*gL(3)*MH2)/SW + 
     -     (2*Alfa2*MW)/SW2*(Abb3*EL*gL(2)*gR(3)*MH2)/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gR(3)*(Abb41 + Abb3*(T + T14)))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*gL(2)*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) + 
     -     1/SW2*(Dval(dd23,iintU7)*
     -     (-(2*Alfa2*MW*(Abb45*EL*gL(2)*gL(3))) - 
     -     2*Alfa2*MW*(Abb41*EL*gL(2)*gR(3))))/SW + 
     -     1/CW2*(Dval(dd22,iintU8)*
     -     (-(4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3))) - 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) - 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) - 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintU8)*
     -     (-(4*Alfa2*MW*(Abb45*EL*gL(1)**3*gL(3))) - 
     -     4*Alfa2*MW*(Abb33*EL*gL(3)*gR(1)**3) - 
     -     4*Alfa2*MW*(Abb41*EL*gL(1)**3*gR(3)) - 
     -     4*Alfa2*MW*(Abb24*EL*gR(1)**3*gR(3))))/SW

      tmp3 = tmp3 + Dval(dd3,iintU8)*((16*Alfa2*MW)/CW2*
     -     (AbbSumU135*EL*gL(1)**3*gL(3)*mass2(1))/SW +(4*Alfa2*MW)/CW2*
     -     (AbbSumU131*EL*gL(3)*gR(1)**3*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU133*EL*gL(1)**3*gR(3)*mass2(1))/SW + 
     -     (16*Alfa2*MW)/CW2*
     -     (AbbSumU130*EL*gR(1)**3*gR(3)*mass2(1))/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb4*EL*gL(1)**3*gL(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb2*EL*gL(3)*gR(1)**3*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb3*EL*gL(1)**3*gR(3)*MH2)/SW + 
     -     (4*Alfa2*MW)/CW2*(Abb1*EL*gR(1)**3*gR(3)*MH2)/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**3*gR(3)*(Abb24 + Abb1*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(3)*gR(1)**3*(Abb33 + Abb2*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gR(3)*(Abb41 + Abb3*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**3*gL(3)*(Abb45 + Abb4*(T + T14)))/SW) 

      tmp4 = Dval(dd0,iintU7)*
     -     ((6*Alfa2*MW)/SW2*(AbbSumU81*EL*mass2(1)*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*(-mass2(2) + MH2 + MW2)*q(2)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW) + 
     -     Dval(dd0,iintU8)*((12*Alfa2*MW)/CW2*
     -     (AbbSumU81*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU77*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU1*EL*gL(1)**2*(-mass2(1) + MH2 + MZ2)*q(1)*q(3))/SW+
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU3*EL*gR(1)**2*(-mass2(1) + MH2 + MZ2)*q(1)*q(3))/SW-
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*q(1)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*(AbbSumU79 + AbbSumU3*(T + T14)))/SW)+
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU1*EL*q(2)*q(3)*Cval(cc0,iintU5))/SW - 
     -     (4*Alfa2*MW)/SW2*(AbbSumU1*EL*q(2)*q(3)*Dval(dd00,iintU7))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (AbbSumU20*EL*q(2)*q(3)*Dval(dd11,iintU7))/SW + 
     -     1/CW2*(Cval(cc0,iintU6)*
     -     (4*Alfa2*MW*(AbbSumU1*EL*gL(1)**2*q(1)*q(3)) + 
     -     4*Alfa2*MW*(AbbSumU3*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd00,iintU8)*
     -     (-(8*Alfa2*MW*(AbbSumU1*EL*gL(1)**2*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumU3*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd11,iintU8)*
     -     (-(4*Alfa2*MW*(AbbSumU20*EL*gL(1)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumU86*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd12,iintU7)*
     -     (-(4*Alfa2*MW*(AbbSumU89*EL*q(2)*q(3)))))/SW + 
     -     Dval(dd1,iintU7)*((2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(2)*
     -     q(3))/SW - (4*Alfa2*MW)/SW2*(EL*q(2)*q(3)*(AbbSumU50 + 
     -     AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW) + 
     -     Dval(dd1,iintU8)*((4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2*MH2*
     -     q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW - 
     -     (8*Alfa2*MW)/CW2*(EL*gL(1)**2*q(1)*q(3)*
     -     (AbbSumU50 + AbbSumU1*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW-
     -     (8*Alfa2*MW)/CW2*(EL*gR(1)**2*q(1)*q(3)*
     -     (AbbSumU79 + AbbSumU3*(-(1/2.D0*S34) + 1/2.D0*T14)))/SW)

      tmp4 = tmp4 + Dval(dd2,iintU8)*
     -     (-((4*Alfa2*MW)/CW2*(dup10*EL*gL(1)**2*q(1)*q(3))/SW) + 
     -     (12*Alfa2*MW)/CW2*
     -     (AbbSumU81*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU77*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*
     -     (2*AbbSumU79 + AbbSumU3*(-S34 + T + T14)))/SW) + 
     -     1/CW2*(Dval(dd13,iintU8)*
     -     (-(4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumU79*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd12,iintU8)*
     -     (-(8*Alfa2*MW*(AbbSumU89*EL*gL(1)**2*q(1)*q(3))) - 
     -     8*Alfa2*MW*(AbbSumU88*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/CW2*(Dval(dd22,iintU8)*
     -     (-(4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumU79*EL*gR(1)**2*q(1)*q(3))))/SW + 
     -     1/SW2*(Dval(dd23,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3)))))/SW + 
     -     1/SW2*(Dval(dd13,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3)))))/SW + 
     -     1/SW2*(Dval(dd22,iintU7)*
     -     (-(2*Alfa2*MW*(AbbSumU50*EL*q(2)*q(3)))))/SW + 
     -     1/SW2*(Dval(dd2,iintU7)*
     -     (-(2*Alfa2*MW*(dup10*EL*q(2)*q(3))) + 
     -     6*Alfa2*MW*(AbbSumU81*EL*mass2(1)*q(2)*q(3)) + 
     -     2*Alfa2*MW*(AbbSumU1*EL*MH2*q(2)*q(3))))/SW + 
     -     1/CW2*(Dval(dd23,iintU8)*
     -     (-(4*Alfa2*MW*(AbbSumU50*EL*gL(1)**2*q(1)*q(3))) - 
     -     4*Alfa2*MW*(AbbSumU79*EL*gR(1)**2*q(1)*q(3))))/SW

      tmp4 = tmp4 + Dval(dd3,iintU7)*
     -     ((8*Alfa2*MW)/SW2*(AbbSumU96*EL*mass2(1)*q(2)*q(3))/SW + 
     -     (2*Alfa2*MW)/SW2*(AbbSumU1*EL*MH2*q(2)*q(3))/SW - 
     -     (2*Alfa2*MW)/SW2*
     -     (EL*q(2)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW) + 
     -     Dval(dd3,iintU8)*((16*Alfa2*MW)/CW2*
     -     (AbbSumU96*EL*gL(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*
     -     (AbbSumU94*EL*gR(1)**2*mass2(1)*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumU1*EL*gL(1)**2*MH2*q(1)*q(3))/SW + 
     -     (4*Alfa2*MW)/CW2*(AbbSumU3*EL*gR(1)**2*MH2*q(1)*q(3))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gL(1)**2*q(1)*q(3)*(AbbSumU50 + AbbSumU1*(T + T14)))/SW - 
     -     (4*Alfa2*MW)/CW2*
     -     (EL*gR(1)**2*q(1)*q(3)*(AbbSumU79 + AbbSumU3*(T + T14)))/SW)

      Cloop(1) = tmp3/(dup7 - MZ2) + tmp4/dup7 + Cloop(1)


      end


********************************************************************************
