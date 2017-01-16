c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa34(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa34)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa34(793:816)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa34(793) = p1mup6*F(9426)+p5mup6*F(9427)+p6mup6*F(9428)+p2mup6
     -   *F(9429)+p3mup6*F(9430)+p4mup6*F(9431)
       Fa34(794) = -2*(p3mup3*p6mup6*F(9222)-mup3mup6*F(9433)+p4mup3*p
     -   6mup6*F(9436)+p5mup3*p6mup6*F(9437)-p6mup3*p6mup6*F(9439))+p
     -   2mup3*p6mup6*F(9442)+p1mup3*p6mup6*F(9447)-8*(p3mup3*p5mup6*
     -   F(9220)+p5mup3*p5mup6*F(9434)+p4mup3*p5mup6*F(9435)-p5mup6*p
     -   6mup3*F(9438)-p2mup3*p5mup6*F(9441)-p4mup6*p5mup3*F(9443)-p1
     -   mup6*p5mup3*F(9445)-p1mup3*p5mup6*F(9446)-p2mup6*(p5mup3*F(9
     -   440)+2*(p4mup3*F(8651)+p6mup3*F(9618)-p1mup3*F(9628)-p2mup3*
     -   F(9635)-p3mup3*F(9637)))+p3mup6*(p5mup3*F(9444)+2*(p4mup3*F(
     -   7109)-p6mup3*F(9624)+p1mup3*F(9630)+p2mup3*F(9637)+p3mup3*F(
     -   9643))))+16*(p1mup6*p4mup3*F(8657)+p4mup3*p4mup6*F(8677)+p4m
     -   up6*p6mup3*F(9620)+p1mup6*p6mup3*F(9626)-p1mup6*p2mup3*F(962
     -   8)-p1mup6*p3mup3*F(9630)-p1mup3*p4mup6*F(9631)-p2mup3*p4mup6
     -   *F(9638)-p3mup3*p4mup6*F(9642)-p1mup3*p1mup6*F(9656))
       Fa34(795) = 16*(p2mup3*p3mup5*F(4087)+p2mup3*p4mup5*F(4089)-p2m
     -   up5*(p4mup3*F(7521)+p5mup3*F(7554))+p3mup5*p4mup3*F(7564)+p3
     -   mup5*p5mup3*F(7569)-p4mup3*p4mup5*F(8733)-p4mup5*p5mup3*F(87
     -   34)-p1mup5*p4mup3*F(8747)-p1mup5*p5mup3*F(8753))+p5mup5*p6mu
     -   p3*F(9450)+p2mup3*p5mup5*F(9453)+2*(p5mup3*p5mup5*F(9448)+p4
     -   mup3*p5mup5*F(9449)+p3mup3*(8*(p3mup5*F(4003)+p1mup5*F(4079)
     -   +p2mup5*F(4087)+p4mup5*F(4091))+p5mup5*F(9458)))+p3mup3*p6mu
     -   p5*F(9459)+8*(2*((p1mup5*p2mup3+p1mup3*p2mup5)*F(4076)+p2mup
     -   3*p2mup5*F(4085))+p6mup3*p6mup5*F(8331)+p4mup3*p6mup5*F(9451
     -   )+p5mup3*p6mup5*F(9452)+p2mup5*p6mup3*F(9454)-p2mup3*p6mup5*
     -   F(9455)+p4mup5*p6mup3*F(9456)+p3mup5*p6mup3*F(9457)+mup3mup5
     -   *F(9460)+p1mup5*p6mup3*F(9462))+p1mup3*(p5mup5*F(9461)-8*p6m
     -   up5*F(9463)+16*(p3mup5*F(4079)+p4mup5*F(5738)+p1mup5*F(9464)
     -   ))
       Fa34(796) = mup5mup6*(p4mup3*F(9466)+p5mup3*F(9467)-4*(p2mup3*F
     -   (9465)+p3mup3*F(9468)))
       Fa34(797) = 4*mup3mup5*(-(p6mup6*F(9473))+2*(p1mup6*F(9471)+p2m
     -   up6*F(9475)+p5mup6*F(9476)))+2*mup3mup6*(p5mup5*F(9472)+4*(p
     -   1mup5*F(9470)-p2mup5*F(9474)-p6mup5*F(9477)))
       Fa34(798) = 4*mup5mup6*p6mup3*F(9478)-8*mup3mup6*p4mup5*F(9479)
       Fa34(799) = p4mup6*(mup3mup5*F(9480)+2*p2mup5*(p4mup3*F(8677)+p
     -   5mup3*F(9481)-p3mup3*F(9642))+p2mup3*(-(p5mup5*F(5687))+2*(p
     -   2mup5*F(8735)+p3mup5*F(9646)+p4mup5*F(9648))))
       Fa34(800) = p4mup6*(p1mup3*p2mup5*F(9631)-p1mup5*p2mup3*F(9633)
     -   )
       Fa34(801) = p2mup6*(p4mup3*p5mup5*F(9483)+4*(p3mup3*p4mup5*F(40
     -   91)+p1mup3*p4mup5*F(5738)+p3mup5*p4mup3*F(7564)-p4mup3*p4mup
     -   5*F(8733)-p4mup5*p5mup3*F(8734)+p2mup3*p4mup5*F(9638)-p2mup5
     -   *p4mup3*F(9640)))+4*(p2mup5*p4mup3*(-(p3mup6*F(7109))+p5mup6
     -   *F(7232))+p2mup3*p4mup5*(p1mup6*F(9633)+p3mup6*F(9646)+p5mup
     -   6*F(9649)))
       Fa34(802) = p2mup6*(p2mup3*p5mup5*F(9484)+p1mup3*(p5mup5*F(9485
     -   )-4*(p3mup5*F(4079)+p1mup5*F(9464)-p2mup5*F(9629)))+4*(p1mup
     -   5*(p4mup3*F(8747)+p5mup3*F(8753))-p2mup3*(p1mup5*F(9628)+p3m
     -   up5*F(9637))-p3mup3*(p3mup5*F(4003)+p1mup5*F(4079)-p2mup5*F(
     -   9639))+p2mup5*p5mup3*F(9641)))+2*(2*p2mup5*(p1mup3*p3mup6*F(
     -   9630)+p1mup3*p5mup6*F(9632)+p3mup3*p3mup6*F(9643)+p3mup3*p5m
     -   up6*F(9644)+p1mup6*(-(p4mup3*F(8657))-p5mup3*F(9486)+p3mup3*
     -   F(9630)+p1mup3*F(9656)))+p2mup3*(p5mup5*(p5mup6*F(5688)+p1mu
     -   p6*F(9487))+p3mup6*(p5mup5*F(5684)-2*(p2mup5*F(8745)+p1mup5*
     -   F(9632)+p3mup5*F(9644)))-2*(p1mup6*p2mup5*F(8748)+p2mup5*p5m
     -   up6*F(9488)+p1mup6*p3mup5*F(9632)+p1mup5*p5mup6*F(9634)+p3mu
     -   p5*p5mup6*F(9647)+p1mup5*p1mup6*F(9657))))
       Fa34(803) = -4*(p2mup5*p3mup6*p5mup3*F(4061)-p2mup5*p5mup3*p5mu
     -   p6*F(7233)-p2mup6*p3mup5*p5mup3*F(7569))+p2mup6*p5mup3*p5mup
     -   5*F(9489)+p2mup6*p3mup3*p5mup5*F(9490)
       Fa34(804) = p6mup3*(-(p2mup6*p4mup5*F(9620))+p2mup5*p4mup6*F(96
     -   20)+p2mup5*p3mup6*F(9624)+(p1mup6*p2mup5-p1mup5*p2mup6)*F(96
     -   26))
       Fa34(805) = -4*p2mup5*p5mup6*p6mup3*F(9621)-p2mup6*(p5mup5*p6mu
     -   p3*F(5723)-4*(p3mup5*p6mup3*F(9624)-p1mup3*p6mup5*F(9627)))+
     -   p2mup3*(4*p6mup5*(p3mup6*F(9621)+p4mup6*F(9622)+p5mup6*F(962
     -   5)+p1mup6*F(9627))+p6mup6*(-(p5mup5*F(9400))-2*(p2mup5*F(430
     -   1)-p3mup5*F(9621)-p4mup5*F(9622)-p1mup5*F(9627))))
       Fa34(806) = p2mup6*p6mup5*(p3mup3*F(9621)+p4mup3*F(9622)+p5mup3
     -   *F(9625))
       Fa34(807) = p6mup6*(p2mup5*p5mup3*F(9491)+p2mup3*p6mup5*F(9617)
     -   )
       Fa34(808) = p2mup5*p6mup6*(p4mup3*F(9492)+p3mup3*F(9624)+p1mup3
     -   *F(9626))
       Fa34(809) = mup5mup6*(p2mup2*p2mup3*F(4110)-p2mup2*p4mup3*F(851
     -   1)+p2mup3*p4mup2*F(8512)+p3mup2*p4mup3*F(8519)-p4mup2*p4mup3
     -   *F(8662)-p5mup2*(p3mup3*F(5680)+p5mup3*F(8736))+p3mup2*p5mup
     -   3*F(8740))
       Fa34(810) = mup5mup6*(p3mup2*p3mup3*F(4324)+p1mup2*p3mup3*F(448
     -   0)+p3mup3*p4mup2*F(5732)+p1mup2*p4mup3*F(8568)+p1mup2*p1mup3
     -   *F(8570)-p1mup3*p4mup2*F(8573)+p1mup3*p3mup2*F(8670)-p4mup3*
     -   p5mup2*F(8724)+p2mup2*(p3mup3*F(4112)+p1mup3*F(8571)+p5mup3*
     -   F(8730))+p2mup3*(p3mup2*F(4322)-p1mup2*F(8566)-p5mup2*F(8731
     -   ))+p1mup2*p5mup3*F(8732)+p4mup2*p5mup3*F(8739)-p1mup3*p5mup2
     -   *F(8743)+mup2mup3*F(9493))
       Fa34(811) = mup5mup6*p6mup2*(p5mup3*F(7580)+p4mup3*F(8675))
       Fa34(812) = mup3mup6*(p4mup2*p5mup5*F(9495)+p5mup2*p5mup5*F(949
     -   7)-p3mup2*p5mup5*F(9499)+p2mup2*(-(p5mup5*F(9494))-4*(p2mup5
     -   *F(4109)-p1mup5*F(4122)-p6mup5*F(9619)))+4*(p2mup5*p3mup2*F(
     -   8176)-p2mup5*p5mup2*F(8194)-p2mup5*p4mup2*F(8757)-p1mup5*p4m
     -   up2*F(8760)-p1mup5*p5mup2*F(8761)+p1mup5*p3mup2*F(8762)+p3mu
     -   p2*p6mup5*F(9621)+p4mup2*p6mup5*F(9622)+p5mup2*p6mup5*F(9625
     -   )))+2*mup3mup5*(p2mup2*p6mup6*F(8377)+p4mup2*p6mup6*F(9496)+
     -   p5mup2*p6mup6*F(9498)+p3mup2*p6mup6*F(9618)-2*(p1mup6*p2mup2
     -   *F(8389)+p2mup2*p5mup6*F(9502)-p5mup2*p5mup6*F(9503)-p4mup2*
     -   p5mup6*F(9504)+mup2mup6*F(9505)-p1mup6*p4mup2*F(9506)-p1mup6
     -   *p5mup2*F(9507)+p1mup6*p3mup2*F(9628)-p2mup6*(p5mup2*F(4105)
     -   +p2mup2*F(9500)+p4mup2*F(9501)-p3mup2*F(9635))+p3mup2*p5mup6
     -   *F(9639)))
       Fa34(813) = 2*mup3mup5*p1mup2*(-(p6mup6*F(9509))+2*(p1mup6*F(86
     -   68)-p2mup6*F(9511)+p5mup6*F(9512)))+4*(mup5mup6*p6mup3*(p1mu
     -   p2*F(8589)+p2mup2*F(8590)-p5mup2*F(9510))-mup3mup5*p5mup6*p6
     -   mup2*F(9621))+mup3mup6*(-(p5mup5*p6mup2*F(5723))+p1mup2*p5mu
     -   p5*F(9508)+4*(p1mup2*p2mup5*F(8738)-p1mup2*p1mup5*F(8742)+p2
     -   mup5*p6mup2*F(9618)+p1mup5*p6mup2*F(9626)-p1mup2*p6mup5*F(96
     -   27)))
       Fa34(814) = 2*mup3mup6*(p2mup2*(p3mup5*F(4111)+p4mup5*F(4114))+
     -   p3mup5*(p5mup2*F(5681)+p4mup2*F(5733)+p3mup2*F(6075)+p1mup2*
     -   F(8729))+p4mup5*(p1mup2*F(8723)-p4mup2*F(8758)-p5mup2*F(8759
     -   )+p3mup2*F(8763))+mup2mup5*F(9505))+mup3mup5*p6mup2*(-(p6mup
     -   6*F(9616))+2*(p2mup6*F(9618)+p1mup6*F(9626)))
       Fa34(815) = mup3mup6*p6mup2*(p4mup5*F(9620)+p3mup5*F(9624))
       Fa34(816) = mup3mup5*(p3mup6*(p5mup2*F(4118)+p2mup2*F(5664)+p1m
     -   up2*F(8676)+p4mup2*F(8678)-p6mup2*F(9624)+p3mup2*F(9637))+p4
     -   mup6*(-(p5mup2*F(4116))+p1mup2*F(8673)+p2mup2*F(8674)-p4mup2
     -   *F(9513)-p6mup2*F(9620)+p3mup2*F(9638)))
       Return
       End
