c************** Calling the Fa functions*************************
       subroutine HexBoxFa10(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa10)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa10(109:120)
       COMMON/HexBoxFaFunctions/Fa
       Fa10(109) = 2*(p1mup3*F(387)+p2mup3*F(388)+p6mup3*F(389))+(p4mu
     &   p3+p5mup3)*F(422)
       Fa10(110) = 2*(p4mup4+p5mup4)*(-(p2mup3*F(291))-p1mup3*F(292)+(
     &   p4mup3+p5mup3)*F(293)+p6mup3*F(294))+mup3mup4*s45*F(396)
       Fa10(111) = -2*mup3mup4*F(397)+(p4mup4+p5mup4)*(4*(p6mup3*F(281
     &   )+p1mup3*F(290)+p2mup3*F(295))+(p4mup3+p5mup3)*F(418))
       Fa10(112) = (p4mup4+p5mup4)*F(289)-p2mup4*F(425)-p6mup4*F(426)+
     &   p1mup4*F(427)
       Fa10(113) = p1mup3*F(428)+8*(p2mup3*F(429)-p6mup3*F(430)-(p4mup
     &   3+p5mup3)*F(431))
       Fa10(114) = -8*p2mup4*s45*F(288)+p6mup4*F(435)+(p4mup4+p5mup4)*
     &   F(436)
       Fa10(115) = p2mup3*F(437)+p1mup3*F(438)-8*((p4mup3+p5mup3)*F(43
     &   9)-p6mup3*F(440))
       Fa10(116) = p2mup4*F(267)+p6mup4*F(441)
       Fa10(117) = p2mup4*F(343)-p1mup4*F(441)
       Fa10(118) = p2mup3*F(274)+p1mup3*F(276)-(p4mup3+p5mup3)*F(277)+
     &   p6mup3*F(283)
       Fa10(119) = p2mup3*F(291)+p1mup3*F(292)-(p4mup3+p5mup3)*F(293)-
     &   p6mup3*F(294)
       Fa10(120) = p1mup4*F(288)+p6mup4*F(334)
       Return
       End
