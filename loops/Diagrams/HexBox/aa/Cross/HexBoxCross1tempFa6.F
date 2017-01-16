c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa6(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa6
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa6(66:78)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa6(66) = mup3mup4*F(255)+p1mup4*p4mup3*F(264)+p1mup4*p5mup3*F(
     &   264)+p4mup3*p6mup4*F(265)+p5mup3*p6mup4*F(265)+p1mup3*(p6mup
     &   4*F(260)-8*((p4mup4+p5mup4)*F(256)-p2mup4*F(261))+p1mup4*F(2
     &   67))-p2mup3*(-(p1mup4*F(262))+8*((p4mup4+p5mup4)*F(257)+p6mu
     &   p4*F(263)+p2mup4*F(269)))+p1mup4*p6mup3*s45*F(270)+8*(p4mup3
     &   *p4mup4*F(258)+p4mup4*p5mup3*F(258)+p4mup3*p5mup4*F(258)+p5m
     &   up3*p5mup4*F(258)+p4mup4*p6mup3*F(259)+p5mup4*p6mup3*F(259)-
     &   p2mup4*p4mup3*F(266)-p2mup4*p5mup3*F(266)+p6mup3*p6mup4*F(26
     &   8)+p2mup4*p6mup3*F(271))
       Fa6(67) = mup3mup4*F(273)+16*(p4mup4+p5mup4)*((p4mup3+p5mup3)*F
     &   (274)+p6mup3*F(275))
       Fa6(68) = (p4mup4+p5mup4)*F(277)+p1mup4*F(278)+8*(p6mup4*F(279)
     &   -p2mup4*F(280))
       Fa6(69) = p1mup4*F(282)+8*((p4mup4+p5mup4)*F(281)-p2mup4*F(283)
     &   -p6mup4*F(284))
       Fa6(70) = p1mup3*F(285)+8*((p4mup3+p5mup3)*F(286)+p6mup3*F(287)
     &   +p2mup3*F(288))
       Fa6(71) = p1mup3*F(289)+8*(p2mup3*F(290)-p6mup3*F(291)+(p4mup3+
     &   p5mup3)*F(292))
       Fa6(72) = (p4mup3+p5mup3)*F(294)+4*(p1mup3*F(293)+p2mup3*F(295)
     &   +p6mup3*F(296))
       Fa6(73) = p2mup3*F(297)+(p4mup3+p5mup3)*F(298)+p1mup3*F(299)+p6
     &   mup3*F(300)
       Fa6(74) = p1mup4*F(301)-p6mup4*F(302)-p2mup4*F(303)+(p4mup4+p5m
     &   up4)*F(304)
       Fa6(75) = p1mup3*F(305)+p6mup3*F(307)+8*((p4mup3+p5mup3)*F(306)
     &   +p2mup3*F(308))
       Fa6(76) = p1mup4*F(310)+p2mup4*F(312)+4*(p6mup4*F(311)-(p4mup4+
     &   p5mup4)*F(313))
       Fa6(77) = -8*(p4mup4+p5mup4)*F(316)+p1mup4*F(317)+p6mup4*F(318)
       Fa6(78) = (p4mup3+p5mup3)*F(319)+p1mup3*F(320)+p6mup3*F(321)+p2
     &   mup3*F(322)
       Return
       End
