c************** Calling the Fa functions*************************
       subroutine HexBox2Fa5(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa5)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa5(77:95)
       COMMON/HexBox2FaFunctions/Fa
       Fa5(77) = p4mup6*F(256)+p5mup6*F(257)+p6mup6*F(258)+p1mup6*F(25
     &   9)
       Fa5(78) = p5mup6*F(260)+p6mup6*F(261)+p4mup6*F(262)-p1mup6*F(26
     &   3)
       Fa5(79) = p5mup6*F(264)-p6mup6*F(265)+2*p4mup6*F(266)+p1mup6*F(
     &   267)
       Fa5(80) = p4mup6*F(230)-p1mup6*F(273)+p5mup6*F(274)+p6mup6*F(16
     &   39)
       Fa5(81) = p5mup6*F(181)+p4mup6*F(275)-p6mup6*F(276)-p1mup6*F(27
     &   7)
       Fa5(82) = p4mup6*F(279)+p5mup6*F(280)+p6mup6*F(281)+p1mup6*F(28
     &   2)
       Fa5(83) = p2mup6*F(289)+p1mup6*F(290)+p6mup6*F(291)+4*p5mup6*F(
     &   292)
       Fa5(84) = p1mup6*p2mup3*F(293)-p1mup3*p2mup6*F(294)-p2mup3*p5mu
     &   p6*F(295)-p2mup3*p6mup6*F(296)-p2mup6*p6mup3*F(297)+p1mup3*p
     &   5mup6*F(298)+p1mup3*p6mup6*F(299)+p1mup6*p6mup3*F(300)+p2mup
     &   3*p2mup6*F(301)-p5mup6*p6mup3*F(302)-p6mup3*p6mup6*F(303)-p1
     &   mup3*p1mup6*F(304)
       Fa5(85) = p2mup3*F(307)-p1mup3*F(308)-p6mup3*F(309)
       Fa5(86) = p2mup6*F(310)+p1mup6*F(311)+p5mup6*F(312)+p6mup6*F(31
     &   3)
       Fa5(87) = -(p1mup6*p6mup3*F(314))+p5mup6*p6mup3*F(315)+p6mup3*p
     &   6mup6*F(316)+p2mup6*p6mup3*F(317)-p1mup3*p5mup6*F(318)+p1mup
     &   3*p2mup6*F(319)+p1mup3*p6mup6*F(320)+p1mup3*p1mup6*F(321)
       Fa5(88) = -(p2mup6*F(322))+p5mup6*F(323)+p6mup6*F(324)+p1mup6*F
     &   (325)
       Fa5(89) = p5mup6*F(326)+p6mup6*F(327)+p2mup6*F(328)+p1mup6*F(32
     &   9)
       Fa5(90) = p2mup6*F(330)+p1mup6*F(331)+p5mup6*F(332)+p6mup6*F(33
     &   3)
       Fa5(91) = p2mup6*F(336)+2*(p5mup6*F(334)+p6mup6*F(335)+p1mup6*F
     &   (337))
       Fa5(92) = p2mup6*F(297)-p1mup6*F(300)+p5mup6*F(302)+p6mup6*F(30
     &   3)
       Fa5(93) = p1mup3*F(339)+p6mup3*F(340)
       Fa5(94) = p6mup3*F(343)-8*(p2mup3*F(342)+p5mup3*F(344))
       Fa5(95) = -8*(p2mup3*F(345)+p6mup3*F(346)-p5mup3*F(347))+p1mup3
     &   *F(348)
       Return
       End
