c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa5(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa5
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa5(77:95)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa5(77) = p6mup6*F(259)+p5mup6*F(260)+p1mup6*F(261)+p4mup6*F(26
     &   2)
       Fa5(78) = p6mup6*F(263)+p5mup6*F(264)+p1mup6*F(265)+p4mup6*F(26
     &   6)
       Fa5(79) = -(p6mup6*F(267))+p1mup6*F(268)+p4mup6*F(269)+p5mup6*F
     &   (270)
       Fa5(80) = p1mup6*F(273)+p5mup6*F(274)+p6mup6*F(275)+8*p2mup6*F(
     &   276)
       Fa5(81) = -(p2mup3*(p5mup6*F(279)+p6mup6*F(280)-p1mup6*F(281)))
     &   +p2mup6*(p5mup3*F(277)+p6mup3*F(278)+p2mup3*F(282))+p5mup3*p
     &   6mup6*F(283)+p6mup3*p6mup6*F(284)+p5mup6*p6mup3*F(285)+2*(p5
     &   mup3*p5mup6*F(286)+p1mup6*p5mup3*F(287))
       Fa5(82) = p1mup6*F(289)+p5mup6*F(290)+p6mup6*F(291)+p2mup6*F(29
     &   2)
       Fa5(83) = -8*(p1mup3*F(293)+p5mup3*F(294))+p6mup3*F(295)
       Fa5(84) = p1mup3*F(296)+p5mup3*F(297)+p6mup3*F(298)+8*p2mup3*F(
     &   299)
       Fa5(85) = p1mup6*p6mup3*F(278)+p2mup6*p6mup3*F(300)+p5mup3*p6mu
     &   p6*F(301)+p6mup3*p6mup6*F(302)-p5mup6*p6mup3*F(303)+2*(p5mup
     &   3*p5mup6*F(304)-p2mup6*p5mup3*F(305)-p1mup6*p5mup3*F(306))
       Fa5(86) = -(p2mup3*(p1mup6*F(309)+p5mup6*F(310)+p6mup6*F(311)))
     &   +p1mup3*(p2mup6*F(308)+p5mup6*F(312)+p6mup6*F(313)+p1mup6*F(
     &   314))+p6mup3*(p1mup6*F(315)+p2mup6*F(316)-p6mup6*F(317)-p5mu
     &   p6*F(318))+p2mup3*p2mup6*F(319)
       Fa5(87) = -(p5mup3*p6mup6*F(317))-p5mup3*p5mup6*F(318)+p2mup3*p
     &   5mup6*F(322)+p2mup3*p6mup6*F(323)-p1mup3*p5mup6*F(324)+p1mup
     &   6*(p5mup3*F(315)+p2mup3*F(321)+p1mup3*F(325))+p1mup3*p6mup6*
     &   F(326)+p2mup6*(p5mup3*F(316)-p1mup3*F(320)+2*p2mup3*F(327))
       Fa5(88) = p6mup3*(p1mup6*F(331)+p6mup6*F(332)-p2mup6*F(333)-p5m
     &   up6*F(335))+p1mup3*(-(p6mup6*F(330))-p5mup6*F(334)+p2mup6*F(
     &   336)+p1mup6*F(337))
       Fa5(89) = p6mup3*p6mup6*F(323)+2*p2mup6*p6mup3*F(327)-p1mup3*p6
     &   mup6*F(338)+p5mup6*(p6mup3*F(322)-p1mup3*F(339)+p5mup3*F(340
     &   ))+p5mup3*p6mup6*F(341)+p2mup6*p5mup3*F(343)+p1mup3*p2mup6*F
     &   (344)+p1mup6*(p6mup3*F(321)+p5mup3*F(342)+p1mup3*F(345))
       Fa5(90) = p2mup6*p5mup3*F(356)+p2mup6*p6mup3*F(357)+p1mup3*(p1m
     &   up6*F(349)+p6mup6*F(350)+p2mup6*F(355)+p5mup6*F(358))+p5mup3
     &   *p5mup6*F(360)-8*(p5mup3*p6mup6*F(353)+p6mup3*p6mup6*F(354)+
     &   p5mup6*p6mup3*F(359)-p2mup3*p6mup6*F(361)-p2mup3*p5mup6*F(36
     &   2)-p2mup3*p2mup6*F(363))+p1mup6*(p5mup3*F(351)+p6mup3*F(352)
     &   +8*p2mup3*F(364))
       Fa5(91) = p6mup3*F(366)+8*(p2mup3*F(367)-p5mup3*F(368))
       Fa5(92) = p6mup3*F(369)+8*(p2mup3*F(370)+p1mup3*F(371))
       Fa5(93) = p1mup3*F(372)+p5mup3*F(373)-8*(p6mup3*F(374)-p2mup3*F
     &   (375))
       Fa5(94) = p1mup6*F(379)+p6mup6*F(380)+p2mup6*F(381)+p5mup6*F(38
     &   2)
       Fa5(95) = p1mup6*F(383)+p6mup6*F(384)+p2mup6*F(385)+p5mup6*F(38
     &   6)
       Return
       End
