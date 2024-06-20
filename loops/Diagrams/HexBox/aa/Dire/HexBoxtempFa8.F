c************** Calling the Fa functions*************************
       subroutine HexBoxFa8(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa8)
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
       Complex*16 Fa8(85:96)
       COMMON/HexBoxFaFunctions/Fa
       Fa8(85) = -4*p2mup4*s45*F(288)+p4mup4*F(355)+p5mup4*F(356)+p6mu
     &   p4*F(357)
       Fa8(86) = p1mup4*F(358)+p2mup4*F(359)-4*((p4mup4+p5mup4)*F(287)
     &   -p6mup4*F(360))
       Fa8(87) = -4*p2mup4*F(346)+(p4mup4+p5mup4)*F(361)+p1mup4*F(362)
       Fa8(88) = p1mup3*F(363)+p2mup3*F(364)-4*(p4mup3*F(365)+p5mup3*F
     &   (366)-p6mup3*F(367))
       Fa8(89) = p1mup3*F(368)+p2mup3*F(369)+p6mup3*F(370)+(p4mup3+p5m
     &   up3)*F(371)
       Fa8(90) = p1mup3*F(372)+p6mup3*F(373)+(p4mup3+p5mup3)*F(374)+p2
     &   mup3*F(375)
       Fa8(91) = -(mup3mup4*F(376))+8*(p4mup4+p5mup4)*((p4mup3+p5mup3)
     &   *F(329)+p1mup3*F(330)-p6mup3*F(332)+p2mup3*F(377))
       Fa8(92) = -4*mup3mup4*F(378)-8*(p4mup4*(p6mup3*F(281)+p1mup3*F(
     &   290)+p2mup3*F(295))+p5mup4*(p6mup3*F(281)+p1mup3*F(290)+p2mu
     &   p3*F(295))-p1mup4*((p4mup3+p5mup3)*F(329)+p1mup3*F(330)-p6mu
     &   p3*F(332)+p2mup3*F(377))-p4mup4*(p4mup3+p5mup3)*F(379)-(p4mu
     &   p3+p5mup3)*p5mup4*F(379)-p2mup4*(p1mup3*F(337)-(p4mup3+p5mup
     &   3)*F(338)+p6mup3*F(339)+p2mup3*F(380)))
       Fa8(93) = p1mup4*F(383)+p6mup4*F(384)-(p4mup4+p5mup4)*F(385)
       Fa8(94) = -(p1mup3*F(243))+p2mup3*F(253)-p6mup3*F(256)+(p4mup3+
     &   p5mup3)*F(386)
       Fa8(95) = 4*(p1mup3*F(387)+p2mup3*F(388)+p6mup3*F(389))+(p4mup3
     &   +p5mup3)*F(390)
       Fa8(96) = 2*(p4mup4+p5mup4)*(-(p1mup3*F(401))-p2mup3*F(402)+(p4
     &   mup3+p5mup3)*F(403)+p6mup3*F(1311))+mup3mup4*F(1308)*P(90)
       Return
       End
