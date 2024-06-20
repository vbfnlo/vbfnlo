c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa7(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa7
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
       Complex*16 Fa7(79:91)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa7(79) = -(p1mup3*F(261))+p2mup3*F(269)-p6mup3*F(271)+(p4mup3+
     &   p5mup3)*F(323)
       Fa7(80) = mup3mup4*F(327)-16*(p4mup4+p5mup4)*(p2mup3*F(328)+p6m
     &   up3*F(329)+(p4mup3+p5mup3)*F(330)+p1mup3*F(331))
       Fa7(81) = 4*mup3mup4*F(313)+8*(p1mup3*p6mup4*F(333)-p1mup3*p1mu
     &   p4*F(334)+p4mup3*p6mup4*F(335)+p5mup3*p6mup4*F(335)-p1mup4*p
     &   2mup3*F(336)+p6mup3*p6mup4*F(337)+p1mup4*p4mup3*F(338)+p1mup
     &   4*p5mup3*F(338)+p1mup4*p6mup3*F(339)+p2mup3*p6mup4*F(340))
       Fa7(82) = 4*mup3mup4*F(283)+8*(p4mup4+p5mup4)*(p1mup3*F(341)-(p
     &   4mup3+p5mup3)*F(342)+p2mup3*F(343)+p6mup3*F(344))
       Fa7(83) = mup3mup4*F(345)+16*(p1mup4*(p2mup3*F(328)+p6mup3*F(32
     &   9)+(p4mup3+p5mup3)*F(330)+p1mup3*F(331))-p4mup3*p6mup4*F(335
     &   )-p5mup3*p6mup4*F(335)-p6mup3*p6mup4*F(337)-p1mup3*p2mup4*F(
     &   341)+p2mup4*p4mup3*F(342)+p2mup4*p5mup3*F(342)-p2mup3*p2mup4
     &   *F(343)-p2mup4*p6mup3*F(344)-p1mup3*p6mup4*F(346)+p2mup3*p6m
     &   up4*F(347))
       Fa7(84) = mup3mup4*F(348)-16*(2*(p4mup3*p6mup4*F(335)+p5mup3*p6
     &   mup4*F(335)+p6mup3*p6mup4*F(337))+p1mup3*p2mup4*F(341)-p2mup
     &   4*p4mup3*F(342)-p2mup4*p5mup3*F(342)+p2mup3*p2mup4*F(343)+p2
     &   mup4*p6mup3*F(344)-p1mup4*p2mup3*F(349)+p1mup3*p6mup4*F(350)
     &   -p1mup4*p6mup3*F(351)-p2mup3*p6mup4*F(352)-p1mup4*(p1mup3+p4
     &   mup3+p5mup3)*F(353))
       Fa7(85) = -8*p2mup4*s45*F(274)+(p4mup4+p5mup4)*F(355)+p6mup4*F(
     &   356)
       Fa7(86) = -(mup3mup4*F(357))+16*(p5mup4*p6mup3*F(351)+(p1mup3+p
     &   4mup3+p5mup3)*p5mup4*F(353)+p4mup4*(p6mup3*F(351)+(p1mup3+p4
     &   mup3+p5mup3)*F(353))-p4mup3*p6mup4*F(359)-p5mup3*p6mup4*F(35
     &   9)+p2mup3*((p4mup4+p5mup4)*F(349)+p6mup4*F(358)-2*p2mup4*F(3
     &   60))+p1mup3*p6mup4*F(361)-p6mup3*p6mup4*F(363)-2*(p1mup3*p2m
     &   up4*F(362)+p2mup4*p4mup3*F(364)+p2mup4*p5mup3*F(364)+p2mup4*
     &   p6mup3*F(365)))
       Fa7(87) = p1mup4*F(369)+8*((p4mup4+p5mup4)*F(368)-p2mup4*F(370)
     &   )
       Fa7(88) = p1mup3*F(371)+8*p2mup3*F(372)+(p4mup3+p5mup3)*F(373)+
     &   p6mup3*F(374)
       Fa7(89) = -(mup3mup4*F(375))+2*(p5mup3*p5mup4*F(335)+p5mup4*p6m
     &   up3*F(337)+p1mup3*p5mup4*F(346)-p2mup3*p5mup4*F(347)+p4mup4*
     &   (p5mup3*F(335)+p6mup3*F(337)+p1mup3*F(346)-p2mup3*F(347))+p1
     &   mup4*p2mup3*F(358)-p1mup4*p5mup3*F(359)+p1mup3*p1mup4*F(361)
     &   -p1mup4*p6mup3*F(363)+p4mup3*((p4mup4+p5mup4)*F(335)-p1mup4*
     &   F(359)+2*p2mup4*F(377))-2*(p1mup3*p2mup4*F(376)-p2mup4*p5mup
     &   3*F(377)+p2mup4*p6mup3*F(378)-p2mup3*p2mup4*F(379)))
       Fa7(90) = mup3mup4*F(380)+2*(p4mup4+p5mup4)*(p1mup3*F(381)-p2mu
     &   p3*F(382))
       Fa7(91) = -(p2mup4*F(275))+(p4mup4+p5mup4)*F(387)
       Return
       End
