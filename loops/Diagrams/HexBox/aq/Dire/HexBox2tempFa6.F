c************** Calling the Fa functions*************************
       subroutine HexBox2Fa6(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa6)
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
       Complex*16 Fa6(96:114)
       COMMON/HexBox2FaFunctions/Fa
       Fa6(96) = p1mup3*(p2mup6*F(350)+p1mup6*F(353)+p5mup6*F(356)+p6m
     &   up6*F(357))+p5mup6*p6mup3*F(360)+p6mup3*p6mup6*F(361)+p2mup6
     &   *(p6mup3*F(352)-8*(p5mup3*F(351)+p2mup3*F(363)))+8*(p1mup6*p
     &   5mup3*F(354)+p1mup6*p6mup3*F(355)-p5mup3*p5mup6*F(358)-p5mup
     &   3*p6mup6*F(359)-p2mup3*p5mup6*F(362)+p2mup3*p6mup6*F(364)-p1
     &   mup6*p2mup3*F(365))
       Fa6(97) = p2mup6*F(366)+p1mup6*F(367)+p5mup6*F(368)+p6mup6*F(36
     &   9)
       Fa6(98) = p2mup6*F(370)-p1mup6*F(371)+p5mup6*F(372)+p6mup6*F(37
     &   3)
       Fa6(99) = p6mup3*F(374)+p1mup3*F(375)-p2mup3*F(376)
       Fa6(100) = -4*mup3mup6*F(377)+8*(-(p2mup6*p6mup3*F(378))+p5mup6
     &   *p6mup3*F(379)+p6mup3*p6mup6*F(380)-p2mup3*p5mup6*F(381)+p1m
     &   up6*p6mup3*F(382)+p2mup3*p2mup6*F(383)+p2mup3*p6mup6*F(384)+
     &   p1mup6*p2mup3*F(385))
       Fa6(101) = 4*mup3mup6*F(386)+8*(p1mup3*p2mup6*F(387)+p1mup3*p6m
     &   up6*F(390)+p1mup3*p1mup6*F(391)+p2mup3*p2mup6*F(392)+p5mup6*
     &   (p6mup3*F(302)+p1mup3*F(389)+p2mup3*F(393))+p2mup3*p6mup6*F(
     &   394)+p1mup6*p2mup3*F(395)+p6mup3*(p6mup6*F(303)+p2mup6*F(388
     &   )+p1mup6*F(396)))
       Fa6(102) = -(p2mup6*F(387))-p1mup6*F(391)+p5mup6*F(397)+p6mup6*
     &   F(398)
       Fa6(103) = p6mup3*F(400)-p1mup3*F(401)+p5mup3*F(402)+p2mup3*F(4
     &   03)
       Fa6(104) = p5mup3*(p2mup6*F(336)+2*(p5mup6*F(334)+p6mup6*F(335)
     &   +p1mup6*F(337)))+p2mup3*(-(p2mup6*F(404))-p1mup6*F(405)+p5mu
     &   p6*F(407)+p6mup6*F(408))+p6mup3*(p2mup6*F(406)-p5mup6*F(409)
     &   -p6mup6*F(410)-2*p1mup6*F(411))
       Fa6(105) = p1mup3*(p2mup6*F(406)-p5mup6*F(409)-p6mup6*F(410)-2*
     &   p1mup6*F(411))-p1mup6*p2mup3*F(412)+p2mup3*p5mup6*F(413)+p2m
     &   up3*p6mup6*F(414)-p2mup6*p5mup3*F(415)+p1mup6*p5mup3*F(416)+
     &   2*p2mup3*p2mup6*F(417)
       Fa6(106) = p6mup3*F(416)+p2mup3*F(418)
       Fa6(107) = -(p1mup3*F(416))+p2mup3*F(419)
       Fa6(108) = p5mup3*F(422)+4*(p6mup3*F(420)+p1mup3*F(421)+p2mup3*
     &   F(423))
       Fa6(109) = p1mup3*F(377)+p6mup3*F(424)
       Fa6(110) = p2mup3*F(418)+p6mup3*F(429)
       Fa6(111) = p2mup3*F(419)-p1mup3*F(429)
       Fa6(112) = p5mup3*F(375)+p2mup3*F(376)-p6mup3*F(430)
       Fa6(113) = mup3mup6*F(431)-2*(p5mup3*(p5mup6*F(389)+p1mup6*F(39
     &   1)-p6mup6*F(398))-p1mup6*p6mup3*F(411)-p5mup6*p6mup3*F(433)-
     &   p6mup3*p6mup6*F(434)+p2mup3*p5mup6*F(435)+p2mup6*(p5mup3*F(3
     &   87)+p6mup3*F(432)-p2mup3*F(436))-p2mup3*p6mup6*F(437)-p1mup6
     &   *p2mup3*F(438))
       Fa6(114) = -(p2mup6*F(406))+p5mup6*F(409)+p6mup6*F(410)+2*p1mup
     &   6*F(411)
       Return
       End
