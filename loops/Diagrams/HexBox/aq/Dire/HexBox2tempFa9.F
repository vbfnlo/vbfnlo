c************** Calling the Fa functions*************************
       subroutine HexBox2Fa9(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa9)
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
       Complex*16 Fa9(153:171)
       COMMON/HexBox2FaFunctions/Fa
       Fa9(153) = p5mup6*p6mup3*F(302)+p6mup3*p6mup6*F(303)+p2mup6*p6m
     &   up3*F(388)+p1mup6*p6mup3*F(396)+p1mup3*p2mup6*F(451)-p1mup3*
     &   p5mup6*F(452)-p1mup3*p6mup6*F(453)-p1mup3*p1mup6*F(454)
       Fa9(154) = p5mup3*(p5mup6*F(302)+p6mup6*F(303)+p2mup6*F(388)+p1
     &   mup6*F(396))+p1mup3*(p1mup6*F(411)-p2mup6*F(432)+p5mup6*F(43
     &   3)+p6mup6*F(434))
       Fa9(155) = p1mup3*F(418)+p6mup3*F(419)
       Fa9(156) = -(p5mup3*F(419))+p1mup3*F(440)
       Fa9(157) = p5mup6*F(302)+p6mup6*F(303)+p2mup6*F(388)+p1mup6*F(3
     &   96)
       Fa9(158) = p2mup6*F(392)+p5mup6*F(393)+p6mup6*F(394)+p1mup6*F(3
     &   95)
       Fa9(159) = 4*p6mup6*F(533)+(p1mup6*F(534))/s16
       Fa9(160) = p1mup6*(p1mup3+p6mup3)*F(535)+p1mup6*p2mup3*F(536)
       Fa9(161) = (p6mup3*F(539))/s16+8*(p4mup3*F(540)+p1mup3*F(541)+p
     &   5mup3*F(542)-p2mup3*F(543))
       Fa9(162) = p1mup6*((p1mup3+p6mup3)*F(543)+p2mup3*F(545)+p5mup3*
     &   F(546)+p4mup3*F(547))
       Fa9(163) = (p1mup3+p6mup3)*F(543)+p2mup3*F(545)+p5mup3*F(546)+p
     &   4mup3*F(547)
       Fa9(164) = (p1mup3+p6mup3)*F(535)+p2mup3*F(536)
       Fa9(165) = p1mup6*(p2mup3*F(556)+8*((p1mup3+p6mup3)*F(555)-p5mu
     &   p3*F(557)-p4mup3*F(558)))
       Fa9(166) = -4*p5mup6*F(533)+(p1mup6*F(559))/s16
       Fa9(167) = p6mup3*F(561)+s16*(p4mup3*F(562)+p1mup3*F(563)+p5mup
     &   3*F(564)+p2mup3*F(565))
       Fa9(168) = 4*mup3mup6*s16*F(533)-8*p1mup6*((p1mup3+p6mup3)*F(56
     &   6)+p5mup3*F(567)+p4mup3*F(568)+p2mup3*F(569))
       Fa9(169) = -4*mup3mup6*s16*F(533)+8*p1mup6*((p1mup3+p6mup3)*F(5
     &   66)+p5mup3*F(567)+p4mup3*F(568)+p2mup3*F(569))
       Fa9(170) = p1mup6*((p1mup3+p6mup3)*F(570)+p5mup3*F(571)+p4mup3*
     &   F(572)+p2mup3*F(573))
       Fa9(171) = p5mup3*F(571)+p4mup3*F(572)+(p6mup3*F(574))/s16+p1mu
     &   p3*F(575)+p2mup3*F(576)
       Return
       End
