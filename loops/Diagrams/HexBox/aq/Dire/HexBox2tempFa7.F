c************** Calling the Fa functions*************************
       subroutine HexBox2Fa7(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa7)
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
       Complex*16 Fa7(115:133)
       COMMON/HexBox2FaFunctions/Fa
       Fa7(115) = p5mup3*F(416)-p2mup3*F(440)
       Fa7(116) = p5mup3*F(377)-p6mup3*F(441)
       Fa7(117) = p5mup3*F(429)-p2mup3*F(440)
       Fa7(118) = p2mup6*F(443)-p5mup6*F(444)+p6mup6*F(445)+p1mup6*F(4
     &   46)
       Fa7(119) = p2mup6*F(447)-p5mup6*F(448)+p6mup6*F(449)+p1mup6*F(4
     &   50)
       Fa7(120) = p5mup6*p6mup3*F(433)+p6mup3*p6mup6*F(434)-p2mup6*(p6
     &   mup3*F(432)+p5mup3*F(451))+p5mup3*p5mup6*F(452)+p5mup3*p6mup
     &   6*F(453)+p1mup6*(p6mup3*F(411)+p5mup3*F(454))
       Fa7(121) = -(p2mup6*F(451))+p5mup6*F(452)+p6mup6*F(453)+p1mup6*
     &   F(454)
       Fa7(122) = p1mup6*F(411)-p2mup6*F(432)+p5mup6*F(433)+p6mup6*F(4
     &   34)
       Fa7(123) = p6mup3*F(455)+p1mup3*F(456)+p5mup3*F(457)
       Fa7(124) = 2*(p5mup3*(p5mup6*F(323)+p6mup6*F(324)+p1mup6*F(325)
     &   )-p2mup6*(p5mup3*F(322)-p6mup3*F(460)))+p6mup3*(p5mup6*F(458
     &   )+p6mup6*F(459)+2*p1mup6*F(461))
       Fa7(125) = -(p6mup3*p6mup6*F(414))-p1mup3*p2mup6*F(462)+p1mup3*
     &   p6mup6*F(465)+p5mup6*(-(p6mup3*F(413))+p1mup3*F(464)+2*p5mup
     &   3*F(467))-2*(p2mup6*p6mup3*F(417)-p2mup6*p5mup3*F(466)-p5mup
     &   3*p6mup6*F(468))+p1mup6*(p6mup3*F(412)-p1mup3*F(463)+p5mup3*
     &   F(469))
       Fa7(126) = p6mup3*F(473)+p1mup3*F(474)+p5mup3*F(475)
       Fa7(127) = p1mup3*p2mup6*F(476)+p1mup3*p1mup6*F(477)-p1mup3*p6m
     &   up6*F(479)+p5mup6*(p6mup3*F(413)-p1mup3*F(478)+p5mup3*F(480)
     &   )+p5mup3*p6mup6*F(481)-p2mup6*p5mup3*F(482)-p1mup6*p5mup3*F(
     &   483)+p6mup3*(p6mup6*F(414)+2*p2mup6*F(417)-p1mup6*F(484))
       Fa7(128) = 4*p5mup3*F(418)+8*p6mup3*F(440)
       Fa7(129) = p5mup6*F(488)+p2mup6*F(489)+p6mup6*F(490)+p1mup6*F(4
     &   91)
       Fa7(130) = p5mup3*(-(p5mup6*F(381))+p2mup6*F(383)+p6mup6*F(384)
     &   +p1mup6*F(385))+p6mup3*(-(p5mup6*F(435))+p2mup6*F(436)+p6mup
     &   6*F(437)+p1mup6*F(438))
       Fa7(131) = p5mup6*F(381)-p2mup6*F(383)-p6mup6*F(384)-p1mup6*F(3
     &   85)
       Fa7(132) = p5mup6*F(435)-p2mup6*F(436)-p6mup6*F(437)-p1mup6*F(4
     &   38)
       Fa7(133) = p2mup6*F(404)+p1mup6*F(405)-p5mup6*F(407)-p6mup6*F(4
     &   08)
       Return
       End
