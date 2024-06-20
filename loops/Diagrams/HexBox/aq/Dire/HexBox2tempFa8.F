c************** Calling the Fa functions*************************
       subroutine HexBox2Fa8(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa8)
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
       Complex*16 Fa8(134:152)
       COMMON/HexBox2FaFunctions/Fa
       Fa8(134) = p5mup6*F(295)+p6mup6*F(296)-p2mup6*F(301)-p1mup6*F(4
     &   92)
       Fa8(135) = 4*(p5mup6*F(413)+p6mup6*F(414))+8*p2mup6*F(417)+p1mu
     &   p6*F(493)
       Fa8(136) = p5mup3*F(418)+p6mup3*F(440)
       Fa8(137) = p2mup6*F(496)+p1mup6*F(497)+p5mup6*F(498)+p6mup6*F(4
     &   99)
       Fa8(138) = p2mup6*F(500)-p1mup6*F(501)+p5mup6*F(502)+p6mup6*F(5
     &   03)
       Fa8(139) = p2mup6*F(504)-p1mup6*F(505)+p5mup6*F(506)+p6mup6*F(5
     &   07)
       Fa8(140) = 16*p5mup3*F(341)-8*p6mup3*F(513)
       Fa8(141) = -(p5mup3*F(415))+p1mup3*F(513)
       Fa8(142) = -(p2mup6*F(517))+p5mup6*F(518)+p6mup6*F(519)-p1mup6*
     &   F(520)
       Fa8(143) = p2mup3*F(376)+p5mup3*F(399)+p1mup3*F(439)+p6mup3*F(5
     &   21)
       Fa8(144) = -(p2mup6*F(462))-p1mup6*F(463)+p5mup6*F(464)+p6mup6*
     &   F(465)
       Fa8(145) = p5mup3*F(386)+p1mup3*F(431)+6*p6mup3*F(1636)
       Fa8(146) = p1mup3*F(418)+2*p6mup3*F(419)
       Fa8(147) = -(p1mup6*F(293))+p5mup6*F(295)+p6mup6*F(296)-p2mup6*
     &   F(301)
       Fa8(148) = -(p1mup6*F(412))+p5mup6*F(413)+p6mup6*F(414)+2*p2mup
     &   6*F(417)
       Fa8(149) = -(p2mup6*F(376))-p1mup6*F(525)+p5mup6*F(526)+p6mup6*
     &   F(1642)
       Fa8(150) = (p1mup3+p5mup3+p6mup3)*(p2mup6*F(376)+p1mup6*F(525)-
     &   p5mup6*F(526)-p6mup6*F(1642))
       Fa8(151) = p1mup3*(p5mup6*F(381)-p2mup6*F(383)-p6mup6*F(384)-p1
     &   mup6*F(385))+p6mup3*(p2mup6*F(392)+p5mup6*F(393)+p6mup6*F(39
     &   4)+p1mup6*F(395))
       Fa8(152) = p2mup6*p5mup3*F(392)+p5mup3*p5mup6*F(393)+p5mup3*p6m
     &   up6*F(394)+p1mup6*p5mup3*F(395)-p1mup3*p5mup6*F(435)+p1mup3*
     &   p2mup6*F(436)+p1mup3*p6mup6*F(437)+p1mup3*p1mup6*F(438)
       Return
       End
