c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa8(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa8
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
       Complex*16 Fa8(134:152)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa8(134) = p1mup3*F(450)+p6mup3*F(471)
       Fa8(135) = p1mup3*F(300)-p5mup3*F(471)
       Fa8(136) = p1mup6*p6mup3*F(325)+p6mup3*p6mup6*F(454)+p5mup6*p6m
     &   up3*F(455)+p2mup6*p6mup3*F(458)+p5mup3*p6mup6*F(508)+p5mup3*
     &   p5mup6*F(509)+p1mup6*p5mup3*F(510)+p2mup6*p5mup3*F(511)
       Fa8(137) = p6mup6*F(403)-p1mup6*F(404)-p2mup6*F(405)+p5mup6*F(5
     &   12)
       Fa8(138) = -(p2mup3*F(450))+p6mup3*F(451)
       Fa8(139) = 8*p6mup3*F(300)+4*p5mup3*F(450)
       Fa8(140) = p2mup3*F(300)+p5mup3*F(451)
       Fa8(141) = p1mup3*F(450)+2*p6mup3*F(471)
       Fa8(142) = p1mup3*F(451)+p2mup3*F(471)
       Fa8(143) = -(p6mup3*(p6mup6*F(317)+p5mup6*F(318)+p1mup6*F(409)+
     &   p2mup6*F(410)))+p1mup3*(p6mup6*F(508)+p5mup6*F(509)+p1mup6*F
     &   (510)+p2mup6*F(511))
       Fa8(144) = p5mup3*(p6mup6*F(317)+p5mup6*F(318)+p1mup6*F(409)+p2
     &   mup6*F(410))+p1mup3*(p1mup6*F(325)+p6mup6*F(454)+p5mup6*F(45
     &   5)+p2mup6*F(458))
       Fa8(145) = p5mup3*(p6mup6*F(394)+p5mup6*F(399)-p2mup6*F(400)-p1
     &   mup6*F(401))+p6mup3*(p6mup6*F(452)+p5mup6*F(456)-p2mup6*F(45
     &   7)-p1mup6*F(459))
       Fa8(146) = p1mup3*(p6mup6*F(394)+p5mup6*F(399)-p2mup6*F(400)-p1
     &   mup6*F(401))+p6mup3*(p6mup6*F(406)+p2mup6*F(408)+p5mup6*F(41
     &   1)+p1mup6*F(412))
       Fa8(147) = p5mup3*(p6mup6*F(406)+p2mup6*F(408)+p5mup6*F(411)+p1
     &   mup6*F(412))+p1mup3*(-(p6mup6*F(452))-p5mup6*F(456)+p2mup6*F
     &   (457)+p1mup6*F(459))
       Fa8(148) = p2mup6*F(466)-p6mup6*F(513)+p1mup6*F(514)-p5mup6*F(5
     &   15)
       Fa8(149) = -(p6mup6*F(394))-p5mup6*F(399)+p2mup6*F(400)+p1mup6*
     &   F(401)
       Fa8(150) = p6mup6*F(508)+p5mup6*F(509)+p1mup6*F(510)+p2mup6*F(5
     &   11)
       Fa8(151) = -(p6mup6*F(452))-p5mup6*F(456)+p2mup6*F(457)+p1mup6*
     &   F(459)
       Fa8(152) = p1mup6*F(325)+p6mup6*F(454)+p5mup6*F(455)+p2mup6*F(4
     &   58)
       Return
       End
