c************** Calling the Fa functions*************************
       subroutine HexBox2Fa10(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa10)
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
       Complex*16 Fa10(172:190)
       COMMON/HexBox2FaFunctions/Fa
       Fa10(172) = p2mup3*F(577)+4*((p1mup3+p6mup3)*F(555)-p4mup3*F(55
     &   8)+p5mup3*F(578))
       Fa10(173) = (p6mup3*F(561))/s16+p4mup3*F(562)+p1mup3*F(563)+p5m
     &   up3*F(564)+p2mup3*F(565)
       Fa10(174) = (p1mup3+p6mup3)*F(566)+p5mup3*F(567)+p4mup3*F(568)+
     &   p2mup3*F(569)
       Fa10(175) = p6mup3*F(584)+s16*(p4mup3*F(570)+p1mup3*F(585)+p5mu
     &   p3*F(586)+p2mup3*F(587))
       Fa10(176) = (p1mup3+p6mup3)*F(570)+p5mup3*F(571)+p4mup3*F(572)+
     &   p2mup3*F(573)
       Fa10(177) = p4mup3*F(570)+(p6mup3*F(584))/s16+p1mup3*F(585)+p5m
     &   up3*F(586)+p2mup3*F(587)
       Fa10(178) = p1mup6*(p5mup3*F(540)+(p1mup3+p6mup3)*F(590)+p4mup3
     &   *F(591)+p2mup3*F(592))
       Fa10(179) = p5mup3*F(540)+(p1mup3+p6mup3)*F(590)+p4mup3*F(591)+
     &   p2mup3*F(592)
       Fa10(180) = s16*(p5mup3*F(565)+p4mup3*F(576)+p1mup3*F(587)+p2mu
     &   p3*F(595))-p6mup3*F(596)
       Fa10(181) = p1mup6*(p2mup3*F(597)-(p1mup3+p6mup3)*F(598)-p5mup3
     &   *F(599)-p4mup3*F(600))
       Fa10(182) = p1mup6*(p5mup3*F(565)+p4mup3*F(576)+(p1mup3+p6mup3)
     &   *F(587)+p2mup3*F(607))
       Fa10(183) = p5mup3*F(565)+p4mup3*F(576)+(p1mup3+p6mup3)*F(587)+
     &   p2mup3*F(607)
       Fa10(184) = p2mup3*F(597)-(p1mup3+p6mup3)*F(598)-p5mup3*F(599)-
     &   p4mup3*F(600)
       Fa10(185) = p1mup6*(-(p2mup3*F(617))+p5mup3*F(618)+p4mup3*F(619
     &   )+(p1mup3+p6mup3)*F(1640))
       Fa10(186) = p1mup6*((p1mup3+p6mup3)*F(620)+p5mup3*F(621)+p4mup3
     &   *F(622)-p2mup3*F(623))
       Fa10(187) = p1mup6*(-(p2mup3*F(547))+p5mup3*F(624)+(p1mup3+p6mu
     &   p3)*F(625)+p4mup3*F(626))
       Fa10(188) = -(p2mup3*F(617))+p5mup3*F(618)+p4mup3*F(619)+(p1mup
     &   3+p6mup3)*F(1640)
       Fa10(189) = (p1mup3+p6mup3)*F(620)+p5mup3*F(621)+p4mup3*F(622)-
     &   p2mup3*F(623)
       Fa10(190) = -(p2mup3*F(547))+p5mup3*F(624)+(p1mup3+p6mup3)*F(62
     &   5)+p4mup3*F(626)
       Return
       End
