c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa10(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   10)
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
       Complex*16 Fa10(172:190)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa10(172) = p2mup3*F(617)-(p1mup3+p6mup3)*F(618)+p4mup3*F(619)+
     &   p5mup3*F(620)
       Fa10(173) = p2mup3*F(621)+(p1mup3+p6mup3)*F(622)-p4mup3*F(623)-
     &   p5mup3*F(624)
       Fa10(174) = -8*((p1mup3+p6mup3)*F(625)-p4mup3*F(626)-p5mup3*F(6
     &   27))+p2mup3*F(628)
       Fa10(175) = (p1mup6+p6mup6)*F(632)-p4mup6*F(633)-p5mup6*F(634)+
     &   p2mup6*F(635)
       Fa10(176) = (p1mup6+p6mup6)*F(636)+p4mup6*F(637)+p2mup6*F(638)
       Fa10(177) = (p1mup6+p6mup6)*F(639)+p2mup6*F(640)+p4mup6*F(641)-
     &   p5mup6*F(642)
       Fa10(178) = (p1mup6+p6mup6)*((p1mup3+p6mup3)*F(533)+p4mup3*F(53
     &   4)+p5mup3*F(535)+p2mup3*F(644))
       Fa10(179) = p1mup6*p6mup3*F(533)+p6mup3*p6mup6*F(533)+p1mup6*p4
     &   mup3*F(534)+p4mup3*p6mup6*F(534)+p1mup6*p5mup3*F(535)+p5mup3
     &   *p6mup6*F(535)+p4mup6*p5mup3*F(567)+p4mup3*p4mup6*F(570)+p4m
     &   up6*p6mup3*F(577)-p2mup3*p4mup6*F(579)+p2mup6*p6mup3*F(592)+
     &   p1mup3*((p1mup6+p6mup6)*F(533)+p4mup6*F(577)+p2mup6*F(592))+
     &   p2mup6*p4mup3*F(593)+p2mup6*p5mup3*F(594)-p2mup3*p2mup6*F(59
     &   5)+p2mup3*(p1mup6+p6mup6)*F(644)+mup3mup6*F(645)
       Fa10(180) = p1mup3*F(648)+p4mup3*F(649)+p5mup3*F(650)+4*(p6mup3
     &   *F(647)+p2mup3*F(651))
       Fa10(181) = p1mup6*F(652)+p6mup6*F(653)-4*(p5mup6*F(634)-p2mup6
     &   *F(654)-p4mup6*F(655))
       Fa10(182) = p1mup6*F(656)+p6mup6*F(657)+4*(p2mup6*F(658)+p4mup6
     &   *F(659))
       Fa10(183) = (p1mup6+p6mup6)*F(653)-4*(p5mup6*F(634)-p2mup6*F(66
     &   0)-p4mup6*F(661))
       Fa10(184) = 4*(p6mup3*F(665)+p2mup3*F(666)+p1mup3*F(667)+p4mup3
     &   *F(668))+p5mup3*F(669)
       Fa10(185) = (p1mup3+p6mup3)*F(670)+p2mup3*F(671)+p4mup3*F(672)+
     &   p5mup3*F(673)
       Fa10(186) = p2mup6*F(555)+p4mup6*F(680)
       Fa10(187) = (p1mup3+p6mup3)*F(533)+p4mup3*F(534)+p5mup3*F(535)+
     &   p2mup3*F(539)
       Fa10(188) = mup3mup6*F(679)+p1mup3*(p2mup6*F(536)+p4mup6*F(682)
     &   )+p2mup6*(p6mup3*F(536)+p2mup3*F(540)+p5mup3*F(683)+p4mup3*F
     &   (685))+p4mup6*(p5mup3*F(681)+p6mup3*F(682)+p4mup3*F(684)+p2m
     &   up3*F(686))
       Fa10(189) = p4mup6*F(554)+p5mup6*F(555)
       Fa10(190) = (p1mup6+p6mup6)*F(555)-p4mup6*F(689)
       Return
       End
