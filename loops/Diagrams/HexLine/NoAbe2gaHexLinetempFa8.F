c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa8(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa8)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa8(169:192)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa8(169) = p2mup3*p5mup2*(p4mup6*p5mup5*F(583)+p5mup5*p5mup6*F(
     -   585)+p4mup5*p4mup6*F(594)-p4mup5*p5mup6*F(595))
       Fa8(170) = p5mup2*(p1mup5*p1mup6*p2mup3*F(586)+p1mup3*(p1mup5*(
     -   -(p4mup6*F(582))+(p2mup6+p3mup6)*F(586)-p5mup6*F(587))+p1mup
     -   6*(p4mup5*F(584)+(p2mup5+p3mup5)*F(586)+p5mup5*F(588))))
       Fa8(171) = (p1mup6*p2mup3*(p2mup5+p3mup5)+p1mup3*p2mup6*(p2mup5
     -   +p3mup5)+p1mup3*p2mup5*p3mup6+p1mup5*p2mup3*(p2mup6+p3mup6))
     -   *p5mup2
       Fa8(172) = p5mup2*(p1mup5*p2mup3*(p4mup6*F(582)+p5mup6*F(587))-
     -   p1mup6*p2mup3*(p4mup5*F(584)+p5mup5*F(588))+p1mup3*(p2mup5*p
     -   4mup6*F(582)+p3mup5*p4mup6*F(582)-(p2mup6+p3mup6)*p4mup5*F(5
     -   84)-p3mup5*p3mup6*F(586)+p2mup5*p5mup6*F(587)-p2mup6*p5mup5*
     -   F(588)))
       Fa8(173) = p5mup2*(-(p2mup3*p2mup5*p4mup6*F(582))+p2mup3*(p2mup
     -   6*p3mup5+p2mup5*(p2mup6+p3mup6))*F(586)+p1mup3*(p4mup6*p5mup
     -   5*F(583)+p5mup5*p5mup6*F(585)-p3mup5*p5mup6*F(587)+p3mup6*p5
     -   mup5*F(588)+p4mup5*p4mup6*F(594)-p4mup5*p5mup6*F(595)))
       Fa8(174) = p2mup3*p5mup2*(-(p2mup5*p5mup6*F(587))-p3mup5*(p4mup
     -   6*F(582)-p3mup6*F(586)+p5mup6*F(587))+(p2mup6+p3mup6)*(p4mup
     -   5*F(584)+p5mup5*F(588)))+p2mup3*(p1mup6*p3mup5+p1mup5*p3mup6
     -   )*p6mup2*F(591)+p6mup2*(p1mup3*((p1mup6+p2mup6)*p3mup5+(p1mu
     -   p5+p2mup5)*p3mup6)*F(591)+p5mup3*p5mup5*p5mup6*F(631)-p3mup5
     -   *p3mup6*p5mup3*F(632))
       Fa8(175) = p6mup2*(p1mup3*(p5mup5*(p5mup6*F(567)+(p1mup6+p2mup6
     -   +p3mup6)*F(568)+p4mup6*F(569))+p1mup6*p4mup5*F(590)+p2mup6*p
     -   4mup5*F(590)+p3mup6*p4mup5*F(590)-p1mup6*(p1mup5+p2mup5)*F(5
     -   91)-p1mup5*p2mup6*F(591)-p2mup5*p2mup6*F(591)-p4mup5*p4mup6*
     -   F(596)-p1mup5*p4mup6*F(634)-p2mup5*p4mup6*F(634)-p3mup5*p4mu
     -   p6*F(634)-p1mup5*p5mup6*F(635)-p2mup5*p5mup6*F(635)-p3mup5*p
     -   5mup6*F(635)-p4mup5*p5mup6*F(640))+p2mup3*(p5mup5*(p5mup6*F(
     -   567)+(p1mup6+p2mup6+p3mup6)*F(568)+p4mup6*F(569))+p1mup6*p4m
     -   up5*F(590)+p2mup6*p4mup5*F(590)+p3mup6*p4mup5*F(590)-p1mup5*
     -   p1mup6*F(591)-p1mup5*p2mup6*F(591)-p2mup6*p3mup5*F(591)-p2mu
     -   p5*p3mup6*F(591)-p4mup5*p4mup6*F(596)-p1mup5*p4mup6*F(634)-p
     -   2mup5*p4mup6*F(634)-p3mup5*p4mup6*F(634)-p1mup5*p5mup6*F(635
     -   )-p2mup5*p5mup6*F(635)-p3mup5*p5mup6*F(635)-p4mup5*p5mup6*F(
     -   640))+p5mup3*(p2mup5*p2mup6*F(632)+p2mup6*p3mup5*F(632)+p2mu
     -   p5*p3mup6*F(632)+p2mup6*p5mup5*F(633)+p3mup6*p5mup5*F(633)-p
     -   4mup6*p5mup5*F(636)+p2mup5*p4mup6*F(637)+p3mup5*p4mup6*F(637
     -   )-p2mup6*p4mup5*F(638)-p3mup6*p4mup5*F(638)+p1mup6*((p2mup5+
     -   p3mup5)*F(632)+p5mup5*F(633)-p4mup5*F(638))+p4mup5*p5mup6*F(
     -   639)+p2mup5*p5mup6*F(641)+p3mup5*p5mup6*F(641)+p1mup5*((p1mu
     -   p6+p2mup6+p3mup6)*F(632)+p4mup6*F(637)+p5mup6*F(641))-p4mup5
     -   *p4mup6*F(642)))
       Fa8(176) = p1mup3*p3mup5*p3mup6*p6mup2*F(591)+p2mup3*(p2mup5*p2
     -   mup6+p3mup5*p3mup6)*p6mup2*F(591)+p2mup6*p3mup5*p5mup2*p6mup
     -   3*F(632)+p1mup5*p3mup6*p5mup2*p6mup3*F(632)+p2mup5*p3mup6*p5
     -   mup2*p6mup3*F(632)+p3mup5*p3mup6*p5mup2*p6mup3*F(632)+p5mup2
     -   *p6mup3*(-(p5mup5*p5mup6*F(631))+(p1mup5+p2mup5)*p2mup6*F(63
     -   2))+p2mup6*p5mup2*p5mup5*p6mup3*F(633)+p3mup6*p5mup2*p5mup5*
     -   p6mup3*F(633)-p4mup6*p5mup2*p5mup5*p6mup3*F(636)+p1mup5*p4mu
     -   p6*p5mup2*p6mup3*F(637)+p2mup5*p4mup6*p5mup2*p6mup3*F(637)+p
     -   3mup5*p4mup6*p5mup2*p6mup3*F(637)-p2mup6*p4mup5*p5mup2*p6mup
     -   3*F(638)-p3mup6*p4mup5*p5mup2*p6mup3*F(638)+p1mup6*(p2mup3*p
     -   2mup5*p6mup2*F(591)+p5mup2*p6mup3*((p1mup5+p2mup5+p3mup5)*F(
     -   632)+p5mup5*F(633)-p4mup5*F(638)))+p4mup5*p5mup2*p5mup6*p6mu
     -   p3*F(639)+p1mup5*p5mup2*p5mup6*p6mup3*F(641)+p2mup5*p5mup2*p
     -   5mup6*p6mup3*F(641)+p3mup5*p5mup2*p5mup6*p6mup3*F(641)-p4mup
     -   5*p4mup6*p5mup2*p6mup3*F(642)
       Fa8(177) = p5mup2*p6mup5*(p2mup3*(p1mup6+p2mup6+p3mup6)*F(571)+
     -   p1mup3*(p5mup6*F(570)+(p1mup6+p2mup6+p3mup6)*F(571)+p4mup6*F
     -   (572)))
       Fa8(178) = p2mup3*p5mup2*p6mup5*(p5mup6*F(570)+p4mup6*F(572))
       Fa8(179) = p6mup5*(p2mup3*p2mup6*p6mup2*F(573)+p2mup3*p3mup6*p6
     -   mup2*F(573)+p2mup3*p4mup6*p6mup2*F(574)+p2mup3*p5mup6*p6mup2
     -   *F(575)+p1mup3*p6mup2*((p1mup6+p2mup6+p3mup6)*F(573)+p4mup6*
     -   F(574)+p5mup6*F(575))-p4mup6*p5mup3*p6mup2*F(643)+p4mup6*p5m
     -   up2*p6mup3*F(643)-p2mup6*p5mup3*p6mup2*F(644)-p3mup6*p5mup3*
     -   p6mup2*F(644)+p2mup6*p5mup2*p6mup3*F(644)+p3mup6*p5mup2*p6mu
     -   p3*F(644)+p1mup6*(p2mup3*p6mup2*F(573)+(-(p5mup3*p6mup2)+p5m
     -   up2*p6mup3)*F(644))+p5mup6*(p5mup3*p6mup2-p5mup2*p6mup3)*F(6
     -   45))
       Fa8(180) = p5mup2*p6mup6*(p1mup5*p2mup3*F(576)+p1mup3*((p1mup5+
     -   p2mup5+p3mup5)*F(576)+p4mup5*F(577)+p5mup5*F(578)))
       Fa8(181) = p6mup6*(p2mup3*(-(p2mup5*p5mup2*F(576))-p3mup5*p5mup
     -   2*F(576)-p4mup5*p5mup2*F(577)-p5mup2*p5mup5*F(578)+p1mup5*p6
     -   mup2*F(579)+p2mup5*p6mup2*F(579)+p3mup5*p6mup2*F(579)+p4mup5
     -   *p6mup2*F(580)+p5mup5*p6mup2*F(647))+p6mup2*((p1mup5+p2mup5+
     -   p3mup5)*p5mup3*F(646)+p1mup3*((p1mup5+p2mup5+p3mup5)*F(579)+
     -   p4mup5*F(580)+p5mup5*F(647))-p5mup3*p5mup5*F(648)))
       Fa8(182) = p6mup6*(p4mup5*p5mup3*p6mup2*F(649)+p5mup2*p6mup3*((
     -   p1mup5+p2mup5+p3mup5)*F(646)-p5mup5*F(648)-p4mup5*F(649))+2*
     -   ((p5mup3*p6mup2-p5mup2*p6mup3)*p6mup5*F(650)+p1mup3*p6mup5*(
     -   p5mup2*F(581)-2*p6mup2*F(9567))+p2mup3*p6mup5*(p5mup2*F(581)
     -   -2*p6mup2*F(9567))))
       Fa8(183) = p5mup3*p6mup2*F(651)-p5mup2*p6mup3*F(652)
       Fa8(184) = -((p1mup2+p2mup2)*(p1mup5*p5mup3*F(653)+p2mup5*p5mup
     -   3*F(653)+p3mup5*p5mup3*F(653)+p1mup5*p6mup3*F(654)+p2mup5*p6
     -   mup3*F(654)+p3mup5*p6mup3*F(654)+p5mup3*p5mup5*F(655)+p5mup5
     -   *p6mup3*F(656)+p5mup3*p6mup5*F(657)+p6mup3*p6mup5*F(658)-p4m
     -   up5*p5mup3*F(659)-p4mup5*p6mup3*F(660)))+8*(p1mup2+p2mup2)*(
     -   p1mup3+p2mup3)*(p1mup5+p2mup5+p3mup5+p4mup5+p5mup5+p6mup5)*F
     -   (9573)
       Fa8(185) = mup5mup6*(p1mup2+p2mup2)*p6mup3
       Fa8(186) = p1mup2*p1mup5*p1mup6*p2mup3+p1mup3*(p1mup2*p1mup6*p2
     -   mup5+p1mup5*(p1mup6*p2mup2+p1mup2*(p2mup6+p3mup6)))
       Fa8(187) = p2mup2*(p1mup5*p1mup6*p2mup3+p1mup3*p1mup6*p2mup5+p1
     -   mup3*p1mup5*p2mup6+p2mup3*p2mup5*p2mup6)*F(586)+p1mup2*(p2mu
     -   p3*(p1mup6*p2mup5+p1mup5*(p2mup6+p3mup6))*F(586)+p1mup3*((p2
     -   mup6*p3mup5+p2mup5*(p2mup6+p3mup6))*F(586)-p1mup5*(p4mup6*F(
     -   582)+p5mup6*F(587))+p1mup6*(p4mup5*F(584)+p3mup5*F(586)+p5mu
     -   p5*F(588))))
       Fa8(188) = p2mup3*(p2mup2*(p1mup6*p2mup5+p1mup5*(p2mup6+p3mup6)
     -   )*F(586)+p1mup2*((p2mup6*p3mup5+p2mup5*(p2mup6+p3mup6))*F(58
     -   6)-p1mup5*(p4mup6*F(582)+p5mup6*F(587))+p1mup6*(p4mup5*F(584
     -   )+p3mup5*F(586)+p5mup5*F(588))))+p1mup3*(p1mup6*p2mup2*p4mup
     -   5*F(584)+p2mup2*((p1mup6+p2mup6)*p3mup5+p2mup5*(p2mup6+p3mup
     -   6))*F(586)-p1mup5*p2mup2*(p4mup6*F(582)-p3mup6*F(586)+p5mup6
     -   *F(587))+p1mup6*p2mup2*p5mup5*F(588)+p1mup2*(-(p3mup5*p4mup6
     -   *F(582))+p4mup6*p5mup5*F(583)+(p2mup6+p3mup6)*p4mup5*F(584)+
     -   p5mup5*p5mup6*F(585)+p3mup5*p3mup6*F(586)-p3mup5*p5mup6*F(58
     -   7)-p2mup5*(p4mup6*F(582)+p5mup6*F(587))+p2mup6*p5mup5*F(588)
     -   +p3mup6*p5mup5*F(588)+p4mup5*p4mup6*F(594)-p4mup5*p5mup6*F(5
     -   95)))
       Fa8(189) = -(p1mup5*p2mup2*p2mup3*(p4mup6*F(582)+p5mup6*F(587))
     -   )+p1mup3*p2mup2*(-(p2mup5*(p4mup6*F(582)+p5mup6*F(587)))-p3m
     -   up5*(p4mup6*F(582)-p3mup6*F(586)+p5mup6*F(587))+(p2mup6+p3mu
     -   p6)*(p4mup5*F(584)+p5mup5*F(588))+p4mup5*p4mup6*F(594))+p2mu
     -   p3*(p1mup6*p2mup2*(p4mup5*F(584)+p3mup5*F(586)+p5mup5*F(588)
     -   )+p1mup2*(-(p3mup5*p4mup6*F(582))+(p2mup6+p3mup6)*p4mup5*F(5
     -   84)+p3mup5*p3mup6*F(586)-p3mup5*p5mup6*F(587)-p2mup5*(p4mup6
     -   *F(582)+p5mup6*F(587))+p2mup6*p5mup5*F(588)+p3mup6*p5mup5*F(
     -   588)+p4mup5*p4mup6*F(594)-p4mup5*p5mup6*F(595)))
       Fa8(190) = (p1mup3*p2mup2+p1mup2*p2mup3)*p5mup5*(p4mup6*F(583)+
     -   p5mup6*F(585))-p1mup3*p2mup2*p4mup5*p5mup6*F(595)
       Fa8(191) = p1mup2*p1mup3*p1mup5*p1mup6*F(586)+p2mup2*p2mup3*(p4
     -   mup6*p5mup5*F(583)+(p2mup6+p3mup6)*p4mup5*F(584)+p3mup5*(-(p
     -   4mup6*F(582))+(p2mup6+p3mup6)*F(586)-p5mup6*F(587))-p2mup5*(
     -   p4mup6*F(582)-p3mup6*F(586)+p5mup6*F(587))+p2mup6*p5mup5*F(5
     -   88)+p3mup6*p5mup5*F(588)+p4mup5*p4mup6*F(594)-p4mup5*p5mup6*
     -   F(595))
       Fa8(192) = -((p1mup2+p2mup2)*p6mup3*(-((p1mup5+p2mup5+p3mup5)*(
     -   p1mup6+p2mup6+p3mup6)*F(632))+p5mup5*(p5mup6*F(631)-(p1mup6+
     -   p2mup6+p3mup6)*F(633)+p4mup6*F(636))+(p1mup6+p2mup6+p3mup6)*
     -   p4mup5*F(638)-p4mup5*p5mup6*F(639)-(p1mup5+p2mup5+p3mup5)*(p
     -   4mup6*F(637)+p5mup6*F(641))+p4mup5*p4mup6*F(642)))
       Return
       End
