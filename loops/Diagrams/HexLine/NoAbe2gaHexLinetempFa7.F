c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa7(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa7)
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
       Complex*16 Fa7(145:168)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa7(145) = p1mup2*((p1mup5+p2mup5)*F(539)+2*(p5mup5*F(434)+p6mu
     -   p5*F(540))+(p3mup5+p4mup5)*F(541))+p2mup2*((p1mup5+p2mup5)*F
     -   (539)+2*(p5mup5*F(434)+p6mup5*F(540))+(p3mup5+p4mup5)*F(541)
     -   )+p5mup2*s12*(2*p5mup5*F(448)-4*((p1mup5+p2mup5)*F(452)+(p3m
     -   up5+p4mup5)*F(454))+p6mup5*F(542))
       Fa7(146) = -((p1mup2+p2mup2)*(mup5mup6*F(516)+2*(p5mup5*p6mup6*
     -   F(522)-(p1mup6+p2mup6)*(p5mup5*F(518)-4*((p3mup5+p4mup5)*F(5
     -   20)-p6mup5*F(524)))-p6mup5*p6mup6*F(525)-2*(p5mup5*p5mup6*F(
     -   521)+p3mup5*p6mup6*F(527)+p4mup5*p6mup6*F(527)+p3mup6*p5mup5
     -   *F(528)+p4mup6*p5mup5*F(528))+p1mup5*(p6mup6*F(517)-4*(-(p5m
     -   up6*F(519))+(p3mup6+p4mup6)*F(523)+(p1mup6+p2mup6)*F(531)))+
     -   p2mup5*(p6mup6*F(517)-4*(-(p5mup6*F(519))+(p3mup6+p4mup6)*F(
     -   523)+(p1mup6+p2mup6)*F(531)))+4*(p3mup5*p5mup6*F(526)+p4mup5
     -   *p5mup6*F(526)-p3mup6*p6mup5*F(529)-p4mup6*p6mup5*F(529)+p5m
     -   up6*p6mup5*F(530)-(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(9576))))
     -   )
       Fa7(147) = mup3mup6*(p1mup2+p2mup2)
       Fa7(148) = mup3mup5*(p1mup2+p2mup2)
       Fa7(149) = p1mup2+p2mup2
       Fa7(150) = 2*(p5mup5*p6mup6*F(522)-p6mup5*p6mup6*F(525)+p1mup6*
     -   (-(p5mup5*F(518))+4*((p3mup5+p4mup5)*F(520)-p6mup5*F(524)-(p
     -   1mup5+p2mup5)*F(531)))+p2mup6*(-(p5mup5*F(518))+4*((p3mup5+p
     -   4mup5)*F(520)-p6mup5*F(524)-(p1mup5+p2mup5)*F(531))))+4*p5mu
     -   p5*p5mup6*F(533)+mup5mup6*F(543)+p1mup5*p6mup6*F(544)+p2mup5
     -   *p6mup6*F(544)+p3mup5*p6mup6*F(545)+p4mup5*p6mup6*F(545)+p3m
     -   up6*p5mup5*F(546)+p4mup6*p5mup5*F(546)+8*((p1mup5+p2mup5)*(p
     -   5mup6*F(519)-(p3mup6+p4mup6)*F(523))+(p3mup5+p4mup5)*p5mup6*
     -   F(526)-p3mup6*p6mup5*F(529)-p4mup6*p6mup5*F(529)+p5mup6*p6mu
     -   p5*F(530)-(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(9576))
       Fa7(151) = (p1mup5+p2mup5)*F(547)+p5mup5*F(548)+p6mup5*F(549)+(
     -   p3mup5+p4mup5)*F(550)
       Fa7(152) = p6mup5*F(540)+(p1mup5+p2mup5)*F(551)+p5mup5*F(552)+(
     -   p3mup5+p4mup5)*F(553)
       Fa7(153) = (p1mup6+p2mup6)*F(554)+p5mup6*F(555)+p6mup6*F(556)+(
     -   p3mup6+p4mup6)*F(557)
       Fa7(154) = (p1mup6+p2mup6)*F(558)+p5mup6*F(559)+p6mup6*F(560)+(
     -   p3mup6+p4mup6)*F(561)
       Fa7(155) = mup3mup6*(p5mup5*F(448)-2*(-(p6mup5*F(449))+(p1mup5+
     -   p2mup5)*F(452)+(p3mup5+p4mup5)*F(454)))+mup3mup5*(-(p5mup6*F
     -   (450))-p6mup6*F(451)+(p1mup6+p2mup6)*F(453)+(p3mup6+p4mup6)*
     -   F(455))
       Fa7(156) = -(p5mup6*F(450))+(p1mup6+p2mup6)*F(453)+(p3mup6+p4mu
     -   p6)*F(455)+p6mup6*F(563)
       Fa7(157) = p5mup5*F(448)+p6mup5*F(564)+(p1mup5+p2mup5)*F(565)+(
     -   p3mup5+p4mup5)*F(566)
       Fa7(158) = p5mup5*(-4*((p1mup6+p2mup6+p3mup6)*p6mup3*F(568)+p4m
     -   up6*p6mup3*F(569)+p5mup3*p6mup6*F(578)-p4mup6*p5mup3*F(583)+
     -   p5mup6*(p6mup3*F(567)-p5mup3*F(585))-(p1mup6+p2mup6+p3mup6)*
     -   p5mup3*F(588))+p6mup3*p6mup6*F(593))-4*p5mup3*(p2mup6*p6mup5
     -   *F(571)+p3mup6*p6mup5*F(571)+p4mup6*p6mup5*F(572)+(p1mup5+p2
     -   mup5+p3mup5)*p6mup6*F(576)+p4mup5*p6mup6*F(577)+2*p6mup5*p6m
     -   up6*F(581)+p1mup5*p4mup6*F(582)+p2mup5*p4mup6*F(582)+p3mup5*
     -   p4mup6*F(582)-p2mup6*p4mup5*F(584)-p3mup6*p4mup5*F(584)-p1mu
     -   p5*p2mup6*F(586)-p2mup5*p2mup6*F(586)-p2mup6*p3mup5*F(586)-p
     -   1mup5*p3mup6*F(586)-p2mup5*p3mup6*F(586)-p3mup5*p3mup6*F(586
     -   )+p1mup6*(p6mup5*F(571)-p4mup5*F(584)-(p1mup5+p2mup5+p3mup5)
     -   *F(586))-p4mup5*p4mup6*F(594)+p5mup6*(p6mup5*F(570)+(p1mup5+
     -   p2mup5+p3mup5)*F(587)+p4mup5*F(595)))+p6mup3*(p1mup5*p4mup6*
     -   F(589)+p2mup5*p4mup6*F(589)+p3mup5*p4mup6*F(589)+p1mup5*p5mu
     -   p6*F(592)+p2mup5*p5mup6*F(592)+p3mup5*p5mup6*F(592)-4*(p3mup
     -   6*p6mup5*F(573)+p4mup6*p6mup5*F(574)+p5mup6*p6mup5*F(575)-p6
     -   mup6*((p1mup5+p2mup5+p3mup5)*F(579)+p4mup5*F(580))+p3mup6*p4
     -   mup5*F(590)+p1mup6*(p6mup5*F(573)+p4mup5*F(590))+p2mup6*(p6m
     -   up5*F(573)+p4mup5*F(590))-p1mup6*(p1mup5+p2mup5+p3mup5)*F(59
     -   1)-p2mup6*(p1mup5+p2mup5+p3mup5)*F(591)-p1mup5*p3mup6*F(591)
     -   -p2mup5*p3mup6*F(591)-p3mup5*p3mup6*F(591)-p4mup5*p4mup6*F(5
     -   96))+p4mup5*p5mup6*F(597)+16*p6mup5*p6mup6*F(9567))
       Fa7(159) = (p1mup2+p2mup2)*(2*(p2mup3*p2mup6*F(598)+p2mup3*p3mu
     -   p6*F(598)+p2mup3*p4mup6*F(598)+p2mup3*p5mup6*F(598)+p2mup3*p
     -   6mup6*F(598)+p1mup3*(p1mup6+p2mup6+p3mup6+p4mup6+p5mup6+p6mu
     -   p6)*F(598))+p2mup6*p5mup3*F(599)+p3mup6*p5mup3*F(599)+p5mup3
     -   *p6mup6*F(600)+p4mup6*p5mup3*F(601)+p4mup6*p6mup3*F(602)+p6m
     -   up3*p6mup6*F(603)+p2mup6*p6mup3*F(604)+p3mup6*p6mup3*F(604)+
     -   p1mup6*(2*p2mup3*F(598)+p5mup3*F(599)+p6mup3*F(604))-p5mup3*
     -   p5mup6*F(605)-p5mup6*p6mup3*F(606))
       Fa7(160) = p5mup2*p5mup6*s12*(p6mup5*F(570)-p5mup5*F(585)+(p1mu
     -   p5+p2mup5)*F(587))+p2mup2*(p1mup5*p1mup6*F(610)+p2mup6*(-(p5
     -   mup5*F(608))+p6mup5*F(609)+(p1mup5+p2mup5)*F(610))+p1mup5*p5
     -   mup6*F(611))
       Fa7(161) = p1mup5*p5mup2*p6mup6*s12*F(576)+p2mup5*p5mup2*p6mup6
     -   *s12*F(576)+p5mup2*p5mup5*p6mup6*s12*F(578)+p2mup6*p5mup2*s1
     -   2*(p6mup5*F(571)-(p1mup5+p2mup5)*F(586)-p5mup5*F(588))+p1mup
     -   6*(p5mup2*s12*(p6mup5*F(571)-(p1mup5+p2mup5)*F(586)-p5mup5*F
     -   (588))+p2mup2*(-(p5mup5*F(608))+p6mup5*F(609)+p2mup5*F(610))
     -   )+p2mup2*p2mup5*p5mup6*F(611)-p1mup5*p2mup2*p6mup6*F(612)-p2
     -   mup2*p2mup5*p6mup6*F(612)+p2mup2*p5mup5*p6mup6*F(613)+p2mup2
     -   *p5mup6*p6mup5*F(614)+2*(p5mup2*p6mup5*p6mup6*s12*F(581)+p2m
     -   up2*p6mup5*p6mup6*F(615))+p2mup2*p5mup5*p5mup6*F(616)
       Fa7(162) = p2mup2*p4mup5*(p4mup6*F(617)+p6mup6*F(618)+p5mup6*F(
     -   619))+p1mup2*(8*p6mup5*p6mup6*F(615)+4*(p1mup6*(-(p5mup5*F(6
     -   08))+p6mup5*F(609)+(p1mup5+p2mup5)*F(610))+p2mup6*(-(p5mup5*
     -   F(608))+p6mup5*F(609)+(p1mup5+p2mup5)*F(610))+p1mup5*p5mup6*
     -   F(611)+p2mup5*p5mup6*F(611)-p1mup5*p6mup6*F(612)-p2mup5*p6mu
     -   p6*F(612)+p5mup5*p6mup6*F(613)+p5mup6*p6mup5*F(614)+p5mup5*p
     -   5mup6*F(616))+p4mup5*p4mup6*F(617)+p4mup5*p6mup6*F(618)+p4mu
     -   p5*p5mup6*F(619))
       Fa7(163) = p3mup5*(p5mup2*s12*(p6mup6*F(576)-(p1mup6+p2mup6)*F(
     -   586)+p5mup6*F(587))+(p1mup2+p2mup2)*((p1mup6+p2mup6)*F(610)+
     -   p5mup6*F(611)-p6mup6*F(612)))-p2mup5*p6mup2*s12*(p6mup6*F(57
     -   9)+p1mup6*F(591)-p5mup6*F(621))+p4mup5*((p1mup2+p2mup2)*(p1m
     -   up6+p2mup6)*F(620)+p5mup2*s12*(p6mup6*F(577)+p5mup6*F(595)+(
     -   p1mup6+p2mup6)*F(622)))
       Fa7(164) = (p1mup2+p2mup2)*p4mup6*((p1mup5+p2mup5+p3mup5)*F(623
     -   )+p5mup5*F(624))
       Fa7(165) = p6mup2*(p5mup5*(p5mup6*F(567)+p2mup6*F(568))-p6mup6*
     -   (p3mup5*F(579)+p4mup5*F(580))+p1mup6*p4mup5*F(590)+p2mup6*p4
     -   mup5*F(590)-p1mup6*p3mup5*F(591)-p2mup6*(p2mup5+p3mup5)*F(59
     -   1)+p3mup5*p5mup6*F(621)+p4mup5*p5mup6*F(625))
       Fa7(166) = p4mup6*p5mup2*p6mup5*s12*F(572)-p1mup5*p6mup2*p6mup6
     -   *s12*F(579)+p1mup5*p4mup6*p5mup2*s12*F(582)+p2mup5*p4mup6*p5
     -   mup2*s12*F(582)-p4mup6*p5mup2*p5mup5*s12*F(583)-p1mup5*p2mup
     -   6*p6mup2*s12*F(591)+p1mup6*p6mup2*s12*(p5mup5*F(568)-p1mup5*
     -   F(591))+p3mup6*(p5mup2*s12*(p6mup5*F(571)-(p1mup5+p2mup5)*F(
     -   586)-p5mup5*F(588))-(p1mup2+p2mup2)*(p5mup5*F(608)-p6mup5*F(
     -   609)-(p1mup5+p2mup5)*F(610)-p4mup5*F(620)))+p1mup5*p5mup6*p6
     -   mup2*s12*F(621)+p1mup2*p4mup6*p6mup5*F(626)+p2mup2*p4mup6*p6
     -   mup5*F(626)+p5mup5*p6mup2*p6mup6*s12*F(627)
       Fa7(167) = p3mup6*(p5mup5*p6mup2*F(568)+p6mup2*p6mup5*F(573)-p3
     -   mup5*p5mup2*F(586)-(p1mup5+p2mup5)*p6mup2*F(591))+p6mup2*((p
     -   1mup6+p2mup6)*p6mup5*F(573)+p5mup6*p6mup5*F(575)+p4mup6*(p5m
     -   up5*F(569)+p6mup5*F(574)+(p1mup5+p2mup5)*F(628)))
       Fa7(168) = p4mup5*s12*(-(p4mup6*(p5mup2*F(594)+p6mup2*F(596)))+
     -   p3mup6*(p6mup2*F(590)+p5mup2*F(622)))+p3mup5*(-(p3mup6*p6mup
     -   2*s12*F(591))+(p1mup2+p2mup2)*p3mup6*F(610)+p4mup6*s12*(p5mu
     -   p2*F(582)+p6mup2*F(628)))
       Return
       End
