c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa10(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa10)
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
       Complex*16 Fa10(217:240)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa10(217) = (p1mup3+p2mup3)*(p5mup5*p5mup6*F(585)+(p1mup5+p2mup
     -   5+p3mup5)*((p1mup6+p2mup6+p3mup6)*F(586)-p5mup6*F(587))+(p1m
     -   up6+p2mup6+p3mup6)*p5mup5*F(588))
       Fa10(218) = p6mup3*(-(p5mup5*p5mup6*F(631))+p1mup6*p2mup5*F(632
     -   )+p2mup5*p2mup6*F(632)+p1mup6*p3mup5*F(632)+p2mup5*p3mup6*F(
     -   632)+p1mup5*(p1mup6+p2mup6+p3mup6)*F(632)+p1mup5*p4mup6*F(63
     -   7)+p2mup5*p4mup6*F(637)+p3mup5*p4mup6*F(637)-(p1mup6+p2mup6+
     -   p3mup6)*p4mup5*F(638))
       Fa10(219) = p1mup3*(p1mup6+p2mup6)*p6mup5*F(571)+p2mup3*(p5mup6
     -   *p6mup5*F(570)+(p1mup6+p2mup6+p3mup6)*p6mup5*F(571)+p4mup6*p
     -   6mup5*F(572)+p6mup6*((p1mup5+p2mup5+p3mup5)*F(576)+p4mup5*F(
     -   577)+p5mup5*F(578)))-p6mup3*(p2mup6*p3mup5*F(632)+p3mup5*p3m
     -   up6*F(632)+p1mup6*p5mup5*F(633)+p2mup6*p5mup5*F(633)+p3mup6*
     -   p5mup5*F(633)-p4mup6*p5mup5*F(636)+p4mup5*p5mup6*F(639)+(p1m
     -   up5+p2mup5)*p5mup6*F(641)+p3mup5*p5mup6*F(641)-p4mup5*p4mup6
     -   *F(642))
       Fa10(220) = p1mup3*p6mup5*(p5mup6*F(570)+p3mup6*F(571)+p4mup6*F
     -   (572))
       Fa10(221) = p6mup3*p6mup5*(p4mup6*F(643)+(p1mup6+p2mup6+p3mup6)
     -   *F(644)-p5mup6*F(645))
       Fa10(222) = p1mup3*p6mup6*((p1mup5+p2mup5+p3mup5)*F(576)+p4mup5
     -   *F(577)+p5mup5*F(578))
       Fa10(223) = p6mup3*p6mup6*((p1mup5+p2mup5+p3mup5)*F(646)-p5mup5
     -   *F(648)-p4mup5*F(649))
       Fa10(224) = p6mup5*p6mup6*(p1mup3*F(581)-p6mup3*F(650))
       Fa10(225) = p4mup6*(p2mup3*(p5mup5*F(569)-p4mup5*F(596)+(p2mup5
     -   +p3mup5)*F(628))+p5mup3*(-(p5mup5*F(636))+(p2mup5+p3mup5)*F(
     -   637)-p4mup5*F(642)))
       Fa10(226) = p4mup6*(p1mup3*(p5mup5*F(569)-p4mup5*F(596)+(p1mup5
     -   +p2mup5+p3mup5)*F(628))+p1mup5*(p2mup3*F(628)+p5mup3*F(637))
     -   )
       Fa10(227) = p4mup5*(p1mup6*p2mup3*F(590)+p2mup3*p2mup6*F(590)+p
     -   2mup3*p3mup6*F(590)+p1mup3*(p1mup6+p2mup6+p3mup6)*F(590)+p1m
     -   up3*p5mup6*F(625)+p2mup3*p5mup6*F(625)-p1mup6*p5mup3*F(638)-
     -   p2mup6*p5mup3*F(638)-p3mup6*p5mup3*F(638)+p5mup3*p5mup6*F(63
     -   9))
       Fa10(228) = p1mup3*(p5mup5*p5mup6*F(567)+(p1mup6+p2mup6+p3mup6)
     -   *p5mup5*F(568)+(p1mup5+p2mup5+p3mup5)*(-((p1mup6+p2mup6+p3mu
     -   p6)*F(591))+p5mup6*F(621)))+p2mup3*(p5mup5*p5mup6*F(567)+(p1
     -   mup6+p2mup6+p3mup6)*p5mup5*F(568)+(p1mup5+p2mup5+p3mup5)*(-(
     -   (p1mup6+p2mup6+p3mup6)*F(591))+p5mup6*F(621)))+p5mup3*(-(p5m
     -   up5*p5mup6*F(631))+(p2mup5*(p1mup6+p2mup6)+p1mup6*p3mup5+p1m
     -   up5*(p1mup6+p2mup6+p3mup6))*F(632)+p1mup6*p5mup5*F(633)+p1mu
     -   p5*p5mup6*F(641))
       Fa10(229) = p5mup3*((p2mup5+p3mup5)*p3mup6*F(632)+p3mup6*p5mup5
     -   *F(633)+p2mup6*(p3mup5*F(632)+p5mup5*F(633))+(p2mup5+p3mup5)
     -   *p5mup6*F(641))
       Fa10(230) = p1mup3*(p1mup6+p2mup6)*p6mup5*F(573)+p2mup3*(p6mup5
     -   *((p1mup6+p2mup6+p3mup6)*F(573)+p4mup6*F(574)+p5mup6*F(575))
     -   -p6mup6*((p1mup5+p2mup5+p3mup5)*F(579)+p4mup5*F(580)-p5mup5*
     -   F(627)))+p5mup3*p5mup6*p6mup5*F(645)
       Fa10(231) = p6mup5*(p1mup3*(p3mup6*F(573)+p4mup6*F(574)+p5mup6*
     -   F(575))-p5mup3*(p4mup6*F(643)+(p1mup6+p2mup6+p3mup6)*F(644))
     -   )
       Fa10(232) = p6mup6*(p2mup5*p5mup3*F(646)+4*p2mup3*p6mup5*F(9567
     -   ))
       Fa10(233) = p6mup6*(p1mup3*((p1mup5+p2mup5+p3mup5)*F(579)+p4mup
     -   5*F(580)-p5mup5*F(627))+p5mup3*((p1mup5+p3mup5)*F(646)-p5mup
     -   5*F(648)-p4mup5*F(649)))
       Fa10(234) = p6mup5*p6mup6*(p5mup3*F(650)-2*p1mup3*F(9567))
       Fa10(235) = p5mup3*(-2*p5mup5*F(655)+(p1mup5+p2mup5+p3mup5)*F(6
     -   86)+p6mup5*F(688)+p4mup5*F(690))+p6mup3*((p1mup5+p2mup5+p3mu
     -   p5)*F(687)-2*(p5mup5*F(656)-p6mup5*F(689))+p4mup5*F(691))+16
     -   *(p1mup3+p2mup3)*(p1mup5+p2mup5+p3mup5+p4mup5+p5mup5+p6mup5)
     -   *F(9573)
       Fa10(236) = -(p2mup6*p6mup5*F(571))-p3mup6*p6mup5*F(571)-p4mup6
     -   *p6mup5*F(572)-p6mup6*((p1mup5+p2mup5+p3mup5)*F(576)+p4mup5*
     -   F(577)+p5mup5*F(578))-2*p6mup5*p6mup6*F(581)-p1mup5*p4mup6*F
     -   (582)-p2mup5*p4mup6*F(582)-p3mup5*p4mup6*F(582)+p4mup6*p5mup
     -   5*F(583)+p2mup6*p4mup5*F(584)+p3mup6*p4mup5*F(584)+p1mup5*p2
     -   mup6*F(586)+p2mup5*p2mup6*F(586)+p2mup6*p3mup5*F(586)+p1mup5
     -   *p3mup6*F(586)+p2mup5*p3mup6*F(586)+p3mup5*p3mup6*F(586)+p2m
     -   up6*p5mup5*F(588)+p3mup6*p5mup5*F(588)+p1mup6*(-(p6mup5*F(57
     -   1))+p4mup5*F(584)+(p1mup5+p2mup5+p3mup5)*F(586)+p5mup5*F(588
     -   ))+p4mup5*p4mup6*F(594)-p5mup6*(p6mup5*F(570)-p5mup5*F(585)+
     -   (p1mup5+p2mup5+p3mup5)*F(587)+p4mup5*F(595))
       Fa10(237) = -2*(p3mup6*p5mup5*F(568)+p4mup6*p5mup5*F(569)+p3mup
     -   6*p6mup5*F(573)+p4mup6*p6mup5*F(574)+p5mup6*p6mup5*F(575)-p6
     -   mup6*((p1mup5+p2mup5+p3mup5)*F(579)+p4mup5*F(580))+p3mup6*p4
     -   mup5*F(590)+p1mup6*(p5mup5*F(568)+p6mup5*F(573)+p4mup5*F(590
     -   ))+p2mup6*(p5mup5*F(568)+p6mup5*F(573)+p4mup5*F(590))-p1mup6
     -   *(p1mup5+p2mup5+p3mup5)*F(591)-p2mup6*(p1mup5+p2mup5+p3mup5)
     -   *F(591)-p1mup5*p3mup6*F(591)-p2mup5*p3mup6*F(591)-p3mup5*p3m
     -   up6*F(591)-p4mup5*p4mup6*F(596)-p5mup5*p5mup6*F(692))+p1mup5
     -   *p5mup6*F(693)+p2mup5*p5mup6*F(693)+p3mup5*p5mup6*F(693)+p5m
     -   up5*p6mup6*F(694)+p1mup5*p4mup6*F(695)+p2mup5*p4mup6*F(695)+
     -   p3mup5*p4mup6*F(695)+p4mup5*p5mup6*F(696)+8*p6mup5*p6mup6*F(
     -   9567)
       Fa10(238) = p1mup6+p2mup6+p3mup6+p4mup6+p5mup6+p6mup6
       Fa10(239) = p1mup6*p5mup2*F(697)+p2mup6*p5mup2*F(697)+p3mup6*p5
     -   mup2*F(697)+p5mup2*p6mup6*F(698)+p4mup6*p5mup2*F(699)+p4mup6
     -   *p6mup2*F(700)+p6mup2*p6mup6*F(701)+p1mup6*p6mup2*F(702)+p2m
     -   up6*p6mup2*F(702)+p3mup6*p6mup2*F(702)+p5mup2*p5mup6*F(703)+
     -   p5mup6*p6mup2*F(704)
       Fa10(240) = (p1mup2+p2mup2)*((p1mup6+p2mup6+p3mup6)*F(674)+p6mu
     -   p6*F(675)+p4mup6*F(677)+p5mup6*F(681))
       Return
       End
