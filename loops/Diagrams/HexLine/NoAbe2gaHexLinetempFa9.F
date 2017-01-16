c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa9(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa9)
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
       Complex*16 Fa9(193:216)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa9(193) = p6mup5*(p2mup2*(p1mup3+p2mup3)*(p5mup6*F(570)+(p1mup
     -   6+p2mup6+p3mup6)*F(571)+p4mup6*F(572))+p1mup2*(p1mup3*(p5mup
     -   6*F(570)+(p2mup6+p3mup6)*F(571)+p4mup6*F(572))+p2mup3*(p5mup
     -   6*F(570)+(p1mup6+p2mup6+p3mup6)*F(571)+p4mup6*F(572))))
       Fa9(194) = (p1mup2+p2mup2)*p6mup3*p6mup5*(p4mup6*F(643)+(p1mup6
     -   +p2mup6+p3mup6)*F(644)-p5mup6*F(645))
       Fa9(195) = p6mup6*((p1mup2+p2mup2)*p2mup3*((p1mup5+p2mup5+p3mup
     -   5)*F(576)+p4mup5*F(577)+p5mup5*F(578))+p1mup3*(p1mup5*p2mup2
     -   *F(576)+(p1mup2+p2mup2)*((p2mup5+p3mup5)*F(576)+p4mup5*F(577
     -   )+p5mup5*F(578))))
       Fa9(196) = (p1mup2+p2mup2)*p6mup6*(p6mup3*((p1mup5+p2mup5+p3mup
     -   5)*F(646)-p5mup5*F(648)-p4mup5*F(649))+2*((p1mup3+p2mup3)*p6
     -   mup5*F(581)-p6mup3*p6mup5*F(650)))
       Fa9(197) = mup5mup6*(p1mup2+p2mup2)*p5mup3
       Fa9(198) = p1mup2*p1mup5*p1mup6*p2mup3*F(591)+p1mup3*(p1mup2*p1
     -   mup6*p2mup5+p1mup5*(p1mup6*p2mup2+p1mup2*(p2mup6+p3mup6)))*F
     -   (591)+p2mup2*p5mup3*(p5mup5*p5mup6*F(631)+p4mup6*p5mup5*F(63
     -   6)-p4mup5*p5mup6*F(639)+p4mup5*p4mup6*F(642))
       Fa9(199) = p2mup2*(p1mup5*p1mup6*p2mup3+p1mup3*p1mup6*p2mup5+p1
     -   mup3*p1mup5*p2mup6+p2mup3*p2mup5*p2mup6)*F(591)+p1mup2*(p2mu
     -   p3*(p1mup6*p2mup5+p1mup5*(p2mup6+p3mup6))*F(591)+p1mup3*(-(p
     -   1mup6*(p5mup5*F(568)+p4mup5*F(590)))+p2mup5*p2mup6*F(591)+((
     -   p1mup6+p2mup6)*p3mup5+p2mup5*p3mup6)*F(591)-p1mup5*(p5mup6*F
     -   (621)+p4mup6*F(628)))-p1mup5*p1mup6*p5mup3*F(632))
       Fa9(200) = p2mup2*(-(p2mup3*(p1mup6*p2mup5+p1mup5*(p2mup6+p3mup
     -   6))*F(591))+p1mup3*(-(p2mup6*p3mup5*F(591))-p1mup5*p3mup6*F(
     -   591)-p2mup5*(p2mup6+p3mup6)*F(591)+p1mup6*(p5mup5*F(568)+p4m
     -   up5*F(590)-p3mup5*F(591))+p1mup5*p5mup6*F(621)+p1mup5*p4mup6
     -   *F(628))+p1mup5*p1mup6*p5mup3*F(632))+p1mup2*(-(p2mup3*p2mup
     -   5*p2mup6*F(591))-p2mup3*p2mup6*p3mup5*F(591)-p2mup3*p2mup5*p
     -   3mup6*F(591)+p1mup5*p2mup3*p5mup6*F(621)+p1mup5*p2mup3*p4mup
     -   6*F(628)+p1mup3*(p5mup5*(p5mup6*F(567)+(p2mup6+p3mup6)*F(568
     -   )+p4mup6*F(569))+p2mup6*p4mup5*F(590)+p3mup6*p4mup5*F(590)-p
     -   3mup5*p3mup6*F(591)-p4mup5*p4mup6*F(596)+p2mup5*p5mup6*F(621
     -   )+p3mup5*p5mup6*F(621)+p4mup5*p5mup6*F(625)+(p2mup5+p3mup5)*
     -   p4mup6*F(628))+p1mup5*p2mup6*p5mup3*F(632)+p1mup5*p3mup6*p5m
     -   up3*F(632)+p1mup5*p4mup6*p5mup3*F(637)+p1mup6*(p2mup3*(p5mup
     -   5*F(568)+p4mup5*F(590)-p3mup5*F(591))+p5mup3*((p2mup5+p3mup5
     -   )*F(632)+p5mup5*F(633)-p4mup5*F(638)))+p1mup5*p5mup3*p5mup6*
     -   F(641))
       Fa9(201) = p1mup3*p2mup2*(p2mup6*(p5mup5*F(568)+p4mup5*F(590))+
     -   p3mup6*(p5mup5*F(568)+p4mup5*F(590)-p3mup5*F(591))-p4mup5*p4
     -   mup6*F(596)+(p2mup5+p3mup5)*p5mup6*F(621)+(p2mup5+p3mup5)*p4
     -   mup6*F(628))+p2mup3*(p1mup5*p2mup2*p5mup6*F(621)+p1mup2*(p2m
     -   up6*(p5mup5*F(568)+p4mup5*F(590))+p3mup6*(p5mup5*F(568)+p4mu
     -   p5*F(590)-p3mup5*F(591))-p4mup5*p4mup6*F(596)+(p2mup5+p3mup5
     -   )*p5mup6*F(621)+p4mup5*p5mup6*F(625))+(p1mup5*p2mup2+p1mup2*
     -   (p2mup5+p3mup5))*p4mup6*F(628))+(p1mup2*p2mup6*(p2mup5+p3mup
     -   5)+p1mup2*p2mup5*p3mup6+p1mup5*p2mup2*(p2mup6+p3mup6))*p5mup
     -   3*F(632)+p1mup6*p2mup2*(p2mup3*(p5mup5*F(568)+p4mup5*F(590)-
     -   p3mup5*F(591))+(p2mup5+p3mup5)*p5mup3*F(632))
       Fa9(202) = p1mup3*p2mup2*(p5mup5*p5mup6*F(567)+p4mup6*p5mup5*F(
     -   569)+p4mup5*p5mup6*F(625))+p2mup2*p5mup3*(p1mup6*p5mup5*F(63
     -   3)+p1mup5*p4mup6*F(637)-p1mup6*p4mup5*F(638)+p1mup5*p5mup6*F
     -   (641))+p1mup2*(p2mup3*p5mup5*(p5mup6*F(567)+p4mup6*F(569))+p
     -   3mup5*p3mup6*p5mup3*F(632)+p5mup3*(p2mup6*p5mup5*F(633)+(p2m
     -   up5+p3mup5)*p4mup6*F(637)-(p2mup6+p3mup6)*p4mup5*F(638)+p2mu
     -   p5*p5mup6*F(641)))
       Fa9(203) = p2mup2*(p2mup3*(p4mup6*p5mup5*F(569)+p2mup6*(p5mup5*
     -   F(568)+p4mup5*F(590)-p3mup5*F(591))+p3mup6*(p5mup5*F(568)+p4
     -   mup5*F(590)-(p2mup5+p3mup5)*F(591))-p4mup5*p4mup6*F(596)+(p2
     -   mup5+p3mup5)*p5mup6*F(621)+p4mup5*p5mup6*F(625)+(p2mup5+p3mu
     -   p5)*p4mup6*F(628))+p5mup3*(p2mup5*(p2mup6+p3mup6)*F(632)+p2m
     -   up5*p4mup6*F(637)+p2mup6*(p3mup5*F(632)-p4mup5*F(638))))-p1m
     -   up2*(p1mup3*p1mup5*p1mup6*F(591)+p5mup3*(p5mup5*(p5mup6*F(63
     -   1)-p3mup6*F(633)+p4mup6*F(636))-p5mup6*(p4mup5*F(639)+p3mup5
     -   *F(641))+p4mup5*p4mup6*F(642)))
       Fa9(204) = p2mup2*(p2mup3*p5mup5*p5mup6*F(567)+p3mup5*p3mup6*p5
     -   mup3*F(632)+p5mup3*((p2mup6+p3mup6)*p5mup5*F(633)+p3mup5*p4m
     -   up6*F(637)-p3mup6*p4mup5*F(638)+(p2mup5+p3mup5)*p5mup6*F(641
     -   )))
       Fa9(205) = p6mup5*(p2mup2*p2mup3*((p1mup6+p2mup6+p3mup6)*F(573)
     -   +p4mup6*F(574)+p5mup6*F(575))+p1mup3*(p1mup6*p2mup2*F(573)+(
     -   p1mup2+p2mup2)*((p2mup6+p3mup6)*F(573)+p4mup6*F(574)+p5mup6*
     -   F(575)))-p2mup2*(p1mup6+p2mup6+p3mup6)*p5mup3*F(644)+p1mup2*
     -   (p2mup3*((p2mup6+p3mup6)*F(573)+p4mup6*F(574)+p5mup6*F(575))
     -   +p1mup6*(p2mup3*F(573)-p5mup3*F(644))-p5mup3*(p4mup6*F(643)+
     -   (p2mup6+p3mup6)*F(644)-p5mup6*F(645))))
       Fa9(206) = p6mup5*(p1mup2*p1mup3*p1mup6*F(573)+p2mup2*p5mup3*(-
     -   (p4mup6*F(643))+p5mup6*F(645)))
       Fa9(207) = p6mup6*(p1mup3*(p1mup5*p2mup2*F(579)+(p1mup2+p2mup2)
     -   *((p2mup5+p3mup5)*F(579)+p4mup5*F(580)+p5mup5*F(647)))+p2mup
     -   2*(p2mup5*p5mup3*F(646)+p1mup5*(p2mup3*F(579)+p5mup3*F(646))
     -   +p2mup3*((p2mup5+p3mup5)*F(579)+p4mup5*F(580)+p5mup5*F(647))
     -   )+p1mup2*(p1mup5*(p2mup3*F(579)+p5mup3*F(646))+p2mup3*((p2mu
     -   p5+p3mup5)*F(579)+p4mup5*F(580)+p5mup5*F(647))+p5mup3*((p2mu
     -   p5+p3mup5)*F(646)-p5mup5*F(648)-p4mup5*F(649))))
       Fa9(208) = p6mup6*(p1mup2*p1mup3*p1mup5*F(579)+p2mup2*p5mup3*(p
     -   3mup5*F(646)-p5mup5*F(648)-p4mup5*F(649)))
       Fa9(209) = -((p1mup2+p2mup2)*p6mup5*p6mup6*(p5mup3*F(650)-2*(p1
     -   mup3+p2mup3)*F(9567)))
       Fa9(210) = (p1mup2+p2mup2)*p6mup3
       Fa9(211) = (p1mup2+p2mup2)*p5mup3
       Fa9(212) = -4*p6mup5*p6mup6*F(615)+mup5mup6*F(661)+p4mup5*p4mup
     -   6*F(662)+(p1mup5+p2mup5+p3mup5)*p4mup6*F(664)+p3mup6*p4mup5*
     -   F(665)+p4mup6*p5mup5*F(666)+p4mup5*p5mup6*F(667)+p3mup6*p6mu
     -   p5*F(668)+p5mup6*p6mup5*F(670)+2*(p3mup6*p5mup5*F(608)+p1mup
     -   5*p6mup6*F(612)+p2mup5*p6mup6*F(612)+p3mup5*p6mup6*F(612)-p4
     -   mup6*p6mup5*F(626)+p4mup5*p6mup6*F(663)+p5mup5*p6mup6*F(669)
     -   +p1mup5*p3mup6*F(671)+p2mup5*p3mup6*F(671)+p3mup5*p3mup6*F(6
     -   71))+p1mup6*(p4mup5*F(665)+p6mup5*F(668)+2*(p5mup5*F(608)+(p
     -   1mup5+p2mup5+p3mup5)*F(671)))+p2mup6*(p4mup5*F(665)+p6mup5*F
     -   (668)+2*(p5mup5*F(608)+(p1mup5+p2mup5+p3mup5)*F(671)))+p1mup
     -   5*p5mup6*F(672)+p2mup5*p5mup6*F(672)+p3mup5*p5mup6*F(672)+p5
     -   mup5*p5mup6*F(673)
       Fa9(213) = p3mup6*p5mup3*F(674)+p5mup3*p6mup6*F(675)-4*(p1mup3*
     -   p3mup6*F(676)+p2mup3*p3mup6*F(676)+p1mup3*p4mup6*F(676)+p2mu
     -   p3*p4mup6*F(676)+p1mup3*p5mup6*F(676)+p2mup3*p5mup6*F(676)+p
     -   1mup3*p6mup6*F(676)+p2mup3*p6mup6*F(676))+p4mup6*p5mup3*F(67
     -   7)+p4mup6*p6mup3*F(678)+p6mup3*p6mup6*F(679)+p3mup6*p6mup3*F
     -   (680)+p1mup6*(p5mup3*F(674)-4*(p1mup3+p2mup3)*F(676)+p6mup3*
     -   F(680))+p2mup6*(p5mup3*F(674)-4*(p1mup3+p2mup3)*F(676)+p6mup
     -   3*F(680))+p5mup3*p5mup6*F(681)+p5mup6*p6mup3*F(682)
       Fa9(214) = p2mup3*p4mup6*((p2mup5+p3mup5)*F(582)-p5mup5*F(583)-
     -   p4mup5*F(594))
       Fa9(215) = p4mup6*(p1mup5*p2mup3*F(582)+p1mup3*((p1mup5+p2mup5+
     -   p3mup5)*F(582)-p5mup5*F(583)-p4mup5*F(594)))
       Fa9(216) = (p1mup3+p2mup3)*p4mup5*(p5mup6*F(595)+(p1mup6+p2mup6
     -   +p3mup6)*F(622))
       Return
       End
