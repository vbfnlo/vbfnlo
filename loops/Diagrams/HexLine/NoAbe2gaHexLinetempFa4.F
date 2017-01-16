c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa4(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa4)
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
       Complex*16 Fa4(73:96)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa4(73) = -16*(p2mup5*p3mup2*p3mup6*F(112)+p2mup5*p3mup6*p4mup2
     -   *F(112)+p2mup5*p3mup2*p4mup6*F(112)+p2mup5*p4mup2*p4mup6*F(1
     -   12)+p2mup5*p3mup6*p5mup2*F(112)+p2mup5*p4mup6*p5mup2*F(112))
     -   +4*(p1mup5*p1mup6*p6mup2*F(114)+p1mup6*p5mup5*p6mup2*F(117)-
     -   p1mup5*p6mup2*p6mup6*F(120)-p5mup5*p6mup2*p6mup6*F(121))+mup
     -   2mup5*p3mup6*F(279)+mup2mup5*p4mup6*F(279)+8*(p3mup2*p3mup6*
     -   p5mup5*F(110)+p3mup6*p4mup2*p5mup5*F(110)+p3mup2*p4mup6*p5mu
     -   p5*F(110)+p4mup2*p4mup6*p5mup5*F(110)+p3mup6*p5mup2*p5mup5*F
     -   (110)+p4mup6*p5mup2*p5mup5*F(110)+p1mup5*p3mup2*p3mup6*F(115
     -   )+p1mup5*p3mup6*p4mup2*F(115)+p1mup5*p3mup2*p4mup6*F(115)+p1
     -   mup5*p4mup2*p4mup6*F(115)+p1mup5*p3mup6*p5mup2*F(115)+p1mup5
     -   *p4mup6*p5mup2*F(115)+p1mup5*p2mup6*p6mup2*F(115)+p1mup2*p3m
     -   up6*(p5mup5*F(110)-2*(p2mup5+p4mup5)*F(112)+p1mup5*F(115))+p
     -   2mup2*p3mup6*(p5mup5*F(110)-2*(p2mup5+p4mup5)*F(112)+p1mup5*
     -   F(115))+p1mup5*p5mup6*p6mup2*F(118)-(p3mup6+p4mup6)*(p1mup2+
     -   p2mup2+p3mup2+p4mup2+p5mup2)*p6mup5*F(9572))
       Fa4(74) = p6mup2*p6mup5*(p1mup6*F(119)+2*(p5mup6*F(9571)+p2mup6
     -   *F(9572)))+p4mup6*(4*p4mup2*(p3mup5+p4mup5)*F(112)-2*p6mup2*
     -   (p5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mup5*F(9572))
     -   )+p3mup6*(4*(p4mup2*p4mup5+p3mup5*p5mup2)*F(112)-2*p6mup2*(p
     -   5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mup5*F(9572)))
       Fa4(75) = p1mup2*p3mup5*p3mup6+p2mup2*p3mup5*p3mup6+p3mup2*p3mu
     -   p5*p3mup6+p3mup5*p3mup6*p4mup2+p3mup2*p3mup6*p4mup5+p3mup2*p
     -   3mup5*p4mup6+p3mup2*p4mup5*p4mup6+p3mup6*p4mup5*p5mup2+p3mup
     -   5*p4mup6*p5mup2+p4mup5*p4mup6*p5mup2+(p3mup5+p4mup5)*(p3mup6
     -   +p4mup6)*p6mup2
       Fa4(76) = -(p5mup2*p5mup5*p5mup6*F(146))+p2mup5*p5mup2*p5mup6*F
     -   (147)-p1mup5*p1mup6*p2mup2*F(149)+p1mup5*p2mup2*p5mup6*F(153
     -   )+p1mup5*p5mup2*p5mup6*F(153)+p5mup2*p5mup6*p6mup5*F(155)+2*
     -   (p1mup5*p2mup2*p2mup6*F(151)-p2mup2*p2mup6*p5mup5*F(282))-4*
     -   (p2mup2*p2mup6*p6mup5*F(9570)-p2mup2*p2mup5*p2mup6*F(9574))+
     -   p3mup2*(-(p1mup5*p1mup6*F(149))+p1mup5*p5mup6*F(153)+p1mup6*
     -   p6mup5*F(154)+p5mup6*p6mup5*F(155)-p1mup5*p6mup6*F(156)-p5mu
     -   p5*(p5mup6*F(146)-p1mup6*F(150)+p6mup6*F(281)+2*p2mup6*F(282
     -   ))+2*(p1mup5*p2mup6*F(151)+p6mup5*p6mup6*F(9568))-4*p2mup6*p
     -   6mup5*F(9570)+p2mup5*(p5mup6*F(147)-p1mup6*F(152)-2*p6mup6*F
     -   (9570)+4*p2mup6*F(9574)))+p4mup2*(-(p1mup5*p1mup6*F(149))+p1
     -   mup5*p5mup6*F(153)+p1mup6*p6mup5*F(154)+p5mup6*p6mup5*F(155)
     -   -p1mup5*p6mup6*F(156)-p5mup5*(p5mup6*F(146)-p1mup6*F(150)+p6
     -   mup6*F(281)+2*p2mup6*F(282))+2*(p1mup5*p2mup6*F(151)+p6mup5*
     -   p6mup6*F(9568))-4*p2mup6*p6mup5*F(9570)+p2mup5*(p5mup6*F(147
     -   )-p1mup6*F(152)-2*p6mup6*F(9570)+4*p2mup6*F(9574)))
       Fa4(77) = -2*mup2mup6*p5mup5*F(188)+mup5mup6*p6mup2*F(280)+mup2
     -   mup6*p1mup5*F(283)+mup2mup6*p6mup5*F(284)+mup2mup6*p2mup5*F(
     -   285)+8*(p2mup6*p5mup2*p5mup5*F(282)-p5mup2*p6mup5*p6mup6*F(9
     -   568)+p2mup5*p5mup2*p6mup6*F(9570))-4*(p1mup6*p5mup2*p5mup5*F
     -   (150)-p1mup6*p2mup5*p5mup2*F(152)+p1mup6*p5mup2*p6mup5*F(154
     -   )-p1mup5*p5mup2*(p1mup6*F(149)-2*p2mup6*F(151)+p6mup6*F(156)
     -   )-p5mup2*p5mup5*p6mup6*F(281)-p2mup2*(p1mup5*p6mup6*F(156)+p
     -   5mup5*(p5mup6*F(146)-p1mup6*F(150)+p6mup6*F(281))-p6mup5*(p1
     -   mup6*F(154)+p5mup6*F(155)+2*p6mup6*F(9568))+p2mup5*(-(p5mup6
     -   *F(147))+p1mup6*F(152)+2*p6mup6*F(9570))))+16*(p2mup6*p5mup2
     -   *p6mup5*F(9570)-p2mup5*p2mup6*p5mup2*F(9574))
       Fa4(78) = p4mup5*(-(mup2mup6*F(285))+4*p2mup2*(p5mup6*F(147)-2*
     -   p6mup6*F(9570)+4*p4mup6*F(9574)))+4*p1mup2*(p2mup5*p5mup6*F(
     -   147)+p4mup5*p5mup6*F(147)-p1mup5*p1mup6*F(149)-p1mup6*p2mup5
     -   *F(152)+p1mup5*p5mup6*F(153)+p1mup6*p6mup5*F(154)+p5mup6*p6m
     -   up5*F(155)-p1mup5*p6mup6*F(156)-p5mup5*(p5mup6*F(146)-p1mup6
     -   *F(150)+p6mup6*F(281)+2*p2mup6*F(282))+2*(p1mup5*p2mup6*F(15
     -   1)+p6mup5*p6mup6*F(9568)-p2mup5*p6mup6*F(9570)-p4mup5*p6mup6
     -   *F(9570))-4*(p2mup6*p6mup5*F(9570)-(p2mup5*p2mup6+p4mup5*p4m
     -   up6)*F(9574)))
       Fa4(79) = -4*(p1mup2*p3mup5*p5mup6*F(147)+p2mup2*p3mup5*p5mup6*
     -   F(147)+p3mup2*p3mup5*p5mup6*F(147)+p3mup5*p4mup2*p5mup6*F(14
     -   7)+p3mup2*p4mup5*p5mup6*F(147)+p4mup2*p4mup5*p5mup6*F(147)+p
     -   3mup5*p5mup2*p5mup6*F(147)+p4mup5*p5mup2*p5mup6*F(147)+p2mup
     -   5*p5mup6*p6mup2*F(147)-p1mup2*p1mup6*p3mup5*F(152)-p1mup6*p2
     -   mup2*p3mup5*F(152)-p1mup6*p3mup2*p3mup5*F(152)-p1mup6*p3mup5
     -   *p4mup2*F(152)-p1mup2*p1mup6*p4mup5*F(152)-p1mup6*p2mup2*p4m
     -   up5*F(152)-p1mup6*p3mup2*p4mup5*F(152)-p1mup6*p4mup2*p4mup5*
     -   F(152)-p1mup6*p3mup5*p5mup2*F(152)-p1mup6*p4mup5*p5mup2*F(15
     -   2)-p1mup6*p2mup5*p6mup2*F(152))+mup2mup6*p3mup5*F(285)+mup2m
     -   up5*p1mup6*F(286)+mup2mup5*p5mup6*F(287)+mup2mup5*p6mup6*F(2
     -   88)+mup2mup5*p2mup6*F(289)+8*(p1mup2*p3mup5*p6mup6*F(9570)+p
     -   2mup2*p3mup5*p6mup6*F(9570)+p3mup2*p3mup5*p6mup6*F(9570)+p3m
     -   up5*p4mup2*p6mup6*F(9570)+p3mup2*p4mup5*p6mup6*F(9570)+p4mup
     -   2*p4mup5*p6mup6*F(9570)+p3mup5*p5mup2*p6mup6*F(9570)+p4mup5*
     -   p5mup2*p6mup6*F(9570)+p2mup5*p6mup2*p6mup6*F(9570))-16*p2mup
     -   6*(p3mup5+p4mup5)*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2)*F(957
     -   4)
       Fa4(80) = (p1mup2+p2mup2)*p4mup6*(p1mup5*F(151)-p5mup5*F(282)+2
     -   *(p2mup5+p3mup5)*F(9574))
       Fa4(81) = p6mup2*(p5mup5*(p5mup6*F(146)+2*p2mup6*F(282))-(p3mup
     -   5+p4mup5)*(p5mup6*F(147)-p1mup6*F(152)-2*p6mup6*F(9570))-4*p
     -   2mup6*(p2mup5+p3mup5+p4mup5)*F(9574))
       Fa4(82) = -4*(p1mup6*p5mup5*p6mup2*F(150)-p1mup5*(-2*(p1mup2*p3
     -   mup6+p2mup2*p3mup6+(p3mup6+p4mup6)*(p3mup2+p4mup2+p5mup2))*F
     -   (151)+p6mup2*(p1mup6*F(149)-2*p2mup6*F(151)-p5mup6*F(153)+p6
     -   mup6*F(156)))-p5mup5*p6mup2*p6mup6*F(281))+8*(p1mup2*p3mup6*
     -   p5mup5*F(282)+p2mup2*p3mup6*p5mup5*F(282)+p3mup2*p3mup6*p5mu
     -   p5*F(282)+p3mup6*p4mup2*p5mup5*F(282)+p3mup2*p4mup6*p5mup5*F
     -   (282)+p4mup2*p4mup6*p5mup5*F(282)+p3mup6*p5mup2*p5mup5*F(282
     -   )+p4mup6*p5mup2*p5mup5*F(282))+mup2mup5*p3mup6*F(289)+mup2mu
     -   p5*p4mup6*F(289)+16*(p1mup2*p3mup6*p6mup5*F(9570)+p2mup2*p3m
     -   up6*p6mup5*F(9570)+p3mup2*p3mup6*p6mup5*F(9570)+p3mup6*p4mup
     -   2*p6mup5*F(9570)+p1mup2*p4mup6*p6mup5*F(9570)+p2mup2*p4mup6*
     -   p6mup5*F(9570)+p3mup2*p4mup6*p6mup5*F(9570)+p4mup2*p4mup6*p6
     -   mup5*F(9570)+p3mup6*p5mup2*p6mup5*F(9570)+p4mup6*p5mup2*p6mu
     -   p5*F(9570)-(p1mup2*p3mup6*(p2mup5+p4mup5)+p2mup2*p3mup6*(p2m
     -   up5+p4mup5)+p2mup5*(p3mup6+p4mup6)*(p3mup2+p4mup2+p5mup2))*F
     -   (9574))
       Fa4(83) = p6mup2*(p1mup6*p6mup5*F(154)+p5mup6*p6mup5*F(155)+2*(
     -   p1mup5*(p3mup6+p4mup6)*F(151)-(p3mup6+p4mup6)*p5mup5*F(282))
     -   -4*(p2mup6+p3mup6+p4mup6)*p6mup5*F(9570))+4*(p4mup2*(p3mup5+
     -   p4mup5)*p4mup6+p2mup5*p4mup6*p6mup2+p3mup6*(p4mup2*p4mup5+p3
     -   mup5*p5mup2+p2mup5*p6mup2))*F(9574)
       Fa4(84) = mup3mup6*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)
       Fa4(85) = (p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*(p1mup5*F
     -   (78)+p5mup5*F(80)+p6mup5*F(82)-(p2mup5+p3mup5+p4mup5)*F(83))
     -   -mup2mup5*F(84)
       Fa4(86) = mup3mup5*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)
       Fa4(87) = p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2
       Fa4(88) = -(p2mup5*(p1mup2+p2mup2+p5mup2+p6mup2)*F(292))-(p3mup
     -   2+p4mup2)*(2*(p1mup5*F(88)+(p3mup5+p4mup5)*F(106)-p5mup5*F(2
     -   90))-p6mup5*F(291)-p2mup5*F(293))+2*(mup2mup5*F(96)+(p1mup2+
     -   p2mup2+p5mup2+p6mup2)*(p5mup5*F(99)+p6mup5*F(294)+(p3mup5+p4
     -   mup5)*F(295)+p1mup5*F(296))*P(59))
       Fa4(89) = p4mup2*p6mup5*F(298)+p2mup5*(p1mup2+p2mup2+p5mup2+p6m
     -   up2)*F(299)+p2mup5*p4mup2*F(300)+p3mup2*(-2*((p3mup5+p4mup5)
     -   *F(142)-p5mup5*F(297))+p6mup5*F(298)+p2mup5*F(300))+p1mup2*p
     -   6mup5*F(302)+p2mup2*p6mup5*F(302)+p5mup2*p6mup5*F(302)+p6mup
     -   2*p6mup5*F(302)+2*(mup2mup5*F(132)-p3mup5*p4mup2*F(142)-p4mu
     -   p2*p4mup5*F(142)+p4mup2*p5mup5*F(297)+p1mup2*p5mup5*F(301)+p
     -   2mup2*p5mup5*F(301)+p5mup2*p5mup5*F(301)+p5mup5*p6mup2*F(301
     -   )+(p3mup5+p4mup5)*(p1mup2+p2mup2+p5mup2+p6mup2)*F(303)+p1mup
     -   5*((p3mup2+p4mup2)*F(124)+(p1mup2+p2mup2+p5mup2+p6mup2)*F(30
     -   4)))
       Fa4(90) = mup5mup6*(p3mup2+p4mup2)
       Fa4(91) = -2*p2mup5*p5mup2*p5mup6*F(100)+(p3mup2+p4mup2)*(p1mup
     -   6*(p5mup5*F(307)-2*p2mup5*F(310)+4*p6mup5*F(311))-p2mup5*p6m
     -   up6*F(313)+p2mup6*p5mup5*F(314)+p2mup6*p6mup5*F(315)-2*(p2mu
     -   p5*p5mup6*F(312)+p5mup5*p5mup6*F(316))-p5mup5*p6mup6*F(317)+
     -   p6mup5*p6mup6*F(318)-4*p5mup6*p6mup5*F(319)+p1mup5*(-(p6mup6
     -   *F(306))+p2mup6*F(308)-4*(p5mup6*F(309)-p1mup6*F(320)))+p2mu
     -   p5*p2mup6*F(321))+p2mup2*p2mup6*(p5mup5*F(99)-p1mup5*F(103)+
     -   p6mup5*F(294)+5*p2mup5*F(9572))
       Fa4(92) = 8*p2mup5*(p2mup2*p5mup6*F(100)+p1mup6*(p2mup2+p5mup2)
     -   *F(104))+p2mup5*(p2mup2+p5mup2)*p6mup6*F(323)+mup2mup6*(-2*(
     -   p5mup5*s234*F(99)-p1mup5*s234*F(103)-p2mup5*F(322))+p6mup5*s
     -   234*F(324))-4*p2mup6*p5mup2*(p5mup5*F(99)-p1mup5*F(103)+p6mu
     -   p5*F(294)+5*p2mup5*F(9572))
       Fa4(93) = 2*mup2mup6*p4mup5*s234*F(107)+p1mup2*(p2mup5*(8*(p5mu
     -   p6*F(100)+p1mup6*F(104))+p6mup6*F(323))-4*p2mup6*(p5mup5*F(9
     -   9)-p1mup5*F(103)+p6mup5*F(294)+5*p2mup5*F(9572)))
       Fa4(94) = 4*p2mup6*(p3mup5+p4mup5)*(p1mup2+p2mup2+p5mup2)*F(107
     -   )+p2mup5*p6mup2*(8*(p5mup6*F(100)+p1mup6*F(104))+p6mup6*F(32
     -   3))+mup2mup5*(4*(p5mup6*s234*F(100)+p1mup6*s234*F(104)+p2mup
     -   6*F(325))+p6mup6*s234*F(326))+2*(mup2mup6*p3mup5*s234*F(107)
     -   +2*(p3mup2+p4mup2)*(p3mup5+p4mup5)*(4*(p1mup6*F(327)+p5mup6*
     -   F(328))+p2mup6*F(329)))+(p3mup2+p4mup2)*(p3mup5+p4mup5)*p6mu
     -   p6*F(330)
       Fa4(95) = (p1mup2+p2mup2)*p2mup5*p4mup6
       Fa4(96) = p2mup6*p6mup2*(-(p5mup5*F(99))+(p3mup5+p4mup5)*F(107)
     -   -5*p2mup5*F(9572))
       Return
       End
