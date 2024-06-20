c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa2(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa2)
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
       Complex*16 Fa2(25:48)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa2(25) = p4mup3*(p1mup2*p1mup6*(-((p2mup5+p3mup5+p4mup5)*F(116
     -   ))+p5mup5*F(117))+p1mup5*(p1mup6*(p2mup2+p3mup2+p4mup2+p5mup
     -   2)*F(114)+2*p1mup2*((p2mup6+p3mup6+p4mup6)*F(115)+p5mup6*F(1
     -   18))))
       Fa2(26) = p4mup3*(4*p1mup2*(p2mup6*(p3mup5+p4mup5)+p2mup5*(p2mu
     -   p6+p3mup6+p4mup6))*F(112)-2*p1mup5*(p2mup6*(p3mup2+p4mup2)+p
     -   2mup2*(p2mup6+p3mup6+p4mup6))*F(115)+p1mup6*(p2mup5*(p3mup2+
     -   p4mup2)+p2mup2*(p2mup5+p3mup5+p4mup5))*F(116))
       Fa2(27) = p4mup3*(p1mup2*(-2*p2mup6*p5mup5*F(110)+p3mup6*(-2*p5
     -   mup5*F(110)+4*(p3mup5+p4mup5)*F(112))+4*(p3mup5*p4mup6*F(112
     -   )+(p2mup5+p3mup5)*p5mup6*F(113)))+p1mup6*(p3mup2*(p3mup5+p4m
     -   up5)+p2mup5*p5mup2+p3mup5*(p4mup2+p5mup2))*F(116)-p1mup6*(p2
     -   mup2+p3mup2)*p5mup5*F(117)-2*(p1mup5*(p3mup2*(p3mup6+p4mup6)
     -   +p2mup6*p5mup2+p3mup6*(p4mup2+p5mup2))*F(115)+p1mup5*(p2mup2
     -   +p3mup2)*p5mup6*F(118)))
       Fa2(28) = p4mup3*(4*p2mup5*p2mup6*(p3mup2+p4mup2)*F(112)+p1mup2
     -   *(-2*(p4mup6*p5mup5*F(110)-p5mup5*p5mup6*F(111))+4*(p4mup5*p
     -   4mup6*F(112)+p4mup5*p5mup6*F(113)))+p2mup2*(-2*p2mup6*p5mup5
     -   *F(110)+4*(p2mup6*(p2mup5+p3mup5+p4mup5)*F(112)+p2mup5*((p3m
     -   up6+p4mup6)*F(112)+p5mup6*F(113))))+(p4mup2+p5mup2)*(p1mup6*
     -   (p4mup5*F(116)-p5mup5*F(117))-2*p1mup5*(p4mup6*F(115)+p5mup6
     -   *F(118))))
       Fa2(29) = p4mup3*(-(p2mup6*p3mup2*p5mup5*F(110))+2*(p2mup6*(p3m
     -   up5*p4mup2+p3mup2*(p3mup5+p4mup5)+(p2mup5+p3mup5)*p5mup2)*F(
     -   112)+(p2mup5*(p3mup2*(p3mup6+p4mup6)+p3mup6*(p4mup2+p5mup2))
     -   +p4mup5*p4mup6*p6mup2)*F(112)+p2mup5*p3mup2*p5mup6*F(113))+p
     -   2mup2*(-(p3mup6*p5mup5*F(110))+2*(p3mup6*(p3mup5+p4mup5)*F(1
     -   12)+p3mup5*(p4mup6*F(112)+p5mup6*F(113)))))
       Fa2(30) = p4mup3*p6mup2*(4*(p2mup5*p4mup6*F(112)+p3mup5*p4mup6*
     -   F(112)+p2mup5*p5mup6*F(113)+p3mup5*p5mup6*F(113)+p4mup5*p5mu
     -   p6*F(113))-p1mup5*p1mup6*F(114)+p2mup6*(4*(p2mup5+p3mup5+p4m
     -   up5)*F(112)-2*(p5mup5*F(110)+p1mup5*F(115)))+p1mup6*p2mup5*F
     -   (116)+p1mup6*p3mup5*F(116)+p1mup6*p4mup5*F(116)-p1mup6*p5mup
     -   5*F(117)-2*(p4mup6*p5mup5*F(110)-p5mup5*p5mup6*F(111)+p1mup5
     -   *p4mup6*F(115)+p3mup6*(p5mup5*F(110)-2*(p2mup5+p4mup5)*F(112
     -   )+p1mup5*F(115))+p1mup5*p5mup6*F(118)))
       Fa2(31) = p4mup3*p6mup5*(p1mup6*(p2mup2+p3mup2+p4mup2+p5mup2)*F
     -   (119)+2*(p2mup6*(p2mup2+p3mup2)+p2mup2*p3mup6)*F(9572)+p1mup
     -   2*(p1mup6*F(119)+2*(p5mup6*F(9571)+(p2mup6+p3mup6+p4mup6)*F(
     -   9572))))
       Fa2(32) = p4mup3*p6mup5*((p2mup2+p3mup2+p4mup2+p5mup2)*p5mup6*F
     -   (9571)+(p2mup2*p4mup6+p3mup2*(p3mup6+p4mup6)+(p2mup6+p3mup6+
     -   p4mup6)*(p4mup2+p5mup2))*F(9572))
       Fa2(33) = p4mup3*p6mup2*p6mup5*(p1mup6*F(119)+2*(p5mup6*F(9571)
     -   +(p2mup6+p3mup6+p4mup6)*F(9572)))
       Fa2(34) = p4mup3*p6mup6*(p1mup5*(p2mup2+p3mup2+p4mup2+p5mup2)*F
     -   (120)-p2mup2*p2mup5*F(122)+p1mup2*(p1mup5*F(120)+p5mup5*F(12
     -   1)-(p2mup5+p3mup5+p4mup5)*F(122)))
       Fa2(35) = p4mup3*p6mup6*(-(p1mup5*p6mup2*F(120))-(p2mup2+p3mup2
     -   +p4mup2+p5mup2)*p5mup5*F(121)+((p3mup5+p4mup5)*(p2mup2+p3mup
     -   2+p4mup2+p5mup2)+p4mup5*p6mup2+p2mup5*(p3mup2+p4mup2+p5mup2+
     -   p6mup2))*F(122))
       Fa2(36) = p4mup3*p6mup6*(p5mup5*p6mup2*F(121)-p3mup5*p6mup2*F(1
     -   22)-2*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*p6mup5*F(9
     -   569))
       Fa2(37) = mup3mup6*(-(p3mup2*p6mup5*F(126))-p4mup2*p6mup5*F(126
     -   )+p2mup5*(p2mup2+p5mup2)*F(130)-4*(p1mup5*(p3mup2+p4mup2)*F(
     -   124)+p3mup2*p5mup5*F(125)+p4mup2*p5mup5*F(125)-p2mup5*p3mup2
     -   *F(131)-p2mup5*p4mup2*F(131)-p2mup2*p5mup5*F(135)-p5mup2*p5m
     -   up5*F(135)-(p2mup2+p5mup2)*p6mup5*F(137)-p1mup5*(p2mup2+p5mu
     -   p2)*F(139)))+mup3mup5*(4*(p1mup6*(p3mup2+p4mup2)*F(127)+mup2
     -   mup6*F(132)+(p3mup2+p4mup2)*(p5mup6*F(128)-p6mup6*F(129)+p2m
     -   up6*F(134)))+(p2mup2+p5mup2)*p6mup6*F(138)-8*(p2mup6*(p2mup2
     -   +p5mup2)*F(133)-(p2mup2+p5mup2)*p5mup6*F(136)-p1mup6*(p2mup2
     -   +p5mup2)*F(140)))
       Fa2(38) = mup3mup6*(p1mup2+p6mup2)*(p2mup5*F(130)+4*(p5mup5*F(1
     -   35)+p6mup5*F(137)+p1mup5*F(139)))+mup3mup5*(8*p5mup6*p6mup2*
     -   F(136)+p1mup2*(p6mup6*F(138)-8*(p2mup6*F(133)-p5mup6*F(136)-
     -   p1mup6*F(140))))
       Fa2(39) = mup3mup5*p6mup2*(p6mup6*F(138)-8*(p2mup6*F(133)-p1mup
     -   6*F(140)))+4*mup3mup6*(mup2mup5*F(141)+(p3mup5+p4mup5)*((p3m
     -   up2+p4mup2)*F(142)+(p1mup2+p2mup2+p5mup2)*F(143)))
       Fa2(40) = mup3mup5*(p3mup6+p4mup6)*((p3mup2+p4mup2)*F(144)+2*(p
     -   1mup2+p2mup2+p5mup2+p6mup2)*F(145))
       Fa2(41) = p4mup3*(p2mup2*(-(p5mup5*p5mup6*F(146))+p4mup5*p5mup6
     -   *F(147)+2*p4mup6*p5mup5*F(148)+4*p4mup5*p4mup6*F(9574))+p3mu
     -   p2*(-(p5mup5*p5mup6*F(146))+p4mup5*p5mup6*F(147)+2*p4mup6*p5
     -   mup5*F(148)+4*(p3mup5*p3mup6+p4mup5*p4mup6)*F(9574))+(p4mup2
     -   +p5mup2)*(-(p5mup5*p5mup6*F(146))+(p2mup5+p3mup5+p4mup5)*p5m
     -   up6*F(147)+2*(p2mup6+p3mup6+p4mup6)*p5mup5*F(148)+4*((p2mup6
     -   +p3mup6)*p4mup5+(p2mup5+p3mup5+p4mup5)*p4mup6)*F(9574)))
       Fa2(42) = p4mup3*(p3mup2*p3mup5*p5mup6*F(147)+2*p3mup2*p3mup6*p
     -   5mup5*F(148)-p1mup2*p1mup5*p1mup6*F(149)+4*(p3mup2*p3mup6*p4
     -   mup5+p3mup5*(p3mup2*p4mup6+p3mup6*(p4mup2+p5mup2)))*F(9574))
       Fa2(43) = p4mup3*(p1mup2*p1mup6*(-(p5mup5*F(150))+(p2mup5+p3mup
     -   5+p4mup5)*F(152))+p1mup5*(p1mup6*(p2mup2+p3mup2+p4mup2+p5mup
     -   2)*F(149)-2*p1mup2*(p2mup6+p3mup6+p4mup6)*F(151)-p1mup2*p5mu
     -   p6*F(153)))
       Fa2(44) = p4mup3*(2*p1mup5*(p2mup6*(p3mup2+p4mup2)+p2mup2*(p2mu
     -   p6+p3mup6+p4mup6))*F(151)-p1mup6*(p2mup5*(p3mup2+p4mup2)+p2m
     -   up2*(p2mup5+p3mup5+p4mup5))*F(152)+4*p1mup2*(p2mup6*(p3mup5+
     -   p4mup5)+p2mup5*(p2mup6+p3mup6+p4mup6))*F(9574))
       Fa2(45) = p4mup3*(p1mup6*(p2mup2+p3mup2)*p5mup5*F(150)+2*p1mup5
     -   *(p3mup6*p4mup2+p3mup2*(p3mup6+p4mup6)+(p2mup6+p3mup6)*p5mup
     -   2)*F(151)-p1mup6*(p3mup5*p4mup2+p3mup2*(p3mup5+p4mup5)+(p2mu
     -   p5+p3mup5)*p5mup2)*F(152)+p1mup5*(p2mup2+p3mup2)*p5mup6*F(15
     -   3)+p1mup2*((p2mup5+p3mup5)*p5mup6*F(147)+2*(p2mup6+p3mup6)*p
     -   5mup5*F(148)+4*(p3mup6*p4mup5+p3mup5*(p3mup6+p4mup6))*F(9574
     -   )))
       Fa2(46) = p4mup3*((p4mup2+p5mup2)*(p1mup6*p5mup5*F(150)+2*p1mup
     -   5*p4mup6*F(151)-p1mup6*p4mup5*F(152)+p1mup5*p5mup6*F(153))+4
     -   *p2mup5*p2mup6*(p3mup2+p4mup2)*F(9574)+p1mup2*(-(p5mup5*p5mu
     -   p6*F(146))+p4mup5*p5mup6*F(147)+2*p4mup6*p5mup5*F(148)+4*p4m
     -   up5*p4mup6*F(9574))+p2mup2*(p2mup5*p5mup6*F(147)+2*p2mup6*p5
     -   mup5*F(148)+4*(p2mup6*(p3mup5+p4mup5)+p2mup5*(p2mup6+p3mup6+
     -   p4mup6))*F(9574)))
       Fa2(47) = p4mup3*(p2mup5*p3mup2*p5mup6*F(147)+p2mup2*p3mup5*p5m
     -   up6*F(147)+2*(p2mup6*p3mup2*p5mup5*F(148)+p2mup2*p3mup6*p5mu
     -   p5*F(148))+4*(p2mup2*p3mup5*p3mup6+p2mup2*p3mup6*p4mup5+p2mu
     -   p2*p3mup5*p4mup6+p2mup6*(p3mup5*p4mup2+p3mup2*(p3mup5+p4mup5
     -   )+(p2mup5+p3mup5)*p5mup2)+p2mup5*(p3mup2*(p3mup6+p4mup6)+p3m
     -   up6*(p4mup2+p5mup2))+p4mup5*p4mup6*p6mup2)*F(9574))
       Fa2(48) = p4mup3*p6mup2*(p2mup5*p5mup6*F(147)+p3mup5*p5mup6*F(1
     -   47)+p4mup5*p5mup6*F(147)-p1mup5*p1mup6*F(149)+p5mup5*(-(p5mu
     -   p6*F(146))+2*(p2mup6+p3mup6+p4mup6)*F(148)+p1mup6*F(150))+2*
     -   p1mup5*(p2mup6+p3mup6+p4mup6)*F(151)-p1mup6*p2mup5*F(152)-p1
     -   mup6*p3mup5*F(152)-p1mup6*p4mup5*F(152)+p1mup5*p5mup6*F(153)
     -   +4*(p3mup6*p4mup5+p2mup6*(p3mup5+p4mup5)+p3mup5*p4mup6+p2mup
     -   5*(p2mup6+p3mup6+p4mup6))*F(9574))
       Return
       End
