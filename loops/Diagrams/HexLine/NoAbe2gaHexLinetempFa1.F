c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa1(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa1)
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
       Complex*16 Fa1(1:24)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa1(1) = mup3mup6*(p1mup5*F(1)+p5mup5*F(3)+p6mup5*F(5)+p2mup5*F
     -   (7)-2*(p3mup5+p4mup5)*F(9))+mup3mup5*(p5mup6*F(4)+2*(p1mup6*
     -   F(2)+p6mup6*F(6))+p2mup6*F(8)+(p3mup6+p4mup6)*F(10))
       Fa1(2) = -(p1mup6*p4mup3*(p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*F
     -   (11))+p2mup2*p4mup3*p5mup6*F(12)+p3mup2*p4mup3*p5mup6*F(12)+
     -   p4mup2*p4mup3*p5mup6*F(12)+p4mup3*p5mup2*p5mup6*F(12)+p4mup3
     -   *p5mup6*p6mup2*F(12)-p2mup2*p4mup3*p6mup6*F(13)-p3mup2*p4mup
     -   3*p6mup6*F(13)-p4mup2*p4mup3*p6mup6*F(13)-p4mup3*p5mup2*p6mu
     -   p6*F(13)-p4mup3*p6mup2*p6mup6*F(13)+p2mup2*p2mup6*p4mup3*F(1
     -   4)+p2mup6*p3mup2*p4mup3*F(14)+p2mup2*p3mup6*p4mup3*F(14)+p3m
     -   up2*p3mup6*p4mup3*F(14)+p2mup6*p4mup2*p4mup3*F(14)+p3mup6*p4
     -   mup2*p4mup3*F(14)+p2mup2*p4mup3*p4mup6*F(14)+p3mup2*p4mup3*p
     -   4mup6*F(14)+p4mup2*p4mup3*p4mup6*F(14)+p2mup6*p4mup3*p5mup2*
     -   F(14)+p3mup6*p4mup3*p5mup2*F(14)+p4mup3*p4mup6*p5mup2*F(14)+
     -   p2mup6*p4mup3*p6mup2*F(14)+p3mup6*p4mup3*p6mup2*F(14)+p4mup3
     -   *p4mup6*p6mup2*F(14)-mup3mup6*p2mup2*F(15)-mup3mup6*p5mup2*F
     -   (15)-mup3mup6*p6mup2*F(15)+p1mup2*(p4mup3*(-(p1mup6*F(11))+p
     -   5mup6*F(12)-p6mup6*F(13)+(p2mup6+p3mup6+p4mup6)*F(14))-mup3m
     -   up6*F(15))-mup3mup6*(p3mup2+p4mup2)*F(16)
       Fa1(3) = 2*mup5mup6*(p2mup2+p5mup2)*F(17)+mup5mup6*(p3mup2+p4mu
     -   p2)*F(18)
       Fa1(4) = p3mup2*p5mup5*p6mup6*F(33)+p4mup2*p5mup5*p6mup6*F(33)+
     -   p3mup2*p5mup6*p6mup5*F(34)+p4mup2*p5mup6*p6mup5*F(34)+p2mup2
     -   *p2mup6*(p6mup5*F(21)+4*(p5mup5*F(20)-p2mup5*F(37)))+4*(-(p1
     -   mup6*p2mup5*(p3mup2+p4mup2)*F(26))+p1mup6*p3mup2*p5mup5*F(27
     -   )+p1mup6*p4mup2*p5mup5*F(27)+p1mup6*p3mup2*p6mup5*F(28)+p1mu
     -   p6*p4mup2*p6mup5*F(28)-p2mup5*p3mup2*p5mup6*F(29)-p2mup5*p4m
     -   up2*p5mup6*F(29)+p2mup5*p3mup2*p6mup6*F(30)+p2mup5*p4mup2*p6
     -   mup6*F(30)+p2mup6*p3mup2*p5mup5*F(31)+p2mup6*p4mup2*p5mup5*F
     -   (31)+p2mup6*p3mup2*p6mup5*F(32)+p2mup6*p4mup2*p6mup5*F(32)+p
     -   3mup2*p6mup5*p6mup6*F(35)+p4mup2*p6mup5*p6mup6*F(35)-p2mup5*
     -   p5mup2*p5mup6*F(36)-p2mup5*p2mup6*p3mup2*F(38)-p2mup5*p2mup6
     -   *p4mup2*F(38)-p3mup2*p5mup5*p5mup6*F(39)-p4mup2*p5mup5*p5mup
     -   6*F(39)+p5mup2*p5mup5*p5mup6*F(41)+p1mup5*((p3mup2+p4mup2)*(
     -   p1mup6*F(22)+p2mup6*F(23)+p5mup6*F(24)+p6mup6*F(25))-p5mup2*
     -   p5mup6*F(42)+p2mup2*(p2mup6*F(19)+p1mup6*F(40)-p5mup6*F(42))
     -   ))+p5mup2*p5mup6*p6mup5*F(43)
       Fa1(5) = 2*mup5mup6*p6mup2*F(17)+p2mup6*p5mup2*(p6mup5*F(21)+4*
     -   (p5mup5*F(20)-p2mup5*F(37)))+p2mup2*p5mup6*p6mup5*F(43)+mup2
     -   mup6*p5mup5*F(47)+mup2mup6*p6mup5*F(48)+mup2mup6*p2mup5*F(49
     -   )+p1mup5*(mup2mup6*F(46)+4*(p2mup6*p5mup2*F(19)+p1mup6*p5mup
     -   2*F(40)+(p2mup2+p5mup2)*p6mup6*F(50)))+4*(-(p2mup2*p2mup5*p5
     -   mup6*F(36))+p2mup2*p5mup5*p5mup6*F(41)+p1mup6*p2mup2*p2mup5*
     -   F(44)+p1mup6*p2mup5*p5mup2*F(44)+p2mup2*p2mup5*p6mup6*F(45)+
     -   p2mup5*p5mup2*p6mup6*F(45)-p1mup6*p2mup2*p5mup5*F(51)-p1mup6
     -   *p5mup2*p5mup5*F(51)-p1mup6*p2mup2*p6mup5*F(52)-p1mup6*p5mup
     -   2*p6mup5*F(52)+p2mup2*p5mup5*p6mup6*F(53)+p5mup2*p5mup5*p6mu
     -   p6*F(53))+8*(p2mup2+p5mup2)*p6mup5*p6mup6*F(54)
       Fa1(6) = p1mup2*(p2mup6*(p6mup5*F(21)+4*(p5mup5*F(20)-p2mup5*F(
     -   37)))+p5mup6*p6mup5*F(43)-4*(p2mup5*p5mup6*F(36)-p5mup5*p5mu
     -   p6*F(41)-p1mup6*p2mup5*F(44)-p2mup5*p6mup6*F(45)-p1mup5*(p2m
     -   up6*F(19)+p1mup6*F(40)-p5mup6*F(42)+p6mup6*F(50))+p1mup6*p5m
     -   up5*F(51)+p1mup6*p6mup5*F(52)-p5mup5*p6mup6*F(53))+8*p6mup5*
     -   p6mup6*F(54))+p4mup5*(2*mup2mup6*F(9)-(p1mup2+p2mup2)*(4*p5m
     -   up6*F(55)-p6mup6*F(56)+16*p4mup6*F(57)))
       Fa1(7) = p1mup2*p3mup5*p6mup6*F(56)+p2mup2*p3mup5*p6mup6*F(56)+
     -   p3mup5*p5mup2*p6mup6*F(56)+p4mup5*p5mup2*p6mup6*F(56)+2*(mup
     -   2mup6*p3mup5*F(9)+mup2mup5*p1mup6*F(58))+mup2mup5*p5mup6*F(5
     -   9)+mup2mup5*p6mup6*F(60)+mup2mup5*p2mup6*F(61)+p3mup2*p3mup5
     -   *p6mup6*F(66)+p3mup5*p4mup2*p6mup6*F(66)+p3mup2*p4mup5*p6mup
     -   6*F(66)+p4mup2*p4mup5*p6mup6*F(66)-4*(-(p2mup5*p6mup2*(-(p5m
     -   up6*F(36))+p1mup6*F(44)+p6mup6*F(45)))+p1mup2*p3mup5*p5mup6*
     -   F(55)+p2mup2*p3mup5*p5mup6*F(55)+p3mup5*p5mup2*p5mup6*F(55)+
     -   p4mup5*p5mup2*p5mup6*F(55)+p1mup2*p2mup6*p3mup5*F(62)+p2mup2
     -   *p2mup6*p3mup5*F(62)+p1mup2*p2mup6*p4mup5*F(62)+p2mup2*p2mup
     -   6*p4mup5*F(62)+p2mup6*p3mup5*p5mup2*F(62)+p2mup6*p4mup5*p5mu
     -   p2*F(62)+p3mup2*p3mup5*p5mup6*F(63)+p3mup5*p4mup2*p5mup6*F(6
     -   3)+p3mup2*p4mup5*p5mup6*F(63)+p4mup2*p4mup5*p5mup6*F(63)+p2m
     -   up6*p3mup2*p3mup5*F(64)+p2mup6*p3mup5*p4mup2*F(64)+p2mup6*p3
     -   mup2*p4mup5*F(64)+p2mup6*p4mup2*p4mup5*F(64)+p1mup6*p3mup2*p
     -   3mup5*F(65)+p1mup6*p3mup5*p4mup2*F(65)+p1mup6*p3mup2*p4mup5*
     -   F(65)+p1mup6*p4mup2*p4mup5*F(65)-p1mup6*(p3mup5+p4mup5)*(p1m
     -   up2+p2mup2+p5mup2)*F(67))
       Fa1(8) = (p1mup2+p2mup2)*p4mup6*(4*p3mup5*F(57)+p2mup5*F(68)+2*
     -   (p1mup5*F(69)-p5mup5*F(70)))
       Fa1(9) = p6mup2*(4*(p5mup5*p5mup6*F(41)+p2mup6*(p5mup5*F(20)-p2
     -   mup5*F(37)-(p3mup5+p4mup5)*F(62)))-(p3mup5+p4mup5)*(-(p6mup6
     -   *F(56))+4*(p5mup6*F(55)-p1mup6*F(67))))
       Fa1(10) = 16*(p1mup2+p2mup2)*p3mup6*p4mup5*F(57)-mup2mup5*p3mup
     -   6*F(71)-mup2mup5*p4mup6*F(71)+4*(p1mup6*p5mup5*p6mup2*F(51)-
     -   p5mup5*p6mup2*p6mup6*F(53)+p2mup5*((p1mup2+p2mup2)*p3mup6+(p
     -   3mup6+p4mup6)*p5mup2)*F(68)+p3mup2*p3mup6*p5mup5*F(72)+p3mup
     -   6*p4mup2*p5mup5*F(72)+p3mup2*p4mup6*p5mup5*F(72)+p4mup2*p4mu
     -   p6*p5mup5*F(72)+p2mup5*p3mup2*p3mup6*F(73)+p2mup5*p3mup6*p4m
     -   up2*F(73)+p2mup5*p3mup2*p4mup6*F(73)+p2mup5*p4mup2*p4mup6*F(
     -   73)+p1mup5*(-(p6mup2*(p2mup6*F(19)+p1mup6*F(40)-p5mup6*F(42)
     -   +p6mup6*F(50)))+2*((p1mup2+p2mup2)*p3mup6+(p3mup6+p4mup6)*p5
     -   mup2)*F(69)+(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(74))+p3mup2*p3
     -   mup6*p6mup5*F(75)+p3mup6*p4mup2*p6mup5*F(75)+p3mup2*p4mup6*p
     -   6mup5*F(75)+p4mup2*p4mup6*p6mup5*F(75))-8*(p1mup2*p3mup6*p5m
     -   up5*F(70)+p2mup2*p3mup6*p5mup5*F(70)+p3mup6*p5mup2*p5mup5*F(
     -   70)+p4mup6*p5mup2*p5mup5*F(70)+(p3mup6+p4mup6)*(p1mup2+p2mup
     -   2+p5mup2)*p6mup5*F(76))
       Fa1(11) = p2mup6*p6mup2*p6mup5*F(21)+p5mup6*p6mup2*p6mup5*F(43)
     -   -4*(p1mup6*p6mup2*p6mup5*F(52)+(p3mup6+p4mup6)*p6mup2*(p2mup
     -   5*F(68)+2*(p1mup5*F(69)-p5mup5*F(70)-p6mup5*F(76)))+4*(p3mup
     -   5*p3mup6*p5mup2*F(57)+p4mup2*(p3mup6*p4mup5+(p3mup5+p4mup5)*
     -   p4mup6)*F(77)))
       Fa1(12) = (p1mup2*p3mup5*p3mup6+p2mup2*p3mup5*p3mup6+(p3mup6*p4
     -   mup5+(p3mup5+p4mup5)*p4mup6)*p5mup2+(p3mup5+p4mup5)*(p3mup6+
     -   p4mup6)*p6mup2)*F(57)+(p3mup5*p3mup6*p4mup2+p3mup2*(p3mup5+p
     -   4mup5)*(p3mup6+p4mup6))*F(77)
       Fa1(13) = (p1mup2+p2mup2+p5mup2+p6mup2)*F(15)+(p3mup2+p4mup2)*F
     -   (16)
       Fa1(14) = mup3mup6*((p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)
     -   *(-(p1mup5*F(78))-p5mup5*F(80)-p6mup5*F(82)+(p2mup5+p3mup5+p
     -   4mup5)*F(83))+mup2mup5*F(84))+mup3mup5*(-((p1mup2+p2mup2+p3m
     -   up2+p4mup2+p5mup2+p6mup2)*(p6mup6*F(13)-(p2mup6+p3mup6+p4mup
     -   6)*F(14)-p1mup6*F(79)-p5mup6*F(81)))-mup2mup6*F(84))
       Fa1(15) = p4mup3*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)*(4
     -   *(p1mup5*F(78)+p5mup5*F(80)-(p2mup5+p3mup5+p4mup5)*F(83))+p6
     -   mup5*F(85))-4*(mup3mup5*(p3mup2+p4mup2)*F(16)-mup3mup5*(p1mu
     -   p2+p2mup2+p5mup2+p6mup2)*F(86))
       Fa1(16) = mup5mup6*(p2mup2+p3mup2+p4mup2)*p4mup3
       Fa1(17) = mup5mup6*p4mup3*(p1mup2+p5mup2)
       Fa1(18) = mup3mup6*(p2mup5*(p2mup2+p5mup2)*F(94)+(p3mup2+p4mup2
     -   )*(p6mup5*F(90)+4*(p1mup5*F(88)+p5mup5*F(89)+p2mup5*F(95)))-
     -   4*(p2mup2+p5mup2)*(p5mup5*F(99)-p6mup5*F(101)-p1mup5*F(103))
     -   *P(59))+mup3mup5*(-8*p2mup6*(p2mup2+p5mup2)*F(97)+(p3mup2+p4
     -   mup2)*(p1mup6*F(91)+p6mup6*F(93)+4*(p5mup6*F(92)+p2mup6*F(98
     -   )))+4*(mup2mup6*F(96)+(p2mup2+p5mup2)*(p6mup6*F(102)+2*(p5mu
     -   p6*F(100)+p1mup6*F(104)))*P(59)))
       Fa1(19) = -8*mup3mup5*p1mup2*p2mup6*F(97)+4*mup3mup5*(2*p5mup6*
     -   p6mup2*F(100)+p1mup2*(p6mup6*F(102)+2*(p5mup6*F(100)+p1mup6*
     -   F(104))))*P(59)+mup3mup6*(p1mup2+p6mup2)*(p2mup5*F(94)+4*(-(
     -   p5mup5*F(99))+p6mup5*F(101)+p1mup5*F(103))*P(59))
       Fa1(20) = mup3mup5*p6mup2*(-2*p2mup6*F(97)+(p6mup6*F(102)+2*p1m
     -   up6*F(104))*P(59))+mup3mup6*(mup2mup5*F(105)+(p3mup5+p4mup5)
     -   *((p3mup2+p4mup2)*F(106)+(p1mup2+p2mup2+p5mup2)*F(107)*P(59)
     -   ))
       Fa1(21) = mup3mup6*(p3mup5+p4mup5)*p6mup2
       Fa1(22) = mup3mup5*(p3mup6+p4mup6)*((p3mup2+p4mup2)*F(108)+2*(p
     -   1mup2+p2mup2+p5mup2+p6mup2)*F(109)*P(59))
       Fa1(23) = p4mup3*(-(p2mup2*p4mup6*p5mup5*F(110))-p3mup2*p4mup6*
     -   p5mup5*F(110)-p4mup2*p4mup6*p5mup5*F(110)-p4mup6*p5mup2*p5mu
     -   p5*F(110)-p3mup6*(p4mup2+p5mup2)*p5mup5*F(110)+p2mup2*p5mup5
     -   *p5mup6*F(111)+p3mup2*p5mup5*p5mup6*F(111)+p4mup2*p5mup5*p5m
     -   up6*F(111)+p5mup2*p5mup5*p5mup6*F(111)+p2mup6*(p4mup2+p5mup2
     -   )*(-(p5mup5*F(110))+2*p4mup5*F(112))+2*(p2mup5*p4mup2*p4mup6
     -   *F(112)+p3mup5*p4mup2*p4mup6*F(112)+p2mup2*p4mup5*p4mup6*F(1
     -   12)+p3mup2*p4mup5*p4mup6*F(112)+p4mup2*p4mup5*p4mup6*F(112)+
     -   p2mup5*p4mup6*p5mup2*F(112)+p3mup5*p4mup6*p5mup2*F(112)+p4mu
     -   p5*p4mup6*p5mup2*F(112)+p3mup6*(p3mup2*p3mup5+p4mup5*(p4mup2
     -   +p5mup2))*F(112)+((p2mup2+p3mup2)*p4mup5+(p2mup5+p3mup5+p4mu
     -   p5)*(p4mup2+p5mup2))*p5mup6*F(113)))
       Fa1(24) = p4mup3*(4*p3mup5*p3mup6*(p4mup2+p5mup2)*F(112)+p3mup2
     -   *(-2*p3mup6*p5mup5*F(110)+4*(p3mup6*p4mup5*F(112)+p3mup5*p4m
     -   up6*F(112)+p3mup5*p5mup6*F(113)))-p1mup2*p1mup5*p1mup6*F(114
     -   ))
       Return
       End
