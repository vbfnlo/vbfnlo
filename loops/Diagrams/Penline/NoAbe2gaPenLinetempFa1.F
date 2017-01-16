c************** Calling the Fa functions*************************
       subroutine NoAbe2gaPenLineFa1(p1mup2,p1mup4,p1mup5,p2mup2,p2mup
     &   4,p2mup5,p3mup2,p3mup4,p3mup5,p4mup2,p4mup4,p4mup5,p5mup2,p5
     &   mup4,p5mup5,mup2mup4,mup2mup5,mup4mup5,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup4, p1mup5, p2mup2, p2mup4, p2mup5, 
     -          p3mup2, p3mup4, p3mup5, p4mup2, p4mup4, p4mup5, 
     -          p5mup2, p5mup4, p5mup5
       Complex*16   mup2mup4, mup2mup5, mup4mup5
        common/NoAbe2gaPenLineFhlFunctions/F
       COMMON/NoAbe2gaPenLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s
     &   23,s34,s45,s15
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45
      Complex*16 Fa(84),F(426)
      Real*8 P(767) 
        COMMON/NoAbe2gaPenLinePFunctions/P
       Complex*16 Fa1(1:28)
       COMMON/NoAbe2gaPenLineFaFunctions/Fa
       Fa1(1) = p3mup5*p5mup4*F(1)+p1mup4*p3mup5*F(2)+p3mup4*p3mup5*F(
     &   3)+p2mup5*p3mup4*F(4)+p1mup4*p2mup5*F(5)+p2mup4*p2mup5*F(6)+
     &   p2mup4*p3mup5*F(7)+p1mup4*p1mup5*F(8)+p1mup5*p3mup4*F(9)+p1m
     &   up5*p5mup4*F(10)+p1mup5*p2mup4*F(11)+p3mup4*p5mup5*F(12)+p2m
     &   up4*p5mup5*F(13)+p2mup5*p5mup4*F(14)+p1mup4*p5mup5*F(15)+p5m
     &   up4*p5mup5*F(16)+mup4mup5*F(17)
       Fa1(2) = p2mup2+p3mup2
       Fa1(3) = p2mup5+p3mup5
       Fa1(4) = p1mup5*F(21)+p5mup5*F(23)+F(20)*Fa1(3)
       Fa1(5) = p5mup2*(p3mup5*F(18)+p1mup5*F(19)+p5mup5*F(22)+p2mup5*
     &   F(24))+mup2mup5*F(25)+Fa1(2)*Fa1(4)
       Fa1(6) = p2mup4+p3mup4
       Fa1(7) = p3mup4*p5mup2*F(26)+p1mup4*p5mup2*F(27)+p1mup4*p2mup2*
     &   F(29)+p1mup4*p3mup2*F(29)+p2mup2*p5mup4*F(30)+p3mup2*p5mup4*
     &   F(30)+p5mup2*p5mup4*F(31)+p2mup4*p5mup2*F(32)+mup2mup4*F(33)
     &   +F(28)*Fa1(2)*Fa1(6)
       Fa1(8) = (p2mup4+p3mup4)*(p2mup5+p3mup5)*p5mup2*F(36)+p1mup5*p5
     &   mup2*p5mup4*F(37)+p1mup4*p1mup5*p5mup2*F(38)+p1mup5*p2mup4*p
     &   5mup2*F(39)+p1mup5*p3mup4*p5mup2*F(39)+p1mup4*p5mup2*p5mup5*
     &   F(40)+p2mup4*p5mup2*p5mup5*F(41)+p3mup4*p5mup2*p5mup5*F(41)+
     &   p5mup2*p5mup4*p5mup5*F(42)+mup2mup4*(p2mup5+p3mup5)*F(43)+mu
     &   p2mup4*p1mup5*F(44)+mup2mup4*p5mup5*F(45)+mup2mup5*p2mup4*F(
     &   46)+mup2mup5*p3mup4*F(46)+mup2mup5*p1mup4*F(47)+mup2mup5*p5m
     &   up4*F(48)+mup4mup5*p5mup2*F(49)+p5mup2*(p5mup4*F(34)+p1mup4*
     &   F(35))*Fa1(3)
       Fa1(9) = (p2mup2+p3mup2)*(p2mup4+p3mup4)*((p2mup5+p3mup5)*F(52)
     &   +p1mup5*F(53))+p1mup5*(p2mup2+p3mup2)*p5mup4*F(54)+p1mup4*p1
     &   mup5*p2mup2*F(55)+p1mup4*p1mup5*p3mup2*F(55)+p1mup4*p2mup2*p
     &   5mup5*F(56)+p1mup4*p3mup2*p5mup5*F(56)+p2mup2*p2mup4*p5mup5*
     &   F(57)+p2mup4*p3mup2*p5mup5*F(57)+p2mup2*p3mup4*p5mup5*F(57)+
     &   p3mup2*p3mup4*p5mup5*F(57)+p2mup2*p5mup4*p5mup5*F(58)+p3mup2
     &   *p5mup4*p5mup5*F(58)+mup2mup4*p3mup5*F(59)+mup2mup4*p1mup5*F
     &   (60)+mup2mup4*p5mup5*F(61)+mup2mup4*p2mup5*F(62)+mup2mup5*p3
     &   mup4*F(63)+mup2mup5*p1mup4*F(64)+mup2mup5*p5mup4*F(65)+mup2m
     &   up5*p2mup4*F(66)+Fa1(2)*(mup4mup5*F(67)+(p5mup4*F(50)+p1mup4
     &   *F(51))*Fa1(3))
       Fa1(10) = p1mup5*F(69)+p5mup5*F(70)+F(68)*Fa1(3)
       Fa1(11) = p3mup5*F(71)+p1mup5*F(72)+p5mup5*F(73)+p2mup5*F(74)
       Fa1(12) = p1mup4*F(79)+p5mup4*F(80)+F(78)*Fa1(6)
       Fa1(13) = p3mup4*F(81)+p1mup4*F(82)+p5mup4*F(83)+p2mup4*F(84)
       Fa1(14) = p1mup5*F(90)+p5mup5*F(92)+F(87)*Fa1(3)
       Fa1(15) = p1mup5*p5mup4*F(88)+p1mup4*p1mup5*F(89)+p1mup4*p5mup5
     &   *F(91)+p5mup4*p5mup5*F(93)+mup4mup5*F(94)+(p5mup4*F(85)+p1mu
     &   p4*F(86))*Fa1(3)+Fa1(6)*Fa1(14)
       Fa1(16) = p1mup5+p2mup5
       Fa1(17) = p1mup4+p2mup4
       Fa1(18) = p3mup5*p5mup4*F(98)+(p1mup5+p2mup5)*p3mup4*F(99)+p3mu
     &   p4*p3mup5*F(101)+p3mup4*p5mup5*F(102)+p1mup4*p5mup5*F(103)+p
     &   2mup4*p5mup5*F(103)+p5mup4*p5mup5*F(104)+mup4mup5*F(105)+p3m
     &   up5*F(100)*Fa1(17)+Fa1(16)*(p5mup4*F(96)+F(97)*Fa1(17))
       Fa1(19) = p3mup5*F(106)+p5mup5*F(108)+F(107)*Fa1(16)
       Fa1(20) = p3mup4*F(109)+p5mup4*F(111)+F(110)*Fa1(17)
       Fa1(21) = p5mup4*F(115)+F(116)*Fa1(17)
       Fa1(22) = p3mup5*F(114)+p5mup5*F(119)+F(117)*Fa1(16)
       Fa1(23) = p3mup5*p5mup4*F(112)+p5mup4*p5mup5*F(120)+mup4mup5*F(
     &   121)+(p3mup5*F(113)+p5mup5*F(118))*Fa1(17)+Fa1(16)*Fa1(21)+p
     &   3mup4*Fa1(22)
       Fa1(24) = p3mup5*F(122)+p5mup5*F(124)+F(123)*Fa1(16)
       Fa1(25) = p3mup5*F(126)+p5mup5*F(127)+F(125)*Fa1(16)
       Fa1(26) = p3mup4*F(130)+p5mup4*F(132)+F(131)*Fa1(17)
       Fa1(27) = p3mup4*F(134)+p5mup4*F(135)+F(133)*Fa1(17)
       Fa1(28) = p5mup4*F(139)+F(140)*Fa1(17)
       Return
       End
