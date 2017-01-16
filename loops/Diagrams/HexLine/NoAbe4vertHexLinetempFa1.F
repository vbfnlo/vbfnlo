c************** Calling the Fa functions*************************
       subroutine NoAbe4vertHexLineFa1(p1mup2,p1mup3,p1mup5,p1mup6,p2m
     -   up2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,
     -   p4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mu
     -   p3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup
     -   6,mup5mup6,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe4vertHexLineFhlFunctions/F
       COMMON/NoAbe4vertHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6s
     -   q,s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(199),F(140)
      Real*8 P(39) 
        COMMON/NoAbe4vertHexLinePFunctions/P
       Complex*16 Fa1(1:19)
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       Fa1(1) = 4*p4mup6*(-(mup3mup5*F(1))+2*p2mup5*(p3mup3*F(2)+p1mup
     -   3*F(3)+p2mup3*F(4)+p4mup3*F(134)))
       Fa1(2) = mup5mup6*(p2mup3*F(5)-p4mup3*F(6)+p3mup3*F(7)+p1mup3*F
     -   (8))
       Fa1(3) = mup3mup6*p4mup5*F(1)+mup5mup6*(p2mup3*F(9)+p3mup3*F(10
     -   )+p1mup3*F(11)-p4mup3*F(12))-2*p2mup6*p4mup5*(p3mup3*F(2)+p1
     -   mup3*F(3)+p2mup3*F(4)+p4mup3*F(134))
       Fa1(4) = p4mup6*(p3mup2*p3mup3*F(18)+p1mup3*p3mup2*F(19)+p2mup3
     -   *p3mup2*F(20)+p1mup3*p4mup2*F(21)+p1mup2*(p2mup3*F(14)+p3mup
     -   3*F(15)+p1mup3*F(17)+p4mup3*F(21))+p2mup3*p4mup2*F(22)+p2mup
     -   2*(p2mup3*F(13)+p1mup3*F(14)+p3mup3*F(16)+p4mup3*F(22))+p3mu
     -   p3*p4mup2*F(23)+p3mup2*p4mup3*F(137)+mup2mup3*F(139)+p4mup2*
     -   p4mup3*F(140))
       Fa1(5) = -(mup2mup6*(p2mup3*F(9)+p3mup3*F(10)+p1mup3*F(11)-p4mu
     -   p3*F(12)))+mup3mup6*p3mup2*F(24)-mup3mup6*p4mup2*F(25)-mup3m
     -   up6*p1mup2*F(39)-mup3mup6*p2mup2*F(40)+2*(p1mup3*p2mup6*p4mu
     -   p2*F(14)+p1mup6*p2mup2*p4mup3*F(14)+p1mup6*p3mup3*p4mup2*F(1
     -   5)+p1mup2*p3mup6*p4mup3*F(15)+p2mup6*p3mup3*p4mup2*F(16)+p2m
     -   up2*p3mup6*p4mup3*F(16)+p1mup3*p1mup6*p4mup2*F(17)+p1mup2*p1
     -   mup6*p4mup3*F(17)+p3mup3*p3mup6*p4mup2*F(18)+p3mup2*p3mup6*p
     -   4mup3*F(18)+p1mup3*p3mup6*p4mup2*F(19)+p1mup6*p3mup2*p4mup3*
     -   F(19)+p1mup6*p4mup2*p4mup3*F(21)+p2mup6*p4mup2*p4mup3*F(22)+
     -   p3mup6*p4mup2*p4mup3*F(23)+p1mup2*p1mup3*p1mup6*F(26)+p1mup6
     -   *p2mup2*p3mup3*F(27)+p1mup3*p1mup6*p2mup2*F(28)+p1mup2*p1mup
     -   6*p3mup3*F(29)+p1mup2*p2mup6*p3mup3*F(31)+p1mup2*p1mup3*p2mu
     -   p6*F(33)+p1mup3*p2mup2*p2mup6*F(34)+p2mup2*p2mup6*p3mup3*F(3
     -   5)+p1mup6*p3mup2*p3mup3*F(36)+p1mup2*p3mup3*p3mup6*F(36)+p2m
     -   up2*p3mup3*p3mup6*F(37)+p1mup3*p2mup2*p3mup6*F(38)+p1mup3*p1
     -   mup6*p3mup2*F(41)+p1mup2*p1mup3*p3mup6*F(41)-p2mup6*p3mup2*p
     -   4mup3*F(43)+mup2mup3*p1mup6*F(44)+mup2mup3*p2mup6*F(45)+p1mu
     -   p3*p3mup2*p3mup6*F(46)+p3mup2*p3mup3*p3mup6*F(48)+p2mup6*p3m
     -   up2*p3mup3*F(49)+p1mup3*p2mup6*p3mup2*F(50)+p2mup3*(p1mup6*(
     -   p4mup2*F(14)+p1mup2*F(28)+p2mup2*F(30)+p3mup2*F(38))+p3mup6*
     -   (p4mup2*F(20)+p1mup2*F(38)+p2mup2*F(42)+p3mup2*F(47))+p2mup6
     -   *(p4mup2*F(13)+p2mup2*F(32)+p1mup2*F(34)+p3mup2*F(51)))+p2mu
     -   p2*p2mup6*p4mup3*F(52)+p1mup2*p2mup6*p4mup3*F(53)+mup2mup3*p
     -   3mup6*F(138))
       Fa1(6) = p2mup3*F(9)+p3mup3*F(10)+p1mup3*F(11)-p4mup3*F(12)
       Fa1(7) = -(mup2mup5*(p2mup3*F(9)+p3mup3*F(10)+p1mup3*F(11)-p4mu
     -   p3*F(12)))+mup3mup5*p3mup2*F(24)-mup3mup5*p4mup2*F(25)-mup3m
     -   up5*p1mup2*F(39)-mup3mup5*p2mup2*F(40)+2*(p1mup3*p2mup5*p4mu
     -   p2*F(14)+p1mup5*p2mup2*p4mup3*F(14)+p1mup5*p3mup3*p4mup2*F(1
     -   5)+p1mup2*p3mup5*p4mup3*F(15)+p2mup5*p3mup3*p4mup2*F(16)+p2m
     -   up2*p3mup5*p4mup3*F(16)+p1mup3*p1mup5*p4mup2*F(17)+p1mup2*p1
     -   mup5*p4mup3*F(17)+p3mup3*p3mup5*p4mup2*F(18)+p3mup2*p3mup5*p
     -   4mup3*F(18)+p1mup3*p3mup5*p4mup2*F(19)+p1mup5*p3mup2*p4mup3*
     -   F(19)+p1mup5*p4mup2*p4mup3*F(21)+p2mup5*p4mup2*p4mup3*F(22)+
     -   p3mup5*p4mup2*p4mup3*F(23)+p1mup2*p1mup3*p1mup5*F(26)+p1mup5
     -   *p2mup2*p3mup3*F(27)+p1mup3*p1mup5*p2mup2*F(28)+p1mup2*p1mup
     -   5*p3mup3*F(29)+p1mup2*p2mup5*p3mup3*F(31)+p1mup2*p1mup3*p2mu
     -   p5*F(33)+p1mup3*p2mup2*p2mup5*F(34)+p2mup2*p2mup5*p3mup3*F(3
     -   5)+p1mup5*p3mup2*p3mup3*F(36)+p1mup2*p3mup3*p3mup5*F(36)+p2m
     -   up2*p3mup3*p3mup5*F(37)+p1mup3*p2mup2*p3mup5*F(38)+p1mup3*p1
     -   mup5*p3mup2*F(41)+p1mup2*p1mup3*p3mup5*F(41)-p2mup5*p3mup2*p
     -   4mup3*F(43)+mup2mup3*p1mup5*F(44)+mup2mup3*p2mup5*F(45)+p1mu
     -   p3*p3mup2*p3mup5*F(46)+p3mup2*p3mup3*p3mup5*F(48)+p2mup5*p3m
     -   up2*p3mup3*F(49)+p1mup3*p2mup5*p3mup2*F(50)+p2mup3*(p1mup5*(
     -   p4mup2*F(14)+p1mup2*F(28)+p2mup2*F(30)+p3mup2*F(38))+p3mup5*
     -   (p4mup2*F(20)+p1mup2*F(38)+p2mup2*F(42)+p3mup2*F(47))+p2mup5
     -   *(p4mup2*F(13)+p2mup2*F(32)+p1mup2*F(34)+p3mup2*F(51)))+p2mu
     -   p2*p2mup5*p4mup3*F(52)+p1mup2*p2mup5*p4mup3*F(53)+mup2mup3*p
     -   3mup5*F(138))
       Fa1(8) = p4mup5*(p3mup2*p3mup3*F(18)+p1mup3*p3mup2*F(19)+p2mup3
     -   *p3mup2*F(20)+p1mup3*p4mup2*F(21)+p1mup2*(p2mup3*F(14)+p3mup
     -   3*F(15)+p1mup3*F(17)+p4mup3*F(21))+p2mup3*p4mup2*F(22)+p2mup
     -   2*(p2mup3*F(13)+p1mup3*F(14)+p3mup3*F(16)+p4mup3*F(22))+p3mu
     -   p3*p4mup2*F(23)+p3mup2*p4mup3*F(137)+mup2mup3*F(139)+p4mup2*
     -   p4mup3*F(140))
       Fa1(9) = 4*p4mup6*(mup2mup5*F(1)+2*(-(p2mup5*p3mup2*F(132))+p3m
     -   up5*(p1mup2*F(3)+p2mup2*F(4)+p3mup2*F(133)+p4mup2*F(134))))
       Fa1(10) = -4*mup5mup6*(p2mup2*F(54)+p3mup2*F(55)+p1mup2*F(56))+
     -   mup5mup6*p4mup2*F(57)
       Fa1(11) = 4*(mup2mup6*p4mup5*F(1)+mup5mup6*p3mup2*F(24)+mup5mup
     -   6*p4mup2*F(58))+mup5mup6*p2mup2*F(59)+mup5mup6*p1mup2*F(60)+
     -   8*(p1mup2*p3mup6*p4mup5*F(3)+p2mup2*p3mup6*p4mup5*F(4)-p2mup
     -   6*p3mup2*p4mup5*F(132)+p3mup2*p3mup6*p4mup5*F(133)+p3mup6*p4
     -   mup2*p4mup5*F(134))
       Fa1(12) = -(p4mup6*(-(mup3mup5*p3mup2*F(132))+mup2mup5*(p3mup3*
     -   F(2)+p1mup3*F(3)+p2mup3*F(4)+p4mup3*F(134))))
       Fa1(13) = 4*mup5mup6*(-(mup2mup3*F(61))-2*(-(p2mup3*p3mup2*F(42
     -   ))+p3mup3*p4mup2*F(43)-p3mup2*p3mup3*F(47)-p2mup3*p4mup2*F(5
     -   2)-p2mup2*(p2mup3*F(32)+p1mup3*F(34)+p3mup3*F(51)+p4mup3*F(5
     -   2))-p1mup3*p4mup2*F(53)-p1mup2*(p4mup3*F(14)+p1mup3*F(28)+p2
     -   mup3*F(30)+p3mup3*F(62))-p1mup3*p3mup2*F(63)-p3mup2*p4mup3*F
     -   (64)-p4mup2*p4mup3*F(136)))
       Fa1(14) = mup5mup6*(p3mup3*p4mup2*F(16)+p4mup2*p4mup3*F(22)+p1m
     -   up2*p3mup3*F(31)+p2mup2*p3mup3*F(35)-p3mup2*p4mup3*F(43)+mup
     -   2mup3*F(45)+p3mup2*p3mup3*F(49)+p1mup3*(p4mup2*F(14)+p1mup2*
     -   F(33)+p2mup2*F(34)+p3mup2*F(50))+p2mup3*(p4mup2*F(13)+p2mup2
     -   *F(32)+p1mup2*F(34)+p3mup2*F(51))+p2mup2*p4mup3*F(52)+p1mup2
     -   *p4mup3*F(53))+mup3mup6*p3mup2*p4mup5*F(132)-mup2mup6*p4mup5
     -   *(p3mup3*F(2)+p1mup3*F(3)+p2mup3*F(4)+p4mup3*F(134))
       Fa1(15) = mup3mup5*p4mup6*(p1mup2*F(3)+p2mup2*F(4)+p3mup2*F(133
     -   )+p4mup2*F(134))
       Fa1(16) = 4*mup5mup6*(-(mup2mup3*F(65))+2*(p1mup2*p4mup3*F(19)+
     -   p3mup3*(p4mup2*F(18)+p1mup2*F(36)+p3mup2*F(48)+p2mup2*F(49))
     -   +p1mup3*(p4mup2*F(15)+p1mup2*F(29)+p3mup2*F(36)+p2mup2*F(66)
     -   )+p2mup3*(p4mup2*F(16)+p2mup2*F(35)+p3mup2*F(37)+p1mup2*F(67
     -   ))+p4mup3*(-(p2mup2*F(43))+p3mup2*F(135)+p4mup2*F(137))))
       Fa1(17) = -(mup3mup6*p4mup5*(p1mup2*F(3)+p2mup2*F(4)+p3mup2*F(1
     -   33)+p4mup2*F(134)))+mup5mup6*(p3mup3*p4mup2*F(18)+p3mup2*p4m
     -   up3*F(18)+p1mup3*p4mup2*F(19)+p2mup3*p4mup2*F(20)+p4mup2*p4m
     -   up3*F(23)+p1mup2*(p4mup3*F(15)+p3mup3*F(36)+p2mup3*F(38)+p1m
     -   up3*F(41))+p2mup2*(p4mup3*F(16)+p3mup3*F(37)+p1mup3*F(38)+p2
     -   mup3*F(42))+p1mup3*p3mup2*F(46)+p2mup3*p3mup2*F(47)+p3mup2*p
     -   3mup3*F(48)+mup2mup3*F(138))
       Fa1(18) = 2*(p3mup2*F(24)+p4mup2*F(58))+p1mup2*F(68)+p2mup2*F(6
     -   9)
       Fa1(19) = -(p3mup2*F(24))+p4mup2*F(25)+p1mup2*F(39)+p2mup2*F(40
     -   )
       Return
       End
