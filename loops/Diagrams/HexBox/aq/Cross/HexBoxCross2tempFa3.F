c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa3(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa3
     &   )
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
       Complex*16 Fa3(39:57)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa3(39) = (p1mup3+p5mup3+p6mup3)*(p1mup6+p6mup6)
       Fa3(40) = 4*p6mup6*s56*F(93)+p5mup6*F(95)
       Fa3(41) = p5mup6*(p5mup3+p6mup3)
       Fa3(42) = p5mup6*(p1mup3*F(108)+(p5mup3+p6mup3)*F(109)+p2mup3*F
     &   (110))
       Fa3(43) = -(p6mup3*F(112))-p2mup3*s56*F(113)+p5mup3*s56*F(1837)
       Fa3(44) = -(p6mup3*F(114))+p2mup3*s56*F(115)+p1mup3*s56*F(1838)
       Fa3(45) = 2*p2mup3*p5mup6*F(118)+mup3mup6*s56*F(119)
       Fa3(46) = p6mup6*s56*F(119)-p5mup6*F(120)
       Fa3(47) = -4*(p5mup3*s56*F(115)-p1mup3*s56*F(125))+p6mup3*F(126
     &   )
       Fa3(48) = p1mup3*F(108)+(p5mup3+p6mup3)*F(109)+p2mup3*F(110)
       Fa3(49) = p5mup3*s56*F(119)+p6mup3*F(130)
       Fa3(50) = p1mup3*s56*F(119)-6*p6mup3*F(1839)
       Fa3(51) = -2*p1mup3*s56*F(119)+12*p6mup3*F(1839)
       Fa3(52) = p5mup3+p6mup3
       Fa3(53) = p5mup6*(p1mup3+p5mup3+p6mup3)
       Fa3(54) = p1mup6*F(132)+p4mup6*F(133)+p5mup6*F(134)+p6mup6*F(13
     &   5)
       Fa3(55) = 8*p1mup6*F(136)+p4mup6*F(137)+p5mup6*F(138)+p6mup6*F(
     &   139)
       Fa3(56) = p1mup6*F(158)+p5mup6*F(159)+p6mup6*F(160)+p4mup6*F(16
     &   1)
       Fa3(57) = p1mup6*F(162)+p5mup6*F(163)+p6mup6*F(164)+p4mup6*F(16
     &   5)
       Return
       End
