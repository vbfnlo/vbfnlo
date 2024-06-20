c************** Calling the Fa functions*************************
       subroutine HexBox2Fa3(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa3(39:57)
       COMMON/HexBox2FaFunctions/Fa
       Fa3(39) = p1mup3*p2mup6*F(99)+3*mup3mup6*F(1631)
       Fa3(40) = (p5mup3+p6mup3)*F(108)-2*p2mup3*F(117)+p1mup3*F(124)
       Fa3(41) = p2mup3*(p5mup6+p6mup6)*F(99)-p2mup6*(p5mup3+p6mup3)*F
     &   (111)
       Fa3(42) = p1mup6*p2mup3*F(99)-p1mup3*p2mup6*F(111)
       Fa3(43) = -4*((p5mup3+p6mup3)*s12*F(96)-p2mup3*F(117))+p1mup3*F
     &   (125)
       Fa3(44) = p5mup3+p6mup3
       Fa3(45) = 8*(-(p2mup6*(p5mup3+p6mup3))+p2mup3*(p5mup6+p6mup6))*
     &   F(99)+4*mup3mup6*F(128)
       Fa3(46) = -(p2mup6*F(117))+6*(p1mup6+p5mup6+p6mup6)*F(1631)
       Fa3(47) = p2mup6*(p5mup3+p6mup3)
       Fa3(48) = p1mup3*s56*F(105)-(p5mup3+p6mup3)*F(130)
       Fa3(49) = (p5mup3+p6mup3)*(p5mup6+p6mup6)
       Fa3(50) = p1mup6*(p5mup3+p6mup3)*F(99)-p1mup3*(p5mup6+p6mup6)*F
     &   (111)
       Fa3(51) = p6mup3*F(123)+4*(p5mup3*F(97)+p1mup3*F(132))
       Fa3(52) = p1mup6*(p5mup3+p6mup3)-p1mup3*(p5mup6+p6mup6)
       Fa3(53) = p2mup6*s56*F(105)+(p5mup6+p6mup6)*F(135)
       Fa3(54) = 2*p1mup6*F(119)-p2mup6*F(130)+(p5mup6+p6mup6)*F(136)
       Fa3(55) = -8*(p5mup6+p6mup6)*F(108)-16*p1mup6*F(128)+p2mup6*F(1
     &   37)
       Fa3(56) = p1mup3*F(128)+6*(p5mup3+p6mup3)*F(1631)
       Fa3(57) = (p1mup3+p5mup3+p6mup3)*(p5mup6+p6mup6)
       Return
       End
