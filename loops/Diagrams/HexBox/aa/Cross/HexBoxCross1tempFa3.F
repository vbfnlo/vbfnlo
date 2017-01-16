c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa3(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa3
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa3(27:39)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa3(27) = p1mup3*F(69)+p6mup3*F(70)+p2mup3*F(73)
       Fa3(28) = 2*p1mup3*(p4mup4+p5mup4)+mup3mup4*P(45)
       Fa3(29) = p4mup4*F(85)+p5mup4*F(88)+p6mup4*F(72)*P(115)
       Fa3(30) = (p4mup4+p5mup4)*F(70)+F(71)*(p2mup4*P(45)+p1mup4*P(11
     &   5))
       Fa3(31) = mup3mup4*F(79)+2*(-(p1mup4*p6mup3)+p1mup3*p6mup4)*F(8
     &   0)
       Fa3(32) = (p4mup4+p5mup4)*F(73)+p6mup4*F(72)*P(45)
       Fa3(33) = 8*(-(p1mup4*p6mup3)+p1mup3*p6mup4)*F(72)-4*mup3mup4*F
     &   (73)
       Fa3(34) = 8*p2mup3*(p4mup4+p5mup4)*F(72)+4*mup3mup4*F(71)*P(115
     &   )
       Fa3(35) = p1mup3+p6mup3
       Fa3(36) = 2*(p4mup4+p5mup4)*(p1mup3+p6mup3)+mup3mup4*P(45)
       Fa3(37) = 4*(-(p1mup4*p6mup3)+p1mup3*p6mup4)*F(72)-2*mup3mup4*F
     &   (73)
       Fa3(38) = -2*p1mup3*(p4mup4+p5mup4)*F(72)+mup3mup4*F(71)*P(45)
       Fa3(39) = p4mup4*F(94)+8*p1mup4*F(95)+p5mup4*F(96)+p6mup4*F(97)
       Return
       End
