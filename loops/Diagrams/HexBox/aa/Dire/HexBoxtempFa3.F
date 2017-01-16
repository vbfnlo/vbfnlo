c************** Calling the Fa functions*************************
       subroutine HexBoxFa3(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa3(25:36)
       COMMON/HexBoxFaFunctions/Fa
       Fa3(25) = -(p2mup4*F(66))+p1mup4*F(81)+p6mup4*F(82)
       Fa3(26) = -((p4mup4+p5mup4)*F(82))+F(67)*(p2mup4*P(50)+p1mup4*P
     &   (116))
       Fa3(27) = -6*p4mup3*F(1305)+p6mup3*F(67)*P(116)
       Fa3(28) = -(p2mup3*F(66))+p1mup3*F(81)+p6mup3*F(82)
       Fa3(29) = 8*p2mup3*(p4mup4+p5mup4)*F(65)+4*mup3mup4*F(67)*P(116
     &   )
       Fa3(30) = 8*(-(p1mup4*p2mup3)+p1mup3*p2mup4)*F(65)-4*mup3mup4*F
     &   (83)
       Fa3(31) = (p4mup4+p5mup4)*F(66)+p6mup4*F(65)*P(50)
       Fa3(32) = 2*(p4mup4+p5mup4)*(p1mup3+p6mup3)+mup3mup4*P(50)
       Fa3(33) = 2*(p1mup4*p6mup3-p1mup3*p6mup4)*F(65)+mup3mup4*F(66)
       Fa3(34) = -2*p1mup3*(p4mup4+p5mup4)*F(65)+mup3mup4*F(67)*P(50)
       Fa3(35) = p1mup3+p6mup3
       Fa3(36) = 2*(p2mup4*p6mup3-p2mup3*p6mup4)*F(65)+mup3mup4*F(80)
       Return
       End
