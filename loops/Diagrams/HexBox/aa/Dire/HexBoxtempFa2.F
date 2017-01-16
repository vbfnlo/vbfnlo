c************** Calling the Fa functions*************************
       subroutine HexBoxFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2)
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
       Complex*16 Fa2(13:24)
       COMMON/HexBoxFaFunctions/Fa
       Fa2(13) = mup3mup4*F(62)+8*(p1mup3*(p4mup4+p5mup4)*F(63)+p4mup4
     &   *p6mup3*s12*F(64)+p5mup4*p6mup3*s12*F(64)-p2mup3*p4mup4*F(66
     &   )-p2mup3*p5mup4*F(66)+p2mup4*p6mup3*F(65)*P(50)+p2mup3*p6mup
     &   4*F(67)*P(50)+p1mup4*p6mup3*F(65)*P(116)+p1mup3*p6mup4*F(67)
     &   *P(116))
       Fa2(14) = (p4mup4+p5mup4)*F(69)+2*p6mup4*F(70)*P(50)
       Fa2(15) = 2*p1mup3*(p4mup4+p5mup4)+mup3mup4*P(50)
       Fa2(16) = (p4mup4+p5mup4)*p6mup3
       Fa2(17) = mup3mup4*F(69)+4*(p1mup4*p6mup3-p1mup3*p6mup4)*F(70)
       Fa2(18) = -2*p2mup4*F(66)+p1mup4*F(75)+p6mup4*F(76)
       Fa2(19) = -2*p2mup3*F(66)+p1mup3*F(75)+p6mup3*F(76)
       Fa2(20) = 4*(p2mup4*p6mup3-p2mup3*p6mup4)*F(65)+mup3mup4*F(77)
       Fa2(21) = 4*(p1mup4*p2mup3-p1mup3*p2mup4)*F(65)+mup3mup4*F(78)
       Fa2(22) = (p4mup4+p5mup4)*F(78)+2*F(65)*(p2mup4*P(50)+p1mup4*P(
     &   116))
       Fa2(23) = p2mup3*P(50)+p1mup3*P(116)
       Fa2(24) = p4mup4*F(80)+p5mup4*F(81)-p6mup4*F(65)*P(116)
       Return
       End
