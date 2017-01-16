c************** Calling the Fa functions*************************
       subroutine HexBoxFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa1)
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
       Complex*16 Fa1(1:12)
       COMMON/HexBoxFaFunctions/Fa
       Fa1(1) = p1mup4+p4mup4+p5mup4+p6mup4
       Fa1(2) = (p4mup4+p5mup4)*F(22)+p1mup4*s45*F(23)-p6mup4*F(24)
       Fa1(3) = p1mup4*F(28)+p6mup4*F(29)+(p4mup4+p5mup4)*F(30)
       Fa1(4) = p4mup4*F(37)+p5mup4*F(38)-p6mup4*F(39)
       Fa1(5) = p1mup4*F(40)+4*((p4mup4+p5mup4)*F(27)+p6mup4*F(41))
       Fa1(6) = (p4mup4+p5mup4)*F(42)+p1mup4*F(43)
       Fa1(7) = (p4mup4+p5mup4)*F(27)+p1mup4*F(45)+p6mup4*F(46)
       Fa1(8) = p1mup4*s45*F(23)+(p4mup4+p5mup4)*F(27)+p6mup4*F(47)
       Fa1(9) = (p4mup4+p5mup4)*F(27)+p1mup4*F(45)-p6mup4*F(50)
       Fa1(10) = p1mup4*s45*F(23)+(p4mup4+p5mup4)*F(27)-p6mup4*F(49)
       Fa1(11) = p4mup4+p5mup4
       Fa1(12) = p1mup4+p6mup4
       Return
       End
