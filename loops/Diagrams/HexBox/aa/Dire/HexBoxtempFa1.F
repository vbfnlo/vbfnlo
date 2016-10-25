c************** Calling the Fa functions*************************
       subroutine HexBoxFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
       Real*8   p1sq, p1p2, p1p3, p1p4, p1p5,p1p6 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq, p5p6 
       Real*8   p6sq 
       Real*8   s12, s13, s14, s15, s16 
       Real*8   s23, s24, s25,s26 
       Real*8   s34, s35,s36 
       Real*8   s45,s46 
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
