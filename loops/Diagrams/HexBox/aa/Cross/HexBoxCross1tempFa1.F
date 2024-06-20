c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa1
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
       Complex*16 Fa1(1:13)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa1(1) = p1mup4+p4mup4+p5mup4+p6mup4
       Fa1(2) = p1mup4*F(23)+8*((p4mup4+p5mup4)*F(22)+p6mup4*F(24))
       Fa1(3) = (p4mup4+p5mup4)*F(25)+p6mup4*F(26)
       Fa1(4) = p1mup4*F(28)-p6mup4*F(29)+4*(p4mup4+p5mup4)*F(30)
       Fa1(5) = (p4mup4+p5mup4)*F(30)+p1mup4*s45*F(32)+p6mup4*F(43)
       Fa1(6) = (p4mup4+p5mup4)*F(30)+p1mup4*F(44)+p6mup4*F(45)
       Fa1(7) = p4mup4*F(48)+p5mup4*F(49)+p6mup4*F(50)
       Fa1(8) = (p4mup4+p5mup4)*F(51)-p1mup4*F(52)
       Fa1(9) = p4mup4+p5mup4
       Fa1(10) = (p4mup4+p5mup4)*F(30)+p1mup4*F(44)-p6mup4*F(58)
       Fa1(11) = (p4mup4+p5mup4)*F(30)+p1mup4*s45*F(32)-p6mup4*F(53)
       Fa1(12) = p1mup4+p6mup4
       Fa1(13) = 4*mup3mup4*F(68)-8*((p4mup4+p5mup4)*(p1mup3*F(69)+p6m
     &   up3*F(70)+p2mup3*F(73))-(p2mup3*p6mup4*F(71)+p2mup4*p6mup3*F
     &   (72))*P(45)-(p1mup3*p6mup4*F(71)+p1mup4*p6mup3*F(72))*P(115)
     &   )
       Return
       End
