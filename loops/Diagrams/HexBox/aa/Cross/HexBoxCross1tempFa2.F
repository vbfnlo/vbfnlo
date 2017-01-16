c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa2
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
       Complex*16 Fa2(14:26)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa2(14) = (p4mup4+p5mup4)*p6mup3
       Fa2(15) = 2*p2mup4*F(73)+p1mup4*F(76)+p6mup4*F(77)
       Fa2(16) = 2*p2mup3*F(73)+p1mup3*F(76)+p6mup3*F(77)
       Fa2(17) = p1mup4*F(69)+p6mup4*F(70)+p2mup4*F(73)
       Fa2(18) = (p4mup4+p5mup4)*F(79)-p6mup4*F(80)*P(45)
       Fa2(19) = (p4mup4+p5mup4)*F(82)+2*p6mup4*F(72)*P(115)
       Fa2(20) = 4*(-(p2mup4*p6mup3)+p2mup3*p6mup4)*F(72)+mup3mup4*F(8
     &   2)
       Fa2(21) = (p4mup4+p5mup4)*F(83)-2*F(72)*(p2mup4*P(45)+p1mup4*P(
     &   115))
       Fa2(22) = p2mup3*P(45)+p1mup3*P(115)
       Fa2(23) = 4*(-(p1mup4*p2mup3)+p1mup3*p2mup4)*F(72)+mup3mup4*F(8
     &   3)
       Fa2(24) = 8*(p2mup4*p6mup3-p2mup3*p6mup4)*F(72)-4*mup3mup4*F(85
     &   )
       Fa2(25) = 8*(-(p1mup4*p2mup3)+p1mup3*p2mup4)*F(72)+4*mup3mup4*F
     &   (86)
       Fa2(26) = -24*p4mup3*F(1452)+4*p6mup3*F(71)*P(115)
       Return
       End
