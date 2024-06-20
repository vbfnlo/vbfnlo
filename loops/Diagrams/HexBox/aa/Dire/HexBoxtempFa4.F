c************** Calling the Fa functions*************************
       subroutine HexBoxFa4(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa4)
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
       Complex*16 Fa4(37:48)
       COMMON/HexBoxFaFunctions/Fa
       Fa4(37) = (p4mup4+p5mup4)*F(77)-2*p6mup4*F(65)*P(116)
       Fa4(38) = 8*p1mup4*F(88)+p6mup4*F(89)+p4mup4*F(90)+p5mup4*F(91)
       Fa4(39) = p1mup4*F(98)+p6mup4*F(99)+p4mup4*F(100)+p5mup4*F(101)
       Fa4(40) = p6mup4*F(109)+p4mup4*F(110)+p5mup4*F(111)+4*p1mup4*F(
     &   112)
       Fa4(41) = p1mup4*F(113)+p6mup4*F(114)+p4mup4*F(115)+p5mup4*F(11
     &   6)
       Fa4(42) = 4*p1mup4*F(117)+p6mup4*F(118)+p4mup4*F(119)+p5mup4*F(
     &   120)
       Fa4(43) = p4mup4*F(127)+4*(p1mup4*F(125)+p6mup4*F(126)+p5mup4*F
     &   (128))
       Fa4(44) = -(p1mup4*F(129))+p6mup4*F(130)+p4mup4*F(131)+p5mup4*F
     &   (132)
       Fa4(45) = p4mup4*F(142)+2*(p1mup4*F(140)-p6mup4*F(141)+p5mup4*F
     &   (143))
       Fa4(46) = p1mup4*F(152)+p6mup4*F(153)+p4mup4*F(154)-p5mup4*F(15
     &   5)
       Fa4(47) = p6mup4*F(163)+p4mup4*F(164)+p5mup4*F(165)+p1mup4*F(16
     &   6)
       Fa4(48) = p6mup4*F(167)+p4mup4*F(168)+p5mup4*F(169)+p1mup4*F(17
     &   0)
       Return
       End
