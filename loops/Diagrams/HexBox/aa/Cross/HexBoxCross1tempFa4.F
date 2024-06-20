c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa4(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa4
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
       Complex*16 Fa4(40:52)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa4(40) = p4mup4*F(98)+8*p1mup4*F(99)+p5mup4*F(100)+p6mup4*F(10
     &   1)
       Fa4(41) = 4*p1mup4*F(103)+p4mup4*F(104)+p5mup4*F(105)+p6mup4*F(
     &   106)
       Fa4(42) = p4mup4*F(126)-p1mup4*F(127)+p5mup4*F(128)+p6mup4*F(12
     &   9)
       Fa4(43) = p4mup4*F(130)+4*(p1mup4*F(131)+p5mup4*F(132)+p6mup4*F
     &   (133))
       Fa4(44) = p4mup4*F(136)+4*p1mup4*F(137)+p5mup4*F(138)+p6mup4*F(
     &   139)
       Fa4(45) = p4mup4*F(140)+p1mup4*F(141)+p5mup4*F(142)+p6mup4*F(14
     &   3)
       Fa4(46) = p4mup4*F(146)-p6mup4*F(147)+p5mup4*F(148)+p1mup4*F(14
     &   9)
       Fa4(47) = p1mup4*F(150)+p4mup4*F(151)+p5mup4*F(152)+p6mup4*F(15
     &   3)
       Fa4(48) = p5mup4*F(154)+p6mup4*F(155)
       Fa4(49) = -(p1mup4*F(156))+p5mup4*F(157)+p4mup4*F(158)+p6mup4*F
     &   (159)
       Fa4(50) = p4mup4*F(162)+p5mup4*F(163)+p6mup4*F(164)+p1mup4*F(16
     &   5)
       Fa4(51) = p4mup4*F(166)+p5mup4*F(167)+p6mup4*F(168)+p1mup4*F(16
     &   9)
       Fa4(52) = p4mup4*F(174)+2*(p1mup4*F(175)+p5mup4*F(176)-p6mup4*F
     &   (177))
       Return
       End
