c************** Calling the Fa functions*************************
       subroutine HexBoxFa5(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa5)
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
       Complex*16 Fa5(49:60)
       COMMON/HexBoxFaFunctions/Fa
       Fa5(49) = p1mup4*F(171)-p6mup4*F(172)-p4mup4*F(173)-p5mup4*F(17
     &   4)
       Fa5(50) = p6mup4*F(177)+p4mup4*F(178)+p5mup4*F(179)+p1mup4*F(18
     &   0)
       Fa5(51) = -4*(p1mup4*F(181)-p6mup4*F(182)-p4mup4*F(183))+8*p5mu
     &   p4*F(184)
       Fa5(52) = p6mup4*F(185)+p4mup4*F(186)+p5mup4*F(187)-p1mup4*F(18
     &   8)
       Fa5(53) = p5mup4*F(191)-p1mup4*F(193)+p4mup4*F(194)
       Fa5(54) = p6mup4*F(195)+p4mup4*F(196)+p5mup4*F(197)+p1mup4*F(19
     &   8)
       Fa5(55) = p1mup4*F(199)+p6mup4*F(200)+p4mup4*F(201)+p5mup4*F(20
     &   2)
       Fa5(56) = p6mup4*F(204)+p4mup4*F(205)+p5mup4*F(206)-p1mup4*F(20
     &   7)
       Fa5(57) = p6mup4*F(208)+p4mup4*F(209)+p5mup4*F(210)+p1mup4*F(21
     &   1)
       Fa5(58) = p1mup4*F(212)-p5mup4*F(213)+p6mup4*F(214)-p4mup4*F(21
     &   5)
       Fa5(59) = p6mup4*F(217)+p4mup4*F(218)+p5mup4*F(219)+p1mup4*F(22
     &   0)
       Fa5(60) = p1mup4*F(221)+p6mup4*F(222)+p5mup4*F(223)+p4mup4*F(22
     &   4)
       Return
       End
